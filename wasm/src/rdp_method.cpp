#include "rdp_method.hpp"

#include "json.hpp"
#include "phylogeny.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iterator>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace rdp {
namespace {

constexpr double kMinimumProbability = 1e-300;
constexpr double kDistanceCorrelationCutoff = 0.05;
constexpr std::size_t kCorrelationFlankInformativeSites = 60;
constexpr std::size_t kEventTreeBootstrapReplicates = 10;
constexpr std::size_t kEventTreeSequenceCap = 100;
constexpr std::size_t kWorkingFragmentSequenceCap = 256;
constexpr std::size_t kFragmentReentryAlignmentLengthLimit = 100000;
constexpr std::uint64_t kNativeCorrectionCap =
    (std::uint64_t{255} * 255 * 255 * 255) / 2;

struct CorrelationRegionProfile {
  std::array<double, 3> category_fraction{};
  std::size_t comparable_sites = 0;
  std::size_t category_sites = 0;
};

using CorrelationRegions = std::array<std::vector<std::size_t>, 5>;
using CorrelationProfiles = std::array<CorrelationRegionProfile, 5>;
using PhylogeneticRegions = std::array<std::vector<std::size_t>, 6>;

struct CorrelationRegionLayout {
  CorrelationRegions regions;
  std::array<std::size_t, 4> boundaries{};
};

std::uint64_t choose_three(std::uint64_t count) {
  if (count < 3) return 0;
  return count * (count - 1) * (count - 2) / 6;
}

std::uint64_t choose_two(std::uint64_t count) {
  if (count < 2) return 0;
  return count * (count - 1) / 2;
}

std::array<std::uint8_t, 2> pair_members(std::uint8_t pair) {
  switch (pair) {
    case 0: return {0, 1};
    case 1: return {0, 2};
    default: return {1, 2};
  }
}

bool source_in_list(
    std::uint8_t outside_pair,
    std::uint8_t inside_pair,
    std::array<std::uint8_t, 3>& in_list) {
  // MakeINList maps the closest representative pair outside and inside the
  // tract into the three rows consumed by MakeACOR. It deliberately has no
  // mapping when the closest pair is unchanged.
  if (outside_pair == 0 && inside_pair == 1) in_list = {1, 0, 2};
  else if (outside_pair == 0 && inside_pair == 2) in_list = {0, 1, 2};
  else if (outside_pair == 1 && inside_pair == 0) in_list = {2, 0, 1};
  else if (outside_pair == 1 && inside_pair == 2) in_list = {0, 2, 1};
  else if (outside_pair == 2 && inside_pair == 0) in_list = {2, 1, 0};
  else if (outside_pair == 2 && inside_pair == 1) in_list = {1, 2, 0};
  else return false;
  return true;
}

std::uint64_t sequence_pair_key(std::uint32_t first, std::uint32_t second) {
  if (first > second) std::swap(first, second);
  return (static_cast<std::uint64_t>(first) << 32U) | second;
}

std::array<std::uint64_t, 3> triplet_pair_keys(
    const std::array<std::uint32_t, 3>& triplet) {
  return {
      sequence_pair_key(triplet[0], triplet[1]),
      sequence_pair_key(triplet[0], triplet[2]),
      sequence_pair_key(triplet[1], triplet[2]),
  };
}

std::array<std::uint32_t, 3> canonical_triplet(
    std::array<std::uint32_t, 3> triplet) {
  std::sort(triplet.begin(), triplet.end());
  return triplet;
}

std::uint64_t signal_signature(const Signal& signal) {
  // FNV-1a is used only to select a small collision bucket. Exact fields are
  // checked before a signal is considered a duplicate.
  std::uint64_t hash = 1469598103934665603ULL;
  const auto mix = [&](std::uint64_t value) {
    hash ^= value;
    hash *= 1099511628211ULL;
  };
  const auto triplet = canonical_triplet(signal.triplet);
  mix(triplet[0]);
  mix(triplet[1]);
  mix(triplet[2]);
  mix(signal.recombinant);
  mix(signal.beginning);
  mix(signal.ending);
  return hash;
}

template <typename T>
void sort_unique(std::vector<T>& values) {
  std::sort(values.begin(), values.end());
  values.erase(std::unique(values.begin(), values.end()), values.end());
}

double binomial_tail(std::size_t trials, std::size_t successes, double probability) {
  if (successes > trials) return 0.0;
  if (probability <= 0.0) return successes == 0 ? 1.0 : 0.0;
  if (probability >= 1.0) return 1.0;

  const double log_p = std::log(probability);
  const double log_q = std::log1p(-probability);
  double term = std::lgamma(static_cast<double>(trials + 1)) -
      std::lgamma(static_cast<double>(successes + 1)) -
      std::lgamma(static_cast<double>(trials - successes + 1)) +
      static_cast<double>(successes) * log_p +
      static_cast<double>(trials - successes) * log_q;
  double log_sum = -std::numeric_limits<double>::infinity();
  for (std::size_t count = successes; count <= trials; ++count) {
    const double high = std::max(log_sum, term);
    const double low = std::min(log_sum, term);
    log_sum = std::isinf(low) ? high : high + std::log1p(std::exp(low - high));
    if (count < trials) {
      term += std::log(static_cast<double>(trials - count)) -
          std::log(static_cast<double>(count + 1)) + log_p - log_q;
    }
  }
  return std::exp(log_sum);
}

double rdp_probability(
    std::size_t region_length,
    std::size_t matching_sites,
    double background_identity,
    std::size_t informative_length) {
  if (region_length == 0 || informative_length == 0) return 1.0;
  double length_factor = 1.0;
  if (region_length >= 170) {
    length_factor = static_cast<double>(region_length) / 169.0;
    const std::size_t different = region_length - matching_sites;
    const std::size_t scaled_different = static_cast<std::size_t>(
        std::llround(static_cast<double>(different) * 169.0 /
                     static_cast<double>(region_length)));
    region_length = 169;
    matching_sites = region_length - std::min(region_length, scaled_different);
  }

  double probability = binomial_tail(region_length, matching_sites, background_identity);
  probability *= static_cast<double>(informative_length) /
      static_cast<double>(std::max<std::size_t>(1, region_length));
  if (length_factor > 1.000001) {
    probability = probability > 0.0 ? std::pow(probability, length_factor) : 0.05;
  }
  if (probability < kMinimumProbability) probability = kMinimumProbability;
  if (!std::isfinite(probability)) probability = 1.0;
  return probability;
}

const char* review_state_name(ReviewState state) {
  switch (state) {
    case ReviewState::accepted: return "accepted";
    case ReviewState::rejected: return "rejected";
    default: return "unreviewed";
  }
}

std::string csv_cell(std::string_view value) {
  std::string result = "\"";
  for (const char character : value) {
    if (character == '"') result += '"';
    result += character;
  }
  result += '"';
  return result;
}

std::size_t wrap_coordinate(std::int64_t coordinate, std::size_t length) {
  const auto signed_length = static_cast<std::int64_t>(length);
  const auto wrapped = (coordinate - 1) % signed_length;
  return static_cast<std::size_t>((wrapped + signed_length) % signed_length + 1);
}

bool triplet_informative_at(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    std::size_t coordinate) {
  const std::size_t position = coordinate - 1;
  const std::uint8_t first = alignment.at(triplet[0], position);
  const std::uint8_t second = alignment.at(triplet[1], position);
  const std::uint8_t third = alignment.at(triplet[2], position);
  if (first == 0 || second == 0 || third == 0) return false;
  return (first == second && first != third) ||
      (first == third && first != second) ||
      (second == third && second != first);
}

std::size_t correlation_boundary(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    std::size_t start,
    std::size_t end_exclusive,
    int direction) {
  if (alignment.length == 0) return 0;
  start = wrap_coordinate(static_cast<std::int64_t>(start), alignment.length);
  end_exclusive = wrap_coordinate(static_cast<std::int64_t>(end_exclusive), alignment.length);
  std::size_t coordinate = start;
  std::size_t informative = 0;
  for (std::size_t visited = 0;
       coordinate != end_exclusive && visited < alignment.length;
       ++visited) {
    if (triplet_informative_at(alignment, triplet, coordinate)) {
      ++informative;
      if (informative >= kCorrelationFlankInformativeSites) return coordinate;
    }
    coordinate = wrap_coordinate(
        static_cast<std::int64_t>(coordinate) + direction,
        alignment.length);
  }
  return coordinate;
}

std::vector<std::size_t> forward_region(
    const Alignment& alignment,
    std::size_t start,
    std::size_t end_exclusive) {
  std::vector<std::size_t> positions;
  if (alignment.length == 0) return positions;
  start = wrap_coordinate(static_cast<std::int64_t>(start), alignment.length);
  end_exclusive = wrap_coordinate(static_cast<std::int64_t>(end_exclusive), alignment.length);
  std::size_t coordinate = start;
  for (std::size_t visited = 0;
       coordinate != end_exclusive && visited < alignment.length;
       ++visited) {
    positions.push_back(coordinate);
    coordinate = wrap_coordinate(
        static_cast<std::int64_t>(coordinate) + 1,
        alignment.length);
  }
  return positions;
}

CorrelationRegionLayout build_correlation_regions(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    std::size_t beginning,
    std::size_t ending) {
  CorrelationRegionLayout layout;
  if (alignment.length == 0) return layout;
  const std::size_t before_beginning = wrap_coordinate(
      static_cast<std::int64_t>(beginning) - 1,
      alignment.length);
  const std::size_t after_beginning = wrap_coordinate(
      static_cast<std::int64_t>(beginning),
      alignment.length);
  const std::size_t before_ending = wrap_coordinate(
      static_cast<std::int64_t>(ending),
      alignment.length);
  const std::size_t after_ending = wrap_coordinate(
      static_cast<std::int64_t>(ending) + 1,
      alignment.length);
  layout.boundaries[0] = correlation_boundary(
      alignment,
      triplet,
      before_beginning,
      after_ending,
      -1);
  layout.boundaries[1] = correlation_boundary(
      alignment,
      triplet,
      after_beginning,
      ending,
      1);
  layout.boundaries[2] = correlation_boundary(
      alignment,
      triplet,
      before_ending,
      beginning,
      -1);
  layout.boundaries[3] = correlation_boundary(
      alignment,
      triplet,
      after_ending,
      before_beginning,
      1);

  // Match MakeBPosLR + MakeSDMP2 exactly: reverse-walk boundaries become
  // inclusive starts, while forward-walk boundaries become exclusive ends.
  layout.regions[0] = forward_region(alignment, layout.boundaries[0], before_beginning);
  layout.regions[1] = forward_region(alignment, beginning, layout.boundaries[1]);
  layout.regions[2] = forward_region(alignment, layout.boundaries[2], ending);
  layout.regions[3] = forward_region(alignment, after_ending, layout.boundaries[3]);
  layout.regions[4] = forward_region(alignment, beginning, ending);
  return layout;
}

std::size_t non_gap_sites_inclusive(
    const Alignment& alignment,
    std::uint32_t sequence,
    std::size_t start,
    std::size_t end_inclusive) {
  if (alignment.length == 0) return 0;
  start = wrap_coordinate(static_cast<std::int64_t>(start), alignment.length);
  end_inclusive = wrap_coordinate(static_cast<std::int64_t>(end_inclusive), alignment.length);
  std::size_t count = 0;
  std::size_t coordinate = start;
  for (std::size_t visited = 0; visited < alignment.length; ++visited) {
    if (alignment.at(sequence, coordinate - 1) != 0) ++count;
    if (coordinate == end_inclusive) break;
    coordinate = wrap_coordinate(
        static_cast<std::int64_t>(coordinate) + 1,
        alignment.length);
  }
  return count;
}

std::array<std::size_t, 2> breakpoint_overlap_sites(
    const Alignment& alignment,
    std::uint32_t sequence,
    const CorrelationRegionLayout& layout) {
  return {
      non_gap_sites_inclusive(
          alignment, sequence, layout.boundaries[0], layout.boundaries[1]),
      non_gap_sites_inclusive(
          alignment, sequence, layout.boundaries[2], layout.boundaries[3]),
  };
}

CorrelationRegionProfile correlation_region_profile(
    const Alignment& alignment,
    std::uint32_t presumed_recombinant,
    std::uint32_t parent_one,
    std::uint32_t parent_two,
    std::uint32_t candidate,
    const std::vector<std::size_t>& positions) {
  CorrelationRegionProfile profile;
  std::array<std::size_t, 3> categories{};
  for (const std::size_t coordinate : positions) {
    const std::size_t position = coordinate - 1;
    const std::uint8_t anchor = alignment.at(presumed_recombinant, position);
    const std::uint8_t first = alignment.at(parent_one, position);
    const std::uint8_t second = alignment.at(parent_two, position);
    const std::uint8_t observed = alignment.at(candidate, position);
    if (anchor == 0 || first == 0 || second == 0 || observed == 0) continue;
    ++profile.comparable_sites;
    if (observed == first && observed != second) {
      ++categories[0];
      ++profile.category_sites;
    } else if (observed == second && observed != first) {
      ++categories[1];
      ++profile.category_sites;
    } else if (first == second && observed != first) {
      ++categories[2];
      ++profile.category_sites;
    }
  }
  if (profile.category_sites > 0) {
    for (std::size_t category = 0; category < 3; ++category) {
      profile.category_fraction[category] =
          static_cast<double>(categories[category]) /
          static_cast<double>(profile.category_sites);
    }
  }
  return profile;
}

CorrelationProfiles correlation_profiles(
    const Alignment& alignment,
    std::uint32_t presumed_recombinant,
    std::uint32_t parent_one,
    std::uint32_t parent_two,
    std::uint32_t candidate,
    const CorrelationRegions& regions) {
  CorrelationProfiles profiles;
  for (std::size_t region = 0; region < regions.size(); ++region) {
    profiles[region] = correlation_region_profile(
        alignment,
        presumed_recombinant,
        parent_one,
        parent_two,
        candidate,
        regions[region]);
  }
  return profiles;
}

std::array<double, 6> correlation_vector(
    const CorrelationProfiles& profiles,
    std::size_t matrix_pair) {
  std::array<double, 6> values{};
  if (matrix_pair == 0) {
    std::copy(profiles[0].category_fraction.begin(), profiles[0].category_fraction.end(), values.begin());
    std::copy(profiles[1].category_fraction.begin(), profiles[1].category_fraction.end(), values.begin() + 3);
  } else if (matrix_pair == 1) {
    std::copy(profiles[3].category_fraction.begin(), profiles[3].category_fraction.end(), values.begin());
    std::copy(profiles[2].category_fraction.begin(), profiles[2].category_fraction.end(), values.begin() + 3);
  } else {
    for (std::size_t category = 0; category < 3; ++category) {
      values[category] =
          (profiles[0].category_fraction[category] +
           profiles[3].category_fraction[category]) /
          2.0;
      values[category + 3] = profiles[4].category_fraction[category];
    }
  }
  return values;
}

std::size_t correlation_minimum_comparable_sites(
    const CorrelationProfiles& profiles,
    std::size_t matrix_pair) {
  if (matrix_pair == 0) {
    return std::min(profiles[0].comparable_sites, profiles[1].comparable_sites);
  }
  if (matrix_pair == 1) {
    return std::min(profiles[3].comparable_sites, profiles[2].comparable_sites);
  }
  return std::min({
      profiles[0].comparable_sites,
      profiles[3].comparable_sites,
      profiles[4].comparable_sites,
  });
}

bool correlation_pair_has_categories(
    const CorrelationProfiles& profiles,
    std::size_t matrix_pair) {
  if (matrix_pair == 0) {
    return profiles[0].category_sites > 0 && profiles[1].category_sites > 0;
  }
  if (matrix_pair == 1) {
    return profiles[3].category_sites > 0 && profiles[2].category_sites > 0;
  }
  return profiles[0].category_sites > 0 && profiles[3].category_sites > 0 &&
      profiles[4].category_sites > 0;
}

bool pearson_six(
    const std::array<double, 6>& first,
    const std::array<double, 6>& second,
    double& correlation) {
  const double first_mean =
      std::accumulate(first.begin(), first.end(), 0.0) / static_cast<double>(first.size());
  const double second_mean =
      std::accumulate(second.begin(), second.end(), 0.0) / static_cast<double>(second.size());
  double covariance = 0.0;
  double first_variance = 0.0;
  double second_variance = 0.0;
  for (std::size_t index = 0; index < first.size(); ++index) {
    const double first_delta = first[index] - first_mean;
    const double second_delta = second[index] - second_mean;
    covariance += first_delta * second_delta;
    first_variance += first_delta * first_delta;
    second_variance += second_delta * second_delta;
  }
  if (first_variance <= 1e-18 || second_variance <= 1e-18) {
    // CalCR assigns one when either six-value vector has zero variance.
    correlation = 1.0;
    return true;
  }
  correlation = covariance / std::sqrt(first_variance * second_variance);
  correlation = std::clamp(correlation, -1.0, 1.0);
  return std::isfinite(correlation);
}

double pearson_six_two_sided_p(double correlation) {
  const double absolute = std::clamp(std::abs(correlation), 0.0, 1.0);
  // With six paired values the Pearson t test has four degrees of freedom.
  // I_(1-r^2)(2, 1/2) simplifies to this exact polynomial.
  return std::clamp(1.0 - 1.5 * absolute + 0.5 * absolute * absolute * absolute, 0.0, 1.0);
}

bool best_category_correlation(
    const std::array<double, 6>& anchor,
    const std::array<double, 6>& candidate,
    double& direct,
    double& selected,
    std::uint8_t& inversion_code) {
  if (!pearson_six(anchor, candidate, direct)) return false;
  selected = direct;
  inversion_code = 0;
  if (direct >= 0.83) return true;

  constexpr std::array<std::array<std::size_t, 3>, 5> permutations{{
      {{1, 0, 2}},
      {{2, 1, 0}},
      {{0, 2, 1}},
      {{1, 2, 0}},
      {{2, 0, 1}},
  }};
  for (std::size_t permutation = 0; permutation < permutations.size(); ++permutation) {
    std::array<double, 6> relabelled{};
    for (std::size_t half = 0; half < 2; ++half) {
      for (std::size_t category = 0; category < 3; ++category) {
        relabelled[half * 3 + category] =
            candidate[half * 3 + permutations[permutation][category]];
      }
    }
    double correlation = 0.0;
    if (!pearson_six(anchor, relabelled, correlation) || correlation <= selected) continue;
    selected = correlation;
    // The supplied CalCR routine stores the two cyclic permutations together
    // as inversion class four.
    inversion_code = static_cast<std::uint8_t>(
        permutation < 3 ? permutation + 1 : 4);
  }
  return true;
}

std::int8_t unique_dominant_category(const CorrelationRegionProfile& profile) {
  std::size_t best = 0;
  for (std::size_t category = 1; category < profile.category_fraction.size(); ++category) {
    if (profile.category_fraction[category] > profile.category_fraction[best]) best = category;
  }
  for (std::size_t category = 0; category < profile.category_fraction.size(); ++category) {
    if (category != best &&
        std::abs(profile.category_fraction[category] - profile.category_fraction[best]) <= 1e-12) {
      return -1;
    }
  }
  return static_cast<std::int8_t>(best);
}

std::array<std::uint8_t, 3> correlation_warnings(
    const CorrelationProfiles& profiles) {
  std::array<std::uint8_t, 3> warnings{};
  const auto same_dominant = [&](std::size_t first, std::size_t second) {
    const std::int8_t first_category = unique_dominant_category(profiles[first]);
    const std::int8_t second_category = unique_dominant_category(profiles[second]);
    return first_category >= 0 && first_category == second_category;
  };
  warnings[0] = same_dominant(0, 1) ? 1 : 0;
  warnings[1] = same_dominant(3, 2) ? 1 : 0;
  // The native caller clears the tract/outside warning when both breakpoint
  // pairs warn, and sets it only when exactly one pair warns.
  warnings[2] = warnings[0] != warnings[1] ? 1 : 0;
  return warnings;
}

PhylogeneticRegions build_phylogenetic_regions(
    const Alignment& alignment,
    const CorrelationRegionLayout& layout) {
  PhylogeneticRegions regions;
  regions[0] = layout.regions[0];
  regions[1] = layout.regions[1];
  regions[2] = layout.regions[3];
  regions[3] = layout.regions[2];
  regions[5] = layout.regions[4];

  std::vector<std::uint8_t> inside(alignment.length + 1, 0);
  for (const std::size_t coordinate : regions[5]) {
    if (coordinate <= alignment.length) inside[coordinate] = 1;
  }
  regions[4].reserve(alignment.length - std::min(alignment.length, regions[5].size()));
  for (std::size_t coordinate = 1; coordinate <= alignment.length; ++coordinate) {
    if (inside[coordinate] == 0) regions[4].push_back(coordinate);
  }
  return regions;
}

std::vector<std::uint32_t> select_tree_sequences(
    const Alignment& alignment,
    const std::vector<std::uint32_t>& active_sequences,
    const std::array<std::uint32_t, 3>& representatives) {
  std::vector<std::uint32_t> selected;
  selected.reserve(std::min(kEventTreeSequenceCap, active_sequences.size()));
  for (const std::uint32_t representative : representatives) {
    if (std::find(active_sequences.begin(), active_sequences.end(), representative) !=
            active_sequences.end() &&
        std::find(selected.begin(), selected.end(), representative) == selected.end()) {
      selected.push_back(representative);
    }
  }

  std::vector<std::pair<double, std::uint32_t>> ranked;
  ranked.reserve(active_sequences.size());
  for (const std::uint32_t sequence : active_sequences) {
    if (std::find(selected.begin(), selected.end(), sequence) != selected.end()) continue;
    double affinity = -1.0;
    for (const std::uint32_t representative : representatives) {
      std::size_t comparable = 0;
      std::size_t identical = 0;
      for (std::size_t position = 0; position < alignment.length; ++position) {
        const std::uint8_t first = alignment.at(sequence, position);
        const std::uint8_t second = alignment.at(representative, position);
        if (first == 0 || second == 0) continue;
        ++comparable;
        if (first == second) ++identical;
      }
      const double similarity = comparable == 0
          ? 0.0
          : static_cast<double>(identical) / static_cast<double>(comparable);
      affinity = std::max(affinity, similarity);
    }
    ranked.emplace_back(affinity, sequence);
  }
  std::stable_sort(ranked.begin(), ranked.end(), [](const auto& left, const auto& right) {
    if (left.first != right.first) return left.first > right.first;
    return left.second < right.second;
  });
  for (const auto& [affinity, sequence] : ranked) {
    (void)affinity;
    if (selected.size() >= kEventTreeSequenceCap) break;
    selected.push_back(sequence);
  }
  std::sort(selected.begin(), selected.end());
  return selected;
}

// MakePhPrScore in the supplied DLL source accepts two observations and uses
// 1 for a constant profile. Keeping that behaviour is important: the native
// consensus explicitly treats an exact +/-1 as non-informative.
bool source_pearson(
    const std::vector<std::pair<double, double>>& values,
    double& correlation) {
  if (values.size() < 2) return false;
  const double count = static_cast<double>(values.size());
  double sum_first = 0.0;
  double sum_second = 0.0;
  double sum_product = 0.0;
  double sum_first_squared = 0.0;
  double sum_second_squared = 0.0;
  for (const auto& [first, second] : values) {
    sum_first += first;
    sum_second += second;
    sum_product += first * second;
    sum_first_squared += first * first;
    sum_second_squared += second * second;
  }
  const double first_variance = count * sum_first_squared - sum_first * sum_first;
  const double second_variance = count * sum_second_squared - sum_second * sum_second;
  if (first_variance <= 1e-12 || second_variance <= 1e-12) {
    correlation = 1.0;
    return true;
  }
  correlation = (count * sum_product - sum_first * sum_second) /
      std::sqrt(first_variance * second_variance);
  correlation = std::clamp(correlation, -1.0, 1.0);
  return std::isfinite(correlation);
}

struct SourcePhylproScores {
  std::array<double, 3> primary{};
  std::array<double, 3> leave_one_out{};
  std::array<double, 3> displacement{};
  std::array<bool, 3> primary_valid{};
  std::array<bool, 3> leave_one_out_valid{};
  std::array<bool, 3> displacement_valid{};
};

template <typename Distance>
SourcePhylproScores source_phylpro_scores(
    const std::array<std::uint32_t, 3>& representatives,
    const std::vector<std::uint32_t>& involved,
    Distance distance) {
  SourcePhylproScores result;
  for (std::size_t role = 0; role < representatives.size(); ++role) {
    std::vector<std::pair<double, double>> values;
    values.reserve(involved.size());
    double displacement = 0.0;
    for (const std::uint32_t sequence : involved) {
      if (sequence == representatives[role]) continue;
      const auto [outside, inside] = distance(representatives[role], sequence);
      if (!std::isfinite(outside) || !std::isfinite(inside)) continue;
      values.emplace_back(outside, inside);
      displacement += std::abs(outside - inside);
    }
    result.primary_valid[role] = source_pearson(values, result.primary[role]);
    result.displacement[role] = displacement;
    result.displacement_valid[role] = values.size() >= 2;
  }

  // The X=1..3 passes in MakePhPrScore remove one representative, correlate
  // the other two anchors, then average those two correlations. The resulting
  // score belongs to the omitted representative.
  for (std::size_t omitted = 0; omitted < representatives.size(); ++omitted) {
    double total = 0.0;
    std::size_t correlations = 0;
    for (std::size_t role = 0; role < representatives.size(); ++role) {
      if (role == omitted) continue;
      std::vector<std::pair<double, double>> values;
      values.reserve(involved.size());
      for (const std::uint32_t sequence : involved) {
        if (sequence == representatives[omitted] ||
            sequence == representatives[role]) {
          continue;
        }
        const auto [outside, inside] = distance(representatives[role], sequence);
        if (!std::isfinite(outside) || !std::isfinite(inside)) continue;
        values.emplace_back(outside, inside);
      }
      double correlation = 0.0;
      if (source_pearson(values, correlation)) {
        total += correlation;
        ++correlations;
      }
    }
    if (correlations == 2) {
      result.leave_one_out[omitted] = total / 2.0;
      result.leave_one_out_valid[omitted] = true;
    }
  }
  return result;
}

template <typename JcDistance, typename TreeDistance>
std::vector<std::uint32_t> source_involved_sequences(
    const std::array<std::uint32_t, 3>& representatives,
    const std::vector<std::uint32_t>& panel,
    JcDistance jc_distance,
    TreeDistance tree_distance) {
  std::vector<std::uint32_t> involved(
      representatives.begin(), representatives.end());
  for (const std::uint32_t candidate : panel) {
    if (std::find(representatives.begin(), representatives.end(), candidate) !=
        representatives.end()) {
      continue;
    }
    bool invalid = false;
    bool outside_too_close = false;
    bool inside_too_close = false;
    for (std::size_t role = 0; role < representatives.size(); ++role) {
      const std::uint32_t anchor = representatives[role];
      const std::uint32_t other_one = representatives[(role + 1) % 3];
      const std::uint32_t other_two = representatives[(role + 2) % 3];
      const auto [jc_outside, jc_inside] = jc_distance(anchor, candidate);
      const auto [tree_outside, tree_inside] = tree_distance(anchor, candidate);
      if (!std::isfinite(jc_outside) || !std::isfinite(jc_inside) ||
          jc_outside >= 10.0 || jc_inside >= 10.0 ||
          !std::isfinite(tree_outside) || !std::isfinite(tree_inside)) {
        invalid = true;
        break;
      }
      const auto [outside_one, inside_one] = tree_distance(anchor, other_one);
      const auto [outside_two, inside_two] = tree_distance(anchor, other_two);
      outside_too_close = outside_too_close ||
          tree_outside < std::min(outside_one, outside_two);
      inside_too_close = inside_too_close ||
          tree_inside < std::min(inside_one, inside_two);
    }
    // MakeDoneThis keeps a candidate when it lies inside the representative
    // subtree in either the outside or inside matrix.
    if (!invalid && (!outside_too_close || !inside_too_close)) {
      involved.push_back(candidate);
    }
  }
  sort_unique(involved);
  return involved;
}

template <typename Distance>
double source_triplet_score(
    std::size_t role,
    const std::array<std::uint32_t, 3>& representatives,
    const std::vector<std::uint32_t>& panel,
    Distance distance) {
  if (panel.size() < 2) return 0.0;
  const std::uint32_t anchor = representatives[role];
  const auto [parent_one_outside, ignored_one] = distance(
      anchor, representatives[(role + 1) % 3]);
  const auto [parent_two_outside, ignored_two] = distance(
      anchor, representatives[(role + 2) % 3]);
  (void)ignored_one;
  (void)ignored_two;
  const double close_threshold = std::min(parent_one_outside, parent_two_outside);

  std::vector<std::size_t> group(panel.size(), 0);
  std::vector<double> outside(panel.size(), 0.0);
  std::vector<double> inside(panel.size(), 0.0);
  std::vector<std::size_t> unresolved;
  unresolved.reserve(panel.size());
  for (std::size_t index = 0; index < panel.size(); ++index) {
    const auto distances = distance(anchor, panel[index]);
    outside[index] = distances.first;
    inside[index] = distances.second;
    if (!(outside[index] < close_threshold)) unresolved.push_back(index);
  }
  std::stable_sort(unresolved.begin(), unresolved.end(), [&](std::size_t left, std::size_t right) {
    if (std::abs(outside[left] - outside[right]) > 1e-9) {
      return outside[left] < outside[right];
    }
    return panel[left] < panel[right];
  });
  std::size_t group_number = 0;
  double prior_distance = -std::numeric_limits<double>::infinity();
  for (const std::size_t index : unresolved) {
    if (group_number == 0 || std::abs(outside[index] - prior_distance) > 1e-9) {
      ++group_number;
      prior_distance = outside[index];
    }
    group[index] = group_number;
  }
  std::vector<std::size_t> group_size(group_number + 1, 0);
  for (const std::size_t value : group) ++group_size[value];

  const auto compare = [](double first, double second) {
    if (std::abs(first - second) <= 1e-9) return 0;
    return first < second ? -1 : 1;
  };
  double score = 0.0;
  for (std::size_t first = 0; first + 1 < panel.size(); ++first) {
    for (std::size_t second = first + 1; second < panel.size(); ++second) {
      if (!std::isfinite(outside[first]) || !std::isfinite(outside[second]) ||
          !std::isfinite(inside[first]) || !std::isfinite(inside[second])) {
        continue;
      }
      const int outside_order = compare(outside[second], outside[first]);
      const int inside_order = compare(inside[second], inside[first]);
      if (outside_order == inside_order) continue;
      const std::size_t denominator =
          group_size[group[first]] * group_size[group[second]];
      if (denominator > 0) score += 1.0 / static_cast<double>(denominator);
    }
  }
  return score;
}

double affinity_margin(
    double anchor_candidate,
    double anchor_parent_one,
    double anchor_parent_two,
    double candidate_parent_one,
    double candidate_parent_two) {
  if (!std::isfinite(anchor_candidate) || anchor_candidate >= 10.0) return -10.0;
  const double competing = std::min({
      anchor_parent_one,
      anchor_parent_two,
      candidate_parent_one,
      candidate_parent_two,
  });
  return competing - anchor_candidate;
}

const char* role_metric_name(RoleMetricKind kind) {
  switch (kind) {
    case RoleMetricKind::phpr: return "PhPr";
    case RoleMetricKind::tree_phpr: return "TreePhPr";
    case RoleMetricKind::collapsed_tree_phpr: return "CollapsedTreePhPr";
    case RoleMetricKind::sub_phpr: return "SubPhPr";
    case RoleMetricKind::tree_sub_phpr: return "TreeSubPhPr";
    case RoleMetricKind::subdist: return "SubDist";
    case RoleMetricKind::tree_subdist: return "TreeSubDist";
    case RoleMetricKind::triplet_score: return "TrpScore";
    case RoleMetricKind::three_set_support: return "ThreeSetSupport";
  }
  return "Unknown";
}

std::string fasta_name(std::string_view value) {
  std::string name;
  name.reserve(value.size());
  for (const char character : value) {
    if (character == '\n' || character == '\r') name.push_back('_');
    else name.push_back(character);
  }
  return name.empty() ? "sequence" : name;
}

void write_fasta_record(
    std::ostringstream& out,
    std::string_view name,
    std::string_view sequence) {
  out << '>' << fasta_name(name) << '\n';
  constexpr std::size_t width = 80;
  for (std::size_t offset = 0; offset < sequence.size(); offset += width) {
    out << sequence.substr(offset, std::min(width, sequence.size() - offset)) << '\n';
  }
}

bool coordinate_in_tract(
    std::size_t coordinate,
    std::size_t beginning,
    std::size_t ending,
    bool wraps_origin) {
  return wraps_origin
      ? coordinate > beginning || coordinate < ending
      : coordinate > beginning && coordinate < ending;
}

}  // namespace

RdpScanner::RdpScanner(const Alignment& alignment)
    : alignment_(alignment) {}

void RdpScanner::reset_working_alignment() {
  // The mutable cyclic alignment needs sequence state, not the original
  // O(N^2) identity matrix or parser diagnostics. Avoid copying those every
  // time a correction/rejection rebuilds the event chain.
  working_alignment_ = {};
  working_alignment_.names = alignment_.names;
  working_alignment_.sequences = alignment_.sequences;
  working_alignment_.states = alignment_.states;
  working_alignment_.sequence_summaries = alignment_.sequence_summaries;
  working_alignment_.length = alignment_.length;
  working_origins_.resize(alignment_.sequence_count());
  std::iota(working_origins_.begin(), working_origins_.end(), 0);
  working_fragment_events_.assign(alignment_.sequence_count(), -1);
  fragment_reentry_capped_ = false;
}

void RdpScanner::rebuild_working_before_event(std::size_t event_index) {
  reset_working_alignment();
  const std::size_t limit = std::min(event_index, events_.size());
  for (std::size_t index = 0; index < limit; ++index) {
    if (events_[index].tract_erased_for_detection) {
      (void)erase_event_tract(events_[index]);
    }
  }
}

std::uint64_t RdpScanner::valid_triplet_count() const {
  const std::uint64_t count = active_sequences_.size();
  std::uint64_t total = choose_three(count);
  if (total == 0) return 0;
  std::vector<std::uint64_t> origin_counts(alignment_.sequence_count(), 0);
  for (const std::uint32_t sequence : active_sequences_) {
    if (sequence < working_origins_.size() &&
        working_origins_[sequence] < origin_counts.size()) {
      ++origin_counts[working_origins_[sequence]];
    }
  }
  for (const std::uint64_t origin_count : origin_counts) {
    if (origin_count < 2) continue;
    const std::uint64_t invalid = choose_two(origin_count) * (count - origin_count) +
        choose_three(origin_count);
    total = invalid > total ? 0 : total - invalid;
  }
  return total;
}

void RdpScanner::refresh_active_sequences() {
  active_sequences_.clear();
  const std::size_t fragment_minimum = std::max<std::size_t>(5, options_.window_sites);
  for (std::size_t sequence = 0; sequence < working_alignment_.sequence_count(); ++sequence) {
    if (sequence >= working_origins_.size()) continue;
    const std::uint32_t origin = working_origins_[sequence];
    if (origin >= options_.mask.size() || options_.mask[origin] != 0) continue;
    std::size_t valid_sites = 0;
    const std::size_t offset = sequence * working_alignment_.length;
    for (std::size_t position = 0; position < working_alignment_.length; ++position) {
      if (working_alignment_.states[offset + position] != 0) ++valid_sites;
    }
    if (valid_sites < fragment_minimum) continue;
    active_sequences_.push_back(static_cast<std::uint32_t>(sequence));
  }
  total_triplets_ = valid_triplet_count();
  correction_tests_ = std::min(total_triplets_, kNativeCorrectionCap);
}

bool RdpScanner::working_triplet_is_valid(
    const std::array<std::uint32_t, 3>& triplet) const {
  if (triplet[0] >= working_origins_.size() || triplet[1] >= working_origins_.size() ||
      triplet[2] >= working_origins_.size()) {
    return false;
  }
  const std::uint32_t first = working_origins_[triplet[0]];
  const std::uint32_t second = working_origins_[triplet[1]];
  const std::uint32_t third = working_origins_[triplet[2]];
  return first != second && first != third && second != third;
}

void RdpScanner::map_signal_to_original(Signal& signal) const {
  for (std::size_t member = 0; member < signal.triplet.size(); ++member) {
    const std::uint32_t working = signal.triplet[member];
    if (working >= working_origins_.size()) continue;
    signal.fragment_event_context[member] = working_fragment_events_[working];
    if (signal.fragment_event_context[member] >= 0) signal.fragment_assisted = true;
    signal.triplet[member] = working_origins_[working];
  }
  if (signal.recombinant < working_origins_.size()) {
    signal.recombinant = working_origins_[signal.recombinant];
  }
  if (signal.major_parent < working_origins_.size()) {
    signal.major_parent = working_origins_[signal.major_parent];
  }
  if (signal.minor_parent < working_origins_.size()) {
    signal.minor_parent = working_origins_[signal.minor_parent];
  }
}

void RdpScanner::reset_round_cursor() {
  cursor_a_ = 0;
  cursor_b_ = 1;
  cursor_c_ = 2;
  processed_triplets_ = 0;
  profile_scratch_.category.clear();
  profile_scratch_.coordinates.clear();
  for (auto& counts : profile_scratch_.rolling_counts) counts.clear();
  round_signal_index_.clear();
}

bool RdpScanner::begin(ScanOptions options, std::string& error) {
  if (alignment_.sequence_count() < 3 || alignment_.length == 0) {
    error = "A valid alignment must be loaded before scanning.";
    return false;
  }
  if (!(options.p_value_cutoff > 0.0 && options.p_value_cutoff <= 1.0)) {
    error = "The highest acceptable p-value must be greater than zero and no greater than one.";
    return false;
  }
  if (options.window_sites < 5) {
    error = "The RDP window must contain at least five variable sites.";
    return false;
  }
  if (options.mask.size() != alignment_.sequence_count()) {
    options.mask.assign(alignment_.sequence_count(), 0);
  }

  options_ = std::move(options);
  reset_working_alignment();
  refresh_active_sequences();
  if (active_sequences_.size() < 3) {
    error = "At least three unmasked sequences are required for an exploratory RDP scan.";
    return false;
  }

  signals_.clear();
  events_.clear();
  reset_round_cursor();
  cumulative_triplets_ = 0;
  cancelled_.store(false);
  running_ = true;
  primary_done_ = false;
  done_ = false;
  scan_round_ = 1;
  round_signal_begin_ = 0;
  fixed_event_count_ = 0;
  cycle_termination_ = "scanning";
  reconciliation_required_after_ = -1;
  return true;
}

int RdpScanner::scan_batch(std::size_t triplet_budget, std::string& error) {
  if (!running_) {
    if (done_) return 1;
    error = "No RDP scan is active.";
    return -1;
  }
  if (active_sequences_.size() < 3 || total_triplets_ == 0) {
    running_ = false;
    primary_done_ = true;
    cycle_termination_ = "fewer-than-three-active-origins";
    return 3;
  }
  const std::size_t budget = std::max<std::size_t>(1, triplet_budget);
  for (std::size_t index = 0; index < budget && !done_; ++index) {
    if (cancelled_.load()) {
      running_ = false;
      return 2;
    }
    const std::array<std::uint32_t, 3> triplet{
        active_sequences_[cursor_a_],
        active_sequences_[cursor_b_],
        active_sequences_[cursor_c_],
    };
    if (working_triplet_is_valid(triplet)) {
      scan_triplet(triplet);
      ++processed_triplets_;
      ++cumulative_triplets_;
    }
    if (!advance_triplet()) {
      if (finish_detection_round(error)) {
        ++scan_round_;
        round_signal_begin_ = signals_.size();
        reset_round_cursor();
        return 4;
      }
      if (!error.empty()) {
        running_ = false;
        return -1;
      }
      running_ = false;
      primary_done_ = true;
      return 3;
    }
  }
  return done_ ? 1 : primary_done_ ? 3 : 0;
}

bool RdpScanner::reconcile(std::string& error) {
  if (!primary_done_) {
    error = "The primary triplet screen must finish before event reconciliation.";
    return false;
  }
  if (cancelled_.load()) {
    error = "The scan was cancelled before event reconciliation.";
    return false;
  }
  reconciliation_required_after_ = -1;
  cycle_termination_ = cycle_termination_ == "scanning"
      ? "no-significant-signals"
      : cycle_termination_;
  done_ = true;
  return true;
}

void RdpScanner::cancel() {
  cancelled_.store(true);
}

bool RdpScanner::advance_triplet() {
  const std::size_t count = active_sequences_.size();
  if (++cursor_c_ < count) return true;
  ++cursor_b_;
  if (cursor_b_ + 1 < count) {
    cursor_c_ = cursor_b_ + 1;
    return true;
  }
  ++cursor_a_;
  if (cursor_a_ + 2 < count) {
    cursor_b_ = cursor_a_ + 1;
    cursor_c_ = cursor_b_ + 1;
    return true;
  }
  return false;
}

bool RdpScanner::build_profile(
    const std::array<std::uint32_t, 3>& triplet,
    TripletProfile& profile) const {
  return build_profile_on(working_alignment_, triplet, profile);
}

bool RdpScanner::build_profile_on(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    TripletProfile& profile) const {
  profile.category.clear();
  profile.coordinates.clear();
  profile.category_counts.fill(0);
  profile.similarities.fill(0.0);
  profile.sequences = triplet;
  if (profile.category.capacity() < alignment.length / 4) {
    profile.category.reserve(alignment.length / 4);
    profile.coordinates.reserve(alignment.length / 4);
  }

  std::array<std::size_t, 3> comparable{};
  std::array<std::size_t, 3> identical{};
  for (std::size_t position = 0; position < alignment.length; ++position) {
    const std::uint8_t first = alignment.at(triplet[0], position);
    const std::uint8_t second = alignment.at(triplet[1], position);
    const std::uint8_t third = alignment.at(triplet[2], position);
    if (first != 0 && second != 0) {
      ++comparable[0];
      if (first == second) ++identical[0];
    }
    if (first != 0 && third != 0) {
      ++comparable[1];
      if (first == third) ++identical[1];
    }
    if (second != 0 && third != 0) {
      ++comparable[2];
      if (second == third) ++identical[2];
    }
    if (first == 0 || second == 0 || third == 0) continue;
    std::uint8_t category = 255;
    if (first == second && first != third) category = 0;
    else if (first == third && first != second) category = 1;
    else if (second == third && second != first) category = 2;
    if (category == 255) continue;
    profile.category.push_back(category);
    profile.coordinates.push_back(position + 1);
    ++profile.category_counts[category];
  }

  const std::size_t half_window = options_.window_sites / 2;
  if (profile.category.size() < half_window * 2) return false;
  for (const std::size_t count : profile.category_counts) {
    if (count < std::max<std::size_t>(1, half_window / 3)) return false;
  }

  for (std::size_t pair = 0; pair < profile.similarities.size(); ++pair) {
    profile.similarities[pair] = comparable[pair] == 0
        ? 0.0
        : static_cast<double>(identical[pair]) / static_cast<double>(comparable[pair]);
  }
  compute_rolling_counts(profile);

  std::uint8_t overall_high_pair = 0;
  if (profile.similarities[0] >= profile.similarities[1] &&
      profile.similarities[0] >= profile.similarities[2]) {
    overall_high_pair = 0;
  } else if (profile.similarities[1] >= profile.similarities[0] &&
             profile.similarities[1] >= profile.similarities[2]) {
    overall_high_pair = 1;
  } else {
    overall_high_pair = 2;
  }
  const bool has_low_window = std::any_of(
      profile.rolling_counts[overall_high_pair].begin(),
      profile.rolling_counts[overall_high_pair].end(),
      [half_window](std::uint32_t count) { return count <= half_window; });
  if (!has_low_window) return false;
  return true;
}

void RdpScanner::compute_rolling_counts(TripletProfile& profile) const {
  const std::size_t length = profile.category.size();
  const std::size_t half = options_.window_sites / 2;
  for (auto& counts : profile.rolling_counts) counts.assign(length, 0);
  if (length == 0) return;

  std::array<std::uint32_t, 3> current{};
  for (std::int64_t offset = -static_cast<std::int64_t>(half);
       offset <= static_cast<std::int64_t>(half);
       ++offset) {
    const auto wrapped = static_cast<std::size_t>(
        (offset % static_cast<std::int64_t>(length) + static_cast<std::int64_t>(length)) %
        static_cast<std::int64_t>(length));
    ++current[profile.category[wrapped]];
  }
  for (std::size_t pair = 0; pair < 3; ++pair) profile.rolling_counts[pair][0] = current[pair];
  for (std::size_t position = 1; position < length; ++position) {
    const std::size_t remove = (position + length - 1 - (half % length)) % length;
    const std::size_t add = (position + half) % length;
    --current[profile.category[remove]];
    ++current[profile.category[add]];
    for (std::size_t pair = 0; pair < 3; ++pair) {
      profile.rolling_counts[pair][position] = current[pair];
    }
  }
}

std::array<std::uint8_t, 3> RdpScanner::ranked_pairs(const TripletProfile& profile) const {
  std::array<std::uint8_t, 3> order{0, 1, 2};
  const double length = static_cast<double>(profile.category.size());
  std::array<double, 3> average{
      static_cast<double>(profile.category_counts[0]) / length,
      static_cast<double>(profile.category_counts[1]) / length,
      static_cast<double>(profile.category_counts[2]) / length,
  };
  std::stable_sort(order.begin(), order.end(), [&](std::uint8_t left, std::uint8_t right) {
    if (std::abs(average[left] - average[right]) > 1e-12) return average[left] > average[right];
    return profile.similarities[left] > profile.similarities[right];
  });
  return order;
}

std::vector<Signal> RdpScanner::candidate_signals(
    const TripletProfile& profile,
    std::uint8_t high_pair,
    std::uint8_t candidate_pair,
    std::uint8_t low_pair,
    bool enforce_cutoff) const {
  std::vector<Signal> result;
  const std::size_t length = profile.category.size();
  std::size_t search = 0;

  const auto candidate_members = pair_members(candidate_pair);
  const auto high_members = pair_members(high_pair);
  std::uint8_t recombinant_local = 255;
  for (const auto member : candidate_members) {
    if (member == high_members[0] || member == high_members[1]) recombinant_local = member;
  }
  if (recombinant_local == 255) return result;
  const std::uint8_t major_local = high_members[0] == recombinant_local
      ? high_members[1]
      : high_members[0];
  const std::uint8_t minor_local = candidate_members[0] == recombinant_local
      ? candidate_members[1]
      : candidate_members[0];

  while (search < length) {
    while (search < length &&
           !(profile.rolling_counts[candidate_pair][search] >
                 profile.rolling_counts[high_pair][search] &&
             profile.rolling_counts[candidate_pair][search] >
                 profile.rolling_counts[low_pair][search])) {
      ++search;
    }
    if (search >= length) break;

    std::size_t beginning = search;
    bool beginning_wrapped = false;
    if (profile.category[beginning] == candidate_pair) {
      while (beginning > 0 && profile.category[beginning - 1] == candidate_pair) --beginning;
    } else {
      std::size_t examined = 0;
      while (profile.category[beginning] != candidate_pair && examined < length) {
        ++beginning;
        ++examined;
        if (beginning == length) {
          if (!options_.circular) break;
          beginning = 0;
          beginning_wrapped = true;
        }
      }
      if (examined == length || (!options_.circular && beginning == length)) break;
    }

    std::size_t ending = beginning;
    std::size_t matching = 0;
    std::size_t region_length = 0;
    bool wrapped = beginning_wrapped;
    for (std::size_t visited = 0; visited < length; ++visited) {
      if (profile.category[ending] == candidate_pair) ++matching;
      ++region_length;
      const std::size_t next = (ending + 1) % length;
      if (next == 0) {
        if (!options_.circular) break;
        wrapped = true;
      }
      const bool candidate_drops =
          profile.rolling_counts[candidate_pair][next] < profile.rolling_counts[high_pair][next] ||
          profile.rolling_counts[candidate_pair][next] < profile.rolling_counts[low_pair][next];
      if (candidate_drops && profile.category[next] != candidate_pair) break;
      ending = next;
      if (ending == beginning) break;
    }

    while (region_length > 1 && profile.category[ending] != candidate_pair) {
      ending = ending == 0 ? length - 1 : ending - 1;
      --region_length;
    }
    const std::size_t different = region_length - matching;
    if (region_length > 2 && matching > static_cast<double>(different) * 0.8) {
      const double background = static_cast<double>(profile.category_counts[candidate_pair]) /
          static_cast<double>(length);
      const double local = rdp_probability(region_length, matching, background, length);
      const double corrected = options_.correction == CorrectionMode::bonferroni
          ? std::min(1.0, local * static_cast<double>(correction_tests_))
          : std::min(1.0, local);
      if (!enforce_cutoff || corrected < options_.p_value_cutoff) {
        Signal signal;
        signal.triplet = profile.sequences;
        signal.recombinant = profile.sequences[recombinant_local];
        signal.major_parent = profile.sequences[major_local];
        signal.minor_parent = profile.sequences[minor_local];
        signal.beginning = profile.coordinates[beginning];
        signal.ending = profile.coordinates[ending];
        signal.wraps_origin = wrapped || signal.beginning > signal.ending;
        signal.informative_beginning = beginning + 1;
        signal.informative_ending = ending + 1;
        signal.local_p_value = local;
        signal.corrected_p_value = corrected;
        signal.correction_tests = correction_tests_;
        signal.pair_similarity = profile.similarities;
        signal.informative_sites = length;
        signal.candidate_pair = candidate_pair;
        result.push_back(signal);
      }
    }

    if (wrapped) break;
    search = std::max(search + 1, ending + 1);
  }
  return result;
}

std::vector<Signal> RdpScanner::triplet_signals(
    const std::array<std::uint32_t, 3>& triplet,
    bool enforce_cutoff) const {
  TripletProfile profile;
  if (!build_profile(triplet, profile)) return {};
  const auto order = ranked_pairs(profile);
  const std::array<double, 3> average{
      static_cast<double>(profile.category_counts[0]) / profile.category.size(),
      static_cast<double>(profile.category_counts[1]) / profile.category.size(),
      static_cast<double>(profile.category_counts[2]) / profile.category.size(),
  };
  std::vector<Signal> candidates;
  auto append = [&](std::uint8_t high, std::uint8_t candidate, std::uint8_t low) {
    auto found = candidate_signals(profile, high, candidate, low, enforce_cutoff);
    candidates.insert(candidates.end(), found.begin(), found.end());
  };
  append(order[0], order[1], order[2]);
  append(order[0], order[2], order[1]);
  if (average[order[0]] < 0.7) append(order[1], order[0], order[2]);
  return candidates;
}

void RdpScanner::scan_triplet(const std::array<std::uint32_t, 3>& triplet) {
  TripletProfile& profile = profile_scratch_;
  if (!build_profile(triplet, profile)) return;
  const auto order = ranked_pairs(profile);
  const std::array<double, 3> average{
      static_cast<double>(profile.category_counts[0]) / profile.category.size(),
      static_cast<double>(profile.category_counts[1]) / profile.category.size(),
      static_cast<double>(profile.category_counts[2]) / profile.category.size(),
  };
  std::vector<Signal> candidates;
  auto append = [&](std::uint8_t high, std::uint8_t candidate, std::uint8_t low) {
    auto found = candidate_signals(profile, high, candidate, low, true);
    candidates.insert(candidates.end(), found.begin(), found.end());
  };
  append(order[0], order[1], order[2]);
  append(order[0], order[2], order[1]);
  if (average[order[0]] < 0.7) append(order[1], order[0], order[2]);

  for (auto& signal : candidates) {
    map_signal_to_original(signal);
    const auto triplet = canonical_triplet(signal.triplet);
    auto& bucket = round_signal_index_[signal_signature(signal)];
    Signal* duplicate = nullptr;
    for (const std::uint32_t signal_id : bucket) {
      if (signal_id < round_signal_begin_ || signal_id >= signals_.size()) continue;
      auto& existing = signals_[signal_id];
      if (canonical_triplet(existing.triplet) == triplet &&
          existing.recombinant == signal.recombinant &&
          existing.beginning == signal.beginning &&
          existing.ending == signal.ending) {
        duplicate = &existing;
        break;
      }
    }
    if (duplicate) {
      if (signal.corrected_p_value < duplicate->corrected_p_value) {
        const auto id = duplicate->id;
        const auto review = duplicate->review_state;
        *duplicate = signal;
        duplicate->id = id;
        duplicate->review_state = review;
      }
      continue;
    }
    signal.id = static_cast<std::uint32_t>(signals_.size());
    bucket.push_back(signal.id);
    signals_.push_back(std::move(signal));
  }
}

double RdpScanner::tract_overlap(
    std::size_t first_beginning,
    std::size_t first_ending,
    std::size_t second_beginning,
    std::size_t second_ending) const {
  if (alignment_.length == 0) return 0.0;
  using Segment = std::array<std::size_t, 2>;
  const auto segments = [&](std::size_t beginning, std::size_t ending) {
    std::vector<Segment> result;
    beginning = std::clamp<std::size_t>(beginning, 1, alignment_.length);
    ending = std::clamp<std::size_t>(ending, 1, alignment_.length);
    if (beginning <= ending) {
      result.push_back({beginning, ending});
    } else {
      result.push_back({beginning, alignment_.length});
      result.push_back({1, ending});
    }
    return result;
  };
  const auto first = segments(first_beginning, first_ending);
  const auto second = segments(second_beginning, second_ending);
  const auto size = [](const std::vector<Segment>& values) {
    std::size_t total = 0;
    for (const auto& segment : values) total += segment[1] - segment[0] + 1;
    return total;
  };
  std::size_t intersection = 0;
  for (const auto& left : first) {
    for (const auto& right : second) {
      const std::size_t beginning = std::max(left[0], right[0]);
      const std::size_t ending = std::min(left[1], right[1]);
      if (beginning <= ending) intersection += ending - beginning + 1;
    }
  }
  const std::size_t combined = size(first) + size(second);
  return combined == 0
      ? 0.0
      : (2.0 * static_cast<double>(intersection)) / static_cast<double>(combined);
}

void RdpScanner::assign_event_support(
    UniqueEvent& event,
    std::vector<std::uint8_t>& assigned,
    const std::unordered_map<std::uint64_t, std::vector<std::uint32_t>>& pair_index) {
  event.support_signal_ids.clear();

  const auto consider = [&](std::uint32_t signal_id) {
    if (signal_id >= signals_.size() || assigned[signal_id] != 0) return false;
    const auto& signal = signals_[signal_id];
    if (signal.review_state == ReviewState::rejected) return false;
    if (tract_overlap(
            event.beginning,
            event.ending,
            signal.beginning,
            signal.ending) <= 0.3) {
      return false;
    }
    assigned[signal_id] = 1;
    signals_[signal_id].event_id = static_cast<std::int32_t>(event.id);
    event.support_signal_ids.push_back(signal_id);
    return true;
  };

  if (event.anchor_signal_id < signals_.size()) consider(event.anchor_signal_id);
  const std::array<std::uint32_t, 3> representatives{
      event.recombinant,
      event.major_parent,
      event.minor_parent,
  };
  const auto seed_pairs = triplet_pair_keys(representatives);
  for (const auto key : seed_pairs) {
    const auto found = pair_index.find(key);
    if (found == pair_index.end()) continue;
    for (const auto signal_id : found->second) consider(signal_id);
  }
  std::stable_sort(
      event.support_signal_ids.begin(),
      event.support_signal_ids.end(),
      [&](std::uint32_t left, std::uint32_t right) {
        const auto& a = signals_[left];
        const auto& b = signals_[right];
        if (a.corrected_p_value != b.corrected_p_value) {
          return a.corrected_p_value < b.corrected_p_value;
        }
        if (a.local_p_value != b.local_p_value) return a.local_p_value < b.local_p_value;
        return left < right;
      });
}

void RdpScanner::build_event_evidence(UniqueEvent& event) {
  for (auto& candidates : event.role_candidates) candidates.clear();
  event.detectable_sequences.clear();
  event.best_local_p_value = 1.0;
  event.best_corrected_p_value = 1.0;
  for (const auto signal_id : event.support_signal_ids) {
    if (signal_id >= signals_.size()) continue;
    const auto& signal = signals_[signal_id];
    event.role_candidates[0].push_back(signal.recombinant);
    event.role_candidates[1].push_back(signal.major_parent);
    event.role_candidates[2].push_back(signal.minor_parent);
    event.detectable_sequences.insert(
        event.detectable_sequences.end(), signal.triplet.begin(), signal.triplet.end());
    event.best_local_p_value = std::min(event.best_local_p_value, signal.local_p_value);
    event.best_corrected_p_value =
        std::min(event.best_corrected_p_value, signal.corrected_p_value);
  }
  event.role_candidates[0].push_back(event.recombinant);
  event.role_candidates[1].push_back(event.major_parent);
  event.role_candidates[2].push_back(event.minor_parent);
  for (auto& candidates : event.role_candidates) sort_unique(candidates);
  sort_unique(event.detectable_sequences);
}

void RdpScanner::refresh_trace_evidence(UniqueEvent& event) {
  event.trace_evidence.clear();
  if (event.major_parent == event.minor_parent) return;
  for (std::size_t sequence = 0; sequence < options_.mask.size(); ++sequence) {
    if (options_.mask[sequence] == 0 || sequence == event.major_parent ||
        sequence == event.minor_parent) {
      continue;
    }
    const std::array<std::uint32_t, 3> triplet{
        static_cast<std::uint32_t>(sequence),
        event.major_parent,
        event.minor_parent,
    };
    const auto candidates = triplet_signals(triplet, false);
    const Signal* best = nullptr;
    for (const auto& candidate : candidates) {
      if (candidate.recombinant != sequence) continue;
      if (tract_overlap(
              event.beginning,
              event.ending,
              candidate.beginning,
              candidate.ending) <= 0.3) {
        continue;
      }
      if (!best || candidate.local_p_value < best->local_p_value) best = &candidate;
    }
    if (!best) continue;
    TraceEvidence evidence;
    evidence.sequence = static_cast<std::uint32_t>(sequence);
    evidence.beginning = best->beginning;
    evidence.ending = best->ending;
    evidence.wraps_origin = best->wraps_origin;
    evidence.local_p_value = best->local_p_value;
    evidence.corrected_p_value = best->corrected_p_value;
    evidence.significant = best->corrected_p_value < options_.p_value_cutoff;
    event.trace_evidence.push_back(evidence);
  }
  std::stable_sort(
      event.trace_evidence.begin(),
      event.trace_evidence.end(),
      [](const TraceEvidence& left, const TraceEvidence& right) {
        if (left.significant != right.significant) return left.significant > right.significant;
        if (left.local_p_value != right.local_p_value) {
          return left.local_p_value < right.local_p_value;
        }
        return left.sequence < right.sequence;
      });
}

void RdpScanner::refresh_role_hypotheses(UniqueEvent& event) {
  const std::array<std::uint32_t, 3> reported_representatives{
      event.recombinant,
      event.major_parent,
      event.minor_parent,
  };
  std::array<std::uint32_t, 3> representatives = reported_representatives;
  if (event.anchor_signal_id < signals_.size()) {
    const Signal& anchor = signals_[event.anchor_signal_id];
    for (std::size_t role = 0; role < representatives.size(); ++role) {
      const auto member = std::find(
          anchor.triplet.begin(), anchor.triplet.end(), reported_representatives[role]);
      if (member == anchor.triplet.end()) continue;
      const std::size_t member_index = static_cast<std::size_t>(member - anchor.triplet.begin());
      const std::int32_t fragment_event = anchor.fragment_event_context[member_index];
      if (fragment_event < 0) continue;
      std::size_t best_sites = 0;
      for (std::size_t candidate = alignment_.sequence_count();
           candidate < working_alignment_.sequence_count();
           ++candidate) {
        if (candidate >= working_origins_.size() ||
            working_origins_[candidate] != reported_representatives[role] ||
            working_fragment_events_[candidate] != fragment_event) {
          continue;
        }
        std::size_t sites = 0;
        for (std::size_t coordinate = 1; coordinate <= working_alignment_.length; ++coordinate) {
          if (!coordinate_in_tract(
                  coordinate, event.beginning, event.ending, event.wraps_origin)) {
            continue;
          }
          if (working_alignment_.at(candidate, coordinate - 1) != 0) ++sites;
        }
        if (sites > best_sites) {
          best_sites = sites;
          representatives[role] = static_cast<std::uint32_t>(candidate);
        }
      }
    }
  }
  const CorrelationRegionLayout layout = build_correlation_regions(
      working_alignment_,
      representatives,
      event.beginning,
      event.ending);
  const CorrelationRegions& regions = layout.regions;
  const PhylogeneticRegions phylogenetic_regions =
      build_phylogenetic_regions(working_alignment_, layout);
  std::vector<std::array<std::size_t, 2>> overlap_sites(alignment_.sequence_count());
  for (std::size_t sequence = 0; sequence < overlap_sites.size(); ++sequence) {
    overlap_sites[sequence] = breakpoint_overlap_sites(
        working_alignment_, static_cast<std::uint32_t>(sequence), layout);
  }

  // FindSets in the supplied source first builds a detectable set for each of
  // the three possible recombinant anchors, then repeatedly closes the sets
  // across the other two anchors until no sequence is added.
  for (std::size_t role = 0; role < representatives.size(); ++role) {
    auto& hypothesis = event.role_hypotheses[role];
    hypothesis = {};
    hypothesis.presumed_recombinant = reported_representatives[role];
    hypothesis.parent_one = reported_representatives[(role + 1) % 3];
    hypothesis.parent_two = reported_representatives[(role + 2) % 3];
    hypothesis.detectable_signal_set.push_back(hypothesis.presumed_recombinant);

    for (const std::uint32_t signal_id : event.support_signal_ids) {
      if (signal_id >= signals_.size()) continue;
      const auto& triplet = signals_[signal_id].triplet;
      const bool has_first = std::find(
          triplet.begin(), triplet.end(), hypothesis.parent_one) != triplet.end();
      const bool has_second = std::find(
          triplet.begin(), triplet.end(), hypothesis.parent_two) != triplet.end();
      if (!has_first || !has_second) continue;
      for (const std::uint32_t sequence : triplet) {
        if (sequence != hypothesis.parent_one && sequence != hypothesis.parent_two) {
          hypothesis.detectable_signal_set.push_back(sequence);
        }
      }
    }
    if (role == 0) {
      for (const auto& trace : event.trace_evidence) {
        hypothesis.detectable_signal_set.push_back(trace.sequence);
      }
    }
    sort_unique(hypothesis.detectable_signal_set);
  }

  bool expanded = true;
  while (expanded) {
    expanded = false;
    for (std::size_t sequence = 0; sequence < alignment_.sequence_count(); ++sequence) {
      std::array<bool, 3> present{};
      for (std::size_t role = 0; role < 3; ++role) {
        const auto& set = event.role_hypotheses[role].detectable_signal_set;
        present[role] = std::binary_search(
            set.begin(), set.end(), static_cast<std::uint32_t>(sequence));
      }
      for (std::size_t role = 0; role < 3; ++role) {
        const std::size_t other_one = (role + 1) % 3;
        const std::size_t other_two = (role + 2) % 3;
        if (present[role] || !present[other_one] || !present[other_two]) continue;
        event.role_hypotheses[role].detectable_signal_set.push_back(
            static_cast<std::uint32_t>(sequence));
        sort_unique(event.role_hypotheses[role].detectable_signal_set);
        expanded = true;
      }
    }
  }

  std::vector<std::uint32_t> tree_candidates;
  tree_candidates.reserve(alignment_.sequence_count() + representatives.size());
  tree_candidates.insert(
      tree_candidates.end(), representatives.begin(), representatives.end());
  for (const std::uint32_t sequence : active_sequences_) {
    if (sequence < alignment_.sequence_count()) tree_candidates.push_back(sequence);
  }
  sort_unique(tree_candidates);
  const std::vector<std::uint32_t> tree_sequences = select_tree_sequences(
      working_alignment_, tree_candidates, representatives);
  event.tree_panel_sequences = tree_sequences.size();
  event.tree_panel_subsampled = tree_sequences.size() < tree_candidates.size();
  std::array<TreeRegionEvidence, 6> tree_evidence;
  const std::uint64_t event_seed = 0x6a09e667f3bcc909ULL ^
      (static_cast<std::uint64_t>(event.id + 1) << 48U) ^
      (static_cast<std::uint64_t>(event.beginning) << 24U) ^
      static_cast<std::uint64_t>(event.ending);
  for (std::size_t region = 0; region < tree_evidence.size(); ++region) {
    tree_evidence[region] = build_tree_region_evidence(
        working_alignment_,
        tree_sequences,
        phylogenetic_regions[region],
        kEventTreeBootstrapReplicates,
        event_seed ^ (0x9e3779b97f4a7c15ULL * (region + 1)));
    event.tree_regions[region] = {
        tree_evidence[region].site_count,
        tree_evidence[region].sequences.size(),
        tree_evidence[region].bootstrap_replicates,
        tree_evidence[region].supported_internal_branches,
        tree_evidence[region].internal_branches,
        tree_evidence[region].usable,
    };
  }

  std::array<std::vector<std::array<double, 3>>, 6> reference_distances;
  for (std::size_t region = 0; region < reference_distances.size(); ++region) {
    reference_distances[region].resize(alignment_.sequence_count());
    for (std::size_t sequence = 0; sequence < alignment_.sequence_count(); ++sequence) {
      reference_distances[region][sequence].fill(10.0);
      if (sequence < options_.mask.size() && options_.mask[sequence] != 0) continue;
      for (std::size_t reference = 0; reference < representatives.size(); ++reference) {
        if (tree_evidence[region].usable &&
            tree_evidence[region].contains(static_cast<std::uint32_t>(sequence))) {
          reference_distances[region][sequence][reference] = tree_evidence[region].jc(
              static_cast<std::uint32_t>(sequence), representatives[reference]);
        } else {
          reference_distances[region][sequence][reference] = jukes_cantor_distance(
              working_alignment_,
              static_cast<std::uint32_t>(sequence),
              representatives[reference],
              phylogenetic_regions[region]);
        }
      }
    }
  }

  const auto reference_distance = [&](
                                      std::size_t region,
                                      std::uint32_t first,
                                      std::uint32_t second) {
    for (std::size_t reference = 0; reference < representatives.size(); ++reference) {
      if (second == representatives[reference] && first < alignment_.sequence_count()) {
        return reference_distances[region][first][reference];
      }
      if (first == representatives[reference] && second < alignment_.sequence_count()) {
        return reference_distances[region][second][reference];
      }
    }
    return jukes_cantor_distance(
        working_alignment_, first, second, phylogenetic_regions[region]);
  };

  // MakeACOR uses the outside/inside phylogenetic distance matrices to reject
  // positive correlations that are incompatible with the representative
  // topology. Prefer the raw patristic matrices used by the source workflow;
  // a direct JC distance keeps the gate defined when the capped tree panel
  // does not contain a candidate.
  const auto affinity_distance = [&](
                                     std::size_t region,
                                     std::uint32_t first,
                                     std::uint32_t second) {
    if (tree_evidence[region].usable && tree_evidence[region].contains(first) &&
        tree_evidence[region].contains(second)) {
      return tree_evidence[region].tree(first, second, false);
    }
    return jukes_cantor_distance(
        working_alignment_, first, second, phylogenetic_regions[region]);
  };

  std::array<std::vector<std::uint8_t>, 3> acceptable_affinity;
  for (auto& role : acceptable_affinity) {
    role.assign(alignment_.sequence_count(), 0);
  }
  std::array<std::uint8_t, 2> closest_pair{};
  for (std::size_t side = 0; side < closest_pair.size(); ++side) {
    const std::size_t region = 4 + side;
    double minimum = std::numeric_limits<double>::infinity();
    for (std::uint8_t pair = 0; pair < 3; ++pair) {
      const auto members = pair_members(pair);
      const double distance = affinity_distance(
          region, representatives[members[0]], representatives[members[1]]);
      if (distance < minimum) {
        minimum = distance;
        closest_pair[side] = pair;
      }
    }
  }
  std::array<std::uint8_t, 3> in_list{};
  if (source_in_list(closest_pair[0], closest_pair[1], in_list)) {
    for (std::size_t sequence = 0; sequence < alignment_.sequence_count(); ++sequence) {
      const auto candidate = static_cast<std::uint32_t>(sequence);
      if (affinity_distance(4, representatives[in_list[0]], candidate) <
              affinity_distance(
                  4, representatives[in_list[2]], representatives[in_list[0]]) ||
          affinity_distance(5, representatives[in_list[1]], candidate) <
              affinity_distance(
                  5, representatives[in_list[1]], representatives[in_list[0]])) {
        acceptable_affinity[in_list[0]][sequence] = 1;
        acceptable_affinity[in_list[1]][sequence] = 1;
      }
      if (affinity_distance(4, representatives[in_list[2]], candidate) <
              affinity_distance(
                  4, representatives[in_list[2]], representatives[in_list[0]]) ||
          affinity_distance(5, representatives[in_list[2]], candidate) <
              affinity_distance(
                  5, representatives[in_list[2]], representatives[in_list[0]])) {
        acceptable_affinity[in_list[2]][sequence] = 1;
      }
    }
  }

  constexpr std::array<std::array<std::size_t, 2>, 3> tree_pairs{{
      {{0, 1}},
      {{2, 3}},
      {{4, 5}},
  }};

  for (std::size_t role = 0; role < representatives.size(); ++role) {
    auto& hypothesis = event.role_hypotheses[role];
    const std::uint32_t analysis_recombinant = representatives[role];
    const std::uint32_t analysis_parent_one = representatives[(role + 1) % 3];
    const std::uint32_t analysis_parent_two = representatives[(role + 2) % 3];

    const CorrelationProfiles anchor_profiles = correlation_profiles(
        working_alignment_,
        analysis_recombinant,
        analysis_parent_one,
        analysis_parent_two,
        analysis_recombinant,
        regions);
    hypothesis.correlation_warnings = correlation_warnings(anchor_profiles);
    for (std::size_t region = 0; region < 4; ++region) {
      const double first_second = reference_distance(
          region, analysis_recombinant, analysis_parent_one);
      const double first_third = reference_distance(
          region, analysis_recombinant, analysis_parent_two);
      const double second_third = reference_distance(
          region, analysis_parent_one, analysis_parent_two);
      const double total = first_second + first_third + second_third;
      const bool weak_triangle = total <= 0.0 ||
          ((1.0 - first_second * 2.0 / total) < 0.4 &&
           (1.0 - first_third * 2.0 / total) < 0.4 &&
           (1.0 - second_third * 2.0 / total) < 0.4);
      if (weak_triangle) hypothesis.correlation_warnings[region < 2 ? 0 : 1] = 1;
    }
    hypothesis.correlation_warnings[2] =
        hypothesis.correlation_warnings[0] != hypothesis.correlation_warnings[1] ? 1 : 0;
    std::vector<std::uint8_t> detectable(alignment_.sequence_count(), 0);
    for (const std::uint32_t sequence : hypothesis.detectable_signal_set) {
      if (sequence < detectable.size()) detectable[sequence] = 1;
    }

    for (std::size_t sequence = 0; sequence < alignment_.sequence_count(); ++sequence) {
      if (sequence == hypothesis.parent_one || sequence == hypothesis.parent_two) continue;
      ++hypothesis.tested_sequences;
      const CorrelationProfiles candidate_profiles = correlation_profiles(
          working_alignment_,
          analysis_recombinant,
          analysis_parent_one,
          analysis_parent_two,
          static_cast<std::uint32_t>(sequence),
          regions);
      DistanceCorrelationEvidence evidence;
      evidence.sequence = static_cast<std::uint32_t>(sequence);
      evidence.detectable_support = detectable[sequence] != 0;
      evidence.acceptable_affinity = acceptable_affinity[role][sequence] != 0;
      evidence.warning_filtered = hypothesis.correlation_warnings;
      evidence.breakpoint_overlap_sites = overlap_sites[sequence];
      evidence.overlap_eligible = evidence.breakpoint_overlap_sites[0] > 10 ||
          evidence.breakpoint_overlap_sites[1] > 10;
      bool has_valid_pair = false;
      double best_p_value = 1.0;

      for (std::size_t pair = 0; pair < 3; ++pair) {
        evidence.minimum_comparable_sites[pair] = std::min(
            correlation_minimum_comparable_sites(anchor_profiles, pair),
            correlation_minimum_comparable_sites(candidate_profiles, pair));
        if (!correlation_pair_has_categories(anchor_profiles, pair) ||
            !correlation_pair_has_categories(candidate_profiles, pair)) {
          continue;
        }
        double correlation = 0.0;
        double direct_correlation = 0.0;
        std::uint8_t inversion_code = 0;
        if (!best_category_correlation(
                correlation_vector(anchor_profiles, pair),
                correlation_vector(candidate_profiles, pair),
                direct_correlation,
                correlation,
                inversion_code)) {
          continue;
        }
        has_valid_pair = true;
        evidence.direct_correlations[pair] = direct_correlation;
        evidence.correlations[pair] = correlation;
        evidence.inversion_codes[pair] = inversion_code;
        evidence.p_values[pair] = pearson_six_two_sided_p(correlation);
        const bool passes_correlation = inversion_code == 0
            ? correlation > 0.0 && evidence.p_values[pair] < kDistanceCorrelationCutoff
            : correlation > 0.83;
        if (hypothesis.correlation_warnings[pair] == 0 && passes_correlation &&
            evidence.p_values[pair] < best_p_value) {
          best_p_value = evidence.p_values[pair];
          evidence.best_matrix_pair = static_cast<std::int8_t>(pair);
        }
      }
      if (has_valid_pair && evidence.overlap_eligible) ++hypothesis.valid_sequences;
      std::size_t strong_correlation_count = 0;
      for (std::size_t pair = 0; pair < 2; ++pair) {
        if (evidence.warning_filtered[pair] != 0 || evidence.correlations[pair] <= 0.95) {
          continue;
        }
        strong_correlation_count += evidence.correlations[pair] > 0.98 ? 2 : 1;
      }
      // MakeRList checks equality rather than >= here; preserve that quirk.
      evidence.strong_correlation_override = strong_correlation_count == 2;
      std::size_t warning_count = 0;
      for (std::size_t pair = 0; pair < 3; ++pair) {
        if (evidence.warning_filtered[pair] != 0) {
          ++warning_count;
          continue;
        }
        if (evidence.inversion_codes[pair] != 0) {
          if (evidence.correlations[pair] > 0.83) evidence.inverse_support = true;
          continue;
        }
        const double probability_score = 1.0 - evidence.p_values[pair];
        if (probability_score >= 0.95) {
          evidence.aggregate_score += (probability_score - 0.95) / 0.05;
        }
      }
      evidence.aggregate_target = warning_count > 0
          ? 0.9 - 0.3 * static_cast<double>(warning_count)
          : 0.9;
      evidence.positive_support = evidence.overlap_eligible &&
          (evidence.acceptable_affinity || evidence.strong_correlation_override) &&
          evidence.aggregate_score > evidence.aggregate_target;
      evidence.stripped_inverse_only = evidence.inverse_support &&
          !evidence.positive_support && sequence != hypothesis.presumed_recombinant;
      // MakeRList temporarily admits inverse-only rows, but the active
      // StripDupInv call removes them before the co-recombinant lists proceed.
      evidence.significant = sequence == hypothesis.presumed_recombinant ||
          evidence.positive_support;
      if (evidence.significant) {
        hypothesis.distance_correlation_set.push_back(evidence.sequence);
      }
      if (evidence.significant || evidence.detectable_support || evidence.inverse_support) {
        hypothesis.distance_evidence.push_back(evidence);
      }
    }

    sort_unique(hypothesis.distance_correlation_set);
    std::stable_sort(
        hypothesis.distance_evidence.begin(),
        hypothesis.distance_evidence.end(),
        [](const DistanceCorrelationEvidence& left, const DistanceCorrelationEvidence& right) {
          if (left.significant != right.significant) return left.significant > right.significant;
          const double left_p = left.best_matrix_pair < 0
              ? 1.0
              : left.p_values[static_cast<std::size_t>(left.best_matrix_pair)];
          const double right_p = right.best_matrix_pair < 0
              ? 1.0
              : right.p_values[static_cast<std::size_t>(right.best_matrix_pair)];
          if (left_p != right_p) return left_p < right_p;
          return left.sequence < right.sequence;
        });

    const auto in_set = [](const std::vector<std::uint32_t>& set, std::uint32_t sequence) {
      return std::binary_search(set.begin(), set.end(), sequence);
    };
    for (std::size_t sequence = 0; sequence < alignment_.sequence_count(); ++sequence) {
      if (sequence == hypothesis.parent_one || sequence == hypothesis.parent_two) continue;
      PhylogeneticCorrelationEvidence evidence;
      evidence.sequence = static_cast<std::uint32_t>(sequence);
      evidence.masked_excluded = sequence < options_.mask.size() && options_.mask[sequence] != 0;
      evidence.distance_fallback = !tree_evidence[0].contains(evidence.sequence) ||
          std::any_of(
              tree_evidence.begin(),
              tree_evidence.end(),
              [](const TreeRegionEvidence& region) { return !region.usable; });
      if (evidence.masked_excluded) {
        evidence.distance_fallback = false;
      } else if (sequence == hypothesis.presumed_recombinant) {
        evidence.included = true;
        evidence.supporting_tree_pairs = 3;
        evidence.best_tree_pair = 2;
        evidence.collapsed_pair_support.fill(1);
        evidence.raw_pair_support.fill(1);
      } else {
        double best_margin = -std::numeric_limits<double>::infinity();
        for (std::size_t pair = 0; pair < tree_pairs.size(); ++pair) {
          bool collapsed_support = true;
          bool raw_support = true;
          double collapsed_pair_margin = std::numeric_limits<double>::infinity();
          double raw_pair_margin = std::numeric_limits<double>::infinity();
          for (const std::size_t region : tree_pairs[pair]) {
            double collapsed_margin = -10.0;
            double raw_margin = -10.0;
            if (!evidence.distance_fallback && tree_evidence[region].usable) {
              collapsed_margin = affinity_margin(
                  tree_evidence[region].tree(analysis_recombinant, evidence.sequence),
                  tree_evidence[region].tree(analysis_recombinant, analysis_parent_one),
                  tree_evidence[region].tree(analysis_recombinant, analysis_parent_two),
                  tree_evidence[region].tree(evidence.sequence, analysis_parent_one),
                  tree_evidence[region].tree(evidence.sequence, analysis_parent_two));
              raw_margin = affinity_margin(
                  tree_evidence[region].tree(
                      analysis_recombinant, evidence.sequence, false),
                  tree_evidence[region].tree(
                      analysis_recombinant, analysis_parent_one, false),
                  tree_evidence[region].tree(
                      analysis_recombinant, analysis_parent_two, false),
                  tree_evidence[region].tree(evidence.sequence, analysis_parent_one, false),
                  tree_evidence[region].tree(evidence.sequence, analysis_parent_two, false));
            } else {
              const auto distance = [&](std::uint32_t first, std::uint32_t second) {
                return reference_distance(region, first, second);
              };
              collapsed_margin = affinity_margin(
                  distance(analysis_recombinant, evidence.sequence),
                  distance(analysis_recombinant, analysis_parent_one),
                  distance(analysis_recombinant, analysis_parent_two),
                  distance(evidence.sequence, analysis_parent_one),
                  distance(evidence.sequence, analysis_parent_two));
              raw_margin = collapsed_margin;
            }
            collapsed_pair_margin = std::min(collapsed_pair_margin, collapsed_margin);
            raw_pair_margin = std::min(raw_pair_margin, raw_margin);
            collapsed_support = collapsed_support && collapsed_margin > 1e-12;
            raw_support = raw_support && raw_margin > 1e-12;
          }
          evidence.collapsed_affinity_margin[pair] = collapsed_pair_margin;
          evidence.raw_affinity_margin[pair] = raw_pair_margin;
          evidence.collapsed_pair_support[pair] = collapsed_support ? 1 : 0;
          evidence.raw_pair_support[pair] = raw_support ? 1 : 0;
          if (collapsed_support) {
            ++evidence.supporting_tree_pairs;
            if (collapsed_pair_margin > best_margin) {
              best_margin = collapsed_pair_margin;
              evidence.best_tree_pair = static_cast<std::int8_t>(pair);
            }
          }
        }
        evidence.included = evidence.supporting_tree_pairs > 0;
      }
      if (evidence.included) {
        hypothesis.phylogenetic_correlation_set.push_back(evidence.sequence);
      }
      if (evidence.included ||
          in_set(hypothesis.detectable_signal_set, evidence.sequence) ||
          in_set(hypothesis.distance_correlation_set, evidence.sequence)) {
        hypothesis.phylogenetic_evidence.push_back(evidence);
      }
    }
    sort_unique(hypothesis.phylogenetic_correlation_set);

    for (std::size_t sequence = 0; sequence < alignment_.sequence_count(); ++sequence) {
      const auto candidate = static_cast<std::uint32_t>(sequence);
      const std::size_t evidence_count =
          static_cast<std::size_t>(in_set(hypothesis.detectable_signal_set, candidate)) +
          static_cast<std::size_t>(in_set(hypothesis.distance_correlation_set, candidate)) +
          static_cast<std::size_t>(in_set(hypothesis.phylogenetic_correlation_set, candidate));
      if (candidate == hypothesis.presumed_recombinant || evidence_count >= 2) {
        hypothesis.complete_two_of_three_set.push_back(candidate);
      }
    }
    sort_unique(hypothesis.complete_two_of_three_set);
    std::stable_sort(
        hypothesis.phylogenetic_evidence.begin(),
        hypothesis.phylogenetic_evidence.end(),
        [](const PhylogeneticCorrelationEvidence& left,
           const PhylogeneticCorrelationEvidence& right) {
          if (left.included != right.included) return left.included > right.included;
          if (left.distance_fallback != right.distance_fallback) {
            return left.distance_fallback < right.distance_fallback;
          }
          if (left.supporting_tree_pairs != right.supporting_tree_pairs) {
            return left.supporting_tree_pairs > right.supporting_tree_pairs;
          }
          return left.sequence < right.sequence;
        });
  }

  // The active first FinalTrim stage copies the direct-polarity correlations,
  // clears warning/inversion rows, counts the same candidate/pair across all
  // three RLists, and suppresses a pair wherever that count exceeds one. Keep
  // this as separately auditable evidence; later unported FinalTrim branches
  // are not yet allowed to prune the browser's two-of-three group.
  for (std::size_t sequence = 0; sequence < alignment_.sequence_count(); ++sequence) {
    for (std::size_t pair = 0; pair < 3; ++pair) {
      std::size_t duplicate_count = 0;
      for (const auto& hypothesis : event.role_hypotheses) {
        for (const auto& evidence : hypothesis.distance_evidence) {
          if (evidence.sequence != sequence || !evidence.significant ||
              evidence.warning_filtered[pair] != 0 ||
              evidence.inversion_codes[pair] != 0) {
            continue;
          }
          if (evidence.correlations[pair] > 0.83) ++duplicate_count;
          break;
        }
      }
      if (duplicate_count <= 1) continue;
      for (auto& hypothesis : event.role_hypotheses) {
        for (auto& evidence : hypothesis.distance_evidence) {
          if (evidence.sequence == sequence && evidence.significant &&
              evidence.warning_filtered[pair] == 0 &&
              evidence.inversion_codes[pair] == 0 &&
              evidence.correlations[pair] > 0.83) {
            evidence.duplicate_filtered[pair] = 1;
          }
        }
      }
    }
  }
  for (auto& hypothesis : event.role_hypotheses) {
    for (auto& evidence : hypothesis.distance_evidence) {
      if (!evidence.significant) continue;
      for (std::size_t pair = 0; pair < 3; ++pair) {
        if (evidence.warning_filtered[pair] == 0 &&
            evidence.inversion_codes[pair] == 0 &&
            evidence.duplicate_filtered[pair] == 0 &&
            evidence.correlations[pair] >= 0.83) {
          evidence.duplicate_cleaned_support = true;
          break;
        }
      }
    }
  }

  event.automatic_co_recombinant_sequences =
      event.role_hypotheses[0].complete_two_of_three_set;
  if (!event.group_manual_adjusted) {
    event.co_recombinant_sequences = event.automatic_co_recombinant_sequences;
  }

  event.role_consensus = {};
  constexpr std::size_t metric_count = 9;
  std::array<std::array<bool, 3>, metric_count> metric_valid{};
  std::array<RoleMetricEvidence, metric_count> metrics;
  metrics[0].kind = RoleMetricKind::phpr;
  metrics[0].weight = 8.0;
  metrics[0].higher_is_recombinant = false;
  metrics[1].kind = RoleMetricKind::tree_phpr;
  metrics[1].weight = 18.0;
  metrics[1].higher_is_recombinant = false;
  metrics[2].kind = RoleMetricKind::collapsed_tree_phpr;
  metrics[2].weight = 20.0;
  metrics[2].higher_is_recombinant = false;
  metrics[3].kind = RoleMetricKind::sub_phpr;
  metrics[3].weight = 10.0;
  metrics[3].higher_is_recombinant = true;
  metrics[4].kind = RoleMetricKind::tree_sub_phpr;
  metrics[4].weight = 8.0;
  metrics[4].higher_is_recombinant = true;
  metrics[5].kind = RoleMetricKind::subdist;
  metrics[5].weight = 2.0;
  metrics[5].higher_is_recombinant = true;
  metrics[6].kind = RoleMetricKind::tree_subdist;
  metrics[6].weight = 10.0;
  metrics[6].higher_is_recombinant = true;
  metrics[7].kind = RoleMetricKind::triplet_score;
  metrics[7].weight = 8.0;
  metrics[7].higher_is_recombinant = true;
  metrics[8].kind = RoleMetricKind::three_set_support;
  // Set membership is reported as context. The supplied decision-tree code
  // records it for later classifiers but does not give it a standalone vote.
  metrics[8].weight = 0.0;
  metrics[8].higher_is_recombinant = true;

  const auto jc_distance = [&](std::uint32_t first, std::uint32_t second) {
    return std::pair{
        tree_evidence[4].jc(first, second),
        tree_evidence[5].jc(first, second),
    };
  };
  const auto raw_tree_distance = [&](std::uint32_t first, std::uint32_t second) {
    return std::pair{
        tree_evidence[4].tree(first, second, false),
        tree_evidence[5].tree(first, second, false),
    };
  };
  const auto collapsed_tree_distance = [&](std::uint32_t first, std::uint32_t second) {
    return std::pair{
        tree_evidence[4].tree(first, second, true),
        tree_evidence[5].tree(first, second, true),
    };
  };
  const std::vector<std::uint32_t> raw_involved = source_involved_sequences(
      representatives, tree_sequences, jc_distance, raw_tree_distance);
  const std::vector<std::uint32_t> collapsed_involved = source_involved_sequences(
      representatives, tree_sequences, jc_distance, collapsed_tree_distance);
  const SourcePhylproScores jc_scores = source_phylpro_scores(
      representatives, raw_involved, jc_distance);
  const SourcePhylproScores raw_tree_scores = source_phylpro_scores(
      representatives, raw_involved, raw_tree_distance);
  const SourcePhylproScores collapsed_tree_scores = source_phylpro_scores(
      representatives, collapsed_involved, collapsed_tree_distance);

  metrics[0].scores = jc_scores.primary;
  metrics[1].scores = raw_tree_scores.primary;
  metrics[2].scores = collapsed_tree_scores.primary;
  metrics[3].scores = jc_scores.leave_one_out;
  metrics[4].scores = raw_tree_scores.leave_one_out;
  metrics[5].scores = jc_scores.displacement;
  metrics[6].scores = raw_tree_scores.displacement;
  for (std::size_t role = 0; role < representatives.size(); ++role) {
    metric_valid[0][role] = jc_scores.primary_valid[role];
    metric_valid[1][role] = raw_tree_scores.primary_valid[role];
    metric_valid[2][role] = collapsed_tree_scores.primary_valid[role];
    metric_valid[3][role] = jc_scores.leave_one_out_valid[role];
    metric_valid[4][role] = raw_tree_scores.leave_one_out_valid[role];
    metric_valid[5][role] = jc_scores.displacement_valid[role];
    metric_valid[6][role] = raw_tree_scores.displacement_valid[role];
    metrics[7].scores[role] = source_triplet_score(
        role, representatives, tree_sequences, raw_tree_distance);
    metric_valid[7][role] = tree_sequences.size() > 11 &&
        tree_evidence[4].usable && tree_evidence[5].usable;
    metrics[8].scores[role] = static_cast<double>(
        event.role_hypotheses[role].complete_two_of_three_set.size());
    metric_valid[8][role] = true;
  }

  // MakeConsensusC disables each PhylPro family when an anchor correlation is
  // exactly +/-1. Apply that family guard to its leave-one-out and displacement
  // derivatives as well.
  const auto family_informative = [](const SourcePhylproScores& scores) {
    return std::all_of(
               scores.primary_valid.begin(), scores.primary_valid.end(),
               [](bool value) { return value; }) &&
        std::none_of(scores.primary.begin(), scores.primary.end(), [](double value) {
          return std::abs(value) >= 0.999999;
        });
  };
  if (!family_informative(jc_scores)) {
    metric_valid[0].fill(false);
    metric_valid[3].fill(false);
    metric_valid[5].fill(false);
  }
  if (!family_informative(raw_tree_scores)) {
    metric_valid[1].fill(false);
    metric_valid[4].fill(false);
    metric_valid[6].fill(false);
  }
  if (!family_informative(collapsed_tree_scores)) metric_valid[2].fill(false);

  for (std::size_t metric_index = 0; metric_index < metrics.size(); ++metric_index) {
    auto& metric = metrics[metric_index];
    if (!std::all_of(
            metric_valid[metric_index].begin(),
            metric_valid[metric_index].end(),
            [](bool value) { return value; })) {
      event.role_consensus.metrics.push_back(metric);
      continue;
    }
    const auto [minimum, maximum] = std::minmax_element(
        metric.scores.begin(), metric.scores.end());
    if (*maximum - *minimum <= 1e-9) {
      event.role_consensus.metrics.push_back(metric);
      continue;
    }
    const double best = metric.higher_is_recombinant ? *maximum : *minimum;
    std::vector<std::size_t> winners;
    for (std::size_t role = 0; role < 3; ++role) {
      if (std::abs(metric.scores[role] - best) <= 1e-9) winners.push_back(role);
    }
    metric.informative = true;
    if (winners.size() == 1) metric.winning_role = static_cast<std::int8_t>(winners.front());

    if (metric_index == 2) {
      // Native special case: the collapsed-tree anchor wins 20 points only
      // when removing that same role also maximises SubPhPr.
      if (metrics[3].informative || std::all_of(
              metric_valid[3].begin(), metric_valid[3].end(),
              [](bool value) { return value; })) {
        for (std::size_t role = 0; role < representatives.size(); ++role) {
          const std::size_t other_one = (role + 1) % 3;
          const std::size_t other_two = (role + 2) % 3;
          if (metric.scores[role] <= metric.scores[other_one] &&
              metric.scores[role] <= metric.scores[other_two] &&
              metrics[3].scores[role] >= metrics[3].scores[other_one] &&
              metrics[3].scores[role] >= metrics[3].scores[other_two]) {
            metric.contributions[role] = metric.weight;
          }
        }
      }
    } else if (metric.weight > 0.0) {
      for (std::size_t role = 0; role < representatives.size(); ++role) {
        const std::size_t other_one = (role + 1) % 3;
        const std::size_t other_two = (role + 2) % 3;
        const auto no_worse = [&](std::size_t other) {
          return metric.higher_is_recombinant
              ? metric.scores[role] >= metric.scores[other]
              : metric.scores[role] <= metric.scores[other];
        };
        const auto strictly_better = [&](std::size_t other) {
          return metric.higher_is_recombinant
              ? metric.scores[role] > metric.scores[other]
              : metric.scores[role] < metric.scores[other];
        };
        if (no_worse(other_one) && no_worse(other_two)) {
          metric.contributions[role] = metric.weight;
        } else if (strictly_better(other_one) || strictly_better(other_two)) {
          metric.contributions[role] = metric.weight / 2.0;
        }
      }
    }
    for (std::size_t role = 0; role < representatives.size(); ++role) {
      event.role_consensus.votes[role] += metric.contributions[role];
    }
    event.role_consensus.metrics.push_back(metric);
  }

  const double total_votes = std::accumulate(
      event.role_consensus.votes.begin(), event.role_consensus.votes.end(), 0.0);
  event.role_consensus.recommended_recombinant = event.recombinant;
  event.role_consensus.recommended_major_parent = event.major_parent;
  event.role_consensus.recommended_minor_parent = event.minor_parent;
  if (total_votes > 0.0) {
    std::array<std::size_t, 3> order{0, 1, 2};
    std::stable_sort(order.begin(), order.end(), [&](std::size_t left, std::size_t right) {
      if (event.role_consensus.votes[left] != event.role_consensus.votes[right]) {
        return event.role_consensus.votes[left] > event.role_consensus.votes[right];
      }
      return left < right;
    });
    if (event.role_consensus.votes[order[0]] -
            event.role_consensus.votes[order[1]] <= 1e-9) {
      return;
    }
    const std::size_t recommended_role = order[0];
    event.role_consensus.informative = true;
    event.role_consensus.recommended_role = static_cast<std::int8_t>(recommended_role);
    event.role_consensus.confidence =
        (event.role_consensus.votes[order[0]] - event.role_consensus.votes[order[1]]) /
        total_votes;
    const auto& hypothesis = event.role_hypotheses[recommended_role];
    event.role_consensus.recommended_recombinant = hypothesis.presumed_recombinant;
    const double outside_first = tree_evidence[4].jc(
        representatives[recommended_role], representatives[(recommended_role + 1) % 3]);
    const double outside_second = tree_evidence[4].jc(
        representatives[recommended_role], representatives[(recommended_role + 2) % 3]);
    if (outside_second < outside_first) {
      event.role_consensus.recommended_major_parent = hypothesis.parent_two;
      event.role_consensus.recommended_minor_parent = hypothesis.parent_one;
    } else {
      event.role_consensus.recommended_major_parent = hypothesis.parent_one;
      event.role_consensus.recommended_minor_parent = hypothesis.parent_two;
    }
  }
}

bool RdpScanner::matches_fixed_event(const Signal& signal) const {
  const std::size_t limit = std::min(fixed_event_count_, events_.size());
  for (std::size_t index = 0; index < limit; ++index) {
    const auto& event = events_[index];
    const std::array<std::uint32_t, 3> representatives{
        event.recombinant,
        event.major_parent,
        event.minor_parent,
    };
    std::size_t shared = 0;
    for (const std::uint32_t sequence : signal.triplet) {
      if (std::find(representatives.begin(), representatives.end(), sequence) !=
          representatives.end()) {
        ++shared;
      }
    }
    if (shared >= 2 && tract_overlap(
            event.beginning,
            event.ending,
            signal.beginning,
            signal.ending) > 0.3) {
      return true;
    }
  }
  return false;
}

RdpScanner::ErasureResult RdpScanner::erase_event_tract(const UniqueEvent& event) {
  std::vector<std::uint32_t> group = event.co_recombinant_sequences;
  if (group.empty()) group.push_back(event.recombinant);
  sort_unique(group);
  ErasureResult result;
  std::vector<std::uint8_t> selected_origins(alignment_.sequence_count(), 0);
  for (const std::uint32_t sequence : group) {
    if (sequence < selected_origins.size()) selected_origins[sequence] = 1;
  }

  const std::size_t existing_sequence_count = working_alignment_.sequence_count();
  const bool retain_fragments =
      working_alignment_.length < kFragmentReentryAlignmentLengthLimit;
  const std::size_t fragment_minimum = std::max({
      std::size_t{5},
      options_.window_sites,
      (working_alignment_.length + 99) / 100,
  });
  for (std::size_t sequence = 0; sequence < existing_sequence_count; ++sequence) {
    if (sequence >= working_origins_.size()) continue;
    const std::uint32_t origin = working_origins_[sequence];
    if (origin >= selected_origins.size() || selected_origins[origin] == 0) continue;
    std::vector<std::uint8_t> fragment;
    if (retain_fragments) fragment.assign(working_alignment_.length, 0);
    std::size_t fragment_sites = 0;
    for (std::size_t coordinate = 1; coordinate <= working_alignment_.length; ++coordinate) {
      if (!coordinate_in_tract(
              coordinate,
              event.beginning,
              event.ending,
              event.wraps_origin)) {
        continue;
      }
      const std::size_t offset =
          static_cast<std::size_t>(sequence) * working_alignment_.length + coordinate - 1;
      if (working_alignment_.states[offset] == 0) continue;
      if (retain_fragments) fragment[coordinate - 1] = working_alignment_.states[offset];
      working_alignment_.states[offset] = 0;
      ++fragment_sites;
      ++result.working_sites;
      if (sequence < alignment_.sequence_count()) ++result.original_sites;
    }
    if (!retain_fragments || fragment_sites < fragment_minimum) continue;

    bool duplicate = false;
    for (std::size_t candidate = alignment_.sequence_count();
         candidate < working_alignment_.sequence_count();
         ++candidate) {
      if (candidate >= working_origins_.size() || working_origins_[candidate] != origin) continue;
      const auto first = working_alignment_.states.begin() +
          static_cast<std::ptrdiff_t>(candidate * working_alignment_.length);
      if (std::equal(fragment.begin(), fragment.end(), first)) {
        duplicate = true;
        break;
      }
    }
    if (duplicate) continue;
    if (working_alignment_.sequence_count() - alignment_.sequence_count() >=
        kWorkingFragmentSequenceCap) {
      fragment_reentry_capped_ = true;
      continue;
    }

    const std::size_t fragment_number =
        working_alignment_.sequence_count() - alignment_.sequence_count() + 1;
    working_alignment_.names.push_back(
        alignment_.names[origin] + "|rdp-fragment-e" + std::to_string(event.id + 1) +
        "-" + std::to_string(fragment_number));
    working_alignment_.sequences.emplace_back(working_alignment_.length, '-');
    working_alignment_.states.insert(
        working_alignment_.states.end(), fragment.begin(), fragment.end());
    working_alignment_.sequence_summaries.push_back(
        {fragment_sites, working_alignment_.length - fragment_sites});
    working_origins_.push_back(origin);
    working_fragment_events_.push_back(static_cast<std::int32_t>(event.id));
    ++result.fragments_added;
  }
  refresh_active_sequences();
  return result;
}

bool RdpScanner::finish_detection_round(std::string& error) {
  std::uint32_t anchor_id = std::numeric_limits<std::uint32_t>::max();
  for (std::size_t index = round_signal_begin_; index < signals_.size(); ++index) {
    const auto& candidate = signals_[index];
    if (candidate.review_state == ReviewState::rejected || matches_fixed_event(candidate)) {
      continue;
    }
    if (anchor_id == std::numeric_limits<std::uint32_t>::max()) {
      anchor_id = static_cast<std::uint32_t>(index);
      continue;
    }
    const auto& anchor = signals_[anchor_id];
    if (candidate.corrected_p_value < anchor.corrected_p_value ||
        (candidate.corrected_p_value == anchor.corrected_p_value &&
         (candidate.local_p_value < anchor.local_p_value ||
          (candidate.local_p_value == anchor.local_p_value && index < anchor_id)))) {
      anchor_id = static_cast<std::uint32_t>(index);
    }
  }
  if (anchor_id == std::numeric_limits<std::uint32_t>::max()) {
    signals_.resize(round_signal_begin_);
    cycle_termination_ = "no-significant-signals";
    return false;
  }

  std::vector<std::uint8_t> assigned(signals_.size(), 1);
  std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> pair_index;
  for (std::size_t index = round_signal_begin_; index < signals_.size(); ++index) {
    if (matches_fixed_event(signals_[index])) continue;
    assigned[index] = 0;
    for (const std::uint64_t key : triplet_pair_keys(signals_[index].triplet)) {
      pair_index[key].push_back(static_cast<std::uint32_t>(index));
    }
  }

  const auto& anchor = signals_[anchor_id];
  UniqueEvent event;
  event.id = static_cast<std::uint32_t>(events_.size());
  event.anchor_signal_id = anchor_id;
  event.recombinant = anchor.recombinant;
  event.major_parent = anchor.major_parent;
  event.minor_parent = anchor.minor_parent;
  event.beginning = anchor.beginning;
  event.ending = anchor.ending;
  event.wraps_origin = anchor.wraps_origin;
  event.detection_round = scan_round_;
  event.fragment_assisted_detection = anchor.fragment_assisted;
  assign_event_support(event, assigned, pair_index);
  if (event.support_signal_ids.empty()) {
    error = "The strongest cyclic RDP signal could not be assigned to an event.";
    cycle_termination_ = "event-assignment-error";
    return false;
  }
  build_event_evidence(event);
  if (signals_[event.anchor_signal_id].correction_tests > 0) {
    correction_tests_ = signals_[event.anchor_signal_id].correction_tests;
  }
  refresh_trace_evidence(event);
  refresh_role_hypotheses(event);

  if (event.role_consensus.informative &&
      (event.role_consensus.recommended_recombinant != event.recombinant ||
       event.role_consensus.recommended_major_parent != event.major_parent ||
       event.role_consensus.recommended_minor_parent != event.minor_parent)) {
    event.recombinant = event.role_consensus.recommended_recombinant;
    event.major_parent = event.role_consensus.recommended_major_parent;
    event.minor_parent = event.role_consensus.recommended_minor_parent;
    build_event_evidence(event);
    refresh_trace_evidence(event);
    refresh_role_hypotheses(event);
  }

  const ErasureResult erasure = erase_event_tract(event);
  event.erased_nucleotide_sites = erasure.original_sites;
  event.erased_working_sites = erasure.working_sites;
  event.fragment_sequences_added = erasure.fragments_added;
  event.tract_erased_for_detection = erasure.working_sites > 0;

  // Signals that were not grouped with the selected event are deliberately
  // discarded: the modified alignment is screened again and only signals
  // that survive that next round can seed a later unique event.
  std::vector<Signal> retained;
  retained.reserve(event.support_signal_ids.size());
  std::uint32_t compact_anchor = 0;
  for (const std::uint32_t signal_id : event.support_signal_ids) {
    Signal signal = signals_[signal_id];
    signal.id = static_cast<std::uint32_t>(round_signal_begin_ + retained.size());
    signal.event_id = static_cast<std::int32_t>(event.id);
    if (signal_id == event.anchor_signal_id) compact_anchor = signal.id;
    retained.push_back(std::move(signal));
  }
  signals_.resize(round_signal_begin_);
  event.support_signal_ids.clear();
  for (auto& signal : retained) {
    event.support_signal_ids.push_back(signal.id);
    signals_.push_back(std::move(signal));
  }
  event.anchor_signal_id = compact_anchor;
  events_.push_back(std::move(event));

  if (events_.back().erased_working_sites == 0) {
    cycle_termination_ = "no-new-tract-sites";
    return false;
  }
  if (active_sequences_.size() < 3 || total_triplets_ == 0) {
    cycle_termination_ = "fewer-than-three-active-origins";
    return false;
  }
  cycle_termination_ = "scanning";
  return true;
}

bool RdpScanner::event_action_allowed(
    std::uint32_t event_id,
    std::string& error) const {
  if (event_id >= events_.size()) {
    error = "The selected reconciled event does not exist.";
    return false;
  }
  if (reconciliation_required_after_ >= 0 &&
      event_id != static_cast<std::uint32_t>(reconciliation_required_after_)) {
    error = "Finish the pending downstream rebuild before changing another event.";
    return false;
  }
  const auto first_unreviewed = std::find_if(
      events_.begin(), events_.end(), [](const UniqueEvent& event) {
        return event.review_state == ReviewState::unreviewed;
      });
  if (first_unreviewed != events_.end() &&
      event_id > static_cast<std::uint32_t>(first_unreviewed - events_.begin())) {
    error = "Record the earlier undecided event before changing this event.";
    return false;
  }
  return true;
}

bool RdpScanner::set_event_review_state(
    std::uint32_t event_id,
    ReviewState state,
    std::string& error) {
  if (!event_action_allowed(event_id, error)) return false;
  auto& event = events_[event_id];
  event.review_state = state;

  // Every detected event influenced the working alignment used to discover
  // later events. Rejecting one therefore invalidates the downstream chain:
  // its erased tract must be restored before screening resumes. Keep the
  // marker even if the decision is toggled again; a redundant rebuild is safer
  // than losing an earlier invalidation marker.
  if (state == ReviewState::rejected && event.tract_erased_for_detection) {
    if (reconciliation_required_after_ < 0 ||
        event_id < static_cast<std::uint32_t>(reconciliation_required_after_)) {
      reconciliation_required_after_ = static_cast<std::int32_t>(event_id);
    }
  }
  return true;
}

bool RdpScanner::update_event(
    std::uint32_t event_id,
    std::uint32_t recombinant,
    std::uint32_t major_parent,
    std::uint32_t minor_parent,
    std::size_t beginning,
      std::size_t ending,
      std::string& error) {
  if (!event_action_allowed(event_id, error)) return false;
  const std::size_t sequence_count = alignment_.sequence_count();
  if (recombinant >= sequence_count || major_parent >= sequence_count ||
      minor_parent >= sequence_count) {
    error = "An edited event refers to a sequence outside the alignment.";
    return false;
  }
  if (recombinant == major_parent || recombinant == minor_parent ||
      major_parent == minor_parent) {
    error = "Recombinant, major-parent-like, and minor-parent-like roles must be distinct.";
    return false;
  }
  if (beginning < 1 || ending < 1 || beginning > alignment_.length ||
      ending > alignment_.length) {
    error = "Edited breakpoints must fall inside the alignment.";
    return false;
  }
  rebuild_working_before_event(event_id);
  auto& event = events_[event_id];
  event.recombinant = recombinant;
  event.major_parent = major_parent;
  event.minor_parent = minor_parent;
  event.beginning = beginning;
  event.ending = ending;
  event.wraps_origin = options_.circular && beginning > ending;
  event.manual_adjusted = true;
  event.review_state = ReviewState::unreviewed;
  if (event.group_manual_adjusted) {
    event.co_recombinant_sequences.erase(
        std::remove_if(
            event.co_recombinant_sequences.begin(),
            event.co_recombinant_sequences.end(),
            [&](std::uint32_t sequence) {
              return sequence == major_parent || sequence == minor_parent;
            }),
        event.co_recombinant_sequences.end());
    event.co_recombinant_sequences.push_back(recombinant);
    sort_unique(event.co_recombinant_sequences);
  }
  build_event_evidence(event);
  if (event.anchor_signal_id < signals_.size() &&
      signals_[event.anchor_signal_id].correction_tests > 0) {
    correction_tests_ = signals_[event.anchor_signal_id].correction_tests;
  }
  refresh_trace_evidence(event);
  refresh_role_hypotheses(event);
  if (reconciliation_required_after_ < 0 ||
      event_id < static_cast<std::uint32_t>(reconciliation_required_after_)) {
    reconciliation_required_after_ = static_cast<std::int32_t>(event_id);
  }
  return true;
}

bool RdpScanner::update_event_group(
    std::uint32_t event_id,
    std::vector<std::uint32_t> sequences,
    bool manual_override,
    std::string& error) {
  if (!event_action_allowed(event_id, error)) return false;
  const auto& current = events_[event_id];
  if (manual_override) {
    if (std::any_of(sequences.begin(), sequences.end(), [&](std::uint32_t sequence) {
          return sequence >= alignment_.sequence_count();
        })) {
      error = "The edited co-recombinant group refers to a sequence outside the alignment.";
      return false;
    }
    if (std::find(sequences.begin(), sequences.end(), current.major_parent) != sequences.end() ||
        std::find(sequences.begin(), sequences.end(), current.minor_parent) != sequences.end()) {
      error = "The current parent representatives cannot be members of the co-recombinant group.";
      return false;
    }
    sequences.push_back(current.recombinant);
    sort_unique(sequences);
  }

  rebuild_working_before_event(event_id);
  auto& event = events_[event_id];
  event.co_recombinant_sequences = manual_override
      ? std::move(sequences)
      : event.automatic_co_recombinant_sequences;
  event.group_manual_adjusted = manual_override;
  event.manual_adjusted = true;
  event.review_state = ReviewState::unreviewed;
  if (reconciliation_required_after_ < 0 ||
      event_id < static_cast<std::uint32_t>(reconciliation_required_after_)) {
    reconciliation_required_after_ = static_cast<std::int32_t>(event_id);
  }
  return true;
}

bool RdpScanner::reconcile_after(std::uint32_t event_id, std::string& error) {
  if (event_id >= events_.size()) {
    error = "The selected reconciled event does not exist.";
    return false;
  }
  if (reconciliation_required_after_ < 0) {
    error = "No corrected or rejected event is waiting for downstream re-identification.";
    return false;
  }
  if (event_id != static_cast<std::uint32_t>(reconciliation_required_after_)) {
    error = "Re-identification must resume from the first changed event.";
    return false;
  }
  if (events_[event_id].review_state == ReviewState::unreviewed) {
    error = "Accept or reject the changed event before re-identifying later events.";
    return false;
  }

  std::vector<UniqueEvent> preserved(
      events_.begin(),
      events_.begin() + static_cast<std::ptrdiff_t>(event_id + 1));
  std::vector<Signal> retained_signals;
  for (std::size_t index = 0; index < preserved.size(); ++index) {
    auto& event = preserved[index];
    std::vector<std::uint32_t> remapped;
    std::uint32_t remapped_anchor = std::numeric_limits<std::uint32_t>::max();
    for (const std::uint32_t old_id : event.support_signal_ids) {
      if (old_id >= signals_.size()) continue;
      Signal signal = signals_[old_id];
      signal.id = static_cast<std::uint32_t>(retained_signals.size());
      signal.event_id = static_cast<std::int32_t>(index);
      if (old_id == event.anchor_signal_id) remapped_anchor = signal.id;
      remapped.push_back(signal.id);
      retained_signals.push_back(std::move(signal));
    }
    if (remapped.empty() || remapped_anchor == std::numeric_limits<std::uint32_t>::max()) {
      error = "A fixed event has lost its representative primary signal.";
      return false;
    }
    event.id = static_cast<std::uint32_t>(index);
    event.anchor_signal_id = remapped_anchor;
    event.support_signal_ids = std::move(remapped);
  }
  events_ = std::move(preserved);
  signals_ = std::move(retained_signals);

  reset_working_alignment();
  for (auto& event : events_) {
    const bool should_erase = event.review_state == ReviewState::accepted ||
        (event.review_state == ReviewState::unreviewed && event.tract_erased_for_detection);
    if (should_erase) {
      const ErasureResult erasure = erase_event_tract(event);
      event.erased_nucleotide_sites = erasure.original_sites;
      event.erased_working_sites = erasure.working_sites;
      event.fragment_sequences_added = erasure.fragments_added;
      event.tract_erased_for_detection = erasure.working_sites > 0;
    } else {
      event.erased_nucleotide_sites = 0;
      event.erased_working_sites = 0;
      event.fragment_sequences_added = 0;
      event.tract_erased_for_detection = false;
    }
  }

  fixed_event_count_ = events_.size();
  round_signal_begin_ = signals_.size();
  scan_round_ = events_.size() + 1;
  cumulative_triplets_ = 0;
  reset_round_cursor();
  running_ = true;
  primary_done_ = false;
  done_ = false;
  cancelled_.store(false);
  cycle_termination_ = "scanning";
  return true;
}

bool RdpScanner::restore(
    ScanOptions options,
    std::vector<Signal> signals,
    std::uint64_t correction_tests,
    std::string& error) {
  if (alignment_.sequence_count() < 3 || alignment_.length == 0) {
    error = "A valid saved alignment is required before restoring an analysis.";
    return false;
  }
  if (options.mask.size() != alignment_.sequence_count()) {
    error = "The saved sequence mask does not match the saved alignment.";
    return false;
  }
  if (!(options.p_value_cutoff > 0.0 && options.p_value_cutoff <= 1.0) ||
      options.window_sites < 5 ||
      std::count(options.mask.begin(), options.mask.end(), 0) < 3) {
    error = "The saved scan settings are outside the supported RDP range.";
    return false;
  }
  for (std::size_t index = 0; index < signals.size(); ++index) {
    auto& signal = signals[index];
    signal.id = static_cast<std::uint32_t>(index);
    signal.correction_tests = signal.correction_tests == 0
        ? std::min(correction_tests, kNativeCorrectionCap)
        : std::min(signal.correction_tests, kNativeCorrectionCap);
    for (const auto sequence : signal.triplet) {
      if (sequence >= alignment_.sequence_count()) {
        error = "A saved RDP signal refers to a sequence outside the alignment.";
        return false;
      }
    }
    if (signal.recombinant >= alignment_.sequence_count() ||
        signal.major_parent >= alignment_.sequence_count() ||
        signal.minor_parent >= alignment_.sequence_count() || signal.beginning < 1 ||
        signal.ending < 1 || signal.beginning > alignment_.length ||
        signal.ending > alignment_.length) {
      error = "A saved RDP signal contains invalid roles or breakpoints.";
      return false;
    }
  }
  options_ = std::move(options);
  signals_ = std::move(signals);
  reset_working_alignment();
  refresh_active_sequences();
  processed_triplets_ = total_triplets_;
  cumulative_triplets_ = processed_triplets_;
  correction_tests_ = correction_tests == 0
      ? std::min(total_triplets_, kNativeCorrectionCap)
      : std::min(correction_tests, kNativeCorrectionCap);
  running_ = false;
  primary_done_ = true;
  done_ = true;
  cancelled_.store(false);
  events_.clear();
  scan_round_ = 1;
  round_signal_begin_ = signals_.size();
  fixed_event_count_ = 0;
  cycle_termination_ = "restored-project";
  reconciliation_required_after_ = -1;
  return true;
}

bool RdpScanner::restore_event_state(
    std::uint32_t event_id,
    std::uint32_t anchor_signal_id,
    std::uint32_t recombinant,
    std::uint32_t major_parent,
    std::uint32_t minor_parent,
    std::size_t beginning,
    std::size_t ending,
    std::size_t detection_round,
    bool tract_erased_for_detection,
    ReviewState review_state,
    bool manual_adjusted,
    std::vector<std::uint32_t> co_recombinant_sequences,
    bool group_manual_adjusted,
    std::string& error) {
  if (event_id != events_.size() || anchor_signal_id >= signals_.size()) {
    error = "Saved events must be restored in their original analysis order.";
    return false;
  }
  if (recombinant >= alignment_.sequence_count() || major_parent >= alignment_.sequence_count() ||
      minor_parent >= alignment_.sequence_count() || recombinant == major_parent ||
      recombinant == minor_parent || major_parent == minor_parent || beginning < 1 ||
      ending < 1 || beginning > alignment_.length || ending > alignment_.length) {
    error = "A saved reconciled event contains invalid roles or breakpoints.";
    return false;
  }
  UniqueEvent event;
  event.id = event_id;
  event.anchor_signal_id = anchor_signal_id;
  event.recombinant = recombinant;
  event.major_parent = major_parent;
  event.minor_parent = minor_parent;
  event.beginning = beginning;
  event.ending = ending;
  event.wraps_origin = options_.circular && beginning > ending;
  event.review_state = review_state;
  event.manual_adjusted = manual_adjusted || group_manual_adjusted;
  event.group_manual_adjusted = group_manual_adjusted;
  event.detection_round = detection_round == 0 ? event_id + 1 : detection_round;
  if (event.group_manual_adjusted) {
    if (std::any_of(
            co_recombinant_sequences.begin(),
            co_recombinant_sequences.end(),
            [&](std::uint32_t sequence) {
              return sequence >= alignment_.sequence_count() || sequence == major_parent ||
                  sequence == minor_parent;
            })) {
      error = "A saved manual co-recombinant group is invalid for its event roles.";
      return false;
    }
    co_recombinant_sequences.push_back(recombinant);
    sort_unique(co_recombinant_sequences);
    event.co_recombinant_sequences = std::move(co_recombinant_sequences);
  }
  for (std::size_t signal_id = 0; signal_id < signals_.size(); ++signal_id) {
    if (signals_[signal_id].event_id == static_cast<std::int32_t>(event_id)) {
      event.support_signal_ids.push_back(static_cast<std::uint32_t>(signal_id));
    }
  }
  if (event.support_signal_ids.empty() ||
      std::find(
          event.support_signal_ids.begin(),
          event.support_signal_ids.end(),
          anchor_signal_id) == event.support_signal_ids.end()) {
    error = "A saved event is missing its grouped primary signal evidence.";
    return false;
  }
  build_event_evidence(event);
  if (signals_[anchor_signal_id].correction_tests > 0) {
    correction_tests_ = signals_[anchor_signal_id].correction_tests;
  }
  refresh_trace_evidence(event);
  refresh_role_hypotheses(event);
  event.fragment_assisted_detection = signals_[anchor_signal_id].fragment_assisted;
  event.tract_erased_for_detection = tract_erased_for_detection;
  if (event.tract_erased_for_detection) {
    const ErasureResult erasure = erase_event_tract(event);
    event.erased_nucleotide_sites = erasure.original_sites;
    event.erased_working_sites = erasure.working_sites;
    event.fragment_sequences_added = erasure.fragments_added;
  }
  events_.push_back(std::move(event));
  scan_round_ = std::max(scan_round_, detection_round + 1);
  cumulative_triplets_ = std::max(
      cumulative_triplets_,
      total_triplets_ * static_cast<std::uint64_t>(events_.size() + 1));
  return true;
}

bool RdpScanner::restore_reconciliation_requirement(
    std::uint32_t event_id,
    std::string& error) {
  if (event_id >= events_.size()) {
    error = "The saved downstream reconciliation marker refers to an unknown event.";
    return false;
  }
  reconciliation_required_after_ = static_cast<std::int32_t>(event_id);
  return true;
}

std::string RdpScanner::progress_json() const {
  double fraction = total_triplets_ == 0
      ? 0.0
      : static_cast<double>(processed_triplets_) / static_cast<double>(total_triplets_);
  if (primary_done_ && !done_) fraction = 0.985;
  if (done_) fraction = 1.0;
  const char* state = cancelled_.load()
      ? "cancelled"
      : done_ ? "done" : (running_ || primary_done_) ? "running" : "idle";
  const char* phase = done_
      ? "complete"
      : primary_done_ ? "reconciliation" : scan_round_ > 1 ? "cyclic-rescan" : "primary";
  std::ostringstream out;
  out << "{\"state\":\"" << state << "\",\"phase\":\"" << phase
      << "\",\"processedTriplets\":"
      << processed_triplets_ << ",\"totalTriplets\":" << total_triplets_
      << ",\"cumulativeTriplets\":" << cumulative_triplets_
      << ",\"scanRound\":" << scan_round_
      << ",\"fixedEventCount\":" << fixed_event_count_
      << ",\"signalCount\":" << signals_.size() << ",\"eventCount\":" << events_.size()
      << ",\"cycleTermination\":\"" << cycle_termination_ << "\",\"fraction\":";
  json::number(out, std::clamp(fraction, 0.0, 1.0));
  out << '}';
  return out.str();
}

std::string RdpScanner::results_json() const {
  const bool final_alignment_ready = done_ && reconciliation_required_after_ < 0 &&
      std::none_of(events_.begin(), events_.end(), [](const UniqueEvent& event) {
        return event.review_state == ReviewState::unreviewed;
      });
  std::ostringstream out;
  out << "{\"engineVersion\":\"0.5.0-session-5\",\"status\":\"cyclic-three-set-reconciled\","
         "\"method\":\"RDP\",\"reconciliationTier\":\"detectable-distance-phylogenetic\","
         "\"cycleMode\":\"strongest-first-tract-erasure-with-bounded-fragment-reentry\","
         "\"finalAlignmentReady\":"
      << (final_alignment_ready ? "true" : "false") << ",\"fragmentReentry\":"
      << (working_alignment_.length < kFragmentReentryAlignmentLengthLimit ? "true" : "false")
      << ",\"fragmentReentryAlignmentLengthLimit\":"
      << kFragmentReentryAlignmentLengthLimit << ",\"fragmentSequenceCap\":"
      << kWorkingFragmentSequenceCap << ",\"fragmentReentryCapped\":"
      << (fragment_reentry_capped_ ? "true" : "false")
      << ",\"workingSequenceCount\":" << working_alignment_.sequence_count()
      << ",\"workingFragmentSequenceCount\":"
      << (working_alignment_.sequence_count() - alignment_.sequence_count())
      << ",\"scanRounds\":"
      << scan_round_ << ",\"cumulativeTriplets\":" << cumulative_triplets_
      << ",\"cycleTermination\":\"" << cycle_termination_ << "\",\"correction\":\""
      << (options_.correction == CorrectionMode::bonferroni ? "bonferroni" : "none")
      << "\",\"correctionTests\":" << correction_tests_
      << ",\"circular\":" << (options_.circular ? "true" : "false")
      << ",\"maskedSequenceIndices\":[";
  bool wrote_mask = false;
  for (std::size_t index = 0; index < options_.mask.size(); ++index) {
    if (options_.mask[index] == 0) continue;
    if (wrote_mask) out << ',';
    wrote_mask = true;
    out << index;
  }
  out << "],\"downstreamReconciliationRequiredAfter\":";
  if (reconciliation_required_after_ < 0) out << "null";
  else out << reconciliation_required_after_;
  out << ",\"pValueCutoff\":";
  json::number(out, options_.p_value_cutoff);
  out << ",\"windowSites\":" << options_.window_sites << ",\"signals\":[";
  for (std::size_t index = 0; index < signals_.size(); ++index) {
    if (index) out << ',';
    const auto& signal = signals_[index];
    out << "{\"id\":" << signal.id << ",\"method\":\"RDP\",\"triplet\":["
        << signal.triplet[0] << ',' << signal.triplet[1] << ',' << signal.triplet[2]
        << "],\"tripletNames\":[";
    for (std::size_t name = 0; name < 3; ++name) {
      if (name) out << ',';
      json::string(out, alignment_.names[signal.triplet[name]]);
    }
    out << "],\"recombinant\":" << signal.recombinant << ",\"recombinantName\":";
    json::string(out, alignment_.names[signal.recombinant]);
    out << ",\"majorParent\":" << signal.major_parent << ",\"majorParentName\":";
    json::string(out, alignment_.names[signal.major_parent]);
    out << ",\"minorParent\":" << signal.minor_parent << ",\"minorParentName\":";
    json::string(out, alignment_.names[signal.minor_parent]);
    out << ",\"beginning\":" << signal.beginning << ",\"ending\":" << signal.ending
        << ",\"wrapsOrigin\":" << (signal.wraps_origin ? "true" : "false")
        << ",\"informativeBeginning\":" << signal.informative_beginning
        << ",\"informativeEnding\":" << signal.informative_ending
        << ",\"localPValue\":";
    json::number(out, signal.local_p_value);
    out << ",\"correctedPValue\":";
    json::number(out, signal.corrected_p_value);
    out << ",\"correctionTests\":" << signal.correction_tests
        << ",\"pairSimilarity\":[";
    for (std::size_t pair = 0; pair < 3; ++pair) {
      if (pair) out << ',';
      json::number(out, signal.pair_similarity[pair]);
    }
    out << "],\"informativeSites\":" << signal.informative_sites
        << ",\"candidatePair\":" << static_cast<unsigned int>(signal.candidate_pair)
        << ",\"fragmentAssisted\":" << (signal.fragment_assisted ? "true" : "false")
        << ",\"fragmentEventContext\":[";
    for (std::size_t member = 0; member < signal.fragment_event_context.size(); ++member) {
      if (member) out << ',';
      if (signal.fragment_event_context[member] < 0) out << "null";
      else out << signal.fragment_event_context[member];
    }
    out
        << "],\"eventId\":";
    if (signal.event_id < 0) out << "null";
    else out << signal.event_id;
    out
        << ",\"reviewState\":\"" << review_state_name(signal.review_state)
        << "\",\"provisionalRoles\":true}";
  }
  out << "],\"events\":[";
  for (std::size_t index = 0; index < events_.size(); ++index) {
    if (index) out << ',';
    const auto& event = events_[index];
    out << "{\"id\":" << event.id << ",\"anchorSignalId\":" << event.anchor_signal_id
        << ",\"detectionRound\":" << event.detection_round
        << ",\"erasedNucleotideSites\":" << event.erased_nucleotide_sites
        << ",\"erasedWorkingSites\":" << event.erased_working_sites
        << ",\"fragmentSequencesAdded\":" << event.fragment_sequences_added
        << ",\"fragmentAssistedDetection\":"
        << (event.fragment_assisted_detection ? "true" : "false")
        << ",\"tractErasedForDetection\":"
        << (event.tract_erased_for_detection ? "true" : "false")
        << ",\"reconciliationBasis\":\"two-shared-sequences-and-30-percent-overlap\""
        << ",\"recombinant\":" << event.recombinant << ",\"recombinantName\":";
    json::string(out, alignment_.names[event.recombinant]);
    out << ",\"majorParent\":" << event.major_parent << ",\"majorParentName\":";
    json::string(out, alignment_.names[event.major_parent]);
    out << ",\"minorParent\":" << event.minor_parent << ",\"minorParentName\":";
    json::string(out, alignment_.names[event.minor_parent]);
    out << ",\"beginning\":" << event.beginning << ",\"ending\":" << event.ending
        << ",\"wrapsOrigin\":" << (event.wraps_origin ? "true" : "false")
        << ",\"bestLocalPValue\":";
    json::number(out, event.best_local_p_value);
    out << ",\"bestCorrectedPValue\":";
    json::number(out, event.best_corrected_p_value);
    out << ",\"supportSignalIds\":[";
    for (std::size_t support = 0; support < event.support_signal_ids.size(); ++support) {
      if (support) out << ',';
      out << event.support_signal_ids[support];
    }
    out << "],\"detectableSequenceIndices\":[";
    for (std::size_t sequence = 0; sequence < event.detectable_sequences.size(); ++sequence) {
      if (sequence) out << ',';
      out << event.detectable_sequences[sequence];
    }
    out << "],\"detectableSequenceNames\":[";
    for (std::size_t sequence = 0; sequence < event.detectable_sequences.size(); ++sequence) {
      if (sequence) out << ',';
      json::string(out, alignment_.names[event.detectable_sequences[sequence]]);
    }
    out << "],\"roleCandidateIndices\":{\"recombinant\":[";
    for (std::size_t role = 0; role < event.role_candidates[0].size(); ++role) {
      if (role) out << ',';
      out << event.role_candidates[0][role];
    }
    out << "],\"majorParent\":[";
    for (std::size_t role = 0; role < event.role_candidates[1].size(); ++role) {
      if (role) out << ',';
      out << event.role_candidates[1][role];
    }
    out << "],\"minorParent\":[";
    for (std::size_t role = 0; role < event.role_candidates[2].size(); ++role) {
      if (role) out << ',';
      out << event.role_candidates[2][role];
    }
    out << "]},\"automaticCoRecombinantSequenceIndices\":[";
    for (std::size_t sequence = 0;
         sequence < event.automatic_co_recombinant_sequences.size();
         ++sequence) {
      if (sequence) out << ',';
      out << event.automatic_co_recombinant_sequences[sequence];
    }
    out << "],\"automaticCoRecombinantSequenceNames\":[";
    for (std::size_t sequence = 0;
         sequence < event.automatic_co_recombinant_sequences.size();
         ++sequence) {
      if (sequence) out << ',';
      json::string(
          out, alignment_.names[event.automatic_co_recombinant_sequences[sequence]]);
    }
    out << "],\"coRecombinantSequenceIndices\":[";
    for (std::size_t sequence = 0; sequence < event.co_recombinant_sequences.size(); ++sequence) {
      if (sequence) out << ',';
      out << event.co_recombinant_sequences[sequence];
    }
    out << "],\"coRecombinantSequenceNames\":[";
    for (std::size_t sequence = 0; sequence < event.co_recombinant_sequences.size(); ++sequence) {
      if (sequence) out << ',';
      json::string(out, alignment_.names[event.co_recombinant_sequences[sequence]]);
    }
    out << "],\"treePanel\":{\"sequenceCount\":" << event.tree_panel_sequences
        << ",\"subsampled\":" << (event.tree_panel_subsampled ? "true" : "false")
        << ",\"sequenceCap\":" << kEventTreeSequenceCap << ",\"regions\":[";
    constexpr std::array<const char*, 6> tree_region_names{
        "5-prime-outside",
        "5-prime-inside",
        "3-prime-outside",
        "3-prime-inside",
        "outside-tract",
        "inside-tract",
    };
    for (std::size_t region = 0; region < event.tree_regions.size(); ++region) {
      if (region) out << ',';
      const auto& summary = event.tree_regions[region];
      out << "{\"name\":\"" << tree_region_names[region] << "\",\"sites\":"
          << summary.sites << ",\"sequences\":" << summary.sequences
          << ",\"bootstrapReplicates\":" << summary.bootstrap_replicates
          << ",\"supportedInternalBranches\":" << summary.supported_internal_branches
          << ",\"internalBranches\":" << summary.internal_branches
          << ",\"usable\":" << (summary.usable ? "true" : "false") << '}';
    }
    out << "]},\"roleConsensus\":{\"method\":\"source-decision-tree-subset\","
           "\"nativeWeightParity\":false,\"informative\":"
        << (event.role_consensus.informative ? "true" : "false")
        << ",\"recommendedRole\":";
    if (event.role_consensus.recommended_role < 0) out << "null";
    else out << static_cast<unsigned int>(event.role_consensus.recommended_role);
    out << ",\"recommendedRecombinant\":"
        << event.role_consensus.recommended_recombinant
        << ",\"recommendedRecombinantName\":";
    json::string(out, alignment_.names[event.role_consensus.recommended_recombinant]);
    out << ",\"recommendedMajorParent\":"
        << event.role_consensus.recommended_major_parent
        << ",\"recommendedMajorParentName\":";
    json::string(out, alignment_.names[event.role_consensus.recommended_major_parent]);
    out << ",\"recommendedMinorParent\":"
        << event.role_consensus.recommended_minor_parent
        << ",\"recommendedMinorParentName\":";
    json::string(out, alignment_.names[event.role_consensus.recommended_minor_parent]);
    out << ",\"confidence\":";
    json::number(out, event.role_consensus.confidence);
    out << ",\"votes\":[";
    for (std::size_t role = 0; role < event.role_consensus.votes.size(); ++role) {
      if (role) out << ',';
      json::number(out, event.role_consensus.votes[role]);
    }
    out << "],\"metrics\":[";
    for (std::size_t metric_index = 0;
         metric_index < event.role_consensus.metrics.size();
         ++metric_index) {
      if (metric_index) out << ',';
      const auto& metric = event.role_consensus.metrics[metric_index];
      out << "{\"method\":\"" << role_metric_name(metric.kind) << "\",\"scores\":[";
      for (std::size_t role = 0; role < metric.scores.size(); ++role) {
        if (role) out << ',';
        json::number(out, metric.scores[role]);
      }
      out << "],\"contributions\":[";
      for (std::size_t role = 0; role < metric.contributions.size(); ++role) {
        if (role) out << ',';
        json::number(out, metric.contributions[role]);
      }
      out << "],\"weight\":";
      json::number(out, metric.weight);
      out << ",\"winningRole\":";
      if (metric.winning_role < 0) out << "null";
      else out << static_cast<unsigned int>(metric.winning_role);
      out << ",\"higherIsRecombinant\":"
          << (metric.higher_is_recombinant ? "true" : "false")
          << ",\"informative\":" << (metric.informative ? "true" : "false") << '}';
    }
    out << "]},\"roleHypotheses\":[";
    for (std::size_t role = 0; role < event.role_hypotheses.size(); ++role) {
      if (role) out << ',';
      const auto& hypothesis = event.role_hypotheses[role];
      out << "{\"presumedRecombinant\":" << hypothesis.presumed_recombinant
          << ",\"presumedRecombinantName\":";
      json::string(out, alignment_.names[hypothesis.presumed_recombinant]);
      out << ",\"parentOne\":" << hypothesis.parent_one << ",\"parentOneName\":";
      json::string(out, alignment_.names[hypothesis.parent_one]);
      out << ",\"parentTwo\":" << hypothesis.parent_two << ",\"parentTwoName\":";
      json::string(out, alignment_.names[hypothesis.parent_two]);
      out << ",\"testedSequences\":" << hypothesis.tested_sequences
          << ",\"validSequences\":" << hypothesis.valid_sequences
          << ",\"detectableSignalSetIndices\":[";
      for (std::size_t sequence = 0; sequence < hypothesis.detectable_signal_set.size(); ++sequence) {
        if (sequence) out << ',';
        out << hypothesis.detectable_signal_set[sequence];
      }
      out << "],\"detectableSignalSetNames\":[";
      for (std::size_t sequence = 0; sequence < hypothesis.detectable_signal_set.size(); ++sequence) {
        if (sequence) out << ',';
        json::string(out, alignment_.names[hypothesis.detectable_signal_set[sequence]]);
      }
      out << "],\"distanceCorrelationSetIndices\":[";
      for (std::size_t sequence = 0; sequence < hypothesis.distance_correlation_set.size(); ++sequence) {
        if (sequence) out << ',';
        out << hypothesis.distance_correlation_set[sequence];
      }
      out << "],\"distanceCorrelationSetNames\":[";
      for (std::size_t sequence = 0; sequence < hypothesis.distance_correlation_set.size(); ++sequence) {
        if (sequence) out << ',';
        json::string(out, alignment_.names[hypothesis.distance_correlation_set[sequence]]);
      }
      out << "],\"phylogeneticCorrelationSetIndices\":[";
      for (std::size_t sequence = 0;
           sequence < hypothesis.phylogenetic_correlation_set.size();
           ++sequence) {
        if (sequence) out << ',';
        out << hypothesis.phylogenetic_correlation_set[sequence];
      }
      out << "],\"phylogeneticCorrelationSetNames\":[";
      for (std::size_t sequence = 0;
           sequence < hypothesis.phylogenetic_correlation_set.size();
           ++sequence) {
        if (sequence) out << ',';
        json::string(out, alignment_.names[hypothesis.phylogenetic_correlation_set[sequence]]);
      }
      out << "],\"completeTwoOfThreeSetIndices\":[";
      for (std::size_t sequence = 0;
           sequence < hypothesis.complete_two_of_three_set.size();
           ++sequence) {
        if (sequence) out << ',';
        out << hypothesis.complete_two_of_three_set[sequence];
      }
      out << "],\"completeTwoOfThreeSetNames\":[";
      for (std::size_t sequence = 0;
           sequence < hypothesis.complete_two_of_three_set.size();
           ++sequence) {
        if (sequence) out << ',';
        json::string(out, alignment_.names[hypothesis.complete_two_of_three_set[sequence]]);
      }
      out << "],\"correlationWarnings\":[";
      for (std::size_t pair = 0; pair < hypothesis.correlation_warnings.size(); ++pair) {
        if (pair) out << ',';
        out << (hypothesis.correlation_warnings[pair] ? "true" : "false");
      }
      out << "],\"distanceCorrelationEvidence\":[";
      for (std::size_t evidence_index = 0;
           evidence_index < hypothesis.distance_evidence.size();
           ++evidence_index) {
        if (evidence_index) out << ',';
        const auto& evidence = hypothesis.distance_evidence[evidence_index];
        out << "{\"sequenceIndex\":" << evidence.sequence << ",\"sequenceName\":";
        json::string(out, alignment_.names[evidence.sequence]);
        out << ",\"correlations\":[";
        for (std::size_t pair = 0; pair < evidence.correlations.size(); ++pair) {
          if (pair) out << ',';
          json::number(out, evidence.correlations[pair]);
        }
        out << "],\"directCorrelations\":[";
        for (std::size_t pair = 0; pair < evidence.direct_correlations.size(); ++pair) {
          if (pair) out << ',';
          json::number(out, evidence.direct_correlations[pair]);
        }
        out << "],\"pValues\":[";
        for (std::size_t pair = 0; pair < evidence.p_values.size(); ++pair) {
          if (pair) out << ',';
          json::number(out, evidence.p_values[pair]);
        }
        out << "],\"inversionCodes\":[";
        for (std::size_t pair = 0; pair < evidence.inversion_codes.size(); ++pair) {
          if (pair) out << ',';
          out << static_cast<unsigned int>(evidence.inversion_codes[pair]);
        }
        out << "],\"warningFiltered\":[";
        for (std::size_t pair = 0; pair < evidence.warning_filtered.size(); ++pair) {
          if (pair) out << ',';
          out << (evidence.warning_filtered[pair] ? "true" : "false");
        }
        out << "],\"duplicateFiltered\":[";
        for (std::size_t pair = 0; pair < evidence.duplicate_filtered.size(); ++pair) {
          if (pair) out << ',';
          out << (evidence.duplicate_filtered[pair] ? "true" : "false");
        }
        out << "],\"minimumComparableSites\":[";
        for (std::size_t pair = 0; pair < evidence.minimum_comparable_sites.size(); ++pair) {
          if (pair) out << ',';
          out << evidence.minimum_comparable_sites[pair];
        }
        out << "],\"breakpointOverlapSites\":["
            << evidence.breakpoint_overlap_sites[0] << ','
            << evidence.breakpoint_overlap_sites[1]
            << "],\"aggregateScore\":";
        json::number(out, evidence.aggregate_score);
        out << ",\"aggregateTarget\":";
        json::number(out, evidence.aggregate_target);
        out
            << ",\"overlapEligible\":"
            << (evidence.overlap_eligible ? "true" : "false")
            << ",\"acceptableAffinity\":"
            << (evidence.acceptable_affinity ? "true" : "false")
            << ",\"strongCorrelationOverride\":"
            << (evidence.strong_correlation_override ? "true" : "false")
            << ",\"bestMatrixPair\":";
        if (evidence.best_matrix_pair < 0) out << "null";
        else out << static_cast<unsigned int>(evidence.best_matrix_pair);
        out << ",\"significant\":" << (evidence.significant ? "true" : "false")
            << ",\"detectableSupport\":"
            << (evidence.detectable_support ? "true" : "false")
            << ",\"positiveSupport\":"
            << (evidence.positive_support ? "true" : "false")
            << ",\"inverseSupport\":"
            << (evidence.inverse_support ? "true" : "false")
            << ",\"strippedInverseOnly\":"
            << (evidence.stripped_inverse_only ? "true" : "false")
            << ",\"duplicateCleanedSupport\":"
            << (evidence.duplicate_cleaned_support ? "true" : "false") << '}';
      }
      out << "],\"phylogeneticCorrelationEvidence\":[";
      for (std::size_t evidence_index = 0;
           evidence_index < hypothesis.phylogenetic_evidence.size();
           ++evidence_index) {
        if (evidence_index) out << ',';
        const auto& evidence = hypothesis.phylogenetic_evidence[evidence_index];
        out << "{\"sequenceIndex\":" << evidence.sequence << ",\"sequenceName\":";
        json::string(out, alignment_.names[evidence.sequence]);
        out << ",\"collapsedAffinityMargins\":[";
        for (std::size_t pair = 0; pair < evidence.collapsed_affinity_margin.size(); ++pair) {
          if (pair) out << ',';
          json::number(out, evidence.collapsed_affinity_margin[pair]);
        }
        out << "],\"rawAffinityMargins\":[";
        for (std::size_t pair = 0; pair < evidence.raw_affinity_margin.size(); ++pair) {
          if (pair) out << ',';
          json::number(out, evidence.raw_affinity_margin[pair]);
        }
        out << "],\"collapsedPairSupport\":[";
        for (std::size_t pair = 0; pair < evidence.collapsed_pair_support.size(); ++pair) {
          if (pair) out << ',';
          out << (evidence.collapsed_pair_support[pair] ? "true" : "false");
        }
        out << "],\"rawPairSupport\":[";
        for (std::size_t pair = 0; pair < evidence.raw_pair_support.size(); ++pair) {
          if (pair) out << ',';
          out << (evidence.raw_pair_support[pair] ? "true" : "false");
        }
        out << "],\"bestTreePair\":";
        if (evidence.best_tree_pair < 0) out << "null";
        else out << static_cast<unsigned int>(evidence.best_tree_pair);
        out << ",\"supportingTreePairs\":"
            << static_cast<unsigned int>(evidence.supporting_tree_pairs)
            << ",\"included\":" << (evidence.included ? "true" : "false")
            << ",\"distanceFallback\":"
            << (evidence.distance_fallback ? "true" : "false")
            << ",\"maskedExcluded\":"
            << (evidence.masked_excluded ? "true" : "false") << '}';
      }
      out << "],\"phylogeneticCorrelationStatus\":\"complete\","
             "\"evidenceSetConsensusComplete\":true,"
             "\"finalTrimDuplicateCorrelationStatus\":\"complete\","
             "\"lateNativeConsensusComplete\":false}";
    }
    out << "],\"traceEvidence\":[";
    for (std::size_t trace = 0; trace < event.trace_evidence.size(); ++trace) {
      if (trace) out << ',';
      const auto& evidence = event.trace_evidence[trace];
      out << "{\"sequenceIndex\":" << evidence.sequence << ",\"sequenceName\":";
      json::string(out, alignment_.names[evidence.sequence]);
      out << ",\"beginning\":" << evidence.beginning << ",\"ending\":" << evidence.ending
          << ",\"wrapsOrigin\":" << (evidence.wraps_origin ? "true" : "false")
          << ",\"localPValue\":";
      json::number(out, evidence.local_p_value);
      out << ",\"correctedPValue\":";
      json::number(out, evidence.corrected_p_value);
      out << ",\"significant\":" << (evidence.significant ? "true" : "false") << '}';
    }
    out << "],\"reviewState\":\"" << review_state_name(event.review_state)
        << "\",\"manualAdjusted\":" << (event.manual_adjusted ? "true" : "false")
        << ",\"groupManualAdjusted\":"
        << (event.group_manual_adjusted ? "true" : "false")
        << ",\"rolesProvisional\":true}";
  }
  out << "],\"notes\":["
         "\"Events are selected cyclically from the strongest unexplained RDP signal; each inferred co-recombinant tract is erased before the next complete triplet screen.\","
         "\"Signal grouping implements the supplied RDP5 detectable-signal rule: two shared triplet sequences and greater than 30% symmetric tract overlap.\","
         "\"Each anchor sequence is treated in turn as the presumed recombinant; three paired six-value correlations use the supplied direct and five category-relabelled Pearson paths.\","
         "\"Six Jukes-Cantor neighbour-joining trees are bootstrapped ten times, branches below 50 percent support are collapsed, and sequences in at least two of the three evidence sets form the co-recombinant group.\","
         "\"Masked-sequence trace checks retain structurally matching RDP profiles even when their corrected p-values are not significant.\","
         "\"Role identification ports MakePhPrScore, leave-one-role-out scores, displacement scores, weighted MakeTrpScore ordering changes, and the corresponding supplied decision-tree contributions.\","
         "\"Manual co-recombinant group edits are preserved separately from the automatic two-of-three set and drive subsequent tract erasure and accepted-event alignment exports.\","
         "\"Below the supplied 100000-site cutoff, erased tracts re-enter subsequent rounds as gap-padded synthetic fragments with original-sequence provenance; same-origin working copies never occupy one triplet and the retained-fragment cap is explicit.\","
         "\"Distance membership applies the supplied MakeACOR topology-affinity gate, MakeRList dual-correlation override, StripDupInv inverse-only removal, and the first FinalTrim duplicate-correlation cleanup; remaining FinalTrim and ConsensusOK decisions are documented parity boundaries.\"]}";
  return out.str();
}

bool RdpScanner::set_review_state(std::uint32_t signal_id, ReviewState state) {
  if (signal_id >= signals_.size()) return false;
  signals_[signal_id].review_state = state;
  return true;
}

SignalPlot RdpScanner::signal_plot(std::uint32_t signal_id, std::string& error) const {
  SignalPlot plot;
  if (signal_id >= signals_.size()) {
    error = "The selected RDP signal does not exist.";
    return plot;
  }
  const auto& signal = signals_[signal_id];
  TripletProfile profile;
  if (!build_profile_on(alignment_, signal.triplet, profile)) {
    error = "The selected triplet no longer passes the RDP informative-site filters.";
    return plot;
  }
  plot.signal_id = signal_id;
  plot.window_sites = options_.window_sites;
  const std::size_t effective_window = options_.window_sites / 2 * 2 + 1;
  const std::size_t stride = std::max<std::size_t>(1, profile.category.size() / 2000);
  for (std::size_t position = 0; position < profile.category.size(); position += stride) {
    PlotPoint point;
    point.alignment_position = profile.coordinates[position];
    for (std::size_t pair = 0; pair < 3; ++pair) {
      point.pair_identity[pair] = static_cast<double>(profile.rolling_counts[pair][position]) /
          static_cast<double>(effective_window);
    }
    plot.points.push_back(point);
  }
  if (!profile.coordinates.empty() &&
      (plot.points.empty() || plot.points.back().alignment_position != profile.coordinates.back())) {
    PlotPoint point;
    point.alignment_position = profile.coordinates.back();
    for (std::size_t pair = 0; pair < 3; ++pair) {
      point.pair_identity[pair] = static_cast<double>(profile.rolling_counts[pair].back()) /
          static_cast<double>(effective_window);
    }
    plot.points.push_back(point);
  }
  return plot;
}

std::string RdpScanner::csv() const {
  std::ostringstream out;
  out << "Event,Method,Detection round,Erased original nucleotide sites,Erased working sites,"
         "Fragments added,Fragment-assisted detection,Tract applied during detection,"
         "Recombinant,Major parent,Minor parent,Beginning,Ending,Wraps origin,"
         "Best local p-value,Best corrected p-value,Supporting signals,Detectable sequences,"
         "Distance-correlation sequences,FinalTrim duplicate-filtered pairs,"
         "Phylogenetic-correlation sequences,"
         "Automatic two-of-three group,Current co-recombinant group,Group manually adjusted,"
         "Masked trace sequences,"
         "Correlation sequences tested,Tree panel sequences,Tree panel subsampled,"
         "Recommended recombinant,Recommended major parent,Recommended minor parent,"
         "Role confidence,Weighted role scores,Review state,Manually adjusted,Native full-consensus parity\n";
  for (const auto& event : events_) {
    const auto& hypothesis = event.role_hypotheses[0];
    std::ostringstream support;
    for (std::size_t index = 0; index < event.support_signal_ids.size(); ++index) {
      if (index) support << ';';
      support << event.support_signal_ids[index] + 1;
    }
    std::ostringstream detectable;
    for (std::size_t index = 0; index < hypothesis.detectable_signal_set.size(); ++index) {
      if (index) detectable << ';';
      detectable << alignment_.names[hypothesis.detectable_signal_set[index]];
    }
    std::ostringstream correlated;
    for (std::size_t index = 0; index < hypothesis.distance_correlation_set.size(); ++index) {
      if (index) correlated << ';';
      correlated << alignment_.names[hypothesis.distance_correlation_set[index]];
    }
    std::ostringstream duplicate_filtered;
    bool first_duplicate = true;
    constexpr std::array<const char*, 3> duplicate_pair_names{
        "5-prime", "3-prime", "tract-outside"};
    for (const auto& evidence : hypothesis.distance_evidence) {
      for (std::size_t pair = 0; pair < evidence.duplicate_filtered.size(); ++pair) {
        if (evidence.duplicate_filtered[pair] == 0) continue;
        if (!first_duplicate) duplicate_filtered << ';';
        first_duplicate = false;
        duplicate_filtered << alignment_.names[evidence.sequence] << ':'
                           << duplicate_pair_names[pair];
      }
    }
    std::ostringstream automatic_consensus;
    for (std::size_t index = 0;
         index < event.automatic_co_recombinant_sequences.size();
         ++index) {
      if (index) automatic_consensus << ';';
      automatic_consensus << alignment_.names[event.automatic_co_recombinant_sequences[index]];
    }
    std::ostringstream consensus;
    for (std::size_t index = 0; index < event.co_recombinant_sequences.size(); ++index) {
      if (index) consensus << ';';
      consensus << alignment_.names[event.co_recombinant_sequences[index]];
    }
    std::ostringstream phylogenetic;
    for (std::size_t index = 0;
         index < hypothesis.phylogenetic_correlation_set.size();
         ++index) {
      if (index) phylogenetic << ';';
      phylogenetic << alignment_.names[hypothesis.phylogenetic_correlation_set[index]];
    }
    std::ostringstream traces;
    for (std::size_t index = 0; index < event.trace_evidence.size(); ++index) {
      if (index) traces << ';';
      traces << alignment_.names[event.trace_evidence[index].sequence];
      if (!event.trace_evidence[index].significant) traces << " (trace)";
    }
    out << event.id + 1 << ",RDP," << event.detection_round << ','
        << event.erased_nucleotide_sites << ','
        << event.erased_working_sites << ',' << event.fragment_sequences_added << ','
        << (event.fragment_assisted_detection ? "yes" : "no") << ','
        << (event.tract_erased_for_detection ? "yes" : "no") << ','
        << csv_cell(alignment_.names[event.recombinant]) << ','
        << csv_cell(alignment_.names[event.major_parent]) << ','
        << csv_cell(alignment_.names[event.minor_parent]) << ',' << event.beginning << ','
        << event.ending << ',' << (event.wraps_origin ? "yes" : "no") << ','
        << std::setprecision(16) << event.best_local_p_value << ','
        << event.best_corrected_p_value << ',' << csv_cell(support.str()) << ','
        << csv_cell(detectable.str()) << ',' << csv_cell(correlated.str()) << ','
        << csv_cell(duplicate_filtered.str()) << ',' << csv_cell(phylogenetic.str()) << ','
        << csv_cell(automatic_consensus.str()) << ','
        << csv_cell(consensus.str()) << ','
        << (event.group_manual_adjusted ? "yes" : "no") << ','
        << csv_cell(traces.str()) << ',' << hypothesis.tested_sequences << ','
        << event.tree_panel_sequences << ',' << (event.tree_panel_subsampled ? "yes" : "no")
        << ',' << csv_cell(alignment_.names[event.role_consensus.recommended_recombinant])
        << ',' << csv_cell(alignment_.names[event.role_consensus.recommended_major_parent])
        << ',' << csv_cell(alignment_.names[event.role_consensus.recommended_minor_parent])
        << ',' << event.role_consensus.confidence << ',';
    std::ostringstream votes;
    for (std::size_t role = 0; role < event.role_consensus.votes.size(); ++role) {
      if (role) votes << ';';
      votes << event.role_consensus.votes[role];
    }
    out << csv_cell(votes.str()) << ',' << review_state_name(event.review_state) << ','
        << (event.manual_adjusted ? "yes" : "no") << ",no\n";
  }
  return out.str();
}

bool RdpScanner::final_alignment_ready(std::string& error) const {
  if (!done_) {
    error = "Finish the cyclic RDP scan before exporting a final alignment.";
    return false;
  }
  if (reconciliation_required_after_ >= 0) {
    error = "Re-identify the downstream event chain before exporting a final alignment.";
    return false;
  }
  if (std::any_of(events_.begin(), events_.end(), [](const UniqueEvent& event) {
        return event.review_state == ReviewState::unreviewed;
      })) {
    error = "Accept or reject every event before exporting a final alignment.";
    return false;
  }
  return true;
}

std::string RdpScanner::recombination_free_fasta() const {
  std::vector<std::string> sequences = alignment_.sequences;
  for (const auto& event : events_) {
    if (event.review_state != ReviewState::accepted) continue;
    std::vector<std::uint32_t> affected = event.co_recombinant_sequences;
    affected.push_back(event.recombinant);
    sort_unique(affected);
    for (const std::uint32_t sequence : affected) {
      if (sequence >= sequences.size()) continue;
      for (std::size_t coordinate = 1; coordinate <= alignment_.length; ++coordinate) {
        if (coordinate_in_tract(
                coordinate, event.beginning, event.ending, event.wraps_origin)) {
          sequences[sequence][coordinate - 1] = '-';
        }
      }
    }
  }

  std::ostringstream out;
  for (std::size_t sequence = 0; sequence < sequences.size(); ++sequence) {
    write_fasta_record(out, alignment_.names[sequence], sequences[sequence]);
  }
  return out.str();
}

std::string RdpScanner::fragmented_fasta() const {
  struct FragmentRecord {
    std::string name;
    std::string sequence;
  };
  std::vector<std::string> remainder = alignment_.sequences;
  std::vector<FragmentRecord> fragments;
  for (const auto& event : events_) {
    if (event.review_state != ReviewState::accepted) continue;
    std::vector<std::uint32_t> affected = event.co_recombinant_sequences;
    affected.push_back(event.recombinant);
    sort_unique(affected);
    for (const std::uint32_t sequence : affected) {
      if (sequence >= remainder.size()) continue;
      std::string fragment(alignment_.length, '-');
      bool copied = false;
      for (std::size_t coordinate = 1; coordinate <= alignment_.length; ++coordinate) {
        if (!coordinate_in_tract(
                coordinate, event.beginning, event.ending, event.wraps_origin)) {
          continue;
        }
        fragment[coordinate - 1] = remainder[sequence][coordinate - 1];
        remainder[sequence][coordinate - 1] = '-';
        copied = true;
      }
      if (!copied) continue;
      std::ostringstream name;
      name << alignment_.names[sequence] << "|RDP_event_" << event.id + 1 << "|bp_"
           << event.beginning << '_' << event.ending;
      if (event.wraps_origin) name << "|circular";
      fragments.push_back({name.str(), std::move(fragment)});
    }
  }

  std::ostringstream out;
  for (std::size_t sequence = 0; sequence < remainder.size(); ++sequence) {
    write_fasta_record(out, alignment_.names[sequence], remainder[sequence]);
  }
  for (const auto& fragment : fragments) {
    write_fasta_record(out, fragment.name, fragment.sequence);
  }
  return out.str();
}

std::string signal_plot_json(const SignalPlot& plot) {
  std::ostringstream out;
  out << "{\"signalId\":" << plot.signal_id << ",\"windowSites\":" << plot.window_sites
      << ",\"points\":[";
  for (std::size_t index = 0; index < plot.points.size(); ++index) {
    if (index) out << ',';
    const auto& point = plot.points[index];
    out << "{\"alignmentPosition\":" << point.alignment_position << ",\"pair12\":";
    json::number(out, point.pair_identity[0]);
    out << ",\"pair13\":";
    json::number(out, point.pair_identity[1]);
    out << ",\"pair23\":";
    json::number(out, point.pair_identity[2]);
    out << '}';
  }
  out << "]}";
  return out.str();
}

}  // namespace rdp
