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
constexpr std::array<std::array<std::size_t, 2>, 3> kSourceCompRoles{{
    {{1, 2}},
    {{0, 2}},
    {{0, 1}},
}};

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
      const std::uint32_t other_one =
          representatives[kSourceCompRoles[role][0]];
      const std::uint32_t other_two =
          representatives[kSourceCompRoles[role][1]];
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
      anchor, representatives[kSourceCompRoles[role][0]]);
  const auto [parent_two_outside, ignored_two] = distance(
      anchor, representatives[kSourceCompRoles[role][1]]);
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

struct FinalTrimDistancePanel {
  double anchor_candidate = 10.0;
  std::array<double, 2> anchor_parents{10.0, 10.0};
  std::array<double, 2> candidate_parents{10.0, 10.0};
  double parent_pair = 10.0;
};

double source_finaltrim_rank_contribution(
    double value,
    double first,
    double second,
    double modifier,
    double partial_tie_reward) {
  if (value < first && value < second) return 4.0 * modifier;
  if (value <= first && value <= second && (value < first || value < second)) {
    return partial_tie_reward * modifier;
  }
  if (value > first && value > second) return -10.0 / modifier;
  if (value > first || value > second) return -2.0 / modifier;
  return 0.0;
}

double source_finaltrim_distance_score(
    const FinalTrimDistancePanel& first,
    const FinalTrimDistancePanel& second,
    bool require_candidate_affinity,
    bool include_anchor_negative_branch,
    double strongest_positive,
    double weakest_positive,
    bool suppress_weaker_positive) {
  const auto closer_to_anchor_parents = [](const FinalTrimDistancePanel& panel) {
    return panel.anchor_candidate < panel.anchor_parents[0] &&
        panel.anchor_candidate < panel.anchor_parents[1];
  };
  const auto closer_to_candidate_parents = [](const FinalTrimDistancePanel& panel) {
    return panel.anchor_candidate < panel.candidate_parents[0] &&
        panel.anchor_candidate < panel.candidate_parents[1];
  };
  if (closer_to_anchor_parents(first) && closer_to_anchor_parents(second) &&
      (!require_candidate_affinity ||
       (closer_to_candidate_parents(first) && closer_to_candidate_parents(second)))) {
    const bool first_parent_pair = first.anchor_candidate < first.parent_pair;
    const bool second_parent_pair = second.anchor_candidate < second.parent_pair;
    if (first_parent_pair && second_parent_pair) return strongest_positive;
    if (suppress_weaker_positive) return 0.0;
    if (first_parent_pair || second_parent_pair) return 2.0;
    return weakest_positive;
  }

  const auto farther_from_candidate_parents = [](const FinalTrimDistancePanel& panel) {
    return panel.anchor_candidate > panel.candidate_parents[0] &&
        panel.anchor_candidate > panel.candidate_parents[1];
  };
  const auto farther_from_anchor_parents = [](const FinalTrimDistancePanel& panel) {
    return panel.anchor_candidate > panel.anchor_parents[0] &&
        panel.anchor_candidate > panel.anchor_parents[1];
  };
  const bool negative =
      (farther_from_candidate_parents(first) &&
       farther_from_candidate_parents(second)) ||
      (include_anchor_negative_branch && farther_from_anchor_parents(first) &&
       farther_from_anchor_parents(second));
  if (!negative) return 0.0;
  const bool first_parent_pair = first.anchor_candidate > first.parent_pair;
  const bool second_parent_pair = second.anchor_candidate > second.parent_pair;
  if (first_parent_pair && second_parent_pair) return -1.0;
  if (first_parent_pair || second_parent_pair) return -0.5;
  return -0.25;
}

struct SourceDetectedTract {
  bool available = false;
  double overlap = 0.0;
  std::size_t beginning = 0;
  std::size_t ending = 0;
  std::int32_t signal_id = -1;
};

using SourceDetectedTractGrid =
    std::array<std::vector<SourceDetectedTract>, 3>;

struct SourceIntervalSegments {
  std::array<std::array<std::size_t, 2>, 2> values{};
  std::size_t count = 0;
};

SourceIntervalSegments source_interval_segments(
    std::size_t beginning,
    std::size_t ending,
    std::size_t length,
    bool preserve_equal_endpoint_revisit) {
  SourceIntervalSegments result;
  if (length == 0 || beginning < 1 || ending < 1 || beginning > length ||
      ending > length) {
    return result;
  }
  if (beginning < ending) {
    result.values[result.count++] = {beginning - 1, ending - 1};
  } else if (beginning > ending) {
    result.values[result.count++] = {beginning - 1, length - 1};
    result.values[result.count++] = {0, ending - 1};
  } else if (preserve_equal_endpoint_revisit) {
    // MakeMatchMatX2P traverses ST..L and 1..EN when ST == EN, so that
    // endpoint is counted twice even though ContainSite itself is Boolean.
    result.values[result.count++] = {beginning - 1, length - 1};
    result.values[result.count++] = {0, ending - 1};
  } else {
    result.values[result.count++] = {0, length - 1};
  }
  return result;
}

std::vector<std::array<double, 3>> source_pattern_scores(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& working_representatives,
    const std::array<std::uint32_t, 3>& reported_representatives,
    std::size_t ordinary_sequence_count,
    const CorrelationRegionLayout& layout,
    std::size_t beginning,
    std::size_t ending) {
  std::vector<std::array<double, 3>> result(ordinary_sequence_count);
  if (alignment.length == 0) return result;
  const std::array<SourceIntervalSegments, 3> regions{{
      source_interval_segments(
          layout.boundaries[0], layout.boundaries[1], alignment.length, true),
      source_interval_segments(
          layout.boundaries[2], layout.boundaries[3], alignment.length, true),
      source_interval_segments(beginning, ending, alignment.length, true),
  }};

  for (std::size_t sequence = 0; sequence < ordinary_sequence_count; ++sequence) {
    std::array<std::array<double, 3>, 3> pattern{};
    for (std::size_t region = 0; region < regions.size(); ++region) {
      for (std::size_t segment = 0; segment < regions[region].count; ++segment) {
        const auto& span = regions[region].values[segment];
        for (std::size_t position = span[0]; position <= span[1]; ++position) {
          const std::uint8_t candidate = alignment.at(sequence, position);
          const std::uint8_t first =
              alignment.at(working_representatives[0], position);
          const std::uint8_t second =
              alignment.at(working_representatives[1], position);
          const std::uint8_t third =
              alignment.at(working_representatives[2], position);
          if (candidate == 0 || first == 0 || second == 0 || third == 0 ||
              (first == second && first == third)) {
            continue;
          }
          if (candidate == first) {
            if (candidate != second) {
              if (second == third) pattern[0][region] += 1.0;
              else if (candidate != third) pattern[0][region] += 0.5;
            }
          } else if (candidate == second) {
            if (first == third) pattern[1][region] += 1.0;
            else if (candidate != third) pattern[1][region] += 0.5;
          } else if (candidate == third) {
            if (first == second) pattern[2][region] += 1.0;
            else if (candidate != second) pattern[2][region] += 0.5;
          }
        }
      }
    }
    for (std::size_t role = 0; role < 3; ++role) {
      if (sequence == reported_representatives[role]) {
        result[sequence][role] = 3.0;
        continue;
      }
      for (std::size_t region = 0; region < regions.size(); ++region) {
        const double total =
            pattern[0][region] + pattern[1][region] + pattern[2][region];
        if (total > 0.0) result[sequence][role] += pattern[role][region] / total;
      }
    }
  }
  return result;
}

double source_match_matrix_distance(
    const std::vector<std::size_t>& valid_prefix,
    const std::vector<std::size_t>& difference_prefix,
    const SourceIntervalSegments& reference_tract,
    const SourceIntervalSegments* candidate_tract) {
  std::size_t valid = 0;
  std::size_t differences = 0;
  const auto add_range = [&](std::size_t beginning, std::size_t ending) {
    if (beginning > ending || ending + 1 >= valid_prefix.size()) return;
    valid += valid_prefix[ending + 1] - valid_prefix[beginning];
    differences +=
        difference_prefix[ending + 1] - difference_prefix[beginning];
  };
  for (std::size_t reference = 0; reference < reference_tract.count; ++reference) {
    const auto& reference_segment = reference_tract.values[reference];
    if (candidate_tract == nullptr) {
      add_range(reference_segment[0], reference_segment[1]);
      continue;
    }
    for (std::size_t candidate = 0; candidate < candidate_tract->count; ++candidate) {
      const auto& candidate_segment = candidate_tract->values[candidate];
      const std::size_t beginning =
          std::max(reference_segment[0], candidate_segment[0]);
      const std::size_t ending =
          std::min(reference_segment[1], candidate_segment[1]);
      if (beginning <= ending) add_range(beginning, ending);
    }
  }

  if (valid <= 30) return 3.0;
  const float difference_fraction = static_cast<float>(differences) /
      static_cast<float>(valid);
  if (difference_fraction >= 0.75F) return 3.0;
  const float identity = 1.0F - difference_fraction;
  const float correction = (4.0F * identity - 1.0F) / 3.0F;
  const float logarithm = static_cast<float>(std::log(correction));
  return static_cast<double>(static_cast<float>(-0.75F * logarithm));
}

double source_finaltrim_detected_region_score(
    const std::array<double, 7>& match,
    bool distinct_closest_pairs,
    std::size_t role,
    const std::array<std::uint8_t, 3>& in_list,
    double repeated_pair_modifier) {
  // MatchMat fields, in active VB order: anchor/candidate, anchor/parent 0,
  // anchor/parent 1 through the source's bare CompMat sequence-index quirk,
  // anchor/parent 1 through ISeqs, candidate/parent 0, candidate/parent 1,
  // and parent 0/parent 1.
  const double anchor_candidate = match[0];
  double score = 0.0;
  if (anchor_candidate != 3.0) {
    if (anchor_candidate < match[1] && anchor_candidate < match[2]) {
      if (anchor_candidate < match[4] && anchor_candidate < match[5]) {
        if (anchor_candidate < match[6]) score = 10.0;
        else if (role != in_list[0]) score = 5.0;
      } else if (anchor_candidate > match[4] && anchor_candidate > match[5]) {
        score = -1.0;
      } else {
        score = -0.5;
      }
    } else if (anchor_candidate > match[1] && anchor_candidate > match[3]) {
      score = anchor_candidate > match[4] && anchor_candidate > match[5]
          ? -50.0
          : -10.0;
    } else {
      score = anchor_candidate > match[4] && anchor_candidate > match[5]
          ? -1.0
          : -0.5;
    }
  } else {
    score = -5.0;
  }
  if (score > 0.0) {
    if (distinct_closest_pairs) {
      if (role == in_list[0]) score /= 2.0;
    } else {
      score *= repeated_pair_modifier;
    }
  }
  return score;
}

using SourceCalcMatchGrid = std::array<std::vector<CalcMatchEvidence>, 3>;

bool source_calc_match_variable_site(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& representatives,
    std::size_t coordinate) {
  if (coordinate < 1 || coordinate > alignment.length) return false;
  const std::size_t position = coordinate - 1;
  const std::uint8_t first = alignment.at(representatives[0], position);
  const std::uint8_t second = alignment.at(representatives[1], position);
  const std::uint8_t third = alignment.at(representatives[2], position);
  return first != 0 && second != 0 && third != 0 &&
      (first != second || first != third);
}

std::size_t source_vb_clng_half(std::size_t value) {
  // VB6 CLng uses banker's rounding. Preserve the half-to-even cases used by
  // CalcMatchY when fewer than 30 variable sites occur inside the tract.
  const std::size_t quotient = value / 2;
  if (value % 2 != 0 && quotient % 2 != 0) return quotient + 1;
  return quotient;
}

std::int64_t source_vb_long(double value) {
  // VB6's narrowing conversion to Long uses round-to-nearest with ties to an
  // even integer. ConsensusOK declares NS As Long, so every assignment back
  // to NS narrows before the following matrix contribution.
  if (!std::isfinite(value)) return 0;
  const double lower = std::floor(value);
  const double fraction = value - lower;
  const auto lower_integer = static_cast<std::int64_t>(lower);
  if (fraction < 0.5) return lower_integer;
  if (fraction > 0.5) return lower_integer + 1;
  return lower_integer % 2 == 0 ? lower_integer : lower_integer + 1;
}

std::int8_t source_calc_match_class(const std::array<double, 6>& match) {
  const double beginning_mean = (match[0] + match[1] + match[4]) / 3.0;
  const double ending_mean = (match[2] + match[3] + match[5]) / 3.0;
  if ((match[0] > 0.8 && match[1] > 0.8 && match[4] > 0.8) ||
      (match[2] > 0.8 && match[3] > 0.8 && match[5] > 0.8)) {
    return 1;
  }
  if (match[4] > 0.9 || match[5] > 0.9) return 1;
  if (match[4] > 0.8 || match[5] > 0.8) return 2;
  if ((match[4] > 0.75 || match[5] > 0.75) &&
      (match[1] > 0.75 || match[2] > 0.75) &&
      (match[0] > 0.75 || match[3] > 0.75)) {
    if (beginning_mean > 0.75 || ending_mean > 0.75) return 2;
    if ((match[0] > 0.75 && match[1] > 0.75 && match[4] > 0.75) ||
        (match[2] > 0.75 && match[3] > 0.75 && match[5] > 0.75)) {
      return 2;
    }
    if (match[0] > 0.7 && match[1] > 0.7 && match[4] > 0.7 &&
        match[2] > 0.7 && match[3] > 0.7 && match[5] > 0.7) {
      return 2;
    }
  } else if ((match[4] < 0.65 && match[5] < 0.65) ||
             (match[1] < 0.65 && match[2] < 0.65) ||
             (match[0] < 0.65 && match[3] < 0.65) ||
             (match[0] < 0.7 && match[1] < 0.7 && match[4] < 0.7 &&
              match[2] < 0.7 && match[3] < 0.7 && match[5] < 0.7) ||
             (beginning_mean < 0.7 && ending_mean < 0.7)) {
    return -1;
  }
  return 0;
}

SourceCalcMatchGrid source_calc_match(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& representatives,
    std::size_t ordinary_sequence_count,
    std::size_t beginning,
    std::size_t ending) {
  SourceCalcMatchGrid result;
  for (auto& role : result) role.resize(ordinary_sequence_count);
  if (alignment.length == 0 || ordinary_sequence_count == 0 || beginning < 1 ||
      ending < 1 || beginning > alignment.length || ending > alignment.length ||
      std::any_of(
          representatives.begin(),
          representatives.end(),
          [&](std::uint32_t sequence) { return sequence >= alignment.sequence_count(); })) {
    return result;
  }

  std::size_t tract_variable_sites = 0;
  std::size_t coordinate = beginning;
  for (std::size_t visited = 0; visited < alignment.length; ++visited) {
    if (source_calc_match_variable_site(alignment, representatives, coordinate)) {
      ++tract_variable_sites;
    }
    coordinate = wrap_coordinate(
        static_cast<std::int64_t>(coordinate) + 1, alignment.length);
    if (coordinate == ending) break;
  }
  const std::size_t target_half_window = tract_variable_sites >= 30
      ? 15
      : source_vb_clng_half(tract_variable_sites);
  if (target_half_window < 2) return result;
  const std::size_t smoothing_half_window = std::max<std::size_t>(5, target_half_window);

  struct SourceFlankScan {
    std::size_t bound = 1;
    std::int64_t cycles = 0;
    bool valid = true;
  };
  const std::size_t source_fragment_limit = alignment.length >
          std::numeric_limits<std::size_t>::max() / 3
      ? std::numeric_limits<std::size_t>::max()
      : alignment.length * 3;
  std::size_t scanned_sites = 0;
  const auto scan_flank = [&](std::size_t start,
                              int direction,
                              std::int64_t initial_cycles) {
    SourceFlankScan scan;
    scan.bound = start;
    scan.cycles = initial_cycles;
    std::size_t variable_sites = 0;
    while (true) {
      ++scanned_sites;
      if (source_calc_match_variable_site(alignment, representatives, scan.bound)) {
        ++variable_sites;
        if (variable_sites == 40) return scan;
      }
      if (scanned_sites > source_fragment_limit) {
        scan.valid = false;
        return scan;
      }
      if (direction < 0) {
        if (scan.bound == 1) {
          scan.bound = alignment.length;
          ++scan.cycles;
          if (scan.cycles > 40 || (scan.cycles > 2 && variable_sites == 0)) return scan;
        } else {
          --scan.bound;
        }
      } else {
        if (scan.bound == alignment.length) {
          scan.bound = 1;
          ++scan.cycles;
          if (scan.cycles > 40 || (scan.cycles > 2 && variable_sites == 0)) return scan;
        } else {
          ++scan.bound;
        }
      }
    }
  };

  const SourceFlankScan beginning_left = scan_flank(beginning, -1, 0);
  const std::size_t pseudo_beginning = scanned_sites;
  if (!beginning_left.valid || scanned_sites > source_fragment_limit) return result;
  std::size_t beginning_right_start = beginning + 1;
  std::int64_t beginning_right_cycles = 0;
  if (beginning_right_start > alignment.length) {
    beginning_right_start = 1;
    ++beginning_right_cycles;
  }
  const SourceFlankScan beginning_right =
      scan_flank(beginning_right_start, 1, beginning_right_cycles);
  const SourceFlankScan ending_left = scan_flank(ending, -1, 0);
  const std::size_t pseudo_ending = scanned_sites;
  if (!beginning_right.valid || !ending_left.valid ||
      scanned_sites > source_fragment_limit) {
    return result;
  }
  std::size_t ending_right_start = ending + 1;
  std::int64_t ending_right_cycles = 0;
  if (ending_right_start > alignment.length) {
    ending_right_start = 1;
    ++ending_right_cycles;
  }
  const SourceFlankScan ending_right =
      scan_flank(ending_right_start, 1, ending_right_cycles);
  if (!ending_right.valid || scanned_sites > source_fragment_limit) return result;

  // MakeLenFrag allocates from the four scan lengths, writes from slot one,
  // and leaves slot zero plus any unused tail as gaps. Retain those sentinel
  // slots because VRPos and the breakpoint samples depend on their offsets.
  const std::size_t fragment_upper_bound = scanned_sites + 1;
  std::vector<std::size_t> fragment_coordinate(fragment_upper_bound + 1, 0);
  std::size_t fragment_length = 0;
  bool reconstruction_valid = true;
  const auto append = [&](std::size_t source_coordinate) {
    if (fragment_length >= fragment_upper_bound) {
      reconstruction_valid = false;
      return;
    }
    fragment_coordinate[++fragment_length] = source_coordinate;
  };
  const auto append_forward_until = [&](std::size_t start,
                                        std::size_t stop,
                                        std::int64_t cycles,
                                        bool nonpositive_stop) {
    std::size_t current = start;
    while (reconstruction_valid) {
      append(current);
      if (!reconstruction_valid) break;
      if (current == alignment.length) {
        current = 1;
        --cycles;
        if (current == stop && (nonpositive_stop ? cycles <= 0 : cycles == 0)) break;
      } else {
        ++current;
        if (current == stop && (nonpositive_stop ? cycles <= 0 : cycles == 0)) break;
      }
    }
    return current;
  };
  coordinate = append_forward_until(
      beginning_left.bound, beginning, beginning_left.cycles, true);
  coordinate = append_forward_until(
      coordinate, beginning_right.bound, beginning_right.cycles, true);
  coordinate = append_forward_until(
      ending_left.bound, ending, ending_left.cycles, false);
  append_forward_until(coordinate, ending_right.bound, ending_right.cycles, false);
  if (!reconstruction_valid || fragment_length > source_fragment_limit ||
      pseudo_beginning > fragment_upper_bound || pseudo_ending > fragment_upper_bound) {
    return result;
  }

  std::vector<std::size_t> variable_coordinates(1, 0);
  std::vector<std::size_t> variable_rank(fragment_upper_bound + 2, 0);
  for (std::size_t slot = 0; slot <= fragment_upper_bound; ++slot) {
    variable_rank[slot] = variable_coordinates.size() - 1;
    const std::size_t source_coordinate = fragment_coordinate[slot];
    if (source_coordinate != 0 &&
        source_calc_match_variable_site(alignment, representatives, source_coordinate)) {
      variable_coordinates.push_back(source_coordinate);
    }
  }
  variable_rank[fragment_upper_bound + 1] = fragment_upper_bound + 1;
  const std::size_t fragment_variable_sites = variable_coordinates.size() - 1;
  if (fragment_variable_sites == 0 || fragment_variable_sites > 160 ||
      fragment_variable_sites < smoothing_half_window || pseudo_beginning < 2 ||
      pseudo_ending + 1 >= variable_rank.size()) {
    return result;
  }

  const auto source_wrapped_variable_index = [&](std::int64_t index) {
    if (index < 1) index += static_cast<std::int64_t>(fragment_variable_sites);
    else if (index > static_cast<std::int64_t>(fragment_variable_sites)) {
      index -= static_cast<std::int64_t>(fragment_variable_sites);
    }
    return index;
  };
  const auto source_checkpoint_index = [&](std::int64_t index) {
    // The KeepTrack pass differs subtly from MakeCntHit2: zero is a valid
    // sentinel sample here and only negative indices wrap to the tail.
    if (index < 0) index += static_cast<std::int64_t>(fragment_variable_sites);
    else if (index > static_cast<std::int64_t>(fragment_variable_sites)) {
      index -= static_cast<std::int64_t>(fragment_variable_sites);
    }
    return index;
  };
  const auto add_score_range = [](const std::vector<float>& smooth,
                                  std::size_t first,
                                  std::size_t last,
                                  double& total,
                                  std::size_t& count) {
    if (first > last) return true;
    if (last >= smooth.size()) return false;
    for (std::size_t index = first; index <= last; ++index) {
      ++count;
      if (smooth[index] > 0.6F) {
        total += (static_cast<double>(smooth[index]) - 0.6) / 0.4;
      }
    }
    return true;
  };

  for (std::size_t role = 0; role < representatives.size(); ++role) {
    const std::size_t other_one = kSourceCompRoles[role][0];
    const std::size_t other_two = kSourceCompRoles[role][1];
    std::vector<std::int16_t> match_map(fragment_variable_sites + 1, 0);
    std::vector<float> smooth(fragment_variable_sites + 1, 0.0F);
    for (std::size_t candidate = 0; candidate < ordinary_sequence_count; ++candidate) {
      CalcMatchEvidence evidence;
      evidence.fragment_variable_sites = fragment_variable_sites;
      evidence.target_half_window = target_half_window;
      evidence.smoothing_half_window = smoothing_half_window;
      std::fill(match_map.begin(), match_map.end(), 0);
      for (std::size_t variable = 1; variable <= fragment_variable_sites; ++variable) {
        const std::size_t position = variable_coordinates[variable] - 1;
        const std::uint8_t anchor = alignment.at(representatives[role], position);
        const std::uint8_t first_other = alignment.at(representatives[other_one], position);
        const std::uint8_t second_other = alignment.at(representatives[other_two], position);
        const std::uint8_t state = alignment.at(candidate, position);
        if (state == anchor) {
          match_map[variable] = 2;
        } else if (state != 0) {
          if ((state == first_other && anchor == second_other) ||
              (anchor == first_other && state == second_other)) {
            match_map[variable] = -1;
          } else if ((state == first_other && state != second_other) ||
                     (state != first_other && state == second_other)) {
            match_map[variable] = 0;
          } else if (state != first_other && state != second_other) {
            match_map[variable] = 1;
          }
        }
      }

      std::fill(smooth.begin(), smooth.end(), 0.0F);
      std::int64_t rolling_total = 0;
      bool candidate_valid = true;
      for (std::int64_t index = 1 - static_cast<std::int64_t>(smoothing_half_window);
           index <= 1 + static_cast<std::int64_t>(smoothing_half_window);
           ++index) {
        const std::int64_t wrapped = source_wrapped_variable_index(index);
        if (wrapped < 0 || wrapped > static_cast<std::int64_t>(fragment_variable_sites)) {
          candidate_valid = false;
          break;
        }
        rolling_total += match_map[static_cast<std::size_t>(wrapped)];
      }
      const double denominator =
          static_cast<double>(smoothing_half_window * 2 + 1) * 2.0;
      if (candidate_valid) {
        smooth[1] = static_cast<float>(static_cast<double>(rolling_total) / denominator);
      }
      for (std::size_t variable = 2;
           candidate_valid && variable <= fragment_variable_sites;
           ++variable) {
        const std::int64_t remove_index = source_wrapped_variable_index(
            static_cast<std::int64_t>(variable) -
            static_cast<std::int64_t>(smoothing_half_window) - 1);
        const std::int64_t add_index = source_wrapped_variable_index(
            static_cast<std::int64_t>(variable) +
            static_cast<std::int64_t>(smoothing_half_window));
        if (remove_index < 0 || add_index < 0 ||
            remove_index > static_cast<std::int64_t>(fragment_variable_sites) ||
            add_index > static_cast<std::int64_t>(fragment_variable_sites)) {
          candidate_valid = false;
          break;
        }
        rolling_total -= match_map[static_cast<std::size_t>(remove_index)];
        rolling_total += match_map[static_cast<std::size_t>(add_index)];
        smooth[variable] =
            static_cast<float>(static_cast<double>(rolling_total) / denominator);
      }

      const std::size_t source_st = pseudo_beginning > 1
          ? pseudo_beginning - 1
          : pseudo_beginning;
      double outside_total = 0.0;
      std::size_t outside_count = 0;
      double inside_total = 0.0;
      std::size_t inside_count = 0;
      candidate_valid = candidate_valid &&
          add_score_range(
              smooth,
              1,
              variable_rank[source_st - 1],
              outside_total,
              outside_count) &&
          add_score_range(
              smooth,
              variable_rank[pseudo_ending + 1],
              variable_rank[fragment_upper_bound],
              outside_total,
              outside_count) &&
          add_score_range(
              smooth,
              variable_rank[pseudo_beginning],
              variable_rank[pseudo_ending],
              inside_total,
              inside_count);
      if (!candidate_valid || outside_count == 0 || inside_count == 0) {
        result[role][candidate] = evidence;
        continue;
      }
      const double outside_match = outside_total / static_cast<double>(outside_count);
      const double inside_match = inside_total / static_cast<double>(inside_count);
      evidence.regional_match_score = outside_match * inside_match;

      constexpr std::array<std::size_t, 6> checkpoint_sources{{0, 0, 1, 1, 0, 1}};
      constexpr std::array<std::int8_t, 6> checkpoint_offsets{{-1, 1, -1, 1, 0, 0}};
      const std::array<std::size_t, 2> checkpoint_centers{
          variable_rank[pseudo_beginning],
          variable_rank[pseudo_ending],
      };
      for (std::size_t checkpoint = 0; checkpoint < evidence.checkpoint_matches.size();
           ++checkpoint) {
        std::int64_t sample = static_cast<std::int64_t>(
            checkpoint_centers[checkpoint_sources[checkpoint]]);
        sample += static_cast<std::int64_t>(checkpoint_offsets[checkpoint]) *
            static_cast<std::int64_t>(smoothing_half_window);
        sample = source_checkpoint_index(sample);
        if (sample < 0 || sample > static_cast<std::int64_t>(fragment_variable_sites)) {
          candidate_valid = false;
          break;
        }
        evidence.checkpoint_matches[checkpoint] = smooth[static_cast<std::size_t>(sample)];
      }
      if (!candidate_valid) {
        result[role][candidate] = evidence;
        continue;
      }
      evidence.raw_breakpoint_match_class =
          source_calc_match_class(evidence.checkpoint_matches);
      evidence.breakpoint_match_class = evidence.raw_breakpoint_match_class;
      evidence.breakpoints_exist[0] =
          evidence.checkpoint_matches[0] > 0.75 ||
              evidence.checkpoint_matches[1] > 0.75 ||
              evidence.checkpoint_matches[4] > 0.75
          ? 1
          : 0;
      evidence.breakpoints_exist[1] =
          evidence.checkpoint_matches[2] > 0.75 ||
              evidence.checkpoint_matches[3] > 0.75 ||
              evidence.checkpoint_matches[5] > 0.75
          ? 1
          : 0;
      evidence.available = true;
      result[role][candidate] = evidence;
    }
  }
  return result;
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
  // Active ModSeqNumY uses its wrapped branch when BPos == EPos, which
  // traverses B..L and 1..E and therefore covers the complete alignment.
  if (beginning == ending) return true;
  return wraps_origin
      ? coordinate >= beginning || coordinate <= ending
      : coordinate >= beginning && coordinate <= ending;
}

void write_breakpoint_confidence_json(
    std::ostringstream& out,
    const BreakpointConfidenceEvidence& confidence) {
  out << "{\"status\":\""
      << (confidence.available ? "complete-active-unvalidated" : "unavailable")
      << "\",\"method\":\"BURT/BenHMM\",\"attempted\":"
      << (confidence.attempted ? "true" : "false")
      << ",\"available\":" << (confidence.available ? "true" : "false")
      << ",\"appliedToEvent\":"
      << (confidence.applied_to_event ? "true" : "false")
      << ",\"informationRichSites\":" << confidence.information_rich_sites
      << ",\"candidateIntervalCount\":" << confidence.candidate_interval_count
      << ",\"bestLogLikelihood\":";
  json::number(out, confidence.best_log_likelihood);
  out << ",\"randomSeed\":" << confidence.random_seed
      << ",\"randomAdapter\":\"msvc-rand-15-bit\""
      << ",\"sourceRandomAdapter\":"
      << (confidence.source_random_adapter ? "true" : "false")
      << ",\"hmmCyclesArgument\":" << confidence.hmm_cycles_argument
      << ",\"serialTrainingStarts\":" << confidence.serial_training_starts
      << ",\"posteriorThresholds\":[0.995,0.999]"
      << ",\"polishedBeginning\":" << confidence.polished_beginning
      << ",\"polishedEnding\":" << confidence.polished_ending
      << ",\"singleTransitionAssignment\":"
      << (confidence.single_transition_assignment ? "true" : "false")
      << ",\"insufficientInsideOrOutsideReverted\":"
      << (confidence.insufficient_inside_or_outside_reverted ? "true" : "false")
      << ",\"unavailableReason\":";
  if (confidence.unavailable_reason.empty()) out << "null";
  else json::string(out, confidence.unavailable_reason);
  out << ",\"boundaries\":[";
  for (std::size_t index = 0; index < confidence.boundaries.size(); ++index) {
    if (index) out << ',';
    const auto& boundary = confidence.boundaries[index];
    out << "{\"name\":\"" << (index == 0 ? "beginning" : "ending")
        << "\",\"inputCoordinate\":" << boundary.input_coordinate
        << ",\"polishedCoordinate\":" << boundary.polished_coordinate
        << ",\"intervalAvailable\":"
        << (boundary.interval_available ? "true" : "false")
        << ",\"sourceIntervalContainsInput\":"
        << (boundary.source_interval_contains_input ? "true" : "false")
        << ",\"confidence99\":{\"beginning\":"
        << boundary.confidence_99_beginning << ",\"ending\":"
        << boundary.confidence_99_ending << ",\"wrapsOrigin\":"
        << (boundary.confidence_99_wraps_origin ? "true" : "false")
        << "},\"hmmCoordinate\":" << boundary.hmm_coordinate
        << ",\"confidence95\":{\"beginning\":"
        << boundary.confidence_95_beginning << ",\"ending\":"
        << boundary.confidence_95_ending << ",\"wrapsOrigin\":"
        << (boundary.confidence_95_wraps_origin ? "true" : "false")
        << "},\"repositioned\":" << (boundary.repositioned ? "true" : "false")
        << ",\"missingDataAdjusted\":"
        << (boundary.missing_data_adjusted ? "true" : "false")
        << ",\"finalGapAdjusted\":"
        << (boundary.final_gap_adjusted ? "true" : "false") << '}';
  }
  out << "]}";
}

const char* maxchi_recheck_status(const MaxChiRecheckEvidence& evidence) {
  if (evidence.representative_skipped) {
    return "representative-skipped";
  }
  if (!evidence.requested) return "not-in-final-distance-list";
  if (!evidence.profile_available) return "profile-unavailable";
  return "complete-active-unvalidated";
}

void write_maxchi_recheck_json(
    std::ostringstream& out,
    const MaxChiRecheckEvidence& evidence) {
  out << "{\"status\":\"" << maxchi_recheck_status(evidence)
      << "\",\"kernel\":\"FastRecCheckMC2-strongest-peak\""
      << ",\"eventDiscoveryApplied\":false"
      << ",\"requested\":" << (evidence.requested ? "true" : "false")
      << ",\"representativeSkipped\":"
      << (evidence.representative_skipped ? "true" : "false")
      << ",\"profileAvailable\":"
      << (evidence.profile_available ? "true" : "false")
      << ",\"missingDataWindowFilterApplied\":"
      << (evidence.missing_data_window_filter_applied ? "true" : "false")
      << ",\"linearEdgeWindowFilterApplied\":"
      << (evidence.linear_edge_window_filter_applied ? "true" : "false")
      << ",\"bonferroniApplied\":"
      << (evidence.bonferroni_applied ? "true" : "false")
      << ",\"correctionTests\":" << evidence.correction_tests
      << ",\"variableSites\":" << evidence.variable_sites
      << ",\"fixedWindowSites\":70"
      << ",\"halfWindow\":" << evidence.half_window
      << ",\"criticalDifference\":" << evidence.critical_difference
      << ",\"grownHalfWindow\":" << evidence.grown_half_window
      << ",\"bestPair\":";
  if (evidence.best_pair < 0) out << "null";
  else out << static_cast<unsigned int>(evidence.best_pair);
  out << ",\"peakAlignmentPosition\":";
  if (evidence.best_pair < 0) out << "null";
  else out << evidence.peak_alignment_position;
  out << ",\"maximumChiSquare\":";
  if (evidence.best_pair < 0) out << "null";
  else json::number(out, evidence.maximum_chi_square);
  out << ",\"localPValue\":";
  if (evidence.best_pair < 0) out << "null";
  else json::number(out, evidence.local_p_value);
  out << ",\"withinTripletPValue\":";
  if (evidence.best_pair < 0) out << "null";
  else json::number(out, evidence.within_triplet_p_value);
  out << ",\"correctedPValue\":";
  if (evidence.best_pair < 0) out << "null";
  else json::number(out, evidence.corrected_p_value);
  out << ",\"sourceRecheckHit\":"
      << (evidence.source_recheck_hit ? "true" : "false") << '}';
}

}  // namespace

std::string curated_sequences_fasta(
    const Alignment& alignment,
    const std::vector<std::uint8_t>& mask,
    const std::vector<std::uint8_t>& disabled,
    bool include_enabled) {
  std::ostringstream out;
  for (std::size_t sequence = 0; sequence < alignment.sequence_count(); ++sequence) {
    const bool masked = sequence < mask.size() && mask[sequence] != 0;
    const bool is_disabled = sequence < disabled.size() && disabled[sequence] != 0;
    const bool enabled = !masked && !is_disabled;
    if (enabled != include_enabled) continue;
    write_fasta_record(out, alignment.names[sequence], alignment.sequences[sequence]);
  }
  return out.str();
}

RdpScanner::RdpScanner(const Alignment& alignment)
    : alignment_(alignment),
      native_input_missing_data_(alignment.sequence_count() * alignment.length, 0) {
  // ModSeqNum marks native MissingData after ten consecutive SeqNum values
  // below ASCII 50, including its literal one-position look-back when the
  // counter first reaches ten. Build that immutable baseline once; cyclic
  // erasures are layered onto it per event in refresh_breakpoint_context.
  for (std::size_t sequence = 0; sequence < alignment_.sequence_count(); ++sequence) {
    std::size_t consecutive_missing = 0;
    for (std::size_t coordinate = 1; coordinate <= alignment_.length; ++coordinate) {
      const auto character = static_cast<unsigned char>(
          alignment_.sequences[sequence][coordinate - 1]);
      if (character < 50) {
        ++consecutive_missing;
        if (consecutive_missing > 10) {
          native_input_missing_data_[
              sequence * alignment_.length + coordinate - 1] = 1;
        } else if (consecutive_missing == 10) {
          const std::size_t first = coordinate > 10 ? coordinate - 10 : 1;
          for (std::size_t marked = first; marked <= coordinate; ++marked) {
            native_input_missing_data_[
                sequence * alignment_.length + marked - 1] = 1;
          }
        }
      } else {
        consecutive_missing = 0;
      }
    }
  }
}

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

bool RdpScanner::sequence_masked(std::uint32_t sequence) const {
  return sequence < options_.mask.size() && options_.mask[sequence] != 0;
}

bool RdpScanner::sequence_disabled(std::uint32_t sequence) const {
  return sequence < options_.disabled.size() && options_.disabled[sequence] != 0;
}

void RdpScanner::refresh_active_sequences() {
  active_sequences_.clear();
  const std::size_t fragment_minimum = std::max<std::size_t>(5, options_.window_sites);
  for (std::size_t sequence = 0; sequence < working_alignment_.sequence_count(); ++sequence) {
    if (sequence >= working_origins_.size()) continue;
    const std::uint32_t origin = working_origins_[sequence];
    if (origin >= options_.mask.size() || origin >= options_.disabled.size() ||
        sequence_masked(origin) || sequence_disabled(origin)) {
      continue;
    }
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
  signal_candidates_scratch_.clear();
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
  if (options.disabled.size() != alignment_.sequence_count()) {
    options.disabled.assign(alignment_.sequence_count(), 0);
  }
  for (std::size_t sequence = 0; sequence < alignment_.sequence_count(); ++sequence) {
    if (options.disabled[sequence] != 0) options.mask[sequence] = 0;
  }

  options_ = std::move(options);
  reset_working_alignment();
  refresh_active_sequences();
  if (active_sequences_.size() < 3) {
    error = "At least three enabled, unmasked sequences are required for an exploratory RDP scan.";
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

void RdpScanner::append_candidate_signals(
    const TripletProfile& profile,
    std::uint8_t high_pair,
    std::uint8_t candidate_pair,
    std::uint8_t low_pair,
    bool enforce_cutoff,
    std::vector<Signal>& output) const {
  const std::size_t length = profile.category.size();
  std::size_t search = 0;

  const auto candidate_members = pair_members(candidate_pair);
  const auto high_members = pair_members(high_pair);
  std::uint8_t recombinant_local = 255;
  for (const auto member : candidate_members) {
    if (member == high_members[0] || member == high_members[1]) recombinant_local = member;
  }
  if (recombinant_local == 255) return;
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
    // FastRecCheckP accepts a detected tract only when EN != BE. Equal
    // endpoints remain meaningful for manual/restored circular events (where
    // ModSeqNumY traverses the full alignment), but primary detection must not
    // manufacture a whole-alignment event from a full circular candidate run.
    if (beginning != ending && region_length > 2 &&
        matching > static_cast<double>(different) * 0.8) {
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
        signal.wraps_origin = wrapped ||
            (options_.circular && signal.beginning >= signal.ending);
        signal.informative_beginning = beginning + 1;
        signal.informative_ending = ending + 1;
        signal.local_p_value = local;
        signal.corrected_p_value = corrected;
        signal.correction_tests = correction_tests_;
        signal.pair_similarity = profile.similarities;
        signal.informative_sites = length;
        signal.candidate_pair = candidate_pair;
        output.push_back(signal);
      }
    }

    if (wrapped) break;
    search = std::max(search + 1, ending + 1);
  }
}

std::vector<Signal> RdpScanner::triplet_signals(
    const std::array<std::uint32_t, 3>& triplet,
    bool enforce_cutoff,
    bool* profile_available,
    TripletProfile* scratch) const {
  if (profile_available) *profile_available = false;
  TripletProfile local_profile;
  TripletProfile& profile = scratch ? *scratch : local_profile;
  if (!build_profile(triplet, profile)) return {};
  if (profile_available) *profile_available = true;
  const auto order = ranked_pairs(profile);
  const std::array<double, 3> average{
      static_cast<double>(profile.category_counts[0]) / profile.category.size(),
      static_cast<double>(profile.category_counts[1]) / profile.category.size(),
      static_cast<double>(profile.category_counts[2]) / profile.category.size(),
  };
  std::vector<Signal> candidates;
  candidates.reserve(4);
  append_candidate_signals(
      profile, order[0], order[1], order[2], enforce_cutoff, candidates);
  append_candidate_signals(
      profile, order[0], order[2], order[1], enforce_cutoff, candidates);
  if (average[order[0]] < 0.7) {
    append_candidate_signals(
        profile, order[1], order[0], order[2], enforce_cutoff, candidates);
  }
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
  auto& candidates = signal_candidates_scratch_;
  candidates.clear();
  if (candidates.capacity() == 0) candidates.reserve(4);
  append_candidate_signals(profile, order[0], order[1], order[2], true, candidates);
  append_candidate_signals(profile, order[0], order[2], order[1], true, candidates);
  if (average[order[0]] < 0.7) {
    append_candidate_signals(profile, order[1], order[0], order[2], true, candidates);
  }

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
  struct SegmentSet {
    std::array<Segment, 2> values{};
    std::size_t count = 0;
  };
  const auto segments = [&](std::size_t beginning, std::size_t ending) {
    SegmentSet result;
    beginning = std::clamp<std::size_t>(beginning, 1, alignment_.length);
    ending = std::clamp<std::size_t>(ending, 1, alignment_.length);
    if (beginning < ending) {
      result.values[result.count++] = {beginning, ending};
    } else if (beginning > ending) {
      result.values[result.count++] = {beginning, alignment_.length};
      result.values[result.count++] = {1, ending};
    } else {
      result.values[result.count++] = {1, alignment_.length};
    }
    return result;
  };
  const auto first = segments(first_beginning, first_ending);
  const auto second = segments(second_beginning, second_ending);
  const auto size = [](const SegmentSet& segments) {
    std::size_t total = 0;
    for (std::size_t index = 0; index < segments.count; ++index) {
      const auto& segment = segments.values[index];
      total += segment[1] - segment[0] + 1;
    }
    return total;
  };
  std::size_t intersection = 0;
  for (std::size_t left_index = 0; left_index < first.count; ++left_index) {
    const auto& left = first.values[left_index];
    for (std::size_t right_index = 0; right_index < second.count; ++right_index) {
      const auto& right = second.values[right_index];
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

void RdpScanner::refresh_breakpoint_confidence(
    UniqueEvent& event,
    bool apply_polished) {
  if (!options_.polish_breakpoints) {
    event.breakpoint_confidence = {};
    event.breakpoint_confidence.unavailable_reason = "disabled";
    event.breakpoint_confidence.polished_beginning = event.beginning;
    event.breakpoint_confidence.polished_ending = event.ending;
    event.breakpoint_confidence.boundaries[0].input_coordinate = event.beginning;
    event.breakpoint_confidence.boundaries[0].polished_coordinate = event.beginning;
    event.breakpoint_confidence.boundaries[1].input_coordinate = event.ending;
    event.breakpoint_confidence.boundaries[1].polished_coordinate = event.ending;
    return;
  }
  const std::array<std::uint32_t, 3> representatives{
      event.recombinant,
      event.major_parent,
      event.minor_parent,
  };
  auto& triplet_missing = breakpoint_polish_missing_scratch_;
  triplet_missing.assign(alignment_.length + 1, 0);
  for (std::size_t coordinate = 1; coordinate <= alignment_.length; ++coordinate) {
    for (const std::uint32_t representative : representatives) {
      const std::size_t offset =
          static_cast<std::size_t>(representative) * alignment_.length + coordinate - 1;
      if (offset < native_input_missing_data_.size() &&
          native_input_missing_data_[offset] != 0) {
        triplet_missing[coordinate] = 1;
        break;
      }
    }
  }
  // PolishBP runs before ModSeqNumY erases the current event, but all earlier
  // accepted cyclic events have already written MissingData=1 over their
  // inclusive recombinant tracts. Reconstruct that triplet union without an
  // N x L copy so BURT's CI repositioning sees the same accumulated history.
  for (std::size_t prior_index = 0;
       prior_index < event.id && prior_index < events_.size();
       ++prior_index) {
    const auto& prior = events_[prior_index];
    if (!prior.tract_erased_for_detection) continue;
    const bool affects_triplet = std::any_of(
        representatives.begin(),
        representatives.end(),
        [&](std::uint32_t sequence) {
          return sequence == prior.recombinant ||
              std::find(
                  prior.co_recombinant_sequences.begin(),
                  prior.co_recombinant_sequences.end(),
                  sequence) != prior.co_recombinant_sequences.end();
        });
    if (!affects_triplet) continue;
    for (std::size_t coordinate = 1; coordinate <= alignment_.length; ++coordinate) {
      if (coordinate_in_tract(
              coordinate,
              prior.beginning,
              prior.ending,
              prior.wraps_origin)) {
        triplet_missing[coordinate] = 1;
      }
    }
  }
  event.breakpoint_confidence = source_polish_breakpoints(
      working_alignment_,
      triplet_missing,
      representatives,
      event.beginning,
      event.ending,
      options_.circular,
      breakpoint_confidence_scratch_);
  if (!apply_polished || !event.breakpoint_confidence.available) return;

  const std::size_t polished_beginning =
      event.breakpoint_confidence.polished_beginning;
  const std::size_t polished_ending =
      event.breakpoint_confidence.polished_ending;
  if (polished_beginning < 1 || polished_beginning > alignment_.length ||
      polished_ending < 1 || polished_ending > alignment_.length) {
    return;
  }
  event.beginning = polished_beginning;
  event.ending = polished_ending;
  event.wraps_origin = options_.circular && event.beginning >= event.ending;
  event.breakpoint_confidence.applied_to_event = true;
}

void RdpScanner::refresh_breakpoint_context(UniqueEvent& event) {
  for (auto& ids : event.adjacent_erasure_event_ids) ids.clear();
  for (auto& ids : event.uncertain_erasure_event_ids) ids.clear();
  event.nearest_erasure_informative_sites = {-1, -1};
  event.breakpoint_uncertainty = {};
  if (alignment_.length == 0) return;

  const auto shifted_coordinate = [&](std::size_t center, std::int64_t offset) {
    const std::int64_t length = static_cast<std::int64_t>(alignment_.length);
    const std::int64_t unwrapped = static_cast<std::int64_t>(center) + offset;
    if (!options_.circular && (unwrapped < 1 || unwrapped > length)) {
      return std::size_t{0};
    }
    return options_.circular
        ? static_cast<std::size_t>(((unwrapped - 1) % length + length) % length + 1)
        : static_cast<std::size_t>(unwrapped);
  };
  const std::array<std::uint32_t, 3> representatives{
      event.recombinant, event.major_parent, event.minor_parent};

  auto& erased_for_representative = breakpoint_erasure_scratch_;
  for (auto& erased : erased_for_representative) {
    erased.assign(alignment_.length + 1, 0);
  }
  auto& input_missing_for_triplet = breakpoint_input_missing_scratch_;
  input_missing_for_triplet.assign(alignment_.length + 1, 0);
  for (std::size_t coordinate = 1; coordinate <= alignment_.length; ++coordinate) {
    for (const std::uint32_t representative : representatives) {
      if (native_input_missing_data_[
              representative * alignment_.length + coordinate - 1] != 0) {
        input_missing_for_triplet[coordinate] = 1;
        break;
      }
    }
  }

  auto& relevant_prior_events = breakpoint_relevant_event_indices_scratch_;
  relevant_prior_events.clear();
  relevant_prior_events.reserve(event.id);
  for (std::size_t prior_index = 0;
       prior_index < event.id && prior_index < events_.size();
       ++prior_index) {
    const auto& prior = events_[prior_index];
    if (!prior.tract_erased_for_detection) continue;
    std::array<bool, 3> affected{};
    for (std::size_t role = 0; role < representatives.size(); ++role) {
      const std::uint32_t sequence = representatives[role];
      affected[role] = sequence == prior.recombinant ||
          std::find(
              prior.co_recombinant_sequences.begin(),
              prior.co_recombinant_sequences.end(),
              sequence) != prior.co_recombinant_sequences.end();
    }
    if (std::none_of(affected.begin(), affected.end(), [](bool value) { return value; })) {
      continue;
    }
    relevant_prior_events.push_back(static_cast<std::uint32_t>(prior_index));
    const auto mark_erased_range = [&](std::size_t first, std::size_t last) {
      if (first > last || first < 1 || last > alignment_.length) return;
      for (std::size_t coordinate = first; coordinate <= last; ++coordinate) {
        for (std::size_t role = 0; role < affected.size(); ++role) {
          if (affected[role]) erased_for_representative[role][coordinate] = 1;
        }
      }
    };
    if (prior.beginning == prior.ending) {
      mark_erased_range(1, alignment_.length);
    } else if (prior.wraps_origin) {
      mark_erased_range(prior.beginning, alignment_.length);
      mark_erased_range(1, prior.ending);
    } else {
      mark_erased_range(prior.beginning, prior.ending);
    }
  }

  const auto working_triplet_informative_at = [&](std::size_t coordinate) {
    std::array<std::uint8_t, 3> states{};
    for (std::size_t role = 0; role < representatives.size(); ++role) {
      states[role] = erased_for_representative[role][coordinate] != 0
          ? 0
          : alignment_.at(representatives[role], coordinate - 1);
    }
    if (states[0] == 0 || states[1] == 0 || states[2] == 0) return false;
    return (states[0] == states[1] && states[0] != states[2]) ||
        (states[0] == states[2] && states[0] != states[1]) ||
        (states[1] == states[2] && states[1] != states[0]);
  };

  // FindSubSeq2/XPosDiff uses the current (already erased) triplet to map
  // alignment coordinates onto information-rich positions. Retain that
  // source convention rather than measuring raw nucleotide distance.
  auto& informative_coordinates = breakpoint_informative_coordinates_scratch_;
  informative_coordinates.clear();
  informative_coordinates.reserve(alignment_.length);
  for (std::size_t coordinate = 1; coordinate <= alignment_.length; ++coordinate) {
    if (working_triplet_informative_at(coordinate)) {
      informative_coordinates.push_back(coordinate);
    }
  }
  const auto information_index_at = [&](std::size_t coordinate) {
    return static_cast<std::size_t>(std::upper_bound(
        informative_coordinates.begin(),
        informative_coordinates.end(),
        coordinate) - informative_coordinates.begin());
  };
  const auto source_xdiffpos = [&](std::int64_t one_based_index) {
    if (one_based_index <= 0 ||
        one_based_index > static_cast<std::int64_t>(informative_coordinates.size())) {
      return std::size_t{0};
    }
    return informative_coordinates[static_cast<std::size_t>(one_based_index - 1)];
  };
  const auto append_coordinates = [&](std::vector<std::size_t>& coordinates,
                                      std::size_t first,
                                      std::size_t last) {
    first = std::max<std::size_t>(first, 1);
    last = std::min(last, alignment_.length);
    if (first > last) return;
    coordinates.reserve(coordinates.size() + last - first + 1);
    for (std::size_t coordinate = first; coordinate <= last; ++coordinate) {
      coordinates.push_back(coordinate);
    }
  };
  const auto make_native_check_range = [&](std::size_t boundary,
                                            std::size_t center,
                                            std::vector<std::size_t>& coordinates) {
    coordinates.clear();
    if (informative_coordinates.empty()) return;

    std::size_t adjusted_center = center;
    std::size_t center_index = information_index_at(adjusted_center);
    // CheckEndsVB searches forward for a beginning before the first mapped
    // site and backward for an ending. Its XPosDiff values at other invariant
    // positions retain the preceding information-rich index.
    if (center_index == 0) {
      adjusted_center = boundary == 0
          ? informative_coordinates.front()
          : informative_coordinates.back();
      center_index = information_index_at(adjusted_center);
    }

    const std::int64_t window = static_cast<std::int64_t>(options_.window_sites);
    const std::int64_t information_sites =
        static_cast<std::int64_t>(informative_coordinates.size());
    const std::int64_t mapped_center = static_cast<std::int64_t>(center_index);
    if (boundary == 0) {
      std::size_t target = 1;
      if (mapped_center - window > 0) {
        target = source_xdiffpos(mapped_center - window);
      } else if (options_.circular && mapped_center - window + information_sites >= 0) {
        target = source_xdiffpos(mapped_center - window + information_sites);
      }
      if (target < adjusted_center) {
        append_coordinates(coordinates, target, adjusted_center);
      } else {
        append_coordinates(coordinates, target, alignment_.length - 1);
        // Preserve CheckEnds' literal CirF comparison. In normal linear
        // ranges this branch is unreachable; in circular wrap ranges the
        // native code omits the two sentinel boundary coordinates.
        if (!options_.circular) {
          coordinates.push_back(alignment_.length);
          coordinates.push_back(1);
        }
        append_coordinates(coordinates, 2, adjusted_center);
      }
    } else {
      std::size_t target = alignment_.length;
      if (mapped_center + window < information_sites) {
        target = source_xdiffpos(mapped_center + window);
      } else if (options_.circular) {
        target = source_xdiffpos(mapped_center + window - information_sites);
      }
      std::size_t start = 1;
      if (mapped_center - window > 0) {
        start = source_xdiffpos(mapped_center - window);
      } else if (options_.circular && mapped_center - window + information_sites >= 0) {
        start = source_xdiffpos(mapped_center - window + information_sites);
      }
      if (start < 1) start = alignment_.length;
      if (target > start) {
        append_coordinates(coordinates, start, target);
      } else {
        append_coordinates(coordinates, start, alignment_.length - 1);
        if (!options_.circular) {
          coordinates.push_back(alignment_.length);
          coordinates.push_back(1);
        }
        append_coordinates(coordinates, 2, target);
      }
    }
  };
  const auto informative_distance_to_tract = [&](std::size_t center,
                                                  const UniqueEvent& prior,
                                                  std::int64_t direction) {
    if (coordinate_in_tract(
            center, prior.beginning, prior.ending, prior.wraps_origin)) {
      return std::size_t{0};
    }
    std::size_t informative_sites = 0;
    for (std::size_t step = 1; step <= alignment_.length; ++step) {
      const std::size_t coordinate = shifted_coordinate(
          center, direction * static_cast<std::int64_t>(step));
      if (coordinate == 0) break;
      if (coordinate_in_tract(
              coordinate, prior.beginning, prior.ending, prior.wraps_origin)) {
        return informative_sites;
      }
      if (working_triplet_informative_at(coordinate)) {
        ++informative_sites;
        if (informative_sites > options_.window_sites) break;
      }
    }
    return std::numeric_limits<std::size_t>::max();
  };
  const std::array<std::size_t, 2> centers{event.beginning, event.ending};
  for (std::size_t boundary = 0; boundary < centers.size(); ++boundary) {
    const std::size_t center = centers[boundary];
    const std::size_t before = shifted_coordinate(center, -1);
    const std::size_t after = shifted_coordinate(center, 1);
    auto& uncertainty = event.breakpoint_uncertainty[boundary];
    uncertainty.native_check_ends_applied = event.id > 0;
    auto& checked_coordinates = breakpoint_check_coordinates_scratch_[boundary];
    checked_coordinates.clear();
    if (uncertainty.native_check_ends_applied) {
      make_native_check_range(boundary, center, checked_coordinates);
    }
    uncertainty.information_profile_available =
        uncertainty.native_check_ends_applied && !informative_coordinates.empty();
    uncertainty.check_coordinate_count = checked_coordinates.size();
    if (!checked_coordinates.empty()) {
      uncertainty.check_range_beginning = checked_coordinates.front();
      uncertainty.check_range_ending = checked_coordinates.back();
      for (std::size_t index = 1; index < checked_coordinates.size(); ++index) {
        if (checked_coordinates[index] < checked_coordinates[index - 1]) {
          uncertainty.check_range_wraps_origin = true;
          break;
        }
      }
      uncertainty.input_missing_data_in_range = std::any_of(
          checked_coordinates.begin(),
          checked_coordinates.end(),
          [&](std::size_t coordinate) {
            return input_missing_for_triplet[coordinate] != 0;
          });
    }
    const std::size_t original_center_index = information_index_at(center);
    uncertainty.linear_edge_within_window = uncertainty.native_check_ends_applied &&
        !options_.circular &&
        (boundary == 0
             ? original_center_index < options_.window_sites
             : original_center_index + options_.window_sites >
                 informative_coordinates.size());

    for (const std::uint32_t prior_index : relevant_prior_events) {
      const auto& prior = events_[prior_index];
      const bool touches =
          (before != 0 && coordinate_in_tract(
              before, prior.beginning, prior.ending, prior.wraps_origin)) ||
          coordinate_in_tract(center, prior.beginning, prior.ending, prior.wraps_origin) ||
          (after != 0 && coordinate_in_tract(
              after, prior.beginning, prior.ending, prior.wraps_origin));
      if (touches) event.adjacent_erasure_event_ids[boundary].push_back(prior.id);

      const bool erased_in_native_range = std::any_of(
          checked_coordinates.begin(),
          checked_coordinates.end(),
          [&](std::size_t coordinate) {
            return coordinate_in_tract(
                coordinate, prior.beginning, prior.ending, prior.wraps_origin);
          });
      if (erased_in_native_range) {
        event.uncertain_erasure_event_ids[boundary].push_back(prior.id);
        const std::size_t negative_distance = informative_distance_to_tract(
            center, prior, -1);
        const std::size_t nearest_distance = boundary == 0
            ? negative_distance
            : std::min(
                negative_distance,
                informative_distance_to_tract(center, prior, 1));
        const auto current = event.nearest_erasure_informative_sites[boundary];
        if (nearest_distance <= options_.window_sites &&
            (current < 0 || nearest_distance < static_cast<std::size_t>(current))) {
          event.nearest_erasure_informative_sites[boundary] =
              static_cast<std::int32_t>(nearest_distance);
        }
      }
    }
    uncertainty.native_check_ends_warning = uncertainty.native_check_ends_applied &&
        (!uncertainty.information_profile_available ||
        uncertainty.input_missing_data_in_range ||
        uncertainty.linear_edge_within_window ||
        !event.uncertain_erasure_event_ids[boundary].empty());
  }
}

void RdpScanner::refresh_trace_evidence(UniqueEvent& event) {
  event.trace_evidence.clear();
  if (event.major_parent == event.minor_parent) return;
  for (std::size_t sequence = 0; sequence < options_.mask.size(); ++sequence) {
    if (!sequence_masked(static_cast<std::uint32_t>(sequence)) ||
        sequence_disabled(static_cast<std::uint32_t>(sequence)) ||
        sequence == event.major_parent ||
        sequence == event.minor_parent) {
      continue;
    }
    const std::array<std::uint32_t, 3> triplet{
        static_cast<std::uint32_t>(sequence),
        event.major_parent,
        event.minor_parent,
    };
    const auto candidates = triplet_signals(
        triplet, false, nullptr, &profile_scratch_);
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
  // A fragment-assisted anchor must use the same working representatives as
  // the distance/tree evidence. Reported identities remain the originals,
  // while the recheck sees the retained fragment rows that actually exposed
  // the current tract.
  event.maxchi_triplet_recheck = maxchi_triplet_recheck(representatives);
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
    hypothesis.parent_one =
        reported_representatives[kSourceCompRoles[role][0]];
    hypothesis.parent_two =
        reported_representatives[kSourceCompRoles[role][1]];
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
      if (sequence_disabled(static_cast<std::uint32_t>(sequence))) continue;
      std::array<bool, 3> present{};
      for (std::size_t role = 0; role < 3; ++role) {
        const auto& set = event.role_hypotheses[role].detectable_signal_set;
        present[role] = std::binary_search(
            set.begin(), set.end(), static_cast<std::uint32_t>(sequence));
      }
      for (std::size_t role = 0; role < 3; ++role) {
        const std::size_t other_one = kSourceCompRoles[role][0];
        const std::size_t other_two = kSourceCompRoles[role][1];
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
  for (std::size_t sequence = 0; sequence < alignment_.sequence_count(); ++sequence) {
    tree_candidates.push_back(static_cast<std::uint32_t>(sequence));
  }
  for (const std::uint32_t sequence : active_sequences_) {
    if (sequence >= alignment_.sequence_count()) tree_candidates.push_back(sequence);
  }
  sort_unique(tree_candidates);
  const std::vector<std::uint32_t> tree_sequences = select_tree_sequences(
      working_alignment_, tree_candidates, representatives);
  event.tree_panel_sequences = tree_sequences.size();
  event.tree_panel_subsampled = tree_sequences.size() < tree_candidates.size();
  event.tree_panel_leaves.clear();
  event.tree_panel_leaves.reserve(tree_sequences.size());
  for (const std::uint32_t working_sequence : tree_sequences) {
    TreePanelLeafSummary leaf;
    leaf.working_sequence = working_sequence;
    leaf.original_sequence = working_sequence < working_origins_.size()
        ? working_origins_[working_sequence]
        : working_sequence;
    leaf.fragment_event = working_sequence < working_fragment_events_.size()
        ? working_fragment_events_[working_sequence]
        : -1;
    event.tree_panel_leaves.push_back(leaf);
  }
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
    auto& summary = event.tree_regions[region];
    summary = {};
    summary.sites = tree_evidence[region].site_count;
    summary.sequences = tree_evidence[region].sequences.size();
    summary.bootstrap_replicates = tree_evidence[region].bootstrap_replicates;
    summary.supported_internal_branches =
        tree_evidence[region].supported_internal_branches;
    summary.internal_branches = tree_evidence[region].internal_branches;
    summary.node_count = tree_evidence[region].topology_node_count;
    summary.root = tree_evidence[region].topology_root;
    summary.topology_edges.reserve(tree_evidence[region].topology_edges.size());
    for (const auto& edge : tree_evidence[region].topology_edges) {
      summary.topology_edges.push_back({
          edge.first,
          edge.second,
          edge.length,
          edge.bootstrap_support,
          edge.internal,
          edge.collapsed,
      });
    }
    summary.usable = tree_evidence[region].usable;
  }

  std::array<std::vector<std::array<double, 3>>, 6> reference_distances;
  for (std::size_t region = 0; region < reference_distances.size(); ++region) {
    reference_distances[region].resize(alignment_.sequence_count());
    for (std::size_t sequence = 0; sequence < alignment_.sequence_count(); ++sequence) {
      reference_distances[region][sequence].fill(10.0);
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
  const bool distinct_closest_pairs = source_in_list(
      closest_pair[0], closest_pair[1], in_list);
  if (distinct_closest_pairs) {
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
  const SourceCalcMatchGrid calc_match_grid = source_calc_match(
      working_alignment_,
      representatives,
      alignment_.sequence_count(),
      event.beginning,
      event.ending);

  for (std::size_t role = 0; role < representatives.size(); ++role) {
    auto& hypothesis = event.role_hypotheses[role];
    const std::uint32_t analysis_recombinant = representatives[role];
    const std::uint32_t analysis_parent_one =
        representatives[kSourceCompRoles[role][0]];
    const std::uint32_t analysis_parent_two =
        representatives[kSourceCompRoles[role][1]];
    double finaltrim_rank_modifier = 1.0;
    if (!distinct_closest_pairs) {
      const auto closest_members = pair_members(closest_pair[0]);
      if (role != closest_members[0] && role != closest_members[1]) {
        finaltrim_rank_modifier = 0.5;
      }
    }
    const bool finaltrim_suppress_weaker_positive =
        distinct_closest_pairs && role == in_list[0];
    const auto finaltrim_positive_modifier = [&](double score) {
      if (score <= 0.0) return score;
      if (distinct_closest_pairs) {
        if (role == in_list[0] || role == in_list[2]) score /= 2.0;
      } else {
        score *= finaltrim_rank_modifier;
      }
      return score;
    };

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
      if (sequence_disabled(static_cast<std::uint32_t>(sequence))) continue;
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
      evidence.calc_match = calc_match_grid[role][sequence];
      evidence.detectable_support = detectable[sequence] != 0;
      evidence.acceptable_affinity = acceptable_affinity[role][sequence] != 0;
      evidence.warning_filtered = hypothesis.correlation_warnings;
      evidence.breakpoint_overlap_sites = overlap_sites[sequence];
      evidence.overlap_eligible = evidence.breakpoint_overlap_sites[0] > 10 ||
          evidence.breakpoint_overlap_sites[1] > 10;
      auto& finaltrim = evidence.final_trim_matrix;
      finaltrim.applies_to_nonrepresentative =
          sequence != reported_representatives[0] &&
          sequence != reported_representatives[1] &&
          sequence != reported_representatives[2];
      const auto tree_distance = [&](
                                     std::size_t region,
                                     std::uint32_t first,
                                     std::uint32_t second,
                                     bool collapsed) {
        if (tree_evidence[region].usable && tree_evidence[region].contains(first) &&
            tree_evidence[region].contains(second)) {
          return tree_evidence[region].tree(first, second, collapsed);
        }
        finaltrim.tree_distance_fallback = true;
        return reference_distance(region, first, second);
      };
      const auto calc_match_tree_distance = [&](std::size_t region,
                                                 std::uint32_t first,
                                                 std::uint32_t second) {
        if (tree_evidence[region].usable && tree_evidence[region].contains(first) &&
            tree_evidence[region].contains(second)) {
          return tree_evidence[region].tree(first, second, false);
        }
        evidence.calc_match.topology_distance_fallback = true;
        return reference_distance(region, first, second);
      };
      if (evidence.calc_match.available &&
          evidence.calc_match.breakpoint_match_class > -1) {
        const std::size_t other_one = kSourceCompRoles[role][0];
        const std::size_t other_two = kSourceCompRoles[role][1];
        const auto representative = [&](std::size_t member) {
          return representatives[member];
        };
        const auto anchor_closest_on_either_side = [&](std::size_t other) {
          return calc_match_tree_distance(
                     4, analysis_recombinant, evidence.sequence) <=
                  calc_match_tree_distance(
                     4, representative(other), evidence.sequence) ||
              calc_match_tree_distance(
                     5, analysis_recombinant, evidence.sequence) <=
                  calc_match_tree_distance(
                     5, representative(other), evidence.sequence);
        };
        bool topology_consistent =
            anchor_closest_on_either_side(other_one) &&
            anchor_closest_on_either_side(other_two);
        constexpr std::array<std::size_t, 2> whole_regions{{4, 5}};
        for (const std::size_t region : whole_regions) {
          if (!topology_consistent) break;
          const double anchor_first = calc_match_tree_distance(
              region, analysis_recombinant, representative(other_one));
          const double anchor_second = calc_match_tree_distance(
              region, analysis_recombinant, representative(other_two));
          const double candidate_first = calc_match_tree_distance(
              region, representative(other_one), evidence.sequence);
          const double candidate_second = calc_match_tree_distance(
              region, representative(other_two), evidence.sequence);
          if ((anchor_first < anchor_second && candidate_first > candidate_second) ||
              (anchor_first > anchor_second && candidate_first < candidate_second)) {
            topology_consistent = false;
          }
        }
        if (!topology_consistent) {
          evidence.calc_match.breakpoint_match_class = -1;
          evidence.calc_match.topology_filtered = true;
        }
      }
      const auto distance_panel = [&](
                                      std::size_t region,
                                      const auto& distance) {
        FinalTrimDistancePanel panel;
        panel.anchor_candidate = distance(
            region, analysis_recombinant, evidence.sequence);
        panel.anchor_parents = {
            distance(region, analysis_recombinant, analysis_parent_one),
            distance(region, analysis_recombinant, analysis_parent_two),
        };
        panel.candidate_parents = {
            distance(region, evidence.sequence, analysis_parent_one),
            distance(region, evidence.sequence, analysis_parent_two),
        };
        panel.parent_pair = distance(region, analysis_parent_one, analysis_parent_two);
        return panel;
      };
      const auto collapsed_distance = [&](std::size_t region, std::uint32_t first,
                                          std::uint32_t second) {
        return tree_distance(region, first, second, true);
      };
      const auto raw_tree_distance_for_finaltrim = [&](std::size_t region, std::uint32_t first,
                                                       std::uint32_t second) {
        return tree_distance(region, first, second, false);
      };
      const auto jc_distance_for_finaltrim = [&](std::size_t region, std::uint32_t first,
                                                 std::uint32_t second) {
        return reference_distance(region, first, second);
      };
      const std::array<FinalTrimDistancePanel, 2> collapsed_panels{
          distance_panel(4, collapsed_distance),
          distance_panel(5, collapsed_distance),
      };
      const std::array<FinalTrimDistancePanel, 2> raw_tree_panels{
          distance_panel(4, raw_tree_distance_for_finaltrim),
          distance_panel(5, raw_tree_distance_for_finaltrim),
      };
      for (const auto& panel : collapsed_panels) {
        finaltrim.collapsed_tree_position_score += source_finaltrim_rank_contribution(
            panel.anchor_candidate,
            panel.anchor_parents[0],
            panel.anchor_parents[1],
            finaltrim_rank_modifier,
            2.0);
        finaltrim.collapsed_tree_position_score += source_finaltrim_rank_contribution(
            panel.anchor_candidate,
            panel.candidate_parents[0],
            panel.candidate_parents[1],
            finaltrim_rank_modifier,
            2.0);
      }
      for (std::size_t side = 0; side < raw_tree_panels.size(); ++side) {
        const auto& panel = raw_tree_panels[side];
        finaltrim.raw_tree_position_score += source_finaltrim_rank_contribution(
            panel.anchor_candidate,
            panel.anchor_parents[0],
            panel.anchor_parents[1],
            finaltrim_rank_modifier,
            2.0);
        finaltrim.raw_tree_position_score += source_finaltrim_rank_contribution(
            panel.anchor_candidate,
            panel.candidate_parents[0],
            panel.candidate_parents[1],
            finaltrim_rank_modifier,
            1.0);
      }
      const std::array<FinalTrimDistancePanel, 2> whole_distance_panels{
          distance_panel(4, jc_distance_for_finaltrim),
          distance_panel(5, jc_distance_for_finaltrim),
      };
      finaltrim.relative_distance_score = finaltrim_positive_modifier(
          source_finaltrim_distance_score(
              whole_distance_panels[0],
              whole_distance_panels[1],
              false,
              false,
              8.0,
              0.5,
              finaltrim_suppress_weaker_positive));
      constexpr std::array<std::array<std::size_t, 2>, 2> breakpoint_regions{{
          {{0, 1}},
          {{2, 3}},
      }};
      for (std::size_t boundary = 0; boundary < breakpoint_regions.size(); ++boundary) {
        const auto first_panel = distance_panel(
            breakpoint_regions[boundary][0], jc_distance_for_finaltrim);
        const auto second_panel = distance_panel(
            breakpoint_regions[boundary][1], jc_distance_for_finaltrim);
        if (hypothesis.correlation_warnings[boundary] != 0 ||
            first_panel.anchor_candidate >= 3.0 ||
            second_panel.anchor_candidate >= 3.0) {
          continue;
        }
        finaltrim.breakpoint_score_available[boundary] = 1;
        finaltrim.breakpoint_distance_scores[boundary] = finaltrim_positive_modifier(
            source_finaltrim_distance_score(
                first_panel,
                second_panel,
                true,
                true,
                6.0,
                1.0,
                finaltrim_suppress_weaker_positive));
      }
      finaltrim.active_consensus_matrix_score =
          (finaltrim.collapsed_tree_position_score +
           finaltrim.raw_tree_position_score) /
              2.0 +
          finaltrim.relative_distance_score +
          finaltrim.breakpoint_distance_scores[0] +
          finaltrim.breakpoint_distance_scores[1];
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
      // Retain every candidate's fixed-size late matrix diagnostics. The
      // source evaluates OKSeq matrix families for all ordinary sequences,
      // including rows that are not ultimately admitted to an RList.
      hypothesis.distance_evidence.push_back(evidence);
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
      evidence.disabled_excluded = sequence_disabled(evidence.sequence);
      evidence.distance_fallback = !tree_evidence[0].contains(evidence.sequence) ||
          std::any_of(
              tree_evidence.begin(),
              tree_evidence.end(),
              [](const TreeRegionEvidence& region) { return !region.usable; });
      if (evidence.disabled_excluded) {
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
      if (sequence_disabled(candidate)) continue;
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
  // this as separately auditable evidence before the later active FinalTrim
  // stages mutate the lists that feed the browser's two-of-three group.
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
          // FinalTrim clears this pair for every occurrence of the sequence in
          // an RList once the pair is duplicated, including an occurrence whose
          // own copied value was already zeroed by a warning or inversion.
          if (evidence.sequence == sequence && evidence.significant) {
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

  // FindActualEvents + MakeMatchMatX2P supply the final active FinalTrim
  // matrix family (OKSeq 14). The temporary native RLists still contain
  // inverse-only rows here, so reconstruct that pre-StripDupInv membership
  // rather than using only the cleaned distance-correlation sets.
  const std::size_t ordinary_sequence_count = alignment_.sequence_count();
  std::array<std::vector<std::uint8_t>, 3> initial_rlist;
  std::array<std::vector<std::uint8_t>, 3> initial_inverse_list;
  for (std::size_t role = 0; role < 3; ++role) {
    initial_rlist[role].assign(ordinary_sequence_count, 0);
    initial_inverse_list[role].assign(ordinary_sequence_count, 0);
    for (const auto& evidence : event.role_hypotheses[role].distance_evidence) {
      if (evidence.sequence >= ordinary_sequence_count) continue;
      if (evidence.significant || evidence.stripped_inverse_only) {
        initial_rlist[role][evidence.sequence] = 1;
      }
      if (evidence.stripped_inverse_only) {
        initial_inverse_list[role][evidence.sequence] = 1;
      }
    }
  }

  SourceDetectedTractGrid detected_tracts;
  for (auto& role : detected_tracts) role.resize(ordinary_sequence_count);
  // Every direct candidate + two-representative signal above this routine's
  // >1/3 overlap threshold also passes event support's >0.3 shared-pair gate,
  // so the retained support catalogue is complete for this lookup and remains
  // reproducible after project restore and later cyclic rounds.
  std::vector<std::uint32_t> catalogue_signal_ids;
  catalogue_signal_ids.reserve(event.support_signal_ids.size());
  for (const std::uint32_t signal_id : event.support_signal_ids) {
    if (signal_id < signals_.size() &&
        signals_[signal_id].review_state != ReviewState::rejected) {
      catalogue_signal_ids.push_back(signal_id);
    }
  }
  std::stable_sort(
      catalogue_signal_ids.begin(),
      catalogue_signal_ids.end(),
      [&](std::uint32_t left, std::uint32_t right) {
        const auto& first = signals_[left];
        const auto& second = signals_[right];
        if (first.recombinant != second.recombinant) {
          return first.recombinant < second.recombinant;
        }
        return left < right;
      });

  for (std::size_t role = 0; role < 3; ++role) {
    const std::size_t comp_zero = kSourceCompRoles[role][0];
    const std::size_t comp_one = kSourceCompRoles[role][1];
    const std::uint32_t parent_zero = reported_representatives[comp_zero];
    const std::uint32_t parent_one = reported_representatives[comp_one];
    // FindActualEvents rejects the whole role permutation if either fixed
    // parental representative entered its other native RList by inversion.
    if (parent_zero >= ordinary_sequence_count ||
        parent_one >= ordinary_sequence_count ||
        initial_inverse_list[comp_zero][parent_zero] != 0 ||
        initial_inverse_list[comp_one][parent_one] != 0) {
      continue;
    }
    for (std::size_t sequence = 0; sequence < ordinary_sequence_count; ++sequence) {
      if (initial_rlist[role][sequence] == 0) continue;
      const auto expected = canonical_triplet({
          static_cast<std::uint32_t>(sequence), parent_zero, parent_one});
      auto& selected = detected_tracts[role][sequence];
      for (const std::uint32_t signal_id : catalogue_signal_ids) {
        if (signal_id >= signals_.size()) continue;
        const auto& signal = signals_[signal_id];
        if (canonical_triplet(signal.triplet) != expected) continue;
        const double overlap = tract_overlap(
            event.beginning,
            event.ending,
            signal.beginning,
            signal.ending);
        if (overlap * 3.0 <= 1.0 || overlap <= selected.overlap) continue;
        selected.available = true;
        selected.overlap = overlap;
        selected.beginning = signal.beginning;
        selected.ending = signal.ending;
        selected.signal_id = static_cast<std::int32_t>(signal.id);
      }
    }
  }

  using SourceMatchMatrix =
      std::array<std::array<std::vector<double>, 3>, 3>;
  SourceMatchMatrix match_matrix;
  for (auto& role : match_matrix) {
    for (auto& reference : role) reference.assign(ordinary_sequence_count, 0.0);
  }
  std::vector<std::size_t> valid_prefix(working_alignment_.length + 1, 0);
  std::vector<std::size_t> difference_prefix(working_alignment_.length + 1, 0);
  for (std::size_t reference_role = 0; reference_role < 3; ++reference_role) {
    const std::uint32_t reported_reference =
        reported_representatives[reference_role];
    const std::uint32_t working_reference = representatives[reference_role];
    for (std::size_t sequence = 0; sequence < ordinary_sequence_count; ++sequence) {
      bool needs_sequence_prefix = false;
      for (std::size_t role = 0; role < 3; ++role) {
        if (reported_reference < ordinary_sequence_count &&
            detected_tracts[role][reported_reference].available) {
          needs_sequence_prefix = true;
          break;
        }
      }
      if (needs_sequence_prefix) {
        valid_prefix[0] = 0;
        difference_prefix[0] = 0;
        for (std::size_t position = 0; position < working_alignment_.length; ++position) {
          const std::uint8_t first = working_alignment_.at(working_reference, position);
          const std::uint8_t second = working_alignment_.at(sequence, position);
          const bool valid = first != 0 && second != 0;
          valid_prefix[position + 1] = valid_prefix[position] + (valid ? 1 : 0);
          difference_prefix[position + 1] =
              difference_prefix[position] + (valid && first != second ? 1 : 0);
        }
      }
      for (std::size_t role = 0; role < 3; ++role) {
        const auto& reference_tract = detected_tracts[role][reported_reference];
        const auto& candidate_tract = detected_tracts[role][sequence];
        if (!reference_tract.available) {
          if (!candidate_tract.available) {
            // The MakeMatchMat fallback is the inside-tract SMatSmall entry.
            match_matrix[role][reference_role][sequence] = reference_distance(
                5,
                static_cast<std::uint32_t>(sequence),
                working_reference);
          }
          continue;
        }
        const SourceIntervalSegments reference_segments = source_interval_segments(
            reference_tract.beginning,
            reference_tract.ending,
            working_alignment_.length,
            true);
        SourceIntervalSegments candidate_segments;
        const SourceIntervalSegments* candidate_segments_pointer = nullptr;
        if (candidate_tract.available) {
          candidate_segments = source_interval_segments(
              candidate_tract.beginning,
              candidate_tract.ending,
              working_alignment_.length,
              false);
          candidate_segments_pointer = &candidate_segments;
        }
        match_matrix[role][reference_role][sequence] = source_match_matrix_distance(
            valid_prefix,
            difference_prefix,
            reference_segments,
            candidate_segments_pointer);
      }
    }
  }

  const auto pattern_scores = source_pattern_scores(
      working_alignment_,
      representatives,
      reported_representatives,
      ordinary_sequence_count,
      layout,
      event.beginning,
      event.ending);

  // FinalTrim OKSeq 6 is the membership snapshot from its first, iterative
  // nearest-nonrecombinant pass. It starts from the post-StripDupInv RLists,
  // clears entries removed by either immediate distance/unfound tests or the
  // paired breakpoint-distance veto, and intentionally does not mark entries
  // added by the later distance-only expansions.
  using SourceCleanedCorrelationGrid =
      std::array<std::vector<std::array<double, 3>>, 3>;
  SourceCleanedCorrelationGrid cleaned_correlations;
  SourceCleanedCorrelationGrid raw_correlations;
  for (auto& role : cleaned_correlations) role.resize(ordinary_sequence_count);
  for (auto& role : raw_correlations) role.resize(ordinary_sequence_count);
  for (std::size_t role = 0; role < 3; ++role) {
    for (const auto& evidence : event.role_hypotheses[role].distance_evidence) {
      if (evidence.sequence >= ordinary_sequence_count) continue;
      for (std::size_t pair = 0; pair < 3; ++pair) {
        raw_correlations[role][evidence.sequence][pair] =
            evidence.correlations[pair];
        if (evidence.warning_filtered[pair] == 0 &&
            evidence.inversion_codes[pair] == 0 &&
            evidence.duplicate_filtered[pair] == 0) {
          cleaned_correlations[role][evidence.sequence][pair] =
              evidence.correlations[pair];
        }
      }
    }
  }

  std::array<std::vector<std::uint32_t>, 3> nearest_lists;
  std::array<std::vector<std::uint8_t>, 3> nearest_membership;
  for (std::size_t role = 0; role < 3; ++role) {
    nearest_membership[role].assign(ordinary_sequence_count, 0);
    // MakeRList appends in ascending sequence order. Reproduce StripDupInv's
    // swap-last removals instead of merely filtering, because FinalTrim later
    // mutates the list in that retained order.
    for (std::size_t sequence = 0; sequence < ordinary_sequence_count; ++sequence) {
      if (initial_rlist[role][sequence] != 0) {
        nearest_lists[role].push_back(static_cast<std::uint32_t>(sequence));
      }
    }
    std::size_t position = 0;
    while (position < nearest_lists[role].size()) {
      const std::uint32_t sequence = nearest_lists[role][position];
      if (initial_inverse_list[role][sequence] == 0) {
        ++position;
        continue;
      }
      nearest_lists[role][position] = nearest_lists[role].back();
      nearest_lists[role].pop_back();
    }
    for (const std::uint32_t sequence : nearest_lists[role]) {
      nearest_membership[role][sequence] = 1;
    }
  }

  const auto finaltrim_tree_distance = [&](std::size_t region,
                                            std::uint32_t first,
                                            std::uint32_t second,
                                            bool collapsed) {
    if (tree_evidence[region].usable && tree_evidence[region].contains(first) &&
        tree_evidence[region].contains(second)) {
      return tree_evidence[region].tree(first, second, collapsed);
    }
    return reference_distance(region, first, second);
  };
  const auto finaltrim_breakpoint_veto = [&](std::size_t role,
                                             std::uint32_t sequence) {
    const std::size_t parent_zero_role = kSourceCompRoles[role][0];
    const std::size_t parent_one_role = kSourceCompRoles[role][1];
    const std::uint32_t anchor = representatives[role];
    const std::uint32_t parent_zero = representatives[parent_zero_role];
    const std::uint32_t parent_one = representatives[parent_one_role];
    const auto region_score = [&](std::size_t region) {
      const double candidate =
          reference_distance(region, anchor, sequence);
      if (candidate > reference_distance(region, anchor, parent_zero) ||
          candidate > reference_distance(region, anchor, parent_one)) {
        if (candidate > reference_distance(region, parent_zero, sequence) ||
            candidate > reference_distance(region, parent_one, sequence)) {
          return 2;
        }
        return 1;
      }
      return -1;
    };

    int go_on = 0;
    if (event.role_hypotheses[role].correlation_warnings[0] == 0) {
      if (cleaned_correlations[role][sequence][0] > 0.95) {
        go_on += region_score(0);
        go_on += region_score(1);
      } else {
        go_on = 1;
      }
    } else {
      go_on = 1;
    }
    if (go_on > 0) {
      if (event.role_hypotheses[role].correlation_warnings[1] == 0) {
        go_on = 0;
        if (cleaned_correlations[role][sequence][1] > 0.95) {
          go_on += region_score(2);
          go_on += region_score(3);
        } else {
          go_on = 1;
        }
      } else {
        go_on = 1;
      }
    }
    return go_on > 0;
  };

  std::array<std::int64_t, 3> previous_last_index{{0, 0, 0}};
  const double source_neighbour_cutoff = ordinary_sequence_count > 0
      ? static_cast<double>(ordinary_sequence_count - 1) * 2.0 / 1000.0
      : 0.0;
  for (std::size_t role = 0; role < 3; ++role) {
    for (;;) {
      const std::int64_t current_last_index = nearest_lists[role].empty()
          ? -1
          : static_cast<std::int64_t>(nearest_lists[role].size() - 1);
      // Preserve the source's oRNum initialization: a one-entry RList has
      // last index zero and skips this pass entirely.
      if (previous_last_index[role] == current_last_index) break;
      previous_last_index[role] = current_last_index;

      std::array<double, 2> farthest_recombinant_distance{};
      std::vector<std::uint8_t> seen(ordinary_sequence_count, 0);
      for (const std::uint32_t sequence : nearest_lists[role]) {
        seen[sequence] = 1;
        farthest_recombinant_distance[0] = std::max(
            farthest_recombinant_distance[0],
            finaltrim_tree_distance(
                4, representatives[role], sequence, false));
        farthest_recombinant_distance[1] = std::max(
            farthest_recombinant_distance[1],
            finaltrim_tree_distance(
                5, representatives[role], sequence, false));
      }

      std::vector<std::uint32_t> nonrecombinant_candidates;
      for (std::size_t sequence = 0; sequence < ordinary_sequence_count; ++sequence) {
        if (sequence_disabled(static_cast<std::uint32_t>(sequence))) continue;
        if (seen[sequence] == 0 &&
            finaltrim_tree_distance(
                4,
                representatives[role],
                static_cast<std::uint32_t>(sequence),
                false) <= farthest_recombinant_distance[0]) {
          nonrecombinant_candidates.push_back(static_cast<std::uint32_t>(sequence));
          seen[sequence] = 1;
        }
        if (seen[sequence] == 0 &&
            finaltrim_tree_distance(
                5,
                representatives[role],
                static_cast<std::uint32_t>(sequence),
                false) <= farthest_recombinant_distance[1]) {
          nonrecombinant_candidates.push_back(static_cast<std::uint32_t>(sequence));
          seen[sequence] = 1;
        }
      }

      std::array<std::vector<std::uint8_t>, 2> mark_remove;
      for (auto& marks : mark_remove) marks.assign(ordinary_sequence_count, 0);
      if (!nonrecombinant_candidates.empty()) {
        const auto prune_side = [&](std::size_t side) {
          double nearest_nonrecombinant_distance = 1000000.0;
          for (const std::uint32_t sequence : nonrecombinant_candidates) {
            const double collapsed_distance = finaltrim_tree_distance(
                4 + side, representatives[role], sequence, true);
            bool below_source_scale = collapsed_distance < source_neighbour_cutoff;
            if (side == 0) {
              below_source_scale = finaltrim_tree_distance(
                  4, representatives[role], sequence, false) <
                  source_neighbour_cutoff;
            }
            if ((side != 0 || ordinary_sequence_count > 3) &&
                collapsed_distance < nearest_nonrecombinant_distance &&
                below_source_scale) {
              nearest_nonrecombinant_distance = collapsed_distance;
            }
          }

          const std::size_t parent_zero_role = kSourceCompRoles[role][0];
          const std::size_t parent_one_role = kSourceCompRoles[role][1];
          std::size_t index = 0;
          while (index < nearest_lists[role].size()) {
            const std::uint32_t sequence = nearest_lists[role][index];
            const double first_correlation =
                cleaned_correlations[role][sequence][0];
            const double second_correlation =
                cleaned_correlations[role][sequence][1];
            if (first_correlation >= 0.99 || second_correlation >= 0.99) {
              ++index;
              continue;
            }
            const bool moderate_correlation =
                (first_correlation < 0.99 && first_correlation > 0.83) ||
                (second_correlation < 0.99 && second_correlation > 0.83);
            if (!moderate_correlation) {
              ++index;
              continue;
            }

            const double candidate_distance = finaltrim_tree_distance(
                4 + side, representatives[role], sequence, true);
            const double first_parent_distance = finaltrim_tree_distance(
                4 + side,
                representatives[role],
                representatives[parent_zero_role],
                true);
            const double second_parent_distance = finaltrim_tree_distance(
                4 + side,
                representatives[role],
                representatives[parent_one_role],
                true);
            if (candidate_distance <= nearest_nonrecombinant_distance &&
                candidate_distance <= first_parent_distance &&
                candidate_distance <= second_parent_distance) {
              ++index;
              continue;
            }

            const bool strong_ticket =
                first_correlation > 0.95 || second_correlation > 0.95;
            const bool unfound = !detected_tracts[role][sequence].available;
            if (unfound || !strong_ticket) {
              nearest_membership[role][sequence] = 0;
              nearest_lists[role][index] = nearest_lists[role].back();
              nearest_lists[role].pop_back();
              continue;
            }
            if (finaltrim_breakpoint_veto(role, sequence)) {
              mark_remove[side][sequence] = 1;
            }
            ++index;
          }
        };

        prune_side(0);
        prune_side(1);
      }

      std::size_t index = 0;
      while (index < nearest_lists[role].size()) {
        const std::uint32_t sequence = nearest_lists[role][index];
        if (mark_remove[0][sequence] != 0 && mark_remove[1][sequence] != 0) {
          nearest_membership[role][sequence] = 0;
          nearest_lists[role][index] = nearest_lists[role].back();
          nearest_lists[role].pop_back();
          continue;
        }
        ++index;
      }
    }
  }

  for (std::size_t role = 0; role < 3; ++role) {
    const std::size_t comp_zero = kSourceCompRoles[role][0];
    const std::size_t comp_one = kSourceCompRoles[role][1];
    double repeated_pair_modifier = 1.0;
    if (!distinct_closest_pairs) {
      const auto closest_members = pair_members(closest_pair[0]);
      if (role != closest_members[0] && role != closest_members[1]) {
        repeated_pair_modifier = 0.5;
      }
    }
    for (auto& evidence : event.role_hypotheses[role].distance_evidence) {
      if (evidence.sequence >= ordinary_sequence_count) continue;
      auto& finaltrim = evidence.final_trim_matrix;
      const auto& detected = detected_tracts[role][evidence.sequence];
      finaltrim.detected_event_match = detected.available;
      finaltrim.detected_event_overlap = detected.overlap;
      finaltrim.detected_event_beginning = detected.beginning;
      finaltrim.detected_event_ending = detected.ending;
      finaltrim.detected_event_signal_id = detected.signal_id;
      finaltrim.detected_region_match_distances = {
          match_matrix[role][role][evidence.sequence],
          match_matrix[role][role][reported_representatives[comp_zero]],
          // Active Module2.bas uses CompMat(WinPP, 1), without ISeqs, in this
          // one positive-branch comparison. Preserve that sequence-index quirk.
          match_matrix[role][role][comp_one],
          match_matrix[role][role][reported_representatives[comp_one]],
          match_matrix[role][comp_zero][evidence.sequence],
          match_matrix[role][comp_one][evidence.sequence],
          match_matrix[role][comp_one][reported_representatives[comp_zero]],
      };
      finaltrim.detected_region_saturated =
          finaltrim.detected_region_match_distances[0] == 3.0;
      finaltrim.detected_region_distance_score =
          source_finaltrim_detected_region_score(
              finaltrim.detected_region_match_distances,
              distinct_closest_pairs,
              role,
              in_list,
              repeated_pair_modifier);
      finaltrim.active_consensus_matrix_score +=
          finaltrim.detected_region_distance_score;

      auto& consensus = evidence.consensus_score;
      const double source_sequence_multiplier =
          ordinary_sequence_count > 0
          ? static_cast<double>(ordinary_sequence_count - 1)
          : 0.0;
      for (std::size_t pair = 0; pair < evidence.p_values.size(); ++pair) {
        if (evidence.warning_filtered[pair] != 0 ||
            evidence.inversion_codes[pair] != 0) {
          continue;
        }
        const double p_value = evidence.p_values[pair];
        if (!(p_value > 0.0)) continue;
        const double corrected = p_value * source_sequence_multiplier;
        if (consensus.corrected_correlation_p_value == 0.0 ||
            corrected < consensus.corrected_correlation_p_value) {
          consensus.corrected_correlation_p_value = corrected;
        }
      }
      consensus.detected_event_overlap = detected.overlap;
      consensus.detectable_set_member = std::binary_search(
          event.role_hypotheses[role].detectable_signal_set.begin(),
          event.role_hypotheses[role].detectable_signal_set.end(),
          evidence.sequence);
      consensus.pattern_score = pattern_scores[evidence.sequence][role];
      consensus.initial_rlist_member = initial_rlist[role][evidence.sequence] != 0;
      consensus.duplicate_cleaned_member = evidence.duplicate_cleaned_support;
      consensus.nearest_nonrecombinant_member =
          nearest_membership[role][evidence.sequence] != 0;
      for (std::size_t pair = 0; pair < evidence.correlations.size(); ++pair) {
        if (evidence.warning_filtered[pair] == 0 &&
            evidence.inversion_codes[pair] == 0) {
          consensus.maximum_direct_correlation = std::max(
              consensus.maximum_direct_correlation,
              evidence.correlations[pair]);
        }
      }
      consensus.representative_sentinel =
          evidence.sequence == reported_representatives[role];
      if (consensus.representative_sentinel) {
        consensus.base_score_before_final_membership = 1000.0;
        consensus.source_long_matrix_multiplier = 1.0;
        continue;
      }
      consensus.other_representative_zero =
          evidence.sequence == reported_representatives[0] ||
          evidence.sequence == reported_representatives[1] ||
          evidence.sequence == reported_representatives[2];
      if (consensus.other_representative_zero) {
        // ConsensusOK leaves the two non-own representative cells at their
        // zero-initialized value; only ISeqs(role) receives the 1000 sentinel.
        continue;
      }

      double base_score = 0.0;
      if (consensus.corrected_correlation_p_value > 0.0 &&
          consensus.corrected_correlation_p_value < 1.0) {
        base_score +=
            -std::log10(consensus.corrected_correlation_p_value) * 20.0;
      }
      base_score += consensus.detected_event_overlap * 5.0;
      base_score += consensus.pattern_score * 4.0;
      if (consensus.detectable_set_member && consensus.initial_rlist_member) {
        base_score *= 1.2;
      } else if (consensus.detectable_set_member || consensus.initial_rlist_member) {
        base_score *= 1.1;
      }
      if (consensus.duplicate_cleaned_member) base_score *= 1.1;
      if (consensus.nearest_nonrecombinant_member) base_score *= 1.1;
      consensus.base_score_before_final_membership = base_score;

      // Preserve the less-obvious source behavior: NS is declared As Long.
      // Each assignment rounds to an even-tied integer, including each added
      // OKSeq component. The three-statement NS<1 branch consequently ends at
      // zero after assigning 2^-1 back into a Long.
      std::int64_t source_ns = source_vb_long(
          (finaltrim.collapsed_tree_position_score +
           finaltrim.raw_tree_position_score) /
          2.0);
      const std::array<double, 6> remaining_matrix_scores{
          finaltrim.relative_distance_score,
          0.0,
          0.0,
          finaltrim.breakpoint_distance_scores[0],
          finaltrim.breakpoint_distance_scores[1],
          finaltrim.detected_region_distance_score,
      };
      for (const double score : remaining_matrix_scores) {
        if (score != 0.0) {
          source_ns = source_vb_long(static_cast<double>(source_ns) + score);
        }
      }
      if (source_ns < 1) {
        source_ns = source_vb_long(-0.5);
        source_ns = source_vb_long(1.0 - static_cast<double>(source_ns));
        source_ns = source_vb_long(
            std::pow(2.0, -static_cast<double>(source_ns)));
      }
      consensus.source_long_matrix_multiplier = static_cast<double>(source_ns);
    }
  }

  // Complete the active RFF=0 FinalTrim path before the browser exposes a
  // co-recombinant group. The supplied routine mutates each RList in-place:
  // two ascending expansion passes surround the OKSeq 7-14 calculations,
  // followed by chosen-role pruning, OKSeq 15, and ConsensusOK's three-pass
  // list rebuild. Preserve that mutation order, including swap-last removals.
  std::array<std::unordered_map<std::uint64_t, double>, 2> whole_direct_cache;
  std::array<std::unordered_map<std::uint64_t, double>, 2> whole_raw_cache;
  std::array<std::unordered_map<std::uint64_t, double>, 2> whole_collapsed_cache;
  const auto source_direct_whole = [&](std::size_t side,
                                       std::uint32_t first,
                                       std::uint32_t second) {
    if (first == second) return 0.0;
    auto& cache = whole_direct_cache[side];
    const std::uint64_t key = sequence_pair_key(first, second);
    const auto found = cache.find(key);
    if (found != cache.end()) return found->second;
    const double value = reference_distance(4 + side, first, second);
    cache.emplace(key, value);
    return value;
  };
  const auto source_raw_whole = [&](std::size_t side,
                                    std::uint32_t first,
                                    std::uint32_t second) {
    if (first == second) return 0.0;
    auto& cache = whole_raw_cache[side];
    const std::uint64_t key = sequence_pair_key(first, second);
    const auto found = cache.find(key);
    if (found != cache.end()) return found->second;
    const std::size_t region = 4 + side;
    const double value = tree_evidence[region].usable &&
            tree_evidence[region].contains(first) &&
            tree_evidence[region].contains(second)
        ? tree_evidence[region].tree(first, second, false)
        : source_direct_whole(side, first, second);
    cache.emplace(key, value);
    return value;
  };
  const auto source_collapsed_whole = [&](std::size_t side,
                                          std::uint32_t first,
                                          std::uint32_t second) {
    if (first == second) return 0.0;
    auto& cache = whole_collapsed_cache[side];
    const std::uint64_t key = sequence_pair_key(first, second);
    const auto found = cache.find(key);
    if (found != cache.end()) return found->second;
    const std::size_t region = 4 + side;
    const double value = tree_evidence[region].usable &&
            tree_evidence[region].contains(first) &&
            tree_evidence[region].contains(second)
        ? tree_evidence[region].tree(first, second, true)
        : source_direct_whole(side, first, second);
    cache.emplace(key, value);
    return value;
  };
  const auto direct_small = [&](std::size_t side,
                                std::size_t role,
                                std::uint32_t sequence) {
    return source_direct_whole(side, representatives[role], sequence);
  };
  const auto raw_small = [&](std::size_t side,
                             std::size_t role,
                             std::uint32_t sequence) {
    return source_raw_whole(side, representatives[role], sequence);
  };
  const auto collapsed_small = [&](std::size_t side,
                                   std::size_t role,
                                   std::uint32_t sequence) {
    return source_collapsed_whole(side, representatives[role], sequence);
  };

  std::array<std::vector<std::uint32_t>, 3> finaltrim_lists = nearest_lists;
  std::array<std::vector<std::uint8_t>, 3> first_expansion_added;
  std::array<std::vector<std::uint8_t>, 3> second_expansion_added;
  std::array<std::vector<std::uint8_t>, 3> selected_role_pruned;
  std::array<std::vector<std::uint8_t>, 3> finaltrim_membership;
  for (std::size_t role = 0; role < 3; ++role) {
    first_expansion_added[role].assign(ordinary_sequence_count, 0);
    second_expansion_added[role].assign(ordinary_sequence_count, 0);
    selected_role_pruned[role].assign(ordinary_sequence_count, 0);
    finaltrim_membership[role].assign(ordinary_sequence_count, 0);

    std::vector<std::uint8_t> already_in(ordinary_sequence_count, 0);
    std::array<double, 2> farthest_recombinant_distance{};
    for (const std::uint32_t sequence : finaltrim_lists[role]) {
      if (sequence >= ordinary_sequence_count) continue;
      already_in[sequence] = 1;
      for (std::size_t side = 0; side < 2; ++side) {
        farthest_recombinant_distance[side] = std::max(
            farthest_recombinant_distance[side],
            raw_small(side, role, sequence));
      }
    }

    std::vector<std::uint32_t> nonrecombinant_candidates;
    std::vector<std::uint8_t> seen = already_in;
    for (std::size_t sequence = 0; sequence < ordinary_sequence_count; ++sequence) {
      const auto candidate = static_cast<std::uint32_t>(sequence);
      if (sequence_disabled(candidate)) continue;
      if (seen[sequence] == 0 &&
          raw_small(0, role, candidate) <= farthest_recombinant_distance[0]) {
        nonrecombinant_candidates.push_back(candidate);
        seen[sequence] = 1;
      }
      if (seen[sequence] == 0 &&
          raw_small(1, role, candidate) <= farthest_recombinant_distance[1]) {
        nonrecombinant_candidates.push_back(candidate);
        seen[sequence] = 1;
      }
    }

    std::array<double, 2> nearest_nonrecombinant_distance{};
    if (!nonrecombinant_candidates.empty()) {
      nearest_nonrecombinant_distance.fill(1000000.0);
      for (const std::uint32_t sequence : nonrecombinant_candidates) {
        const double outside_collapsed = collapsed_small(0, role, sequence);
        // The source guards only the outside-tree search with NextNo > 2 and
        // applies its scale test to the uncollapsed outside matrix.
        if (ordinary_sequence_count > 3 &&
            outside_collapsed < nearest_nonrecombinant_distance[0] &&
            raw_small(0, role, sequence) < source_neighbour_cutoff) {
          nearest_nonrecombinant_distance[0] = outside_collapsed;
        }
        const double inside_collapsed = collapsed_small(1, role, sequence);
        if (inside_collapsed < nearest_nonrecombinant_distance[1] &&
            inside_collapsed < source_neighbour_cutoff) {
          nearest_nonrecombinant_distance[1] = inside_collapsed;
        }
      }
    }

    for (const std::uint32_t representative : reported_representatives) {
      if (representative < ordinary_sequence_count) already_in[representative] = 1;
    }
    for (std::size_t sequence = 0; sequence < ordinary_sequence_count; ++sequence) {
      if (already_in[sequence] != 0) continue;
      const auto candidate = static_cast<std::uint32_t>(sequence);
      if (sequence_disabled(candidate)) continue;
      if (collapsed_small(0, role, candidate) <=
              nearest_nonrecombinant_distance[0] &&
          collapsed_small(1, role, candidate) <=
              nearest_nonrecombinant_distance[1] &&
          (cleaned_correlations[role][sequence][0] > 0.83 ||
           cleaned_correlations[role][sequence][1] > 0.83 ||
           cleaned_correlations[role][sequence][2] > 0.83)) {
        finaltrim_lists[role].push_back(candidate);
        already_in[sequence] = 1;
        first_expansion_added[role][sequence] = 1;
      }
    }

    const std::size_t parent_zero_role = kSourceCompRoles[role][0];
    const std::size_t parent_one_role = kSourceCompRoles[role][1];
    const double outside_parent_zero = raw_small(
        0, role, representatives[parent_zero_role]);
    const double outside_parent_one = raw_small(
        0, role, representatives[parent_one_role]);
    const double inside_parent_zero = raw_small(
        1, role, representatives[parent_zero_role]);
    const double inside_parent_one = raw_small(
        1, role, representatives[parent_one_role]);
    for (std::size_t sequence = 0; sequence < ordinary_sequence_count; ++sequence) {
      if (already_in[sequence] != 0) continue;
      const auto candidate = static_cast<std::uint32_t>(sequence);
      if (sequence_disabled(candidate)) continue;
      if (raw_small(0, role, candidate) < outside_parent_zero &&
          raw_small(0, role, candidate) < outside_parent_one &&
          raw_small(1, role, candidate) < inside_parent_zero &&
          raw_small(1, role, candidate) < inside_parent_one) {
        finaltrim_lists[role].push_back(candidate);
        already_in[sequence] = 1;
        second_expansion_added[role][sequence] = 1;
      }
    }
  }

  // RWinPP is the selected recombinant role in the active call path. At this
  // point role zero is the event's current selected role; a later role vote
  // can reorder the event and rerun this routine. UBFC is represented by two
  // usable collapsed whole-tract trees, with the documented JC fallback for
  // capped-out candidates.
  const bool collapsed_trim_available =
      tree_evidence[4].usable && tree_evidence[5].usable;
  const auto selected_trim_distance = [&](std::size_t side,
                                          std::size_t role,
                                          std::uint32_t sequence) {
    return collapsed_trim_available
        ? collapsed_small(side, role, sequence)
        : raw_small(side, role, sequence);
  };
  const auto lacks_free_correlation_ticket = [&](std::size_t role,
                                                  std::uint32_t sequence) {
    return raw_correlations[role][sequence][0] < 0.99 &&
        raw_correlations[role][sequence][1] < 0.99 &&
        raw_correlations[role][sequence][2] < 0.99;
  };
  if (distinct_closest_pairs) {
    std::size_t position = 0;
    const std::size_t first_role = in_list[0];
    while (position < finaltrim_lists[first_role].size()) {
      const std::uint32_t sequence = finaltrim_lists[first_role][position];
      const double outside_value = selected_trim_distance(0, first_role, sequence);
      const double outside_limit = selected_trim_distance(
          0, first_role, representatives[in_list[1]]);
      const double inside_value = selected_trim_distance(1, first_role, sequence);
      // Preserve the active collapsed-path index spelling from Module2.bas:
      // SCMatSmall(INList(1), ISeqs(INList(0))). The matrix is symmetric, but
      // retaining it makes the source correspondence explicit.
      const double inside_limit = collapsed_trim_available
          ? selected_trim_distance(1, in_list[1], representatives[in_list[0]])
          : selected_trim_distance(1, first_role, representatives[in_list[1]]);
      const bool keep = outside_value <= outside_limit &&
          inside_value <= inside_limit;
      if (!keep && lacks_free_correlation_ticket(first_role, sequence)) {
        selected_role_pruned[first_role][sequence] = 1;
        finaltrim_lists[first_role][position] = finaltrim_lists[first_role].back();
        finaltrim_lists[first_role].pop_back();
      } else {
        ++position;
      }
    }

    position = 0;
    const std::size_t middle_role = in_list[1];
    while (position < finaltrim_lists[middle_role].size()) {
      const std::uint32_t sequence = finaltrim_lists[middle_role][position];
      const bool keep = selected_trim_distance(0, middle_role, sequence) <=
              selected_trim_distance(
                  0, middle_role, representatives[in_list[0]]) ||
          selected_trim_distance(1, middle_role, sequence) <=
              selected_trim_distance(
                  1, middle_role, representatives[in_list[2]]);
      if (!keep && lacks_free_correlation_ticket(middle_role, sequence)) {
        selected_role_pruned[middle_role][sequence] = 1;
        finaltrim_lists[middle_role][position] = finaltrim_lists[middle_role].back();
        finaltrim_lists[middle_role].pop_back();
      } else {
        ++position;
      }
    }

    // Do not reset position here. The supplied third Do While inherits the
    // final index from INList(1), an observable ordering quirk of FinalTrim.
    const std::size_t last_role = in_list[2];
    while (position < finaltrim_lists[last_role].size()) {
      const std::uint32_t sequence = finaltrim_lists[last_role][position];
      const bool keep = selected_trim_distance(0, last_role, sequence) <=
              selected_trim_distance(
                  0, last_role, representatives[in_list[1]]) ||
          selected_trim_distance(1, last_role, sequence) <=
              selected_trim_distance(
                  1, last_role, representatives[in_list[1]]);
      if (!keep && lacks_free_correlation_ticket(last_role, sequence)) {
        selected_role_pruned[last_role][sequence] = 1;
        finaltrim_lists[last_role][position] = finaltrim_lists[last_role].back();
        finaltrim_lists[last_role].pop_back();
      } else {
        ++position;
      }
    }
  }

  for (std::size_t role = 0; role < 3; ++role) {
    for (const std::uint32_t sequence : finaltrim_lists[role]) {
      if (sequence < ordinary_sequence_count) {
        finaltrim_membership[role][sequence] = 1;
      }
    }
  }

  // OKSeq 15 is now known, so finish the exact ConsensusOK CScore sequence.
  std::array<std::vector<double>, 3> consensus_scores;
  for (std::size_t role = 0; role < 3; ++role) {
    consensus_scores[role].assign(ordinary_sequence_count, 0.0);
    for (auto& evidence : event.role_hypotheses[role].distance_evidence) {
      if (evidence.sequence >= ordinary_sequence_count) continue;
      auto& score = evidence.consensus_score;
      score.final_trim_first_expansion_added =
          first_expansion_added[role][evidence.sequence] != 0;
      score.final_trim_second_expansion_added =
          second_expansion_added[role][evidence.sequence] != 0;
      score.selected_role_pruned_out =
          selected_role_pruned[role][evidence.sequence] != 0;
      score.final_trim_member =
          finaltrim_membership[role][evidence.sequence] != 0;
      if (score.representative_sentinel) {
        score.score_after_final_membership = 1000.0;
        score.score_after_rcorrx = 1000.0;
        score.final_score = 1000.0;
      } else if (!score.other_representative_zero) {
        score.score_after_final_membership =
            score.base_score_before_final_membership *
            (score.final_trim_member ? 2.0 : 1.0);
        score.score_after_rcorrx =
            (score.score_after_final_membership +
             score.maximum_direct_correlation) *
            score.maximum_direct_correlation;
        score.final_score = score.score_after_rcorrx *
            score.source_long_matrix_multiplier;
      }
      score.complete = true;
      consensus_scores[role][evidence.sequence] = score.final_score;
    }
    const std::uint32_t own_representative = reported_representatives[role];
    if (own_representative < ordinary_sequence_count) {
      consensus_scores[role][own_representative] = 1000.0;
    }
  }

  // ConsensusOK first applies its raw-tree topology filter to OKSeq 18 for
  // every ordinary sequence, including representative cells that are absent
  // from a role's candidate-evidence vector.
  std::array<std::vector<std::int8_t>, 3> consensus_match_class;
  for (std::size_t role = 0; role < 3; ++role) {
    consensus_match_class[role].resize(ordinary_sequence_count);
    const std::size_t other_one = kSourceCompRoles[role][0];
    const std::size_t other_two = kSourceCompRoles[role][1];
    for (std::size_t sequence = 0; sequence < ordinary_sequence_count; ++sequence) {
      const auto candidate = static_cast<std::uint32_t>(sequence);
      std::int8_t match_class =
          calc_match_grid[role][sequence].breakpoint_match_class;
      if (match_class > -1) {
        bool topology_consistent =
            (raw_small(0, role, candidate) <= raw_small(0, other_one, candidate) ||
             raw_small(1, role, candidate) <= raw_small(1, other_one, candidate)) &&
            (raw_small(0, role, candidate) <= raw_small(0, other_two, candidate) ||
             raw_small(1, role, candidate) <= raw_small(1, other_two, candidate));
        for (std::size_t side = 0; side < 2 && topology_consistent; ++side) {
          const double anchor_first = raw_small(
              side, role, representatives[other_one]);
          const double anchor_second = raw_small(
              side, role, representatives[other_two]);
          const double candidate_first = raw_small(side, other_one, candidate);
          const double candidate_second = raw_small(side, other_two, candidate);
          if ((anchor_first < anchor_second &&
               candidate_first > candidate_second) ||
              (anchor_first > anchor_second &&
               candidate_first < candidate_second)) {
            topology_consistent = false;
          }
        }
        if (!topology_consistent) match_class = -1;
      }
      consensus_match_class[role][sequence] = match_class;
    }
    for (auto& evidence : event.role_hypotheses[role].distance_evidence) {
      if (evidence.sequence >= ordinary_sequence_count) continue;
      const std::int8_t filtered = consensus_match_class[role][evidence.sequence];
      if (filtered != evidence.calc_match.raw_breakpoint_match_class) {
        evidence.calc_match.topology_filtered = true;
      }
      evidence.calc_match.breakpoint_match_class = filtered;
    }
  }

  const auto consensus_score_gate = [&](std::size_t role, std::size_t sequence) {
    const double regional =
        calc_match_grid[role][sequence].regional_match_score;
    const double score = consensus_scores[role][sequence];
    return (score > 30.0 && regional > 0.05) ||
        score * regional > 2.0;
  };
  std::array<std::vector<std::uint32_t>, 3> consensus_lists;
  std::array<std::vector<std::uint8_t>, 3> primary_membership;
  std::array<std::vector<std::uint8_t>, 3> equivalent_membership;
  std::array<std::vector<std::uint8_t>, 3> straggler_membership;
  for (std::size_t role = 0; role < 3; ++role) {
    primary_membership[role].assign(ordinary_sequence_count, 0);
    equivalent_membership[role].assign(ordinary_sequence_count, 0);
    straggler_membership[role].assign(ordinary_sequence_count, 0);
    const std::size_t other_one = kSourceCompRoles[role][0];
    const std::size_t other_two = kSourceCompRoles[role][1];

    for (std::size_t sequence = 0; sequence < ordinary_sequence_count; ++sequence) {
      const auto candidate = static_cast<std::uint32_t>(sequence);
      if (sequence_disabled(candidate)) continue;
      const std::int8_t match_class = consensus_match_class[role][sequence];
      if (match_class <= -1) continue;
      const double regional =
          calc_match_grid[role][sequence].regional_match_score;
      bool include = false;
      if (consensus_score_gate(role, sequence) && match_class == 1) {
        include = true;
      } else if (regional > 0.1 && match_class == 1) {
        include = true;
      } else {
        bool grouped_with_exact_match = false;
        if ((consensus_score_gate(role, sequence) && match_class > 0) ||
            match_class > 1) {
          for (std::size_t other = 0; other < ordinary_sequence_count; ++other) {
            if (other == sequence || consensus_match_class[role][other] != 1) continue;
            const auto comparison = static_cast<std::uint32_t>(other);
            if (raw_small(1, role, comparison) >=
                    raw_small(1, role, candidate) &&
                raw_small(0, role, comparison) >=
                    raw_small(0, role, candidate)) {
              grouped_with_exact_match = true;
              break;
            }
          }
        }
        include = grouped_with_exact_match;

        if (!include && collapsed_trim_available && match_class > 0) {
          const bool closest_on_either_side =
              (collapsed_small(0, role, candidate) <=
                       collapsed_small(0, other_one, candidate) &&
               collapsed_small(0, role, candidate) <=
                       collapsed_small(0, other_two, candidate)) ||
              (collapsed_small(1, role, candidate) <=
                       collapsed_small(1, other_one, candidate) &&
               collapsed_small(1, role, candidate) <=
                       collapsed_small(1, other_two, candidate));
          const bool within_parent_bounds =
              collapsed_small(0, role, representatives[other_one]) >=
                      collapsed_small(0, other_one, candidate) &&
              collapsed_small(0, role, representatives[other_two]) >=
                      collapsed_small(0, other_two, candidate) &&
              collapsed_small(1, role, representatives[other_one]) >=
                      collapsed_small(1, other_one, candidate) &&
              collapsed_small(1, role, representatives[other_two]) >=
                      collapsed_small(1, other_two, candidate);
          if (closest_on_either_side && within_parent_bounds) {
            const bool zero_with_first_parent =
                direct_small(0, other_one, candidate) == 0.0 &&
                direct_small(1, other_one, candidate) == 0.0;
            const bool zero_with_second_parent =
                direct_small(0, other_two, candidate) == 0.0 &&
                direct_small(1, other_two, candidate) == 0.0;
            include = !zero_with_first_parent && !zero_with_second_parent;
          } else if (direct_small(0, role, candidate) == 0.0 &&
                     direct_small(1, role, candidate) == 0.0) {
            include = true;
          }
        }

        if (!include) {
          const bool exact_raw_position =
              raw_small(0, role, representatives[other_one]) ==
                      raw_small(0, other_one, candidate) &&
              raw_small(0, role, representatives[other_two]) ==
                      raw_small(0, other_two, candidate) &&
              raw_small(1, role, representatives[other_one]) ==
                      raw_small(1, other_one, candidate) &&
              raw_small(1, role, representatives[other_two]) ==
                      raw_small(1, other_two, candidate);
          // The x=0 fallback in the supplied source spells the full direct
          // matrix row as literal sequence zero rather than ISeqs(0).
          const std::uint32_t direct_zero_anchor = role == 0
              ? 0U
              : representatives[role];
          const bool direct_zero =
              source_direct_whole(0, direct_zero_anchor, candidate) == 0.0 &&
              source_direct_whole(1, direct_zero_anchor, candidate) == 0.0;
          include = exact_raw_position || direct_zero;
        }
      }
      if (include) {
        consensus_lists[role].push_back(candidate);
        primary_membership[role][sequence] = 1;
      }
    }

    // ConsensusOK snapshots membership before each of its two widening
    // passes. Entries appended during a pass cannot seed another entry in the
    // same pass.
    const std::vector<std::uint8_t> primary_snapshot = primary_membership[role];
    for (std::size_t sequence = 0; sequence < ordinary_sequence_count; ++sequence) {
      if (primary_snapshot[sequence] != 0 ||
          consensus_match_class[role][sequence] <= 0) {
        continue;
      }
      const auto candidate = static_cast<std::uint32_t>(sequence);
      if (sequence_disabled(candidate)) continue;
      bool include = false;
      if (consensus_score_gate(role, sequence)) {
        for (std::size_t other = 0; other < ordinary_sequence_count; ++other) {
          if (other == sequence || primary_snapshot[other] == 0 ||
              consensus_match_class[role][other] <= 0) {
            continue;
          }
          const auto comparison = static_cast<std::uint32_t>(other);
          bool exact_six_distances = true;
          for (std::size_t reference_role = 0; reference_role < 3;
               ++reference_role) {
            exact_six_distances = exact_six_distances &&
                raw_small(0, reference_role, comparison) ==
                    raw_small(0, reference_role, candidate) &&
                raw_small(1, reference_role, comparison) ==
                    raw_small(1, reference_role, candidate);
          }
          if (exact_six_distances) {
            include = true;
            break;
          }
        }
      } else {
        for (std::size_t other = 0; other < ordinary_sequence_count; ++other) {
          if (other == sequence || primary_snapshot[other] == 0) continue;
          const auto comparison = static_cast<std::uint32_t>(other);
          if (direct_small(1, role, comparison) >
                  source_direct_whole(1, comparison, candidate) &&
              direct_small(0, role, comparison) >
                  source_direct_whole(0, comparison, candidate)) {
            include = true;
            break;
          }
        }
      }
      if (include) {
        consensus_lists[role].push_back(candidate);
        equivalent_membership[role][sequence] = 1;
      }
    }

    std::vector<std::uint8_t> widened_snapshot = primary_snapshot;
    for (std::size_t sequence = 0; sequence < ordinary_sequence_count; ++sequence) {
      if (equivalent_membership[role][sequence] != 0) {
        widened_snapshot[sequence] = 1;
      }
    }
    for (std::size_t sequence = 0; sequence < ordinary_sequence_count; ++sequence) {
      if (widened_snapshot[sequence] != 0 ||
          consensus_match_class[role][sequence] <= -1) {
        continue;
      }
      const auto candidate = static_cast<std::uint32_t>(sequence);
      if (sequence_disabled(candidate)) continue;
      bool include = false;
      for (std::size_t other = 0; other < ordinary_sequence_count; ++other) {
        if (other == sequence || widened_snapshot[other] == 0) continue;
        const auto comparison = static_cast<std::uint32_t>(other);
        if (raw_small(1, role, comparison) >=
                source_raw_whole(1, comparison, candidate) &&
            raw_small(0, role, comparison) >=
                source_raw_whole(0, comparison, candidate) &&
            direct_small(1, role, comparison) >=
                source_direct_whole(1, comparison, candidate) &&
            direct_small(0, role, comparison) >=
                source_direct_whole(0, comparison, candidate)) {
          include = true;
          break;
        }
      }
      if (include) {
        consensus_lists[role].push_back(candidate);
        straggler_membership[role][sequence] = 1;
      }
    }
  }

  const bool consensus_fallback_restored = std::any_of(
      consensus_lists.begin(), consensus_lists.end(),
      [](const auto& list) { return list.empty(); });
  if (consensus_fallback_restored) consensus_lists = finaltrim_lists;

  const auto consensus_rebuilt_lists = consensus_lists;
  std::array<std::vector<std::uint8_t>, 3> selected_tree_cleanup_pruned;
  std::array<std::vector<std::uint8_t>, 3> selected_tree_cleanup_added;
  for (std::size_t role = 0; role < 3; ++role) {
    selected_tree_cleanup_pruned[role].assign(ordinary_sequence_count, 0);
    selected_tree_cleanup_added[role].assign(ordinary_sequence_count, 0);
  }

  // The selected-role portion after the RFF guard is shared by both native
  // call shapes. With ConservativeGroup at its active default zero it forces
  // every rebuilt list to move through the two whole-region trees like its
  // representative, then admits strict four-matrix inliers. This block runs
  // once after role selection whether the earlier list build used RFF=0 in
  // the same call or a prior call.
  std::array<std::vector<float>, 2> movement_distance;
  for (auto& side : movement_distance) side.assign(ordinary_sequence_count, 0.0F);
  for (std::size_t first = 0; first < ordinary_sequence_count; ++first) {
    for (std::size_t second = 0; second < ordinary_sequence_count; ++second) {
      movement_distance[0][first] = static_cast<float>(
          movement_distance[0][first] + static_cast<float>(source_direct_whole(
              0,
              static_cast<std::uint32_t>(first),
              static_cast<std::uint32_t>(second))));
      movement_distance[1][first] = static_cast<float>(
          movement_distance[1][first] + static_cast<float>(source_direct_whole(
              1,
              static_cast<std::uint32_t>(first),
              static_cast<std::uint32_t>(second))));
    }
  }
  std::array<std::array<std::size_t, 3>, 2> representative_rank{};
  for (std::size_t role = 0; role < 3; ++role) {
    std::array<float, 2> representative_movement{};
    for (std::size_t side = 0; side < 2; ++side) {
      for (std::size_t sequence = 0; sequence < ordinary_sequence_count; ++sequence) {
        representative_movement[side] = static_cast<float>(
            representative_movement[side] + static_cast<float>(source_direct_whole(
                side,
                representatives[role],
                static_cast<std::uint32_t>(sequence))));
      }
      for (std::size_t sequence = 0; sequence < ordinary_sequence_count; ++sequence) {
        if (representative_movement[side] > movement_distance[side][sequence]) {
          ++representative_rank[side][role];
        }
      }
    }
  }
  const double source_rank_denominator = ordinary_sequence_count > 1
      ? static_cast<double>(ordinary_sequence_count - 1)
      : 1.0;

  for (std::size_t role = 0; role < 3; ++role) {
    const std::size_t comp_zero = kSourceCompRoles[role][0];
    const std::size_t comp_one = kSourceCompRoles[role][1];
    const std::uint32_t parent_zero = representatives[comp_zero];
    const std::uint32_t parent_one = representatives[comp_one];
    std::vector<std::uint8_t> remove(ordinary_sequence_count, 0);
    bool inlier_on_both_sides = true;

    const bool inside_anchor_outlier =
        raw_small(1, role, parent_zero) > raw_small(1, comp_zero, parent_one) ||
        (direct_small(1, role, parent_zero) >
                 direct_small(1, comp_zero, parent_one) &&
         direct_small(1, role, parent_one) >
                 direct_small(1, comp_zero, parent_one));
    if (inside_anchor_outlier) {
      inlier_on_both_sides = false;
      const double between_parents = raw_small(1, comp_zero, parent_one);
      for (const std::uint32_t sequence : consensus_lists[role]) {
        if (sequence == representatives[role] ||
            direct_small(1, role, sequence) <= 0.0) {
          continue;
        }
        if (raw_small(1, role, parent_zero) !=
                raw_small(1, comp_zero, sequence) ||
            raw_small(1, role, parent_one) !=
                raw_small(1, comp_one, sequence) ||
            raw_small(1, role, sequence) > between_parents ||
            direct_small(1, role, sequence) >
                direct_small(1, comp_zero, sequence) ||
            direct_small(1, role, sequence) >
                direct_small(1, comp_one, sequence)) {
          remove[sequence] = 1;
        }
      }
      if (raw_small(0, role, parent_zero) < raw_small(0, role, parent_one)) {
        for (const std::uint32_t sequence : consensus_lists[role]) {
          if (sequence != representatives[role] &&
              raw_small(0, role, parent_zero) !=
                  raw_small(0, comp_zero, sequence)) {
            remove[sequence] = 1;
          }
        }
      } else {
        for (const std::uint32_t sequence : consensus_lists[role]) {
          if (sequence != representatives[role] &&
              raw_small(0, role, parent_one) !=
                  raw_small(0, comp_one, sequence)) {
            remove[sequence] = 1;
          }
        }
      }
    }

    const bool outside_anchor_outlier =
        raw_small(0, role, parent_zero) > raw_small(0, comp_zero, parent_one) ||
        (direct_small(0, role, parent_zero) >
                 direct_small(0, comp_zero, parent_one) &&
         direct_small(0, role, parent_one) >
                 direct_small(0, comp_zero, parent_one));
    if (outside_anchor_outlier) {
      inlier_on_both_sides = false;
      const double between_parents = raw_small(0, comp_zero, parent_one);
      for (const std::uint32_t sequence : consensus_lists[role]) {
        if (sequence == representatives[role] ||
            direct_small(0, role, sequence) <= 0.0) {
          continue;
        }
        if (raw_small(0, role, parent_zero) !=
                raw_small(0, comp_zero, sequence) ||
            raw_small(0, role, parent_one) !=
                raw_small(0, comp_one, sequence) ||
            raw_small(0, role, sequence) > between_parents ||
            direct_small(0, role, sequence) >
                direct_small(0, comp_zero, sequence) ||
            direct_small(0, role, sequence) >
                direct_small(0, comp_one, sequence)) {
          remove[sequence] = 1;
        }
      }
      if (raw_small(1, role, parent_zero) < raw_small(1, role, parent_one)) {
        for (const std::uint32_t sequence : consensus_lists[role]) {
          if (sequence != representatives[role] &&
              raw_small(1, role, parent_zero) !=
                  raw_small(1, comp_zero, sequence)) {
            remove[sequence] = 1;
          }
        }
      } else {
        for (const std::uint32_t sequence : consensus_lists[role]) {
          if (sequence != representatives[role] &&
              raw_small(1, role, parent_one) !=
                  raw_small(1, comp_one, sequence)) {
            remove[sequence] = 1;
          }
        }
      }
    }

    if (inlier_on_both_sides) {
      const bool parent_zero_closer_inside =
          raw_small(1, role, parent_zero) < raw_small(1, role, parent_one);
      const std::size_t inside_outlier =
          parent_zero_closer_inside ? comp_one : comp_zero;
      const std::size_t inside_inlier =
          parent_zero_closer_inside ? comp_zero : comp_one;
      const std::uint32_t inside_inlier_sequence = representatives[inside_inlier];
      for (const std::uint32_t sequence : consensus_lists[role]) {
        if (sequence == representatives[role]) continue;
        if (raw_small(1, role, sequence) >=
                raw_small(1, inside_outlier, sequence) ||
            direct_small(1, role, sequence) >
                direct_small(1, role, inside_inlier_sequence) * 6.0) {
          remove[sequence] = 1;
        }
        if (raw_small(1, role, sequence) >
                raw_small(1, role, inside_inlier_sequence) &&
            direct_small(1, role, inside_inlier_sequence) <
                direct_small(0, role, inside_inlier_sequence) &&
            direct_small(1, inside_inlier, sequence) >
                direct_small(0, inside_inlier, sequence)) {
          remove[sequence] = 1;
        }
      }
      const double inside_rank = static_cast<double>(
          representative_rank[1][inside_outlier]);
      const double outside_rank = static_cast<double>(
          representative_rank[0][inside_outlier]);
      if ((inside_rank / source_rank_denominator > 0.95 &&
           outside_rank / source_rank_denominator < 0.75) ||
          (inside_rank - outside_rank) / source_rank_denominator > 0.5) {
        for (const std::uint32_t sequence : consensus_lists[role]) {
          if (sequence != representatives[role] &&
              raw_small(1, inside_inlier, sequence) >
                  raw_small(1, role, inside_inlier_sequence)) {
            remove[sequence] = 1;
          }
        }
      }

      const bool parent_zero_closer_outside =
          raw_small(0, role, parent_zero) < raw_small(0, role, parent_one);
      const std::size_t outside_outlier =
          parent_zero_closer_outside ? comp_one : comp_zero;
      const std::size_t outside_inlier =
          parent_zero_closer_outside ? comp_zero : comp_one;
      const std::uint32_t outside_inlier_sequence = representatives[outside_inlier];
      for (const std::uint32_t sequence : consensus_lists[role]) {
        if (sequence == representatives[role]) continue;
        if (raw_small(0, role, sequence) >=
                raw_small(0, outside_outlier, sequence) ||
            direct_small(0, role, sequence) >
                direct_small(0, role, outside_inlier_sequence) * 6.0) {
          remove[sequence] = 1;
        }
        if (raw_small(0, role, sequence) >
                raw_small(0, role, outside_inlier_sequence) &&
            direct_small(0, role, outside_inlier_sequence) <
                direct_small(1, role, outside_inlier_sequence) &&
            direct_small(0, outside_inlier, sequence) >
                direct_small(1, outside_inlier, sequence)) {
          remove[sequence] = 1;
        }
      }
      const double outside_rank = static_cast<double>(
          representative_rank[0][outside_outlier]);
      const double inside_rank = static_cast<double>(
          representative_rank[1][outside_outlier]);
      if ((outside_rank / source_rank_denominator > 0.95 &&
           inside_rank / source_rank_denominator < 0.75) ||
          (outside_rank - inside_rank) / source_rank_denominator > 0.5) {
        for (const std::uint32_t sequence : consensus_lists[role]) {
          if (sequence != representatives[role] &&
              raw_small(0, outside_inlier, sequence) >
                  raw_small(0, role, outside_inlier_sequence)) {
            remove[sequence] = 1;
          }
        }
      }
    }

    std::size_t position = 0;
    while (position < consensus_lists[role].size()) {
      const std::uint32_t sequence = consensus_lists[role][position];
      if (remove[sequence] != 0) {
        selected_tree_cleanup_pruned[role][sequence] = 1;
        consensus_lists[role][position] = consensus_lists[role].back();
        consensus_lists[role].pop_back();
      } else {
        ++position;
      }
    }

    std::array<double, 4> highest_distance{};
    for (const std::uint32_t sequence : consensus_lists[role]) {
      highest_distance[0] = std::max(
          highest_distance[0], raw_small(1, role, sequence));
      highest_distance[1] = std::max(
          highest_distance[1], raw_small(0, role, sequence));
      highest_distance[2] = std::max(
          highest_distance[2], direct_small(1, role, sequence));
      highest_distance[3] = std::max(
          highest_distance[3], direct_small(0, role, sequence));
    }
    std::vector<std::uint8_t> present(ordinary_sequence_count, 0);
    for (const std::uint32_t sequence : consensus_lists[role]) {
      if (sequence < ordinary_sequence_count) present[sequence] = 1;
    }
    for (std::size_t sequence = 0; sequence < ordinary_sequence_count; ++sequence) {
      const auto candidate = static_cast<std::uint32_t>(sequence);
      if (sequence_disabled(candidate)) continue;
      bool admitted = false;
      // Preserve line 24895's fourth comparison: FAMatSmall is compared with
      // the maximum direct outside distance HDF rather than FMatSmall.
      if (raw_small(1, role, candidate) <= highest_distance[0] &&
          raw_small(0, role, candidate) <= highest_distance[1] &&
          direct_small(1, role, candidate) <= highest_distance[2] &&
          raw_small(0, role, candidate) <= highest_distance[3] &&
          raw_small(1, role, candidate) < raw_small(1, comp_zero, candidate) &&
          raw_small(1, role, candidate) < raw_small(1, comp_one, candidate) &&
          raw_small(0, role, candidate) < raw_small(0, comp_zero, candidate) &&
          raw_small(0, role, candidate) < raw_small(0, comp_one, candidate)) {
        admitted = true;
      }
      // `ILP(0) = 0 And ILP(1) = 0 Or x = x` is always true in the supplied
      // source, so this strict direct+raw four-matrix gate is unconditional
      // whenever the bounded admission above fails.
      if (!admitted &&
          raw_small(1, role, candidate) < raw_small(1, comp_zero, candidate) &&
          raw_small(1, role, candidate) < raw_small(1, comp_one, candidate) &&
          direct_small(1, role, candidate) < direct_small(1, comp_zero, candidate) &&
          direct_small(1, role, candidate) < direct_small(1, comp_one, candidate) &&
          raw_small(0, role, candidate) < raw_small(0, comp_zero, candidate) &&
          raw_small(0, role, candidate) < raw_small(0, comp_one, candidate) &&
          direct_small(0, role, candidate) < direct_small(0, comp_zero, candidate) &&
          direct_small(0, role, candidate) < direct_small(0, comp_one, candidate)) {
        admitted = true;
      }
      if (admitted && present[sequence] == 0) {
        consensus_lists[role].push_back(candidate);
        present[sequence] = 1;
        selected_tree_cleanup_added[role][sequence] = 1;
      }
    }
  }

  for (auto& list : consensus_lists) {
    list.erase(
        std::remove_if(
            list.begin(),
            list.end(),
            [&](std::uint32_t sequence) { return sequence_disabled(sequence); }),
        list.end());
  }

  // The next active block in native FinalTrim rechecks every surviving RList
  // candidate against the other two role representatives. For an RDP-origin
  // event ProgramFlag=0 dispatches XOver here; the other method families are
  // deliberately not inferred. Native temporarily widens LowestProb to at
  // least LowP*100000 (and its corrected project threshold). Generate the
  // structurally valid browser RDP candidates without the ordinary scan gate,
  // then retain only rows inside that explicit widened local threshold. The
  // normal corrected-cutoff result remains separate and auditable.
  std::array<std::vector<PostGroupRdpRecheckEvidence>, 3>
      post_group_rdp_rechecks;
  std::array<std::vector<MaxChiRecheckEvidence>, 3>
      post_group_maxchi_rechecks;
  const double post_group_local_cutoff = std::min(
      1.0,
      std::max(
          options_.p_value_cutoff,
          event.best_local_p_value * 100000.0));
  for (std::size_t role = 0; role < 3; ++role) {
    post_group_rdp_rechecks[role].resize(ordinary_sequence_count);
    post_group_maxchi_rechecks[role].resize(ordinary_sequence_count);
    const std::size_t comp_zero = kSourceCompRoles[role][0];
    const std::size_t comp_one = kSourceCompRoles[role][1];
    for (const std::uint32_t sequence : consensus_lists[role]) {
      if (sequence >= ordinary_sequence_count) continue;
      const std::array<std::uint32_t, 3> triplet{
          sequence,
          representatives[comp_zero],
          representatives[comp_one],
      };
      auto& recheck = post_group_rdp_rechecks[role][sequence];
      recheck.local_p_value_cutoff = post_group_local_cutoff;
      if (sequence == representatives[role]) {
        // Module2.bas line 25099 skips ISeqs(WinPP) itself.
        recheck.representative_skipped = true;
        post_group_maxchi_rechecks[role][sequence].representative_skipped = true;
        continue;
      }
      post_group_maxchi_rechecks[role][sequence] =
          maxchi_triplet_recheck(triplet);
      recheck.requested = true;
      bool profile_available = false;
      const auto candidates = triplet_signals(
          triplet, false, &profile_available, &profile_scratch_);
      recheck.profile_available = profile_available;
      for (const auto& candidate : candidates) {
        if (candidate.local_p_value > post_group_local_cutoff) continue;
        ++recheck.emitted_signal_count;
        if (candidate.recombinant != sequence) continue;
        ++recheck.candidate_signal_count;
        const double overlap = tract_overlap(
            event.beginning,
            event.ending,
            candidate.beginning,
            candidate.ending);
        if (overlap <= 0.3) continue;
        ++recheck.overlapping_signal_count;
        if (!recheck.event_redetected ||
            candidate.local_p_value < recheck.best_local_p_value ||
            (candidate.local_p_value == recheck.best_local_p_value &&
             overlap > recheck.best_overlap)) {
          recheck.event_redetected = true;
          recheck.best_beginning = candidate.beginning;
          recheck.best_ending = candidate.ending;
          recheck.best_wraps_origin = candidate.wraps_origin;
          recheck.best_overlap = overlap;
          recheck.best_local_p_value = candidate.local_p_value;
          recheck.best_corrected_p_value = candidate.corrected_p_value;
        }
      }
      recheck.significant = recheck.event_redetected &&
          recheck.best_corrected_p_value < options_.p_value_cutoff;
    }
  }

  for (std::size_t role = 0; role < 3; ++role) {
    std::vector<std::uint8_t> rebuilt_membership(ordinary_sequence_count, 0);
    for (const std::uint32_t sequence : consensus_rebuilt_lists[role]) {
      if (sequence < ordinary_sequence_count) rebuilt_membership[sequence] = 1;
    }
    std::vector<std::uint8_t> final_distance_membership(
        ordinary_sequence_count, 0);
    for (const std::uint32_t sequence : consensus_lists[role]) {
      if (sequence < ordinary_sequence_count) {
        final_distance_membership[sequence] = 1;
      }
    }
    auto& hypothesis = event.role_hypotheses[role];
    hypothesis.distance_correlation_set = consensus_lists[role];
    sort_unique(hypothesis.distance_correlation_set);
    for (auto& evidence : hypothesis.distance_evidence) {
      if (evidence.sequence >= ordinary_sequence_count) continue;
      auto& score = evidence.consensus_score;
      score.consensus_primary_member =
          primary_membership[role][evidence.sequence] != 0;
      score.consensus_equivalent_member =
          equivalent_membership[role][evidence.sequence] != 0;
      score.consensus_straggler_member =
          straggler_membership[role][evidence.sequence] != 0;
      score.consensus_rebuilt_member =
          rebuilt_membership[evidence.sequence] != 0;
      score.consensus_fallback_restored = consensus_fallback_restored;
      score.selected_tree_cleanup_pruned_out =
          selected_tree_cleanup_pruned[role][evidence.sequence] != 0;
      score.selected_tree_cleanup_added =
          selected_tree_cleanup_added[role][evidence.sequence] != 0;
      score.final_distance_member =
          final_distance_membership[evidence.sequence] != 0;
      evidence.post_group_rdp_recheck =
          post_group_rdp_rechecks[role][evidence.sequence];
      evidence.post_group_maxchi_recheck =
          post_group_maxchi_rechecks[role][evidence.sequence];
      evidence.significant = score.final_distance_member;
    }
    std::stable_sort(
        hypothesis.distance_evidence.begin(),
        hypothesis.distance_evidence.end(),
        [](const DistanceCorrelationEvidence& left,
           const DistanceCorrelationEvidence& right) {
          if (left.significant != right.significant) {
            return left.significant > right.significant;
          }
          if (left.consensus_score.final_score !=
              right.consensus_score.final_score) {
            return left.consensus_score.final_score >
                right.consensus_score.final_score;
          }
          return left.sequence < right.sequence;
        });

    hypothesis.complete_two_of_three_set.clear();
    for (std::size_t sequence = 0; sequence < ordinary_sequence_count; ++sequence) {
      const auto candidate = static_cast<std::uint32_t>(sequence);
      if (sequence_disabled(candidate)) continue;
      const std::size_t evidence_count =
          static_cast<std::size_t>(std::binary_search(
              hypothesis.detectable_signal_set.begin(),
              hypothesis.detectable_signal_set.end(), candidate)) +
          static_cast<std::size_t>(std::binary_search(
              hypothesis.distance_correlation_set.begin(),
              hypothesis.distance_correlation_set.end(), candidate)) +
          static_cast<std::size_t>(std::binary_search(
              hypothesis.phylogenetic_correlation_set.begin(),
              hypothesis.phylogenetic_correlation_set.end(), candidate));
      if (candidate == hypothesis.presumed_recombinant || evidence_count >= 2) {
        hypothesis.complete_two_of_three_set.push_back(candidate);
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

MaxChiRecheckEvidence RdpScanner::maxchi_triplet_recheck(
    const std::array<std::uint32_t, 3>& triplet) {
  MaxChiRecheckEvidence evidence;
  evidence.requested = true;
  if (working_alignment_.length == 0 ||
      std::any_of(triplet.begin(), triplet.end(), [&](std::uint32_t sequence) {
        return sequence >= working_alignment_.sequence_count() ||
            sequence >= working_origins_.size();
      })) {
    return evidence;
  }

  auto& triplet_missing_data = maxchi_workspace_.triplet_missing_data;
  triplet_missing_data.assign(working_alignment_.length, 0);
  for (std::size_t coordinate = 0; coordinate < working_alignment_.length; ++coordinate) {
    for (const std::uint32_t working_sequence : triplet) {
      const std::uint32_t origin = working_origins_[working_sequence];
      if (origin >= alignment_.sequence_count()) continue;
      const std::size_t original_offset = origin * alignment_.length + coordinate;
      const bool input_missing = original_offset < native_input_missing_data_.size() &&
          native_input_missing_data_[original_offset] != 0;
      const bool erased_or_fragment_gap = alignment_.at(origin, coordinate) != 0 &&
          working_alignment_.at(working_sequence, coordinate) == 0;
      if (input_missing || erased_or_fragment_gap) {
        triplet_missing_data[coordinate] = 1;
        break;
      }
    }
  }

  MaxChiRecheckOptions recheck_options;
  recheck_options.circular = options_.circular;
  recheck_options.bonferroni =
      options_.correction == CorrectionMode::bonferroni;
  recheck_options.p_value_cutoff = options_.p_value_cutoff;
  recheck_options.correction_tests =
      std::max<std::uint64_t>(1, correction_tests_);
  return maxchi_recheck(
      working_alignment_,
      triplet,
      triplet_missing_data,
      recheck_options,
      maxchi_workspace_);
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
  event.wraps_origin = anchor.wraps_origin ||
      (options_.circular && event.beginning >= event.ending);
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

  const std::size_t pre_polish_beginning = event.beginning;
  const std::size_t pre_polish_ending = event.ending;
  refresh_breakpoint_confidence(event, true);
  if (event.beginning != pre_polish_beginning || event.ending != pre_polish_ending) {
    // PolishBP runs after representative refinement. Its revised boundaries
    // become the stored event coordinates, so rebuild boundary-dependent
    // evidence once before cyclic erasure and ordered review.
    build_event_evidence(event);
    refresh_trace_evidence(event);
    refresh_role_hypotheses(event);
  }
  refresh_breakpoint_context(event);

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
  const bool roles_changed =
      events_[event_id].recombinant != recombinant ||
      events_[event_id].major_parent != major_parent ||
      events_[event_id].minor_parent != minor_parent;
  const bool breakpoints_changed =
      events_[event_id].beginning != beginning ||
      events_[event_id].ending != ending;
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
  if (sequence_disabled(recombinant) || sequence_disabled(major_parent) ||
      sequence_disabled(minor_parent)) {
    error = "Disabled sequences cannot be assigned an analysed event role.";
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
  event.wraps_origin = options_.circular && beginning >= ending;
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
  const std::size_t pre_polish_beginning = event.beginning;
  const std::size_t pre_polish_ending = event.ending;
  // The supplied role-replacement path reruns PolishBP. Preserve an explicit
  // breakpoint edit, but apply the native result when only representatives
  // changed (including the UI's weighted role recommendation).
  refresh_breakpoint_confidence(
      event,
      roles_changed && !breakpoints_changed);
  if (event.beginning != pre_polish_beginning || event.ending != pre_polish_ending) {
    build_event_evidence(event);
    refresh_trace_evidence(event);
    refresh_role_hypotheses(event);
  }
  refresh_breakpoint_context(event);
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
          return sequence >= alignment_.sequence_count() || sequence_disabled(sequence);
        })) {
      error = "The edited co-recombinant group contains an unavailable or disabled sequence.";
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
  if (options.mask.size() != alignment_.sequence_count() ||
      options.disabled.size() != alignment_.sequence_count()) {
    error = "The saved sequence curation state does not match the saved alignment.";
    return false;
  }
  std::size_t primary_sequence_count = 0;
  for (std::size_t sequence = 0; sequence < alignment_.sequence_count(); ++sequence) {
    if (options.disabled[sequence] != 0) options.mask[sequence] = 0;
    if (options.mask[sequence] == 0 && options.disabled[sequence] == 0) {
      ++primary_sequence_count;
    }
  }
  if (!(options.p_value_cutoff > 0.0 && options.p_value_cutoff <= 1.0) ||
      options.window_sites < 5 ||
      primary_sequence_count < 3) {
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
  if (sequence_disabled(recombinant) || sequence_disabled(major_parent) ||
      sequence_disabled(minor_parent)) {
    error = "A saved event assigns an analysed role to a disabled sequence.";
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
  event.wraps_origin = options_.circular && beginning >= ending;
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
                  sequence == minor_parent || sequence_disabled(sequence);
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
  // Saved event coordinates are authoritative. Recompute statistical
  // evidence for the restored triplet without silently rewriting the user's
  // checkpoint during import.
  refresh_breakpoint_confidence(event, false);
  refresh_breakpoint_context(event);
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
  out << "{\"engineVersion\":\"0.9.0-session-9\",\"status\":\"cyclic-three-set-reconciled\","
         "\"method\":\"RDP\",\"reconciliationTier\":\"detectable-distance-phylogenetic\","
         "\"cycleMode\":\"strongest-first-tract-erasure-with-bounded-fragment-reentry\","
         "\"lateConsensus\":{\"status\":\"active-rdp-maxchi-post-group-recheck\","
         "\"groupPruningApplied\":true,\"nativeGroupMembershipComplete\":true,"
         "\"primaryRdpPostGroupRecheckApplied\":true,"
         "\"nativePrimaryRdpRecheckComplete\":true,"
         "\"maxChiTripletRecheckApplied\":true,"
         "\"maxChiPostGroupRecheckApplied\":true,"
         "\"maxChiKernelStatus\":\"source-shaped-strongest-peak-unvalidated\","
         "\"maxChiEventDiscoveryApplied\":false,"
         "\"nativeMaxChiFullRecheckComplete\":false,"
         "\"implementedStages\":["
         "\"MakeINList/MakeACOR\",\"MakeRList\",\"StripDupInv\","
         "\"FinalTrim duplicate-correlation cleanup\","
         "\"FinalTrim OKSeq 6 nearest-nonrecombinant fixed-point membership\","
         "\"FinalTrim active OKSeq 7-14 matrix scoring (10/11 source-zero)\","
         "\"FinalTrim first and second ascending final-list expansions\","
         "\"FinalTrim selected-role pruning with source list-order semantics\","
         "\"FinalTrim OKSeq 15 final membership\","
         "\"CalcMatchY active OKSeq 17/18 grouping inputs\","
         "\"ConsensusOK OKSeq 18 topology-consistency filter\","
         "\"ConsensusOK complete CScore with OKSeq 15 RCorrX and Long NS semantics\","
         "\"ConsensusOK primary equivalence and straggler list rebuild passes\","
         "\"ConsensusOK empty-role fallback\","
         "\"Selected-role conservative raw/direct tree cleanup and strict inlier admission\","
         "\"Active two-of-three group update\","
         "\"FinalTrim primary-RDP post-group signal and probability recheck\","
         "\"FastRecCheckMC2 strongest-peak MaxChi representative and finalized-list recheck\"],"
         "\"pendingStages\":["
         "\"MaxChi exploratory event discovery and native multi-peak destroy/retry scheduler\","
         "\"remaining non-RDP method signal rechecks\"]},"
         "\"breakpointInspection\":{\"available\":true,"
         "\"source\":\"original-alignment\",\"maxFlankSites\":100,\"maxRows\":64,"
         "\"nativeCheckEndsStatus\":\"complete-active-unvalidated\","
         "\"nativeCheckEndsAfterFirstEvent\":true,\"inputMissingRunLength\":10,"
         "\"uncertaintyReasons\":[\"prior-erasure\",\"input-missing-data\","
         "\"linear-edge\",\"profile-unavailable\"],"
         "\"statisticalConfidence\":{\"available\":true,"
         "\"status\":\"complete-active-unvalidated\","
         "\"method\":\"BURT/BenHMM\",\"hmmCyclesArgument\":20,"
         "\"serialTrainingStarts\":21,"
         "\"posteriorThresholds\":[0.995,0.999],"
         "\"randomSeed\":3,\"randomAdapter\":\"msvc-rand-15-bit\","
         "\"optionAvailable\":true,\"enabledByDefault\":true,"
         "\"canRepositionDetectedEvents\":true}},"
         "\"treeInspection\":{\"available\":true,"
         "\"source\":\"reconciliation-tree-panel\",\"regionCount\":6,"
         "\"bootstrapCollapseCutoff\":0.5,\"payload\":\"on-demand-edge-lists\"},"
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
      << ",\"polishBreakpoints\":"
      << (options_.polish_breakpoints ? "true" : "false")
      << ",\"maskedSequenceIndices\":[";
  bool wrote_mask = false;
  for (std::size_t index = 0; index < options_.mask.size(); ++index) {
    if (options_.mask[index] == 0) continue;
    if (wrote_mask) out << ',';
    wrote_mask = true;
    out << index;
  }
  out << "],\"disabledSequenceIndices\":[";
  bool wrote_disabled = false;
  for (std::size_t index = 0; index < options_.disabled.size(); ++index) {
    if (options_.disabled[index] == 0) continue;
    if (wrote_disabled) out << ',';
    wrote_disabled = true;
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
        << ",\"breakpointConfidence\":";
    write_breakpoint_confidence_json(out, event.breakpoint_confidence);
    out << ",\"breakpointContext\":{\"source\":\"cyclic-erasure-history\","
           "\"beginning\":{\"erasureAdjacent\":"
        << (event.adjacent_erasure_event_ids[0].empty() ? "false" : "true")
        << ",\"erasureWithinRdpWindow\":"
        << (event.uncertain_erasure_event_ids[0].empty() ? "false" : "true")
        << ",\"uncertainDueToErasure\":"
        << (event.uncertain_erasure_event_ids[0].empty() ? "false" : "true")
        << ",\"nativeCheckEndsApplied\":"
        << (event.breakpoint_uncertainty[0].native_check_ends_applied
                ? "true"
                : "false")
        << ",\"nativeCheckEndsWarning\":"
        << (event.breakpoint_uncertainty[0].native_check_ends_warning
                ? "true"
                : "false")
        << ",\"informationProfileAvailable\":"
        << (event.breakpoint_uncertainty[0].information_profile_available
                ? "true"
                : "false")
        << ",\"inputMissingDataInCheckRange\":"
        << (event.breakpoint_uncertainty[0].input_missing_data_in_range
                ? "true"
                : "false")
        << ",\"linearEdgeWithinRdpWindow\":"
        << (event.breakpoint_uncertainty[0].linear_edge_within_window
                ? "true"
                : "false")
        << ",\"nativeCheckRange\":{\"beginning\":"
        << event.breakpoint_uncertainty[0].check_range_beginning
        << ",\"ending\":"
        << event.breakpoint_uncertainty[0].check_range_ending
        << ",\"wrapsOrigin\":"
        << (event.breakpoint_uncertainty[0].check_range_wraps_origin
                ? "true"
                : "false")
        << ",\"coordinateCount\":"
        << event.breakpoint_uncertainty[0].check_coordinate_count << '}'
        << ",\"rdpWindowInformativeSites\":" << options_.window_sites
        << ",\"nearestErasureInformativeSites\":";
    if (event.nearest_erasure_informative_sites[0] < 0) {
      out << "null";
    } else {
      out << event.nearest_erasure_informative_sites[0];
    }
    out << ",\"uncertainPriorEventIds\":[";
    for (std::size_t prior = 0;
         prior < event.uncertain_erasure_event_ids[0].size();
         ++prior) {
      if (prior) out << ',';
      out << event.uncertain_erasure_event_ids[0][prior];
    }
    out << "],\"priorEventIds\":[";
    for (std::size_t prior = 0;
         prior < event.adjacent_erasure_event_ids[0].size();
         ++prior) {
      if (prior) out << ',';
      out << event.adjacent_erasure_event_ids[0][prior];
    }
    out << "]},\"ending\":{\"erasureAdjacent\":"
        << (event.adjacent_erasure_event_ids[1].empty() ? "false" : "true")
        << ",\"erasureWithinRdpWindow\":"
        << (event.uncertain_erasure_event_ids[1].empty() ? "false" : "true")
        << ",\"uncertainDueToErasure\":"
        << (event.uncertain_erasure_event_ids[1].empty() ? "false" : "true")
        << ",\"nativeCheckEndsApplied\":"
        << (event.breakpoint_uncertainty[1].native_check_ends_applied
                ? "true"
                : "false")
        << ",\"nativeCheckEndsWarning\":"
        << (event.breakpoint_uncertainty[1].native_check_ends_warning
                ? "true"
                : "false")
        << ",\"informationProfileAvailable\":"
        << (event.breakpoint_uncertainty[1].information_profile_available
                ? "true"
                : "false")
        << ",\"inputMissingDataInCheckRange\":"
        << (event.breakpoint_uncertainty[1].input_missing_data_in_range
                ? "true"
                : "false")
        << ",\"linearEdgeWithinRdpWindow\":"
        << (event.breakpoint_uncertainty[1].linear_edge_within_window
                ? "true"
                : "false")
        << ",\"nativeCheckRange\":{\"beginning\":"
        << event.breakpoint_uncertainty[1].check_range_beginning
        << ",\"ending\":"
        << event.breakpoint_uncertainty[1].check_range_ending
        << ",\"wrapsOrigin\":"
        << (event.breakpoint_uncertainty[1].check_range_wraps_origin
                ? "true"
                : "false")
        << ",\"coordinateCount\":"
        << event.breakpoint_uncertainty[1].check_coordinate_count << '}'
        << ",\"rdpWindowInformativeSites\":" << options_.window_sites
        << ",\"nearestErasureInformativeSites\":";
    if (event.nearest_erasure_informative_sites[1] < 0) {
      out << "null";
    } else {
      out << event.nearest_erasure_informative_sites[1];
    }
    out << ",\"uncertainPriorEventIds\":[";
    for (std::size_t prior = 0;
         prior < event.uncertain_erasure_event_ids[1].size();
         ++prior) {
      if (prior) out << ',';
      out << event.uncertain_erasure_event_ids[1][prior];
    }
    out << "],\"priorEventIds\":[";
    for (std::size_t prior = 0;
         prior < event.adjacent_erasure_event_ids[1].size();
         ++prior) {
      if (prior) out << ',';
      out << event.adjacent_erasure_event_ids[1][prior];
    }
    out << "]}}"
        << ",\"maxChiTripletRecheck\":";
    write_maxchi_recheck_json(out, event.maxchi_triplet_recheck);
    out << ",\"bestLocalPValue\":";
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
            << (evidence.duplicate_cleaned_support ? "true" : "false")
            << ",\"finalTrimMatrix\":{\"status\":\"complete-active-rff0\","
               "\"implementedOKSeqElements\":[7,8,9,10,11,12,13,14],"
               "\"inactiveZeroOKSeqElements\":[10,11],"
               "\"collapsedTreePositionScore\":";
        json::number(out, evidence.final_trim_matrix.collapsed_tree_position_score);
        out << ",\"rawTreePositionScore\":";
        json::number(out, evidence.final_trim_matrix.raw_tree_position_score);
        out << ",\"relativeDistanceScore\":";
        json::number(out, evidence.final_trim_matrix.relative_distance_score);
        out << ",\"breakpointDistanceScores\":[";
        json::number(out, evidence.final_trim_matrix.breakpoint_distance_scores[0]);
        out << ',';
        json::number(out, evidence.final_trim_matrix.breakpoint_distance_scores[1]);
        out << "],\"breakpointScoreAvailable\":["
            << (evidence.final_trim_matrix.breakpoint_score_available[0]
                    ? "true"
                    : "false")
            << ','
            << (evidence.final_trim_matrix.breakpoint_score_available[1]
                    ? "true"
                    : "false")
            << "],\"detectedRegionDistanceScore\":";
        json::number(out, evidence.final_trim_matrix.detected_region_distance_score);
        out << ",\"detectedRegionMatchDistances\":[";
        for (std::size_t distance = 0;
             distance < evidence.final_trim_matrix.detected_region_match_distances.size();
             ++distance) {
          if (distance) out << ',';
          json::number(
              out,
              evidence.final_trim_matrix.detected_region_match_distances[distance]);
        }
        out << "],\"detectedEventMatch\":"
            << (evidence.final_trim_matrix.detected_event_match ? "true" : "false")
            << ",\"detectedEventOverlap\":";
        json::number(out, evidence.final_trim_matrix.detected_event_overlap);
        out << ",\"detectedEventBeginning\":";
        if (evidence.final_trim_matrix.detected_event_match) {
          out << evidence.final_trim_matrix.detected_event_beginning;
        } else {
          out << "null";
        }
        out << ",\"detectedEventEnding\":";
        if (evidence.final_trim_matrix.detected_event_match) {
          out << evidence.final_trim_matrix.detected_event_ending;
        } else {
          out << "null";
        }
        out << ",\"detectedEventSignalId\":";
        if (evidence.final_trim_matrix.detected_event_signal_id >= 0) {
          out << evidence.final_trim_matrix.detected_event_signal_id;
        } else {
          out << "null";
        }
        out << ",\"detectedRegionSaturated\":"
            << (evidence.final_trim_matrix.detected_region_saturated ? "true" : "false")
            << ",\"sourceSequenceIndexQuirkApplied\":"
            << (evidence.final_trim_matrix.source_sequence_index_quirk_applied
                    ? "true"
                    : "false")
            << ",\"activeConsensusMatrixScore\":";
        json::number(out, evidence.final_trim_matrix.active_consensus_matrix_score);
        out << ",\"treeDistanceFallback\":"
            << (evidence.final_trim_matrix.tree_distance_fallback ? "true" : "false")
            << ",\"appliesToNonrepresentative\":"
            << (evidence.final_trim_matrix.applies_to_nonrepresentative
                    ? "true"
                    : "false")
            << "},\"calcMatch\":{\"status\":\""
            << (evidence.calc_match.available ? "complete-active-rff0" : "unavailable")
            << "\",\"implementedOKSeqElements\":[17,18],"
               "\"standardGroupingThresholds\":true,\"fragmentVariableSites\":"
            << evidence.calc_match.fragment_variable_sites
            << ",\"targetHalfWindow\":" << evidence.calc_match.target_half_window
            << ",\"smoothingHalfWindow\":"
            << evidence.calc_match.smoothing_half_window
            << ",\"regionalMatchScore\":";
        json::number(out, evidence.calc_match.regional_match_score);
        out << ",\"rawBreakpointMatchClass\":"
            << static_cast<int>(evidence.calc_match.raw_breakpoint_match_class)
            << ",\"breakpointMatchClass\":"
            << static_cast<int>(evidence.calc_match.breakpoint_match_class)
            << ",\"checkpointMatches\":[";
        for (std::size_t checkpoint = 0;
             checkpoint < evidence.calc_match.checkpoint_matches.size();
             ++checkpoint) {
          if (checkpoint) out << ',';
          json::number(out, evidence.calc_match.checkpoint_matches[checkpoint]);
        }
        out << "],\"breakpointsExist\":["
            << (evidence.calc_match.breakpoints_exist[0] ? "true" : "false") << ','
            << (evidence.calc_match.breakpoints_exist[1] ? "true" : "false")
            << "],\"topologyFiltered\":"
            << (evidence.calc_match.topology_filtered ? "true" : "false")
            << ",\"topologyDistanceFallback\":"
            << (evidence.calc_match.topology_distance_fallback ? "true" : "false");
        if (!evidence.calc_match.available) {
          out << ",\"unavailableReason\":"
                 "\"insufficient-variable-sites-or-source-fragment-bound\"";
        }
        out << "},\"consensusScore\":{"
               "\"status\":\"complete-active-rff0\","
               "\"implementedOKSeqElements\":[0,1,2,3,4,5,6,15],"
               "\"correctedCorrelationPValue\":";
        json::number(out, evidence.consensus_score.corrected_correlation_p_value);
        out << ",\"detectedEventOverlap\":";
        json::number(out, evidence.consensus_score.detected_event_overlap);
        out << ",\"patternScore\":";
        json::number(out, evidence.consensus_score.pattern_score);
        out << ",\"detectableSetMember\":"
            << (evidence.consensus_score.detectable_set_member ? "true" : "false")
            << ",\"initialRListMember\":"
            << (evidence.consensus_score.initial_rlist_member ? "true" : "false")
            << ",\"duplicateCleanedMember\":"
            << (evidence.consensus_score.duplicate_cleaned_member ? "true" : "false")
            << ",\"nearestNonrecombinantMember\":"
            << (evidence.consensus_score.nearest_nonrecombinant_member
                    ? "true"
                    : "false")
            << ",\"finalTrimFirstExpansionAdded\":"
            << (evidence.consensus_score.final_trim_first_expansion_added
                    ? "true"
                    : "false")
            << ",\"finalTrimSecondExpansionAdded\":"
            << (evidence.consensus_score.final_trim_second_expansion_added
                    ? "true"
                    : "false")
            << ",\"selectedRolePrunedOut\":"
            << (evidence.consensus_score.selected_role_pruned_out
                    ? "true"
                    : "false")
            << ",\"finalTrimMember\":"
            << (evidence.consensus_score.final_trim_member ? "true" : "false")
            << ",\"maximumDirectCorrelation\":";
        json::number(out, evidence.consensus_score.maximum_direct_correlation);
        out << ",\"baseScoreBeforeFinalMembership\":";
        json::number(
            out,
            evidence.consensus_score.base_score_before_final_membership);
        out << ",\"scoreAfterFinalMembership\":";
        json::number(out, evidence.consensus_score.score_after_final_membership);
        out << ",\"scoreAfterRcorrx\":";
        json::number(out, evidence.consensus_score.score_after_rcorrx);
        out << ",\"sourceLongMatrixMultiplier\":";
        json::number(out, evidence.consensus_score.source_long_matrix_multiplier);
        out << ",\"sourceLongNsSemantics\":true,\"finalScore\":";
        json::number(out, evidence.consensus_score.final_score);
        out << ",\"consensusPrimaryMember\":"
            << (evidence.consensus_score.consensus_primary_member
                    ? "true"
                    : "false")
            << ",\"consensusEquivalentMember\":"
            << (evidence.consensus_score.consensus_equivalent_member
                    ? "true"
                    : "false")
            << ",\"consensusStragglerMember\":"
            << (evidence.consensus_score.consensus_straggler_member
                    ? "true"
                    : "false")
            << ",\"consensusRebuiltMember\":"
            << (evidence.consensus_score.consensus_rebuilt_member
                    ? "true"
                    : "false")
            << ",\"consensusFallbackRestored\":"
            << (evidence.consensus_score.consensus_fallback_restored
                    ? "true"
                    : "false")
            << ",\"selectedTreeCleanupPrunedOut\":"
            << (evidence.consensus_score.selected_tree_cleanup_pruned_out
                    ? "true"
                    : "false")
            << ",\"selectedTreeCleanupAdded\":"
            << (evidence.consensus_score.selected_tree_cleanup_added
                    ? "true"
                    : "false")
            << ",\"finalDistanceMember\":"
            << (evidence.consensus_score.final_distance_member
                    ? "true"
                    : "false");
        out << ",\"representativeSentinel\":"
            << (evidence.consensus_score.representative_sentinel ? "true" : "false")
            << ",\"otherRepresentativeZero\":"
            << (evidence.consensus_score.other_representative_zero ? "true" : "false")
            << ",\"complete\":"
            << (evidence.consensus_score.complete ? "true" : "false")
            << "},\"postGroupRdpRecheck\":{\"status\":\"";
        if (evidence.post_group_rdp_recheck.representative_skipped) {
          out << "representative-skipped";
        } else if (!evidence.post_group_rdp_recheck.requested) {
          out << "not-in-final-distance-list";
        } else if (!evidence.post_group_rdp_recheck.profile_available) {
          out << "profile-unavailable";
        } else {
          out << "complete";
        }
        out << "\",\"requested\":"
            << (evidence.post_group_rdp_recheck.requested ? "true" : "false")
            << ",\"representativeSkipped\":"
            << (evidence.post_group_rdp_recheck.representative_skipped
                    ? "true"
                    : "false")
            << ",\"profileAvailable\":"
            << (evidence.post_group_rdp_recheck.profile_available
                    ? "true"
                    : "false")
            << ",\"thresholdMode\":\"native-lowp-times-100000-lift\","
               "\"localPValueCutoff\":";
        json::number(
            out, evidence.post_group_rdp_recheck.local_p_value_cutoff);
        out << ",\"emittedSignalCount\":"
            << evidence.post_group_rdp_recheck.emitted_signal_count
            << ",\"candidateSignalCount\":"
            << evidence.post_group_rdp_recheck.candidate_signal_count
            << ",\"overlappingSignalCount\":"
            << evidence.post_group_rdp_recheck.overlapping_signal_count
            << ",\"eventRedetected\":"
            << (evidence.post_group_rdp_recheck.event_redetected
                    ? "true"
                    : "false")
            << ",\"significant\":"
            << (evidence.post_group_rdp_recheck.significant ? "true" : "false")
            << ",\"bestBeginning\":";
        if (evidence.post_group_rdp_recheck.event_redetected) {
          out << evidence.post_group_rdp_recheck.best_beginning;
        } else {
          out << "null";
        }
        out << ",\"bestEnding\":";
        if (evidence.post_group_rdp_recheck.event_redetected) {
          out << evidence.post_group_rdp_recheck.best_ending;
        } else {
          out << "null";
        }
        out << ",\"bestWrapsOrigin\":"
            << (evidence.post_group_rdp_recheck.best_wraps_origin
                    ? "true"
                    : "false")
            << ",\"bestOverlap\":";
        json::number(out, evidence.post_group_rdp_recheck.best_overlap);
        out << ",\"bestLocalPValue\":";
        if (evidence.post_group_rdp_recheck.event_redetected) {
          json::number(out, evidence.post_group_rdp_recheck.best_local_p_value);
        } else {
          out << "null";
        }
        out << ",\"bestCorrectedPValue\":";
        if (evidence.post_group_rdp_recheck.event_redetected) {
          json::number(
              out, evidence.post_group_rdp_recheck.best_corrected_p_value);
        } else {
          out << "null";
        }
        out << "},\"postGroupMaxChiRecheck\":";
        write_maxchi_recheck_json(out, evidence.post_group_maxchi_recheck);
        out << '}';
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
            << ",\"disabledExcluded\":"
            << (evidence.disabled_excluded ? "true" : "false") << '}';
      }
      out << "],\"phylogeneticCorrelationStatus\":\"complete\","
             "\"evidenceSetConsensusComplete\":true,"
             "\"finalTrimDuplicateCorrelationStatus\":\"complete\","
             "\"finalTrimMatrixStatus\":\"complete-active-rff0\","
             "\"finalTrimMembershipStatus\":\"complete-active-rff0\","
             "\"calcMatchStatus\":\"complete-active-rff0-with-bounded-unavailable-cases\","
             "\"consensusScoreStatus\":\"complete-active-rff0\","
             "\"selectedTreeCleanupStatus\":\"complete-active\","
             "\"nativeGroupMembershipComplete\":true,"
             "\"primaryRdpPostGroupRecheckStatus\":\"complete-active\","
             "\"nativePrimaryRdpRecheckComplete\":true,"
             "\"maxChiPostGroupRecheckStatus\":\"source-shaped-strongest-peak-unvalidated\","
             "\"nativeMaxChiFullRecheckComplete\":false,"
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
         "\"Masked rows skip primary triplets but retain secondary RDP and grouping evidence; disabled rows skip event evidence and remain only as bounded phylogenetic context. Trace checks retain structurally matching masked-row profiles even when their corrected p-values are not significant.\","
         "\"Role identification ports MakePhPrScore, leave-one-role-out scores, displacement scores, weighted MakeTrpScore ordering changes, and the corresponding supplied decision-tree contributions.\","
         "\"Manual co-recombinant group edits are preserved separately from the automatic two-of-three set and drive subsequent tract erasure and accepted-event alignment exports.\","
         "\"Breakpoint context ports the supplied RDP CheckEnds path: the current triplet information map after prior erasures, source-shaped input MissingData runs, distinct beginning/ending native ranges, strict linear-edge gates, and the literal wrap comparison. Aggregate warning, erased-event attribution, immediate contact, and nearest information-rich distance remain separate; the original-alignment view brackets expected parent states for manual review.\","
         "\"BURT statistical confidence ports the supplied non-segmented PolishBP/BenHMM path: sorted three-symbol recoding, source circular expansion, 21 seeded Viterbi starts, forward/reverse posteriors, strict 0.995/0.999 range searches, signed interval matching, missing-data repositioning, three-usable-site reversion, and final gap relocation. These source-labelled 99/95 percent ranges remain distinct from parent-state brackets and CheckEnds warnings.\","
         "\"The on-demand graphical tree view reuses compact edge lists from the six reconciliation topologies; arbitrary rooting and weak-branch expansion are display-only and never change event evidence.\","
         "\"Below the supplied 100000-site cutoff, erased tracts re-enter subsequent rounds as gap-padded synthetic fragments with original-sequence provenance; same-origin working copies never occupy one triplet and the retained-fragment cap is explicit.\","
         "\"Distance membership applies the supplied MakeACOR topology-affinity gate, MakeRList dual-correlation override, StripDupInv inverse-only removal, FinalTrim fixed-point and two expansion passes, selected-role pruning, OKSeq 15, all three ConsensusOK list-rebuild passes, shared selected-tree cleanup, and the primary-RDP post-group recheck. Post-group rechecks for the unported non-RDP method families remain a documented parity boundary.\"]}";
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

std::string RdpScanner::event_alignment_json(
    std::uint32_t event_id,
    std::size_t flank_sites,
    std::size_t row_limit,
    std::string& error) const {
  if (event_id >= events_.size()) {
    error = "The selected RDP event does not exist.";
    return {};
  }
  if (alignment_.length == 0 || alignment_.sequence_count() == 0) {
    error = "The original alignment is not available for breakpoint inspection.";
    return {};
  }

  const auto& event = events_[event_id];
  flank_sites = std::clamp<std::size_t>(flank_sites, 5, 100);
  if (options_.circular) {
    // Do not repeat coordinates when a requested window is longer than a
    // short circular alignment.
    flank_sites = std::min(flank_sites, (alignment_.length - 1) / 2);
  }
  row_limit = std::clamp<std::size_t>(row_limit, 3, 64);

  std::vector<std::uint32_t> candidates;
  candidates.reserve(std::min<std::size_t>(alignment_.sequence_count(), row_limit * 2));
  const auto add_candidate = [&](std::uint32_t sequence) {
    if (sequence >= alignment_.sequence_count() ||
        std::find(candidates.begin(), candidates.end(), sequence) != candidates.end()) {
      return;
    }
    candidates.push_back(sequence);
  };
  const auto add_candidates = [&](const std::vector<std::uint32_t>& sequences) {
    for (const std::uint32_t sequence : sequences) add_candidate(sequence);
  };

  // The manual's graphical breakpoint workflow starts with the current
  // recombinant and parents, then the inferred recombinant group and nearby
  // evidence. Preserve that order so row caps never hide a representative.
  add_candidate(event.recombinant);
  add_candidate(event.major_parent);
  add_candidate(event.minor_parent);
  add_candidates(event.co_recombinant_sequences);
  add_candidates(event.automatic_co_recombinant_sequences);
  for (const auto& trace : event.trace_evidence) add_candidate(trace.sequence);
  add_candidates(event.detectable_sequences);
  if (!event.role_hypotheses.empty()) {
    add_candidates(event.role_hypotheses[0].distance_correlation_set);
    add_candidates(event.role_hypotheses[0].phylogenetic_correlation_set);
  }

  const std::size_t candidate_row_count = candidates.size();
  if (candidates.size() > row_limit) candidates.resize(row_limit);

  struct AlignmentPanel {
    const char* name;
    std::size_t center;
    std::vector<std::size_t> coordinates;
    std::size_t center_index = 0;
    std::size_t left_informative_coordinate = 0;
    std::size_t right_informative_coordinate = 0;
    std::size_t transition_span_sites = 0;
    std::vector<std::uint32_t> adjacent_erasure_event_ids;
    std::vector<std::uint32_t> uncertain_erasure_event_ids;
    std::int32_t nearest_erasure_informative_sites = -1;
    BreakpointUncertaintyEvidence uncertainty;
  };
  const auto shifted_coordinate = [&](std::size_t center, std::int64_t offset) {
    const std::int64_t length = static_cast<std::int64_t>(alignment_.length);
    const std::int64_t unwrapped = static_cast<std::int64_t>(center) + offset;
    if (!options_.circular && (unwrapped < 1 || unwrapped > length)) {
      return std::size_t{0};
    }
    return options_.circular
        ? static_cast<std::size_t>(((unwrapped - 1) % length + length) % length + 1)
        : static_cast<std::size_t>(unwrapped);
  };
  const auto parent_pattern = [&](std::size_t coordinate) {
    const std::uint8_t recombinant = alignment_.at(event.recombinant, coordinate - 1);
    const std::uint8_t major = alignment_.at(event.major_parent, coordinate - 1);
    const std::uint8_t minor = alignment_.at(event.minor_parent, coordinate - 1);
    if (recombinant == 0 || major == 0 || minor == 0) return std::uint8_t{0};
    if (recombinant == major && recombinant != minor) return std::uint8_t{1};
    if (recombinant == minor && recombinant != major) return std::uint8_t{2};
    return std::uint8_t{0};
  };
  const auto make_panel = [&](
      const char* name,
      std::size_t center,
      std::uint8_t expected_left,
      std::uint8_t expected_right,
      const std::vector<std::uint32_t>& adjacent_erasure_event_ids,
      const std::vector<std::uint32_t>& uncertain_erasure_event_ids,
      std::int32_t nearest_erasure_informative_sites,
      const BreakpointUncertaintyEvidence& uncertainty) {
    AlignmentPanel panel{name, center, {}, 0};
    panel.adjacent_erasure_event_ids = adjacent_erasure_event_ids;
    panel.uncertain_erasure_event_ids = uncertain_erasure_event_ids;
    panel.nearest_erasure_informative_sites = nearest_erasure_informative_sites;
    panel.uncertainty = uncertainty;
    panel.coordinates.reserve(flank_sites * 2 + 1);
    for (std::int64_t offset = -static_cast<std::int64_t>(flank_sites);
         offset <= static_cast<std::int64_t>(flank_sites);
         ++offset) {
      const std::size_t coordinate = shifted_coordinate(center, offset);
      if (coordinate == 0) continue;
      if (offset == 0) panel.center_index = panel.coordinates.size();
      panel.coordinates.push_back(coordinate);
    }

    // The manual recommends placing each breakpoint between consecutive
    // parent-informative states. Report the closest expected state on either
    // side as a transparent review interval, not as a statistical CI.
    for (std::size_t step = 0; step <= flank_sites; ++step) {
      const std::size_t coordinate = shifted_coordinate(
          center, -static_cast<std::int64_t>(step));
      if (coordinate != 0 && parent_pattern(coordinate) == expected_left) {
        panel.left_informative_coordinate = coordinate;
        break;
      }
    }
    for (std::size_t step = 0; step <= flank_sites; ++step) {
      const std::size_t coordinate = shifted_coordinate(
          center, static_cast<std::int64_t>(step));
      if (coordinate != 0 && parent_pattern(coordinate) == expected_right) {
        panel.right_informative_coordinate = coordinate;
        break;
      }
    }
    if (panel.left_informative_coordinate != 0 &&
        panel.right_informative_coordinate != 0) {
      panel.transition_span_sites = panel.right_informative_coordinate >=
              panel.left_informative_coordinate
          ? panel.right_informative_coordinate - panel.left_informative_coordinate
          : alignment_.length - panel.left_informative_coordinate +
              panel.right_informative_coordinate;
    }

    return panel;
  };
  const std::array<AlignmentPanel, 2> panels{
      make_panel(
          "beginning",
          event.beginning,
          1,
          2,
          event.adjacent_erasure_event_ids[0],
          event.uncertain_erasure_event_ids[0],
          event.nearest_erasure_informative_sites[0],
          event.breakpoint_uncertainty[0]),
      make_panel(
          "ending",
          event.ending,
          2,
          1,
          event.adjacent_erasure_event_ids[1],
          event.uncertain_erasure_event_ids[1],
          event.nearest_erasure_informative_sites[1],
          event.breakpoint_uncertainty[1]),
  };

  const auto contains = [](const std::vector<std::uint32_t>& values, std::uint32_t sequence) {
    return std::find(values.begin(), values.end(), sequence) != values.end();
  };
  const auto is_trace = [&](std::uint32_t sequence) {
    return std::any_of(
        event.trace_evidence.begin(),
        event.trace_evidence.end(),
        [sequence](const TraceEvidence& trace) { return trace.sequence == sequence; });
  };
  const auto role_name = [&](std::uint32_t sequence) -> const char* {
    if (sequence == event.recombinant) return "recombinant";
    if (sequence == event.major_parent) return "major-parent";
    if (sequence == event.minor_parent) return "minor-parent";
    if (contains(event.co_recombinant_sequences, sequence)) return "co-recombinant";
    return "evidence";
  };

  std::ostringstream out;
  out << "{\"eventId\":" << event.id << ",\"alignmentLength\":" << alignment_.length
      << ",\"circular\":" << (options_.circular ? "true" : "false")
      << ",\"fragmentAssisted\":"
      << (event.fragment_assisted_detection ? "true" : "false")
      << ",\"requestedFlankSites\":" << flank_sites
      << ",\"candidateRowCount\":" << candidate_row_count
      << ",\"omittedRowCount\":" << candidate_row_count - candidates.size()
      << ",\"rows\":[";
  for (std::size_t row = 0; row < candidates.size(); ++row) {
    if (row) out << ',';
    const std::uint32_t sequence = candidates[row];
    out << "{\"sequenceIndex\":" << sequence << ",\"sequenceName\":";
    json::string(out, alignment_.names[sequence]);
    out << ",\"role\":\"" << role_name(sequence) << "\",\"masked\":"
        << (sequence_masked(sequence) ? "true" : "false")
        << ",\"disabled\":" << (sequence_disabled(sequence) ? "true" : "false")
        << ",\"currentGroupMember\":"
        << (contains(event.co_recombinant_sequences, sequence) ? "true" : "false")
        << ",\"automaticGroupMember\":"
        << (contains(event.automatic_co_recombinant_sequences, sequence) ? "true" : "false")
        << ",\"trace\":" << (is_trace(sequence) ? "true" : "false")
        << ",\"panels\":[";
    for (std::size_t panel_index = 0; panel_index < panels.size(); ++panel_index) {
      if (panel_index) out << ',';
      std::string bases;
      bases.reserve(panels[panel_index].coordinates.size());
      for (const std::size_t coordinate : panels[panel_index].coordinates) {
        bases.push_back(alignment_.sequences[sequence][coordinate - 1]);
      }
      json::string(out, bases);
    }
    out << "]}";
  }
  out << "],\"panels\":[";
  for (std::size_t panel_index = 0; panel_index < panels.size(); ++panel_index) {
    if (panel_index) out << ',';
    const auto& panel = panels[panel_index];
    const auto& confidence = event.breakpoint_confidence.boundaries[panel_index];
    out << "{\"name\":\"" << panel.name << "\",\"center\":" << panel.center
        << ",\"centerIndex\":" << panel.center_index
        << ",\"statisticalConfidence\":{\"name\":\"" << panel.name
        << "\",\"inputCoordinate\":" << confidence.input_coordinate
        << ",\"polishedCoordinate\":" << confidence.polished_coordinate
        << ",\"intervalAvailable\":"
        << (confidence.interval_available ? "true" : "false")
        << ",\"sourceIntervalContainsInput\":"
        << (confidence.source_interval_contains_input ? "true" : "false")
        << ",\"confidence99\":{\"beginning\":"
        << confidence.confidence_99_beginning << ",\"ending\":"
        << confidence.confidence_99_ending << ",\"wrapsOrigin\":"
        << (confidence.confidence_99_wraps_origin ? "true" : "false")
        << "},\"hmmCoordinate\":" << confidence.hmm_coordinate
        << ",\"confidence95\":{\"beginning\":"
        << confidence.confidence_95_beginning << ",\"ending\":"
        << confidence.confidence_95_ending << ",\"wrapsOrigin\":"
        << (confidence.confidence_95_wraps_origin ? "true" : "false")
        << "},\"repositioned\":"
        << (confidence.repositioned ? "true" : "false")
        << ",\"missingDataAdjusted\":"
        << (confidence.missing_data_adjusted ? "true" : "false")
        << ",\"finalGapAdjusted\":"
        << (confidence.final_gap_adjusted ? "true" : "false")
        << "},\"parentTransition\":{\"expectedLeft\":\""
        << (panel_index == 0 ? "major" : "minor")
        << "\",\"expectedRight\":\"" << (panel_index == 0 ? "minor" : "major")
        << "\",\"leftInformativeCoordinate\":";
    if (panel.left_informative_coordinate == 0) out << "null";
    else out << panel.left_informative_coordinate;
    out << ",\"rightInformativeCoordinate\":";
    if (panel.right_informative_coordinate == 0) out << "null";
    else out << panel.right_informative_coordinate;
    out << ",\"spanSites\":";
    if (panel.left_informative_coordinate == 0 ||
        panel.right_informative_coordinate == 0) {
      out << "null";
    } else {
      out << panel.transition_span_sites;
    }
    out << ",\"supported\":"
        << (panel.left_informative_coordinate != 0 &&
                    panel.right_informative_coordinate != 0
                ? "true"
                : "false")
        << "},\"adjacentErasureEventIds\":[";
    for (std::size_t index = 0; index < panel.adjacent_erasure_event_ids.size(); ++index) {
      if (index) out << ',';
      out << panel.adjacent_erasure_event_ids[index];
    }
    out << "],\"erasureAdjacent\":"
        << (panel.adjacent_erasure_event_ids.empty() ? "false" : "true")
        << ",\"uncertainErasureEventIds\":[";
    for (std::size_t index = 0; index < panel.uncertain_erasure_event_ids.size(); ++index) {
      if (index) out << ',';
      out << panel.uncertain_erasure_event_ids[index];
    }
    out << "],\"erasureWithinRdpWindow\":"
        << (panel.uncertain_erasure_event_ids.empty() ? "false" : "true")
        << ",\"uncertainDueToErasure\":"
        << (panel.uncertain_erasure_event_ids.empty() ? "false" : "true")
        << ",\"nativeCheckEndsApplied\":"
        << (panel.uncertainty.native_check_ends_applied ? "true" : "false")
        << ",\"nativeCheckEndsWarning\":"
        << (panel.uncertainty.native_check_ends_warning ? "true" : "false")
        << ",\"informationProfileAvailable\":"
        << (panel.uncertainty.information_profile_available ? "true" : "false")
        << ",\"inputMissingDataInCheckRange\":"
        << (panel.uncertainty.input_missing_data_in_range ? "true" : "false")
        << ",\"linearEdgeWithinRdpWindow\":"
        << (panel.uncertainty.linear_edge_within_window ? "true" : "false")
        << ",\"nativeCheckRange\":{\"beginning\":"
        << panel.uncertainty.check_range_beginning
        << ",\"ending\":" << panel.uncertainty.check_range_ending
        << ",\"wrapsOrigin\":"
        << (panel.uncertainty.check_range_wraps_origin ? "true" : "false")
        << ",\"coordinateCount\":" << panel.uncertainty.check_coordinate_count
        << '}'
        << ",\"rdpWindowInformativeSites\":" << options_.window_sites
        << ",\"nearestErasureInformativeSites\":";
    if (panel.nearest_erasure_informative_sites < 0) {
      out << "null";
    } else {
      out << panel.nearest_erasure_informative_sites;
    }
    out
        << ",\"coordinates\":[";
    for (std::size_t index = 0; index < panel.coordinates.size(); ++index) {
      if (index) out << ',';
      out << panel.coordinates[index];
    }
    out << "]}";
  }
  out << "]}";
  return out.str();
}

std::string RdpScanner::event_trees_json(
    std::uint32_t event_id,
    std::string& error) const {
  if (event_id >= events_.size()) {
    error = "The selected RDP event does not exist.";
    return {};
  }

  const auto& event = events_[event_id];
  const auto contains = [](const std::vector<std::uint32_t>& values, std::uint32_t sequence) {
    return std::find(values.begin(), values.end(), sequence) != values.end();
  };
  const auto is_trace = [&](std::uint32_t sequence) {
    return std::any_of(
        event.trace_evidence.begin(),
        event.trace_evidence.end(),
        [sequence](const TraceEvidence& trace) { return trace.sequence == sequence; });
  };
  const auto role_name = [&](std::uint32_t sequence) -> const char* {
    if (sequence == event.recombinant) return "recombinant";
    if (sequence == event.major_parent) return "major-parent";
    if (sequence == event.minor_parent) return "minor-parent";
    if (contains(event.co_recombinant_sequences, sequence)) return "co-recombinant";
    return "evidence";
  };
  constexpr std::array<const char*, 6> tree_region_names{
      "5-prime-outside",
      "5-prime-inside",
      "3-prime-outside",
      "3-prime-inside",
      "outside-tract",
      "inside-tract",
  };

  std::ostringstream out;
  out << "{\"eventId\":" << event.id
      << ",\"method\":\"neighbour-joining\""
         ",\"distance\":\"Jukes-Cantor\""
         ",\"displayRooting\":\"arbitrary-internal-node\""
         ",\"bootstrapCollapseCutoff\":0.5"
      << ",\"bootstrapReplicates\":" << kEventTreeBootstrapReplicates
      << ",\"subsampled\":" << (event.tree_panel_subsampled ? "true" : "false")
      << ",\"sequenceCap\":" << kEventTreeSequenceCap
      << ",\"fragmentAssisted\":"
      << (event.fragment_assisted_detection ? "true" : "false")
      << ",\"leaves\":[";
  for (std::size_t index = 0; index < event.tree_panel_leaves.size(); ++index) {
    if (index) out << ',';
    const auto& leaf = event.tree_panel_leaves[index];
    const std::uint32_t original = leaf.original_sequence;
    out << "{\"node\":" << index
        << ",\"workingSequenceIndex\":" << leaf.working_sequence
        << ",\"sequenceIndex\":" << original
        << ",\"sequenceName\":";
    if (original < alignment_.names.size()) json::string(out, alignment_.names[original]);
    else json::string(out, "unmapped sequence");
    out << ",\"fragmentEventId\":";
    if (leaf.fragment_event < 0) out << "null";
    else out << leaf.fragment_event;
    out << ",\"role\":\"" << role_name(original) << "\",\"masked\":"
        << (sequence_masked(original) ? "true" : "false")
        << ",\"disabled\":" << (sequence_disabled(original) ? "true" : "false")
        << ",\"currentGroupMember\":"
        << (contains(event.co_recombinant_sequences, original) ? "true" : "false")
        << ",\"automaticGroupMember\":"
        << (contains(event.automatic_co_recombinant_sequences, original)
                ? "true"
                : "false")
        << ",\"trace\":" << (is_trace(original) ? "true" : "false") << '}';
  }
  out << "],\"regions\":[";
  for (std::size_t region_index = 0; region_index < event.tree_regions.size(); ++region_index) {
    if (region_index) out << ',';
    const auto& region = event.tree_regions[region_index];
    out << "{\"name\":\"" << tree_region_names[region_index]
        << "\",\"sites\":" << region.sites
        << ",\"sequences\":" << region.sequences
        << ",\"usable\":" << (region.usable ? "true" : "false")
        << ",\"nodeCount\":" << region.node_count
        << ",\"root\":" << region.root
        << ",\"bootstrapReplicates\":" << region.bootstrap_replicates
        << ",\"supportedInternalBranches\":"
        << region.supported_internal_branches
        << ",\"internalBranches\":" << region.internal_branches
        << ",\"edges\":[";
    for (std::size_t edge_index = 0; edge_index < region.topology_edges.size(); ++edge_index) {
      if (edge_index) out << ',';
      const auto& edge = region.topology_edges[edge_index];
      out << "{\"from\":" << edge.first << ",\"to\":" << edge.second
          << ",\"length\":";
      json::number(out, edge.length);
      out << ",\"bootstrapSupport\":";
      if (edge.internal) json::number(out, edge.bootstrap_support);
      else out << "null";
      out << ",\"internal\":" << (edge.internal ? "true" : "false")
          << ",\"collapsed\":" << (edge.collapsed ? "true" : "false") << '}';
    }
    out << "]}";
  }
  out << "]}";
  return out.str();
}

std::string RdpScanner::csv() const {
  std::ostringstream out;
  out << "Event,Method,Detection round,Erased original nucleotide sites,Erased working sites,"
         "Fragments added,Fragment-assisted detection,Tract applied during detection,"
         "Recombinant,Major parent,Minor parent,Beginning,Ending,Wraps origin,"
         "BURT confidence status,BURT applied to event,Beginning BURT confidence,"
         "Ending BURT confidence,"
         "Beginning adjacent erased events,Ending adjacent erased events,"
         "Beginning uncertain erased events,Ending uncertain erased events,"
         "Beginning native CheckEnds applied,Ending native CheckEnds applied,"
         "Beginning native CheckEnds warning,Ending native CheckEnds warning,"
         "Beginning input MissingData in check range,Ending input MissingData in check range,"
         "Beginning linear edge within RDP window,Ending linear edge within RDP window,"
         "Beginning native CheckEnds range,Ending native CheckEnds range,"
         "Beginning nearest erased informative sites,Ending nearest erased informative sites,"
         "Best local p-value,Best corrected p-value,"
         "MaxChi triplet recheck status,MaxChi variable sites,MaxChi half-window,"
         "MaxChi grown half-window,MaxChi best pair,MaxChi peak alignment position,"
         "MaxChi maximum chi-square,MaxChi local p-value,MaxChi within-triplet p-value,"
         "MaxChi correction tests,"
         "MaxChi corrected p-value,MaxChi source recheck hit,"
         "Supporting signals,Detectable sequences,"
         "Distance-correlation sequences,FinalTrim duplicate-filtered pairs,"
         "Phylogenetic-correlation sequences,"
         "Automatic two-of-three group,Current co-recombinant group,Group manually adjusted,"
         "Masked trace sequences,Curated masked sequences,Curated disabled sequences,"
         "Correlation sequences tested,Tree panel sequences,Tree panel subsampled,"
         "Recommended recombinant,Recommended major parent,Recommended minor parent,"
         "Role confidence,Weighted role scores,Late-consensus diagnostics,"
         "Review state,Manually adjusted,Native full-consensus parity\n";
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
    std::ostringstream curated_masked;
    std::ostringstream curated_disabled;
    bool first_masked = true;
    bool first_disabled = true;
    for (std::size_t sequence = 0; sequence < alignment_.sequence_count(); ++sequence) {
      if (sequence_masked(static_cast<std::uint32_t>(sequence))) {
        if (!first_masked) curated_masked << ';';
        first_masked = false;
        curated_masked << alignment_.names[sequence];
      }
      if (sequence_disabled(static_cast<std::uint32_t>(sequence))) {
        if (!first_disabled) curated_disabled << ';';
        first_disabled = false;
        curated_disabled << alignment_.names[sequence];
      }
    }
    std::array<std::ostringstream, 2> adjacent_erasure_events;
    std::array<std::ostringstream, 2> uncertain_erasure_events;
    std::array<std::ostringstream, 2> native_check_ranges;
    std::array<std::ostringstream, 2> burt_confidence;
    for (std::size_t boundary = 0;
         boundary < event.adjacent_erasure_event_ids.size();
         ++boundary) {
      for (std::size_t index = 0;
           index < event.adjacent_erasure_event_ids[boundary].size();
           ++index) {
        if (index) adjacent_erasure_events[boundary] << ';';
        adjacent_erasure_events[boundary]
            << event.adjacent_erasure_event_ids[boundary][index] + 1;
      }
      for (std::size_t index = 0;
           index < event.uncertain_erasure_event_ids[boundary].size();
           ++index) {
        if (index) uncertain_erasure_events[boundary] << ';';
        uncertain_erasure_events[boundary]
            << event.uncertain_erasure_event_ids[boundary][index] + 1;
      }
      const auto& uncertainty = event.breakpoint_uncertainty[boundary];
      if (uncertainty.check_coordinate_count > 0) {
        native_check_ranges[boundary]
            << uncertainty.check_range_beginning << '-'
            << uncertainty.check_range_ending;
        if (uncertainty.check_range_wraps_origin) {
          native_check_ranges[boundary] << " (wraps origin)";
        }
        native_check_ranges[boundary]
            << "; " << uncertainty.check_coordinate_count << " coordinates";
      }
      const auto& confidence = event.breakpoint_confidence.boundaries[boundary];
      if (confidence.interval_available) {
        burt_confidence[boundary]
            << "99=" << confidence.confidence_99_beginning << '-'
            << confidence.confidence_99_ending
            << "; HMM=" << confidence.hmm_coordinate
            << "; 95=" << confidence.confidence_95_beginning << '-'
            << confidence.confidence_95_ending
            << "; input=" << confidence.input_coordinate
            << "; polished=" << confidence.polished_coordinate
            << "; input-inside="
            << (confidence.source_interval_contains_input ? "yes" : "no");
        if (confidence.confidence_99_wraps_origin) {
          burt_confidence[boundary] << "; 99-wraps";
        }
        if (confidence.confidence_95_wraps_origin) {
          burt_confidence[boundary] << "; 95-wraps";
        }
        if (confidence.missing_data_adjusted) {
          burt_confidence[boundary] << "; missing-adjusted";
        }
        if (confidence.final_gap_adjusted) {
          burt_confidence[boundary] << "; gap-adjusted";
        }
      }
    }
    std::ostringstream finaltrim_matrix_diagnostics;
    bool first_finaltrim_matrix = true;
    for (const auto& evidence : hypothesis.distance_evidence) {
      const auto& matrix = evidence.final_trim_matrix;
      if (!matrix.applies_to_nonrepresentative) continue;
      if (!first_finaltrim_matrix) finaltrim_matrix_diagnostics << ';';
      first_finaltrim_matrix = false;
      finaltrim_matrix_diagnostics << alignment_.names[evidence.sequence]
          << ":C=" << matrix.collapsed_tree_position_score
          << "/R=" << matrix.raw_tree_position_score
          << "/D=" << matrix.relative_distance_score
          << "/5=";
      if (matrix.breakpoint_score_available[0]) {
        finaltrim_matrix_diagnostics << matrix.breakpoint_distance_scores[0];
      } else {
        finaltrim_matrix_diagnostics << "warning";
      }
      finaltrim_matrix_diagnostics << "/3=";
      if (matrix.breakpoint_score_available[1]) {
        finaltrim_matrix_diagnostics << matrix.breakpoint_distance_scores[1];
      } else {
        finaltrim_matrix_diagnostics << "warning";
      }
      finaltrim_matrix_diagnostics << "/M14="
                                   << matrix.detected_region_distance_score
                                   << "/matrix="
                                   << matrix.active_consensus_matrix_score;
      if (matrix.detected_event_match) {
        finaltrim_matrix_diagnostics << "/event="
                                     << matrix.detected_event_signal_id + 1
                                     << '@' << matrix.detected_event_beginning
                                     << '-' << matrix.detected_event_ending
                                     << "/event-overlap="
                                     << matrix.detected_event_overlap;
      } else {
        finaltrim_matrix_diagnostics << "/event=none";
      }
      if (matrix.detected_region_saturated) {
        finaltrim_matrix_diagnostics << "/M14-saturated";
      }
      const auto& score_evidence = evidence.consensus_score;
      finaltrim_matrix_diagnostics
          << "/CScore=" << score_evidence.final_score
          << "/P0=" << score_evidence.corrected_correlation_p_value
          << "/Pat=" << score_evidence.pattern_score
          << "/RCorrX=" << score_evidence.maximum_direct_correlation
          << "/NSLong=" << score_evidence.source_long_matrix_multiplier
          << "/OK2=" << static_cast<unsigned int>(score_evidence.detectable_set_member)
          << "/OK4=" << static_cast<unsigned int>(score_evidence.initial_rlist_member)
          << "/OK5=" << static_cast<unsigned int>(score_evidence.duplicate_cleaned_member)
          << "/OK6="
          << static_cast<unsigned int>(score_evidence.nearest_nonrecombinant_member)
          << "/OK15=" << static_cast<unsigned int>(score_evidence.final_trim_member)
          << "/consensus="
          << static_cast<unsigned int>(score_evidence.consensus_rebuilt_member)
          << "/final="
          << static_cast<unsigned int>(score_evidence.final_distance_member);
      if (score_evidence.final_trim_first_expansion_added) {
        finaltrim_matrix_diagnostics << "/expand-1";
      }
      if (score_evidence.final_trim_second_expansion_added) {
        finaltrim_matrix_diagnostics << "/expand-2";
      }
      if (score_evidence.selected_role_pruned_out) {
        finaltrim_matrix_diagnostics << "/selected-pruned";
      }
      if (score_evidence.consensus_primary_member) {
        finaltrim_matrix_diagnostics << "/consensus-primary";
      } else if (score_evidence.consensus_equivalent_member) {
        finaltrim_matrix_diagnostics << "/consensus-equivalent";
      } else if (score_evidence.consensus_straggler_member) {
        finaltrim_matrix_diagnostics << "/consensus-straggler";
      }
      if (score_evidence.consensus_fallback_restored) {
        finaltrim_matrix_diagnostics << "/fallback-restored";
      }
      if (score_evidence.selected_tree_cleanup_pruned_out) {
        finaltrim_matrix_diagnostics << "/tree-cleanup-pruned";
      }
      if (score_evidence.selected_tree_cleanup_added) {
        finaltrim_matrix_diagnostics << "/tree-cleanup-added";
      }
      const auto& recheck = evidence.post_group_rdp_recheck;
      if (recheck.representative_skipped) {
        finaltrim_matrix_diagnostics << "/RDP-recheck=representative-skipped";
      } else if (recheck.requested) {
        finaltrim_matrix_diagnostics
            << "/RDP-recheck="
            << (recheck.profile_available ? "complete" : "profile-unavailable")
            << "/RDP-emitted=" << recheck.emitted_signal_count
            << "/RDP-candidate=" << recheck.candidate_signal_count
            << "/RDP-overlap=" << recheck.overlapping_signal_count
            << "/RDP-cutoff=" << recheck.local_p_value_cutoff;
        if (recheck.event_redetected) {
          finaltrim_matrix_diagnostics
              << "/RDP-best=" << recheck.best_beginning << '-'
              << recheck.best_ending
              << "/RDP-best-overlap=" << recheck.best_overlap
              << "/RDP-local-P=" << recheck.best_local_p_value
              << "/RDP-corrected-P=" << recheck.best_corrected_p_value
              << (recheck.significant ? "/RDP-significant" : "/RDP-trace");
        } else {
          finaltrim_matrix_diagnostics << "/RDP-event=not-redetected";
        }
      } else {
        finaltrim_matrix_diagnostics << "/RDP-recheck=not-final-list";
      }
      const auto& maxchi = evidence.post_group_maxchi_recheck;
      finaltrim_matrix_diagnostics
          << "/MaxChi-recheck=" << maxchi_recheck_status(maxchi);
      if (maxchi.profile_available) {
        finaltrim_matrix_diagnostics
            << "/MaxChi-variable-sites=" << maxchi.variable_sites
            << "/MaxChi-half-window=" << maxchi.half_window
            << "/MaxChi-grown-half-window=" << maxchi.grown_half_window;
        if (maxchi.best_pair >= 0) {
          finaltrim_matrix_diagnostics
              << "/MaxChi-pair=" << static_cast<int>(maxchi.best_pair)
              << "/MaxChi-peak=" << maxchi.peak_alignment_position
              << "/MaxChi-chi2=" << maxchi.maximum_chi_square
              << "/MaxChi-local-P=" << maxchi.local_p_value
              << "/MaxChi-within-triplet-P=" << maxchi.within_triplet_p_value
              << "/MaxChi-corrected-P=" << maxchi.corrected_p_value
              << (maxchi.source_recheck_hit ? "/MaxChi-hit" : "/MaxChi-no-hit");
        } else {
          finaltrim_matrix_diagnostics << "/MaxChi-peak=none";
        }
      }
      if (matrix.tree_distance_fallback) finaltrim_matrix_diagnostics << "/JC-fallback";
      const auto& match = evidence.calc_match;
      if (match.available) {
        finaltrim_matrix_diagnostics << "/M17=" << match.regional_match_score
                                     << "/M18="
                                     << static_cast<int>(match.breakpoint_match_class)
                                     << "/BP="
                                     << static_cast<unsigned int>(match.breakpoints_exist[0])
                                     << static_cast<unsigned int>(match.breakpoints_exist[1])
                                     << "/SWin=" << match.smoothing_half_window;
        if (match.topology_filtered) {
          finaltrim_matrix_diagnostics << "/M18raw="
                                       << static_cast<int>(match.raw_breakpoint_match_class)
                                       << "/topology-filtered";
        }
        if (match.topology_distance_fallback) {
          finaltrim_matrix_diagnostics << "/topology-JC-fallback";
        }
      } else {
        finaltrim_matrix_diagnostics << "/M17-M18=unavailable";
      }
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
        << csv_cell(event.breakpoint_confidence.available
                        ? "complete-active-unvalidated"
                        : event.breakpoint_confidence.unavailable_reason)
        << ',' << (event.breakpoint_confidence.applied_to_event ? "yes" : "no") << ','
        << csv_cell(burt_confidence[0].str()) << ','
        << csv_cell(burt_confidence[1].str()) << ','
        << csv_cell(adjacent_erasure_events[0].str()) << ','
        << csv_cell(adjacent_erasure_events[1].str()) << ','
        << csv_cell(uncertain_erasure_events[0].str()) << ','
        << csv_cell(uncertain_erasure_events[1].str()) << ',';
    for (std::size_t boundary = 0; boundary < 2; ++boundary) {
      out << (event.breakpoint_uncertainty[boundary].native_check_ends_applied
                  ? "yes"
                  : "no")
          << ',';
    }
    for (std::size_t boundary = 0; boundary < 2; ++boundary) {
      out << (event.breakpoint_uncertainty[boundary].native_check_ends_warning
                  ? "yes"
                  : "no")
          << ',';
    }
    for (std::size_t boundary = 0; boundary < 2; ++boundary) {
      out << (event.breakpoint_uncertainty[boundary].input_missing_data_in_range
                  ? "yes"
                  : "no")
          << ',';
    }
    for (std::size_t boundary = 0; boundary < 2; ++boundary) {
      out << (event.breakpoint_uncertainty[boundary].linear_edge_within_window
                  ? "yes"
                  : "no")
          << ',';
    }
    out << csv_cell(native_check_ranges[0].str()) << ','
        << csv_cell(native_check_ranges[1].str()) << ',';
    if (event.nearest_erasure_informative_sites[0] < 0) {
      out << ',';
    } else {
      out << event.nearest_erasure_informative_sites[0] << ',';
    }
    if (event.nearest_erasure_informative_sites[1] >= 0) {
      out << event.nearest_erasure_informative_sites[1];
    }
    out << ','
        << std::setprecision(16) << event.best_local_p_value << ','
        << event.best_corrected_p_value << ','
        << maxchi_recheck_status(event.maxchi_triplet_recheck) << ','
        << event.maxchi_triplet_recheck.variable_sites << ','
        << event.maxchi_triplet_recheck.half_window << ','
        << event.maxchi_triplet_recheck.grown_half_window << ',';
    if (event.maxchi_triplet_recheck.best_pair >= 0) {
      out << static_cast<int>(event.maxchi_triplet_recheck.best_pair) << ','
          << event.maxchi_triplet_recheck.peak_alignment_position << ','
          << event.maxchi_triplet_recheck.maximum_chi_square << ','
          << event.maxchi_triplet_recheck.local_p_value << ','
          << event.maxchi_triplet_recheck.within_triplet_p_value << ','
          << event.maxchi_triplet_recheck.correction_tests << ','
          << event.maxchi_triplet_recheck.corrected_p_value << ',';
    } else {
      out << ",,,,,,,";
    }
    out << (event.maxchi_triplet_recheck.source_recheck_hit ? "yes" : "no") << ','
        << csv_cell(support.str()) << ','
        << csv_cell(detectable.str()) << ',' << csv_cell(correlated.str()) << ','
        << csv_cell(duplicate_filtered.str()) << ',' << csv_cell(phylogenetic.str()) << ','
        << csv_cell(automatic_consensus.str()) << ','
        << csv_cell(consensus.str()) << ','
        << (event.group_manual_adjusted ? "yes" : "no") << ','
        << csv_cell(traces.str()) << ','
        << csv_cell(curated_masked.str()) << ','
        << csv_cell(curated_disabled.str()) << ','
        << hypothesis.tested_sequences << ','
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
    out << csv_cell(votes.str()) << ','
        << csv_cell(finaltrim_matrix_diagnostics.str()) << ','
        << review_state_name(event.review_state) << ','
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

std::string RdpScanner::enabled_sequences_fasta() const {
  return curated_sequences_fasta(alignment_, options_.mask, options_.disabled, true);
}

std::string RdpScanner::masked_or_disabled_sequences_fasta() const {
  return curated_sequences_fasta(alignment_, options_.mask, options_.disabled, false);
}

std::string RdpScanner::recombinant_sequences_removed_fasta() const {
  std::vector<std::uint8_t> removed(alignment_.sequence_count(), 0);
  for (const auto& event : events_) {
    if (event.review_state != ReviewState::accepted) continue;
    if (event.recombinant < removed.size()) removed[event.recombinant] = 1;
    for (const std::uint32_t sequence : event.co_recombinant_sequences) {
      if (sequence < removed.size()) removed[sequence] = 1;
    }
  }

  std::ostringstream out;
  for (std::size_t sequence = 0; sequence < alignment_.sequence_count(); ++sequence) {
    if (removed[sequence] != 0) continue;
    write_fasta_record(out, alignment_.names[sequence], alignment_.sequences[sequence]);
  }
  return out.str();
}

std::string RdpScanner::recombinant_columns_removed_fasta() const {
  std::vector<std::int32_t> removal_delta(alignment_.length + 2, 0);
  const auto mark_range = [&](std::size_t first, std::size_t last) {
    if (first < 1 || last < first || last > alignment_.length) return;
    ++removal_delta[first];
    --removal_delta[last + 1];
  };
  for (const auto& event : events_) {
    if (event.review_state != ReviewState::accepted) continue;
    if (event.beginning == event.ending) {
      mark_range(1, alignment_.length);
    } else if (event.wraps_origin) {
      mark_range(event.beginning, alignment_.length);
      mark_range(1, event.ending);
    } else {
      mark_range(event.beginning, event.ending);
    }
  }

  std::vector<std::uint8_t> removed(alignment_.length + 1, 0);
  std::int32_t active = 0;
  std::size_t retained_length = 0;
  for (std::size_t coordinate = 1; coordinate <= alignment_.length; ++coordinate) {
    active += removal_delta[coordinate];
    removed[coordinate] = active > 0 ? 1 : 0;
    if (removed[coordinate] == 0) ++retained_length;
  }
  if (retained_length == 0) return {};

  std::ostringstream out;
  for (std::size_t sequence = 0; sequence < alignment_.sequence_count(); ++sequence) {
    std::string retained;
    retained.reserve(retained_length);
    for (std::size_t coordinate = 1; coordinate <= alignment_.length; ++coordinate) {
      if (removed[coordinate] == 0) {
        retained.push_back(alignment_.sequences[sequence][coordinate - 1]);
      }
    }
    write_fasta_record(out, alignment_.names[sequence], retained);
  }
  return out.str();
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
