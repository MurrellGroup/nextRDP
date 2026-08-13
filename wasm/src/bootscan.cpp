#include "bootscan.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace rdp {
namespace {

constexpr double kMinimumProbability = 1e-300;
constexpr std::array<std::array<std::size_t, 2>, 3> kPairs{{
    {{0, 1}},
    {{0, 2}},
    {{1, 2}},
}};

class MicrosoftCRand {
 public:
  explicit MicrosoftCRand(std::uint32_t seed) : state_(seed) {}

  std::uint32_t next() {
    state_ = state_ * 214013U + 2531011U;
    return (state_ >> 16U) & 0x7fffU;
  }

 private:
  std::uint32_t state_;
};

std::size_t wrapped_coordinate(long long coordinate, std::size_t length) {
  const long long modulus = static_cast<long long>(length);
  coordinate %= modulus;
  if (coordinate <= 0) coordinate += modulus;
  return static_cast<std::size_t>(coordinate);
}

bool coordinate_in_tract(
    std::size_t coordinate,
    std::size_t beginning,
    std::size_t ending) {
  if (beginning < ending) return coordinate >= beginning && coordinate <= ending;
  return coordinate >= beginning || coordinate <= ending;
}

float source_jukes_cantor_distance(
    std::size_t valid,
    std::size_t differences) {
  if (valid == 0) return 10.0F;
  const float identity = static_cast<float>(valid - differences) /
      static_cast<float>(valid);
  if (!(identity > 0.25F)) return 10.0F;
  const float transformed = (4.0F * identity - 1.0F) / 3.0F;
  if (!(transformed > 0.0F)) return 10.0F;
  const float distance = -0.75F * std::log(transformed);
  return std::isfinite(distance) ? distance : 10.0F;
}

double binomial_tail(
    std::size_t trials,
    std::size_t successes,
    double probability) {
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

std::size_t source_vb_round_nonnegative(double value) {
  if (!(value > 0.0)) return 0;
  const double lower_value = std::floor(value);
  const auto lower = static_cast<std::size_t>(lower_value);
  const double fraction = value - lower_value;
  if (fraction > 0.5 || (fraction == 0.5 && (lower & 1U) != 0)) {
    return lower + 1;
  }
  return lower;
}

void source_seqboot2(
    std::size_t window_sites,
    std::size_t bootstrap_replicates,
    std::uint32_t seed,
    std::vector<std::uint32_t>& weights) {
  const std::size_t stride = bootstrap_replicates;
  weights.assign(window_sites * stride, 0);
  for (std::size_t site = 0; site < window_sites; ++site) {
    weights[site * stride] = 1;
  }

  MicrosoftCRand random(seed == 0 ? 3U : seed);
  (void)random.next();
  (void)random.next();
  constexpr double kRandMaximum = 32767.0;
  const double upper = static_cast<double>(window_sites - 1);
  for (std::size_t draw = 0; draw < window_sites; ++draw) {
    for (std::size_t replicate = 1;
         replicate < bootstrap_replicates;
         ++replicate) {
      const std::size_t sampled_site = std::min(
          window_sites - 1,
          static_cast<std::size_t>(
              static_cast<double>(random.next()) / kRandMaximum * upper));
      ++weights[sampled_site * stride + replicate];
    }
    // SEQBOOT2 allocates indices 0..BSBootReps and generates the otherwise
    // unused final replicate too. FastBootDist/DubToInt later retain only
    // 0..BSBootReps-1, but that discarded draw still advances rand() before
    // the next site. Consume it to keep every retained replicate on the exact
    // supplied Microsoft-CRT random stream.
    (void)random.next();
  }
}

std::size_t rounded_window_index(std::size_t coordinate, std::size_t step) {
  return static_cast<std::size_t>(std::llround(
      static_cast<double>(coordinate) / static_cast<double>(step)));
}

std::uint64_t pair_profile_key(std::uint32_t first, std::uint32_t second) {
  if (second < first) std::swap(first, second);
  return (static_cast<std::uint64_t>(first) << 32U) |
      static_cast<std::uint64_t>(second);
}

void clear_pair_profile_cache(BootscanWorkspace& workspace) {
  workspace.pair_profile_cache.clear();
  workspace.pair_profile_fifo.clear();
  workspace.pair_profile_cache_bytes = 0;
}

bool effective_bootscan_geometry(
    std::size_t length,
    std::size_t requested_window,
    std::size_t requested_step,
    std::size_t& window_sites,
    std::size_t& step_sites,
    std::size_t& source_num_windows) {
  if (length < 10) return false;
  window_sites = requested_window;
  if (window_sites > length / 2) window_sites = length / 2;
  if (window_sites < 5) return false;
  step_sites = requested_step;
  if (step_sites > length / 4) step_sites = length / 4;
  step_sites = std::max<std::size_t>(1, step_sites);
  if (step_sites > window_sites / 2) {
    step_sites = std::max<std::size_t>(1, window_sites / 2);
  }
  source_num_windows = length / step_sites + 2;
  return source_num_windows > 2;
}

void configure_discovery_workspace(
    const Alignment& alignment,
    std::size_t window_sites,
    std::size_t step_sites,
    std::size_t bootstrap_replicates,
    std::uint32_t random_seed,
    std::size_t cache_limit_bytes,
    BootscanWorkspace& workspace) {
  random_seed = random_seed == 0 ? 3U : random_seed;
  const bool changed =
      workspace.discovery_alignment_length != alignment.length ||
      workspace.discovery_window_sites != window_sites ||
      workspace.discovery_step_sites != step_sites ||
      workspace.discovery_bootstrap_replicates != bootstrap_replicates ||
      workspace.discovery_random_seed != random_seed;
  if (changed) {
    clear_pair_profile_cache(workspace);
    source_seqboot2(
        window_sites,
        bootstrap_replicates,
        random_seed,
        workspace.bootstrap_weights);
    workspace.discovery_alignment_length = alignment.length;
    workspace.discovery_window_sites = window_sites;
    workspace.discovery_step_sites = step_sites;
    workspace.discovery_bootstrap_replicates = bootstrap_replicates;
    workspace.discovery_random_seed = random_seed;
  }
  workspace.pair_profile_cache_limit_bytes = cache_limit_bytes;
}

std::shared_ptr<BootscanPairDistanceProfile> pair_distance_profile(
    const Alignment& alignment,
    std::uint32_t first,
    std::uint32_t second,
    std::size_t window_sites,
    std::size_t step_sites,
    std::size_t bootstrap_replicates,
    std::size_t source_num_windows,
    BootscanWorkspace& workspace) {
  const std::uint64_t key = pair_profile_key(first, second);
  const auto cached = workspace.pair_profile_cache.find(key);
  if (cached != workspace.pair_profile_cache.end()) {
    ++workspace.pair_profile_cache_hits;
    return cached->second;
  }
  ++workspace.pair_profile_cache_misses;

  auto profile = std::make_shared<BootscanPairDistanceProfile>();
  profile->distances.assign(
      (source_num_windows + 1) * bootstrap_replicates,
      10.0F);
  workspace.pair_valid_scratch.resize(bootstrap_replicates);
  workspace.pair_difference_scratch.resize(bootstrap_replicates);
  for (std::size_t window = 0; window <= source_num_windows; ++window) {
    std::fill(
        workspace.pair_valid_scratch.begin(),
        workspace.pair_valid_scratch.end(),
        0);
    std::fill(
        workspace.pair_difference_scratch.begin(),
        workspace.pair_difference_scratch.end(),
        0);
    const long long offset = static_cast<long long>(window * step_sites) -
        static_cast<long long>(window_sites / 2);
    for (std::size_t local = 1; local <= window_sites; ++local) {
      const std::size_t coordinate = wrapped_coordinate(
          offset + static_cast<long long>(local), alignment.length);
      const std::uint8_t first_state = alignment.at(first, coordinate - 1);
      const std::uint8_t second_state = alignment.at(second, coordinate - 1);
      if (first_state == 0 || second_state == 0) continue;
      const bool differs = first_state != second_state;
      const std::size_t weight_offset = (local - 1) * bootstrap_replicates;
      for (std::size_t replicate = 0;
           replicate < bootstrap_replicates;
           ++replicate) {
        const std::uint32_t weight =
            workspace.bootstrap_weights[weight_offset + replicate];
        workspace.pair_valid_scratch[replicate] += weight;
        if (differs) workspace.pair_difference_scratch[replicate] += weight;
      }
    }
    const std::size_t distance_offset = window * bootstrap_replicates;
    for (std::size_t replicate = 0;
         replicate < bootstrap_replicates;
         ++replicate) {
      profile->distances[distance_offset + replicate] =
          source_jukes_cantor_distance(
              workspace.pair_valid_scratch[replicate],
              workspace.pair_difference_scratch[replicate]);
    }
  }

  const std::size_t profile_bytes =
      profile->distances.size() * sizeof(profile->distances.front());
  const std::size_t cache_limit = workspace.pair_profile_cache_limit_bytes;
  if (profile_bytes <= cache_limit) {
    while (!workspace.pair_profile_fifo.empty() &&
           workspace.pair_profile_cache_bytes + profile_bytes > cache_limit) {
      const std::uint64_t oldest = workspace.pair_profile_fifo.front();
      workspace.pair_profile_fifo.pop_front();
      const auto entry = workspace.pair_profile_cache.find(oldest);
      if (entry == workspace.pair_profile_cache.end()) continue;
      workspace.pair_profile_cache_bytes -=
          entry->second->distances.size() * sizeof(float);
      workspace.pair_profile_cache.erase(entry);
      ++workspace.pair_profile_cache_evictions;
    }
    workspace.pair_profile_cache.emplace(key, profile);
    workspace.pair_profile_fifo.push_back(key);
    workspace.pair_profile_cache_bytes += profile_bytes;
    workspace.pair_profile_cache_peak_bytes = std::max(
        workspace.pair_profile_cache_peak_bytes,
        workspace.pair_profile_cache_bytes);
  }
  return profile;
}

void combine_pair_profiles(
    const std::array<std::shared_ptr<BootscanPairDistanceProfile>, 3>& profiles,
    const std::vector<std::uint8_t>& triplet_missing_data,
    std::size_t length,
    std::size_t step_sites,
    std::size_t bootstrap_replicates,
    std::size_t source_num_windows,
    BootscanWorkspace& workspace,
    bool& erased_window_filter_applied) {
  workspace.support_counts.assign(source_num_windows + 1, {0, 0, 0});
  workspace.usable_windows.assign(source_num_windows + 1, 0);
  for (std::size_t window = 0; window <= source_num_windows; ++window) {
    const std::size_t offset = window * bootstrap_replicates;
    for (std::size_t replicate = 0;
         replicate < bootstrap_replicates;
         ++replicate) {
      const float first = profiles[0]->distances[offset + replicate];
      const float second = profiles[1]->distances[offset + replicate];
      const float third = profiles[2]->distances[offset + replicate];
      if (!(first < 2.0F && second < 2.0F && third < 2.0F)) continue;
      workspace.usable_windows[window] = 1;
      if (first < second && first < third) {
        ++workspace.support_counts[window][0];
      } else if (second < first && second < third) {
        ++workspace.support_counts[window][1];
      } else if (third < first && third < second) {
        ++workspace.support_counts[window][2];
      }
    }
  }
  erased_window_filter_applied = false;
  if (triplet_missing_data.size() != length) return;
  for (std::size_t coordinate = 1; coordinate <= length; ++coordinate) {
    if (triplet_missing_data[coordinate - 1] == 0) continue;
    const std::size_t window = coordinate / step_sites;
    if (window >= workspace.support_counts.size()) continue;
    workspace.support_counts[window] = {0, 0, 0};
    workspace.usable_windows[window] = 0;
    erased_window_filter_applied = true;
  }
}

void prepare_informative_scores(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    BootscanWorkspace& workspace) {
  workspace.position_to_informative.assign(alignment.length + 1, 0);
  for (auto& scores : workspace.pair_scores) {
    scores.assign(1, 0);
    scores.reserve(alignment.length + 1);
  }
  std::size_t informative_sites = 0;
  for (std::size_t coordinate = 1; coordinate <= alignment.length; ++coordinate) {
    workspace.position_to_informative[coordinate] = informative_sites;
    const std::array<std::uint8_t, 3> states{
        alignment.at(triplet[0], coordinate - 1),
        alignment.at(triplet[1], coordinate - 1),
        alignment.at(triplet[2], coordinate - 1),
    };
    if (states[0] == 0 || states[1] == 0 || states[2] == 0 ||
        (states[0] == states[1] && states[0] == states[2])) {
      continue;
    }
    ++informative_sites;
    for (std::size_t pair = 0; pair < 3; ++pair) {
      const auto members = kPairs[pair];
      workspace.pair_scores[pair].push_back(
          states[members[0]] == states[members[1]] ? 1 : 0);
    }
    workspace.position_to_informative[coordinate] = informative_sites;
  }
}

std::uint8_t local_pair_index(std::uint8_t first, std::uint8_t second) {
  if (second < first) std::swap(first, second);
  if (first == 0 && second == 1) return 0;
  if (first == 0 && second == 2) return 1;
  return 2;
}

std::size_t window_coordinate(
    std::size_t window,
    std::size_t step_sites,
    std::size_t length,
    bool circular) {
  const std::size_t raw = window * step_sites;
  if (!circular) return std::clamp<std::size_t>(raw, 1, length);
  return raw == 0 ? length : ((raw - 1) % length) + 1;
}

bool coordinate_in_candidate_tract(
    std::size_t coordinate,
    std::size_t beginning,
    std::size_t ending,
    bool wraps_origin) {
  if (!wraps_origin) return coordinate >= beginning && coordinate <= ending;
  return coordinate >= beginning || coordinate <= ending;
}

void assign_candidate_roles(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    std::uint8_t supported_pair,
    std::size_t beginning,
    std::size_t ending,
    bool wraps_origin,
    const std::array<double, 3>& pair_similarity,
    BootscanDiscoveryCandidate& candidate) {
  const auto members = kPairs[supported_pair];
  const std::uint8_t outside = static_cast<std::uint8_t>(
      3U - members[0] - members[1]);
  std::array<std::size_t, 2> valid{};
  std::array<std::size_t, 2> identical{};
  for (std::size_t coordinate = 1; coordinate <= alignment.length; ++coordinate) {
    if (coordinate_in_candidate_tract(
            coordinate, beginning, ending, wraps_origin)) {
      continue;
    }
    const std::uint8_t outside_state =
        alignment.at(triplet[outside], coordinate - 1);
    if (outside_state == 0) continue;
    for (std::size_t member = 0; member < 2; ++member) {
      const std::uint8_t state =
          alignment.at(triplet[members[member]], coordinate - 1);
      if (state == 0) continue;
      ++valid[member];
      if (state == outside_state) ++identical[member];
    }
  }
  std::array<double, 2> outside_identity{};
  for (std::size_t member = 0; member < 2; ++member) {
    outside_identity[member] = valid[member] == 0
        ? pair_similarity[local_pair_index(members[member], outside)]
        : static_cast<double>(identical[member]) /
            static_cast<double>(valid[member]);
  }
  const std::size_t recombinant_member =
      outside_identity[1] > outside_identity[0] ? 1 : 0;
  const std::size_t minor_member = 1 - recombinant_member;
  candidate.recombinant_local = members[recombinant_member];
  candidate.major_parent_local = outside;
  candidate.minor_parent_local = members[minor_member];
}

void score_candidate_probability(
    std::uint8_t supported_pair,
    std::size_t beginning,
    std::size_t ending,
    bool wraps_origin,
    bool bonferroni,
    std::uint64_t correction_tests,
    BootscanWorkspace& workspace,
    BootscanDiscoveryCandidate& candidate) {
  const std::size_t informative_sites =
      workspace.pair_scores[supported_pair].empty()
      ? 0
      : workspace.pair_scores[supported_pair].size() - 1;
  candidate.informative_sites = informative_sites;
  if (informative_sites == 0) return;
  const std::size_t beginning_index =
      workspace.position_to_informative[beginning];
  const std::size_t ending_index = workspace.position_to_informative[ending];
  std::size_t tract_sites = 0;
  std::size_t tract_matches = 0;
  std::size_t outside_matches = 0;
  const auto add_scores = [&](std::size_t first,
                              std::size_t last,
                              std::size_t& total) {
    if (first > last || first >= workspace.pair_scores[supported_pair].size()) {
      return;
    }
    last = std::min(last, informative_sites);
    for (std::size_t index = first; index <= last; ++index) {
      total += workspace.pair_scores[supported_pair][index];
    }
  };

  // MakeScoresBS works on XPosDiff indices rather than directly on alignment
  // coordinates. That distinction is observable when a boundary falls on an
  // invariant site: the preceding information-rich score is deliberately
  // included in the tract. Preserve that source convention (including score
  // slot zero's sentinel) instead of silently shifting the probability.
  if (!wraps_origin) {
    tract_sites = ending_index >= beginning_index
        ? ending_index - beginning_index + 1
        : 0;
    add_scores(beginning_index, ending_index, tract_matches);
    if (beginning_index > 0) {
      add_scores(1, beginning_index - 1, outside_matches);
    }
    add_scores(ending_index + 1, informative_sites, outside_matches);
  } else {
    tract_sites = ending_index + (informative_sites - beginning_index) + 1;
    add_scores(1, ending_index, tract_matches);
    add_scores(beginning_index, informative_sites, tract_matches);
    if (ending_index + 1 < beginning_index) {
      add_scores(ending_index + 1, beginning_index - 1, outside_matches);
    }
  }
  candidate.tract_informative_sites = tract_sites;
  candidate.tract_pair_matches = tract_matches;
  candidate.outside_pair_matches = outside_matches;
  if (tract_sites <= 2) return;

  const double independent_probability = tract_sites >= informative_sites
      ? 1.0
      : std::clamp(
          static_cast<double>(tract_matches + outside_matches) /
              static_cast<double>(informative_sites),
          0.0,
          1.0);
  std::size_t probability_length = tract_sites;
  std::size_t probability_matches = tract_matches;
  double exponent = 1.0;
  if (probability_length >= 170) {
    probability_matches = source_vb_round_nonnegative(
        static_cast<double>(probability_matches) * 169.0 /
        static_cast<double>(probability_length));
    exponent = static_cast<double>(probability_length) / 169.0;
    probability_length = 169;
  }
  double probability = binomial_tail(
      probability_length,
      probability_matches,
      independent_probability);
  probability *= static_cast<double>(informative_sites) /
      static_cast<double>(probability_length);
  probability = std::min(1.0, probability);
  probability = probability > 0.0
      ? std::pow(probability, exponent)
      : 0.0;
  if (!std::isfinite(probability)) probability = 1.0;
  candidate.raw_p_value = probability > 0.0
      ? std::max(kMinimumProbability, probability)
      : kMinimumProbability;
  candidate.corrected_p_value = bonferroni
      ? std::min(
          1.0,
          candidate.raw_p_value *
              static_cast<double>(std::max<std::uint64_t>(1, correction_tests)))
      : candidate.raw_p_value;
}

}  // namespace

void bootscan_reset_discovery_cache(BootscanWorkspace& workspace) {
  clear_pair_profile_cache(workspace);
  workspace.discovery_alignment_length = 0;
  workspace.discovery_window_sites = 0;
  workspace.discovery_step_sites = 0;
  workspace.discovery_bootstrap_replicates = 0;
  workspace.discovery_random_seed = 0;
}

BootscanDiscoverySummary bootscan_discover(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    const std::vector<std::uint8_t>& triplet_missing_data,
    const std::array<double, 3>& pair_similarity,
    const BootscanDiscoveryOptions& options,
    BootscanWorkspace& workspace,
    std::vector<BootscanDiscoveryCandidate>& output) {
  output.clear();
  BootscanDiscoverySummary summary;
  summary.bonferroni_applied = options.bonferroni;
  summary.correction_tests = std::max<std::uint64_t>(1, options.correction_tests);
  summary.bootstrap_replicates = options.bootstrap_replicates;
  summary.support_cutoff = options.support_cutoff;
  if (alignment.length == 0 || options.bootstrap_replicates == 0 ||
      std::any_of(triplet.begin(), triplet.end(), [&](std::uint32_t sequence) {
        return sequence >= alignment.sequence_count();
      })) {
    return summary;
  }

  std::size_t window_sites = 0;
  std::size_t step_sites = 0;
  std::size_t source_num_windows = 0;
  if (!effective_bootscan_geometry(
          alignment.length,
          options.window_sites,
          options.step_sites,
          window_sites,
          step_sites,
          source_num_windows)) {
    return summary;
  }
  summary.window_sites = window_sites;
  summary.step_sites = step_sites;
  summary.windows_scanned = source_num_windows + 1;
  configure_discovery_workspace(
      alignment,
      window_sites,
      step_sites,
      options.bootstrap_replicates,
      options.random_seed,
      options.pair_cache_limit_bytes,
      workspace);

  const std::uint64_t hits_before = workspace.pair_profile_cache_hits;
  const std::uint64_t misses_before = workspace.pair_profile_cache_misses;
  const std::uint64_t evictions_before = workspace.pair_profile_cache_evictions;
  std::array<std::shared_ptr<BootscanPairDistanceProfile>, 3> pair_profiles;
  for (std::size_t pair = 0; pair < 3; ++pair) {
    const auto members = kPairs[pair];
    pair_profiles[pair] = pair_distance_profile(
        alignment,
        triplet[members[0]],
        triplet[members[1]],
        window_sites,
        step_sites,
        options.bootstrap_replicates,
        source_num_windows,
        workspace);
  }
  summary.pair_profiles_requested = 3;
  summary.pair_profile_cache_hits = static_cast<std::size_t>(
      workspace.pair_profile_cache_hits - hits_before);
  summary.pair_profile_cache_misses = static_cast<std::size_t>(
      workspace.pair_profile_cache_misses - misses_before);
  summary.pair_profile_cache_evictions = static_cast<std::size_t>(
      workspace.pair_profile_cache_evictions - evictions_before);
  summary.pair_profile_cache_bytes = workspace.pair_profile_cache_bytes;
  summary.pair_profile_cache_peak_bytes = workspace.pair_profile_cache_peak_bytes;
  combine_pair_profiles(
      pair_profiles,
      triplet_missing_data,
      alignment.length,
      step_sites,
      options.bootstrap_replicates,
      source_num_windows,
      workspace,
      summary.erased_window_filter_applied);
  prepare_informative_scores(alignment, triplet, workspace);
  const std::size_t informative_sites = workspace.pair_scores[0].empty()
      ? 0
      : workspace.pair_scores[0].size() - 1;
  // BSXoverR serializes NumWins windows (0..NumWins-1) after creating the
  // Scan indices 1..NumWins-1, including the final retained centered window;
  // that last centered sample is important for joining a circular-origin signal.
  const std::size_t last_window = source_num_windows - 1;
  const std::size_t active_window_count = last_window;
  summary.profile_available = informative_sites > 2 &&
      std::any_of(
          workspace.usable_windows.begin() + 1,
          workspace.usable_windows.begin() + last_window + 1,
          [](std::uint8_t usable) { return usable != 0; });
  if (!summary.profile_available) return summary;

  const double seed_threshold =
      options.support_cutoff * static_cast<double>(options.bootstrap_replicates);
  const std::size_t allowed_gap = std::max<std::size_t>(
      1,
      source_vb_round_nonnegative(
          static_cast<double>(window_sites) /
          static_cast<double>(step_sites)));
  const auto next_window = [&](std::size_t window, int direction) {
    if (direction > 0) {
      if (window < last_window) return window + 1;
      return options.circular ? std::size_t{1} : std::size_t{0};
    }
    if (window > 1) return window - 1;
    return options.circular ? last_window : std::size_t{0};
  };
  const auto dominant = [&](std::size_t window, std::size_t pair, double floor) {
    const auto& support = workspace.support_counts[window];
    return workspace.usable_windows[window] != 0 &&
        support[pair] > support[(pair + 1) % 3] &&
        support[pair] > support[(pair + 2) % 3] &&
        static_cast<double>(support[pair]) >
            floor * static_cast<double>(options.bootstrap_replicates);
  };

  for (std::size_t pair = 0; pair < 3; ++pair) {
    std::vector<std::uint8_t> consumed(last_window + 1, 0);
    for (std::size_t seed = 1; seed <= last_window; ++seed) {
      if (consumed[seed] != 0 ||
          static_cast<double>(workspace.support_counts[seed][pair]) <
              seed_threshold) {
        continue;
      }
      std::size_t first = seed;
      std::size_t last = seed;
      std::size_t region_windows = 1;
      for (const int direction : {-1, 1}) {
        std::size_t cursor = seed;
        std::size_t gap = 0;
        for (std::size_t visited = 1; visited < active_window_count; ++visited) {
          const std::size_t next = next_window(cursor, direction);
          if (next == 0 || next == seed) break;
          const bool full =
              static_cast<double>(workspace.support_counts[next][pair]) >=
              seed_threshold;
          if (full) {
            gap = 0;
          } else {
            if (!dominant(next, pair, 0.4)) break;
            ++gap;
            if (gap > allowed_gap) break;
          }
          cursor = next;
          ++region_windows;
          if (direction < 0) first = cursor;
          else last = cursor;
        }
      }
      if (region_windows >= active_window_count) continue;

      std::size_t cursor = first;
      std::size_t usable_windows = 0;
      double support_total = 0.0;
      double maximum_support = 0.0;
      std::size_t guarded = 0;
      while (guarded++ < active_window_count) {
        consumed[cursor] = 1;
        if (workspace.usable_windows[cursor] != 0) ++usable_windows;
        const double support =
            static_cast<double>(workspace.support_counts[cursor][pair]) /
            static_cast<double>(options.bootstrap_replicates);
        support_total += support;
        maximum_support = std::max(maximum_support, support);
        if (cursor == last) break;
        cursor = next_window(cursor, 1);
        if (cursor == 0) break;
      }
      if (guarded > active_window_count || usable_windows == 0) continue;

      BootscanDiscoveryCandidate candidate;
      candidate.supported_pair = static_cast<std::uint8_t>(pair);
      candidate.candidate_pair = static_cast<std::uint8_t>(pair);
      candidate.beginning = window_coordinate(
          first, step_sites, alignment.length, options.circular);
      candidate.ending = window_coordinate(
          last, step_sites, alignment.length, options.circular);
      candidate.wraps_origin = options.circular && first > last;
      if (!options.circular) {
        if (first == 1) candidate.beginning = 1;
        if (last == last_window) candidate.ending = alignment.length;
      }
      if (candidate.beginning == candidate.ending) continue;
      candidate.windows_scored = guarded;
      candidate.usable_windows = usable_windows;
      candidate.maximum_pair_support = maximum_support;
      candidate.mean_pair_support = support_total /
          static_cast<double>(guarded);
      candidate.bootstrap_p_value = std::max(
          1.0 / (static_cast<double>(options.bootstrap_replicates) * 10.0),
          1.0 - candidate.mean_pair_support);
      candidate.pair_similarity = pair_similarity;
      candidate.erased_window_filter_applied =
          summary.erased_window_filter_applied;
      candidate.informative_beginning =
          workspace.position_to_informative[candidate.beginning];
      candidate.informative_ending =
          workspace.position_to_informative[candidate.ending];
      assign_candidate_roles(
          alignment,
          triplet,
          static_cast<std::uint8_t>(pair),
          candidate.beginning,
          candidate.ending,
          candidate.wraps_origin,
          pair_similarity,
          candidate);
      score_candidate_probability(
          static_cast<std::uint8_t>(pair),
          candidate.beginning,
          candidate.ending,
          candidate.wraps_origin,
          options.bonferroni,
          summary.correction_tests,
          workspace,
          candidate);
      ++summary.candidate_regions_scored;
      if (candidate.tract_informative_sites > 2 &&
          candidate.corrected_p_value > 0.0 &&
          candidate.corrected_p_value < options.p_value_cutoff) {
        output.push_back(std::move(candidate));
      }
    }
  }
  std::stable_sort(
      output.begin(), output.end(),
      [](const BootscanDiscoveryCandidate& first,
         const BootscanDiscoveryCandidate& second) {
        if (first.corrected_p_value != second.corrected_p_value) {
          return first.corrected_p_value < second.corrected_p_value;
        }
        if (first.raw_p_value != second.raw_p_value) {
          return first.raw_p_value < second.raw_p_value;
        }
        if (first.beginning != second.beginning) {
          return first.beginning < second.beginning;
        }
        return first.supported_pair < second.supported_pair;
      });
  summary.emitted_candidates = output.size();
  return summary;
}

BootscanPlotProfile bootscan_plot_profile(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    const std::vector<std::uint8_t>& triplet_missing_data,
    const BootscanDiscoveryOptions& options,
    BootscanWorkspace& workspace) {
  BootscanPlotProfile plot;
  std::size_t window_sites = 0;
  std::size_t step_sites = 0;
  std::size_t source_num_windows = 0;
  if (alignment.length == 0 || options.bootstrap_replicates == 0 ||
      std::any_of(triplet.begin(), triplet.end(), [&](std::uint32_t sequence) {
        return sequence >= alignment.sequence_count();
      }) ||
      !effective_bootscan_geometry(
          alignment.length,
          options.window_sites,
          options.step_sites,
          window_sites,
          step_sites,
          source_num_windows)) {
    return plot;
  }
  configure_discovery_workspace(
      alignment,
      window_sites,
      step_sites,
      options.bootstrap_replicates,
      options.random_seed,
      options.pair_cache_limit_bytes,
      workspace);
  std::array<std::shared_ptr<BootscanPairDistanceProfile>, 3> pair_profiles;
  for (std::size_t pair = 0; pair < 3; ++pair) {
    const auto members = kPairs[pair];
    pair_profiles[pair] = pair_distance_profile(
        alignment,
        triplet[members[0]],
        triplet[members[1]],
        window_sites,
        step_sites,
        options.bootstrap_replicates,
        source_num_windows,
        workspace);
  }
  bool erased = false;
  combine_pair_profiles(
      pair_profiles,
      triplet_missing_data,
      alignment.length,
      step_sites,
      options.bootstrap_replicates,
      source_num_windows,
      workspace,
      erased);
  const std::size_t last_window = source_num_windows - 1;
  plot.window_sites = window_sites;
  plot.step_sites = step_sites;
  plot.bootstrap_replicates = options.bootstrap_replicates;
  plot.coordinates.reserve(last_window);
  for (auto& support : plot.pair_support) support.reserve(last_window);
  for (std::size_t window = 1; window <= last_window; ++window) {
    plot.coordinates.push_back(window_coordinate(
        window, step_sites, alignment.length, options.circular));
    for (std::size_t pair = 0; pair < 3; ++pair) {
      plot.pair_support[pair].push_back(
          static_cast<double>(workspace.support_counts[window][pair]) /
          static_cast<double>(options.bootstrap_replicates));
    }
    if (workspace.usable_windows[window] != 0) plot.available = true;
  }
  return plot;
}

BootscanRecheckEvidence bootscan_recheck(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    const std::vector<std::uint8_t>& triplet_missing_data,
    std::size_t event_beginning,
    std::size_t event_ending,
    const BootscanRecheckOptions& options,
    BootscanWorkspace& workspace) {
  BootscanRecheckEvidence evidence;
  evidence.requested = true;
  evidence.bonferroni_applied = options.bonferroni;
  evidence.correction_tests = std::max<std::uint64_t>(1, options.correction_tests);
  evidence.window_sites = options.window_sites;
  evidence.step_sites = options.step_sites;
  evidence.bootstrap_replicates = options.bootstrap_replicates;
  evidence.random_seed = options.random_seed == 0 ? 3U : options.random_seed;
  evidence.support_cutoff = options.support_cutoff;

  const std::size_t length = alignment.length;
  if (length == 0 || options.bootstrap_replicates == 0 ||
      std::any_of(triplet.begin(), triplet.end(), [&](std::uint32_t sequence) {
        return sequence >= alignment.sequence_count();
      })) {
    return evidence;
  }

  std::size_t window_sites = options.window_sites;
  if (window_sites > length / 2) window_sites = length / 2;
  if (window_sites < 5) return evidence;
  std::size_t step_sites = options.step_sites;
  if (step_sites > length / 4) step_sites = length / 4;
  step_sites = std::max<std::size_t>(1, step_sites);
  if (step_sites > window_sites / 2) {
    step_sites = std::max<std::size_t>(1, window_sites / 2);
  }
  const std::size_t bootstrap_replicates = options.bootstrap_replicates;
  evidence.window_sites = window_sites;
  evidence.step_sites = step_sites;

  event_beginning = std::clamp<std::size_t>(event_beginning, 1, length);
  event_ending = std::clamp<std::size_t>(event_ending, 1, length);

  source_seqboot2(
      window_sites,
      bootstrap_replicates,
      evidence.random_seed,
      workspace.bootstrap_weights);
  const std::size_t replicate_stride = bootstrap_replicates;
  const std::size_t source_num_windows = length / step_sites + 2;
  workspace.support_counts.assign(source_num_windows + 1, {0, 0, 0});
  workspace.usable_windows.assign(source_num_windows + 1, 0);

  std::array<std::vector<std::uint32_t>, 3> valid;
  std::array<std::vector<std::uint32_t>, 3> differences;
  for (std::size_t pair = 0; pair < 3; ++pair) {
    valid[pair].resize(bootstrap_replicates);
    differences[pair].resize(bootstrap_replicates);
  }

  for (std::size_t window = 0; window <= source_num_windows; ++window) {
    for (std::size_t pair = 0; pair < 3; ++pair) {
      std::fill(valid[pair].begin(), valid[pair].end(), 0);
      std::fill(differences[pair].begin(), differences[pair].end(), 0);
    }
    const long long offset = static_cast<long long>(window * step_sites) -
        static_cast<long long>(window_sites / 2);
    for (std::size_t local = 1; local <= window_sites; ++local) {
      const std::size_t coordinate = wrapped_coordinate(
          offset + static_cast<long long>(local), length);
      std::array<std::uint8_t, 3> states{};
      for (std::size_t member = 0; member < 3; ++member) {
        states[member] = alignment.at(triplet[member], coordinate - 1);
      }
      for (std::size_t pair = 0; pair < 3; ++pair) {
        const auto members = kPairs[pair];
        if (states[members[0]] == 0 || states[members[1]] == 0) continue;
        const bool differs = states[members[0]] != states[members[1]];
        const std::size_t weight_offset = (local - 1) * replicate_stride;
        for (std::size_t replicate = 0;
             replicate < bootstrap_replicates;
             ++replicate) {
          const std::uint32_t weight =
              workspace.bootstrap_weights[weight_offset + replicate];
          valid[pair][replicate] += weight;
          if (differs) differences[pair][replicate] += weight;
        }
      }
    }

    for (std::size_t replicate = 0;
         replicate < bootstrap_replicates;
         ++replicate) {
      std::array<float, 3> distances{};
      for (std::size_t pair = 0; pair < 3; ++pair) {
        distances[pair] = source_jukes_cantor_distance(
            valid[pair][replicate], differences[pair][replicate]);
      }
      if (!(distances[0] < 2.0F && distances[1] < 2.0F && distances[2] < 2.0F)) {
        continue;
      }
      workspace.usable_windows[window] = 1;
      if (distances[0] < distances[1] && distances[0] < distances[2]) {
        ++workspace.support_counts[window][0];
      } else if (distances[1] < distances[0] && distances[1] < distances[2]) {
        ++workspace.support_counts[window][1];
      } else if (distances[2] < distances[0] && distances[2] < distances[1]) {
        ++workspace.support_counts[window][2];
      }
    }
  }
  evidence.windows_scanned = source_num_windows + 1;

  if (triplet_missing_data.size() == length) {
    for (std::size_t coordinate = 1; coordinate <= length; ++coordinate) {
      if (triplet_missing_data[coordinate - 1] == 0) continue;
      const std::size_t window = coordinate / step_sites;
      if (window >= workspace.support_counts.size()) continue;
      workspace.support_counts[window] = {0, 0, 0};
      workspace.usable_windows[window] = 0;
      evidence.erased_window_filter_applied = true;
    }
  }

  workspace.event_window_indices.clear();
  const std::size_t last_source_window = length / step_sites > 0
      ? length / step_sites - 1
      : 0;
  if (last_source_window == 0) return evidence;
  std::size_t first_window = rounded_window_index(event_beginning, step_sites);
  std::size_t last_window = rounded_window_index(event_ending, step_sites);
  first_window = first_window > 0 ? first_window - 1 : 1;
  last_window = std::min(last_source_window, last_window + 1);
  first_window = std::clamp<std::size_t>(first_window, 1, last_source_window);
  last_window = std::clamp<std::size_t>(last_window, 1, last_source_window);
  if (event_beginning < event_ending) {
    if (first_window <= last_window) {
      for (std::size_t window = first_window; window <= last_window; ++window) {
        workspace.event_window_indices.push_back(window);
      }
    }
  } else {
    for (std::size_t window = 1; window <= last_window; ++window) {
      workspace.event_window_indices.push_back(window);
    }
    for (std::size_t window = first_window; window <= last_source_window; ++window) {
      workspace.event_window_indices.push_back(window);
    }
  }
  std::sort(workspace.event_window_indices.begin(), workspace.event_window_indices.end());
  workspace.event_window_indices.erase(
      std::unique(
          workspace.event_window_indices.begin(),
          workspace.event_window_indices.end()),
      workspace.event_window_indices.end());
  evidence.event_windows_scored = workspace.event_window_indices.size();

  std::array<std::size_t, 3> region_valid{};
  std::array<std::size_t, 3> region_differences{};
  for (std::size_t coordinate = 1; coordinate <= length; ++coordinate) {
    if (!coordinate_in_tract(coordinate, event_beginning, event_ending)) continue;
    for (std::size_t pair = 0; pair < 3; ++pair) {
      const auto members = kPairs[pair];
      const std::uint8_t first = alignment.at(triplet[members[0]], coordinate - 1);
      const std::uint8_t second = alignment.at(triplet[members[1]], coordinate - 1);
      if (first == 0 || second == 0) continue;
      ++region_valid[pair];
      if (first != second) ++region_differences[pair];
    }
  }
  std::array<float, 3> region_distances{};
  for (std::size_t pair = 0; pair < 3; ++pair) {
    region_distances[pair] = source_jukes_cantor_distance(
        region_valid[pair], region_differences[pair]);
  }
  std::size_t scored_pair = 2;
  if (region_distances[0] > region_distances[1] &&
      region_distances[0] > region_distances[2]) {
    scored_pair = 0;
  } else if (region_distances[1] > region_distances[0] &&
             region_distances[1] > region_distances[2]) {
    scored_pair = 1;
  }
  evidence.scored_pair = static_cast<std::int8_t>(scored_pair);

  double scored_support_total = 0.0;
  for (const std::size_t window : workspace.event_window_indices) {
    if (window >= workspace.support_counts.size()) continue;
    if (workspace.usable_windows[window] != 0) ++evidence.usable_event_windows;
    for (std::size_t pair = 0; pair < 3; ++pair) {
      const double support = static_cast<double>(workspace.support_counts[window][pair]) /
          static_cast<double>(bootstrap_replicates);
      evidence.maximum_pair_support = std::max(evidence.maximum_pair_support, support);
      if (support >= options.support_cutoff) evidence.support_gate_passed = true;
    }
    scored_support_total +=
        static_cast<double>(workspace.support_counts[window][scored_pair]) /
        static_cast<double>(bootstrap_replicates);
  }
  if (!workspace.event_window_indices.empty()) {
    evidence.mean_scored_pair_support = scored_support_total /
        static_cast<double>(workspace.event_window_indices.size());
  }

  workspace.position_to_informative.assign(length + 1, 0);
  for (auto& scores : workspace.pair_scores) {
    scores.assign(1, 0);
    scores.reserve(length + 1);
  }
  std::size_t informative_sites = 0;
  for (std::size_t coordinate = 1; coordinate <= length; ++coordinate) {
    workspace.position_to_informative[coordinate] = informative_sites;
    const std::array<std::uint8_t, 3> states{
        alignment.at(triplet[0], coordinate - 1),
        alignment.at(triplet[1], coordinate - 1),
        alignment.at(triplet[2], coordinate - 1),
    };
    if (states[0] == 0 || states[1] == 0 || states[2] == 0 ||
        (states[0] == states[1] && states[0] == states[2])) {
      continue;
    }
    ++informative_sites;
    for (std::size_t pair = 0; pair < 3; ++pair) {
      const auto members = kPairs[pair];
      workspace.pair_scores[pair].push_back(
          states[members[0]] == states[members[1]] ? 1 : 0);
    }
    workspace.position_to_informative[coordinate] = informative_sites;
  }
  evidence.informative_sites = informative_sites;
  if (informative_sites == 0 || workspace.event_window_indices.empty()) return evidence;

  const std::size_t beginning_index = workspace.position_to_informative[event_beginning];
  const std::size_t ending_index = workspace.position_to_informative[event_ending];
  std::size_t tract_length = 0;
  std::size_t tract_matches = 0;
  std::size_t outside_matches = 0;
  const auto add_scores = [&](std::size_t first, std::size_t last, std::size_t& target) {
    if (first > last || last > informative_sites) return;
    for (std::size_t index = std::max<std::size_t>(1, first);
         index <= last;
         ++index) {
      target += workspace.pair_scores[scored_pair][index];
    }
  };
  if (event_beginning < event_ending) {
    tract_length = ending_index >= beginning_index
        ? ending_index - beginning_index
        : 0;
    add_scores(beginning_index, ending_index, tract_matches);
    if (beginning_index > 1) add_scores(1, beginning_index - 1, outside_matches);
    if (ending_index < informative_sites) {
      add_scores(ending_index + 1, informative_sites, outside_matches);
    }
  } else {
    tract_length = ending_index + (informative_sites - beginning_index);
    if (ending_index + 1 < beginning_index) {
      add_scores(ending_index + 1, beginning_index - 1, outside_matches);
    }
    add_scores(1, ending_index, tract_matches);
    add_scores(beginning_index, informative_sites, tract_matches);
  }
  evidence.tract_informative_sites = tract_length;
  evidence.tract_pair_matches = tract_matches;
  evidence.outside_pair_matches = outside_matches;

  if (tract_length > 2) {
    const double independent_probability = std::clamp(
        static_cast<double>(tract_matches + outside_matches) /
            static_cast<double>(informative_sites),
        0.0,
        1.0);
    std::size_t probability_length = tract_length;
    std::size_t probability_matches = tract_matches;
    double exponent = 1.0;
    if (probability_length >= 170) {
      probability_matches = static_cast<std::size_t>(std::llround(
          static_cast<double>(probability_matches) * 169.0 /
          static_cast<double>(probability_length)));
      exponent = static_cast<double>(probability_length) / 169.0;
      probability_length = 169;
    }
    double probability = binomial_tail(
        probability_length,
        probability_matches,
        independent_probability);
    probability *= static_cast<double>(informative_sites) /
        static_cast<double>(probability_length);
    probability = std::min(1.0, probability);
    probability = probability > 0.0
        ? std::pow(probability, exponent)
        : 0.0;
    if (!std::isfinite(probability)) probability = 1.0;
    evidence.local_p_value = probability > 0.0
        ? std::max(kMinimumProbability, probability)
        : 0.0;
    evidence.corrected_p_value = options.bonferroni
        ? std::min(
            1.0,
            evidence.local_p_value *
                static_cast<double>(evidence.correction_tests))
        : evidence.local_p_value;
  }
  evidence.profile_available = evidence.usable_event_windows > 0 &&
      evidence.tract_informative_sites > 2;
  evidence.source_recheck_hit = evidence.profile_available &&
      evidence.support_gate_passed &&
      evidence.corrected_p_value > 0.0 &&
      evidence.corrected_p_value < options.p_value_cutoff;
  return evidence;
}

}  // namespace rdp
