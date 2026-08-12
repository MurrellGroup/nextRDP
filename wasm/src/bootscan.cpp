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
  }
}

std::size_t rounded_window_index(std::size_t coordinate, std::size_t step) {
  return static_cast<std::size_t>(std::llround(
      static_cast<double>(coordinate) / static_cast<double>(step)));
}

}  // namespace

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
