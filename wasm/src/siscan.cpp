#include "siscan.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <queue>
#include <tuple>
#include <utility>
#include <vector>

namespace rdp {
namespace {

constexpr double kMinimumProbability = 1.0e-300;
constexpr std::array<std::array<std::uint8_t, 2>, 3> kPairs{{
    {{0, 1}},
    {{0, 2}},
    {{1, 2}},
}};
constexpr std::array<std::array<std::uint8_t, 3>, 2> kPartitionGroups{{
    {{2, 3, 5}},
    {{8, 9, 10}},
}};
constexpr std::array<std::array<std::uint8_t, 3>, 2> kSummedGroups{{
    {{1, 2, 3}},
    {{4, 5, 7}},
}};

// SetUpSiScan's Seq34Conv table. Row zero and column zero are the source's
// unused sentinel cells. Categories 11 and 15 are converted to zero below for
// the default SSVarPFlag=2 profile.
constexpr std::array<std::array<std::uint8_t, 6>, 6> kSeq34Conv{{
    {{0, 0, 0, 0, 0, 0}},
    {{0, 4, 6, 7, 1, 1}},
    {{0, 12, 2, 8, 2, 2}},
    {{0, 13, 9, 3, 3, 3}},
    {{0, 10, 14, 5, 5, 5}},
    {{0, 15, 11, 11, 11, 11}},
}};

class MicrosoftCRand {
 public:
  explicit MicrosoftCRand(std::uint32_t state) : state_(state) {}

  std::uint32_t next() {
    state_ = state_ * 214013U + 2531011U;
    return (state_ >> 16U) & 0x7fffU;
  }

  [[nodiscard]] std::uint32_t state() const { return state_; }

 private:
  std::uint32_t state_;
};

struct TreeHeapEntry {
  float distance = 0.0F;
  std::uint32_t first = 0;
  std::uint32_t second = 0;
  std::uint32_t first_generation = 0;
  std::uint32_t second_generation = 0;
};

struct TreeHeapLater {
  bool operator()(const TreeHeapEntry& left, const TreeHeapEntry& right) const {
    return std::tie(left.distance, left.first, left.second) >
        std::tie(right.distance, right.first, right.second);
  }
};

struct FinalScore {
  std::int8_t pair = -1;
  std::uint8_t score = 0;
  SiscanScoreFamily family = SiscanScoreFamily::unavailable;
  double maximum_z = 0.0;
};

std::uint32_t effective_origin(
    std::size_t sequence,
    const std::vector<std::uint32_t>& origins) {
  return sequence < origins.size()
      ? origins[sequence]
      : static_cast<std::uint32_t>(sequence);
}

bool origin_disabled(
    std::uint32_t origin,
    const std::vector<std::uint8_t>& disabled_origins) {
  return origin < disabled_origins.size() && disabled_origins[origin] != 0;
}

float source_direct_similarity(
    const Alignment& alignment,
    std::size_t first,
    std::size_t second) {
  std::size_t valid = 0;
  std::size_t matches = 0;
  for (std::size_t position = 0; position < alignment.length; ++position) {
    const std::uint8_t a = alignment.at(first, position);
    const std::uint8_t b = alignment.at(second, position);
    if (a == 0 || b == 0) continue;
    ++valid;
    if (a == b) ++matches;
  }
  return valid == 0
      ? 0.0F
      : static_cast<float>(matches) / static_cast<float>(valid);
}

void build_source_wpgma_context(
    const Alignment& alignment,
    const std::vector<std::uint32_t>& origins,
    const std::vector<std::uint8_t>& disabled_origins,
    SiscanWorkspace& workspace) {
  const std::size_t count = alignment.sequence_count();
  std::vector<std::uint32_t> normalized_origins(count);
  std::vector<std::uint8_t> eligible(count, 1);
  for (std::size_t sequence = 0; sequence < count; ++sequence) {
    normalized_origins[sequence] = effective_origin(sequence, origins);
    eligible[sequence] = origin_disabled(
        normalized_origins[sequence], disabled_origins) ? 0 : 1;
  }
  const bool reusable = workspace.context_ready &&
      workspace.context_sequence_count == count &&
      workspace.context_alignment_length == alignment.length &&
      workspace.context_origins == normalized_origins &&
      workspace.context_eligible == eligible;
  if (reusable) return;

  workspace.context_sequence_count = count;
  workspace.context_alignment_length = alignment.length;
  workspace.context_origins = std::move(normalized_origins);
  workspace.context_eligible = std::move(eligible);
  workspace.direct_similarity.assign(count * count, 0.0F);
  workspace.tree_similarity.assign(count * count, 0.0F);
  ++workspace.context_builds;

  const bool supplied_similarity =
      alignment.pair_similarity.size() == count * count;
  for (std::size_t first = 0; first < count; ++first) {
    workspace.direct_similarity[first * count + first] = 1.0F;
    workspace.tree_similarity[first * count + first] = 1.0F;
    for (std::size_t second = first + 1; second < count; ++second) {
      if (workspace.context_eligible[first] == 0 ||
          workspace.context_eligible[second] == 0) {
        continue;
      }
      const float similarity = supplied_similarity
          ? alignment.pair_similarity[first * count + second]
          : source_direct_similarity(alignment, first, second);
      workspace.direct_similarity[first * count + second] = similarity;
      workspace.direct_similarity[second * count + first] = similarity;
      ++workspace.context_pair_comparisons;
    }
  }

  std::vector<float> distances(count * count, 100.0F);
  std::vector<std::uint8_t> active(count, 0);
  std::vector<std::uint32_t> generations(count, 0);
  std::vector<std::vector<std::uint32_t>> clusters(count);
  std::priority_queue<
      TreeHeapEntry,
      std::vector<TreeHeapEntry>,
      TreeHeapLater> heap;
  std::size_t active_count = 0;
  for (std::size_t sequence = 0; sequence < count; ++sequence) {
    if (workspace.context_eligible[sequence] == 0) continue;
    active[sequence] = 1;
    clusters[sequence].push_back(static_cast<std::uint32_t>(sequence));
    ++active_count;
  }
  for (std::size_t first = 0; first < count; ++first) {
    if (active[first] == 0) continue;
    for (std::size_t second = first + 1; second < count; ++second) {
      if (active[second] == 0) continue;
      const float similarity =
          workspace.direct_similarity[first * count + second];
      const float distance = similarity > 0.0F ? 1.0F - similarity : 0.999F;
      distances[first * count + second] = distance;
      distances[second * count + first] = distance;
      heap.push({
          distance,
          static_cast<std::uint32_t>(first),
          static_cast<std::uint32_t>(second),
          0,
          0,
      });
    }
  }

  while (active_count > 1 && !heap.empty()) {
    TreeHeapEntry chosen;
    bool found = false;
    while (!heap.empty()) {
      chosen = heap.top();
      heap.pop();
      const std::size_t first = chosen.first;
      const std::size_t second = chosen.second;
      if (active[first] == 0 || active[second] == 0 ||
          generations[first] != chosen.first_generation ||
          generations[second] != chosen.second_generation ||
          distances[first * count + second] != chosen.distance) {
        continue;
      }
      found = true;
      break;
    }
    if (!found) break;

    const std::size_t first = chosen.first;
    const std::size_t second = chosen.second;
    const float cophenetic_similarity =
        1.0F - chosen.distance / 2.0F;
    for (const std::uint32_t left : clusters[first]) {
      for (const std::uint32_t right : clusters[second]) {
        workspace.tree_similarity[left * count + right] =
            cophenetic_similarity;
        workspace.tree_similarity[right * count + left] =
            cophenetic_similarity;
      }
    }

    std::vector<std::pair<std::size_t, float>> updates;
    updates.reserve(active_count > 2 ? active_count - 2 : 0);
    for (std::size_t other = 0; other < count; ++other) {
      if (active[other] == 0 || other == first || other == second) continue;
      const float first_distance = distances[first * count + other];
      const float second_distance = distances[second * count + other];
      float merged_distance = 0.999999F;
      if (first_distance < 1.0F && second_distance < 1.0F) {
        merged_distance = (first_distance + second_distance) / 2.0F;
      } else if (second_distance < 1.0F) {
        merged_distance = second_distance;
      }
      updates.emplace_back(other, merged_distance);
    }

    clusters[first].insert(
        clusters[first].end(),
        clusters[second].begin(),
        clusters[second].end());
    clusters[second].clear();
    active[second] = 0;
    --active_count;
    ++generations[first];
    ++generations[second];
    ++workspace.context_tree_merges;
    for (const auto& [other, distance] : updates) {
      distances[first * count + other] = distance;
      distances[other * count + first] = distance;
      const std::uint32_t low = static_cast<std::uint32_t>(
          std::min(first, other));
      const std::uint32_t high = static_cast<std::uint32_t>(
          std::max(first, other));
      heap.push({
          distance,
          low,
          high,
          generations[low],
          generations[high],
      });
    }
  }
  workspace.context_ready = true;
}

bool fourth_sequence_compatible(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    std::uint32_t candidate) {
  for (std::size_t position = 0; position < alignment.length; ++position) {
    if (alignment.at(candidate, position) != 0) continue;
    if (alignment.at(triplet[0], position) != 0 &&
        alignment.at(triplet[1], position) != 0 &&
        alignment.at(triplet[2], position) != 0) {
      return false;
    }
  }
  return true;
}

std::int64_t nearest_source_outlier(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    const std::vector<std::uint32_t>& origins,
    const std::vector<std::uint8_t>& disabled_origins,
    SiscanWorkspace& workspace) {
  build_source_wpgma_context(
      alignment, origins, disabled_origins, workspace);
  const std::size_t count = alignment.sequence_count();
  if (count < 4 || std::any_of(
          triplet.begin(), triplet.end(),
          [&](std::uint32_t sequence) { return sequence >= count; })) {
    return -1;
  }
  const auto tree = [&](std::uint32_t first, std::uint32_t second) {
    return workspace.tree_similarity[first * count + second];
  };
  const auto direct = [&](std::uint32_t first, std::uint32_t second) {
    return workspace.direct_similarity[first * count + second];
  };

  std::uint32_t inside = 0;
  std::uint32_t outside = 0;
  if (tree(triplet[0], triplet[1]) > tree(triplet[0], triplet[2])) {
    inside = triplet[0];
    outside = triplet[2];
  } else if (tree(triplet[0], triplet[2]) >
             tree(triplet[0], triplet[1])) {
    inside = triplet[0];
    outside = triplet[1];
  } else {
    inside = triplet[1];
    outside = triplet[0];
  }

  const std::array<std::uint32_t, 3> triplet_origins{
      effective_origin(triplet[0], origins),
      effective_origin(triplet[1], origins),
      effective_origin(triplet[2], origins),
  };
  const auto eligible = [&](std::uint32_t candidate) {
    if (candidate >= count || workspace.context_eligible[candidate] == 0) {
      return false;
    }
    const std::uint32_t origin = effective_origin(candidate, origins);
    if (std::find(
            triplet_origins.begin(), triplet_origins.end(), origin) !=
        triplet_origins.end()) {
      return false;
    }
    return fourth_sequence_compatible(alignment, triplet, candidate);
  };

  const float inside_similarity = tree(outside, inside);
  float current_tree_similarity = 0.0F;
  float current_direct_similarity = -1.0F;
  std::int64_t winner = -1;
  for (std::uint32_t candidate = 0; candidate < count; ++candidate) {
    if (!eligible(candidate) || candidate == outside) continue;
    const float candidate_tree_similarity = tree(outside, candidate);
    if (!(candidate_tree_similarity < inside_similarity) ||
        candidate_tree_similarity < current_tree_similarity) {
      continue;
    }
    const float candidate_direct_similarity = direct(outside, candidate);
    if (candidate_tree_similarity > current_tree_similarity ||
        (candidate_tree_similarity == current_tree_similarity &&
         candidate_direct_similarity > current_direct_similarity)) {
      winner = candidate;
      current_tree_similarity = candidate_tree_similarity;
      current_direct_similarity = candidate_direct_similarity;
    }
  }
  // GetSSOL treats CTD==0 as failure and enters its direct-distance fallback.
  if (winner >= 0 && current_tree_similarity > 0.0F) return winner;

  float minimum_mean_similarity = 1.0F;
  for (std::uint32_t candidate = 0; candidate < count; ++candidate) {
    if (!eligible(candidate)) continue;
    const float mean_similarity =
        (direct(candidate, triplet[0]) + direct(candidate, triplet[1]) +
         direct(candidate, triplet[2])) /
        3.0F;
    if (mean_similarity < minimum_mean_similarity) {
      minimum_mean_similarity = mean_similarity;
      winner = candidate;
    }
  }
  if (winner >= 0) return winner;

  // Source GetSSOL has a final missing-data rescue. The browser makes its
  // implicit intent deterministic: least newly missing sites, then the most
  // distant direct mean, then the earliest working row.
  std::size_t minimum_missing = std::numeric_limits<std::size_t>::max();
  minimum_mean_similarity = 1.0F;
  for (std::uint32_t candidate = 0; candidate < count; ++candidate) {
    if (workspace.context_eligible[candidate] == 0) continue;
    const std::uint32_t origin = effective_origin(candidate, origins);
    if (std::find(
            triplet_origins.begin(), triplet_origins.end(), origin) !=
        triplet_origins.end()) {
      continue;
    }
    std::size_t missing = 0;
    for (std::size_t position = 0; position < alignment.length; ++position) {
      if (alignment.at(candidate, position) == 0 &&
          alignment.at(triplet[0], position) != 0 &&
          alignment.at(triplet[1], position) != 0 &&
          alignment.at(triplet[2], position) != 0) {
        ++missing;
      }
    }
    const float mean_similarity =
        (direct(candidate, triplet[0]) + direct(candidate, triplet[1]) +
         direct(candidate, triplet[2])) /
        3.0F;
    if (missing < minimum_missing ||
        (missing == minimum_missing &&
         mean_similarity < minimum_mean_similarity)) {
      minimum_missing = missing;
      minimum_mean_similarity = mean_similarity;
      winner = candidate;
    }
  }
  return winner;
}

void ensure_vertical_random_prefix(
    std::size_t count,
    std::uint32_t seed,
    SiscanWorkspace& workspace) {
  seed = seed == 0 ? 3U : seed;
  if (workspace.random_seed != seed) {
    workspace.random_seed = seed;
    workspace.random_state = seed;
    workspace.random_started = false;
    workspace.vertical_random_prefix.clear();
  }
  MicrosoftCRand random(workspace.random_state);
  if (!workspace.random_started) {
    (void)random.next();
    workspace.random_started = true;
  }
  if (workspace.vertical_random_prefix.size() < count) {
    ++workspace.random_prefix_extensions;
  }
  while (workspace.vertical_random_prefix.size() < count) {
    const double fraction =
        static_cast<double>(random.next()) / 32767.0;
    const auto value = static_cast<std::uint8_t>(
        std::min(12.0, std::floor(fraction * 11.999 + 1.0)));
    workspace.vertical_random_prefix.push_back(value);
    ++workspace.random_values_generated;
  }
  workspace.random_state = random.state();
}

std::uint8_t default_filtered_category(std::uint8_t category) {
  return category == 11 || category == 15 ? 0 : category;
}

std::uint8_t vertical_category(
    std::uint8_t category,
    std::uint8_t random_value) {
  std::uint8_t mapped = 0;
  if (category == 1 || category == 15) {
    mapped = category;
  } else if (category >= 2 && category <= 7) {
    mapped = static_cast<std::uint8_t>(2 + (random_value - 1) % 6);
  } else if (category >= 8 && category <= 10) {
    mapped = static_cast<std::uint8_t>(8 + (random_value - 1) % 3);
  } else if (category >= 11 && category <= 14) {
    mapped = static_cast<std::uint8_t>(11 + (random_value - 1) % 4);
  }
  return default_filtered_category(mapped);
}

std::uint8_t triplet_category(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    std::size_t position) {
  const std::uint8_t first = alignment.at(triplet[0], position);
  const std::uint8_t second = alignment.at(triplet[1], position);
  const std::uint8_t third = alignment.at(triplet[2], position);
  if (first == 0 || second == 0 || third == 0) return 0;
  if (first == second) return first == third ? 5 : 2;
  if (first == third) return 3;
  if (second == third) return 4;
  return 1;
}

void build_triplet_categories(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    SiscanWorkspace& workspace) {
  workspace.triplet_categories.resize(alignment.length);
  for (std::size_t position = 0; position < alignment.length; ++position) {
    workspace.triplet_categories[position] =
        triplet_category(alignment, triplet, position);
  }
}

std::uint8_t pattern_category(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    std::uint32_t outlier,
    const std::vector<std::uint8_t>& triplet_categories,
    std::size_t position) {
  const std::uint8_t relation = triplet_categories[position];
  if (relation == 0 || relation > 5) return 0;
  const std::uint8_t outlier_state = alignment.at(outlier, position);
  std::uint8_t column = 4;
  if (alignment.at(triplet[0], position) == outlier_state) {
    column = 1;
  } else if (relation == 1) {
    if (alignment.at(triplet[1], position) == outlier_state) column = 2;
    else if (alignment.at(triplet[2], position) == outlier_state) column = 3;
  } else if (relation == 2) {
    if (alignment.at(triplet[2], position) == outlier_state) column = 3;
  } else if (relation == 3 || relation == 4) {
    if (alignment.at(triplet[1], position) == outlier_state) column = 2;
  }
  return default_filtered_category(kSeq34Conv[relation][column]);
}

void count_patterns(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    std::uint32_t outlier,
    std::size_t beginning,
    std::size_t span,
    SiscanWorkspace& workspace) {
  workspace.pattern_counts.fill(0);
  if (alignment.length == 0) return;
  for (std::size_t offset = 0; offset < span; ++offset) {
    const std::size_t position = (beginning + offset) % alignment.length;
    const std::uint8_t category = pattern_category(
        alignment,
        triplet,
        outlier,
        workspace.triplet_categories,
        position);
    ++workspace.pattern_counts[category];
  }
}

double source_z_score(
    const std::vector<std::uint32_t>& scores,
    std::size_t category,
    std::size_t stride,
    std::size_t permutations) {
  const std::size_t offset = category * stride;
  if (permutations == 0 || scores[offset] == 0) return 0.0;
  long double mean = 0.0L;
  long double square_mean = 0.0L;
  for (std::size_t permutation = 1;
       permutation <= permutations;
       ++permutation) {
    const long double value = scores[offset + permutation];
    mean += value;
    square_mean += value * value;
  }
  mean /= static_cast<long double>(permutations);
  square_mean /= static_cast<long double>(permutations);
  long double variance = square_mean - mean * mean;
  if (variance < 0.0L && variance > -1.0e-12L) variance = 0.0L;
  if (!(variance > 0.0L)) return 0.0;
  return static_cast<double>(
      (static_cast<long double>(scores[offset]) - mean) /
      std::sqrt(variance));
}

std::size_t permute_patterns(
    std::size_t permutations,
    std::uint32_t seed,
    SiscanWorkspace& workspace,
    std::array<double, 16>& partition_z,
    std::array<double, 15>& summed_z) {
  const std::size_t stride = permutations + 1;
  std::size_t randomized_patterns = 0;
  for (std::size_t category = 2; category <= 14; ++category) {
    randomized_patterns += workspace.pattern_counts[category];
  }
  const std::size_t draws = randomized_patterns * permutations;
  ensure_vertical_random_prefix(draws, seed, workspace);
  workspace.permutation_scores.assign(16 * stride, 0);
  workspace.summed_scores.assign(15 * stride, 0);
  for (std::size_t category = 0; category <= 15; ++category) {
    workspace.permutation_scores[category * stride] =
        workspace.pattern_counts[category];
  }

  std::size_t cursor = 0;
  for (std::uint8_t category = 2; category <= 14; ++category) {
    for (std::size_t occurrence = 0;
         occurrence < workspace.pattern_counts[category];
         ++occurrence) {
      for (std::size_t permutation = 1;
           permutation <= permutations;
           ++permutation) {
        const std::uint8_t mapped = vertical_category(
            category, workspace.vertical_random_prefix[cursor++]);
        ++workspace.permutation_scores[mapped * stride + permutation];
      }
    }
  }

  for (std::size_t permutation = 0;
       permutation <= permutations;
       ++permutation) {
    const auto p = [&](std::size_t category) {
      return workspace.permutation_scores[category * stride + permutation];
    };
    workspace.summed_scores[1 * stride + permutation] = p(2) + p(7) + p(8);
    workspace.summed_scores[2 * stride + permutation] = p(3) + p(6) + p(9);
    workspace.summed_scores[3 * stride + permutation] = p(4) + p(5) + p(10);
    workspace.summed_scores[4 * stride + permutation] =
        p(2) + p(8) + p(12) + p(11);
    workspace.summed_scores[5 * stride + permutation] =
        p(3) + p(9) + p(13) + p(11);
    workspace.summed_scores[7 * stride + permutation] =
        p(5) + p(10) + p(14) + p(11);
  }

  partition_z.fill(0.0);
  summed_z.fill(0.0);
  for (std::size_t category = 1; category <= 15; ++category) {
    if (category == 11 || category == 15) continue;
    partition_z[category] = source_z_score(
        workspace.permutation_scores, category, stride, permutations);
  }
  for (std::size_t category = 1; category <= 12; ++category) {
    if (category == 6 || category == 8 || category == 9 ||
        category == 10 || category == 11 || category == 12) {
      continue;
    }
    summed_z[category] = source_z_score(
        workspace.summed_scores, category, stride, permutations);
  }
  return draws;
}

std::uint8_t group_pair(
    SiscanScoreFamily family,
    std::uint8_t category) {
  const auto& groups = family == SiscanScoreFamily::partition
      ? kPartitionGroups
      : kSummedGroups;
  for (std::uint8_t group = 0; group < groups.size(); ++group) {
    for (std::uint8_t pair = 0; pair < 3; ++pair) {
      if (groups[group][pair] == category) return pair;
    }
  }
  return 3;
}

std::uint8_t strongest_window_pair(
    const std::array<double, 16>& partition_z,
    const std::array<double, 15>& summed_z,
    std::uint8_t fallback) {
  double maximum = 0.0;
  std::uint8_t score = 0;
  SiscanScoreFamily family = SiscanScoreFamily::unavailable;
  for (std::uint8_t category = 0; category <= 15; ++category) {
    if (partition_z[category] > maximum) {
      maximum = partition_z[category];
      score = category;
      family = SiscanScoreFamily::partition;
    }
  }
  for (std::uint8_t category = 0; category <= 14; ++category) {
    if (summed_z[category] > maximum) {
      maximum = summed_z[category];
      score = category;
      family = SiscanScoreFamily::summed;
    }
  }
  const std::uint8_t pair = group_pair(family, score);
  return pair < 3 ? pair : fallback;
}

FinalScore strongest_region_score(
    std::uint8_t high_pair,
    const std::array<double, 16>& partition_z,
    const std::array<double, 15>& summed_z) {
  std::array<std::uint8_t, 2> low_pairs{};
  std::size_t low_index = 0;
  for (std::uint8_t pair = 0; pair < 3; ++pair) {
    if (pair != high_pair) low_pairs[low_index++] = pair;
  }
  FinalScore best;
  for (std::size_t group = 0; group < 2; ++group) {
    const double high = std::fabs(partition_z[kPartitionGroups[group][high_pair]]);
    const double low_one =
        std::fabs(partition_z[kPartitionGroups[group][low_pairs[0]]]);
    const double low_two =
        std::fabs(partition_z[kPartitionGroups[group][low_pairs[1]]]);
    if (!(high < low_one || high < low_two)) continue;
    for (std::uint8_t pair = 0; pair < 3; ++pair) {
      if (pair == high_pair) continue;
      const std::uint8_t category = kPartitionGroups[group][pair];
      const double value = std::fabs(partition_z[category]);
      if (value > best.maximum_z) {
        best.pair = pair;
        best.score = category;
        best.family = SiscanScoreFamily::partition;
        best.maximum_z = value;
      }
    }
  }
  for (std::size_t group = 0; group < 2; ++group) {
    const double high = std::fabs(summed_z[kSummedGroups[group][high_pair]]);
    const double low_one =
        std::fabs(summed_z[kSummedGroups[group][low_pairs[0]]]);
    const double low_two =
        std::fabs(summed_z[kSummedGroups[group][low_pairs[1]]]);
    if (!(high < low_one || high < low_two)) continue;
    for (std::uint8_t pair = 0; pair < 3; ++pair) {
      if (pair == high_pair) continue;
      const std::uint8_t category = kSummedGroups[group][pair];
      const double value = std::fabs(summed_z[category]);
      if (value > best.maximum_z) {
        best.pair = pair;
        best.score = category;
        best.family = SiscanScoreFamily::summed;
        best.maximum_z = value;
      }
    }
  }
  return best;
}

double source_normal_z(double z) {
  if (!(z > 0.0)) return 1.0;
  long double x = 0.0L;
  long double y = 0.0L;
  if (std::fabs(z) < 5.9999999) {
    y = 0.5L * std::fabs(static_cast<long double>(z));
    if (y >= 3.0L) {
      x = 1.0L;
    } else if (y < 1.0L) {
      const long double w = y * y;
      x = ((((((((0.000124818987L * w - 0.001075204047L) * w +
                       0.005198775019L) *
                          w -
                      0.019198292004L) *
                         w +
                     0.059054035642L) *
                        w -
                    0.151968751364L) *
                       w +
                   0.319152932694L) *
                      w -
                  0.5319230073L) *
                     w +
                 0.797884560593L) *
          y * 2.0L;
    } else {
      y -= 2.0L;
      x = (((((((((((((-0.000045255659L * y + 0.00015252929L) * y -
                            0.000019538132L) *
                               y -
                           0.000676904986L) *
                              y +
                          0.001390604284L) *
                             y -
                         0.00079462082L) *
                            y -
                        0.002034254874L) *
                           y +
                       0.006549791214L) *
                          y -
                      0.010557625006L) *
                         y +
                     0.011630447319L) *
                        y -
                    0.009279453341L) *
                       y +
                   0.005353579108L) *
                      y -
                  0.002141268741L) *
                     y +
                 0.000535310849L) *
                    y +
                0.999936657524L;
    }
    return std::clamp(
        static_cast<double>(std::min(x + 1.0L, 1.0L - x)),
        kMinimumProbability,
        1.0);
  }
  const long double exponent =
      (std::fabs(static_cast<long double>(z)) - 5.999999L) * 10.0L;
  return std::clamp(
      static_cast<double>(1.0e-9L / std::pow(1.6L, exponent)),
      kMinimumProbability,
      1.0);
}

std::uint8_t global_closest_pair(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet) {
  std::array<std::size_t, 3> differences{};
  for (std::size_t position = 0; position < alignment.length; ++position) {
    const std::uint8_t first = alignment.at(triplet[0], position);
    const std::uint8_t second = alignment.at(triplet[1], position);
    const std::uint8_t third = alignment.at(triplet[2], position);
    if (first == 0 || second == 0 || third == 0) continue;
    if (first != second) ++differences[0];
    if (first != third) ++differences[1];
    if (second != third) ++differences[2];
  }
  if (differences[0] < differences[1] &&
      differences[0] < differences[2]) {
    return 0;
  }
  if (differences[1] < differences[0] &&
      differences[1] < differences[2]) {
    return 1;
  }
  return 2;
}

std::array<std::uint8_t, 3> source_relation_categories(
    std::uint8_t high_pair) {
  if (high_pair == 0) return {{2, 3, 4}};
  if (high_pair == 1) return {{3, 2, 4}};
  return {{4, 2, 3}};
}

std::size_t quick_check_window(
    const std::vector<std::uint8_t>& categories,
    std::size_t starting_window,
    std::size_t window_sites,
    std::size_t step_sites,
    std::uint8_t high_pair) {
  const auto relation = source_relation_categories(high_pair);
  const std::size_t length = categories.size();
  const std::size_t window_count =
      length > window_sites
      ? (length - window_sites - 1) / step_sites + 1
      : 0;
  std::size_t total = 0;
  for (std::size_t window = starting_window;
       window < window_count;
       ++window) {
    std::array<std::size_t, 6> tally{};
    const std::size_t beginning = window * step_sites;
    const std::size_t ending = std::min(
        length - 1, beginning + window_sites);
    for (std::size_t position = beginning; position <= ending; ++position) {
      ++tally[categories[position]];
    }
    if (tally[relation[0]] < tally[relation[1]] ||
        tally[relation[0]] < tally[relation[2]]) {
      total = tally[relation[0]] + tally[relation[1]] + tally[relation[2]];
    }
    // Deliberately preserve QuickCheckB's missing braces: these two tests run
    // even when `total` retains its preceding value.
    if (tally[relation[1]] * 2 > total ||
        tally[relation[2]] * 2 > total) {
      return window;
    }
  }
  return window_count;
}

std::array<std::uint8_t, 3> roles_for_pair_change(
    std::uint8_t global_pair,
    std::uint8_t candidate_pair) {
  std::uint8_t recombinant = 0;
  for (std::uint8_t member = 0; member < 3; ++member) {
    const bool global_member =
        kPairs[global_pair][0] == member || kPairs[global_pair][1] == member;
    const bool candidate_member =
        kPairs[candidate_pair][0] == member ||
        kPairs[candidate_pair][1] == member;
    if (global_member && candidate_member) {
      recombinant = member;
      break;
    }
  }
  const std::uint8_t major = kPairs[global_pair][0] == recombinant
      ? kPairs[global_pair][1]
      : kPairs[global_pair][0];
  const std::uint8_t minor = kPairs[candidate_pair][0] == recombinant
      ? kPairs[candidate_pair][1]
      : kPairs[candidate_pair][0];
  return {{recombinant, major, minor}};
}

bool candidate_pair_site(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    std::uint8_t pair,
    std::size_t position) {
  const std::uint8_t first = alignment.at(triplet[kPairs[pair][0]], position);
  const std::uint8_t second = alignment.at(triplet[kPairs[pair][1]], position);
  const std::uint8_t other_local = static_cast<std::uint8_t>(
      3 - kPairs[pair][0] - kPairs[pair][1]);
  const std::uint8_t other = alignment.at(triplet[other_local], position);
  return first != 0 && second != 0 && other != 0 &&
      first == second && first != other;
}

bool shrink_region(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    std::uint8_t pair,
    std::size_t first_window,
    std::size_t last_window,
    std::size_t window_sites,
    std::size_t step_sites,
    std::size_t& beginning,
    std::size_t& ending) {
  if (alignment.length == 0) return false;
  const std::size_t left = std::min(
      alignment.length - 1, first_window * step_sites);
  const std::size_t right = std::min(
      alignment.length - 1,
      last_window * step_sites + window_sites - 1);
  bool found = false;
  for (std::size_t position = left; position <= right; ++position) {
    if (candidate_pair_site(alignment, triplet, pair, position)) {
      beginning = position + 1;
      found = true;
      break;
    }
  }
  if (!found) return false;
  for (std::size_t position = right + 1; position-- > left;) {
    if (candidate_pair_site(alignment, triplet, pair, position)) {
      ending = position + 1;
      break;
    }
  }
  return beginning != ending;
}

std::size_t interval_span(
    std::size_t length,
    std::size_t beginning,
    std::size_t ending) {
  if (length == 0 || beginning == 0 || ending == 0) return 0;
  return beginning <= ending
      ? ending - beginning + 1
      : (length - beginning + 1) + ending;
}

std::size_t interval_informative_sites(
    const std::vector<std::uint8_t>& categories,
    std::size_t beginning,
    std::size_t ending) {
  if (categories.empty()) return 0;
  const std::size_t span = interval_span(categories.size(), beginning, ending);
  std::size_t informative = 0;
  for (std::size_t offset = 0; offset < span; ++offset) {
    const std::size_t position = (beginning - 1 + offset) % categories.size();
    if (categories[position] >= 1 && categories[position] <= 4) ++informative;
  }
  return informative;
}

void probability_chain(
    double maximum_z,
    std::size_t alignment_length,
    std::size_t region_span,
    const SiscanOptions& options,
    double& normal_tail,
    double& region_adjusted,
    double& window_adjusted,
    double& corrected) {
  normal_tail = source_normal_z(maximum_z);
  const std::size_t source_region_length =
      region_span > 1 ? region_span - 1 : 1;
  region_adjusted = std::clamp(
      normal_tail * static_cast<double>(alignment_length) /
          static_cast<double>(source_region_length),
      kMinimumProbability,
      1.0);
  window_adjusted = std::clamp(
      region_adjusted * static_cast<double>(alignment_length) /
          static_cast<double>(std::max<std::size_t>(1, options.window_sites)),
      kMinimumProbability,
      1.0);
  corrected = options.bonferroni
      ? std::clamp(
          window_adjusted * static_cast<double>(
              std::max<std::uint64_t>(1, options.correction_tests)),
          kMinimumProbability,
          1.0)
      : region_adjusted;
}

bool validate_common(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    const SiscanOptions& options) {
  if (alignment.length < 6 || options.window_sites < 5 ||
      options.window_sites >= alignment.length || options.step_sites == 0 ||
      options.step_sites > options.window_sites / 2 ||
      options.scan_permutations == 0 || options.p_value_permutations == 0) {
    return false;
  }
  return std::none_of(
      triplet.begin(), triplet.end(),
      [&](std::uint32_t sequence) {
        return sequence >= alignment.sequence_count();
      });
}

}  // namespace

void siscan_reset_round_context(SiscanWorkspace& workspace) {
  workspace.context_ready = false;
  workspace.context_sequence_count = 0;
  workspace.context_alignment_length = 0;
  workspace.direct_similarity.clear();
  workspace.tree_similarity.clear();
  workspace.context_eligible.clear();
  workspace.context_origins.clear();
}

SiscanDiscoverySummary siscan_discover(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    const std::vector<std::uint8_t>& triplet_missing_data,
    const std::array<double, 3>& pair_similarity,
    const std::vector<std::uint32_t>& origins,
    const std::vector<std::uint8_t>& disabled_origins,
    const SiscanOptions& options,
    SiscanWorkspace& workspace,
    std::vector<SiscanDiscoveryCandidate>& output) {
  (void)triplet_missing_data;
  output.clear();
  SiscanDiscoverySummary summary;
  summary.bonferroni_applied = options.bonferroni;
  summary.correction_tests = std::max<std::uint64_t>(1, options.correction_tests);
  summary.window_sites = options.window_sites;
  summary.step_sites = options.step_sites;
  summary.scan_permutations = options.scan_permutations;
  summary.p_value_permutations = options.p_value_permutations;
  if (!validate_common(alignment, triplet, options)) return summary;

  build_triplet_categories(alignment, triplet, workspace);
  const std::int64_t outlier = nearest_source_outlier(
      alignment,
      triplet,
      origins,
      disabled_origins,
      workspace);
  if (outlier < 0) return summary;
  summary.outlier_available = true;
  summary.outlier_sequence = static_cast<std::uint32_t>(outlier);
  const std::uint8_t high_pair = global_closest_pair(alignment, triplet);
  const std::size_t window_count =
      (alignment.length - options.window_sites - 1) / options.step_sites + 1;
  summary.windows_considered = window_count;
  workspace.window_map.assign(window_count, high_pair);
  for (auto& values : workspace.partition_z) values.assign(window_count, 0.0);
  for (auto& values : workspace.summed_z) values.assign(window_count, 0.0);

  std::size_t next_window = 0;
  while (next_window < window_count) {
    const std::size_t window = options.fast_scan
        ? quick_check_window(
            workspace.triplet_categories,
            next_window,
            options.window_sites,
            options.step_sites,
            high_pair)
        : next_window;
    if (window >= window_count) break;
    count_patterns(
        alignment,
        triplet,
        static_cast<std::uint32_t>(outlier),
        window * options.step_sites,
        options.window_sites,
        workspace);
    std::array<double, 16> partition_z{};
    std::array<double, 15> summed_z{};
    summary.permutation_draws += permute_patterns(
        options.scan_permutations,
        options.random_seed,
        workspace,
        partition_z,
        summed_z);
    for (std::size_t category = 0; category < partition_z.size(); ++category) {
      workspace.partition_z[category][window] = partition_z[category];
    }
    for (std::size_t category = 0; category < summed_z.size(); ++category) {
      workspace.summed_z[category][window] = summed_z[category];
    }
    workspace.window_map[window] = strongest_window_pair(
        partition_z, summed_z, high_pair);
    ++summary.windows_scored;
    next_window = window + 1;
  }
  summary.windows_fast_skipped =
      summary.windows_considered - summary.windows_scored;
  summary.profile_available = summary.windows_scored > 0;
  if (!summary.profile_available) return summary;

  std::size_t window = 0;
  while (window < window_count) {
    if (workspace.window_map[window] == high_pair ||
        workspace.window_map[window] > 2) {
      ++window;
      continue;
    }
    const std::uint8_t candidate_pair = workspace.window_map[window];
    const std::size_t first_window = window;
    while (window + 1 < window_count &&
           workspace.window_map[window + 1] == candidate_pair) {
      ++window;
    }
    const std::size_t last_window = window;
    ++summary.candidate_regions_scored;
    std::size_t beginning = 0;
    std::size_t ending = 0;
    if (!shrink_region(
            alignment,
            triplet,
            candidate_pair,
            first_window,
            last_window,
            options.window_sites,
            options.step_sites,
            beginning,
            ending)) {
      ++window;
      continue;
    }

    const std::size_t span = interval_span(
        alignment.length, beginning, ending);
    count_patterns(
        alignment,
        triplet,
        static_cast<std::uint32_t>(outlier),
        beginning - 1,
        span,
        workspace);
    std::array<double, 16> partition_z{};
    std::array<double, 15> summed_z{};
    const std::size_t draws = permute_patterns(
        options.p_value_permutations,
        options.random_seed,
        workspace,
        partition_z,
        summed_z);
    summary.permutation_draws += draws;
    const FinalScore score = strongest_region_score(
        high_pair, partition_z, summed_z);
    if (score.pair < 0 || !(score.maximum_z > 0.0)) {
      ++window;
      continue;
    }

    SiscanDiscoveryCandidate candidate;
    candidate.beginning = beginning;
    candidate.ending = ending;
    candidate.wraps_origin = beginning > ending;
    candidate.informative_beginning = beginning;
    candidate.informative_ending = ending;
    candidate.global_pair = high_pair;
    candidate.candidate_pair = candidate_pair;
    candidate.outlier_sequence = static_cast<std::uint32_t>(outlier);
    candidate.windows_in_region = last_window - first_window + 1;
    candidate.informative_sites = interval_informative_sites(
        workspace.triplet_categories, beginning, ending);
    candidate.permutation_draws = draws;
    candidate.selected_score = score.score;
    candidate.selected_score_family = score.family;
    candidate.maximum_z = score.maximum_z;
    candidate.pair_similarity = pair_similarity;
    probability_chain(
        score.maximum_z,
        alignment.length,
        span,
        options,
        candidate.normal_tail_p_value,
        candidate.region_length_adjusted_p_value,
        candidate.window_adjusted_p_value,
        candidate.corrected_p_value);
    const auto roles = roles_for_pair_change(high_pair, candidate_pair);
    candidate.recombinant_local = roles[0];
    candidate.major_parent_local = roles[1];
    candidate.minor_parent_local = roles[2];
    if (candidate.corrected_p_value < options.p_value_cutoff) {
      output.push_back(candidate);
      ++summary.emitted_candidates;
    }
    ++window;
  }
  return summary;
}

SiscanRecheckEvidence siscan_recheck(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    const std::vector<std::uint8_t>& triplet_missing_data,
    const std::vector<std::uint32_t>& origins,
    const std::vector<std::uint8_t>& disabled_origins,
    std::size_t event_beginning,
    std::size_t event_ending,
    const SiscanOptions& options,
    SiscanWorkspace& workspace) {
  (void)triplet_missing_data;
  SiscanRecheckEvidence evidence;
  evidence.requested = true;
  evidence.bonferroni_applied = options.bonferroni;
  evidence.correction_tests = std::max<std::uint64_t>(1, options.correction_tests);
  if (!validate_common(alignment, triplet, options) ||
      event_beginning == 0 || event_ending == 0 ||
      event_beginning > alignment.length || event_ending > alignment.length) {
    return evidence;
  }
  build_triplet_categories(alignment, triplet, workspace);
  const std::int64_t outlier = nearest_source_outlier(
      alignment,
      triplet,
      origins,
      disabled_origins,
      workspace);
  if (outlier < 0) return evidence;
  evidence.outlier_available = true;
  evidence.outlier_sequence = static_cast<std::uint32_t>(outlier);
  evidence.global_pair = global_closest_pair(alignment, triplet);
  const std::size_t span = interval_span(
      alignment.length, event_beginning, event_ending);
  count_patterns(
      alignment,
      triplet,
      static_cast<std::uint32_t>(outlier),
      event_beginning - 1,
      span,
      workspace);
  std::array<double, 16> partition_z{};
  std::array<double, 15> summed_z{};
  evidence.permutation_draws = permute_patterns(
      options.p_value_permutations,
      options.random_seed,
      workspace,
      partition_z,
      summed_z);
  const FinalScore score = strongest_region_score(
      evidence.global_pair, partition_z, summed_z);
  evidence.profile_available = true;
  evidence.scored_pair = score.pair;
  evidence.selected_score = score.score;
  evidence.selected_score_family = score.family;
  evidence.maximum_z = score.maximum_z;
  evidence.informative_sites = interval_informative_sites(
      workspace.triplet_categories, event_beginning, event_ending);
  probability_chain(
      score.maximum_z,
      alignment.length,
      span,
      options,
      evidence.normal_tail_p_value,
      evidence.region_length_adjusted_p_value,
      evidence.window_adjusted_p_value,
      evidence.corrected_p_value);
  evidence.source_recheck_hit =
      score.pair >= 0 && evidence.corrected_p_value < options.p_value_cutoff;
  return evidence;
}

SiscanPlotProfile siscan_plot_profile(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    const std::vector<std::uint8_t>& triplet_missing_data,
    const std::vector<std::uint32_t>& origins,
    const std::vector<std::uint8_t>& disabled_origins,
    const SiscanOptions& options,
    SiscanWorkspace& workspace,
    std::int64_t fixed_outlier) {
  (void)triplet_missing_data;
  SiscanPlotProfile plot;
  plot.window_sites = options.window_sites;
  plot.step_sites = options.step_sites;
  plot.permutations = options.scan_permutations;
  if (!validate_common(alignment, triplet, options)) return plot;
  build_triplet_categories(alignment, triplet, workspace);
  std::int64_t outlier = fixed_outlier;
  if (outlier < 0 ||
      static_cast<std::size_t>(outlier) >= alignment.sequence_count()) {
    outlier = nearest_source_outlier(
        alignment,
        triplet,
        origins,
        disabled_origins,
        workspace);
  }
  if (outlier < 0) return plot;
  plot.outlier_sequence = static_cast<std::uint32_t>(outlier);
  const std::size_t window_count =
      (alignment.length - options.window_sites - 1) / options.step_sites + 1;
  plot.coordinates.reserve(window_count);
  for (auto& values : plot.pair_z) values.reserve(window_count);
  for (std::size_t window = 0; window < window_count; ++window) {
    count_patterns(
        alignment,
        triplet,
        static_cast<std::uint32_t>(outlier),
        window * options.step_sites,
        options.window_sites,
        workspace);
    std::array<double, 16> partition_z{};
    std::array<double, 15> summed_z{};
    (void)permute_patterns(
        options.scan_permutations,
        options.random_seed,
        workspace,
        partition_z,
        summed_z);
    plot.coordinates.push_back(
        window * options.step_sites + options.window_sites / 2 + 1);
    for (std::uint8_t pair = 0; pair < 3; ++pair) {
      double strongest = 0.0;
      for (std::size_t group = 0; group < 2; ++group) {
        const double partition =
            partition_z[kPartitionGroups[group][pair]];
        const double summed = summed_z[kSummedGroups[group][pair]];
        if (std::fabs(partition) > std::fabs(strongest)) {
          strongest = partition;
        }
        if (std::fabs(summed) > std::fabs(strongest)) strongest = summed;
      }
      // The manual plot is signed. The compact browser plot keeps one curve
      // per sister pair, selecting the eligible P/S category with the
      // greatest absolute deviation while retaining its original sign.
      plot.pair_z[pair].push_back(strongest);
    }
  }
  plot.available = !plot.coordinates.empty();
  return plot;
}

}  // namespace rdp
