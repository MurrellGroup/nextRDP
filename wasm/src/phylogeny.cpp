#include "phylogeny.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace rdp {
namespace {

constexpr double kSaturatedDistance = 10.0;
constexpr double kBootstrapCollapseCutoff = 0.5;

struct Edge {
  std::size_t to = 0;
  double length = 0.0;
};

struct NjTree {
  std::size_t leaf_count = 0;
  std::vector<std::vector<Edge>> adjacency;
};

class DeterministicRandom {
 public:
  explicit DeterministicRandom(std::uint64_t seed)
      : state_(seed == 0 ? 0x9e3779b97f4a7c15ULL : seed) {}

  std::size_t index(std::size_t upper_bound) {
    state_ ^= state_ >> 12U;
    state_ ^= state_ << 25U;
    state_ ^= state_ >> 27U;
    const std::uint64_t value = state_ * 0x2545f4914f6cdd1dULL;
    return upper_bound == 0 ? 0 : static_cast<std::size_t>(value % upper_bound);
  }

 private:
  std::uint64_t state_;
};

std::vector<double> distance_matrix(
    const Alignment& alignment,
    const std::vector<std::uint32_t>& sequences,
    const std::vector<std::size_t>& positions) {
  const std::size_t count = sequences.size();
  std::vector<double> distances(count * count, 0.0);
  std::vector<std::uint32_t> comparable(count * count, 0);
  std::vector<std::uint32_t> differences(count * count, 0);
  std::vector<std::uint8_t> states(count, 0);
  for (const std::size_t coordinate : positions) {
    if (coordinate < 1 || coordinate > alignment.length) continue;
    for (std::size_t sequence = 0; sequence < count; ++sequence) {
      states[sequence] = alignment.at(sequences[sequence], coordinate - 1);
    }
    for (std::size_t first = 0; first + 1 < count; ++first) {
      if (states[first] == 0) continue;
      for (std::size_t second = first + 1; second < count; ++second) {
        if (states[second] == 0) continue;
        ++comparable[first * count + second];
        if (states[first] != states[second]) ++differences[first * count + second];
      }
    }
  }
  for (std::size_t first = 0; first < count; ++first) {
    for (std::size_t second = first + 1; second < count; ++second) {
      const std::size_t index = first * count + second;
      double distance = kSaturatedDistance;
      if (comparable[index] >= 10) {
        const double p = static_cast<double>(differences[index]) /
            static_cast<double>(comparable[index]);
        const double argument = 1.0 - (4.0 * p / 3.0);
        if (argument > 0.0) {
          const double corrected = -0.75 * std::log(argument);
          if (std::isfinite(corrected)) distance = std::min(corrected, kSaturatedDistance);
        }
      }
      distances[first * count + second] = distance;
      distances[second * count + first] = distance;
    }
  }
  return distances;
}

void add_edge(NjTree& tree, std::size_t first, std::size_t second, double length) {
  length = std::isfinite(length) ? std::max(0.0, length) : kSaturatedDistance;
  tree.adjacency[first].push_back({second, length});
  tree.adjacency[second].push_back({first, length});
}

NjTree neighbor_joining(const std::vector<double>& leaf_distances, std::size_t leaf_count) {
  NjTree tree;
  tree.leaf_count = leaf_count;
  if (leaf_count == 0) return tree;
  const std::size_t node_capacity = leaf_count == 1 ? 1 : leaf_count * 2 - 2;
  tree.adjacency.resize(node_capacity);
  if (leaf_count == 1) return tree;

  std::vector<double> distances(node_capacity * node_capacity, 0.0);
  for (std::size_t first = 0; first < leaf_count; ++first) {
    for (std::size_t second = 0; second < leaf_count; ++second) {
      distances[first * node_capacity + second] =
          leaf_distances[first * leaf_count + second];
    }
  }

  std::vector<std::size_t> active(leaf_count);
  std::iota(active.begin(), active.end(), 0);
  std::size_t next_node = leaf_count;
  while (active.size() > 2) {
    const std::size_t count = active.size();
    std::vector<double> row_sums(count, 0.0);
    for (std::size_t first = 0; first < count; ++first) {
      for (std::size_t second = 0; second < count; ++second) {
        row_sums[first] += distances[active[first] * node_capacity + active[second]];
      }
    }

    std::size_t best_first = 0;
    std::size_t best_second = 1;
    double best_q = std::numeric_limits<double>::infinity();
    for (std::size_t first = 0; first + 1 < count; ++first) {
      for (std::size_t second = first + 1; second < count; ++second) {
        const double q = static_cast<double>(count - 2) *
                distances[active[first] * node_capacity + active[second]] -
            row_sums[first] - row_sums[second];
        if (q < best_q - 1e-12 ||
            (std::abs(q - best_q) <= 1e-12 &&
             std::pair{active[first], active[second]} <
                 std::pair{active[best_first], active[best_second]})) {
          best_q = q;
          best_first = first;
          best_second = second;
        }
      }
    }

    const std::size_t first_node = active[best_first];
    const std::size_t second_node = active[best_second];
    const double joined_distance = distances[first_node * node_capacity + second_node];
    const double first_length = 0.5 * joined_distance +
        (row_sums[best_first] - row_sums[best_second]) /
            (2.0 * static_cast<double>(count - 2));
    const double second_length = joined_distance - first_length;
    const std::size_t joined_node = next_node++;
    add_edge(tree, joined_node, first_node, first_length);
    add_edge(tree, joined_node, second_node, second_length);

    for (std::size_t index = 0; index < count; ++index) {
      if (index == best_first || index == best_second) continue;
      const std::size_t other = active[index];
      const double distance = 0.5 *
          (distances[first_node * node_capacity + other] +
           distances[second_node * node_capacity + other] - joined_distance);
      distances[joined_node * node_capacity + other] = std::max(0.0, distance);
      distances[other * node_capacity + joined_node] = std::max(0.0, distance);
    }

    std::vector<std::size_t> next_active;
    next_active.reserve(count - 1);
    for (std::size_t index = 0; index < count; ++index) {
      if (index != best_first && index != best_second) next_active.push_back(active[index]);
    }
    next_active.push_back(joined_node);
    std::sort(next_active.begin(), next_active.end());
    active = std::move(next_active);
  }

  add_edge(
      tree,
      active[0],
      active[1],
      distances[active[0] * node_capacity + active[1]]);
  tree.adjacency.resize(next_node);
  return tree;
}

void collect_leaves(
    const NjTree& tree,
    std::size_t node,
    std::size_t parent,
    std::vector<std::size_t>& leaves) {
  if (node < tree.leaf_count) leaves.push_back(node);
  for (const Edge& edge : tree.adjacency[node]) {
    if (edge.to == parent) continue;
    collect_leaves(tree, edge.to, node, leaves);
  }
}

std::string canonical_split_key(std::vector<std::size_t> leaves, std::size_t leaf_count) {
  std::sort(leaves.begin(), leaves.end());
  std::vector<std::size_t> complement;
  complement.reserve(leaf_count - leaves.size());
  std::size_t cursor = 0;
  for (std::size_t leaf = 0; leaf < leaf_count; ++leaf) {
    if (cursor < leaves.size() && leaves[cursor] == leaf) ++cursor;
    else complement.push_back(leaf);
  }
  if (complement.size() < leaves.size() ||
      (complement.size() == leaves.size() && complement < leaves)) {
    leaves = std::move(complement);
  }
  std::string key;
  key.reserve(leaves.size() * 5);
  for (const std::size_t leaf : leaves) {
    key += std::to_string(leaf);
    key.push_back(',');
  }
  return key;
}

std::unordered_map<std::string, std::pair<std::size_t, std::size_t>> internal_splits(
    const NjTree& tree) {
  std::unordered_map<std::string, std::pair<std::size_t, std::size_t>> splits;
  if (tree.leaf_count < 4) return splits;
  for (std::size_t first = 0; first < tree.adjacency.size(); ++first) {
    for (const Edge& edge : tree.adjacency[first]) {
      if (first >= edge.to || first < tree.leaf_count || edge.to < tree.leaf_count) continue;
      std::vector<std::size_t> leaves;
      collect_leaves(tree, first, edge.to, leaves);
      if (leaves.size() < 2 || leaves.size() > tree.leaf_count - 2) continue;
      splits.emplace(canonical_split_key(std::move(leaves), tree.leaf_count),
                     std::pair{first, edge.to});
    }
  }
  return splits;
}

std::vector<double> patristic_distances(
    const NjTree& tree,
    const std::unordered_set<std::string>& collapsed_splits) {
  const std::size_t count = tree.leaf_count;
  std::vector<double> result(count * count, 0.0);
  std::unordered_map<std::uint64_t, std::string> edge_keys;
  if (!collapsed_splits.empty()) {
    for (const auto& [key, endpoints] : internal_splits(tree)) {
      const std::uint64_t edge_key =
          (static_cast<std::uint64_t>(std::min(endpoints.first, endpoints.second)) << 32U) |
          static_cast<std::uint64_t>(std::max(endpoints.first, endpoints.second));
      edge_keys.emplace(edge_key, key);
    }
  }

  for (std::size_t source = 0; source < count; ++source) {
    std::vector<std::pair<std::size_t, double>> stack{{source, 0.0}};
    std::vector<std::size_t> parent(tree.adjacency.size(), tree.adjacency.size());
    parent[source] = source;
    while (!stack.empty()) {
      const auto [node, distance] = stack.back();
      stack.pop_back();
      if (node < count) result[source * count + node] = distance;
      for (const Edge& edge : tree.adjacency[node]) {
        if (parent[edge.to] != tree.adjacency.size()) continue;
        parent[edge.to] = node;
        double length = edge.length;
        const std::uint64_t edge_key =
            (static_cast<std::uint64_t>(std::min(node, edge.to)) << 32U) |
            static_cast<std::uint64_t>(std::max(node, edge.to));
        const auto found = edge_keys.find(edge_key);
        if (found != edge_keys.end() && collapsed_splits.contains(found->second)) length = 0.0;
        stack.emplace_back(edge.to, distance + length);
      }
    }
  }
  return result;
}

}  // namespace

bool TreeRegionEvidence::contains(std::uint32_t sequence) const {
  return sequence < alignment_to_tree.size() && alignment_to_tree[sequence] >= 0;
}

double TreeRegionEvidence::jc(std::uint32_t first, std::uint32_t second) const {
  if (!usable || !contains(first) || !contains(second)) return kSaturatedDistance;
  const std::size_t count = sequences.size();
  if (jukes_cantor.size() != count * count) return kSaturatedDistance;
  return jukes_cantor[
      static_cast<std::size_t>(alignment_to_tree[first]) * count +
      static_cast<std::size_t>(alignment_to_tree[second])];
}

double TreeRegionEvidence::tree(
    std::uint32_t first,
    std::uint32_t second,
    bool collapsed) const {
  if (!usable || !contains(first) || !contains(second)) return kSaturatedDistance;
  const std::size_t count = sequences.size();
  const std::size_t index = static_cast<std::size_t>(alignment_to_tree[first]) * count +
      static_cast<std::size_t>(alignment_to_tree[second]);
  if (raw_patristic.size() != count * count ||
      collapsed_patristic.size() != count * count) {
    return kSaturatedDistance;
  }
  return collapsed ? collapsed_patristic[index] : raw_patristic[index];
}

double jukes_cantor_distance(
    const Alignment& alignment,
    std::uint32_t first,
    std::uint32_t second,
    const std::vector<std::size_t>& one_based_positions,
    std::size_t minimum_comparable_sites) {
  if (first >= alignment.sequence_count() || second >= alignment.sequence_count()) {
    return kSaturatedDistance;
  }
  std::size_t comparable = 0;
  std::size_t differences = 0;
  for (const std::size_t coordinate : one_based_positions) {
    if (coordinate < 1 || coordinate > alignment.length) continue;
    const std::uint8_t first_state = alignment.at(first, coordinate - 1);
    const std::uint8_t second_state = alignment.at(second, coordinate - 1);
    if (first_state == 0 || second_state == 0) continue;
    ++comparable;
    if (first_state != second_state) ++differences;
  }
  if (comparable < minimum_comparable_sites) return kSaturatedDistance;
  const double p = static_cast<double>(differences) / static_cast<double>(comparable);
  const double argument = 1.0 - (4.0 * p / 3.0);
  if (argument <= 0.0) return kSaturatedDistance;
  const double distance = -0.75 * std::log(argument);
  return std::isfinite(distance) ? std::min(distance, kSaturatedDistance) : kSaturatedDistance;
}

TreeRegionEvidence build_tree_region_evidence(
    const Alignment& alignment,
    const std::vector<std::uint32_t>& sequences,
    const std::vector<std::size_t>& one_based_positions,
    std::size_t bootstrap_replicates,
    std::uint64_t seed) {
  TreeRegionEvidence evidence;
  evidence.sequences = sequences;
  evidence.site_count = one_based_positions.size();
  evidence.alignment_to_tree.assign(alignment.sequence_count(), -1);
  for (std::size_t index = 0; index < sequences.size(); ++index) {
    if (sequences[index] < evidence.alignment_to_tree.size()) {
      evidence.alignment_to_tree[sequences[index]] = static_cast<std::int32_t>(index);
    }
  }
  if (sequences.size() < 3 || one_based_positions.size() < 10) return evidence;

  evidence.jukes_cantor = distance_matrix(alignment, sequences, one_based_positions);
  const NjTree base_tree = neighbor_joining(evidence.jukes_cantor, sequences.size());
  const auto base_splits = internal_splits(base_tree);
  evidence.internal_branches = base_splits.size();

  std::unordered_map<std::string, std::size_t> support;
  support.reserve(base_splits.size());
  for (const auto& [key, endpoints] : base_splits) {
    (void)endpoints;
    support.emplace(key, 0);
  }

  DeterministicRandom random(seed);
  std::vector<std::size_t> sampled(one_based_positions.size());
  for (std::size_t replicate = 0; replicate < bootstrap_replicates; ++replicate) {
    for (std::size_t& coordinate : sampled) {
      coordinate = one_based_positions[random.index(one_based_positions.size())];
    }
    const auto replicate_distances = distance_matrix(alignment, sequences, sampled);
    const NjTree replicate_tree = neighbor_joining(replicate_distances, sequences.size());
    const auto replicate_splits = internal_splits(replicate_tree);
    for (auto& [key, count] : support) {
      if (replicate_splits.contains(key)) ++count;
    }
  }
  evidence.bootstrap_replicates = bootstrap_replicates;

  std::unordered_set<std::string> collapsed;
  std::unordered_map<std::uint64_t, double> internal_edge_support;
  std::unordered_set<std::uint64_t> collapsed_edges;
  for (const auto& [key, count] : support) {
    const double fraction = bootstrap_replicates == 0
        ? 1.0
        : static_cast<double>(count) / static_cast<double>(bootstrap_replicates);
    const auto split = base_splits.find(key);
    if (split != base_splits.end()) {
      const auto [first, second] = split->second;
      const std::uint64_t edge_key =
          (static_cast<std::uint64_t>(std::min(first, second)) << 32U) |
          static_cast<std::uint64_t>(std::max(first, second));
      internal_edge_support.emplace(edge_key, fraction);
      if (fraction < kBootstrapCollapseCutoff) collapsed_edges.insert(edge_key);
    }
    if (fraction < kBootstrapCollapseCutoff) {
      collapsed.insert(key);
    } else {
      ++evidence.supported_internal_branches;
    }
  }
  evidence.raw_patristic = patristic_distances(base_tree, {});
  evidence.collapsed_patristic = patristic_distances(base_tree, collapsed);
  evidence.topology_node_count = base_tree.adjacency.size();
  evidence.topology_root = base_tree.leaf_count < base_tree.adjacency.size()
      ? base_tree.leaf_count
      : 0;
  evidence.topology_edges.reserve(
      base_tree.adjacency.empty() ? 0 : base_tree.adjacency.size() - 1);
  for (std::size_t first = 0; first < base_tree.adjacency.size(); ++first) {
    for (const Edge& edge : base_tree.adjacency[first]) {
      if (first >= edge.to) continue;
      const std::uint64_t edge_key =
          (static_cast<std::uint64_t>(first) << 32U) |
          static_cast<std::uint64_t>(edge.to);
      const bool internal = first >= base_tree.leaf_count &&
          edge.to >= base_tree.leaf_count;
      const auto support_value = internal_edge_support.find(edge_key);
      evidence.topology_edges.push_back({
          static_cast<std::uint32_t>(first),
          static_cast<std::uint32_t>(edge.to),
          edge.length,
          support_value == internal_edge_support.end() ? 1.0 : support_value->second,
          internal,
          internal && collapsed_edges.contains(edge_key),
      });
    }
  }
  evidence.usable = true;
  return evidence;
}

}  // namespace rdp
