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
  std::size_t root = 0;
  std::size_t negative_branches_normalized = 0;
  std::vector<std::vector<Edge>> adjacency;
};

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

float source_jukes_cantor_distance(
    std::uint32_t comparable,
    std::uint32_t differences) {
  if (comparable == 0) return static_cast<float>(kSaturatedDistance);
  const float valid = static_cast<float>(comparable);
  const float identity = (valid - static_cast<float>(differences)) / valid;
  if (!(identity > 0.25F)) return static_cast<float>(kSaturatedDistance);
  const float argument = (4.0F * identity - 1.0F) / 3.0F;
  if (!(argument > 0.0F)) return static_cast<float>(kSaturatedDistance);
  const float corrected = -0.75F * std::log(argument);
  return std::isfinite(corrected)
      ? corrected
      : static_cast<float>(kSaturatedDistance);
}

std::vector<std::uint32_t> source_event_bootstrap_weights_impl(
    std::size_t site_count,
    std::size_t bootstrap_replicates,
    std::uint32_t seed) {
  const std::size_t stride = bootstrap_replicates + 1;
  std::vector<std::uint32_t> weights(site_count * stride, 0);
  for (std::size_t site = 0; site < site_count; ++site) {
    weights[site * stride] = 1;
  }
  if (site_count == 0 || bootstrap_replicates == 0) return weights;

  MicrosoftCRand random(seed == 0 ? 3U : seed);
  // Supplied SEQBOOT2 calls rand() once immediately after srand(), then once
  // again before its site/replicate draw loop. Both discarded values are part
  // of the stream observed by every retained event-tree replicate.
  (void)random.next();
  (void)random.next();
  constexpr double kRandMaximum = 32767.0;
  const double upper = static_cast<double>(site_count - 1);
  for (std::size_t draw = 0; draw < site_count; ++draw) {
    for (std::size_t replicate = 1;
         replicate <= bootstrap_replicates;
         ++replicate) {
      const std::size_t sampled_site = std::min(
          site_count - 1,
          static_cast<std::size_t>(
              static_cast<double>(random.next()) / kRandMaximum * upper));
      ++weights[sampled_site * stride + replicate];
    }
  }
  return weights;
}

std::vector<float> weighted_distance_matrix(
    const Alignment& alignment,
    const std::vector<std::uint32_t>& sequences,
    const std::vector<std::size_t>& positions,
    const std::vector<std::uint32_t>& weights,
    std::size_t weight_stride,
    std::size_t replicate) {
  const std::size_t count = sequences.size();
  std::vector<float> distances(count * count, 0.0F);
  std::vector<std::uint32_t> comparable(count * count, 0);
  std::vector<std::uint32_t> differences(count * count, 0);
  std::vector<std::uint8_t> states(count, 0);
  for (std::size_t site = 0; site < positions.size(); ++site) {
    const std::size_t coordinate = positions[site];
    if (coordinate < 1 || coordinate > alignment.length) continue;
    const std::uint32_t weight = weights.empty()
        ? 1
        : weights[site * weight_stride + replicate];
    if (weight == 0) continue;
    for (std::size_t sequence = 0; sequence < count; ++sequence) {
      states[sequence] = alignment.at(sequences[sequence], coordinate - 1);
    }
    for (std::size_t first = 0; first + 1 < count; ++first) {
      if (states[first] == 0) continue;
      for (std::size_t second = first + 1; second < count; ++second) {
        if (states[second] == 0) continue;
        comparable[first * count + second] += weight;
        if (states[first] != states[second]) {
          differences[first * count + second] += weight;
        }
      }
    }
  }
  for (std::size_t first = 0; first < count; ++first) {
    for (std::size_t second = first + 1; second < count; ++second) {
      const std::size_t index = first * count + second;
      const float distance = source_jukes_cantor_distance(
          comparable[index], differences[index]);
      distances[first * count + second] = distance;
      distances[second * count + first] = distance;
    }
  }
  return distances;
}

double source_serialized_branch_length(float source_length, NjTree& tree) {
  if (!std::isfinite(source_length)) return kSaturatedDistance;
  if (source_length < 0.0F) ++tree.negative_branches_normalized;
  // Clearcut's Newick writer suppresses the minus sign and truncates at five
  // decimal places before Tree2Array reads the branch back.
  const double magnitude = std::abs(static_cast<double>(source_length));
  return std::trunc(magnitude * 100000.0) / 100000.0;
}

void add_source_edge(
    NjTree& tree,
    std::size_t first,
    std::size_t second,
    float source_length) {
  const double length = source_serialized_branch_length(source_length, tree);
  tree.adjacency[first].push_back({second, length});
  tree.adjacency[second].push_back({first, length});
}

NjTree source_clearcut_neighbor_joining(
    const std::vector<float>& leaf_distances,
    std::size_t leaf_count) {
  NjTree tree;
  tree.leaf_count = leaf_count;
  if (leaf_count == 0) return tree;
  // Clearcut emits a rooted binary Newick tree: n leaves plus n-1 internal
  // nodes. The last two active clades each receive the full remaining
  // distance, preserving NJ_decompose(..., NJ_LAST)'s observable quirk.
  const std::size_t node_capacity = leaf_count == 1 ? 1 : leaf_count * 2 - 1;
  tree.adjacency.resize(node_capacity);
  if (leaf_count == 1) {
    tree.root = 0;
    return tree;
  }

  std::vector<float> distances = leaf_distances;
  std::vector<std::size_t> active_nodes(leaf_count);
  std::iota(active_nodes.begin(), active_nodes.end(), 0);
  std::size_t next_node = leaf_count;
  // Build Clearcut's packed-matrix r/r2 vectors in the same pair visitation
  // order. NJ_init_r assigns r2 after completing each row, leaving the final
  // calloc-initialised r2 entry at zero.
  std::vector<float> row_sums(leaf_count, 0.0F);
  std::vector<float> transformed_sums(leaf_count, 0.0F);
  for (std::size_t first = 0; first + 1 < leaf_count; ++first) {
    for (std::size_t second = first + 1; second < leaf_count; ++second) {
      const float value = distances[first * leaf_count + second];
      row_sums[first] += value;
      row_sums[second] += value;
    }
    transformed_sums[first] = row_sums[first] /
        static_cast<float>(leaf_count - 2);
  }
  while (active_nodes.size() > 2) {
    const std::size_t count = active_nodes.size();
    std::size_t best_first = 0;
    std::size_t best_second = 1;
    float best_q = std::numeric_limits<float>::infinity();
    for (std::size_t first = 0; first + 1 < count; ++first) {
      for (std::size_t second = first + 1; second < count; ++second) {
        const float q = distances[first * count + second] -
            transformed_sums[first] - transformed_sums[second];
        // NJ_min_transform uses a strict comparison while scanning rows and
        // columns in order, so the first exactly tied pair remains selected.
        if (q < best_q) {
          best_q = q;
          best_first = first;
          best_second = second;
        }
      }
    }

    const std::size_t first_node = active_nodes[best_first];
    const std::size_t second_node = active_nodes[best_second];
    const float joined_distance = distances[best_first * count + best_second];
    const float first_length = 0.5F *
        (joined_distance + transformed_sums[best_first] -
         transformed_sums[best_second]);
    const float second_length = 0.5F *
        (joined_distance + transformed_sums[best_second] -
         transformed_sums[best_first]);
    const std::size_t joined_node = next_node++;
    add_source_edge(tree, joined_node, first_node, first_length);
    add_source_edge(tree, joined_node, second_node, second_length);

    // Mirror NJ_compute_r followed by NJ_collapse. Preserving these update and
    // accumulation orders matters because the supplied implementation retains
    // every intermediate in single precision rather than recomputing row sums.
    for (std::size_t index = best_first + 1; index < count; ++index) {
      row_sums[index] -= distances[best_first * count + index];
      if (index > best_second) {
        row_sums[index] -= distances[best_second * count + index];
      }
    }
    for (std::size_t index = 0; index < best_second; ++index) {
      if (index < best_first) {
        row_sums[index] -= distances[index * count + best_first];
      }
      row_sums[index] -= distances[index * count + best_second];
    }

    std::vector<float> mutated = distances;
    row_sums[best_first] = 0.0F;
    const float next_divisor = static_cast<float>(count - 3);
    for (std::size_t index = best_first + 1; index < count; ++index) {
      const float joined = 0.5F *
          ((distances[best_first * count + index] - first_length) +
           (distances[best_second * count + index] - second_length));
      mutated[best_first * count + index] = joined;
      mutated[index * count + best_first] = joined;
      row_sums[best_first] += joined;
      row_sums[index] += joined;
      transformed_sums[index] = row_sums[index] / next_divisor;
    }
    for (std::size_t index = 0; index < best_first; ++index) {
      const float joined = 0.5F *
          ((distances[index * count + best_first] - first_length) +
           (distances[index * count + best_second] - second_length));
      mutated[index * count + best_first] = joined;
      mutated[best_first * count + index] = joined;
      row_sums[best_first] += joined;
      row_sums[index] += joined;
      transformed_sums[index] = row_sums[index] / next_divisor;
    }
    transformed_sums[best_first] = row_sums[best_first] / next_divisor;

    // Put the joined clade at row a, copy updated row zero and its node into
    // row b, then drop row zero exactly as the packed matrix does.
    active_nodes[best_first] = joined_node;
    active_nodes[best_second] = active_nodes[0];
    for (std::size_t index = 0; index < count; ++index) {
      const float copied = mutated[index];
      mutated[best_second * count + index] = copied;
      mutated[index * count + best_second] = copied;
    }
    mutated[best_second * count + best_second] = 0.0F;
    row_sums[best_second] = row_sums[0];
    transformed_sums[best_second] = transformed_sums[0];

    std::vector<float> compacted((count - 1) * (count - 1), 0.0F);
    for (std::size_t first = 1; first < count; ++first) {
      for (std::size_t second = 1; second < count; ++second) {
        compacted[(first - 1) * (count - 1) + second - 1] =
            mutated[first * count + second];
      }
    }
    active_nodes.erase(active_nodes.begin());
    row_sums.erase(row_sums.begin());
    transformed_sums.erase(transformed_sums.begin());
    distances = std::move(compacted);
  }

  const std::size_t root = next_node++;
  const float remaining = distances[1];
  add_source_edge(tree, root, active_nodes[0], remaining);
  add_source_edge(tree, root, active_nodes[1], remaining);
  tree.root = root;
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

std::unordered_map<std::uint64_t, std::string> internal_edge_splits(
    const NjTree& tree) {
  std::unordered_map<std::uint64_t, std::string> splits;
  if (tree.leaf_count < 4) return splits;
  for (std::size_t first = 0; first < tree.adjacency.size(); ++first) {
    for (const Edge& edge : tree.adjacency[first]) {
      if (first >= edge.to || first < tree.leaf_count || edge.to < tree.leaf_count) continue;
      std::vector<std::size_t> leaves;
      collect_leaves(tree, first, edge.to, leaves);
      if (leaves.size() < 2 || leaves.size() > tree.leaf_count - 2) continue;
      const std::uint64_t edge_key =
          (static_cast<std::uint64_t>(first) << 32U) |
          static_cast<std::uint64_t>(edge.to);
      splits.emplace(
          edge_key,
          canonical_split_key(std::move(leaves), tree.leaf_count));
    }
  }
  return splits;
}

std::unordered_map<std::string, std::pair<std::size_t, std::size_t>> internal_splits(
    const NjTree& tree) {
  std::unordered_map<std::string, std::pair<std::size_t, std::size_t>> splits;
  for (const auto& [edge_key, split] : internal_edge_splits(tree)) {
    const std::size_t first = static_cast<std::size_t>(edge_key >> 32U);
    const std::size_t second = static_cast<std::size_t>(edge_key & 0xffffffffULL);
    splits.emplace(split, std::pair{first, second});
  }
  return splits;
}

std::uint64_t topology_edge_key(std::size_t first, std::size_t second) {
  return (static_cast<std::uint64_t>(std::min(first, second)) << 32U) |
      static_cast<std::uint64_t>(std::max(first, second));
}

double source_parsed_branch_length(double serialized_length) {
  // NJ_output_tree2 writes five fractional digits, but TreeToArrayP's seven-
  // character parser consumes only four. It then clamps every parsed branch
  // to [0, 1] before midpoint rooting.
  return std::min(
      1.0,
      std::trunc(std::max(0.0, serialized_length) * 10000.0) / 10000.0);
}

std::vector<double> source_node_distances(
    const NjTree& tree,
    std::size_t source) {
  std::vector<double> distances(
      tree.adjacency.size(), std::numeric_limits<double>::infinity());
  if (source >= tree.adjacency.size()) return distances;
  std::vector<std::size_t> stack{source};
  distances[source] = 0.0;
  while (!stack.empty()) {
    const std::size_t node = stack.back();
    stack.pop_back();
    for (const Edge& edge : tree.adjacency[node]) {
      if (std::isfinite(distances[edge.to])) continue;
      distances[edge.to] = distances[node] +
          source_parsed_branch_length(edge.length);
      stack.push_back(edge.to);
    }
  }
  return distances;
}

struct SourceMidpointRankTree {
  std::vector<double> ranked_distances;
  std::size_t rank_levels = 0;
  std::size_t root_node = std::numeric_limits<std::size_t>::max();
  std::size_t root_edge_first = 0;
  std::size_t root_edge_second = 0;
  double root_offset_from_first = 0.0;
  double root_edge_length = 0.0;
};

std::vector<double> source_rank_tree_distances(
    const std::vector<double>& tree_distances,
    std::size_t count,
    std::size_t& rank_levels) {
  std::vector<float> rounded(count * count, 0.0F);
  std::vector<float> levels;
  levels.reserve(count * (count - 1) / 2);
  for (std::size_t first = 0; first < count; ++first) {
    for (std::size_t second = first + 1; second < count; ++second) {
      // Tree2ArrayP rounds path lengths to five decimals; MakeTreeArrayXP2
      // then rounds them to four before assigning ascending topology ranks.
      const float input = static_cast<float>(
          tree_distances[first * count + second]);
      const float value = std::round(input * 10000.0F) / 10000.0F;
      rounded[first * count + second] = value;
      rounded[second * count + first] = value;
      levels.push_back(value);
    }
  }
  std::sort(levels.begin(), levels.end());
  levels.erase(std::unique(levels.begin(), levels.end()), levels.end());
  rank_levels = levels.size();

  std::vector<double> ranked(count * count, 0.0);
  for (std::size_t first = 0; first < count; ++first) {
    for (std::size_t second = first + 1; second < count; ++second) {
      const auto found = std::lower_bound(
          levels.begin(), levels.end(), rounded[first * count + second]);
      const double rank = static_cast<double>(found - levels.begin() + 1) / 1000.0;
      ranked[first * count + second] = rank;
      ranked[second * count + first] = rank;
    }
  }
  return ranked;
}

SourceMidpointRankTree source_midpoint_rank_tree_distances(const NjTree& tree) {
  SourceMidpointRankTree result;
  const std::size_t count = tree.leaf_count;
  if (count == 0) return result;

  std::vector<double> paths(count * count, 0.0);
  for (std::size_t source = 0; source < count; ++source) {
    const auto distances = source_node_distances(tree, source);
    for (std::size_t target = source + 1; target < count; ++target) {
      const float input = static_cast<float>(distances[target]);
      const float rounded = std::round(input * 100000.0F) / 100000.0F;
      paths[source * count + target] = rounded;
      paths[target * count + source] = rounded;
    }
  }

  std::size_t diameter_first = 0;
  std::size_t diameter_second = count > 1 ? 1 : 0;
  double diameter = 0.0;
  for (std::size_t first = 0; first + 1 < count; ++first) {
    for (std::size_t second = first + 1; second < count; ++second) {
      const double distance = paths[first * count + second];
      // TreeMidP scans in this order and retains the first exactly tied pair.
      if (diameter < distance) {
        diameter = distance;
        diameter_first = first;
        diameter_second = second;
      }
    }
  }
  const double midpoint_distance = diameter / 2.0;

  std::vector<std::size_t> parent(
      tree.adjacency.size(), tree.adjacency.size());
  std::vector<double> parent_length(tree.adjacency.size(), 0.0);
  std::vector<std::size_t> stack{diameter_first};
  parent[diameter_first] = diameter_first;
  while (!stack.empty()) {
    const std::size_t node = stack.back();
    stack.pop_back();
    for (const Edge& edge : tree.adjacency[node]) {
      if (parent[edge.to] != tree.adjacency.size()) continue;
      parent[edge.to] = node;
      parent_length[edge.to] = source_parsed_branch_length(edge.length);
      stack.push_back(edge.to);
    }
  }
  std::vector<std::size_t> path;
  for (std::size_t node = diameter_second;; node = parent[node]) {
    path.push_back(node);
    if (node == diameter_first) break;
  }
  std::reverse(path.begin(), path.end());

  result.root_node = diameter_first;
  double traversed = 0.0;
  for (std::size_t index = 0; index + 1 < path.size(); ++index) {
    const std::size_t first = path[index];
    const std::size_t second = path[index + 1];
    const double length = parent_length[second];
    if (traversed + length < midpoint_distance) {
      traversed += length;
      continue;
    }
    const double offset = std::clamp(
        midpoint_distance - traversed, 0.0, length);
    if (offset <= 1e-12) {
      result.root_node = first;
    } else if (length - offset <= 1e-12) {
      result.root_node = second;
    } else {
      result.root_node = std::numeric_limits<std::size_t>::max();
      result.root_edge_first = first;
      result.root_edge_second = second;
      result.root_offset_from_first = offset;
      result.root_edge_length = length;
    }
    break;
  }

  std::vector<double> root_distances(count, 0.0);
  if (result.root_node != std::numeric_limits<std::size_t>::max()) {
    const auto distances = source_node_distances(tree, result.root_node);
    std::copy_n(distances.begin(), count, root_distances.begin());
  } else {
    const auto from_first = source_node_distances(tree, result.root_edge_first);
    const auto from_second = source_node_distances(tree, result.root_edge_second);
    for (std::size_t leaf = 0; leaf < count; ++leaf) {
      root_distances[leaf] = std::min(
          from_first[leaf] + result.root_offset_from_first,
          from_second[leaf] +
              (result.root_edge_length - result.root_offset_from_first));
    }
  }

  std::vector<double> ultrametric = paths;
  std::vector<double> extensions(count, 0.0);
  for (std::size_t leaf = 0; leaf < count; ++leaf) {
    extensions[leaf] = midpoint_distance - root_distances[leaf];
  }
  for (std::size_t first = 0; first < count; ++first) {
    for (std::size_t second = first + 1; second < count; ++second) {
      float adjusted = static_cast<float>(paths[first * count + second]);
      adjusted += static_cast<float>(extensions[first]);
      adjusted += static_cast<float>(extensions[second]);
      adjusted = std::round(adjusted * 100000.0F) / 100000.0F;
      ultrametric[first * count + second] = adjusted;
      ultrametric[second * count + first] = adjusted;
    }
  }
  result.ranked_distances = source_rank_tree_distances(
      ultrametric, count, result.rank_levels);
  return result;
}

std::vector<double> source_collapsed_rank_tree_distances(
    const NjTree& tree,
    const SourceMidpointRankTree& raw,
    const std::unordered_set<std::string>& weak_splits,
    std::unordered_set<std::string>& effective_collapsed_splits,
    std::size_t& rank_levels) {
  const std::size_t count = tree.leaf_count;
  if (weak_splits.empty() || raw.ranked_distances.size() != count * count) {
    rank_levels = raw.rank_levels;
    return raw.ranked_distances;
  }

  struct TopologyNeighbor {
    std::size_t node = 0;
    std::uint64_t original_edge = 0;
    bool midpoint_half = false;
  };
  const bool virtual_root =
      raw.root_node == std::numeric_limits<std::size_t>::max();
  const std::size_t root = virtual_root ? tree.adjacency.size() : raw.root_node;
  std::vector<std::vector<TopologyNeighbor>> adjacency(
      tree.adjacency.size() + (virtual_root ? 1 : 0));
  const std::uint64_t midpoint_edge = virtual_root
      ? topology_edge_key(raw.root_edge_first, raw.root_edge_second)
      : std::numeric_limits<std::uint64_t>::max();
  for (std::size_t first = 0; first < tree.adjacency.size(); ++first) {
    for (const Edge& edge : tree.adjacency[first]) {
      if (first >= edge.to) continue;
      const std::uint64_t key = topology_edge_key(first, edge.to);
      if (virtual_root && key == midpoint_edge) continue;
      adjacency[first].push_back({edge.to, key, false});
      adjacency[edge.to].push_back({first, key, false});
    }
  }
  if (virtual_root) {
    adjacency[root].push_back({raw.root_edge_first, midpoint_edge, true});
    adjacency[root].push_back({raw.root_edge_second, midpoint_edge, true});
    adjacency[raw.root_edge_first].push_back({root, midpoint_edge, true});
    adjacency[raw.root_edge_second].push_back({root, midpoint_edge, true});
  }

  std::vector<std::size_t> parent(adjacency.size(), adjacency.size());
  std::vector<std::uint64_t> parent_edge(
      adjacency.size(), std::numeric_limits<std::uint64_t>::max());
  std::vector<bool> parent_midpoint_half(adjacency.size(), false);
  std::vector<std::size_t> depth(adjacency.size(), 0);
  std::vector<std::size_t> order{root};
  parent[root] = root;
  for (std::size_t cursor = 0; cursor < order.size(); ++cursor) {
    const std::size_t node = order[cursor];
    for (const TopologyNeighbor& edge : adjacency[node]) {
      if (parent[edge.node] != adjacency.size()) continue;
      parent[edge.node] = node;
      parent_edge[edge.node] = edge.original_edge;
      parent_midpoint_half[edge.node] = edge.midpoint_half;
      depth[edge.node] = depth[node] + 1;
      order.push_back(edge.node);
    }
  }

  std::vector<std::size_t> descendant_leaf(
      adjacency.size(), std::numeric_limits<std::size_t>::max());
  for (auto cursor = order.rbegin(); cursor != order.rend(); ++cursor) {
    const std::size_t node = *cursor;
    if (node < count) descendant_leaf[node] = node;
    for (const TopologyNeighbor& edge : adjacency[node]) {
      if (parent[edge.node] != node) continue;
      descendant_leaf[node] = std::min(
          descendant_leaf[node], descendant_leaf[edge.node]);
    }
  }
  std::vector<double> node_rank(adjacency.size(), 0.0);
  for (const std::size_t node : order) {
    std::size_t first_leaf = std::numeric_limits<std::size_t>::max();
    std::size_t second_leaf = std::numeric_limits<std::size_t>::max();
    for (const TopologyNeighbor& edge : adjacency[node]) {
      if (parent[edge.node] != node ||
          descendant_leaf[edge.node] == std::numeric_limits<std::size_t>::max()) {
        continue;
      }
      if (first_leaf == std::numeric_limits<std::size_t>::max()) {
        first_leaf = descendant_leaf[edge.node];
      } else {
        second_leaf = descendant_leaf[edge.node];
        break;
      }
    }
    if (second_leaf != std::numeric_limits<std::size_t>::max()) {
      node_rank[node] = raw.ranked_distances[first_leaf * count + second_leaf];
    }
  }

  const auto edge_splits = internal_edge_splits(tree);
  std::unordered_map<std::string, std::size_t> split_occurrences;
  for (const auto& [key, split] : edge_splits) {
    (void)key;
    ++split_occurrences[split];
  }
  const auto collapses_parent_edge = [&](std::size_t node) {
    if (node == root || parent_midpoint_half[node]) return false;
    const auto split = edge_splits.find(parent_edge[node]);
    if (split == edge_splits.end() ||
        split_occurrences[split->second] != 1 ||
        !weak_splits.contains(split->second)) {
      return false;
    }
    effective_collapsed_splits.insert(split->second);
    return true;
  };

  const auto lca = [&](std::size_t first, std::size_t second) {
    while (depth[first] > depth[second]) first = parent[first];
    while (depth[second] > depth[first]) second = parent[second];
    while (first != second) {
      first = parent[first];
      second = parent[second];
    }
    return first;
  };

  std::vector<double> collapsed(count * count, 0.0);
  for (std::size_t first = 0; first < count; ++first) {
    for (std::size_t second = first + 1; second < count; ++second) {
      std::size_t node = lca(first, second);
      while (node != root && collapses_parent_edge(node)) node = parent[node];
      while (node != root && node_rank[node] == 0.0) node = parent[node];
      const double rank = node_rank[node] > 0.0
          ? node_rank[node]
          : raw.ranked_distances[first * count + second];
      collapsed[first * count + second] = rank;
      collapsed[second * count + first] = rank;
    }
  }
  std::vector<double> levels;
  for (std::size_t first = 0; first < count; ++first) {
    for (std::size_t second = first + 1; second < count; ++second) {
      levels.push_back(collapsed[first * count + second]);
    }
  }
  std::sort(levels.begin(), levels.end());
  levels.erase(std::unique(levels.begin(), levels.end()), levels.end());
  rank_levels = levels.size();
  return collapsed;
}

std::size_t source_vb_round_nonnegative(double value) {
  if (!(value > 0.0)) return 0;
  const double floor_value = std::floor(value);
  const auto lower = static_cast<std::size_t>(floor_value);
  const double fraction = value - floor_value;
  if (fraction > 0.5 || (fraction == 0.5 && (lower & 1U) != 0)) {
    return lower + 1;
  }
  return lower;
}

}  // namespace

std::vector<std::uint32_t> source_event_bootstrap_weights(
    std::size_t site_count,
    std::size_t bootstrap_replicates,
    std::uint32_t seed) {
  return source_event_bootstrap_weights_impl(
      site_count, bootstrap_replicates, seed);
}

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
  if (raw_tree_distances.size() != count * count ||
      collapsed_tree_distances.size() != count * count) {
    return kSaturatedDistance;
  }
  return collapsed ? collapsed_tree_distances[index] : raw_tree_distances[index];
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
  evidence.bootstrap_random_seed = static_cast<std::uint32_t>(seed);
  evidence.source_clearcut_float_nj = true;
  evidence.source_ranked_tree_distances = true;
  evidence.source_midpoint_ultrametric = true;
  evidence.source_parent_rank_collapse = true;
  evidence.source_seqboot2_bootstrap = true;
  evidence.source_bootstrap_pseudocount = true;
  if (sequences.size() < 3 || one_based_positions.empty()) return evidence;

  const std::size_t stride = bootstrap_replicates + 1;
  const auto weights = source_event_bootstrap_weights_impl(
      one_based_positions.size(),
      bootstrap_replicates,
      evidence.bootstrap_random_seed);
  const auto base_distances = weighted_distance_matrix(
      alignment,
      sequences,
      one_based_positions,
      weights,
      stride,
      0);
  evidence.jukes_cantor.assign(base_distances.begin(), base_distances.end());
  const NjTree base_tree = source_clearcut_neighbor_joining(
      base_distances, sequences.size());
  evidence.negative_branches_normalized =
      base_tree.negative_branches_normalized;
  const auto base_splits = internal_splits(base_tree);
  evidence.internal_branches = base_splits.size();

  std::unordered_map<std::string, std::size_t> support;
  support.reserve(base_splits.size());
  for (const auto& [key, endpoints] : base_splits) {
    (void)endpoints;
    support.emplace(key, 0);
  }

  for (std::size_t replicate = 1;
       replicate <= bootstrap_replicates;
       ++replicate) {
    const auto replicate_distances = weighted_distance_matrix(
        alignment,
        sequences,
        one_based_positions,
        weights,
        stride,
        replicate);
    const NjTree replicate_tree = source_clearcut_neighbor_joining(
        replicate_distances, sequences.size());
    const auto replicate_splits = internal_splits(replicate_tree);
    for (auto& [key, count] : support) {
      if (replicate_splits.contains(key)) ++count;
    }
  }
  evidence.bootstrap_replicates = bootstrap_replicates;

  std::unordered_set<std::string> weak_splits;
  std::unordered_map<std::string, double> split_support;
  for (const auto& [key, count] : support) {
    // TreeRepsP starts every base-tree split with its retained replicate-zero
    // observation, divides by reps + 1, converts to an integer percentage
    // with VB6's CLng rounding, then CollapseNodes compares that integer with
    // the 50% cutoff.
    const std::size_t percent = bootstrap_replicates == 0
        ? 100
        : source_vb_round_nonnegative(
              (static_cast<double>(count) + 1.0) /
              (static_cast<double>(bootstrap_replicates) + 1.0) * 100.0);
    const double fraction = static_cast<double>(percent) / 100.0;
    split_support.emplace(key, fraction);
    if (percent < static_cast<std::size_t>(kBootstrapCollapseCutoff * 100.0)) {
      weak_splits.insert(key);
    } else {
      ++evidence.supported_internal_branches;
    }
  }
  const SourceMidpointRankTree raw_ranked =
      source_midpoint_rank_tree_distances(base_tree);
  evidence.raw_tree_distances = raw_ranked.ranked_distances;
  evidence.raw_distance_rank_levels = raw_ranked.rank_levels;
  std::unordered_set<std::string> effective_collapsed_splits;
  evidence.collapsed_tree_distances = source_collapsed_rank_tree_distances(
      base_tree,
      raw_ranked,
      weak_splits,
      effective_collapsed_splits,
      evidence.collapsed_distance_rank_levels);
  evidence.topology_node_count = base_tree.adjacency.size();
  evidence.topology_root = base_tree.root;
  const auto base_edge_splits = internal_edge_splits(base_tree);
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
      const auto edge_split = base_edge_splits.find(edge_key);
      const auto support_value = edge_split == base_edge_splits.end()
          ? split_support.end()
          : split_support.find(edge_split->second);
      const bool edge_collapsed = internal &&
          edge_split != base_edge_splits.end() &&
          effective_collapsed_splits.contains(edge_split->second);
      evidence.topology_edges.push_back({
          static_cast<std::uint32_t>(first),
          static_cast<std::uint32_t>(edge.to),
          edge.length,
          support_value == split_support.end() ? 1.0 : support_value->second,
          internal,
          edge_collapsed,
      });
    }
  }
  evidence.usable = true;
  return evidence;
}

}  // namespace rdp
