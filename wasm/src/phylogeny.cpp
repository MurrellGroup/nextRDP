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
  // The installed 32-bit FinishDists path retains these intermediates in its
  // x87 register stack despite the source-level float temporaries.  Computing
  // in double and rounding only the stored matrix entry reproduces every one
  // of the 600 captured Dataset0 outside/inside distances exactly.  Forcing
  // source-level float rounding changes many entries by one to three ULPs and
  // can change Clearcut's choice between near-tied joins.
  const double valid = static_cast<double>(comparable);
  const double identity =
      (valid - static_cast<double>(differences)) / valid;
  if (identity == 1.0) return 0.0F;
  if (!(identity > 0.25)) return static_cast<float>(kSaturatedDistance);
  const double argument = (4.0 * identity - 1.0) / 3.0;
  if (!(argument > 0.0)) return static_cast<float>(kSaturatedDistance);
  const float corrected = static_cast<float>(-0.75 * std::log(argument));
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
  const std::size_t node_capacity = leaf_count == 1 ? 1 : leaf_count * 2 - 1;
  tree.adjacency.resize(node_capacity);
  if (leaf_count == 1) {
    tree.root = 0;
    return tree;
  }

  // Literal packed-matrix Clearcut NJ used by DNA5.dll. Keeping the moving
  // value/r/r2 bases and mutation order matters even where a recomputed
  // full-matrix implementation produces the same unrooted splits: the final
  // two active clades determine Newick rooting, and Tree2ArrayP2 observes it.
  const auto packed_cells = [](std::size_t size) {
    return size * (size + 1) / 2;
  };
  const auto packed_index = [](std::size_t first,
                               std::size_t second,
                               std::size_t size) {
    return first * (2 * size - first - 1) / 2 + second;
  };
  std::vector<float> distances(packed_cells(leaf_count), 0.0F);
  for (std::size_t first = 0; first + 1 < leaf_count; ++first) {
    for (std::size_t second = first + 1; second < leaf_count; ++second) {
      distances[packed_index(first, second, leaf_count)] =
          leaf_distances[first + second * leaf_count];
    }
  }
  std::vector<float> row_sums(leaf_count, 0.0F);
  std::vector<float> transformed_sums(leaf_count, 0.0F);
  std::vector<std::size_t> active_nodes(leaf_count);
  std::iota(active_nodes.begin(), active_nodes.end(), 0);
  std::size_t value_base = 0;
  std::size_t sum_base = 0;
  std::size_t size = leaf_count;
  std::size_t next_node = leaf_count;

  std::size_t scan = 0;
  for (std::size_t first = 0; first + 1 < size; ++first) {
    ++scan;
    for (std::size_t second = first + 1; second < size; ++second) {
      const float value = distances[value_base + scan++];
      row_sums[sum_base + first] += value;
      row_sums[sum_base + second] += value;
    }
    transformed_sums[sum_base + first] =
        row_sums[sum_base + first] / static_cast<float>(size - 2);
  }

  while (size > 2) {
    std::size_t best_first = 0;
    std::size_t best_second = 1;
    float best = std::numeric_limits<float>::infinity();
    std::size_t position = 0;
    for (std::size_t first = 0; first < size; ++first) {
      ++position;
      for (std::size_t second = first + 1; second < size; ++second) {
        const float transformed =
            distances[value_base + position++] -
            (transformed_sums[sum_base + first] +
             transformed_sums[sum_base + second]);
        if (transformed < best) {
          best = transformed;
          best_first = first;
          best_second = second;
        }
      }
    }

    const auto at = [&](std::size_t first, std::size_t second) -> float& {
      if (second < first) std::swap(first, second);
      return distances[
          value_base + packed_index(first, second, size)];
    };
    const float joined_distance = at(best_first, best_second);
    const float first_length = joined_distance * 0.5F +
        (transformed_sums[sum_base + best_first] -
         transformed_sums[sum_base + best_second]) * 0.5F;
    const float second_length = joined_distance * 0.5F +
        (transformed_sums[sum_base + best_second] -
         transformed_sums[sum_base + best_first]) * 0.5F;
    const std::size_t joined_node = next_node++;
    add_source_edge(
        tree, joined_node, active_nodes[best_first], first_length);
    add_source_edge(
        tree, joined_node, active_nodes[best_second], second_length);
    active_nodes[best_first] = joined_node;
    active_nodes[best_second] = active_nodes[0];

    for (std::size_t index = best_first + 1; index < size; ++index) {
      row_sums[sum_base + index] -= at(best_first, index);
      if (index > best_second) {
        row_sums[sum_base + index] -= at(best_second, index);
      }
    }
    for (std::size_t index = 0; index < best_second; ++index) {
      if (index < best_first) {
        row_sums[sum_base + index] -= at(index, best_first);
      }
      row_sums[sum_base + index] -= at(index, best_second);
    }

    row_sums[sum_base + best_first] = 0.0F;
    for (std::size_t index = best_first + 1; index < size; ++index) {
      const float joined =
          ((at(best_first, index) - first_length) +
           (at(best_second, index) - second_length)) * 0.5F;
      at(best_first, index) = joined;
      row_sums[sum_base + best_first] += joined;
      row_sums[sum_base + index] += joined;
      transformed_sums[sum_base + index] =
          row_sums[sum_base + index] / static_cast<float>(size - 3);
    }
    for (std::size_t index = 0; index < best_first; ++index) {
      const float joined =
          ((at(index, best_first) - first_length) +
           (at(index, best_second) - second_length)) * 0.5F;
      at(index, best_first) = joined;
      row_sums[sum_base + best_first] += joined;
      row_sums[sum_base + index] += joined;
      transformed_sums[sum_base + index] =
          row_sums[sum_base + index] / static_cast<float>(size - 3);
    }
    transformed_sums[sum_base + best_first] =
        row_sums[sum_base + best_first] / static_cast<float>(size - 3);

    for (std::size_t index = 0; index < best_second; ++index) {
      at(index, best_second) = distances[value_base + index];
    }
    std::size_t source = best_second + 1;
    for (std::size_t index = best_second + 1; index < size; ++index) {
      at(best_second, index) = distances[value_base + source++];
    }
    row_sums[sum_base + best_second] = row_sums[sum_base];
    transformed_sums[sum_base + best_second] = transformed_sums[sum_base];

    active_nodes.erase(active_nodes.begin());
    value_base += size;
    ++sum_base;
    --size;
  }

  const float remaining =
      distances[value_base + packed_index(0, 1, size)];
  const std::size_t root = next_node++;
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

void append_source_branch_length(std::string& output, double branch_length) {
  output.push_back(':');
  double value = std::abs(branch_length);
  if (value < 1.0) {
    output += "0.";
  } else {
    output.push_back(static_cast<char>(static_cast<int>(value) + '0'));
    output.push_back('.');
  }
  value -= static_cast<int>(value);
  int modulus = 100000;
  int remainder = static_cast<int>(value * static_cast<double>(modulus));
  int digit = remainder / (modulus / 10);
  while (modulus > 10) {
    output.push_back(static_cast<char>(digit + '0'));
    modulus /= 10;
    remainder -= digit * modulus;
    digit = remainder / (modulus / 10);
  }
}

void append_source_taxon_name(
    std::string& output,
    std::size_t taxon,
    int name_modulus) {
  output.push_back('S');
  int remainder = static_cast<int>(taxon);
  int modulus = name_modulus;
  int digit = remainder / modulus;
  output.push_back(static_cast<char>(digit + '0'));
  while (modulus > 1) {
    remainder -= digit * modulus;
    digit = remainder / (modulus / 10);
    output.push_back(static_cast<char>(digit + '0'));
    modulus /= 10;
  }
}

void append_source_clearcut_newick_node(
    const NjTree& tree,
    std::size_t node,
    std::size_t parent,
    std::size_t root_left,
    double parent_length,
    int name_modulus,
    std::string& output) {
  if (node < tree.leaf_count) {
    append_source_taxon_name(output, node, name_modulus);
    append_source_branch_length(output, parent_length);
    return;
  }
  output.push_back('(');
  bool first_child = true;
  for (const Edge& edge : tree.adjacency[node]) {
    if (edge.to == parent) continue;
    if (!first_child) output.push_back(',');
    append_source_clearcut_newick_node(
        tree, edge.to, node, root_left, edge.length, name_modulus, output);
    first_child = false;
  }
  output.push_back(')');
  // NJ_output_tree2 omits the length of the root and, unusually, the root's
  // left internal child. TreeToArrayP2 then scans forward to the next decimal
  // for that child. Preserve that parser-visible omission exactly.
  if (node != tree.root && node != root_left) {
    append_source_branch_length(output, parent_length);
  }
}

std::vector<std::uint8_t> source_clearcut_newick_holder(const NjTree& tree) {
  if (tree.adjacency.empty()) return {};
  int name_modulus = 10;
  if (tree.leaf_count > 100) name_modulus = 100;
  if (tree.leaf_count > 1000) name_modulus = 1000;
  const std::size_t root_left = tree.adjacency[tree.root].empty()
      ? tree.root
      : tree.adjacency[tree.root].front().to;
  std::string serialized;
  serialized.reserve(tree.leaf_count * 24);
  append_source_clearcut_newick_node(
      tree,
      tree.root,
      tree.adjacency.size(),
      root_left,
      0.0,
      name_modulus,
      serialized);
  serialized.push_back(';');
  std::vector<std::uint8_t> holder(serialized.size() + 1, 0);
  std::copy(serialized.begin(), serialized.end(), holder.begin() + 1);
  return holder;
}

std::vector<double> source_tree2arrayp2_rank_distances(
    const NjTree& tree,
    std::size_t& rank_levels) {
  const std::size_t count = tree.leaf_count;
  rank_levels = 0;
  if (count == 0) return {};
  const int upper = static_cast<int>(count) - 1;
  const int name_length = std::max(
      2, static_cast<int>(std::to_string(upper).size()));
  const int max_position = upper * 3 + 100;
  const auto holder = source_clearcut_newick_holder(tree);
  if (holder.empty()) return std::vector<double>(count * count, 0.0);
  const auto semicolon = std::find(holder.begin(), holder.end(), ';');
  const int tree_length = static_cast<int>(semicolon - holder.begin());
  std::vector<int> node_order(max_position + 1, 0);
  std::vector<int> done_node(max_position + 1, 0);
  std::vector<double> node_length(max_position + 1, 0.0);
  std::vector<double> num_done(max_position + 1, 0.0);
  std::vector<float> tree_matrix(count * count, 0.0F);

  const auto parse_length = [&](int decimal_position) {
    if (decimal_position < 2 ||
        decimal_position >= static_cast<int>(holder.size())) {
      return 0.0;
    }
    int adjustment = 2;
    double value = 0.0;
    for (int offset = 0; offset <= 6; ++offset) {
      const int position = decimal_position - 2 + offset;
      const int character = position >= 0 &&
              position < static_cast<int>(holder.size())
          ? holder[position]
          : 0;
      if (character > '0' - 1 && character < '9' + 1) {
        value += static_cast<double>(character - '0') *
            std::pow(10.0, 6 - offset - adjustment);
      } else {
        --adjustment;
      }
    }
    return holder[decimal_position - 2] == '-'
        ? 0.0
        : value / 10000.0;
  };

  int last_position = 0;
  int current_position = 0;
  int current_node = upper;
  while (last_position < tree_length && current_position <= max_position) {
    ++last_position;
    if (holder[last_position] == 'S' && current_position <= max_position) {
      for (int offset = 1; offset <= name_length; ++offset) {
        node_order[current_position] += static_cast<int>(
            0.1 + (holder[last_position + offset] - '0') *
                std::pow(10.0, name_length - offset));
      }
      int tree_position = last_position + 2;
      while (tree_position < tree_length) {
        if (holder[tree_position] == '.') {
          const int node = node_order[current_position];
          node_length[node] += parse_length(tree_position);
          done_node[node] = 1;
          break;
        }
        ++tree_position;
      }
      while (current_position <= max_position) {
        ++current_position;
        if (node_order[current_position] == 0) break;
      }
    } else if (holder[last_position] == '(' && current_position <= max_position) {
      ++current_node;
      node_order[current_position] = current_node;
      while (node_order[current_position] != 0) {
        ++current_position;
        if (current_position > max_position) break;
      }
      if (current_node != upper + 1) {
        int tree_node = current_node;
        int temporary_position = current_position;
        int tree_position = last_position;
        do {
          ++tree_position;
          if (holder[tree_position] == '(') {
            ++tree_node;
            ++temporary_position;
            if (temporary_position > max_position) break;
          } else if (holder[tree_position] == 'S') {
            tree_position += name_length + 7;
            ++temporary_position;
            if (temporary_position > max_position) break;
          } else if (holder[tree_position] == ')') {
            ++temporary_position;
            if (temporary_position > max_position) break;
            --tree_node;
            if (tree_node == current_node - 1) {
              while (tree_position < tree_length) {
                ++tree_position;
                if (holder[tree_position] == '.') break;
              }
              node_length[current_node] += parse_length(tree_position);
              done_node[current_node] = 1;
              node_order[temporary_position - 1] = current_node;
            }
          } else if (holder[tree_position] == ';' ||
                     holder[tree_position] == 0) {
            break;
          }
        } while (done_node[current_node] == 0);
      }
    } else if (holder[last_position] == ';' || holder[last_position] == 0) {
      break;
    }
  }

  for (double& length : node_length) {
    length = std::clamp(length, 0.0, 1.0);
  }
  int saw_empty = 0;
  for (int& node : node_order) {
    if (node != 0) continue;
    if (saw_empty == 0) {
      saw_empty = 1;
    } else {
      node = upper + 1;
      break;
    }
  }
  for (int first_position = 0;
       first_position <= max_position;
       ++first_position) {
    std::fill(num_done.begin(), num_done.end(), 1.0);
    if (node_order[first_position] > upper) continue;
    double distance = node_length[node_order[first_position]];
    for (int second_position = first_position + 1;
         second_position <= max_position;
         ++second_position) {
      if (node_order[second_position] == upper + 1) break;
      if (node_order[second_position] > upper) {
        distance += node_length[node_order[second_position]] *
            num_done[node_order[second_position]];
        num_done[node_order[second_position]] *= -1.0;
      } else {
        const float value = static_cast<float>(
            distance + node_length[node_order[second_position]]);
        tree_matrix[
            node_order[first_position] + node_order[second_position] * count] =
            value;
        tree_matrix[
            node_order[second_position] + node_order[first_position] * count] =
            value;
      }
    }
  }
  for (std::size_t first = 0; first < count; ++first) {
    for (std::size_t second = first + 1; second < count; ++second) {
      const float value = static_cast<float>(
          std::round(tree_matrix[first + second * count] * 100000.0) /
          100000.0);
      tree_matrix[first + second * count] = value;
      tree_matrix[second + first * count] = value;
    }
  }
  for (double& length : node_length) {
    length = std::round(length * 100000.0) / 100000.0;
  }

  double diameter = 0.0;
  int diameter_first = 0;
  int diameter_second = 0;
  for (int first = 0; first < upper; ++first) {
    for (int second = first + 1; second <= upper; ++second) {
      if (diameter < tree_matrix[first + second * count]) {
        diameter = tree_matrix[first + second * count];
        diameter_first = first;
        diameter_second = second;
      }
    }
  }
  const double midpoint = diameter / 2.0;
  std::fill(num_done.begin(), num_done.end(), 1.0);
  int midpoint_position = 0;
  double traversed = 0.0;
  for (int position = 0; position <= max_position; ++position) {
    int other = -1;
    if (node_order[position] == diameter_first) other = diameter_second;
    else if (node_order[position] == diameter_second) other = diameter_first;
    if (other < 0) continue;
    int end_position = position + 1;
    do {
      if (node_order[end_position] > upper) {
        num_done[node_order[end_position]] *= -1.0;
      } else if (node_order[end_position] == other) {
        break;
      }
      ++end_position;
    } while (node_order[end_position] != other);

    traversed += node_length[node_order[position]];
    if (traversed < midpoint) {
      for (int path_position = position + 1;
           path_position <= end_position;
           ++path_position) {
        if (node_order[path_position] > upper ||
            node_order[path_position] == node_order[end_position]) {
          if (num_done[node_order[path_position]] == -1.0 ||
              node_order[path_position] == node_order[end_position]) {
            num_done[node_order[path_position]] *= -1.0;
            if (traversed + node_length[node_order[path_position]] < midpoint) {
              traversed += node_length[node_order[path_position]];
            } else {
              midpoint_position = path_position;
              break;
            }
          }
        }
      }
    } else {
      midpoint_position = position;
    }
    break;
  }

  int maximum_nonzero = max_position;
  for (int position = max_position; position >= 0; --position) {
    if (node_order[position] != 0) {
      maximum_nonzero = position;
      break;
    }
  }
  int left_root = 0;
  for (; left_root <= maximum_nonzero; ++left_root) {
    if (node_order[left_root] == node_order[midpoint_position]) break;
  }
  int right_root = left_root;
  if (node_order[midpoint_position] > upper) {
    for (int position = maximum_nonzero; position >= 0; --position) {
      if (node_order[position] == node_order[midpoint_position]) {
        right_root = position;
        break;
      }
    }
  }

  std::vector<int> reordered;
  reordered.reserve(maximum_nonzero + 5);
  reordered.push_back(upper * 2 + 3);
  for (int position = left_root; position <= right_root; ++position) {
    reordered.push_back(node_order[position]);
  }
  reordered.push_back(upper * 2 + 2);
  for (int position = left_root - 1; position > 0; --position) {
    reordered.push_back(node_order[position]);
  }
  for (int position = maximum_nonzero - 1;
       position >= right_root + 1;
       --position) {
    reordered.push_back(node_order[position]);
  }
  reordered.push_back(upper * 2 + 2);
  reordered.push_back(upper * 2 + 3);

  std::vector<int> done(max_position + 1, 0);
  std::vector<double> result(count * count, 0.0);
  int rank = upper + 1;
  for (int outer_close = static_cast<int>(reordered.size()) - 1;
       outer_close > 0;
       --outer_close) {
    if (reordered[outer_close] <= upper ||
        done[reordered[outer_close]] != 0) {
      continue;
    }
    done[reordered[outer_close]] = 1;
    const int outer_label = reordered[outer_close];
    int inner_close = -1;
    int inner_label = -1;
    for (int position = outer_close - 1; position > 0; --position) {
      if (reordered[position] > upper && done[reordered[position]] == 0) {
        inner_label = reordered[position];
        inner_close = position;
        break;
      }
    }
    int inner_open = -1;
    for (int position = inner_close - 1; position > 0; --position) {
      if (reordered[position] == inner_label) {
        inner_open = position;
        break;
      }
    }
    if (inner_open < 0) continue;
    int outer_open = -1;
    for (int position = inner_open - 1; position >= 0; --position) {
      if (reordered[position] == outer_label) {
        outer_open = position;
        break;
      }
    }
    if (outer_open < 0) continue;
    std::vector<int> first_list;
    std::vector<int> second_list;
    for (int position = outer_open; position < inner_open; ++position) {
      if (reordered[position] <= upper) first_list.push_back(reordered[position]);
    }
    for (int position = inner_close + 1;
         position < outer_close;
         ++position) {
      if (reordered[position] <= upper) first_list.push_back(reordered[position]);
    }
    for (int position = inner_open + 1;
         position < inner_close;
         ++position) {
      if (reordered[position] <= upper) second_list.push_back(reordered[position]);
    }
    if (first_list.empty() || second_list.empty()) continue;
    --rank;
    for (const int first : first_list) {
      for (const int second : second_list) {
        result[first * count + second] = rank / 1000.0;
        result[second * count + first] = rank / 1000.0;
      }
    }
  }
  for (std::size_t position = 0; position + 1 < reordered.size(); ++position) {
    if (reordered[position] <= upper && reordered[position + 1] <= upper) {
      --rank;
      result[reordered[position] * count + reordered[position + 1]] =
          rank / 1000.0;
      result[reordered[position + 1] * count + reordered[position]] =
          rank / 1000.0;
    }
  }
  std::vector<double> levels;
  levels.reserve(count > 0 ? count - 1 : 0);
  for (std::size_t first = 0; first < count; ++first) {
    for (std::size_t second = first + 1; second < count; ++second) {
      levels.push_back(result[first * count + second]);
    }
  }
  std::sort(levels.begin(), levels.end());
  levels.erase(std::unique(levels.begin(), levels.end()), levels.end());
  rank_levels = levels.size();
  return result;
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
  // RDP 5.93 uses Tree2ArrayP2 here. Unlike the older Tree2ArrayP path, it
  // assigns ranks from a midpoint-rerooted Newick token stream rather than
  // ranking branch-length-derived ultrametric distances. Its parser quirks
  // are observable in role consensus, so retain them in the conversion.
  result.ranked_distances =
      source_tree2arrayp2_rank_distances(tree, result.rank_levels);
  return result;
}

struct SourceNodeGroups {
  std::vector<float> trace_values;
  std::vector<std::vector<std::uint8_t>> members;
  std::size_t populated = 0;
};

SourceNodeGroups source_midpoint_node_groups(
    const SourceMidpointRankTree& raw,
    std::size_t count) {
  SourceNodeGroups groups;
  groups.trace_values.assign(count, 0.0F);
  groups.members.assign(count, std::vector<std::uint8_t>(count, 0));
  if (raw.ranked_distances.size() != count * count || count == 0) {
    return groups;
  }

  // MakeNodeDepth first converts each matrix value back to its integer
  // thousandth, then scans those levels in ascending order.  A node's match
  // row is the union of all leaves participating in a pair at that level.
  const int maximum_level = static_cast<int>((count - 1) * 2);
  std::vector<std::vector<std::uint8_t>> by_level(
      static_cast<std::size_t>(maximum_level + 1),
      std::vector<std::uint8_t>(count, 0));
  std::vector<std::uint8_t> present(
      static_cast<std::size_t>(maximum_level + 1), 0);
  for (std::size_t first = 0; first + 1 < count; ++first) {
    for (std::size_t second = first + 1; second < count; ++second) {
      const int level = static_cast<int>(std::nearbyint(
          raw.ranked_distances[first * count + second] * 1000.0));
      if (level < 0 || level > maximum_level) continue;
      present[static_cast<std::size_t>(level)] = 1;
      by_level[static_cast<std::size_t>(level)][first] = 1;
      by_level[static_cast<std::size_t>(level)][second] = 1;
    }
  }
  for (int level = 0;
       level <= maximum_level && groups.populated < count;
       ++level) {
    if (present[static_cast<std::size_t>(level)] == 0) continue;
    groups.trace_values[groups.populated] =
        static_cast<float>(level) / 1000.0F;
    groups.members[groups.populated] =
        std::move(by_level[static_cast<std::size_t>(level)]);
    ++groups.populated;
  }
  return groups;
}

std::vector<std::vector<std::uint8_t>> source_newick_node_groups(
    const NjTree& tree) {
  const std::size_t count = tree.leaf_count;
  std::vector<std::vector<std::uint8_t>> groups(
      count, std::vector<std::uint8_t>(count, 0));
  if (count == 0) return groups;
  const int upper = static_cast<int>(count) - 1;
  const int name_length = std::max(
      2, static_cast<int>(std::to_string(upper).size()));
  const auto holder = source_clearcut_newick_holder(tree);
  std::size_t node_count = 0;
  for (std::size_t position = 1;
       position < holder.size() && node_count < count;
       ++position) {
    if (holder[position] != '(') continue;
    int depth = 1;
    for (std::size_t cursor = position + 1;
         cursor < holder.size() && depth > 0;
         ++cursor) {
      if (holder[cursor] == '(') {
        ++depth;
      } else if (holder[cursor] == ')') {
        --depth;
      } else if (holder[cursor] == 'S') {
        int taxon = 0;
        for (int digit = 1; digit <= name_length; ++digit) {
          taxon += static_cast<int>(
              0.1 + (holder[cursor + static_cast<std::size_t>(digit)] - '0') *
                  std::pow(10.0, name_length - digit));
        }
        if (taxon >= 0 && taxon <= upper) {
          groups[node_count][static_cast<std::size_t>(taxon)] = 1;
        }
      }
    }
    ++node_count;
  }
  return groups;
}

void source_add_tree_group_support(
    const NjTree& replicate_tree,
    const SourceNodeGroups& base_groups,
    std::vector<std::size_t>& support) {
  const std::size_t count = replicate_tree.leaf_count;
  const auto replicate_groups = source_newick_node_groups(replicate_tree);
  std::vector<std::uint8_t> done(count, 0);
  // This is TreeGroupsXP's loop order, including its all-zero unused node row
  // and its one-use rule for each replicate clade.  Both the exact group and
  // its complement count as a match.
  for (std::size_t base = 0; base < count; ++base) {
    for (std::size_t candidate = 0; candidate < count; ++candidate) {
      if (done[candidate] != 0) continue;
      std::size_t misses = 0;
      std::size_t hits = 0;
      for (std::size_t leaf = 0; leaf < count; ++leaf) {
        if (replicate_groups[candidate][leaf] ==
            base_groups.members[base][leaf]) {
          ++hits;
        } else {
          ++misses;
        }
        if (misses > 0 && hits > 0) break;
      }
      if (misses == 0 || hits == 0) {
        ++support[base];
        done[candidate] = 1;
      }
    }
  }
}

std::vector<double> source_collapsed_rank_tree_distances(
    const SourceMidpointRankTree& raw,
    const SourceNodeGroups& groups,
    const std::vector<std::size_t>& support_percent,
    std::size_t count,
    std::size_t& rank_levels) {
  if (raw.ranked_distances.size() != count * count) {
    rank_levels = 0;
    return {};
  }
  std::vector<float> working(raw.ranked_distances.begin(),
                             raw.ranked_distances.end());
  std::vector<float> collapsed(count * count, 0.0F);
  for (std::size_t node = 0; node < count; ++node) {
    const float trace = groups.trace_values[node];
    const bool weak = node >= support_percent.size() ||
        support_percent[node] <
            static_cast<std::size_t>(kBootstrapCollapseCutoff * 100.0);
    if (weak) {
      bool found_pair = false;
      float next_distance = 100000.0F;
      for (std::size_t first = 0; first + 1 < count && !found_pair; ++first) {
        for (std::size_t second = first + 1; second < count; ++second) {
          if (working[first * count + second] != trace) continue;
          found_pair = true;
          for (std::size_t other = 0; other < count; ++other) {
            const float first_distance = working[other * count + first];
            const float second_distance = working[other * count + second];
            if (first_distance == second_distance &&
                first_distance > trace && first_distance < next_distance) {
              next_distance = first_distance;
            }
          }
          break;
        }
      }
      if (!found_pair) continue;
      if (next_distance < 100000.0F) {
        for (std::size_t first = 0; first + 1 < count; ++first) {
          for (std::size_t second = first + 1; second < count; ++second) {
            if (working[first * count + second] != trace) continue;
            working[first * count + second] = next_distance;
            working[second * count + first] = next_distance;
          }
        }
        continue;
      }
    }
    for (std::size_t first = 0; first + 1 < count; ++first) {
      for (std::size_t second = first + 1; second < count; ++second) {
        if (working[first * count + second] != trace) continue;
        collapsed[first * count + second] = trace;
        collapsed[second * count + first] = trace;
      }
    }
  }

  std::vector<double> result(collapsed.begin(), collapsed.end());
  std::vector<double> levels;
  for (std::size_t first = 0; first + 1 < count; ++first) {
    for (std::size_t second = first + 1; second < count; ++second) {
      levels.push_back(result[first * count + second]);
    }
  }
  std::sort(levels.begin(), levels.end());
  levels.erase(std::unique(levels.begin(), levels.end()), levels.end());
  rank_levels = levels.size();
  return result;
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
  const SourceMidpointRankTree raw_ranked =
      source_midpoint_rank_tree_distances(base_tree);
  const SourceNodeGroups base_node_groups =
      source_midpoint_node_groups(raw_ranked, sequences.size());
  std::vector<std::size_t> node_group_support(sequences.size(), 0);
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
    source_add_tree_group_support(
        replicate_tree, base_node_groups, node_group_support);
    const auto replicate_splits = internal_splits(replicate_tree);
    for (auto& [key, count] : support) {
      if (replicate_splits.contains(key)) ++count;
    }
  }
  evidence.bootstrap_replicates = bootstrap_replicates;

  std::vector<std::size_t> node_group_support_percent(sequences.size(), 0);
  for (std::size_t node = 0; node < sequences.size(); ++node) {
    node_group_support_percent[node] = bootstrap_replicates == 0
        ? 100
        : source_vb_round_nonnegative(
              (static_cast<double>(node_group_support[node]) + 1.0) /
              (static_cast<double>(bootstrap_replicates) + 1.0) * 100.0);
  }

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
  evidence.raw_tree_distances = raw_ranked.ranked_distances;
  evidence.raw_distance_rank_levels = raw_ranked.rank_levels;
  std::unordered_set<std::string> effective_collapsed_splits;
  if (bootstrap_replicates == 0) {
    // TestMoveInTreeAlt's active Reps=0 branch assigns tFAMat/tSAMat directly
    // into FCMat/SCMat.  Preserve identity at the stored-value level: passing
    // through the general float collapse buffer introduces sub-ULP changes
    // that MakePhPrScore can amplify even though no node was collapsed.
    evidence.collapsed_tree_distances = evidence.raw_tree_distances;
    evidence.collapsed_distance_rank_levels =
        evidence.raw_distance_rank_levels;
  } else {
    evidence.collapsed_tree_distances = source_collapsed_rank_tree_distances(
        raw_ranked,
        base_node_groups,
        node_group_support_percent,
        sequences.size(),
        evidence.collapsed_distance_rank_levels);
  }
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
