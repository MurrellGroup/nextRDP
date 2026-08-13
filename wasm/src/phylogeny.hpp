#pragma once

#include "alignment.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace rdp {

struct TreeTopologyEdge {
  std::uint32_t first = 0;
  std::uint32_t second = 0;
  double length = 0.0;
  double bootstrap_support = 1.0;
  bool internal = false;
  bool collapsed = false;
};

struct TreeRegionEvidence {
  std::vector<std::uint32_t> sequences;
  std::vector<std::int32_t> alignment_to_tree;
  std::vector<double> jukes_cantor;
  std::vector<double> raw_tree_distances;
  std::vector<double> collapsed_tree_distances;
  std::vector<TreeTopologyEdge> topology_edges;
  std::size_t topology_node_count = 0;
  std::size_t topology_root = 0;
  std::size_t site_count = 0;
  std::size_t bootstrap_replicates = 0;
  std::size_t supported_internal_branches = 0;
  std::size_t internal_branches = 0;
  std::size_t raw_distance_rank_levels = 0;
  std::size_t collapsed_distance_rank_levels = 0;
  std::size_t negative_branches_normalized = 0;
  std::uint32_t bootstrap_random_seed = 0;
  bool source_clearcut_float_nj = false;
  bool source_ranked_tree_distances = false;
  bool source_midpoint_ultrametric = false;
  bool source_parent_rank_collapse = false;
  bool source_seqboot2_bootstrap = false;
  bool source_bootstrap_pseudocount = false;
  bool usable = false;

  [[nodiscard]] bool contains(std::uint32_t sequence) const;
  [[nodiscard]] double jc(std::uint32_t first, std::uint32_t second) const;
  [[nodiscard]] double tree(
      std::uint32_t first,
      std::uint32_t second,
      bool collapsed = true) const;
};

// The supplied event-tree path calls SEQBOOT2 with a retained unresampled
// replicate at index zero and bootstrap replicates at indices 1..reps. The
// returned vector is site-major with a stride of bootstrap_replicates + 1.
std::vector<std::uint32_t> source_event_bootstrap_weights(
    std::size_t site_count,
    std::size_t bootstrap_replicates,
    std::uint32_t seed);

double jukes_cantor_distance(
    const Alignment& alignment,
    std::uint32_t first,
    std::uint32_t second,
    const std::vector<std::size_t>& one_based_positions,
    std::size_t minimum_comparable_sites = 10);

TreeRegionEvidence build_tree_region_evidence(
    const Alignment& alignment,
    const std::vector<std::uint32_t>& sequences,
    const std::vector<std::size_t>& one_based_positions,
    std::size_t bootstrap_replicates,
    std::uint64_t seed);

}  // namespace rdp
