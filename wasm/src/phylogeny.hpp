#pragma once

#include "alignment.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace rdp {

struct TreeRegionEvidence {
  std::vector<std::uint32_t> sequences;
  std::vector<std::int32_t> alignment_to_tree;
  std::vector<double> jukes_cantor;
  std::vector<double> raw_patristic;
  std::vector<double> collapsed_patristic;
  std::size_t site_count = 0;
  std::size_t bootstrap_replicates = 0;
  std::size_t supported_internal_branches = 0;
  std::size_t internal_branches = 0;
  bool usable = false;

  [[nodiscard]] bool contains(std::uint32_t sequence) const;
  [[nodiscard]] double jc(std::uint32_t first, std::uint32_t second) const;
  [[nodiscard]] double tree(
      std::uint32_t first,
      std::uint32_t second,
      bool collapsed = true) const;
};

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
