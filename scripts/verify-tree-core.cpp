#include "alignment.hpp"
#include "phylogeny.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {

int fail(const std::string& message) {
  std::cerr << "Event-tree core check failed: " << message << '\n';
  return 1;
}

std::string pseudo_random_sequence(std::size_t length, std::uint32_t seed) {
  std::string sequence(length, 'A');
  constexpr std::array<char, 4> bases{'A', 'C', 'G', 'T'};
  for (char& base : sequence) {
    seed = seed * 1664525U + 1013904223U;
    base = bases[(seed >> 16U) & 3U];
  }
  return sequence;
}

std::string periodic_mutant(const std::string& source, std::size_t period) {
  std::string sequence = source;
  constexpr std::string_view bases = "ACGT";
  for (std::size_t site = period / 2; site < sequence.size(); site += period) {
    const std::size_t state = bases.find(sequence[site]);
    sequence[site] = bases[(state + 1 + (site & 1U)) & 3U];
  }
  return sequence;
}

}  // namespace

int main() {
  const std::vector<std::uint32_t> expected_weights{
      1, 2, 0,
      1, 1, 2,
      1, 1, 2,
      1, 1, 1,
      1, 0, 0,
  };
  const auto weights = rdp::source_event_bootstrap_weights(5, 2, 3);
  if (weights != expected_weights) {
    return fail("the Microsoft-CRT SEQBOOT2 stream changed");
  }
  for (std::size_t replicate = 0; replicate < 3; ++replicate) {
    std::size_t sum = 0;
    for (std::size_t site = 0; site < 5; ++site) {
      sum += weights[site * 3 + replicate];
    }
    if (sum != 5) return fail("a retained/bootstrap replicate has the wrong site total");
  }

  constexpr std::size_t length = 240;
  const std::string first = pseudo_random_sequence(length, 0x12345678U);
  const std::string first_sibling = periodic_mutant(first, 31);
  const std::string second = periodic_mutant(first, 4);
  const std::string second_sibling = periodic_mutant(second, 29);
  const std::string third = pseudo_random_sequence(length, 0xa5a5f00dU);
  const std::string third_sibling = periodic_mutant(third, 27);
  auto parsed = rdp::build_alignment(
      "FASTA",
      {"first", "first-sibling", "second", "second-sibling", "third", "third-sibling"},
      {first, first_sibling, second, second_sibling, third, third_sibling});
  if (!parsed.ok()) return fail(parsed.error);

  std::vector<std::size_t> positions(length);
  std::iota(positions.begin(), positions.end(), 1);
  const std::vector<std::uint32_t> sequences{0, 1, 2, 3, 4, 5};
  const auto evidence = rdp::build_tree_region_evidence(
      parsed.alignment, sequences, positions, 0, 3);
  const auto repeated = rdp::build_tree_region_evidence(
      parsed.alignment, sequences, positions, 0, 3);

  if (!evidence.usable ||
      !evidence.source_clearcut_float_nj ||
      !evidence.source_ranked_tree_distances ||
      !evidence.source_midpoint_ultrametric ||
      !evidence.source_parent_rank_collapse ||
      !evidence.source_seqboot2_bootstrap ||
      !evidence.source_bootstrap_pseudocount ||
      evidence.bootstrap_random_seed != 3 || evidence.bootstrap_replicates != 0) {
    return fail("the supplied-source event-tree provenance is missing");
  }
  if (evidence.topology_node_count != 11 || evidence.topology_root != 10 ||
      evidence.topology_edges.size() != 10) {
    return fail("Clearcut-shaped rooted topology dimensions changed");
  }
  if (evidence.raw_tree_distances.size() != 36 ||
      evidence.collapsed_tree_distances.size() != 36 ||
      evidence.raw_distance_rank_levels == 0 ||
      evidence.collapsed_distance_rank_levels == 0 ||
      evidence.raw_distance_rank_levels > 15 ||
      evidence.collapsed_distance_rank_levels > 15) {
    return fail("rank-coded tree matrices have invalid dimensions or levels");
  }
  bool found_positive_rank = false;
  for (std::size_t first_index = 0; first_index < sequences.size(); ++first_index) {
    for (std::size_t second_index = 0; second_index < sequences.size(); ++second_index) {
      const std::size_t index = first_index * sequences.size() + second_index;
      const double value = evidence.collapsed_tree_distances[index];
      if (first_index == second_index) {
        if (value != 0.0) return fail("a ranked tree-matrix diagonal is nonzero");
      } else {
        const double thousandths = value * 1000.0;
        if (value < 0.0 || std::abs(thousandths - std::round(thousandths)) > 1e-12) {
          return fail("an analytical tree distance is not a non-negative 1/1000 rank");
        }
        found_positive_rank |= value > 0.0;
      }
      if (value != evidence.raw_tree_distances[index]) {
        return fail("the active zero-replicate path did not copy the raw rank matrix");
      }
      if (value != repeated.collapsed_tree_distances[index]) {
        return fail("identical source seeds did not reproduce the tree matrix");
      }
    }
  }
  if (!found_positive_rank) {
    return fail("Tree2ArrayP2 emitted no positive analytical rank");
  }
  if (!(evidence.jc(0, 1) < evidence.jc(0, 2)) ||
      !(evidence.jc(2, 3) < evidence.jc(2, 4))) {
    return fail("the base Jukes-Cantor matrix lost the fixture's close pairs");
  }
  if (!(evidence.tree(0, 1, false) < evidence.tree(0, 2, false)) ||
      !(evidence.tree(2, 3, false) < evidence.tree(2, 4, false))) {
    return fail("the source-ranked base topology lost the fixture's close clades");
  }
  for (const auto& edge : evidence.topology_edges) {
    if (std::abs(edge.length * 100000.0 - std::round(edge.length * 100000.0)) > 1e-8) {
      return fail("a displayed branch bypassed five-decimal source serialization");
    }
    if (edge.internal && (edge.bootstrap_support != 1.0 || edge.collapsed)) {
      return fail("the zero-replicate path applied bootstrap collapse metadata");
    }
  }

  const auto one_site = rdp::build_tree_region_evidence(
      parsed.alignment,
      {0, 1, 2},
      {1},
      0,
      3);
  if (!one_site.usable || one_site.site_count != 1) {
    return fail("the supplied FastBootDist valid-site rule regressed to a ten-site minimum");
  }

  std::cout
      << "Event-tree core verified: Microsoft-CRT SEQBOOT2 weights, supplied Clearcut-shaped "
         "float NJ, zero-replicate raw-tree copy, five-decimal branch serialization, "
         "and Tree2ArrayP2 midpoint analytical ranks ("
      << evidence.supported_internal_branches << '/' << evidence.internal_branches
      << " internal branches supported).\n";
  return 0;
}
