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

bool is_source_support(double value) {
  constexpr std::array<double, 11> allowed{
      0.09, 0.18, 0.27, 0.36, 0.45, 0.55,
      0.64, 0.73, 0.82, 0.91, 1.00,
  };
  return std::any_of(allowed.begin(), allowed.end(), [value](double expected) {
    return std::abs(value - expected) < 1e-12;
  });
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
      parsed.alignment, sequences, positions, 10, 3);
  const auto repeated = rdp::build_tree_region_evidence(
      parsed.alignment, sequences, positions, 10, 3);

  if (!evidence.usable ||
      !evidence.source_clearcut_float_nj ||
      !evidence.source_ranked_tree_distances ||
      !evidence.source_midpoint_ultrametric ||
      !evidence.source_parent_rank_collapse ||
      !evidence.source_seqboot2_bootstrap ||
      !evidence.source_bootstrap_pseudocount ||
      evidence.bootstrap_random_seed != 3) {
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
  for (std::size_t first_index = 0; first_index < sequences.size(); ++first_index) {
    for (std::size_t second_index = 0; second_index < sequences.size(); ++second_index) {
      const std::size_t index = first_index * sequences.size() + second_index;
      const double value = evidence.collapsed_tree_distances[index];
      if (first_index == second_index) {
        if (value != 0.0) return fail("a ranked tree-matrix diagonal is nonzero");
      } else {
        const double thousandths = value * 1000.0;
        if (value <= 0.0 || std::abs(thousandths - std::round(thousandths)) > 1e-12) {
          return fail("an analytical tree distance is not a positive 1/1000 rank");
        }
      }
      if (value != repeated.collapsed_tree_distances[index]) {
        return fail("identical source seeds did not reproduce the tree matrix");
      }
    }
  }
  for (std::size_t first_index = 0; first_index + 2 < sequences.size(); ++first_index) {
    for (std::size_t second_index = first_index + 1;
         second_index + 1 < sequences.size();
         ++second_index) {
      for (std::size_t third_index = second_index + 1;
           third_index < sequences.size();
           ++third_index) {
        std::array<double, 3> raw_distances{
            evidence.tree(first_index, second_index, false),
            evidence.tree(first_index, third_index, false),
            evidence.tree(second_index, third_index, false),
        };
        std::sort(raw_distances.begin(), raw_distances.end());
        if (raw_distances[1] != raw_distances[2]) {
          return fail("the TreeMidP/UltraTreeDistP matrix is not ultrametric");
        }
      }
    }
  }
  for (std::size_t first_index = 0; first_index < sequences.size(); ++first_index) {
    for (std::size_t second_index = first_index + 1;
         second_index < sequences.size();
         ++second_index) {
      const double collapsed = evidence.tree(first_index, second_index, true);
      bool found_original_rank = false;
      for (std::size_t raw_first = 0; raw_first < sequences.size(); ++raw_first) {
        for (std::size_t raw_second = raw_first + 1;
             raw_second < sequences.size();
             ++raw_second) {
          found_original_rank |=
              collapsed == evidence.tree(raw_first, raw_second, false);
        }
      }
      if (!found_original_rank) {
        return fail("CollapseNodes-style promotion compressed or invented a rank");
      }
    }
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
    if (edge.internal && !is_source_support(edge.bootstrap_support)) {
      return fail("an internal branch bypassed the replicate-zero pseudocount percentage");
    }
  }

  bool exercised_parent_rank_collapse = false;
  for (std::uint32_t fixture = 1; fixture <= 32 && !exercised_parent_rank_collapse; ++fixture) {
    const std::string ancestor = pseudo_random_sequence(120, 0x9e3779b9U ^ fixture);
    std::vector<std::string> star_sequences(6, ancestor);
    std::uint32_t mutation_state = 0x85ebca6bU * fixture;
    constexpr std::string_view bases = "ACGT";
    for (std::size_t sequence = 0; sequence < star_sequences.size(); ++sequence) {
      for (char& base : star_sequences[sequence]) {
        mutation_state = mutation_state * 1664525U + 1013904223U;
        if ((mutation_state >> 24U) >= 48U) continue;
        const std::size_t current = bases.find(base);
        base = bases[(current + 1U + ((mutation_state >> 16U) % 3U)) & 3U];
      }
    }
    auto star = rdp::build_alignment(
        "FASTA",
        {"star-0", "star-1", "star-2", "star-3", "star-4", "star-5"},
        std::move(star_sequences));
    if (!star.ok()) return fail(star.error);
    std::vector<std::size_t> star_positions(120);
    std::iota(star_positions.begin(), star_positions.end(), 1);
    const auto weak = rdp::build_tree_region_evidence(
        star.alignment, sequences, star_positions, 10, 3);
    const bool edge_collapsed = std::any_of(
        weak.topology_edges.begin(),
        weak.topology_edges.end(),
        [](const rdp::TreeTopologyEdge& edge) { return edge.collapsed; });
    if (!edge_collapsed) continue;
    if (weak.raw_tree_distances == weak.collapsed_tree_distances) {
      return fail("a weak source node was labelled collapsed without parent-rank promotion");
    }
    for (std::size_t first_index = 0; first_index < sequences.size(); ++first_index) {
      for (std::size_t second_index = first_index + 1;
           second_index < sequences.size();
           ++second_index) {
        const double collapsed = weak.tree(first_index, second_index, true);
        bool retained_rank = false;
        for (std::size_t raw_first = 0; raw_first < sequences.size(); ++raw_first) {
          for (std::size_t raw_second = raw_first + 1;
               raw_second < sequences.size();
               ++raw_second) {
            retained_rank |= collapsed == weak.tree(raw_first, raw_second, false);
          }
        }
        if (!retained_rank) {
          return fail("CollapseNodes promotion invented or recompressed a weak-node rank");
        }
      }
    }
    exercised_parent_rank_collapse = true;
  }
  if (!exercised_parent_rank_collapse) {
    return fail("the deterministic weak-tree corpus did not exercise parent-rank promotion");
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
         "float NJ, replicate-zero support, 50% collapse, five-decimal branch serialization, "
         "and midpoint-rooted ultrametric analytical ranks ("
      << evidence.supported_internal_branches << '/' << evidence.internal_branches
      << " internal branches supported).\n";
  return 0;
}
