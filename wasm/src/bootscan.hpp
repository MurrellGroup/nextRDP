#pragma once

#include "alignment.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace rdp {

struct BootscanRecheckOptions {
  bool circular = true;
  bool bonferroni = true;
  double p_value_cutoff = 0.05;
  std::uint64_t correction_tests = 1;
  std::size_t window_sites = 200;
  std::size_t step_sites = 20;
  std::size_t bootstrap_replicates = 100;
  double support_cutoff = 0.70;
  std::uint32_t random_seed = 3;
};

struct BootscanRecheckEvidence {
  bool requested = false;
  bool representative_skipped = false;
  bool profile_available = false;
  bool source_distance_mode = true;
  bool source_binomial_probability = true;
  bool source_circular_windows = true;
  bool erased_window_filter_applied = false;
  bool bonferroni_applied = false;
  std::uint64_t correction_tests = 1;
  std::size_t window_sites = 200;
  std::size_t step_sites = 20;
  std::size_t bootstrap_replicates = 100;
  std::uint32_t random_seed = 3;
  double support_cutoff = 0.70;
  std::size_t windows_scanned = 0;
  std::size_t event_windows_scored = 0;
  std::size_t usable_event_windows = 0;
  std::size_t informative_sites = 0;
  std::size_t tract_informative_sites = 0;
  std::size_t tract_pair_matches = 0;
  std::size_t outside_pair_matches = 0;
  std::int8_t scored_pair = -1;
  double maximum_pair_support = 0.0;
  double mean_scored_pair_support = 0.0;
  double local_p_value = 1.0;
  double corrected_p_value = 1.0;
  bool support_gate_passed = false;
  bool source_recheck_hit = false;
};

struct BootscanWorkspace {
  // SEQBOOT2 stores weights site-major with replicate zero holding the
  // unpermuted window. Reuse these buffers across event/list rechecks.
  std::vector<std::uint32_t> bootstrap_weights;
  std::vector<std::array<std::uint32_t, 3>> support_counts;
  std::vector<std::uint8_t> usable_windows;
  std::vector<std::size_t> event_window_indices;
  std::array<std::vector<std::uint8_t>, 3> pair_scores;
  std::vector<std::size_t> position_to_informative;
};

// Implements the supplied distance-mode BSXoverM/DrawBSPlotsIII confirmation
// path. It regenerates the seeded SEQBOOT2 weights, uses the source's strict
// closest-pair voting and event-window support gate, then calculates the
// ordinary binomial BOOTSCAN probability without moving reconciled bounds.
[[nodiscard]] BootscanRecheckEvidence bootscan_recheck(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    const std::vector<std::uint8_t>& triplet_missing_data,
    std::size_t event_beginning,
    std::size_t event_ending,
    const BootscanRecheckOptions& options,
    BootscanWorkspace& workspace);

}  // namespace rdp
