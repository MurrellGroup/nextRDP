#pragma once

#include "alignment.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <memory>
#include <unordered_map>
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

struct BootscanDiscoveryOptions {
  bool circular = true;
  bool bonferroni = true;
  double p_value_cutoff = 0.05;
  std::uint64_t correction_tests = 1;
  std::size_t window_sites = 200;
  std::size_t step_sites = 20;
  std::size_t bootstrap_replicates = 100;
  double support_cutoff = 0.70;
  std::uint32_t random_seed = 3;
  // The desktop BSXoverR path writes every pair/window bootstrap profile to
  // temporary files. WebAssembly has no persistent scratch filesystem, so a
  // bounded FIFO retains the same reusable pair summaries in memory.
  std::size_t pair_cache_limit_bytes = 64U * 1024U * 1024U;
};

struct BootscanDiscoveryCandidate {
  std::size_t beginning = 0;
  std::size_t ending = 0;
  bool wraps_origin = false;
  std::size_t informative_beginning = 0;
  std::size_t informative_ending = 0;
  std::uint8_t supported_pair = 0;
  std::uint8_t recombinant_local = 0;
  std::uint8_t major_parent_local = 2;
  std::uint8_t minor_parent_local = 1;
  std::uint8_t candidate_pair = 0;
  std::size_t windows_scored = 0;
  std::size_t usable_windows = 0;
  std::size_t informative_sites = 0;
  std::size_t tract_informative_sites = 0;
  std::size_t tract_pair_matches = 0;
  std::size_t outside_pair_matches = 0;
  double maximum_pair_support = 0.0;
  double mean_pair_support = 0.0;
  double bootstrap_p_value = 1.0;
  double raw_p_value = 1.0;
  double corrected_p_value = 1.0;
  std::array<double, 3> pair_similarity{};
  bool source_distance_mode = true;
  bool source_binomial_probability = true;
  bool source_strict_closest_pair_voting = true;
  bool erased_window_filter_applied = false;
};

struct BootscanDiscoverySummary {
  bool profile_available = false;
  bool bonferroni_applied = false;
  bool source_distance_mode = true;
  bool source_binomial_probability = true;
  bool source_strict_closest_pair_voting = true;
  bool erased_window_filter_applied = false;
  std::uint64_t correction_tests = 1;
  std::size_t window_sites = 200;
  std::size_t step_sites = 20;
  std::size_t bootstrap_replicates = 100;
  double support_cutoff = 0.70;
  std::size_t windows_scanned = 0;
  std::size_t candidate_regions_scored = 0;
  std::size_t emitted_candidates = 0;
  std::size_t pair_profiles_requested = 0;
  std::size_t pair_profile_cache_hits = 0;
  std::size_t pair_profile_cache_misses = 0;
  std::size_t pair_profile_cache_evictions = 0;
  std::size_t pair_profile_cache_bytes = 0;
  std::size_t pair_profile_cache_peak_bytes = 0;
};

struct BootscanPlotProfile {
  bool available = false;
  std::size_t window_sites = 0;
  std::size_t step_sites = 0;
  std::size_t bootstrap_replicates = 0;
  std::vector<std::size_t> coordinates;
  std::array<std::vector<double>, 3> pair_support;
};

struct BootscanPairDistanceProfile {
  std::vector<float> distances;
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

  // Automated BSXoverR pair/window distance summaries. Entries are shared so
  // a bounded-cache eviction cannot invalidate the three profiles currently
  // being combined into one triplet plot.
  std::unordered_map<
      std::uint64_t,
      std::shared_ptr<BootscanPairDistanceProfile>> pair_profile_cache;
  std::deque<std::uint64_t> pair_profile_fifo;
  std::size_t pair_profile_cache_bytes = 0;
  std::size_t pair_profile_cache_peak_bytes = 0;
  std::size_t pair_profile_cache_limit_bytes = 64U * 1024U * 1024U;
  std::size_t discovery_alignment_length = 0;
  std::size_t discovery_window_sites = 0;
  std::size_t discovery_step_sites = 0;
  std::size_t discovery_bootstrap_replicates = 0;
  std::uint32_t discovery_random_seed = 0;
  std::uint64_t pair_profile_cache_hits = 0;
  std::uint64_t pair_profile_cache_misses = 0;
  std::uint64_t pair_profile_cache_evictions = 0;
  std::vector<std::uint32_t> pair_valid_scratch;
  std::vector<std::uint32_t> pair_difference_scratch;
};

// Invalidates pair/window summaries between cyclic rounds. Unchanged triplets
// are already replayed from the XOverList-shaped shortlist; changed rows must
// never see distances calculated before the preceding tract erasure.
void bootscan_reset_discovery_cache(BootscanWorkspace& workspace);

// Implements the supplied automated BSXoverR -> GetPltVal -> ScanBSPlots ->
// MakeBSEvent distance-mode path. Seeded SEQBOOT2 weights are shared by every
// pair, pair/window distance profiles are cached across triplets, strict
// closest-pair votes define candidate regions, and MakeScoresBS-shaped
// binomial probabilities decide which regions enter the common signal list.
[[nodiscard]] BootscanDiscoverySummary bootscan_discover(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    const std::vector<std::uint8_t>& triplet_missing_data,
    const std::array<double, 3>& pair_similarity,
    const BootscanDiscoveryOptions& options,
    BootscanWorkspace& workspace,
    std::vector<BootscanDiscoveryCandidate>& output);

// Reconstructs the three strict closest-pair bootstrap-support curves for the
// review plot. It uses the same primary kernel but does not discover events.
[[nodiscard]] BootscanPlotProfile bootscan_plot_profile(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    const std::vector<std::uint8_t>& triplet_missing_data,
    const BootscanDiscoveryOptions& options,
    BootscanWorkspace& workspace);

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
