#pragma once

#include "alignment.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace rdp {

struct SiscanOptions {
  bool circular = true;
  bool bonferroni = true;
  double p_value_cutoff = 0.05;
  std::uint64_t correction_tests = 1;
  std::size_t window_sites = 200;
  std::size_t step_sites = 20;
  std::size_t scan_permutations = 100;
  std::size_t p_value_permutations = 1000;
  std::uint32_t random_seed = 3;
  bool fast_scan = true;
};

enum class SiscanScoreFamily : std::uint8_t {
  unavailable = 0,
  partition = 1,
  summed = 2,
};

struct SiscanDiscoveryCandidate {
  std::size_t beginning = 0;
  std::size_t ending = 0;
  bool wraps_origin = false;
  std::size_t informative_beginning = 0;
  std::size_t informative_ending = 0;
  std::uint8_t recombinant_local = 0;
  std::uint8_t major_parent_local = 1;
  std::uint8_t minor_parent_local = 2;
  std::uint8_t global_pair = 0;
  std::uint8_t candidate_pair = 0;
  std::uint32_t outlier_sequence = 0;
  std::size_t windows_in_region = 0;
  std::size_t informative_sites = 0;
  std::size_t permutation_draws = 0;
  std::uint8_t selected_score = 0;
  SiscanScoreFamily selected_score_family = SiscanScoreFamily::unavailable;
  double maximum_z = 0.0;
  double normal_tail_p_value = 1.0;
  double region_length_adjusted_p_value = 1.0;
  double window_adjusted_p_value = 1.0;
  double corrected_p_value = 1.0;
  std::array<double, 3> pair_similarity{};
  bool source_nearest_outlier = true;
  bool source_gap_stripping = true;
  bool source_variable_patterns = true;
  bool source_fast_window_quirk = true;
};

struct SiscanDiscoverySummary {
  bool profile_available = false;
  bool outlier_available = false;
  bool bonferroni_applied = false;
  bool source_nearest_outlier = true;
  bool source_gap_stripping = true;
  bool source_variable_patterns = true;
  bool source_fast_window_quirk = true;
  bool source_linear_window_scan = true;
  std::uint64_t correction_tests = 1;
  std::uint32_t outlier_sequence = 0;
  std::size_t window_sites = 200;
  std::size_t step_sites = 20;
  std::size_t scan_permutations = 100;
  std::size_t p_value_permutations = 1000;
  std::size_t windows_considered = 0;
  std::size_t windows_scored = 0;
  std::size_t windows_fast_skipped = 0;
  std::size_t candidate_regions_scored = 0;
  std::size_t emitted_candidates = 0;
  std::size_t permutation_draws = 0;
};

struct SiscanRecheckEvidence {
  bool requested = false;
  bool representative_skipped = false;
  bool profile_available = false;
  bool outlier_available = false;
  bool bonferroni_applied = false;
  bool source_nearest_outlier = true;
  bool source_gap_stripping = true;
  bool source_variable_patterns = true;
  std::uint64_t correction_tests = 1;
  std::uint32_t outlier_sequence = 0;
  std::size_t informative_sites = 0;
  std::size_t permutation_draws = 0;
  std::uint8_t global_pair = 0;
  std::int8_t scored_pair = -1;
  std::uint8_t selected_score = 0;
  SiscanScoreFamily selected_score_family = SiscanScoreFamily::unavailable;
  double maximum_z = 0.0;
  double normal_tail_p_value = 1.0;
  double region_length_adjusted_p_value = 1.0;
  double window_adjusted_p_value = 1.0;
  double corrected_p_value = 1.0;
  bool source_recheck_hit = false;
};

struct SiscanPlotProfile {
  bool available = false;
  std::uint32_t outlier_sequence = 0;
  std::size_t window_sites = 0;
  std::size_t step_sites = 0;
  std::size_t permutations = 0;
  std::vector<std::size_t> coordinates;
  std::array<std::vector<double>, 3> pair_z;
};

struct SiscanWorkspace {
  // Round-wide MakeDistanceBakB/WPGMA context. The desktop constructs this
  // once and GetSSOL reuses it for every triplet; the browser does the same.
  std::size_t context_sequence_count = 0;
  std::size_t context_alignment_length = 0;
  bool context_ready = false;
  std::vector<float> direct_similarity;
  std::vector<float> tree_similarity;
  std::vector<std::uint8_t> context_eligible;
  std::vector<std::uint32_t> context_origins;

  // MakeVRand is one seeded flat byte matrix. DoPerms3 consumes only a prefix
  // of it, so retaining an extensible prefix is both source-faithful and much
  // cheaper than regenerating the same random values for each window.
  std::uint32_t random_seed = 0;
  std::uint32_t random_state = 0;
  bool random_started = false;
  std::vector<std::uint8_t> vertical_random_prefix;

  std::vector<std::uint8_t> triplet_categories;
  std::vector<std::uint8_t> window_map;
  std::array<std::vector<double>, 16> partition_z;
  std::array<std::vector<double>, 15> summed_z;
  std::array<std::uint32_t, 16> pattern_counts{};
  std::vector<std::uint32_t> permutation_scores;
  std::vector<std::uint32_t> summed_scores;
  std::array<double, 16> final_partition_z{};
  std::array<double, 15> final_summed_z{};
  std::size_t context_builds = 0;
  std::size_t context_pair_comparisons = 0;
  std::size_t context_tree_merges = 0;
  std::size_t random_prefix_extensions = 0;
  std::size_t random_values_generated = 0;
};

// Invalidates only the state-dependent distance/tree context. The deterministic
// MakeVRand prefix survives cyclic rounds and remains reusable for the same
// seed, mirroring the desktop template file without requiring a filesystem.
void siscan_reset_round_context(SiscanWorkspace& workspace);

// Implements the supplied SSXoverC/GetSSOL/Get3Score/GetPScores2/DoPerms3/
// MakeZValue2/DoSums/FindMaxZ/ShrinkRegionC path. `origins` is the browser's
// TraceSub equivalent and prevents a cyclic fragment from becoming its own
// fourth sequence. `disabled_origins` is indexed by original sequence.
[[nodiscard]] SiscanDiscoverySummary siscan_discover(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    const std::vector<std::uint8_t>& triplet_missing_data,
    const std::array<double, 3>& pair_similarity,
    const std::vector<std::uint32_t>& origins,
    const std::vector<std::uint8_t>& disabled_origins,
    const SiscanOptions& options,
    SiscanWorkspace& workspace,
    std::vector<SiscanDiscoveryCandidate>& output);

// Fixed-bound confirmation for the supplied late SISCAN role. It applies the
// same fourth-sequence choice and full-region 1,000-permutation score without
// moving already reconciled event coordinates.
[[nodiscard]] SiscanRecheckEvidence siscan_recheck(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    const std::vector<std::uint8_t>& triplet_missing_data,
    const std::vector<std::uint32_t>& origins,
    const std::vector<std::uint8_t>& disabled_origins,
    std::size_t event_beginning,
    std::size_t event_ending,
    const SiscanOptions& options,
    SiscanWorkspace& workspace);

// Reconstructs the all-window SISCAN Z traces used when ShowPlotFlag disables
// the fast screen in the supplied UI. Each compact curve is the signed
// eligible P or S category with greatest absolute Z for that triplet pair at
// the window; the full 15-plus-9 desktop display remains a validation boundary.
[[nodiscard]] SiscanPlotProfile siscan_plot_profile(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    const std::vector<std::uint8_t>& triplet_missing_data,
    const std::vector<std::uint32_t>& origins,
    const std::vector<std::uint8_t>& disabled_origins,
    const SiscanOptions& options,
    SiscanWorkspace& workspace,
    std::int64_t fixed_outlier = -1);

}  // namespace rdp
