#pragma once

#include "alignment.hpp"
#include "maxchi.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace rdp {

struct ChimaeraDiscoveryOptions {
  bool circular = true;
  bool bonferroni = true;
  double p_value_cutoff = 0.05;
  std::uint64_t correction_tests = 1;
  std::size_t fixed_window_sites = 60;
  std::size_t maximum_peak_attempts = 100;
};

struct ChimaeraRecheckOptions {
  bool circular = true;
  bool bonferroni = true;
  double p_value_cutoff = 0.05;
  std::uint64_t correction_tests = 1;
  std::size_t fixed_window_sites = 60;
};

struct ChimaeraDiscoveryCandidate {
  std::size_t beginning = 0;
  std::size_t ending = 0;
  bool wraps_origin = false;
  std::size_t informative_beginning = 0;
  std::size_t informative_ending = 0;
  std::uint8_t target_local = 0;
  std::uint8_t recombinant_local = 0;
  std::uint8_t major_parent_local = 1;
  std::uint8_t minor_parent_local = 2;
  std::uint8_t candidate_pair = 0;
  MaxChiTractSide tract_side = MaxChiTractSide::unavailable;
  std::size_t peak_attempt = 0;
  std::size_t peak_alignment_position = 0;
  std::size_t information_rich_sites = 0;
  std::size_t initial_half_window = 0;
  std::size_t grown_half_window = 0;
  std::size_t critical_difference = 0;
  double maximum_chi_square = 0.0;
  double raw_p_value = 1.0;
  double within_triplet_p_value = 1.0;
  double corrected_p_value = 1.0;
  double left_flank_chi_square = 0.0;
  double right_flank_chi_square = 0.0;
  double inside_parent_one_match_rate = 0.0;
  double outside_parent_one_match_rate = 0.0;
  std::array<double, 3> pair_similarity{};
  bool missing_data_window_filter_applied = false;
  bool linear_edge_window_filter_applied = false;
};

struct ChimaeraDiscoverySummary {
  bool bonferroni_applied = false;
  std::uint64_t correction_tests = 1;
  std::size_t fixed_window_sites = 60;
  std::size_t target_profiles_scanned = 0;
  std::size_t peak_attempts = 0;
  std::size_t destroyed_peak_regions = 0;
  std::size_t emitted_candidates = 0;
  std::size_t peak_limit_targets = 0;
};

struct ChimaeraRecheckEvidence {
  bool requested = false;
  bool representative_skipped = false;
  bool profile_available = false;
  bool missing_data_window_filter_applied = false;
  bool linear_edge_window_filter_applied = false;
  bool bonferroni_applied = false;
  std::uint64_t correction_tests = 1;
  std::size_t fixed_window_sites = 60;
  std::size_t target_profiles_scanned = 0;
  std::int8_t best_target = -1;
  std::size_t information_rich_sites = 0;
  std::size_t half_window = 0;
  std::size_t critical_difference = 0;
  std::size_t grown_half_window = 0;
  std::size_t peak_alignment_position = 0;
  double maximum_chi_square = 0.0;
  double local_p_value = 1.0;
  double within_triplet_p_value = 1.0;
  double corrected_p_value = 1.0;
  bool source_recheck_hit = false;
};

struct ChimaeraPlotProfile {
  bool available = false;
  std::uint8_t target_local = 0;
  std::size_t half_window = 0;
  std::vector<std::size_t> coordinates;
  std::vector<double> chi_square;
};

// Implements the supplied AlistChi -> FastRecCheckChim -> CXoverA path. The
// full variable-site profile is shared with MaxChi, then each triplet member is
// treated in turn as the candidate recombinant. Sites where neither proposed
// parent matches that target are discarded before the binary chi-square scan.
[[nodiscard]] ChimaeraDiscoverySummary chimaera_discover_prepared(
    const MaxChiWorkspace& variable_workspace,
    const std::vector<std::uint8_t>& triplet_missing_data,
    const std::array<double, 3>& pair_similarity,
    const ChimaeraDiscoveryOptions& options,
    MaxChiWorkspace& target_workspace,
    std::vector<ChimaeraDiscoveryCandidate>& output);

// Source-shaped FastRecCheckChim secondary check. The caller supplies the
// already prepared three-pair variable profile so MaxChi and CHIMAERA late
// corroboration share one alignment-byte pass. Each possible target is
// screened, and the strongest corrected target statistic is retained without
// changing the reconciled event coordinates.
[[nodiscard]] ChimaeraRecheckEvidence chimaera_recheck_prepared(
    const MaxChiWorkspace& variable_workspace,
    const std::vector<std::uint8_t>& triplet_missing_data,
    const ChimaeraRecheckOptions& options,
    MaxChiWorkspace& target_workspace);

// Reconstructs the target-specific CHIMAERA chi-square trace for manual review.
[[nodiscard]] ChimaeraPlotProfile chimaera_plot_profile(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    std::uint8_t target_local,
    const std::vector<std::uint8_t>& triplet_missing_data,
    bool circular,
    std::size_t fixed_window_sites,
    double p_value_cutoff,
    MaxChiWorkspace& variable_workspace,
    MaxChiWorkspace& target_workspace);

}  // namespace rdp
