#pragma once

#include "alignment.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace rdp {

struct MaxChiRecheckOptions {
  bool circular = true;
  bool bonferroni = true;
  double p_value_cutoff = 0.05;
  std::uint64_t correction_tests = 1;
  std::size_t fixed_window_sites = 70;
};

struct MaxChiDiscoveryOptions {
  bool circular = true;
  bool bonferroni = true;
  double p_value_cutoff = 0.05;
  std::uint64_t correction_tests = 1;
  std::size_t fixed_window_sites = 70;
  std::size_t maximum_peak_attempts = 100;
};

enum class MaxChiTractSide : std::int8_t {
  left = -1,
  unavailable = 0,
  right = 1,
};

struct MaxChiDiscoveryCandidate {
  std::size_t beginning = 0;
  std::size_t ending = 0;
  bool wraps_origin = false;
  std::size_t informative_beginning = 0;
  std::size_t informative_ending = 0;
  std::uint8_t recombinant_local = 0;
  std::uint8_t major_parent_local = 1;
  std::uint8_t minor_parent_local = 2;
  std::uint8_t candidate_pair = 0;
  std::int8_t peak_pair = -1;
  MaxChiTractSide tract_side = MaxChiTractSide::unavailable;
  std::size_t peak_attempt = 0;
  std::size_t peak_alignment_position = 0;
  std::size_t variable_sites = 0;
  std::size_t initial_half_window = 0;
  std::size_t grown_half_window = 0;
  std::size_t critical_difference = 0;
  double maximum_chi_square = 0.0;
  double raw_p_value = 1.0;
  double within_triplet_p_value = 1.0;
  double corrected_p_value = 1.0;
  double left_flank_chi_square = 0.0;
  double right_flank_chi_square = 0.0;
  std::array<double, 3> pair_similarity{};
  bool missing_data_window_filter_applied = false;
  bool linear_edge_window_filter_applied = false;
};

struct MaxChiDiscoverySummary {
  bool profile_available = false;
  bool missing_data_window_filter_applied = false;
  bool linear_edge_window_filter_applied = false;
  bool bonferroni_applied = false;
  bool peak_attempt_limit_reached = false;
  std::uint64_t correction_tests = 1;
  std::size_t variable_sites = 0;
  std::size_t fixed_window_sites = 70;
  std::size_t half_window = 0;
  std::size_t critical_difference = 0;
  std::size_t peak_attempts = 0;
  std::size_t destroyed_peak_regions = 0;
  std::size_t emitted_candidates = 0;
  double initial_maximum_chi_square = 0.0;
};

struct MaxChiPlotProfile {
  bool available = false;
  std::size_t half_window = 0;
  std::vector<std::size_t> coordinates;
  std::array<std::vector<double>, 3> chi_square;
};

struct MaxChiRecheckEvidence {
  bool requested = false;
  bool representative_skipped = false;
  bool profile_available = false;
  bool missing_data_window_filter_applied = false;
  bool linear_edge_window_filter_applied = false;
  bool bonferroni_applied = false;
  std::uint64_t correction_tests = 1;
  std::size_t variable_sites = 0;
  std::size_t fixed_window_sites = 70;
  std::size_t half_window = 0;
  std::size_t critical_difference = 0;
  std::size_t grown_half_window = 0;
  std::int8_t best_pair = -1;
  std::size_t peak_alignment_position = 0;
  double maximum_chi_square = 0.0;
  double local_p_value = 1.0;
  double within_triplet_p_value = 1.0;
  double corrected_p_value = 1.0;
  bool source_recheck_hit = false;
};

struct MaxChiWorkspace {
  std::vector<std::size_t> coordinates;
  std::array<std::vector<std::uint8_t>, 3> matches;
  std::array<std::vector<double>, 3> chi_values;
  std::array<std::vector<double>, 3> smooth_chi;
  std::vector<std::size_t> variable_prefix;
  std::vector<std::uint8_t> banned_windows;
  std::vector<std::uint8_t> missing_boundaries;
  std::vector<std::uint8_t> triplet_missing_data;
};

// Source-shaped exploratory MCXoverF path. Raw chi-square peaks are ordered
// across all three pair profiles, grown symmetrically, assigned to the left or
// right tract by FindSide, and removed with the supplied smoothed-peak retry
// lifecycle. The output vector is cleared but retains capacity between calls.
[[nodiscard]] MaxChiDiscoverySummary maxchi_discover(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    const std::vector<std::uint8_t>& triplet_missing_data,
    const MaxChiDiscoveryOptions& options,
    MaxChiWorkspace& workspace,
    std::vector<MaxChiDiscoveryCandidate>& output);

// Runs discovery from coordinates/matches/variable_prefix already populated in
// workspace. The combined RDP/MaxChi triplet path uses this to scan alignment
// bytes only once; maxchi_discover remains the standalone entry point.
[[nodiscard]] MaxChiDiscoverySummary maxchi_discover_prepared(
    const std::vector<std::uint8_t>& triplet_missing_data,
    const MaxChiDiscoveryOptions& options,
    MaxChiWorkspace& workspace,
    std::vector<MaxChiDiscoveryCandidate>& output);

// Reconstructs the three manual-review chi-square traces without running the
// event scheduler. Intended for the on-demand browser plot only.
[[nodiscard]] MaxChiPlotProfile maxchi_plot_profile(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    const std::vector<std::uint8_t>& triplet_missing_data,
    bool circular,
    std::size_t fixed_window_sites,
    double p_value_cutoff,
    MaxChiWorkspace& workspace);

// Source-shaped MaxChi recheck used by the supplied FastRecCheckMC2/AlistMC3
// path. It reports the strongest triplet statistic for review and late-list
// confirmation without repositioning an already reconciled event. The separate
// maxchi_discover entry point implements exploratory MCXoverF discovery.
[[nodiscard]] MaxChiRecheckEvidence maxchi_recheck(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    const std::vector<std::uint8_t>& triplet_missing_data,
    const MaxChiRecheckOptions& options,
    MaxChiWorkspace& workspace);

}  // namespace rdp
