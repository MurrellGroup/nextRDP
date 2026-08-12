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

struct MaxChiRecheckEvidence {
  bool requested = false;
  bool representative_skipped = false;
  bool profile_available = false;
  bool missing_data_window_filter_applied = false;
  bool linear_edge_window_filter_applied = false;
  bool bonferroni_applied = false;
  std::uint64_t correction_tests = 1;
  std::size_t variable_sites = 0;
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
  std::vector<std::size_t> variable_prefix;
  std::vector<std::uint8_t> banned_windows;
  std::vector<std::uint8_t> missing_boundaries;
  std::vector<std::uint8_t> triplet_missing_data;
};

// Source-shaped MaxChi recheck used by the supplied FastRecCheckMC2/AlistMC3
// path. This deliberately does not discover or authoritatively reposition an
// event: it reports the strongest triplet statistic for review and late-list
// confirmation while primary RDP remains the event-discovery method.
[[nodiscard]] MaxChiRecheckEvidence maxchi_recheck(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    const std::vector<std::uint8_t>& triplet_missing_data,
    const MaxChiRecheckOptions& options,
    MaxChiWorkspace& workspace);

}  // namespace rdp
