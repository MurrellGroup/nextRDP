#pragma once

#include "alignment.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace rdp {

enum class PhylproGapMode : std::uint8_t {
  ignore_pairwise = 0,
  strip_columns = 1,
};

struct PhylproOptions {
  std::size_t window_sites = 60;
  PhylproGapMode gap_mode = PhylproGapMode::ignore_pairwise;
  bool include_self = false;
  bool circular = true;
};

struct PhylproPoint {
  std::size_t alignment_position = 0;
  std::array<float, 3> correlations{1.0F, 1.0F, 1.0F};
};

struct PhylproProfile {
  std::size_t requested_window_sites = 0;
  std::size_t half_window_sites = 0;
  std::size_t eligible_columns = 0;
  std::size_t context_sequences = 0;
  std::size_t target_context_comparisons = 0;
  std::size_t rolling_updates = 0;
  bool window_capped = false;
  bool circular = true;
  bool include_self = false;
  PhylproGapMode gap_mode = PhylproGapMode::ignore_pairwise;
  std::array<float, 3> minimum_correlations{1.0F, 1.0F, 1.0F};
  std::array<std::size_t, 3> minimum_positions{};
  float minimum_value = 1.0F;
  float maximum_value = 1.0F;
  std::vector<std::size_t> eligible_alignment_positions;
  std::vector<PhylproPoint> points;
};

// Source-shaped PHYLPRO review profile. The supplied desktop path constructs
// every pairwise distance row even though its ordinary plot displays only the
// three selected roles. This implementation updates only those three rows;
// their Pearson coefficients are identical while work falls from O(L*N^2) to
// O(L*N).
[[nodiscard]] PhylproProfile phylpro_profile(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& target_sequences,
    const std::vector<std::uint8_t>& disabled,
    const PhylproOptions& options,
    std::string& error);

}  // namespace rdp
