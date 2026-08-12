#pragma once

#include "alignment.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace rdp {

struct BreakpointConfidenceBoundary {
  std::int64_t confidence_99_beginning = -1;
  std::int64_t confidence_99_ending = -1;
  std::int64_t hmm_coordinate = -1;
  std::int64_t confidence_95_beginning = -1;
  std::int64_t confidence_95_ending = -1;
  std::size_t input_coordinate = 0;
  std::size_t polished_coordinate = 0;
  bool interval_available = false;
  bool source_interval_contains_input = false;
  bool confidence_99_wraps_origin = false;
  bool confidence_95_wraps_origin = false;
  bool repositioned = false;
  bool missing_data_adjusted = false;
  bool final_gap_adjusted = false;
};

struct BreakpointConfidenceEvidence {
  std::array<BreakpointConfidenceBoundary, 2> boundaries;
  std::size_t information_rich_sites = 0;
  std::size_t candidate_interval_count = 0;
  std::size_t polished_beginning = 0;
  std::size_t polished_ending = 0;
  double best_log_likelihood = 0.0;
  std::uint32_t random_seed = 3;
  std::size_t hmm_cycles_argument = 20;
  std::size_t serial_training_starts = 21;
  bool attempted = false;
  bool available = false;
  bool applied_to_event = false;
  bool single_transition_assignment = false;
  bool insufficient_inside_or_outside_reverted = false;
  bool source_random_adapter = true;
  std::string unavailable_reason = "not-run";
};

// The supplied DLL repeatedly allocates alignment-sized HMM matrices. The
// browser keeps the same flattened state order but retains vector capacity
// across events, which removes allocator churn without changing the path.
struct BurtConfidenceWorkspace {
  std::vector<std::uint8_t> symbols;
  std::vector<std::uint8_t> expanded_symbols;
  std::vector<std::size_t> information_coordinates;
  std::vector<std::size_t> position_information_counts;
  std::vector<double> lattice;
  std::vector<double> backpointers;
  std::vector<std::int32_t> path;
  std::vector<std::int32_t> best_path;
  std::vector<double> forward;
  std::vector<double> reverse;
  std::vector<double> posterior;
};

[[nodiscard]] BreakpointConfidenceEvidence source_polish_breakpoints(
    const Alignment& current_alignment,
    const std::vector<std::uint8_t>& triplet_missing_data,
    const std::array<std::uint32_t, 3>& representatives,
    std::size_t beginning,
    std::size_t ending,
    bool circular,
    BurtConfidenceWorkspace& workspace);

}  // namespace rdp
