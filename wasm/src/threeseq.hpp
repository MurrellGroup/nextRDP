#pragma once

#include "maxchi.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace rdp {

struct ThreeSeqDiscoveryOptions {
  bool circular = true;
  bool correction_enabled = true;
  // TSXOver only invokes CheckSplit3Seq after at least one event has already
  // been erased from the working alignment.
  bool post_erasure_split_enabled = false;
  double p_value_cutoff = 0.05;
  std::uint64_t correction_tests = 1;
  // The supplied desktop path reads a precomputed exact table and falls back
  // to SiegmundDiscrete outside it. The browser computes the same finite walk
  // distribution on demand while it remains bounded, then uses the supplied
  // approximation. This cap limits one probability evaluation, not a scan.
  std::uint64_t maximum_exact_state_transitions = 8'000'000;
};

enum class ThreeSeqWalkDirection : std::int8_t {
  descent = -1,
  ascent = 1,
};

struct ThreeSeqDiscoveryCandidate {
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
  ThreeSeqWalkDirection direction = ThreeSeqWalkDirection::descent;
  std::size_t information_rich_sites = 0;
  std::size_t parent_one_matches = 0;
  std::size_t parent_two_matches = 0;
  std::size_t probability_excursion = 0;
  std::size_t maximum_excursion = 0;
  double raw_p_value = 1.0;
  double corrected_p_value = 1.0;
  std::array<double, 3> pair_similarity{};
  bool exact_probability = false;
  bool siegmund_fallback = false;
  bool missing_data_split_applied = false;
};

struct ThreeSeqDiscoverySummary {
  bool correction_applied = false;
  std::uint64_t correction_tests = 1;
  std::size_t target_profiles_scanned = 0;
  std::size_t exact_probability_evaluations = 0;
  std::size_t approximate_probability_evaluations = 0;
  std::size_t emitted_candidates = 0;
};

struct ThreeSeqRecheckEvidence {
  bool requested = false;
  bool representative_skipped = false;
  bool profile_available = false;
  bool correction_applied = false;
  std::uint64_t correction_tests = 1;
  std::size_t target_profiles_scanned = 0;
  std::size_t exact_probability_evaluations = 0;
  std::size_t approximate_probability_evaluations = 0;
  std::size_t qualifying_orientations = 0;
  std::size_t source_list_entries = 0;
  std::int8_t best_target = -1;
  ThreeSeqWalkDirection best_direction = ThreeSeqWalkDirection::descent;
  std::size_t information_rich_sites = 0;
  std::size_t parent_one_matches = 0;
  std::size_t parent_two_matches = 0;
  std::size_t probability_excursion = 0;
  std::size_t maximum_excursion = 0;
  std::size_t beginning = 0;
  std::size_t ending = 0;
  bool wraps_origin = false;
  double raw_p_value = 1.0;
  double corrected_p_value = 1.0;
  bool exact_probability = false;
  bool siegmund_fallback = false;
  bool missing_data_split_applied = false;
  bool source_recheck_hit = false;
};

struct ThreeSeqPlotProfile {
  bool available = false;
  std::vector<std::size_t> coordinates;
  std::array<std::vector<double>, 3> target_walks;
};

struct ThreeSeqProbabilityKey {
  std::size_t plus = 0;
  std::size_t minus = 0;
  std::size_t excursion = 0;

  bool operator==(const ThreeSeqProbabilityKey&) const = default;
};

struct ThreeSeqProbabilityKeyHash {
  std::size_t operator()(const ThreeSeqProbabilityKey& key) const noexcept;
};

struct ThreeSeqWorkspace {
  std::vector<std::size_t> coordinates;
  std::vector<std::int8_t> steps;
  std::vector<std::int64_t> heights;
  // Seq3PVals' supplied lookup table is Single. Keeping the compact DP in
  // float preserves its intended precision while minimizing WASM memory and
  // avoiding slow extended-precision arithmetic in large bounded profiles.
  std::vector<float> probability_state;
  std::vector<float> probability_next;
  std::unordered_map<ThreeSeqProbabilityKey, float, ThreeSeqProbabilityKeyHash>
      exact_probability_cache;
};

// Implements the supplied automated FindSubSeqTS/TS2 -> GetTSPVal/Seq3PVals
// -> CheckwrapC -> TSXOver route. All three candidate-recombinant
// rotations reuse the variable-site equality profile already prepared for
// MaxChi and CHIMAERA, so 3SEQ adds no alignment-byte pass.
[[nodiscard]] ThreeSeqDiscoverySummary threeseq_discover_prepared(
    const MaxChiWorkspace& variable_workspace,
    const std::vector<std::uint8_t>& triplet_missing_data,
    const std::array<double, 3>& pair_similarity,
    const ThreeSeqDiscoveryOptions& options,
    ThreeSeqWorkspace& workspace,
    std::vector<ThreeSeqDiscoveryCandidate>& output);

// Implements the supplied late TSXOver(1) shape. After the ordinary initial
// gate it evaluates both split orientations and accounts for the inverse-
// interval list copy emitted for each, without moving reconciled event bounds.
[[nodiscard]] ThreeSeqRecheckEvidence threeseq_recheck_prepared(
    const MaxChiWorkspace& variable_workspace,
    const std::vector<std::uint8_t>& triplet_missing_data,
    const ThreeSeqDiscoveryOptions& options,
    ThreeSeqWorkspace& workspace);

// Reconstructs the three target-specific hypergeometric random walks on the
// original alignment. The traces share a union coordinate grid for the web
// plot; probability-envelope permutations remain a separate manual mode.
[[nodiscard]] ThreeSeqPlotProfile threeseq_plot_profile(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    ThreeSeqWorkspace& workspace,
    MaxChiWorkspace& variable_workspace);

}  // namespace rdp
