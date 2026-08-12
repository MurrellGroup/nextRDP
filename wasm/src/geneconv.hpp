#pragma once

#include "alignment.hpp"
#include "maxchi.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace rdp {

struct GeneconvDiscoveryOptions {
  bool circular = true;
  bool bonferroni = true;
  double p_value_cutoff = 0.05;
  std::uint64_t correction_tests = 1;
  std::size_t mismatch_scale = 1;
  std::size_t maximum_overlapping_fragments = 1;
};

struct GeneconvDiscoveryCandidate {
  std::size_t beginning = 0;
  std::size_t ending = 0;
  bool wraps_origin = false;
  std::size_t informative_beginning = 0;
  std::size_t informative_ending = 0;
  std::uint8_t track = 0;
  std::uint8_t recombinant_local = 0;
  std::uint8_t major_parent_local = 2;
  std::uint8_t minor_parent_local = 1;
  std::uint8_t candidate_pair = 0;
  std::size_t polymorphic_sites = 0;
  std::size_t positive_sites = 0;
  std::size_t discordant_sites = 0;
  std::size_t mismatch_penalty = 0;
  std::size_t fragment_score = 0;
  std::size_t critical_score = 0;
  double lambda = 0.0;
  double karlin_altschul_k = 0.0;
  double raw_p_value = 1.0;
  double corrected_p_value = 1.0;
  std::array<double, 3> pair_similarity{};
  bool karlin_altschul_probability = true;
  bool ignored_indels = true;
  bool overlap_filter_applied = true;
  bool minimum_fragment_filters_applied = false;
};

struct GeneconvDiscoverySummary {
  bool profile_available = false;
  bool bonferroni_applied = false;
  bool ignored_indels = true;
  bool overlap_filter_applied = true;
  bool minimum_fragment_filters_applied = false;
  bool source_skew_filter_rejected = false;
  std::uint64_t correction_tests = 1;
  std::size_t polymorphic_sites = 0;
  std::size_t tracks_screened = 0;
  std::size_t fragments_scored = 0;
  std::size_t qualified_fragments = 0;
  std::size_t overlap_rejected_fragments = 0;
  std::size_t numerical_fallback_tracks = 0;
  std::size_t emitted_candidates = 0;
};

struct GeneconvRecheckEvidence {
  bool requested = false;
  bool representative_skipped = false;
  bool profile_available = false;
  bool bonferroni_applied = false;
  bool ignored_indels = true;
  bool overlap_filter_applied = true;
  bool minimum_fragment_filters_applied = false;
  bool source_skew_filter_rejected = false;
  std::uint64_t correction_tests = 1;
  std::size_t polymorphic_sites = 0;
  std::size_t tracks_screened = 0;
  std::size_t fragments_scored = 0;
  std::size_t qualified_fragments = 0;
  std::size_t overlap_rejected_fragments = 0;
  std::size_t numerical_fallback_tracks = 0;
  std::int8_t best_track = -1;
  std::int8_t recombinant_local = -1;
  std::size_t beginning = 0;
  std::size_t ending = 0;
  bool wraps_origin = false;
  std::size_t fragment_score = 0;
  std::size_t critical_score = 0;
  double raw_p_value = 1.0;
  double corrected_p_value = 1.0;
  bool source_recheck_hit = false;
};

struct GeneconvRun {
  std::size_t beginning = 0;
  std::size_t ending = 0;
  std::size_t length = 0;
  bool positive = false;
  bool wraps_origin = false;
};

struct GeneconvCategoryRun {
  std::uint8_t category = 3;
  std::size_t beginning = 0;
  std::size_t ending = 0;
  std::size_t length = 0;
  bool wraps_origin = false;
};

struct GeneconvScoredFragment {
  std::uint8_t track = 0;
  std::size_t beginning = 0;
  std::size_t ending = 0;
  bool wraps_origin = false;
  std::size_t positive_sites = 0;
  std::size_t discordant_sites = 0;
  std::size_t mismatch_penalty = 0;
  std::size_t fragment_score = 0;
  std::size_t critical_score = 0;
  double lambda = 0.0;
  double karlin_altschul_k = 0.0;
  double raw_p_value = 1.0;
};

struct GeneconvWorkspace {
  std::vector<std::uint8_t> categories;
  std::vector<GeneconvCategoryRun> category_runs;
  std::array<std::vector<GeneconvRun>, 6> runs;
  std::vector<GeneconvScoredFragment> scored_fragments;
  std::vector<std::uint32_t> overlap_tree_max;
  std::vector<std::uint32_t> overlap_tree_lazy;
  std::vector<std::int64_t> prefix_scores;
  std::vector<std::size_t> prefix_positive_sites;
  std::vector<std::size_t> prefix_discordant_sites;
  std::vector<std::size_t> next_lower_prefix;
  std::vector<std::int64_t> range_max_values;
  std::vector<std::size_t> range_max_indices;
  std::vector<std::size_t> monotonic_stack;
};

struct GeneconvPlotProfile {
  bool available = false;
  std::vector<std::size_t> coordinates;
  std::array<std::vector<double>, 3> negative_log10_p_value;
};

// Implements the supplied FindSubSeqGCAP6 -> GetFragsP ->
// GetMaxFragScoreP -> CalcKMaxP -> GCCalcPValP2 -> GCXoverD ordinary
// automated triplet path. The MaxChi workspace already contains the exact
// non-monomorphic coordinates and three pair-match bits, so the combined
// primary-method scheduler does not rescan alignment bytes.
[[nodiscard]] GeneconvDiscoverySummary geneconv_discover_prepared(
    const MaxChiWorkspace& variable_workspace,
    const std::array<double, 3>& pair_similarity,
    const GeneconvDiscoveryOptions& options,
    GeneconvWorkspace& workspace,
    std::vector<GeneconvDiscoveryCandidate>& output);

// Runs the same supplied ordinary six-track kernel as a non-coordinate-
// changing representative/final-list corroboration pass. The caller provides
// an already prepared MaxChi workspace so MaxChi, CHIMAERA, and GENECONV late
// checks share one alignment-byte scan.
[[nodiscard]] GeneconvRecheckEvidence geneconv_recheck_prepared(
    const MaxChiWorkspace& variable_workspace,
    const std::array<double, 3>& pair_similarity,
    const GeneconvDiscoveryOptions& options,
    GeneconvWorkspace& workspace,
    std::vector<GeneconvDiscoveryCandidate>& candidates);

// Reconstructs a compact three-colour -log10(KA P) fragment envelope for
// manual review. Inner/outer tracks share the colours used by the supplied
// plot: 0/5, 1/4, and 2/3.
[[nodiscard]] GeneconvPlotProfile geneconv_plot_profile(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    const GeneconvDiscoveryOptions& options,
    GeneconvWorkspace& workspace);

}  // namespace rdp
