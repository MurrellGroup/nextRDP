#pragma once

#include "alignment.hpp"
#include "burt_confidence.hpp"
#include "maxchi.hpp"

#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace rdp {

enum class CorrectionMode : std::uint8_t {
  bonferroni = 0,
  none = 1,
};

enum class ReviewState : std::uint8_t {
  unreviewed = 0,
  accepted = 1,
  rejected = 2,
};

struct ScanOptions {
  bool circular = true;
  CorrectionMode correction = CorrectionMode::bonferroni;
  double p_value_cutoff = 0.05;
  std::size_t window_sites = 30;
  bool polish_breakpoints = true;
  std::vector<std::uint8_t> mask;
  std::vector<std::uint8_t> disabled;
};

[[nodiscard]] std::string curated_sequences_fasta(
    const Alignment& alignment,
    const std::vector<std::uint8_t>& mask,
    const std::vector<std::uint8_t>& disabled,
    bool include_enabled);

struct Signal {
  std::uint32_t id = 0;
  std::array<std::uint32_t, 3> triplet{};
  std::uint32_t recombinant = 0;
  std::uint32_t major_parent = 0;
  std::uint32_t minor_parent = 0;
  std::size_t beginning = 0;
  std::size_t ending = 0;
  bool wraps_origin = false;
  std::size_t informative_beginning = 0;
  std::size_t informative_ending = 0;
  double local_p_value = 1.0;
  double corrected_p_value = 1.0;
  std::uint64_t correction_tests = 0;
  std::array<double, 3> pair_similarity{};
  std::size_t informative_sites = 0;
  std::uint8_t candidate_pair = 0;
  bool fragment_assisted = false;
  std::array<std::int32_t, 3> fragment_event_context{-1, -1, -1};
  ReviewState review_state = ReviewState::unreviewed;
  std::int32_t event_id = -1;
};

struct TraceEvidence {
  std::uint32_t sequence = 0;
  std::size_t beginning = 0;
  std::size_t ending = 0;
  bool wraps_origin = false;
  double local_p_value = 1.0;
  double corrected_p_value = 1.0;
  bool significant = false;
};

struct FinalTrimMatrixEvidence {
  double collapsed_tree_position_score = 0.0;
  double raw_tree_position_score = 0.0;
  double relative_distance_score = 0.0;
  std::array<double, 2> breakpoint_distance_scores{};
  std::array<std::uint8_t, 2> breakpoint_score_available{};
  double detected_region_distance_score = 0.0;
  std::array<double, 7> detected_region_match_distances{};
  double detected_event_overlap = 0.0;
  std::size_t detected_event_beginning = 0;
  std::size_t detected_event_ending = 0;
  std::int32_t detected_event_signal_id = -1;
  double active_consensus_matrix_score = 0.0;
  bool detected_event_match = false;
  bool detected_region_saturated = false;
  bool source_sequence_index_quirk_applied = true;
  bool tree_distance_fallback = false;
  bool applies_to_nonrepresentative = false;
};

struct CalcMatchEvidence {
  bool available = false;
  std::size_t fragment_variable_sites = 0;
  std::size_t target_half_window = 0;
  std::size_t smoothing_half_window = 0;
  double regional_match_score = 0.0;
  std::int8_t raw_breakpoint_match_class = 0;
  std::int8_t breakpoint_match_class = 0;
  std::array<double, 6> checkpoint_matches{};
  std::array<std::uint8_t, 2> breakpoints_exist{};
  bool topology_filtered = false;
  bool topology_distance_fallback = false;
};

struct PostGroupRdpRecheckEvidence {
  bool requested = false;
  bool representative_skipped = false;
  bool profile_available = false;
  double local_p_value_cutoff = 1.0;
  std::size_t emitted_signal_count = 0;
  std::size_t candidate_signal_count = 0;
  std::size_t overlapping_signal_count = 0;
  bool event_redetected = false;
  bool significant = false;
  std::size_t best_beginning = 0;
  std::size_t best_ending = 0;
  bool best_wraps_origin = false;
  double best_overlap = 0.0;
  double best_local_p_value = 1.0;
  double best_corrected_p_value = 1.0;
};

struct ConsensusScoreEvidence {
  double corrected_correlation_p_value = 0.0;
  double detected_event_overlap = 0.0;
  double pattern_score = 0.0;
  double maximum_direct_correlation = 0.0;
  double base_score_before_final_membership = 0.0;
  double source_long_matrix_multiplier = 0.0;
  double score_after_final_membership = 0.0;
  double score_after_rcorrx = 0.0;
  double final_score = 0.0;
  bool detectable_set_member = false;
  bool initial_rlist_member = false;
  bool duplicate_cleaned_member = false;
  bool nearest_nonrecombinant_member = false;
  bool final_trim_first_expansion_added = false;
  bool final_trim_second_expansion_added = false;
  bool selected_role_pruned_out = false;
  bool final_trim_member = false;
  bool consensus_primary_member = false;
  bool consensus_equivalent_member = false;
  bool consensus_straggler_member = false;
  bool consensus_rebuilt_member = false;
  bool consensus_fallback_restored = false;
  bool selected_tree_cleanup_pruned_out = false;
  bool selected_tree_cleanup_added = false;
  bool final_distance_member = false;
  bool representative_sentinel = false;
  bool other_representative_zero = false;
  bool complete = false;
};

struct DistanceCorrelationEvidence {
  std::uint32_t sequence = 0;
  std::array<double, 3> correlations{};
  std::array<double, 3> direct_correlations{};
  std::array<double, 3> p_values{1.0, 1.0, 1.0};
  std::array<std::uint8_t, 3> inversion_codes{};
  std::array<std::uint8_t, 3> warning_filtered{};
  std::array<std::uint8_t, 3> duplicate_filtered{};
  std::array<std::size_t, 3> minimum_comparable_sites{};
  std::array<std::size_t, 2> breakpoint_overlap_sites{};
  double aggregate_score = 0.0;
  double aggregate_target = 0.9;
  std::int8_t best_matrix_pair = -1;
  bool overlap_eligible = false;
  bool acceptable_affinity = false;
  bool strong_correlation_override = false;
  bool significant = false;
  bool detectable_support = false;
  bool positive_support = false;
  bool inverse_support = false;
  bool stripped_inverse_only = false;
  bool duplicate_cleaned_support = false;
  FinalTrimMatrixEvidence final_trim_matrix;
  CalcMatchEvidence calc_match;
  ConsensusScoreEvidence consensus_score;
  PostGroupRdpRecheckEvidence post_group_rdp_recheck;
  MaxChiRecheckEvidence post_group_maxchi_recheck;
};

struct PhylogeneticCorrelationEvidence {
  std::uint32_t sequence = 0;
  std::array<double, 3> collapsed_affinity_margin{};
  std::array<double, 3> raw_affinity_margin{};
  std::array<std::uint8_t, 3> collapsed_pair_support{};
  std::array<std::uint8_t, 3> raw_pair_support{};
  std::int8_t best_tree_pair = -1;
  std::uint8_t supporting_tree_pairs = 0;
  bool included = false;
  bool distance_fallback = false;
  bool disabled_excluded = false;
};

struct TreeTopologyEdgeSummary {
  std::uint32_t first = 0;
  std::uint32_t second = 0;
  double length = 0.0;
  double bootstrap_support = 1.0;
  bool internal = false;
  bool collapsed = false;
};

struct TreeRegionSummary {
  std::size_t sites = 0;
  std::size_t sequences = 0;
  std::size_t bootstrap_replicates = 0;
  std::size_t supported_internal_branches = 0;
  std::size_t internal_branches = 0;
  std::size_t node_count = 0;
  std::size_t root = 0;
  std::vector<TreeTopologyEdgeSummary> topology_edges;
  bool usable = false;
};

struct TreePanelLeafSummary {
  std::uint32_t working_sequence = 0;
  std::uint32_t original_sequence = 0;
  std::int32_t fragment_event = -1;
};

enum class RoleMetricKind : std::uint8_t {
  phpr = 0,
  tree_phpr = 1,
  collapsed_tree_phpr = 2,
  sub_phpr = 3,
  tree_sub_phpr = 4,
  subdist = 5,
  tree_subdist = 6,
  triplet_score = 7,
  three_set_support = 8,
};

struct RoleMetricEvidence {
  RoleMetricKind kind = RoleMetricKind::phpr;
  std::array<double, 3> scores{};
  std::array<double, 3> contributions{};
  double weight = 0.0;
  std::int8_t winning_role = -1;
  bool higher_is_recombinant = false;
  bool informative = false;
};

struct RoleConsensusEvidence {
  std::vector<RoleMetricEvidence> metrics;
  std::array<double, 3> votes{};
  std::int8_t recommended_role = -1;
  std::uint32_t recommended_recombinant = 0;
  std::uint32_t recommended_major_parent = 0;
  std::uint32_t recommended_minor_parent = 0;
  double confidence = 0.0;
  bool informative = false;
};

struct RoleHypothesisEvidence {
  std::uint32_t presumed_recombinant = 0;
  std::uint32_t parent_one = 0;
  std::uint32_t parent_two = 0;
  std::vector<std::uint32_t> detectable_signal_set;
  std::vector<std::uint32_t> distance_correlation_set;
  std::vector<std::uint32_t> phylogenetic_correlation_set;
  std::vector<std::uint32_t> complete_two_of_three_set;
  std::vector<DistanceCorrelationEvidence> distance_evidence;
  std::vector<PhylogeneticCorrelationEvidence> phylogenetic_evidence;
  std::array<std::uint8_t, 3> correlation_warnings{};
  std::size_t tested_sequences = 0;
  std::size_t valid_sequences = 0;
};

struct BreakpointUncertaintyEvidence {
  bool native_check_ends_applied = false;
  bool information_profile_available = false;
  bool native_check_ends_warning = false;
  bool input_missing_data_in_range = false;
  bool linear_edge_within_window = false;
  std::size_t check_range_beginning = 0;
  std::size_t check_range_ending = 0;
  std::size_t check_coordinate_count = 0;
  bool check_range_wraps_origin = false;
};

struct UniqueEvent {
  std::uint32_t id = 0;
  std::uint32_t anchor_signal_id = 0;
  std::uint32_t recombinant = 0;
  std::uint32_t major_parent = 0;
  std::uint32_t minor_parent = 0;
  std::size_t beginning = 0;
  std::size_t ending = 0;
  bool wraps_origin = false;
  double best_local_p_value = 1.0;
  double best_corrected_p_value = 1.0;
  std::vector<std::uint32_t> support_signal_ids;
  std::array<std::vector<std::uint32_t>, 3> role_candidates;
  std::vector<std::uint32_t> detectable_sequences;
  std::vector<TraceEvidence> trace_evidence;
  std::array<RoleHypothesisEvidence, 3> role_hypotheses;
  std::vector<std::uint32_t> automatic_co_recombinant_sequences;
  std::vector<std::uint32_t> co_recombinant_sequences;
  std::array<TreeRegionSummary, 6> tree_regions;
  std::size_t tree_panel_sequences = 0;
  bool tree_panel_subsampled = false;
  std::vector<TreePanelLeafSummary> tree_panel_leaves;
  RoleConsensusEvidence role_consensus;
  std::size_t detection_round = 0;
  std::size_t erased_nucleotide_sites = 0;
  std::size_t erased_working_sites = 0;
  std::size_t fragment_sequences_added = 0;
  bool fragment_assisted_detection = false;
  bool tract_erased_for_detection = false;
  ReviewState review_state = ReviewState::unreviewed;
  bool manual_adjusted = false;
  bool group_manual_adjusted = false;
  std::array<std::vector<std::uint32_t>, 2> adjacent_erasure_event_ids;
  std::array<std::vector<std::uint32_t>, 2> uncertain_erasure_event_ids;
  std::array<std::int32_t, 2> nearest_erasure_informative_sites{-1, -1};
  std::array<BreakpointUncertaintyEvidence, 2> breakpoint_uncertainty;
  BreakpointConfidenceEvidence breakpoint_confidence;
  MaxChiRecheckEvidence maxchi_triplet_recheck;
};

struct PlotPoint {
  std::size_t alignment_position = 0;
  std::array<double, 3> pair_identity{};
};

struct SignalPlot {
  std::uint32_t signal_id = 0;
  std::size_t window_sites = 0;
  std::vector<PlotPoint> points;
};

class RdpScanner {
 public:
  explicit RdpScanner(const Alignment& alignment);

  bool begin(ScanOptions options, std::string& error);
  int scan_batch(std::size_t triplet_budget, std::string& error);
  bool reconcile(std::string& error);
  void cancel();

  [[nodiscard]] std::string progress_json() const;
  [[nodiscard]] std::string results_json() const;
  [[nodiscard]] std::string csv() const;
  [[nodiscard]] std::string enabled_sequences_fasta() const;
  [[nodiscard]] std::string masked_or_disabled_sequences_fasta() const;
  [[nodiscard]] std::string recombinant_sequences_removed_fasta() const;
  [[nodiscard]] std::string recombinant_columns_removed_fasta() const;
  [[nodiscard]] std::string recombination_free_fasta() const;
  [[nodiscard]] std::string fragmented_fasta() const;
  [[nodiscard]] bool final_alignment_ready(std::string& error) const;
  [[nodiscard]] SignalPlot signal_plot(std::uint32_t signal_id, std::string& error) const;
  [[nodiscard]] std::string event_alignment_json(
      std::uint32_t event_id,
      std::size_t flank_sites,
      std::size_t row_limit,
      std::string& error) const;
  [[nodiscard]] std::string event_trees_json(
      std::uint32_t event_id,
      std::string& error) const;
  bool set_review_state(std::uint32_t signal_id, ReviewState state);
  bool set_event_review_state(
      std::uint32_t event_id,
      ReviewState state,
      std::string& error);
  bool update_event(
      std::uint32_t event_id,
      std::uint32_t recombinant,
      std::uint32_t major_parent,
      std::uint32_t minor_parent,
      std::size_t beginning,
      std::size_t ending,
      std::string& error);
  bool update_event_group(
      std::uint32_t event_id,
      std::vector<std::uint32_t> sequences,
      bool manual_override,
      std::string& error);
  bool reconcile_after(std::uint32_t event_id, std::string& error);
  bool restore(
      ScanOptions options,
      std::vector<Signal> signals,
      std::uint64_t correction_tests,
      std::string& error);
  bool restore_event_state(
      std::uint32_t event_id,
      std::uint32_t anchor_signal_id,
      std::uint32_t recombinant,
      std::uint32_t major_parent,
      std::uint32_t minor_parent,
      std::size_t beginning,
      std::size_t ending,
      std::size_t detection_round,
      bool tract_erased_for_detection,
      ReviewState review_state,
      bool manual_adjusted,
      std::vector<std::uint32_t> co_recombinant_sequences,
      bool group_manual_adjusted,
      std::string& error);
  bool restore_reconciliation_requirement(std::uint32_t event_id, std::string& error);

  [[nodiscard]] const ScanOptions& options() const { return options_; }
  [[nodiscard]] const std::vector<Signal>& signals() const { return signals_; }
  [[nodiscard]] const std::vector<UniqueEvent>& events() const { return events_; }
  [[nodiscard]] std::uint64_t processed_triplets() const { return processed_triplets_; }
  [[nodiscard]] std::uint64_t total_triplets() const { return total_triplets_; }
  [[nodiscard]] bool done() const { return done_; }
  [[nodiscard]] bool cancelled() const { return cancelled_.load(); }

 private:
  struct TripletProfile {
    std::array<std::uint32_t, 3> sequences{};
    std::vector<std::uint8_t> category;
    std::vector<std::size_t> coordinates;
    std::array<std::size_t, 3> category_counts{};
    std::array<double, 3> similarities{};
    std::array<std::vector<std::uint32_t>, 3> rolling_counts;
  };

  struct ErasureResult {
    std::size_t original_sites = 0;
    std::size_t working_sites = 0;
    std::size_t fragments_added = 0;
  };

  const Alignment& alignment_;
  Alignment working_alignment_;
  std::vector<std::uint8_t> native_input_missing_data_;
  std::vector<std::uint32_t> working_origins_;
  std::vector<std::int32_t> working_fragment_events_;
  ScanOptions options_;
  std::vector<std::uint32_t> active_sequences_;
  std::vector<Signal> signals_;
  std::vector<UniqueEvent> events_;
  std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> round_signal_index_;
  TripletProfile profile_scratch_;
  MaxChiWorkspace maxchi_workspace_;
  std::vector<Signal> signal_candidates_scratch_;
  std::array<std::vector<std::uint8_t>, 3> breakpoint_erasure_scratch_;
  std::vector<std::uint8_t> breakpoint_input_missing_scratch_;
  std::vector<std::uint8_t> breakpoint_polish_missing_scratch_;
  std::vector<std::size_t> breakpoint_informative_coordinates_scratch_;
  std::array<std::vector<std::size_t>, 2> breakpoint_check_coordinates_scratch_;
  std::vector<std::uint32_t> breakpoint_relevant_event_indices_scratch_;
  BurtConfidenceWorkspace breakpoint_confidence_scratch_;
  std::size_t cursor_a_ = 0;
  std::size_t cursor_b_ = 1;
  std::size_t cursor_c_ = 2;
  std::uint64_t processed_triplets_ = 0;
  std::uint64_t cumulative_triplets_ = 0;
  std::uint64_t total_triplets_ = 0;
  std::uint64_t correction_tests_ = 0;
  bool running_ = false;
  bool primary_done_ = false;
  bool done_ = false;
  std::size_t scan_round_ = 1;
  std::size_t round_signal_begin_ = 0;
  std::size_t fixed_event_count_ = 0;
  bool fragment_reentry_capped_ = false;
  std::string cycle_termination_ = "not-started";
  std::int32_t reconciliation_required_after_ = -1;
  std::atomic_bool cancelled_{false};

  [[nodiscard]] bool build_profile(
      const std::array<std::uint32_t, 3>& triplet,
      TripletProfile& profile) const;
  [[nodiscard]] bool sequence_masked(std::uint32_t sequence) const;
  [[nodiscard]] bool sequence_disabled(std::uint32_t sequence) const;
  [[nodiscard]] bool build_profile_on(
      const Alignment& alignment,
      const std::array<std::uint32_t, 3>& triplet,
      TripletProfile& profile) const;
  void scan_triplet(const std::array<std::uint32_t, 3>& triplet);
  void compute_rolling_counts(TripletProfile& profile) const;
  [[nodiscard]] std::array<std::uint8_t, 3> ranked_pairs(const TripletProfile& profile) const;
  void append_candidate_signals(
      const TripletProfile& profile,
      std::uint8_t high_pair,
      std::uint8_t candidate_pair,
      std::uint8_t low_pair,
      bool enforce_cutoff,
      std::vector<Signal>& output) const;
  [[nodiscard]] std::vector<Signal> triplet_signals(
      const std::array<std::uint32_t, 3>& triplet,
      bool enforce_cutoff,
      bool* profile_available = nullptr,
      TripletProfile* scratch = nullptr) const;
  [[nodiscard]] double tract_overlap(
      std::size_t first_beginning,
      std::size_t first_ending,
      std::size_t second_beginning,
      std::size_t second_ending) const;
  void build_event_evidence(UniqueEvent& event);
  void refresh_breakpoint_confidence(UniqueEvent& event, bool apply_polished);
  void refresh_breakpoint_context(UniqueEvent& event);
  void refresh_trace_evidence(UniqueEvent& event);
  void refresh_role_hypotheses(UniqueEvent& event);
  [[nodiscard]] MaxChiRecheckEvidence maxchi_triplet_recheck(
      const std::array<std::uint32_t, 3>& triplet);
  [[nodiscard]] bool finish_detection_round(std::string& error);
  [[nodiscard]] ErasureResult erase_event_tract(const UniqueEvent& event);
  void refresh_active_sequences();
  [[nodiscard]] std::uint64_t valid_triplet_count() const;
  [[nodiscard]] bool working_triplet_is_valid(
      const std::array<std::uint32_t, 3>& triplet) const;
  void map_signal_to_original(Signal& signal) const;
  void reset_round_cursor();
  void reset_working_alignment();
  void rebuild_working_before_event(std::size_t event_index);
  [[nodiscard]] bool event_action_allowed(
      std::uint32_t event_id,
      std::string& error) const;
  [[nodiscard]] bool matches_fixed_event(const Signal& signal) const;
  void assign_event_support(
      UniqueEvent& event,
      std::vector<std::uint8_t>& assigned,
      const std::unordered_map<std::uint64_t, std::vector<std::uint32_t>>& pair_index);
  [[nodiscard]] bool advance_triplet();
};

std::string signal_plot_json(const SignalPlot& plot);

}  // namespace rdp
