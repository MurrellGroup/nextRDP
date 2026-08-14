#pragma once

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

uint32_t rdp_create(void);
void rdp_destroy(uint32_t handle);
const char* rdp_version(void);

int rdp_load_alignment(uint32_t handle, const uint8_t* bytes, size_t length);
const char* rdp_get_summary_json(uint32_t handle);

int rdp_scan_begin(
    uint32_t handle,
    int circular,
    int correction_mode,
    double p_value_cutoff,
    uint32_t window_sites,
    int maxchi_enabled,
    uint32_t maxchi_window_sites,
    int chimaera_enabled,
    uint32_t chimaera_window_sites,
    int geneconv_enabled,
    uint32_t geneconv_mismatch_scale,
    uint32_t geneconv_max_overlaps,
    int threeseq_enabled,
    int bootscan_primary_enabled,
    int bootscan_secondary_enabled,
    uint32_t bootscan_window_sites,
    uint32_t bootscan_step_sites,
    uint32_t bootscan_bootstrap_replicates,
    double bootscan_support_cutoff,
    uint32_t bootscan_random_seed,
    int siscan_primary_enabled,
    int siscan_secondary_enabled,
    uint32_t siscan_window_sites,
    uint32_t siscan_step_sites,
    uint32_t siscan_scan_permutations,
    uint32_t siscan_p_value_permutations,
    uint32_t siscan_random_seed,
    int polish_breakpoints,
    int query_reference_mode,
    const uint32_t* reference_groups,
    size_t reference_group_count,
    const uint8_t* masked_sequences,
    size_t mask_length,
    const uint8_t* disabled_sequences,
    size_t disabled_length);
int rdp_scan_batch(uint32_t handle, uint32_t triplet_budget);
int rdp_reconcile(uint32_t handle);
void rdp_cancel(uint32_t handle);

const char* rdp_get_progress_json(uint32_t handle);
const char* rdp_get_results_json(uint32_t handle);
const char* rdp_get_signal_plot_json(uint32_t handle, uint32_t signal_id);
const char* rdp_get_event_alignment_json(
    uint32_t handle,
    uint32_t event_id,
    uint32_t flank_sites,
    uint32_t row_limit);
const char* rdp_get_event_trees_json(uint32_t handle, uint32_t event_id);
const char* rdp_get_event_phylpro_json(
    uint32_t handle,
    uint32_t event_id,
    uint32_t window_sites,
    int gap_mode,
    int include_self);
int rdp_set_review_state(uint32_t handle, uint32_t signal_id, int state);
int rdp_set_event_review_state(uint32_t handle, uint32_t event_id, int state);
int rdp_update_event(
    uint32_t handle,
    uint32_t event_id,
    uint32_t recombinant,
    uint32_t major_parent,
    uint32_t minor_parent,
    uint32_t beginning,
    uint32_t ending);
int rdp_update_event_group(
    uint32_t handle,
    uint32_t event_id,
    const uint32_t* sequence_indices,
    size_t sequence_count,
    int manual_override);
int rdp_reconcile_after(uint32_t handle, uint32_t event_id);

int rdp_restore_alignment_begin(uint32_t handle, uint32_t sequence_count);
int rdp_restore_alignment_record(
    uint32_t handle,
    uint32_t index,
    const uint8_t* name,
    size_t name_length,
    const uint8_t* sequence,
    size_t sequence_length);
int rdp_restore_alignment_finish(
    uint32_t handle,
    const uint8_t* format,
    size_t format_length);
int rdp_restore_scan_begin(
    uint32_t handle,
    int circular,
    int correction_mode,
    double p_value_cutoff,
    uint32_t window_sites,
    int maxchi_enabled,
    uint32_t maxchi_window_sites,
    int chimaera_enabled,
    uint32_t chimaera_window_sites,
    int geneconv_enabled,
    uint32_t geneconv_mismatch_scale,
    uint32_t geneconv_max_overlaps,
    int threeseq_enabled,
    int bootscan_primary_enabled,
    int bootscan_secondary_enabled,
    uint32_t bootscan_window_sites,
    uint32_t bootscan_step_sites,
    uint32_t bootscan_bootstrap_replicates,
    double bootscan_support_cutoff,
    uint32_t bootscan_random_seed,
    int siscan_primary_enabled,
    int siscan_secondary_enabled,
    uint32_t siscan_window_sites,
    uint32_t siscan_step_sites,
    uint32_t siscan_scan_permutations,
    uint32_t siscan_p_value_permutations,
    uint32_t siscan_random_seed,
    int polish_breakpoints,
    int query_reference_mode,
    const uint32_t* reference_groups,
    size_t reference_group_count,
    const uint8_t* masked_sequences,
    size_t mask_length,
    const uint8_t* disabled_sequences,
    size_t disabled_length);
int rdp_restore_signal(
    uint32_t handle,
    uint32_t triplet_0,
    uint32_t triplet_1,
    uint32_t triplet_2,
    uint32_t recombinant,
    uint32_t major_parent,
    uint32_t minor_parent,
    uint32_t beginning,
    uint32_t ending,
    int wraps_origin,
    uint32_t informative_beginning,
    uint32_t informative_ending,
    double local_p_value,
    double corrected_p_value,
    uint32_t correction_tests,
    double pair_similarity_0,
    double pair_similarity_1,
    double pair_similarity_2,
    uint32_t informative_sites,
    uint32_t candidate_pair,
    int fragment_assisted,
    int fragment_event_0,
    int fragment_event_1,
    int fragment_event_2,
    int review_state,
    int event_id,
    int method);
int rdp_restore_maxchi_discovery(
    uint32_t handle,
    uint32_t signal_id,
    int peak_pair,
    int tract_side,
    uint32_t peak_attempt,
    uint32_t peak_alignment_position,
    uint32_t variable_sites,
    uint32_t initial_half_window,
    uint32_t grown_half_window,
    uint32_t critical_difference,
    double maximum_chi_square,
    double raw_p_value,
    double within_triplet_p_value,
    double left_flank_chi_square,
    double right_flank_chi_square,
    int missing_data_window_filter_applied,
    int linear_edge_window_filter_applied);
int rdp_restore_chimaera_discovery(
    uint32_t handle,
    uint32_t signal_id,
    uint32_t target_local,
    int tract_side,
    uint32_t peak_attempt,
    uint32_t peak_alignment_position,
    uint32_t information_rich_sites,
    uint32_t initial_half_window,
    uint32_t grown_half_window,
    uint32_t critical_difference,
    double maximum_chi_square,
    double raw_p_value,
    double within_triplet_p_value,
    double left_flank_chi_square,
    double right_flank_chi_square,
    double inside_parent_one_match_rate,
    double outside_parent_one_match_rate,
    int missing_data_window_filter_applied,
    int linear_edge_window_filter_applied);
int rdp_restore_geneconv_discovery(
    uint32_t handle,
    uint32_t signal_id,
    uint32_t track,
    uint32_t polymorphic_sites,
    uint32_t positive_sites,
    uint32_t discordant_sites,
    uint32_t mismatch_penalty,
    uint32_t fragment_score,
    uint32_t critical_score,
    double lambda,
    double karlin_altschul_k,
    double raw_p_value);
int rdp_restore_threeseq_discovery(
    uint32_t handle,
    uint32_t signal_id,
    uint32_t target_local,
    int walk_direction,
    uint32_t information_rich_sites,
    uint32_t parent_one_matches,
    uint32_t parent_two_matches,
    uint32_t probability_excursion,
    uint32_t maximum_excursion,
    double raw_p_value,
    int exact_probability,
    int siegmund_fallback,
    int missing_data_split_applied);
int rdp_restore_bootscan_discovery(
    uint32_t handle,
    uint32_t signal_id,
    uint32_t supported_pair,
    uint32_t windows_scored,
    uint32_t usable_windows,
    uint32_t informative_sites,
    uint32_t tract_informative_sites,
    uint32_t tract_pair_matches,
    uint32_t outside_pair_matches,
    double maximum_pair_support,
    double mean_pair_support,
    double bootstrap_p_value,
    double raw_p_value,
    int erased_window_filter_applied);
int rdp_restore_siscan_discovery(
    uint32_t handle,
    uint32_t signal_id,
    uint32_t global_pair,
    uint32_t candidate_pair,
    uint32_t outlier_sequence,
    uint32_t windows_in_region,
    uint32_t informative_sites,
    double permutation_draws,
    uint32_t selected_score,
    uint32_t selected_score_family,
    double maximum_z,
    double normal_tail_p_value,
    double region_length_adjusted_p_value,
    double window_adjusted_p_value);
int rdp_restore_scan_finish(
    uint32_t handle,
    uint32_t correction_tests,
    double cumulative_triplets,
    uint32_t scan_rounds,
    double maxchi_profiles_scanned,
    double maxchi_peak_attempts,
    double maxchi_candidates_found,
    double maxchi_peak_limit_triplets,
    double chimaera_profiles_scanned,
    double chimaera_peak_attempts,
    double chimaera_candidates_found,
    double chimaera_peak_limit_targets,
    double geneconv_fragments_scored,
    double geneconv_qualified_fragments,
    double geneconv_candidates_found,
    double geneconv_overlap_rejections,
    double geneconv_numerical_fallback_tracks,
    double threeseq_profiles_scanned,
    double threeseq_exact_evaluations,
    double threeseq_approximate_evaluations,
    double threeseq_candidates_found,
    double bootscan_profiles_scanned,
    double bootscan_candidate_regions_scored,
    double bootscan_candidates_found,
    double bootscan_pair_profiles_requested,
    double bootscan_pair_profile_cache_hits,
    double bootscan_pair_profile_cache_misses,
    double bootscan_pair_profile_cache_evictions,
    double bootscan_pair_profile_cache_peak_bytes,
    double siscan_profiles_scanned,
    double siscan_windows_scored,
    double siscan_candidate_regions_scored,
    double siscan_candidates_found,
    double siscan_permutation_draws,
    double siscan_context_builds,
    double siscan_context_pair_comparisons,
    double siscan_context_tree_merges,
    double siscan_random_values_generated,
    const uint8_t* cycle_termination,
    size_t cycle_termination_length);
int rdp_restore_event_state(
    uint32_t handle,
    uint32_t event_id,
    uint32_t anchor_signal_id,
    uint32_t recombinant,
    uint32_t major_parent,
    uint32_t minor_parent,
    uint32_t beginning,
    uint32_t ending,
    uint32_t detection_round,
    int tract_erased_for_detection,
    int review_state,
    int manual_adjusted,
    const uint32_t* co_recombinant_sequences,
    size_t co_recombinant_sequence_count,
    int group_manual_adjusted);
int rdp_restore_reconciliation_required_after(uint32_t handle, uint32_t event_id);

const char* rdp_export_csv(uint32_t handle);
const char* rdp_export_enabled_sequences_fasta(
    uint32_t handle,
    const uint8_t* masked_sequences,
    size_t mask_length,
    const uint8_t* disabled_sequences,
    size_t disabled_length);
const char* rdp_export_masked_or_disabled_sequences_fasta(
    uint32_t handle,
    const uint8_t* masked_sequences,
    size_t mask_length,
    const uint8_t* disabled_sequences,
    size_t disabled_length);
const char* rdp_export_recombinant_sequences_removed_fasta(uint32_t handle);
const char* rdp_export_recombinant_columns_removed_fasta(uint32_t handle);
const char* rdp_export_recombination_free_fasta(uint32_t handle);
const char* rdp_export_fragmented_fasta(uint32_t handle);
const char* rdp_export_project_json(uint32_t handle);
const char* rdp_get_error(uint32_t handle);

#ifdef __cplusplus
}
#endif
