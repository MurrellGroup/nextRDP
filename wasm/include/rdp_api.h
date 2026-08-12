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
    int polish_breakpoints,
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
    int polish_breakpoints,
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
    int event_id);
int rdp_restore_scan_finish(uint32_t handle, uint32_t correction_tests);
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
