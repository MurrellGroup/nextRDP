#include "rdp_api.h"

#include "alignment.hpp"
#include "json.hpp"
#include "rdp_method.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#ifdef __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#define RDP_KEEPALIVE EMSCRIPTEN_KEEPALIVE
#else
#define RDP_KEEPALIVE
#endif

namespace {

struct Context {
  rdp::Alignment alignment;
  std::unique_ptr<rdp::RdpScanner> scanner;
  std::vector<std::string> restore_names;
  std::vector<std::string> restore_sequences;
  rdp::ScanOptions restore_options;
  std::vector<rdp::Signal> restore_signals;
  std::string cache;
  std::string error;
  bool loaded = false;
  bool restoring_alignment = false;
  bool restoring_scan = false;
};

std::vector<std::unique_ptr<Context>>& contexts() {
  static std::vector<std::unique_ptr<Context>> value;
  return value;
}

Context* context_for(std::uint32_t handle) {
  if (handle == 0 || handle > contexts().size()) return nullptr;
  return contexts()[handle - 1].get();
}

const char* invalid_handle_error() {
  static const std::string value = "The RDP analysis context is invalid.";
  return value.c_str();
}

const char* cached(Context& context, std::string value) {
  context.cache = std::move(value);
  return context.cache.c_str();
}

rdp::ReviewState review_state_from_int(int state) {
  return state == 1
      ? rdp::ReviewState::accepted
      : state == 2 ? rdp::ReviewState::rejected : rdp::ReviewState::unreviewed;
}

std::string bytes_to_string(const std::uint8_t* bytes, std::size_t length) {
  if (!bytes || length == 0) return {};
  return std::string(reinterpret_cast<const char*>(bytes), length);
}

std::string project_json(const Context& context) {
  std::ostringstream out;
  out << "{\"schema\":\"org.rdp-web.project/v1alpha9\","
         "\"engineVersion\":\"0.9.0-session-9\",\"dataset\":{";
  out << "\"format\":";
  rdp::json::string(out, context.alignment.format);
  out << ",\"alignmentLength\":" << context.alignment.length << ",\"sequences\":[";
  for (std::size_t index = 0; index < context.alignment.sequence_count(); ++index) {
    if (index) out << ',';
    out << "{\"name\":";
    rdp::json::string(out, context.alignment.names[index]);
    out << ",\"sequence\":";
    rdp::json::string(out, context.alignment.sequences[index]);
    out << '}';
  }
  out << "]},\"analysis\":";
  if (context.scanner) out << context.scanner->results_json();
  else out << "null";
  out << '}';
  return out.str();
}

}  // namespace

extern "C" {

RDP_KEEPALIVE std::uint32_t rdp_create(void) {
  auto& storage = contexts();
  for (std::size_t index = 0; index < storage.size(); ++index) {
    if (!storage[index]) {
      storage[index] = std::make_unique<Context>();
      return static_cast<std::uint32_t>(index + 1);
    }
  }
  storage.push_back(std::make_unique<Context>());
  return static_cast<std::uint32_t>(storage.size());
}

RDP_KEEPALIVE void rdp_destroy(std::uint32_t handle) {
  if (handle == 0 || handle > contexts().size()) return;
  contexts()[handle - 1].reset();
}

RDP_KEEPALIVE const char* rdp_version(void) {
  return "0.9.0-session-9";
}

RDP_KEEPALIVE int rdp_load_alignment(
    std::uint32_t handle,
    const std::uint8_t* bytes,
    std::size_t length) {
  Context* context = context_for(handle);
  if (!context) return 0;
  context->error.clear();
  context->cache.clear();
  if (!bytes || length == 0) {
    context->error = "The selected alignment file is empty.";
    return 0;
  }
  auto parsed = rdp::parse_alignment(
      std::string_view(reinterpret_cast<const char*>(bytes), length));
  if (!parsed.ok()) {
    context->error = parsed.error;
    return 0;
  }
  context->scanner.reset();
  context->restore_names.clear();
  context->restore_sequences.clear();
  context->restore_signals.clear();
  context->restoring_alignment = false;
  context->restoring_scan = false;
  context->alignment = std::move(parsed.alignment);
  context->loaded = true;
  return 1;
}

RDP_KEEPALIVE const char* rdp_get_summary_json(std::uint32_t handle) {
  Context* context = context_for(handle);
  if (!context) return invalid_handle_error();
  if (!context->loaded) {
    context->error = "No alignment is loaded.";
    return "";
  }
  return cached(*context, rdp::alignment_summary_json(context->alignment));
}

RDP_KEEPALIVE int rdp_scan_begin(
    std::uint32_t handle,
    int circular,
    int correction_mode,
    double p_value_cutoff,
    std::uint32_t window_sites,
    int polish_breakpoints,
    const std::uint8_t* masked_sequences,
    std::size_t mask_length,
    const std::uint8_t* disabled_sequences,
    std::size_t disabled_length) {
  Context* context = context_for(handle);
  if (!context || !context->loaded) return 0;
  context->error.clear();
  rdp::ScanOptions options;
  options.circular = circular != 0;
  options.correction = correction_mode == 0
      ? rdp::CorrectionMode::bonferroni
      : rdp::CorrectionMode::none;
  options.p_value_cutoff = p_value_cutoff;
  options.window_sites = window_sites;
  options.polish_breakpoints = polish_breakpoints != 0;
  options.mask.assign(context->alignment.sequence_count(), 0);
  if (masked_sequences && mask_length == options.mask.size()) {
    options.mask.assign(masked_sequences, masked_sequences + mask_length);
  }
  options.disabled.assign(context->alignment.sequence_count(), 0);
  if (disabled_sequences && disabled_length == options.disabled.size()) {
    options.disabled.assign(disabled_sequences, disabled_sequences + disabled_length);
  }

  context->scanner = std::make_unique<rdp::RdpScanner>(context->alignment);
  if (!context->scanner->begin(std::move(options), context->error)) {
    context->scanner.reset();
    return 0;
  }
  return 1;
}

RDP_KEEPALIVE int rdp_scan_batch(std::uint32_t handle, std::uint32_t triplet_budget) {
  Context* context = context_for(handle);
  if (!context || !context->scanner) return -1;
  context->error.clear();
  return context->scanner->scan_batch(triplet_budget, context->error);
}

RDP_KEEPALIVE int rdp_reconcile(std::uint32_t handle) {
  Context* context = context_for(handle);
  if (!context || !context->scanner) return 0;
  context->error.clear();
  return context->scanner->reconcile(context->error) ? 1 : 0;
}

RDP_KEEPALIVE void rdp_cancel(std::uint32_t handle) {
  Context* context = context_for(handle);
  if (context && context->scanner) context->scanner->cancel();
}

RDP_KEEPALIVE const char* rdp_get_progress_json(std::uint32_t handle) {
  Context* context = context_for(handle);
  if (!context || !context->scanner) return "";
  return cached(*context, context->scanner->progress_json());
}

RDP_KEEPALIVE const char* rdp_get_results_json(std::uint32_t handle) {
  Context* context = context_for(handle);
  if (!context || !context->scanner) return "";
  return cached(*context, context->scanner->results_json());
}

RDP_KEEPALIVE const char* rdp_get_signal_plot_json(
    std::uint32_t handle,
    std::uint32_t signal_id) {
  Context* context = context_for(handle);
  if (!context || !context->scanner) return "";
  context->error.clear();
  const auto plot = context->scanner->signal_plot(signal_id, context->error);
  if (!context->error.empty()) return "";
  return cached(*context, rdp::signal_plot_json(plot));
}

RDP_KEEPALIVE const char* rdp_get_event_alignment_json(
    std::uint32_t handle,
    std::uint32_t event_id,
    std::uint32_t flank_sites,
    std::uint32_t row_limit) {
  Context* context = context_for(handle);
  if (!context || !context->scanner) return "";
  context->error.clear();
  const std::string view = context->scanner->event_alignment_json(
      event_id, flank_sites, row_limit, context->error);
  if (!context->error.empty()) return "";
  return cached(*context, view);
}

RDP_KEEPALIVE const char* rdp_get_event_trees_json(
    std::uint32_t handle,
    std::uint32_t event_id) {
  Context* context = context_for(handle);
  if (!context || !context->scanner) return "";
  context->error.clear();
  const std::string view = context->scanner->event_trees_json(event_id, context->error);
  if (!context->error.empty()) return "";
  return cached(*context, view);
}

RDP_KEEPALIVE int rdp_set_review_state(
    std::uint32_t handle,
    std::uint32_t signal_id,
    int state) {
  Context* context = context_for(handle);
  if (!context || !context->scanner) return 0;
  context->error.clear();
  const auto review_state = review_state_from_int(state);
  if (!context->scanner->set_review_state(signal_id, review_state)) {
    context->error = "The selected RDP signal does not exist.";
    return 0;
  }
  return 1;
}

RDP_KEEPALIVE int rdp_set_event_review_state(
    std::uint32_t handle,
    std::uint32_t event_id,
    int state) {
  Context* context = context_for(handle);
  if (!context || !context->scanner) return 0;
  context->error.clear();
  return context->scanner->set_event_review_state(
             event_id, review_state_from_int(state), context->error)
      ? 1
      : 0;
}

RDP_KEEPALIVE int rdp_update_event(
    std::uint32_t handle,
    std::uint32_t event_id,
    std::uint32_t recombinant,
    std::uint32_t major_parent,
    std::uint32_t minor_parent,
    std::uint32_t beginning,
    std::uint32_t ending) {
  Context* context = context_for(handle);
  if (!context || !context->scanner) return 0;
  context->error.clear();
  return context->scanner->update_event(
             event_id,
             recombinant,
             major_parent,
             minor_parent,
             beginning,
             ending,
             context->error)
      ? 1
      : 0;
}

RDP_KEEPALIVE int rdp_update_event_group(
    std::uint32_t handle,
    std::uint32_t event_id,
    const std::uint32_t* sequence_indices,
    std::size_t sequence_count,
    int manual_override) {
  Context* context = context_for(handle);
  if (!context || !context->scanner || (!sequence_indices && sequence_count > 0)) return 0;
  context->error.clear();
  std::vector<std::uint32_t> group;
  if (sequence_count > 0) {
    group.assign(sequence_indices, sequence_indices + sequence_count);
  }
  return context->scanner->update_event_group(
             event_id, std::move(group), manual_override != 0, context->error)
      ? 1
      : 0;
}

RDP_KEEPALIVE int rdp_reconcile_after(std::uint32_t handle, std::uint32_t event_id) {
  Context* context = context_for(handle);
  if (!context || !context->scanner) return 0;
  context->error.clear();
  return context->scanner->reconcile_after(event_id, context->error) ? 1 : 0;
}

RDP_KEEPALIVE int rdp_restore_alignment_begin(
    std::uint32_t handle,
    std::uint32_t sequence_count) {
  Context* context = context_for(handle);
  if (!context || sequence_count < 3) return 0;
  context->error.clear();
  context->scanner.reset();
  context->loaded = false;
  context->restore_names.assign(sequence_count, std::string{});
  context->restore_sequences.assign(sequence_count, std::string{});
  context->restore_signals.clear();
  context->restoring_alignment = true;
  context->restoring_scan = false;
  return 1;
}

RDP_KEEPALIVE int rdp_restore_alignment_record(
    std::uint32_t handle,
    std::uint32_t index,
    const std::uint8_t* name,
    std::size_t name_length,
    const std::uint8_t* sequence,
    std::size_t sequence_length) {
  Context* context = context_for(handle);
  if (!context || !context->restoring_alignment || index >= context->restore_names.size() ||
      !sequence || sequence_length == 0) {
    return 0;
  }
  context->restore_names[index] = bytes_to_string(name, name_length);
  context->restore_sequences[index] = bytes_to_string(sequence, sequence_length);
  return 1;
}

RDP_KEEPALIVE int rdp_restore_alignment_finish(
    std::uint32_t handle,
    const std::uint8_t* format,
    std::size_t format_length) {
  Context* context = context_for(handle);
  if (!context || !context->restoring_alignment) return 0;
  context->error.clear();
  auto parsed = rdp::build_alignment(
      bytes_to_string(format, format_length),
      std::move(context->restore_names),
      std::move(context->restore_sequences));
  context->restoring_alignment = false;
  if (!parsed.ok()) {
    context->error = parsed.error;
    return 0;
  }
  context->alignment = std::move(parsed.alignment);
  context->loaded = true;
  return 1;
}

RDP_KEEPALIVE int rdp_restore_scan_begin(
    std::uint32_t handle,
    int circular,
    int correction_mode,
    double p_value_cutoff,
    std::uint32_t window_sites,
    int polish_breakpoints,
    const std::uint8_t* masked_sequences,
    std::size_t mask_length,
    const std::uint8_t* disabled_sequences,
    std::size_t disabled_length) {
  Context* context = context_for(handle);
  if (!context || !context->loaded || !masked_sequences ||
      mask_length != context->alignment.sequence_count() || !disabled_sequences ||
      disabled_length != context->alignment.sequence_count()) {
    return 0;
  }
  context->error.clear();
  context->restore_options = {};
  context->restore_options.circular = circular != 0;
  context->restore_options.correction = correction_mode == 0
      ? rdp::CorrectionMode::bonferroni
      : rdp::CorrectionMode::none;
  context->restore_options.p_value_cutoff = p_value_cutoff;
  context->restore_options.window_sites = window_sites;
  context->restore_options.polish_breakpoints = polish_breakpoints != 0;
  context->restore_options.mask.assign(masked_sequences, masked_sequences + mask_length);
  context->restore_options.disabled.assign(
      disabled_sequences, disabled_sequences + disabled_length);
  context->restore_signals.clear();
  context->restoring_scan = true;
  return 1;
}

RDP_KEEPALIVE int rdp_restore_signal(
    std::uint32_t handle,
    std::uint32_t triplet_0,
    std::uint32_t triplet_1,
    std::uint32_t triplet_2,
    std::uint32_t recombinant,
    std::uint32_t major_parent,
    std::uint32_t minor_parent,
    std::uint32_t beginning,
    std::uint32_t ending,
    int wraps_origin,
    std::uint32_t informative_beginning,
    std::uint32_t informative_ending,
    double local_p_value,
    double corrected_p_value,
    std::uint32_t correction_tests,
    double pair_similarity_0,
    double pair_similarity_1,
    double pair_similarity_2,
    std::uint32_t informative_sites,
    std::uint32_t candidate_pair,
    int fragment_assisted,
    int fragment_event_0,
    int fragment_event_1,
    int fragment_event_2,
    int review_state,
    int event_id) {
  Context* context = context_for(handle);
  if (!context || !context->restoring_scan) return 0;
  rdp::Signal signal;
  signal.id = static_cast<std::uint32_t>(context->restore_signals.size());
  signal.triplet = {triplet_0, triplet_1, triplet_2};
  signal.recombinant = recombinant;
  signal.major_parent = major_parent;
  signal.minor_parent = minor_parent;
  signal.beginning = beginning;
  signal.ending = ending;
  signal.wraps_origin = wraps_origin != 0;
  signal.informative_beginning = informative_beginning;
  signal.informative_ending = informative_ending;
  signal.local_p_value = local_p_value;
  signal.corrected_p_value = corrected_p_value;
  signal.correction_tests = correction_tests;
  signal.pair_similarity = {pair_similarity_0, pair_similarity_1, pair_similarity_2};
  signal.informative_sites = informative_sites;
  signal.candidate_pair = static_cast<std::uint8_t>(candidate_pair);
  signal.fragment_assisted = fragment_assisted != 0;
  signal.fragment_event_context = {fragment_event_0, fragment_event_1, fragment_event_2};
  signal.review_state = review_state_from_int(review_state);
  signal.event_id = event_id;
  context->restore_signals.push_back(signal);
  return 1;
}

RDP_KEEPALIVE int rdp_restore_scan_finish(
    std::uint32_t handle,
    std::uint32_t correction_tests) {
  Context* context = context_for(handle);
  if (!context || !context->loaded || !context->restoring_scan) return 0;
  context->error.clear();
  context->scanner = std::make_unique<rdp::RdpScanner>(context->alignment);
  const bool restored = context->scanner->restore(
      std::move(context->restore_options),
      std::move(context->restore_signals),
      correction_tests,
      context->error);
  context->restoring_scan = false;
  if (!restored) {
    context->scanner.reset();
    return 0;
  }
  return 1;
}

RDP_KEEPALIVE int rdp_restore_event_state(
    std::uint32_t handle,
    std::uint32_t event_id,
    std::uint32_t anchor_signal_id,
    std::uint32_t recombinant,
    std::uint32_t major_parent,
    std::uint32_t minor_parent,
    std::uint32_t beginning,
    std::uint32_t ending,
    std::uint32_t detection_round,
    int tract_erased_for_detection,
    int review_state,
    int manual_adjusted,
    const std::uint32_t* co_recombinant_sequences,
    std::size_t co_recombinant_sequence_count,
    int group_manual_adjusted) {
  Context* context = context_for(handle);
  if (!context || !context->scanner ||
      (!co_recombinant_sequences && co_recombinant_sequence_count > 0)) {
    return 0;
  }
  context->error.clear();
  std::vector<std::uint32_t> group;
  if (co_recombinant_sequence_count > 0) {
    group.assign(
        co_recombinant_sequences,
        co_recombinant_sequences + co_recombinant_sequence_count);
  }
  return context->scanner->restore_event_state(
             event_id,
             anchor_signal_id,
             recombinant,
             major_parent,
             minor_parent,
             beginning,
             ending,
             detection_round,
             tract_erased_for_detection != 0,
             review_state_from_int(review_state),
             manual_adjusted != 0,
             std::move(group),
             group_manual_adjusted != 0,
             context->error)
      ? 1
      : 0;
}

RDP_KEEPALIVE int rdp_restore_reconciliation_required_after(
    std::uint32_t handle,
    std::uint32_t event_id) {
  Context* context = context_for(handle);
  if (!context || !context->scanner) return 0;
  context->error.clear();
  return context->scanner->restore_reconciliation_requirement(event_id, context->error) ? 1 : 0;
}

RDP_KEEPALIVE const char* rdp_export_csv(std::uint32_t handle) {
  Context* context = context_for(handle);
  if (!context || !context->scanner) return "";
  return cached(*context, context->scanner->csv());
}

RDP_KEEPALIVE const char* rdp_export_enabled_sequences_fasta(
    std::uint32_t handle,
    const std::uint8_t* masked_sequences,
    std::size_t mask_length,
    const std::uint8_t* disabled_sequences,
    std::size_t disabled_length) {
  Context* context = context_for(handle);
  if (!context || !context->loaded) return "";
  context->error.clear();
  if (!masked_sequences || !disabled_sequences ||
      mask_length != context->alignment.sequence_count() ||
      disabled_length != context->alignment.sequence_count()) {
    context->error = "The sequence curation state does not match the loaded alignment.";
    return "";
  }
  std::vector<std::uint8_t> mask(masked_sequences, masked_sequences + mask_length);
  std::vector<std::uint8_t> disabled(
      disabled_sequences, disabled_sequences + disabled_length);
  for (std::size_t sequence = 0; sequence < mask.size(); ++sequence) {
    if (disabled[sequence] != 0) mask[sequence] = 0;
  }
  std::string fasta = rdp::curated_sequences_fasta(
      context->alignment, mask, disabled, true);
  if (fasta.empty()) {
    context->error =
        "No enabled, unmasked sequence is available for the enabled-only alignment.";
    return "";
  }
  return cached(*context, std::move(fasta));
}

RDP_KEEPALIVE const char* rdp_export_masked_or_disabled_sequences_fasta(
    std::uint32_t handle,
    const std::uint8_t* masked_sequences,
    std::size_t mask_length,
    const std::uint8_t* disabled_sequences,
    std::size_t disabled_length) {
  Context* context = context_for(handle);
  if (!context || !context->loaded) return "";
  context->error.clear();
  if (!masked_sequences || !disabled_sequences ||
      mask_length != context->alignment.sequence_count() ||
      disabled_length != context->alignment.sequence_count()) {
    context->error = "The sequence curation state does not match the loaded alignment.";
    return "";
  }
  std::vector<std::uint8_t> mask(masked_sequences, masked_sequences + mask_length);
  std::vector<std::uint8_t> disabled(
      disabled_sequences, disabled_sequences + disabled_length);
  for (std::size_t sequence = 0; sequence < mask.size(); ++sequence) {
    if (disabled[sequence] != 0) mask[sequence] = 0;
  }
  std::string fasta = rdp::curated_sequences_fasta(
      context->alignment, mask, disabled, false);
  if (fasta.empty()) {
    context->error =
        "No masked or disabled sequence is available for the excluded-row alignment.";
    return "";
  }
  return cached(*context, std::move(fasta));
}

RDP_KEEPALIVE const char* rdp_export_recombinant_sequences_removed_fasta(
    std::uint32_t handle) {
  Context* context = context_for(handle);
  if (!context || !context->scanner) return "";
  context->error.clear();
  if (!context->scanner->final_alignment_ready(context->error)) return "";
  std::string fasta = context->scanner->recombinant_sequences_removed_fasta();
  if (fasta.empty()) {
    context->error =
        "Every sequence belongs to an accepted recombinant group; the sequence-removed alignment would be empty.";
    return "";
  }
  return cached(*context, std::move(fasta));
}

RDP_KEEPALIVE const char* rdp_export_recombinant_columns_removed_fasta(
    std::uint32_t handle) {
  Context* context = context_for(handle);
  if (!context || !context->scanner) return "";
  context->error.clear();
  if (!context->scanner->final_alignment_ready(context->error)) return "";
  std::string fasta = context->scanner->recombinant_columns_removed_fasta();
  if (fasta.empty()) {
    context->error =
        "Accepted event tracts cover every alignment column; the column-removed alignment would be empty.";
    return "";
  }
  return cached(*context, std::move(fasta));
}

RDP_KEEPALIVE const char* rdp_export_recombination_free_fasta(std::uint32_t handle) {
  Context* context = context_for(handle);
  if (!context || !context->scanner) return "";
  context->error.clear();
  if (!context->scanner->final_alignment_ready(context->error)) return "";
  return cached(*context, context->scanner->recombination_free_fasta());
}

RDP_KEEPALIVE const char* rdp_export_fragmented_fasta(std::uint32_t handle) {
  Context* context = context_for(handle);
  if (!context || !context->scanner) return "";
  context->error.clear();
  if (!context->scanner->final_alignment_ready(context->error)) return "";
  return cached(*context, context->scanner->fragmented_fasta());
}

RDP_KEEPALIVE const char* rdp_export_project_json(std::uint32_t handle) {
  Context* context = context_for(handle);
  if (!context || !context->loaded) return "";
  return cached(*context, project_json(*context));
}

RDP_KEEPALIVE const char* rdp_get_error(std::uint32_t handle) {
  Context* context = context_for(handle);
  return context ? context->error.c_str() : invalid_handle_error();
}

}  // extern "C"
