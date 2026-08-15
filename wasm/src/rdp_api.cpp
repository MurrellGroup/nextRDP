#include "rdp_api.h"

#include "alignment.hpp"
#include "json.hpp"
#include "rdp_method.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
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
  std::size_t requested_worker_threads = 1;
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

std::uint64_t restored_counter(double value) {
  if (!std::isfinite(value) || value <= 0.0) return 0;
  const double maximum = static_cast<double>(
      std::numeric_limits<std::uint64_t>::max());
  if (value >= maximum) return std::numeric_limits<std::uint64_t>::max();
  return static_cast<std::uint64_t>(value);
}

std::string project_json(const Context& context) {
  std::ostringstream out;
  out << "{\"schema\":\"org.rdp-web.project/v1alpha19\","
         "\"engineVersion\":\"0.26.0-session-26\",\"dataset\":{";
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
  return "0.26.0-session-26";
}

RDP_KEEPALIVE std::uint32_t rdp_set_worker_threads(
    std::uint32_t handle,
    std::uint32_t requested) {
  Context* context = context_for(handle);
  if (!context) return 0;
  context->requested_worker_threads =
      std::clamp<std::size_t>(requested, 1, 6);
  if (context->scanner) {
    context->scanner->set_worker_threads(context->requested_worker_threads);
    return static_cast<std::uint32_t>(context->scanner->worker_threads());
  }
#if defined(RDP_ENABLE_THREADS)
  return static_cast<std::uint32_t>(context->requested_worker_threads);
#else
  return 1;
#endif
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
    int maxchi_enabled,
    std::uint32_t maxchi_window_sites,
    int chimaera_enabled,
    std::uint32_t chimaera_window_sites,
    int geneconv_enabled,
    std::uint32_t geneconv_mismatch_scale,
    std::uint32_t geneconv_max_overlaps,
    int threeseq_enabled,
    int bootscan_primary_enabled,
    int bootscan_secondary_enabled,
    std::uint32_t bootscan_window_sites,
    std::uint32_t bootscan_step_sites,
    std::uint32_t bootscan_bootstrap_replicates,
    double bootscan_support_cutoff,
    std::uint32_t bootscan_random_seed,
    int siscan_primary_enabled,
    int siscan_secondary_enabled,
    std::uint32_t siscan_window_sites,
    std::uint32_t siscan_step_sites,
    std::uint32_t siscan_scan_permutations,
    std::uint32_t siscan_p_value_permutations,
    std::uint32_t siscan_random_seed,
    int polish_breakpoints,
    int query_reference_mode,
    const std::uint32_t* reference_groups,
    std::size_t reference_group_count,
    const std::uint8_t* masked_sequences,
    std::size_t mask_length,
    const std::uint8_t* disabled_sequences,
    std::size_t disabled_length) {
  Context* context = context_for(handle);
  if (!context || !context->loaded) return 0;
  context->error.clear();
  const std::size_t sequence_count = context->alignment.sequence_count();
  if (!reference_groups || reference_group_count != sequence_count ||
      !masked_sequences || mask_length != sequence_count ||
      !disabled_sequences || disabled_length != sequence_count) {
    context->error = "The scan's sequence-role buffers do not match the loaded alignment.";
    return 0;
  }
  rdp::ScanOptions options;
  options.circular = circular != 0;
  options.correction = correction_mode == 0
      ? rdp::CorrectionMode::bonferroni
      : rdp::CorrectionMode::none;
  options.p_value_cutoff = p_value_cutoff;
  options.window_sites = window_sites;
  options.maxchi_enabled = maxchi_enabled != 0;
  options.maxchi_window_sites = maxchi_window_sites;
  options.chimaera_enabled = chimaera_enabled != 0;
  options.chimaera_window_sites = chimaera_window_sites;
  options.geneconv_enabled = geneconv_enabled != 0;
  options.geneconv_mismatch_scale = geneconv_mismatch_scale;
  options.geneconv_max_overlaps = geneconv_max_overlaps;
  options.threeseq_enabled = threeseq_enabled != 0;
  options.bootscan_primary_enabled = bootscan_primary_enabled != 0;
  options.bootscan_secondary_enabled = bootscan_secondary_enabled != 0;
  options.bootscan_window_sites = bootscan_window_sites;
  options.bootscan_step_sites = bootscan_step_sites;
  options.bootscan_bootstrap_replicates = bootscan_bootstrap_replicates;
  options.bootscan_support_cutoff = bootscan_support_cutoff;
  options.bootscan_random_seed = bootscan_random_seed;
  options.siscan_primary_enabled = siscan_primary_enabled != 0;
  options.siscan_secondary_enabled = siscan_secondary_enabled != 0;
  options.siscan_window_sites = siscan_window_sites;
  options.siscan_step_sites = siscan_step_sites;
  options.siscan_scan_permutations = siscan_scan_permutations;
  options.siscan_p_value_permutations = siscan_p_value_permutations;
  options.siscan_random_seed = siscan_random_seed;
  options.polish_breakpoints = polish_breakpoints != 0;
  options.analysis_mode = query_reference_mode != 0
      ? rdp::AnalysisMode::query_reference
      : rdp::AnalysisMode::exploratory;
  options.reference_groups.assign(
      reference_groups, reference_groups + reference_group_count);
  options.mask.assign(masked_sequences, masked_sequences + mask_length);
  options.disabled.assign(disabled_sequences, disabled_sequences + disabled_length);

  context->scanner = std::make_unique<rdp::RdpScanner>(context->alignment);
  context->scanner->set_worker_threads(context->requested_worker_threads);
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

RDP_KEEPALIVE const char* rdp_get_event_phylpro_json(
    std::uint32_t handle,
    std::uint32_t event_id,
    std::uint32_t window_sites,
    int gap_mode,
    int include_self) {
  Context* context = context_for(handle);
  if (!context || !context->scanner) return "";
  context->error.clear();
  if (gap_mode != 0 && gap_mode != 1) {
    context->error = "The selected PHYLPRO gap mode is unavailable.";
    return "";
  }
  const std::string view = context->scanner->event_phylpro_json(
      event_id,
      window_sites,
      gap_mode == 1
          ? rdp::PhylproGapMode::strip_columns
          : rdp::PhylproGapMode::ignore_pairwise,
      include_self != 0,
      context->error);
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
    context->error = "The selected recombination signal does not exist.";
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
    int maxchi_enabled,
    std::uint32_t maxchi_window_sites,
    int chimaera_enabled,
    std::uint32_t chimaera_window_sites,
    int geneconv_enabled,
    std::uint32_t geneconv_mismatch_scale,
    std::uint32_t geneconv_max_overlaps,
    int threeseq_enabled,
    int bootscan_primary_enabled,
    int bootscan_secondary_enabled,
    std::uint32_t bootscan_window_sites,
    std::uint32_t bootscan_step_sites,
    std::uint32_t bootscan_bootstrap_replicates,
    double bootscan_support_cutoff,
    std::uint32_t bootscan_random_seed,
    int siscan_primary_enabled,
    int siscan_secondary_enabled,
    std::uint32_t siscan_window_sites,
    std::uint32_t siscan_step_sites,
    std::uint32_t siscan_scan_permutations,
    std::uint32_t siscan_p_value_permutations,
    std::uint32_t siscan_random_seed,
    int polish_breakpoints,
    int query_reference_mode,
    const std::uint32_t* reference_groups,
    std::size_t reference_group_count,
    const std::uint8_t* masked_sequences,
    std::size_t mask_length,
    const std::uint8_t* disabled_sequences,
    std::size_t disabled_length) {
  Context* context = context_for(handle);
  if (!context || !context->loaded) {
    return 0;
  }
  context->error.clear();
  const std::size_t sequence_count = context->alignment.sequence_count();
  if (!reference_groups || reference_group_count != sequence_count || !masked_sequences ||
      mask_length != sequence_count || !disabled_sequences ||
      disabled_length != sequence_count) {
    context->error = "The saved scan's sequence-role buffers do not match its alignment.";
    return 0;
  }
  context->restore_options = {};
  context->restore_options.circular = circular != 0;
  context->restore_options.correction = correction_mode == 0
      ? rdp::CorrectionMode::bonferroni
      : rdp::CorrectionMode::none;
  context->restore_options.p_value_cutoff = p_value_cutoff;
  context->restore_options.window_sites = window_sites;
  context->restore_options.maxchi_enabled = maxchi_enabled != 0;
  context->restore_options.maxchi_window_sites = maxchi_window_sites;
  context->restore_options.chimaera_enabled = chimaera_enabled != 0;
  context->restore_options.chimaera_window_sites = chimaera_window_sites;
  context->restore_options.geneconv_enabled = geneconv_enabled != 0;
  context->restore_options.geneconv_mismatch_scale = geneconv_mismatch_scale;
  context->restore_options.geneconv_max_overlaps = geneconv_max_overlaps;
  context->restore_options.threeseq_enabled = threeseq_enabled != 0;
  context->restore_options.bootscan_primary_enabled =
      bootscan_primary_enabled != 0;
  context->restore_options.bootscan_secondary_enabled =
      bootscan_secondary_enabled != 0;
  context->restore_options.bootscan_window_sites = bootscan_window_sites;
  context->restore_options.bootscan_step_sites = bootscan_step_sites;
  context->restore_options.bootscan_bootstrap_replicates =
      bootscan_bootstrap_replicates;
  context->restore_options.bootscan_support_cutoff = bootscan_support_cutoff;
  context->restore_options.bootscan_random_seed = bootscan_random_seed;
  context->restore_options.siscan_primary_enabled =
      siscan_primary_enabled != 0;
  context->restore_options.siscan_secondary_enabled =
      siscan_secondary_enabled != 0;
  context->restore_options.siscan_window_sites = siscan_window_sites;
  context->restore_options.siscan_step_sites = siscan_step_sites;
  context->restore_options.siscan_scan_permutations =
      siscan_scan_permutations;
  context->restore_options.siscan_p_value_permutations =
      siscan_p_value_permutations;
  context->restore_options.siscan_random_seed = siscan_random_seed;
  context->restore_options.polish_breakpoints = polish_breakpoints != 0;
  context->restore_options.analysis_mode = query_reference_mode != 0
      ? rdp::AnalysisMode::query_reference
      : rdp::AnalysisMode::exploratory;
  context->restore_options.reference_groups.assign(
      reference_groups, reference_groups + reference_group_count);
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
    int event_id,
    int method) {
  Context* context = context_for(handle);
  if (!context || !context->restoring_scan) return 0;
  rdp::Signal signal;
  signal.id = static_cast<std::uint32_t>(context->restore_signals.size());
  signal.method = method == 1
      ? rdp::SignalMethod::maxchi
      : method == 2
          ? rdp::SignalMethod::chimaera
          : method == 3
              ? rdp::SignalMethod::geneconv
              : method == 4
                  ? rdp::SignalMethod::threeseq
                  : method == 5
                      ? rdp::SignalMethod::bootscan
                  : method == 6
                      ? rdp::SignalMethod::siscan
                  : rdp::SignalMethod::rdp;
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

RDP_KEEPALIVE int rdp_restore_maxchi_discovery(
    std::uint32_t handle,
    std::uint32_t signal_id,
    int peak_pair,
    int tract_side,
    std::uint32_t peak_attempt,
    std::uint32_t peak_alignment_position,
    std::uint32_t variable_sites,
    std::uint32_t initial_half_window,
    std::uint32_t grown_half_window,
    std::uint32_t critical_difference,
    double maximum_chi_square,
    double raw_p_value,
    double within_triplet_p_value,
    double left_flank_chi_square,
    double right_flank_chi_square,
    int missing_data_window_filter_applied,
    int linear_edge_window_filter_applied) {
  Context* context = context_for(handle);
  if (!context || !context->restoring_scan ||
      signal_id >= context->restore_signals.size()) {
    return 0;
  }
  auto& signal = context->restore_signals[signal_id];
  if (signal.method != rdp::SignalMethod::maxchi) return 0;
  auto& discovery = signal.maxchi_discovery;
  discovery.beginning = signal.beginning;
  discovery.ending = signal.ending;
  discovery.wraps_origin = signal.wraps_origin;
  discovery.informative_beginning = signal.informative_beginning;
  discovery.informative_ending = signal.informative_ending;
  const auto local_member = [&](std::uint32_t sequence) {
    for (std::uint8_t member = 0; member < signal.triplet.size(); ++member) {
      if (signal.triplet[member] == sequence) return member;
    }
    return static_cast<std::uint8_t>(0);
  };
  discovery.recombinant_local = local_member(signal.recombinant);
  discovery.major_parent_local = local_member(signal.major_parent);
  discovery.minor_parent_local = local_member(signal.minor_parent);
  discovery.candidate_pair = signal.candidate_pair;
  discovery.peak_pair = static_cast<std::int8_t>(peak_pair);
  discovery.tract_side = tract_side < 0
      ? rdp::MaxChiTractSide::left
      : tract_side > 0
          ? rdp::MaxChiTractSide::right
          : rdp::MaxChiTractSide::unavailable;
  discovery.peak_attempt = peak_attempt;
  discovery.peak_alignment_position = peak_alignment_position;
  discovery.variable_sites = variable_sites;
  discovery.initial_half_window = initial_half_window;
  discovery.grown_half_window = grown_half_window;
  discovery.critical_difference = critical_difference;
  discovery.maximum_chi_square = maximum_chi_square;
  discovery.raw_p_value = raw_p_value;
  discovery.within_triplet_p_value = within_triplet_p_value;
  discovery.corrected_p_value = signal.corrected_p_value;
  discovery.left_flank_chi_square = left_flank_chi_square;
  discovery.right_flank_chi_square = right_flank_chi_square;
  discovery.pair_similarity = signal.pair_similarity;
  discovery.missing_data_window_filter_applied =
      missing_data_window_filter_applied != 0;
  discovery.linear_edge_window_filter_applied =
      linear_edge_window_filter_applied != 0;
  return 1;
}

RDP_KEEPALIVE int rdp_restore_chimaera_discovery(
    std::uint32_t handle,
    std::uint32_t signal_id,
    std::uint32_t target_local,
    int tract_side,
    std::uint32_t peak_attempt,
    std::uint32_t peak_alignment_position,
    std::uint32_t information_rich_sites,
    std::uint32_t initial_half_window,
    std::uint32_t grown_half_window,
    std::uint32_t critical_difference,
    double maximum_chi_square,
    double raw_p_value,
    double within_triplet_p_value,
    double left_flank_chi_square,
    double right_flank_chi_square,
    double inside_parent_one_match_rate,
    double outside_parent_one_match_rate,
    int missing_data_window_filter_applied,
    int linear_edge_window_filter_applied) {
  Context* context = context_for(handle);
  if (!context || !context->restoring_scan ||
      signal_id >= context->restore_signals.size() || target_local > 2) {
    return 0;
  }
  auto& signal = context->restore_signals[signal_id];
  if (signal.method != rdp::SignalMethod::chimaera) return 0;
  auto& discovery = signal.chimaera_discovery;
  discovery.beginning = signal.beginning;
  discovery.ending = signal.ending;
  discovery.wraps_origin = signal.wraps_origin;
  discovery.informative_beginning = signal.informative_beginning;
  discovery.informative_ending = signal.informative_ending;
  const auto local_member = [&](std::uint32_t sequence) {
    for (std::uint8_t member = 0; member < signal.triplet.size(); ++member) {
      if (signal.triplet[member] == sequence) return member;
    }
    return static_cast<std::uint8_t>(0);
  };
  discovery.target_local = static_cast<std::uint8_t>(target_local);
  discovery.recombinant_local = local_member(signal.recombinant);
  discovery.major_parent_local = local_member(signal.major_parent);
  discovery.minor_parent_local = local_member(signal.minor_parent);
  discovery.candidate_pair = signal.candidate_pair;
  discovery.tract_side = tract_side < 0
      ? rdp::MaxChiTractSide::left
      : tract_side > 0
          ? rdp::MaxChiTractSide::right
          : rdp::MaxChiTractSide::unavailable;
  discovery.peak_attempt = peak_attempt;
  discovery.peak_alignment_position = peak_alignment_position;
  discovery.information_rich_sites = information_rich_sites;
  discovery.initial_half_window = initial_half_window;
  discovery.grown_half_window = grown_half_window;
  discovery.critical_difference = critical_difference;
  discovery.maximum_chi_square = maximum_chi_square;
  discovery.raw_p_value = raw_p_value;
  discovery.within_triplet_p_value = within_triplet_p_value;
  discovery.corrected_p_value = signal.corrected_p_value;
  discovery.left_flank_chi_square = left_flank_chi_square;
  discovery.right_flank_chi_square = right_flank_chi_square;
  discovery.inside_parent_one_match_rate = inside_parent_one_match_rate;
  discovery.outside_parent_one_match_rate = outside_parent_one_match_rate;
  discovery.pair_similarity = signal.pair_similarity;
  discovery.missing_data_window_filter_applied =
      missing_data_window_filter_applied != 0;
  discovery.linear_edge_window_filter_applied =
      linear_edge_window_filter_applied != 0;
  return discovery.recombinant_local == discovery.target_local ? 1 : 0;
}

RDP_KEEPALIVE int rdp_restore_geneconv_discovery(
    std::uint32_t handle,
    std::uint32_t signal_id,
    std::uint32_t track,
    std::uint32_t polymorphic_sites,
    std::uint32_t positive_sites,
    std::uint32_t discordant_sites,
    std::uint32_t mismatch_penalty,
    std::uint32_t fragment_score,
    std::uint32_t critical_score,
    double lambda,
    double karlin_altschul_k,
    double raw_p_value) {
  Context* context = context_for(handle);
  if (!context || !context->restoring_scan ||
      signal_id >= context->restore_signals.size() || track > 5) {
    return 0;
  }
  auto& signal = context->restore_signals[signal_id];
  if (signal.method != rdp::SignalMethod::geneconv) return 0;
  auto& discovery = signal.geneconv_discovery;
  discovery.beginning = signal.beginning;
  discovery.ending = signal.ending;
  discovery.wraps_origin = signal.wraps_origin;
  discovery.informative_beginning = signal.informative_beginning;
  discovery.informative_ending = signal.informative_ending;
  const auto local_member = [&](std::uint32_t sequence) {
    for (std::uint8_t member = 0; member < signal.triplet.size(); ++member) {
      if (signal.triplet[member] == sequence) return member;
    }
    return static_cast<std::uint8_t>(0);
  };
  discovery.track = static_cast<std::uint8_t>(track);
  discovery.recombinant_local = local_member(signal.recombinant);
  discovery.major_parent_local = local_member(signal.major_parent);
  discovery.minor_parent_local = local_member(signal.minor_parent);
  discovery.candidate_pair = signal.candidate_pair;
  discovery.polymorphic_sites = polymorphic_sites;
  discovery.positive_sites = positive_sites;
  discovery.discordant_sites = discordant_sites;
  discovery.mismatch_penalty = mismatch_penalty;
  discovery.fragment_score = fragment_score;
  discovery.critical_score = critical_score;
  discovery.lambda = lambda;
  discovery.karlin_altschul_k = karlin_altschul_k;
  discovery.raw_p_value = raw_p_value;
  discovery.corrected_p_value = signal.corrected_p_value;
  discovery.pair_similarity = signal.pair_similarity;
  return 1;
}

RDP_KEEPALIVE int rdp_restore_threeseq_discovery(
    std::uint32_t handle,
    std::uint32_t signal_id,
    std::uint32_t target_local,
    int walk_direction,
    std::uint32_t information_rich_sites,
    std::uint32_t parent_one_matches,
    std::uint32_t parent_two_matches,
    std::uint32_t probability_excursion,
    std::uint32_t maximum_excursion,
    double raw_p_value,
    int exact_probability,
    int siegmund_fallback,
    int missing_data_split_applied) {
  Context* context = context_for(handle);
  if (!context || !context->restoring_scan ||
      signal_id >= context->restore_signals.size() || target_local > 2 ||
      (walk_direction != -1 && walk_direction != 1) ||
      information_rich_sites < 4 ||
      static_cast<std::uint64_t>(parent_one_matches) + parent_two_matches !=
          information_rich_sites ||
      probability_excursion > information_rich_sites ||
      maximum_excursion > parent_two_matches ||
      (missing_data_split_applied == 0 &&
       probability_excursion > maximum_excursion) ||
      (missing_data_split_applied == 0 && parent_two_matches > 0 &&
       maximum_excursion == 1) ||
      (missing_data_split_applied == 0 &&
       parent_two_matches >= parent_one_matches &&
       parent_two_matches - parent_one_matches == maximum_excursion) ||
      !std::isfinite(raw_p_value) || !(raw_p_value > 0.0) ||
      raw_p_value > 1.0 ||
      (exact_probability != 0 && exact_probability != 1) ||
      (siegmund_fallback != 0 && siegmund_fallback != 1) ||
      (missing_data_split_applied != 0 && missing_data_split_applied != 1) ||
      ((exact_probability != 0) == (siegmund_fallback != 0))) {
    return 0;
  }
  auto& signal = context->restore_signals[signal_id];
  if (signal.method != rdp::SignalMethod::threeseq ||
      signal.informative_sites != information_rich_sites ||
      signal.local_p_value != raw_p_value) {
    return 0;
  }
  auto& discovery = signal.threeseq_discovery;
  discovery.beginning = signal.beginning;
  discovery.ending = signal.ending;
  discovery.wraps_origin = signal.wraps_origin;
  discovery.informative_beginning = signal.informative_beginning;
  discovery.informative_ending = signal.informative_ending;
  const auto local_member = [&](std::uint32_t sequence) {
    for (std::uint8_t member = 0; member < signal.triplet.size(); ++member) {
      if (signal.triplet[member] == sequence) return member;
    }
    return static_cast<std::uint8_t>(3);
  };
  discovery.target_local = static_cast<std::uint8_t>(target_local);
  discovery.recombinant_local = local_member(signal.recombinant);
  discovery.major_parent_local = local_member(signal.major_parent);
  discovery.minor_parent_local = local_member(signal.minor_parent);
  discovery.candidate_pair = signal.candidate_pair;
  discovery.direction = walk_direction > 0
      ? rdp::ThreeSeqWalkDirection::ascent
      : rdp::ThreeSeqWalkDirection::descent;
  discovery.information_rich_sites = information_rich_sites;
  discovery.parent_one_matches = parent_one_matches;
  discovery.parent_two_matches = parent_two_matches;
  discovery.probability_excursion = probability_excursion;
  discovery.maximum_excursion = maximum_excursion;
  discovery.raw_p_value = raw_p_value;
  discovery.corrected_p_value = signal.corrected_p_value;
  discovery.pair_similarity = signal.pair_similarity;
  discovery.exact_probability = exact_probability != 0;
  discovery.siegmund_fallback = siegmund_fallback != 0;
  discovery.missing_data_split_applied = missing_data_split_applied != 0;
  const std::uint8_t parent_one = static_cast<std::uint8_t>((target_local + 1) % 3);
  const std::uint8_t parent_two = static_cast<std::uint8_t>((target_local + 2) % 3);
  const std::uint8_t expected_major = walk_direction > 0 ? parent_two : parent_one;
  const std::uint8_t expected_minor = walk_direction > 0 ? parent_one : parent_two;
  const std::uint8_t pair_first = std::min(discovery.target_local, expected_minor);
  const std::uint8_t pair_second = std::max(discovery.target_local, expected_minor);
  const std::uint8_t expected_pair = pair_first == 0 && pair_second == 1
      ? 0
      : pair_first == 0 && pair_second == 2 ? 1 : 2;
  return discovery.recombinant_local == discovery.target_local &&
          discovery.major_parent_local == expected_major &&
          discovery.minor_parent_local == expected_minor &&
          discovery.candidate_pair == expected_pair
      ? 1
      : 0;
}

RDP_KEEPALIVE int rdp_restore_bootscan_discovery(
    std::uint32_t handle,
    std::uint32_t signal_id,
    std::uint32_t supported_pair,
    std::uint32_t windows_scored,
    std::uint32_t usable_windows,
    std::uint32_t informative_sites,
    std::uint32_t tract_informative_sites,
    std::uint32_t tract_pair_matches,
    std::uint32_t outside_pair_matches,
    double maximum_pair_support,
    double mean_pair_support,
    double bootstrap_p_value,
    double raw_p_value,
    int erased_window_filter_applied) {
  Context* context = context_for(handle);
  if (!context || !context->restoring_scan ||
      signal_id >= context->restore_signals.size() || supported_pair > 2 ||
      windows_scored == 0 || usable_windows > windows_scored ||
      tract_informative_sites > informative_sites ||
      tract_pair_matches > tract_informative_sites ||
      outside_pair_matches > informative_sites - tract_informative_sites ||
      !std::isfinite(maximum_pair_support) || maximum_pair_support < 0.0 ||
      maximum_pair_support > 1.0 || !std::isfinite(mean_pair_support) ||
      mean_pair_support < 0.0 || mean_pair_support > maximum_pair_support ||
      !std::isfinite(bootstrap_p_value) || bootstrap_p_value <= 0.0 ||
      bootstrap_p_value > 1.0 || !std::isfinite(raw_p_value) ||
      raw_p_value <= 0.0 || raw_p_value > 1.0 ||
      (erased_window_filter_applied != 0 &&
       erased_window_filter_applied != 1)) {
    return 0;
  }
  auto& signal = context->restore_signals[signal_id];
  if (signal.method != rdp::SignalMethod::bootscan ||
      signal.informative_sites != informative_sites ||
      signal.candidate_pair != supported_pair || signal.local_p_value != raw_p_value) {
    return 0;
  }
  const auto local_member = [&](std::uint32_t sequence) {
    for (std::uint8_t member = 0; member < signal.triplet.size(); ++member) {
      if (signal.triplet[member] == sequence) return member;
    }
    return static_cast<std::uint8_t>(3);
  };
  auto& discovery = signal.bootscan_discovery;
  discovery.beginning = signal.beginning;
  discovery.ending = signal.ending;
  discovery.wraps_origin = signal.wraps_origin;
  discovery.informative_beginning = signal.informative_beginning;
  discovery.informative_ending = signal.informative_ending;
  discovery.supported_pair = static_cast<std::uint8_t>(supported_pair);
  discovery.recombinant_local = local_member(signal.recombinant);
  discovery.major_parent_local = local_member(signal.major_parent);
  discovery.minor_parent_local = local_member(signal.minor_parent);
  discovery.candidate_pair = signal.candidate_pair;
  discovery.windows_scored = windows_scored;
  discovery.usable_windows = usable_windows;
  discovery.informative_sites = informative_sites;
  discovery.tract_informative_sites = tract_informative_sites;
  discovery.tract_pair_matches = tract_pair_matches;
  discovery.outside_pair_matches = outside_pair_matches;
  discovery.maximum_pair_support = maximum_pair_support;
  discovery.mean_pair_support = mean_pair_support;
  discovery.bootstrap_p_value = bootstrap_p_value;
  discovery.raw_p_value = raw_p_value;
  discovery.corrected_p_value = signal.corrected_p_value;
  discovery.pair_similarity = signal.pair_similarity;
  discovery.erased_window_filter_applied = erased_window_filter_applied != 0;
  return discovery.recombinant_local < 3 &&
          discovery.major_parent_local < 3 &&
          discovery.minor_parent_local < 3
      ? 1
      : 0;
}

RDP_KEEPALIVE int rdp_restore_siscan_discovery(
    std::uint32_t handle,
    std::uint32_t signal_id,
    std::uint32_t global_pair,
    std::uint32_t candidate_pair,
    std::uint32_t outlier_sequence,
    std::uint32_t windows_in_region,
    std::uint32_t informative_sites,
    double permutation_draws,
    std::uint32_t selected_score,
    std::uint32_t selected_score_family,
    double maximum_z,
    double normal_tail_p_value,
    double region_length_adjusted_p_value,
    double window_adjusted_p_value) {
  Context* context = context_for(handle);
  if (!context || !context->restoring_scan ||
      signal_id >= context->restore_signals.size() || global_pair > 2 ||
      candidate_pair > 2 || candidate_pair == global_pair ||
      outlier_sequence >= context->alignment.sequence_count() ||
      windows_in_region == 0 || selected_score > 15 ||
      (selected_score_family != 1 && selected_score_family != 2) ||
      !std::isfinite(permutation_draws) || permutation_draws < 0.0 ||
      !std::isfinite(maximum_z) || maximum_z <= 0.0 ||
      !std::isfinite(normal_tail_p_value) || normal_tail_p_value <= 0.0 ||
      normal_tail_p_value > 1.0 ||
      !std::isfinite(region_length_adjusted_p_value) ||
      region_length_adjusted_p_value <= 0.0 ||
      region_length_adjusted_p_value > 1.0 ||
      !std::isfinite(window_adjusted_p_value) ||
      window_adjusted_p_value <= 0.0 || window_adjusted_p_value > 1.0) {
    return 0;
  }
  auto& signal = context->restore_signals[signal_id];
  if (signal.method != rdp::SignalMethod::siscan ||
      signal.informative_sites != informative_sites ||
      signal.candidate_pair != candidate_pair ||
      signal.local_p_value != window_adjusted_p_value) {
    return 0;
  }
  const auto local_member = [&](std::uint32_t sequence) {
    for (std::uint8_t member = 0; member < signal.triplet.size(); ++member) {
      if (signal.triplet[member] == sequence) return member;
    }
    return static_cast<std::uint8_t>(3);
  };
  auto& discovery = signal.siscan_discovery;
  discovery.beginning = signal.beginning;
  discovery.ending = signal.ending;
  discovery.wraps_origin = signal.wraps_origin;
  discovery.informative_beginning = signal.informative_beginning;
  discovery.informative_ending = signal.informative_ending;
  discovery.recombinant_local = local_member(signal.recombinant);
  discovery.major_parent_local = local_member(signal.major_parent);
  discovery.minor_parent_local = local_member(signal.minor_parent);
  discovery.global_pair = static_cast<std::uint8_t>(global_pair);
  discovery.candidate_pair = static_cast<std::uint8_t>(candidate_pair);
  discovery.outlier_sequence = outlier_sequence;
  discovery.windows_in_region = windows_in_region;
  discovery.informative_sites = informative_sites;
  discovery.permutation_draws = static_cast<std::size_t>(
      restored_counter(permutation_draws));
  discovery.selected_score = static_cast<std::uint8_t>(selected_score);
  discovery.selected_score_family = selected_score_family == 1
      ? rdp::SiscanScoreFamily::partition
      : rdp::SiscanScoreFamily::summed;
  discovery.maximum_z = maximum_z;
  discovery.normal_tail_p_value = normal_tail_p_value;
  discovery.region_length_adjusted_p_value =
      region_length_adjusted_p_value;
  discovery.window_adjusted_p_value = window_adjusted_p_value;
  discovery.corrected_p_value = signal.corrected_p_value;
  discovery.pair_similarity = signal.pair_similarity;
  return discovery.recombinant_local < 3 &&
          discovery.major_parent_local < 3 &&
          discovery.minor_parent_local < 3
      ? 1
      : 0;
}

RDP_KEEPALIVE int rdp_restore_scan_finish(
    std::uint32_t handle,
    std::uint32_t correction_tests,
    double cumulative_triplets,
    std::uint32_t scan_rounds,
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
    const std::uint8_t* cycle_termination,
    std::size_t cycle_termination_length) {
  Context* context = context_for(handle);
  if (!context || !context->loaded || !context->restoring_scan ||
      (!cycle_termination && cycle_termination_length > 0)) {
    return 0;
  }
  context->error.clear();
  context->scanner = std::make_unique<rdp::RdpScanner>(context->alignment);
  context->scanner->set_worker_threads(context->requested_worker_threads);
  const bool restored = context->scanner->restore(
      std::move(context->restore_options),
      std::move(context->restore_signals),
      correction_tests,
      restored_counter(cumulative_triplets),
      scan_rounds,
      restored_counter(maxchi_profiles_scanned),
      restored_counter(maxchi_peak_attempts),
      restored_counter(maxchi_candidates_found),
      restored_counter(maxchi_peak_limit_triplets),
      restored_counter(chimaera_profiles_scanned),
      restored_counter(chimaera_peak_attempts),
      restored_counter(chimaera_candidates_found),
      restored_counter(chimaera_peak_limit_targets),
      restored_counter(geneconv_fragments_scored),
      restored_counter(geneconv_qualified_fragments),
      restored_counter(geneconv_candidates_found),
      restored_counter(geneconv_overlap_rejections),
      restored_counter(geneconv_numerical_fallback_tracks),
      restored_counter(threeseq_profiles_scanned),
      restored_counter(threeseq_exact_evaluations),
      restored_counter(threeseq_approximate_evaluations),
      restored_counter(threeseq_candidates_found),
      restored_counter(bootscan_profiles_scanned),
      restored_counter(bootscan_candidate_regions_scored),
      restored_counter(bootscan_candidates_found),
      restored_counter(bootscan_pair_profiles_requested),
      restored_counter(bootscan_pair_profile_cache_hits),
      restored_counter(bootscan_pair_profile_cache_misses),
      restored_counter(bootscan_pair_profile_cache_evictions),
      restored_counter(bootscan_pair_profile_cache_peak_bytes),
      restored_counter(siscan_profiles_scanned),
      restored_counter(siscan_windows_scored),
      restored_counter(siscan_candidate_regions_scored),
      restored_counter(siscan_candidates_found),
      restored_counter(siscan_permutation_draws),
      restored_counter(siscan_context_builds),
      restored_counter(siscan_context_pair_comparisons),
      restored_counter(siscan_context_tree_merges),
      restored_counter(siscan_random_values_generated),
      bytes_to_string(cycle_termination, cycle_termination_length),
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
