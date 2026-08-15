#include "alignment.hpp"
#include "rdp_api.h"
#include "rdp_method.hpp"
#include "siscan.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

[[noreturn]] void fail(const std::string& message) {
  std::cerr << "SISCAN core verification failed: " << message << '\n';
  std::exit(1);
}

void require(bool condition, const std::string& message) {
  if (!condition) fail(message);
}

bool positive_json_integer(const std::string& json, const std::string& field) {
  const std::string marker = "\"" + field + "\":";
  const std::size_t position = json.find(marker);
  if (position == std::string::npos) return false;
  const char* beginning = json.c_str() + position + marker.size();
  char* ending = nullptr;
  return std::strtoull(beginning, &ending, 10) > 0 && ending != beginning;
}

}  // namespace

int main() {
  constexpr std::size_t length = 240;
  const std::string alphabet = "ACGT";
  std::string parent_one(length, 'A');
  std::string parent_two(length, 'A');
  std::string outlier(length, 'A');
  for (std::size_t position = 0; position < length; ++position) {
    parent_one[position] = alphabet[position % 4];
    parent_two[position] = alphabet[(position + 1) % 4];
    outlier[position] = alphabet[(position + 2) % 4];
  }
  std::string recombinant = parent_one;
  for (std::size_t position = 80; position <= 139; ++position) {
    recombinant[position] = parent_two[position];
  }
  auto parsed = rdp::build_alignment(
      "synthetic-siscan",
      {"recombinant", "parent-one", "parent-two", "outlier"},
      {recombinant, parent_one, parent_two, outlier});
  require(parsed.ok(), parsed.error);

  rdp::SiscanOptions options;
  options.circular = false;
  options.bonferroni = false;
  options.p_value_cutoff = 1.0;
  options.window_sites = 40;
  options.step_sites = 5;
  options.scan_permutations = 20;
  options.p_value_permutations = 100;
  options.random_seed = 3;
  options.fast_scan = true;
  const std::array<std::uint32_t, 3> triplet{{0, 1, 2}};
  const std::array<double, 3> similarities{
      parsed.alignment.similarity(0, 1),
      parsed.alignment.similarity(0, 2),
      parsed.alignment.similarity(1, 2),
  };
  std::vector<std::uint32_t> origins(parsed.alignment.sequence_count());
  std::iota(origins.begin(), origins.end(), 0U);
  std::vector<std::uint8_t> disabled(parsed.alignment.sequence_count(), 0);
  std::vector<std::uint8_t> missing(length, 0);
  rdp::SiscanWorkspace workspace;
  std::vector<rdp::SiscanDiscoveryCandidate> candidates;
  const rdp::SiscanDiscoverySummary summary = rdp::siscan_discover(
      parsed.alignment,
      triplet,
      missing,
      similarities,
      origins,
      disabled,
      options,
      workspace,
      candidates);
  require(summary.profile_available, "synthetic profile was unavailable");
  require(summary.outlier_available && summary.outlier_sequence == 3,
          "nearest fourth sequence was not selected deterministically");
  require(summary.windows_scored > 0, "fast window screen scored no windows");
  require(summary.candidate_regions_scored > 0,
          "pair-switch windows did not form a candidate region");
  require(!candidates.empty(), "the planted sister-scan tract was not emitted");
  const auto& candidate = candidates.front();
  require(candidate.global_pair == 0 && candidate.candidate_pair == 1,
          "the global/local sister-pair switch was misclassified");
  require(candidate.recombinant_local == 0 &&
              candidate.major_parent_local == 1 &&
              candidate.minor_parent_local == 2,
          "pair geometry did not recover the planted roles");
  require(candidate.beginning <= 90 && candidate.ending >= 130,
          "the planted tract was not covered by the shrunken region");
  require(candidate.maximum_z > 0.0 &&
              candidate.normal_tail_p_value > 0.0 &&
              candidate.normal_tail_p_value <= 1.0,
          "the native NormalZ probability chain is invalid");
  require(workspace.context_builds == 1 && workspace.context_tree_merges == 3,
          "the round-wide WPGMA context was not constructed once");
  constexpr std::array<std::uint8_t, 20> microsoft_seed_three_prefix{{
      3, 4, 4, 3, 9, 7, 9, 10, 5, 2,
      11, 1, 7, 8, 1, 9, 6, 4, 8, 9,
  }};
  require(
      workspace.vertical_random_prefix.size() >=
          microsoft_seed_three_prefix.size() &&
          std::equal(
              microsoft_seed_three_prefix.begin(),
              microsoft_seed_three_prefix.end(),
              workspace.vertical_random_prefix.begin()),
      "the seeded Microsoft-CRT random prefix drifted from MakeVRand");

  rdp::SiscanWorkspace same_origin_workspace;
  std::vector<std::uint32_t> same_origin = origins;
  same_origin[3] = 0;
  std::vector<rdp::SiscanDiscoveryCandidate> blocked_candidates;
  const rdp::SiscanDiscoverySummary same_origin_summary = rdp::siscan_discover(
      parsed.alignment,
      triplet,
      missing,
      similarities,
      same_origin,
      disabled,
      options,
      same_origin_workspace,
      blocked_candidates);
  require(!same_origin_summary.outlier_available && blocked_candidates.empty(),
          "a same-origin fragment was accepted as its own fourth sequence");

  rdp::SiscanWorkspace disabled_workspace;
  std::vector<std::uint8_t> outlier_disabled = disabled;
  outlier_disabled[3] = 1;
  const rdp::SiscanDiscoverySummary disabled_summary = rdp::siscan_discover(
      parsed.alignment,
      triplet,
      missing,
      similarities,
      origins,
      outlier_disabled,
      options,
      disabled_workspace,
      blocked_candidates);
  require(!disabled_summary.outlier_available && blocked_candidates.empty(),
          "a disabled original sequence was accepted as the fourth sequence");
  require(disabled_workspace.context_pair_comparisons == 3 &&
              disabled_workspace.context_tree_merges == 2,
          "disabled rows were not pruned from the WPGMA context work");

  const std::size_t generated = workspace.random_values_generated;
  std::vector<rdp::SiscanDiscoveryCandidate> replay;
  const rdp::SiscanDiscoverySummary replay_summary = rdp::siscan_discover(
      parsed.alignment,
      triplet,
      missing,
      similarities,
      origins,
      disabled,
      options,
      workspace,
      replay);
  require(workspace.context_builds == 1,
          "unchanged triplet replay rebuilt the cached distance tree");
  require(workspace.random_values_generated == generated,
          "unchanged triplet replay regenerated the flat random prefix");
  require(replay_summary.emitted_candidates == summary.emitted_candidates &&
              replay.size() == candidates.size(),
          "seeded discovery replay changed the candidate count");
  require(replay.front().beginning == candidate.beginning &&
              replay.front().ending == candidate.ending &&
              replay.front().corrected_p_value == candidate.corrected_p_value,
          "seeded discovery replay was not bit-deterministic");

  const rdp::SiscanRecheckEvidence recheck = rdp::siscan_recheck(
      parsed.alignment,
      triplet,
      missing,
      origins,
      disabled,
      candidate.beginning,
      candidate.ending,
      options,
      workspace);
  require(recheck.profile_available && recheck.outlier_sequence == 3,
          "fixed-bound SISCAN confirmation was unavailable");
  require(recheck.maximum_z == candidate.maximum_z &&
              recheck.corrected_p_value == candidate.corrected_p_value,
          "fixed-bound confirmation diverged from discovery scoring");

  options.fast_scan = false;
  const rdp::SiscanPlotProfile plot = rdp::siscan_plot_profile(
      parsed.alignment,
      triplet,
      missing,
      origins,
      disabled,
      options,
      workspace,
      3);
  require(plot.available && !plot.coordinates.empty(),
          "all-window review plot was unavailable");
  require(plot.coordinates.size() == plot.pair_z[0].size() &&
              plot.coordinates.size() == plot.pair_z[1].size() &&
              plot.coordinates.size() == plot.pair_z[2].size(),
          "review plot pair traces do not share a coordinate grid");
  bool plot_has_negative = false;
  bool plot_has_positive = false;
  for (const auto& trace : plot.pair_z) {
    plot_has_negative = plot_has_negative || std::any_of(
        trace.begin(), trace.end(), [](double value) { return value < 0.0; });
    plot_has_positive = plot_has_positive || std::any_of(
        trace.begin(), trace.end(), [](double value) { return value > 0.0; });
  }
  require(plot_has_negative && plot_has_positive,
          "the compact SISCAN plot did not preserve signed Z extrema");

  rdp::siscan_reset_round_context(workspace);
  std::vector<rdp::SiscanDiscoveryCandidate> next_round;
  (void)rdp::siscan_discover(
      parsed.alignment,
      triplet,
      missing,
      similarities,
      origins,
      disabled,
      options,
      workspace,
      next_round);
  require(workspace.context_builds == 2,
          "cyclic round invalidation did not rebuild state-dependent context");
  require(workspace.random_values_generated == generated,
          "cyclic round invalidation discarded the reusable random template");

  rdp::ScanOptions restored_options;
  restored_options.mask.assign(parsed.alignment.sequence_count(), 0);
  restored_options.disabled.assign(parsed.alignment.sequence_count(), 0);
  restored_options.reference_groups.assign(parsed.alignment.sequence_count(), 0);
  rdp::RdpScanner restored_scanner(parsed.alignment);
  std::string restored_error;
  const bool restored = restored_scanner.restore(
      restored_options,
      {},
      4, 4, 1,
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0,
      7, 11, 13, 17,
      "restored-project",
      restored_error);
  require(restored, "SISCAN metric restore failed: " + restored_error);
  const std::string restored_progress = restored_scanner.progress_json();
  require(
      restored_progress.find("\"siscanContextBuilds\":7") !=
              std::string::npos &&
          restored_progress.find("\"siscanContextPairComparisons\":11") !=
              std::string::npos &&
          restored_progress.find("\"siscanContextTreeMerges\":13") !=
              std::string::npos &&
          restored_progress.find("\"siscanRandomValuesGenerated\":17") !=
              std::string::npos,
      "SISCAN cache metrics did not survive project restoration");

  std::string fasta;
  const std::array<std::pair<const char*, const std::string*>, 4> records{{
      {"siscan-recombinant", &recombinant},
      {"siscan-parent-one", &parent_one},
      {"siscan-parent-two", &parent_two},
      {"siscan-outlier", &outlier},
  }};
  for (const auto& [name, sequence] : records) {
    fasta += '>';
    fasta += name;
    fasta += '\n';
    fasta += *sequence;
    fasta += '\n';
  }
  const std::uint32_t handle = rdp_create();
  require(handle != 0, "the public API could not create a context");
  const auto destroy_context = [&]() { rdp_destroy(handle); };
  if (rdp_load_alignment(
          handle,
          reinterpret_cast<const std::uint8_t*>(fasta.data()),
          fasta.size()) != 1) {
    const std::string message = rdp_get_error(handle);
    destroy_context();
    fail("the public API could not load the SISCAN fixture: " + message);
  }
  const std::array<std::uint32_t, 4> reference_groups{};
  const std::array<std::uint8_t, 4> flags{};
  const auto begin_public_siscan = [&]() {
    return rdp_scan_begin(
          handle,
          0,
          1,
          0.05,
          30,
          0,
          70,
          0,
          60,
          0,
          1,
          1,
          0,
          0,
          0,
          200,
          20,
          100,
          0.70,
          3,
          1,
          0,
          40,
          5,
          20,
          100,
          3,
          0,
          0,
          reference_groups.data(),
          reference_groups.size(),
          flags.data(),
          flags.size(),
          flags.data(),
          flags.size());
  };
  if (begin_public_siscan() != 1) {
    const std::string message = rdp_get_error(handle);
    destroy_context();
    fail("the public API could not start SISCAN: " + message);
  }
  int status = 0;
  for (std::size_t batch = 0;
       batch < 1000 && status != 4 && status != 3;
       ++batch) {
    status = rdp_scan_batch(handle, 64);
    if (status < 0) {
      const std::string message = rdp_get_error(handle);
      destroy_context();
      fail("the public SISCAN scan failed: " + message);
    }
  }
  const std::string progress = rdp_get_progress_json(handle);
  if (status != 4 ||
      !positive_json_integer(progress, "siscanProfilesScanned") ||
      !positive_json_integer(progress, "siscanWindowsScored") ||
      !positive_json_integer(progress, "siscanCandidateRegionsScored") ||
      !positive_json_integer(progress, "siscanCandidatesFound") ||
      !positive_json_integer(progress, "siscanPermutationDraws") ||
      progress.find("\"siscanContextBuilds\":1") == std::string::npos ||
      progress.find("\"siscanContextTreeMerges\":3") == std::string::npos ||
      !positive_json_integer(progress, "siscanRandomValuesGenerated")) {
    destroy_context();
    fail("SISCAN did not traverse the public scheduler/context caches");
  }
  rdp_cancel(handle);
  status = rdp_scan_batch(handle, 1);
  if (status != 3 || rdp_reconcile(handle) != 1) {
    const std::string message = rdp_get_error(handle);
    destroy_context();
    fail("the completed public SISCAN event could not be reconciled: " + message);
  }
  const std::string results = rdp_get_results_json(handle);
  const std::size_t siscan_method = results.find("\"method\":\"SISCAN\"");
  if (results.find("\"siscanPrimaryEnabled\":true") == std::string::npos ||
      siscan_method == std::string::npos ||
      results.find("\"siscanDiscovery\":{\"status\":\"source-shaped-active-unvalidated\"") ==
          std::string::npos ||
      results.find("\"permutationGenerator\":\"microsoft-crt-flat-prefix\"") ==
          std::string::npos) {
    destroy_context();
    fail("SISCAN evidence did not survive public result reconciliation");
  }
  const std::size_t signal_id_field = results.rfind("\"id\":", siscan_method);
  if (signal_id_field == std::string::npos) {
    destroy_context();
    fail("the reconciled SISCAN signal has no public ID");
  }
  const char* signal_id_beginning =
      results.c_str() + signal_id_field + std::string_view("\"id\":").size();
  char* signal_id_ending = nullptr;
  const unsigned long signal_id =
      std::strtoul(signal_id_beginning, &signal_id_ending, 10);
  if (signal_id_ending == signal_id_beginning) {
    destroy_context();
    fail("the reconciled SISCAN signal ID is malformed");
  }
  const std::string plot_json = rdp_get_signal_plot_json(
      handle, static_cast<std::uint32_t>(signal_id));
  if (plot_json.find("\"method\":\"SISCAN\"") == std::string::npos ||
      plot_json.find("\"metric\":\"sister-scan-z-score\"") == std::string::npos ||
      plot_json.find("\"minimumValue\":-") == std::string::npos ||
      plot_json.find("\"points\":[") == std::string::npos) {
    destroy_context();
    fail("the public SISCAN review plot was not reconstructed");
  }
  const std::string csv = rdp_export_csv(handle);
  if (csv.find("SISCAN triplet recheck status") == std::string::npos ||
      csv.find("SSXoverC/GetSSOL/Get3Score/GetPScores2") == std::string::npos) {
    destroy_context();
    fail("the public CSV omitted the SISCAN audit trail");
  }
  const std::string project = rdp_export_project_json(handle);
  if (project.find("\"schema\":\"org.rdp-web.project/v1alpha19\"") ==
          std::string::npos ||
      project.find("\"engineVersion\":\"0.25.0-session-25\"") ==
          std::string::npos) {
    destroy_context();
    fail("the Session 25 SISCAN project version/v1alpha19 schema contract is missing");
  }

  // A browser context is intentionally reusable. Its second scan must retain
  // the deterministic random template without leaking the first scan's cache
  // instrumentation into the new progress/result contract.
  if (begin_public_siscan() != 1) {
    const std::string message = rdp_get_error(handle);
    destroy_context();
    fail("the public API could not restart SISCAN: " + message);
  }
  const int rerun_status = rdp_scan_batch(handle, 1);
  if (rerun_status < 0) {
    const std::string message = rdp_get_error(handle);
    destroy_context();
    fail("the restarted SISCAN scan failed: " + message);
  }
  const std::string rerun_progress = rdp_get_progress_json(handle);
  if (rerun_progress.find("\"siscanContextBuilds\":1") == std::string::npos ||
      rerun_progress.find("\"siscanContextTreeMerges\":3") ==
          std::string::npos) {
    destroy_context();
    fail("a restarted analysis inherited stale SISCAN cache telemetry");
  }
  rdp_cancel(handle);
  (void)rdp_scan_batch(handle, 1);
  destroy_context();

  std::cout << "SISCAN core verification passed: "
            << summary.windows_scored << " windows, "
            << summary.candidate_regions_scored << " regions, "
            << summary.permutation_draws << " draws, Z="
            << candidate.maximum_z << ", p="
            << candidate.corrected_p_value
            << "; same-origin/disabled outlier gates and public scheduler/reconciliation/"
               "signed plot/CSV/project/context-restart passed\n";
  return 0;
}
