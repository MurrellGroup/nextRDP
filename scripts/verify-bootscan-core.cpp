#include "alignment.hpp"
#include "bootscan.hpp"
#include "rdp_api.h"

#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace {

std::uint32_t random_state = 0x51a7c3d9U;

char next_base() {
  random_state = random_state * 1664525U + 1013904223U;
  constexpr char bases[] = "ACGT";
  return bases[(random_state >> 16U) & 3U];
}

std::string random_sequence(std::size_t length) {
  std::string sequence(length, 'A');
  for (char& base : sequence) base = next_base();
  return sequence;
}

std::string periodic_mutant(const std::string& source, std::size_t period) {
  std::string result = source;
  constexpr char bases[] = "ACGT";
  for (std::size_t index = 0; index < result.size(); index += period) {
    const std::size_t current = std::string_view(bases).find(result[index]);
    result[index] = bases[(current + 1U + (index & 1U)) & 3U];
  }
  return result;
}

int fail(const std::string& message) {
  std::cerr << "BootScan core check failed: " << message << '\n';
  return 1;
}

bool positive_json_integer(std::string_view json, std::string_view field) {
  const std::string marker = "\"" + std::string(field) + "\":";
  const std::size_t position = json.find(marker);
  if (position == std::string_view::npos) return false;
  const char* beginning = json.data() + position + marker.size();
  char* ending = nullptr;
  const unsigned long long value = std::strtoull(beginning, &ending, 10);
  return ending != beginning && value > 0;
}

}  // namespace

int main() {
  constexpr std::size_t length = 600;
  const std::string parent_one = random_sequence(length);
  const std::string parent_two = periodic_mutant(parent_one, 3);
  std::string recombinant = parent_one;
  recombinant.replace(length / 2, length / 2, parent_two, length / 2, length / 2);
  const std::string background = random_sequence(length);

  auto parsed = rdp::build_alignment(
      "FASTA",
      {"recombinant", "parent-one", "parent-two", "background"},
      {recombinant, parent_one, parent_two, background});
  if (!parsed.ok()) return fail(parsed.error);
  const rdp::Alignment& alignment = parsed.alignment;

  rdp::BootscanDiscoveryOptions options;
  options.circular = false;
  options.bonferroni = false;
  options.p_value_cutoff = 1.0;
  options.window_sites = 100;
  options.step_sites = 20;
  options.bootstrap_replicates = 60;
  options.support_cutoff = 0.70;
  options.random_seed = 3;
  options.pair_cache_limit_bytes = 4U * 1024U * 1024U;

  rdp::BootscanWorkspace workspace;
  std::vector<rdp::BootscanDiscoveryCandidate> candidates;
  const std::array<std::uint32_t, 3> first_triplet{0, 1, 2};
  const std::array<double, 3> first_similarity{
      alignment.similarity(0, 1),
      alignment.similarity(0, 2),
      alignment.similarity(1, 2),
  };
  const auto first = rdp::bootscan_discover(
      alignment,
      first_triplet,
      {},
      first_similarity,
      options,
      workspace,
      candidates);
  if (!first.profile_available || first.pair_profile_cache_misses != 3) {
    return fail("the opening triplet did not build its three distance profiles");
  }
  if (candidates.size() < 2) {
    return fail("the deterministic two-parent mosaic did not produce both supported regions");
  }
  bool saw_first_parent_region = false;
  bool saw_second_parent_region = false;
  for (const auto& candidate : candidates) {
    if (candidate.recombinant_local != 0 ||
        !(candidate.raw_p_value > 0.0 && candidate.raw_p_value < 1.0)) {
      return fail("a mosaic candidate has invalid roles or binomial probability");
    }
    saw_first_parent_region |= candidate.supported_pair == 0;
    saw_second_parent_region |= candidate.supported_pair == 1;
  }
  if (!saw_first_parent_region || !saw_second_parent_region) {
    return fail("strict pair voting did not recover both mosaic halves");
  }

  const std::array<std::uint32_t, 3> shared_pair_triplet{0, 1, 3};
  const std::array<double, 3> shared_pair_similarity{
      alignment.similarity(0, 1),
      alignment.similarity(0, 3),
      alignment.similarity(1, 3),
  };
  const auto second = rdp::bootscan_discover(
      alignment,
      shared_pair_triplet,
      {},
      shared_pair_similarity,
      options,
      workspace,
      candidates);
  if (second.pair_profile_cache_hits != 1 ||
      second.pair_profile_cache_misses != 2 ||
      workspace.pair_profile_cache_hits != 1) {
    return fail("a shared exact pair was not reused across triplets");
  }
  if (workspace.pair_profile_cache_bytes == 0 ||
      workspace.pair_profile_cache_peak_bytes > options.pair_cache_limit_bytes) {
    return fail("the bounded cache byte accounting is invalid");
  }

  rdp::bootscan_reset_discovery_cache(workspace);
  if (!workspace.pair_profile_cache.empty() ||
      workspace.pair_profile_cache_bytes != 0) {
    return fail("round-boundary invalidation retained stale profiles");
  }

  std::string fasta;
  const std::array<std::pair<std::string_view, const std::string*>, 4> records{{
      {"recombinant", &recombinant},
      {"parent-one", &parent_one},
      {"parent-two", &parent_two},
      {"background", &background},
  }};
  for (const auto& [name, sequence] : records) {
    fasta += '>';
    fasta += name;
    fasta += '\n';
    fasta += *sequence;
    fasta += '\n';
  }
  const std::uint32_t handle = rdp_create();
  if (handle == 0) return fail("the public API could not create a context");
  const auto destroy_context = [&]() { rdp_destroy(handle); };
  if (rdp_load_alignment(
          handle,
          reinterpret_cast<const std::uint8_t*>(fasta.data()),
          fasta.size()) != 1) {
    const std::string message = rdp_get_error(handle);
    destroy_context();
    return fail("the public API could not load the mosaic: " + message);
  }
  const std::array<std::uint32_t, 4> reference_groups{};
  const std::array<std::uint8_t, 4> flags{};
  if (rdp_scan_begin(
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
          1,
          0,
          100,
          20,
          60,
          0.70,
          3,
          0,
          0,
          200,
          20,
          100,
          1000,
          3,
          0,
          0,
          reference_groups.data(),
          reference_groups.size(),
          flags.data(),
          flags.size(),
          flags.data(),
          flags.size()) != 1) {
    const std::string message = rdp_get_error(handle);
    destroy_context();
    return fail("the public API could not start primary BootScan: " + message);
  }
  int status = 0;
  for (std::size_t batch = 0; batch < 1000 && status != 4 && status != 3; ++batch) {
    status = rdp_scan_batch(handle, 64);
    if (status < 0) {
      const std::string message = rdp_get_error(handle);
      destroy_context();
      return fail("the primary-BootScan API scan failed: " + message);
    }
  }
  const std::string progress = rdp_get_progress_json(handle);
  if (status != 4 ||
      !positive_json_integer(progress, "bootscanCandidatesFound") ||
      !positive_json_integer(progress, "bootscanPairProfileCacheHits") ||
      !positive_json_integer(progress, "bootscanPairProfileCacheMisses")) {
    destroy_context();
    return fail("primary BootScan did not reach the public cyclic scheduler/cache");
  }
  rdp_cancel(handle);
  if (rdp_scan_batch(handle, 1) != 3 || rdp_reconcile(handle) != 1) {
    const std::string message = rdp_get_error(handle);
    destroy_context();
    return fail("the completed primary-BootScan event could not be reconciled: " + message);
  }
  const std::string results = rdp_get_results_json(handle);
  const std::size_t bootscan_method = results.find("\"method\":\"BOOTSCAN\"");
  if (results.find("\"bootscanPrimaryEnabled\":true") == std::string::npos ||
      bootscan_method == std::string::npos ||
      results.find("\"bootscanDiscovery\":{\"status\":\"source-shaped-active-unvalidated\"") ==
          std::string::npos) {
    destroy_context();
    return fail("primary BootScan evidence did not survive public result reconciliation");
  }
  if (results.find("\"distanceEncoding\":\"source-midpoint-ultrametric-ranks\"") ==
          std::string::npos ||
      results.find("\"flankVariableSiteTarget\":20") == std::string::npos ||
      results.find("\"collapseEncoding\":\"parent-rank-promotion-no-recompression\"") ==
          std::string::npos) {
    destroy_context();
    return fail("supplied event-tree provenance did not survive public reconciliation");
  }
  const std::string event_trees = rdp_get_event_trees_json(handle, 0);
  if (event_trees.find("\"njKernel\":\"supplied-clearcut-float\"") ==
          std::string::npos ||
      event_trees.find("\"analyticalBranchParsing\":"
                       "\"four-decimal-clamped-complete-edge-repair\"") ==
          std::string::npos ||
      event_trees.find("\"regions\":[") == std::string::npos) {
    destroy_context();
    return fail("the lazy six-region tree endpoint lost its Session 21 contract");
  }
  const std::string invalid_phylpro = rdp_get_event_phylpro_json(handle, 0, 40, 2, 0);
  if (!invalid_phylpro.empty() ||
      std::string(rdp_get_error(handle)).find("gap mode is unavailable") == std::string::npos) {
    destroy_context();
    return fail("the public PHYLPRO ABI accepted an unknown gap policy");
  }
  const std::string event_phylpro = rdp_get_event_phylpro_json(handle, 0, 40, 0, 0);
  if (event_phylpro.find(
          "\"kernel\":\"FindSubSeqPP-MakePDstMat-UpdatePDstMat-PPRegression\"") ==
          std::string::npos ||
      event_phylpro.find("\"optimization\":\"three-target-rows-linear-in-context\"") ==
          std::string::npos ||
      event_phylpro.find("\"significanceTest\":\"not-implemented-in-supplied-rdp5\"") ==
          std::string::npos ||
      event_phylpro.find("\"points\":[]") != std::string::npos) {
    destroy_context();
    return fail("the lazy PHYLPRO event-review endpoint lost its Session 22 contract");
  }
  const std::string project = rdp_export_project_json(handle);
  if (project.find("\"schema\":\"org.rdp-web.project/v1alpha19\"") ==
          std::string::npos ||
      project.find("\"engineVersion\":\"0.25.0-session-25\"") ==
          std::string::npos) {
    destroy_context();
    return fail("the Session 25 project version/v1alpha19 schema contract is missing");
  }
  const std::size_t signal_id_field = results.rfind("\"id\":", bootscan_method);
  if (signal_id_field == std::string::npos) {
    destroy_context();
    return fail("the reconciled BootScan signal has no public ID");
  }
  const char* signal_id_beginning =
      results.c_str() + signal_id_field + std::string_view("\"id\":").size();
  char* signal_id_ending = nullptr;
  const unsigned long signal_id =
      std::strtoul(signal_id_beginning, &signal_id_ending, 10);
  if (signal_id_ending == signal_id_beginning) {
    destroy_context();
    return fail("the reconciled BootScan signal ID is malformed");
  }
  const std::string plot = rdp_get_signal_plot_json(
      handle, static_cast<std::uint32_t>(signal_id));
  if (plot.find("\"method\":\"BOOTSCAN\"") == std::string::npos ||
      plot.find("\"metric\":\"bootstrap-support\"") == std::string::npos ||
      plot.find("\"points\":[]") != std::string::npos) {
    destroy_context();
    return fail("the public API did not reconstruct the BootScan support plot");
  }
  destroy_context();

  std::cout << "BootScan core verified: two mosaic regions, MakeScoresBS probability, "
               "shared-pair cache hit, round invalidation, reconciled public-API evidence, "
               "support plot, Session 21 tree endpoint, Session 22 PHYLPRO endpoint, Session 23 scheduler, and "
               "v1alpha19 project export.\n";
  return 0;
}
