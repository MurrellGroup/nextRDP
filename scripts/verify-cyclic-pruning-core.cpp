#include "rdp_api.h"

#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

namespace {

std::uint32_t random_state = 123456789U;

char next_base() {
  random_state = random_state * 1664525U + 1013904223U;
  constexpr char bases[] = "ACGT";
  return bases[random_state % 4U];
}

std::string random_sequence(std::size_t length) {
  std::string sequence(length, 'A');
  for (char& base : sequence) base = next_base();
  return sequence;
}

std::string periodic_mutant(const std::string& source, std::size_t period) {
  std::string result = source;
  constexpr std::string_view bases = "ACGT";
  for (std::size_t index = 0; index < result.size(); ++index) {
    if (index % period != 0) continue;
    const std::size_t current = bases.find(result[index]);
    result[index] = bases[(current + 1U + (index % 2U)) % bases.size()];
  }
  return result;
}

std::string mosaic(
    const std::string& first,
    const std::string& second,
    std::size_t split) {
  std::string result = first;
  constexpr std::string_view bases = "ACGT";
  for (std::size_t index = 0; index < result.size(); ++index) {
    result[index] = index < split ? first[index] : second[index];
    if (index % 47 == 0 && first[index] == second[index]) {
      result[index] = bases[(bases.find(result[index]) + 2U) % bases.size()];
    }
  }
  return result;
}

int fail(const std::string& message) {
  std::cerr << "Cyclic pruning core check failed: " << message << '\n';
  return 1;
}

unsigned long long json_integer(
    std::string_view json,
    std::string_view field) {
  const std::string marker = "\"" + std::string(field) + "\":";
  const std::size_t position = json.find(marker);
  if (position == std::string_view::npos) return 0;
  const char* beginning = json.data() + position + marker.size();
  char* ending = nullptr;
  const unsigned long long value = std::strtoull(beginning, &ending, 10);
  return ending == beginning ? 0 : value;
}

}  // namespace

int main() {
  constexpr std::size_t length = 900;
  const std::string parent_a = random_sequence(length);
  const std::string parent_b = periodic_mutant(parent_a, 3);
  const std::string recombinant_a = mosaic(parent_a, parent_b, 450);
  const std::string parent_c = random_sequence(length);
  const std::string parent_d = periodic_mutant(parent_c, 4);
  const std::string recombinant_b = mosaic(parent_c, parent_d, 300);
  std::vector<std::string> sequences{
      parent_a,
      parent_b,
      recombinant_a,
      parent_c,
      parent_d,
      recombinant_b,
  };
  for (std::size_t index = 0; index < 4; ++index) {
    sequences.push_back(random_sequence(length));
  }

  std::string fasta;
  for (std::size_t index = 0; index < sequences.size(); ++index) {
    fasta += ">cyclic-" + std::to_string(index) + "\n";
    fasta += sequences[index] + "\n";
  }

  const std::uint32_t handle = rdp_create();
  if (handle == 0) return fail("the public API could not create a context");
  const auto destroy = [&]() { rdp_destroy(handle); };
  if (rdp_load_alignment(
          handle,
          reinterpret_cast<const std::uint8_t*>(fasta.data()),
          fasta.size()) != 1) {
    const std::string error = rdp_get_error(handle);
    destroy();
    return fail("alignment load failed: " + error);
  }
  const std::vector<std::uint32_t> reference_groups(sequences.size(), 0);
  const std::vector<std::uint8_t> flags(sequences.size(), 0);
  if (rdp_scan_begin(
          handle,
          1,
          1,
          0.05,
          30,
          0, 70,
          0, 60,
          0, 1, 1,
          0,
          0, 0, 200, 20, 100, 0.7, 3,
          0, 0, 200, 20, 100, 1000, 3,
          0,
          0,
          reference_groups.data(),
          reference_groups.size(),
          flags.data(),
          flags.size(),
          flags.data(),
          flags.size()) != 1) {
    const std::string error = rdp_get_error(handle);
    destroy();
    return fail("scan start failed: " + error);
  }

  int status = 0;
  for (std::size_t batch = 0; batch < 10000 && status != 3; ++batch) {
    status = rdp_scan_batch(handle, 10000);
    if (status < 0) {
      const std::string error = rdp_get_error(handle);
      destroy();
      return fail("scan failed: " + error);
    }
  }
  const std::string progress = rdp_get_progress_json(handle);
  if (status != 3 || json_integer(progress, "eventCount") < 2 ||
      json_integer(progress, "processedTriplets") !=
          json_integer(progress, "totalTriplets") ||
      json_integer(progress, "cleanTripletsPruned") == 0 ||
      json_integer(progress, "methodScansSkipped") == 0 ||
      json_integer(progress, "invalidScheduleTripletsSkipped") == 0 ||
      json_integer(progress, "fragmentSequencesPruned") == 0) {
    std::cerr << progress << '\n';
    destroy();
    return fail("the clean-triplet or event-free-fragment shortcut was not exercised");
  }
  if (rdp_reconcile(handle) != 1) {
    const std::string error = rdp_get_error(handle);
    destroy();
    return fail("reconciliation failed: " + error);
  }
  const std::string results = rdp_get_results_json(handle);
  if (results.find("\"cycleMode\":\"strongest-first-tract-erasure-with-bounded-fragment-reentry\"") ==
      std::string::npos) {
    destroy();
    return fail("cyclic result provenance is missing");
  }
  if (std::getenv("RDP_DUMP_RESULTS") != nullptr) {
    std::cout << "RDP_RESULTS_JSON=" << results << '\n';
  }
  destroy();
  std::cout << "Cyclic pruning verified: "
            << json_integer(progress, "cleanTripletsPruned")
            << " clean triplets pruned, "
            << json_integer(progress, "methodScansSkipped")
            << " method scans skipped, "
            << json_integer(progress, "invalidScheduleTripletsSkipped")
            << " same-origin combinations bypassed, and "
            << json_integer(progress, "fragmentSequencesPruned")
            << " event-free fragments removed with swap/reindex compaction.\n";
  return 0;
}
