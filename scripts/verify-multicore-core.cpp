#include "rdp_api.h"

#include <chrono>
#include <cstdint>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

namespace {

std::uint32_t random_state = 0x5a17c9e3U;

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
    const std::size_t base = bases.find(result[index]);
    result[index] = bases[(base + 1U + index % 2U) % bases.size()];
  }
  return result;
}

std::string mosaic(
    const std::string& first,
    const std::string& second,
    std::size_t beginning,
    std::size_t ending) {
  std::string result = first;
  for (std::size_t index = beginning; index < ending; ++index) {
    result[index] = second[index];
  }
  return result;
}

struct Run {
  std::string results;
  std::string progress;
  double milliseconds = 0.0;
};

bool run_scan(
    const std::string& fasta,
    std::size_t sequence_count,
    std::uint32_t threads,
    Run& output,
    std::string& error) {
  const std::uint32_t handle = rdp_create();
  if (handle == 0) {
    error = "the public API could not create a context";
    return false;
  }
  const auto destroy = [&] { rdp_destroy(handle); };
  if (rdp_set_worker_threads(handle, threads) != threads) {
    error = "the threaded API did not accept " + std::to_string(threads) + " CPUs";
    destroy();
    return false;
  }
  if (rdp_load_alignment(
          handle,
          reinterpret_cast<const std::uint8_t*>(fasta.data()),
          fasta.size()) != 1) {
    error = std::string("alignment load failed: ") + rdp_get_error(handle);
    destroy();
    return false;
  }
  const std::vector<std::uint32_t> reference_groups(sequence_count, 0);
  const std::vector<std::uint8_t> flags(sequence_count, 0);
  const auto started = std::chrono::steady_clock::now();
  if (rdp_scan_begin(
          handle,
          1, 1, 0.05, 30,
          1, 70,
          1, 60,
          1, 1, 1,
          1,
          0, 0, 120, 20, 30, 0.7, 3,
          0, 0, 120, 20, 20, 100, 3,
          0,
          0,
          reference_groups.data(), reference_groups.size(),
          flags.data(), flags.size(),
          flags.data(), flags.size()) != 1) {
    error = std::string("scan start failed: ") + rdp_get_error(handle);
    destroy();
    return false;
  }
  int status = 0;
  for (std::size_t batch = 0; batch < 10000 && status != 3; ++batch) {
    status = rdp_scan_batch(handle, 4096);
    if (status < 0) {
      error = std::string("scan failed: ") + rdp_get_error(handle);
      destroy();
      return false;
    }
  }
  if (status != 3 || rdp_reconcile(handle) != 1) {
    error = std::string("scan did not reconcile: ") + rdp_get_error(handle);
    destroy();
    return false;
  }
  const auto finished = std::chrono::steady_clock::now();
  output.milliseconds = std::chrono::duration<double, std::milli>(finished - started).count();
  output.progress = rdp_get_progress_json(handle);
  output.results = rdp_get_results_json(handle);
  destroy();
  return true;
}

}  // namespace

int main() {
  constexpr std::size_t length = 30000;
  const std::string parent_one = random_sequence(length);
  const std::string parent_two = periodic_mutant(parent_one, 3);
  const std::string recombinant_one = mosaic(parent_one, parent_two, 9000, 21500);
  const std::string parent_three = random_sequence(length);
  const std::string parent_four = periodic_mutant(parent_three, 4);
  const std::string recombinant_two = mosaic(parent_three, parent_four, 4500, 16500);
  std::vector<std::string> sequences{
      parent_one,
      parent_two,
      recombinant_one,
      parent_three,
      parent_four,
      recombinant_two,
      random_sequence(length),
      random_sequence(length),
  };
  std::string fasta;
  for (std::size_t index = 0; index < sequences.size(); ++index) {
    fasta += ">multicore-" + std::to_string(index) + "\n";
    fasta += sequences[index] + "\n";
  }

  Run serial;
  Run parallel;
  std::string error;
  if (!run_scan(fasta, sequences.size(), 1, serial, error) ||
      !run_scan(fasta, sequences.size(), 4, parallel, error)) {
    std::cerr << "Multicore core check failed: " << error << '\n';
    return 1;
  }
  if (serial.results != parallel.results || serial.progress != parallel.progress) {
    std::cerr << "Multicore core check failed: 1-CPU and 4-CPU results differ\n";
    return 1;
  }
  if (parallel.results.find("\"discoveryMethods\":[\"RDP\",\"GENECONV\",\"MAXCHI\",\"CHIMAERA\",\"3SEQ\"]") ==
      std::string::npos) {
    std::cerr << "Multicore core check failed: the heavy-method fixture was not exercised\n";
    return 1;
  }
  std::cout << "Deterministic multicore dispatch verified: exact 1-CPU/4-CPU progress and results; "
            << "serial " << serial.milliseconds << " ms, parallel "
            << parallel.milliseconds << " ms (timings informational).\n";
  return 0;
}
