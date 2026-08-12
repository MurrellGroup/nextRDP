#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace rdp {

struct SequenceSummary {
  std::size_t valid_sites = 0;
  std::size_t missing_sites = 0;
};

struct Alignment {
  std::string format;
  std::vector<std::string> names;
  std::vector<std::string> sequences;
  std::vector<std::uint8_t> states;
  std::vector<float> pair_similarity;
  std::vector<SequenceSummary> sequence_summaries;
  std::vector<std::uint8_t> suggested_mask;
  std::vector<std::size_t> partition_boundaries;
  std::vector<std::string> warnings;
  std::size_t length = 0;
  std::size_t variable_site_count = 0;
  std::size_t informative_site_count = 0;
  double minimum_pair_identity = -1.0;
  double mean_pair_identity = -1.0;

  [[nodiscard]] std::size_t sequence_count() const { return names.size(); }
  [[nodiscard]] std::uint8_t at(std::size_t sequence, std::size_t position) const {
    return states[sequence * length + position];
  }
  [[nodiscard]] float similarity(std::size_t first, std::size_t second) const {
    return pair_similarity[first * sequence_count() + second];
  }
};

struct AlignmentParseResult {
  Alignment alignment;
  std::string error;
  [[nodiscard]] bool ok() const { return error.empty(); }
};

AlignmentParseResult parse_alignment(std::string_view input);
AlignmentParseResult build_alignment(
    std::string format,
    std::vector<std::string> names,
    std::vector<std::string> sequences);
std::string alignment_summary_json(const Alignment& alignment, const std::vector<std::uint8_t>& mask = {});

}  // namespace rdp
