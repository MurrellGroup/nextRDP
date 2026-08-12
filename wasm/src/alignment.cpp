#include "alignment.hpp"

#include "json.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>

namespace rdp {
namespace {

using Records = std::vector<std::pair<std::string, std::string>>;

std::string trim(std::string_view value) {
  const auto first = value.find_first_not_of(" \t\n\r");
  if (first == std::string_view::npos) return {};
  const auto last = value.find_last_not_of(" \t\n\r");
  return std::string(value.substr(first, last - first + 1));
}

std::string lowercase(std::string_view value) {
  std::string result(value);
  std::transform(result.begin(), result.end(), result.begin(), [](unsigned char character) {
    return static_cast<char>(std::tolower(character));
  });
  return result;
}

bool starts_with_case_insensitive(std::string_view value, std::string_view prefix) {
  if (value.size() < prefix.size()) return false;
  for (std::size_t index = 0; index < prefix.size(); ++index) {
    if (std::tolower(static_cast<unsigned char>(value[index])) !=
        std::tolower(static_cast<unsigned char>(prefix[index]))) {
      return false;
    }
  }
  return true;
}

std::vector<std::string> lines(std::string_view input) {
  std::vector<std::string> result;
  std::string line;
  line.reserve(256);
  for (const char character : input) {
    if (character == '\r') continue;
    if (character == '\n') {
      result.push_back(std::move(line));
      line.clear();
    } else {
      line.push_back(character);
    }
  }
  result.push_back(std::move(line));
  return result;
}

std::string sequence_fragment(std::string_view value) {
  std::string result;
  result.reserve(value.size());
  for (const unsigned char character : value) {
    if (std::isalpha(character)) {
      result.push_back(static_cast<char>(std::toupper(character)));
    } else if (character == '-' || character == '.' || character == '?' ||
               character == '!' || character == '*') {
      result.push_back(static_cast<char>(character));
    }
  }
  return result;
}

std::vector<std::string> words(std::string_view line) {
  std::istringstream input{std::string(line)};
  std::vector<std::string> result;
  for (std::string word; input >> word;) result.push_back(std::move(word));
  return result;
}

void append_record(
    Records& records,
    std::unordered_map<std::string, std::size_t>& indices,
    std::string name,
    std::string fragment) {
  if (name.empty() || fragment.empty()) return;
  const auto found = indices.find(name);
  if (found == indices.end()) {
    indices.emplace(name, records.size());
    records.emplace_back(std::move(name), std::move(fragment));
  } else {
    records[found->second].second += fragment;
  }
}

Records parse_fasta_or_gde(const std::vector<std::string>& input_lines) {
  Records records;
  std::string name;
  std::string sequence;
  auto flush = [&]() {
    if (!name.empty()) records.emplace_back(std::move(name), std::move(sequence));
    name.clear();
    sequence.clear();
  };

  for (const auto& raw_line : input_lines) {
    const std::string line = trim(raw_line);
    if (line.empty()) continue;
    if (line.front() == '>' || line.front() == '%') {
      flush();
      name = trim(std::string_view(line).substr(1));
      if (const auto whitespace = name.find_first_of(" \t"); whitespace != std::string::npos) {
        name = name.substr(0, whitespace);
      }
    } else if (!name.empty() && line.front() != ';') {
      sequence += sequence_fragment(line);
    }
  }
  flush();
  return records;
}

Records parse_clustal(const std::vector<std::string>& input_lines) {
  Records records;
  std::unordered_map<std::string, std::size_t> indices;
  for (std::size_t line_index = 1; line_index < input_lines.size(); ++line_index) {
    const auto& raw = input_lines[line_index];
    const std::string line = trim(raw);
    if (line.empty()) continue;
    if (raw.size() > 0 && std::isspace(static_cast<unsigned char>(raw.front())) &&
        line.find_first_not_of("*:. \t") == std::string::npos) {
      continue;
    }
    const auto tokens = words(line);
    if (tokens.size() < 2) continue;
    const std::string fragment = sequence_fragment(tokens[1]);
    append_record(records, indices, tokens[0], fragment);
  }
  return records;
}

Records parse_mega(const std::vector<std::string>& input_lines) {
  Records records;
  std::unordered_map<std::string, std::size_t> indices;
  std::string current_name;
  for (const auto& raw : input_lines) {
    const std::string line = trim(raw);
    if (line.empty() || line.front() == '!') continue;
    if (line.front() == '#') {
      if (starts_with_case_insensitive(line, "#mega") ||
          starts_with_case_insensitive(line, "#title") ||
          starts_with_case_insensitive(line, "#format")) {
        continue;
      }
      const auto space = line.find_first_of(" \t");
      current_name = trim(std::string_view(line).substr(1, space == std::string::npos
                                                              ? std::string::npos
                                                              : space - 1));
      if (space != std::string::npos) {
        append_record(
            records,
            indices,
            current_name,
            sequence_fragment(std::string_view(line).substr(space + 1)));
      }
    } else if (!current_name.empty()) {
      append_record(records, indices, current_name, sequence_fragment(line));
    }
  }
  return records;
}

std::pair<std::string, std::string> nexus_record(std::string_view line) {
  const std::string clean = trim(line);
  if (clean.empty()) return {};
  if (clean.front() == '\'' || clean.front() == '"') {
    const char quote = clean.front();
    const auto end = clean.find(quote, 1);
    if (end == std::string::npos) return {};
    return {
        clean.substr(1, end - 1),
        sequence_fragment(std::string_view(clean).substr(end + 1)),
    };
  }
  const auto tokens = words(clean);
  if (tokens.size() < 2) return {};
  return {tokens[0], sequence_fragment(tokens[1])};
}

Records parse_nexus(const std::vector<std::string>& input_lines) {
  Records records;
  std::unordered_map<std::string, std::size_t> indices;
  bool in_matrix = false;
  for (const auto& raw : input_lines) {
    std::string line = trim(raw);
    if (line.empty() || starts_with_case_insensitive(line, "[")) continue;
    if (!in_matrix) {
      const std::string low = lowercase(line);
      const auto matrix = low.find("matrix");
      if (matrix == std::string::npos) continue;
      in_matrix = true;
      line = trim(std::string_view(line).substr(matrix + 6));
    }
    const auto semicolon = line.find(';');
    const std::string row = semicolon == std::string::npos ? line : line.substr(0, semicolon);
    const auto [name, fragment] = nexus_record(row);
    append_record(records, indices, name, fragment);
    if (semicolon != std::string::npos) break;
  }
  return records;
}

std::pair<std::string, std::string> phylip_named_row(std::string_view raw) {
  const std::string line = trim(raw);
  const auto tokens = words(line);
  if (tokens.size() >= 2) {
    std::string fragment;
    for (std::size_t index = 1; index < tokens.size(); ++index) {
      fragment += sequence_fragment(tokens[index]);
    }
    return {tokens[0], std::move(fragment)};
  }
  // Strict PHYLIP reserves the first ten columns for the sequence name.
  if (raw.size() > 10) {
    const std::string name = trim(raw.substr(0, 10));
    const std::string fragment = sequence_fragment(raw.substr(10));
    if (!name.empty() && !fragment.empty()) return {name, fragment};
  }
  return {};
}

bool phylip_complete(
    const Records& records,
    std::size_t expected_sequences,
    std::size_t expected_length) {
  if (records.size() != expected_sequences) return false;
  return std::all_of(records.begin(), records.end(), [&](const auto& record) {
    return record.second.size() == expected_length;
  });
}

Records parse_phylip_interleaved(
    const std::vector<std::string>& input_lines,
    std::size_t expected_sequences,
    std::size_t expected_length) {
  Records records;
  records.reserve(expected_sequences);
  std::unordered_map<std::string, std::size_t> indices;
  std::size_t line_index = 1;
  while (line_index < input_lines.size() && records.size() < expected_sequences) {
    const auto& raw = input_lines[line_index++];
    if (trim(raw).empty()) continue;
    auto [name, fragment] = phylip_named_row(raw);
    if (name.empty() || fragment.empty() || fragment.size() > expected_length) return {};
    if (!indices.emplace(name, records.size()).second) return {};
    records.emplace_back(std::move(name), std::move(fragment));
  }
  if (records.size() != expected_sequences) return {};

  std::size_t continuation_index = 0;
  for (; line_index < input_lines.size(); ++line_index) {
    const std::string line = trim(input_lines[line_index]);
    if (line.empty()) {
      continuation_index = 0;
      continue;
    }
    const auto tokens = words(line);
    if (tokens.empty()) continue;
    auto named = tokens.size() > 1 ? indices.find(tokens[0]) : indices.end();
    std::size_t target = continuation_index % expected_sequences;
    std::string fragment;
    if (named != indices.end()) {
      target = named->second;
      for (std::size_t index = 1; index < tokens.size(); ++index) {
        fragment += sequence_fragment(tokens[index]);
      }
    } else {
      fragment = sequence_fragment(line);
    }
    records[target].second += fragment;
    if (records[target].second.size() > expected_length) return {};
    ++continuation_index;
  }
  return phylip_complete(records, expected_sequences, expected_length)
      ? records
      : Records{};
}

Records parse_phylip_sequential(
    const std::vector<std::string>& input_lines,
    std::size_t expected_sequences,
    std::size_t expected_length) {
  Records records;
  records.reserve(expected_sequences);
  std::size_t line_index = 1;
  while (records.size() < expected_sequences) {
    while (line_index < input_lines.size() && trim(input_lines[line_index]).empty()) ++line_index;
    if (line_index >= input_lines.size()) return {};
    auto [name, sequence] = phylip_named_row(input_lines[line_index++]);
    if (name.empty() || sequence.empty() || sequence.size() > expected_length) return {};
    while (sequence.size() < expected_length) {
      while (line_index < input_lines.size() && trim(input_lines[line_index]).empty()) ++line_index;
      if (line_index >= input_lines.size()) return {};
      sequence += sequence_fragment(trim(input_lines[line_index++]));
      if (sequence.size() > expected_length) return {};
    }
    records.emplace_back(std::move(name), std::move(sequence));
  }
  return phylip_complete(records, expected_sequences, expected_length)
      ? records
      : Records{};
}

Records parse_phylip(const std::vector<std::string>& input_lines) {
  if (input_lines.empty()) return {};
  const auto header = words(input_lines.front());
  if (header.size() < 2) return {};
  std::size_t expected_sequences = 0;
  std::size_t expected_length = 0;
  try {
    expected_sequences = static_cast<std::size_t>(std::stoull(header[0]));
    expected_length = static_cast<std::size_t>(std::stoull(header[1]));
  } catch (...) {
    return {};
  }
  if (expected_sequences == 0 || expected_length == 0) return {};

  Records records = parse_phylip_interleaved(
      input_lines, expected_sequences, expected_length);
  if (!records.empty()) return records;
  return parse_phylip_sequential(input_lines, expected_sequences, expected_length);
}

std::uint8_t encode_base(char base) {
  switch (static_cast<char>(std::toupper(static_cast<unsigned char>(base)))) {
    case 'A': return 1;
    case 'C': return 2;
    case 'G': return 3;
    case 'T':
    case 'U': return 4;
    default: return 0;
  }
}

std::uint64_t choose_three(std::uint64_t count) {
  if (count < 3) return 0;
  return count * (count - 1) * (count - 2) / 6;
}

void make_suggested_mask(Alignment& alignment) {
  const std::size_t count = alignment.sequence_count();
  alignment.suggested_mask.assign(count, 0);
  if (count <= 3) return;

  // The supplied VB workflow keeps three representatives when no diversity exists.
  if (alignment.minimum_pair_identity >= 1.0 - 1e-7) {
    std::fill(alignment.suggested_mask.begin() + 3, alignment.suggested_mask.end(), 1);
    return;
  }

  std::vector<float> closest_identity(count, -1.0F);
  std::vector<std::size_t> closest_sequence(count, count);
  const auto active = [&](std::size_t index) {
    return alignment.suggested_mask[index] == 0;
  };
  const auto refresh_closest = [&](std::size_t first) {
    closest_identity[first] = -1.0F;
    closest_sequence[first] = count;
    if (!active(first)) return;
    for (std::size_t second = 0; second < count; ++second) {
      if (first == second || !active(second)) continue;
      const float identity = alignment.similarity(first, second);
      if (identity > closest_identity[first]) {
        closest_identity[first] = identity;
        closest_sequence[first] = second;
      }
    }
  };
  for (std::size_t index = 0; index < count; ++index) refresh_closest(index);

  std::size_t active_count = count;
  while (active_count >= 4) {
    // This is the threshold calculation used by AutoMaskmnu_Click in the supplied
    // VB source, including the RDP window opportunity factor.
    const double correction = static_cast<double>(choose_three(active_count)) *
        (static_cast<double>(alignment.length) / 30.0);
    if (!(correction > 0.05)) break;
    const double minimum_distance =
        std::log(correction / 0.05) / std::log(4.0) /
        static_cast<double>(alignment.length);

    float best_identity = -1.0F;
    std::size_t first = count;
    std::size_t second = count;
    for (std::size_t index = 0; index < count; ++index) {
      if (!active(index) || closest_identity[index] <= best_identity) continue;
      best_identity = closest_identity[index];
      first = index;
      second = closest_sequence[index];
    }
    if (first == count || second == count || !active(second) ||
        1.0 - static_cast<double>(best_identity) >= minimum_distance) {
      break;
    }

    std::size_t remove = first;
    const double first_length =
        static_cast<double>(alignment.sequence_summaries[first].valid_sites);
    const double second_length =
        static_cast<double>(alignment.sequence_summaries[second].valid_sites);
    if (second_length > first_length * 0.95 && first_length > second_length * 0.95) {
      double first_total_distance = 0.0;
      double second_total_distance = 0.0;
      for (std::size_t other = 0; other < count; ++other) {
        first_total_distance += 1.0 - alignment.similarity(first, other);
        second_total_distance += 1.0 - alignment.similarity(second, other);
      }
      // Preserve the representative that is, on average, more distinct.
      remove = first_total_distance > second_total_distance ? second : first;
    } else {
      // Preserve the representative with more usable alignment columns.
      remove = first_length > second_length ? second : first;
    }

    alignment.suggested_mask[remove] = 1;
    --active_count;
    closest_identity[remove] = -1.0F;
    closest_sequence[remove] = count;
    for (std::size_t index = 0; index < count; ++index) {
      if (active(index) && closest_sequence[index] == remove) refresh_closest(index);
    }
  }
}

AlignmentParseResult finalise(std::string format, Records records) {
  AlignmentParseResult result;
  if (records.size() < 3) {
    result.error = "RDP requires at least three aligned nucleotide sequences.";
    return result;
  }

  const std::size_t length = records.front().second.size();
  if (length == 0) {
    result.error = "The alignment contains no nucleotide columns.";
    return result;
  }
  for (const auto& [name, sequence] : records) {
    if (sequence.size() != length) {
      std::ostringstream message;
      message << "The sequences are not aligned: '" << name << "' has " << sequence.size()
              << " columns; expected " << length << ".";
      result.error = message.str();
      return result;
    }
  }

  Alignment alignment;
  alignment.format = std::move(format);
  alignment.length = length;
  alignment.names.reserve(records.size());
  alignment.sequences.reserve(records.size());
  alignment.states.resize(records.size() * length);
  alignment.sequence_summaries.resize(records.size());

  std::unordered_map<std::string, std::size_t> name_counts;
  for (std::size_t sequence_index = 0; sequence_index < records.size(); ++sequence_index) {
    auto [name, sequence] = std::move(records[sequence_index]);
    if (name.empty()) name = "Sequence_" + std::to_string(sequence_index + 1);
    const std::size_t occurrence = ++name_counts[name];
    if (occurrence > 1) {
      alignment.warnings.push_back("Duplicate sequence name retained: " + name);
    }
    alignment.names.push_back(std::move(name));
    alignment.sequences.push_back(std::move(sequence));
    auto& summary = alignment.sequence_summaries[sequence_index];
    for (std::size_t position = 0; position < length; ++position) {
      const char original = alignment.sequences.back()[position];
      const std::uint8_t state = encode_base(original);
      alignment.states[sequence_index * length + position] = state;
      if (state == 0) ++summary.missing_sites;
      else ++summary.valid_sites;
      if (sequence_index == 0 && original == '!') {
        alignment.partition_boundaries.push_back(position + 1);
      }
    }
  }

  if (!alignment.partition_boundaries.empty()) {
    alignment.warnings.push_back(
        "Concatenation markers were retained as missing states; partition-aware scans are not yet enabled in this snapshot.");
  }

  std::array<std::size_t, 5> counts{};
  for (std::size_t position = 0; position < length; ++position) {
    counts.fill(0);
    for (std::size_t sequence = 0; sequence < alignment.sequence_count(); ++sequence) {
      ++counts[alignment.at(sequence, position)];
    }
    std::size_t observed_states = 0;
    std::size_t repeated_states = 0;
    for (std::size_t state = 1; state <= 4; ++state) {
      if (counts[state] > 0) ++observed_states;
      if (counts[state] >= 2) ++repeated_states;
    }
    if (observed_states > 1) ++alignment.variable_site_count;
    if (repeated_states >= 2) ++alignment.informative_site_count;
  }

  const std::size_t sequence_count = alignment.sequence_count();
  alignment.pair_similarity.assign(sequence_count * sequence_count, 1.0F);
  double identity_sum = 0.0;
  std::size_t pair_count = 0;
  alignment.minimum_pair_identity = 1.0;
  for (std::size_t first = 0; first < sequence_count; ++first) {
    for (std::size_t second = first + 1; second < sequence_count; ++second) {
      std::size_t valid = 0;
      std::size_t same = 0;
      for (std::size_t position = 0; position < length; ++position) {
        const std::uint8_t a = alignment.at(first, position);
        const std::uint8_t b = alignment.at(second, position);
        if (a == 0 || b == 0) continue;
        ++valid;
        if (a == b) ++same;
      }
      const float identity = valid == 0 ? 0.0F : static_cast<float>(same) / static_cast<float>(valid);
      alignment.pair_similarity[first * sequence_count + second] = identity;
      alignment.pair_similarity[second * sequence_count + first] = identity;
      alignment.minimum_pair_identity = std::min(alignment.minimum_pair_identity, static_cast<double>(identity));
      identity_sum += identity;
      ++pair_count;
    }
  }
  alignment.mean_pair_identity = pair_count == 0 ? -1.0 : identity_sum / static_cast<double>(pair_count);
  if (alignment.minimum_pair_identity >= 0.0 && alignment.minimum_pair_identity < 0.6) {
    alignment.warnings.push_back(
        "At least one sequence pair is below 60% identity; the RDP5 manual warns that such alignments are especially vulnerable to false signals.");
  }
  for (std::size_t index = 0; index < sequence_count; ++index) {
    if (alignment.sequence_summaries[index].valid_sites == 0) {
      alignment.warnings.push_back("Sequence has no unambiguous nucleotide sites: " + alignment.names[index]);
    }
  }

  make_suggested_mask(alignment);
  const auto auto_masked = static_cast<std::size_t>(std::count(
      alignment.suggested_mask.begin(), alignment.suggested_mask.end(), 1));
  if (auto_masked > 0) {
    alignment.warnings.push_back(
        std::to_string(auto_masked) +
        " near-identical sequence(s) were auto-masked using the supplied RDP5 optimisation workflow; you can change the selection before scanning.");
  }

  result.alignment = std::move(alignment);
  return result;
}

}  // namespace

AlignmentParseResult parse_alignment(std::string_view input) {
  if (input.empty()) return {{}, "The selected file is empty."};
  const auto input_lines = lines(input);
  std::string first;
  std::size_t first_line_index = 0;
  for (; first_line_index < input_lines.size(); ++first_line_index) {
    first = trim(input_lines[first_line_index]);
    if (!first.empty()) break;
  }
  if (first.empty()) return {{}, "The selected file contains no readable text."};
  const std::vector<std::string> content_lines(
      input_lines.begin() + static_cast<std::ptrdiff_t>(first_line_index), input_lines.end());

  Records records;
  std::string format;
  if (first.front() == '>' || first.front() == '%') {
    format = first.front() == '>' ? "FASTA" : "GDE";
    records = parse_fasta_or_gde(content_lines);
  } else if (starts_with_case_insensitive(first, "clustal") ||
             starts_with_case_insensitive(first, "muscle")) {
    format = "CLUSTAL";
    records = parse_clustal(content_lines);
  } else if (starts_with_case_insensitive(first, "#nexus")) {
    format = "NEXUS";
    records = parse_nexus(content_lines);
  } else if (starts_with_case_insensitive(first, "#mega")) {
    format = "MEGA";
    records = parse_mega(content_lines);
  } else {
    const auto header = words(first);
    const auto is_decimal = [](const std::string& value) {
      return !value.empty() &&
          std::all_of(value.begin(), value.end(), [](unsigned char character) {
            return std::isdigit(character) != 0;
          });
    };
    const bool looks_phylip = header.size() >= 2 &&
        is_decimal(header[0]) && is_decimal(header[1]);
    if (looks_phylip) {
      format = "PHYLIP";
      records = parse_phylip(content_lines);
    }
  }

  if (records.empty()) {
    return {
        {},
        "The alignment format was not recognised. This snapshot accepts FASTA, GDE, CLUSTAL, PHYLIP, NEXUS, and MEGA text alignments.",
    };
  }
  return finalise(std::move(format), std::move(records));
}

AlignmentParseResult build_alignment(
    std::string format,
    std::vector<std::string> names,
    std::vector<std::string> sequences) {
  if (names.size() != sequences.size()) {
    return {{}, "The saved project has different numbers of sequence names and sequences."};
  }
  Records records;
  records.reserve(names.size());
  for (std::size_t index = 0; index < names.size(); ++index) {
    records.emplace_back(std::move(names[index]), std::move(sequences[index]));
  }
  return finalise(std::move(format), std::move(records));
}

std::string alignment_summary_json(const Alignment& alignment, const std::vector<std::uint8_t>& mask) {
  const std::size_t count = alignment.sequence_count();
  const bool explicit_mask = mask.size() == count;
  const auto is_masked = [&](std::size_t index) {
    return explicit_mask
        ? mask[index] != 0
        : index < alignment.suggested_mask.size() && alignment.suggested_mask[index] != 0;
  };
  std::size_t active = 0;
  for (std::size_t index = 0; index < count; ++index) {
    if (!is_masked(index)) ++active;
  }
  const double recommended_distance = alignment.length == 0
      ? 0.0
      : (2.0 * std::log(4.0 * static_cast<double>(std::max<std::size_t>(1, active)))) /
          static_cast<double>(alignment.length);

  std::ostringstream out;
  out << '{';
  out << "\"format\":";
  json::string(out, alignment.format);
  out << ",\"sequenceCount\":" << count;
  out << ",\"alignmentLength\":" << alignment.length;
  out << ",\"activeSequenceCount\":" << active;
  out << ",\"tripletCount\":" << choose_three(active);
  out << ",\"variableSiteCount\":" << alignment.variable_site_count;
  out << ",\"informativeSiteCount\":" << alignment.informative_site_count;
  out << ",\"minimumPairIdentity\":";
  json::number(out, alignment.minimum_pair_identity);
  out << ",\"meanPairIdentity\":";
  json::number(out, alignment.mean_pair_identity);
  out << ",\"recommendedMinimumDistance\":";
  json::number(out, recommended_distance);
  out << ",\"partitionBoundaries\":[";
  for (std::size_t index = 0; index < alignment.partition_boundaries.size(); ++index) {
    if (index) out << ',';
    out << alignment.partition_boundaries[index];
  }
  out << "],\"sequences\":[";
  for (std::size_t index = 0; index < count; ++index) {
    if (index) out << ',';
    const auto& summary = alignment.sequence_summaries[index];
    out << "{\"index\":" << index << ",\"name\":";
    json::string(out, alignment.names[index]);
    out << ",\"validSites\":" << summary.valid_sites;
    out << ",\"missingSites\":" << summary.missing_sites;
    out << ",\"missingFraction\":";
    json::number(out, alignment.length == 0
                          ? 0.0
                          : static_cast<double>(summary.missing_sites) /
                                static_cast<double>(alignment.length));
    out << ",\"masked\":" << (is_masked(index) ? "true" : "false") << '}';
  }
  out << "],\"warnings\":[";
  for (std::size_t index = 0; index < alignment.warnings.size(); ++index) {
    if (index) out << ',';
    json::string(out, alignment.warnings[index]);
  }
  out << "]}";
  return out.str();
}

}  // namespace rdp
