#include "burt_confidence.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

namespace rdp {
namespace {

constexpr std::size_t kStates = 3;
constexpr std::size_t kSymbols = 3;
constexpr std::size_t kHmmCyclesArgument = 20;
constexpr std::size_t kHmmTrainingStarts = kHmmCyclesArgument + 1;
constexpr std::size_t kMaximumTrainingIterations = 100;
constexpr double kPseudocount = 0.01;
constexpr std::uint32_t kSourceDefaultSeed = 3;
constexpr double kNegativeLatticeSentinel = -10000000000000000.0;

struct CandidateInterval {
  std::array<std::int64_t, 5> values{};
};

struct HmmOutput {
  std::vector<CandidateInterval> intervals;
  double best_log_likelihood = 0.0;
  bool trained = false;
};

// The supplied DNA5 DLL is a Windows build and calls the Microsoft C runtime
// srand/rand pair. Keep its 15-bit sequence local so Emscripten's libc does
// not silently change BURT results across browsers or toolchain releases.
class SourceMsvcRandom {
 public:
  explicit SourceMsvcRandom(std::uint32_t seed) : state_(seed) {}

  double unit() {
    state_ = state_ * 214013U + 2531011U;
    const std::uint32_t value = (state_ >> 16U) & 0x7fffU;
    return static_cast<double>(value) / 32767.0;
  }

 private:
  std::uint32_t state_;
};

std::int64_t source_vb_long(double value) {
  if (!std::isfinite(value)) return 0;
  const double lower = std::floor(value);
  const double fraction = value - lower;
  const auto lower_integer = static_cast<std::int64_t>(lower);
  if (fraction < 0.5) return lower_integer;
  if (fraction > 0.5) return lower_integer + 1;
  return lower_integer % 2 == 0 ? lower_integer : lower_integer + 1;
}

std::int64_t absolute_source_coordinate(std::int64_t value) {
  return value < 0 ? -value : value;
}

std::size_t bounded_coordinate(std::int64_t value, std::size_t length) {
  if (length == 0) return 0;
  value = absolute_source_coordinate(value);
  if (value < 1) return 1;
  if (value > static_cast<std::int64_t>(length)) return length;
  return static_cast<std::size_t>(value);
}

bool information_rich_state(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& sorted,
    std::size_t coordinate,
    std::uint8_t& symbol) {
  const std::size_t position = coordinate - 1;
  const std::uint8_t first = alignment.at(sorted[0], position);
  const std::uint8_t second = alignment.at(sorted[1], position);
  const std::uint8_t third = alignment.at(sorted[2], position);
  if (first == 0 || second == 0 || third == 0) return false;
  if (first == second && first != third) {
    symbol = 0;
    return true;
  }
  if (second == third && first != second) {
    symbol = 1;
    return true;
  }
  if (first == third && first != second) {
    symbol = 2;
    return true;
  }
  return false;
}

double source_log_sum(const std::array<double, kStates>& values) {
  double maximum = kNegativeLatticeSentinel;
  for (const double value : values) {
    if (maximum < value) maximum = value;
  }
  double total = 0.0;
  for (const double value : values) total += std::exp(value - maximum);
  return std::log(total) + maximum;
}

double source_viterbi(
    const std::vector<std::uint8_t>& symbols,
    std::size_t sequence_length,
    const std::array<double, 9>& transition,
    const std::array<double, 9>& emission,
    const std::array<double, 3>& initial,
    BurtConfidenceWorkspace& workspace) {
  const std::size_t stride = sequence_length + 1;

  for (std::size_t state = 0; state < kStates; ++state) {
    workspace.lattice[state * stride] =
        emission[symbols[0] + state * kSymbols] + initial[state];
  }
  std::array<double, 9> alternatives{};
  for (std::size_t position = 1; position <= sequence_length; ++position) {
    for (std::size_t previous = 0; previous < kStates; ++previous) {
      for (std::size_t next = 0; next < kStates; ++next) {
        alternatives[previous + next * kStates] =
            workspace.lattice[position - 1 + previous * stride] +
            transition[previous + next * kStates] +
            emission[symbols[position] + next * kSymbols];
      }
    }
    for (std::size_t next = 0; next < kStates; ++next) {
      const std::size_t offset = position + next * stride;
      workspace.lattice[offset] = kNegativeLatticeSentinel;
      for (std::size_t previous = 0; previous < kStates; ++previous) {
        const double candidate = alternatives[previous + next * kStates];
        if (workspace.lattice[offset] < candidate) {
          workspace.lattice[offset] = candidate;
          workspace.backpointers[offset] = static_cast<double>(previous);
        }
      }
    }
  }

  // Preserve GetLaticePathP's terminal assignment: it stores the selected
  // terminal state's predecessor rather than the terminal state itself.
  double maximum_likelihood = kNegativeLatticeSentinel;
  for (std::size_t state = 0; state < kStates; ++state) {
    const std::size_t offset = sequence_length + state * stride;
    if (workspace.lattice[offset] > maximum_likelihood) {
      maximum_likelihood = workspace.lattice[offset];
      workspace.path[sequence_length] =
          static_cast<std::int32_t>(workspace.backpointers[offset]);
    }
  }
  for (std::size_t position = sequence_length; position-- > 0;) {
    const std::size_t state = static_cast<std::size_t>(
        std::clamp<std::int32_t>(workspace.path[position + 1], 0, 2));
    workspace.path[position] = static_cast<std::int32_t>(
        workspace.backpointers[position + state * stride]);
  }
  return maximum_likelihood;
}

HmmOutput source_train_hmm(
    const std::vector<std::uint8_t>& symbols,
    std::size_t sequence_length,
    std::size_t alignment_length,
    BurtConfidenceWorkspace& workspace) {
  HmmOutput output;
  if (symbols.size() <= sequence_length || alignment_length == 0) return output;

  std::array<double, 9> transition{};
  std::array<double, 9> emission{};
  std::array<double, 3> initial{};
  std::array<double, 9> best_transition{};
  std::array<double, 9> best_emission{};
  std::array<double, 3> best_initial{};
  std::array<double, 9> transition_counts{};
  std::array<double, 9> state_counts{};
  SourceMsvcRandom random(kSourceDefaultSeed);
  workspace.best_path.clear();
  // DoHMMCyclesSerial callocates these buffers once, outside both its start
  // and training-iteration loops. Clear once per event to preserve that
  // lifetime while avoiding an O(profile length) fill on every Viterbi pass.
  const std::size_t training_stride = sequence_length + 1;
  workspace.lattice.assign(training_stride * kStates, 0.0);
  workspace.backpointers.assign(training_stride * kStates, 0.0);
  workspace.path.assign(sequence_length + 4, 0);

  double best_likelihood = -1000000.0;
  double path_maximum = 0.0;
  for (std::size_t start = 0; start < kHmmTrainingStarts; ++start) {
    const float transition_mass =
        static_cast<float>(5.0 / static_cast<double>(alignment_length));
    const float self_transition_mass = 1.0F - transition_mass;
    for (std::size_t previous = 0; previous < kStates; ++previous) {
      for (std::size_t next = 0; next < kStates; ++next) {
        transition[previous + next * kStates] = previous == next
            ? std::log(static_cast<double>(self_transition_mass))
            : std::log(
                static_cast<double>(transition_mass) /
                static_cast<double>(kStates - 1));
      }
    }

    const double imbalance =
        static_cast<double>(static_cast<int>(kSymbols * random.unit() + 1.0)) /
        10.0;
    std::array<std::uint8_t, kSymbols + 1> used_imbalances{};
    for (std::size_t state = 0; state < kStates; ++state) {
      while (true) {
        const int selected =
            static_cast<int>((kStates + 3) * random.unit()) - 2;
        if (selected < 0 || selected >= static_cast<int>(kStates) ||
            used_imbalances[static_cast<std::size_t>(selected)] != 0) {
          continue;
        }
        for (std::size_t symbol = 0; symbol < kSymbols; ++symbol) {
          emission[symbol + state * kSymbols] =
              symbol == static_cast<std::size_t>(selected)
              ? 1.0 / static_cast<double>(kSymbols) + imbalance * 2.0
              : 1.0 / static_cast<double>(kSymbols) - imbalance;
        }
        used_imbalances[static_cast<std::size_t>(selected)] = 1;
        break;
      }
    }
    for (double& value : emission) value = std::log(value);
    for (double& value : initial) value = std::log(1.0 / kStates);

    double likelihood = kNegativeLatticeSentinel;
    for (std::size_t iteration = 1;
         iteration <= kMaximumTrainingIterations;
         ++iteration) {
      likelihood = source_viterbi(
          symbols,
          sequence_length,
          transition,
          emission,
          initial,
          workspace);
      if (path_maximum == likelihood) break;
      path_maximum = likelihood;

      transition_counts.fill(0.0);
      state_counts.fill(0.0);
      for (std::size_t position = 0; position < sequence_length; ++position) {
        const std::size_t previous = static_cast<std::size_t>(
            std::clamp<std::int32_t>(workspace.path[position], 0, 2));
        const std::size_t next = static_cast<std::size_t>(
            std::clamp<std::int32_t>(workspace.path[position + 1], 0, 2));
        transition_counts[previous + next * kStates] += 1.0;
      }
      for (std::size_t position = 0; position <= sequence_length; ++position) {
        const std::size_t state = static_cast<std::size_t>(
            std::clamp<std::int32_t>(workspace.path[position], 0, 2));
        state_counts[symbols[position] + state * kSymbols] += 1.0;
      }
      for (std::size_t previous = 0; previous < kStates; ++previous) {
        double total = 0.0;
        for (std::size_t next = 0; next < kStates; ++next) {
          total += transition_counts[previous + next * kStates];
        }
        total += kPseudocount * kStates;
        for (std::size_t next = 0; next < kStates; ++next) {
          const std::size_t index = previous + next * kStates;
          transition_counts[index] += kPseudocount;
          transition[index] = std::log(transition_counts[index] / total);
        }
      }
      for (std::size_t state = 0; state < kStates; ++state) {
        double total = 0.0;
        for (std::size_t symbol = 0; symbol < kSymbols; ++symbol) {
          total += state_counts[symbol + state * kSymbols];
        }
        total += kPseudocount * kSymbols;
        for (std::size_t symbol = 0; symbol < kSymbols; ++symbol) {
          emission[symbol + state * kSymbols] = std::log(
              (state_counts[symbol + state * kSymbols] + kPseudocount) /
              total);
        }
      }
    }

    if (likelihood > best_likelihood) {
      best_likelihood = likelihood;
      best_transition = transition;
      best_emission = emission;
      best_initial = initial;
      workspace.best_path.assign(
          workspace.path.begin(),
          workspace.path.begin() + static_cast<std::ptrdiff_t>(sequence_length + 1));
    }
  }

  if (workspace.best_path.size() < sequence_length + 1 ||
      best_likelihood <= -1000000.0) {
    return output;
  }
  workspace.path = workspace.best_path;
  output.best_log_likelihood = static_cast<double>(
      static_cast<float>(best_likelihood));
  output.trained = true;

  const std::size_t stride = sequence_length + 1;
  workspace.forward.assign(stride * kStates, 0.0);
  workspace.reverse.assign(stride * kStates, 0.0);
  for (std::size_t state = 0; state < kStates; ++state) {
    workspace.forward[state * stride] =
        best_emission[symbols[0] + state * kSymbols] + best_initial[state];
  }
  for (std::size_t position = 1; position <= sequence_length; ++position) {
    for (std::size_t next = 0; next < kStates; ++next) {
      std::array<double, kStates> alternatives{};
      for (std::size_t previous = 0; previous < kStates; ++previous) {
        alternatives[previous] =
            workspace.forward[position - 1 + previous * stride] +
            best_transition[previous + next * kStates] +
            best_emission[symbols[position] + next * kSymbols];
      }
      workspace.forward[position + next * stride] = source_log_sum(alternatives);
    }
  }
  for (std::size_t state = 0; state < kStates; ++state) {
    workspace.reverse[sequence_length + state * stride] = 0.0;
  }
  for (std::size_t position = sequence_length; position-- > 0;) {
    for (std::size_t previous = 0; previous < kStates; ++previous) {
      std::array<double, kStates> alternatives{};
      for (std::size_t next = 0; next < kStates; ++next) {
        alternatives[next] =
            workspace.reverse[position + 1 + next * stride] +
            best_transition[previous + next * kStates] +
            best_emission[symbols[position + 1] + next * kSymbols];
      }
      workspace.reverse[position + previous * stride] = source_log_sum(alternatives);
    }
  }

  workspace.posterior.assign(stride * kStates, 0.0);
  for (std::size_t position = 0; position <= sequence_length; ++position) {
    std::array<double, kStates> combined{};
    double maximum = -1e22;
    for (std::size_t state = 0; state < kStates; ++state) {
      combined[state] = workspace.forward[position + state * stride] +
          workspace.reverse[position + state * stride];
      if (maximum < combined[state]) maximum = combined[state];
    }
    double total = 0.0;
    for (double& value : combined) {
      value = std::exp(value - maximum);
      total += value;
    }
    for (std::size_t state = 0; state < kStates; ++state) {
      workspace.posterior[position * kStates + state] = combined[state] / total;
    }
  }
  return output;
}

std::int64_t source_xdiffpos(
    const BurtConfidenceWorkspace& workspace,
    std::int64_t index) {
  if (index <= 0 ||
      index > static_cast<std::int64_t>(workspace.information_coordinates.size())) {
    return 0;
  }
  return static_cast<std::int64_t>(
      workspace.information_coordinates[static_cast<std::size_t>(index - 1)]);
}

bool any_state_exceeds(
    const BurtConfidenceWorkspace& workspace,
    std::size_t position,
    double threshold) {
  for (std::size_t state = 0; state < kStates; ++state) {
    if (workspace.posterior[position * kStates + state] > threshold) return true;
  }
  return false;
}

void repair_interval_midpoint(
    CandidateInterval& interval,
    std::size_t alignment_length) {
  const std::int64_t left = interval.values[3];
  const std::int64_t right = interval.values[4];
  if (left < right) {
    if (!(interval.values[2] > left && interval.values[2] < right)) {
      interval.values[2] = left + source_vb_long(
          static_cast<double>(right - left) / 2.0);
    }
  } else if (!(interval.values[2] > left || interval.values[2] < right)) {
    interval.values[2] = source_vb_long(
        static_cast<double>(
            left + static_cast<std::int64_t>(alignment_length) - right) /
        2.0);
    if (interval.values[2] > static_cast<std::int64_t>(alignment_length)) {
      interval.values[2] -= static_cast<std::int64_t>(alignment_length);
    }
  }
}

HmmOutput source_ben_hmm(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& representatives,
    bool circular,
    BurtConfidenceWorkspace& workspace) {
  HmmOutput output;
  std::array<std::uint32_t, 3> sorted = representatives;
  std::sort(sorted.begin(), sorted.end());
  workspace.symbols.clear();
  workspace.information_coordinates.clear();
  workspace.position_information_counts.assign(alignment.length + 1, 0);
  workspace.symbols.reserve(alignment.length + 1);
  workspace.information_coordinates.reserve(alignment.length);

  std::size_t information_count = 0;
  for (std::size_t coordinate = 1; coordinate <= alignment.length; ++coordinate) {
    std::uint8_t symbol = 0;
    if (information_rich_state(alignment, sorted, coordinate, symbol)) {
      workspace.symbols.push_back(symbol);
      workspace.information_coordinates.push_back(coordinate);
      ++information_count;
    }
    workspace.position_information_counts[coordinate] = information_count;
  }

  // BenHMM passes SLen=Y while RecodeB's real symbols occupy 0..Y-1. Its
  // zero-initialized element Y is consequently a deliberate sentinel symbol.
  workspace.symbols.resize(information_count + 1, 0);
  if (information_count == 0) return output;
  std::size_t hmm_length = information_count;
  std::size_t circular_offset = 0;
  const std::vector<std::uint8_t>* training_symbols = &workspace.symbols;
  if (circular) {
    circular_offset = static_cast<std::size_t>(source_vb_long(
        static_cast<double>(information_count) / 2.0));
    hmm_length = information_count + circular_offset * 2;
    workspace.expanded_symbols.assign(hmm_length + 1, 0);
    for (std::size_t index = 0; index < circular_offset; ++index) {
      workspace.expanded_symbols[
          information_count + index + circular_offset + 1] =
          workspace.symbols[index];
    }
    for (std::size_t index = circular_offset;
         index < information_count;
         ++index) {
      workspace.expanded_symbols[index - circular_offset] = workspace.symbols[index];
    }
    for (std::size_t index = 0; index < information_count; ++index) {
      workspace.expanded_symbols[index + circular_offset + 1] = workspace.symbols[index];
    }
    training_symbols = &workspace.expanded_symbols;
  }

  output = source_train_hmm(
      *training_symbols,
      hmm_length,
      alignment.length,
      workspace);
  if (!output.trained) return output;

  std::size_t selected_length = hmm_length;
  if (circular && hmm_length > 0) {
    selected_length = static_cast<std::size_t>(source_vb_long(
        static_cast<double>(hmm_length) / 2.0));
    for (std::size_t position = 0; position <= selected_length; ++position) {
      const std::size_t source_position = position + circular_offset + 1;
      if (source_position >= workspace.path.size()) break;
      workspace.path[position] = workspace.path[source_position];
      for (std::size_t state = 0; state < kStates; ++state) {
        workspace.posterior[position * kStates + state] =
            workspace.posterior[source_position * kStates + state];
      }
    }
  }

  if (selected_length == 0) return output;
  if (selected_length >= 3 && workspace.path.size() > selected_length) {
    workspace.path[selected_length - 1] = circular
        ? workspace.path[2]
        : workspace.path[selected_length - 2];
  }

  if (selected_length < 4) return output;
  for (std::size_t position = 2; position <= selected_length - 2; ++position) {
    if (workspace.path[position] == workspace.path[position + 1]) continue;
    CandidateInterval interval;
    if (position < selected_length - 2) {
      const std::int64_t current = source_xdiffpos(
          workspace, static_cast<std::int64_t>(position));
      const std::int64_t next = source_xdiffpos(
          workspace, static_cast<std::int64_t>(position + 1));
      interval.values[2] = source_vb_long(
          static_cast<double>(current) +
          static_cast<double>(next - current) / 2.0);
    } else {
      const std::int64_t current = source_xdiffpos(
          workspace, static_cast<std::int64_t>(position));
      interval.values[2] = source_vb_long(
          static_cast<double>(current) +
          static_cast<double>(
              static_cast<std::int64_t>(alignment.length) - current) /
              2.0);
    }
    if (interval.values[2] < 0) interval.values[2] = 0;

    for (std::int64_t search = static_cast<std::int64_t>(position);
         search >= -static_cast<std::int64_t>(selected_length);
         --search) {
      std::int64_t mapped = search;
      if (search <= 0) {
        if (circular) mapped = search + static_cast<std::int64_t>(selected_length);
        else {
          interval.values[0] = 1;
          break;
        }
      }
      if (mapped < 0 || mapped > static_cast<std::int64_t>(selected_length)) continue;
      const std::size_t state_position = static_cast<std::size_t>(mapped);
      if (any_state_exceeds(workspace, state_position, 0.995) &&
          interval.values[3] == 0) {
        interval.values[3] = mapped > 1
            ? source_xdiffpos(workspace, mapped - 1) + 1
            : 1;
      }
      if (any_state_exceeds(workspace, state_position, 0.999)) {
        interval.values[0] = mapped > 1
            ? source_xdiffpos(workspace, mapped - 1) + 1
            : 1;
        break;
      }
    }

    for (std::int64_t search = static_cast<std::int64_t>(position + 1);
         search <= static_cast<std::int64_t>(selected_length * 2);
         ++search) {
      std::int64_t mapped = search;
      if (search > static_cast<std::int64_t>(selected_length)) {
        if (circular) mapped = search - static_cast<std::int64_t>(selected_length);
        else {
          interval.values[1] = static_cast<std::int64_t>(alignment.length);
          interval.values[4] = static_cast<std::int64_t>(alignment.length);
          break;
        }
      }
      if (mapped < 0 || mapped > static_cast<std::int64_t>(selected_length)) continue;
      const std::size_t state_position = static_cast<std::size_t>(mapped);
      if (any_state_exceeds(workspace, state_position, 0.995) &&
          interval.values[4] == 0) {
        interval.values[4] = mapped < static_cast<std::int64_t>(selected_length - 1)
            ? source_xdiffpos(workspace, mapped + 1) - 1
            : static_cast<std::int64_t>(alignment.length);
      }
      if (any_state_exceeds(workspace, state_position, 0.999)) {
        interval.values[1] = mapped < static_cast<std::int64_t>(selected_length - 1)
            ? source_xdiffpos(workspace, mapped + 1) - 1
            : static_cast<std::int64_t>(alignment.length);
        break;
      }
    }

    if (circular) {
      for (const std::size_t field : {std::size_t{0}, std::size_t{1},
                                      std::size_t{3}, std::size_t{4}}) {
        if (interval.values[field] != 1 &&
            interval.values[field] != static_cast<std::int64_t>(alignment.length)) {
          continue;
        }
        if (field == 1 || field == 4) {
          interval.values[field] = source_xdiffpos(workspace, 1) - 1;
          if (interval.values[field] < 1) interval.values[field] = 1;
        } else {
          interval.values[field] = source_xdiffpos(
              workspace, static_cast<std::int64_t>(selected_length)) + 1;
          if (interval.values[field] > static_cast<std::int64_t>(alignment.length)) {
            interval.values[field] = static_cast<std::int64_t>(alignment.length);
          }
        }
      }
    }
    repair_interval_midpoint(interval, alignment.length);
    output.intervals.push_back(interval);
  }
  return output;
}

std::array<std::int64_t, 5> source_match_breakpoint(
    std::vector<CandidateInterval>& intervals,
    std::size_t breakpoint,
    std::size_t alignment_length,
    bool circular,
    const BurtConfidenceWorkspace& workspace) {
  std::array<std::int64_t, 5> selected{};
  if (intervals.empty() || breakpoint < 1 || breakpoint > alignment_length) return selected;

  std::vector<std::int64_t> distances(intervals.size(), 0);
  const std::int64_t total_information = static_cast<std::int64_t>(
      workspace.position_information_counts[alignment_length]);
  for (std::size_t index = 0; index < intervals.size(); ++index) {
    auto& interval = intervals[index];
    const std::int64_t factor = alignment_length == 0
        ? 0
        : interval.values[0] / static_cast<std::int64_t>(alignment_length);
    if (factor > 0) {
      for (std::int64_t& value : interval.values) {
        if (value > static_cast<std::int64_t>(alignment_length)) {
          value -= static_cast<std::int64_t>(alignment_length) * factor;
        }
      }
    }
    const std::int64_t left = interval.values[0];
    const std::int64_t right = interval.values[1];
    const std::int64_t point = interval.values[2];
    const std::int64_t breakpoint_information = static_cast<std::int64_t>(
        workspace.position_information_counts[breakpoint]);
    const std::int64_t point_information = point <= 0
        ? 0
        : static_cast<std::int64_t>(
            workspace.position_information_counts[
                bounded_coordinate(point, alignment_length)]);
    std::int64_t distance = absolute_source_coordinate(
        breakpoint_information - point_information);
    if (circular && total_information - distance < distance) {
      distance = total_information - distance;
    }
    const std::int64_t breakpoint_signed = static_cast<std::int64_t>(breakpoint);
    const bool contains = left < right
        ? breakpoint_signed >= left && breakpoint_signed <= right
        : breakpoint_signed >= left || breakpoint_signed <= right;
    if (contains) {
      if (total_information - distance < distance) {
        distance = -(total_information - distance);
      } else {
        distance = -distance;
      }
    }
    distances[index] = distance;
  }

  std::int64_t smallest = static_cast<std::int64_t>(alignment_length);
  std::size_t selected_index = 0;
  bool found = false;
  for (std::size_t index = 0; index < distances.size(); ++index) {
    if (smallest == 0) break;
    if (smallest > 0) {
      if (smallest > distances[index]) {
        smallest = distances[index];
        selected_index = index;
        found = true;
      }
    } else if (absolute_source_coordinate(smallest) >
                   absolute_source_coordinate(distances[index]) &&
               distances[index] <= 0) {
      smallest = distances[index];
      selected_index = index;
      found = true;
    }
  }
  if (!found) return selected;
  selected = intervals[selected_index].values;
  if (smallest > 0) {
    for (std::int64_t& value : selected) value = -value;
  }
  return selected;
}

bool triplet_missing_at(
    const std::vector<std::uint8_t>& triplet_missing,
    std::size_t alignment_length,
    std::size_t coordinate) {
  if (coordinate < 1 || coordinate > alignment_length) return true;
  return coordinate < triplet_missing.size() && triplet_missing[coordinate] != 0;
}

bool triplet_difference_at(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& representatives,
    std::size_t coordinate) {
  if (coordinate < 1 || coordinate > alignment.length) return false;
  const std::uint8_t first = alignment.at(representatives[0], coordinate - 1);
  const std::uint8_t second = alignment.at(representatives[1], coordinate - 1);
  const std::uint8_t third = alignment.at(representatives[2], coordinate - 1);
  return first != 0 && second != 0 && third != 0 &&
      (first != second || first != third);
}

void source_missing_data_reposition(
    const Alignment& alignment,
    const std::vector<std::uint8_t>& triplet_missing,
    bool circular,
    std::array<std::int64_t, 10>& confidence,
    std::size_t& beginning,
    std::size_t& ending,
    std::array<bool, 2>& adjusted) {
  const std::size_t length = alignment.length;
  if (confidence[0] > -1) {
    std::int64_t coordinate = confidence[1];
    for (std::size_t guard = 0; guard <= length * 2 + 2; ++guard) {
      --coordinate;
      if (coordinate < 1) {
        if (circular) coordinate = static_cast<std::int64_t>(length);
        else {
          coordinate = 2;
          beginning = bounded_coordinate(coordinate, length);
          adjusted[0] = true;
          break;
        }
      }
      if (coordinate == confidence[2]) {
        beginning = bounded_coordinate(confidence[2], length);
        break;
      }
      if (!triplet_missing_at(
              triplet_missing,
              length,
              static_cast<std::size_t>(coordinate))) {
        continue;
      }
      if (confidence[2] < confidence[5]) {
        if (coordinate < confidence[5]) {
          confidence[2] = coordinate + 1;
          adjusted[0] = true;
          break;
        }
        beginning = bounded_coordinate(confidence[2], length);
        break;
      }
      confidence[2] = coordinate + 1;
      adjusted[0] = true;
      break;
    }
  }

  if (confidence[3] > -1) {
    std::int64_t coordinate = confidence[3];
    for (std::size_t guard = 0; guard <= length * 2 + 2; ++guard) {
      ++coordinate;
      if (coordinate > static_cast<std::int64_t>(length)) coordinate = 1;
      if (coordinate == confidence[5]) {
        ending = bounded_coordinate(confidence[5], length);
        break;
      }
      if (!triplet_missing_at(
              triplet_missing,
              length,
              static_cast<std::size_t>(coordinate))) {
        continue;
      }
      if (confidence[2] < confidence[5]) {
        if (coordinate > confidence[2]) {
          confidence[5] = coordinate - 1;
          ending = bounded_coordinate(confidence[5], length);
          adjusted[1] = true;
          break;
        }
        ending = bounded_coordinate(confidence[5], length);
        break;
      }
      confidence[5] = coordinate - 1;
      ending = bounded_coordinate(confidence[5], length);
      adjusted[1] = true;
      break;
    }
  }

  if (confidence[1] > -1 &&
      triplet_missing_at(triplet_missing, length, beginning)) {
    bool found = false;
    for (std::size_t coordinate = beginning; coordinate <= length; ++coordinate) {
      if (!triplet_missing_at(triplet_missing, length, coordinate)) {
        beginning = coordinate;
        confidence[0] = static_cast<std::int64_t>(coordinate);
        adjusted[0] = true;
        found = true;
        break;
      }
    }
    if (!found) {
      for (std::size_t coordinate = 1; coordinate <= beginning; ++coordinate) {
        if (!triplet_missing_at(triplet_missing, length, coordinate)) {
          beginning = coordinate;
          confidence[0] = static_cast<std::int64_t>(coordinate);
          adjusted[0] = true;
          break;
        }
      }
    }
  }

  if (confidence[5] > -1 &&
      triplet_missing_at(triplet_missing, length, ending)) {
    bool found = false;
    for (std::size_t coordinate = ending; coordinate >= 1; --coordinate) {
      if (!triplet_missing_at(triplet_missing, length, coordinate)) {
        ending = coordinate;
        confidence[4] = static_cast<std::int64_t>(coordinate);
        adjusted[1] = true;
        found = true;
        break;
      }
      if (coordinate == 1) break;
    }
    if (!found) {
      for (std::size_t coordinate = length; coordinate >= 1; --coordinate) {
        if (!triplet_missing_at(triplet_missing, length, coordinate)) {
          ending = coordinate;
          confidence[4] = static_cast<std::int64_t>(coordinate);
          adjusted[1] = true;
          break;
        }
        if (coordinate == 1) break;
      }
    }
  }
}

void source_final_gap_reposition(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& representatives,
    bool circular,
    std::size_t& beginning,
    std::size_t& ending,
    std::array<bool, 2>& adjusted) {
  const auto any_gap = [&](std::size_t coordinate) {
    return alignment.at(representatives[0], coordinate - 1) == 0 ||
        alignment.at(representatives[1], coordinate - 1) == 0 ||
        alignment.at(representatives[2], coordinate - 1) == 0;
  };
  const auto source_repeated_second_valid = [&](std::size_t coordinate) {
    return alignment.at(representatives[0], coordinate - 1) != 0 &&
        alignment.at(representatives[1], coordinate - 1) != 0 &&
        alignment.at(representatives[1], coordinate - 1) != 0;
  };
  if (any_gap(beginning)) {
    std::size_t coordinate = beginning + 1;
    std::size_t cycles = 0;
    for (std::size_t guard = 0; guard <= alignment.length * 2 + 2; ++guard) {
      if (coordinate > alignment.length) {
        if (!circular) break;
        coordinate = 1;
        if (++cycles > 1) break;
      }
      if (source_repeated_second_valid(coordinate)) {
        beginning = coordinate;
        adjusted[0] = true;
        break;
      }
      ++coordinate;
    }
  }

  if (any_gap(ending)) {
    std::int64_t coordinate = static_cast<std::int64_t>(ending) - 1;
    std::size_t cycles = 0;
    for (std::size_t guard = 0; guard <= alignment.length * 2 + 2; ++guard) {
      // Preserve the active source's literal X > 0 branch. In a linear
      // analysis it exits immediately; circular data tests the final site.
      if (coordinate > 0) {
        if (circular) {
          coordinate = static_cast<std::int64_t>(alignment.length);
          if (++cycles > 1) break;
        } else {
          break;
        }
      }
      if (coordinate < 1) {
        if (circular) coordinate += static_cast<std::int64_t>(alignment.length);
        else {
          ending = 1;
          adjusted[1] = true;
          break;
        }
      }
      if (source_repeated_second_valid(static_cast<std::size_t>(coordinate))) {
        ending = static_cast<std::size_t>(coordinate);
        adjusted[1] = true;
        break;
      }
      --coordinate;
    }
  }
}

}  // namespace

BreakpointConfidenceEvidence source_polish_breakpoints(
    const Alignment& current_alignment,
    const std::vector<std::uint8_t>& triplet_missing_data,
    const std::array<std::uint32_t, 3>& representatives,
    std::size_t beginning,
    std::size_t ending,
    bool circular,
    BurtConfidenceWorkspace& workspace) {
  BreakpointConfidenceEvidence result;
  result.attempted = true;
  result.random_seed = kSourceDefaultSeed;
  result.hmm_cycles_argument = kHmmCyclesArgument;
  result.serial_training_starts = kHmmTrainingStarts;
  result.polished_beginning = beginning;
  result.polished_ending = ending;
  result.boundaries[0].input_coordinate = beginning;
  result.boundaries[1].input_coordinate = ending;
  result.boundaries[0].polished_coordinate = beginning;
  result.boundaries[1].polished_coordinate = ending;

  if (current_alignment.length == 0 || beginning < 1 || ending < 1 ||
      beginning > current_alignment.length || ending > current_alignment.length ||
      std::any_of(
          representatives.begin(),
          representatives.end(),
          [&](std::uint32_t sequence) {
            return sequence >= current_alignment.sequence_count();
          })) {
    result.unavailable_reason = "invalid-triplet-or-breakpoint";
    return result;
  }

  HmmOutput hmm = source_ben_hmm(
      current_alignment,
      representatives,
      circular,
      workspace);
  result.information_rich_sites = workspace.information_coordinates.size();
  result.best_log_likelihood = hmm.best_log_likelihood;
  result.candidate_interval_count = hmm.intervals.size();
  if (!hmm.trained) {
    result.unavailable_reason = result.information_rich_sites == 0
        ? "no-information-rich-sites"
        : "hmm-training-unavailable";
    return result;
  }
  if (hmm.intervals.empty()) {
    result.unavailable_reason = "no-hmm-state-transition";
    return result;
  }

  std::array<std::int64_t, 10> confidence{};
  const auto beginning_match = source_match_breakpoint(
      hmm.intervals,
      beginning,
      current_alignment.length,
      circular,
      workspace);
  confidence[0] = beginning_match[0];
  confidence[1] = beginning_match[1];
  confidence[2] = beginning_match[2];
  confidence[6] = beginning_match[3];
  confidence[7] = beginning_match[4];
  const auto ending_match = source_match_breakpoint(
      hmm.intervals,
      ending,
      current_alignment.length,
      circular,
      workspace);
  confidence[3] = ending_match[0];
  confidence[4] = ending_match[1];
  confidence[5] = ending_match[2];
  confidence[8] = ending_match[3];
  confidence[9] = ending_match[4];
  const std::array<bool, 2> source_interval_contains_input{
      std::any_of(
          beginning_match.begin(), beginning_match.end(),
          [](std::int64_t value) { return value != 0; }) &&
          beginning_match[2] >= 0,
      std::any_of(
          ending_match.begin(), ending_match.end(),
          [](std::int64_t value) { return value != 0; }) &&
          ending_match[2] >= 0,
  };

  const std::size_t original_beginning = beginning;
  const std::size_t original_ending = ending;
  const double event_span = beginning < ending
      ? static_cast<double>(ending - beginning)
      : static_cast<double>(current_alignment.length + ending - beginning);
  const auto source_near_enough = [&](std::size_t input, std::int64_t candidate) {
    const std::int64_t difference = absolute_source_coordinate(
        static_cast<std::int64_t>(input) - absolute_source_coordinate(candidate));
    return static_cast<double>(difference) < event_span / 2.0 ||
        (circular &&
         static_cast<double>(
             static_cast<std::int64_t>(current_alignment.length) - difference) <
             event_span / 2.0);
  };

  if (absolute_source_coordinate(
          absolute_source_coordinate(confidence[2]) -
          absolute_source_coordinate(confidence[5])) > 2) {
    if (confidence[2] >= 0) {
      beginning = bounded_coordinate(confidence[2], current_alignment.length);
    } else if (source_near_enough(beginning, confidence[2])) {
      beginning = bounded_coordinate(confidence[2], current_alignment.length);
      for (const std::size_t index : {std::size_t{0}, std::size_t{1},
                                      std::size_t{2}, std::size_t{6},
                                      std::size_t{7}}) {
        confidence[index] = absolute_source_coordinate(confidence[index]);
      }
    }
    if (confidence[5] >= 0) {
      ending = bounded_coordinate(confidence[5], current_alignment.length);
    } else if (source_near_enough(ending, confidence[5])) {
      ending = bounded_coordinate(confidence[5], current_alignment.length);
      for (const std::size_t index : {std::size_t{3}, std::size_t{4},
                                      std::size_t{5}, std::size_t{8},
                                      std::size_t{9}}) {
        confidence[index] = absolute_source_coordinate(confidence[index]);
      }
    }
  } else {
    result.single_transition_assignment = true;
    const std::int64_t point = absolute_source_coordinate(confidence[2]);
    const std::int64_t beginning_distance = absolute_source_coordinate(
        static_cast<std::int64_t>(beginning) - point);
    const std::int64_t ending_distance = absolute_source_coordinate(
        static_cast<std::int64_t>(ending) - point);
    const bool beginning_is_closer =
        (beginning_distance <= ending_distance ||
         (circular &&
          static_cast<std::int64_t>(current_alignment.length) - beginning_distance <=
              ending_distance)) &&
        (!circular ||
         beginning_distance <=
             static_cast<std::int64_t>(current_alignment.length) - ending_distance ||
         static_cast<std::int64_t>(current_alignment.length) - beginning_distance <=
             static_cast<std::int64_t>(current_alignment.length) - ending_distance);
    if (beginning_is_closer) {
      if (source_near_enough(beginning, confidence[2])) {
        beginning = bounded_coordinate(confidence[2], current_alignment.length);
      } else {
        confidence[0] = -1;
        confidence[1] = -1;
        confidence[2] = -1;
      }
      confidence[3] = -1;
      confidence[4] = -1;
      confidence[5] = -1;
    } else {
      if (source_near_enough(ending, confidence[5])) {
        ending = bounded_coordinate(confidence[5], current_alignment.length);
      } else {
        confidence[3] = -1;
        confidence[4] = -1;
        confidence[5] = -1;
      }
      confidence[0] = -1;
      confidence[1] = -1;
      confidence[2] = -1;
    }
  }

  if (confidence[0] == 0 && confidence[1] == 0 && confidence[2] == 0) {
    confidence[0] = -1;
    confidence[1] = -1;
    confidence[2] = -1;
    beginning = original_beginning;
  }
  if (confidence[3] == 0 && confidence[4] == 0 && confidence[5] == 0) {
    confidence[3] = -1;
    confidence[4] = -1;
    confidence[5] = -1;
    ending = original_ending;
  }

  std::array<bool, 2> missing_adjusted{};
  source_missing_data_reposition(
      current_alignment,
      triplet_missing_data,
      circular,
      confidence,
      beginning,
      ending,
      missing_adjusted);

  std::size_t outside_differences = 0;
  std::size_t inside_differences = 0;
  if (beginning < ending) {
    for (std::size_t coordinate = 1; coordinate < beginning; ++coordinate) {
      if (triplet_difference_at(current_alignment, representatives, coordinate)) {
        ++outside_differences;
      }
    }
    for (std::size_t coordinate = ending + 1;
         coordinate <= current_alignment.length;
         ++coordinate) {
      if (triplet_difference_at(current_alignment, representatives, coordinate)) {
        ++outside_differences;
      }
    }
    for (std::size_t coordinate = beginning; coordinate <= ending; ++coordinate) {
      if (triplet_difference_at(current_alignment, representatives, coordinate)) {
        ++inside_differences;
      }
    }
  } else {
    for (std::size_t coordinate = 1; coordinate <= ending; ++coordinate) {
      if (triplet_difference_at(current_alignment, representatives, coordinate)) {
        ++inside_differences;
      }
    }
    for (std::size_t coordinate = ending + 1; coordinate < beginning; ++coordinate) {
      if (triplet_difference_at(current_alignment, representatives, coordinate)) {
        ++outside_differences;
      }
    }
    for (std::size_t coordinate = beginning;
         coordinate <= current_alignment.length;
         ++coordinate) {
      if (triplet_difference_at(current_alignment, representatives, coordinate)) {
        ++inside_differences;
      }
    }
  }
  if (outside_differences < 3 || inside_differences < 3) {
    beginning = original_beginning;
    ending = original_ending;
    result.insufficient_inside_or_outside_reverted = true;
  }

  std::array<bool, 2> final_gap_adjusted{};
  source_final_gap_reposition(
      current_alignment,
      representatives,
      circular,
      beginning,
      ending,
      final_gap_adjusted);

  const std::array<std::array<std::size_t, 5>, 2> mapping{{
      {{0, 1, 2, 6, 7}},
      {{3, 4, 5, 8, 9}},
  }};
  for (std::size_t boundary_index = 0; boundary_index < 2; ++boundary_index) {
    auto& boundary = result.boundaries[boundary_index];
    const auto& fields = mapping[boundary_index];
    boundary.confidence_99_beginning = confidence[fields[0]];
    boundary.confidence_99_ending = confidence[fields[1]];
    boundary.hmm_coordinate = confidence[fields[2]];
    boundary.confidence_95_beginning = confidence[fields[3]];
    boundary.confidence_95_ending = confidence[fields[4]];
    boundary.polished_coordinate = boundary_index == 0 ? beginning : ending;
    boundary.interval_available = confidence[fields[2]] != -1 &&
        !(confidence[fields[0]] == -1 && confidence[fields[1]] == -1);
    boundary.source_interval_contains_input =
        boundary.interval_available && source_interval_contains_input[boundary_index];
    boundary.confidence_99_wraps_origin = boundary.interval_available &&
        absolute_source_coordinate(confidence[fields[0]]) >
            absolute_source_coordinate(confidence[fields[1]]);
    boundary.confidence_95_wraps_origin = boundary.interval_available &&
        absolute_source_coordinate(confidence[fields[3]]) >
            absolute_source_coordinate(confidence[fields[4]]);
    boundary.repositioned =
        boundary.polished_coordinate != boundary.input_coordinate;
    boundary.missing_data_adjusted = missing_adjusted[boundary_index];
    boundary.final_gap_adjusted = final_gap_adjusted[boundary_index];
  }
  result.polished_beginning = beginning;
  result.polished_ending = ending;
  result.available = result.boundaries[0].interval_available ||
      result.boundaries[1].interval_available;
  result.unavailable_reason = result.available ? "" : "no-matched-breakpoint-interval";
  return result;
}

}  // namespace rdp
