#include "phylpro.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace rdp {
namespace {

std::size_t vb_round_half(std::size_t value) {
  const std::size_t quotient = value / 2;
  if ((value & 1U) == 0U) return quotient;
  // VB6 CInt uses round-to-nearest-even. For an exact x.5 half, increment
  // only when the truncated integer is odd.
  return quotient + (quotient & 1U);
}

bool disabled_at(const std::vector<std::uint8_t>& disabled, std::size_t sequence) {
  return sequence < disabled.size() && disabled[sequence] != 0;
}

bool mismatch(
    const Alignment& alignment,
    std::size_t first,
    std::size_t second,
    std::size_t position) {
  const std::uint8_t a = alignment.at(first, position);
  const std::uint8_t b = alignment.at(second, position);
  return a != 0 && b != 0 && a != b;
}

float source_pearson(
    const std::vector<float>& left,
    const std::vector<float>& right,
    std::size_t target_context_index,
    bool include_self) {
  double sum_x = 0.0;
  double sum_y = 0.0;
  double sum_xy = 0.0;
  double sum_x2 = 0.0;
  double sum_y2 = 0.0;
  std::size_t observations = 0;
  for (std::size_t index = 0; index < left.size(); ++index) {
    if (!include_self && index == target_context_index) continue;
    const double x = static_cast<double>(left[index]);
    const double y = static_cast<double>(right[index]);
    sum_x += x;
    sum_y += y;
    sum_xy += x * y;
    sum_x2 += x * x;
    sum_y2 += y * y;
    ++observations;
  }

  if (observations == 0 || sum_x2 <= 0.0 || sum_y2 <= 0.0) return 1.0F;
  const double count = static_cast<double>(observations);
  const double numerator = count * sum_xy - sum_x * sum_y;
  const double left_variance = count * sum_x2 - sum_x * sum_x;
  const double right_variance = count * sum_y2 - sum_y * sum_y;
  if (left_variance <= 0.0 || right_variance <= 0.0) return 1.0F;
  const double denominator = std::sqrt(left_variance) * std::sqrt(right_variance);
  if (!(denominator > 0.0) || !std::isfinite(denominator)) return 1.0F;
  const double correlation = numerator / denominator;
  if (!std::isfinite(correlation)) return 1.0F;
  // PPRegression narrows AA / BB directly to Single; do not add a clamp that
  // could hide a source-rounding excursion just outside [-1, 1].
  return static_cast<float>(correlation);
}

std::vector<std::size_t> eligible_columns(
    const Alignment& alignment,
    const std::vector<std::size_t>& context,
    PhylproGapMode gap_mode) {
  std::vector<std::size_t> columns;
  columns.reserve(alignment.length);
  for (std::size_t position = 0; position < alignment.length; ++position) {
    std::array<bool, 5> seen{};
    bool missing = false;
    std::size_t observed_states = 0;
    for (const std::size_t sequence : context) {
      const std::uint8_t state = alignment.at(sequence, position);
      if (state == 0) {
        missing = true;
        continue;
      }
      if (!seen[state]) {
        seen[state] = true;
        ++observed_states;
      }
    }
    if (gap_mode == PhylproGapMode::strip_columns && missing) continue;
    // The supplied active FindSubSeqPP mapping retains polymorphic columns
    // after its selected gap policy. This differs from the manual's prose
    // saying that the conceptual method queries all alignment columns; the
    // literal calculation path is retained here and reported in the API.
    if (observed_states >= 2) columns.push_back(position);
  }
  return columns;
}

struct TargetWindows {
  std::size_t context_index = 0;
  std::vector<float> left;
  std::vector<float> right;
};

void add_position(
    const Alignment& alignment,
    std::size_t target,
    const std::vector<std::size_t>& context,
    std::size_t position,
    float direction,
    std::vector<float>& distances) {
  for (std::size_t index = 0; index < context.size(); ++index) {
    if (mismatch(alignment, target, context[index], position)) {
      distances[index] += direction;
    }
  }
}

}  // namespace

PhylproProfile phylpro_profile(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& target_sequences,
    const std::vector<std::uint8_t>& disabled,
    const PhylproOptions& options,
    std::string& error) {
  error.clear();
  PhylproProfile result;
  result.requested_window_sites = options.window_sites;
  result.circular = options.circular;
  result.include_self = options.include_self;
  result.gap_mode = options.gap_mode;

  if (alignment.length == 0 || alignment.sequence_count() < 3) {
    error = "PHYLPRO requires at least three aligned sequences.";
    return result;
  }
  if (options.window_sites < 2) {
    error = "The PHYLPRO window must contain at least two sites.";
    return result;
  }
  std::array<std::uint32_t, 3> sorted_targets = target_sequences;
  std::sort(sorted_targets.begin(), sorted_targets.end());
  if (std::adjacent_find(sorted_targets.begin(), sorted_targets.end()) != sorted_targets.end()) {
    error = "PHYLPRO requires three distinct role sequences.";
    return result;
  }
  for (const std::uint32_t target : target_sequences) {
    if (target >= alignment.sequence_count()) {
      error = "A PHYLPRO role sequence is outside the alignment.";
      return result;
    }
    if (disabled_at(disabled, target)) {
      error = "Disabled sequences cannot be PHYLPRO role representatives.";
      return result;
    }
  }

  std::vector<std::size_t> context;
  context.reserve(alignment.sequence_count());
  for (std::size_t sequence = 0; sequence < alignment.sequence_count(); ++sequence) {
    if (!disabled_at(disabled, sequence)) context.push_back(sequence);
  }
  if (context.size() < 3) {
    error = "PHYLPRO requires at least three non-disabled context sequences.";
    return result;
  }
  result.context_sequences = context.size();

  std::array<TargetWindows, 3> targets;
  for (std::size_t role = 0; role < target_sequences.size(); ++role) {
    const auto found = std::find(context.begin(), context.end(), target_sequences[role]);
    if (found == context.end()) {
      error = "A PHYLPRO role sequence is missing from the context.";
      return result;
    }
    targets[role].context_index = static_cast<std::size_t>(found - context.begin());
    targets[role].left.assign(context.size(), 0.0F);
    targets[role].right.assign(context.size(), 0.0F);
  }

  result.eligible_alignment_positions = eligible_columns(alignment, context, options.gap_mode);
  result.eligible_columns = result.eligible_alignment_positions.size();
  if (result.eligible_columns < 2) {
    error = "PHYLPRO found fewer than two eligible polymorphic columns.";
    return result;
  }

  const std::size_t requested_half = std::max<std::size_t>(1, vb_round_half(options.window_sites));
  const std::size_t maximum_half = options.circular
      ? std::max<std::size_t>(1, vb_round_half(result.eligible_columns))
      : std::max<std::size_t>(1, result.eligible_columns / 2);
  result.half_window_sites = std::min(requested_half, maximum_half);
  result.window_capped = result.half_window_sites != requested_half;
  const std::size_t half = result.half_window_sites;
  const std::size_t column_count = result.eligible_columns;

  auto populate = [&](std::size_t left_begin, std::size_t right_begin) {
    for (std::size_t offset = 0; offset < half; ++offset) {
      const std::size_t left_position = result.eligible_alignment_positions[left_begin + offset];
      const std::size_t right_position = result.eligible_alignment_positions[right_begin + offset];
      for (std::size_t role = 0; role < targets.size(); ++role) {
        add_position(
            alignment,
            target_sequences[role],
            context,
            left_position,
            1.0F,
            targets[role].left);
        add_position(
            alignment,
            target_sequences[role],
            context,
            right_position,
            1.0F,
            targets[role].right);
      }
    }
  };

  result.target_context_comparisons =
      3 * context.size() * (2 * half);
  if (options.circular) {
    populate(column_count - half, 0);
  } else {
    populate(0, half);
  }

  result.minimum_value = 1.0F;
  result.maximum_value = -1.0F;
  result.minimum_correlations.fill(1.0F);
  const auto emit_point = [&](std::size_t column_index, PhylproProfile& profile) {
    PhylproPoint point;
    point.alignment_position = profile.eligible_alignment_positions[column_index] + 1;
    for (std::size_t role = 0; role < targets.size(); ++role) {
      point.correlations[role] = source_pearson(
          targets[role].left,
          targets[role].right,
          targets[role].context_index,
          options.include_self);
      if (profile.points.empty() ||
          point.correlations[role] < profile.minimum_correlations[role]) {
        profile.minimum_correlations[role] = point.correlations[role];
        profile.minimum_positions[role] = point.alignment_position;
      }
      profile.minimum_value = std::min(profile.minimum_value, point.correlations[role]);
      profile.maximum_value = std::max(profile.maximum_value, point.correlations[role]);
    }
    profile.points.push_back(point);
  };

  if (options.circular) {
    result.points.reserve(column_count);
    for (std::size_t center = 0; center < column_count; ++center) {
      emit_point(center, result);
      if (center + 1 == column_count) break;
      const std::size_t old_left = (center + column_count - half) % column_count;
      const std::size_t boundary = center;
      const std::size_t new_right = (center + half) % column_count;
      for (std::size_t role = 0; role < targets.size(); ++role) {
        add_position(
            alignment,
            target_sequences[role],
            context,
            result.eligible_alignment_positions[old_left],
            -1.0F,
            targets[role].left);
        add_position(
            alignment,
            target_sequences[role],
            context,
            result.eligible_alignment_positions[boundary],
            1.0F,
            targets[role].left);
        add_position(
            alignment,
            target_sequences[role],
            context,
            result.eligible_alignment_positions[boundary],
            -1.0F,
            targets[role].right);
        add_position(
            alignment,
            target_sequences[role],
            context,
            result.eligible_alignment_positions[new_right],
            1.0F,
            targets[role].right);
      }
      result.rolling_updates += 3 * context.size() * 4;
    }
  } else {
    const std::size_t first_center = half;
    const std::size_t last_center = column_count - half;
    result.points.reserve(last_center >= first_center ? last_center - first_center + 1 : 0);
    for (std::size_t center = first_center; center <= last_center; ++center) {
      const std::size_t plotted_column = std::min(center, column_count - 1);
      emit_point(plotted_column, result);
      if (center == last_center) break;
      const std::size_t old_left = center - half;
      const std::size_t boundary = center;
      const std::size_t new_right = center + half;
      for (std::size_t role = 0; role < targets.size(); ++role) {
        add_position(
            alignment,
            target_sequences[role],
            context,
            result.eligible_alignment_positions[old_left],
            -1.0F,
            targets[role].left);
        add_position(
            alignment,
            target_sequences[role],
            context,
            result.eligible_alignment_positions[boundary],
            1.0F,
            targets[role].left);
        add_position(
            alignment,
            target_sequences[role],
            context,
            result.eligible_alignment_positions[boundary],
            -1.0F,
            targets[role].right);
        add_position(
            alignment,
            target_sequences[role],
            context,
            result.eligible_alignment_positions[new_right],
            1.0F,
            targets[role].right);
      }
      result.rolling_updates += 3 * context.size() * 4;
    }
  }

  if (result.points.empty()) {
    error = "PHYLPRO could not form a complete pair of half-windows.";
    return result;
  }
  return result;
}

}  // namespace rdp
