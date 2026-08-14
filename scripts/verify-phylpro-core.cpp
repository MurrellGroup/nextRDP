#include "alignment.hpp"
#include "phylpro.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {

[[noreturn]] void fail(const std::string& message) {
  std::cerr << "PHYLPRO core verification failed: " << message << '\n';
  std::exit(1);
}

void require(bool condition, const std::string& message) {
  if (!condition) fail(message);
}

std::size_t vb_round_half(std::size_t value) {
  const std::size_t quotient = value / 2;
  return quotient + ((value & 1U) != 0U && (quotient & 1U) != 0U ? 1U : 0U);
}

bool disabled_at(const std::vector<std::uint8_t>& disabled, std::size_t sequence) {
  return sequence < disabled.size() && disabled[sequence] != 0;
}

std::vector<std::size_t> reference_columns(
    const rdp::Alignment& alignment,
    const std::vector<std::size_t>& context,
    rdp::PhylproGapMode gap_mode) {
  std::vector<std::size_t> columns;
  for (std::size_t position = 0; position < alignment.length; ++position) {
    std::array<bool, 5> seen{};
    bool missing = false;
    std::size_t states = 0;
    for (const std::size_t sequence : context) {
      const std::uint8_t state = alignment.at(sequence, position);
      if (state == 0) {
        missing = true;
      } else if (!seen[state]) {
        seen[state] = true;
        ++states;
      }
    }
    if (gap_mode == rdp::PhylproGapMode::strip_columns && missing) continue;
    if (states >= 2) columns.push_back(position);
  }
  return columns;
}

float reference_pearson(
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
    const double x = left[index];
    const double y = right[index];
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
  return static_cast<float>(correlation);
}

std::vector<rdp::PhylproPoint> brute_force_profile(
    const rdp::Alignment& alignment,
    const std::array<std::uint32_t, 3>& targets,
    const std::vector<std::uint8_t>& disabled,
    const rdp::PhylproOptions& options) {
  std::vector<std::size_t> context;
  for (std::size_t sequence = 0; sequence < alignment.sequence_count(); ++sequence) {
    if (!disabled_at(disabled, sequence)) context.push_back(sequence);
  }
  const auto columns = reference_columns(alignment, context, options.gap_mode);
  const std::size_t requested_half = std::max<std::size_t>(1, vb_round_half(options.window_sites));
  const std::size_t maximum_half = options.circular
      ? std::max<std::size_t>(1, vb_round_half(columns.size()))
      : std::max<std::size_t>(1, columns.size() / 2);
  const std::size_t half = std::min(requested_half, maximum_half);

  std::vector<std::size_t> centers;
  if (options.circular) {
    centers.resize(columns.size());
    std::iota(centers.begin(), centers.end(), 0);
  } else {
    for (std::size_t center = half; center <= columns.size() - half; ++center) {
      centers.push_back(center);
    }
  }

  std::vector<rdp::PhylproPoint> points;
  points.reserve(centers.size());
  for (const std::size_t center : centers) {
    rdp::PhylproPoint point;
    point.alignment_position = columns[std::min(center, columns.size() - 1)] + 1;
    for (std::size_t role = 0; role < targets.size(); ++role) {
      const auto target_it = std::find(context.begin(), context.end(), targets[role]);
      require(target_it != context.end(), "brute-force role disappeared from context");
      std::vector<float> left(context.size(), 0.0F);
      std::vector<float> right(context.size(), 0.0F);
      for (std::size_t offset = 0; offset < half; ++offset) {
        const std::size_t left_index = options.circular
            ? (center + columns.size() - half + offset) % columns.size()
            : center - half + offset;
        const std::size_t right_index = options.circular
            ? (center + offset) % columns.size()
            : center + offset;
        for (std::size_t context_index = 0; context_index < context.size(); ++context_index) {
          const std::uint8_t target_left = alignment.at(targets[role], columns[left_index]);
          const std::uint8_t context_left = alignment.at(
              context[context_index], columns[left_index]);
          const std::uint8_t target_right = alignment.at(targets[role], columns[right_index]);
          const std::uint8_t context_right = alignment.at(
              context[context_index], columns[right_index]);
          if (target_left != 0 && context_left != 0 && target_left != context_left) {
            left[context_index] += 1.0F;
          }
          if (target_right != 0 && context_right != 0 && target_right != context_right) {
            right[context_index] += 1.0F;
          }
        }
      }
      point.correlations[role] = reference_pearson(
          left,
          right,
          static_cast<std::size_t>(target_it - context.begin()),
          options.include_self);
    }
    points.push_back(point);
  }
  return points;
}

void compare_with_brute_force(
    const rdp::Alignment& alignment,
    const std::array<std::uint32_t, 3>& targets,
    const std::vector<std::uint8_t>& disabled,
    const rdp::PhylproOptions& options,
    const std::string& label) {
  std::string error;
  const auto optimized = rdp::phylpro_profile(
      alignment, targets, disabled, options, error);
  require(error.empty(), label + ": " + error);
  const auto reference = brute_force_profile(
      alignment, targets, disabled, options);
  require(optimized.points.size() == reference.size(),
          label + ": optimized/brute-force point counts differ");
  for (std::size_t point = 0; point < reference.size(); ++point) {
    require(optimized.points[point].alignment_position == reference[point].alignment_position,
            label + ": alignment-position mapping changed");
    for (std::size_t role = 0; role < targets.size(); ++role) {
      require(std::abs(
                  optimized.points[point].correlations[role] -
                  reference[point].correlations[role]) < 1e-6F,
              label + ": rolling three-target row differs from full recomputation");
    }
  }
}

char next_base(char base, std::size_t shift = 1) {
  constexpr std::string_view bases = "ACGT";
  const std::size_t index = bases.find(base);
  return bases[(index + shift) % bases.size()];
}

}  // namespace

int main() {
  constexpr std::size_t length = 240;
  constexpr std::string_view bases = "ACGT";
  std::string parent_one(length, 'A');
  for (std::size_t position = 0; position < length; ++position) {
    parent_one[position] = bases[(position * 7 + position / 5) % bases.size()];
  }
  std::string parent_two = parent_one;
  for (std::size_t position = 0; position < length; ++position) {
    parent_two[position] = next_base(parent_two[position], 1 + (position & 1U));
  }
  std::string recombinant = parent_one;
  for (std::size_t position = 80; position < 160; ++position) {
    recombinant[position] = parent_two[position];
  }
  std::string parent_one_sibling = parent_one;
  std::string parent_two_sibling = parent_two;
  for (std::size_t position = 5; position < length; position += 17) {
    parent_one_sibling[position] = next_base(parent_one_sibling[position], 2);
  }
  for (std::size_t position = 8; position < length; position += 19) {
    parent_two_sibling[position] = next_base(parent_two_sibling[position], 2);
  }
  std::string outlier(length, 'A');
  for (std::size_t position = 0; position < length; ++position) {
    outlier[position] = bases[(position * 11 + 3) % bases.size()];
  }
  std::string gapped = parent_one_sibling;
  for (std::size_t position = 3; position < length; position += 17) {
    gapped[position] = '-';
  }

  auto parsed = rdp::build_alignment(
      "synthetic-phylpro",
      {"recombinant", "parent-one", "parent-two", "parent-one-sibling",
       "parent-two-sibling", "outlier", "gapped-context"},
      {recombinant, parent_one, parent_two, parent_one_sibling,
       parent_two_sibling, outlier, gapped});
  require(parsed.ok(), parsed.error);
  const std::array<std::uint32_t, 3> targets{{0, 1, 2}};
  std::vector<std::uint8_t> enabled(parsed.alignment.sequence_count(), 0);

  rdp::PhylproOptions circular;
  circular.window_sites = 40;
  circular.circular = true;
  circular.gap_mode = rdp::PhylproGapMode::ignore_pairwise;
  compare_with_brute_force(
      parsed.alignment, targets, enabled, circular, "circular source windows");

  std::string error;
  const auto profile = rdp::phylpro_profile(
      parsed.alignment, targets, enabled, circular, error);
  require(error.empty(), error);
  require(profile.eligible_columns == length && profile.points.size() == length,
          "ignore-missing-pairwise did not retain all polymorphic columns");
  require(profile.context_sequences == 7 &&
              profile.target_context_comparisons == 3 * 7 * 40 &&
              profile.rolling_updates == (length - 1) * 3 * 7 * 4,
          "O(L*N) target-row work telemetry changed");
  const auto recombinant_minimum = profile.minimum_positions[0];
  const auto boundary_distance = [](std::size_t position, std::size_t boundary) {
    const std::size_t direct = position > boundary ? position - boundary : boundary - position;
    return std::min(direct, length - direct);
  };
  require(
      std::min(boundary_distance(recombinant_minimum, 81),
               boundary_distance(recombinant_minimum, 161)) <= 30,
      "the planted mosaic did not create a recombinant PHYLPRO minimum near a boundary");
  require(profile.minimum_correlations[0] < 0.5F,
          "the planted mosaic did not create a strong recombinant correlation trough");

  rdp::PhylproOptions stripped = circular;
  stripped.gap_mode = rdp::PhylproGapMode::strip_columns;
  compare_with_brute_force(
      parsed.alignment, targets, enabled, stripped, "strip-any-missing source windows");
  const auto stripped_profile = rdp::phylpro_profile(
      parsed.alignment, targets, enabled, stripped, error);
  require(error.empty(), error);
  require(stripped_profile.eligible_columns < profile.eligible_columns,
          "strip-any-missing did not remove gapped columns");

  rdp::PhylproOptions include_self = circular;
  include_self.include_self = true;
  compare_with_brute_force(
      parsed.alignment, targets, enabled, include_self, "include-self regression");
  const auto self_profile = rdp::phylpro_profile(
      parsed.alignment, targets, enabled, include_self, error);
  require(error.empty(), error);
  bool self_changed = false;
  for (std::size_t point = 0; point < profile.points.size(); ++point) {
    for (std::size_t role = 0; role < targets.size(); ++role) {
      self_changed |= std::abs(
          profile.points[point].correlations[role] -
          self_profile.points[point].correlations[role]) > 1e-6F;
    }
  }
  require(self_changed, "including the zero-distance self observation changed no coefficient");

  std::vector<std::uint8_t> one_disabled = enabled;
  one_disabled[5] = 1;
  compare_with_brute_force(
      parsed.alignment, targets, one_disabled, circular, "disabled-context exclusion");
  const auto disabled_profile = rdp::phylpro_profile(
      parsed.alignment, targets, one_disabled, circular, error);
  require(error.empty() && disabled_profile.context_sequences == 6,
          "disabled context was not excluded exactly once");

  rdp::PhylproOptions linear = circular;
  linear.circular = false;
  compare_with_brute_force(
      parsed.alignment, targets, enabled, linear, "linear complete half-windows");
  const auto linear_profile = rdp::phylpro_profile(
      parsed.alignment, targets, enabled, linear, error);
  require(error.empty() &&
              linear_profile.points.size() == length - 2 * linear_profile.half_window_sites + 1 &&
              linear_profile.points.front().alignment_position ==
                  linear_profile.half_window_sites + 1,
          "linear topology emitted an incomplete edge window");

  rdp::PhylproOptions capped = linear;
  capped.window_sites = 5000;
  compare_with_brute_force(
      parsed.alignment, targets, enabled, capped, "linear half-window cap");
  const auto capped_profile = rdp::phylpro_profile(
      parsed.alignment, targets, enabled, capped, error);
  require(error.empty() && capped_profile.window_capped &&
              capped_profile.half_window_sites == length / 2 &&
              capped_profile.points.size() == 1,
          "linear oversized window was not capped to one complete split");

  std::vector<std::uint8_t> disabled_role = enabled;
  disabled_role[0] = 1;
  const auto disabled_role_profile = rdp::phylpro_profile(
      parsed.alignment, targets, disabled_role, circular, error);
  (void)disabled_role_profile;
  require(error.find("Disabled sequences") != std::string::npos,
          "disabled role validation was lost");
  const auto duplicate_role_profile = rdp::phylpro_profile(
      parsed.alignment, std::array<std::uint32_t, 3>{{0, 0, 2}}, enabled, circular, error);
  (void)duplicate_role_profile;
  require(error.find("three distinct") != std::string::npos,
          "distinct-role validation was lost");

  auto flat_parsed = rdp::build_alignment(
      "flat-phylpro",
      {"flat-a", "flat-c", "flat-g"},
      {std::string(20, 'A'), std::string(20, 'C'), std::string(20, 'G')});
  require(flat_parsed.ok(), flat_parsed.error);
  rdp::PhylproOptions flat_options;
  flat_options.window_sites = 10;
  flat_options.circular = true;
  const auto flat_profile = rdp::phylpro_profile(
      flat_parsed.alignment,
      std::array<std::uint32_t, 3>{{0, 1, 2}},
      std::vector<std::uint8_t>(3, 0),
      flat_options,
      error);
  require(error.empty() && flat_profile.minimum_positions ==
              std::array<std::size_t, 3>{{1, 1, 1}} &&
              flat_profile.minimum_correlations ==
              std::array<float, 3>{{1.0F, 1.0F, 1.0F}},
          "zero-variance source fallback did not retain the first minimum position");

  std::cout
      << "PHYLPRO core verified: source-shaped polymorphic/gap mapping, circular and linear "
         "windows, rolling O(L*N) target rows against brute-force recomputation, Pearson "
         "self policy, disabled context, planted mosaic trough, and window caps passed.\n";
  return 0;
}
