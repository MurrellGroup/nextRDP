#include "threeseq.hpp"

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

constexpr std::array<std::uint8_t, 3> kParentOne{{1, 2, 0}};
constexpr std::array<std::uint8_t, 3> kParentTwo{{2, 0, 1}};
// Module5 invokes TSXOver as (0,1,2), (1,2,0), then (2,0,1). Seq3 is the
// selected recombinant in each call, giving this stable target order.
constexpr std::array<std::uint8_t, 3> kSourceTargetOrder{{2, 0, 1}};
// Pair-profile slots are 0:1, 0:2, and 1:2. These mappings mirror the
// supplied Seq1/Seq2 parent rotation rather than the parent member indices.
constexpr std::array<std::uint8_t, 3> kTargetParentOnePair{{0, 2, 1}};
constexpr std::array<std::uint8_t, 3> kTargetParentTwoPair{{1, 0, 2}};

std::uint8_t pair_index_for_members(std::uint8_t first, std::uint8_t second) {
  if (first > second) std::swap(first, second);
  if (first == 0 && second == 1) return 0;
  if (first == 0 && second == 2) return 1;
  return 2;
}

bool build_target_walk(
    const MaxChiWorkspace& variable_workspace,
    std::uint8_t target_local,
    ThreeSeqWorkspace& workspace) {
  workspace.coordinates.clear();
  workspace.steps.clear();
  workspace.heights.clear();
  if (target_local > 2 || variable_workspace.coordinates.empty()) return false;
  for (const auto& matches : variable_workspace.matches) {
    if (matches.size() != variable_workspace.coordinates.size()) return false;
  }

  const std::size_t parent_one_pair = kTargetParentOnePair[target_local];
  const std::size_t parent_two_pair = kTargetParentTwoPair[target_local];
  workspace.coordinates.reserve(variable_workspace.coordinates.size());
  workspace.steps.reserve(variable_workspace.coordinates.size());
  workspace.heights.reserve(variable_workspace.coordinates.size());
  std::int64_t height = 0;
  for (std::size_t index = 0;
       index < variable_workspace.coordinates.size();
       ++index) {
    const bool matches_parent_one =
        variable_workspace.matches[parent_one_pair][index] != 0;
    const bool matches_parent_two =
        variable_workspace.matches[parent_two_pair][index] != 0;
    // FindSubSeqTS retains a site only when the proposed parents differ and
    // the target matches exactly one of them. At a non-monomorphic triplet
    // site both equality flags cannot be true.
    if (matches_parent_one == matches_parent_two) continue;
    const std::int8_t step = matches_parent_one ? 1 : -1;
    height += step;
    workspace.coordinates.push_back(variable_workspace.coordinates[index]);
    workspace.steps.push_back(step);
    workspace.heights.push_back(height);
  }
  // FindSubSeqTS returns the last populated index (Y - 1); TSXOver exits when
  // that value is below three, so the automated path requires four sites.
  return workspace.coordinates.size() >= 4;
}

double approximate_normal_pdf(double value) {
  constexpr double kPi = 3.14159265359;
  return std::exp(-0.5 * value * value) / std::sqrt(2.0 * kPi);
}

double approximate_normal_cdf(double value) {
  constexpr double b1 = 0.31938153;
  constexpr double b2 = -0.356563782;
  constexpr double b3 = 1.781477937;
  constexpr double b4 = -1.821255978;
  constexpr double b5 = 1.330274429;
  constexpr double p = 0.2316419;
  constexpr double c = 0.39894228;
  if (value > 6.0) return 1.0;
  if (value < -6.0) return 0.0;
  const double t = 1.0 / (1.0 + p * std::fabs(value));
  const double tail = c * std::exp(-value * value / 2.0) * t *
      (t * (t * (t * (t * b5 + b4) + b3) + b2) + b1);
  return value >= 0.0 ? 1.0 - tail : tail;
}

double approximate_nu(double value) {
  if (!(value > 0.0)) return 0.0;
  const double half = value / 2.0;
  const double cdf = approximate_normal_cdf(half);
  const double denominator = value *
      (approximate_normal_pdf(half) + value * cdf / 2.0);
  return denominator > 0.0 ? ((cdf - 0.5) * 2.0) / denominator : 0.0;
}

double source_siegmund_discrete(
    std::size_t plus,
    std::size_t minus,
    std::size_t excursion) {
  const double total = static_cast<double>(plus + minus);
  if (!(total > 0.0)) return -1.0;
  const double boundary = static_cast<double>(excursion) - 0.5;
  const double endpoint = static_cast<double>(minus) - static_cast<double>(plus);
  double exponent = -2.0 * boundary * (boundary - endpoint) / total;
  if (exponent > 700.0) exponent = 700.0;
  double first = std::exp(exponent);
  if (first > 1.0e200) first = 1.0e200;
  const double second = first *
      (2.0 * (2.0 * boundary - endpoint) * (boundary - endpoint) / total + 1.0);
  const double argument = 2.0 * (2.0 * boundary - endpoint) / total;
  const double nu = approximate_nu(argument);
  double third = nu * nu * second;
  if (third < -1.0) third = -1.0;
  const double result = 1.0 - std::exp(-third);
  return std::isfinite(result) ? result : -1.0;
}

bool exact_transition_budget_ok(
    std::size_t plus,
    std::size_t minus,
    std::size_t excursion,
    std::uint64_t maximum) {
  if (excursion == 0) return true;
  const long double estimate =
      static_cast<long double>(plus + minus + 1) *
      static_cast<long double>(plus + 1) *
      static_cast<long double>(excursion);
  return estimate <= static_cast<long double>(maximum) &&
      excursion <= 4096 && plus <= 4096 && minus <= 4096;
}

float exact_excursion_probability(
    std::size_t plus,
    std::size_t minus,
    std::size_t excursion,
    ThreeSeqWorkspace& workspace) {
  // This evaluates the same finite hypergeometric random-walk distribution
  // represented by Seq3PVals/Get3SeqPvalC without allocating its four-
  // dimensional desktop table. State d is the current drop below the running
  // maximum. Paths reaching d == excursion have crossed the test boundary.
  if (excursion == 0 || plus == 0 || minus == 0) return 1.0F;
  const std::size_t stride = excursion;
  const std::size_t state_count = (plus + 1) * stride;
  workspace.probability_state.assign(state_count, 0.0L);
  workspace.probability_next.assign(state_count, 0.0L);
  workspace.probability_state[0] = 1.0F;
  const std::size_t total = plus + minus;

  for (std::size_t used = 0; used < total; ++used) {
    std::fill(
        workspace.probability_next.begin(),
        workspace.probability_next.end(),
        0.0F);
    const std::size_t first_plus = used > minus ? used - minus : 0;
    const std::size_t last_plus = std::min(plus, used);
    const float remaining = static_cast<float>(total - used);
    for (std::size_t used_plus = first_plus;
         used_plus <= last_plus;
         ++used_plus) {
      const std::size_t used_minus = used - used_plus;
      const std::size_t base = used_plus * stride;
      for (std::size_t drop = 0; drop < excursion; ++drop) {
        const float probability = workspace.probability_state[base + drop];
        if (!(probability > 0.0F)) continue;
        if (used_plus < plus) {
          const std::size_t next_drop = drop == 0 ? 0 : drop - 1;
          workspace.probability_next[(used_plus + 1) * stride + next_drop] +=
              probability * static_cast<float>(plus - used_plus) / remaining;
        }
        if (used_minus < minus && drop + 1 < excursion) {
          workspace.probability_next[base + drop + 1] +=
              probability * static_cast<float>(minus - used_minus) / remaining;
        }
      }
    }
    workspace.probability_state.swap(workspace.probability_next);
  }

  float survival = 0.0F;
  const std::size_t final_base = plus * stride;
  for (std::size_t drop = 0; drop < excursion; ++drop) {
    survival += workspace.probability_state[final_base + drop];
  }
  return std::clamp(1.0F - survival, 0.0F, 1.0F);
}

struct ProbabilityResult {
  double probability = 1.0;
  bool exact = false;
  bool approximation = false;
};

ProbabilityResult excursion_probability(
    std::size_t plus,
    std::size_t minus,
    std::size_t excursion,
    const ThreeSeqDiscoveryOptions& options,
    ThreeSeqWorkspace& workspace) {
  ProbabilityResult result;
  if (excursion == 0 || plus == 0 || minus == 0) return result;
  const ThreeSeqProbabilityKey key{plus, minus, excursion};
  if (const auto found = workspace.exact_probability_cache.find(key);
      found != workspace.exact_probability_cache.end()) {
    result.probability = found->second;
    result.exact = true;
    return result;
  }
  if (exact_transition_budget_ok(
          plus,
          minus,
          excursion,
          options.maximum_exact_state_transitions)) {
    const float probability =
        exact_excursion_probability(plus, minus, excursion, workspace);
    if (workspace.exact_probability_cache.size() >= 8192) {
      workspace.exact_probability_cache.clear();
    }
    workspace.exact_probability_cache.emplace(key, probability);
    result.probability = probability;
    result.exact = true;
    return result;
  }
  result.probability = source_siegmund_discrete(plus, minus, excursion);
  result.approximation = true;
  if (!(result.probability > 0.0) || result.probability >= 1.0) {
    // GetTSPVal rescales into the available exact table if the supplied
    // approximation leaves the open probability interval. A bounded exact
    // evaluation is preferable when scaling can make one available.
    const std::size_t maximum_dimension = std::max({plus, minus, excursion});
    if (maximum_dimension > 1) {
      const long double target = std::cbrt(
          static_cast<long double>(options.maximum_exact_state_transitions));
      const long double initial_scale = std::max<long double>(
          1.0L,
          static_cast<long double>(maximum_dimension) /
              std::max<long double>(2.0L, target));
      const std::size_t scaled_plus =
          static_cast<std::size_t>(plus / initial_scale);
      const std::size_t scaled_minus =
          static_cast<std::size_t>(minus / initial_scale);
      const std::size_t scaled_excursion =
          static_cast<std::size_t>(excursion / initial_scale);
      if (exact_transition_budget_ok(
              scaled_plus,
              scaled_minus,
              scaled_excursion,
              options.maximum_exact_state_transitions)) {
        const double scaled = exact_excursion_probability(
            scaled_plus, scaled_minus, scaled_excursion, workspace);
        // GetTSPVal first scales all three table coordinates by PVM, then
        // replaces PVM with onM/nM (or onN/nN when nM rounded to zero).
        // Retaining that second ratio matters at the exact-table boundary.
        const long double source_scale = scaled_plus > 0
            ? static_cast<long double>(plus) /
                static_cast<long double>(scaled_plus)
            : scaled_minus > 0
                ? static_cast<long double>(minus) /
                    static_cast<long double>(scaled_minus)
                : initial_scale;
        result.probability = scaled > 0.0
            ? std::pow(scaled, static_cast<double>(source_scale))
            : 0.0;
        // The supplied path floors a positive pre-exponent probability that
        // underflows during the PVM exponentiation to 1e-300.
        if (result.probability == 0.0 && scaled > 0.0) {
          result.probability = 1.0e-300;
        }
      }
    }
  }
  result.probability = std::clamp(result.probability, 0.0, 1.0);
  return result;
}

struct ExcursionResult {
  std::size_t beginning = 0;
  std::size_t ending = 0;
  std::size_t beginning_index = 0;
  std::size_t ending_index = 0;
  std::size_t probability_magnitude = 0;
  std::size_t magnitude = 0;
};

struct SubProbabilityResult {
  ProbabilityResult probability;
  std::size_t excursion = 0;
};

struct SplitResult {
  std::size_t beginning = 0;
  std::size_t ending = 0;
  std::size_t probability_excursion = 0;
  ProbabilityResult probability;
  std::size_t exact_evaluations = 0;
  std::size_t approximate_evaluations = 0;
  bool missing_data_found = false;
};

ExcursionResult source_excursion(
    const std::vector<std::size_t>& coordinates,
    const std::vector<std::int8_t>& steps,
    std::int8_t multiplier,
    bool circular,
    std::size_t alignment_length,
    std::vector<std::int64_t>& heights) {
  ExcursionResult result;
  if (coordinates.empty() || coordinates.size() != steps.size()) return result;
  heights.resize(steps.size());
  std::int64_t current = 0;
  std::int64_t maximum = 0;
  for (std::size_t index = 0; index < steps.size(); ++index) {
    current += static_cast<std::int64_t>(steps[index]) * multiplier;
    heights[index] = current;
    if (current > maximum) {
      maximum = current;
      result.beginning = coordinates[index];
      result.beginning_index = index;
    }
    const std::int64_t excursion = maximum - current;
    if (excursion > static_cast<std::int64_t>(result.magnitude)) {
      result.magnitude = static_cast<std::size_t>(excursion);
      result.ending = coordinates[index];
      result.ending_index = index;
    }
  }
  // TSXOver reads the exact/fallback probability before calling CheckwrapC.
  // The latter may increase nK while extending the selected tract through the
  // origin, but the source does not recalculate the probability afterwards.
  result.probability_magnitude = result.magnitude;
  if (result.beginning == 0) {
    result.beginning = coordinates.front();
    result.beginning_index = 0;
  }

  // Literal CheckwrapC prefix extension. It runs before the linear/circular
  // branch, then advances BE to the first informative site within the tract.
  maximum = heights[result.beginning_index];
  const std::int64_t terminal_height = heights.back();
  const std::size_t loop_end = result.beginning < result.ending
      ? result.beginning_index
      : result.ending_index;
  for (std::size_t index = 0; index <= loop_end && index < heights.size(); ++index) {
    const std::int64_t extended = terminal_height + heights[index];
    if (extended > maximum) {
      maximum = extended;
      result.beginning = coordinates[index];
      result.beginning_index = index;
    }
    const std::int64_t excursion = maximum - extended;
    if (excursion > static_cast<std::int64_t>(result.magnitude)) {
      result.magnitude = static_cast<std::size_t>(excursion);
      result.ending = coordinates[index];
      result.ending_index = index;
    }
  }
  result.beginning_index = (result.beginning_index + 1) % coordinates.size();
  result.beginning = coordinates[result.beginning_index];

  if (!circular && result.beginning > result.ending) {
    const std::size_t original_beginning_index = result.beginning_index;
    const std::size_t shifted_beginning = result.ending_index + 1 < coordinates.size()
        ? coordinates[result.ending_index + 1]
        : 1;
    result.ending = original_beginning_index > 0
        ? coordinates[original_beginning_index - 1]
        : alignment_length;
    result.beginning = shifted_beginning;
    result.beginning_index = result.ending_index + 1 < coordinates.size()
        ? result.ending_index + 1
        : 0;
    result.ending_index = original_beginning_index > 0
        ? original_beginning_index - 1
        : coordinates.size() - 1;
  }
  return result;
}

SubProbabilityResult source_sub_probability(
    const std::vector<std::size_t>& coordinates,
    const std::vector<std::int8_t>& steps,
    std::int8_t multiplier,
    std::size_t beginning,
    std::size_t ending,
    std::size_t plus,
    std::size_t minus,
    const ThreeSeqDiscoveryOptions& options,
    ThreeSeqWorkspace& workspace) {
  SubProbabilityResult result;
  if (coordinates.empty() || coordinates.size() != steps.size()) return result;

  workspace.heights.resize(steps.size());
  std::int64_t height = 0;
  for (std::size_t index = 0; index < steps.size(); ++index) {
    height += static_cast<std::int64_t>(steps[index]) * multiplier;
    workspace.heights[index] = height;
  }

  // FindSubSeqTS2 leaves XPosDiff[x] as the number of retained sites at or
  // before x. This is exactly upper_bound on the compact coordinate vector.
  const auto mapped_position = [&](std::size_t coordinate) {
    return static_cast<std::ptrdiff_t>(std::upper_bound(
        coordinates.begin(), coordinates.end(), coordinate) - coordinates.begin());
  };
  const std::ptrdiff_t last =
      static_cast<std::ptrdiff_t>(coordinates.size() - 1);
  const std::ptrdiff_t mapped_beginning = mapped_position(beginning);
  std::ptrdiff_t cursor = mapped_beginning > 0 ? mapped_beginning - 1 : last;
  const std::ptrdiff_t mapped_ending = mapped_position(ending);
  std::ptrdiff_t endpoint = mapped_ending + 1 > last
      ? last
      : mapped_ending + 1;
  if (cursor == endpoint) --endpoint;
  if (endpoint == 0) endpoint = last;
  if (endpoint < 0) endpoint = 0;

  std::int64_t maximum = std::numeric_limits<std::int64_t>::lowest();
  std::int64_t minimum = std::numeric_limits<std::int64_t>::max();
  std::int64_t wrap_offset = 0;
  std::size_t guard = coordinates.size() * 2 + 2;
  while (cursor != endpoint && guard-- > 0) {
    if (cursor < 0) cursor = last + cursor;
    if (cursor > last) {
      wrap_offset = workspace.heights.back();
      cursor = cursor - last - 1;
    }
    if (cursor < 0 || cursor > last) return result;
    const std::int64_t value =
        workspace.heights[static_cast<std::size_t>(cursor)] + wrap_offset;
    maximum = std::max(maximum, value);
    minimum = std::min(minimum, value);
    // SubPVal tests this again after an in-loop wrap, so the wrapped endpoint
    // itself is included even though the ordinary loop condition excludes it.
    if (cursor == endpoint) break;
    ++cursor;
  }
  if (maximum == std::numeric_limits<std::int64_t>::lowest() ||
      minimum == std::numeric_limits<std::int64_t>::max() || maximum < minimum) {
    return result;
  }
  result.excursion = static_cast<std::size_t>(maximum - minimum);
  result.probability = excursion_probability(
      plus, minus, result.excursion, options, workspace);
  return result;
}

SplitResult source_check_split(
    const ExcursionResult& excursion,
    const ProbabilityResult& probability,
    std::size_t plus,
    std::size_t minus,
    std::int8_t multiplier,
    const std::vector<std::size_t>& coordinates,
    const std::vector<std::int8_t>& steps,
    const std::vector<std::uint8_t>& triplet_missing_data,
    std::size_t alignment_length,
    const ThreeSeqDiscoveryOptions& options,
    ThreeSeqWorkspace& workspace) {
  SplitResult result;
  result.beginning = excursion.beginning;
  result.ending = excursion.ending;
  result.probability_excursion = excursion.probability_magnitude;
  result.probability = probability;
  if (alignment_length == 0 ||
      triplet_missing_data.size() != alignment_length) {
    return result;
  }

  const auto is_missing = [&](std::ptrdiff_t coordinate) {
    return coordinate > 0 &&
        coordinate <= static_cast<std::ptrdiff_t>(alignment_length) &&
        triplet_missing_data[static_cast<std::size_t>(coordinate - 1)] != 0;
  };
  const auto normalize_boundary = [&](std::ptrdiff_t coordinate) {
    if (coordinate <= 0) return alignment_length;
    if (coordinate > static_cast<std::ptrdiff_t>(alignment_length)) {
      return std::size_t{1};
    }
    return static_cast<std::size_t>(coordinate);
  };
  const auto account = [&](const SubProbabilityResult& sub) {
    result.exact_evaluations += static_cast<std::size_t>(sub.probability.exact);
    result.approximate_evaluations +=
        static_cast<std::size_t>(sub.probability.approximation);
  };

  const std::size_t original_beginning = excursion.beginning;
  const std::size_t original_ending = excursion.ending;
  std::ptrdiff_t coordinate =
      static_cast<std::ptrdiff_t>(original_beginning) - 1;
  bool forward_missing_found = false;
  std::size_t guard = alignment_length + 2;
  while (coordinate != static_cast<std::ptrdiff_t>(original_ending) && guard-- > 0) {
    if (coordinate > static_cast<std::ptrdiff_t>(alignment_length)) {
      coordinate = 1;
      // Literal CheckSplit3Seq early exit after wrapping to an ending at the
      // first source coordinate.
      if (original_ending <= 1) break;
    }
    if (coordinate < 0) coordinate = static_cast<std::ptrdiff_t>(alignment_length);
    if (is_missing(coordinate)) {
      result.missing_data_found = true;
      forward_missing_found = true;
      // CheckSplit3Seq discards the incoming probability as soon as the first
      // missing state is found, then admits only a sub-tract value below one.
      result.probability = {};
      const SubProbabilityResult sub = source_sub_probability(
          coordinates,
          steps,
          multiplier,
          original_beginning,
          static_cast<std::size_t>(coordinate),
          plus,
          minus,
          options,
          workspace);
      account(sub);
      if (sub.probability.probability < result.probability.probability) {
        result.probability = sub.probability;
        result.probability_excursion = sub.excursion;
        result.beginning = original_beginning;
        result.ending = normalize_boundary(coordinate - 1);
      }
      break;
    }
    if (coordinate == static_cast<std::ptrdiff_t>(original_ending)) break;
    ++coordinate;
    if (coordinate == static_cast<std::ptrdiff_t>(original_ending)) break;
  }

  if (forward_missing_found) {
    coordinate = static_cast<std::ptrdiff_t>(original_ending) + 1;
    guard = alignment_length + 2;
    while (coordinate != static_cast<std::ptrdiff_t>(original_beginning) && guard-- > 0) {
      if (coordinate < 0) coordinate = static_cast<std::ptrdiff_t>(alignment_length);
      if (coordinate > static_cast<std::ptrdiff_t>(alignment_length)) coordinate = 1;
      if (is_missing(coordinate)) {
        result.missing_data_found = true;
        const SubProbabilityResult sub = source_sub_probability(
            coordinates,
            steps,
            multiplier,
            static_cast<std::size_t>(coordinate),
            original_ending,
            plus,
            minus,
            options,
            workspace);
        account(sub);
        if (sub.probability.probability < result.probability.probability) {
          result.probability = sub.probability;
          result.probability_excursion = sub.excursion;
          result.beginning = normalize_boundary(coordinate + 1);
          result.ending = original_ending;
        }
        break;
      }
      if (coordinate == static_cast<std::ptrdiff_t>(original_beginning)) break;
      --coordinate;
    }
  }
  return result;
}

double source_corrected_probability(
    double raw,
    bool correction_enabled,
    std::uint64_t tests) {
  if (!correction_enabled) return std::clamp(raw, 0.0, 1.0);
  tests = std::max<std::uint64_t>(1, tests);
  const double product = raw * static_cast<double>(tests);
  // TSXOver uses the Dunn–Šidák expression only above 1e-15 and uses the
  // unbounded product for smaller tails. expm1/log1p retain the same branch
  // while avoiding cancellation in the browser's double-precision runtime.
  double corrected = 1.0;
  if (raw < 1.0 && product < 1.0) {
    corrected = raw > 1.0e-15
        ? -std::expm1(static_cast<double>(tests) * std::log1p(-raw))
        : product;
  }
  // The later TSXOver guard also replaces a literal zero with xpValue.
  if (corrected == 0.0 && product > 0.0) corrected = product;
  return std::clamp(corrected, 0.0, 1.0);
}

double source_findall_corrected_probability(
    double raw,
    bool correction_enabled,
    std::uint64_t tests) {
  if (!correction_enabled) return std::clamp(raw, 0.0, 1.0);
  tests = std::max<std::uint64_t>(1, tests);
  const double product = raw * static_cast<double>(tests);
  // After TSXOver(1) swaps to the second orientation it uses Dunn–Šidák for
  // every open-interval P, without the earlier 1e-15 shortcut.
  double corrected = 1.0;
  if (raw > 0.0 && raw < 1.0) {
    corrected = -std::expm1(
        static_cast<double>(tests) * std::log1p(-raw));
  }
  if (corrected == 0.0 && product > 0.0) corrected = product;
  return std::clamp(corrected, 0.0, 1.0);
}

bool source_threshold_passes(
    double raw,
    double corrected,
    bool correction_enabled,
    std::uint64_t tests,
    double cutoff) {
  const double xp_value = correction_enabled
      ? raw * static_cast<double>(std::max<std::uint64_t>(1, tests))
      : raw;
  // Literal TSXOver call gate. The dPValue == 1 branch exists for a rounded
  // Dunn–Šidák result, but still compares the unbounded product strictly.
  return xp_value > 0.0 &&
      ((corrected < 1.0 && corrected <= cutoff) ||
       (corrected == 1.0 && xp_value < cutoff));
}

}  // namespace

std::size_t ThreeSeqProbabilityKeyHash::operator()(
    const ThreeSeqProbabilityKey& key) const noexcept {
  std::size_t seed = key.plus + 0x9e3779b97f4a7c15ULL;
  seed ^= key.minus + 0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
  seed ^= key.excursion + 0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
  return seed;
}

ThreeSeqDiscoverySummary threeseq_discover_prepared(
    const MaxChiWorkspace& variable_workspace,
    const std::vector<std::uint8_t>& triplet_missing_data,
    const std::array<double, 3>& pair_similarity,
    const ThreeSeqDiscoveryOptions& options,
    ThreeSeqWorkspace& workspace,
    std::vector<ThreeSeqDiscoveryCandidate>& output) {
  output.clear();
  ThreeSeqDiscoverySummary summary;
  summary.correction_applied = options.correction_enabled;
  summary.correction_tests = std::max<std::uint64_t>(1, options.correction_tests);
  output.reserve(3);
  const std::size_t alignment_length = variable_workspace.variable_prefix.empty()
      ? 0
      : variable_workspace.variable_prefix.size() - 1;

  for (const std::uint8_t target_local : kSourceTargetOrder) {
    if (!build_target_walk(variable_workspace, target_local, workspace)) continue;
    ++summary.target_profiles_scanned;
    const std::size_t plus = static_cast<std::size_t>(std::count(
        workspace.steps.begin(), workspace.steps.end(), std::int8_t{1}));
    const std::size_t minus = workspace.steps.size() - plus;

    ExcursionResult descent = source_excursion(
        workspace.coordinates,
        workspace.steps,
        1,
        options.circular,
        alignment_length,
        workspace.heights);
    const ProbabilityResult descent_probability = excursion_probability(
        plus, minus, descent.probability_magnitude, options, workspace);
    ExcursionResult ascent = source_excursion(
        workspace.coordinates,
        workspace.steps,
        -1,
        options.circular,
        alignment_length,
        workspace.heights);
    const ProbabilityResult ascent_probability = excursion_probability(
        minus, plus, ascent.probability_magnitude, options, workspace);
    summary.exact_probability_evaluations +=
        static_cast<std::size_t>(descent_probability.exact) +
        static_cast<std::size_t>(ascent_probability.exact);
    summary.approximate_probability_evaluations +=
        static_cast<std::size_t>(descent_probability.approximation) +
        static_cast<std::size_t>(ascent_probability.approximation);

    bool use_ascent =
        ascent_probability.probability < descent_probability.probability;
    ExcursionResult selected = use_ascent ? ascent : descent;
    ProbabilityResult probability =
        use_ascent ? ascent_probability : descent_probability;
    std::size_t selected_plus = use_ascent ? minus : plus;
    std::size_t selected_minus = use_ascent ? plus : minus;
    std::size_t probability_excursion = selected.probability_magnitude;
    bool missing_data_split_applied = false;
    // Literal TSXOver low-information exits, applied after parent swapping.
    if ((selected_minus > 0 && selected.magnitude == 1) ||
        (selected_minus >= selected_plus &&
         selected_minus - selected_plus == selected.magnitude)) {
      continue;
    }
    // TSXOver maps a zero table value to 10 before comparing orientations and
    // requires xpvalue > 0, so a literal underflow is not a discovery call.
    if (!(probability.probability > 0.0) || probability.probability > 1.0) continue;
    double corrected = source_corrected_probability(
        probability.probability,
        options.correction_enabled,
        summary.correction_tests);
    if (!source_threshold_passes(
            probability.probability,
            corrected,
            options.correction_enabled,
            summary.correction_tests,
            options.p_value_cutoff)) {
      continue;
    }

    // TSXOver enters CheckSplit3Seq only after the original candidate has
    // passed its first corrected threshold and at least one previous event is
    // present. The selected orientation is tried first. If trimming makes it
    // worse than the unsplit reverse orientation, the reverse walk is split
    // too and wins only on a strictly lower probability.
    if (options.post_erasure_split_enabled) {
      SplitResult split = source_check_split(
          selected,
          probability,
          selected_plus,
          selected_minus,
          use_ascent ? -1 : 1,
          workspace.coordinates,
          workspace.steps,
          triplet_missing_data,
          alignment_length,
          options,
          workspace);
      summary.exact_probability_evaluations += split.exact_evaluations;
      summary.approximate_probability_evaluations +=
          split.approximate_evaluations;

      const ExcursionResult alternative = use_ascent ? descent : ascent;
      const ProbabilityResult alternative_probability =
          use_ascent ? descent_probability : ascent_probability;
      const std::size_t alternative_plus = use_ascent ? plus : minus;
      const std::size_t alternative_minus = use_ascent ? minus : plus;
      if (split.probability.probability > alternative_probability.probability) {
        SplitResult alternative_split = source_check_split(
            alternative,
            alternative_probability,
            alternative_plus,
            alternative_minus,
            use_ascent ? 1 : -1,
            workspace.coordinates,
            workspace.steps,
            triplet_missing_data,
            alignment_length,
            options,
            workspace);
        summary.exact_probability_evaluations +=
            alternative_split.exact_evaluations;
        summary.approximate_probability_evaluations +=
            alternative_split.approximate_evaluations;
        if (split.probability.probability >
            alternative_split.probability.probability) {
          use_ascent = !use_ascent;
          selected = alternative;
          selected_plus = alternative_plus;
          selected_minus = alternative_minus;
          split = std::move(alternative_split);
        }
      }
      selected.beginning = split.beginning;
      selected.ending = split.ending;
      probability = split.probability;
      probability_excursion = split.probability_excursion;
      missing_data_split_applied = split.missing_data_found;
      corrected = source_corrected_probability(
          probability.probability,
          options.correction_enabled,
          summary.correction_tests);
      if (!source_threshold_passes(
              probability.probability,
              corrected,
              options.correction_enabled,
              summary.correction_tests,
              options.p_value_cutoff)) {
        continue;
      }
    }

    ThreeSeqDiscoveryCandidate candidate;
    candidate.beginning = selected.beginning;
    candidate.ending = selected.ending;
    candidate.wraps_origin = options.circular &&
        candidate.beginning >= candidate.ending;
    candidate.informative_beginning = static_cast<std::size_t>(std::upper_bound(
        workspace.coordinates.begin(),
        workspace.coordinates.end(),
        candidate.beginning) - workspace.coordinates.begin());
    candidate.informative_ending = static_cast<std::size_t>(std::upper_bound(
        workspace.coordinates.begin(),
        workspace.coordinates.end(),
        candidate.ending) - workspace.coordinates.begin());
    candidate.target_local = target_local;
    candidate.recombinant_local = target_local;
    // The selected descent is the tract matching parent two. When the ascent
    // wins, SwapRound exchanges parents/counts before the same assignment.
    candidate.major_parent_local = use_ascent
        ? kParentTwo[target_local]
        : kParentOne[target_local];
    candidate.minor_parent_local = use_ascent
        ? kParentOne[target_local]
        : kParentTwo[target_local];
    candidate.candidate_pair = pair_index_for_members(
        candidate.recombinant_local, candidate.minor_parent_local);
    candidate.direction = use_ascent
        ? ThreeSeqWalkDirection::ascent
        : ThreeSeqWalkDirection::descent;
    candidate.information_rich_sites = workspace.steps.size();
    candidate.parent_one_matches = selected_plus;
    candidate.parent_two_matches = selected_minus;
    candidate.probability_excursion = probability_excursion;
    candidate.maximum_excursion = selected.magnitude;
    candidate.raw_p_value = probability.probability;
    candidate.corrected_p_value = corrected;
    candidate.pair_similarity = pair_similarity;
    candidate.exact_probability = probability.exact;
    candidate.siegmund_fallback = probability.approximation;
    candidate.missing_data_split_applied = missing_data_split_applied;
    output.push_back(candidate);
  }
  summary.emitted_candidates = output.size();
  return summary;
}

ThreeSeqRecheckEvidence threeseq_recheck_prepared(
    const MaxChiWorkspace& variable_workspace,
    const std::vector<std::uint8_t>& triplet_missing_data,
    const ThreeSeqDiscoveryOptions& options,
    ThreeSeqWorkspace& workspace) {
  ThreeSeqRecheckEvidence evidence;
  evidence.requested = true;
  evidence.correction_applied = options.correction_enabled;
  evidence.correction_tests =
      std::max<std::uint64_t>(1, options.correction_tests);
  const std::size_t alignment_length = variable_workspace.variable_prefix.empty()
      ? 0
      : variable_workspace.variable_prefix.size() - 1;

  const auto account_probability = [&](const ProbabilityResult& probability) {
    evidence.exact_probability_evaluations +=
        static_cast<std::size_t>(probability.exact);
    evidence.approximate_probability_evaluations +=
        static_cast<std::size_t>(probability.approximation);
  };
  const auto account_split = [&](const SplitResult& split) {
    evidence.exact_probability_evaluations += split.exact_evaluations;
    evidence.approximate_probability_evaluations +=
        split.approximate_evaluations;
  };

  for (const std::uint8_t target_local : kSourceTargetOrder) {
    if (!build_target_walk(variable_workspace, target_local, workspace)) continue;
    evidence.profile_available = true;
    ++evidence.target_profiles_scanned;
    const std::size_t plus = static_cast<std::size_t>(std::count(
        workspace.steps.begin(), workspace.steps.end(), std::int8_t{1}));
    const std::size_t minus = workspace.steps.size() - plus;
    const ExcursionResult descent = source_excursion(
        workspace.coordinates,
        workspace.steps,
        1,
        options.circular,
        alignment_length,
        workspace.heights);
    const ProbabilityResult descent_probability = excursion_probability(
        plus, minus, descent.probability_magnitude, options, workspace);
    const ExcursionResult ascent = source_excursion(
        workspace.coordinates,
        workspace.steps,
        -1,
        options.circular,
        alignment_length,
        workspace.heights);
    const ProbabilityResult ascent_probability = excursion_probability(
        minus, plus, ascent.probability_magnitude, options, workspace);
    account_probability(descent_probability);
    account_probability(ascent_probability);

    const bool first_is_ascent =
        ascent_probability.probability < descent_probability.probability;
    const ExcursionResult& first_excursion =
        first_is_ascent ? ascent : descent;
    const ProbabilityResult& first_probability =
        first_is_ascent ? ascent_probability : descent_probability;
    const std::size_t first_plus = first_is_ascent ? minus : plus;
    const std::size_t first_minus = first_is_ascent ? plus : minus;
    // TSXOver exits before entering its Findall loop when the initially
    // selected orientation fails either low-information or corrected gate.
    if ((first_minus > 0 && first_excursion.magnitude == 1) ||
        (first_minus >= first_plus &&
         first_minus - first_plus == first_excursion.magnitude) ||
        !(first_probability.probability > 0.0) ||
        first_probability.probability > 1.0) {
      continue;
    }
    const double initial_corrected = source_corrected_probability(
        first_probability.probability,
        options.correction_enabled,
        evidence.correction_tests);
    if (!source_threshold_passes(
            first_probability.probability,
            initial_corrected,
            options.correction_enabled,
            evidence.correction_tests,
            options.p_value_cutoff)) {
      continue;
    }

    SplitResult first_split;
    first_split.beginning = first_excursion.beginning;
    first_split.ending = first_excursion.ending;
    first_split.probability_excursion =
        first_excursion.probability_magnitude;
    first_split.probability = first_probability;
    const bool second_is_ascent = !first_is_ascent;
    const ExcursionResult& second_excursion =
        second_is_ascent ? ascent : descent;
    const ProbabilityResult& second_probability =
        second_is_ascent ? ascent_probability : descent_probability;
    const std::size_t second_plus = second_is_ascent ? minus : plus;
    const std::size_t second_minus = second_is_ascent ? plus : minus;
    SplitResult second_split;
    second_split.beginning = second_excursion.beginning;
    second_split.ending = second_excursion.ending;
    second_split.probability_excursion =
        second_excursion.probability_magnitude;
    second_split.probability = second_probability;
    if (options.post_erasure_split_enabled) {
      first_split = source_check_split(
          first_excursion,
          first_probability,
          first_plus,
          first_minus,
          first_is_ascent ? -1 : 1,
          workspace.coordinates,
          workspace.steps,
          triplet_missing_data,
          alignment_length,
          options,
          workspace);
      // FindallFlag=1 forces the second CheckSplit3Seq call even when the
      // first trimmed probability remains lower.
      second_split = source_check_split(
          second_excursion,
          second_probability,
          second_plus,
          second_minus,
          second_is_ascent ? -1 : 1,
          workspace.coordinates,
          workspace.steps,
          triplet_missing_data,
          alignment_length,
          options,
          workspace);
      account_split(first_split);
      account_split(second_split);
    }

    const auto consider = [&](
        bool use_ascent,
        const ExcursionResult& excursion,
        const SplitResult& split,
        std::size_t selected_plus,
        std::size_t selected_minus,
        bool second_findall_orientation) {
      const double raw = split.probability.probability;
      const double corrected = second_findall_orientation
          ? source_findall_corrected_probability(
                raw, options.correction_enabled, evidence.correction_tests)
          : source_corrected_probability(
                raw, options.correction_enabled, evidence.correction_tests);
      if (!source_threshold_passes(
              raw,
              corrected,
              options.correction_enabled,
              evidence.correction_tests,
              options.p_value_cutoff)) {
        return;
      }
      ++evidence.qualifying_orientations;
      // Findall copies each call with parents and interval reversed.
      evidence.source_list_entries += 2;
      if (evidence.best_target >= 0 &&
          (corrected > evidence.corrected_p_value ||
           (corrected == evidence.corrected_p_value &&
            raw >= evidence.raw_p_value))) {
        return;
      }
      evidence.best_target = static_cast<std::int8_t>(target_local);
      evidence.best_direction = use_ascent
          ? ThreeSeqWalkDirection::ascent
          : ThreeSeqWalkDirection::descent;
      evidence.information_rich_sites = workspace.steps.size();
      evidence.parent_one_matches = selected_plus;
      evidence.parent_two_matches = selected_minus;
      evidence.probability_excursion = split.probability_excursion;
      evidence.maximum_excursion = excursion.magnitude;
      evidence.beginning = split.beginning;
      evidence.ending = split.ending;
      evidence.wraps_origin = options.circular &&
          split.beginning >= split.ending;
      evidence.raw_p_value = raw;
      evidence.corrected_p_value = corrected;
      evidence.exact_probability = split.probability.exact;
      evidence.siegmund_fallback = split.probability.approximation;
      evidence.missing_data_split_applied = split.missing_data_found;
    };
    consider(
        first_is_ascent,
        first_excursion,
        first_split,
        first_plus,
        first_minus,
        false);
    consider(
        second_is_ascent,
        second_excursion,
        second_split,
        second_plus,
        second_minus,
        true);
  }
  evidence.source_recheck_hit = evidence.best_target >= 0;
  return evidence;
}

ThreeSeqPlotProfile threeseq_plot_profile(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    ThreeSeqWorkspace& workspace,
    MaxChiWorkspace& variable_workspace) {
  ThreeSeqPlotProfile plot;
  if (alignment.length == 0) return plot;
  variable_workspace.coordinates.clear();
  for (auto& matches : variable_workspace.matches) matches.clear();
  variable_workspace.variable_prefix.assign(alignment.length + 1, 0);
  for (std::size_t coordinate = 1; coordinate <= alignment.length; ++coordinate) {
    const std::uint8_t first = alignment.at(triplet[0], coordinate - 1);
    const std::uint8_t second = alignment.at(triplet[1], coordinate - 1);
    const std::uint8_t third = alignment.at(triplet[2], coordinate - 1);
    if (first != 0 && second != 0 && third != 0 &&
        (first != second || first != third)) {
      variable_workspace.coordinates.push_back(coordinate);
      variable_workspace.matches[0].push_back(first == second ? 1 : 0);
      variable_workspace.matches[1].push_back(first == third ? 1 : 0);
      variable_workspace.matches[2].push_back(second == third ? 1 : 0);
    }
    variable_workspace.variable_prefix[coordinate] =
        variable_workspace.coordinates.size();
  }
  if (variable_workspace.coordinates.empty()) return plot;

  plot.coordinates = variable_workspace.coordinates;
  for (std::uint8_t target = 0; target < 3; ++target) {
    if (!build_target_walk(variable_workspace, target, workspace)) {
      plot.target_walks[target].assign(plot.coordinates.size(), 0.0);
      continue;
    }
    auto& trace = plot.target_walks[target];
    trace.resize(plot.coordinates.size());
    std::size_t source = 0;
    double height = 0.0;
    for (std::size_t output_index = 0;
         output_index < plot.coordinates.size();
         ++output_index) {
      while (source < workspace.coordinates.size() &&
             workspace.coordinates[source] <= plot.coordinates[output_index]) {
        height += workspace.steps[source];
        ++source;
      }
      trace[output_index] = height;
    }
  }
  plot.available = true;
  return plot;
}

}  // namespace rdp
