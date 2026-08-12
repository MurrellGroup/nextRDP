#include "geneconv.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <queue>
#include <utility>
#include <vector>

namespace rdp {
namespace {

constexpr double kMinimumProbability = 1.0e-300;
constexpr std::array<std::uint8_t, 6> kRecombinantLocal{{0, 0, 1, 0, 1, 2}};
constexpr std::array<std::uint8_t, 6> kMinorParentLocal{{1, 2, 2, 1, 0, 1}};
constexpr std::array<std::uint8_t, 6> kMajorParentLocal{{2, 1, 0, 2, 2, 0}};
constexpr std::array<std::uint8_t, 6> kCandidatePair{{0, 1, 2, 0, 0, 2}};
constexpr std::array<std::uint8_t, 6> kPlotPair{{0, 1, 2, 2, 1, 0}};

bool track_positive(std::uint8_t track, std::uint8_t category) {
  if (track < 3) return category == track;
  if (track == 3) return category == 2 || category == 3;
  if (track == 4) return category == 1 || category == 3;
  return category == 0 || category == 3;
}

void build_categories(
    const MaxChiWorkspace& variable_workspace,
    GeneconvWorkspace& workspace,
    std::array<std::size_t, 4>& category_counts) {
  const std::size_t length = variable_workspace.coordinates.size();
  workspace.categories.resize(length);
  category_counts.fill(0);
  for (std::size_t position = 0; position < length; ++position) {
    std::uint8_t category = 3;
    if (position < variable_workspace.matches[0].size() &&
        variable_workspace.matches[0][position] != 0) {
      category = 0;
    } else if (position < variable_workspace.matches[1].size() &&
               variable_workspace.matches[1][position] != 0) {
      category = 1;
    } else if (position < variable_workspace.matches[2].size() &&
               variable_workspace.matches[2][position] != 0) {
      category = 2;
    }
    workspace.categories[position] = category;
    ++category_counts[category];
  }
}

void append_source_run(
    std::vector<GeneconvRun>& runs,
    const GeneconvCategoryRun& category_run,
    bool positive,
    bool merge_adjacent_positive) {
  // GetFragsP coalesces adjacent positive outer fragments (not negative
  // fragments), while each inner category fragment retains its own slot.
  if (merge_adjacent_positive && positive && !runs.empty() &&
      runs.back().positive) {
    runs.back().ending = category_run.ending;
    runs.back().length += category_run.length;
    runs.back().wraps_origin = runs.back().wraps_origin ||
        category_run.wraps_origin || runs.back().beginning > runs.back().ending;
    return;
  }
  runs.push_back({
      category_run.beginning,
      category_run.ending,
      category_run.length,
      positive,
      category_run.wraps_origin,
  });
}

bool build_source_runs(bool circular, GeneconvWorkspace& workspace) {
  workspace.category_runs.clear();
  for (auto& runs : workspace.runs) runs.clear();
  if (workspace.categories.empty()) return false;
  if (workspace.category_runs.capacity() < workspace.categories.size()) {
    workspace.category_runs.reserve(workspace.categories.size());
  }
  std::size_t beginning = 0;
  std::uint8_t category = workspace.categories[0];
  for (std::size_t position = 1; position <= workspace.categories.size(); ++position) {
    const bool at_end = position == workspace.categories.size();
    if (!at_end && workspace.categories[position] == category) continue;
    workspace.category_runs.push_back({
        category,
        beginning,
        position - 1,
        position - beginning,
        false,
    });
    if (!at_end) {
      beginning = position;
      category = workspace.categories[position];
    }
  }

  if (circular && workspace.category_runs.size() == 1) {
    // GetFragsP returns zero when its circular scan reaches the starting
    // position without encountering another category.
    return false;
  }
  if (circular && workspace.category_runs.size() > 1 &&
      workspace.category_runs.front().category ==
          workspace.category_runs.back().category) {
    // Literal source ordering: keep the first run, and extend the terminal
    // run through that same first run. GetMaxFragScoreP then walks this list
    // once; it does not duplicate every signed run for circular data.
    auto& terminal = workspace.category_runs.back();
    terminal.ending = workspace.category_runs.front().ending;
    terminal.length += workspace.category_runs.front().length;
    terminal.wraps_origin = true;
  }

  for (auto& runs : workspace.runs) {
    if (runs.capacity() < workspace.category_runs.size()) {
      runs.reserve(workspace.category_runs.size());
    }
  }
  for (const auto& category_run : workspace.category_runs) {
    for (std::uint8_t track = 0; track < 6; ++track) {
      append_source_run(
          workspace.runs[track],
          category_run,
          track_positive(track, category_run.category),
          track >= 3);
    }
  }
  return true;
}

struct RangeMaximum {
  std::int64_t value = std::numeric_limits<std::int64_t>::lowest();
  std::size_t index = 0;
  bool available = false;
};

RangeMaximum better_maximum(RangeMaximum left, RangeMaximum right) {
  if (!left.available) return right;
  if (!right.available) return left;
  if (right.value > left.value ||
      (right.value == left.value && right.index > left.index)) {
    return right;
  }
  return left;
}

void build_range_maximum_tree(
    const std::vector<std::int64_t>& values,
    GeneconvWorkspace& workspace,
    std::size_t& base) {
  base = 1;
  while (base < values.size()) base *= 2;
  workspace.range_max_values.assign(
      base * 2, std::numeric_limits<std::int64_t>::lowest());
  workspace.range_max_indices.assign(base * 2, 0);
  for (std::size_t index = 0; index < values.size(); ++index) {
    workspace.range_max_values[base + index] = values[index];
    workspace.range_max_indices[base + index] = index;
  }
  for (std::size_t node = base; node-- > 1;) {
    const std::size_t left = node * 2;
    const std::size_t right = left + 1;
    const bool take_right =
        workspace.range_max_values[right] > workspace.range_max_values[left] ||
        (workspace.range_max_values[right] == workspace.range_max_values[left] &&
         workspace.range_max_indices[right] > workspace.range_max_indices[left]);
    const std::size_t selected = take_right ? right : left;
    workspace.range_max_values[node] = workspace.range_max_values[selected];
    workspace.range_max_indices[node] = workspace.range_max_indices[selected];
  }
}

RangeMaximum range_maximum(
    const GeneconvWorkspace& workspace,
    std::size_t base,
    std::size_t beginning,
    std::size_t ending) {
  if (beginning > ending) return {};
  beginning += base;
  ending += base;
  RangeMaximum left_result;
  RangeMaximum right_result;
  while (beginning <= ending) {
    if ((beginning & 1U) != 0) {
      left_result = better_maximum(
          left_result,
          {workspace.range_max_values[beginning],
           workspace.range_max_indices[beginning], true});
      ++beginning;
    }
    if ((ending & 1U) == 0) {
      right_result = better_maximum(
          {workspace.range_max_values[ending],
           workspace.range_max_indices[ending], true},
          right_result);
      if (ending == 0) break;
      --ending;
    }
    beginning /= 2;
    ending /= 2;
  }
  return better_maximum(left_result, right_result);
}

bool solve_source_lambda_k(
    std::size_t discordant_sites,
    std::size_t polymorphic_sites,
    std::size_t mismatch_penalty,
    double& lambda,
    double& k,
    bool& numerical_fallback) {
  lambda = 0.0;
  k = 0.0;
  numerical_fallback = false;
  if (discordant_sites == 0 || discordant_sites >= polymorphic_sites ||
      mismatch_penalty == 0) {
    return false;
  }
  const double p = static_cast<double>(discordant_sites) /
      static_cast<double>(polymorphic_sites);
  const double q = 1.0 - p;
  const double penalty = static_cast<double>(mismatch_penalty);
  const double weighted_p = penalty * p;
  const double initial_log = std::log(weighted_p / q) / (penalty + 1.0);
  double z = std::exp(2.0 * initial_log);
  bool converged = std::isfinite(z) && z > 1.0;
  if (converged) {
    for (std::size_t iteration = 0; iteration < 256; ++iteration) {
      const double inverse = std::pow(z, -penalty);
      const double residual = q * z + p * inverse - 1.0;
      const double derivative = q - weighted_p * inverse / z;
      if (!std::isfinite(residual) || !std::isfinite(derivative) ||
          std::fabs(derivative) < 1.0e-15) {
        converged = false;
        break;
      }
      const double delta = residual / derivative;
      const double next = z - delta;
      if (!std::isfinite(next) || next <= 1.0) {
        converged = false;
        break;
      }
      z = next;
      if (std::fabs(delta) <= 1.0e-6 && std::fabs(residual) <= 1.0e-6) {
        converged = true;
        break;
      }
      converged = iteration + 1 < 256;
    }
  }
  if (!converged || !(z > 1.0)) {
    // The supplied routine is an unbounded Newton loop. Keep its starting
    // point above, but use a bounded bracket only when that path is unstable
    // so a hostile browser dataset cannot hang a worker.
    numerical_fallback = true;
    const auto residual = [&](double candidate) {
      return q * candidate + p * std::pow(candidate, -penalty) - 1.0;
    };
    double low = 1.0 + 1.0e-10;
    double high = 2.0;
    while (residual(high) <= 0.0 && high < 1.0e12) high *= 2.0;
    if (!(residual(high) > 0.0)) return false;
    for (std::size_t iteration = 0; iteration < 160; ++iteration) {
      const double middle = (low + high) / 2.0;
      if (residual(middle) > 0.0) high = middle;
      else low = middle;
    }
    z = (low + high) / 2.0;
  }
  lambda = std::log(z);
  const double first = std::exp(lambda) - 1.0;
  const double tail = std::exp(-(penalty + 1.0) * lambda);
  k = first * (q - weighted_p * tail);
  return std::isfinite(lambda) && lambda > 0.0 &&
      std::isfinite(k) && k > 0.0;
}

std::size_t source_critical_score(
    double p_value_cutoff,
    double lambda,
    double k,
    std::size_t polymorphic_sites) {
  double critical = 4.0;
  if (p_value_cutoff > 0.0 && p_value_cutoff < 1.0 &&
      lambda > 0.0 && k > 0.0) {
    // CalcKMaxP performs these two transforms directly rather than through
    // log1p; retain its rounding route at the critical-score boundary.
    const double first = -std::log(1.0 - p_value_cutoff);
    if (first > 0.0) {
      const double transformed = -std::log(first);
      const double candidate =
          (std::log(k * static_cast<double>(polymorphic_sites)) + transformed) /
          lambda;
      if (std::isfinite(candidate)) critical = std::max(4.0, candidate);
    }
  }
  if (critical >= static_cast<double>(std::numeric_limits<std::size_t>::max())) {
    return std::numeric_limits<std::size_t>::max();
  }
  return static_cast<std::size_t>(critical);
}

double source_ka_probability(
    std::size_t score,
    double lambda,
    double k,
    std::size_t polymorphic_sites) {
  const float log_k_length = static_cast<float>(
      std::log(k * static_cast<double>(polymorphic_sites)));
  float ka_score = static_cast<float>(
      lambda * static_cast<double>(score) - static_cast<double>(log_k_length));
  if (!(ka_score > 0.0F)) return 1.0;
  double probability = 1.0;
  if (ka_score < 32.0F) {
    const double tail = std::exp(-static_cast<double>(ka_score));
    probability = 1.0 - std::exp(-tail);
  } else {
    const float original = ka_score;
    if (ka_score > 700.0F) ka_score = 701.0F;
    probability = std::exp(-static_cast<double>(ka_score));
    if (original > 700.0F) {
      probability /= static_cast<double>(original - 700.0F);
    }
  }
  if (!std::isfinite(probability)) return 1.0;
  // GCCalcPValP2 can underflow to zero. GCGetHiPValP and the normal GCXoverD
  // branch retain that literal value, so do not replace it with a display floor.
  return std::clamp(probability, 0.0, 1.0);
}

void prepare_prefix_queries(
    const std::vector<GeneconvRun>& runs,
    std::size_t mismatch_penalty,
    GeneconvWorkspace& workspace,
    std::size_t& tree_base) {
  workspace.prefix_scores.assign(runs.size() + 1, 0);
  workspace.prefix_positive_sites.assign(runs.size() + 1, 0);
  workspace.prefix_discordant_sites.assign(runs.size() + 1, 0);
  for (std::size_t index = 0; index < runs.size(); ++index) {
    const auto& run = runs[index];
    const std::int64_t weight = run.positive
        ? static_cast<std::int64_t>(run.length)
        : -static_cast<std::int64_t>(mismatch_penalty) *
            static_cast<std::int64_t>(run.length);
    workspace.prefix_scores[index + 1] = workspace.prefix_scores[index] + weight;
    workspace.prefix_positive_sites[index + 1] =
        workspace.prefix_positive_sites[index] + (run.positive ? run.length : 0);
    workspace.prefix_discordant_sites[index + 1] =
        workspace.prefix_discordant_sites[index] + (run.positive ? 0 : run.length);
  }
  workspace.next_lower_prefix.assign(
      workspace.prefix_scores.size(), workspace.prefix_scores.size());
  workspace.monotonic_stack.clear();
  workspace.monotonic_stack.reserve(workspace.prefix_scores.size());
  for (std::size_t index = workspace.prefix_scores.size(); index-- > 0;) {
    while (!workspace.monotonic_stack.empty() &&
           workspace.prefix_scores[workspace.monotonic_stack.back()] >=
               workspace.prefix_scores[index]) {
      workspace.monotonic_stack.pop_back();
    }
    if (!workspace.monotonic_stack.empty()) {
      workspace.next_lower_prefix[index] = workspace.monotonic_stack.back();
    }
    workspace.monotonic_stack.push_back(index);
  }
  build_range_maximum_tree(workspace.prefix_scores, workspace, tree_base);
}

void score_track_fragments(
    std::uint8_t track,
    std::size_t discordant_sites,
    std::size_t mismatch_penalty,
    std::size_t critical_score,
    double lambda,
    double k,
    GeneconvWorkspace& workspace,
    GeneconvDiscoverySummary& summary) {
  const auto& runs = workspace.runs[track];
  if (runs.empty()) return;
  std::size_t tree_base = 0;
  prepare_prefix_queries(
      runs, mismatch_penalty, workspace, tree_base);
  for (std::size_t run_index = 0; run_index < runs.size(); ++run_index) {
    if (!runs[run_index].positive) continue;
    ++summary.fragments_scored;
    const std::size_t maximum_prefix = runs.size();
    const std::size_t first_lower = workspace.next_lower_prefix[run_index];
    const std::size_t final_prefix = first_lower <= maximum_prefix
        ? first_lower - 1
        : maximum_prefix;
    const RangeMaximum maximum = range_maximum(
        workspace, tree_base, run_index + 1, final_prefix);
    if (!maximum.available || maximum.value < workspace.prefix_scores[run_index]) {
      continue;
    }
    const std::size_t ending_run_index = maximum.index - 1;
    const GeneconvRun& beginning_run = runs[run_index];
    const GeneconvRun& ending_run = runs[ending_run_index];
    const std::size_t positive_sites =
        workspace.prefix_positive_sites[maximum.index] -
        workspace.prefix_positive_sites[run_index];
    const std::size_t discordant_fragment_sites =
        workspace.prefix_discordant_sites[maximum.index] -
        workspace.prefix_discordant_sites[run_index];
    // GetMaxFragScoreP narrows MissPen to float and evaluates the cumulative
    // positive/mismatch totals before truncating to an integer score.
    const float source_score = static_cast<float>(positive_sites) -
        static_cast<float>(discordant_fragment_sites) *
            static_cast<float>(mismatch_penalty);
    const std::size_t score = source_score <= 0.0F
        ? 0
        : source_score >=
                static_cast<float>(std::numeric_limits<std::size_t>::max())
            ? std::numeric_limits<std::size_t>::max()
            : static_cast<std::size_t>(source_score);
    if (score <= 3 || score <= critical_score) continue;
    GeneconvScoredFragment fragment;
    fragment.track = track;
    fragment.beginning = beginning_run.beginning;
    fragment.ending = ending_run.ending;
    fragment.wraps_origin = beginning_run.wraps_origin || ending_run.wraps_origin ||
        fragment.beginning > fragment.ending;
    fragment.positive_sites = positive_sites;
    fragment.discordant_sites = discordant_fragment_sites;
    fragment.mismatch_penalty = mismatch_penalty;
    fragment.fragment_score = score;
    fragment.critical_score = critical_score;
    fragment.lambda = lambda;
    fragment.karlin_altschul_k = k;
    fragment.raw_p_value = source_ka_probability(
        score, lambda, k, workspace.categories.size());
    workspace.scored_fragments.push_back(fragment);
    ++summary.qualified_fragments;
  }
  (void)discordant_sites;
}

void overlap_push(GeneconvWorkspace& workspace, std::size_t node) {
  const std::uint32_t lazy = workspace.overlap_tree_lazy[node];
  if (lazy == 0) return;
  for (const std::size_t child : {node * 2, node * 2 + 1}) {
    workspace.overlap_tree_max[child] += lazy;
    workspace.overlap_tree_lazy[child] += lazy;
  }
  workspace.overlap_tree_lazy[node] = 0;
}

std::uint32_t overlap_range_maximum(
    GeneconvWorkspace& workspace,
    std::size_t node,
    std::size_t beginning,
    std::size_t ending,
    std::size_t query_beginning,
    std::size_t query_ending) {
  if (query_beginning <= beginning && ending <= query_ending) {
    return workspace.overlap_tree_max[node];
  }
  overlap_push(workspace, node);
  const std::size_t middle = beginning + (ending - beginning) / 2;
  std::uint32_t result = 0;
  if (query_beginning <= middle) {
    result = overlap_range_maximum(
        workspace, node * 2, beginning, middle,
        query_beginning, query_ending);
  }
  if (query_ending > middle) {
    result = std::max(
        result,
        overlap_range_maximum(
            workspace, node * 2 + 1, middle + 1, ending,
            query_beginning, query_ending));
  }
  return result;
}

void overlap_range_add(
    GeneconvWorkspace& workspace,
    std::size_t node,
    std::size_t beginning,
    std::size_t ending,
    std::size_t query_beginning,
    std::size_t query_ending) {
  if (query_beginning <= beginning && ending <= query_ending) {
    ++workspace.overlap_tree_max[node];
    ++workspace.overlap_tree_lazy[node];
    return;
  }
  overlap_push(workspace, node);
  const std::size_t middle = beginning + (ending - beginning) / 2;
  if (query_beginning <= middle) {
    overlap_range_add(
        workspace, node * 2, beginning, middle,
        query_beginning, query_ending);
  }
  if (query_ending > middle) {
    overlap_range_add(
        workspace, node * 2 + 1, middle + 1, ending,
        query_beginning, query_ending);
  }
  workspace.overlap_tree_max[node] = std::max(
      workspace.overlap_tree_max[node * 2],
      workspace.overlap_tree_max[node * 2 + 1]);
}

std::uint32_t fragment_maximum_coverage(
    const GeneconvScoredFragment& fragment,
    std::size_t length,
    GeneconvWorkspace& workspace) {
  if (length == 0) return std::numeric_limits<std::uint32_t>::max();
  if (fragment.wraps_origin || fragment.beginning > fragment.ending) {
    return std::max(
        overlap_range_maximum(
            workspace, 1, 0, length - 1, fragment.beginning, length - 1),
        overlap_range_maximum(
            workspace, 1, 0, length - 1, 0, fragment.ending));
  }
  return overlap_range_maximum(
      workspace, 1, 0, length - 1, fragment.beginning, fragment.ending);
}

void add_fragment_coverage(
    const GeneconvScoredFragment& fragment,
    std::size_t length,
    GeneconvWorkspace& workspace) {
  if (fragment.wraps_origin || fragment.beginning > fragment.ending) {
    overlap_range_add(
        workspace, 1, 0, length - 1, fragment.beginning, length - 1);
    overlap_range_add(
        workspace, 1, 0, length - 1, 0, fragment.ending);
  } else {
    overlap_range_add(
        workspace, 1, 0, length - 1, fragment.beginning, fragment.ending);
  }
}

std::size_t beginning_coordinate(
    const MaxChiWorkspace& variable_workspace,
    std::size_t informative_beginning,
    std::size_t alignment_length) {
  if (informative_beginning == 0) return 1;
  const std::size_t previous = variable_workspace.coordinates[informative_beginning - 1];
  return std::min(alignment_length, previous + 1);
}

std::size_t ending_coordinate(
    const MaxChiWorkspace& variable_workspace,
    std::size_t informative_ending,
    std::size_t alignment_length) {
  if (informative_ending + 1 >= variable_workspace.coordinates.size()) {
    return alignment_length > 1 ? alignment_length - 1 : alignment_length;
  }
  const std::size_t next = variable_workspace.coordinates[informative_ending + 1];
  return next > 1 ? next - 1 : 1;
}

void prepare_variable_workspace(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    MaxChiWorkspace& variable_workspace,
    std::array<double, 3>& pair_similarity) {
  variable_workspace.coordinates.clear();
  for (auto& matches : variable_workspace.matches) matches.clear();
  variable_workspace.coordinates.reserve(alignment.length);
  for (auto& matches : variable_workspace.matches) matches.reserve(alignment.length);
  std::array<std::size_t, 3> comparable{};
  std::array<std::size_t, 3> identical{};
  for (std::size_t position = 0; position < alignment.length; ++position) {
    const std::uint8_t first = alignment.at(triplet[0], position);
    const std::uint8_t second = alignment.at(triplet[1], position);
    const std::uint8_t third = alignment.at(triplet[2], position);
    if (first != 0 && second != 0) {
      ++comparable[0];
      if (first == second) ++identical[0];
    }
    if (first != 0 && third != 0) {
      ++comparable[1];
      if (first == third) ++identical[1];
    }
    if (second != 0 && third != 0) {
      ++comparable[2];
      if (second == third) ++identical[2];
    }
    if (first == 0 || second == 0 || third == 0 ||
        (first == second && first == third)) {
      continue;
    }
    variable_workspace.coordinates.push_back(position + 1);
    variable_workspace.matches[0].push_back(first == second ? 1 : 0);
    variable_workspace.matches[1].push_back(first == third ? 1 : 0);
    variable_workspace.matches[2].push_back(second == third ? 1 : 0);
  }
  for (std::size_t pair = 0; pair < 3; ++pair) {
    pair_similarity[pair] = comparable[pair] == 0
        ? 0.0
        : static_cast<double>(identical[pair]) /
            static_cast<double>(comparable[pair]);
  }
}

}  // namespace

GeneconvDiscoverySummary geneconv_discover_prepared(
    const MaxChiWorkspace& variable_workspace,
    const std::array<double, 3>& pair_similarity,
    const GeneconvDiscoveryOptions& options,
    GeneconvWorkspace& workspace,
    std::vector<GeneconvDiscoveryCandidate>& output) {
  GeneconvDiscoverySummary summary;
  summary.bonferroni_applied = options.bonferroni;
  summary.correction_tests = std::max<std::uint64_t>(1, options.correction_tests);
  summary.polymorphic_sites = variable_workspace.coordinates.size();
  output.clear();
  workspace.scored_fragments.clear();
  workspace.category_runs.clear();
  workspace.overlap_tree_max.clear();
  workspace.overlap_tree_lazy.clear();
  for (auto& runs : workspace.runs) runs.clear();
  if (variable_workspace.coordinates.empty() ||
      variable_workspace.matches[0].size() != variable_workspace.coordinates.size() ||
      variable_workspace.matches[1].size() != variable_workspace.coordinates.size() ||
      variable_workspace.matches[2].size() != variable_workspace.coordinates.size()) {
    return summary;
  }

  std::array<std::size_t, 4> category_counts{};
  build_categories(variable_workspace, workspace, category_counts);
  const std::size_t length = workspace.categories.size();
  if (category_counts[0] == length || category_counts[1] == length ||
      category_counts[2] == length) {
    return summary;
  }
  std::array<std::size_t, 6> discordant_sites{{
      length - category_counts[0],
      length - category_counts[1],
      length - category_counts[2],
      category_counts[0] + category_counts[1],
      category_counts[0] + category_counts[2],
      category_counts[1] + category_counts[2],
  }};
  const auto [minimum_difference, maximum_difference] = std::minmax_element(
      discordant_sites.begin(), discordant_sites.end());
  if (*minimum_difference < 3 && *maximum_difference > *minimum_difference * 10) {
    summary.source_skew_filter_rejected = true;
    return summary;
  }
  summary.profile_available = true;
  if (!build_source_runs(options.circular, workspace)) return summary;
  const float scaled_length = static_cast<float>(length) *
      static_cast<float>(std::max<std::size_t>(1, options.mismatch_scale));
  const double raw_cutoff = options.bonferroni
      ? options.p_value_cutoff / static_cast<double>(summary.correction_tests)
      : options.p_value_cutoff;

  for (std::uint8_t track = 0; track < 6; ++track) {
    if (track >= 3 && discordant_sites[track] == 0) discordant_sites[track] = 1;
    ++summary.tracks_screened;
    const std::size_t mismatch_penalty = discordant_sites[track] == 0
        ? 0
        : static_cast<std::size_t>(
              scaled_length / static_cast<float>(discordant_sites[track])) + 1;
    double lambda = 0.0;
    double k = 0.0;
    bool numerical_fallback = false;
    if (!solve_source_lambda_k(
            discordant_sites[track], length, mismatch_penalty,
            lambda, k, numerical_fallback)) {
      continue;
    }
    if (numerical_fallback) ++summary.numerical_fallback_tracks;
    const std::size_t critical_score = source_critical_score(
        raw_cutoff, lambda, k, length);
    score_track_fragments(
        track,
        discordant_sites[track],
        mismatch_penalty,
        critical_score,
        lambda,
        k,
        workspace,
        summary);
  }

  std::stable_sort(
      workspace.scored_fragments.begin(),
      workspace.scored_fragments.end(),
      [](const GeneconvScoredFragment& left, const GeneconvScoredFragment& right) {
        return left.raw_p_value < right.raw_p_value;
      });
  workspace.overlap_tree_max.assign(length * 4 + 4, 0);
  workspace.overlap_tree_lazy.assign(length * 4 + 4, 0);
  const std::size_t overlap_limit = std::max<std::size_t>(
      1, options.maximum_overlapping_fragments);
  const std::size_t alignment_length = variable_workspace.variable_prefix.empty()
      ? variable_workspace.coordinates.back()
      : variable_workspace.variable_prefix.size() - 1;
  for (const auto& fragment : workspace.scored_fragments) {
    if (!(fragment.raw_p_value < raw_cutoff)) continue;
    if (fragment_maximum_coverage(fragment, length, workspace) >= overlap_limit) {
      ++summary.overlap_rejected_fragments;
      continue;
    }
    add_fragment_coverage(fragment, length, workspace);
    GeneconvDiscoveryCandidate candidate;
    candidate.track = fragment.track;
    candidate.recombinant_local = kRecombinantLocal[fragment.track];
    candidate.major_parent_local = kMajorParentLocal[fragment.track];
    candidate.minor_parent_local = kMinorParentLocal[fragment.track];
    candidate.candidate_pair = kCandidatePair[fragment.track];
    candidate.beginning = beginning_coordinate(
        variable_workspace, fragment.beginning, alignment_length);
    candidate.ending = ending_coordinate(
        variable_workspace, fragment.ending, alignment_length);
    candidate.wraps_origin = fragment.wraps_origin ||
        (options.circular && candidate.beginning >= candidate.ending);
    candidate.informative_beginning = fragment.beginning + 1;
    candidate.informative_ending = fragment.ending + 1;
    candidate.polymorphic_sites = length;
    candidate.positive_sites = fragment.positive_sites;
    candidate.discordant_sites = fragment.discordant_sites;
    candidate.mismatch_penalty = fragment.mismatch_penalty;
    candidate.fragment_score = fragment.fragment_score;
    candidate.critical_score = fragment.critical_score;
    candidate.lambda = fragment.lambda;
    candidate.karlin_altschul_k = fragment.karlin_altschul_k;
    candidate.raw_p_value = fragment.raw_p_value;
    candidate.corrected_p_value = options.bonferroni
        ? std::min(1.0, fragment.raw_p_value *
              static_cast<double>(summary.correction_tests))
        : fragment.raw_p_value;
    candidate.pair_similarity = pair_similarity;
    output.push_back(candidate);
  }
  summary.emitted_candidates = output.size();
  return summary;
}

GeneconvRecheckEvidence geneconv_recheck_prepared(
    const MaxChiWorkspace& variable_workspace,
    const std::array<double, 3>& pair_similarity,
    const GeneconvDiscoveryOptions& options,
    GeneconvWorkspace& workspace,
    std::vector<GeneconvDiscoveryCandidate>& candidates) {
  GeneconvRecheckEvidence evidence;
  evidence.requested = true;
  const GeneconvDiscoverySummary summary = geneconv_discover_prepared(
      variable_workspace, pair_similarity, options, workspace, candidates);
  evidence.profile_available = summary.profile_available;
  evidence.bonferroni_applied = summary.bonferroni_applied;
  evidence.ignored_indels = summary.ignored_indels;
  evidence.overlap_filter_applied = summary.overlap_filter_applied;
  evidence.minimum_fragment_filters_applied =
      summary.minimum_fragment_filters_applied;
  evidence.source_skew_filter_rejected = summary.source_skew_filter_rejected;
  evidence.correction_tests = summary.correction_tests;
  evidence.polymorphic_sites = summary.polymorphic_sites;
  evidence.tracks_screened = summary.tracks_screened;
  evidence.fragments_scored = summary.fragments_scored;
  evidence.qualified_fragments = summary.qualified_fragments;
  evidence.overlap_rejected_fragments = summary.overlap_rejected_fragments;
  evidence.numerical_fallback_tracks = summary.numerical_fallback_tracks;

  if (!candidates.empty()) {
    // geneconv_discover_prepared preserves the supplied stable lowest-raw-P
    // ordering, so the first surviving overlap-eligible fragment is the
    // source-shaped recheck representative.
    const auto& best = candidates.front();
    evidence.best_track = static_cast<std::int8_t>(best.track);
    evidence.recombinant_local =
        static_cast<std::int8_t>(best.recombinant_local);
    evidence.beginning = best.beginning;
    evidence.ending = best.ending;
    evidence.wraps_origin = best.wraps_origin;
    evidence.fragment_score = best.fragment_score;
    evidence.critical_score = best.critical_score;
    evidence.raw_p_value = best.raw_p_value;
    evidence.corrected_p_value = best.corrected_p_value;
    evidence.source_recheck_hit =
        best.corrected_p_value < options.p_value_cutoff;
  }
  return evidence;
}

GeneconvPlotProfile geneconv_plot_profile(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    const GeneconvDiscoveryOptions& options,
    GeneconvWorkspace& workspace) {
  GeneconvPlotProfile plot;
  MaxChiWorkspace variable_workspace;
  std::array<double, 3> pair_similarity{};
  prepare_variable_workspace(
      alignment, triplet, variable_workspace, pair_similarity);
  variable_workspace.variable_prefix.assign(alignment.length + 1, 0);
  std::size_t variable_index = 0;
  for (std::size_t coordinate = 1; coordinate <= alignment.length; ++coordinate) {
    while (variable_index < variable_workspace.coordinates.size() &&
           variable_workspace.coordinates[variable_index] <= coordinate) {
      ++variable_index;
    }
    variable_workspace.variable_prefix[coordinate] = variable_index;
  }
  std::vector<GeneconvDiscoveryCandidate> ignored;
  const auto summary = geneconv_discover_prepared(
      variable_workspace, pair_similarity, options, workspace, ignored);
  if (!summary.profile_available) return plot;
  plot.available = true;
  plot.coordinates = variable_workspace.coordinates;
  for (auto& values : plot.negative_log10_p_value) {
    values.assign(plot.coordinates.size(), 0.0);
  }
  struct PlotInterval {
    std::size_t beginning = 0;
    std::size_t ending = 0;
    double value = 0.0;
  };
  std::array<std::vector<PlotInterval>, 3> intervals;
  for (const auto& fragment : workspace.scored_fragments) {
    const double value = -std::log10(std::max(kMinimumProbability, fragment.raw_p_value));
    const std::size_t pair = kPlotPair[fragment.track];
    if (fragment.wraps_origin || fragment.beginning > fragment.ending) {
      intervals[pair].push_back({
          fragment.beginning, plot.coordinates.size() - 1, value});
      intervals[pair].push_back({0, fragment.ending, value});
    } else {
      intervals[pair].push_back({fragment.beginning, fragment.ending, value});
    }
  }
  for (std::size_t pair = 0; pair < intervals.size(); ++pair) {
    auto& pair_intervals = intervals[pair];
    std::sort(
        pair_intervals.begin(), pair_intervals.end(),
        [](const PlotInterval& left, const PlotInterval& right) {
          return left.beginning < right.beginning;
        });
    std::priority_queue<std::pair<double, std::size_t>> active;
    std::size_t interval = 0;
    for (std::size_t position = 0; position < plot.coordinates.size(); ++position) {
      while (interval < pair_intervals.size() &&
             pair_intervals[interval].beginning <= position) {
        active.emplace(
            pair_intervals[interval].value, pair_intervals[interval].ending);
        ++interval;
      }
      while (!active.empty() && active.top().second < position) active.pop();
      if (!active.empty()) {
        plot.negative_log10_p_value[pair][position] = active.top().first;
      }
    }
  }
  return plot;
}

}  // namespace rdp
