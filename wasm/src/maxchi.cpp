#include "maxchi.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace rdp {
namespace {

double source_normal_z(double z) {
  // NormalZ in the supplied DNA5 DLL. Retaining the polynomial matters here:
  // ChiPVal2P is used both for the initial MaxChi gate and the grown window.
  long double x = 0.0L;
  long double y = 0.0L;
  if (std::fabs(z) < 5.9999999) {
    if (z == 0.0) return 0.0;
    y = 0.5L * std::fabs(static_cast<long double>(z));
    if (y >= 3.0L) {
      x = 1.0L;
    } else if (y < 1.0L) {
      const long double w = y * y;
      x = ((((((((0.000124818987L * w - 0.001075204047L) * w +
                       0.005198775019L) *
                          w -
                      0.019198292004L) *
                         w +
                     0.059054035642L) *
                        w -
                    0.151968751364L) *
                       w +
                   0.319152932694L) *
                      w -
                  0.5319230073L) *
                     w +
                 0.797884560593L) *
          y * 2.0L;
    } else {
      y -= 2.0L;
      x = (((((((((((((-0.000045255659L * y + 0.00015252929L) * y -
                            0.000019538132L) *
                               y -
                           0.000676904986L) *
                              y +
                          0.001390604284L) *
                             y -
                         0.00079462082L) *
                            y -
                        0.002034254874L) *
                           y +
                       0.006549791214L) *
                          y -
                      0.010557625006L) *
                         y +
                     0.011630447319L) *
                        y -
                    0.009279453341L) *
                       y +
                   0.005353579108L) *
                      y -
                  0.002141268741L) *
                     y +
                 0.000535310849L) *
                    y +
                0.999936657524L;
    }
    return static_cast<double>(std::min(x + 1.0L, 1.0L - x));
  }
  const long double exponent = (std::fabs(static_cast<long double>(z)) - 5.999999L) * 10.0L;
  return static_cast<double>(1.0e-9L / std::pow(1.6L, exponent));
}

double source_chi_p_value(double chi_square) {
  if (!(chi_square > 0.0)) return 1.0;
  long double probability = source_normal_z(-std::sqrt(chi_square));
  if (probability == 0.0L) {
    probability = 1.0e-10L / (static_cast<long double>(chi_square) - 34.0L);
  }
  return std::clamp(static_cast<double>(probability), 0.0, 1.0);
}

double chi_square(
    std::size_t half_window,
    std::size_t left_matches,
    std::size_t right_matches) {
  const double a = static_cast<double>(left_matches);
  const double c = static_cast<double>(right_matches);
  const double h = static_cast<double>(half_window);
  const double b = h - a;
  const double d = h - c;
  if (!(a + c > 0.0) || !(b + d > 0.0)) return 0.0;
  const double cross = a * d - b * c;
  return (cross * cross * 2.0) / (h * (a + c) * (b + d));
}

std::size_t source_critical_difference(
    std::size_t half_window,
    double screening_probability) {
  // GetCriticalDiff derives a cheap absolute match-count screen from the chi
  // cutoff, then subtracts one because CalcChiVals uses a strict comparison.
  screening_probability = std::max(0.0001, screening_probability);
  double low = 0.0;
  double high = std::max<double>(5.0, static_cast<double>(half_window) * 2.0);
  while (source_chi_p_value(high) > screening_probability && high < 1.0e6) {
    high *= 2.0;
  }
  for (std::size_t iteration = 0; iteration < 80; ++iteration) {
    const double middle = (low + high) / 2.0;
    if (source_chi_p_value(middle) < screening_probability) high = middle;
    else low = middle;
  }
  const double critical_chi = (low + high) / 2.0;
  for (std::size_t difference = 1; difference <= half_window; ++difference) {
    if (chi_square(half_window, 0, difference) > critical_chi) {
      return difference - 1;
    }
  }
  return half_window > 0 ? half_window - 1 : 0;
}

void build_variable_profile(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    MaxChiWorkspace& workspace) {
  workspace.coordinates.clear();
  workspace.coordinates.reserve(alignment.length);
  for (auto& pair : workspace.matches) {
    pair.clear();
    pair.reserve(alignment.length);
  }
  workspace.variable_prefix.assign(alignment.length + 1, 0);
  for (std::size_t coordinate = 1; coordinate <= alignment.length; ++coordinate) {
    const std::uint8_t first = alignment.at(triplet[0], coordinate - 1);
    const std::uint8_t second = alignment.at(triplet[1], coordinate - 1);
    const std::uint8_t third = alignment.at(triplet[2], coordinate - 1);
    if (first != 0 && second != 0 && third != 0 &&
        (first != second || first != third)) {
      workspace.coordinates.push_back(coordinate);
      workspace.matches[0].push_back(first == second ? 1 : 0);
      workspace.matches[1].push_back(first == third ? 1 : 0);
      workspace.matches[2].push_back(second == third ? 1 : 0);
    }
    workspace.variable_prefix[coordinate] = workspace.coordinates.size();
  }
}

void make_banned_windows(
    MaxChiWorkspace& workspace,
    const std::vector<std::uint8_t>& triplet_missing_data,
    std::size_t half_window,
    bool circular,
    bool& filter_applied) {
  const std::size_t variable_sites = workspace.coordinates.size();
  workspace.banned_windows.assign(variable_sites + 1, 0);
  workspace.missing_boundaries.assign(variable_sites + 1, 0);
  auto& banned = workspace.banned_windows;
  auto& missing_boundaries = workspace.missing_boundaries;
  filter_applied = false;
  if (triplet_missing_data.size() == workspace.variable_prefix.size() - 1) {
    for (std::size_t coordinate = 1;
         coordinate < workspace.variable_prefix.size();
         ++coordinate) {
      if (triplet_missing_data[coordinate - 1] == 0) continue;
      filter_applied = true;
      const std::size_t mapped = workspace.variable_prefix[coordinate];
      if (variable_sites == 0) continue;
      missing_boundaries[mapped] = 1;
      if (mapped + half_window - 1 <= variable_sites) {
        for (std::size_t position = mapped;
             position < mapped + half_window;
             ++position) {
          banned[position] = 1;
        }
      } else {
        for (std::size_t position = mapped; position <= variable_sites; ++position) {
          banned[position] = 1;
        }
        const std::size_t wrapped_end = mapped + half_window - 1 - variable_sites;
        for (std::size_t position = 0; position < wrapped_end; ++position) {
          banned[position] = 1;
        }
      }
      if (mapped < variable_sites) {
        missing_boundaries[mapped + 1] = 1;
        if (mapped + 2 > half_window) {
          for (std::size_t position = mapped + 2 - half_window;
               position < mapped + 2;
               ++position) {
            banned[position] = 1;
          }
        } else {
          for (std::size_t position = 0; position < mapped + 2; ++position) {
            banned[position] = 1;
          }
          for (std::size_t position = mapped + 2 + variable_sites - half_window;
               position <= variable_sites;
               ++position) {
            banned[position] = 1;
          }
        }
      } else {
        missing_boundaries[1] = 1;
      }
    }
  }
  if (variable_sites > 0 &&
      (missing_boundaries[variable_sites] != 0 || missing_boundaries[1] != 0)) {
    const std::size_t beginning = variable_sites >= half_window
        ? variable_sites - half_window + 2
        : 1;
    for (std::size_t position = beginning; position <= variable_sites; ++position) {
      banned[position] = 1;
    }
  }
  if (!circular && variable_sites > 0) {
    missing_boundaries[1] = 1;
    missing_boundaries[variable_sites] = 1;
    const std::size_t beginning = variable_sites >= half_window
        ? variable_sites - half_window + 2
        : 1;
    for (std::size_t position = beginning; position <= variable_sites; ++position) {
      banned[position] = 1;
    }
  }
}

std::size_t circular_sum(
    const std::vector<std::uint8_t>& values,
    std::ptrdiff_t beginning,
    std::size_t count) {
  std::size_t total = 0;
  const std::size_t length = values.size();
  for (std::size_t offset = 0; offset < count; ++offset) {
    std::ptrdiff_t index = beginning + static_cast<std::ptrdiff_t>(offset);
    index %= static_cast<std::ptrdiff_t>(length);
    if (index < 0) index += static_cast<std::ptrdiff_t>(length);
    total += values[static_cast<std::size_t>(index)];
  }
  return total;
}

bool boundary_banned(
    const std::vector<std::uint8_t>& banned,
    std::size_t boundary,
    std::size_t half_window) {
  const std::size_t length = banned.size() - 1;
  const std::size_t native_boundary = boundary;
  std::ptrdiff_t other = static_cast<std::ptrdiff_t>(native_boundary) -
      static_cast<std::ptrdiff_t>(half_window);
  if (other < 1) other += static_cast<std::ptrdiff_t>(length);
  const std::size_t other_index = static_cast<std::size_t>(other);
  return banned[native_boundary] != 0 || banned[other_index] != 0;
}

struct Peak {
  std::size_t boundary = 0;
  std::size_t pair = 0;
  std::size_t half_window = 0;
  double chi = 0.0;
};

Peak strongest_peak(
    const MaxChiWorkspace& profile,
    const std::vector<std::uint8_t>& banned,
    std::size_t half_window,
    std::size_t critical_difference) {
  Peak best;
  const std::size_t length = profile.coordinates.size();
  std::array<std::size_t, 3> left_matches{};
  std::array<std::size_t, 3> right_matches{};
  for (std::size_t pair = 0; pair < profile.matches.size(); ++pair) {
    left_matches[pair] = circular_sum(
        profile.matches[pair],
        -static_cast<std::ptrdiff_t>(half_window),
        half_window);
    right_matches[pair] = circular_sum(profile.matches[pair], 0, half_window);
  }
  for (std::size_t boundary = 0; boundary < length; ++boundary) {
    if (!boundary_banned(banned, boundary, half_window)) {
      for (std::size_t pair = 0; pair < profile.matches.size(); ++pair) {
        const std::size_t left = left_matches[pair];
        const std::size_t right = right_matches[pair];
        const std::ptrdiff_t difference =
            static_cast<std::ptrdiff_t>(left) - static_cast<std::ptrdiff_t>(right);
        if (difference > static_cast<std::ptrdiff_t>(critical_difference) ||
            difference < -static_cast<std::ptrdiff_t>(critical_difference)) {
          const double value = chi_square(half_window, left, right);
          if (value > best.chi) {
            best = {boundary, pair, half_window, value};
          }
        }
      }
    }
    for (std::size_t pair = 0; pair < profile.matches.size(); ++pair) {
      const auto& matches = profile.matches[pair];
      const std::size_t leaving_left =
          (boundary + length - half_window) % length;
      const std::size_t entering_left = boundary;
      const std::size_t leaving_right = boundary;
      const std::size_t entering_right = (boundary + half_window) % length;
      left_matches[pair] = left_matches[pair] - matches[leaving_left] +
          matches[entering_left];
      right_matches[pair] = right_matches[pair] - matches[leaving_right] +
          matches[entering_right];
    }
  }
  return best;
}

void grow_peak(
    const MaxChiWorkspace& profile,
    const std::vector<std::uint8_t>& missing_boundaries,
    Peak& peak) {
  if (!(peak.chi > 0.0)) return;
  const std::size_t length = profile.coordinates.size();
  const auto& scores = profile.matches[peak.pair];
  const std::size_t source_boundary = peak.boundary == 0 ? 1 : peak.boundary;
  std::size_t test_window = static_cast<std::size_t>(
      static_cast<double>(peak.half_window) / 4.0 + 0.51);
  test_window = std::clamp<std::size_t>(test_window, 6, peak.half_window);
  test_window = std::min(test_window, length / 2);
  if (test_window == 0) return;

  std::size_t left_matches = circular_sum(
      scores,
      static_cast<std::ptrdiff_t>(source_boundary) -
          static_cast<std::ptrdiff_t>(test_window),
      test_window);
  std::size_t right_matches = circular_sum(scores, source_boundary, test_window);
  std::size_t maximum_failures = peak.half_window * 2;
  maximum_failures = std::min(maximum_failures, (length - test_window * 2) / 2);
  if (maximum_failures == 0) maximum_failures = 1;

  std::size_t failures = 0;
  std::size_t current_window = test_window + 1;
  while (failures <= maximum_failures && current_window * 2 <= length) {
    const std::size_t left_index = static_cast<std::size_t>(
        (static_cast<std::ptrdiff_t>(source_boundary) -
         static_cast<std::ptrdiff_t>(current_window) +
         static_cast<std::ptrdiff_t>(length)) %
        static_cast<std::ptrdiff_t>(length));
    const std::size_t right_index =
        (source_boundary + current_window - 1) % length;
    const std::size_t left_native = left_index + 1;
    const std::size_t right_native = right_index + 1;
    left_matches += scores[left_index];
    right_matches += scores[right_index];
    const double value = chi_square(current_window, left_matches, right_matches);
    if (value >= peak.chi) {
      peak.chi = value;
      peak.half_window = current_window;
      failures = 0;
    } else {
      ++failures;
    }
    // GrowMChiWin2P2 evaluates the newly enlarged window before consulting
    // MDMap at its new edges. Preserve that ordering so a half-window may end
    // on a missing-data boundary but cannot traverse it.
    if (missing_boundaries[left_native] != 0 ||
        missing_boundaries[right_native] != 0) {
      break;
    }
    ++current_window;
    if (current_window * 2 <= length) {
      const std::size_t next_left = static_cast<std::size_t>(
          (static_cast<std::ptrdiff_t>(source_boundary) -
           static_cast<std::ptrdiff_t>(current_window) +
           static_cast<std::ptrdiff_t>(length)) %
          static_cast<std::ptrdiff_t>(length));
      const std::size_t next_right =
          (source_boundary + current_window - 1) % length;
      if (missing_boundaries[next_left + 1] != 0 ||
          missing_boundaries[next_right + 1] != 0) {
        break;
      }
    }
  }
}

}  // namespace

MaxChiRecheckEvidence maxchi_recheck(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    const std::vector<std::uint8_t>& triplet_missing_data,
    const MaxChiRecheckOptions& options,
    MaxChiWorkspace& workspace) {
  MaxChiRecheckEvidence evidence;
  evidence.requested = true;
  evidence.bonferroni_applied = options.bonferroni;
  evidence.correction_tests = options.bonferroni
      ? std::max<std::uint64_t>(1, options.correction_tests)
      : 1;
  if (alignment.length == 0 ||
      std::any_of(triplet.begin(), triplet.end(), [&](std::uint32_t sequence) {
        return sequence >= alignment.sequence_count();
      })) {
    return evidence;
  }

  build_variable_profile(alignment, triplet, workspace);
  evidence.variable_sites = workspace.coordinates.size();
  if (evidence.variable_sites < 7) return evidence;

  std::size_t half_window = static_cast<std::size_t>(
      static_cast<double>(options.fixed_window_sites) / 2.0 + 0.51);
  // GetCriticalDiff in the supplied RDP5 VB workflow uses LowestProb / 6 for
  // the MaxChi match-difference pre-screen. Project-wide correction is applied
  // later to the source xMPV value, not to this cheap gate.
  const double screening_probability = options.p_value_cutoff / 6.0;
  std::size_t critical_difference =
      source_critical_difference(half_window, screening_probability);
  if (evidence.variable_sites < critical_difference * 2) return evidence;
  if (half_window * 2 > evidence.variable_sites) {
    half_window = static_cast<std::size_t>(
        static_cast<double>(evidence.variable_sites) * 0.75 / 2.0 + 0.51);
    if (half_window > 0) --half_window;
  }
  if (half_window <= critical_difference) {
    half_window = static_cast<std::size_t>(
        static_cast<double>(evidence.variable_sites) / 2.0 + 0.51);
    if (half_window > 0) --half_window;
  }
  if (half_window < 6 || half_window * 2 > evidence.variable_sites) return evidence;

  evidence.half_window = half_window;
  evidence.critical_difference = critical_difference;
  bool filter_applied = false;
  make_banned_windows(
      workspace,
      triplet_missing_data,
      half_window,
      options.circular,
      filter_applied);
  evidence.missing_data_window_filter_applied = filter_applied;
  evidence.linear_edge_window_filter_applied = !options.circular;
  Peak peak = strongest_peak(
      workspace,
      workspace.banned_windows,
      half_window,
      critical_difference);
  evidence.profile_available = true;
  if (!(peak.chi > 0.0)) return evidence;

  const double initial_tail = source_chi_p_value(peak.chi);
  const double initial_within = std::min(
      1.0,
      initial_tail * static_cast<double>(evidence.variable_sites) /
          static_cast<double>(half_window) * 3.0);
  if (initial_within <= options.p_value_cutoff && initial_tail < 1.0) {
    grow_peak(workspace, workspace.missing_boundaries, peak);
  }

  evidence.best_pair = static_cast<std::int8_t>(peak.pair);
  evidence.grown_half_window = peak.half_window;
  evidence.maximum_chi_square = peak.chi;
  const std::size_t peak_coordinate_index = peak.boundary == 0
      ? workspace.coordinates.size() - 1
      : peak.boundary - 1;
  evidence.peak_alignment_position = workspace.coordinates[peak_coordinate_index];
  evidence.local_p_value = source_chi_p_value(peak.chi);
  const std::size_t probability_window = std::min(half_window, peak.half_window);
  evidence.within_triplet_p_value = std::min(
      1.0,
      evidence.local_p_value * static_cast<double>(evidence.variable_sites) /
          static_cast<double>(probability_window) * 3.0);
  evidence.corrected_p_value = options.bonferroni
      ? std::min(
          1.0,
          evidence.within_triplet_p_value *
              static_cast<double>(evidence.correction_tests))
      : evidence.within_triplet_p_value;
  evidence.source_recheck_hit =
      evidence.corrected_p_value < options.p_value_cutoff;
  return evidence;
}

}  // namespace rdp
