#include "maxchi.hpp"
#include "chimaera.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <queue>
#include <tuple>
#include <utility>
#include <vector>

namespace rdp {
namespace {

constexpr double kMinimumProbability = 1.0e-300;

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

void calculate_chi_profiles(
    const MaxChiWorkspace& profile,
    const std::vector<std::uint8_t>& banned,
    std::size_t half_window,
    std::size_t critical_difference,
    std::array<std::vector<double>, 3>& chi_values) {
  const std::size_t length = profile.coordinates.size();
  // ChiVals is backed by LenSeq + 1 cells per pair in the DLL. Keep the
  // LenXoverSeq cell as zero padding because SmoothChiValsP reads it.
  for (auto& values : chi_values) values.assign(length + 1, 0.0);
  if (length == 0 || half_window == 0 || half_window * 2 > length) return;

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
        const std::ptrdiff_t difference =
            static_cast<std::ptrdiff_t>(left_matches[pair]) -
            static_cast<std::ptrdiff_t>(right_matches[pair]);
        if (difference > static_cast<std::ptrdiff_t>(critical_difference) ||
            difference < -static_cast<std::ptrdiff_t>(critical_difference)) {
          chi_values[pair][boundary] = chi_square(
              half_window, left_matches[pair], right_matches[pair]);
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
}

double source_profile_value(const std::vector<double>& values, std::ptrdiff_t index) {
  // SmoothChiValsP operates on a LenSeq-sized backing array while only the
  // first LenXoverSeq elements are populated. Its literal index LenXoverSeq
  // therefore reads the zero-filled padding cell rather than wrapping to 0.
  const std::ptrdiff_t length = static_cast<std::ptrdiff_t>(values.size());
  if (index < 0) index += length;
  if (index < 0 || index >= length) return 0.0;
  return values[static_cast<std::size_t>(index)];
}

void source_smooth_chi(
    const std::array<std::vector<double>, 3>& chi_values,
    std::array<std::vector<double>, 3>& smooth_chi) {
  constexpr std::ptrdiff_t q_window = 5;
  constexpr std::ptrdiff_t divisor = q_window * 2 + 1;
  for (std::size_t pair = 0; pair < chi_values.size(); ++pair) {
    const auto& source = chi_values[pair];
    auto& smooth = smooth_chi[pair];
    smooth.assign(source.size(), 0.0);
    const std::ptrdiff_t length =
        static_cast<std::ptrdiff_t>(source.size()) - 1;
    if (length == 0) continue;
    double running = 0.0;
    // Literal SmoothChiValsP bounds: -5 through +6, divided by 11. The
    // off-by-one is part of the supplied implementation and is intentional.
    for (std::ptrdiff_t position = -q_window;
         position <= 1 + q_window;
         ++position) {
      running += position < 1
          ? source_profile_value(source, length + position)
          : source_profile_value(source, position);
    }
    smooth[0] = running / static_cast<double>(divisor);
    for (std::ptrdiff_t position = 1 - q_window;
         position < length - q_window;
         ++position) {
      if (position > 0) {
        const std::ptrdiff_t incoming = position + divisor <= length
            ? position + divisor
            : position + divisor - length;
        running = running - source_profile_value(source, position) +
            source_profile_value(source, incoming);
      } else {
        running = running - source_profile_value(source, length + position) +
            source_profile_value(source, position + divisor);
      }
      const std::ptrdiff_t target = position + q_window;
      if (target >= 0 && target < length) {
        smooth[static_cast<std::size_t>(target)] =
            running / static_cast<double>(divisor);
      }
    }
  }
}

std::size_t native_wrap(std::ptrdiff_t position, std::size_t length) {
  if (length == 0) return 0;
  const std::ptrdiff_t native_length = static_cast<std::ptrdiff_t>(length);
  position = (position - 1) % native_length;
  if (position < 0) position += native_length;
  return static_cast<std::size_t>(position + 1);
}

std::uint8_t native_match(
    const std::vector<std::uint8_t>& scores,
    std::ptrdiff_t position) {
  if (scores.empty()) return 0;
  return scores[native_wrap(position, scores.size()) - 1];
}

bool native_missing(
    const std::vector<std::uint8_t>& missing_boundaries,
    std::ptrdiff_t position) {
  if (missing_boundaries.size() <= 1) return false;
  const std::size_t length = missing_boundaries.size() - 1;
  return missing_boundaries[native_wrap(position, length)] != 0;
}

double native_smooth_value(
    const std::vector<double>& smooth,
    std::ptrdiff_t position) {
  if (smooth.size() <= 1) return 0.0;
  const std::size_t length = smooth.size() - 1;
  return smooth[native_wrap(position, length)];
}

struct GrownDiscoveryPeak {
  std::size_t boundary = 1;
  std::size_t pair = 0;
  std::size_t window = 0;
  std::size_t top_left = 0;
  std::size_t top_right = 0;
  double chi = 0.0;
};

GrownDiscoveryPeak grow_discovery_peak(
    const MaxChiWorkspace& profile,
    const std::vector<std::uint8_t>& missing_boundaries,
    const Peak& peak,
    std::size_t initial_half_window) {
  GrownDiscoveryPeak result;
  const std::size_t length = profile.coordinates.size();
  result.boundary = peak.boundary == 0 ? 1 : peak.boundary;
  result.pair = peak.pair;
  result.window = initial_half_window;
  result.chi = peak.chi;
  if (length == 0 || peak.pair >= profile.matches.size()) return result;
  const auto& scores = profile.matches[peak.pair];

  std::size_t test_window = static_cast<std::size_t>(
      static_cast<double>(initial_half_window) / 4.0 + 0.51);
  test_window = std::max<std::size_t>(6, test_window);
  test_window = std::min(test_window, initial_half_window);
  test_window = std::min(test_window, length / 2);
  if (test_window == 0) return result;

  std::size_t left_matches = 0;
  std::size_t right_matches = 0;
  for (std::ptrdiff_t position =
           static_cast<std::ptrdiff_t>(result.boundary + 1) -
               static_cast<std::ptrdiff_t>(test_window);
       position <= static_cast<std::ptrdiff_t>(result.boundary);
       ++position) {
    left_matches += native_match(scores, position);
  }
  for (std::size_t offset = 1; offset <= test_window; ++offset) {
    right_matches += native_match(
        scores,
        static_cast<std::ptrdiff_t>(result.boundary + offset));
  }
  if (test_window >= initial_half_window) {
    result.top_left = left_matches;
    result.top_right = right_matches;
  }

  std::size_t maximum_failures = initial_half_window * 2;
  maximum_failures = std::min(
      maximum_failures,
      (length - test_window * 2) / 2);
  if (maximum_failures == 0) maximum_failures = 1;
  std::size_t failures = 0;
  std::size_t current_window = test_window + 1;
  std::ptrdiff_t left_position =
      static_cast<std::ptrdiff_t>(result.boundary) -
      static_cast<std::ptrdiff_t>(current_window) + 1;
  std::ptrdiff_t right_position =
      static_cast<std::ptrdiff_t>(result.boundary + current_window);

  while (failures <= maximum_failures && current_window * 2 <= length) {
    left_matches += native_match(scores, left_position);
    right_matches += native_match(scores, right_position);
    const double value = chi_square(
        current_window, left_matches, right_matches);
    if (value >= result.chi) {
      result.chi = value;
      result.window = current_window;
      result.top_left = left_matches;
      result.top_right = right_matches;
      failures = 0;
    } else {
      if (current_window == initial_half_window &&
          result.top_left == 0 && result.top_right == 0) {
        result.top_left = left_matches;
        result.top_right = right_matches;
      }
      ++failures;
    }
    // GrowMChiWin2P2 evaluates the enlarged window before consulting MDMap.
    if (native_missing(missing_boundaries, right_position) ||
        native_missing(missing_boundaries, left_position)) {
      break;
    }
    ++right_position;
    --left_position;
    if (native_missing(missing_boundaries, right_position) ||
        native_missing(missing_boundaries, left_position)) {
      break;
    }
    ++current_window;
  }
  if (result.top_left == 0 && result.top_right == 0) {
    result.top_left = circular_sum(
        scores,
        static_cast<std::ptrdiff_t>(result.boundary) -
            static_cast<std::ptrdiff_t>(result.window),
        result.window);
    result.top_right = circular_sum(
        scores,
        static_cast<std::ptrdiff_t>(result.boundary),
        result.window);
  }
  return result;
}

std::array<double, 2> find_side_chi(
    const std::vector<std::uint8_t>& scores,
    const GrownDiscoveryPeak& peak,
    std::size_t left_boundary,
    std::size_t right_boundary) {
  std::size_t outer_left = 0;
  for (std::ptrdiff_t position =
           static_cast<std::ptrdiff_t>(left_boundary) -
               static_cast<std::ptrdiff_t>(peak.window);
       position < static_cast<std::ptrdiff_t>(left_boundary);
       ++position) {
    outer_left += native_match(scores, position);
  }
  std::size_t outer_right = 0;
  for (std::size_t offset = 1; offset <= peak.window; ++offset) {
    outer_right += native_match(
        scores,
        static_cast<std::ptrdiff_t>(right_boundary + offset));
  }
  return {
      chi_square(peak.window, peak.top_left, outer_left),
      chi_square(peak.window, peak.top_right, outer_right),
  };
}

std::size_t optimize_left_breakpoint(
    std::size_t left_boundary,
    double initial_chi,
    std::size_t top_left,
    std::size_t peak_boundary,
    std::size_t window,
    const std::vector<std::uint8_t>& scores,
    const std::vector<std::uint8_t>& missing_boundaries) {
  const std::size_t length = scores.size();
  double outside = 0.0;
  for (std::ptrdiff_t position =
           static_cast<std::ptrdiff_t>(left_boundary) -
               static_cast<std::ptrdiff_t>(window);
       position < static_cast<std::ptrdiff_t>(left_boundary);
       ++position) {
    outside += native_match(scores, position);
  }
  std::ptrdiff_t peak_cursor = static_cast<std::ptrdiff_t>(peak_boundary) - 1;
  std::ptrdiff_t middle = static_cast<std::ptrdiff_t>(left_boundary) - 1;
  std::ptrdiff_t outside_cursor =
      static_cast<std::ptrdiff_t>(left_boundary) -
      static_cast<std::ptrdiff_t>(window);
  peak_cursor = static_cast<std::ptrdiff_t>(native_wrap(peak_cursor, length));
  middle = static_cast<std::ptrdiff_t>(native_wrap(middle, length));
  outside_cursor = static_cast<std::ptrdiff_t>(native_wrap(outside_cursor, length));
  double inside = static_cast<double>(top_left);
  std::int64_t best_scaled = static_cast<std::int64_t>(initial_chi * 10000.0);
  std::size_t best_position = static_cast<std::size_t>(middle);
  std::size_t failures = 0;
  const std::size_t failure_limit = window / 10;
  for (std::size_t step = 0; step < length * 2; ++step) {
    outside -= native_match(scores, middle);
    inside -= native_match(scores, peak_cursor);
    inside += native_match(scores, middle);
    peak_cursor = static_cast<std::ptrdiff_t>(native_wrap(peak_cursor - 1, length));
    outside_cursor = static_cast<std::ptrdiff_t>(native_wrap(outside_cursor - 1, length));
    middle = static_cast<std::ptrdiff_t>(native_wrap(middle - 1, length));
    if (native_missing(missing_boundaries, peak_cursor) ||
        native_missing(missing_boundaries, outside_cursor) ||
        native_missing(missing_boundaries, middle)) {
      break;
    }
    outside += native_match(scores, outside_cursor);
    const double value = chi_square(
        window,
        static_cast<std::size_t>(std::max(0.0, outside)),
        static_cast<std::size_t>(std::max(0.0, inside)));
    const std::int64_t scaled = static_cast<std::int64_t>(value * 10000.0);
    if (scaled > best_scaled) {
      best_scaled = scaled;
      best_position = static_cast<std::size_t>(middle);
      failures = 0;
    } else if (++failures > failure_limit) {
      break;
    }
  }
  return native_wrap(static_cast<std::ptrdiff_t>(best_position) + 1, length);
}

std::size_t optimize_right_breakpoint(
    std::size_t right_boundary,
    double initial_chi,
    std::size_t top_right,
    std::size_t peak_boundary,
    std::size_t window,
    const std::vector<std::uint8_t>& scores,
    const std::vector<std::uint8_t>& missing_boundaries) {
  const std::size_t length = scores.size();
  double outside = 0.0;
  for (std::size_t offset = 1; offset <= window; ++offset) {
    outside += native_match(
        scores,
        static_cast<std::ptrdiff_t>(right_boundary + offset));
  }
  std::ptrdiff_t outside_cursor =
      static_cast<std::ptrdiff_t>(right_boundary + window);
  std::ptrdiff_t middle = static_cast<std::ptrdiff_t>(right_boundary + 1);
  std::ptrdiff_t peak_cursor = static_cast<std::ptrdiff_t>(peak_boundary);
  outside_cursor = static_cast<std::ptrdiff_t>(native_wrap(outside_cursor, length));
  middle = static_cast<std::ptrdiff_t>(native_wrap(middle, length));
  peak_cursor = static_cast<std::ptrdiff_t>(native_wrap(peak_cursor, length));
  double inside = static_cast<double>(top_right);
  std::int64_t best_scaled = static_cast<std::int64_t>(initial_chi * 10000.0);
  std::size_t best_position = static_cast<std::size_t>(middle);
  std::size_t failures = 0;
  const std::size_t failure_limit = window / 10;
  for (std::size_t step = 0; step < length * 2; ++step) {
    inside -= native_match(scores, peak_cursor);
    outside -= native_match(scores, middle);
    inside += native_match(scores, middle);
    outside_cursor = static_cast<std::ptrdiff_t>(native_wrap(outside_cursor + 1, length));
    peak_cursor = static_cast<std::ptrdiff_t>(native_wrap(peak_cursor + 1, length));
    middle = static_cast<std::ptrdiff_t>(native_wrap(middle + 1, length));
    if (native_missing(missing_boundaries, outside_cursor) ||
        native_missing(missing_boundaries, peak_cursor) ||
        native_missing(missing_boundaries, middle)) {
      break;
    }
    outside += native_match(scores, outside_cursor);
    const double value = chi_square(
        window,
        static_cast<std::size_t>(std::max(0.0, inside)),
        static_cast<std::size_t>(std::max(0.0, outside)));
    const std::int64_t scaled = static_cast<std::int64_t>(value * 10000.0);
    if (scaled > best_scaled) {
      best_scaled = scaled;
      best_position = static_cast<std::size_t>(middle);
      failures = 0;
    } else if (++failures > failure_limit) {
      break;
    }
  }
  return native_wrap(static_cast<std::ptrdiff_t>(best_position) - 1, length);
}

void zero_circular_range(
    std::vector<double>& values,
    std::size_t beginning,
    std::size_t ending) {
  if (values.size() <= 1) return;
  const std::size_t length = values.size() - 1;
  beginning = native_wrap(beginning, length);
  ending = native_wrap(ending, length);
  std::size_t position = beginning;
  for (std::size_t count = 0; count < length; ++count) {
    values[position] = 0.0;
    if (position == ending) return;
    position = native_wrap(static_cast<std::ptrdiff_t>(position) + 1, length);
  }
  std::fill(values.begin(), values.end(), 0.0);
}

void destroy_rejected_peak(
    std::size_t peak_boundary,
    const std::vector<double>& smooth,
    std::vector<double>& chi_values) {
  if (chi_values.size() <= 1) return;
  const std::size_t length = chi_values.size() - 1;
  std::size_t right = native_wrap(peak_boundary, length);
  std::size_t left = right;
  bool erase_all = false;
  for (std::size_t count = 0; count < length; ++count) {
    const double current = native_smooth_value(smooth, right);
    if (!(current > 0.0 &&
          (current >= native_smooth_value(smooth, right + 1) ||
           current >= native_smooth_value(smooth, right + 2)))) {
      break;
    }
    right = native_wrap(static_cast<std::ptrdiff_t>(right) + 1, length);
    if (right == left) {
      erase_all = true;
      break;
    }
  }
  if (!erase_all) {
    for (std::size_t count = 0; count < length; ++count) {
      const double current = native_smooth_value(smooth, left);
      if (!(current > 0.0 &&
            (current >= native_smooth_value(smooth, left - 1) ||
             current >= native_smooth_value(smooth, left - 2)))) {
        break;
      }
      left = native_wrap(static_cast<std::ptrdiff_t>(left) - 1, length);
      if (left == right) {
        erase_all = true;
        break;
      }
    }
  }
  if (erase_all) std::fill(chi_values.begin(), chi_values.end(), 0.0);
  else zero_circular_range(chi_values, left, right);
}

void destroy_completed_peak_region(
    std::size_t left_boundary,
    std::size_t right_boundary,
    const std::vector<double>& smooth,
    std::vector<double>& chi_values) {
  if (chi_values.size() <= 1) return;
  const std::size_t length = chi_values.size() - 1;
  std::size_t left = native_wrap(left_boundary, length);
  std::size_t right = native_wrap(right_boundary, length);
  const std::size_t source_left = left;
  bool erase_all = false;

  for (std::size_t count = 0; count < length * 2; ++count) {
    if (!(native_smooth_value(smooth, left) >= native_smooth_value(smooth, left + 1) ||
          native_smooth_value(smooth, left) >= native_smooth_value(smooth, left + 2))) {
      break;
    }
    left = native_wrap(static_cast<std::ptrdiff_t>(left) - 1, length);
    if (left == right) {
      erase_all = true;
      break;
    }
  }
  if (!erase_all) {
    for (std::size_t count = 0; count < length * 2; ++count) {
      if (!((native_smooth_value(smooth, left) <= native_smooth_value(smooth, left + 1) ||
             native_smooth_value(smooth, left) <= native_smooth_value(smooth, left + 2)) &&
            native_smooth_value(smooth, left) > 1.0)) {
        break;
      }
      left = native_wrap(static_cast<std::ptrdiff_t>(left) - 1, length);
      if (left == right) {
        erase_all = true;
        break;
      }
    }
  }
  if (!erase_all) {
    for (std::size_t count = 0; count < length; ++count) {
      right = native_wrap(static_cast<std::ptrdiff_t>(right) + 1, length);
      if (right == left || right == source_left) {
        erase_all = true;
        break;
      }
      if (native_smooth_value(smooth, right) <
          native_smooth_value(smooth, right - 1)) {
        break;
      }
    }
  }
  if (!erase_all) {
    for (std::size_t count = 0; count < length; ++count) {
      right = native_wrap(static_cast<std::ptrdiff_t>(right) + 1, length);
      if (right == left || right == source_left) {
        erase_all = true;
        break;
      }
      if (native_smooth_value(smooth, right) == 0.0 ||
          native_smooth_value(smooth, right) >
              native_smooth_value(smooth, right - 1)) {
        break;
      }
    }
  }
  if (erase_all) std::fill(chi_values.begin(), chi_values.end(), 0.0);
  else zero_circular_range(chi_values, left, right);
}

void include_source_peak_in_destroy_region(
    std::size_t peak_boundary,
    std::size_t& left_boundary,
    std::size_t& right_boundary) {
  // MCXoverF expands whichever optimized interval does not contain pMaxX to
  // the numerically closer endpoint before calling DestroyPeaks. This is a
  // plain source-index comparison, not a shortest circular-distance test.
  const bool contains_peak = left_boundary < right_boundary
      ? peak_boundary >= left_boundary && peak_boundary <= right_boundary
      : peak_boundary > left_boundary || peak_boundary < right_boundary;
  if (contains_peak) return;
  const std::size_t left_distance = peak_boundary >= left_boundary
      ? peak_boundary - left_boundary
      : left_boundary - peak_boundary;
  const std::size_t right_distance = peak_boundary >= right_boundary
      ? peak_boundary - right_boundary
      : right_boundary - peak_boundary;
  if (left_distance > right_distance) right_boundary = peak_boundary;
  else left_boundary = peak_boundary;
}

constexpr std::array<std::array<std::uint8_t, 2>, 3> kPairMembers{{
    {{0, 1}},
    {{0, 2}},
    {{1, 2}},
}};

bool native_position_in_tract(
    std::size_t position,
    std::size_t beginning,
    std::size_t ending) {
  if (beginning < ending) return position >= beginning && position <= ending;
  if (beginning > ending) return position >= beginning || position <= ending;
  return true;
}

void assign_discovery_roles(
    const MaxChiWorkspace& profile,
    std::size_t beginning,
    std::size_t ending,
    std::size_t peak_pair,
    MaxChiDiscoveryCandidate& candidate) {
  const std::size_t length = profile.coordinates.size();
  std::array<double, 3> inside_rate{};
  std::array<double, 3> outside_rate{};
  std::array<double, 3> delta{};
  std::size_t inside_sites = 0;
  for (std::size_t position = 1; position <= length; ++position) {
    if (native_position_in_tract(position, beginning, ending)) ++inside_sites;
  }
  const std::size_t outside_sites = length - inside_sites;
  for (std::size_t pair = 0; pair < profile.matches.size(); ++pair) {
    const auto& matches = profile.matches[pair];
    std::size_t inside_matches = 0;
    std::size_t total_matches = 0;
    for (std::size_t position = 1; position <= length; ++position) {
      total_matches += matches[position - 1];
      if (native_position_in_tract(position, beginning, ending)) {
        inside_matches += matches[position - 1];
      }
    }
    inside_rate[pair] = inside_sites == 0
        ? 0.0
        : static_cast<double>(inside_matches) / static_cast<double>(inside_sites);
    outside_rate[pair] = outside_sites == 0
        ? 0.0
        : static_cast<double>(total_matches - inside_matches) /
            static_cast<double>(outside_sites);
    delta[pair] = inside_rate[pair] - outside_rate[pair];
    candidate.pair_similarity[pair] = length == 0
        ? 0.0
        : static_cast<double>(total_matches) / static_cast<double>(length);
  }

  std::size_t high_pair = 0;
  for (std::size_t pair = 1; pair < delta.size(); ++pair) {
    if (delta[pair] > delta[high_pair]) high_pair = pair;
  }
  std::size_t low_pair = high_pair == 0 ? 1 : 0;
  for (std::size_t pair = 0; pair < delta.size(); ++pair) {
    if (pair == high_pair) continue;
    if (delta[pair] < delta[low_pair]) low_pair = pair;
  }
  if (outside_sites == 0 || std::abs(delta[high_pair] - delta[low_pair]) < 1.0e-12) {
    // A whole-profile/tied tract has no directional contrast. Retain the
    // selected source pair as one role edge and choose the most contrasting
    // remaining edge, leaving the late role consensus to arbitrate.
    high_pair = std::min<std::size_t>(peak_pair, 2);
    low_pair = high_pair == 0 ? 1 : 0;
    for (std::size_t pair = 0; pair < delta.size(); ++pair) {
      if (pair != high_pair && delta[pair] < delta[low_pair]) low_pair = pair;
    }
  }

  std::uint8_t recombinant = 0;
  for (const std::uint8_t member : kPairMembers[high_pair]) {
    if (member == kPairMembers[low_pair][0] ||
        member == kPairMembers[low_pair][1]) {
      recombinant = member;
      break;
    }
  }
  const std::uint8_t minor = kPairMembers[high_pair][0] == recombinant
      ? kPairMembers[high_pair][1]
      : kPairMembers[high_pair][0];
  const std::uint8_t major = kPairMembers[low_pair][0] == recombinant
      ? kPairMembers[low_pair][1]
      : kPairMembers[low_pair][0];
  candidate.recombinant_local = recombinant;
  candidate.major_parent_local = major;
  candidate.minor_parent_local = minor;
  candidate.candidate_pair = static_cast<std::uint8_t>(high_pair);
}

bool choose_discovery_window(
    std::size_t variable_sites,
    std::size_t fixed_window_sites,
    double screening_probability,
    std::size_t& half_window,
    std::size_t& critical_difference) {
  half_window = static_cast<std::size_t>(
      static_cast<double>(fixed_window_sites) / 2.0 + 0.51);
  critical_difference = source_critical_difference(
      half_window, screening_probability);
  if (variable_sites < critical_difference * 2 || variable_sites < 7) return false;
  if (half_window * 2 > variable_sites) {
    half_window = static_cast<std::size_t>(
        static_cast<double>(variable_sites) * 0.75 / 2.0 + 0.51);
    if (half_window > 0) --half_window;
  }
  if (half_window <= critical_difference) {
    half_window = static_cast<std::size_t>(
        static_cast<double>(variable_sites) / 2.0 + 0.51);
    if (half_window > 0) --half_window;
  }
  return half_window >= 6 && half_window * 2 <= variable_sites;
}

struct HeapPeak {
  double chi = 0.0;
  std::size_t boundary = 0;
  std::size_t pair = 0;
};

struct HeapPeakLowerPriority {
  bool operator()(const HeapPeak& left, const HeapPeak& right) const {
    if (left.chi != right.chi) return left.chi < right.chi;
    if (left.boundary != right.boundary) return left.boundary > right.boundary;
    return left.pair > right.pair;
  }
};

}  // namespace

MaxChiDiscoverySummary maxchi_discover(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    const std::vector<std::uint8_t>& triplet_missing_data,
    const MaxChiDiscoveryOptions& options,
    MaxChiWorkspace& workspace,
    std::vector<MaxChiDiscoveryCandidate>& output) {
  if (alignment.length == 0 ||
      std::any_of(triplet.begin(), triplet.end(), [&](std::uint32_t sequence) {
        return sequence >= alignment.sequence_count();
      })) {
    output.clear();
    MaxChiDiscoverySummary summary;
    summary.bonferroni_applied = options.bonferroni;
    summary.fixed_window_sites = options.fixed_window_sites;
    summary.correction_tests = options.bonferroni
        ? std::max<std::uint64_t>(1, options.correction_tests)
        : 1;
    return summary;
  }

  build_variable_profile(alignment, triplet, workspace);
  return maxchi_discover_prepared(
      triplet_missing_data, options, workspace, output);
}

MaxChiDiscoverySummary maxchi_discover_prepared(
    const std::vector<std::uint8_t>& triplet_missing_data,
    const MaxChiDiscoveryOptions& options,
    MaxChiWorkspace& workspace,
    std::vector<MaxChiDiscoveryCandidate>& output) {
  output.clear();
  MaxChiDiscoverySummary summary;
  summary.bonferroni_applied = options.bonferroni;
  summary.fixed_window_sites = options.fixed_window_sites;
  summary.correction_tests = options.bonferroni
      ? std::max<std::uint64_t>(1, options.correction_tests)
      : 1;
  if (workspace.variable_prefix.empty() ||
      std::any_of(
          workspace.matches.begin(),
          workspace.matches.end(),
          [&](const std::vector<std::uint8_t>& matches) {
            return matches.size() != workspace.coordinates.size();
          })) {
    return summary;
  }
  summary.variable_sites = workspace.coordinates.size();
  const double screening_probability = options.p_value_cutoff / 6.0;
  if (!choose_discovery_window(
          summary.variable_sites,
          options.fixed_window_sites,
          screening_probability,
          summary.half_window,
          summary.critical_difference)) {
    return summary;
  }

  bool filter_applied = false;
  make_banned_windows(
      workspace,
      triplet_missing_data,
      summary.half_window,
      options.circular,
      filter_applied);
  summary.missing_data_window_filter_applied = filter_applied;
  summary.linear_edge_window_filter_applied = !options.circular;
  calculate_chi_profiles(
      workspace,
      workspace.banned_windows,
      summary.half_window,
      summary.critical_difference,
      workspace.chi_values);
  summary.profile_available = true;

  std::vector<HeapPeak> peak_storage;
  peak_storage.reserve(summary.variable_sites * workspace.chi_values.size());
  for (std::size_t boundary = 0; boundary < summary.variable_sites; ++boundary) {
    for (std::size_t pair = 0; pair < workspace.chi_values.size(); ++pair) {
      const double value = workspace.chi_values[pair][boundary];
      if (value > 0.0) peak_storage.push_back({value, boundary, pair});
      summary.initial_maximum_chi_square =
          std::max(summary.initial_maximum_chi_square, value);
    }
  }
  if (!(summary.initial_maximum_chi_square > 0.0)) return summary;
  // MCXoverF applies the global raw-maximum gate before SmoothChiValsP and
  // before entering its Redox retry loop. Preserve that ordering: an
  // insignificant triplet performs no smoothing work and records no attempted
  // peak, while later surviving peaks retain the same early-stop gate.
  const double initial_tail = source_chi_p_value(
      summary.initial_maximum_chi_square);
  const double initial_within = std::min(
      1.0,
      initial_tail * static_cast<double>(summary.variable_sites) /
          static_cast<double>(summary.half_window) * 3.0);
  if (initial_tail >= 1.0 || initial_within > options.p_value_cutoff) {
    return summary;
  }
  source_smooth_chi(workspace.chi_values, workspace.smooth_chi);
  // The supplied code repeatedly rescans every surviving raw peak. A single
  // linear-time heapify preserves the same chi/boundary/pair ordering while
  // making subsequent peak extraction logarithmic and destruction lazy.
  std::priority_queue<HeapPeak, std::vector<HeapPeak>, HeapPeakLowerPriority> peaks(
      HeapPeakLowerPriority{}, std::move(peak_storage));

  std::size_t wasted_attempts = 0;
  const std::size_t attempt_limit =
      std::max<std::size_t>(1, options.maximum_peak_attempts);
  while (!peaks.empty() && summary.peak_attempts < attempt_limit) {
    HeapPeak heap_peak = peaks.top();
    peaks.pop();
    if (heap_peak.pair >= workspace.chi_values.size() ||
        heap_peak.boundary >= workspace.chi_values[heap_peak.pair].size() ||
        workspace.chi_values[heap_peak.pair][heap_peak.boundary] <= 0.0 ||
        workspace.chi_values[heap_peak.pair][heap_peak.boundary] != heap_peak.chi) {
      continue;
    }
    ++summary.peak_attempts;
    Peak peak{
        heap_peak.boundary,
        heap_peak.pair,
        summary.half_window,
        heap_peak.chi,
    };
    const double raw_peak_probability = source_chi_p_value(peak.chi);
    const double initial_within = std::min(
        1.0,
        raw_peak_probability * static_cast<double>(summary.variable_sites) /
            static_cast<double>(summary.half_window) * 3.0);
    // FindMChiP returns the global maximum. Once it fails this uncorrected
    // within-triplet gate every remaining peak must fail as well.
    if (raw_peak_probability >= 1.0 || initial_within > options.p_value_cutoff) break;

    GrownDiscoveryPeak grown = grow_discovery_peak(
        workspace, workspace.missing_boundaries, peak, summary.half_window);
    const std::size_t probability_window = std::max<std::size_t>(
        1, std::min(summary.half_window, grown.window));
    const double raw_probability = source_chi_p_value(grown.chi);
    const double within_triplet = std::min(
        1.0,
        raw_probability * static_cast<double>(summary.variable_sites) /
            static_cast<double>(probability_window) * 3.0);
    const double corrected = options.bonferroni
        ? std::min(
            1.0,
            within_triplet * static_cast<double>(summary.correction_tests))
        : within_triplet;
    const std::size_t source_peak = grown.boundary;

    if (corrected < options.p_value_cutoff) {
      const std::size_t provisional_left = native_wrap(
          static_cast<std::ptrdiff_t>(source_peak) -
              static_cast<std::ptrdiff_t>(grown.window),
          summary.variable_sites);
      const std::size_t provisional_right = native_wrap(
          static_cast<std::ptrdiff_t>(source_peak + grown.window - 1),
          summary.variable_sites);
      const auto& scores = workspace.matches[grown.pair];
      const auto side_chi = find_side_chi(
          scores, grown, provisional_left, provisional_right);
      const bool use_left = side_chi[0] >= side_chi[1];
      std::size_t left_boundary = provisional_left;
      std::size_t right_boundary = provisional_right;
      std::size_t tract_beginning = 0;
      std::size_t tract_ending = 0;
      if (use_left) {
        left_boundary = optimize_left_breakpoint(
            provisional_left,
            side_chi[0],
            grown.top_left,
            source_peak,
            grown.window,
            scores,
            workspace.missing_boundaries);
        std::ptrdiff_t native_end = static_cast<std::ptrdiff_t>(source_peak) - 1;
        native_end = static_cast<std::ptrdiff_t>(
            native_wrap(native_end, summary.variable_sites));
        for (std::size_t count = 0;
             count < summary.variable_sites &&
             native_missing(workspace.missing_boundaries, native_end);
             ++count) {
          native_end = static_cast<std::ptrdiff_t>(
              native_wrap(native_end - 1, summary.variable_sites));
        }
        native_end = static_cast<std::ptrdiff_t>(
            native_wrap(native_end + 1, summary.variable_sites));
        tract_beginning = left_boundary + 1 > summary.variable_sites
            ? left_boundary
            : left_boundary + 1;
        tract_ending = static_cast<std::size_t>(native_end);
        right_boundary = tract_ending;
      } else {
        right_boundary = optimize_right_breakpoint(
            provisional_right,
            side_chi[1],
            grown.top_right,
            source_peak,
            grown.window,
            scores,
            workspace.missing_boundaries);
        std::ptrdiff_t native_begin = static_cast<std::ptrdiff_t>(source_peak);
        if (native_begin < static_cast<std::ptrdiff_t>(summary.variable_sites)) {
          ++native_begin;
        }
        native_begin = static_cast<std::ptrdiff_t>(
            native_wrap(native_begin, summary.variable_sites));
        for (std::size_t count = 0;
             count < summary.variable_sites &&
             native_missing(workspace.missing_boundaries, native_begin);
             ++count) {
          native_begin = static_cast<std::ptrdiff_t>(
              native_wrap(native_begin + 1, summary.variable_sites));
        }
        native_begin = static_cast<std::ptrdiff_t>(
            native_wrap(native_begin - 1, summary.variable_sites));
        tract_beginning = static_cast<std::size_t>(native_begin) + 1 >
                summary.variable_sites
            ? static_cast<std::size_t>(native_begin)
            : static_cast<std::size_t>(native_begin) + 1;
        tract_ending = right_boundary;
        left_boundary = native_wrap(source_peak, summary.variable_sites);
      }

      bool emitted = false;
      if (tract_beginning >= 1 && tract_beginning <= summary.variable_sites &&
          tract_ending >= 1 && tract_ending <= summary.variable_sites &&
          (options.circular || tract_beginning < tract_ending)) {
        MaxChiDiscoveryCandidate candidate;
        candidate.beginning = workspace.coordinates[tract_beginning - 1];
        candidate.ending = workspace.coordinates[tract_ending - 1];
        candidate.wraps_origin = options.circular &&
            candidate.beginning >= candidate.ending;
        candidate.informative_beginning = tract_beginning;
        candidate.informative_ending = tract_ending;
        candidate.peak_pair = static_cast<std::int8_t>(grown.pair);
        candidate.tract_side = use_left
            ? MaxChiTractSide::left
            : MaxChiTractSide::right;
        candidate.peak_attempt = summary.peak_attempts;
        const std::size_t peak_coordinate_index = peak.boundary == 0
            ? summary.variable_sites - 1
            : peak.boundary - 1;
        candidate.peak_alignment_position =
            workspace.coordinates[peak_coordinate_index];
        candidate.variable_sites = summary.variable_sites;
        candidate.initial_half_window = summary.half_window;
        candidate.grown_half_window = grown.window;
        candidate.critical_difference = summary.critical_difference;
        candidate.maximum_chi_square = grown.chi;
        candidate.raw_p_value = raw_probability;
        candidate.within_triplet_p_value = within_triplet;
        candidate.corrected_p_value = corrected;
        candidate.left_flank_chi_square = side_chi[0];
        candidate.right_flank_chi_square = side_chi[1];
        candidate.missing_data_window_filter_applied =
            summary.missing_data_window_filter_applied;
        candidate.linear_edge_window_filter_applied =
            summary.linear_edge_window_filter_applied;
        assign_discovery_roles(
            workspace,
            tract_beginning,
            tract_ending,
            grown.pair,
            candidate);
        output.push_back(candidate);
        emitted = true;
        // WasteOfTime is reset only after MCXoverF stores a usable hit.
        wasted_attempts = 0;
      }
      include_source_peak_in_destroy_region(
          source_peak, left_boundary, right_boundary);
      destroy_completed_peak_region(
          left_boundary,
          right_boundary,
          workspace.smooth_chi[grown.pair],
          workspace.chi_values[grown.pair]);
      ++summary.destroyed_peak_regions;
      if (!emitted && ++wasted_attempts >= 3) break;
    } else {
      ++wasted_attempts;
      destroy_rejected_peak(
          source_peak,
          workspace.smooth_chi[grown.pair],
          workspace.chi_values[grown.pair]);
      ++summary.destroyed_peak_regions;
      if (wasted_attempts >= 3) break;
    }
    workspace.chi_values[grown.pair][heap_peak.boundary] = 0.0;
  }
  bool surviving_peak = false;
  if (summary.peak_attempts >= attempt_limit) {
    for (std::size_t pair = 0;
         pair < workspace.chi_values.size() && !surviving_peak;
         ++pair) {
      surviving_peak = std::any_of(
          workspace.chi_values[pair].begin(),
          workspace.chi_values[pair].begin() +
              static_cast<std::ptrdiff_t>(summary.variable_sites),
          [](double value) { return value > 0.0; });
    }
  }
  summary.peak_attempt_limit_reached = surviving_peak;
  summary.emitted_candidates = output.size();
  return summary;
}

MaxChiPlotProfile maxchi_plot_profile(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    const std::vector<std::uint8_t>& triplet_missing_data,
    bool circular,
    std::size_t fixed_window_sites,
    double p_value_cutoff,
    MaxChiWorkspace& workspace) {
  MaxChiPlotProfile plot;
  if (alignment.length == 0 ||
      std::any_of(triplet.begin(), triplet.end(), [&](std::uint32_t sequence) {
        return sequence >= alignment.sequence_count();
      })) {
    return plot;
  }
  build_variable_profile(alignment, triplet, workspace);
  std::size_t critical_difference = 0;
  if (!choose_discovery_window(
          workspace.coordinates.size(),
          fixed_window_sites,
          std::clamp(p_value_cutoff / 6.0, kMinimumProbability, 1.0),
          plot.half_window,
          critical_difference)) {
    return plot;
  }
  bool ignored_filter = false;
  make_banned_windows(
      workspace,
      triplet_missing_data,
      plot.half_window,
      circular,
      ignored_filter);
  calculate_chi_profiles(
      workspace,
      workspace.banned_windows,
      plot.half_window,
      critical_difference,
      workspace.chi_values);
  plot.coordinates.reserve(workspace.coordinates.size());
  for (auto& values : plot.chi_square) {
    values.reserve(workspace.coordinates.size());
  }
  // Boundary 0 maps to the final variable coordinate. Rotate it to the end so
  // the browser receives monotonically increasing alignment positions.
  for (std::size_t output_index = 0;
       output_index < workspace.coordinates.size();
       ++output_index) {
    const std::size_t boundary =
        (output_index + 1) % workspace.coordinates.size();
    const std::size_t coordinate_index = boundary == 0
        ? workspace.coordinates.size() - 1
        : boundary - 1;
    plot.coordinates.push_back(workspace.coordinates[coordinate_index]);
    for (std::size_t pair = 0; pair < plot.chi_square.size(); ++pair) {
      plot.chi_square[pair].push_back(workspace.chi_values[pair][boundary]);
    }
  }
  plot.available = true;
  return plot;
}

MaxChiRecheckEvidence maxchi_recheck(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    const std::vector<std::uint8_t>& triplet_missing_data,
    const MaxChiRecheckOptions& options,
    MaxChiWorkspace& workspace) {
  MaxChiRecheckEvidence evidence;
  evidence.requested = true;
  evidence.fixed_window_sites = options.fixed_window_sites;
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

namespace {

constexpr std::array<std::uint8_t, 3> kChimaeraParentOne{{1, 2, 0}};
constexpr std::array<std::uint8_t, 3> kChimaeraParentTwo{{2, 0, 1}};
constexpr std::array<std::uint8_t, 3> kChimaeraScorePair{{0, 2, 1}};
constexpr std::array<std::uint8_t, 3> kChimaeraOtherPair{{1, 0, 2}};

std::uint8_t pair_index_for_members(std::uint8_t first, std::uint8_t second) {
  if (first > second) std::swap(first, second);
  if (first == 0 && second == 1) return 0;
  if (first == 0 && second == 2) return 1;
  return 2;
}

bool build_chimaera_target_profile(
    const MaxChiWorkspace& variable_workspace,
    std::uint8_t target_local,
    MaxChiWorkspace& target_workspace) {
  target_workspace.coordinates.clear();
  for (auto& matches : target_workspace.matches) matches.clear();
  target_workspace.variable_prefix.assign(
      variable_workspace.variable_prefix.size(), 0);
  if (target_local > 2 || variable_workspace.variable_prefix.empty()) return false;
  for (const auto& matches : variable_workspace.matches) {
    if (matches.size() != variable_workspace.coordinates.size()) return false;
  }

  const std::size_t score_pair = kChimaeraScorePair[target_local];
  const std::size_t other_pair = kChimaeraOtherPair[target_local];
  target_workspace.coordinates.reserve(variable_workspace.coordinates.size());
  target_workspace.matches[0].reserve(variable_workspace.coordinates.size());
  for (std::size_t index = 0;
       index < variable_workspace.coordinates.size();
       ++index) {
    const bool target_matches_parent_one =
        variable_workspace.matches[score_pair][index] != 0;
    const bool target_matches_parent_two =
        variable_workspace.matches[other_pair][index] != 0;
    if (!target_matches_parent_one && !target_matches_parent_two) continue;
    target_workspace.coordinates.push_back(variable_workspace.coordinates[index]);
    target_workspace.matches[0].push_back(target_matches_parent_one ? 1 : 0);
  }
  target_workspace.matches[1].assign(target_workspace.coordinates.size(), 0);
  target_workspace.matches[2].assign(target_workspace.coordinates.size(), 0);

  std::size_t included = 0;
  for (std::size_t coordinate = 1;
       coordinate < target_workspace.variable_prefix.size();
       ++coordinate) {
    while (included < target_workspace.coordinates.size() &&
           target_workspace.coordinates[included] <= coordinate) {
      ++included;
    }
    target_workspace.variable_prefix[coordinate] = included;
  }
  return true;
}

struct ChimaeraTargetSummary {
  bool profile_available = false;
  bool peak_limit_reached = false;
  std::size_t peak_attempts = 0;
  std::size_t destroyed_peak_regions = 0;
};

ChimaeraTargetSummary discover_chimaera_target(
    std::uint8_t target_local,
    const std::vector<std::uint8_t>& triplet_missing_data,
    const std::array<double, 3>& pair_similarity,
    const ChimaeraDiscoveryOptions& options,
    MaxChiWorkspace& workspace,
    std::vector<ChimaeraDiscoveryCandidate>& output) {
  ChimaeraTargetSummary summary;
  const std::size_t information_rich_sites = workspace.coordinates.size();
  std::size_t half_window = 0;
  std::size_t critical_difference = 0;
  if (!choose_discovery_window(
          information_rich_sites,
          options.fixed_window_sites,
          options.p_value_cutoff / 6.0,
          half_window,
          critical_difference)) {
    return summary;
  }

  bool missing_filter_applied = false;
  make_banned_windows(
      workspace,
      triplet_missing_data,
      half_window,
      options.circular,
      missing_filter_applied);
  calculate_chi_profiles(
      workspace,
      workspace.banned_windows,
      half_window,
      critical_difference,
      workspace.chi_values);
  summary.profile_available = true;

  std::vector<HeapPeak> peak_storage;
  peak_storage.reserve(information_rich_sites);
  double initial_maximum = 0.0;
  for (std::size_t boundary = 0; boundary < information_rich_sites; ++boundary) {
    const double value = workspace.chi_values[0][boundary];
    if (value > 0.0) peak_storage.push_back({value, boundary, 0});
    initial_maximum = std::max(initial_maximum, value);
  }
  if (!(initial_maximum > 0.0)) return summary;
  const double initial_tail = source_chi_p_value(initial_maximum);
  const double initial_within = std::min(
      1.0,
      initial_tail * static_cast<double>(information_rich_sites) /
          static_cast<double>(half_window) * 3.0);
  if (initial_tail >= 1.0 || initial_within > options.p_value_cutoff) {
    return summary;
  }
  source_smooth_chi(workspace.chi_values, workspace.smooth_chi);
  std::priority_queue<HeapPeak, std::vector<HeapPeak>, HeapPeakLowerPriority> peaks(
      HeapPeakLowerPriority{}, std::move(peak_storage));

  std::size_t wasted_attempts = 0;
  const std::size_t attempt_limit =
      std::max<std::size_t>(1, options.maximum_peak_attempts);
  while (!peaks.empty() && summary.peak_attempts < attempt_limit) {
    const HeapPeak heap_peak = peaks.top();
    peaks.pop();
    if (heap_peak.boundary >= workspace.chi_values[0].size() ||
        workspace.chi_values[0][heap_peak.boundary] <= 0.0 ||
        workspace.chi_values[0][heap_peak.boundary] != heap_peak.chi) {
      continue;
    }
    ++summary.peak_attempts;
    Peak peak{heap_peak.boundary, 0, half_window, heap_peak.chi};
    const double peak_tail = source_chi_p_value(peak.chi);
    const double peak_within = std::min(
        1.0,
        peak_tail * static_cast<double>(information_rich_sites) /
            static_cast<double>(half_window) * 3.0);
    if (peak_tail >= 1.0 || peak_within > options.p_value_cutoff) break;

    GrownDiscoveryPeak grown = grow_discovery_peak(
        workspace, workspace.missing_boundaries, peak, half_window);
    const std::size_t probability_window = std::max<std::size_t>(
        1, std::min(half_window, grown.window));
    const double raw_probability = source_chi_p_value(grown.chi);
    const double within_triplet = std::min(
        1.0,
        raw_probability * static_cast<double>(information_rich_sites) /
            static_cast<double>(probability_window) * 3.0);
    const double corrected = options.bonferroni
        ? std::min(
            1.0,
            within_triplet * static_cast<double>(
                std::max<std::uint64_t>(1, options.correction_tests)))
        : within_triplet;
    const std::size_t source_peak = grown.boundary;

    if (corrected < options.p_value_cutoff) {
      const std::size_t provisional_left = native_wrap(
          static_cast<std::ptrdiff_t>(source_peak) -
              static_cast<std::ptrdiff_t>(grown.window),
          information_rich_sites);
      const std::size_t provisional_right = native_wrap(
          static_cast<std::ptrdiff_t>(source_peak + grown.window - 1),
          information_rich_sites);
      const auto& scores = workspace.matches[0];
      const auto side_chi = find_side_chi(
          scores, grown, provisional_left, provisional_right);
      const bool use_left = side_chi[0] >= side_chi[1];
      std::size_t left_boundary = provisional_left;
      std::size_t right_boundary = provisional_right;
      std::size_t tract_beginning = 0;
      std::size_t tract_ending = 0;
      if (use_left) {
        left_boundary = optimize_left_breakpoint(
            provisional_left,
            side_chi[0],
            grown.top_left,
            source_peak,
            grown.window,
            scores,
            workspace.missing_boundaries);
        std::ptrdiff_t native_end = static_cast<std::ptrdiff_t>(source_peak) - 1;
        native_end = static_cast<std::ptrdiff_t>(
            native_wrap(native_end, information_rich_sites));
        for (std::size_t count = 0;
             count < information_rich_sites &&
             native_missing(workspace.missing_boundaries, native_end);
             ++count) {
          native_end = static_cast<std::ptrdiff_t>(
              native_wrap(native_end - 1, information_rich_sites));
        }
        native_end = static_cast<std::ptrdiff_t>(
            native_wrap(native_end + 1, information_rich_sites));
        tract_beginning = left_boundary + 1 > information_rich_sites
            ? left_boundary
            : left_boundary + 1;
        tract_ending = static_cast<std::size_t>(native_end);
        right_boundary = tract_ending;
      } else {
        right_boundary = optimize_right_breakpoint(
            provisional_right,
            side_chi[1],
            grown.top_right,
            source_peak,
            grown.window,
            scores,
            workspace.missing_boundaries);
        std::ptrdiff_t native_begin = static_cast<std::ptrdiff_t>(source_peak);
        if (native_begin < static_cast<std::ptrdiff_t>(information_rich_sites)) {
          ++native_begin;
        }
        native_begin = static_cast<std::ptrdiff_t>(
            native_wrap(native_begin, information_rich_sites));
        for (std::size_t count = 0;
             count < information_rich_sites &&
             native_missing(workspace.missing_boundaries, native_begin);
             ++count) {
          native_begin = static_cast<std::ptrdiff_t>(
              native_wrap(native_begin + 1, information_rich_sites));
        }
        native_begin = static_cast<std::ptrdiff_t>(
            native_wrap(native_begin - 1, information_rich_sites));
        tract_beginning = static_cast<std::size_t>(native_begin) + 1 >
                information_rich_sites
            ? static_cast<std::size_t>(native_begin)
            : static_cast<std::size_t>(native_begin) + 1;
        tract_ending = right_boundary;
        left_boundary = native_wrap(source_peak, information_rich_sites);
      }

      bool emitted = false;
      if (tract_beginning >= 1 && tract_beginning <= information_rich_sites &&
          tract_ending >= 1 && tract_ending <= information_rich_sites &&
          (options.circular || tract_beginning < tract_ending)) {
        std::size_t inside_sites = 0;
        std::size_t inside_parent_one_matches = 0;
        std::size_t total_parent_one_matches = 0;
        for (std::size_t position = 1;
             position <= information_rich_sites;
             ++position) {
          total_parent_one_matches += scores[position - 1];
          if (native_position_in_tract(position, tract_beginning, tract_ending)) {
            ++inside_sites;
            inside_parent_one_matches += scores[position - 1];
          }
        }
        const std::size_t outside_sites = information_rich_sites - inside_sites;
        ChimaeraDiscoveryCandidate candidate;
        candidate.beginning = workspace.coordinates[tract_beginning - 1];
        candidate.ending = workspace.coordinates[tract_ending - 1];
        candidate.wraps_origin = options.circular &&
            candidate.beginning >= candidate.ending;
        candidate.informative_beginning = tract_beginning;
        candidate.informative_ending = tract_ending;
        candidate.target_local = target_local;
        candidate.recombinant_local = target_local;
        candidate.inside_parent_one_match_rate = inside_sites == 0
            ? 0.0
            : static_cast<double>(inside_parent_one_matches) /
                static_cast<double>(inside_sites);
        candidate.outside_parent_one_match_rate = outside_sites == 0
            ? 0.0
            : static_cast<double>(
                total_parent_one_matches - inside_parent_one_matches) /
                static_cast<double>(outside_sites);
        const bool parent_one_is_minor =
            candidate.inside_parent_one_match_rate >=
            candidate.outside_parent_one_match_rate;
        candidate.minor_parent_local = parent_one_is_minor
            ? kChimaeraParentOne[target_local]
            : kChimaeraParentTwo[target_local];
        candidate.major_parent_local = parent_one_is_minor
            ? kChimaeraParentTwo[target_local]
            : kChimaeraParentOne[target_local];
        candidate.candidate_pair = pair_index_for_members(
            candidate.recombinant_local, candidate.minor_parent_local);
        candidate.tract_side = use_left
            ? MaxChiTractSide::left
            : MaxChiTractSide::right;
        candidate.peak_attempt = summary.peak_attempts;
        const std::size_t peak_coordinate_index = peak.boundary == 0
            ? information_rich_sites - 1
            : peak.boundary - 1;
        candidate.peak_alignment_position =
            workspace.coordinates[peak_coordinate_index];
        candidate.information_rich_sites = information_rich_sites;
        candidate.initial_half_window = half_window;
        candidate.grown_half_window = grown.window;
        candidate.critical_difference = critical_difference;
        candidate.maximum_chi_square = grown.chi;
        candidate.raw_p_value = raw_probability;
        candidate.within_triplet_p_value = within_triplet;
        candidate.corrected_p_value = corrected;
        candidate.left_flank_chi_square = side_chi[0];
        candidate.right_flank_chi_square = side_chi[1];
        candidate.pair_similarity = pair_similarity;
        candidate.missing_data_window_filter_applied = missing_filter_applied;
        candidate.linear_edge_window_filter_applied = !options.circular;
        output.push_back(candidate);
        emitted = true;
        wasted_attempts = 0;
      }
      include_source_peak_in_destroy_region(
          source_peak, left_boundary, right_boundary);
      destroy_completed_peak_region(
          left_boundary,
          right_boundary,
          workspace.smooth_chi[0],
          workspace.chi_values[0]);
      ++summary.destroyed_peak_regions;
      if (!emitted && ++wasted_attempts >= 3) break;
    } else {
      ++wasted_attempts;
      destroy_rejected_peak(
          source_peak, workspace.smooth_chi[0], workspace.chi_values[0]);
      ++summary.destroyed_peak_regions;
      if (wasted_attempts >= 3) break;
    }
    workspace.chi_values[0][heap_peak.boundary] = 0.0;
  }
  if (summary.peak_attempts >= attempt_limit) {
    summary.peak_limit_reached = std::any_of(
        workspace.chi_values[0].begin(),
        workspace.chi_values[0].begin() +
            static_cast<std::ptrdiff_t>(information_rich_sites),
        [](double value) { return value > 0.0; });
  }
  return summary;
}

}  // namespace

ChimaeraDiscoverySummary chimaera_discover_prepared(
    const MaxChiWorkspace& variable_workspace,
    const std::vector<std::uint8_t>& triplet_missing_data,
    const std::array<double, 3>& pair_similarity,
    const ChimaeraDiscoveryOptions& options,
    MaxChiWorkspace& target_workspace,
    std::vector<ChimaeraDiscoveryCandidate>& output) {
  output.clear();
  ChimaeraDiscoverySummary summary;
  summary.bonferroni_applied = options.bonferroni;
  summary.correction_tests = options.bonferroni
      ? std::max<std::uint64_t>(1, options.correction_tests)
      : 1;
  summary.fixed_window_sites = options.fixed_window_sites;
  for (std::uint8_t target_local = 0; target_local < 3; ++target_local) {
    if (!build_chimaera_target_profile(
            variable_workspace, target_local, target_workspace)) {
      continue;
    }
    const ChimaeraTargetSummary target_summary = discover_chimaera_target(
        target_local,
        triplet_missing_data,
        pair_similarity,
        options,
        target_workspace,
        output);
    if (target_summary.profile_available) ++summary.target_profiles_scanned;
    summary.peak_attempts += target_summary.peak_attempts;
    summary.destroyed_peak_regions += target_summary.destroyed_peak_regions;
    if (target_summary.peak_limit_reached) ++summary.peak_limit_targets;
  }
  summary.emitted_candidates = output.size();
  return summary;
}

ChimaeraRecheckEvidence chimaera_recheck_prepared(
    const MaxChiWorkspace& variable_workspace,
    const std::vector<std::uint8_t>& triplet_missing_data,
    const ChimaeraRecheckOptions& options,
    MaxChiWorkspace& target_workspace) {
  ChimaeraRecheckEvidence evidence;
  evidence.requested = true;
  evidence.fixed_window_sites = options.fixed_window_sites;
  evidence.bonferroni_applied = options.bonferroni;
  evidence.correction_tests = options.bonferroni
      ? std::max<std::uint64_t>(1, options.correction_tests)
      : 1;
  bool selected = false;
  for (std::uint8_t target_local = 0; target_local < 3; ++target_local) {
    if (!build_chimaera_target_profile(
            variable_workspace, target_local, target_workspace)) {
      continue;
    }
    const std::size_t information_rich_sites =
        target_workspace.coordinates.size();
    std::size_t half_window = 0;
    std::size_t critical_difference = 0;
    if (!choose_discovery_window(
            information_rich_sites,
            options.fixed_window_sites,
            options.p_value_cutoff / 6.0,
            half_window,
            critical_difference)) {
      continue;
    }
    ++evidence.target_profiles_scanned;
    evidence.profile_available = true;
    bool missing_filter_applied = false;
    make_banned_windows(
        target_workspace,
        triplet_missing_data,
        half_window,
        options.circular,
        missing_filter_applied);
    evidence.missing_data_window_filter_applied =
        evidence.missing_data_window_filter_applied || missing_filter_applied;
    evidence.linear_edge_window_filter_applied = !options.circular;
    Peak peak = strongest_peak(
        target_workspace,
        target_workspace.banned_windows,
        half_window,
        critical_difference);
    if (!(peak.chi > 0.0)) continue;

    const double initial_tail = source_chi_p_value(peak.chi);
    const double initial_within = std::min(
        1.0,
        initial_tail * static_cast<double>(information_rich_sites) /
            static_cast<double>(half_window) * 3.0);
    if (initial_within <= options.p_value_cutoff && initial_tail < 1.0) {
      grow_peak(target_workspace, target_workspace.missing_boundaries, peak);
    }
    const double local_p_value = source_chi_p_value(peak.chi);
    const std::size_t probability_window =
        std::max<std::size_t>(1, std::min(half_window, peak.half_window));
    const double within_triplet_p_value = std::min(
        1.0,
        local_p_value * static_cast<double>(information_rich_sites) /
            static_cast<double>(probability_window) * 3.0);
    const double corrected_p_value = options.bonferroni
        ? std::min(
            1.0,
            within_triplet_p_value *
                static_cast<double>(evidence.correction_tests))
        : within_triplet_p_value;
    if (selected &&
        std::make_tuple(
            corrected_p_value,
            within_triplet_p_value,
            local_p_value,
            target_local) >=
            std::make_tuple(
                evidence.corrected_p_value,
                evidence.within_triplet_p_value,
                evidence.local_p_value,
                static_cast<std::uint8_t>(evidence.best_target))) {
      continue;
    }

    selected = true;
    evidence.best_target = static_cast<std::int8_t>(target_local);
    evidence.information_rich_sites = information_rich_sites;
    evidence.half_window = half_window;
    evidence.critical_difference = critical_difference;
    evidence.grown_half_window = peak.half_window;
    const std::size_t peak_coordinate_index = peak.boundary == 0
        ? information_rich_sites - 1
        : peak.boundary - 1;
    evidence.peak_alignment_position =
        target_workspace.coordinates[peak_coordinate_index];
    evidence.maximum_chi_square = peak.chi;
    evidence.local_p_value = local_p_value;
    evidence.within_triplet_p_value = within_triplet_p_value;
    evidence.corrected_p_value = corrected_p_value;
  }
  evidence.source_recheck_hit = selected &&
      evidence.corrected_p_value < options.p_value_cutoff;
  return evidence;
}

ChimaeraPlotProfile chimaera_plot_profile(
    const Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet,
    std::uint8_t target_local,
    const std::vector<std::uint8_t>& triplet_missing_data,
    bool circular,
    std::size_t fixed_window_sites,
    double p_value_cutoff,
    MaxChiWorkspace& variable_workspace,
    MaxChiWorkspace& target_workspace) {
  ChimaeraPlotProfile plot;
  plot.target_local = target_local;
  if (alignment.length == 0 || target_local > 2 ||
      std::any_of(triplet.begin(), triplet.end(), [&](std::uint32_t sequence) {
        return sequence >= alignment.sequence_count();
      })) {
    return plot;
  }
  build_variable_profile(alignment, triplet, variable_workspace);
  if (!build_chimaera_target_profile(
          variable_workspace, target_local, target_workspace)) {
    return plot;
  }
  std::size_t critical_difference = 0;
  if (!choose_discovery_window(
          target_workspace.coordinates.size(),
          fixed_window_sites,
          std::clamp(p_value_cutoff / 6.0, kMinimumProbability, 1.0),
          plot.half_window,
          critical_difference)) {
    return plot;
  }
  bool ignored_filter = false;
  make_banned_windows(
      target_workspace,
      triplet_missing_data,
      plot.half_window,
      circular,
      ignored_filter);
  calculate_chi_profiles(
      target_workspace,
      target_workspace.banned_windows,
      plot.half_window,
      critical_difference,
      target_workspace.chi_values);
  plot.coordinates.reserve(target_workspace.coordinates.size());
  plot.chi_square.reserve(target_workspace.coordinates.size());
  for (std::size_t output_index = 0;
       output_index < target_workspace.coordinates.size();
       ++output_index) {
    const std::size_t boundary =
        (output_index + 1) % target_workspace.coordinates.size();
    const std::size_t coordinate_index = boundary == 0
        ? target_workspace.coordinates.size() - 1
        : boundary - 1;
    plot.coordinates.push_back(target_workspace.coordinates[coordinate_index]);
    plot.chi_square.push_back(target_workspace.chi_values[0][boundary]);
  }
  plot.available = true;
  return plot;
}

}  // namespace rdp
