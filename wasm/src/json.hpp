#pragma once

#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <string_view>

namespace rdp::json {

inline std::string escape(std::string_view value) {
  std::ostringstream out;
  for (const unsigned char character : value) {
    switch (character) {
      case '"': out << "\\\""; break;
      case '\\': out << "\\\\"; break;
      case '\b': out << "\\b"; break;
      case '\f': out << "\\f"; break;
      case '\n': out << "\\n"; break;
      case '\r': out << "\\r"; break;
      case '\t': out << "\\t"; break;
      default:
        if (character < 0x20) {
          out << "\\u" << std::hex << std::setw(4) << std::setfill('0')
              << static_cast<int>(character) << std::dec;
        } else {
          out << static_cast<char>(character);
        }
    }
  }
  return out.str();
}

inline void string(std::ostringstream& out, std::string_view value) {
  out << '"' << escape(value) << '"';
}

inline void number(std::ostringstream& out, double value) {
  if (!std::isfinite(value)) {
    out << "null";
    return;
  }
  out << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
}

}  // namespace rdp::json
