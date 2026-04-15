#include <cmath>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "ipinstance.h"
#include "timer.h"

namespace {

std::string EscapeJson(const std::string& value) {
  std::string escaped;
  escaped.reserve(value.size());
  for (char c : value) {
    switch (c) {
      case '\\':
        escaped += "\\\\";
        break;
      case '"':
        escaped += "\\\"";
        break;
      default:
        escaped += c;
        break;
    }
  }
  return escaped;
}

std::string FormatTime(double seconds) {
  std::ostringstream out;
  out << std::fixed << std::setprecision(2) << seconds;
  return out.str();
}

std::string FormatObjective(double value) {
  std::ostringstream out;
  if (std::abs(value - std::round(value)) <= 1e-6) {
    out << static_cast<long long>(std::llround(value));
  } else {
    out << std::fixed << std::setprecision(6) << value;
  }
  return out.str();
}

}  // namespace

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <input>\n";
    return 65;
  }

  const std::string input_file = argv[1];
  const std::string filename = std::filesystem::path(input_file).filename().string();

  Timer timer;
  timer.Start();

  SolveResult result;
  try {
    IPInstance instance(input_file);
    result = instance.Solve();
  } catch (const std::exception& ex) {
    timer.Stop();
    std::cerr << ex.what() << '\n';
    std::cout << "{\"Instance\":\"" << EscapeJson(filename)
              << "\",\"Time\":\"" << FormatTime(timer.GetTime())
              << "\",\"Result\":\"--\",\"Solution\":\"--\"}\n";
    return 1;
  }

  timer.Stop();

  std::cout << "{\"Instance\":\"" << EscapeJson(filename)
            << "\",\"Time\":\"" << FormatTime(timer.GetTime()) << "\",\"Result\":";
  if (result.has_solution) {
    std::cout << FormatObjective(result.objective_value);
  } else {
    std::cout << "\"--\"";
  }
  std::cout << ",\"Solution\":\"" << (result.is_optimal ? "OPT" : "--") << "\"}\n";

  return 0;
}
