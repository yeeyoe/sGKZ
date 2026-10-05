#pragma once

#include <atomic>
#include <array>
#include <cstddef>
#include <functional>
#include <string>
#include <utility>
#include <vector>

namespace gkz_gui {

struct GuiPoint {
  long long x = 0;
  long long y = 0;
};

struct GuiSolveResult {
  bool certified = false;
  bool cancelled = false;
  std::string message;
  std::vector<std::string> sigma_exact;
  std::vector<std::string> sigma_vee_exact;
  std::vector<std::array<std::size_t, 3>> triangulation_faces;
  std::vector<std::vector<std::pair<double, double>>> subdivision_cells;
  int iterations = 0;
  std::size_t active_size = 0;
};

using ProgressCallback = std::function<void(int, std::size_t, const std::string&)>;

GuiSolveResult solve_points(const std::vector<GuiPoint>& points,
                            std::atomic_bool* cancel = nullptr,
                            ProgressCallback progress = {});

}  // namespace gkz_gui
