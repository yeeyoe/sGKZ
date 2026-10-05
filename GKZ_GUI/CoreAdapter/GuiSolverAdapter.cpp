#include "GuiSolverAdapter.hpp"

#include "gkz/gkz.hpp"

#include <stdexcept>
#include <utility>

namespace gkz_gui {
namespace {

std::string rational_string(const CGAL::Gmpq& value) {
  std::ostringstream out;
  if (value.denominator() == 1) out << value.numerator();
  else out << value;
  return out.str();
}

}  // namespace

GuiSolveResult solve_points(const std::vector<GuiPoint>& points,
                            std::atomic_bool* cancel,
                            ProgressCallback progress) {
  GuiSolveResult out;
  try {
    if (cancel && cancel->load()) { out.cancelled = true; out.message = "cancelled"; return out; }
    std::vector<gkz::IntPoint> input;
    input.reserve(points.size());
    for (const auto& point : points) input.push_back({point.x, point.y});
    const gkz::PointConfiguration configuration = gkz::PointConfiguration::from_points(std::move(input));
    gkz::SolverOptions options;
    options.exact_certification = true;
    options.verbose = false;
    options.cancellation = cancel;
    options.progress_callback = [progress](int iteration, std::size_t active) {
      if (progress) progress(iteration, active, "Solving");
    };
    const gkz::ShortestGkzSolver solver(options);
    const auto result = solver.solve(configuration);
    out.iterations = result.iterations;
    out.active_size = result.active_vectors.size();
    if (result.cancelled || (cancel && cancel->load())) {
      out.cancelled = true;
      out.message = "Cancelled";
      return out;
    }
    if (cancel && cancel->load()) { out.cancelled = true; out.message = "cancelled"; return out; }
    if (!result.exact.certified) {
      out.message = result.exact.message.empty() ? "Exact certification failed" : result.exact.message;
      return out;
    }
    out.certified = true;
    out.message = "Exact verification passed";
    for (const auto& value : result.exact.sigma) out.sigma_exact.push_back(rational_string(value));
    if (progress) progress(result.iterations, out.active_size, "Building subdivision");
    const auto plot = gkz::compute_certified_plot_data(configuration, result);
    for (const auto& value : plot.sigma_vee)
      out.sigma_vee_exact.push_back(rational_string(value));
    out.triangulation_faces = plot.triangulation_faces;
    for (const auto& indices : plot.subdivision_cells) {
      std::vector<std::pair<double, double>> cell;
      cell.reserve(indices.size());
      for (const auto index : indices) {
        const auto& point = configuration.points().at(index);
        cell.push_back({static_cast<double>(point.x), static_cast<double>(point.y)});
      }
      out.subdivision_cells.push_back(std::move(cell));
    }
    if (progress) progress(result.iterations, out.active_size, "Complete");
    return out;
  } catch (const std::exception& error) {
    out.certified = false;
    out.message = error.what();
    return out;
  }
}

}  // namespace gkz_gui
