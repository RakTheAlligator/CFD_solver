#pragma once

#include <filesystem>

namespace cfd::app
{

/// Starts the repository live-convergence plotter as a detached process.
///
/// Launch failures are reported to stderr and never propagated to the
/// numerical application.
void launch_live_convergence_plotter(const std::filesystem::path &convergence_file) noexcept;

} // namespace cfd::app
