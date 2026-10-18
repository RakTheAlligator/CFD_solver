#pragma once

#include "cfd/compressible/EulerSteadyMonitor.hpp"

#include <filesystem>
#include <functional>
#include <string_view>

namespace cfd
{
/// Configuration borrowed only for an I/O call. Full mesh metadata are compared,
/// not a cell count/hash. Time/step budgets may be extended on restart; the
/// physical model, CFL, retry policy and convergence criteria must match exactly.
struct EulerCheckpointConfiguration
{
    std::reference_wrapper<const Mesh> mesh;
    CaloricallyPerfectGas gas;
    std::reference_wrapper<const EulerBoundaryConditions> boundaries;
    AxisymmetricEulerOptions solver;
    EulerSteadyOptions monitor;
    std::string_view physical_context;
};

struct EulerCheckpointHistory
{
    EulerIntegralQuantities initial;
    std::array<long double, 4> integrated_boundary{}, balance_scales{};
    long double integrated_source{};
    double maximum_step_error{};
    std::vector<double> output_times{};
    Index next_output{};
};

/// Explicit, owning snapshot. No VTU conversion and no allocation per Euler step.
struct EulerCheckpoint
{
    AxisymmetricEulerRestartState solver;
    EulerSteadyRestartState monitor;
    EulerCheckpointHistory history;
};

/// Versioned exact text round-trip, including extended-precision ledger values.
/// Atomic no-replace publication on Linux (renameat2); existing checkpoints are
/// never overwritten. Parent must exist. No fsync/power-loss durability guarantee.
/// @throws std::runtime_error / std::filesystem::filesystem_error On invalid
///         data, incompatibility, corruption, truncation or I/O failure.
void write_euler_checkpoint(const std::filesystem::path &path, const EulerCheckpointConfiguration &configuration,
                            const EulerCheckpoint &checkpoint);
[[nodiscard]] EulerCheckpoint read_euler_checkpoint(const std::filesystem::path &path,
                                                    const EulerCheckpointConfiguration &configuration);
} // namespace cfd
