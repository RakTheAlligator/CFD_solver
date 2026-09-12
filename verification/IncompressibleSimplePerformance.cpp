#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVelocityField.hpp"
#include "cfd/field/FaceFluxField.hpp"
#include "cfd/field/PressureCorrectionBoundaryConditions.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/mesh/Boundary.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/mesh/Types.hpp"
#include "cfd/meshing/GmshMesher.hpp"
#include "cfd/meshing/RectangleGeometry.hpp"
#include "cfd/numerics/IncompressibleSimpleSolver.hpp"
#include "cfd/numerics/ScalarConvectionOperator.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace
{

using BenchmarkClock = std::chrono::steady_clock;

constexpr double domain_length{4.0};
constexpr double domain_height{1.0};
constexpr double density{1.0};
constexpr double dynamic_viscosity{0.1};
constexpr double inlet_pressure{0.04};
constexpr double outlet_pressure{0.0};

constexpr std::array target_mesh_sizes{
    0.20,
    0.14,
    0.10,
};

constexpr std::size_t measured_run_count{5};

struct RunMeasurement
{
    double solver_setup_seconds{};
    cfd::Index iteration_count{};

    double velocity_relative_change{};
    double rhie_chow_flux_relative_residual{};
    double provisional_continuity{};
    double corrected_continuity{};

    cfd::SimpleTimingBreakdown timings;
};

struct LinearIterationStatistics
{
    std::uint64_t simple_iteration_count{};
    std::uint64_t u_momentum_iteration_count{};
    std::uint64_t v_momentum_iteration_count{};
    std::uint64_t pressure_correction_iteration_count{};
};

struct PhaseTiming
{
    std::string_view name;
    double seconds{};
};

struct TimingStatistics
{
    double minimum_seconds{};
    double median_seconds{};
    double mean_seconds{};
    double maximum_seconds{};
};

struct MeshBenchmarkResult
{
    double target_mesh_size{};

    cfd::Index node_count{};
    cfd::Index cell_count{};
    cfd::Index face_count{};

    double gmsh_seconds{};
    double mesh_build_seconds{};
    cfd::MeshBuildTimings mesh_build_timings;

    RunMeasurement warm_up;
    LinearIterationStatistics warm_up_linear_iterations;

    std::vector<RunMeasurement> measurements;
};

struct LevelSummary
{
    double target_mesh_size{};

    cfd::Index cell_count{};
    cfd::Index face_count{};
    cfd::Index simple_iteration_count{};

    double median_total_seconds{};
    double seconds_per_simple_iteration{};

    double average_u_linear_iterations{};
    double average_v_linear_iterations{};
    double average_pressure_correction_linear_iterations{};
};

[[nodiscard]]
double elapsed_seconds_since(const BenchmarkClock::time_point start) noexcept
{
    return std::chrono::duration<double>(BenchmarkClock::now() - start).count();
}

void accumulate_linear_iterations(LinearIterationStatistics &statistics, const cfd::SimpleIterationInfo &info) noexcept
{
    ++statistics.simple_iteration_count;

    statistics.u_momentum_iteration_count += static_cast<std::uint64_t>(info.u_solve.iteration_count);

    statistics.v_momentum_iteration_count += static_cast<std::uint64_t>(info.v_solve.iteration_count);

    statistics.pressure_correction_iteration_count +=
        static_cast<std::uint64_t>(info.pressure_correction_solve.iteration_count);
}

[[nodiscard]]
double average_linear_iterations(const std::uint64_t linear_iteration_count,
                                 const std::uint64_t simple_iteration_count) noexcept
{
    if (simple_iteration_count == 0)
    {
        return 0.0;
    }

    return static_cast<double>(linear_iteration_count) / static_cast<double>(simple_iteration_count);
}

[[nodiscard]]
cfd::ScalarBoundaryConditions velocity_boundary_conditions(const cfd::Mesh &mesh)
{
    std::vector<cfd::ScalarBoundaryCondition> conditions;
    conditions.reserve(mesh.boundary_groups().size());

    for (const cfd::BoundaryGroup &group : mesh.boundary_groups())
    {
        if (group.name == "wall")
        {
            conditions.emplace_back(cfd::ScalarBoundaryConditionType::Dirichlet, 0.0);
        }
        else if (group.name == "inlet" || group.name == "outlet")
        {
            conditions.emplace_back(cfd::ScalarBoundaryConditionType::Neumann, 0.0);
        }
        else
        {
            throw std::runtime_error("SIMPLE performance mesh has an unexpected boundary group.");
        }
    }

    return {
        mesh.boundary_groups().size(),
        std::move(conditions),
    };
}

[[nodiscard]]
cfd::ScalarBoundaryConditions pressure_boundary_conditions(const cfd::Mesh &mesh)
{
    std::vector<cfd::ScalarBoundaryCondition> conditions;
    conditions.reserve(mesh.boundary_groups().size());

    for (const cfd::BoundaryGroup &group : mesh.boundary_groups())
    {
        if (group.name == "inlet")
        {
            conditions.emplace_back(cfd::ScalarBoundaryConditionType::Dirichlet, inlet_pressure);
        }
        else if (group.name == "outlet")
        {
            conditions.emplace_back(cfd::ScalarBoundaryConditionType::Dirichlet, outlet_pressure);
        }
        else if (group.name == "wall")
        {
            conditions.emplace_back(cfd::ScalarBoundaryConditionType::Neumann, 0.0);
        }
        else
        {
            throw std::runtime_error("SIMPLE performance mesh has an unexpected boundary group.");
        }
    }

    return {
        mesh.boundary_groups().size(),
        std::move(conditions),
    };
}

[[nodiscard]]
cfd::PressureCorrectionBoundaryConditions pressure_correction_boundary_conditions(const cfd::Mesh &mesh)
{
    std::vector<cfd::PressureCorrectionBoundaryConditionType> conditions;
    conditions.reserve(mesh.boundary_groups().size());

    for (const cfd::BoundaryGroup &group : mesh.boundary_groups())
    {
        if (group.name == "inlet" || group.name == "outlet")
        {
            conditions.push_back(cfd::PressureCorrectionBoundaryConditionType::FixedPressure);
        }
        else if (group.name == "wall")
        {
            conditions.push_back(cfd::PressureCorrectionBoundaryConditionType::FixedMassFlux);
        }
        else
        {
            throw std::runtime_error("SIMPLE performance mesh has an unexpected boundary group.");
        }
    }

    return {
        mesh.boundary_groups().size(),
        std::move(conditions),
    };
}

[[nodiscard]]
cfd::IncompressibleSimpleOptions simple_options()
{
    return {
        .maximum_iterations = 4000,
        .momentum_relaxation_factor = 1.0,
        .pressure_relaxation_factor = 0.1,
        .rhie_chow_flux_relaxation_factor = 0.1,
        .velocity_relative_tolerance = 1.0e-10,
        .rhie_chow_flux_relative_tolerance = 1.0e-8,
        .continuity_relative_tolerance = 1.0e-10,
        .momentum_linear_solver =
            {
                .relative_tolerance = 1.0e-12,
                .maximum_iterations = 5000,
            },
        .pressure_correction_linear_solver =
            {
                .relative_tolerance = 1.0e-12,
                .maximum_iterations = 5000,
            },
    };
}

[[nodiscard]]
RunMeasurement run_benchmark(const cfd::Mesh &mesh, const cfd::ScalarBoundaryConditions &velocity_conditions,
                             const cfd::ScalarBoundaryConditions &pressure_conditions,
                             const cfd::PressureCorrectionBoundaryConditions &pressure_correction_conditions,
                             const cfd::SimpleIterationCallback &iteration_callback = {})
{
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    cfd::FaceFluxField mass_flux{mesh.face_count()};

    const BenchmarkClock::time_point setup_start{BenchmarkClock::now()};

    cfd::IncompressibleSimpleSolver solver{
        mesh, density, dynamic_viscosity, cfd::ScalarConvectionScheme::Linear, simple_options(),
    };

    const double setup_seconds{elapsed_seconds_since(setup_start)};

    const cfd::IncompressibleSimpleResult result{solver.solve(velocity_conditions, velocity_conditions,
                                                              pressure_conditions, pressure_correction_conditions,
                                                              velocity, pressure, mass_flux, iteration_callback)};

    if (!result.converged)
    {
        throw std::runtime_error("SIMPLE performance run did not converge.");
    }

    return {
        .solver_setup_seconds = setup_seconds,
        .iteration_count = result.iteration_count,
        .velocity_relative_change = result.velocity_relative_change,
        .rhie_chow_flux_relative_residual = result.rhie_chow_flux_relative_residual,
        .provisional_continuity = result.provisional_continuity_relative_residual,
        .corrected_continuity = result.continuity_relative_residual,
        .timings = result.timings,
    };
}

[[nodiscard]]
TimingStatistics timing_statistics(const std::vector<RunMeasurement> &measurements)
{
    if (measurements.empty())
    {
        throw std::runtime_error("SIMPLE performance campaign has no measured runs.");
    }

    std::vector<double> total_times;
    total_times.reserve(measurements.size());

    double sum{};

    for (const RunMeasurement &measurement : measurements)
    {
        const double total_seconds{measurement.timings.total_seconds};

        total_times.push_back(total_seconds);
        sum += total_seconds;
    }

    std::ranges::sort(total_times);

    return {
        .minimum_seconds = total_times.front(),
        .median_seconds = total_times[total_times.size() / 2],
        .mean_seconds = sum / static_cast<double>(total_times.size()),
        .maximum_seconds = total_times.back(),
    };
}

[[nodiscard]]
const RunMeasurement &representative_measurement(const std::vector<RunMeasurement> &measurements,
                                                 const double median_seconds)
{
    if (measurements.empty())
    {
        throw std::runtime_error("Cannot select a representative SIMPLE performance run.");
    }

    const auto iterator{std::ranges::min_element(measurements, {}, [median_seconds](const RunMeasurement &run) {
        return std::abs(run.timings.total_seconds - median_seconds);
    })};

    return *iterator;
}

[[nodiscard]]
std::array<PhaseTiming, 12> representative_phases(const cfd::SimpleTimingBreakdown &timings)
{
    return {
        PhaseTiming{"gradient reconstruction", timings.velocity_gradient_reconstruction_seconds +
                                                   timings.pressure_gradient_reconstruction_seconds +
                                                   timings.pressure_correction_gradient_reconstruction_seconds},
        PhaseTiming{"momentum assembly", timings.momentum_assembly_seconds},
        PhaseTiming{"momentum sparse preparation", timings.momentum_matrix_preparation_seconds},
        PhaseTiming{"u-momentum linear solve", timings.u_momentum_linear_solve_seconds},
        PhaseTiming{"v-momentum linear solve", timings.v_momentum_linear_solve_seconds},
        PhaseTiming{"momentum pressure response", timings.momentum_pressure_response_seconds},
        PhaseTiming{"Rhie-Chow interpolation", timings.rhie_chow_interpolation_seconds},
        PhaseTiming{"pressure-correction assembly", timings.pressure_correction_assembly_seconds},
        PhaseTiming{"pressure sparse preparation", timings.pressure_correction_matrix_preparation_seconds},
        PhaseTiming{"pressure-correction solve", timings.pressure_correction_linear_solve_seconds},
        PhaseTiming{"field corrections", timings.field_correction_seconds},
        PhaseTiming{"convergence/diagnostics", timings.momentum_residual_diagnostics_seconds +
                                                   timings.provisional_continuity_diagnostics_seconds +
                                                   timings.convergence_diagnostics_seconds},
    };
}

void print_phase_breakdown(const RunMeasurement &measurement)
{
    std::array phases{representative_phases(measurement.timings)};

    std::ranges::sort(phases, {}, &PhaseTiming::seconds);

    double measured_phase_sum{};

    for (const PhaseTiming &phase : phases)
    {
        measured_phase_sum += phase.seconds;
    }

    std::cout << "\nRepresentative run phase breakdown "
                 "(closest to median total):\n";

    for (auto iterator = phases.rbegin(); iterator != phases.rend(); ++iterator)
    {
        std::cout << "  " << std::left << std::setw(32) << iterator->name << std::right << std::setw(11)
                  << iterator->seconds << " s  " << std::setw(7)
                  << 100.0 * iterator->seconds / measurement.timings.total_seconds << " %\n";
    }

    const double other_seconds{std::max(0.0, measurement.timings.total_seconds - measured_phase_sum)};

    std::cout << "  " << std::left << std::setw(32) << "other orchestration" << std::right << std::setw(11)
              << other_seconds << " s  " << std::setw(7) << 100.0 * other_seconds / measurement.timings.total_seconds
              << " %\n"
              << "\nGradient detail:\n"
              << "  velocity gradients: " << measurement.timings.velocity_gradient_reconstruction_seconds << " s\n"
              << "  pressure gradient:  " << measurement.timings.pressure_gradient_reconstruction_seconds << " s\n"
              << "  p-prime gradient:   " << measurement.timings.pressure_correction_gradient_reconstruction_seconds
              << " s\n"
              << "\nDiagnostic detail:\n"
              << "  momentum residual SpMV: " << measurement.timings.momentum_residual_diagnostics_seconds << " s  "
              << 100.0 * measurement.timings.momentum_residual_diagnostics_seconds / measurement.timings.total_seconds
              << " %\n"
              << "  provisional continuity: " << measurement.timings.provisional_continuity_diagnostics_seconds
              << " s\n"
              << "  corrected/convergence:  " << measurement.timings.convergence_diagnostics_seconds << " s\n";
}

[[nodiscard]]
MeshBenchmarkResult benchmark_mesh_level(const double target_mesh_size)
{
    const BenchmarkClock::time_point gmsh_start{BenchmarkClock::now()};

    cfd::RawMeshData raw_mesh{cfd::generate_mesh(
        {
            .length = domain_length,
            .height = domain_height,
        },
        {
            .mesh_size = target_mesh_size,
            .cell_type = cfd::CellType::Triangle,
        })};

    const double gmsh_seconds{elapsed_seconds_since(gmsh_start)};

    const BenchmarkClock::time_point mesh_build_start{BenchmarkClock::now()};

    cfd::MeshBuildResult build_result{cfd::build_mesh(std::move(raw_mesh))};

    const double mesh_build_seconds{elapsed_seconds_since(mesh_build_start)};

    const cfd::Mesh &mesh{build_result.mesh};

    const cfd::ScalarBoundaryConditions velocity_conditions{velocity_boundary_conditions(mesh)};

    const cfd::ScalarBoundaryConditions pressure_conditions{pressure_boundary_conditions(mesh)};

    const cfd::PressureCorrectionBoundaryConditions pressure_correction_conditions{
        pressure_correction_boundary_conditions(mesh)};

    LinearIterationStatistics linear_iterations;

    const cfd::SimpleIterationCallback iteration_callback{[&linear_iterations](const cfd::SimpleIterationInfo &info) {
        accumulate_linear_iterations(linear_iterations, info);
    }};

    const RunMeasurement warm_up{run_benchmark(mesh, velocity_conditions, pressure_conditions,
                                               pressure_correction_conditions, iteration_callback)};

    if (linear_iterations.simple_iteration_count != static_cast<std::uint64_t>(warm_up.iteration_count))
    {
        throw std::runtime_error("Warm-up callback did not observe every SIMPLE iteration.");
    }

    std::vector<RunMeasurement> measurements;
    measurements.reserve(measured_run_count);

    for (std::size_t run_id = 0; run_id < measured_run_count; ++run_id)
    {
        RunMeasurement measurement{
            run_benchmark(mesh, velocity_conditions, pressure_conditions, pressure_correction_conditions)};

        if (measurement.iteration_count != warm_up.iteration_count)
        {
            throw std::runtime_error("Measured SIMPLE iteration count differs from the warm-up.");
        }

        measurements.push_back(measurement);
    }

    return {
        .target_mesh_size = target_mesh_size,
        .node_count = mesh.node_count(),
        .cell_count = mesh.cell_count(),
        .face_count = mesh.face_count(),
        .gmsh_seconds = gmsh_seconds,
        .mesh_build_seconds = mesh_build_seconds,
        .mesh_build_timings = build_result.timings,
        .warm_up = warm_up,
        .warm_up_linear_iterations = linear_iterations,
        .measurements = std::move(measurements),
    };
}

[[nodiscard]]
LevelSummary make_level_summary(const MeshBenchmarkResult &result)
{
    const TimingStatistics statistics{timing_statistics(result.measurements)};

    if (result.warm_up.iteration_count == 0)
    {
        throw std::runtime_error("SIMPLE performance run completed zero outer iterations.");
    }

    const double simple_iteration_count{static_cast<double>(result.warm_up.iteration_count)};

    const LinearIterationStatistics &linear{result.warm_up_linear_iterations};

    return {
        .target_mesh_size = result.target_mesh_size,
        .cell_count = result.cell_count,
        .face_count = result.face_count,
        .simple_iteration_count = result.warm_up.iteration_count,
        .median_total_seconds = statistics.median_seconds,
        .seconds_per_simple_iteration = statistics.median_seconds / simple_iteration_count,
        .average_u_linear_iterations =
            average_linear_iterations(linear.u_momentum_iteration_count, linear.simple_iteration_count),
        .average_v_linear_iterations =
            average_linear_iterations(linear.v_momentum_iteration_count, linear.simple_iteration_count),
        .average_pressure_correction_linear_iterations =
            average_linear_iterations(linear.pressure_correction_iteration_count, linear.simple_iteration_count),
    };
}

void print_level_measurements(const MeshBenchmarkResult &result)
{
    const TimingStatistics statistics{timing_statistics(result.measurements)};

    const RunMeasurement &representative{representative_measurement(result.measurements, statistics.median_seconds)};

    const double topology_seconds{result.mesh_build_timings.topology.count() / 1000.0};

    const double geometry_seconds{result.mesh_build_timings.geometry.count() / 1000.0};

    const double other_build_seconds{std::max(0.0, result.mesh_build_seconds - topology_seconds - geometry_seconds)};

    const LinearIterationStatistics &linear{result.warm_up_linear_iterations};

    std::cout << std::scientific << std::setprecision(6)
              << "\n============================================================\n"
              << "Mesh target size: " << result.target_mesh_size << '\n'
              << "  nodes: " << result.node_count << '\n'
              << "  cells: " << result.cell_count << '\n'
              << "  faces: " << result.face_count << "\n\n"
              << "Preprocessing "
                 "(not included in SIMPLE totals):\n"
              << "  Gmsh generation/extraction:       " << result.gmsh_seconds << " s\n"
              << "  topology build/validation:        " << topology_seconds << " s\n"
              << "  geometry build/validation:        " << geometry_seconds << " s\n"
              << "  raw validation/final transfer remainder: " << other_build_seconds << " s\n"
              << "  total Mesh build:                 " << result.mesh_build_seconds << " s\n"
              << "  VTU export:                       "
                 "not performed\n\n"
              << "Warm-up "
                 "(excluded from measured timings):\n"
              << "  SIMPLE iterations: " << result.warm_up.iteration_count << '\n'
              << "  solver setup: " << result.warm_up.solver_setup_seconds << " s\n"
              << "  SIMPLE total: " << result.warm_up.timings.total_seconds << " s\n"
              << "  diagnostics: r_U=" << result.warm_up.velocity_relative_change
              << ", r_RC=" << result.warm_up.rhie_chow_flux_relative_residual
              << ", r_cont*=" << result.warm_up.provisional_continuity
              << ", r_cont=" << result.warm_up.corrected_continuity << '\n'
              << "  mean u-BiCGSTAB iterations/SIMPLE: "
              << average_linear_iterations(linear.u_momentum_iteration_count, linear.simple_iteration_count) << '\n'
              << "  mean v-BiCGSTAB iterations/SIMPLE: "
              << average_linear_iterations(linear.v_momentum_iteration_count, linear.simple_iteration_count) << '\n'
              << "  mean p'-CG iterations/SIMPLE:      "
              << average_linear_iterations(linear.pressure_correction_iteration_count, linear.simple_iteration_count)
              << "\n\n"
              << "Measured SIMPLE totals:\n";

    for (std::size_t run_id = 0; run_id < result.measurements.size(); ++run_id)
    {
        std::cout << "  run " << run_id + 1 << ": " << result.measurements[run_id].timings.total_seconds
                  << " s, iterations: " << result.measurements[run_id].iteration_count << '\n';
    }

    std::cout << "  minimum: " << statistics.minimum_seconds << " s\n"
              << "  median:  " << statistics.median_seconds << " s\n"
              << "  mean:    " << statistics.mean_seconds << " s\n"
              << "  maximum: " << statistics.maximum_seconds << " s\n"
              << "  representative solver setup: " << representative.solver_setup_seconds << " s\n"
              << "  representative diagnostics: r_U=" << representative.velocity_relative_change
              << ", r_RC=" << representative.rhie_chow_flux_relative_residual
              << ", r_cont*=" << representative.provisional_continuity
              << ", r_cont=" << representative.corrected_continuity << '\n';

    print_phase_breakdown(representative);
}

void print_campaign_summary(const std::vector<MeshBenchmarkResult> &results)
{
    if (results.empty())
    {
        throw std::runtime_error("SIMPLE performance campaign has no mesh levels.");
    }

    std::vector<LevelSummary> summaries;
    summaries.reserve(results.size());

    for (const MeshBenchmarkResult &result : results)
    {
        summaries.push_back(make_level_summary(result));
    }

    std::cout << "\n============================================================\n"
              << "SIMPLE multi-mesh scaling summary\n\n"
              << std::left << std::setw(12) << "h_target" << std::setw(12) << "cells" << std::setw(12) << "faces"
              << std::setw(12) << "SIMPLE" << std::setw(16) << "median[s]" << std::setw(16) << "s/SIMPLE"
              << std::setw(16) << "u_it/SIMPLE" << std::setw(16) << "v_it/SIMPLE" << std::setw(16) << "p'_it/SIMPLE"
              << '\n';

    for (const LevelSummary &summary : summaries)
    {
        std::cout << std::scientific << std::setprecision(5) << std::left << std::setw(12) << summary.target_mesh_size
                  << std::setw(12) << summary.cell_count << std::setw(12) << summary.face_count << std::setw(12)
                  << summary.simple_iteration_count << std::setw(16) << summary.median_total_seconds << std::setw(16)
                  << summary.seconds_per_simple_iteration << std::setw(16) << summary.average_u_linear_iterations
                  << std::setw(16) << summary.average_v_linear_iterations << std::setw(16)
                  << summary.average_pressure_correction_linear_iterations << '\n';
    }

    if (summaries.size() < 2)
    {
        return;
    }

    std::cout << "\nAdjacent-mesh empirical scaling\n\n"
              << std::left << std::setw(14) << "cell_ratio" << std::setw(14) << "time_ratio" << std::setw(14)
              << "time_exp" << std::setw(18) << "iter_cost_ratio" << std::setw(18) << "iter_cost_exp" << '\n';

    for (std::size_t index = 1; index < summaries.size(); ++index)
    {
        const LevelSummary &coarse{summaries[index - 1]};

        const LevelSummary &fine{summaries[index]};

        const double cell_ratio{static_cast<double>(fine.cell_count) / static_cast<double>(coarse.cell_count)};

        const double time_ratio{fine.median_total_seconds / coarse.median_total_seconds};

        const double iteration_cost_ratio{fine.seconds_per_simple_iteration / coarse.seconds_per_simple_iteration};

        if (!(cell_ratio > 1.0) || !(time_ratio > 0.0) || !(iteration_cost_ratio > 0.0))
        {
            throw std::runtime_error("Invalid values encountered while computing performance scaling.");
        }

        const double logarithmic_cell_ratio{std::log(cell_ratio)};

        const double total_time_exponent{std::log(time_ratio) / logarithmic_cell_ratio};

        const double iteration_cost_exponent{std::log(iteration_cost_ratio) / logarithmic_cell_ratio};

        std::cout << std::scientific << std::setprecision(5) << std::left << std::setw(14) << cell_ratio
                  << std::setw(14) << time_ratio << std::setw(14) << total_time_exponent << std::setw(18)
                  << iteration_cost_ratio << std::setw(18) << iteration_cost_exponent << '\n';
    }
}

} // namespace

int main()
{
    try
    {
        std::cout << "SIMPLE Release multi-mesh performance campaign\n"
                  << "Triangle meshes, cold-start solve, " << measured_run_count << " measured runs per mesh\n";

        std::vector<MeshBenchmarkResult> results;
        results.reserve(target_mesh_sizes.size());

        for (const double target_mesh_size : target_mesh_sizes)
        {
            MeshBenchmarkResult result{benchmark_mesh_level(target_mesh_size)};

            print_level_measurements(result);

            results.push_back(std::move(result));
        }

        print_campaign_summary(results);

        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "SIMPLE performance campaign failed: " << error.what() << '\n';

        return 1;
    }
}
