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
constexpr double target_mesh_size{0.1};
constexpr std::size_t measured_run_count{5};

struct RunMeasurement
{
    double solver_setup_seconds{};
    cfd::Index iteration_count{};
    double velocity_relative_change{};
    double provisional_continuity{};
    double corrected_continuity{};
    cfd::SimpleTimingBreakdown timings;
};

struct PhaseTiming
{
    std::string_view name;
    double seconds{};
};

[[nodiscard]]
double elapsed_seconds_since(const BenchmarkClock::time_point start) noexcept
{
    return std::chrono::duration<double>(BenchmarkClock::now() - start).count();
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
    return {mesh.boundary_groups().size(), std::move(conditions)};
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
    return {mesh.boundary_groups().size(), std::move(conditions)};
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
    return {mesh.boundary_groups().size(), std::move(conditions)};
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
        .continuity_relative_tolerance = 1.0e-10,
        .momentum_linear_solver = {.relative_tolerance = 1.0e-12, .maximum_iterations = 5000},
        .pressure_correction_linear_solver = {.relative_tolerance = 1.0e-12, .maximum_iterations = 5000},
    };
}

[[nodiscard]]
RunMeasurement run_benchmark(const cfd::Mesh &mesh, const cfd::ScalarBoundaryConditions &velocity_conditions,
                             const cfd::ScalarBoundaryConditions &pressure_conditions,
                             const cfd::PressureCorrectionBoundaryConditions &pressure_correction_conditions)
{
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    cfd::FaceFluxField mass_flux{mesh.face_count()};

    const BenchmarkClock::time_point setup_start{BenchmarkClock::now()};
    cfd::IncompressibleSimpleSolver solver{mesh, density, dynamic_viscosity, cfd::ScalarConvectionScheme::Linear,
                                           simple_options()};
    const double setup_seconds{elapsed_seconds_since(setup_start)};

    const cfd::IncompressibleSimpleResult result{solver.solve(velocity_conditions, velocity_conditions,
                                                              pressure_conditions, pressure_correction_conditions,
                                                              velocity, pressure, mass_flux)};
    if (!result.converged)
    {
        throw std::runtime_error("SIMPLE performance run did not converge.");
    }
    return {
        .solver_setup_seconds = setup_seconds,
        .iteration_count = result.iteration_count,
        .velocity_relative_change = result.velocity_relative_change,
        .provisional_continuity = result.provisional_continuity_relative_residual,
        .corrected_continuity = result.continuity_relative_residual,
        .timings = result.timings,
    };
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

    std::cout << "\nRepresentative run phase breakdown (closest to median total):\n";
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

void print_measurements(const cfd::Mesh &mesh, const double gmsh_seconds, const double mesh_build_seconds,
                        const cfd::MeshBuildTimings &mesh_build_timings, const RunMeasurement &warm_up,
                        const std::vector<RunMeasurement> &measurements)
{
    std::vector<double> total_times;
    total_times.reserve(measurements.size());
    for (const RunMeasurement &measurement : measurements)
    {
        total_times.push_back(measurement.timings.total_seconds);
    }
    std::ranges::sort(total_times);
    const double minimum{total_times.front()};
    const double median{total_times[total_times.size() / 2]};
    const double maximum{total_times.back()};
    double sum{};
    for (const double total : total_times)
    {
        sum += total;
    }
    const double mean{sum / static_cast<double>(total_times.size())};

    const auto representative_iterator{std::ranges::min_element(measurements, {}, [median](const RunMeasurement &run) {
        return std::abs(run.timings.total_seconds - median);
    })};
    const RunMeasurement &representative{*representative_iterator};
    const double topology_seconds{mesh_build_timings.topology.count() / 1000.0};
    const double geometry_seconds{mesh_build_timings.geometry.count() / 1000.0};
    const double other_build_seconds{std::max(0.0, mesh_build_seconds - topology_seconds - geometry_seconds)};

    std::cout << std::scientific << std::setprecision(6) << "SIMPLE Release performance campaign\n"
              << "Mesh:\n"
              << "  target size: " << target_mesh_size << '\n'
              << "  nodes:       " << mesh.node_count() << '\n'
              << "  cells:       " << mesh.cell_count() << '\n'
              << "  faces:       " << mesh.face_count() << "\n\n"
              << "Preprocessing (not included in SIMPLE totals):\n"
              << "  Gmsh generation/extraction:       " << gmsh_seconds << " s\n"
              << "  topology build/validation:        " << topology_seconds << " s\n"
              << "  geometry build/validation:        " << geometry_seconds << " s\n"
              << "  raw validation/final transfer remainder: " << other_build_seconds << " s\n"
              << "  total Mesh build:                 " << mesh_build_seconds << " s\n"
              << "  VTU export:                       not performed\n\n"
              << "Warm-up:\n"
              << "  iterations: " << warm_up.iteration_count << ", setup: " << warm_up.solver_setup_seconds
              << " s, SIMPLE: " << warm_up.timings.total_seconds << " s\n\n"
              << "Measured SIMPLE totals:\n";
    for (std::size_t run_id = 0; run_id < measurements.size(); ++run_id)
    {
        std::cout << "  run " << run_id + 1 << ": " << measurements[run_id].timings.total_seconds
                  << " s, iterations: " << measurements[run_id].iteration_count << '\n';
    }
    std::cout << "  minimum: " << minimum << " s\n"
              << "  median:  " << median << " s\n"
              << "  mean:    " << mean << " s\n"
              << "  maximum: " << maximum << " s\n"
              << "  representative solver setup: " << representative.solver_setup_seconds << " s\n"
              << "  representative diagnostics: r_U=" << representative.velocity_relative_change
              << ", r_cont*=" << representative.provisional_continuity
              << ", r_cont=" << representative.corrected_continuity << '\n';
    print_phase_breakdown(representative);
}

} // namespace

int main()
{
    try
    {
        const BenchmarkClock::time_point gmsh_start{BenchmarkClock::now()};
        cfd::RawMeshData raw_mesh{
            cfd::generate_mesh({.length = domain_length, .height = domain_height},
                               {.mesh_size = target_mesh_size, .cell_type = cfd::CellType::Triangle})};
        const double gmsh_seconds{elapsed_seconds_since(gmsh_start)};

        const BenchmarkClock::time_point mesh_build_start{BenchmarkClock::now()};
        cfd::MeshBuildResult build_result{cfd::build_mesh(std::move(raw_mesh))};
        const double mesh_build_seconds{elapsed_seconds_since(mesh_build_start)};
        const cfd::Mesh &mesh{build_result.mesh};
        const cfd::ScalarBoundaryConditions velocity_conditions{velocity_boundary_conditions(mesh)};
        const cfd::ScalarBoundaryConditions pressure_conditions{pressure_boundary_conditions(mesh)};
        const cfd::PressureCorrectionBoundaryConditions pressure_correction_conditions{
            pressure_correction_boundary_conditions(mesh)};

        const RunMeasurement warm_up{
            run_benchmark(mesh, velocity_conditions, pressure_conditions, pressure_correction_conditions)};
        std::vector<RunMeasurement> measurements;
        measurements.reserve(measured_run_count);
        for (std::size_t run_id = 0; run_id < measured_run_count; ++run_id)
        {
            measurements.push_back(
                run_benchmark(mesh, velocity_conditions, pressure_conditions, pressure_correction_conditions));
        }

        print_measurements(mesh, gmsh_seconds, mesh_build_seconds, build_result.timings, warm_up, measurements);
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "SIMPLE performance campaign failed: " << error.what() << '\n';
        return 1;
    }
}
