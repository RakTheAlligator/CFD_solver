#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVelocityField.hpp"
#include "cfd/field/FaceFluxField.hpp"
#include "cfd/field/PressureCorrectionBoundaryConditions.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/mesh/Boundary.hpp"
#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/mesh/Types.hpp"
#include "cfd/meshing/RawMeshData.hpp"
#include "cfd/numerics/IncompressibleSimpleSolver.hpp"
#include "cfd/numerics/ScalarConvectionOperator.hpp"

#include "support/VerificationStatistics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace
{

using cfd::verification::ErrorAccumulator;
using cfd::verification::observed_order;

constexpr double domain_length{4.0};
constexpr double domain_height{1.0};
constexpr double density{1.0};
constexpr double dynamic_viscosity{0.1};
constexpr double inlet_pressure{0.04};
constexpr double outlet_pressure{0.0};

constexpr cfd::BoundaryId bottom_boundary_id{0};
constexpr cfd::BoundaryId right_boundary_id{1};
constexpr cfd::BoundaryId top_boundary_id{2};
constexpr cfd::BoundaryId left_boundary_id{3};

struct GridLevel
{
    cfd::Index x_cell_count;
    cfd::Index y_cell_count;
};

constexpr std::array grid_levels{
    GridLevel{16, 4},
    GridLevel{32, 8},
    GridLevel{64, 16},
};

struct LevelResult
{
    cfd::Index cell_count{};
    double mesh_spacing{};
    cfd::Index simple_iterations{};
    double velocity_relative_change{};
    double continuity_relative_residual{};
    double maximum_mass_imbalance{};
    double relative_u_error{};
    double normalized_pressure_error{};
    double maximum_absolute_v{};
    double inlet_mass_flow{};
    double outlet_mass_flow{};
    double relative_mass_flow_mismatch{};
    std::optional<double> observed_u_order;
};

[[nodiscard]]
double exact_u(const double y) noexcept
{
    const double pressure_gradient{(outlet_pressure - inlet_pressure) / domain_length};
    return -pressure_gradient * y * (domain_height - y) / (2.0 * dynamic_viscosity);
}

[[nodiscard]]
double exact_pressure(const double x) noexcept
{
    return inlet_pressure + (outlet_pressure - inlet_pressure) * x / domain_length;
}

[[nodiscard]]
cfd::RawMeshData make_channel_raw_mesh(const GridLevel &level)
{
    const cfd::Index x_node_count{level.x_cell_count + 1};
    const cfd::Index y_node_count{level.y_cell_count + 1};
    cfd::RawMeshData raw_mesh;
    raw_mesh.nodes.reserve(x_node_count * y_node_count);
    for (cfd::Index j = 0; j < y_node_count; ++j)
    {
        for (cfd::Index i = 0; i < x_node_count; ++i)
        {
            raw_mesh.nodes.push_back(
                {domain_length * static_cast<double>(i) / static_cast<double>(level.x_cell_count),
                 domain_height * static_cast<double>(j) / static_cast<double>(level.y_cell_count)});
        }
    }

    const cfd::Index cell_count{level.x_cell_count * level.y_cell_count};
    raw_mesh.cell_types.reserve(cell_count);
    raw_mesh.cell_nodes.reserve(4 * cell_count);
    raw_mesh.cell_node_offsets.reserve(cell_count + 1);
    raw_mesh.cell_node_offsets.push_back(0);
    for (cfd::Index j = 0; j < level.y_cell_count; ++j)
    {
        for (cfd::Index i = 0; i < level.x_cell_count; ++i)
        {
            const cfd::Index lower_left{j * x_node_count + i};
            raw_mesh.cell_types.push_back(cfd::CellType::Quadrilateral);
            raw_mesh.cell_nodes.insert(
                raw_mesh.cell_nodes.end(),
                {lower_left, lower_left + 1, lower_left + x_node_count + 1, lower_left + x_node_count});
            raw_mesh.cell_node_offsets.push_back(raw_mesh.cell_nodes.size());
        }
    }

    raw_mesh.boundary_groups = {
        {bottom_boundary_id, "bottom"},
        {right_boundary_id, "right"},
        {top_boundary_id, "top"},
        {left_boundary_id, "left"},
    };
    raw_mesh.boundary_edges.reserve(2 * (level.x_cell_count + level.y_cell_count));
    for (cfd::Index i = 0; i < level.x_cell_count; ++i)
    {
        raw_mesh.boundary_edges.push_back({{i, i + 1}, bottom_boundary_id});
    }
    for (cfd::Index j = 0; j < level.y_cell_count; ++j)
    {
        const cfd::Index lower_right{j * x_node_count + level.x_cell_count};
        raw_mesh.boundary_edges.push_back({{lower_right, lower_right + x_node_count}, right_boundary_id});
    }
    const cfd::Index top_row{level.y_cell_count * x_node_count};
    for (cfd::Index i = level.x_cell_count; i > 0; --i)
    {
        raw_mesh.boundary_edges.push_back({{top_row + i, top_row + i - 1}, top_boundary_id});
    }
    for (cfd::Index j = level.y_cell_count; j > 0; --j)
    {
        raw_mesh.boundary_edges.push_back({{j * x_node_count, (j - 1) * x_node_count}, left_boundary_id});
    }
    return raw_mesh;
}

[[nodiscard]]
cfd::ScalarBoundaryConditions velocity_boundary_conditions(const cfd::Mesh &mesh)
{
    return {mesh.boundary_groups().size(),
            {
                {cfd::ScalarBoundaryConditionType::Dirichlet, 0.0},
                {cfd::ScalarBoundaryConditionType::Neumann, 0.0},
                {cfd::ScalarBoundaryConditionType::Dirichlet, 0.0},
                {cfd::ScalarBoundaryConditionType::Neumann, 0.0},
            }};
}

[[nodiscard]]
cfd::ScalarBoundaryConditions pressure_boundary_conditions(const cfd::Mesh &mesh)
{
    return {mesh.boundary_groups().size(),
            {
                {cfd::ScalarBoundaryConditionType::Neumann, 0.0},
                {cfd::ScalarBoundaryConditionType::Dirichlet, outlet_pressure},
                {cfd::ScalarBoundaryConditionType::Neumann, 0.0},
                {cfd::ScalarBoundaryConditionType::Dirichlet, inlet_pressure},
            }};
}

[[nodiscard]]
cfd::PressureCorrectionBoundaryConditions pressure_correction_boundary_conditions(const cfd::Mesh &mesh)
{
    return {mesh.boundary_groups().size(),
            {
                cfd::PressureCorrectionBoundaryConditionType::FixedMassFlux,
                cfd::PressureCorrectionBoundaryConditionType::FixedPressure,
                cfd::PressureCorrectionBoundaryConditionType::FixedMassFlux,
                cfd::PressureCorrectionBoundaryConditionType::FixedPressure,
            }};
}

[[nodiscard]]
cfd::IncompressibleSimpleOptions simple_options()
{
    return {
        .maximum_iterations = 300,
        .momentum_relaxation_factor = 1.0,
        .pressure_relaxation_factor = 0.3,
        .velocity_relative_tolerance = 1.0e-10,
        .continuity_relative_tolerance = 1.0e-10,
        .momentum_linear_solver = {.relative_tolerance = 1.0e-12, .maximum_iterations = 5000},
        .pressure_correction_linear_solver = {.relative_tolerance = 1.0e-12, .maximum_iterations = 5000},
    };
}

[[nodiscard]]
double relative_mass_flow_mismatch(const double inlet_mass_flow, const double outlet_mass_flow) noexcept
{
    const double scale{std::max(std::abs(inlet_mass_flow), std::abs(outlet_mass_flow))};
    const double difference{std::abs(inlet_mass_flow - outlet_mass_flow)};
    if (scale == 0.0)
    {
        return difference == 0.0 ? 0.0 : std::numeric_limits<double>::infinity();
    }
    return difference / scale;
}

[[nodiscard]]
LevelResult solve_level(const GridLevel &level)
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_channel_raw_mesh(level))};
    const cfd::Mesh &mesh{build_result.mesh};
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        pressure[cell_id] = exact_pressure(mesh.cell_centers()[cell_id].x);
    }

    const cfd::ScalarBoundaryConditions velocity_conditions{velocity_boundary_conditions(mesh)};
    const cfd::ScalarBoundaryConditions pressure_conditions{pressure_boundary_conditions(mesh)};
    const cfd::PressureCorrectionBoundaryConditions pressure_correction_conditions{
        pressure_correction_boundary_conditions(mesh)};
    cfd::IncompressibleSimpleSolver solver{mesh, density, dynamic_viscosity, cfd::ScalarConvectionScheme::Linear,
                                           simple_options()};
    const cfd::IncompressibleSimpleResult result{solver.solve(velocity_conditions, velocity_conditions,
                                                              pressure_conditions, pressure_correction_conditions,
                                                              velocity, pressure, mass_flux)};
    if (!result.converged)
    {
        throw std::runtime_error("Poiseuille SIMPLE solve did not converge.");
    }

    ErrorAccumulator u_error;
    ErrorAccumulator u_reference;
    ErrorAccumulator pressure_error;
    ErrorAccumulator pressure_reference;
    double maximum_absolute_v{};
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        const double area{mesh.cell_areas()[cell_id]};
        const double exact_u_value{exact_u(mesh.cell_centers()[cell_id].y)};
        const double exact_pressure_value{exact_pressure(mesh.cell_centers()[cell_id].x)};
        u_error.add(area, velocity.u()[cell_id] - exact_u_value);
        u_reference.add(area, exact_u_value);
        pressure_error.add(area, pressure[cell_id] - exact_pressure_value);
        pressure_reference.add(area, exact_pressure_value);
        maximum_absolute_v = std::max(maximum_absolute_v, std::abs(velocity.v()[cell_id]));
    }

    double inlet_mass_flow{};
    double outlet_mass_flow{};
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (!mesh.face_adjacencies()[face_id].is_boundary())
        {
            continue;
        }
        switch (mesh.face_boundary_ids()[face_id])
        {
        case left_boundary_id:
            inlet_mass_flow -= mass_flux[face_id];
            break;
        case right_boundary_id:
            outlet_mass_flow += mass_flux[face_id];
            break;
        case bottom_boundary_id:
        case top_boundary_id:
            if (mass_flux[face_id] != 0.0)
            {
                throw std::runtime_error("Poiseuille wall FixedMassFlux value changed from zero.");
            }
            break;
        default:
            throw std::runtime_error("Poiseuille mesh contains an unexpected boundary ID.");
        }
    }

    const double u_reference_rms{u_reference.finish().area_weighted_rms_error};
    const double pressure_reference_rms{pressure_reference.finish().area_weighted_rms_error};
    return {
        .cell_count = mesh.cell_count(),
        .mesh_spacing = domain_height / static_cast<double>(level.y_cell_count),
        .simple_iterations = result.iteration_count,
        .velocity_relative_change = result.velocity_relative_change,
        .continuity_relative_residual = result.continuity_relative_residual,
        .maximum_mass_imbalance = result.maximum_mass_imbalance,
        .relative_u_error = u_error.finish().area_weighted_rms_error / u_reference_rms,
        .normalized_pressure_error = pressure_error.finish().area_weighted_rms_error / pressure_reference_rms,
        .maximum_absolute_v = maximum_absolute_v,
        .inlet_mass_flow = inlet_mass_flow,
        .outlet_mass_flow = outlet_mass_flow,
        .relative_mass_flow_mismatch = relative_mass_flow_mismatch(inlet_mass_flow, outlet_mass_flow),
        .observed_u_order = std::nullopt,
    };
}

void validate_results(const std::vector<LevelResult> &results)
{
    for (std::size_t index = 1; index < results.size(); ++index)
    {
        if (!(results[index].relative_u_error < results[index - 1].relative_u_error))
        {
            throw std::runtime_error("Poiseuille relative u error did not decrease under refinement.");
        }
    }
    for (const LevelResult &result : results)
    {
        if (result.maximum_absolute_v > 1.0e-10 || result.continuity_relative_residual > 1.0e-10 ||
            result.maximum_mass_imbalance > 1.0e-12 || result.relative_mass_flow_mismatch > 1.0e-10)
        {
            throw std::runtime_error("Poiseuille velocity or continuity diagnostics exceed their tolerances.");
        }
    }
    if (results.back().relative_u_error > 1.0e-2)
    {
        throw std::runtime_error("Poiseuille fine-grid relative u error is excessive.");
    }
}

void print_results(const std::vector<LevelResult> &results)
{
    const double pressure_gradient{(outlet_pressure - inlet_pressure) / domain_length};
    const double mean_velocity{domain_height * domain_height * (-pressure_gradient) / (12.0 * dynamic_viscosity)};
    const double maximum_velocity{domain_height * domain_height * (-pressure_gradient) / (8.0 * dynamic_viscosity)};
    std::cout << "Steady incompressible SIMPLE v1 - plane Poiseuille flow - Cartesian QUAD\n"
              << "L=" << domain_length << ", H=" << domain_height << ", rho=" << density << ", mu=" << dynamic_viscosity
              << ", p_in=" << inlet_pressure << ", p_out=" << outlet_pressure << ", U_mean=" << mean_velocity
              << ", u_max=" << maximum_velocity << "\n\n"
              << std::left << std::setw(9) << "cells" << std::setw(8) << "iters" << std::setw(14) << "rel_L2(u)"
              << std::setw(10) << "order" << std::setw(14) << "max|v|" << std::setw(14) << "r_U" << std::setw(14)
              << "r_cont" << std::setw(14) << "max|Rm|" << std::setw(14) << "flow_mismatch" << std::setw(14)
              << "p_error" << '\n';
    for (const LevelResult &result : results)
    {
        std::cout << std::left << std::setw(9) << result.cell_count << std::setw(8) << result.simple_iterations
                  << std::scientific << std::setprecision(5) << std::setw(14) << result.relative_u_error;
        if (result.observed_u_order.has_value())
        {
            std::cout << std::fixed << std::setprecision(3) << std::setw(10) << *result.observed_u_order;
        }
        else
        {
            std::cout << std::setw(10) << "-";
        }
        std::cout << std::scientific << std::setprecision(5) << std::setw(14) << result.maximum_absolute_v
                  << std::setw(14) << result.velocity_relative_change << std::setw(14)
                  << result.continuity_relative_residual << std::setw(14) << result.maximum_mass_imbalance
                  << std::setw(14) << result.relative_mass_flow_mismatch << std::setw(14)
                  << result.normalized_pressure_error << '\n';
    }

    const LevelResult &finest{results.back()};
    std::cout << "\nFinest-grid inlet mass flow:  " << std::scientific << finest.inlet_mass_flow
              << "\nFinest-grid outlet mass flow: " << finest.outlet_mass_flow
              << "\nFinal velocity relative change: " << finest.velocity_relative_change << '\n';
}

} // namespace

int main()
{
    try
    {
        std::vector<LevelResult> results;
        results.reserve(grid_levels.size());
        for (const GridLevel &level : grid_levels)
        {
            results.push_back(solve_level(level));
        }
        for (std::size_t index = 1; index < results.size(); ++index)
        {
            results[index].observed_u_order =
                observed_order(results[index - 1].relative_u_error, results[index].relative_u_error,
                               results[index - 1].mesh_spacing, results[index].mesh_spacing);
        }
        validate_results(results);
        print_results(results);
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "Incompressible SIMPLE Poiseuille verification failed: " << error.what() << '\n';
        return 1;
    }
}
