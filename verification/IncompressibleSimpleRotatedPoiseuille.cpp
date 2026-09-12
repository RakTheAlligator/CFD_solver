#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVelocityField.hpp"
#include "cfd/field/FaceFluxField.hpp"
#include "cfd/field/PressureCorrectionBoundaryConditions.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/math/Point2.hpp"
#include "cfd/math/Vector2.hpp"
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
#include <string>
#include <vector>

namespace
{

using cfd::verification::ErrorAccumulator;
using cfd::verification::ErrorStatistics;
using cfd::verification::observed_order;

constexpr double domain_length{4.0};
constexpr double domain_height{1.0};
constexpr double density{1.0};
constexpr double dynamic_viscosity{0.1};
constexpr double inlet_pressure{0.04};
constexpr double outlet_pressure{0.0};

constexpr double corrected_continuity_tolerance{1.0e-10};
constexpr double maximum_transverse_velocity_tolerance{1.0e-10};
constexpr double maximum_mass_imbalance_tolerance{1.0e-12};
constexpr double mass_flow_mismatch_tolerance{1.0e-10};
constexpr double cold_start_pressure_error_tolerance{1.0e-8};
constexpr double fine_grid_velocity_error_tolerance{1.0e-2};
constexpr double fine_grid_mass_flow_error_tolerance{2.0e-2};
constexpr double minimum_observed_velocity_order{1.8};

// Orthonormal basis of the rotated channel.
//
// The 3-4-5 direction avoids a special orientation such as 45 degrees while
// retaining simple exact constants for the basis components.
constexpr cfd::Vector2 streamwise_direction{0.8, 0.6};
constexpr cfd::Vector2 transverse_direction{-0.6, 0.8};

constexpr cfd::BoundaryId bottom_boundary_id{0};
constexpr cfd::BoundaryId right_boundary_id{1};
constexpr cfd::BoundaryId top_boundary_id{2};
constexpr cfd::BoundaryId left_boundary_id{3};

struct GridLevel
{
    cfd::Index streamwise_cell_count;
    cfd::Index transverse_cell_count;
};

constexpr std::array grid_levels{
    GridLevel{16, 4},
    GridLevel{32, 8},
    GridLevel{64, 16},
};

constexpr GridLevel cold_start_grid{32, 8};

struct LevelResult
{
    cfd::Index cell_count{};
    double mesh_spacing{};

    bool converged{};
    cfd::Index simple_iterations{};

    double velocity_relative_change{};
    double rhie_chow_flux_relative_residual{};
    double provisional_continuity_relative_residual{};
    double continuity_relative_residual{};
    double maximum_mass_imbalance{};

    double relative_velocity_error{};
    double normalized_pressure_error{};
    double maximum_absolute_transverse_velocity{};

    double inlet_mass_flow{};
    double outlet_mass_flow{};
    double relative_mass_flow_mismatch{};
    double relative_mass_flow_error{};

    std::optional<double> observed_velocity_order;
};

[[nodiscard]]
cfd::Point2 physical_position(const double streamwise_coordinate_value,
                              const double transverse_coordinate_value) noexcept
{
    return {
        streamwise_coordinate_value * streamwise_direction.x + transverse_coordinate_value * transverse_direction.x,
        streamwise_coordinate_value * streamwise_direction.y + transverse_coordinate_value * transverse_direction.y,
    };
}

[[nodiscard]]
double streamwise_coordinate(const cfd::Point2 &position) noexcept
{
    return position.x * streamwise_direction.x + position.y * streamwise_direction.y;
}

[[nodiscard]]
double transverse_coordinate(const cfd::Point2 &position) noexcept
{
    return position.x * transverse_direction.x + position.y * transverse_direction.y;
}

[[nodiscard]]
double pressure_gradient() noexcept
{
    return (outlet_pressure - inlet_pressure) / domain_length;
}

[[nodiscard]]
double exact_streamwise_velocity(const double transverse_position) noexcept
{
    return -pressure_gradient() * transverse_position * (domain_height - transverse_position) /
           (2.0 * dynamic_viscosity);
}

[[nodiscard]]
cfd::Vector2 exact_velocity(const cfd::Point2 &position) noexcept
{
    const double streamwise_velocity{exact_streamwise_velocity(transverse_coordinate(position))};

    return {
        streamwise_velocity * streamwise_direction.x,
        streamwise_velocity * streamwise_direction.y,
    };
}

[[nodiscard]]
double exact_pressure(const cfd::Point2 &position) noexcept
{
    return inlet_pressure + pressure_gradient() * streamwise_coordinate(position);
}

[[nodiscard]]
double exact_mass_flow() noexcept
{
    const double volume_flow{-pressure_gradient() * domain_height * domain_height * domain_height /
                             (12.0 * dynamic_viscosity)};

    return density * volume_flow;
}

[[nodiscard]]
cfd::RawMeshData make_rotated_channel_raw_mesh(const GridLevel &level)
{
    const cfd::Index streamwise_node_count{level.streamwise_cell_count + 1};

    const cfd::Index transverse_node_count{level.transverse_cell_count + 1};

    cfd::RawMeshData raw_mesh;

    raw_mesh.nodes.reserve(streamwise_node_count * transverse_node_count);

    for (cfd::Index j = 0; j < transverse_node_count; ++j)
    {
        const double n{domain_height * static_cast<double>(j) / static_cast<double>(level.transverse_cell_count)};

        for (cfd::Index i = 0; i < streamwise_node_count; ++i)
        {
            const double s{domain_length * static_cast<double>(i) / static_cast<double>(level.streamwise_cell_count)};

            raw_mesh.nodes.push_back(physical_position(s, n));
        }
    }

    const cfd::Index cell_count{level.streamwise_cell_count * level.transverse_cell_count};

    raw_mesh.cell_types.reserve(cell_count);
    raw_mesh.cell_nodes.reserve(4 * cell_count);
    raw_mesh.cell_node_offsets.reserve(cell_count + 1);
    raw_mesh.cell_node_offsets.push_back(0);

    for (cfd::Index j = 0; j < level.transverse_cell_count; ++j)
    {
        for (cfd::Index i = 0; i < level.streamwise_cell_count; ++i)
        {
            const cfd::Index lower_left{j * streamwise_node_count + i};

            raw_mesh.cell_types.push_back(cfd::CellType::Quadrilateral);

            raw_mesh.cell_nodes.insert(raw_mesh.cell_nodes.end(), {
                                                                      lower_left,
                                                                      lower_left + 1,
                                                                      lower_left + streamwise_node_count + 1,
                                                                      lower_left + streamwise_node_count,
                                                                  });

            raw_mesh.cell_node_offsets.push_back(raw_mesh.cell_nodes.size());
        }
    }

    raw_mesh.boundary_groups = {
        {bottom_boundary_id, "bottom"},
        {right_boundary_id, "right"},
        {top_boundary_id, "top"},
        {left_boundary_id, "left"},
    };

    raw_mesh.boundary_edges.reserve(2 * (level.streamwise_cell_count + level.transverse_cell_count));

    for (cfd::Index i = 0; i < level.streamwise_cell_count; ++i)
    {
        raw_mesh.boundary_edges.push_back({{i, i + 1}, bottom_boundary_id});
    }

    for (cfd::Index j = 0; j < level.transverse_cell_count; ++j)
    {
        const cfd::Index lower_right{j * streamwise_node_count + level.streamwise_cell_count};

        raw_mesh.boundary_edges.push_back({
            {
                lower_right,
                lower_right + streamwise_node_count,
            },
            right_boundary_id,
        });
    }

    const cfd::Index top_row{level.transverse_cell_count * streamwise_node_count};

    for (cfd::Index i = level.streamwise_cell_count; i > 0; --i)
    {
        raw_mesh.boundary_edges.push_back({
            {
                top_row + i,
                top_row + i - 1,
            },
            top_boundary_id,
        });
    }

    for (cfd::Index j = level.transverse_cell_count; j > 0; --j)
    {
        raw_mesh.boundary_edges.push_back({
            {
                j * streamwise_node_count,
                (j - 1) * streamwise_node_count,
            },
            left_boundary_id,
        });
    }

    return raw_mesh;
}

[[nodiscard]]
cfd::ScalarBoundaryConditions velocity_boundary_conditions(const cfd::Mesh &mesh)
{
    return {
        mesh.boundary_groups().size(),
        {
            {cfd::ScalarBoundaryConditionType::Dirichlet, 0.0},
            {cfd::ScalarBoundaryConditionType::Neumann, 0.0},
            {cfd::ScalarBoundaryConditionType::Dirichlet, 0.0},
            {cfd::ScalarBoundaryConditionType::Neumann, 0.0},
        },
    };
}

[[nodiscard]]
cfd::ScalarBoundaryConditions pressure_boundary_conditions(const cfd::Mesh &mesh)
{
    return {
        mesh.boundary_groups().size(),
        {
            {cfd::ScalarBoundaryConditionType::Neumann, 0.0},
            {cfd::ScalarBoundaryConditionType::Dirichlet, outlet_pressure},
            {cfd::ScalarBoundaryConditionType::Neumann, 0.0},
            {cfd::ScalarBoundaryConditionType::Dirichlet, inlet_pressure},
        },
    };
}

[[nodiscard]]
cfd::PressureCorrectionBoundaryConditions pressure_correction_boundary_conditions(const cfd::Mesh &mesh)
{
    return {
        mesh.boundary_groups().size(),
        {
            cfd::PressureCorrectionBoundaryConditionType::FixedMassFlux,
            cfd::PressureCorrectionBoundaryConditionType::FixedPressure,
            cfd::PressureCorrectionBoundaryConditionType::FixedMassFlux,
            cfd::PressureCorrectionBoundaryConditionType::FixedPressure,
        },
    };
}

[[nodiscard]]
cfd::IncompressibleSimpleOptions simple_options()
{
    return {
        .maximum_iterations = 300,
        .momentum_relaxation_factor = 1.0,
        .pressure_relaxation_factor = 0.3,
        .rhie_chow_flux_relaxation_factor = 1.0,
        .velocity_relative_tolerance = 1.0e-10,
        .rhie_chow_flux_relative_tolerance = 1.0e-10,
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
cfd::IncompressibleSimpleOptions cold_start_options()
{
    cfd::IncompressibleSimpleOptions options{simple_options()};

    options.maximum_iterations = 2000;
    options.pressure_relaxation_factor = 0.1;
    options.rhie_chow_flux_relaxation_factor = 0.3;

    return options;
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
double relative_mass_flow_error(const double inlet_mass_flow, const double outlet_mass_flow)
{
    const double reference_mass_flow{exact_mass_flow()};

    if (!(reference_mass_flow > 0.0) || !std::isfinite(reference_mass_flow))
    {
        throw std::runtime_error("Rotated Poiseuille reference mass flow is invalid.");
    }

    return std::max(std::abs(inlet_mass_flow - reference_mass_flow), std::abs(outlet_mass_flow - reference_mass_flow)) /
           reference_mass_flow;
}

[[nodiscard]]
double relative_vector_error(const ErrorStatistics &x_error, const ErrorStatistics &y_error,
                             const ErrorStatistics &x_reference, const ErrorStatistics &y_reference)
{
    const double squared_error{x_error.area_weighted_squared_error + y_error.area_weighted_squared_error};

    const double squared_reference{x_reference.area_weighted_squared_error + y_reference.area_weighted_squared_error};

    if (!std::isfinite(squared_error) || squared_error < 0.0 || !std::isfinite(squared_reference) ||
        !(squared_reference > 0.0))
    {
        throw std::runtime_error("Rotated Poiseuille velocity error normalization is invalid.");
    }

    return std::sqrt(squared_error / squared_reference);
}

[[nodiscard]]
LevelResult solve_level(const GridLevel &level, const cfd::IncompressibleSimpleOptions &options,
                        const bool initialize_with_exact_pressure)
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_rotated_channel_raw_mesh(level))};

    const cfd::Mesh &mesh{build_result.mesh};

    cfd::CellVelocityField velocity{mesh.cell_count()};

    cfd::CellScalarField pressure{mesh.cell_count()};

    cfd::FaceFluxField mass_flux{mesh.face_count()};

    if (initialize_with_exact_pressure)
    {
        for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
        {
            pressure[cell_id] = exact_pressure(mesh.cell_centers()[cell_id]);
        }
    }

    const cfd::ScalarBoundaryConditions velocity_conditions{velocity_boundary_conditions(mesh)};

    const cfd::ScalarBoundaryConditions pressure_conditions{pressure_boundary_conditions(mesh)};

    const cfd::PressureCorrectionBoundaryConditions pressure_correction_conditions{
        pressure_correction_boundary_conditions(mesh)};

    cfd::IncompressibleSimpleSolver solver{
        mesh, density, dynamic_viscosity, cfd::ScalarConvectionScheme::Linear, options,
    };

    const cfd::IncompressibleSimpleResult result{solver.solve(velocity_conditions, velocity_conditions,
                                                              pressure_conditions, pressure_correction_conditions,
                                                              velocity, pressure, mass_flux)};

    ErrorAccumulator x_velocity_error;
    ErrorAccumulator y_velocity_error;
    ErrorAccumulator x_velocity_reference;
    ErrorAccumulator y_velocity_reference;
    ErrorAccumulator pressure_error;
    ErrorAccumulator pressure_reference;

    double maximum_absolute_transverse_velocity{};

    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        const double area{mesh.cell_areas()[cell_id]};

        const cfd::Point2 &center{mesh.cell_centers()[cell_id]};

        const cfd::Vector2 exact_velocity_value{exact_velocity(center)};

        const double exact_pressure_value{exact_pressure(center)};

        x_velocity_error.add(area, velocity.u()[cell_id] - exact_velocity_value.x);

        y_velocity_error.add(area, velocity.v()[cell_id] - exact_velocity_value.y);

        x_velocity_reference.add(area, exact_velocity_value.x);

        y_velocity_reference.add(area, exact_velocity_value.y);

        pressure_error.add(area, pressure[cell_id] - exact_pressure_value);

        pressure_reference.add(area, exact_pressure_value);

        const double transverse_velocity{velocity.u()[cell_id] * transverse_direction.x +
                                         velocity.v()[cell_id] * transverse_direction.y};

        maximum_absolute_transverse_velocity =
            std::max(maximum_absolute_transverse_velocity, std::abs(transverse_velocity));
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
                throw std::runtime_error("Rotated Poiseuille wall FixedMassFlux value changed from zero.");
            }
            break;

        default:
            throw std::runtime_error("Rotated Poiseuille mesh contains an unexpected boundary ID.");
        }
    }

    const ErrorStatistics x_velocity_error_statistics{x_velocity_error.finish()};

    const ErrorStatistics y_velocity_error_statistics{y_velocity_error.finish()};

    const ErrorStatistics x_velocity_reference_statistics{x_velocity_reference.finish()};

    const ErrorStatistics y_velocity_reference_statistics{y_velocity_reference.finish()};

    const ErrorStatistics pressure_error_statistics{pressure_error.finish()};

    const ErrorStatistics pressure_reference_statistics{pressure_reference.finish()};

    if (!(pressure_reference_statistics.area_weighted_rms_error > 0.0))
    {
        throw std::runtime_error("Rotated Poiseuille pressure reference norm is invalid.");
    }

    return {
        .cell_count = mesh.cell_count(),
        .mesh_spacing = domain_height / static_cast<double>(level.transverse_cell_count),

        .converged = result.converged,
        .simple_iterations = result.iteration_count,

        .velocity_relative_change = result.velocity_relative_change,

        .rhie_chow_flux_relative_residual = result.rhie_chow_flux_relative_residual,

        .provisional_continuity_relative_residual = result.provisional_continuity_relative_residual,

        .continuity_relative_residual = result.continuity_relative_residual,

        .maximum_mass_imbalance = result.maximum_mass_imbalance,

        .relative_velocity_error =
            relative_vector_error(x_velocity_error_statistics, y_velocity_error_statistics,
                                  x_velocity_reference_statistics, y_velocity_reference_statistics),

        .normalized_pressure_error =
            pressure_error_statistics.area_weighted_rms_error / pressure_reference_statistics.area_weighted_rms_error,

        .maximum_absolute_transverse_velocity = maximum_absolute_transverse_velocity,

        .inlet_mass_flow = inlet_mass_flow,

        .outlet_mass_flow = outlet_mass_flow,

        .relative_mass_flow_mismatch = relative_mass_flow_mismatch(inlet_mass_flow, outlet_mass_flow),

        .relative_mass_flow_error = relative_mass_flow_error(inlet_mass_flow, outlet_mass_flow),

        .observed_velocity_order = std::nullopt,
    };
}

[[nodiscard]]
LevelResult solve_level(const GridLevel &level)
{
    return solve_level(level, simple_options(), true);
}

void validate_convergence_diagnostics(const LevelResult &result, const cfd::IncompressibleSimpleOptions &options,
                                      const char *context)
{
    if (!result.converged)
    {
        throw std::runtime_error(std::string{context} + " SIMPLE solve did not converge.");
    }

    if (result.velocity_relative_change > options.velocity_relative_tolerance)
    {
        throw std::runtime_error(std::string{context} + " velocity relative change exceeds its tolerance.");
    }

    if (result.rhie_chow_flux_relative_residual > options.rhie_chow_flux_relative_tolerance)
    {
        throw std::runtime_error(std::string{context} + " Rhie-Chow flux residual exceeds its tolerance.");
    }

    if (result.provisional_continuity_relative_residual > options.continuity_relative_tolerance)
    {
        throw std::runtime_error(std::string{context} + " provisional continuity residual exceeds its tolerance.");
    }

    if (result.continuity_relative_residual > corrected_continuity_tolerance)
    {
        throw std::runtime_error(std::string{context} + " corrected continuity residual is too large.");
    }

    if (result.maximum_mass_imbalance > maximum_mass_imbalance_tolerance)
    {
        throw std::runtime_error(std::string{context} + " maximum mass imbalance is too large.");
    }

    if (result.relative_mass_flow_mismatch > mass_flow_mismatch_tolerance)
    {
        throw std::runtime_error(std::string{context} + " inlet/outlet mass-flow mismatch is too large.");
    }
}

void validate_cold_start_result(const LevelResult &result, const cfd::IncompressibleSimpleOptions &options)
{
    validate_convergence_diagnostics(result, options, "Cold-start rotated Poiseuille");

    if (result.maximum_absolute_transverse_velocity > maximum_transverse_velocity_tolerance)
    {
        throw std::runtime_error("Cold-start rotated Poiseuille transverse velocity is too large.");
    }

    if (result.normalized_pressure_error > cold_start_pressure_error_tolerance)
    {
        throw std::runtime_error("Cold-start rotated Poiseuille pressure error is too large.");
    }
}

void validate_results(const std::vector<LevelResult> &results)
{
    if (results.empty())
    {
        throw std::runtime_error("Rotated Poiseuille verification has no grid levels.");
    }

    const cfd::IncompressibleSimpleOptions options{simple_options()};

    for (std::size_t index = 1; index < results.size(); ++index)
    {
        if (!(results[index].relative_velocity_error < results[index - 1].relative_velocity_error))
        {
            throw std::runtime_error("Rotated Poiseuille relative velocity error did not decrease under refinement.");
        }

        if (!(results[index].relative_mass_flow_error < results[index - 1].relative_mass_flow_error))
        {
            throw std::runtime_error("Rotated Poiseuille mass-flow error did not decrease under refinement.");
        }

        const std::optional<double> &observed_order{results[index].observed_velocity_order};

        if (!observed_order.has_value())
        {
            throw std::runtime_error("Rotated Poiseuille observed velocity order is unavailable.");
        }

        if (observed_order.value() < minimum_observed_velocity_order)
        {
            throw std::runtime_error("Rotated Poiseuille observed velocity order is too low.");
        }
    }

    for (const LevelResult &result : results)
    {
        validate_convergence_diagnostics(result, options, "Rotated Poiseuille");

        if (result.maximum_absolute_transverse_velocity > maximum_transverse_velocity_tolerance)
        {
            throw std::runtime_error("Rotated Poiseuille transverse velocity is too large.");
        }
    }

    const LevelResult &finest{results.back()};

    if (finest.relative_velocity_error > fine_grid_velocity_error_tolerance)
    {
        throw std::runtime_error("Rotated Poiseuille fine-grid relative velocity error is excessive.");
    }

    if (finest.relative_mass_flow_error > fine_grid_mass_flow_error_tolerance)
    {
        throw std::runtime_error("Rotated Poiseuille fine-grid mass-flow error is excessive.");
    }
}

void print_results(const std::vector<LevelResult> &results)
{
    const double mean_velocity{domain_height * domain_height * (-pressure_gradient()) / (12.0 * dynamic_viscosity)};

    const double maximum_velocity{domain_height * domain_height * (-pressure_gradient()) / (8.0 * dynamic_viscosity)};

    std::cout << "Steady incompressible SIMPLE v1 - rotated plane Poiseuille flow - orthogonal QUAD\n"
              << "L=" << domain_length << ", H=" << domain_height << ", rho=" << density << ", mu=" << dynamic_viscosity
              << ", p_in=" << inlet_pressure << ", p_out=" << outlet_pressure << ", e_s=(" << streamwise_direction.x
              << ", " << streamwise_direction.y << ')' << ", U_mean=" << mean_velocity << ", U_max=" << maximum_velocity
              << "\n\n"

              << std::left << std::setw(9) << "cells" << std::setw(8) << "iters" << std::setw(16) << "rel_L2(U)"
              << std::setw(10) << "order" << std::setw(16) << "max|Un|" << std::setw(16) << "r_U" << std::setw(16)
              << "r_RC" << std::setw(16) << "r_cont(prov)" << std::setw(16) << "r_cont(corr)" << std::setw(16)
              << "max|Rm|" << std::setw(16) << "flow_mismatch" << std::setw(16) << "flow_error" << std::setw(16)
              << "p_error" << '\n';

    for (const LevelResult &result : results)
    {
        std::cout << std::left << std::setw(9) << result.cell_count << std::setw(8) << result.simple_iterations
                  << std::scientific << std::setprecision(5) << std::setw(16) << result.relative_velocity_error;

        if (result.observed_velocity_order.has_value())
        {
            std::cout << std::fixed << std::setprecision(3) << std::setw(10) << *result.observed_velocity_order;
        }
        else
        {
            std::cout << std::setw(10) << "-";
        }

        std::cout << std::scientific << std::setprecision(5) << std::setw(16)
                  << result.maximum_absolute_transverse_velocity << std::setw(16) << result.velocity_relative_change
                  << std::setw(16) << result.rhie_chow_flux_relative_residual << std::setw(16)
                  << result.provisional_continuity_relative_residual << std::setw(16)
                  << result.continuity_relative_residual << std::setw(16) << result.maximum_mass_imbalance
                  << std::setw(16) << result.relative_mass_flow_mismatch << std::setw(16)
                  << result.relative_mass_flow_error << std::setw(16) << result.normalized_pressure_error << '\n';
    }

    const LevelResult &finest{results.back()};

    std::cout << "\nExact mass flow:              " << std::scientific << exact_mass_flow()
              << "\nFinest-grid inlet mass flow:  " << finest.inlet_mass_flow
              << "\nFinest-grid outlet mass flow: " << finest.outlet_mass_flow
              << "\nFinal velocity relative change: " << finest.velocity_relative_change << '\n';
}

void print_cold_start_result(const LevelResult &result)
{
    std::cout << "\nCold-start coupled rotated solve (" << cold_start_grid.streamwise_cell_count << 'x'
              << cold_start_grid.transverse_cell_count << ", " << result.cell_count << " cells)\n"
              << std::scientific << std::setprecision(6) << "  iterations: " << result.simple_iterations
              << ", r_U: " << result.velocity_relative_change << ", r_RC: " << result.rhie_chow_flux_relative_residual
              << ", r_cont(provisional): " << result.provisional_continuity_relative_residual
              << ", r_cont(corrected): " << result.continuity_relative_residual
              << "\n  rel_L2(U): " << result.relative_velocity_error
              << ", p_error: " << result.normalized_pressure_error
              << ", max|Un|: " << result.maximum_absolute_transverse_velocity
              << ", max|Rm|: " << result.maximum_mass_imbalance
              << ", flow mismatch: " << result.relative_mass_flow_mismatch
              << ", flow error: " << result.relative_mass_flow_error << '\n';
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
            results[index].observed_velocity_order =
                observed_order(results[index - 1].relative_velocity_error, results[index].relative_velocity_error,
                               results[index - 1].mesh_spacing, results[index].mesh_spacing);
        }

        validate_results(results);
        print_results(results);

        const cfd::IncompressibleSimpleOptions options{cold_start_options()};

        const LevelResult cold_start_result{solve_level(cold_start_grid, options, false)};

        print_cold_start_result(cold_start_result);

        validate_cold_start_result(cold_start_result, options);

        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "Rotated incompressible SIMPLE Poiseuille verification failed: " << error.what() << '\n';

        return 1;
    }
}
