#include "cfd/numerics/IncompressibleSimpleSolver.hpp"

#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVelocityField.hpp"
#include "cfd/field/FaceFluxField.hpp"
#include "cfd/field/PressureCorrectionBoundaryConditions.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/math/Vector2.hpp"
#include "cfd/mesh/Boundary.hpp"
#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/mesh/Types.hpp"
#include "cfd/meshing/RawMeshData.hpp"

#include "support/TestUtils.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace
{

using cfd::test::require;
using cfd::test::require_throws;
using cfd::test::require_throws_with_message;

constexpr cfd::BoundaryId bottom_boundary_id{0};
constexpr cfd::BoundaryId right_boundary_id{1};
constexpr cfd::BoundaryId top_boundary_id{2};
constexpr cfd::BoundaryId left_boundary_id{3};

static_assert(!std::is_copy_constructible_v<cfd::IncompressibleSimpleSolver>);
static_assert(!std::is_copy_assignable_v<cfd::IncompressibleSimpleSolver>);
static_assert(!std::is_move_constructible_v<cfd::IncompressibleSimpleSolver>);
static_assert(!std::is_move_assignable_v<cfd::IncompressibleSimpleSolver>);
static_assert(std::is_trivially_copyable_v<cfd::SimpleTimingBreakdown>);

[[nodiscard]]
cfd::RawMeshData make_channel_raw_mesh(const cfd::Index x_cell_count, const cfd::Index y_cell_count,
                                       const double length = 2.0, const double height = 1.0)
{
    cfd::RawMeshData raw_mesh;
    const cfd::Index x_node_count{x_cell_count + 1};
    const cfd::Index y_node_count{y_cell_count + 1};
    raw_mesh.nodes.reserve(x_node_count * y_node_count);
    for (cfd::Index j = 0; j < y_node_count; ++j)
    {
        for (cfd::Index i = 0; i < x_node_count; ++i)
        {
            raw_mesh.nodes.push_back({length * static_cast<double>(i) / static_cast<double>(x_cell_count),
                                      height * static_cast<double>(j) / static_cast<double>(y_cell_count)});
        }
    }

    raw_mesh.cell_types.reserve(x_cell_count * y_cell_count);
    raw_mesh.cell_nodes.reserve(4 * x_cell_count * y_cell_count);
    raw_mesh.cell_node_offsets.reserve(x_cell_count * y_cell_count + 1);
    raw_mesh.cell_node_offsets.push_back(0);
    for (cfd::Index j = 0; j < y_cell_count; ++j)
    {
        for (cfd::Index i = 0; i < x_cell_count; ++i)
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
    for (cfd::Index i = 0; i < x_cell_count; ++i)
    {
        raw_mesh.boundary_edges.push_back({{i, i + 1}, bottom_boundary_id});
    }
    for (cfd::Index j = 0; j < y_cell_count; ++j)
    {
        const cfd::Index lower_right{j * x_node_count + x_cell_count};
        raw_mesh.boundary_edges.push_back({{lower_right, lower_right + x_node_count}, right_boundary_id});
    }
    const cfd::Index top_row{y_cell_count * x_node_count};
    for (cfd::Index i = x_cell_count; i > 0; --i)
    {
        raw_mesh.boundary_edges.push_back({{top_row + i, top_row + i - 1}, top_boundary_id});
    }
    for (cfd::Index j = y_cell_count; j > 0; --j)
    {
        raw_mesh.boundary_edges.push_back({{j * x_node_count, (j - 1) * x_node_count}, left_boundary_id});
    }
    return raw_mesh;
}

[[nodiscard]]
cfd::RawMeshData make_disconnected_two_cell_raw_mesh()
{
    constexpr cfd::BoundaryId wall_boundary_id{0};

    cfd::RawMeshData raw_mesh;
    raw_mesh.nodes = {
        {0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}, {2.0, 0.0}, {3.0, 0.0}, {3.0, 1.0}, {2.0, 1.0},
    };
    raw_mesh.cell_types = {
        cfd::CellType::Quadrilateral,
        cfd::CellType::Quadrilateral,
    };
    raw_mesh.cell_nodes = {
        0, 1, 2, 3, 4, 5, 6, 7,
    };
    raw_mesh.cell_node_offsets = {0, 4, 8};
    raw_mesh.boundary_groups = {{wall_boundary_id, "wall"}};
    raw_mesh.boundary_edges = {
        {{0, 1}, wall_boundary_id}, {{1, 2}, wall_boundary_id}, {{2, 3}, wall_boundary_id}, {{3, 0}, wall_boundary_id},
        {{4, 5}, wall_boundary_id}, {{5, 6}, wall_boundary_id}, {{6, 7}, wall_boundary_id}, {{7, 4}, wall_boundary_id},
    };
    return raw_mesh;
}

[[nodiscard]]
cfd::ScalarBoundaryConditions uniform_scalar_conditions(const cfd::Mesh &mesh,
                                                        const cfd::ScalarBoundaryConditionType type,
                                                        const double value = 0.0)
{
    return {mesh.boundary_groups().size(),
            std::vector<cfd::ScalarBoundaryCondition>(mesh.boundary_groups().size(), {type, value})};
}

[[nodiscard]]
cfd::PressureCorrectionBoundaryConditions uniform_pressure_correction_conditions(
    const cfd::Mesh &mesh, const cfd::PressureCorrectionBoundaryConditionType type)
{
    return {mesh.boundary_groups().size(),
            std::vector<cfd::PressureCorrectionBoundaryConditionType>(mesh.boundary_groups().size(), type)};
}

[[nodiscard]]
cfd::ScalarBoundaryConditions channel_velocity_conditions(const cfd::Mesh &mesh)
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
cfd::ScalarBoundaryConditions channel_pressure_conditions(const cfd::Mesh &mesh, const double inlet_pressure,
                                                          const double outlet_pressure)
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
cfd::PressureCorrectionBoundaryConditions channel_pressure_correction_conditions(const cfd::Mesh &mesh)
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
cfd::IncompressibleSimpleOptions test_options()
{
    return {
        .maximum_iterations = 200,
        .momentum_relaxation_factor = 1.0,
        .pressure_relaxation_factor = 0.3,
        .velocity_relative_tolerance = 1.0e-9,
        .continuity_relative_tolerance = 1.0e-10,
        .momentum_linear_solver = {.relative_tolerance = 1.0e-12, .maximum_iterations = 1000},
        .pressure_correction_linear_solver = {.relative_tolerance = 1.0e-12, .maximum_iterations = 1000},
    };
}

void test_constructor_and_options_validation()
{
    require(cfd::IncompressibleSimpleOptions{}.rhie_chow_flux_relaxation_factor == 1.0,
            "SIMPLE changed the default unrelaxed Rhie-Chow flux path.");
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_channel_raw_mesh(2, 2))};
    const cfd::Mesh &mesh{build_result.mesh};
    require_throws<std::invalid_argument>(
        [&]() { cfd::IncompressibleSimpleSolver solver{mesh, 0.0, 0.1, cfd::ScalarConvectionScheme::Linear}; },
        "SIMPLE accepted zero density.");
    require_throws<std::invalid_argument>(
        [&]() { cfd::IncompressibleSimpleSolver solver{mesh, 1.0, 0.0, cfd::ScalarConvectionScheme::Linear}; },
        "SIMPLE accepted zero dynamic viscosity.");

    cfd::IncompressibleSimpleOptions options{test_options()};
    options.maximum_iterations = 0;
    require_throws<std::invalid_argument>(
        [&]() { cfd::IncompressibleSimpleSolver solver{mesh, 1.0, 0.1, cfd::ScalarConvectionScheme::Linear, options}; },
        "SIMPLE accepted zero outer iterations.");
    options = test_options();
    options.momentum_relaxation_factor = 0.7;
    require_throws<std::invalid_argument>(
        [&]() { cfd::IncompressibleSimpleSolver solver{mesh, 1.0, 0.1, cfd::ScalarConvectionScheme::Linear, options}; },
        "SIMPLE v1 accepted momentum under-relaxation.");

    for (const double invalid_factor :
         {0.0, -1.0, 1.1, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()})
    {
        options = test_options();
        options.pressure_relaxation_factor = invalid_factor;
        require_throws<std::invalid_argument>(
            [&]() {
                cfd::IncompressibleSimpleSolver solver{mesh, 1.0, 0.1, cfd::ScalarConvectionScheme::Linear, options};
            },
            "SIMPLE accepted an invalid pressure relaxation factor.");

        options = test_options();
        options.rhie_chow_flux_relaxation_factor = invalid_factor;
        require_throws<std::invalid_argument>(
            [&]() {
                cfd::IncompressibleSimpleSolver solver{mesh, 1.0, 0.1, cfd::ScalarConvectionScheme::Linear, options};
            },
            "SIMPLE accepted an invalid Rhie-Chow flux relaxation factor.");
    }
    for (const double invalid_tolerance :
         {0.0, -1.0, 1.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()})
    {
        options = test_options();
        options.velocity_relative_tolerance = invalid_tolerance;
        require_throws<std::invalid_argument>(
            [&]() {
                cfd::IncompressibleSimpleSolver solver{mesh, 1.0, 0.1, cfd::ScalarConvectionScheme::Linear, options};
            },
            "SIMPLE accepted an invalid velocity tolerance.");
        options = test_options();
        options.continuity_relative_tolerance = invalid_tolerance;
        require_throws<std::invalid_argument>(
            [&]() {
                cfd::IncompressibleSimpleSolver solver{mesh, 1.0, 0.1, cfd::ScalarConvectionScheme::Linear, options};
            },
            "SIMPLE accepted an invalid continuity tolerance.");
        options = test_options();
        options.rhie_chow_flux_relative_tolerance = invalid_tolerance;
        require_throws<std::invalid_argument>(
            [&]() {
                cfd::IncompressibleSimpleSolver solver{mesh, 1.0, 0.1, cfd::ScalarConvectionScheme::Linear, options};
            },
            "SIMPLE accepted an invalid Rhie-Chow flux relative tolerance.");
    }
}

void test_rejects_disconnected_cell_domain()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_disconnected_two_cell_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};

    require(mesh.cell_count() == 2, "Disconnected SIMPLE fixture did not build both valid cells.");
    require_throws_with_message<std::invalid_argument>(
        [&mesh]() {
            const cfd::IncompressibleSimpleSolver solver{mesh, 1.0, 0.1, cfd::ScalarConvectionScheme::Linear,
                                                         test_options()};
        },
        "single connected cell domain", "SIMPLE accepted a disconnected cell domain.");
}

void test_pressure_boundary_semantic_consistency()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_channel_raw_mesh(2, 2))};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarBoundaryConditions velocity_conditions{
        uniform_scalar_conditions(mesh, cfd::ScalarBoundaryConditionType::Dirichlet)};
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    cfd::FaceFluxField mass_flux{mesh.face_count()};

    {
        const cfd::ScalarBoundaryConditions pressure_conditions{
            uniform_scalar_conditions(mesh, cfd::ScalarBoundaryConditionType::Neumann)};
        const cfd::PressureCorrectionBoundaryConditions pressure_correction_conditions{
            uniform_pressure_correction_conditions(mesh, cfd::PressureCorrectionBoundaryConditionType::FixedPressure)};
        cfd::IncompressibleSimpleSolver solver{mesh, 1.0, 0.1, cfd::ScalarConvectionScheme::Linear, test_options()};
        require_throws<std::invalid_argument>(
            [&]() {
                static_cast<void>(solver.solve(velocity_conditions, velocity_conditions, pressure_conditions,
                                               pressure_correction_conditions, velocity, pressure, mass_flux));
            },
            "SIMPLE accepted FixedPressure with physical pressure Neumann data.");
    }
    {
        const cfd::ScalarBoundaryConditions pressure_conditions{
            uniform_scalar_conditions(mesh, cfd::ScalarBoundaryConditionType::Dirichlet)};
        const cfd::PressureCorrectionBoundaryConditions pressure_correction_conditions{
            uniform_pressure_correction_conditions(mesh, cfd::PressureCorrectionBoundaryConditionType::FixedMassFlux)};
        cfd::IncompressibleSimpleSolver solver{mesh, 1.0, 0.1, cfd::ScalarConvectionScheme::Linear, test_options()};
        require_throws<std::invalid_argument>(
            [&]() {
                static_cast<void>(solver.solve(velocity_conditions, velocity_conditions, pressure_conditions,
                                               pressure_correction_conditions, velocity, pressure, mass_flux));
            },
            "SIMPLE accepted FixedMassFlux with physical pressure Dirichlet data.");
    }
}

void test_all_fixed_mass_flux_compatibility_and_gauge_path()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_channel_raw_mesh(2, 2))};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarBoundaryConditions velocity_conditions{
        uniform_scalar_conditions(mesh, cfd::ScalarBoundaryConditionType::Dirichlet)};
    const cfd::ScalarBoundaryConditions pressure_conditions{
        uniform_scalar_conditions(mesh, cfd::ScalarBoundaryConditionType::Neumann)};
    const cfd::PressureCorrectionBoundaryConditions pressure_correction_conditions{
        uniform_pressure_correction_conditions(mesh, cfd::PressureCorrectionBoundaryConditionType::FixedMassFlux)};

    {
        cfd::CellVelocityField velocity{mesh.cell_count()};
        cfd::CellScalarField pressure{mesh.cell_count()};
        cfd::FaceFluxField mass_flux{mesh.face_count()};
        mass_flux[0] = 1.0;
        cfd::IncompressibleSimpleSolver solver{mesh, 1.0, 0.1, cfd::ScalarConvectionScheme::Linear, test_options()};
        require_throws<std::runtime_error>(
            [&]() {
                static_cast<void>(solver.solve(velocity_conditions, velocity_conditions, pressure_conditions,
                                               pressure_correction_conditions, velocity, pressure, mass_flux));
            },
            "SIMPLE accepted globally incompatible all-FixedMassFlux boundary flow.");
    }
    {
        cfd::CellVelocityField velocity{mesh.cell_count()};
        cfd::CellScalarField pressure{mesh.cell_count()};
        cfd::FaceFluxField mass_flux{mesh.face_count()};
        cfd::IncompressibleSimpleOptions options{test_options()};
        options.maximum_iterations = 2;
        cfd::IncompressibleSimpleSolver solver{mesh, 1.0, 0.1, cfd::ScalarConvectionScheme::Linear, options};
        const cfd::IncompressibleSimpleResult result{solver.solve(velocity_conditions, velocity_conditions,
                                                                  pressure_conditions, pressure_correction_conditions,
                                                                  velocity, pressure, mass_flux)};
        require(result.converged, "Compatible all-FixedMassFlux zero-flow case did not converge through gauge fixing.");
        require(result.maximum_mass_imbalance == 0.0, "Gauge-fixed zero-flow case created mass imbalance.");
    }
}

void test_fixed_pressure_anchored_zero_flow_path()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_channel_raw_mesh(2, 2))};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarBoundaryConditions velocity_conditions{
        uniform_scalar_conditions(mesh, cfd::ScalarBoundaryConditionType::Dirichlet)};
    std::vector<cfd::ScalarBoundaryCondition> pressure_values(mesh.boundary_groups().size(),
                                                              {cfd::ScalarBoundaryConditionType::Neumann, 0.0});
    pressure_values.at(left_boundary_id) = {cfd::ScalarBoundaryConditionType::Dirichlet, 0.0};
    const cfd::ScalarBoundaryConditions pressure_conditions{mesh.boundary_groups().size(), std::move(pressure_values)};
    std::vector<cfd::PressureCorrectionBoundaryConditionType> correction_values(
        mesh.boundary_groups().size(), cfd::PressureCorrectionBoundaryConditionType::FixedMassFlux);
    correction_values.at(left_boundary_id) = cfd::PressureCorrectionBoundaryConditionType::FixedPressure;
    const cfd::PressureCorrectionBoundaryConditions pressure_correction_conditions{mesh.boundary_groups().size(),
                                                                                   std::move(correction_values)};
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    cfd::IncompressibleSimpleOptions options{test_options()};
    options.maximum_iterations = 2;
    cfd::IncompressibleSimpleSolver solver{mesh, 1.0, 0.1, cfd::ScalarConvectionScheme::Linear, options};

    const cfd::IncompressibleSimpleResult result{solver.solve(velocity_conditions, velocity_conditions,
                                                              pressure_conditions, pressure_correction_conditions,
                                                              velocity, pressure, mass_flux)};

    require(result.converged, "FixedPressure-anchored zero-flow case did not converge without an artificial gauge.");
    require(result.maximum_mass_imbalance == 0.0, "FixedPressure-anchored zero flow created mass imbalance.");
}

void test_input_validation_precedes_iterations()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_channel_raw_mesh(2, 2))};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarBoundaryConditions velocity_conditions{channel_velocity_conditions(mesh)};
    const cfd::ScalarBoundaryConditions pressure_conditions{channel_pressure_conditions(mesh, 0.02, 0.0)};
    const cfd::PressureCorrectionBoundaryConditions pressure_correction_conditions{
        channel_pressure_correction_conditions(mesh)};
    cfd::IncompressibleSimpleSolver solver{mesh, 1.0, 0.1, cfd::ScalarConvectionScheme::Linear, test_options()};
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    cfd::CellVelocityField wrong_velocity{mesh.cell_count() - 1};
    require_throws<std::invalid_argument>(
        [&]() {
            static_cast<void>(solver.solve(velocity_conditions, velocity_conditions, pressure_conditions,
                                           pressure_correction_conditions, wrong_velocity, pressure, mass_flux));
        },
        "SIMPLE accepted the wrong velocity cardinality.");

    cfd::CellScalarField wrong_pressure{mesh.cell_count() - 1};
    require_throws<std::invalid_argument>(
        [&]() {
            static_cast<void>(solver.solve(velocity_conditions, velocity_conditions, pressure_conditions,
                                           pressure_correction_conditions, velocity, wrong_pressure, mass_flux));
        },
        "SIMPLE accepted the wrong pressure cardinality.");
    cfd::FaceFluxField wrong_mass_flux{mesh.face_count() - 1};
    require_throws<std::invalid_argument>(
        [&]() {
            static_cast<void>(solver.solve(velocity_conditions, velocity_conditions, pressure_conditions,
                                           pressure_correction_conditions, velocity, pressure, wrong_mass_flux));
        },
        "SIMPLE accepted the wrong mass-flux cardinality.");

    const cfd::ScalarBoundaryConditions wrong_boundary_conditions{1,
                                                                  {{cfd::ScalarBoundaryConditionType::Dirichlet, 0.0}}};
    require_throws<std::invalid_argument>(
        [&]() {
            static_cast<void>(solver.solve(wrong_boundary_conditions, velocity_conditions, pressure_conditions,
                                           pressure_correction_conditions, velocity, pressure, mass_flux));
        },
        "SIMPLE accepted the wrong boundary-condition cardinality.");

    velocity.u()[1] = std::numeric_limits<double>::quiet_NaN();
    require_throws<std::invalid_argument>(
        [&]() {
            static_cast<void>(solver.solve(velocity_conditions, velocity_conditions, pressure_conditions,
                                           pressure_correction_conditions, velocity, pressure, mass_flux));
        },
        "SIMPLE accepted non-finite initial velocity.");
    velocity.u()[1] = 0.0;
    pressure[1] = std::numeric_limits<double>::infinity();
    require_throws<std::invalid_argument>(
        [&]() {
            static_cast<void>(solver.solve(velocity_conditions, velocity_conditions, pressure_conditions,
                                           pressure_correction_conditions, velocity, pressure, mass_flux));
        },
        "SIMPLE accepted non-finite initial pressure.");
    pressure[1] = 0.0;
    mass_flux[1] = std::numeric_limits<double>::quiet_NaN();
    require_throws<std::invalid_argument>(
        [&]() {
            static_cast<void>(solver.solve(velocity_conditions, velocity_conditions, pressure_conditions,
                                           pressure_correction_conditions, velocity, pressure, mass_flux));
        },
        "SIMPLE accepted non-finite initial mass flux.");
}

[[nodiscard]]
cfd::IncompressibleSimpleResult solve_pressure_driven_channel(
    const cfd::Mesh &mesh, cfd::IncompressibleSimpleOptions options, cfd::CellVelocityField &velocity,
    cfd::CellScalarField &pressure, cfd::FaceFluxField &mass_flux, const bool initialize_exact_pressure = true,
    const cfd::SimpleIterationCallback &iteration_callback = {})
{
    constexpr double length{2.0};
    constexpr double inlet_pressure{0.02};
    constexpr double outlet_pressure{0.0};
    if (initialize_exact_pressure)
    {
        for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
        {
            pressure[cell_id] =
                inlet_pressure + (outlet_pressure - inlet_pressure) * mesh.cell_centers()[cell_id].x / length;
        }
    }
    const cfd::ScalarBoundaryConditions velocity_conditions{channel_velocity_conditions(mesh)};
    const cfd::ScalarBoundaryConditions pressure_conditions{
        channel_pressure_conditions(mesh, inlet_pressure, outlet_pressure)};
    const cfd::PressureCorrectionBoundaryConditions pressure_correction_conditions{
        channel_pressure_correction_conditions(mesh)};
    cfd::IncompressibleSimpleSolver solver{mesh, 1.0, 0.1, cfd::ScalarConvectionScheme::Linear, options};
    return solver.solve(velocity_conditions, velocity_conditions, pressure_conditions, pressure_correction_conditions,
                        velocity, pressure, mass_flux, iteration_callback);
}

struct OneIterationFluxResult
{
    cfd::IncompressibleSimpleResult solve_result;
    std::vector<double> mass_flux;
};

[[nodiscard]]
OneIterationFluxResult run_one_iteration_with_flux_relaxation(const cfd::Mesh &mesh, const double relaxation_factor,
                                                              const std::vector<double> &initial_mass_flux)
{
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    std::copy(initial_mass_flux.begin(), initial_mass_flux.end(), mass_flux.values().begin());
    cfd::IncompressibleSimpleOptions options{test_options()};
    options.maximum_iterations = 1;
    options.velocity_relative_tolerance = 1.0e-14;
    options.rhie_chow_flux_relaxation_factor = relaxation_factor;
    const cfd::IncompressibleSimpleResult result{
        solve_pressure_driven_channel(mesh, options, velocity, pressure, mass_flux)};
    return {result, {mass_flux.values().begin(), mass_flux.values().end()}};
}

void test_reports_each_completed_iteration()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_channel_raw_mesh(2, 2))};
    const cfd::Mesh &mesh{build_result.mesh};
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    cfd::IncompressibleSimpleOptions options{test_options()};
    options.maximum_iterations = 10;
    std::vector<cfd::SimpleIterationInfo> iteration_infos;
    iteration_infos.reserve(options.maximum_iterations);

    const cfd::IncompressibleSimpleResult result{solve_pressure_driven_channel(
        mesh, options, velocity, pressure, mass_flux, true,
        [&iteration_infos](const cfd::SimpleIterationInfo &info) { iteration_infos.push_back(info); })};

    require(result.converged, "SIMPLE callback fixture did not reach its converged final iteration.");
    require(iteration_infos.size() == result.iteration_count,
            "SIMPLE callback count differs from the completed outer-iteration count.");
    require(iteration_infos.size() > 1, "SIMPLE callback fixture did not complete multiple outer iterations.");
    constexpr double residual_upper_bound{1.0 + 64.0 * std::numeric_limits<double>::epsilon()};
    for (std::size_t index = 0; index < iteration_infos.size(); ++index)
    {
        const cfd::SimpleIterationInfo &info{iteration_infos[index]};
        require(info.iteration == index + 1, "SIMPLE callback iteration numbers are not consecutive and 1-based.");
        require(info.u_solve.converged && info.v_solve.converged && info.pressure_correction_solve.converged,
                "SIMPLE callback reported a non-converged inner linear solve.");
        require(std::isfinite(info.u_solve.estimated_relative_error) &&
                    std::isfinite(info.v_solve.estimated_relative_error) &&
                    std::isfinite(info.pressure_correction_solve.estimated_relative_error),
                "SIMPLE callback reported a non-finite inner linear-solve error.");
        require(std::isfinite(info.x_velocity_equation_residual) && info.x_velocity_equation_residual >= 0.0 &&
                    info.x_velocity_equation_residual <= residual_upper_bound,
                "SIMPLE callback reported an invalid x-velocity equation residual.");
        require(std::isfinite(info.y_velocity_equation_residual) && info.y_velocity_equation_residual >= 0.0 &&
                    info.y_velocity_equation_residual <= residual_upper_bound,
                "SIMPLE callback reported an invalid y-velocity equation residual.");
        require(std::isfinite(info.rhie_chow_flux_relative_residual) && info.rhie_chow_flux_relative_residual >= 0.0,
                "SIMPLE callback reported an invalid Rhie-Chow flux residual.");
    }

    const cfd::SimpleIterationInfo &final_info{iteration_infos.back()};
    require(final_info.iteration == result.iteration_count,
            "Final SIMPLE callback iteration differs from the returned iteration count.");
    require(final_info.velocity_relative_change == result.velocity_relative_change,
            "Final SIMPLE callback velocity diagnostic differs from the returned result.");
    require(final_info.provisional_continuity_relative_residual == result.provisional_continuity_relative_residual,
            "Final SIMPLE callback provisional-continuity diagnostic differs from the returned result.");
    require(final_info.corrected_continuity_relative_residual == result.continuity_relative_residual,
            "Final SIMPLE callback corrected-continuity diagnostic differs from the returned result.");
    require(final_info.maximum_pressure_correction == result.maximum_pressure_correction,
            "Final SIMPLE callback pressure-correction diagnostic differs from the returned result.");
    require(final_info.rhie_chow_flux_relative_residual == result.rhie_chow_flux_relative_residual,
            "Final SIMPLE callback Rhie-Chow flux diagnostic differs from the returned result.");
    const cfd::SimpleTimingBreakdown &timings{result.timings};
    const std::array phase_durations{
        timings.velocity_gradient_reconstruction_seconds,
        timings.pressure_gradient_reconstruction_seconds,
        timings.momentum_assembly_seconds,
        timings.momentum_residual_diagnostics_seconds,
        timings.momentum_matrix_preparation_seconds,
        timings.u_momentum_linear_solve_seconds,
        timings.v_momentum_linear_solve_seconds,
        timings.momentum_pressure_response_seconds,
        timings.rhie_chow_interpolation_seconds,
        timings.pressure_correction_assembly_seconds,
        timings.provisional_continuity_diagnostics_seconds,
        timings.pressure_correction_matrix_preparation_seconds,
        timings.pressure_correction_linear_solve_seconds,
        timings.pressure_correction_gradient_reconstruction_seconds,
        timings.field_correction_seconds,
        timings.convergence_diagnostics_seconds,
    };
    for (const double duration : phase_durations)
    {
        require(std::isfinite(duration) && duration >= 0.0, "SIMPLE returned an invalid cumulative phase duration.");
    }
    require(std::isfinite(timings.total_seconds) && timings.total_seconds > 0.0,
            "SIMPLE returned an invalid total solve duration.");
}

void test_rhie_chow_flux_relaxation_path()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_channel_raw_mesh(2, 2))};
    const cfd::Mesh &mesh{build_result.mesh};
    constexpr cfd::Vector2 uniform_flux_velocity{0.01, 0.001};
    std::vector<double> initial_mass_flux(mesh.face_count());
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        const cfd::Vector2 &area_vector{mesh.face_area_vectors()[face_id]};
        initial_mass_flux[face_id] = uniform_flux_velocity.x * area_vector.x + uniform_flux_velocity.y * area_vector.y;
    }

    // A uniform vector dotted with every owner-oriented area vector is
    // divergence-free. The linear pressure-correction projection therefore
    // preserves the same blend that was applied to the provisional flux.
    const OneIterationFluxResult unrelaxed{run_one_iteration_with_flux_relaxation(mesh, 1.0, initial_mass_flux)};
    constexpr double relaxation_factor{0.25};
    const OneIterationFluxResult relaxed{
        run_one_iteration_with_flux_relaxation(mesh, relaxation_factor, initial_mass_flux)};
    require(unrelaxed.solve_result.maximum_mass_imbalance < 1.0e-12 &&
                relaxed.solve_result.maximum_mass_imbalance < 1.0e-12,
            "Flux relaxation prevented pressure correction from restoring continuity.");

    const auto face_adjacencies{mesh.face_adjacencies()};
    const auto face_boundary_ids{mesh.face_boundary_ids()};
    const cfd::PressureCorrectionBoundaryConditions pressure_correction_conditions{
        channel_pressure_correction_conditions(mesh)};
    bool checked_internal_face{};
    bool checked_fixed_pressure_face{};
    bool checked_nonzero_fixed_mass_flux_face{};
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (face_adjacencies[face_id].is_boundary() && pressure_correction_conditions[face_boundary_ids[face_id]] ==
                                                           cfd::PressureCorrectionBoundaryConditionType::FixedMassFlux)
        {
            require(relaxed.mass_flux[face_id] == initial_mass_flux[face_id],
                    "Flux relaxation changed a FixedMassFlux boundary value.");
            checked_nonzero_fixed_mass_flux_face =
                checked_nonzero_fixed_mass_flux_face || initial_mass_flux[face_id] != 0.0;
            continue;
        }

        checked_internal_face = checked_internal_face || !face_adjacencies[face_id].is_boundary();
        checked_fixed_pressure_face = checked_fixed_pressure_face || face_adjacencies[face_id].is_boundary();

        const double expected_flux{relaxation_factor * unrelaxed.mass_flux[face_id] +
                                   (1.0 - relaxation_factor) * initial_mass_flux[face_id]};
        require(std::abs(relaxed.mass_flux[face_id] - expected_flux) < 1.0e-11,
                "Computed-face flux does not follow the requested Rhie-Chow relaxation blend.");
    }
    require(checked_internal_face && checked_fixed_pressure_face && checked_nonzero_fixed_mass_flux_face,
            "Flux-relaxation fixture did not exercise every required face category.");
}

void test_pressure_correction_gradient_boundary_mapping()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_channel_raw_mesh(8, 4))};
    const cfd::Mesh &mesh{build_result.mesh};
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    cfd::FaceFluxField mass_flux{mesh.face_count()};

    const cfd::IncompressibleSimpleResult result{
        solve_pressure_driven_channel(mesh, test_options(), velocity, pressure, mass_flux, false)};

    require(result.converged, "Channel with initially uniform pressure did not converge.");
    double maximum_pressure_error{};
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        const double exact_pressure{0.02 * (1.0 - mesh.cell_centers()[cell_id].x / 2.0)};
        maximum_pressure_error = std::max(maximum_pressure_error, std::abs(pressure[cell_id] - exact_pressure));
    }
    require(maximum_pressure_error < 1.0e-8,
            "FixedPressure/FixedMassFlux p-prime reconstruction mapping produced an incorrect pressure field.");
}

void test_maximum_iteration_nonconvergence_and_channel_convergence()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_channel_raw_mesh(8, 4))};
    const cfd::Mesh &mesh{build_result.mesh};
    {
        cfd::CellVelocityField velocity{mesh.cell_count()};
        cfd::CellScalarField pressure{mesh.cell_count()};
        cfd::FaceFluxField mass_flux{mesh.face_count()};
        cfd::IncompressibleSimpleOptions options{test_options()};
        options.maximum_iterations = 1;
        options.velocity_relative_tolerance = 1.0e-14;
        const cfd::IncompressibleSimpleResult result{
            solve_pressure_driven_channel(mesh, options, velocity, pressure, mass_flux)};
        require(!result.converged && result.iteration_count == 1,
                "SIMPLE did not report maximum-iteration non-convergence.");
    }
    {
        cfd::CellVelocityField velocity{mesh.cell_count()};
        cfd::CellScalarField pressure{mesh.cell_count()};
        cfd::FaceFluxField mass_flux{mesh.face_count()};
        const cfd::IncompressibleSimpleResult result{
            solve_pressure_driven_channel(mesh, test_options(), velocity, pressure, mass_flux)};
        require(result.converged, "Small deterministic pressure-driven channel did not converge.");
        require(result.velocity_relative_change <= test_options().velocity_relative_tolerance,
                "Converged channel exceeds the velocity-change tolerance.");
        require(result.continuity_relative_residual <= test_options().continuity_relative_tolerance,
                "Converged channel exceeds the continuity tolerance.");
        require(result.provisional_continuity_relative_residual <= test_options().continuity_relative_tolerance,
                "Converged channel exceeds the provisional continuity tolerance.");
        require(result.maximum_mass_imbalance < 1.0e-12,
                "Converged channel retains excessive absolute mass imbalance.");
    }
}

void test_does_not_converge_while_relaxed_face_flux_is_still_changing()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_channel_raw_mesh(1, 1, 1.0, 1.0))};
    const cfd::Mesh &mesh{build_result.mesh};

    const cfd::ScalarBoundaryConditions velocity_conditions{channel_velocity_conditions(mesh)};
    const cfd::ScalarBoundaryConditions pressure_conditions{channel_pressure_conditions(mesh, 1.0, 0.0)};
    const cfd::PressureCorrectionBoundaryConditions pressure_correction_conditions{
        channel_pressure_correction_conditions(mesh)};

    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count(), 0.5};
    cfd::FaceFluxField mass_flux{mesh.face_count()};

    cfd::IncompressibleSimpleOptions options{test_options()};
    options.maximum_iterations = 2;
    options.pressure_relaxation_factor = 0.3;
    options.rhie_chow_flux_relaxation_factor = 0.3;
    options.velocity_relative_tolerance = 1.0e-12;
    options.rhie_chow_flux_relative_tolerance = 1.0e-12;
    options.continuity_relative_tolerance = 1.0e-12;

    cfd::IncompressibleSimpleSolver solver{
        mesh, 1.0, 1.0, cfd::ScalarConvectionScheme::Linear, options,
    };

    const cfd::IncompressibleSimpleResult result{solver.solve(velocity_conditions, velocity_conditions,
                                                              pressure_conditions, pressure_correction_conditions,
                                                              velocity, pressure, mass_flux)};

    require(result.velocity_relative_change <= options.velocity_relative_tolerance,
            "Relaxed-flux false-convergence fixture does not satisfy the velocity-change criterion.");

    require(result.provisional_continuity_relative_residual <= options.continuity_relative_tolerance,
            "Relaxed-flux false-convergence fixture does not satisfy the provisional-continuity criterion.");
    require(result.rhie_chow_flux_relative_residual > options.rhie_chow_flux_relative_tolerance,
            "Relaxed-flux false-convergence fixture does not retain a significant "
            "Rhie-Chow flux fixed-point residual.");
    require(!result.converged, "SIMPLE reported convergence while the relaxed face flux was still changing.");
}

void test_does_not_converge_while_pressure_correction_remains_large()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_channel_raw_mesh(1, 1, 1.0, 1.0))};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarBoundaryConditions velocity_conditions{channel_velocity_conditions(mesh)};
    const cfd::ScalarBoundaryConditions pressure_conditions{channel_pressure_conditions(mesh, 1.0, 0.0)};
    const cfd::PressureCorrectionBoundaryConditions pressure_correction_conditions{
        channel_pressure_correction_conditions(mesh)};
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    cfd::IncompressibleSimpleOptions options{test_options()};
    options.maximum_iterations = 2;
    options.pressure_relaxation_factor = 0.3;
    options.rhie_chow_flux_relaxation_factor = 1.0;
    options.velocity_relative_tolerance = 1.0e-12;
    options.continuity_relative_tolerance = 1.0e-12;
    cfd::IncompressibleSimpleSolver solver{mesh, 1.0, 1.0, cfd::ScalarConvectionScheme::Linear, options};

    const cfd::IncompressibleSimpleResult result{solver.solve(velocity_conditions, velocity_conditions,
                                                              pressure_conditions, pressure_correction_conditions,
                                                              velocity, pressure, mass_flux)};

    require(result.velocity_relative_change <= options.velocity_relative_tolerance,
            "False-convergence fixture does not satisfy the old velocity-change criterion.");
    require(result.continuity_relative_residual <= options.continuity_relative_tolerance,
            "False-convergence fixture does not satisfy the old corrected-continuity criterion.");
    require(result.provisional_continuity_relative_residual > options.continuity_relative_tolerance,
            "False-convergence fixture does not retain a significant provisional continuity residual.");
    require(result.maximum_pressure_correction > 0.1,
            "False-convergence fixture no longer retains a substantial pressure correction.");
    require(!result.converged, "SIMPLE reported convergence while pressure correction remains substantial.");
}

} // namespace

int main()
{
    int failure_count{};

    failure_count +=
        cfd::test::run_test("SIMPLE constructor and options validation", test_constructor_and_options_validation);
    failure_count +=
        cfd::test::run_test("SIMPLE connected cell-domain precondition", test_rejects_disconnected_cell_domain);
    failure_count +=
        cfd::test::run_test("SIMPLE pressure boundary semantics", test_pressure_boundary_semantic_consistency);
    failure_count += cfd::test::run_test("SIMPLE all-FixedMassFlux compatibility and gauge",
                                         test_all_fixed_mass_flux_compatibility_and_gauge_path);
    failure_count +=
        cfd::test::run_test("SIMPLE FixedPressure anchored path", test_fixed_pressure_anchored_zero_flow_path);
    failure_count += cfd::test::run_test("SIMPLE solve input validation", test_input_validation_precedes_iterations);
    failure_count += cfd::test::run_test("SIMPLE p-prime gradient boundary mapping",
                                         test_pressure_correction_gradient_boundary_mapping);
    failure_count += cfd::test::run_test("SIMPLE completed-iteration callback", test_reports_each_completed_iteration);
    failure_count += cfd::test::run_test("SIMPLE Rhie-Chow flux relaxation", test_rhie_chow_flux_relaxation_path);
    failure_count += cfd::test::run_test("SIMPLE channel convergence",
                                         test_maximum_iteration_nonconvergence_and_channel_convergence);
    failure_count += cfd::test::run_test("SIMPLE rejects convergence while relaxed face flux is changing",
                                         test_does_not_converge_while_relaxed_face_flux_is_still_changing);
    failure_count += cfd::test::run_test("SIMPLE rejects convergence with a large pressure correction",
                                         test_does_not_converge_while_pressure_correction_remains_large);
    return cfd::test::finish_tests(failure_count, "incompressible SIMPLE solver");
}
