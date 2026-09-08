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

constexpr cfd::BoundaryId bottom_boundary_id{0};
constexpr cfd::BoundaryId right_boundary_id{1};
constexpr cfd::BoundaryId top_boundary_id{2};
constexpr cfd::BoundaryId left_boundary_id{3};

static_assert(!std::is_copy_constructible_v<cfd::IncompressibleSimpleSolver>);
static_assert(!std::is_copy_assignable_v<cfd::IncompressibleSimpleSolver>);
static_assert(!std::is_move_constructible_v<cfd::IncompressibleSimpleSolver>);
static_assert(!std::is_move_assignable_v<cfd::IncompressibleSimpleSolver>);

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
    }
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
    cfd::CellScalarField &pressure, cfd::FaceFluxField &mass_flux, const bool initialize_exact_pressure = true)
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
                        velocity, pressure, mass_flux);
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
        require(result.maximum_mass_imbalance < 1.0e-12,
                "Converged channel retains excessive absolute mass imbalance.");
    }
}

} // namespace

int main()
{
    int failure_count{};

    failure_count +=
        cfd::test::run_test("SIMPLE constructor and options validation", test_constructor_and_options_validation);
    failure_count +=
        cfd::test::run_test("SIMPLE pressure boundary semantics", test_pressure_boundary_semantic_consistency);
    failure_count += cfd::test::run_test("SIMPLE all-FixedMassFlux compatibility and gauge",
                                         test_all_fixed_mass_flux_compatibility_and_gauge_path);
    failure_count +=
        cfd::test::run_test("SIMPLE FixedPressure anchored path", test_fixed_pressure_anchored_zero_flow_path);
    failure_count += cfd::test::run_test("SIMPLE solve input validation", test_input_validation_precedes_iterations);
    failure_count += cfd::test::run_test("SIMPLE p-prime gradient boundary mapping",
                                         test_pressure_correction_gradient_boundary_mapping);
    failure_count += cfd::test::run_test("SIMPLE Rhie-Chow flux relaxation", test_rhie_chow_flux_relaxation_path);
    failure_count += cfd::test::run_test("SIMPLE channel convergence",
                                         test_maximum_iteration_nonconvergence_and_channel_convergence);

    return cfd::test::finish_tests(failure_count, "incompressible SIMPLE solver");
}
