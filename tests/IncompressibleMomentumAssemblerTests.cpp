#include "cfd/numerics/IncompressibleMomentumAssembler.hpp"

#include "cfd/field/CellVectorField.hpp"
#include "cfd/field/CellVelocityField.hpp"
#include "cfd/field/FaceFluxField.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/linear_algebra/ScalarLinearSystem.hpp"
#include "cfd/math/Vector2.hpp"
#include "cfd/mesh/Boundary.hpp"
#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/mesh/Types.hpp"
#include "cfd/meshing/RawMeshData.hpp"
#include "cfd/numerics/ScalarConvectionOperator.hpp"
#include "cfd/numerics/ScalarDiffusionOperator.hpp"

#include "support/MeshFixtures.hpp"
#include "support/TestUtils.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{

using cfd::test::make_single_quadrilateral_raw_mesh;
using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws;
using cfd::test::test_tolerance;

[[nodiscard]]
cfd::RawMeshData make_two_cell_sheared_raw_mesh()
{
    cfd::RawMeshData raw_mesh;
    raw_mesh.nodes = {
        {0.0, 0.0}, {1.0, 0.0}, {2.0, 0.0}, {0.4, 1.0}, {1.4, 1.0}, {2.4, 1.0},
    };
    raw_mesh.cell_types = {
        cfd::CellType::Quadrilateral,
        cfd::CellType::Quadrilateral,
    };
    raw_mesh.cell_nodes = {
        0, 1, 4, 3, 1, 2, 5, 4,
    };
    raw_mesh.cell_node_offsets = {0, 4, 8};
    raw_mesh.boundary_groups = {
        {0, "bottom_0"}, {1, "bottom_1"}, {2, "right"}, {3, "top_1"}, {4, "top_0"}, {5, "left"},
    };
    raw_mesh.boundary_edges = {
        {{0, 1}, 0}, {{1, 2}, 1}, {{2, 5}, 2}, {{5, 4}, 3}, {{4, 3}, 4}, {{3, 0}, 5},
    };
    return raw_mesh;
}

[[nodiscard]]
cfd::ScalarBoundaryConditions make_uniform_boundary_conditions(const cfd::Index boundary_count,
                                                               const cfd::ScalarBoundaryConditionType type,
                                                               const double value)
{
    std::vector<cfd::ScalarBoundaryCondition> conditions;
    conditions.reserve(boundary_count);
    for (cfd::Index boundary_id = 0; boundary_id < boundary_count; ++boundary_id)
    {
        conditions.emplace_back(type, value);
    }
    return {boundary_count, std::move(conditions)};
}

void require_system_near(const cfd::ScalarLinearSystem &actual, const cfd::ScalarLinearSystem &expected,
                         const std::string &context)
{
    require(actual.cell_count() == expected.cell_count(), context + " cell counts differ.");
    require(actual.face_count() == expected.face_count(), context + " face counts differ.");

    for (cfd::Index cell_id = 0; cell_id < actual.cell_count(); ++cell_id)
    {
        require_near(actual.diagonal()[cell_id], expected.diagonal()[cell_id], test_tolerance,
                     context + " diagonal differs.");
        require_near(actual.rhs()[cell_id], expected.rhs()[cell_id], test_tolerance, context + " RHS differs.");
    }

    for (cfd::Index face_id = 0; face_id < actual.face_count(); ++face_id)
    {
        require_near(actual.owner_neighbor_coefficients()[face_id], expected.owner_neighbor_coefficients()[face_id],
                     test_tolerance, context + " owner-neighbor coefficient differs.");
        require_near(actual.neighbor_owner_coefficients()[face_id], expected.neighbor_owner_coefficients()[face_id],
                     test_tolerance, context + " neighbor-owner coefficient differs.");
    }
}

void seed_system(cfd::ScalarLinearSystem &system, const double value)
{
    std::fill(system.diagonal().begin(), system.diagonal().end(), value);
    std::fill(system.owner_neighbor_coefficients().begin(), system.owner_neighbor_coefficients().end(), value + 1.0);
    std::fill(system.neighbor_owner_coefficients().begin(), system.neighbor_owner_coefficients().end(), value + 2.0);
    std::fill(system.rhs().begin(), system.rhs().end(), value + 3.0);
}

void test_pressure_source_sign_and_component_separation()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_quadrilateral_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressibleMomentumAssembler assembler{mesh, 1.0};
    const cfd::CellVelocityField previous_velocity{mesh.cell_count()};
    const cfd::CellVectorField zero_gradient{mesh.cell_count()};
    const cfd::CellVectorField pressure_gradient{mesh.cell_count(), {2.0, -3.0}};
    const cfd::ScalarBoundaryConditions boundary_conditions{make_uniform_boundary_conditions(
        mesh.boundary_groups().size(), cfd::ScalarBoundaryConditionType::Neumann, 0.0)};
    const cfd::FaceFluxField mass_flux{mesh.face_count()};
    cfd::ScalarLinearSystem u_system{mesh};
    cfd::ScalarLinearSystem v_system{mesh};

    assembler.assemble(previous_velocity, zero_gradient, zero_gradient, pressure_gradient, boundary_conditions,
                       boundary_conditions, mass_flux, 1.0, u_system, v_system);

    const double cell_area{mesh.cell_areas()[0]};
    require_near(u_system.rhs()[0], -2.0 * cell_area, test_tolerance,
                 "The x-pressure gradient has the wrong u-momentum source sign.");
    require_near(v_system.rhs()[0], 3.0 * cell_area, test_tolerance,
                 "The y-pressure gradient has the wrong v-momentum source sign.");
}

void test_reuses_scalar_operators_and_alpha_one_preserves_unrelaxed_assembly()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_sheared_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    constexpr double dynamic_viscosity{0.7};
    constexpr cfd::ScalarConvectionScheme scheme{cfd::ScalarConvectionScheme::Linear};
    const cfd::IncompressibleMomentumAssembler assembler{mesh, dynamic_viscosity, scheme};
    cfd::CellVelocityField previous_velocity{mesh.cell_count()};
    previous_velocity.u()[0] = 1.2;
    previous_velocity.u()[1] = -0.4;
    previous_velocity.v()[0] = -2.0;
    previous_velocity.v()[1] = 0.8;
    cfd::CellVectorField u_gradient{mesh.cell_count()};
    u_gradient[0] = {0.8, -2.1};
    u_gradient[1] = {-1.3, 4.2};
    cfd::CellVectorField v_gradient{mesh.cell_count()};
    v_gradient[0] = {-0.2, 1.7};
    v_gradient[1] = {2.4, -0.6};
    cfd::CellVectorField pressure_gradient{mesh.cell_count()};
    pressure_gradient[0] = {2.0, -3.0};
    pressure_gradient[1] = {-1.5, 0.25};
    std::vector<cfd::ScalarBoundaryCondition> u_conditions{
        {cfd::ScalarBoundaryConditionType::Dirichlet, 1.2},  {cfd::ScalarBoundaryConditionType::Neumann, -0.4},
        {cfd::ScalarBoundaryConditionType::Dirichlet, 0.7},  {cfd::ScalarBoundaryConditionType::Neumann, 0.3},
        {cfd::ScalarBoundaryConditionType::Dirichlet, -1.0}, {cfd::ScalarBoundaryConditionType::Neumann, 0.2},
    };
    std::vector<cfd::ScalarBoundaryCondition> v_conditions{
        {cfd::ScalarBoundaryConditionType::Neumann, 0.1},  {cfd::ScalarBoundaryConditionType::Dirichlet, -2.0},
        {cfd::ScalarBoundaryConditionType::Neumann, 0.5},  {cfd::ScalarBoundaryConditionType::Dirichlet, 0.9},
        {cfd::ScalarBoundaryConditionType::Neumann, -0.2}, {cfd::ScalarBoundaryConditionType::Dirichlet, 1.4},
    };
    const cfd::ScalarBoundaryConditions u_boundary_conditions{mesh.boundary_groups().size(), std::move(u_conditions)};
    const cfd::ScalarBoundaryConditions v_boundary_conditions{mesh.boundary_groups().size(), std::move(v_conditions)};
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    for (cfd::Index face_id = 0; face_id < mass_flux.size(); ++face_id)
    {
        const double magnitude{0.2 * static_cast<double>(face_id + 1)};
        mass_flux[face_id] = face_id % 2 == 0 ? magnitude : -magnitude;
    }
    cfd::ScalarLinearSystem actual_u_system{mesh};
    cfd::ScalarLinearSystem actual_v_system{mesh};

    assembler.assemble(previous_velocity, u_gradient, v_gradient, pressure_gradient, u_boundary_conditions,
                       v_boundary_conditions, mass_flux, 1.0, actual_u_system, actual_v_system);

    const cfd::ScalarDiffusionOperator diffusion{mesh, dynamic_viscosity};
    const cfd::ScalarConvectionOperator convection{mesh, scheme};
    cfd::ScalarLinearSystem expected_u_system{mesh};
    cfd::ScalarLinearSystem expected_v_system{mesh};
    diffusion.add_matrix_contributions(u_boundary_conditions, expected_u_system);
    convection.add_matrix_contributions(u_boundary_conditions, mass_flux, expected_u_system);
    diffusion.add_boundary_rhs(u_boundary_conditions, expected_u_system.rhs());
    convection.add_boundary_rhs(u_boundary_conditions, mass_flux, expected_u_system.rhs());
    diffusion.add_non_orthogonal_rhs(u_boundary_conditions, u_gradient, expected_u_system.rhs());
    diffusion.add_matrix_contributions(v_boundary_conditions, expected_v_system);
    convection.add_matrix_contributions(v_boundary_conditions, mass_flux, expected_v_system);
    diffusion.add_boundary_rhs(v_boundary_conditions, expected_v_system.rhs());
    convection.add_boundary_rhs(v_boundary_conditions, mass_flux, expected_v_system.rhs());
    diffusion.add_non_orthogonal_rhs(v_boundary_conditions, v_gradient, expected_v_system.rhs());
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        expected_u_system.rhs()[cell_id] -= mesh.cell_areas()[cell_id] * pressure_gradient[cell_id].x;
        expected_v_system.rhs()[cell_id] -= mesh.cell_areas()[cell_id] * pressure_gradient[cell_id].y;
    }

    require_system_near(actual_u_system, expected_u_system, "The alpha-one u-momentum assembly");
    require_system_near(actual_v_system, expected_v_system, "The alpha-one v-momentum assembly");
}

void test_equation_under_relaxation()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_quadrilateral_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressibleMomentumAssembler assembler{mesh, 2.0};
    const cfd::CellVelocityField previous_velocity{mesh.cell_count(), {3.0, -4.0}};
    const cfd::CellVectorField zero_gradient{mesh.cell_count()};
    const cfd::ScalarBoundaryConditions u_boundary_conditions{make_uniform_boundary_conditions(
        mesh.boundary_groups().size(), cfd::ScalarBoundaryConditionType::Dirichlet, 1.0)};
    const cfd::ScalarBoundaryConditions v_boundary_conditions{make_uniform_boundary_conditions(
        mesh.boundary_groups().size(), cfd::ScalarBoundaryConditionType::Dirichlet, 2.0)};
    const cfd::FaceFluxField mass_flux{mesh.face_count()};
    cfd::ScalarLinearSystem unrelaxed_u_system{mesh};
    cfd::ScalarLinearSystem unrelaxed_v_system{mesh};
    cfd::ScalarLinearSystem relaxed_u_system{mesh};
    cfd::ScalarLinearSystem relaxed_v_system{mesh};

    assembler.assemble(previous_velocity, zero_gradient, zero_gradient, zero_gradient, u_boundary_conditions,
                       v_boundary_conditions, mass_flux, 1.0, unrelaxed_u_system, unrelaxed_v_system);
    constexpr double relaxation_factor{0.5};
    assembler.assemble(previous_velocity, zero_gradient, zero_gradient, zero_gradient, u_boundary_conditions,
                       v_boundary_conditions, mass_flux, relaxation_factor, relaxed_u_system, relaxed_v_system);

    const double relaxation_rhs_factor{(1.0 - relaxation_factor) / relaxation_factor};
    require_near(unrelaxed_u_system.diagonal()[0], 16.0, test_tolerance,
                 "The relaxation fixture has an unexpected unrelaxed u diagonal.");
    require_near(unrelaxed_v_system.diagonal()[0], 16.0, test_tolerance,
                 "The relaxation fixture has an unexpected unrelaxed v diagonal.");
    require_near(unrelaxed_u_system.rhs()[0], 16.0, test_tolerance,
                 "The relaxation fixture has an unexpected unrelaxed u RHS.");
    require_near(unrelaxed_v_system.rhs()[0], 32.0, test_tolerance,
                 "The relaxation fixture has an unexpected unrelaxed v RHS.");
    require_near(relaxed_u_system.diagonal()[0], unrelaxed_u_system.diagonal()[0] / relaxation_factor, test_tolerance,
                 "The u-momentum diagonal was not equation-under-relaxed.");
    require_near(relaxed_v_system.diagonal()[0], unrelaxed_v_system.diagonal()[0] / relaxation_factor, test_tolerance,
                 "The v-momentum diagonal was not equation-under-relaxed.");
    require_near(relaxed_u_system.rhs()[0],
                 unrelaxed_u_system.rhs()[0] + relaxation_rhs_factor * unrelaxed_u_system.diagonal()[0] * 3.0,
                 test_tolerance, "The u-momentum relaxation source is incorrect.");
    require_near(relaxed_v_system.rhs()[0],
                 unrelaxed_v_system.rhs()[0] + relaxation_rhs_factor * unrelaxed_v_system.diagonal()[0] * -4.0,
                 test_tolerance, "The v-momentum relaxation source is incorrect.");
}

void test_non_orthogonal_correction_reuses_diffusion_operator()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_sheared_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressibleMomentumAssembler assembler{mesh, 1.0};
    const cfd::CellVelocityField previous_velocity{mesh.cell_count()};
    const cfd::CellVectorField u_gradient{mesh.cell_count(), {0.0, 1.0}};
    const cfd::CellVectorField v_gradient{mesh.cell_count(), {1.0, -0.5}};
    const cfd::CellVectorField zero_pressure_gradient{mesh.cell_count()};
    const cfd::ScalarBoundaryConditions boundary_conditions{make_uniform_boundary_conditions(
        mesh.boundary_groups().size(), cfd::ScalarBoundaryConditionType::Neumann, 0.0)};
    const cfd::FaceFluxField mass_flux{mesh.face_count()};
    cfd::ScalarLinearSystem u_system{mesh};
    cfd::ScalarLinearSystem v_system{mesh};

    assembler.assemble(previous_velocity, u_gradient, v_gradient, zero_pressure_gradient, boundary_conditions,
                       boundary_conditions, mass_flux, 1.0, u_system, v_system);

    const cfd::ScalarDiffusionOperator diffusion{mesh, 1.0};
    std::vector<double> expected_u_rhs(mesh.cell_count());
    std::vector<double> expected_v_rhs(mesh.cell_count());
    diffusion.add_non_orthogonal_rhs(boundary_conditions, u_gradient, expected_u_rhs);
    diffusion.add_non_orthogonal_rhs(boundary_conditions, v_gradient, expected_v_rhs);
    require(std::abs(expected_u_rhs[0]) > test_tolerance || std::abs(expected_u_rhs[1]) > test_tolerance,
            "The non-orthogonal fixture produced a trivial u correction.");
    require(std::abs(expected_v_rhs[0]) > test_tolerance || std::abs(expected_v_rhs[1]) > test_tolerance,
            "The non-orthogonal fixture produced a trivial v correction.");
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        require_near(u_system.rhs()[cell_id], expected_u_rhs[cell_id], test_tolerance,
                     "The u non-orthogonal correction differs from ScalarDiffusionOperator.");
        require_near(v_system.rhs()[cell_id], expected_v_rhs[cell_id], test_tolerance,
                     "The v non-orthogonal correction differs from ScalarDiffusionOperator.");
    }
}

void test_output_systems_are_cleared()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_quadrilateral_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressibleMomentumAssembler assembler{mesh, 1.0};
    const cfd::CellVelocityField previous_velocity{mesh.cell_count()};
    const cfd::CellVectorField zero_gradient{mesh.cell_count()};
    const cfd::ScalarBoundaryConditions u_boundary_conditions{make_uniform_boundary_conditions(
        mesh.boundary_groups().size(), cfd::ScalarBoundaryConditionType::Dirichlet, 1.5)};
    const cfd::ScalarBoundaryConditions v_boundary_conditions{make_uniform_boundary_conditions(
        mesh.boundary_groups().size(), cfd::ScalarBoundaryConditionType::Dirichlet, -0.75)};
    const cfd::FaceFluxField mass_flux{mesh.face_count()};
    cfd::ScalarLinearSystem expected_u_system{mesh};
    cfd::ScalarLinearSystem expected_v_system{mesh};
    cfd::ScalarLinearSystem seeded_u_system{mesh};
    cfd::ScalarLinearSystem seeded_v_system{mesh};
    seed_system(seeded_u_system, 17.0);
    seed_system(seeded_v_system, -23.0);

    assembler.assemble(previous_velocity, zero_gradient, zero_gradient, zero_gradient, u_boundary_conditions,
                       v_boundary_conditions, mass_flux, 1.0, expected_u_system, expected_v_system);
    assembler.assemble(previous_velocity, zero_gradient, zero_gradient, zero_gradient, u_boundary_conditions,
                       v_boundary_conditions, mass_flux, 1.0, seeded_u_system, seeded_v_system);

    require_system_near(seeded_u_system, expected_u_system, "The rebuilt u-momentum system");
    require_system_near(seeded_v_system, expected_v_system, "The rebuilt v-momentum system");
}

void test_different_component_boundary_data_remain_independent()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_quadrilateral_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressibleMomentumAssembler assembler{mesh, 2.0};
    const cfd::CellVelocityField previous_velocity{mesh.cell_count()};
    const cfd::CellVectorField zero_gradient{mesh.cell_count()};
    const cfd::ScalarBoundaryConditions u_boundary_conditions{make_uniform_boundary_conditions(
        mesh.boundary_groups().size(), cfd::ScalarBoundaryConditionType::Dirichlet, 1.25)};
    const cfd::ScalarBoundaryConditions v_boundary_conditions{make_uniform_boundary_conditions(
        mesh.boundary_groups().size(), cfd::ScalarBoundaryConditionType::Dirichlet, -0.5)};
    const cfd::FaceFluxField mass_flux{mesh.face_count()};
    cfd::ScalarLinearSystem u_system{mesh};
    cfd::ScalarLinearSystem v_system{mesh};

    assembler.assemble(previous_velocity, zero_gradient, zero_gradient, zero_gradient, u_boundary_conditions,
                       v_boundary_conditions, mass_flux, 1.0, u_system, v_system);

    require_near(u_system.rhs()[0], 20.0, test_tolerance, "The u boundary value produced an incorrect RHS.");
    require_near(v_system.rhs()[0], -8.0, test_tolerance, "The v boundary value produced an incorrect RHS.");
}

void test_constructor_rejects_invalid_dynamic_viscosity()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_quadrilateral_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const auto require_rejected = [&mesh](const double dynamic_viscosity) {
        require_throws<std::invalid_argument>(
            [&mesh, dynamic_viscosity]() {
                const cfd::IncompressibleMomentumAssembler assembler{mesh, dynamic_viscosity};
            },
            "Momentum assembly accepted an invalid dynamic viscosity.");
    };

    require_rejected(0.0);
    require_rejected(-1.0);
    require_rejected(std::numeric_limits<double>::quiet_NaN());
    require_rejected(std::numeric_limits<double>::infinity());
}

void test_rejects_invalid_assembly_inputs_before_mutating_outputs()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_quadrilateral_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressibleMomentumAssembler assembler{mesh, 1.0};
    const cfd::CellVelocityField previous_velocity{mesh.cell_count()};
    const cfd::CellVectorField gradient{mesh.cell_count()};
    const cfd::ScalarBoundaryConditions boundary_conditions{make_uniform_boundary_conditions(
        mesh.boundary_groups().size(), cfd::ScalarBoundaryConditionType::Neumann, 0.0)};
    const cfd::FaceFluxField mass_flux{mesh.face_count()};
    cfd::ScalarLinearSystem u_system{mesh};
    cfd::ScalarLinearSystem v_system{mesh};
    seed_system(u_system, 31.0);
    seed_system(v_system, -37.0);

    const auto require_invalid_relaxation = [&](const double relaxation_factor) {
        require_throws<std::invalid_argument>(
            [&]() {
                assembler.assemble(previous_velocity, gradient, gradient, gradient, boundary_conditions,
                                   boundary_conditions, mass_flux, relaxation_factor, u_system, v_system);
            },
            "Momentum assembly accepted an invalid relaxation factor.");
    };
    require_invalid_relaxation(0.0);
    require_invalid_relaxation(-0.1);
    require_invalid_relaxation(1.1);
    require_invalid_relaxation(std::numeric_limits<double>::quiet_NaN());
    require_invalid_relaxation(std::numeric_limits<double>::infinity());

    require_throws<std::invalid_argument>(
        [&]() {
            const cfd::CellVelocityField wrong_velocity{mesh.cell_count() + 1};
            assembler.assemble(wrong_velocity, gradient, gradient, gradient, boundary_conditions, boundary_conditions,
                               mass_flux, 1.0, u_system, v_system);
        },
        "Momentum assembly accepted a velocity field with incorrect cardinality.");
    require_throws<std::invalid_argument>(
        [&]() {
            const cfd::CellVectorField wrong_gradient{mesh.cell_count() + 1};
            assembler.assemble(previous_velocity, wrong_gradient, gradient, gradient, boundary_conditions,
                               boundary_conditions, mass_flux, 1.0, u_system, v_system);
        },
        "Momentum assembly accepted a u gradient with incorrect cardinality.");
    require_throws<std::invalid_argument>(
        [&]() {
            const cfd::CellVectorField wrong_gradient{mesh.cell_count() + 1};
            assembler.assemble(previous_velocity, gradient, wrong_gradient, gradient, boundary_conditions,
                               boundary_conditions, mass_flux, 1.0, u_system, v_system);
        },
        "Momentum assembly accepted a v gradient with incorrect cardinality.");
    require_throws<std::invalid_argument>(
        [&]() {
            const cfd::CellVectorField wrong_gradient{mesh.cell_count() + 1};
            assembler.assemble(previous_velocity, gradient, gradient, wrong_gradient, boundary_conditions,
                               boundary_conditions, mass_flux, 1.0, u_system, v_system);
        },
        "Momentum assembly accepted a pressure gradient with incorrect cardinality.");
    require_throws<std::invalid_argument>(
        [&]() {
            const cfd::FaceFluxField wrong_mass_flux{mesh.face_count() + 1};
            assembler.assemble(previous_velocity, gradient, gradient, gradient, boundary_conditions,
                               boundary_conditions, wrong_mass_flux, 1.0, u_system, v_system);
        },
        "Momentum assembly accepted a mass flux with incorrect cardinality.");

    const cfd::ScalarBoundaryConditions wrong_boundary_conditions{make_uniform_boundary_conditions(
        mesh.boundary_groups().size() + 1, cfd::ScalarBoundaryConditionType::Neumann, 0.0)};
    require_throws<std::invalid_argument>(
        [&]() {
            assembler.assemble(previous_velocity, gradient, gradient, gradient, wrong_boundary_conditions,
                               boundary_conditions, mass_flux, 1.0, u_system, v_system);
        },
        "Momentum assembly accepted an incorrect u boundary-condition count.");
    require_throws<std::invalid_argument>(
        [&]() {
            assembler.assemble(previous_velocity, gradient, gradient, gradient, boundary_conditions,
                               wrong_boundary_conditions, mass_flux, 1.0, u_system, v_system);
        },
        "Momentum assembly accepted an incorrect v boundary-condition count.");

    cfd::MeshBuildResult other_build_result{cfd::build_mesh(make_single_quadrilateral_raw_mesh())};
    cfd::ScalarLinearSystem other_system{other_build_result.mesh};
    require_throws<std::invalid_argument>(
        [&]() {
            assembler.assemble(previous_velocity, gradient, gradient, gradient, boundary_conditions,
                               boundary_conditions, mass_flux, 1.0, other_system, v_system);
        },
        "Momentum assembly accepted a u system referencing a different Mesh.");
    require_throws<std::invalid_argument>(
        [&]() {
            assembler.assemble(previous_velocity, gradient, gradient, gradient, boundary_conditions,
                               boundary_conditions, mass_flux, 1.0, u_system, other_system);
        },
        "Momentum assembly accepted a v system referencing a different Mesh.");
    require_throws<std::invalid_argument>(
        [&]() {
            assembler.assemble(previous_velocity, gradient, gradient, gradient, boundary_conditions,
                               boundary_conditions, mass_flux, 1.0, u_system, u_system);
        },
        "Momentum assembly accepted the same system for both outputs.");

    require_near(u_system.diagonal()[0], 31.0, 0.0, "Invalid input mutated the u-system diagonal.");
    require_near(u_system.rhs()[0], 34.0, 0.0, "Invalid input mutated the u-system RHS.");
    require_near(v_system.diagonal()[0], -37.0, 0.0, "Invalid input mutated the v-system diagonal.");
    require_near(v_system.rhs()[0], -34.0, 0.0, "Invalid input mutated the v-system RHS.");
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        require_near(u_system.owner_neighbor_coefficients()[face_id], 32.0, 0.0,
                     "Invalid input mutated a u-system owner-neighbor coefficient.");
        require_near(u_system.neighbor_owner_coefficients()[face_id], 33.0, 0.0,
                     "Invalid input mutated a u-system neighbor-owner coefficient.");
        require_near(v_system.owner_neighbor_coefficients()[face_id], -36.0, 0.0,
                     "Invalid input mutated a v-system owner-neighbor coefficient.");
        require_near(v_system.neighbor_owner_coefficients()[face_id], -35.0, 0.0,
                     "Invalid input mutated a v-system neighbor-owner coefficient.");
    }
}

void test_does_not_mutate_inputs()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_sheared_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressibleMomentumAssembler assembler{mesh, 0.9, cfd::ScalarConvectionScheme::Linear};
    cfd::CellVelocityField previous_velocity{mesh.cell_count()};
    previous_velocity.u()[0] = 1.0;
    previous_velocity.u()[1] = 2.0;
    previous_velocity.v()[0] = -3.0;
    previous_velocity.v()[1] = -4.0;
    cfd::CellVectorField u_gradient{mesh.cell_count()};
    u_gradient[0] = {0.1, 0.2};
    u_gradient[1] = {0.3, 0.4};
    cfd::CellVectorField v_gradient{mesh.cell_count()};
    v_gradient[0] = {-0.5, -0.6};
    v_gradient[1] = {-0.7, -0.8};
    cfd::CellVectorField pressure_gradient{mesh.cell_count()};
    pressure_gradient[0] = {2.0, -3.0};
    pressure_gradient[1] = {4.0, -5.0};
    const cfd::ScalarBoundaryConditions u_boundary_conditions{make_uniform_boundary_conditions(
        mesh.boundary_groups().size(), cfd::ScalarBoundaryConditionType::Dirichlet, 1.25)};
    const cfd::ScalarBoundaryConditions v_boundary_conditions{make_uniform_boundary_conditions(
        mesh.boundary_groups().size(), cfd::ScalarBoundaryConditionType::Neumann, -0.75)};
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        mass_flux[face_id] = static_cast<double>(face_id) - 2.5;
    }
    const cfd::CellVelocityField previous_velocity_before{previous_velocity};
    const cfd::CellVectorField u_gradient_before{u_gradient};
    const cfd::CellVectorField v_gradient_before{v_gradient};
    const cfd::CellVectorField pressure_gradient_before{pressure_gradient};
    const cfd::FaceFluxField mass_flux_before{mass_flux};
    cfd::ScalarLinearSystem u_system{mesh};
    cfd::ScalarLinearSystem v_system{mesh};

    assembler.assemble(previous_velocity, u_gradient, v_gradient, pressure_gradient, u_boundary_conditions,
                       v_boundary_conditions, mass_flux, 0.6, u_system, v_system);

    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        require_near(previous_velocity.u()[cell_id], previous_velocity_before.u()[cell_id], 0.0,
                     "Momentum assembly modified the previous u velocity.");
        require_near(previous_velocity.v()[cell_id], previous_velocity_before.v()[cell_id], 0.0,
                     "Momentum assembly modified the previous v velocity.");
        require_near(u_gradient[cell_id].x, u_gradient_before[cell_id].x, 0.0,
                     "Momentum assembly modified a u-gradient x component.");
        require_near(u_gradient[cell_id].y, u_gradient_before[cell_id].y, 0.0,
                     "Momentum assembly modified a u-gradient y component.");
        require_near(v_gradient[cell_id].x, v_gradient_before[cell_id].x, 0.0,
                     "Momentum assembly modified a v-gradient x component.");
        require_near(v_gradient[cell_id].y, v_gradient_before[cell_id].y, 0.0,
                     "Momentum assembly modified a v-gradient y component.");
        require_near(pressure_gradient[cell_id].x, pressure_gradient_before[cell_id].x, 0.0,
                     "Momentum assembly modified a pressure-gradient x component.");
        require_near(pressure_gradient[cell_id].y, pressure_gradient_before[cell_id].y, 0.0,
                     "Momentum assembly modified a pressure-gradient y component.");
    }
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        require_near(mass_flux[face_id], mass_flux_before[face_id], 0.0,
                     "Momentum assembly modified an input mass flux.");
    }
    for (cfd::BoundaryId boundary_id = 0; boundary_id < mesh.boundary_groups().size(); ++boundary_id)
    {
        require(u_boundary_conditions[boundary_id].type == cfd::ScalarBoundaryConditionType::Dirichlet,
                "Momentum assembly modified a u boundary-condition type.");
        require_near(u_boundary_conditions[boundary_id].value, 1.25, 0.0,
                     "Momentum assembly modified a u boundary-condition value.");
        require(v_boundary_conditions[boundary_id].type == cfd::ScalarBoundaryConditionType::Neumann,
                "Momentum assembly modified a v boundary-condition type.");
        require_near(v_boundary_conditions[boundary_id].value, -0.75, 0.0,
                     "Momentum assembly modified a v boundary-condition value.");
    }
}

} // namespace

int main()
{
    int failure_count{};

    failure_count += cfd::test::run_test("momentum pressure source sign and component separation",
                                         test_pressure_source_sign_and_component_separation);
    failure_count += cfd::test::run_test("momentum scalar-operator reuse and alpha-one assembly",
                                         test_reuses_scalar_operators_and_alpha_one_preserves_unrelaxed_assembly);
    failure_count += cfd::test::run_test("momentum equation under-relaxation", test_equation_under_relaxation);
    failure_count += cfd::test::run_test("momentum non-orthogonal correction",
                                         test_non_orthogonal_correction_reuses_diffusion_operator);
    failure_count += cfd::test::run_test("momentum output clearing", test_output_systems_are_cleared);
    failure_count += cfd::test::run_test("momentum component boundary independence",
                                         test_different_component_boundary_data_remain_independent);
    failure_count +=
        cfd::test::run_test("momentum viscosity validation", test_constructor_rejects_invalid_dynamic_viscosity);
    failure_count += cfd::test::run_test("momentum assembly input validation",
                                         test_rejects_invalid_assembly_inputs_before_mutating_outputs);
    failure_count += cfd::test::run_test("momentum input immutability", test_does_not_mutate_inputs);

    return cfd::test::finish_tests(failure_count, "incompressible momentum assembler");
}
