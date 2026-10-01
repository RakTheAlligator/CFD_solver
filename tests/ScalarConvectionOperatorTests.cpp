#include "cfd/numerics/ScalarConvectionOperator.hpp"

#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVectorField.hpp"
#include "cfd/field/FaceFluxField.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/linear_algebra/ScalarLinearSystem.hpp"
#include "cfd/mesh/Boundary.hpp"
#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/meshing/RawMeshData.hpp"

#include "support/MeshFixtures.hpp"
#include "support/TestUtils.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace
{

using cfd::test::make_single_quadrilateral_raw_mesh;
using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws;
using cfd::test::test_tolerance;

static_assert(!std::is_copy_constructible_v<cfd::ScalarConvectionOperator>);
static_assert(!std::is_copy_assignable_v<cfd::ScalarConvectionOperator>);
static_assert(std::is_nothrow_move_constructible_v<cfd::ScalarConvectionOperator>);
static_assert(!std::is_move_assignable_v<cfd::ScalarConvectionOperator>);
static_assert(std::is_nothrow_constructible_v<cfd::ScalarConvectionOperator, const cfd::Mesh &>);

[[nodiscard]]
cfd::RawMeshData make_two_cell_rectangle_raw_mesh()
{
    constexpr cfd::BoundaryId wall_boundary_id{0};

    cfd::RawMeshData raw_mesh;
    raw_mesh.nodes = {
        {0.0, 0.0}, {1.0, 0.0}, {2.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}, {2.0, 1.0},
    };
    raw_mesh.cell_types = {
        cfd::CellType::Quadrilateral,
        cfd::CellType::Quadrilateral,
    };
    raw_mesh.cell_nodes = {
        0, 1, 4, 3, 1, 2, 5, 4,
    };
    raw_mesh.cell_node_offsets = {0, 4, 8};
    raw_mesh.boundary_groups = {{wall_boundary_id, "wall"}};
    raw_mesh.boundary_edges = {
        {{0, 1}, wall_boundary_id}, {{1, 2}, wall_boundary_id}, {{2, 5}, wall_boundary_id},
        {{5, 4}, wall_boundary_id}, {{4, 3}, wall_boundary_id}, {{3, 0}, wall_boundary_id},
    };
    return raw_mesh;
}

[[nodiscard]]
cfd::RawMeshData make_nonuniform_two_cell_rectangle_raw_mesh()
{
    cfd::RawMeshData raw_mesh{make_two_cell_rectangle_raw_mesh()};
    raw_mesh.nodes[2].x = 3.0;
    raw_mesh.nodes[5].x = 3.0;
    return raw_mesh;
}

[[nodiscard]]
cfd::RawMeshData make_single_cell_four_boundary_raw_mesh()
{
    cfd::RawMeshData raw_mesh{make_single_quadrilateral_raw_mesh()};
    raw_mesh.boundary_groups = {
        {0, "left"},
        {1, "right"},
        {2, "bottom"},
        {3, "top"},
    };
    raw_mesh.boundary_edges = {
        {{3, 0}, 0},
        {{1, 2}, 1},
        {{0, 1}, 2},
        {{2, 3}, 3},
    };
    return raw_mesh;
}

[[nodiscard]]
cfd::RawMeshData make_two_cell_six_boundary_raw_mesh()
{
    cfd::RawMeshData raw_mesh{make_two_cell_rectangle_raw_mesh()};
    raw_mesh.boundary_groups = {
        {0, "bottom_0"}, {1, "bottom_1"}, {2, "right"}, {3, "top_1"}, {4, "top_0"}, {5, "left"},
    };
    raw_mesh.boundary_edges = {
        {{0, 1}, 0}, {{1, 2}, 1}, {{2, 5}, 2}, {{5, 4}, 3}, {{4, 3}, 4}, {{3, 0}, 5},
    };
    return raw_mesh;
}

[[nodiscard]]
cfd::ScalarBoundaryConditions make_uniform_conditions(const cfd::Index boundary_count,
                                                      const cfd::ScalarBoundaryConditionType type, const double value)
{
    std::vector<cfd::ScalarBoundaryCondition> conditions;
    conditions.reserve(boundary_count);
    for (cfd::Index boundary_id = 0; boundary_id < boundary_count; ++boundary_id)
    {
        conditions.emplace_back(type, value);
    }
    return {boundary_count, std::move(conditions)};
}

[[nodiscard]]
cfd::Index internal_face_id(const cfd::Mesh &mesh)
{
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (!mesh.face_adjacencies()[face_id].is_boundary())
        {
            return face_id;
        }
    }
    throw std::runtime_error("Two-cell convection fixture has no internal face.");
}

[[nodiscard]]
cfd::Index first_boundary_face_id(const cfd::Mesh &mesh)
{
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (mesh.face_adjacencies()[face_id].is_boundary())
        {
            return face_id;
        }
    }
    throw std::runtime_error("Convection fixture has no boundary face.");
}

[[nodiscard]]
cfd::BoundaryId find_boundary_id(const cfd::Mesh &mesh, const std::string_view name)
{
    for (const cfd::BoundaryGroup &group : mesh.boundary_groups())
    {
        if (group.name == name)
        {
            return group.id;
        }
    }
    throw std::runtime_error("Convection fixture is missing boundary '" + std::string(name) + "'.");
}

void set_boundary_flux(const cfd::Mesh &mesh, cfd::FaceFluxField &face_flux, const cfd::BoundaryId boundary_id,
                       const double value)
{
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (mesh.face_boundary_ids()[face_id] == boundary_id)
        {
            face_flux[face_id] = value;
        }
    }
}

void test_positive_internal_flux()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator convection{mesh};
    const cfd::ScalarBoundaryConditions conditions{
        make_uniform_conditions(1, cfd::ScalarBoundaryConditionType::Neumann, 0.0)};
    cfd::CellScalarField field{mesh.cell_count()};
    field[0] = 2.0;
    field[1] = 5.0;
    cfd::FaceFluxField face_flux{mesh.face_count()};
    const cfd::Index face_id{internal_face_id(mesh)};
    face_flux[face_id] = 3.0;
    cfd::CellScalarField balance{mesh.cell_count()};
    cfd::ScalarLinearSystem system{mesh};

    convection.compute_flux_balance(field, conditions, face_flux, balance);
    convection.add_matrix_contributions(conditions, face_flux, system);

    require_near(balance[0], 6.0, 0.0, "Positive internal flux gave an incorrect owner balance.");
    require_near(balance[1], -6.0, 0.0, "Positive internal flux gave an incorrect neighbor balance.");
    require_near(system.diagonal()[0], 3.0, 0.0, "Positive internal flux gave an incorrect owner diagonal.");
    require_near(system.diagonal()[1], 0.0, 0.0, "Positive internal flux changed the neighbor diagonal.");
    require_near(system.owner_neighbor_coefficients()[face_id], 0.0, 0.0,
                 "Positive internal flux changed A(owner,neighbor).");
    require_near(system.neighbor_owner_coefficients()[face_id], -3.0, 0.0,
                 "Positive internal flux gave an incorrect A(neighbor,owner).");
}

void test_negative_internal_flux()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator convection{mesh};
    const cfd::ScalarBoundaryConditions conditions{
        make_uniform_conditions(1, cfd::ScalarBoundaryConditionType::Neumann, 0.0)};
    cfd::CellScalarField field{mesh.cell_count()};
    field[0] = 2.0;
    field[1] = 5.0;
    cfd::FaceFluxField face_flux{mesh.face_count()};
    const cfd::Index face_id{internal_face_id(mesh)};
    face_flux[face_id] = -4.0;
    cfd::CellScalarField balance{mesh.cell_count()};
    cfd::ScalarLinearSystem system{mesh};

    convection.compute_flux_balance(field, conditions, face_flux, balance);
    convection.add_matrix_contributions(conditions, face_flux, system);

    require_near(balance[0], -20.0, 0.0, "Negative internal flux gave an incorrect owner balance.");
    require_near(balance[1], 20.0, 0.0, "Negative internal flux gave an incorrect neighbor balance.");
    require_near(system.diagonal()[0], 0.0, 0.0, "Negative internal flux changed the owner diagonal.");
    require_near(system.diagonal()[1], 4.0, 0.0, "Negative internal flux gave an incorrect neighbor diagonal.");
    require_near(system.owner_neighbor_coefficients()[face_id], -4.0, 0.0,
                 "Negative internal flux gave an incorrect A(owner,neighbor).");
    require_near(system.neighbor_owner_coefficients()[face_id], 0.0, 0.0,
                 "Negative internal flux changed A(neighbor,owner).");
}

void test_internal_face_conservation()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator convection{mesh};
    const cfd::ScalarBoundaryConditions conditions{
        make_uniform_conditions(1, cfd::ScalarBoundaryConditionType::Neumann, 0.0)};
    cfd::CellScalarField field{mesh.cell_count()};
    field[0] = -1.25;
    field[1] = 3.5;
    cfd::FaceFluxField face_flux{mesh.face_count()};
    face_flux[internal_face_id(mesh)] = -2.75;
    cfd::CellScalarField balance{mesh.cell_count()};

    convection.compute_flux_balance(field, conditions, face_flux, balance);

    require_near(balance[0] + balance[1], 0.0, test_tolerance,
                 "Internal convective face contributions are not conservative.");
}

void test_dirichlet_inflow()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_quadrilateral_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator convection{mesh};
    const cfd::ScalarBoundaryConditions conditions{
        make_uniform_conditions(1, cfd::ScalarBoundaryConditionType::Dirichlet, 7.0)};
    const cfd::CellScalarField field{mesh.cell_count(), 3.0};
    cfd::FaceFluxField face_flux{mesh.face_count()};
    face_flux[first_boundary_face_id(mesh)] = -2.0;
    cfd::CellScalarField balance{mesh.cell_count()};
    cfd::ScalarLinearSystem system{mesh};

    convection.compute_flux_balance(field, conditions, face_flux, balance);
    convection.add_matrix_contributions(conditions, face_flux, system);
    convection.add_boundary_rhs(conditions, face_flux, system.rhs());

    require_near(balance[0], -14.0, 0.0, "Dirichlet inflow gave an incorrect direct balance.");
    require_near(system.diagonal()[0], 0.0, 0.0, "Dirichlet inflow changed the matrix diagonal.");
    require_near(system.rhs()[0], 14.0, 0.0, "Dirichlet inflow gave an incorrect boundary RHS.");
}

void test_dirichlet_outflow()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_quadrilateral_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator convection{mesh};
    const cfd::ScalarBoundaryConditions conditions{
        make_uniform_conditions(1, cfd::ScalarBoundaryConditionType::Dirichlet, 1000.0)};
    const cfd::CellScalarField field{mesh.cell_count(), 3.0};
    cfd::FaceFluxField face_flux{mesh.face_count()};
    face_flux[first_boundary_face_id(mesh)] = 2.0;
    cfd::CellScalarField balance{mesh.cell_count()};
    cfd::ScalarLinearSystem system{mesh};

    convection.compute_flux_balance(field, conditions, face_flux, balance);
    convection.add_matrix_contributions(conditions, face_flux, system);
    convection.add_boundary_rhs(conditions, face_flux, system.rhs());

    require_near(balance[0], 6.0, 0.0, "Dirichlet outflow did not use the owner value.");
    require_near(system.diagonal()[0], 2.0, 0.0, "Dirichlet outflow gave an incorrect diagonal.");
    require_near(system.rhs()[0], 0.0, 0.0, "Dirichlet outflow incorrectly used the boundary value.");
}

void test_zero_gradient_inflow_and_outflow()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_quadrilateral_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator convection{mesh};
    const cfd::ScalarBoundaryConditions conditions{
        make_uniform_conditions(1, cfd::ScalarBoundaryConditionType::Neumann, 0.0)};
    const cfd::CellScalarField field{mesh.cell_count(), 3.0};
    const cfd::Index face_id{first_boundary_face_id(mesh)};
    cfd::FaceFluxField face_flux{mesh.face_count()};
    cfd::CellScalarField balance{mesh.cell_count()};
    cfd::ScalarLinearSystem system{mesh};

    face_flux[face_id] = -2.0;
    convection.compute_flux_balance(field, conditions, face_flux, balance);
    convection.add_matrix_contributions(conditions, face_flux, system);
    convection.add_boundary_rhs(conditions, face_flux, system.rhs());
    require_near(balance[0], -6.0, 0.0, "zeroGradient inflow gave an incorrect balance.");
    require_near(system.diagonal()[0], -2.0, 0.0, "zeroGradient inflow gave an incorrect diagonal.");
    require_near(system.rhs()[0], 0.0, 0.0, "zeroGradient inflow changed the RHS.");

    face_flux[face_id] = 2.0;
    system.clear();
    convection.compute_flux_balance(field, conditions, face_flux, balance);
    convection.add_matrix_contributions(conditions, face_flux, system);
    convection.add_boundary_rhs(conditions, face_flux, system.rhs());
    require_near(balance[0], 6.0, 0.0, "zeroGradient outflow gave an incorrect balance.");
    require_near(system.diagonal()[0], 2.0, 0.0, "zeroGradient outflow gave an incorrect diagonal.");
    require_near(system.rhs()[0], 0.0, 0.0, "zeroGradient outflow changed the RHS.");
}

void test_nonzero_neumann_inflow()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_quadrilateral_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator convection{mesh};
    const cfd::ScalarBoundaryConditions conditions{
        make_uniform_conditions(1, cfd::ScalarBoundaryConditionType::Neumann, 4.0)};
    const cfd::CellScalarField field{mesh.cell_count(), 3.0};
    cfd::FaceFluxField face_flux{mesh.face_count()};
    face_flux[first_boundary_face_id(mesh)] = -2.0;
    cfd::CellScalarField balance{mesh.cell_count()};
    cfd::ScalarLinearSystem system{mesh};

    convection.compute_flux_balance(field, conditions, face_flux, balance);
    convection.add_matrix_contributions(conditions, face_flux, system);
    convection.add_boundary_rhs(conditions, face_flux, system.rhs());

    require_near(balance[0], -10.0, test_tolerance, "Non-zero Neumann inflow reconstructed an incorrect face value.");
    require_near(system.diagonal()[0], -2.0, 0.0, "Neumann inflow gave an incorrect diagonal.");
    require_near(system.rhs()[0], 4.0, test_tolerance, "Neumann inflow gave an incorrect RHS sign or distance.");
}

void test_nonzero_neumann_outflow()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_quadrilateral_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator convection{mesh};
    const cfd::ScalarBoundaryConditions conditions{
        make_uniform_conditions(1, cfd::ScalarBoundaryConditionType::Neumann, 1000.0)};
    const cfd::CellScalarField field{mesh.cell_count(), 3.0};
    cfd::FaceFluxField face_flux{mesh.face_count()};
    face_flux[first_boundary_face_id(mesh)] = 2.0;
    cfd::CellScalarField balance{mesh.cell_count()};
    cfd::ScalarLinearSystem system{mesh};

    convection.compute_flux_balance(field, conditions, face_flux, balance);
    convection.add_matrix_contributions(conditions, face_flux, system);
    convection.add_boundary_rhs(conditions, face_flux, system.rhs());

    require_near(balance[0], 6.0, 0.0, "Neumann gradient altered the upwind outflow value.");
    require_near(system.diagonal()[0], 2.0, 0.0, "Neumann outflow gave an incorrect diagonal.");
    require_near(system.rhs()[0], 0.0, 0.0, "Neumann outflow incorrectly changed the RHS.");
}

void test_zero_face_flux()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator convection{mesh};
    const cfd::ScalarBoundaryConditions conditions{
        make_uniform_conditions(1, cfd::ScalarBoundaryConditionType::Dirichlet, 9.0)};
    const cfd::CellScalarField field{mesh.cell_count(), 3.0};
    const cfd::FaceFluxField face_flux{mesh.face_count()};
    cfd::CellScalarField balance{mesh.cell_count(), 100.0};
    cfd::ScalarLinearSystem system{mesh};

    convection.compute_flux_balance(field, conditions, face_flux, balance);
    convection.add_matrix_contributions(conditions, face_flux, system);
    convection.add_boundary_rhs(conditions, face_flux, system.rhs());

    require(std::ranges::all_of(balance.values(), [](const double value) { return value == 0.0; }),
            "Zero face flux produced a non-zero direct balance.");
    require(std::ranges::all_of(system.diagonal(), [](const double value) { return value == 0.0; }) &&
                std::ranges::all_of(system.owner_neighbor_coefficients(),
                                    [](const double value) { return value == 0.0; }) &&
                std::ranges::all_of(system.neighbor_owner_coefficients(),
                                    [](const double value) { return value == 0.0; }) &&
                std::ranges::all_of(system.rhs(), [](const double value) { return value == 0.0; }),
            "Zero face flux produced non-zero assembled contributions.");
}

void test_constant_field_preservation()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_cell_four_boundary_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator convection{mesh};
    std::vector<cfd::ScalarBoundaryCondition> condition_values{
        {cfd::ScalarBoundaryConditionType::Dirichlet, 4.0},
        {cfd::ScalarBoundaryConditionType::Neumann, 0.0},
        {cfd::ScalarBoundaryConditionType::Dirichlet, 4.0},
        {cfd::ScalarBoundaryConditionType::Neumann, 0.0},
    };
    const cfd::ScalarBoundaryConditions conditions{4, std::move(condition_values)};
    const cfd::CellScalarField field{mesh.cell_count(), 4.0};
    cfd::FaceFluxField face_flux{mesh.face_count()};
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "left"), -1.0);
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "right"), 1.0);
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "bottom"), -2.0);
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "top"), 2.0);
    cfd::CellScalarField balance{mesh.cell_count()};

    convection.compute_flux_balance(field, conditions, face_flux, balance);

    require_near(balance[0], 0.0, test_tolerance,
                 "Divergence-free face flux did not preserve a compatible constant field.");
}

void test_assembly_matches_direct_balance()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_six_boundary_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator convection{mesh};
    std::vector<cfd::ScalarBoundaryCondition> condition_values{
        {cfd::ScalarBoundaryConditionType::Dirichlet, 1.2},  {cfd::ScalarBoundaryConditionType::Neumann, -0.4},
        {cfd::ScalarBoundaryConditionType::Dirichlet, -2.0}, {cfd::ScalarBoundaryConditionType::Neumann, 0.7},
        {cfd::ScalarBoundaryConditionType::Dirichlet, 3.0},  {cfd::ScalarBoundaryConditionType::Neumann, -1.1},
    };
    const cfd::ScalarBoundaryConditions conditions{6, std::move(condition_values)};
    cfd::CellScalarField field{mesh.cell_count()};
    field[0] = 2.3;
    field[1] = -0.8;
    cfd::FaceFluxField face_flux{mesh.face_count()};
    face_flux[internal_face_id(mesh)] = -1.7;
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "bottom_0"), -0.7);
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "bottom_1"), -1.1);
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "right"), 0.6);
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "top_1"), 0.9);
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "left"), -0.3);
    cfd::CellScalarField balance{mesh.cell_count()};
    cfd::ScalarLinearSystem system{mesh};

    convection.compute_flux_balance(field, conditions, face_flux, balance);
    convection.add_matrix_contributions(conditions, face_flux, system);
    convection.add_boundary_rhs(conditions, face_flux, system.rhs());
    std::vector<double> matrix_product(mesh.cell_count());
    system.apply_matrix(field.values(), matrix_product);

    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        require_near(matrix_product[cell_id] - system.rhs()[cell_id], balance[cell_id], test_tolerance,
                     "Assembled convection matrix and RHS do not reconstruct the direct balance.");
    }
}

void test_assembly_is_additive()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator convection{mesh};
    const cfd::ScalarBoundaryConditions conditions{
        make_uniform_conditions(1, cfd::ScalarBoundaryConditionType::Dirichlet, 7.0)};
    cfd::FaceFluxField face_flux{mesh.face_count()};
    const cfd::Index internal_id{internal_face_id(mesh)};
    const cfd::Index boundary_id{first_boundary_face_id(mesh)};
    const cfd::FaceAdjacency &internal_adjacency{mesh.face_adjacencies()[internal_id]};
    const cfd::Index boundary_owner{mesh.face_adjacencies()[boundary_id].owner};
    face_flux[internal_id] = 3.0;
    face_flux[boundary_id] = -2.0;

    cfd::ScalarLinearSystem system{mesh};
    system.diagonal()[internal_adjacency.owner] = 10.0;
    system.diagonal()[internal_adjacency.neighbor] = 20.0;
    system.owner_neighbor_coefficients()[internal_id] = 30.0;
    system.neighbor_owner_coefficients()[internal_id] = 40.0;
    system.rhs()[boundary_owner] = 50.0;

    convection.add_matrix_contributions(conditions, face_flux, system);
    convection.add_boundary_rhs(conditions, face_flux, system.rhs());

    require_near(system.diagonal()[internal_adjacency.owner], 13.0, 0.0,
                 "Convection assembly overwrote the existing owner diagonal.");
    require_near(system.diagonal()[internal_adjacency.neighbor], 20.0, 0.0,
                 "Convection assembly changed the neighbor diagonal for positive flux.");
    require_near(system.owner_neighbor_coefficients()[internal_id], 30.0, 0.0,
                 "Convection assembly changed A(owner,neighbor) for positive flux.");
    require_near(system.neighbor_owner_coefficients()[internal_id], 37.0, 0.0,
                 "Convection assembly overwrote A(neighbor,owner).");
    require_near(system.rhs()[boundary_owner], 64.0, 0.0, "Convection boundary assembly overwrote the existing RHS.");
}

void test_linear_internal_interpolation_is_geometry_aware()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_nonuniform_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator convection{mesh, cfd::ScalarConvectionScheme::Linear};
    const cfd::ScalarBoundaryConditions conditions{
        make_uniform_conditions(1, cfd::ScalarBoundaryConditionType::Neumann, 0.0)};
    const cfd::Index face_id{internal_face_id(mesh)};
    const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};
    cfd::CellScalarField field{mesh.cell_count()};
    field[adjacency.owner] = 2.0;
    field[adjacency.neighbor] = 8.0;
    cfd::FaceFluxField face_flux{mesh.face_count()};
    face_flux[face_id] = 3.0;
    cfd::CellScalarField balance{mesh.cell_count()};
    cfd::ScalarLinearSystem system{mesh};

    convection.compute_flux_balance(field, conditions, face_flux, balance);
    convection.add_matrix_contributions(conditions, face_flux, system);

    require_near(balance[adjacency.owner], 12.0, test_tolerance,
                 "Linear interpolation gave an incorrect owner balance for lambda=1/3.");
    require_near(balance[adjacency.neighbor], -12.0, test_tolerance,
                 "Linear interpolation gave an incorrect neighbor balance for lambda=1/3.");
    require_near(balance[adjacency.owner] + balance[adjacency.neighbor], 0.0, test_tolerance,
                 "Linear internal-face contributions are not conservative.");
    require_near(system.diagonal()[adjacency.owner], 2.0, test_tolerance,
                 "Linear interpolation gave an incorrect owner diagonal for lambda=1/3.");
    require_near(system.owner_neighbor_coefficients()[face_id], 1.0, test_tolerance,
                 "Linear interpolation gave an incorrect A(owner,neighbor) for lambda=1/3.");
    require_near(system.neighbor_owner_coefficients()[face_id], -2.0, test_tolerance,
                 "Linear interpolation gave an incorrect A(neighbor,owner) for lambda=1/3.");
    require_near(system.diagonal()[adjacency.neighbor], -1.0, test_tolerance,
                 "Linear interpolation gave an incorrect neighbor diagonal for lambda=1/3.");
}

void test_linear_dirichlet_boundary_for_both_flux_directions()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_quadrilateral_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator convection{mesh, cfd::ScalarConvectionScheme::Linear};
    const cfd::ScalarBoundaryConditions conditions{
        make_uniform_conditions(1, cfd::ScalarBoundaryConditionType::Dirichlet, 7.0)};
    const cfd::CellScalarField field{mesh.cell_count(), 3.0};
    const cfd::Index face_id{first_boundary_face_id(mesh)};

    for (const double carrier_flux : {-2.0, 2.0})
    {
        cfd::FaceFluxField face_flux{mesh.face_count()};
        face_flux[face_id] = carrier_flux;
        cfd::CellScalarField balance{mesh.cell_count()};
        cfd::ScalarLinearSystem system{mesh};
        convection.compute_flux_balance(field, conditions, face_flux, balance);
        convection.add_matrix_contributions(conditions, face_flux, system);
        convection.add_boundary_rhs(conditions, face_flux, system.rhs());

        require_near(balance[0], carrier_flux * 7.0, 0.0,
                     "Linear Dirichlet boundary did not use the prescribed value.");
        require_near(system.diagonal()[0], 0.0, 0.0, "Linear Dirichlet boundary incorrectly changed the diagonal.");
        require_near(system.rhs()[0], -carrier_flux * 7.0, 0.0,
                     "Linear Dirichlet boundary gave an incorrect RHS contribution.");
        require_near(system.diagonal()[0] * field[0] - system.rhs()[0], balance[0], 0.0,
                     "Linear Dirichlet boundary assembly does not match its direct balance.");
    }
}

void test_linear_zero_gradient_boundary_for_both_flux_directions()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_quadrilateral_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator convection{mesh, cfd::ScalarConvectionScheme::Linear};
    const cfd::ScalarBoundaryConditions conditions{
        make_uniform_conditions(1, cfd::ScalarBoundaryConditionType::Neumann, 0.0)};
    const cfd::CellScalarField field{mesh.cell_count(), 3.0};
    const cfd::Index face_id{first_boundary_face_id(mesh)};

    for (const double carrier_flux : {-2.0, 2.0})
    {
        cfd::FaceFluxField face_flux{mesh.face_count()};
        face_flux[face_id] = carrier_flux;
        cfd::CellScalarField balance{mesh.cell_count()};
        cfd::ScalarLinearSystem system{mesh};
        convection.compute_flux_balance(field, conditions, face_flux, balance);
        convection.add_matrix_contributions(conditions, face_flux, system);
        convection.add_boundary_rhs(conditions, face_flux, system.rhs());

        require_near(balance[0], carrier_flux * 3.0, 0.0, "Linear zeroGradient boundary did not use the owner value.");
        require_near(system.diagonal()[0], carrier_flux, 0.0,
                     "Linear zeroGradient boundary gave an incorrect diagonal.");
        require_near(system.rhs()[0], 0.0, 0.0, "Linear zeroGradient boundary changed the RHS.");
    }
}

void test_linear_nonzero_neumann_boundary_for_both_flux_directions()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_quadrilateral_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator convection{mesh, cfd::ScalarConvectionScheme::Linear};
    const cfd::ScalarBoundaryConditions conditions{
        make_uniform_conditions(1, cfd::ScalarBoundaryConditionType::Neumann, 4.0)};
    const cfd::CellScalarField field{mesh.cell_count(), 3.0};
    const cfd::Index face_id{first_boundary_face_id(mesh)};

    for (const double carrier_flux : {-2.0, 2.0})
    {
        cfd::FaceFluxField face_flux{mesh.face_count()};
        face_flux[face_id] = carrier_flux;
        cfd::CellScalarField balance{mesh.cell_count()};
        cfd::ScalarLinearSystem system{mesh};
        convection.compute_flux_balance(field, conditions, face_flux, balance);
        convection.add_matrix_contributions(conditions, face_flux, system);
        convection.add_boundary_rhs(conditions, face_flux, system.rhs());

        require_near(balance[0], carrier_flux * 5.0, test_tolerance,
                     "Linear Neumann boundary reconstructed an incorrect face value.");
        require_near(system.diagonal()[0], carrier_flux, 0.0, "Linear Neumann boundary gave an incorrect diagonal.");
        require_near(system.rhs()[0], -carrier_flux * 2.0, test_tolerance,
                     "Linear Neumann boundary gave an incorrect RHS contribution.");
        require_near(system.diagonal()[0] * field[0] - system.rhs()[0], balance[0], test_tolerance,
                     "Linear Neumann boundary assembly does not match its direct balance.");
    }
}

void test_linear_constant_field_preservation()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_cell_four_boundary_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator convection{mesh, cfd::ScalarConvectionScheme::Linear};
    std::vector<cfd::ScalarBoundaryCondition> condition_values{
        {cfd::ScalarBoundaryConditionType::Dirichlet, 4.0},
        {cfd::ScalarBoundaryConditionType::Neumann, 0.0},
        {cfd::ScalarBoundaryConditionType::Dirichlet, 4.0},
        {cfd::ScalarBoundaryConditionType::Neumann, 0.0},
    };
    const cfd::ScalarBoundaryConditions conditions{4, std::move(condition_values)};
    const cfd::CellScalarField field{mesh.cell_count(), 4.0};
    cfd::FaceFluxField face_flux{mesh.face_count()};
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "left"), -1.0);
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "right"), 1.0);
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "bottom"), -2.0);
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "top"), 2.0);
    cfd::CellScalarField balance{mesh.cell_count()};

    convection.compute_flux_balance(field, conditions, face_flux, balance);

    require_near(balance[0], 0.0, test_tolerance, "Linear convection did not preserve a compatible constant field.");
}

void test_linear_assembly_matches_direct_balance()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_six_boundary_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator convection{mesh, cfd::ScalarConvectionScheme::Linear};
    std::vector<cfd::ScalarBoundaryCondition> condition_values{
        {cfd::ScalarBoundaryConditionType::Dirichlet, 1.2},  {cfd::ScalarBoundaryConditionType::Neumann, -0.4},
        {cfd::ScalarBoundaryConditionType::Dirichlet, -2.0}, {cfd::ScalarBoundaryConditionType::Neumann, 0.7},
        {cfd::ScalarBoundaryConditionType::Dirichlet, 3.0},  {cfd::ScalarBoundaryConditionType::Neumann, -1.1},
    };
    const cfd::ScalarBoundaryConditions conditions{6, std::move(condition_values)};
    cfd::CellScalarField field{mesh.cell_count()};
    field[0] = 2.3;
    field[1] = -0.8;
    cfd::FaceFluxField face_flux{mesh.face_count()};
    face_flux[internal_face_id(mesh)] = -1.7;
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "bottom_0"), -0.7);
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "bottom_1"), -1.1);
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "right"), 0.6);
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "top_1"), 0.9);
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "left"), -0.3);
    cfd::CellScalarField balance{mesh.cell_count()};
    cfd::ScalarLinearSystem system{mesh};

    convection.compute_flux_balance(field, conditions, face_flux, balance);
    convection.add_matrix_contributions(conditions, face_flux, system);
    convection.add_boundary_rhs(conditions, face_flux, system.rhs());
    std::vector<double> matrix_product(mesh.cell_count());
    system.apply_matrix(field.values(), matrix_product);

    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        require_near(matrix_product[cell_id] - system.rhs()[cell_id], balance[cell_id], test_tolerance,
                     "Linear convection assembly does not reconstruct the direct balance.");
    }
}

void test_linear_assembly_is_additive()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_nonuniform_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator convection{mesh, cfd::ScalarConvectionScheme::Linear};
    const cfd::ScalarBoundaryConditions conditions{
        make_uniform_conditions(1, cfd::ScalarBoundaryConditionType::Dirichlet, 7.0)};
    cfd::FaceFluxField face_flux{mesh.face_count()};
    const cfd::Index internal_id{internal_face_id(mesh)};
    const cfd::Index boundary_id{first_boundary_face_id(mesh)};
    const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[internal_id]};
    const cfd::Index boundary_owner{mesh.face_adjacencies()[boundary_id].owner};
    face_flux[internal_id] = 3.0;
    face_flux[boundary_id] = 2.0;

    cfd::ScalarLinearSystem system{mesh};
    system.diagonal()[adjacency.owner] = 10.0;
    system.diagonal()[adjacency.neighbor] = 20.0;
    system.owner_neighbor_coefficients()[internal_id] = 30.0;
    system.neighbor_owner_coefficients()[internal_id] = 40.0;
    system.rhs()[boundary_owner] = 50.0;

    convection.add_matrix_contributions(conditions, face_flux, system);
    convection.add_boundary_rhs(conditions, face_flux, system.rhs());

    require_near(system.diagonal()[adjacency.owner], 12.0, test_tolerance,
                 "Linear convection overwrote the existing owner diagonal.");
    require_near(system.diagonal()[adjacency.neighbor], 19.0, test_tolerance,
                 "Linear convection overwrote the existing neighbor diagonal.");
    require_near(system.owner_neighbor_coefficients()[internal_id], 31.0, test_tolerance,
                 "Linear convection overwrote A(owner,neighbor).");
    require_near(system.neighbor_owner_coefficients()[internal_id], 38.0, test_tolerance,
                 "Linear convection overwrote A(neighbor,owner).");
    require_near(system.rhs()[boundary_owner], 36.0, test_tolerance, "Linear convection overwrote the existing RHS.");
}

void test_linear_validation_and_input_immutability()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator convection{mesh, cfd::ScalarConvectionScheme::Linear};
    cfd::CellScalarField field{mesh.cell_count()};
    field[0] = 1.25;
    field[1] = -0.75;
    cfd::FaceFluxField face_flux{mesh.face_count()};
    face_flux[internal_face_id(mesh)] = -2.5;
    const cfd::ScalarBoundaryConditions conditions{
        make_uniform_conditions(1, cfd::ScalarBoundaryConditionType::Neumann, 0.5)};
    const cfd::CellScalarField original_field{field};
    const cfd::FaceFluxField original_flux{face_flux};
    cfd::CellScalarField balance{mesh.cell_count()};
    cfd::ScalarLinearSystem system{mesh};

    require_throws<std::invalid_argument>(
        [&convection, &field, &conditions, &balance, &mesh]() {
            const cfd::FaceFluxField wrong_flux{mesh.face_count() + 1};
            convection.compute_flux_balance(field, conditions, wrong_flux, balance);
        },
        "Linear convection accepted face fluxes with incorrect cardinality.");
    require_throws<std::invalid_argument>(
        [&convection, &field, &conditions, &face_flux]() {
            convection.compute_flux_balance(field, conditions, face_flux, field);
        },
        "Linear convection accepted an output alias of its input field.");

    convection.compute_flux_balance(field, conditions, face_flux, balance);
    convection.add_matrix_contributions(conditions, face_flux, system);
    convection.add_boundary_rhs(conditions, face_flux, system.rhs());

    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        require_near(field[cell_id], original_field[cell_id], 0.0, "Linear convection modified its input field.");
    }
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        require_near(face_flux[face_id], original_flux[face_id], 0.0,
                     "Linear convection modified its input face flux.");
    }
    require(conditions[0].type == cfd::ScalarBoundaryConditionType::Neumann && conditions[0].value == 0.5,
            "Linear convection modified its boundary condition.");
}

void test_hybrid_uniform_internal_switch()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator convection{mesh, cfd::ScalarConvectionScheme::Hybrid};
    const cfd::ScalarBoundaryConditions conditions{
        make_uniform_conditions(1, cfd::ScalarBoundaryConditionType::Neumann, 0.0)};
    const cfd::Index face_id{internal_face_id(mesh)};
    const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};
    cfd::CellScalarField field{mesh.cell_count()};
    field[adjacency.owner] = 2.0;
    field[adjacency.neighbor] = 8.0;
    const std::vector<double> conductances(mesh.face_count(), 1.0);

    const auto require_face_value = [&](const double carrier_flux, const double expected_face_value,
                                        const std::string &context) {
        cfd::FaceFluxField face_flux{mesh.face_count()};
        face_flux[face_id] = carrier_flux;
        cfd::CellScalarField balance{mesh.cell_count()};

        convection.compute_flux_balance(field, conditions, face_flux, balance, conductances);

        require_near(balance[adjacency.owner], carrier_flux * expected_face_value, test_tolerance,
                     context + " owner balance is incorrect.");
        require_near(balance[adjacency.neighbor], -carrier_flux * expected_face_value, test_tolerance,
                     context + " neighbor balance is incorrect.");
        require_near(balance[adjacency.owner] + balance[adjacency.neighbor], 0.0, test_tolerance,
                     context + " is not conservative.");
    };

    require_face_value(1.5, 5.0, "Positive sub-threshold Hybrid flux");
    require_face_value(-1.5, 5.0, "Negative sub-threshold Hybrid flux");
    require_face_value(2.0, 5.0, "Positive threshold-equality Hybrid flux");
    require_face_value(-2.0, 5.0, "Negative threshold-equality Hybrid flux");
    require_face_value(2.5, 2.0, "Positive super-threshold Hybrid flux");
    require_face_value(-2.5, 8.0, "Negative super-threshold Hybrid flux");
}

void test_hybrid_nonuniform_internal_switch()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_nonuniform_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator convection{mesh, cfd::ScalarConvectionScheme::Hybrid};
    const cfd::ScalarBoundaryConditions conditions{
        make_uniform_conditions(1, cfd::ScalarBoundaryConditionType::Neumann, 0.0)};
    const cfd::Index face_id{internal_face_id(mesh)};
    const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};
    cfd::CellScalarField field{mesh.cell_count()};
    field[adjacency.owner] = 2.0;
    field[adjacency.neighbor] = 8.0;
    const std::vector<double> conductances(mesh.face_count(), 1.0);

    cfd::FaceFluxField positive_flux{mesh.face_count()};
    positive_flux[face_id] = 2.5;
    cfd::CellScalarField positive_balance{mesh.cell_count()};
    convection.compute_flux_balance(field, conditions, positive_flux, positive_balance, conductances);
    require_near(positive_balance[adjacency.owner], 10.0, test_tolerance,
                 "Hybrid did not use Linear for F*lambda <= D on a nonuniform face.");

    cfd::FaceFluxField negative_flux{mesh.face_count()};
    negative_flux[face_id] = -2.0;
    cfd::CellScalarField negative_balance{mesh.cell_count()};
    convection.compute_flux_balance(field, conditions, negative_flux, negative_balance, conductances);
    require_near(negative_balance[adjacency.owner], -16.0, test_tolerance,
                 "Hybrid did not use Upwind for (-F)*(1-lambda) > D on a nonuniform face.");
    require_near(negative_balance[adjacency.owner] + negative_balance[adjacency.neighbor], 0.0, test_tolerance,
                 "Nonuniform Hybrid internal-face contributions are not conservative.");
}

void test_hybrid_boundary_switch()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_quadrilateral_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator convection{mesh, cfd::ScalarConvectionScheme::Hybrid};
    const cfd::CellScalarField field{mesh.cell_count(), 3.0};
    const cfd::Index face_id{first_boundary_face_id(mesh)};
    const std::vector<double> conductances(mesh.face_count(), 1.0);

    const auto require_boundary_case = [&](const cfd::ScalarBoundaryConditionType type, const double boundary_value,
                                           const double carrier_flux, const double expected_face_value,
                                           const std::string &context) {
        const cfd::ScalarBoundaryConditions conditions{make_uniform_conditions(1, type, boundary_value)};
        cfd::FaceFluxField face_flux{mesh.face_count()};
        face_flux[face_id] = carrier_flux;
        cfd::CellScalarField balance{mesh.cell_count()};
        cfd::ScalarLinearSystem system{mesh};
        std::vector<double> matrix_product(mesh.cell_count());

        convection.compute_flux_balance(field, conditions, face_flux, balance, conductances);
        convection.add_matrix_contributions(conditions, face_flux, system, conductances);
        convection.add_boundary_rhs(conditions, face_flux, system.rhs(), conductances);
        system.apply_matrix(field.values(), matrix_product);

        require_near(balance[0], carrier_flux * expected_face_value, test_tolerance,
                     context + " direct balance is incorrect.");
        require_near(matrix_product[0] - system.rhs()[0], balance[0], test_tolerance,
                     context + " assembly does not reconstruct the direct balance.");
    };

    require_boundary_case(cfd::ScalarBoundaryConditionType::Dirichlet, 7.0, -2.0, 7.0, "Hybrid Dirichlet inflow");
    require_boundary_case(cfd::ScalarBoundaryConditionType::Dirichlet, 7.0, 0.5, 7.0,
                          "Hybrid Dirichlet sub-threshold outflow");
    require_boundary_case(cfd::ScalarBoundaryConditionType::Dirichlet, 7.0, 1.0, 7.0,
                          "Hybrid Dirichlet threshold-equality outflow");
    require_boundary_case(cfd::ScalarBoundaryConditionType::Dirichlet, 7.0, 2.0, 3.0,
                          "Hybrid Dirichlet super-threshold outflow");
    require_boundary_case(cfd::ScalarBoundaryConditionType::Neumann, 4.0, -2.0, 5.0, "Hybrid Neumann inflow");
    require_boundary_case(cfd::ScalarBoundaryConditionType::Neumann, 4.0, 0.5, 5.0,
                          "Hybrid Neumann sub-threshold outflow");
    require_boundary_case(cfd::ScalarBoundaryConditionType::Neumann, 4.0, 1.0, 5.0,
                          "Hybrid Neumann threshold-equality outflow");
    require_boundary_case(cfd::ScalarBoundaryConditionType::Neumann, 4.0, 2.0, 3.0,
                          "Hybrid Neumann super-threshold outflow");
    require_boundary_case(cfd::ScalarBoundaryConditionType::Dirichlet, 7.0, 0.0, 7.0, "Hybrid zero boundary flux");
}

void test_hybrid_assembly_matches_direct_balance()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_six_boundary_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator convection{mesh, cfd::ScalarConvectionScheme::Hybrid};
    std::vector<cfd::ScalarBoundaryCondition> condition_values{
        {cfd::ScalarBoundaryConditionType::Dirichlet, 1.2},  {cfd::ScalarBoundaryConditionType::Neumann, -0.4},
        {cfd::ScalarBoundaryConditionType::Dirichlet, -2.0}, {cfd::ScalarBoundaryConditionType::Neumann, 0.7},
        {cfd::ScalarBoundaryConditionType::Dirichlet, 3.0},  {cfd::ScalarBoundaryConditionType::Neumann, -1.1},
    };
    const cfd::ScalarBoundaryConditions conditions{6, std::move(condition_values)};
    cfd::CellScalarField field{mesh.cell_count()};
    field[0] = 2.3;
    field[1] = -0.8;
    cfd::FaceFluxField face_flux{mesh.face_count()};
    const cfd::Index internal_id{internal_face_id(mesh)};
    face_flux[internal_id] = -1.7;
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "bottom_0"), -0.7);
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "bottom_1"), -1.1);
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "right"), 1.5);
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "top_1"), 0.9);
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "left"), -0.3);
    std::vector<double> conductances(mesh.face_count(), 1.0);
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (mesh.face_boundary_ids()[face_id] == find_boundary_id(mesh, "right"))
        {
            conductances[face_id] = 0.5;
        }
    }
    conductances[internal_id] = 2.0;
    cfd::CellScalarField balance{mesh.cell_count()};
    cfd::ScalarLinearSystem system{mesh};
    std::vector<double> matrix_product(mesh.cell_count());

    convection.compute_flux_balance(field, conditions, face_flux, balance, conductances);
    convection.add_matrix_contributions(conditions, face_flux, system, conductances);
    convection.add_boundary_rhs(conditions, face_flux, system.rhs(), conductances);
    system.apply_matrix(field.values(), matrix_product);

    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        require_near(matrix_product[cell_id] - system.rhs()[cell_id], balance[cell_id], test_tolerance,
                     "Hybrid convection assembly does not reconstruct the direct balance.");
    }
}

void test_hybrid_face_classification_diagnostic()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_six_boundary_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator hybrid{mesh, cfd::ScalarConvectionScheme::Hybrid};
    cfd::FaceFluxField face_flux{mesh.face_count()};
    const cfd::Index internal_id{internal_face_id(mesh)};
    face_flux[internal_id] = 3.0;
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "bottom_0"), -2.0);
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "bottom_1"), 0.0);
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "right"), 2.0);
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "top_1"), 0.5);
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "top_0"), 2.0);
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "left"), 1.0);
    const std::vector<double> conductances(mesh.face_count(), 1.0);

    const cfd::HybridConvectionFaceCounts upwind_internal{hybrid.classify_hybrid_faces(face_flux, conductances)};
    require(upwind_internal.internal_linear_faces == 0 && upwind_internal.internal_upwind_faces == 1,
            "Hybrid face classification gave incorrect internal Upwind counts.");
    require(upwind_internal.boundary_linear_faces == 4 && upwind_internal.boundary_upwind_faces == 2,
            "Hybrid face classification gave incorrect boundary counts.");

    face_flux[internal_id] = 2.0;
    const cfd::HybridConvectionFaceCounts linear_internal{hybrid.classify_hybrid_faces(face_flux, conductances)};
    require(linear_internal.internal_linear_faces == 1 && linear_internal.internal_upwind_faces == 0,
            "Hybrid face classification did not classify threshold equality as Linear.");
    require(linear_internal.boundary_linear_faces == upwind_internal.boundary_linear_faces &&
                linear_internal.boundary_upwind_faces == upwind_internal.boundary_upwind_faces,
            "Changing an internal flux changed Hybrid boundary counts.");

    require_throws<std::invalid_argument>(
        [&]() {
            const cfd::FaceFluxField wrong_flux{mesh.face_count() + 1};
            static_cast<void>(hybrid.classify_hybrid_faces(wrong_flux, conductances));
        },
        "Hybrid face classification accepted an incorrect flux cardinality.");
    require_throws<std::invalid_argument>(
        [&]() {
            const std::vector<double> wrong_conductances(mesh.face_count() - 1, 1.0);
            static_cast<void>(hybrid.classify_hybrid_faces(face_flux, wrong_conductances));
        },
        "Hybrid face classification accepted an incorrect conductance cardinality.");
    require_throws<std::invalid_argument>(
        [&]() {
            std::vector<double> invalid_conductances(mesh.face_count(), 1.0);
            invalid_conductances.back() = 0.0;
            static_cast<void>(hybrid.classify_hybrid_faces(face_flux, invalid_conductances));
        },
        "Hybrid face classification accepted a non-positive conductance.");

    const cfd::ScalarConvectionOperator linear{mesh, cfd::ScalarConvectionScheme::Linear};
    require_throws<std::invalid_argument>(
        [&]() { static_cast<void>(linear.classify_hybrid_faces(face_flux, conductances)); },
        "A non-Hybrid convection operator accepted Hybrid face classification.");
}

void test_hybrid_conductance_validation()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator convection{mesh, cfd::ScalarConvectionScheme::Hybrid};
    const cfd::ScalarBoundaryConditions conditions{
        make_uniform_conditions(1, cfd::ScalarBoundaryConditionType::Neumann, 0.0)};
    const cfd::CellScalarField field{mesh.cell_count(), 1.0};
    const cfd::FaceFluxField face_flux{mesh.face_count()};
    cfd::CellScalarField balance{mesh.cell_count()};
    cfd::ScalarLinearSystem system{mesh};

    require_throws<std::invalid_argument>(
        [&]() { convection.compute_flux_balance(field, conditions, face_flux, balance); },
        "Hybrid convection accepted absent diffusion conductances.");
    require_throws<std::invalid_argument>(
        [&]() {
            const std::vector<double> wrong_conductances(mesh.face_count() - 1, 1.0);
            convection.add_matrix_contributions(conditions, face_flux, system, wrong_conductances);
        },
        "Hybrid convection accepted an incorrect conductance cardinality.");

    const auto require_invalid_conductance = [&](const double invalid_conductance) {
        std::vector<double> conductances(mesh.face_count(), 1.0);
        conductances[first_boundary_face_id(mesh)] = invalid_conductance;
        require_throws<std::invalid_argument>(
            [&]() { convection.add_boundary_rhs(conditions, face_flux, system.rhs(), conductances); },
            "Hybrid convection accepted an invalid diffusion conductance.");
    };
    require_invalid_conductance(0.0);
    require_invalid_conductance(-1.0);
    require_invalid_conductance(std::numeric_limits<double>::quiet_NaN());
    require_invalid_conductance(std::numeric_limits<double>::infinity());
    require_invalid_conductance(-std::numeric_limits<double>::infinity());
}

void test_rejects_unsupported_scheme()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_quadrilateral_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};

    require_throws<std::invalid_argument>(
        [&mesh]() {
            static_cast<void>(cfd::ScalarConvectionOperator{
                mesh,
                // NOLINTNEXTLINE(clang-analyzer-optin.core.EnumCastOutOfRange)
                static_cast<cfd::ScalarConvectionScheme>(255),
            });
        },
        "Scalar convection accepted an unsupported scheme.");
}

void test_rejects_incompatible_inputs_and_aliasing()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator convection{mesh};
    cfd::CellScalarField field{mesh.cell_count()};
    const cfd::ScalarBoundaryConditions conditions{
        make_uniform_conditions(1, cfd::ScalarBoundaryConditionType::Neumann, 0.0)};
    const cfd::FaceFluxField face_flux{mesh.face_count()};
    cfd::CellScalarField balance{mesh.cell_count()};

    require_throws<std::invalid_argument>(
        [&convection, &conditions, &face_flux, &balance]() {
            const cfd::CellScalarField wrong_field{3};
            convection.compute_flux_balance(wrong_field, conditions, face_flux, balance);
        },
        "Scalar convection accepted a field with incorrect cardinality.");
    require_throws<std::invalid_argument>(
        [&convection, &field, &conditions, &face_flux]() {
            cfd::CellScalarField wrong_output{3};
            convection.compute_flux_balance(field, conditions, face_flux, wrong_output);
        },
        "Scalar convection accepted an output with incorrect cardinality.");
    require_throws<std::invalid_argument>(
        [&convection, &field, &conditions, &balance, &mesh]() {
            const cfd::FaceFluxField wrong_flux{mesh.face_count() + 1};
            convection.compute_flux_balance(field, conditions, wrong_flux, balance);
        },
        "Scalar convection accepted face fluxes with incorrect cardinality.");
    require_throws<std::invalid_argument>(
        [&convection, &field, &face_flux, &balance]() {
            const cfd::ScalarBoundaryConditions wrong_conditions{
                make_uniform_conditions(2, cfd::ScalarBoundaryConditionType::Neumann, 0.0)};
            convection.compute_flux_balance(field, wrong_conditions, face_flux, balance);
        },
        "Scalar convection accepted boundary conditions with incorrect cardinality.");
    require_throws<std::invalid_argument>(
        [&convection, &field, &conditions, &face_flux]() {
            convection.compute_flux_balance(field, conditions, face_flux, field);
        },
        "Scalar convection accepted an output alias of its input field.");
    require_throws<std::invalid_argument>(
        [&convection, &conditions, &face_flux]() {
            cfd::MeshBuildResult other_build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
            cfd::ScalarLinearSystem other_system{other_build_result.mesh};
            convection.add_matrix_contributions(conditions, face_flux, other_system);
        },
        "Scalar convection accepted a system referencing a different Mesh instance.");
    require_throws<std::invalid_argument>(
        [&convection, &conditions, &face_flux, &mesh]() {
            std::vector<double> wrong_rhs(mesh.cell_count() + 1);
            convection.add_boundary_rhs(conditions, face_flux, wrong_rhs);
        },
        "Scalar convection accepted a boundary RHS with incorrect cardinality.");
}

void test_does_not_mutate_inputs()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator convection{mesh};
    cfd::CellScalarField field{mesh.cell_count()};
    field[0] = 1.25;
    field[1] = -0.75;
    cfd::FaceFluxField face_flux{mesh.face_count()};
    face_flux[internal_face_id(mesh)] = -2.5;
    const cfd::ScalarBoundaryConditions conditions{
        make_uniform_conditions(1, cfd::ScalarBoundaryConditionType::Neumann, 0.5)};
    const cfd::CellScalarField original_field{field};
    const cfd::FaceFluxField original_flux{face_flux};
    cfd::CellScalarField balance{mesh.cell_count()};
    cfd::ScalarLinearSystem system{mesh};

    convection.compute_flux_balance(field, conditions, face_flux, balance);
    convection.add_matrix_contributions(conditions, face_flux, system);
    convection.add_boundary_rhs(conditions, face_flux, system.rhs());

    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        require_near(field[cell_id], original_field[cell_id], 0.0, "Scalar convection modified its input field.");
    }
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        require_near(face_flux[face_id], original_flux[face_id], 0.0,
                     "Scalar convection modified its input face flux.");
    }
    require(conditions[0].type == cfd::ScalarBoundaryConditionType::Neumann && conditions[0].value == 0.5,
            "Scalar convection modified its boundary condition.");
}

void test_linear_upwind_exact_linear_internal_reconstruction()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_nonuniform_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator convection{mesh, cfd::ScalarConvectionScheme::LinearUpwind};
    const cfd::ScalarBoundaryConditions conditions{
        make_uniform_conditions(1, cfd::ScalarBoundaryConditionType::Neumann, 0.0)};
    const auto cell_centers{mesh.cell_centers()};
    const cfd::Index face_id{internal_face_id(mesh)};
    const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};
    const cfd::Point2 &face_center{mesh.face_centers()[face_id]};
    cfd::CellScalarField field{mesh.cell_count()};
    cfd::CellVectorField gradient{mesh.cell_count(), {2.0, -3.0}};
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        field[cell_id] = 2.0 * cell_centers[cell_id].x - 3.0 * cell_centers[cell_id].y + 4.0;
    }
    const double exact_face_value{2.0 * face_center.x - 3.0 * face_center.y + 4.0};

    for (const double carrier_flux : {3.0, -4.0})
    {
        cfd::FaceFluxField face_flux{mesh.face_count()};
        face_flux[face_id] = carrier_flux;
        cfd::CellScalarField balance{mesh.cell_count()};

        convection.compute_flux_balance(field, conditions, face_flux, gradient, balance);

        require_near(balance[adjacency.owner], carrier_flux * exact_face_value, test_tolerance,
                     "LinearUpwind did not exactly reconstruct a linear internal-face value.");
        require_near(balance[adjacency.neighbor], -carrier_flux * exact_face_value, test_tolerance,
                     "LinearUpwind gave an incorrect conservative neighbor balance.");
    }
}

void test_linear_upwind_selects_upwind_gradient_and_adds_conservatively()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_nonuniform_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator convection{mesh, cfd::ScalarConvectionScheme::LinearUpwind};
    const cfd::Index face_id{internal_face_id(mesh)};
    const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};
    const auto cell_centers{mesh.cell_centers()};
    const cfd::Point2 &face_center{mesh.face_centers()[face_id]};
    cfd::CellVectorField gradient{mesh.cell_count()};
    gradient[adjacency.owner] = {2.0, 0.0};
    gradient[adjacency.neighbor] = {-5.0, 0.0};

    for (const double carrier_flux : {3.0, -4.0})
    {
        cfd::FaceFluxField face_flux{mesh.face_count()};
        face_flux[face_id] = carrier_flux;
        std::vector<double> rhs{10.0, -20.0};
        const cfd::Index upwind_id{carrier_flux >= 0.0 ? adjacency.owner : adjacency.neighbor};
        const double correction{gradient[upwind_id].x * (face_center.x - cell_centers[upwind_id].x) +
                                gradient[upwind_id].y * (face_center.y - cell_centers[upwind_id].y)};
        const double correction_flux{carrier_flux * correction};

        convection.add_deferred_correction_rhs(gradient, face_flux, rhs);

        require_near(rhs[adjacency.owner], 10.0 - correction_flux, test_tolerance,
                     "LinearUpwind selected the wrong upwind gradient for the owner RHS.");
        require_near(rhs[adjacency.neighbor], -20.0 + correction_flux, test_tolerance,
                     "LinearUpwind selected the wrong upwind gradient for the neighbor RHS.");
        require_near((rhs[adjacency.owner] - 10.0) + (rhs[adjacency.neighbor] + 20.0), 0.0, test_tolerance,
                     "LinearUpwind deferred correction is not conservative.");
    }

    const cfd::CellVectorField zero_gradient{mesh.cell_count()};
    cfd::FaceFluxField face_flux{mesh.face_count()};
    face_flux[face_id] = 3.0;
    std::vector<double> rhs{10.0, -20.0};
    convection.add_deferred_correction_rhs(zero_gradient, face_flux, rhs);
    require_near(rhs[0], 10.0, 0.0, "A zero gradient changed the LinearUpwind owner RHS.");
    require_near(rhs[1], -20.0, 0.0, "A zero gradient changed the LinearUpwind neighbor RHS.");
}

void test_linear_upwind_matrix_and_base_boundary_rhs_match_upwind()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_six_boundary_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator upwind{mesh, cfd::ScalarConvectionScheme::FirstOrderUpwind};
    const cfd::ScalarConvectionOperator linear_upwind{mesh, cfd::ScalarConvectionScheme::LinearUpwind};
    const cfd::ScalarBoundaryConditions conditions{
        make_uniform_conditions(mesh.boundary_groups().size(), cfd::ScalarBoundaryConditionType::Dirichlet, 7.0)};
    cfd::FaceFluxField face_flux{mesh.face_count()};
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        const double magnitude{0.5 * static_cast<double>(face_id + 1)};
        face_flux[face_id] = face_id % 2 == 0 ? magnitude : -magnitude;
    }
    cfd::ScalarLinearSystem upwind_system{mesh};
    cfd::ScalarLinearSystem linear_upwind_system{mesh};

    upwind.add_matrix_contributions(conditions, face_flux, upwind_system);
    upwind.add_boundary_rhs(conditions, face_flux, upwind_system.rhs());
    linear_upwind.add_matrix_contributions(conditions, face_flux, linear_upwind_system);
    linear_upwind.add_boundary_rhs(conditions, face_flux, linear_upwind_system.rhs());

    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        require_near(linear_upwind_system.diagonal()[cell_id], upwind_system.diagonal()[cell_id], 0.0,
                     "LinearUpwind matrix diagonal differs from FirstOrderUpwind.");
        require_near(linear_upwind_system.rhs()[cell_id], upwind_system.rhs()[cell_id], 0.0,
                     "LinearUpwind base boundary RHS differs from FirstOrderUpwind.");
    }
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        require_near(linear_upwind_system.owner_neighbor_coefficients()[face_id],
                     upwind_system.owner_neighbor_coefficients()[face_id], 0.0,
                     "LinearUpwind owner-neighbor coefficient differs from FirstOrderUpwind.");
        require_near(linear_upwind_system.neighbor_owner_coefficients()[face_id],
                     upwind_system.neighbor_owner_coefficients()[face_id], 0.0,
                     "LinearUpwind neighbor-owner coefficient differs from FirstOrderUpwind.");
    }
}

void test_linear_upwind_assembly_matches_direct_balance()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_six_boundary_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator convection{mesh, cfd::ScalarConvectionScheme::LinearUpwind};
    std::vector<cfd::ScalarBoundaryCondition> condition_values{
        {cfd::ScalarBoundaryConditionType::Dirichlet, 1.2},  {cfd::ScalarBoundaryConditionType::Neumann, -0.4},
        {cfd::ScalarBoundaryConditionType::Dirichlet, -2.0}, {cfd::ScalarBoundaryConditionType::Neumann, 0.7},
        {cfd::ScalarBoundaryConditionType::Dirichlet, 3.0},  {cfd::ScalarBoundaryConditionType::Neumann, -1.1},
    };
    const cfd::ScalarBoundaryConditions conditions{6, std::move(condition_values)};
    cfd::CellScalarField field{mesh.cell_count()};
    field[0] = 2.3;
    field[1] = -0.8;
    cfd::CellVectorField gradient{mesh.cell_count()};
    gradient[0] = {0.7, -1.1};
    gradient[1] = {-0.4, 1.3};
    cfd::FaceFluxField face_flux{mesh.face_count()};
    face_flux[internal_face_id(mesh)] = -1.7;
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "bottom_0"), -0.7);
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "bottom_1"), -1.1);
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "right"), 0.6);
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "top_1"), 0.9);
    set_boundary_flux(mesh, face_flux, find_boundary_id(mesh, "left"), -0.3);
    cfd::CellScalarField balance{mesh.cell_count()};
    cfd::ScalarLinearSystem system{mesh};

    convection.compute_flux_balance(field, conditions, face_flux, gradient, balance);
    convection.add_matrix_contributions(conditions, face_flux, system);
    convection.add_boundary_rhs(conditions, face_flux, system.rhs());
    convection.add_deferred_correction_rhs(gradient, face_flux, system.rhs());
    std::vector<double> matrix_product(mesh.cell_count());
    system.apply_matrix(field.values(), matrix_product);

    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        require_near(matrix_product[cell_id] - system.rhs()[cell_id], balance[cell_id], test_tolerance,
                     "LinearUpwind assembly does not reconstruct its direct flux balance.");
    }
}

void test_linear_upwind_boundary_treatments()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_quadrilateral_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator convection{mesh, cfd::ScalarConvectionScheme::LinearUpwind};
    const cfd::CellScalarField field{mesh.cell_count(), 3.0};
    const cfd::CellVectorField gradient{mesh.cell_count(), {2.0, 4.0}};
    const cfd::Index face_id{first_boundary_face_id(mesh)};
    const cfd::Index owner_id{mesh.face_adjacencies()[face_id].owner};
    const cfd::Point2 &cell_center{mesh.cell_centers()[owner_id]};
    const cfd::Point2 &face_center{mesh.face_centers()[face_id]};
    const cfd::Vector2 &area_vector{mesh.face_area_vectors()[face_id]};
    const double normal_distance{
        ((face_center.x - cell_center.x) * area_vector.x + (face_center.y - cell_center.y) * area_vector.y) /
        mesh.face_lengths()[face_id]};
    const double reconstruction{gradient[owner_id].x * (face_center.x - cell_center.x) +
                                gradient[owner_id].y * (face_center.y - cell_center.y)};

    const auto require_case = [&](const cfd::ScalarBoundaryConditionType type, const double condition_value,
                                  const double carrier_flux, const double expected_face_value,
                                  const double expected_diagonal, const double expected_rhs,
                                  const std::string &context) {
        const cfd::ScalarBoundaryConditions conditions{make_uniform_conditions(1, type, condition_value)};
        cfd::FaceFluxField face_flux{mesh.face_count()};
        face_flux[face_id] = carrier_flux;
        cfd::CellScalarField balance{mesh.cell_count(), 99.0};
        cfd::ScalarLinearSystem system{mesh};

        convection.compute_flux_balance(field, conditions, face_flux, gradient, balance);
        convection.add_matrix_contributions(conditions, face_flux, system);
        convection.add_boundary_rhs(conditions, face_flux, system.rhs());
        convection.add_deferred_correction_rhs(gradient, face_flux, system.rhs());

        require_near(balance[owner_id], carrier_flux * expected_face_value, test_tolerance,
                     context + " direct balance is incorrect.");
        require_near(system.diagonal()[owner_id], expected_diagonal, test_tolerance,
                     context + " diagonal is incorrect.");
        require_near(system.rhs()[owner_id], expected_rhs, test_tolerance, context + " RHS is incorrect.");
    };

    require_case(cfd::ScalarBoundaryConditionType::Dirichlet, 7.0, -2.0, 7.0, 0.0, 14.0,
                 "LinearUpwind Dirichlet inflow");
    require_case(cfd::ScalarBoundaryConditionType::Dirichlet, 1000.0, 2.0, 3.0 + reconstruction, 2.0,
                 -2.0 * reconstruction, "LinearUpwind Dirichlet outflow");
    require_case(cfd::ScalarBoundaryConditionType::Neumann, 4.0, -2.0, 3.0 + 4.0 * normal_distance, -2.0,
                 8.0 * normal_distance, "LinearUpwind Neumann inflow");
    require_case(cfd::ScalarBoundaryConditionType::Neumann, 1000.0, 2.0, 3.0 + reconstruction, 2.0,
                 -2.0 * reconstruction, "LinearUpwind Neumann outflow");
    require_case(cfd::ScalarBoundaryConditionType::Dirichlet, 1000.0, 0.0, 0.0, 0.0, 0.0, "LinearUpwind zero flux");
}

void test_barth_jespersen_constant_monotone_and_bounded_reconstruction()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_cell_four_boundary_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator unlimited{mesh, cfd::ScalarConvectionScheme::LinearUpwind,
                                                  cfd::ScalarConvectionLimiter::None};
    const cfd::ScalarConvectionOperator limited{mesh, cfd::ScalarConvectionScheme::LinearUpwind,
                                                cfd::ScalarConvectionLimiter::BarthJespersen};
    const cfd::ScalarBoundaryConditions linear_conditions{
        4,
        {
            {cfd::ScalarBoundaryConditionType::Dirichlet, 0.0},
            {cfd::ScalarBoundaryConditionType::Dirichlet, 1.0},
            {cfd::ScalarBoundaryConditionType::Dirichlet, 0.5},
            {cfd::ScalarBoundaryConditionType::Dirichlet, 0.5},
        },
    };
    const cfd::BoundaryId left_id{find_boundary_id(mesh, "left")};
    const cfd::BoundaryId right_id{find_boundary_id(mesh, "right")};
    cfd::CellScalarField field{mesh.cell_count(), 0.5};
    cfd::CellVectorField gradient{mesh.cell_count(), {4.0, 0.0}};
    cfd::FaceFluxField face_flux{mesh.face_count()};
    cfd::CellScalarField unlimited_balance{mesh.cell_count()};
    cfd::CellScalarField limited_balance{mesh.cell_count()};
    std::vector<double> limiter(mesh.cell_count(), -1.0);

    set_boundary_flux(mesh, face_flux, right_id, 2.0);
    unlimited.compute_flux_balance(field, linear_conditions, face_flux, gradient, unlimited_balance);
    limited.compute_flux_balance(field, linear_conditions, face_flux, gradient, limiter, limited_balance);
    require_near(limiter[0], 0.25, test_tolerance,
                 "Barth-Jespersen computed an incorrect positive-overshoot coefficient.");
    require_near(unlimited_balance[0], 5.0, test_tolerance,
                 "The unlimited reconstruction did not expose the positive overshoot fixture.");
    require_near(limited_balance[0], 2.0, test_tolerance,
                 "Barth-Jespersen did not clamp the outflow reconstruction to phi_max.");

    std::ranges::fill(face_flux.values(), 0.0);
    set_boundary_flux(mesh, face_flux, left_id, 2.0);
    unlimited.compute_flux_balance(field, linear_conditions, face_flux, gradient, unlimited_balance);
    limited.compute_flux_balance(field, linear_conditions, face_flux, gradient, limiter, limited_balance);
    require_near(limiter[0], 0.25, test_tolerance,
                 "Barth-Jespersen computed an incorrect negative-undershoot coefficient.");
    require_near(unlimited_balance[0], -3.0, test_tolerance,
                 "The unlimited reconstruction did not expose the negative undershoot fixture.");
    require_near(limited_balance[0], 0.0, test_tolerance,
                 "Barth-Jespersen did not clamp the outflow reconstruction to phi_min.");

    gradient[0] = {1.0, 0.0};
    std::ranges::fill(face_flux.values(), 0.0);
    set_boundary_flux(mesh, face_flux, right_id, 2.0);
    limited.compute_flux_balance(field, linear_conditions, face_flux, gradient, limiter, limited_balance);
    require_near(limiter[0], 1.0, test_tolerance, "Barth-Jespersen limited a monotone linear reconstruction.");
    require_near(limited_balance[0], 2.0, test_tolerance, "Barth-Jespersen lost exact monotone linear reconstruction.");

    const cfd::ScalarBoundaryConditions constant_conditions{
        make_uniform_conditions(4, cfd::ScalarBoundaryConditionType::Neumann, 0.0)};
    field[0] = 3.0;
    gradient[0] = {};
    std::ranges::fill(face_flux.values(), 0.0);
    limited.compute_flux_balance(field, constant_conditions, face_flux, gradient, limiter, limited_balance);
    require_near(limiter[0], 1.0, 0.0, "A constant field did not retain a unit limiter coefficient.");
    require_near(limited_balance[0], 0.0, 0.0, "A constant zero-flux field produced a convective balance.");

    const cfd::ScalarBoundaryConditions extremum_conditions{
        make_uniform_conditions(4, cfd::ScalarBoundaryConditionType::Dirichlet, 0.0)};
    field[0] = 1.0;
    gradient[0] = {2.0, 0.0};
    set_boundary_flux(mesh, face_flux, right_id, 2.0);
    limited.compute_flux_balance(field, extremum_conditions, face_flux, gradient, limiter, limited_balance);
    require_near(limiter[0], 0.0, 0.0, "Barth-Jespersen did not suppress reconstruction at a local extremum.");
    require_near(limited_balance[0], 2.0, test_tolerance,
                 "A local-extremum reconstruction did not fall back to the owner value.");
}

void test_barth_jespersen_internal_flux_selection_conservation_and_assembly()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_six_boundary_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator limited{mesh, cfd::ScalarConvectionScheme::LinearUpwind,
                                                cfd::ScalarConvectionLimiter::BarthJespersen};
    const cfd::ScalarConvectionOperator upwind{mesh, cfd::ScalarConvectionScheme::FirstOrderUpwind};
    const cfd::ScalarBoundaryConditions conditions{
        6,
        {
            {cfd::ScalarBoundaryConditionType::Dirichlet, 0.5},
            {cfd::ScalarBoundaryConditionType::Dirichlet, 1.5},
            {cfd::ScalarBoundaryConditionType::Dirichlet, 2.0},
            {cfd::ScalarBoundaryConditionType::Dirichlet, 1.5},
            {cfd::ScalarBoundaryConditionType::Dirichlet, 0.5},
            {cfd::ScalarBoundaryConditionType::Dirichlet, 0.0},
        },
    };
    cfd::CellScalarField field{mesh.cell_count()};
    field[0] = 0.5;
    field[1] = 1.5;
    const cfd::CellVectorField gradient{mesh.cell_count(), {4.0, 0.0}};
    const cfd::Index face_id{internal_face_id(mesh)};

    for (const double carrier_flux : {3.0, -4.0})
    {
        cfd::FaceFluxField face_flux{mesh.face_count()};
        face_flux[face_id] = carrier_flux;
        cfd::CellScalarField balance{mesh.cell_count()};
        std::vector<double> limiter(mesh.cell_count(), -1.0);
        std::vector<double> rhs{10.0, -20.0};

        limited.compute_flux_balance(field, conditions, face_flux, gradient, limiter, balance);
        limited.add_deferred_correction_rhs(field, conditions, gradient, face_flux, limiter, rhs);

        require_near(limiter[0], 0.25, test_tolerance, "The owner received an incorrect limiter coefficient.");
        require_near(limiter[1], 0.25, test_tolerance, "The neighbor received an incorrect limiter coefficient.");
        require_near(balance[0], carrier_flux, test_tolerance,
                     "The limited internal face used an incorrect upwind reconstruction.");
        require_near(balance[1], -carrier_flux, test_tolerance, "The limited internal face is not conservative.");
        require_near((rhs[0] - 10.0) + (rhs[1] + 20.0), 0.0, test_tolerance,
                     "The limited deferred correction is not conservative or additive.");

        cfd::ScalarLinearSystem limited_system{mesh};
        cfd::ScalarLinearSystem upwind_system{mesh};
        limited.add_matrix_contributions(conditions, face_flux, limited_system);
        upwind.add_matrix_contributions(conditions, face_flux, upwind_system);
        limited.add_boundary_rhs(conditions, face_flux, limited_system.rhs());
        limited.add_deferred_correction_rhs(field, conditions, gradient, face_flux, limiter, limited_system.rhs());
        std::vector<double> matrix_product(mesh.cell_count());
        limited_system.apply_matrix(field.values(), matrix_product);

        for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
        {
            require_near(limited_system.diagonal()[cell_id], upwind_system.diagonal()[cell_id], 0.0,
                         "Barth-Jespersen changed the FirstOrderUpwind matrix diagonal.");
            require_near(matrix_product[cell_id] - limited_system.rhs()[cell_id], balance[cell_id], test_tolerance,
                         "Limited direct balance differs from A*phi-rhs.");
        }
        for (cfd::Index current_face = 0; current_face < mesh.face_count(); ++current_face)
        {
            require_near(limited_system.owner_neighbor_coefficients()[current_face],
                         upwind_system.owner_neighbor_coefficients()[current_face], 0.0,
                         "Barth-Jespersen changed an owner-neighbor matrix coefficient.");
            require_near(limited_system.neighbor_owner_coefficients()[current_face],
                         upwind_system.neighbor_owner_coefficients()[current_face], 0.0,
                         "Barth-Jespersen changed a neighbor-owner matrix coefficient.");
        }
    }
}

void test_barth_jespersen_boundary_bounds_and_flow_directions()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_cell_four_boundary_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator limited{mesh, cfd::ScalarConvectionScheme::LinearUpwind,
                                                cfd::ScalarConvectionLimiter::BarthJespersen};
    const cfd::CellScalarField field{mesh.cell_count(), 0.5};
    const cfd::CellVectorField gradient{mesh.cell_count(), {4.0, 0.0}};
    std::vector<double> limiter(mesh.cell_count(), -1.0);
    cfd::CellScalarField balance{mesh.cell_count()};
    cfd::FaceFluxField face_flux{mesh.face_count()};
    const cfd::BoundaryId left_id{find_boundary_id(mesh, "left")};
    const cfd::BoundaryId right_id{find_boundary_id(mesh, "right")};

    const cfd::ScalarBoundaryConditions dirichlet_conditions{
        4,
        {
            {cfd::ScalarBoundaryConditionType::Dirichlet, 0.0},
            {cfd::ScalarBoundaryConditionType::Dirichlet, 1.0},
            {cfd::ScalarBoundaryConditionType::Dirichlet, 0.5},
            {cfd::ScalarBoundaryConditionType::Dirichlet, 0.5},
        },
    };
    set_boundary_flux(mesh, face_flux, left_id, -2.0);
    limited.compute_flux_balance(field, dirichlet_conditions, face_flux, gradient, limiter, balance);
    require_near(balance[0], 0.0, 0.0, "Barth-Jespersen changed Dirichlet inflow treatment.");
    require_near(limiter[0], 0.25, test_tolerance, "Dirichlet face values did not participate in limiter bounds.");

    const cfd::ScalarBoundaryConditions neumann_conditions{
        4,
        {
            {cfd::ScalarBoundaryConditionType::Neumann, -1.0},
            {cfd::ScalarBoundaryConditionType::Neumann, 1.0},
            {cfd::ScalarBoundaryConditionType::Neumann, 0.0},
            {cfd::ScalarBoundaryConditionType::Neumann, 0.0},
        },
    };
    std::ranges::fill(face_flux.values(), 0.0);
    set_boundary_flux(mesh, face_flux, right_id, 2.0);
    limited.compute_flux_balance(field, neumann_conditions, face_flux, gradient, limiter, balance);
    require_near(limiter[0], 0.25, test_tolerance, "Neumann closure values did not participate in limiter bounds.");
    require_near(balance[0], 2.0, test_tolerance, "Barth-Jespersen gave an incorrect Neumann outflow.");

    std::ranges::fill(face_flux.values(), 0.0);
    set_boundary_flux(mesh, face_flux, left_id, -2.0);
    limited.compute_flux_balance(field, neumann_conditions, face_flux, gradient, limiter, balance);
    require_near(balance[0], 0.0, test_tolerance, "Barth-Jespersen changed Neumann inflow treatment.");

    std::ranges::fill(face_flux.values(), 0.0);
    limited.compute_flux_balance(field, neumann_conditions, face_flux, gradient, limiter, balance);
    require_near(limiter[0], 0.25, test_tolerance, "Zero-flux faces did not participate in spatial limiting.");
    require_near(balance[0], 0.0, 0.0, "Zero face flux produced a limited convective contribution.");
}

void test_linear_upwind_validation_noop_and_input_immutability()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator convection{mesh, cfd::ScalarConvectionScheme::LinearUpwind};
    const cfd::ScalarBoundaryConditions conditions{
        make_uniform_conditions(1, cfd::ScalarBoundaryConditionType::Neumann, 0.0)};
    cfd::CellScalarField field{mesh.cell_count(), 2.0};
    cfd::CellVectorField gradient{mesh.cell_count(), {0.5, -0.25}};
    cfd::FaceFluxField face_flux{mesh.face_count()};
    face_flux[internal_face_id(mesh)] = 2.0;
    cfd::CellScalarField balance{mesh.cell_count()};
    std::vector<double> rhs(mesh.cell_count());
    const cfd::CellScalarField field_before{field};
    const cfd::CellVectorField gradient_before{gradient};
    const cfd::FaceFluxField face_flux_before{face_flux};

    require_throws<std::invalid_argument>(
        [&]() { convection.compute_flux_balance(field, conditions, face_flux, balance); },
        "LinearUpwind accepted the flux-balance overload without a gradient.");
    require_throws<std::invalid_argument>(
        [&]() {
            const cfd::CellVectorField wrong_gradient{mesh.cell_count() + 1};
            convection.add_deferred_correction_rhs(wrong_gradient, face_flux, rhs);
        },
        "LinearUpwind accepted an incorrect gradient cardinality.");
    require_throws<std::invalid_argument>(
        [&]() {
            const cfd::FaceFluxField wrong_flux{mesh.face_count() + 1};
            convection.add_deferred_correction_rhs(gradient, wrong_flux, rhs);
        },
        "LinearUpwind accepted an incorrect face-flux cardinality.");
    require_throws<std::invalid_argument>(
        [&]() {
            std::vector<double> wrong_rhs(mesh.cell_count() + 1);
            convection.add_deferred_correction_rhs(gradient, face_flux, wrong_rhs);
        },
        "LinearUpwind accepted an incorrect RHS cardinality.");
    require_throws<std::invalid_argument>(
        [&]() {
            const cfd::CellVectorField wrong_gradient{mesh.cell_count() + 1};
            convection.compute_flux_balance(field, conditions, face_flux, wrong_gradient, balance);
        },
        "LinearUpwind direct balance accepted an incorrect gradient cardinality.");
    require_throws<std::invalid_argument>(
        [&]() { convection.compute_flux_balance(field, conditions, face_flux, gradient, field); },
        "LinearUpwind direct balance accepted an aliased output.");

    convection.compute_flux_balance(field, conditions, face_flux, gradient, balance);
    convection.add_deferred_correction_rhs(gradient, face_flux, rhs);
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        require_near(field[cell_id], field_before[cell_id], 0.0, "LinearUpwind modified its scalar input.");
        require_near(gradient[cell_id].x, gradient_before[cell_id].x, 0.0,
                     "LinearUpwind modified an input gradient x component.");
        require_near(gradient[cell_id].y, gradient_before[cell_id].y, 0.0,
                     "LinearUpwind modified an input gradient y component.");
    }
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        require_near(face_flux[face_id], face_flux_before[face_id], 0.0, "LinearUpwind modified an input face flux.");
    }

    for (const cfd::ScalarConvectionScheme scheme :
         {cfd::ScalarConvectionScheme::FirstOrderUpwind, cfd::ScalarConvectionScheme::Linear,
          cfd::ScalarConvectionScheme::Hybrid})
    {
        const cfd::ScalarConvectionOperator other_scheme{mesh, scheme};
        std::vector<double> unchanged_rhs(mesh.cell_count(), 17.0);
        other_scheme.add_deferred_correction_rhs(gradient, face_flux, unchanged_rhs);
        require(std::ranges::all_of(unchanged_rhs, [](const double value) { return value == 17.0; }),
                "A non-LinearUpwind scheme changed the deferred-correction RHS.");
    }
}

void test_barth_jespersen_validation_none_compatibility_and_input_immutability()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_cell_four_boundary_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarConvectionOperator limited{mesh, cfd::ScalarConvectionScheme::LinearUpwind,
                                                cfd::ScalarConvectionLimiter::BarthJespersen};
    const cfd::ScalarConvectionOperator default_unlimited{mesh, cfd::ScalarConvectionScheme::LinearUpwind};
    const cfd::ScalarConvectionOperator explicit_unlimited{mesh, cfd::ScalarConvectionScheme::LinearUpwind,
                                                           cfd::ScalarConvectionLimiter::None};
    const cfd::ScalarBoundaryConditions conditions{
        make_uniform_conditions(4, cfd::ScalarBoundaryConditionType::Dirichlet, 1.0)};
    cfd::CellScalarField field{mesh.cell_count(), 0.5};
    cfd::CellVectorField gradient{mesh.cell_count(), {1.0, -0.5}};
    cfd::FaceFluxField face_flux{mesh.face_count()};
    face_flux[first_boundary_face_id(mesh)] = 2.0;
    cfd::CellScalarField default_balance{mesh.cell_count()};
    cfd::CellScalarField explicit_balance{mesh.cell_count()};
    cfd::CellScalarField limited_balance{mesh.cell_count()};
    std::vector<double> limiter(mesh.cell_count());
    std::vector<double> rhs(mesh.cell_count());

    default_unlimited.compute_flux_balance(field, conditions, face_flux, gradient, default_balance);
    explicit_unlimited.compute_flux_balance(field, conditions, face_flux, gradient, explicit_balance);
    require_near(explicit_balance[0], default_balance[0], 0.0,
                 "Explicit ScalarConvectionLimiter::None changed LinearUpwind.");

    for (const cfd::ScalarConvectionScheme scheme :
         {cfd::ScalarConvectionScheme::FirstOrderUpwind, cfd::ScalarConvectionScheme::Linear,
          cfd::ScalarConvectionScheme::Hybrid})
    {
        require_throws<std::invalid_argument>(
            [&mesh, scheme]() {
                const cfd::ScalarConvectionOperator invalid{mesh, scheme, cfd::ScalarConvectionLimiter::BarthJespersen};
            },
            "Barth-Jespersen was accepted with a non-LinearUpwind scheme.");
    }
    require_throws<std::invalid_argument>(
        [&mesh]() {
            static_cast<void>(cfd::ScalarConvectionOperator{
                mesh,
                cfd::ScalarConvectionScheme::LinearUpwind,
                // NOLINTNEXTLINE(clang-analyzer-optin.core.EnumCastOutOfRange)
                static_cast<cfd::ScalarConvectionLimiter>(255),
            });
        },
        "An unsupported scalar convection limiter was accepted.");
    require_throws<std::invalid_argument>(
        [&]() { limited.compute_flux_balance(field, conditions, face_flux, gradient, limited_balance); },
        "Limited LinearUpwind accepted the overload without a limiter workspace.");
    require_throws<std::invalid_argument>(
        [&]() {
            std::vector<double> wrong_workspace(mesh.cell_count() + 1);
            limited.compute_flux_balance(field, conditions, face_flux, gradient, wrong_workspace, limited_balance);
        },
        "Barth-Jespersen accepted an incorrect workspace cardinality.");
    require_throws<std::invalid_argument>(
        [&]() {
            cfd::CellScalarField non_finite_field{mesh.cell_count(), std::numeric_limits<double>::quiet_NaN()};
            limited.compute_flux_balance(non_finite_field, conditions, face_flux, gradient, limiter, limited_balance);
        },
        "Barth-Jespersen accepted a non-finite scalar value.");
    require_throws<std::invalid_argument>(
        [&]() {
            cfd::CellVectorField non_finite_gradient{mesh.cell_count(), {std::numeric_limits<double>::infinity(), 0.0}};
            limited.compute_flux_balance(field, conditions, face_flux, non_finite_gradient, limiter, limited_balance);
        },
        "Barth-Jespersen accepted a non-finite gradient.");
    require_throws<std::invalid_argument>(
        [&]() {
            cfd::FaceFluxField non_finite_flux{mesh.face_count()};
            non_finite_flux[first_boundary_face_id(mesh)] = std::numeric_limits<double>::infinity();
            limited.compute_flux_balance(field, conditions, non_finite_flux, gradient, limiter, limited_balance);
        },
        "Barth-Jespersen accepted a non-finite face flux.");
    require_throws<std::invalid_argument>(
        [&]() {
            limited.compute_flux_balance(field, conditions, face_flux, gradient, field.values(), limited_balance);
        },
        "Barth-Jespersen accepted a workspace aliasing its scalar input.");

    const auto offset_flux_workspace{face_flux.values().subspan(1, mesh.cell_count())};
    require(offset_flux_workspace.data() != face_flux.values().data(),
            "The partial-overlap fixture did not use a nonzero face-flux offset.");
    require_throws<std::invalid_argument>(
        [&]() {
            limited.compute_flux_balance(field, conditions, face_flux, gradient, offset_flux_workspace,
                                         limited_balance);
        },
        "Barth-Jespersen accepted a workspace partially overlapping face-flux storage.");
    require_throws<std::invalid_argument>(
        [&]() {
            limited.add_deferred_correction_rhs(field, conditions, gradient, face_flux, offset_flux_workspace, rhs);
        },
        "Barth-Jespersen assembly accepted a workspace partially overlapping face-flux storage.");
    const auto offset_flux_rhs{face_flux.values().subspan(1, mesh.cell_count())};
    require_throws<std::invalid_argument>(
        [&]() { default_unlimited.add_deferred_correction_rhs(gradient, face_flux, offset_flux_rhs); },
        "Unlimited LinearUpwind accepted an RHS partially overlapping face-flux storage.");
    require_throws<std::invalid_argument>(
        [&]() {
            limited.add_deferred_correction_rhs(field, conditions, gradient, face_flux, limiter, offset_flux_rhs);
        },
        "Barth-Jespersen accepted an RHS partially overlapping face-flux storage.");
    require_throws<std::invalid_argument>(
        [&]() {
            limited.compute_flux_balance(field, conditions, face_flux, gradient, limited_balance.values(),
                                         limited_balance);
        },
        "Barth-Jespersen accepted a workspace overlapping its output.");
    require_throws<std::invalid_argument>(
        [&]() { limited.add_deferred_correction_rhs(field, conditions, gradient, face_flux, rhs, rhs); },
        "Barth-Jespersen accepted a workspace aliasing its RHS.");
    require_throws<std::invalid_argument>(
        [&]() { limited.add_deferred_correction_rhs(field, conditions, gradient, face_flux, limiter, field.values()); },
        "Barth-Jespersen accepted an RHS overlapping its scalar input.");

    const cfd::CellScalarField field_before{field};
    const cfd::CellVectorField gradient_before{gradient};
    const cfd::FaceFluxField face_flux_before{face_flux};
    limited.compute_flux_balance(field, conditions, face_flux, gradient, limiter, limited_balance);
    limited.add_deferred_correction_rhs(field, conditions, gradient, face_flux, limiter, rhs);
    require_near(field[0], field_before[0], 0.0, "Barth-Jespersen modified its scalar input.");
    require_near(gradient[0].x, gradient_before[0].x, 0.0, "Barth-Jespersen modified its gradient input.");
    require_near(gradient[0].y, gradient_before[0].y, 0.0, "Barth-Jespersen modified its gradient input.");
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        require_near(face_flux[face_id], face_flux_before[face_id], 0.0,
                     "Barth-Jespersen modified its face-flux input.");
    }
}

} // namespace

int main()
{
    int failure_count{};

    failure_count += cfd::test::run_test("scalar convection positive internal flux", test_positive_internal_flux);
    failure_count += cfd::test::run_test("scalar convection negative internal flux", test_negative_internal_flux);
    failure_count +=
        cfd::test::run_test("scalar convection internal-face conservation", test_internal_face_conservation);
    failure_count += cfd::test::run_test("scalar convection Dirichlet inflow", test_dirichlet_inflow);
    failure_count += cfd::test::run_test("scalar convection Dirichlet outflow", test_dirichlet_outflow);
    failure_count +=
        cfd::test::run_test("scalar convection zeroGradient inflow and outflow", test_zero_gradient_inflow_and_outflow);
    failure_count += cfd::test::run_test("scalar convection non-zero Neumann inflow", test_nonzero_neumann_inflow);
    failure_count += cfd::test::run_test("scalar convection non-zero Neumann outflow", test_nonzero_neumann_outflow);
    failure_count += cfd::test::run_test("scalar convection zero face flux", test_zero_face_flux);
    failure_count +=
        cfd::test::run_test("scalar convection constant-field preservation", test_constant_field_preservation);
    failure_count += cfd::test::run_test("scalar convection assembly identity", test_assembly_matches_direct_balance);
    failure_count += cfd::test::run_test("scalar convection additive assembly", test_assembly_is_additive);
    failure_count += cfd::test::run_test("linear convection geometry-aware interpolation",
                                         test_linear_internal_interpolation_is_geometry_aware);
    failure_count += cfd::test::run_test("linear convection Dirichlet boundary",
                                         test_linear_dirichlet_boundary_for_both_flux_directions);
    failure_count += cfd::test::run_test("linear convection zeroGradient boundary",
                                         test_linear_zero_gradient_boundary_for_both_flux_directions);
    failure_count += cfd::test::run_test("linear convection non-zero Neumann boundary",
                                         test_linear_nonzero_neumann_boundary_for_both_flux_directions);
    failure_count +=
        cfd::test::run_test("linear convection constant-field preservation", test_linear_constant_field_preservation);
    failure_count +=
        cfd::test::run_test("linear convection assembly identity", test_linear_assembly_matches_direct_balance);
    failure_count += cfd::test::run_test("linear convection additive assembly", test_linear_assembly_is_additive);
    failure_count += cfd::test::run_test("linear convection validation and input immutability",
                                         test_linear_validation_and_input_immutability);

    failure_count += cfd::test::run_test("LinearUpwind exact linear reconstruction",
                                         test_linear_upwind_exact_linear_internal_reconstruction);
    failure_count += cfd::test::run_test("LinearUpwind upwind-gradient selection and conservation",
                                         test_linear_upwind_selects_upwind_gradient_and_adds_conservatively);
    failure_count += cfd::test::run_test("LinearUpwind matrix equals FirstOrderUpwind",
                                         test_linear_upwind_matrix_and_base_boundary_rhs_match_upwind);
    failure_count +=
        cfd::test::run_test("LinearUpwind assembly identity", test_linear_upwind_assembly_matches_direct_balance);
    failure_count += cfd::test::run_test("LinearUpwind boundary treatments", test_linear_upwind_boundary_treatments);
    failure_count += cfd::test::run_test("LinearUpwind validation, no-op, and input immutability",
                                         test_linear_upwind_validation_noop_and_input_immutability);
    failure_count += cfd::test::run_test("Barth-Jespersen constant, monotone, and bounded reconstruction",
                                         test_barth_jespersen_constant_monotone_and_bounded_reconstruction);
    failure_count += cfd::test::run_test("Barth-Jespersen internal selection, conservation, and assembly",
                                         test_barth_jespersen_internal_flux_selection_conservation_and_assembly);
    failure_count += cfd::test::run_test("Barth-Jespersen boundary bounds and flow directions",
                                         test_barth_jespersen_boundary_bounds_and_flow_directions);
    failure_count += cfd::test::run_test("Barth-Jespersen validation, None compatibility, and input immutability",
                                         test_barth_jespersen_validation_none_compatibility_and_input_immutability);
    failure_count +=
        cfd::test::run_test("hybrid convection uniform internal switch", test_hybrid_uniform_internal_switch);
    failure_count +=
        cfd::test::run_test("hybrid convection nonuniform internal switch", test_hybrid_nonuniform_internal_switch);
    failure_count += cfd::test::run_test("hybrid convection boundary switch", test_hybrid_boundary_switch);
    failure_count +=
        cfd::test::run_test("hybrid convection assembly identity", test_hybrid_assembly_matches_direct_balance);
    failure_count += cfd::test::run_test("hybrid convection face classification diagnostic",
                                         test_hybrid_face_classification_diagnostic);
    failure_count +=
        cfd::test::run_test("hybrid convection conductance validation", test_hybrid_conductance_validation);
    failure_count +=
        cfd::test::run_test("scalar convection unsupported scheme rejection", test_rejects_unsupported_scheme);
    failure_count += cfd::test::run_test("scalar convection input validation and alias rejection",
                                         test_rejects_incompatible_inputs_and_aliasing);
    failure_count += cfd::test::run_test("scalar convection input immutability", test_does_not_mutate_inputs);

    return cfd::test::finish_tests(failure_count, "scalar convection operator");
}
