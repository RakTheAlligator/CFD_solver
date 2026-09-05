#include "cfd/numerics/ScalarConvectionOperator.hpp"

#include "cfd/field/CellScalarField.hpp"
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
    failure_count +=
        cfd::test::run_test("scalar convection unsupported scheme rejection", test_rejects_unsupported_scheme);
    failure_count += cfd::test::run_test("scalar convection input validation and alias rejection",
                                         test_rejects_incompatible_inputs_and_aliasing);
    failure_count += cfd::test::run_test("scalar convection input immutability", test_does_not_mutate_inputs);

    return cfd::test::finish_tests(failure_count, "scalar convection operator");
}
