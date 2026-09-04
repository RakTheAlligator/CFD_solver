#include "cfd/linear_algebra/ScalarLinearSystem.hpp"

#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/meshing/RawMeshData.hpp"

#include "support/TestUtils.hpp"

#include <algorithm>
#include <array>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace
{

using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws;
using cfd::test::test_tolerance;

static_assert(!std::is_copy_constructible_v<cfd::ScalarLinearSystem>);
static_assert(!std::is_copy_assignable_v<cfd::ScalarLinearSystem>);
static_assert(std::is_nothrow_move_constructible_v<cfd::ScalarLinearSystem>);
static_assert(!std::is_move_assignable_v<cfd::ScalarLinearSystem>);

[[nodiscard]]
cfd::RawMeshData make_two_cell_mesh()
{
    constexpr cfd::BoundaryId wall_boundary_id{0};

    cfd::RawMeshData raw_mesh;
    raw_mesh.nodes = {
        {0.0, 0.0}, {1.0, 0.0}, {2.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}, {2.0, 1.0},
    };
    raw_mesh.cell_types = {cfd::CellType::Quadrilateral, cfd::CellType::Quadrilateral};
    raw_mesh.cell_nodes = {0, 1, 4, 3, 1, 2, 5, 4};
    raw_mesh.cell_node_offsets = {0, 4, 8};
    raw_mesh.boundary_groups = {{wall_boundary_id, "wall"}};
    raw_mesh.boundary_edges = {
        {{0, 1}, wall_boundary_id}, {{1, 2}, wall_boundary_id}, {{2, 5}, wall_boundary_id},
        {{5, 4}, wall_boundary_id}, {{4, 3}, wall_boundary_id}, {{3, 0}, wall_boundary_id},
    };
    return raw_mesh;
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
    throw std::runtime_error("Two-cell scalar-system fixture has no internal face.");
}

void test_construction_cardinalities_and_initialization()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarLinearSystem system{mesh};

    require(system.cell_count() == mesh.cell_count(), "Scalar system has an incorrect cell count.");
    require(system.face_count() == mesh.face_count(), "Scalar system has an incorrect face count.");
    require(&system.mesh() == &mesh, "Scalar system does not reference its construction Mesh.");
    require(std::ranges::all_of(system.diagonal(), [](const double value) { return value == 0.0; }),
            "Scalar system diagonal is not zero-initialized.");
    require(std::ranges::all_of(system.owner_neighbor_coefficients(), [](const double value) { return value == 0.0; }),
            "Scalar system owner-neighbor coefficients are not zero-initialized.");
    require(std::ranges::all_of(system.neighbor_owner_coefficients(), [](const double value) { return value == 0.0; }),
            "Scalar system neighbor-owner coefficients are not zero-initialized.");
    require(std::ranges::all_of(system.rhs(), [](const double value) { return value == 0.0; }),
            "Scalar system right-hand side is not zero-initialized.");
}

void test_clear_operations_are_separated()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    cfd::ScalarLinearSystem system{mesh};

    std::ranges::fill(system.diagonal(), 1.0);
    std::ranges::fill(system.owner_neighbor_coefficients(), 2.0);
    std::ranges::fill(system.neighbor_owner_coefficients(), 3.0);
    std::ranges::fill(system.rhs(), 4.0);
    system.clear_matrix();

    require(std::ranges::all_of(system.diagonal(), [](const double value) { return value == 0.0; }),
            "clear_matrix() did not clear the diagonal.");
    require(std::ranges::all_of(system.owner_neighbor_coefficients(), [](const double value) { return value == 0.0; }),
            "clear_matrix() did not clear owner-neighbor coefficients.");
    require(std::ranges::all_of(system.neighbor_owner_coefficients(), [](const double value) { return value == 0.0; }),
            "clear_matrix() did not clear neighbor-owner coefficients.");
    require(std::ranges::all_of(system.rhs(), [](const double value) { return value == 4.0; }),
            "clear_matrix() modified the right-hand side.");

    std::ranges::fill(system.diagonal(), 5.0);
    system.clear_rhs();
    require(std::ranges::all_of(system.diagonal(), [](const double value) { return value == 5.0; }),
            "clear_rhs() modified matrix coefficients.");
    require(std::ranges::all_of(system.rhs(), [](const double value) { return value == 0.0; }),
            "clear_rhs() did not clear the right-hand side.");

    std::ranges::fill(system.owner_neighbor_coefficients(), 6.0);
    std::ranges::fill(system.neighbor_owner_coefficients(), 7.0);
    std::ranges::fill(system.rhs(), 8.0);
    system.clear();
    require(std::ranges::all_of(system.diagonal(), [](const double value) { return value == 0.0; }) &&
                std::ranges::all_of(system.owner_neighbor_coefficients(),
                                    [](const double value) { return value == 0.0; }) &&
                std::ranges::all_of(system.neighbor_owner_coefficients(),
                                    [](const double value) { return value == 0.0; }) &&
                std::ranges::all_of(system.rhs(), [](const double value) { return value == 0.0; }),
            "clear() did not clear the complete scalar system.");
}

void test_apply_matrix_supports_non_symmetric_face_coefficients()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    cfd::ScalarLinearSystem system{mesh};
    const cfd::Index face_id{internal_face_id(mesh)};

    system.diagonal()[0] = 2.0;
    system.diagonal()[1] = 3.0;
    system.owner_neighbor_coefficients()[face_id] = 4.0;
    system.neighbor_owner_coefficients()[face_id] = -5.0;
    for (cfd::Index boundary_face_id = 0; boundary_face_id < mesh.face_count(); ++boundary_face_id)
    {
        if (mesh.face_adjacencies()[boundary_face_id].is_boundary())
        {
            system.owner_neighbor_coefficients()[boundary_face_id] = 1000.0;
            system.neighbor_owner_coefficients()[boundary_face_id] = -1000.0;
        }
    }

    const std::array input{7.0, 11.0};
    std::array output{99.0, 99.0};
    system.apply_matrix(input, output);

    require_near(output[0], 58.0, test_tolerance, "A(owner,neighbor) was applied incorrectly.");
    require_near(output[1], -2.0, test_tolerance, "A(neighbor,owner) was applied incorrectly.");
}

void test_apply_matrix_validates_cardinality_and_aliasing()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_mesh())};
    const cfd::ScalarLinearSystem system{build_result.mesh};
    const std::array valid_input{1.0, 2.0};
    std::array valid_output{0.0, 0.0};

    require_throws<std::invalid_argument>(
        [&system, &valid_output]() {
            const std::array wrong_input{1.0, 2.0, 3.0};
            system.apply_matrix(wrong_input, valid_output);
        },
        "Scalar system accepted an input with incorrect cardinality.");
    require_throws<std::invalid_argument>(
        [&system, &valid_input]() {
            std::array wrong_output{0.0, 0.0, 0.0};
            system.apply_matrix(valid_input, wrong_output);
        },
        "Scalar system accepted an output with incorrect cardinality.");
    require_throws<std::invalid_argument>(
        [&system]() {
            std::array aliased{1.0, 2.0};
            system.apply_matrix(aliased, aliased);
        },
        "Scalar system accepted overlapping input and output spans.");
}

} // namespace

int main()
{
    int failure_count{};

    failure_count +=
        cfd::test::run_test("scalar linear system construction", test_construction_cardinalities_and_initialization);
    failure_count += cfd::test::run_test("scalar linear system clearing", test_clear_operations_are_separated);
    failure_count += cfd::test::run_test("scalar linear system non-symmetric matrix application",
                                         test_apply_matrix_supports_non_symmetric_face_coefficients);
    failure_count += cfd::test::run_test("scalar linear system cardinality and alias validation",
                                         test_apply_matrix_validates_cardinality_and_aliasing);

    return cfd::test::finish_tests(failure_count, "scalar linear system");
}
