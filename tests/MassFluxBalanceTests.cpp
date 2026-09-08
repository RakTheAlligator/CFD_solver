#include "cfd/numerics/MassFluxBalance.hpp"

#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/FaceFluxField.hpp"
#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/mesh/Types.hpp"
#include "cfd/meshing/RawMeshData.hpp"

#include "support/TestUtils.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace
{

using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws;

constexpr cfd::BoundaryId bottom_boundary_id{0};
constexpr cfd::BoundaryId right_boundary_id{1};
constexpr cfd::BoundaryId top_boundary_id{2};
constexpr cfd::BoundaryId left_boundary_id{3};

[[nodiscard]]
cfd::RawMeshData make_two_cell_rectangle_raw_mesh()
{
    cfd::RawMeshData raw_mesh;
    raw_mesh.nodes = {
        {0.0, 0.0}, {1.0, 0.0}, {2.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}, {2.0, 1.0},
    };
    raw_mesh.cell_types = {cfd::CellType::Quadrilateral, cfd::CellType::Quadrilateral};
    raw_mesh.cell_nodes = {0, 1, 4, 3, 1, 2, 5, 4};
    raw_mesh.cell_node_offsets = {0, 4, 8};
    raw_mesh.boundary_groups = {
        {bottom_boundary_id, "bottom"},
        {right_boundary_id, "right"},
        {top_boundary_id, "top"},
        {left_boundary_id, "left"},
    };
    raw_mesh.boundary_edges = {
        {{0, 1}, bottom_boundary_id}, {{1, 2}, bottom_boundary_id}, {{2, 5}, right_boundary_id},
        {{5, 4}, top_boundary_id},    {{4, 3}, top_boundary_id},    {{3, 0}, left_boundary_id},
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
    throw std::runtime_error("Mass-flux balance fixture has no internal face.");
}

[[nodiscard]]
cfd::Index boundary_face_id(const cfd::Mesh &mesh, const cfd::BoundaryId boundary_id)
{
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (mesh.face_boundary_ids()[face_id] == boundary_id)
        {
            return face_id;
        }
    }
    throw std::runtime_error("Mass-flux balance fixture has no requested boundary face.");
}

void test_owner_neighbor_and_boundary_signs()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::Index internal_id{internal_face_id(mesh)};
    const cfd::FaceAdjacency &internal_adjacency{mesh.face_adjacencies()[internal_id]};
    const cfd::Index left_face_id{boundary_face_id(mesh, left_boundary_id)};
    const cfd::Index right_face_id{boundary_face_id(mesh, right_boundary_id)};
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    mass_flux[internal_id] = 2.0;
    mass_flux[left_face_id] = 3.0;
    mass_flux[right_face_id] = -4.0;
    cfd::CellScalarField imbalance{mesh.cell_count(), 99.0};

    cfd::compute_cell_mass_imbalance(mesh, mass_flux, imbalance);

    require_near(imbalance[internal_adjacency.owner], 5.0, 0.0, "Owner or owner-boundary mass-flux sign is incorrect.");
    require_near(imbalance[internal_adjacency.neighbor], -6.0, 0.0,
                 "Neighbor or neighbor-boundary mass-flux sign is incorrect.");
}

void test_multiple_faces_accumulate_and_internal_flux_is_globally_conservative()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::Index internal_id{internal_face_id(mesh)};
    const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[internal_id]};
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    mass_flux[internal_id] = 7.0;
    cfd::CellScalarField imbalance{mesh.cell_count(), 99.0};
    cfd::compute_cell_mass_imbalance(mesh, mass_flux, imbalance);
    require_near(imbalance[adjacency.owner], 7.0, 0.0, "Internal owner accumulation is incorrect.");
    require_near(imbalance[adjacency.neighbor], -7.0, 0.0, "Internal neighbor accumulation is incorrect.");
    require_near(imbalance[0] + imbalance[1], 0.0, 0.0, "Internal face contribution is not globally conservative.");

    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        mass_flux[face_id] = 0.0;
    }
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        const cfd::FaceAdjacency &face_adjacency{mesh.face_adjacencies()[face_id]};
        if (face_adjacency.is_boundary() && face_adjacency.owner == adjacency.owner)
        {
            mass_flux[face_id] = 1.0;
        }
    }
    cfd::compute_cell_mass_imbalance(mesh, mass_flux, imbalance);
    require_near(imbalance[adjacency.owner], 3.0, 0.0,
                 "Multiple boundary faces did not accumulate on their owner cell.");
    require_near(imbalance[adjacency.neighbor], 0.0, 0.0, "Boundary faces contributed to a non-owner cell.");
}

void test_cardinality_validation_preserves_output()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    cfd::FaceFluxField wrong_mass_flux{mesh.face_count() - 1};
    cfd::CellScalarField imbalance{mesh.cell_count(), 11.0};
    cfd::CellScalarField wrong_imbalance{mesh.cell_count() - 1, 12.0};

    require_throws<std::invalid_argument>([&]() { cfd::compute_cell_mass_imbalance(mesh, wrong_mass_flux, imbalance); },
                                          "Mass-flux balance accepted the wrong face cardinality.");
    require(imbalance[0] == 11.0 && imbalance[1] == 11.0, "Wrong face cardinality modified the mass-imbalance output.");
    require_throws<std::invalid_argument>([&]() { cfd::compute_cell_mass_imbalance(mesh, mass_flux, wrong_imbalance); },
                                          "Mass-flux balance accepted the wrong cell cardinality.");
    require(wrong_imbalance[0] == 12.0, "Wrong cell cardinality modified the mass-imbalance output.");
}

void test_nonfinite_flux_rejection_is_transactional()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::Index later_face_id{mesh.face_count() - 1};
    for (const double invalid_flux : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
                                      -std::numeric_limits<double>::infinity()})
    {
        cfd::FaceFluxField mass_flux{mesh.face_count(), 1.0};
        mass_flux[later_face_id] = invalid_flux;
        cfd::CellScalarField imbalance{mesh.cell_count()};
        imbalance[0] = 31.0;
        imbalance[1] = 32.0;

        require_throws<std::runtime_error>([&]() { cfd::compute_cell_mass_imbalance(mesh, mass_flux, imbalance); },
                                           "Mass-flux balance accepted a non-finite face flux.");
        require(imbalance[0] == 31.0 && imbalance[1] == 32.0,
                "Later non-finite face flux partially modified the output.");
    }
}

} // namespace

int main()
{
    int failure_count{};

    failure_count += cfd::test::run_test("mass-flux balance signs", test_owner_neighbor_and_boundary_signs);
    failure_count += cfd::test::run_test("mass-flux balance accumulation and conservation",
                                         test_multiple_faces_accumulate_and_internal_flux_is_globally_conservative);
    failure_count +=
        cfd::test::run_test("mass-flux balance cardinality validation", test_cardinality_validation_preserves_output);
    failure_count += cfd::test::run_test("mass-flux balance transactional validation",
                                         test_nonfinite_flux_rejection_is_transactional);

    return cfd::test::finish_tests(failure_count, "mass-flux balance");
}
