#include "cfd/mesh/Boundary.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/mesh/Types.hpp"
#include "cfd/meshing/RawMeshData.hpp"

#include "support/MeshFixtures.hpp"
#include "support/TestUtils.hpp"

#include <utility>

namespace
{

using cfd::test::make_two_triangle_raw_mesh;
using cfd::test::require;

[[nodiscard]]
bool face_matches_nodes(const cfd::Face &face, const cfd::Index node_0_id, const cfd::Index node_1_id)
{
    return (face.node_ids[0] == node_0_id && face.node_ids[1] == node_1_id) ||
           (face.node_ids[0] == node_1_id && face.node_ids[1] == node_0_id);
}

[[nodiscard]]
cfd::Index find_face_id(const cfd::Mesh &mesh, const cfd::Index node_0_id, const cfd::Index node_1_id)
{
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (face_matches_nodes(mesh.faces()[face_id], node_0_id, node_1_id))
        {
            return face_id;
        }
    }

    return cfd::invalid_index;
}

[[nodiscard]]
bool cell_contains_face(const cfd::Mesh &mesh, const cfd::Index cell_id, const cfd::Index face_id)
{
    const cfd::Index cell_face_begin_offset{mesh.cell_node_offsets()[cell_id]};
    const cfd::Index cell_face_end_offset{mesh.cell_node_offsets()[cell_id + 1]};

    for (cfd::Index cell_face_position = cell_face_begin_offset; cell_face_position < cell_face_end_offset;
         ++cell_face_position)
    {
        if (mesh.cell_faces()[cell_face_position] == face_id)
        {
            return true;
        }
    }

    return false;
}

void require_boundary_face(const cfd::Mesh &mesh, const cfd::Index node_0_id, const cfd::Index node_1_id,
                           const cfd::Index expected_owner_id, const cfd::BoundaryId expected_boundary_id)
{
    const cfd::Index face_id{find_face_id(mesh, node_0_id, node_1_id)};

    require(face_id != cfd::invalid_index, "Expected boundary face was not constructed.");

    const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};

    require(adjacency.owner == expected_owner_id, "Boundary face has an incorrect owner.");
    require(adjacency.is_boundary(), "Expected boundary face was classified as internal.");
    require(mesh.face_boundary_ids()[face_id] == expected_boundary_id, "Boundary face has an incorrect boundary ID.");
}

void test_two_triangle_topology()
{
    constexpr cfd::BoundaryId boundary_id{0};

    cfd::RawMeshData raw_mesh{make_two_triangle_raw_mesh()};
    cfd::MeshBuildResult build_result{cfd::build_mesh(std::move(raw_mesh))};

    const cfd::Mesh &mesh{build_result.mesh};

    require(mesh.node_count() == 4, "Two-triangle mesh must contain four nodes.");
    require(mesh.cell_count() == 2, "Two-triangle mesh must contain two cells.");
    require(mesh.face_count() == 5, "Two triangles sharing one edge must produce five unique faces.");
    require(mesh.cell_faces().size() == 6, "Two triangles must contain six cell-face incidences.");

    cfd::Index internal_face_count{};
    cfd::Index boundary_face_count{};

    for (const cfd::FaceAdjacency &adjacency : mesh.face_adjacencies())
    {
        if (adjacency.is_boundary())
        {
            ++boundary_face_count;
        }
        else
        {
            ++internal_face_count;
        }
    }

    require(internal_face_count == 1, "Two-triangle mesh must contain one internal face.");
    require(boundary_face_count == 4, "Two-triangle mesh must contain four boundary faces.");

    const cfd::Index shared_face_id{find_face_id(mesh, 0, 2)};

    require(shared_face_id != cfd::invalid_index, "Shared triangle face was not constructed.");

    const cfd::FaceAdjacency &shared_adjacency{mesh.face_adjacencies()[shared_face_id]};

    // Cell 0 encounters the shared edge first in the fixture, so topology
    // construction assigns cell 0 as owner and cell 1 as neighbor.
    require(shared_adjacency.owner == 0, "Shared face must have cell 0 as owner.");
    require(shared_adjacency.neighbor == 1, "Shared face must have cell 1 as neighbor.");

    require(!shared_adjacency.is_boundary(), "Shared face must be internal.");
    require(mesh.face_boundary_ids()[shared_face_id] == cfd::invalid_boundary_id,
            "Internal face must not have a boundary ID.");

    require(cell_contains_face(mesh, 0, shared_face_id), "Cell 0 does not reference the shared face.");
    require(cell_contains_face(mesh, 1, shared_face_id), "Cell 1 does not reference the shared face.");

    require_boundary_face(mesh, 0, 1, 0, boundary_id);
    require_boundary_face(mesh, 1, 2, 0, boundary_id);
    require_boundary_face(mesh, 2, 3, 1, boundary_id);
    require_boundary_face(mesh, 3, 0, 1, boundary_id);
}

} // namespace

int main()
{
    int failure_count{};

    failure_count += cfd::test::run_test("two triangle topology", test_two_triangle_topology);

    return cfd::test::finish_tests(failure_count, "mesh topology");
}