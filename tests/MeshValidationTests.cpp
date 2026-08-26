#include "cfd/mesh/Boundary.hpp"
#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/meshing/RawMeshData.hpp"
#include "cfd/meshing/RawMeshValidation.hpp"

#include "support/MeshFixtures.hpp"
#include "support/TestUtils.hpp"

#include <bit>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace
{

using cfd::test::make_non_manifold_raw_mesh;
using cfd::test::make_single_quadrilateral_raw_mesh;
using cfd::test::make_single_triangle_raw_mesh;
using cfd::test::make_two_triangle_raw_mesh;
using cfd::test::require_throws_with_message;

void require_raw_mesh_rejected(const cfd::RawMeshData &raw_mesh, const std::string_view expected_message,
                               const std::string &failure_message)
{
    require_throws_with_message<std::runtime_error>([&raw_mesh]() { cfd::validate_raw_mesh(raw_mesh); },
                                                    expected_message, failure_message);
}

void test_rejects_empty_nodes()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};
    raw_mesh.nodes.clear();

    require_raw_mesh_rejected(raw_mesh, "nodes array is empty.", "Raw mesh validation accepted an empty node array.");
}

void test_rejects_non_finite_node_x()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};
    raw_mesh.nodes[0].x = std::numeric_limits<double>::quiet_NaN();

    require_raw_mesh_rejected(raw_mesh, "node 0 contains a non-finite coordinate.",
                              "Raw mesh validation accepted a NaN node coordinate.");
}

void test_rejects_non_finite_node_y()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};
    raw_mesh.nodes[0].y = std::numeric_limits<double>::infinity();

    require_raw_mesh_rejected(raw_mesh, "node 0 contains a non-finite coordinate.",
                              "Raw mesh validation accepted an infinite node coordinate.");
}

void test_rejects_empty_cell_types()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};
    raw_mesh.cell_types.clear();

    require_raw_mesh_rejected(raw_mesh, "cell_types array is empty.",
                              "Raw mesh validation accepted an empty cell-type array.");
}

void test_rejects_empty_cell_nodes()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};
    raw_mesh.cell_nodes.clear();

    require_raw_mesh_rejected(raw_mesh, "cell_nodes array is empty.",
                              "Raw mesh validation accepted an empty cell-connectivity array.");
}

void test_rejects_empty_cell_node_offsets()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};
    raw_mesh.cell_node_offsets.clear();

    require_raw_mesh_rejected(raw_mesh, "cell_node_offsets array is empty.",
                              "Raw mesh validation accepted an empty cell-offset array.");
}

void test_rejects_incorrect_cell_offset_count()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};
    raw_mesh.cell_node_offsets = {0};

    require_raw_mesh_rejected(raw_mesh, "cell_node_offsets size must equal the number of cells + 1.",
                              "Raw mesh validation accepted an incorrect cell-offset count.");
}

void test_rejects_non_zero_first_cell_offset()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};
    raw_mesh.cell_node_offsets.front() = 1;

    require_raw_mesh_rejected(raw_mesh, "cell_node_offsets must start at 0.",
                              "Raw mesh validation accepted a non-zero first cell offset.");
}

void test_rejects_non_increasing_cell_offsets()
{
    cfd::RawMeshData raw_mesh{make_two_triangle_raw_mesh()};
    raw_mesh.cell_node_offsets[2] = raw_mesh.cell_node_offsets[1];

    require_raw_mesh_rejected(raw_mesh, "cell_node_offsets must be strictly increasing.",
                              "Raw mesh validation accepted non-increasing cell offsets.");
}

void test_rejects_cell_offset_outside_connectivity()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};
    raw_mesh.cell_node_offsets.back() = raw_mesh.cell_nodes.size() + 1;

    require_raw_mesh_rejected(raw_mesh, "cell_node_offsets contains an offset outside cell_nodes.",
                              "Raw mesh validation accepted a cell offset outside the connectivity array.");
}

void test_rejects_incorrect_last_cell_offset()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};
    raw_mesh.cell_node_offsets.back() = raw_mesh.cell_nodes.size() - 1;

    require_raw_mesh_rejected(raw_mesh, "last cell_node_offset must equal cell_nodes.size().",
                              "Raw mesh validation accepted an incorrect final cell offset.");
}

void test_rejects_wrong_triangle_node_count()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};
    raw_mesh.cell_nodes.push_back(0);
    raw_mesh.cell_node_offsets.back() = raw_mesh.cell_nodes.size();

    require_raw_mesh_rejected(raw_mesh, "cell 0 has 4 nodes, but its type requires 3.",
                              "Raw mesh validation accepted a triangle with four nodes.");
}

void test_rejects_unsupported_cell_type()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};
    raw_mesh.cell_types[0] = std::bit_cast<cfd::CellType>(std::uint8_t{255});

    require_raw_mesh_rejected(raw_mesh, "cell 0 has an unsupported CellType.",
                              "Raw mesh validation accepted an unsupported cell type.");
}

void test_rejects_duplicate_node_in_cell()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};
    raw_mesh.cell_nodes[2] = raw_mesh.cell_nodes[1];

    require_raw_mesh_rejected(raw_mesh, "cell 0 contains duplicate node indices.",
                              "Raw mesh validation accepted a duplicated node inside a cell.");
}

void test_rejects_invalid_cell_node_index()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};
    raw_mesh.cell_nodes[2] = raw_mesh.nodes.size();

    require_raw_mesh_rejected(raw_mesh, "which is outside the nodes array.",
                              "Raw mesh validation accepted an out-of-range cell node index.");
}

void test_rejects_empty_boundary_groups()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};
    raw_mesh.boundary_groups.clear();

    require_raw_mesh_rejected(raw_mesh, "boundary_groups array is empty.",
                              "Raw mesh validation accepted an empty boundary-group array.");
}

void test_rejects_reserved_boundary_group_id()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};

    raw_mesh.boundary_groups[0].id = cfd::invalid_boundary_id;
    for (cfd::BoundaryEdge &edge : raw_mesh.boundary_edges)
    {
        edge.boundary_id = cfd::invalid_boundary_id;
    }

    require_raw_mesh_rejected(raw_mesh, "uses the reserved invalid boundary ID.",
                              "Raw mesh validation accepted invalid_boundary_id as a physical boundary ID.");
}

void test_rejects_empty_boundary_group_name()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};
    raw_mesh.boundary_groups[0].name.clear();

    require_raw_mesh_rejected(raw_mesh, "a boundary group has an empty name.",
                              "Raw mesh validation accepted an empty boundary-group name.");
}

void test_rejects_duplicate_boundary_group_id()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};
    raw_mesh.boundary_groups.push_back({
        raw_mesh.boundary_groups[0].id,
        "duplicate-id",
    });

    require_raw_mesh_rejected(raw_mesh, "boundary group ID 0 is duplicated.",
                              "Raw mesh validation accepted a duplicated boundary-group ID.");
}

void test_rejects_duplicate_boundary_group_name()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};
    raw_mesh.boundary_groups.push_back({
        1,
        raw_mesh.boundary_groups[0].name,
    });

    require_raw_mesh_rejected(raw_mesh, "boundary group name \"wall\" is duplicated.",
                              "Raw mesh validation accepted a duplicated boundary-group name.");
}

void test_rejects_empty_boundary_edges()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};
    raw_mesh.boundary_edges.clear();

    require_raw_mesh_rejected(raw_mesh, "boundary_edges array is empty.",
                              "Raw mesh validation accepted an empty boundary-edge array.");
}

void test_rejects_invalid_boundary_node_index()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};
    raw_mesh.boundary_edges[0].node_ids[0] = raw_mesh.nodes.size();

    require_raw_mesh_rejected(raw_mesh, "references a node outside the nodes array.",
                              "Raw mesh validation accepted an out-of-range boundary node index.");
}

void test_rejects_zero_length_topological_boundary_edge()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};
    raw_mesh.boundary_edges[0].node_ids[1] = raw_mesh.boundary_edges[0].node_ids[0];

    require_raw_mesh_rejected(raw_mesh, "references the same node twice.",
                              "Raw mesh validation accepted a zero-length topological boundary edge.");
}

void test_rejects_unknown_boundary_group_id()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};
    raw_mesh.boundary_edges[0].boundary_id = 42;

    require_raw_mesh_rejected(raw_mesh, "references unknown boundary group ID 42.",
                              "Raw mesh validation accepted an unknown boundary-group ID.");
}

void test_rejects_duplicate_boundary_edge()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};
    raw_mesh.boundary_edges.push_back({
        {1, 0},
        raw_mesh.boundary_groups[0].id,
    });

    require_raw_mesh_rejected(raw_mesh, "the same boundary edge appears more than once.",
                              "Raw mesh validation accepted a duplicated boundary edge.");
}

void test_rejects_unused_boundary_group()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};
    raw_mesh.boundary_groups.push_back({
        1,
        "unused",
    });

    require_raw_mesh_rejected(raw_mesh, "boundary group \"unused\" contains no boundary edge.",
                              "Raw mesh validation accepted an unused boundary group.");
}

void test_rejects_open_boundary()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};
    raw_mesh.boundary_edges.pop_back();

    require_throws_with_message<std::runtime_error>(
        [&raw_mesh]() { static_cast<void>(cfd::build_mesh(std::move(raw_mesh))); },
        "Topology validation failed:", "Mesh construction accepted an incomplete physical boundary.");
}

void test_rejects_zero_area_triangle()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};

    raw_mesh.nodes[0] = {0.0, 0.0};
    raw_mesh.nodes[1] = {1.0, 0.0};
    raw_mesh.nodes[2] = {2.0, 0.0};

    cfd::validate_raw_mesh(raw_mesh);

    require_throws_with_message<std::runtime_error>(
        [&raw_mesh]() { static_cast<void>(cfd::build_mesh(std::move(raw_mesh))); },
        "Geometry construction failed:", "Mesh construction accepted a zero-area triangle.");
}

void test_rejects_non_manifold_face()
{
    cfd::RawMeshData raw_mesh{make_non_manifold_raw_mesh()};

    cfd::validate_raw_mesh(raw_mesh);

    require_throws_with_message<std::runtime_error>(
        [&raw_mesh]() { static_cast<void>(cfd::build_mesh(std::move(raw_mesh))); },
        "Topology construction failed:", "Mesh construction accepted a face shared by more than two cells.");
}

void test_rejects_concave_quadrilateral()
{
    cfd::RawMeshData raw_mesh{make_single_quadrilateral_raw_mesh()};

    raw_mesh.nodes[2] = {
        0.25,
        0.25,
    };

    cfd::validate_raw_mesh(raw_mesh);

    require_throws_with_message<std::runtime_error>(
        [&raw_mesh]() { static_cast<void>(cfd::build_mesh(std::move(raw_mesh))); },
        "non-convex or self-intersecting quadrilateral", "Mesh construction accepted a concave quadrilateral.");
}

} // namespace

int main()
{
    int failure_count{};

    failure_count += cfd::test::run_test("reject empty nodes", test_rejects_empty_nodes);
    failure_count += cfd::test::run_test("reject non-finite node x", test_rejects_non_finite_node_x);
    failure_count += cfd::test::run_test("reject non-finite node y", test_rejects_non_finite_node_y);

    failure_count += cfd::test::run_test("reject empty cell types", test_rejects_empty_cell_types);
    failure_count += cfd::test::run_test("reject empty cell nodes", test_rejects_empty_cell_nodes);
    failure_count += cfd::test::run_test("reject empty cell offsets", test_rejects_empty_cell_node_offsets);
    failure_count +=
        cfd::test::run_test("reject incorrect cell-offset count", test_rejects_incorrect_cell_offset_count);
    failure_count += cfd::test::run_test("reject non-zero first cell offset", test_rejects_non_zero_first_cell_offset);
    failure_count +=
        cfd::test::run_test("reject non-increasing cell offsets", test_rejects_non_increasing_cell_offsets);
    failure_count +=
        cfd::test::run_test("reject cell offset outside connectivity", test_rejects_cell_offset_outside_connectivity);
    failure_count += cfd::test::run_test("reject incorrect last cell offset", test_rejects_incorrect_last_cell_offset);
    failure_count += cfd::test::run_test("reject wrong triangle node count", test_rejects_wrong_triangle_node_count);
    failure_count += cfd::test::run_test("reject unsupported cell type", test_rejects_unsupported_cell_type);
    failure_count += cfd::test::run_test("reject duplicate node in cell", test_rejects_duplicate_node_in_cell);
    failure_count += cfd::test::run_test("reject invalid cell node index", test_rejects_invalid_cell_node_index);

    failure_count += cfd::test::run_test("reject empty boundary groups", test_rejects_empty_boundary_groups);
    failure_count += cfd::test::run_test("reject reserved boundary group ID", test_rejects_reserved_boundary_group_id);
    failure_count += cfd::test::run_test("reject empty boundary group name", test_rejects_empty_boundary_group_name);
    failure_count +=
        cfd::test::run_test("reject duplicate boundary group ID", test_rejects_duplicate_boundary_group_id);
    failure_count +=
        cfd::test::run_test("reject duplicate boundary group name", test_rejects_duplicate_boundary_group_name);

    failure_count += cfd::test::run_test("reject empty boundary edges", test_rejects_empty_boundary_edges);
    failure_count +=
        cfd::test::run_test("reject invalid boundary node index", test_rejects_invalid_boundary_node_index);
    failure_count += cfd::test::run_test("reject zero-length topological boundary edge",
                                         test_rejects_zero_length_topological_boundary_edge);
    failure_count += cfd::test::run_test("reject unknown boundary group ID", test_rejects_unknown_boundary_group_id);
    failure_count += cfd::test::run_test("reject duplicate boundary edge", test_rejects_duplicate_boundary_edge);
    failure_count += cfd::test::run_test("reject unused boundary group", test_rejects_unused_boundary_group);

    failure_count += cfd::test::run_test("reject open boundary", test_rejects_open_boundary);
    failure_count += cfd::test::run_test("reject zero-area triangle", test_rejects_zero_area_triangle);
    failure_count += cfd::test::run_test("reject non-manifold face", test_rejects_non_manifold_face);
    failure_count += cfd::test::run_test("reject concave quadrilateral", test_rejects_concave_quadrilateral);

    return cfd::test::finish_tests(failure_count, "mesh validation");
}
