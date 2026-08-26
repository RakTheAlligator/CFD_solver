#include "cfd/mesh/Boundary.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/meshing/RawMeshData.hpp"
#include "cfd/meshing/RawMeshValidation.hpp"

#include "support/MeshFixtures.hpp"
#include "support/TestUtils.hpp"

#include <stdexcept>
#include <utility>

namespace
{

using cfd::test::make_non_manifold_raw_mesh;
using cfd::test::make_single_triangle_raw_mesh;
using cfd::test::require_throws_with_message;

void test_rejects_inconsistent_cell_offsets()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};

    raw_mesh.cell_node_offsets.back() = 2;

    require_throws_with_message<std::runtime_error>(
        [&raw_mesh]() { cfd::validate_raw_mesh(raw_mesh); },
        "Raw mesh validation failed:", "Raw mesh validation accepted inconsistent cell offsets.");
}

void test_rejects_duplicate_node_in_cell()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};

    raw_mesh.cell_nodes[2] = raw_mesh.cell_nodes[1];

    require_throws_with_message<std::runtime_error>(
        [&raw_mesh]() { cfd::validate_raw_mesh(raw_mesh); },
        "Raw mesh validation failed:", "Raw mesh validation accepted a duplicated node inside a cell.");
}

void test_rejects_invalid_node_index()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};

    raw_mesh.cell_nodes[2] = raw_mesh.nodes.size();

    require_throws_with_message<std::runtime_error>(
        [&raw_mesh]() { cfd::validate_raw_mesh(raw_mesh); },
        "Raw mesh validation failed:", "Raw mesh validation accepted an out-of-range node index.");
}

void test_rejects_duplicate_boundary_edge()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};

    raw_mesh.boundary_edges.push_back({
        {1, 0},
        raw_mesh.boundary_groups[0].id,
    });

    require_throws_with_message<std::runtime_error>(
        [&raw_mesh]() { cfd::validate_raw_mesh(raw_mesh); },
        "Raw mesh validation failed:", "Raw mesh validation accepted a duplicated boundary edge.");
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

void test_rejects_invalid_boundary_group_id()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};

    raw_mesh.boundary_groups[0].id = cfd::invalid_boundary_id;

    for (cfd::BoundaryEdge &edge : raw_mesh.boundary_edges)
    {
        edge.boundary_id = cfd::invalid_boundary_id;
    }

    require_throws_with_message<std::runtime_error>(
        [&raw_mesh]() { cfd::validate_raw_mesh(raw_mesh); },
        "Raw mesh validation failed:", "Raw mesh validation accepted invalid_boundary_id as a physical boundary ID.");
}

} // namespace

int main()
{
    int failure_count{};

    failure_count += cfd::test::run_test("reject inconsistent cell offsets", test_rejects_inconsistent_cell_offsets);

    failure_count += cfd::test::run_test("reject duplicate node in cell", test_rejects_duplicate_node_in_cell);

    failure_count += cfd::test::run_test("reject invalid node index", test_rejects_invalid_node_index);

    failure_count += cfd::test::run_test("reject duplicate boundary edge", test_rejects_duplicate_boundary_edge);

    failure_count += cfd::test::run_test("reject open boundary", test_rejects_open_boundary);

    failure_count += cfd::test::run_test("reject zero-area triangle", test_rejects_zero_area_triangle);

    failure_count += cfd::test::run_test("reject non-manifold face", test_rejects_non_manifold_face);

    failure_count += cfd::test::run_test("reject invalid boundary group ID", test_rejects_invalid_boundary_group_id);

    return cfd::test::finish_tests(failure_count, "mesh validation");
}