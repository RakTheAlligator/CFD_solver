#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/mesh/MeshStatistics.hpp"
#include "cfd/meshing/GmshMesher.hpp"
#include "cfd/meshing/RawMeshData.hpp"

#include "support/MeshFixtures.hpp"
#include "support/TestUtils.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <type_traits>
#include <utility>
#include <vector>

namespace
{

using cfd::test::fail;
using cfd::test::make_single_triangle_raw_mesh;
using cfd::test::require;
using cfd::test::require_near;
using cfd::test::test_tolerance;

static_assert(!std::is_copy_constructible_v<cfd::Mesh>);
static_assert(!std::is_copy_assignable_v<cfd::Mesh>);

static_assert(std::is_nothrow_move_constructible_v<cfd::Mesh>);
static_assert(std::is_nothrow_move_assignable_v<cfd::Mesh>);

double compute_single_closed_boundary_area(const cfd::RawMeshData &raw_mesh)
{
    using cfd::Index;
    using cfd::invalid_index;

    if (raw_mesh.boundary_edges.empty())
    {
        fail("Cannot compute boundary area: boundary has no edges.");
    }

    const std::array<Index, 2> empty_neighbors{invalid_index, invalid_index};

    std::vector<std::array<Index, 2>> neighbors(raw_mesh.nodes.size(), empty_neighbors);
    std::vector<std::uint8_t> degrees(raw_mesh.nodes.size(), 0);

    const auto add_neighbor = [&](const Index node_id, const Index neighbor_id) {
        switch (degrees[node_id])
        {
        case 0:
            neighbors[node_id][0] = neighbor_id;
            ++degrees[node_id];
            return;

        case 1:
            if (neighbors[node_id][0] == neighbor_id)
            {
                fail("Boundary contains a duplicated edge.");
            }

            neighbors[node_id][1] = neighbor_id;
            ++degrees[node_id];
            return;

        default:
            fail("Boundary is not a simple closed contour: a boundary node has more than two neighbors.");
        }
    };

    for (const cfd::BoundaryEdge &edge : raw_mesh.boundary_edges)
    {
        const Index node_0{edge.node_ids[0]};
        const Index node_1{edge.node_ids[1]};

        if (node_0 >= raw_mesh.nodes.size() || node_1 >= raw_mesh.nodes.size())
        {
            fail("Boundary contains an invalid node index.");
        }

        if (node_0 == node_1)
        {
            fail("Boundary contains a zero-length topological edge.");
        }

        add_neighbor(node_0, node_1);
        add_neighbor(node_1, node_0);
    }

    Index start_node{invalid_index};
    Index boundary_node_count{};

    for (Index node_id = 0; node_id < degrees.size(); ++node_id)
    {
        if (degrees[node_id] == 0)
        {
            continue;
        }

        if (degrees[node_id] != 2)
        {
            fail("Boundary is open: every boundary node must have exactly two neighbors.");
        }

        if (start_node == invalid_index)
        {
            start_node = node_id;
        }

        ++boundary_node_count;
    }

    if (start_node == invalid_index)
    {
        fail("Cannot compute boundary area: no boundary node was found.");
    }

    std::vector<bool> visited(raw_mesh.nodes.size(), false);

    Index previous_node{invalid_index};
    Index current_node{start_node};
    Index traversed_edge_count{};

    double twice_signed_area{};

    while (true)
    {
        if (visited[current_node])
        {
            fail("Boundary revisits a node before closing the contour.");
        }

        visited[current_node] = true;

        const std::array<Index, 2> &current_neighbors{neighbors[current_node]};

        Index next_node{current_neighbors[0]};

        if (previous_node != invalid_index)
        {
            if (current_neighbors[0] == previous_node)
            {
                next_node = current_neighbors[1];
            }
            else if (current_neighbors[1] != previous_node)
            {
                fail("Boundary connectivity is inconsistent.");
            }
        }

        const cfd::Node &current{raw_mesh.nodes[current_node]};
        const cfd::Node &next{raw_mesh.nodes[next_node]};

        twice_signed_area += current.x * next.y - next.x * current.y;

        ++traversed_edge_count;

        previous_node = current_node;
        current_node = next_node;

        if (current_node == start_node)
        {
            break;
        }

        if (traversed_edge_count >= raw_mesh.boundary_edges.size())
        {
            fail("Boundary traversal did not close.");
        }
    }

    if (traversed_edge_count != raw_mesh.boundary_edges.size())
    {
        fail("Boundary contains more than one disconnected closed contour.");
    }

    if (boundary_node_count != traversed_edge_count)
    {
        fail("Boundary node and edge counts are inconsistent for a simple closed contour.");
    }

    const double area{0.5 * std::abs(twice_signed_area)};

    if (!std::isfinite(area) || !(area > 0.0))
    {
        fail("Boundary has an invalid enclosed area.");
    }

    return area;
}

double compute_total_cell_area(const cfd::Mesh &mesh)
{
    return std::accumulate(mesh.cell_areas().begin(), mesh.cell_areas().end(), 0.0);
}

void check_mesh_geometry_invariants(const cfd::Mesh &mesh)
{
    require(mesh.cell_areas().size() == mesh.cell_count(), "Cell-area storage size is inconsistent.");
    require(mesh.cell_centers().size() == mesh.cell_count(), "Cell-center storage size is inconsistent.");

    require(mesh.face_lengths().size() == mesh.face_count(), "Face-length storage size is inconsistent.");
    require(mesh.face_centers().size() == mesh.face_count(), "Face-center storage size is inconsistent.");
    require(mesh.face_area_vectors().size() == mesh.face_count(), "Face-area-vector storage size is inconsistent.");
    require(mesh.cell_qualities().size() == mesh.cell_count(), "Cell-quality storage size is inconsistent.");

    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        const double length{mesh.face_lengths()[face_id]};
        const cfd::Vector2 &area_vector{mesh.face_area_vectors()[face_id]};

        const double area_vector_norm{std::hypot(area_vector.x, area_vector.y)};

        require_near(area_vector_norm, length, test_tolerance, "Face area-vector norm does not match face length.");

        const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};
        const cfd::Vector2 &face_center{mesh.face_centers()[face_id]};
        const cfd::Vector2 &owner_center{mesh.cell_centers()[adjacency.owner]};

        const double owner_orientation{area_vector.x * (face_center.x - owner_center.x) +
                                       area_vector.y * (face_center.y - owner_center.y)};

        require(owner_orientation > 0.0, "Face area vector is not oriented outward from its owner.");

        if (!adjacency.is_boundary())
        {
            const cfd::Vector2 &neighbor_center{mesh.cell_centers()[adjacency.neighbor]};

            const double owner_to_neighbor_orientation{area_vector.x * (neighbor_center.x - owner_center.x) +
                                                       area_vector.y * (neighbor_center.y - owner_center.y)};

            require(owner_to_neighbor_orientation > 0.0,
                    "Internal face orientation is inconsistent with owner and neighbor cells.");
        }
    }
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        if (mesh.cell_types()[cell_id] == cfd::CellType::Triangle)
        {
            const double quality{mesh.cell_qualities()[cell_id]};

            require(std::isfinite(quality), "Triangle quality is not finite.");

            require(quality > 0.0, "Triangle quality is not positive.");

            require(quality <= 1.0 + test_tolerance, "Triangle quality is greater than 1.");
        }
    }
}

void test_single_triangle()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};

    const double boundary_area{compute_single_closed_boundary_area(raw_mesh)};

    require_near(boundary_area, 0.5, test_tolerance, "Single-triangle boundary area is incorrect.");

    cfd::MeshBuildResult result{cfd::build_mesh(std::move(raw_mesh))};

    const cfd::Mesh &mesh{result.mesh};

    require(mesh.node_count() == 3, "Single triangle must contain three nodes.");
    require(mesh.cell_count() == 1, "Single triangle must contain one cell.");
    require(mesh.face_count() == 3, "Single triangle must contain three faces.");

    require_near(mesh.cell_areas()[0], 0.5, test_tolerance, "Single-triangle cell area is incorrect.");

    require_near(mesh.cell_centers()[0].x, 1.0 / 3.0, test_tolerance,
                 "Single-triangle centroid x-coordinate is incorrect.");

    require_near(mesh.cell_centers()[0].y, 1.0 / 3.0, test_tolerance,
                 "Single-triangle centroid y-coordinate is incorrect.");

    const double total_cell_area{compute_total_cell_area(mesh)};

    require_near(total_cell_area, boundary_area, test_tolerance,
                 "Single-triangle boundary area and cell-area sum differ.");

    check_mesh_geometry_invariants(mesh);

    const double expected_quality{std::sqrt(3.0) / 2.0};

    require_near(mesh.cell_qualities()[0], expected_quality, test_tolerance, "Single-triangle quality is incorrect.");
}
void test_single_triangle_statistics()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};
    cfd::MeshBuildResult result{cfd::build_mesh(std::move(raw_mesh))};

    const cfd::Mesh &mesh{result.mesh};

    const cfd::MeshStatistics statistics{cfd::compute_mesh_statistics(mesh)};

    require(statistics.internal_face_count == 0, "Single triangle must have no internal faces.");

    require(statistics.boundary_face_count == 3, "Single triangle must have three boundary faces.");

    require_near(statistics.total_cell_area, 0.5, test_tolerance, "Single-triangle total area statistic is incorrect.");

    require_near(statistics.cell_areas.minimum, 0.5, test_tolerance,
                 "Single-triangle minimum cell area statistic is incorrect.");

    require_near(statistics.cell_areas.mean, 0.5, test_tolerance,
                 "Single-triangle mean cell area statistic is incorrect.");

    require_near(statistics.cell_areas.maximum, 0.5, test_tolerance,
                 "Single-triangle maximum cell area statistic is incorrect.");

    const double expected_quality{std::sqrt(3.0) / 2.0};

    require_near(statistics.triangle_quality.minimum, expected_quality, test_tolerance,
                 "Single-triangle minimum quality statistic is incorrect.");

    require_near(statistics.triangle_quality.mean, expected_quality, test_tolerance,
                 "Single-triangle mean quality statistic is incorrect.");

    require_near(statistics.triangle_quality.maximum, expected_quality, test_tolerance,
                 "Single-triangle maximum quality statistic is incorrect.");

    require(statistics.worst_quality_cell == 0, "Single triangle must be its own worst-quality cell.");
}

void test_small_translated_equilateral_triangle()
{
    constexpr double origin_x{5.0};
    constexpr double origin_y{1.0};
    constexpr double side{0.01};

    const double height{std::sqrt(3.0) * side / 2.0};
    const double expected_area{std::sqrt(3.0) * side * side / 4.0};

    const double expected_centroid_x{origin_x + side / 2.0};
    const double expected_centroid_y{origin_y + height / 3.0};

    cfd::RawMeshData raw_mesh{cfd::test::make_equilateral_triangle_raw_mesh(origin_x, origin_y, side)};

    cfd::MeshBuildResult result{cfd::build_mesh(std::move(raw_mesh))};

    const cfd::Mesh &mesh{result.mesh};

    require_near(mesh.cell_areas()[0], expected_area, test_tolerance, "Translated small triangle area is incorrect.");

    require_near(mesh.cell_centers()[0].x, expected_centroid_x, test_tolerance,
                 "Translated small triangle centroid x-coordinate is incorrect.");

    require_near(mesh.cell_centers()[0].y, expected_centroid_y, test_tolerance,
                 "Translated small triangle centroid y-coordinate is incorrect.");

    require_near(mesh.cell_qualities()[0], 1.0, test_tolerance,
                 "Translated small equilateral triangle quality is incorrect.");

    check_mesh_geometry_invariants(mesh);
}
void test_reference_rectangle()
{
    const cfd::RectangleGeometry geometry{
        .length = 5.0,
        .height = 1.0,
    };

    const cfd::MeshGenerationOptions options{
        .mesh_size = 0.2,
        .cell_type = cfd::CellType::Triangle,
    };

    cfd::RawMeshData raw_mesh{cfd::generate_mesh(geometry, options)};

    const cfd::Index boundary_edge_count{raw_mesh.boundary_edges.size()};
    const double boundary_area{compute_single_closed_boundary_area(raw_mesh)};
    const double analytical_area{geometry.length * geometry.height};

    require_near(boundary_area, analytical_area, test_tolerance,
                 "Rectangle boundary area does not match analytical area.");

    cfd::MeshBuildResult result{cfd::build_mesh(std::move(raw_mesh))};

    const cfd::Mesh &mesh{result.mesh};

    const double total_cell_area{compute_total_cell_area(mesh)};

    require_near(total_cell_area, analytical_area, test_tolerance,
                 "Rectangle cell-area sum does not match analytical area.");

    require_near(total_cell_area, boundary_area, test_tolerance,
                 "Rectangle cell-area sum does not match Shoelace boundary area.");

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

    require(boundary_face_count == boundary_edge_count,
            "Constructed boundary-face count does not match imported boundary-edge count.");

    require(mesh.cell_faces().size() == 2 * internal_face_count + boundary_face_count,
            "Global face-incidence relation is not satisfied.");

    check_mesh_geometry_invariants(mesh);
}

} // namespace

int main()
{
    int failure_count{};

    failure_count += cfd::test::run_test("single triangle geometry", test_single_triangle);

    failure_count +=
        cfd::test::run_test("small translated equilateral triangle", test_small_translated_equilateral_triangle);

    failure_count += cfd::test::run_test("reference rectangle preprocessing", test_reference_rectangle);

    failure_count += cfd::test::run_test("single triangle statistics", test_single_triangle_statistics);

    return cfd::test::finish_tests(failure_count, "mesh geometry");
}