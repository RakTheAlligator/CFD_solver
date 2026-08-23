#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/meshing/GmshMesher.hpp"
#include "cfd/meshing/RawMeshData.hpp"
#include "cfd/meshing/RawMeshValidation.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace
{

constexpr double test_tolerance{1.0e-12};

[[noreturn]]
void fail(const std::string &message)
{
    throw std::runtime_error(message);
}

void require(const bool condition, const std::string &message)
{
    if (!condition)
    {
        fail(message);
    }
}

void require_near(const double actual, const double expected, const double tolerance, const std::string &message)
{
    const double scale{std::max({1.0, std::abs(actual), std::abs(expected)})};

    if (std::abs(actual - expected) > tolerance * scale)
    {
        fail(message + " (actual=" + std::to_string(actual) + ", expected=" + std::to_string(expected) + ")");
    }
}

template <typename Exception = std::exception, typename Function>
void require_throws(Function &&function, const std::string &message)
{
    try
    {
        std::forward<Function>(function)();
    }
    catch (const Exception &)
    {
        return;
    }
    catch (...)
    {
        fail(message + " (unexpected exception type).");
    }

    fail(message + " (no exception was thrown).");
}

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
}

cfd::RawMeshData make_single_triangle_raw_mesh()
{
    constexpr cfd::BoundaryId boundary_id{0};

    cfd::RawMeshData raw_mesh;

    raw_mesh.nodes = {
        {0.0, 0.0},
        {1.0, 0.0},
        {0.0, 1.0},
    };

    raw_mesh.cell_types = {
        cfd::CellType::Triangle,
    };

    raw_mesh.cell_nodes = {
        0,
        1,
        2,
    };

    raw_mesh.cell_node_offsets = {
        0,
        3,
    };

    raw_mesh.boundary_groups = {
        {boundary_id, "wall"},
    };

    raw_mesh.boundary_edges = {
        {{0, 1}, boundary_id},
        {{1, 2}, boundary_id},
        {{2, 0}, boundary_id},
    };

    return raw_mesh;
}

cfd::RawMeshData make_non_manifold_raw_mesh()
{
    constexpr cfd::BoundaryId boundary_id{0};

    cfd::RawMeshData raw_mesh;

    raw_mesh.nodes = {
        {0.0, 0.0}, {1.0, 0.0}, {0.5, 1.0}, {0.5, -1.0}, {0.5, 2.0},
    };

    raw_mesh.cell_types = {
        cfd::CellType::Triangle,
        cfd::CellType::Triangle,
        cfd::CellType::Triangle,
    };

    raw_mesh.cell_nodes = {
        0, 1, 2, 1, 0, 3, 0, 1, 4,
    };

    raw_mesh.cell_node_offsets = {
        0,
        3,
        6,
        9,
    };

    raw_mesh.boundary_groups = {
        {boundary_id, "wall"},
    };

    raw_mesh.boundary_edges = {
        {{1, 2}, boundary_id}, {{2, 0}, boundary_id}, {{0, 3}, boundary_id},
        {{3, 1}, boundary_id}, {{1, 4}, boundary_id}, {{4, 0}, boundary_id},
    };

    return raw_mesh;
}

void test_single_triangle()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};

    const double boundary_area{compute_single_closed_boundary_area(raw_mesh)};

    require_near(boundary_area, 0.5, test_tolerance, "Single-triangle boundary area is incorrect.");

    cfd::Mesh mesh{cfd::build_mesh(std::move(raw_mesh))};

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

    cfd::Mesh mesh{cfd::build_mesh(std::move(raw_mesh))};

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

    require(mesh.cell_nodes().size() == 2 * internal_face_count + boundary_face_count,
            "Global face-incidence relation is not satisfied.");

    check_mesh_geometry_invariants(mesh);
}

void test_rejects_inconsistent_cell_offsets()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};

    raw_mesh.cell_node_offsets.back() = 2;

    require_throws<std::runtime_error>([&raw_mesh]() { cfd::validate_raw_mesh(raw_mesh); },
                                       "Raw mesh validation accepted inconsistent cell offsets.");
}

void test_rejects_duplicate_node_in_cell()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};

    raw_mesh.cell_nodes[2] = raw_mesh.cell_nodes[1];

    require_throws<std::runtime_error>([&raw_mesh]() { cfd::validate_raw_mesh(raw_mesh); },
                                       "Raw mesh validation accepted a duplicated node inside a cell.");
}

void test_rejects_invalid_node_index()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};

    raw_mesh.cell_nodes[2] = raw_mesh.nodes.size();

    require_throws<std::runtime_error>([&raw_mesh]() { cfd::validate_raw_mesh(raw_mesh); },
                                       "Raw mesh validation accepted an out-of-range node index.");
}

void test_rejects_duplicate_boundary_edge()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};

    raw_mesh.boundary_edges.push_back({
        {1, 0},
        raw_mesh.boundary_groups[0].id,
    });

    require_throws<std::runtime_error>([&raw_mesh]() { cfd::validate_raw_mesh(raw_mesh); },
                                       "Raw mesh validation accepted a duplicated boundary edge.");
}

void test_rejects_open_boundary()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};

    raw_mesh.boundary_edges.pop_back();

    require_throws<std::runtime_error>(
        [&raw_mesh]() { static_cast<void>(compute_single_closed_boundary_area(raw_mesh)); },
        "Boundary-area reconstruction accepted an open contour.");

    require_throws<std::runtime_error>([&raw_mesh]() { static_cast<void>(cfd::build_mesh(std::move(raw_mesh))); },
                                       "Mesh construction accepted an incomplete physical boundary.");
}

void test_rejects_zero_area_triangle()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};

    raw_mesh.nodes[0] = {0.0, 0.0};
    raw_mesh.nodes[1] = {1.0, 0.0};
    raw_mesh.nodes[2] = {2.0, 0.0};

    cfd::validate_raw_mesh(raw_mesh);

    require_throws<std::runtime_error>([&raw_mesh]() { static_cast<void>(cfd::build_mesh(std::move(raw_mesh))); },
                                       "Mesh construction accepted a zero-area triangle.");
}

void test_rejects_non_manifold_face()
{
    cfd::RawMeshData raw_mesh{make_non_manifold_raw_mesh()};

    cfd::validate_raw_mesh(raw_mesh);

    require_throws<std::runtime_error>([&raw_mesh]() { static_cast<void>(cfd::build_mesh(std::move(raw_mesh))); },
                                       "Mesh construction accepted a face shared by more than two cells.");
}
void test_rejects_invalid_boundary_group_id()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};

    raw_mesh.boundary_groups[0].id = cfd::invalid_boundary_id;

    for (cfd::BoundaryEdge &edge : raw_mesh.boundary_edges)
    {
        edge.boundary_id = cfd::invalid_boundary_id;
    }

    require_throws<std::runtime_error>([&raw_mesh]() { cfd::validate_raw_mesh(raw_mesh); },
                                       "Raw mesh validation accepted invalid_boundary_id as a physical boundary ID.");
}

template <typename TestFunction> int run_test(const std::string_view name, TestFunction test_function)
{
    try
    {
        test_function();
        std::cout << "[PASS] " << name << '\n';
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "[FAIL] " << name << ": " << error.what() << '\n';
        return 1;
    }
}

} // namespace

int main()
{
    int failure_count{};

    failure_count += run_test("single triangle geometry", test_single_triangle);
    failure_count += run_test("reference rectangle preprocessing", test_reference_rectangle);

    failure_count += run_test("reject inconsistent cell offsets", test_rejects_inconsistent_cell_offsets);
    failure_count += run_test("reject duplicate node in cell", test_rejects_duplicate_node_in_cell);
    failure_count += run_test("reject invalid node index", test_rejects_invalid_node_index);
    failure_count += run_test("reject duplicate boundary edge", test_rejects_duplicate_boundary_edge);
    failure_count += run_test("reject open boundary", test_rejects_open_boundary);
    failure_count += run_test("reject zero-area triangle", test_rejects_zero_area_triangle);
    failure_count += run_test("reject non-manifold face", test_rejects_non_manifold_face);
    failure_count += run_test("reject invalid boundary group ID", test_rejects_invalid_boundary_group_id);

    if (failure_count == 0)
    {
        std::cout << "[PASS] All mesh preprocessing tests passed.\n";
        return 0;
    }

    std::cerr << "[FAIL] " << failure_count << " mesh preprocessing test(s) failed.\n";
    return 1;
}