#include "cfd/mesh/Boundary.hpp"
#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/mesh/Types.hpp"
#include "cfd/meshing/BackwardFacingStepGeometry.hpp"
#include "cfd/meshing/GeometryInput.hpp"
#include "cfd/meshing/GmshMesher.hpp"
#include "cfd/meshing/RawMeshData.hpp"
#include "cfd/meshing/RawMeshValidation.hpp"
#include "cfd/meshing/RectangleGeometry.hpp"

#include "support/TestUtils.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace
{

using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws_with_message;

[[nodiscard]]
cfd::BoundaryId find_boundary_id(const cfd::Mesh &mesh, const std::string_view name)
{
    for (const cfd::BoundaryGroup &group : mesh.boundary_groups())
    {
        if (std::string_view{group.name} == name)
        {
            return group.id;
        }
    }

    return cfd::invalid_boundary_id;
}

[[nodiscard]]
double compute_boundary_length(const cfd::Mesh &mesh, const cfd::BoundaryId boundary_id)
{
    double total_length{};

    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (mesh.face_boundary_ids()[face_id] == boundary_id)
        {
            total_length += mesh.face_lengths()[face_id];
        }
    }

    return total_length;
}

[[nodiscard]]
cfd::BoundaryId find_raw_boundary_id(const cfd::RawMeshData &raw_mesh, const std::string_view name)
{
    for (const cfd::BoundaryGroup &group : raw_mesh.boundary_groups)
    {
        if (std::string_view{group.name} == name)
        {
            return group.id;
        }
    }
    return cfd::invalid_boundary_id;
}

[[nodiscard]]
bool contains_node(const cfd::RawMeshData &raw_mesh, const double x, const double y, const double tolerance)
{
    for (const cfd::Node &node : raw_mesh.nodes)
    {
        if (std::abs(node.x - x) <= tolerance && std::abs(node.y - y) <= tolerance)
        {
            return true;
        }
    }
    return false;
}

[[nodiscard]]
bool point_on_vertical_segment(const cfd::Node &point, const double x, const double y_min, const double y_max,
                               const double tolerance)
{
    return std::abs(point.x - x) <= tolerance && point.y >= y_min - tolerance && point.y <= y_max + tolerance;
}

[[nodiscard]]
bool point_on_horizontal_segment(const cfd::Node &point, const double y, const double x_min, const double x_max,
                                 const double tolerance)
{
    return std::abs(point.y - y) <= tolerance && point.x >= x_min - tolerance && point.x <= x_max + tolerance;
}

void require_invalid_meshing_input(const cfd::RectangleGeometry &geometry, const cfd::MeshGenerationOptions &options,
                                   const std::string_view expected_message, const std::string &failure_message)
{
    require_throws_with_message<std::invalid_argument>(
        [&geometry, &options]() { static_cast<void>(cfd::generate_mesh(geometry, options)); }, expected_message,
        failure_message);
}

void test_rejects_non_finite_rectangle_length()
{
    const cfd::RectangleGeometry geometry{
        .length = std::numeric_limits<double>::quiet_NaN(),
        .height = 1.0,
    };

    const cfd::MeshGenerationOptions options{
        .mesh_size = 0.2,
        .cell_type = cfd::CellType::Triangle,
    };

    require_invalid_meshing_input(geometry, options, "Rectangle length must be finite and positive.",
                                  "Gmsh mesher accepted a non-finite rectangle length.");
}

void test_rejects_infinite_rectangle_height()
{
    const cfd::RectangleGeometry geometry{
        .length = 5.0,
        .height = std::numeric_limits<double>::infinity(),
    };

    const cfd::MeshGenerationOptions options{
        .mesh_size = 0.2,
        .cell_type = cfd::CellType::Triangle,
    };

    require_invalid_meshing_input(geometry, options, "Rectangle height must be finite and positive.",
                                  "Gmsh mesher accepted an infinite rectangle height.");
}

void test_rejects_non_finite_mesh_size()
{
    const cfd::RectangleGeometry geometry{
        .length = 5.0,
        .height = 1.0,
    };

    const cfd::MeshGenerationOptions options{
        .mesh_size = std::numeric_limits<double>::quiet_NaN(),
        .cell_type = cfd::CellType::Triangle,
    };

    require_invalid_meshing_input(geometry, options, "Mesh size must be finite and positive.",
                                  "Gmsh mesher accepted a non-finite mesh size.");
}

void test_rejects_non_positive_dimensions()
{
    const cfd::MeshGenerationOptions options{
        .mesh_size = 0.2,
        .cell_type = cfd::CellType::Triangle,
    };

    require_invalid_meshing_input({.length = 0.0, .height = 1.0}, options,
                                  "Rectangle length must be finite and positive.",
                                  "Gmsh mesher accepted a zero rectangle length.");

    require_invalid_meshing_input({.length = 5.0, .height = -1.0}, options,
                                  "Rectangle height must be finite and positive.",
                                  "Gmsh mesher accepted a negative rectangle height.");
}

void test_rejects_non_positive_mesh_size()
{
    const cfd::RectangleGeometry geometry{
        .length = 5.0,
        .height = 1.0,
    };

    require_invalid_meshing_input(geometry,
                                  {
                                      .mesh_size = 0.0,
                                      .cell_type = cfd::CellType::Triangle,
                                  },
                                  "Mesh size must be finite and positive.", "Gmsh mesher accepted a zero mesh size.");

    require_invalid_meshing_input(geometry,
                                  {
                                      .mesh_size = -0.2,
                                      .cell_type = cfd::CellType::Triangle,
                                  },
                                  "Mesh size must be finite and positive.",
                                  "Gmsh mesher accepted a negative mesh size.");
}

void test_rejects_invalid_backward_facing_step_dimensions()
{
    const cfd::MeshGenerationOptions options{.mesh_size = 0.2, .cell_type = cfd::CellType::Triangle};
    require_throws_with_message<std::invalid_argument>(
        [&options]() {
            static_cast<void>(cfd::generate_mesh(
                cfd::BackwardFacingStepGeometry{
                    .upstream_length = 0.0, .downstream_length = 4.0, .channel_height = 1.0, .step_height = 0.5},
                options));
        },
        "upstream length", "Gmsh mesher accepted a zero backward-facing-step upstream length.");
    require_throws_with_message<std::invalid_argument>(
        [&options]() {
            static_cast<void>(cfd::generate_mesh(
                cfd::BackwardFacingStepGeometry{
                    .upstream_length = 1.0, .downstream_length = 4.0, .channel_height = 1.0, .step_height = 1.0},
                options));
        },
        "less than the channel height", "Gmsh mesher accepted an invalid backward-facing-step height.");
}

void test_generates_triangular_rectangle()
{
    const cfd::RectangleGeometry geometry{
        .length = 1.0,
        .height = 1.0,
    };

    const cfd::MeshGenerationOptions options{
        .mesh_size = 0.2,
        .cell_type = cfd::CellType::Triangle,
    };

    const cfd::RawMeshData raw_mesh{cfd::generate_mesh(geometry, options)};

    cfd::validate_raw_mesh(raw_mesh);

    require(!raw_mesh.cell_types.empty(), "Triangular Gmsh mesh contains no cells.");

    for (cfd::Index cell_id = 0; cell_id < raw_mesh.cell_types.size(); ++cell_id)
    {
        require(raw_mesh.cell_types[cell_id] == cfd::CellType::Triangle,
                "Triangular Gmsh mesh contains a non-triangular cell.");

        const cfd::Index node_count{raw_mesh.cell_node_offsets[cell_id + 1] - raw_mesh.cell_node_offsets[cell_id]};

        require(node_count == 3, "Triangular Gmsh cell does not contain three nodes.");
    }
}

void test_generates_quadrilateral_rectangle()
{
    const cfd::RectangleGeometry geometry{
        .length = 1.0,
        .height = 1.0,
    };

    const cfd::MeshGenerationOptions options{
        .mesh_size = 0.2,
        .cell_type = cfd::CellType::Quadrilateral,
    };

    const cfd::RawMeshData raw_mesh{cfd::generate_mesh(geometry, options)};

    cfd::validate_raw_mesh(raw_mesh);

    require(!raw_mesh.cell_types.empty(), "Quadrilateral Gmsh mesh contains no cells.");

    for (cfd::Index cell_id = 0; cell_id < raw_mesh.cell_types.size(); ++cell_id)
    {
        require(raw_mesh.cell_types[cell_id] == cfd::CellType::Quadrilateral,
                "Quadrilateral Gmsh mesh contains a non-quadrilateral cell.");

        const cfd::Index node_count{raw_mesh.cell_node_offsets[cell_id + 1] - raw_mesh.cell_node_offsets[cell_id]};

        require(node_count == 4, "Quadrilateral Gmsh cell does not contain four nodes.");
    }
}

void test_builds_quadrilateral_rectangle_end_to_end()
{
    constexpr double length{5.0};
    constexpr double height{1.0};
    constexpr double mesh_size{0.2};
    constexpr double tolerance{1.0e-10};

    const cfd::RectangleGeometry geometry{
        .length = length,
        .height = height,
    };

    const cfd::MeshGenerationOptions options{
        .mesh_size = mesh_size,
        .cell_type = cfd::CellType::Quadrilateral,
    };

    cfd::RawMeshData raw_mesh{cfd::generate_mesh(geometry, options)};
    cfd::MeshBuildResult build_result{cfd::build_mesh(std::move(raw_mesh))};

    const cfd::Mesh &mesh{build_result.mesh};

    require(mesh.node_count() > 0, "Quadrilateral end-to-end mesh contains no nodes.");
    require(mesh.cell_count() > 0, "Quadrilateral end-to-end mesh contains no cells.");
    require(mesh.face_count() > 0, "Quadrilateral end-to-end mesh contains no faces.");

    require(mesh.cell_faces().size() == 4 * mesh.cell_count(),
            "Quadrilateral mesh has an incorrect number of cell-face incidences.");

    double total_area{};

    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        require(mesh.cell_types()[cell_id] == cfd::CellType::Quadrilateral,
                "Quadrilateral end-to-end mesh contains a non-quadrilateral cell.");

        const double area{mesh.cell_areas()[cell_id]};

        require(std::isfinite(area) && area > 0.0, "Quadrilateral end-to-end mesh contains an invalid cell area.");

        total_area += area;

        const double quality{mesh.cell_qualities()[cell_id]};

        require(std::isfinite(quality) && quality > 0.0,
                "Quadrilateral end-to-end mesh contains an invalid cell quality.");
        require(quality <= 1.0 + tolerance, "Quadrilateral end-to-end mesh contains a cell quality greater than 1.");
    }

    require_near(total_area, length * height, tolerance, "Quadrilateral end-to-end mesh has an incorrect total area.");

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

    // Each quadrilateral contributes four cell-face incidences. Internal faces
    // are counted by two cells and boundary faces by one.
    require(4 * mesh.cell_count() == 2 * internal_face_count + boundary_face_count,
            "Quadrilateral end-to-end mesh violates the face-incidence identity.");

    const cfd::BoundaryId inlet_id{find_boundary_id(mesh, "inlet")};
    const cfd::BoundaryId wall_id{find_boundary_id(mesh, "wall")};
    const cfd::BoundaryId outlet_id{find_boundary_id(mesh, "outlet")};

    require(inlet_id != cfd::invalid_boundary_id, "Quadrilateral mesh does not contain an inlet boundary.");
    require(wall_id != cfd::invalid_boundary_id, "Quadrilateral mesh does not contain a wall boundary.");
    require(outlet_id != cfd::invalid_boundary_id, "Quadrilateral mesh does not contain an outlet boundary.");

    require_near(compute_boundary_length(mesh, inlet_id), height, tolerance,
                 "Quadrilateral inlet boundary length is incorrect.");
    require_near(compute_boundary_length(mesh, outlet_id), height, tolerance,
                 "Quadrilateral outlet boundary length is incorrect.");
    require_near(compute_boundary_length(mesh, wall_id), 2.0 * length, tolerance,
                 "Quadrilateral wall boundary length is incorrect.");
}

void test_rectangle_boundary_groups()
{
    constexpr double length{5.0};
    constexpr double height{1.0};
    constexpr double mesh_size{0.2};
    constexpr double tolerance{1.0e-10};

    const cfd::RectangleGeometry geometry{
        .length = length,
        .height = height,
    };

    const cfd::MeshGenerationOptions options{
        .mesh_size = mesh_size,
        .cell_type = cfd::CellType::Triangle,
    };

    cfd::RawMeshData raw_mesh{cfd::generate_mesh(geometry, options)};
    cfd::MeshBuildResult build_result{cfd::build_mesh(std::move(raw_mesh))};

    const cfd::Mesh &mesh{build_result.mesh};

    const cfd::BoundaryId inlet_id{find_boundary_id(mesh, "inlet")};
    const cfd::BoundaryId wall_id{find_boundary_id(mesh, "wall")};
    const cfd::BoundaryId outlet_id{find_boundary_id(mesh, "outlet")};

    require(inlet_id != cfd::invalid_boundary_id, "Rectangle mesh does not contain an inlet boundary.");
    require(wall_id != cfd::invalid_boundary_id, "Rectangle mesh does not contain a wall boundary.");
    require(outlet_id != cfd::invalid_boundary_id, "Rectangle mesh does not contain an outlet boundary.");

    require_near(compute_boundary_length(mesh, inlet_id), height, tolerance, "Inlet boundary length is incorrect.");
    require_near(compute_boundary_length(mesh, outlet_id), height, tolerance, "Outlet boundary length is incorrect.");
    require_near(compute_boundary_length(mesh, wall_id), 2.0 * length, tolerance, "Wall boundary length is incorrect.");
}

void verify_backward_facing_step_mesh(const cfd::CellType requested_cell_type, const bool use_variant_dispatch)
{
    constexpr cfd::BackwardFacingStepGeometry geometry{
        .upstream_length = 1.0,
        .downstream_length = 3.0,
        .channel_height = 1.0,
        .step_height = 0.5,
    };
    constexpr cfd::MeshGenerationOptions options{.mesh_size = 0.2, .cell_type = cfd::CellType::Triangle};
    constexpr double tolerance{1.0e-10};
    const cfd::MeshGenerationOptions selected_options{.mesh_size = options.mesh_size, .cell_type = requested_cell_type};
    cfd::RawMeshData raw_mesh{use_variant_dispatch ? cfd::generate_mesh(cfd::GeometryInput{geometry}, selected_options)
                                                   : cfd::generate_mesh(geometry, selected_options)};

    cfd::validate_raw_mesh(raw_mesh);
    require(!raw_mesh.cell_types.empty(), "Backward-facing-step mesh contains no cells.");
    for (const cfd::CellType cell_type : raw_mesh.cell_types)
    {
        require(cell_type == requested_cell_type,
                "Backward-facing-step mesh contains a cell of an unexpected topology.");
    }

    require(raw_mesh.boundary_groups.size() == 3,
            "Backward-facing-step mesh does not contain exactly three boundary groups.");
    for (cfd::BoundaryId boundary_id = 0; boundary_id < raw_mesh.boundary_groups.size(); ++boundary_id)
    {
        require(raw_mesh.boundary_groups[boundary_id].id == boundary_id,
                "Backward-facing-step boundary IDs are not compact.");
    }
    const cfd::BoundaryId inlet_id{find_raw_boundary_id(raw_mesh, "inlet")};
    const cfd::BoundaryId wall_id{find_raw_boundary_id(raw_mesh, "wall")};
    const cfd::BoundaryId outlet_id{find_raw_boundary_id(raw_mesh, "outlet")};
    require(inlet_id != cfd::invalid_boundary_id, "Backward-facing-step mesh has no inlet group.");
    require(wall_id != cfd::invalid_boundary_id, "Backward-facing-step mesh has no wall group.");
    require(outlet_id != cfd::invalid_boundary_id, "Backward-facing-step mesh has no outlet group.");

    const double outlet_x{geometry.upstream_length + geometry.downstream_length};
    for (const cfd::BoundaryEdge &edge : raw_mesh.boundary_edges)
    {
        const cfd::Node &first{raw_mesh.nodes[edge.node_ids[0]]};
        const cfd::Node &second{raw_mesh.nodes[edge.node_ids[1]]};
        if (edge.boundary_id == inlet_id)
        {
            require(
                point_on_vertical_segment(first, 0.0, geometry.step_height, geometry.channel_height, tolerance) &&
                    point_on_vertical_segment(second, 0.0, geometry.step_height, geometry.channel_height, tolerance),
                "Backward-facing-step inlet edge is not on the inlet segment.");
        }
        else if (edge.boundary_id == outlet_id)
        {
            require(point_on_vertical_segment(first, outlet_x, 0.0, geometry.channel_height, tolerance) &&
                        point_on_vertical_segment(second, outlet_x, 0.0, geometry.channel_height, tolerance),
                    "Backward-facing-step outlet edge is not on the outlet segment.");
        }
        else
        {
            require(edge.boundary_id == wall_id, "Backward-facing-step edge has an unknown boundary ID.");
            const bool on_top{point_on_horizontal_segment(first, geometry.channel_height, 0.0, outlet_x, tolerance) &&
                              point_on_horizontal_segment(second, geometry.channel_height, 0.0, outlet_x, tolerance)};
            const bool on_upstream_bottom{
                point_on_horizontal_segment(first, geometry.step_height, 0.0, geometry.upstream_length, tolerance) &&
                point_on_horizontal_segment(second, geometry.step_height, 0.0, geometry.upstream_length, tolerance)};
            const bool on_step{
                point_on_vertical_segment(first, geometry.upstream_length, 0.0, geometry.step_height, tolerance) &&
                point_on_vertical_segment(second, geometry.upstream_length, 0.0, geometry.step_height, tolerance)};
            const bool on_downstream_bottom{
                point_on_horizontal_segment(first, 0.0, geometry.upstream_length, outlet_x, tolerance) &&
                point_on_horizontal_segment(second, 0.0, geometry.upstream_length, outlet_x, tolerance)};
            require(on_top || on_upstream_bottom || on_step || on_downstream_bottom,
                    "Backward-facing-step wall edge is not on the exterior L-shaped contour.");
        }
    }

    constexpr std::array characteristic_points{
        cfd::Point2{0.0, geometry.step_height},
        cfd::Point2{geometry.upstream_length, geometry.step_height},
        cfd::Point2{geometry.upstream_length, 0.0},
        cfd::Point2{geometry.upstream_length + geometry.downstream_length, 0.0},
        cfd::Point2{geometry.upstream_length + geometry.downstream_length, geometry.channel_height},
        cfd::Point2{0.0, geometry.channel_height},
    };
    for (const cfd::Point2 &point : characteristic_points)
    {
        require(contains_node(raw_mesh, point.x, point.y, tolerance),
                "Backward-facing-step mesh is missing a characteristic contour point.");
    }

    cfd::MeshBuildResult build_result{cfd::build_mesh(std::move(raw_mesh))};
    const cfd::Mesh &mesh{build_result.mesh};
    double total_area{};
    for (const double area : mesh.cell_areas())
    {
        total_area += area;
    }
    const double analytical_area{geometry.upstream_length * (geometry.channel_height - geometry.step_height) +
                                 geometry.downstream_length * geometry.channel_height};
    require_near(total_area, analytical_area, tolerance, "Backward-facing-step mesh has an incorrect total area.");
    require_near(compute_boundary_length(mesh, inlet_id), geometry.channel_height - geometry.step_height, tolerance,
                 "Backward-facing-step inlet length is incorrect.");
    require_near(compute_boundary_length(mesh, outlet_id), geometry.channel_height, tolerance,
                 "Backward-facing-step outlet length is incorrect.");
    require_near(compute_boundary_length(mesh, wall_id),
                 2.0 * (geometry.upstream_length + geometry.downstream_length) + geometry.step_height, tolerance,
                 "Backward-facing-step wall length is incorrect.");
}

void test_generates_triangular_backward_facing_step()
{
    verify_backward_facing_step_mesh(cfd::CellType::Triangle, false);
}

void test_generates_quadrilateral_backward_facing_step()
{
    verify_backward_facing_step_mesh(cfd::CellType::Quadrilateral, true);
}

} // namespace

int main()
{
    int failure_count{};

    failure_count +=
        cfd::test::run_test("reject non-finite rectangle length", test_rejects_non_finite_rectangle_length);
    failure_count += cfd::test::run_test("reject infinite rectangle height", test_rejects_infinite_rectangle_height);
    failure_count += cfd::test::run_test("reject non-finite mesh size", test_rejects_non_finite_mesh_size);
    failure_count += cfd::test::run_test("reject non-positive dimensions", test_rejects_non_positive_dimensions);
    failure_count += cfd::test::run_test("reject non-positive mesh size", test_rejects_non_positive_mesh_size);
    failure_count += cfd::test::run_test("reject invalid backward-facing-step dimensions",
                                         test_rejects_invalid_backward_facing_step_dimensions);
    failure_count += cfd::test::run_test("generate triangular rectangle", test_generates_triangular_rectangle);
    failure_count += cfd::test::run_test("generate quadrilateral rectangle", test_generates_quadrilateral_rectangle);
    failure_count +=
        cfd::test::run_test("build quadrilateral rectangle end-to-end", test_builds_quadrilateral_rectangle_end_to_end);
    failure_count += cfd::test::run_test("rectangle boundary groups", test_rectangle_boundary_groups);
    failure_count +=
        cfd::test::run_test("generate triangular backward-facing step", test_generates_triangular_backward_facing_step);
    failure_count += cfd::test::run_test("generate quadrilateral backward-facing step",
                                         test_generates_quadrilateral_backward_facing_step);

    return cfd::test::finish_tests(failure_count, "Gmsh mesher");
}
