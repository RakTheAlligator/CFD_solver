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

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

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

[[nodiscard]]
std::vector<double> x_coordinates_on_line(const cfd::RawMeshData &raw_mesh, const double y, const double tolerance)
{
    std::vector<double> coordinates;
    for (const cfd::Node &node : raw_mesh.nodes)
    {
        if (std::abs(node.y - y) <= tolerance)
        {
            coordinates.push_back(node.x);
        }
    }
    std::ranges::sort(coordinates);
    const auto duplicate_begin{std::ranges::unique(coordinates, [tolerance](const double left, const double right) {
                                   return std::abs(left - right) <= tolerance;
                               }).begin()};
    coordinates.erase(duplicate_begin, coordinates.end());
    return coordinates;
}

[[nodiscard]]
std::vector<double> y_coordinates_on_line(const cfd::RawMeshData &raw_mesh, const double x, const double tolerance)
{
    std::vector<double> coordinates;
    for (const cfd::Node &node : raw_mesh.nodes)
    {
        if (std::abs(node.x - x) <= tolerance)
        {
            coordinates.push_back(node.y);
        }
    }
    std::ranges::sort(coordinates);
    const auto duplicate_begin{std::ranges::unique(coordinates, [tolerance](const double left, const double right) {
                                   return std::abs(left - right) <= tolerance;
                               }).begin()};
    coordinates.erase(duplicate_begin, coordinates.end());
    return coordinates;
}

[[nodiscard]]
double minimum_adjacent_spacing(const std::vector<double> &coordinates)
{
    require(coordinates.size() >= 2, "Not enough line nodes to measure mesh spacing.");
    double minimum{std::numeric_limits<double>::max()};
    for (std::size_t index = 1; index < coordinates.size(); ++index)
    {
        minimum = std::min(minimum, coordinates[index] - coordinates[index - 1]);
    }
    return minimum;
}

[[nodiscard]]
double maximum_adjacent_spacing(const std::vector<double> &coordinates)
{
    require(coordinates.size() >= 2, "Not enough line nodes to measure mesh spacing.");
    double maximum{};
    for (std::size_t index = 1; index < coordinates.size(); ++index)
    {
        maximum = std::max(maximum, coordinates[index] - coordinates[index - 1]);
    }
    return maximum;
}

[[nodiscard]]
double spacing_adjacent_to(const std::vector<double> &coordinates, const double location, const double tolerance)
{
    const auto iterator{std::ranges::find_if(coordinates, [location, tolerance](const double coordinate) {
        return std::abs(coordinate - location) <= tolerance;
    })};
    require(iterator != coordinates.end(), "Refinement location is absent from the mesh line.");
    const std::size_t index{static_cast<std::size_t>(std::distance(coordinates.begin(), iterator))};
    double spacing{std::numeric_limits<double>::max()};
    if (index > 0)
    {
        spacing = location - coordinates[index - 1];
    }
    if (index + 1 < coordinates.size())
    {
        spacing = std::min(spacing, coordinates[index + 1] - location);
    }
    return spacing;
}

[[nodiscard]]
double next_coordinate_after(const std::vector<double> &coordinates, const double location, const double tolerance)
{
    const auto iterator{std::ranges::find_if(
        coordinates, [location, tolerance](const double coordinate) { return coordinate > location + tolerance; })};
    require(iterator != coordinates.end(), "No structured line exists after the requested location.");
    return *iterator;
}

[[nodiscard]]
double maximum_successive_spacing_ratio(const std::vector<double> &coordinates)
{
    require(coordinates.size() >= 3, "Not enough line nodes to measure consecutive mesh-spacing ratios.");
    double maximum_ratio{1.0};
    double previous_spacing{coordinates[1] - coordinates[0]};
    for (std::size_t index = 2; index < coordinates.size(); ++index)
    {
        const double spacing{coordinates[index] - coordinates[index - 1]};
        maximum_ratio = std::max(maximum_ratio, std::max(spacing / previous_spacing, previous_spacing / spacing));
        previous_spacing = spacing;
    }
    return maximum_ratio;
}

struct StepTransitionMeasurements
{
    double adjacent_spacing{};
    double transition_end_spacing{};
    double bulk_spacing{};
    double transition_width{};
    std::size_t transition_cell_count{};
};

[[nodiscard]]
StepTransitionMeasurements measure_downstream_step_transition(const cfd::RawMeshData &raw_mesh,
                                                              const cfd::BackwardFacingStepGeometry &geometry,
                                                              const double tolerance)
{
    const std::vector<double> coordinates{x_coordinates_on_line(raw_mesh, geometry.step_height, tolerance)};
    const auto step_iterator{std::ranges::find_if(coordinates, [&geometry, tolerance](const double coordinate) {
        return std::abs(coordinate - geometry.upstream_length) <= tolerance;
    })};
    require(step_iterator != coordinates.end(), "Backward-facing-step location is absent from the structured line.");
    const std::size_t step_index{static_cast<std::size_t>(std::distance(coordinates.begin(), step_iterator))};
    require(step_index + 3 < coordinates.size(), "Downstream structured line is too short to identify a bulk region.");

    std::vector<double> downstream_spacings;
    downstream_spacings.reserve(coordinates.size() - step_index - 1);
    for (std::size_t index = step_index + 1; index < coordinates.size(); ++index)
    {
        downstream_spacings.push_back(coordinates[index] - coordinates[index - 1]);
    }

    const double bulk_spacing{downstream_spacings.back()};
    const double uniform_tolerance{std::sqrt(std::numeric_limits<double>::epsilon()) * std::max(1.0, bulk_spacing)};
    std::size_t first_bulk_spacing{downstream_spacings.size()};
    for (std::size_t candidate = 1; candidate + 1 < downstream_spacings.size(); ++candidate)
    {
        const bool uniform_tail{
            std::ranges::all_of(downstream_spacings.begin() + static_cast<std::ptrdiff_t>(candidate),
                                downstream_spacings.end(), [bulk_spacing, uniform_tolerance](const double spacing) {
                                    return std::abs(spacing - bulk_spacing) <= uniform_tolerance;
                                })};
        if (uniform_tail)
        {
            first_bulk_spacing = candidate;
            break;
        }
    }
    require(first_bulk_spacing < downstream_spacings.size(),
            "No uniform downstream bulk was found after the local step transition.");

    return {
        .adjacent_spacing = downstream_spacings.front(),
        .transition_end_spacing = downstream_spacings[first_bulk_spacing - 1],
        .bulk_spacing = bulk_spacing,
        .transition_width = coordinates[step_index + first_bulk_spacing] - geometry.upstream_length,
        .transition_cell_count = first_bulk_spacing,
    };
}

void require_representative_lines_respect_growth_rate(const cfd::RawMeshData &raw_mesh,
                                                      const cfd::BackwardFacingStepGeometry &geometry,
                                                      const double maximum_growth_rate, const double tolerance)
{
    const double outlet_x{geometry.upstream_length + geometry.downstream_length};
    const double ratio_tolerance{std::sqrt(std::numeric_limits<double>::epsilon()) * maximum_growth_rate};
    const std::vector<double> outlet_coordinates{y_coordinates_on_line(raw_mesh, outlet_x, tolerance)};
    std::vector<double> lower_outlet_coordinates;
    std::vector<double> upper_outlet_coordinates;
    std::ranges::copy_if(
        outlet_coordinates, std::back_inserter(lower_outlet_coordinates),
        [&geometry, tolerance](const double coordinate) { return coordinate <= geometry.step_height + tolerance; });
    std::ranges::copy_if(
        outlet_coordinates, std::back_inserter(upper_outlet_coordinates),
        [&geometry, tolerance](const double coordinate) { return coordinate >= geometry.step_height - tolerance; });
    const std::array line_coordinates{
        x_coordinates_on_line(raw_mesh, geometry.step_height, tolerance),
        x_coordinates_on_line(raw_mesh, geometry.channel_height, tolerance),
        std::move(lower_outlet_coordinates),
        std::move(upper_outlet_coordinates),
    };
    for (std::size_t line_index = 0; line_index < line_coordinates.size(); ++line_index)
    {
        const double ratio{maximum_successive_spacing_ratio(line_coordinates.at(line_index))};
        require(ratio <= maximum_growth_rate + ratio_tolerance,
                "Generated structured line " + std::to_string(line_index) + " has spacing ratio " +
                    std::to_string(ratio) + " above maximumGrowthRate " + std::to_string(maximum_growth_rate) + ".");
    }
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
        const bool on_upstream_interface{
            point_on_vertical_segment(first, geometry.upstream_length, geometry.step_height, geometry.channel_height,
                                      tolerance) &&
            point_on_vertical_segment(second, geometry.upstream_length, geometry.step_height, geometry.channel_height,
                                      tolerance)};
        const bool on_downstream_interface{
            point_on_horizontal_segment(first, geometry.step_height, geometry.upstream_length, outlet_x, tolerance) &&
            point_on_horizontal_segment(second, geometry.step_height, geometry.upstream_length, outlet_x, tolerance)};
        require(!on_upstream_interface && !on_downstream_interface,
                "A backward-facing-step block interface was classified as a CFD boundary.");
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

    if (requested_cell_type == cfd::CellType::Quadrilateral)
    {
        for (cfd::Index first_id = 0; first_id < raw_mesh.nodes.size(); ++first_id)
        {
            for (cfd::Index second_id = first_id + 1; second_id < raw_mesh.nodes.size(); ++second_id)
            {
                const cfd::Node &first{raw_mesh.nodes[first_id]};
                const cfd::Node &second{raw_mesh.nodes[second_id]};
                require(std::abs(first.x - second.x) > tolerance || std::abs(first.y - second.y) > tolerance,
                        "Structured backward-facing-step interfaces contain duplicate nodes.");
            }
        }

        const std::vector<double> streamwise_coordinates{
            x_coordinates_on_line(raw_mesh, geometry.step_height, tolerance)};
        const std::vector<double> wall_normal_coordinates{y_coordinates_on_line(raw_mesh, outlet_x, tolerance)};
        require(spacing_adjacent_to(streamwise_coordinates, geometry.upstream_length, tolerance) <
                    maximum_adjacent_spacing(streamwise_coordinates),
                "Automatic mesh is not refined near the backward-facing step.");
        require(minimum_adjacent_spacing(wall_normal_coordinates) < maximum_adjacent_spacing(wall_normal_coordinates),
                "Automatic mesh is not refined toward solid walls.");
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

void test_backward_facing_step_refinement_options_change_spacing()
{
    constexpr cfd::BackwardFacingStepGeometry geometry{
        .upstream_length = 1.0,
        .downstream_length = 3.0,
        .channel_height = 1.0,
        .step_height = 0.5,
    };
    constexpr cfd::MeshGenerationOptions generation{.mesh_size = 0.2, .cell_type = cfd::CellType::Quadrilateral};
    constexpr cfd::AutomaticMeshingOptions automatic;
    constexpr double tolerance{1.0e-10};

    cfd::BackwardFacingStepMeshingOptions no_step_refinement;
    no_step_refinement.step_refinement_factor = 1.0;
    const cfd::RawMeshData uniform_step_mesh{cfd::generate_mesh(geometry, generation, automatic, no_step_refinement)};
    cfd::BackwardFacingStepMeshingOptions stronger_step_refinement;
    stronger_step_refinement.step_refinement_factor = 0.25;
    const cfd::RawMeshData refined_step_mesh{
        cfd::generate_mesh(geometry, generation, automatic, stronger_step_refinement)};

    const double uniform_step_spacing{
        spacing_adjacent_to(x_coordinates_on_line(uniform_step_mesh, geometry.step_height, tolerance),
                            geometry.upstream_length, tolerance)};
    const double refined_step_spacing{
        spacing_adjacent_to(x_coordinates_on_line(refined_step_mesh, geometry.step_height, tolerance),
                            geometry.upstream_length, tolerance)};
    require(refined_step_spacing < uniform_step_spacing,
            "Changing stepRefinementFactor did not strengthen step refinement.");

    cfd::AutomaticMeshingOptions no_wall_refinement{automatic};
    no_wall_refinement.wall_refinement_factor = 1.0;
    const cfd::RawMeshData uniform_wall_mesh{
        cfd::generate_mesh(geometry, generation, no_wall_refinement, cfd::BackwardFacingStepMeshingOptions{})};
    cfd::AutomaticMeshingOptions stronger_wall_refinement{automatic};
    stronger_wall_refinement.wall_refinement_factor = 0.25;
    const cfd::RawMeshData refined_wall_mesh{
        cfd::generate_mesh(geometry, generation, stronger_wall_refinement, cfd::BackwardFacingStepMeshingOptions{})};
    const double outlet_x{geometry.upstream_length + geometry.downstream_length};
    const double uniform_wall_spacing{
        minimum_adjacent_spacing(y_coordinates_on_line(uniform_wall_mesh, outlet_x, tolerance))};
    const double refined_wall_spacing{
        minimum_adjacent_spacing(y_coordinates_on_line(refined_wall_mesh, outlet_x, tolerance))};
    require(refined_wall_spacing < uniform_wall_spacing,
            "Changing wallRefinementFactor did not strengthen wall refinement.");
}

void test_backward_facing_step_internal_interface_refinement_relaxes_downstream()
{
    constexpr cfd::BackwardFacingStepGeometry geometry{
        .upstream_length = 1.0,
        .downstream_length = 3.0,
        .channel_height = 1.0,
        .step_height = 0.5,
    };
    constexpr cfd::MeshGenerationOptions generation{.mesh_size = 0.2, .cell_type = cfd::CellType::Quadrilateral};
    constexpr double tolerance{1.0e-10};
    const cfd::RawMeshData raw_mesh{cfd::generate_mesh(geometry, generation)};

    const std::vector<double> streamwise_coordinates{x_coordinates_on_line(raw_mesh, geometry.step_height, tolerance)};
    const double immediate_downstream_x{
        next_coordinate_after(streamwise_coordinates, geometry.upstream_length, tolerance)};
    const double outlet_x{geometry.upstream_length + geometry.downstream_length};
    const double immediate_interface_spacing{spacing_adjacent_to(
        y_coordinates_on_line(raw_mesh, immediate_downstream_x, tolerance), geometry.step_height, tolerance)};
    const double far_interface_spacing{
        spacing_adjacent_to(y_coordinates_on_line(raw_mesh, outlet_x, tolerance), geometry.step_height, tolerance)};

    require(far_interface_spacing > immediate_interface_spacing,
            "Wall refinement inherited at the backward-facing step did not relax downstream.");
}

void test_backward_facing_step_refinement_is_local_downstream()
{
    constexpr cfd::MeshGenerationOptions generation{.mesh_size = 0.2, .cell_type = cfd::CellType::Quadrilateral};
    constexpr cfd::AutomaticMeshingOptions automatic;
    constexpr cfd::BackwardFacingStepMeshingOptions step_options;
    constexpr double tolerance{1.0e-10};
    constexpr cfd::BackwardFacingStepGeometry shorter_geometry{
        .upstream_length = 1.0,
        .downstream_length = 4.0,
        .channel_height = 1.0,
        .step_height = 0.5,
    };
    constexpr cfd::BackwardFacingStepGeometry longer_geometry{
        .upstream_length = 1.0,
        .downstream_length = 8.0,
        .channel_height = 1.0,
        .step_height = 0.5,
    };

    const cfd::RawMeshData shorter_mesh{cfd::generate_mesh(shorter_geometry, generation, automatic, step_options)};
    const cfd::RawMeshData longer_mesh{cfd::generate_mesh(longer_geometry, generation, automatic, step_options)};
    const StepTransitionMeasurements shorter{
        measure_downstream_step_transition(shorter_mesh, shorter_geometry, tolerance)};
    const StepTransitionMeasurements longer{
        measure_downstream_step_transition(longer_mesh, longer_geometry, tolerance)};

    const double target_step_spacing{generation.mesh_size *
                                     std::min(automatic.wall_refinement_factor, step_options.step_refinement_factor)};
    require(shorter.adjacent_spacing <= target_step_spacing + tolerance,
            "Step-adjacent spacing exceeds the effective refinement target.");
    require(shorter.adjacent_spacing < shorter.bulk_spacing,
            "Downstream step refinement does not return to a coarser bulk spacing.");
    require(shorter.bulk_spacing > 0.8 * generation.mesh_size && shorter.bulk_spacing <= generation.mesh_size,
            "Downstream bulk spacing is not close to the base mesh size.");
    require(longer.transition_width < 1.25 * shorter.transition_width,
            "Increasing downstreamLength extended the local step transition toward the outlet.");
}

void test_backward_facing_step_refinement_layers_control_transition()
{
    constexpr cfd::BackwardFacingStepGeometry geometry{
        .upstream_length = 1.0,
        .downstream_length = 8.0,
        .channel_height = 1.0,
        .step_height = 0.5,
    };
    constexpr cfd::MeshGenerationOptions generation{.mesh_size = 0.2, .cell_type = cfd::CellType::Quadrilateral};
    constexpr cfd::AutomaticMeshingOptions automatic;
    constexpr double tolerance{1.0e-10};

    cfd::BackwardFacingStepMeshingOptions short_transition;
    short_transition.step_refinement_layers = 3;
    const cfd::RawMeshData short_mesh{cfd::generate_mesh(geometry, generation, automatic, short_transition)};

    cfd::BackwardFacingStepMeshingOptions long_transition;
    long_transition.step_refinement_layers = 8;
    const cfd::RawMeshData long_mesh{cfd::generate_mesh(geometry, generation, automatic, long_transition)};

    const StepTransitionMeasurements short_measurements{
        measure_downstream_step_transition(short_mesh, geometry, tolerance)};
    const StepTransitionMeasurements long_measurements{
        measure_downstream_step_transition(long_mesh, geometry, tolerance)};
    require(long_measurements.transition_cell_count > short_measurements.transition_cell_count,
            "Increasing stepRefinementLayers did not add cells to the local step transition.");
    require(long_measurements.transition_width > short_measurements.transition_width,
            "Increasing stepRefinementLayers did not widen the local step transition.");
}

void test_backward_facing_step_respects_actual_maximum_growth_rate()
{
    constexpr cfd::BackwardFacingStepGeometry geometry{
        .upstream_length = 1.0,
        .downstream_length = 6.0,
        .channel_height = 1.0,
        .step_height = 0.5,
    };
    constexpr cfd::MeshGenerationOptions generation{.mesh_size = 0.2, .cell_type = cfd::CellType::Quadrilateral};
    constexpr cfd::BackwardFacingStepMeshingOptions step_options;
    constexpr double tolerance{1.0e-10};

    cfd::AutomaticMeshingOptions restrictive;
    restrictive.maximum_growth_rate = 1.05;
    const cfd::RawMeshData restrictive_mesh{cfd::generate_mesh(geometry, generation, restrictive, step_options)};

    cfd::AutomaticMeshingOptions permissive;
    permissive.maximum_growth_rate = 1.4;
    const cfd::RawMeshData permissive_mesh{cfd::generate_mesh(geometry, generation, permissive, step_options)};

    require_representative_lines_respect_growth_rate(restrictive_mesh, geometry, restrictive.maximum_growth_rate,
                                                     tolerance);
    require_representative_lines_respect_growth_rate(permissive_mesh, geometry, permissive.maximum_growth_rate,
                                                     tolerance);
    const StepTransitionMeasurements restrictive_transition{
        measure_downstream_step_transition(restrictive_mesh, geometry, tolerance)};
    const StepTransitionMeasurements permissive_transition{
        measure_downstream_step_transition(permissive_mesh, geometry, tolerance)};
    const auto require_transition_bulk_join = [tolerance](const StepTransitionMeasurements &measurements,
                                                          const double maximum_growth_rate) {
        const double ratio{std::max(measurements.transition_end_spacing / measurements.bulk_spacing,
                                    measurements.bulk_spacing / measurements.transition_end_spacing)};
        require(ratio <= maximum_growth_rate + tolerance,
                "The generated step-transition-to-bulk join exceeds maximumGrowthRate.");
    };
    require_transition_bulk_join(restrictive_transition, restrictive.maximum_growth_rate);
    require_transition_bulk_join(permissive_transition, permissive.maximum_growth_rate);
    require(restrictive_mesh.cell_types.size() > permissive_mesh.cell_types.size(),
            "A more restrictive maximumGrowthRate did not increase the structured transition resolution.");
}

void test_disabled_automatic_meshing_uses_legacy_recombination()
{
    constexpr cfd::BackwardFacingStepGeometry geometry{
        .upstream_length = 1.0,
        .downstream_length = 3.0,
        .channel_height = 1.0,
        .step_height = 0.5,
    };
    constexpr cfd::MeshGenerationOptions generation{.mesh_size = 0.2, .cell_type = cfd::CellType::Quadrilateral};
    constexpr cfd::BackwardFacingStepMeshingOptions step_options;
    cfd::AutomaticMeshingOptions legacy_options;
    legacy_options.enabled = false;

    cfd::RawMeshData legacy_mesh{cfd::generate_mesh(geometry, generation, legacy_options, step_options)};
    cfd::validate_raw_mesh(legacy_mesh);
    for (const cfd::CellType cell_type : legacy_mesh.cell_types)
    {
        require(cell_type == cfd::CellType::Quadrilateral,
                "Legacy backward-facing-step recombination produced a non-quadrilateral cell.");
    }

    const cfd::RawMeshData structured_mesh{
        cfd::generate_mesh(geometry, generation, cfd::AutomaticMeshingOptions{}, step_options)};
    require(legacy_mesh.nodes.size() != structured_mesh.nodes.size() ||
                legacy_mesh.cell_types.size() != structured_mesh.cell_types.size(),
            "Disabling automatic meshing did not select the legacy recombination path.");

    cfd::MeshBuildResult build_result{cfd::build_mesh(std::move(legacy_mesh))};
    require(build_result.mesh.cell_count() > 0, "Legacy backward-facing-step mesh could not be built.");
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
    failure_count += cfd::test::run_test("backward-facing-step refinement options change spacing",
                                         test_backward_facing_step_refinement_options_change_spacing);
    failure_count += cfd::test::run_test("backward-facing-step internal-interface refinement relaxes downstream",
                                         test_backward_facing_step_internal_interface_refinement_relaxes_downstream);
    failure_count += cfd::test::run_test("backward-facing-step refinement is local downstream",
                                         test_backward_facing_step_refinement_is_local_downstream);
    failure_count += cfd::test::run_test("backward-facing-step refinement layers control transition",
                                         test_backward_facing_step_refinement_layers_control_transition);
    failure_count += cfd::test::run_test("backward-facing-step respects actual maximum growth rate",
                                         test_backward_facing_step_respects_actual_maximum_growth_rate);
    failure_count += cfd::test::run_test("disabled automatic meshing uses legacy recombination",
                                         test_disabled_automatic_meshing_uses_legacy_recombination);

    return cfd::test::finish_tests(failure_count, "Gmsh mesher");
}
