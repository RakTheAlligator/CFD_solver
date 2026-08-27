#include "mesh_build/GeometryBuilder.hpp"

#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Node.hpp"
#include "cfd/mesh/Types.hpp"
#include "cfd/mesh/Vector2.hpp"
#include "cfd/meshing/RawMeshData.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace cfd::detail
{

namespace
{

[[noreturn]]
void throw_geometry_build_error(const std::string &message)
{
    throw std::runtime_error("Geometry construction failed: " + message);
}

[[nodiscard]]
double compute_triangle_quality(const RawMeshData &raw_mesh, const Index cell_id, const double area)
{
    const Index cell_node_begin_offset{raw_mesh.cell_node_offsets[cell_id]};

    const Node &node_0{raw_mesh.nodes[raw_mesh.cell_nodes[cell_node_begin_offset]]};
    const Node &node_1{raw_mesh.nodes[raw_mesh.cell_nodes[cell_node_begin_offset + 1]]};
    const Node &node_2{raw_mesh.nodes[raw_mesh.cell_nodes[cell_node_begin_offset + 2]]};

    const double dx_01{node_1.x - node_0.x};
    const double dy_01{node_1.y - node_0.y};

    const double dx_12{node_2.x - node_1.x};
    const double dy_12{node_2.y - node_1.y};

    const double dx_20{node_0.x - node_2.x};
    const double dy_20{node_0.y - node_2.y};

    const double squared_length_sum{dx_01 * dx_01 + dy_01 * dy_01 + dx_12 * dx_12 + dy_12 * dy_12 + dx_20 * dx_20 +
                                    dy_20 * dy_20};

    if (!std::isfinite(squared_length_sum) || !(squared_length_sum > 0.0))
    {
        throw_geometry_build_error("cell " + std::to_string(cell_id) + " has invalid edge lengths.");
    }

    // Normalized triangle quality:
    //
    //     q = 4 sqrt(3) A / (l_01^2 + l_12^2 + l_20^2)
    //
    // An equilateral triangle gives q = 1; distorted or degenerate triangles
    // approach zero.
    return 4.0 * std::sqrt(3.0) * area / squared_length_sum;
}

[[nodiscard]]
double compute_quadrilateral_quality(const RawMeshData &raw_mesh, const Index cell_id)
{
    constexpr Index node_count{4};

    const Index cell_node_begin_offset{raw_mesh.cell_node_offsets[cell_id]};

    double minimum_corner_quality{std::numeric_limits<double>::infinity()};

    for (Index local_node_index = 0; local_node_index < node_count; ++local_node_index)
    {
        const Index previous_local_node_index{(local_node_index + node_count - 1) % node_count};
        const Index next_local_node_index{(local_node_index + 1) % node_count};

        const Node &previous{raw_mesh.nodes[raw_mesh.cell_nodes[cell_node_begin_offset + previous_local_node_index]]};
        const Node &current{raw_mesh.nodes[raw_mesh.cell_nodes[cell_node_begin_offset + local_node_index]]};
        const Node &next{raw_mesh.nodes[raw_mesh.cell_nodes[cell_node_begin_offset + next_local_node_index]]};

        const double previous_dx{previous.x - current.x};
        const double previous_dy{previous.y - current.y};

        const double next_dx{next.x - current.x};
        const double next_dy{next.y - current.y};

        const double corner_cross_magnitude{std::abs(previous_dx * next_dy - previous_dy * next_dx)};
        const double squared_length_sum{previous_dx * previous_dx + previous_dy * previous_dy + next_dx * next_dx +
                                        next_dy * next_dy};

        if (!std::isfinite(corner_cross_magnitude) || !std::isfinite(squared_length_sum) || !(squared_length_sum > 0.0))
        {
            throw_geometry_build_error("cell " + std::to_string(cell_id) + " has invalid quadrilateral edge geometry.");
        }

        // Local corner quality:
        //
        //     q_i = 2 |a x b| / (|a|^2 + |b|^2)
        //
        // q_i = 1 for two perpendicular edges of equal length. Both angular
        // distortion and unequal adjacent edge lengths reduce the metric.
        const double corner_quality{2.0 * corner_cross_magnitude / squared_length_sum};

        // A cell is only as good as its most distorted corner.
        minimum_corner_quality = std::min(minimum_corner_quality, corner_quality);
    }

    return minimum_corner_quality;
}

void ensure_valid_quadrilateral_shape(const RawMeshData &raw_mesh, const Index cell_id)
{
    constexpr Index node_count{4};

    const Index cell_node_begin_offset{raw_mesh.cell_node_offsets[cell_id]};

    bool is_orientation_initialized{};
    bool is_reference_orientation_positive{};

    // For an ordered convex quadrilateral, every consecutive edge pair must
    // produce a non-zero cross product with the same sign. A sign change
    // indicates a concave or self-intersecting cell.
    for (Index local_node_index = 0; local_node_index < node_count; ++local_node_index)
    {
        const Index previous_local_node_index{(local_node_index + node_count - 1) % node_count};
        const Index next_local_node_index{(local_node_index + 1) % node_count};

        const Node &previous{raw_mesh.nodes[raw_mesh.cell_nodes[cell_node_begin_offset + previous_local_node_index]]};
        const Node &current{raw_mesh.nodes[raw_mesh.cell_nodes[cell_node_begin_offset + local_node_index]]};
        const Node &next{raw_mesh.nodes[raw_mesh.cell_nodes[cell_node_begin_offset + next_local_node_index]]};

        const double incoming_x{current.x - previous.x};
        const double incoming_y{current.y - previous.y};

        const double outgoing_x{next.x - current.x};
        const double outgoing_y{next.y - current.y};

        const double signed_corner_cross{incoming_x * outgoing_y - incoming_y * outgoing_x};

        if (!std::isfinite(signed_corner_cross) || signed_corner_cross == 0.0)
        {
            throw_geometry_build_error("cell " + std::to_string(cell_id) + " has a degenerate quadrilateral corner.");
        }

        const bool is_current_orientation_positive{signed_corner_cross > 0.0};

        if (!is_orientation_initialized)
        {
            is_reference_orientation_positive = is_current_orientation_positive;
            is_orientation_initialized = true;
            continue;
        }

        if (is_current_orientation_positive != is_reference_orientation_positive)
        {
            throw_geometry_build_error("cell " + std::to_string(cell_id) +
                                       " has a non-convex or self-intersecting quadrilateral.");
        }
    }
}

void build_cell_geometry(const RawMeshData &raw_mesh, GeometryBuildData &geometry)
{
    const Index cell_count{raw_mesh.cell_types.size()};

    geometry.cell_areas.resize(cell_count);
    geometry.cell_centers.resize(cell_count);
    geometry.cell_qualities.assign(cell_count, std::numeric_limits<double>::quiet_NaN());

    for (Index cell_id = 0; cell_id < cell_count; ++cell_id)
    {
        const Index cell_node_begin_offset{raw_mesh.cell_node_offsets[cell_id]};
        const Index cell_node_end_offset{raw_mesh.cell_node_offsets[cell_id + 1]};
        const Index node_count{cell_node_end_offset - cell_node_begin_offset};

        // Use the first cell node as a local origin.
        //
        // The Shoelace formula is translation invariant, but evaluating it
        // with local coordinates reduces cancellation when small cells are
        // located far from the global origin.
        const Node &reference_node{raw_mesh.nodes[raw_mesh.cell_nodes[cell_node_begin_offset]]};

        double twice_signed_area{};
        double centroid_x_numerator{};
        double centroid_y_numerator{};

        // Accumulate signed polygon area and centroid numerators in one
        // traversal. The sign preserves the polygon orientation; the physical
        // area stored below is its absolute value.
        for (Index local_node_index = 0; local_node_index < node_count; ++local_node_index)
        {
            const Index current_node_id{raw_mesh.cell_nodes[cell_node_begin_offset + local_node_index]};
            const Index next_node_id{raw_mesh.cell_nodes[cell_node_begin_offset + (local_node_index + 1) % node_count]};

            const Node &current{raw_mesh.nodes[current_node_id]};
            const Node &next{raw_mesh.nodes[next_node_id]};

            const double current_x{current.x - reference_node.x};
            const double current_y{current.y - reference_node.y};

            const double next_x{next.x - reference_node.x};
            const double next_y{next.y - reference_node.y};

            const double shoelace_cross{current_x * next_y - next_x * current_y};

            twice_signed_area += shoelace_cross;

            centroid_x_numerator += (current_x + next_x) * shoelace_cross;
            centroid_y_numerator += (current_y + next_y) * shoelace_cross;
        }

        if (twice_signed_area == 0.0)
        {
            throw_geometry_build_error("cell " + std::to_string(cell_id) + " has zero area.");
        }

        const double area{0.5 * std::abs(twice_signed_area)};

        // For a polygon with signed double area S, the centroid denominator is
        // 3S. Add the local origin back after evaluating the translated formula.
        const double centroid_x{reference_node.x + centroid_x_numerator / (3.0 * twice_signed_area)};
        const double centroid_y{reference_node.y + centroid_y_numerator / (3.0 * twice_signed_area)};

        geometry.cell_areas[cell_id] = area;
        geometry.cell_centers[cell_id] = {centroid_x, centroid_y};

        switch (raw_mesh.cell_types[cell_id])
        {
        case CellType::Triangle:
            geometry.cell_qualities[cell_id] = compute_triangle_quality(raw_mesh, cell_id, area);
            break;

        case CellType::Quadrilateral:
            ensure_valid_quadrilateral_shape(raw_mesh, cell_id);
            geometry.cell_qualities[cell_id] = compute_quadrilateral_quality(raw_mesh, cell_id);
            break;
        }
    }
}

void build_face_geometry(const RawMeshData &raw_mesh, const TopologyBuildData &topology, GeometryBuildData &geometry)
{
    const Index face_count{topology.faces.size()};

    geometry.face_centers.resize(face_count);
    geometry.face_lengths.resize(face_count);
    geometry.face_area_vectors.resize(face_count);

    for (Index face_id = 0; face_id < face_count; ++face_id)
    {
        const Face &face{topology.faces[face_id]};

        const Node &node_0{raw_mesh.nodes[face.node_ids[0]]};
        const Node &node_1{raw_mesh.nodes[face.node_ids[1]]};

        const double dx{node_1.x - node_0.x};
        const double dy{node_1.y - node_0.y};

        const double face_length{std::hypot(dx, dy)};

        if (face_length == 0.0)
        {
            throw_geometry_build_error("face " + std::to_string(face_id) + " has zero length.");
        }

        const Vector2 face_center{0.5 * (node_0.x + node_1.x), 0.5 * (node_0.y + node_1.y)};

        // Rotating the edge vector (dx, dy) clockwise gives (dy, -dx), a
        // perpendicular vector whose magnitude is exactly the edge length.
        // This is the two-dimensional finite-volume face-area vector before
        // owner-based orientation is enforced.
        Vector2 face_area_vector{dy, -dx};

        const FaceAdjacency &adjacency{topology.face_adjacencies[face_id]};
        const Vector2 &owner_center{geometry.cell_centers[adjacency.owner]};

        const double owner_to_face_x{face_center.x - owner_center.x};
        const double owner_to_face_y{face_center.y - owner_center.y};

        const double owner_to_face_dot{face_area_vector.x * owner_to_face_x + face_area_vector.y * owner_to_face_y};

        if (owner_to_face_dot == 0.0)
        {
            throw_geometry_build_error("face " + std::to_string(face_id) +
                                       " has ambiguous orientation relative to its owner cell.");
        }

        // Face node ordering is not a geometric orientation convention. Flip
        // the normal when necessary so Sf always points outward from the owner.
        if (owner_to_face_dot < 0.0)
        {
            face_area_vector.x = -face_area_vector.x;
            face_area_vector.y = -face_area_vector.y;
        }

        geometry.face_centers[face_id] = face_center;
        geometry.face_lengths[face_id] = face_length;
        geometry.face_area_vectors[face_id] = face_area_vector;
    }
}

} // namespace

GeometryBuildData build_geometry(const RawMeshData &raw_mesh, const TopologyBuildData &topology)
{
    GeometryBuildData geometry;

    // Cell centers must be available before face-area vectors can be oriented
    // relative to their owner cells.
    build_cell_geometry(raw_mesh, geometry);
    build_face_geometry(raw_mesh, topology, geometry);

    return geometry;
}

} // namespace cfd::detail