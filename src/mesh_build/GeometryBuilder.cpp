#include "mesh_build/GeometryBuilder.hpp"

#include "cfd/mesh/Node.hpp"
#include "cfd/mesh/Types.hpp"
#include "cfd/mesh/Cell.hpp"

#include "cfd/meshing/RawMeshData.hpp"

#include <cmath>
#include <stdexcept>
#include <string>
#include <limits>
namespace cfd::detail
{

namespace
{
[[noreturn]]
void throw_geometry_build_error(const std::string &message)
{
    throw std::runtime_error("Geometry construction failed: " + message);
} // end of throw_geometry_build_error

[[nodiscard]]
double compute_triangle_quality(const RawMeshData &raw_mesh, const Index cell_id, const double area)
{
    const Index begin{raw_mesh.cell_node_offsets[cell_id]};

    const Node &node_0{raw_mesh.nodes[raw_mesh.cell_nodes[begin]]};
    const Node &node_1{raw_mesh.nodes[raw_mesh.cell_nodes[begin + 1]]};
    const Node &node_2{raw_mesh.nodes[raw_mesh.cell_nodes[begin + 2]]};

    const double dx_01{node_1.x - node_0.x};
    const double dy_01{node_1.y - node_0.y};

    const double dx_12{node_2.x - node_1.x};
    const double dy_12{node_2.y - node_1.y};

    const double dx_20{node_0.x - node_2.x};
    const double dy_20{node_0.y - node_2.y};

    const double squared_length_sum{
        dx_01 * dx_01 + dy_01 * dy_01 +
        dx_12 * dx_12 + dy_12 * dy_12 +
        dx_20 * dx_20 + dy_20 * dy_20};

    if (!std::isfinite(squared_length_sum) || !(squared_length_sum > 0.0))
    {
        throw_geometry_build_error(
            "cell " + std::to_string(cell_id) + " has invalid edge lengths.");
    }

    return 4.0 * std::sqrt(3.0) * area / squared_length_sum;
}

void build_cell_geometry(const RawMeshData &raw_mesh, GeometryBuildData &geometry)
{
    const Index cell_count{raw_mesh.cell_types.size()};

    geometry.cell_areas.resize(cell_count);
    geometry.cell_centers.resize(cell_count);

    geometry.cell_qualities.assign(
        cell_count,
        std::numeric_limits<double>::quiet_NaN());

    for (Index cell_id = 0; cell_id < cell_count; ++cell_id)
    {
        const Index begin{raw_mesh.cell_node_offsets[cell_id]};
        const Index end{raw_mesh.cell_node_offsets[cell_id + 1]};
        const Index node_count{end - begin};

        double twice_signed_area{};
        double centroid_x_numerator{};
        double centroid_y_numerator{};

        for (Index i = 0; i < node_count; ++i)
        {
            const Index current_node_id{raw_mesh.cell_nodes[begin + i]};
            const Index next_node_id{raw_mesh.cell_nodes[begin + (i + 1) % node_count]};

            const Node &current{raw_mesh.nodes[current_node_id]};
            const Node &next{raw_mesh.nodes[next_node_id]};

            const double cross{current.x * next.y - next.x * current.y};

            twice_signed_area += cross;

            centroid_x_numerator += (current.x + next.x) * cross;
            centroid_y_numerator += (current.y + next.y) * cross;
        }
        if (twice_signed_area == 0.0)
        {
            throw_geometry_build_error("cell " + std::to_string(cell_id) + " has zero area.");
        }
        geometry.cell_areas[cell_id] = 0.5 * std::abs(twice_signed_area);
        geometry.cell_centers[cell_id] = {centroid_x_numerator / (3.0 * twice_signed_area),
                                          centroid_y_numerator / (3.0 * twice_signed_area)};
        if (raw_mesh.cell_types[cell_id] == CellType::Triangle)
        {
            geometry.cell_qualities[cell_id] =
                compute_triangle_quality(raw_mesh, cell_id, geometry.cell_areas[cell_id]);
        }
    }
} // end of build_cell_geometry

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

        Vector2 area_vector{dy, -dx};

        const FaceAdjacency &adjacency{topology.face_adjacencies[face_id]};
        const Vector2 &owner_center{geometry.cell_centers[adjacency.owner]};

        const double owner_to_face_x{face_center.x - owner_center.x};
        const double owner_to_face_y{face_center.y - owner_center.y};

        const double orientation{area_vector.x * owner_to_face_x + area_vector.y * owner_to_face_y};

        if (orientation == 0.0)
        {
            throw_geometry_build_error("face " + std::to_string(face_id) +
                                       " has ambiguous orientation relative to its owner cell.");
        }
        if (orientation < 0.0)
        {
            area_vector.x = -area_vector.x;
            area_vector.y = -area_vector.y;
        }
        geometry.face_centers[face_id] = face_center;
        geometry.face_lengths[face_id] = face_length;
        geometry.face_area_vectors[face_id] = area_vector;
    }
} // end of build_face_geometry
} // namespace

GeometryBuildData build_geometry(const RawMeshData &raw_mesh, const TopologyBuildData &topology)
{
    GeometryBuildData geometry;

    build_cell_geometry(raw_mesh, geometry);

    build_face_geometry(raw_mesh, topology, geometry);

    return geometry;
} // end of build_geometry
} // namespace cfd::detail
