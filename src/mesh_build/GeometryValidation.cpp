#include "mesh_build/GeometryValidation.hpp"

#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Node.hpp"
#include "cfd/mesh/Types.hpp"

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

struct ScalarAccumulator
{
    double minimum{std::numeric_limits<double>::infinity()};
    double maximum{-std::numeric_limits<double>::infinity()};
    double sum{};
};

void update_stats(ScalarAccumulator &stats, const double value)
{
    stats.minimum = std::min(stats.minimum, value);
    stats.maximum = std::max(stats.maximum, value);
    stats.sum += value;
}

[[noreturn]]
void throw_geometry_validation_error(const std::string &message)
{
    throw std::runtime_error("Geometry validation failed: " + message);
}

[[nodiscard]]
bool nearly_equal(const double a, const double b) noexcept
{
    constexpr double tolerance{64.0 * std::numeric_limits<double>::epsilon()};

    const double scale{std::max(std::abs(a), std::abs(b))};

    return std::abs(a - b) <= tolerance * scale;
}

void validate_geometry_storage(const RawMeshData &raw_mesh, const TopologyBuildData &topology,
                               const GeometryBuildData &geometry)
{
    if (geometry.cell_areas.size() != raw_mesh.cell_types.size())
    {
        throw_geometry_validation_error("cell areas size mismatch.");
    }

    if (geometry.cell_centers.size() != raw_mesh.cell_types.size())
    {
        throw_geometry_validation_error("cell centers size mismatch.");
    }

    if (geometry.face_centers.size() != topology.faces.size())
    {
        throw_geometry_validation_error("face centers size mismatch.");
    }

    if (geometry.face_lengths.size() != topology.faces.size())
    {
        throw_geometry_validation_error("face lengths size mismatch.");
    }

    if (geometry.face_area_vectors.size() != topology.faces.size())
    {
        throw_geometry_validation_error("face area vectors size mismatch.");
    }
}

void validate_cell_geometry(const RawMeshData &raw_mesh, const GeometryBuildData &geometry, GeometryStats &stats)
{
    ScalarAccumulator area_stats;
    ScalarAccumulator size_stats;
    ScalarAccumulator triangle_quality_stats;

    const Index cell_count{raw_mesh.cell_types.size()};

    double total_area{};
    Index triangle_count{};
    Index worst_quality_cell{invalid_index};
    double worst_quality{std::numeric_limits<double>::infinity()};

    for (Index cell_id = 0; cell_id < cell_count; ++cell_id)
    {
        const double area{geometry.cell_areas[cell_id]};
        const Vector2 &center{geometry.cell_centers[cell_id]};

        if (!std::isfinite(area))
        {
            throw_geometry_validation_error("cell " + std::to_string(cell_id) + " has non-finite area.");
        }

        if (!std::isfinite(center.x) || !std::isfinite(center.y))
        {
            throw_geometry_validation_error("cell " + std::to_string(cell_id) + " has non-finite center coordinates.");
        }

        if (area <= 0.0)
        {
            throw_geometry_validation_error("cell " + std::to_string(cell_id) + " has non-positive area.");
        }

        const double cell_size{std::sqrt(area)};

        update_stats(area_stats, area);
        update_stats(size_stats, cell_size);

        total_area += area;

        if (raw_mesh.cell_types[cell_id] == CellType::Triangle)
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

            const double squared_length_sum{dx_01 * dx_01 + dy_01 * dy_01 + dx_12 * dx_12 + dy_12 * dy_12 +
                                            dx_20 * dx_20 + dy_20 * dy_20};

            if (!std::isfinite(squared_length_sum) || !(squared_length_sum > 0.0))
            {
                throw_geometry_validation_error("cell " + std::to_string(cell_id) + " has invalid edge lengths.");
            }

            const double quality{4.0 * std::sqrt(3.0) * area / squared_length_sum};

            if (!std::isfinite(quality) || !(quality > 0.0))
            {
                throw_geometry_validation_error("cell " + std::to_string(cell_id) + " has invalid triangle quality.");
            }

            ++triangle_count;

            update_stats(triangle_quality_stats, quality);

            if (quality < worst_quality)
            {
                worst_quality = quality;
                worst_quality_cell = cell_id;
            }
        }
    }

    stats.total_cell_area = total_area;

    stats.cell_areas = {.minimum = area_stats.minimum,
                        .maximum = area_stats.maximum,
                        .mean = area_stats.sum / static_cast<double>(cell_count)};

    stats.cell_sizes = {.minimum = size_stats.minimum,
                        .maximum = size_stats.maximum,
                        .mean = size_stats.sum / static_cast<double>(cell_count)};

    if (triangle_count > 0)
    {
        stats.triangle_quality = {.minimum = triangle_quality_stats.minimum,
                                  .maximum = triangle_quality_stats.maximum,
                                  .mean = triangle_quality_stats.sum / static_cast<double>(triangle_count)};

        stats.worst_quality_cell = worst_quality_cell;
    }
}

void validate_face_geometry(const TopologyBuildData &topology, const GeometryBuildData &geometry, GeometryStats &stats)
{
    ScalarAccumulator length_stats;

    const Index face_count{topology.faces.size()};

    for (Index face_id = 0; face_id < face_count; ++face_id)
    {
        const double length{geometry.face_lengths[face_id]};
        const Vector2 &center{geometry.face_centers[face_id]};
        const Vector2 &area_vector{geometry.face_area_vectors[face_id]};

        if (!std::isfinite(length))
        {
            throw_geometry_validation_error("face " + std::to_string(face_id) + " has non-finite length.");
        }

        if (!std::isfinite(center.x) || !std::isfinite(center.y))
        {
            throw_geometry_validation_error("face " + std::to_string(face_id) + " has non-finite center coordinates.");
        }

        if (!std::isfinite(area_vector.x) || !std::isfinite(area_vector.y))
        {
            throw_geometry_validation_error("face " + std::to_string(face_id) +
                                            " has non-finite area vector components.");
        }

        if (length <= 0.0)
        {
            throw_geometry_validation_error("face " + std::to_string(face_id) + " has non-positive length.");
        }

        const double area_vector_norm{std::hypot(area_vector.x, area_vector.y)};

        if (!nearly_equal(area_vector_norm, length))
        {
            throw_geometry_validation_error("face " + std::to_string(face_id) +
                                            " has inconsistent length and area-vector norm.");
        }

        const FaceAdjacency &adjacency{topology.face_adjacencies[face_id]};
        const Vector2 &owner_center{geometry.cell_centers[adjacency.owner]};

        const double owner_orientation{area_vector.x * (center.x - owner_center.x) +
                                       area_vector.y * (center.y - owner_center.y)};

        if (!(owner_orientation > 0.0))
        {
            throw_geometry_validation_error("face " + std::to_string(face_id) +
                                            " has an area vector not oriented outward from its owner.");
        }

        if (adjacency.neighbor != invalid_index)
        {
            const Vector2 &neighbor_center{geometry.cell_centers[adjacency.neighbor]};

            const double owner_to_neighbor_orientation{area_vector.x * (neighbor_center.x - owner_center.x) +
                                                       area_vector.y * (neighbor_center.y - owner_center.y)};

            if (!(owner_to_neighbor_orientation > 0.0))
            {
                throw_geometry_validation_error("internal face " + std::to_string(face_id) +
                                                " has inconsistent owner-neighbor orientation.");
            }
        }

        update_stats(length_stats, length);
    }

    stats.face_lengths = {.minimum = length_stats.minimum,
                          .maximum = length_stats.maximum,
                          .mean = length_stats.sum / static_cast<double>(face_count)};
}

} // namespace

GeometryStats validate_geometry(const RawMeshData &raw_mesh, const TopologyBuildData &topology,
                                const GeometryBuildData &geometry)
{
    validate_geometry_storage(raw_mesh, topology, geometry);

    GeometryStats stats;

    validate_cell_geometry(raw_mesh, geometry, stats);
    validate_face_geometry(topology, geometry, stats);

    return stats;
}

} // namespace cfd::detail