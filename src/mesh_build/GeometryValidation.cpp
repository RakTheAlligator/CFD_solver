#include "mesh_build/GeometryValidation.hpp"

#include "cfd/mesh/Face.hpp"
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
void throw_geometry_validation_error(const std::string &message)
{
    throw std::runtime_error("Geometry validation failed: " + message);
}

[[nodiscard]]
bool nearly_equal(const double a, const double b) noexcept
{
    constexpr double tolerance{64.0 * std::numeric_limits<double>::epsilon()};

    // Use a scale-relative comparison because this helper validates positive
    // geometric quantities whose magnitudes may vary with mesh resolution.
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

    if (geometry.cell_qualities.size() != raw_mesh.cell_types.size())
    {
        throw_geometry_validation_error("cell qualities size mismatch.");
    }
}

void validate_cell_geometry(const RawMeshData &raw_mesh, const GeometryBuildData &geometry)
{
    const Index cell_count{raw_mesh.cell_types.size()};

    for (Index cell_id = 0; cell_id < cell_count; ++cell_id)
    {
        const double area{geometry.cell_areas[cell_id]};
        const Vector2 &cell_center{geometry.cell_centers[cell_id]};

        if (!std::isfinite(area))
        {
            throw_geometry_validation_error("cell " + std::to_string(cell_id) + " has non-finite area.");
        }

        if (!std::isfinite(cell_center.x) || !std::isfinite(cell_center.y))
        {
            throw_geometry_validation_error("cell " + std::to_string(cell_id) + " has non-finite center coordinates.");
        }

        if (!(area > 0.0))
        {
            throw_geometry_validation_error("cell " + std::to_string(cell_id) + " has non-positive area.");
        }

        const double quality{geometry.cell_qualities[cell_id]};

        if (!std::isfinite(quality) || !(quality > 0.0))
        {
            throw_geometry_validation_error("cell " + std::to_string(cell_id) + " has invalid cell quality.");
        }

        // The analytical quality metrics are bounded by one. Allow only the
        // small overshoot that can result from floating-point roundoff.
        if (quality > 1.0 && !nearly_equal(quality, 1.0))
        {
            throw_geometry_validation_error("cell " + std::to_string(cell_id) + " has cell quality greater than 1.");
        }
    }
}

void validate_face_geometry(const TopologyBuildData &topology, const GeometryBuildData &geometry)
{
    const Index face_count{topology.faces.size()};

    for (Index face_id = 0; face_id < face_count; ++face_id)
    {
        const double length{geometry.face_lengths[face_id]};
        const Vector2 &face_center{geometry.face_centers[face_id]};
        const Vector2 &face_area_vector{geometry.face_area_vectors[face_id]};

        if (!std::isfinite(length))
        {
            throw_geometry_validation_error("face " + std::to_string(face_id) + " has non-finite length.");
        }

        if (!std::isfinite(face_center.x) || !std::isfinite(face_center.y))
        {
            throw_geometry_validation_error("face " + std::to_string(face_id) + " has non-finite center coordinates.");
        }

        if (!std::isfinite(face_area_vector.x) || !std::isfinite(face_area_vector.y))
        {
            throw_geometry_validation_error("face " + std::to_string(face_id) +
                                            " has non-finite area vector components.");
        }

        if (!(length > 0.0))
        {
            throw_geometry_validation_error("face " + std::to_string(face_id) + " has non-positive length.");
        }

        // In two dimensions, the magnitude of the face-area vector is the
        // physical face length.
        const double area_vector_norm{std::hypot(face_area_vector.x, face_area_vector.y)};

        if (!nearly_equal(area_vector_norm, length))
        {
            throw_geometry_validation_error("face " + std::to_string(face_id) +
                                            " has inconsistent length and area-vector norm.");
        }

        const FaceAdjacency &adjacency{topology.face_adjacencies[face_id]};
        const Vector2 &owner_center{geometry.cell_centers[adjacency.owner]};

        // Sf must point outward from the owner. The owner-to-face vector
        // provides an orientation test independent of the stored face-node
        // ordering.
        const double owner_to_face_dot{face_area_vector.x * (face_center.x - owner_center.x) +
                                       face_area_vector.y * (face_center.y - owner_center.y)};

        if (!(owner_to_face_dot > 0.0))
        {
            throw_geometry_validation_error("face " + std::to_string(face_id) +
                                            " has an area vector not oriented outward from its owner.");
        }

        if (adjacency.neighbor == invalid_index)
        {
            continue;
        }

        const Vector2 &neighbor_center{geometry.cell_centers[adjacency.neighbor]};

        // For an internal face, the owner-oriented Sf must have a positive
        // projection along the owner-to-neighbor direction. It is not generally
        // parallel to that direction on a non-orthogonal mesh.
        const double owner_to_neighbor_dot{face_area_vector.x * (neighbor_center.x - owner_center.x) +
                                           face_area_vector.y * (neighbor_center.y - owner_center.y)};

        if (!(owner_to_neighbor_dot > 0.0))
        {
            throw_geometry_validation_error("internal face " + std::to_string(face_id) +
                                            " has inconsistent owner-neighbor orientation.");
        }
    }
}

void validate_cell_face_closure(const RawMeshData &raw_mesh, const TopologyBuildData &topology,
                                const GeometryBuildData &geometry)
{
    constexpr double tolerance_factor{256.0 * std::numeric_limits<double>::epsilon()};

    const Index cell_count{raw_mesh.cell_types.size()};

    for (Index cell_id = 0; cell_id < cell_count; ++cell_id)
    {
        const Index cell_face_begin_offset{raw_mesh.cell_node_offsets[cell_id]};
        const Index cell_face_end_offset{raw_mesh.cell_node_offsets[cell_id + 1]};

        double area_vector_sum_x{};
        double area_vector_sum_y{};
        double perimeter{};

        // face_area_vectors are stored outward from each face owner. Reverse
        // the sign when the current cell is the neighbor to recover the
        // outward vector relative to this cell.
        for (Index cell_face_position = cell_face_begin_offset; cell_face_position < cell_face_end_offset;
             ++cell_face_position)
        {
            const Index face_id{topology.cell_faces[cell_face_position]};
            const FaceAdjacency &adjacency{topology.face_adjacencies[face_id]};
            const Vector2 &face_area_vector{geometry.face_area_vectors[face_id]};

            perimeter += geometry.face_lengths[face_id];

            if (adjacency.owner == cell_id)
            {
                area_vector_sum_x += face_area_vector.x;
                area_vector_sum_y += face_area_vector.y;
            }
            else if (adjacency.neighbor == cell_id)
            {
                area_vector_sum_x -= face_area_vector.x;
                area_vector_sum_y -= face_area_vector.y;
            }
            else
            {
                throw_geometry_validation_error("cell " + std::to_string(cell_id) +
                                                " references a face to which it does not belong.");
            }
        }

        // A closed control volume satisfies sum(Sf) = 0. Scale the floating-
        // point tolerance with the cell perimeter so both sides retain the
        // dimensions of length and the test adapts to cell size.
        const double closure_norm{std::hypot(area_vector_sum_x, area_vector_sum_y)};
        const double closure_tolerance{tolerance_factor * perimeter};

        if (!std::isfinite(closure_norm) || closure_norm > closure_tolerance)
        {
            throw_geometry_validation_error("cell " + std::to_string(cell_id) +
                                            " does not satisfy face-area-vector closure.");
        }
    }
}

} // namespace

void validate_geometry(const RawMeshData &raw_mesh, const TopologyBuildData &topology,
                       const GeometryBuildData &geometry)
{
    validate_geometry_storage(raw_mesh, topology, geometry);
    validate_cell_geometry(raw_mesh, geometry);
    validate_face_geometry(topology, geometry);
    validate_cell_face_closure(raw_mesh, topology, geometry);
}

} // namespace cfd::detail