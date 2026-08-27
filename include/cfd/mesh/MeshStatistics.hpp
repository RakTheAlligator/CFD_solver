#pragma once

#include "cfd/mesh/Types.hpp"

namespace cfd
{

class Mesh;

/// Minimum, maximum, and arithmetic mean of a scalar mesh quantity.
struct ScalarStatistics
{
    double minimum{};
    double maximum{};
    double mean{};
};

/// Aggregate statistics derived from a validated Mesh.
///
/// Geometric quantities use SI units. Cell quality is dimensionless.
struct MeshStatistics
{
    Index internal_face_count{};
    Index boundary_face_count{};

    double total_cell_area{};

    ScalarStatistics cell_areas;
    ScalarStatistics cell_sizes;
    ScalarStatistics face_lengths;
    ScalarStatistics cell_quality;

    /// Internal ID of the cell with the minimum quality.
    ///
    /// `invalid_index` indicates that no quality value was available.
    Index worst_quality_cell_id{invalid_index};
};

/// Computes descriptive statistics from mesh topology and geometry.
///
/// @param mesh Validated mesh to inspect.
/// @return Statistics derived from the mesh without modifying it.
[[nodiscard]]
MeshStatistics compute_mesh_statistics(const Mesh &mesh);

} // namespace cfd