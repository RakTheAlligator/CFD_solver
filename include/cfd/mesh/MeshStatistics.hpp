#pragma once

#include "cfd/mesh/Types.hpp"
#include <array>

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
/// Geometric quantities use SI units. Cell quality and neighboring-cell size
/// ratios are dimensionless; face non-orthogonality is reported in degrees.
struct MeshStatistics
{
    Index internal_face_count{};
    Index boundary_face_count{};

    double total_cell_area{};

    ScalarStatistics cell_areas;
    ScalarStatistics cell_sizes;
    ScalarStatistics face_lengths;
    ScalarStatistics cell_quality;
    ScalarStatistics internal_face_non_orthogonality_degrees;
    ScalarStatistics internal_face_neighbor_cell_size_ratios;

    /// Internal ID of the cell with the minimum quality.
    ///
    /// `invalid_index` indicates that no quality value was available.
    Index worst_quality_cell_id{invalid_index};

    /// Internal face with the maximum non-orthogonality angle.
    ///
    /// `invalid_index` indicates that the mesh has no internal face.
    Index maximum_non_orthogonality_face_id{invalid_index};

    /// Internal face with the maximum neighboring-cell size ratio.
    ///
    /// `invalid_index` indicates that the mesh has no internal face.
    Index maximum_neighbor_cell_size_ratio_face_id{invalid_index};
    Index triangle_cell_count{};
    Index quadrilateral_cell_count{};
    /// Lmax^2/(2*area) for TRI, Lmax^2/area for QUAD, where Lmax is
    /// the longest edge. Captures thin triangles as well as stretched quads;
    /// informational, not a validity test for anisotropic wall cells.
    ScalarStatistics cell_aspect_ratios;
    /// Tangential face-centre offset from the centre-line intersection,
    /// normalized by face length (boundary: owner-normal projection).
    ScalarStatistics finite_volume_face_skewness;
    /// Quality bins [0,.1), ... [.9,1]; no cell-sized diagnostic buffer.
    std::array<Index, 10> quality_histogram{};
};

/// Computes descriptive statistics from mesh topology and geometry.
///
/// @param mesh Validated mesh to inspect.
/// @return Statistics derived from the mesh without modifying it.
[[nodiscard]]
MeshStatistics compute_mesh_statistics(const Mesh &mesh);

} // namespace cfd
