#pragma once

#include "cfd/meshing/GmshMesher.hpp"
#include "cfd/meshing/RpaRegionLayout.hpp"

#include <span>
#include <vector>

namespace cfd::detail
{
/// Column-major structured grid; generation-only data, not a Mesh cache.
struct RpaStructuredGrid
{
    std::vector<double> x;
    Index radial_cells{};
    std::vector<Point2> nodes;

    [[nodiscard]] const Point2 &at(Index column, Index row) const
    {
        return nodes[column * (radial_cells + 1) + row];
    }
};

struct RpaMeshDiscretization
{
    RpaStructuredGrid engine, band, plume, farfield;
    std::vector<double> lip_y;
    Index transition_upstream_cells{};
    double transition_size{};
};

/// Angular quality: minimum signed corner sine for TRI and QUAD.
/// Every corner must be positive. Orthogonal QUAD anisotropy is not penalized;
/// a thin TRI still has small angles, even when its signed area is positive.
[[nodiscard]] double rpa_angular_quality(std::span<const Point2> vertices);

/// Maximum vertical deviation of a linear wall face from the imported profile.
/// Includes its endpoints and every profile vertex inside the face interval.
/// Diagnostic only: this error does not refine the mean-height columns.
[[nodiscard]] double rpa_wall_chord_error(std::span<const Point2> profile, Point2 a, Point2 b);

/// Pure geometry/numerics: no Gmsh calls. Fixed minimum row counts are 8/4/2
/// for ENGINE/BAND/FARFIELD. CAD junctions are mandatory column positions.
[[nodiscard]] RpaMeshDiscretization discretize_rpa(const RpaNozzleDomain &domain, const RpaRegionLayout &layout,
                                                   const MeshGenerationOptions &generation,
                                                   const RpaNozzleMeshingOptions &options);

/// Checks the actual generated cells, not only the candidate structured grids.
void require_rpa_mesh_quality(const RawMeshData &mesh, const RpaNozzleMeshingOptions &options);
} // namespace cfd::detail
