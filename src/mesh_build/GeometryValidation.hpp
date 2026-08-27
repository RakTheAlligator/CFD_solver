#pragma once

#include "mesh_build/MeshBuildData.hpp"

namespace cfd
{

struct RawMeshData;

namespace detail
{

// Validates geometric data before it is transferred into the final Mesh.
//
// Validation checks storage consistency, finite and positive geometric
// quantities, face-area-vector orientation, and finite-volume closure for each
// cell.
void validate_geometry(const RawMeshData &raw_mesh, const TopologyBuildData &topology,
                       const GeometryBuildData &geometry);

} // namespace detail

} // namespace cfd