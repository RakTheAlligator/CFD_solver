#pragma once

#include "mesh_build/MeshBuildData.hpp"

namespace cfd
{

struct RawMeshData;

namespace detail
{

// Builds temporary geometric data from validated raw mesh and topology.
//
// The returned arrays are indexed directly by internal cell and face IDs.
// Geometric consistency is checked separately by validate_geometry() before
// the data is transferred into the final Mesh.
[[nodiscard]]
GeometryBuildData build_geometry(const RawMeshData &raw_mesh, const TopologyBuildData &topology);

} // namespace detail

} // namespace cfd