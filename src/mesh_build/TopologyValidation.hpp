#pragma once

#include "mesh_build/MeshBuildData.hpp"

namespace cfd
{

struct RawMeshData;

namespace detail
{

// Validates the topology produced from RawMeshData before it is transferred
// into the final Mesh.
//
// Validation checks connectivity sizes and references, owner/neighbor
// consistency, boundary assignments, and the one-to-one relationship between
// local cell edges and constructed mesh faces.
void validate_topology(const RawMeshData &raw_mesh, const TopologyBuildData &topology);

} // namespace detail

} // namespace cfd