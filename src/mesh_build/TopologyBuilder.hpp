#pragma once

#include "mesh_build/MeshBuildData.hpp"

namespace cfd
{

struct RawMeshData;

namespace detail
{

// Builds temporary face topology from structurally validated raw mesh data.
//
// Faces are deduplicated independently of local cell orientation. The first
// cell encountering a face becomes its owner; a second adjacent cell becomes
// its neighbor. Physical boundary IDs are then attached from RawMeshData.
[[nodiscard]]
TopologyBuildData build_topology(const RawMeshData &raw_mesh);

} // namespace detail

} // namespace cfd