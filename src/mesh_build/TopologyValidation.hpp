#pragma once

#include "mesh_build/MeshBuildData.hpp"

namespace cfd
{

struct RawMeshData;

namespace detail
{

void validate_topology(const RawMeshData &raw_mesh, const TopologyBuildData &topology);

} // namespace detail

} // namespace cfd