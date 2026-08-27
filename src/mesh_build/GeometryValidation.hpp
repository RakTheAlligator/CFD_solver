#pragma once

#include "mesh_build/MeshBuildData.hpp"

namespace cfd
{

struct RawMeshData;

namespace detail
{

void validate_geometry(const RawMeshData &raw_mesh, const TopologyBuildData &topology,
                       const GeometryBuildData &geometry);

} // namespace detail

} // namespace cfd