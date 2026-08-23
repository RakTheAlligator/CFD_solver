#pragma once

#include "mesh_build/BuildData.hpp"

namespace cfd
{

struct RawMeshData;

namespace detail
{

[[nodiscard]]
GeometryBuildData build_geometry(const RawMeshData &raw_mesh, const TopologyBuildData &topology);

} // namespace detail

} // namespace cfd