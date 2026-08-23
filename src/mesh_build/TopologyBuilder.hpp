#pragma once

#include "mesh_build/BuildData.hpp"

namespace cfd
{

struct RawMeshData;

namespace detail
{

[[nodiscard]]
TopologyBuildData build_topology(const RawMeshData &raw_mesh);

} // namespace detail

} // namespace cfd