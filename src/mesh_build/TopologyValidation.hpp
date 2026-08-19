#pragma once

#include "mesh_build/BuildData.hpp"

namespace cfd {

struct RawMeshData;

namespace detail {

[[nodiscard]]
TopologyStats validate_topology(
    const RawMeshData& raw_mesh,
    const TopologyBuildData& topology);

} // namespace detail

} // namespace cfd