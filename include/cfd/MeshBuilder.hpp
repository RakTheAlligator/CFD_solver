#pragma once

#include "cfd/Mesh.hpp"

namespace cfd {

struct RawMeshData;

[[nodiscard]]
Mesh build_mesh(RawMeshData&& raw_mesh);

} // namespace cfd