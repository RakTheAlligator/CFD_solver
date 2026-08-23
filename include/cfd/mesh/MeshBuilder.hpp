#pragma once

#include "cfd/mesh/Mesh.hpp"

namespace cfd
{

struct RawMeshData;

[[nodiscard]]
Mesh build_mesh(RawMeshData &&raw_mesh);

} // namespace cfd