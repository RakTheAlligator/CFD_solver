#pragma once

#include "cfd/Cell.hpp"
#include "cfd/Geometry.hpp"
#include "cfd/RawMeshData.hpp"

namespace cfd {

struct MeshGenerationOptions {
    double mesh_size{};
    CellType cell_type{CellType::Triangle};
};

[[nodiscard]]
RawMeshData generate_mesh(
    const RectangleGeometry& geometry,
    const MeshGenerationOptions& options
);

} // namespace cfd