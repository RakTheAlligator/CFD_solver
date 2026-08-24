#pragma once

#include "cfd/mesh/Cell.hpp"
#include "cfd/meshing/DomainGeometry.hpp"
#include "cfd/meshing/RawMeshData.hpp"

namespace cfd
{

struct MeshGenerationOptions
{
    // Target characteristic mesh length in metres.
    double mesh_size{};

    CellType cell_type{CellType::Triangle};
};
[[nodiscard]]
RawMeshData generate_mesh(const RectangleGeometry &geometry, const MeshGenerationOptions &options);

} // namespace cfd