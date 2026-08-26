#pragma once

#include "cfd/mesh/Mesh.hpp"

#include <chrono>

namespace cfd
{

struct RawMeshData;

struct MeshBuildTimings
{
    std::chrono::duration<double, std::milli> topology{};
    std::chrono::duration<double, std::milli> geometry{};
};

struct MeshBuildResult
{
    Mesh mesh;
    MeshBuildTimings timings;
};

[[nodiscard]]
MeshBuildResult build_mesh(RawMeshData &&raw_mesh);

} // namespace cfd