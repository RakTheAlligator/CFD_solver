#pragma once

#include "cfd/mesh/Mesh.hpp"

#include <chrono>

namespace cfd
{

struct RawMeshData;

/// Timings of the main mesh-construction stages.
///
/// The durations measure topology construction and geometry construction,
/// including the validation associated with each stage.
struct MeshBuildTimings
{
    std::chrono::duration<double, std::milli> topology{};
    std::chrono::duration<double, std::milli> geometry{};
};

/// Result of validated mesh construction.
struct MeshBuildResult
{
    Mesh mesh;
    MeshBuildTimings timings;
};

/// Builds the validated internal Mesh representation from raw mesh data.
///
/// Raw mesh storage is consumed so that large data arrays can be transferred
/// into the final Mesh without unnecessary copies.
///
/// The preprocessing sequence is:
/// raw-mesh validation, topology construction and validation, then geometry
/// construction and validation.
///
/// @param raw_mesh Raw mesh representation whose storage may be moved into the
///        resulting Mesh. The object must not be relied upon after this call.
/// @return The constructed Mesh together with topology and geometry timings.
/// @throws std::runtime_error If raw-mesh, topology, or geometry validation
///         fails.
[[nodiscard]]
MeshBuildResult build_mesh(RawMeshData &&raw_mesh);

} // namespace cfd