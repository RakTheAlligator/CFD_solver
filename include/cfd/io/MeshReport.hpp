#pragma once

#include <iosfwd>

namespace cfd
{

class Mesh;
struct MeshBuildTimings;
struct MeshStatistics;

/// Writes a human-readable summary of mesh topology, geometry, quality, and
/// construction timings to an output stream.
///
/// The output stream remains owned by the caller.
///
/// @param output Destination stream.
/// @param mesh Mesh whose topology is reported.
/// @param statistics Precomputed statistics associated with the mesh.
/// @param timings Timings recorded during mesh construction.
void write_mesh_report(std::ostream &output, const Mesh &mesh, const MeshStatistics &statistics,
                       const MeshBuildTimings &timings);

} // namespace cfd