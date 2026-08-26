#pragma once

#include <iosfwd>

namespace cfd
{

class Mesh;
struct MeshBuildTimings;
struct MeshStatistics;

void write_mesh_report(std::ostream &output, const Mesh &mesh, const MeshStatistics &statistics,
                       const MeshBuildTimings &timings);

} // namespace cfd