#pragma once

#include "cfd/meshing/GmshMesher.hpp"

#include <iosfwd>
#include <span>

namespace cfd
{
class Mesh;
/// Reports mutually exclusive geometry-based regions, not CAD entity IDs.
/// Optional output receives logical region IDs 1..5 for VTU coloring without
/// copying Mesh arrays: engine interior, outer band, plume, transition, upper farfield.
/// Also reports ENGINE/BAND wall-face chord errors and any import-tolerance
/// exceedance, independently of topological/angular validity. No refinement
/// or geometric-accuracy certification is performed by this report.
/// @throws std::runtime_error If a cell centre lies outside the fluid domain
///         (including inside the fictitious wall).
void write_rpa_mesh_diagnostics(std::ostream &output, const Mesh &mesh, const RpaNozzleDomain &domain,
                                const MeshGenerationOptions &generation, const RpaNozzleMeshingOptions &options,
                                std::span<double> region_ids = {});
} // namespace cfd
