#pragma once

#include "cfd/mesh/Cell.hpp"
#include "cfd/meshing/RawMeshData.hpp"
#include "cfd/meshing/RectangleGeometry.hpp"

namespace cfd
{

/// Options controlling two-dimensional mesh generation.
struct MeshGenerationOptions
{
    /// Target characteristic mesh length in metres.
    ///
    /// This value controls the requested spatial resolution; generated element
    /// sizes are not required to match it exactly.
    double mesh_size{};

    /// Requested topology of the generated two-dimensional cells.
    CellType cell_type{CellType::Triangle};
};

/// Generates a rectangular two-dimensional mesh using Gmsh.
///
/// Gmsh node and physical-group tags are converted to the solver's internal
/// zero-based IDs before the data is returned. The generated rectangle uses
/// the logical boundary groups `inlet` (left), `wall` (top and bottom), and
/// `outlet` (right).
///
/// @param geometry Dimensions of the rectangular domain.
/// @param options Meshing resolution and requested cell topology.
/// @return Raw mesh data ready for preprocessing and validation.
/// @throws std::invalid_argument If the geometry, mesh size, or requested cell
///         type is invalid.
/// @throws std::runtime_error If the generated Gmsh data is inconsistent with
///         the representation expected by the preprocessing pipeline.
[[nodiscard]]
RawMeshData generate_mesh(const RectangleGeometry &geometry, const MeshGenerationOptions &options);

} // namespace cfd