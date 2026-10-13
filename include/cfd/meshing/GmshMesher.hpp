#pragma once

#include "cfd/mesh/Cell.hpp"
#include "cfd/meshing/BackwardFacingStepGeometry.hpp"
#include "cfd/meshing/GeometryInput.hpp"
#include "cfd/meshing/RawMeshData.hpp"
#include "cfd/meshing/RectangleGeometry.hpp"

#include <optional>

namespace cfd
{

/// Generic policy for built-in automatic mesh strategies.
struct AutomaticMeshingOptions
{
    /// Enables geometry-specific automatic meshing when available.
    bool enabled{true};

    /// Maximum ratio between consecutive cell sizes on a structured line.
    double maximum_growth_rate{1.2};

    /// Target local wall-adjacent size relative to the base mesh size.
    ///
    /// Geometry, conformity, or growth-rate constraints may require a smaller
    /// actual size.
    double wall_refinement_factor{0.5};

    /// Nominal number of cells in each wall transition.
    Index wall_refinement_layers{4};
};

/// Automatic-meshing policy specific to a backward-facing-step domain.
///
/// The default factors request a moderate twofold refinement relative to the
/// base size. Layer counts scale transition widths with the base size; they are
/// policy defaults rather than physical constants.
struct BackwardFacingStepMeshingOptions
{
    /// Target local step-adjacent size relative to the base mesh size.
    ///
    /// Geometry, conformity, or growth-rate constraints may require a smaller
    /// actual size.
    double step_refinement_factor{0.5};

    /// Nominal number of cells returning from step refinement to the bulk.
    Index step_refinement_layers{6};
};

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

/// RPA-specific size policy. Extents are metres; absent extents use detected
/// throat/exit scales. maximumGrowthRate controls smooth size-field ramps,
/// not an exact cell-to-cell ratio on an unstructured mesh.
struct RpaNozzleMeshingOptions
{
    double convergent_refinement_factor{0.75};
    double throat_refinement_factor{0.25};
    double divergent_refinement_factor{0.75};
    double exit_refinement_factor{0.25};
    double plume_refinement_factor{0.5};
    std::optional<double> throat_refinement_length{};
    std::optional<double> exit_refinement_length{};
    std::optional<double> near_plume_length{};
};

/// Generates only fluid cells in the Cartesian RPA meridian/plume domain.
/// CAD import is automatic; this does not implement axisymmetric physics.
/// @throws std::invalid_argument If domain or meshing options are inconsistent.
/// @throws std::runtime_error If CAD import or exclusive requested topology fails.
[[nodiscard]]
RawMeshData generate_mesh(const RpaNozzleGeometry &geometry, const MeshGenerationOptions &options,
                          const AutomaticMeshingOptions &automatic_options = {},
                          const RpaNozzleMeshingOptions &rpa_options = {});

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

/// Generates a backward-facing-step mesh using Gmsh.
///
/// The logical boundary groups are `inlet` (upstream vertical face), `wall`
/// (the four solid-wall segments), and `outlet` (downstream vertical face).
///
/// @throws std::invalid_argument If the geometry, mesh size, or requested cell
///         type is invalid.
/// @throws std::runtime_error If Gmsh returns an unsupported mesh.
[[nodiscard]]
RawMeshData generate_mesh(const BackwardFacingStepGeometry &geometry, const MeshGenerationOptions &options);

/// Generates a backward-facing-step mesh with an explicit geometry-specific
/// automatic-meshing policy.
[[nodiscard]]
RawMeshData generate_mesh(const BackwardFacingStepGeometry &geometry, const MeshGenerationOptions &options,
                          const AutomaticMeshingOptions &automatic_options,
                          const BackwardFacingStepMeshingOptions &step_options);

/// Dispatches mesh generation for the selected closed geometry alternative.
[[nodiscard]]
RawMeshData generate_mesh(const GeometryInput &geometry, const MeshGenerationOptions &options);

/// Dispatches mesh generation with explicit geometry-specific policies.
[[nodiscard]]
RawMeshData generate_mesh(const GeometryInput &geometry, const MeshGenerationOptions &options,
                          const AutomaticMeshingOptions &automatic_options,
                          const BackwardFacingStepMeshingOptions &step_options,
                          const RpaNozzleMeshingOptions &rpa_options = {});

} // namespace cfd
