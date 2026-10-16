#pragma once

#include "cfd/mesh/Cell.hpp"
#include "cfd/meshing/BackwardFacingStepGeometry.hpp"
#include "cfd/meshing/GeometryInput.hpp"
#include "cfd/meshing/MeshSizeField.hpp"
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

    /// Mixed TRI/QUAD extraction; actual cells retain their existing CellType.
    /// RPA only: requires cell_type == Quadrilateral, with QUAD blocks and
    /// triangular transition patches. Rectangle/BFS do not support this mode.
    bool hybrid{};
};

/// Explicit wall-normal layer prescription, independent of any y+ model.
/// Physical prescriptions are absent by default. Enabled layers require first_cell_height or thickness.
/// If both are supplied they must describe the same geometric series.
struct WallLayerOptions
{
    /// Absent/false: no physical prescription; true: explicit physical layers
    /// with a height/thickness prescription. Automatic core grading is separate.
    std::optional<bool> enabled{};
    std::optional<double> first_cell_height{};
    Index layer_count{4};
    double growth_rate{1.2};
    std::optional<double> thickness{};

    [[nodiscard]] bool active() const noexcept
    {
        return enabled.value_or(false);
    }
};

/// RPA discretization input. Region geometry is separate from the numerical
/// settings; no legacy regional refinement factors are used.
struct RpaNozzleMeshingOptions
{
    RpaRegionLayoutOptions regions{};
    /// Minimum signed corner sine for TRI/QUAD (scaled corner Jacobian).
    /// Orthogonal QUAD anisotropy does not change its corner angles; thin TRI
    /// cells do have small angles and can fail this threshold.
    double min_cell_quality{0.05};
    WallLayerOptions inner_wall_layers{};
    WallLayerOptions outer_wall_layers{};
    std::optional<Index> target_cell_budget{};
    std::optional<std::filesystem::path> external_size_field_file{};
};

/// Generates only fluid cells in the Cartesian RPA meridian/plume domain.
/// CAD import is automatic; this does not implement axisymmetric physics.
/// @throws std::invalid_argument If domain or meshing options are inconsistent.
/// @throws std::runtime_error If CAD import or exclusive requested topology fails.
[[nodiscard]]
RawMeshData generate_mesh(const RpaNozzleGeometry &geometry, const MeshGenerationOptions &options,
                          const RpaNozzleMeshingOptions &rpa_options = {}, const MeshSizeField &size_field = {});

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
/// RPA requires default automatic_options; its dedicated options replace the
/// former generic automatic-meshing policy.
[[nodiscard]]
RawMeshData generate_mesh(const GeometryInput &geometry, const MeshGenerationOptions &options,
                          const AutomaticMeshingOptions &automatic_options,
                          const BackwardFacingStepMeshingOptions &step_options,
                          const RpaNozzleMeshingOptions &rpa_options = {}, const MeshSizeField &size_field = {});

} // namespace cfd
