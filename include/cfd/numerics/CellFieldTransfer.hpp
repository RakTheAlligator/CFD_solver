#pragma once

#include "cfd/math/Point2.hpp"
#include "cfd/math/Vector2.hpp"
#include "cfd/mesh/Types.hpp"

#include <vector>

namespace cfd
{

class CellScalarField;
class CellVectorField;
class Mesh;
class ScalarBoundaryConditions;

/// Reusable cell-centered initial-guess transfer between independent 2D meshes.
///
/// Construction locates each target center in a source triangle or convex quad
/// using a temporary bounding-box tree. Inclusion allows a distance tolerance
/// of `64 * epsilon * scale`, where scale is the maximum absolute source-cell
/// vertex coordinate or bounding-box extent (no unit-sized floor).
/// On shared edges/vertices, the smallest admissible source cell ID wins.
/// No nearest-cell fallback or integral-conservation guarantee is provided.
///
/// @note Both meshes are non-owning dependencies: they must outlive this object
///       and must not be moved from or replaced. Fields do not carry Mesh
///       identity; callers must supply values in the associated cell ordering.
/// @note The mapping owns one source ID and displacement per target cell.
///       Valid applications allocate nothing and copy no field-sized arrays.
///       Velocity components use the same scalar API; face fluxes are not transferred.
class CellFieldTransfer
{
  public:
    /// Builds a fixed geometric mapping; the spatial search storage is then released.
    /// @throws std::invalid_argument If either mesh has no cells.
    /// @throws std::runtime_error If a target center is not covered by the source mesh.
    CellFieldTransfer(const Mesh &source_mesh, const Mesh &target_mesh);

    CellFieldTransfer(const CellFieldTransfer &) = delete;
    CellFieldTransfer &operator=(const CellFieldTransfer &) = delete;
    CellFieldTransfer(CellFieldTransfer &&) noexcept = default;
    CellFieldTransfer &operator=(CellFieldTransfer &&) noexcept = delete;
    ~CellFieldTransfer() = default;

    /// Overwrites each target value with its containing source-cell value.
    ///
    /// Used source values are validated before any target value is written.
    /// @throws std::invalid_argument If mesh storage/cardinalities changed,
    ///         field cardinalities differ, or source and target are the same field.
    /// @throws std::runtime_error If a used source value is non-finite.
    void apply_piecewise_constant(const CellScalarField &source, CellScalarField &target) const;

    /// Computes source WLS gradients, then writes `phi_S + grad(phi_S) . (x_T - x_S)`.
    ///
    /// Boundary data and the caller-owned gradient workspace belong to the source
    /// mesh. Reuse this workspace across components/fields. Target values remain
    /// unchanged on failure; the workspace may have been overwritten by WLS.
    /// This reconstruction is unlimited and need not preserve local extrema.
    /// @throws std::invalid_argument For incompatible meshes, fields, workspace,
    ///         boundary cardinality, or in-place transfer.
    /// @throws std::runtime_error For non-finite source/reconstructed values or
    ///         an unusable source WLS stencil.
    void apply_linear_reconstruction(const CellScalarField &source,
                                     const ScalarBoundaryConditions &source_boundary_conditions,
                                     CellVectorField &source_gradient_workspace, CellScalarField &target) const;

  private:
    struct Entry
    {
        Index source_cell{};
        Vector2 displacement{};
    };

    void validate_fields(const CellScalarField &source, const CellScalarField &target) const;

    const Mesh *source_mesh_;
    const Mesh *target_mesh_;
    const Point2 *source_centers_;
    const Point2 *target_centers_;
    Index source_cell_count_;
    std::vector<Entry> entries_;
};

} // namespace cfd
