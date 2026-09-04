#pragma once

#include "cfd/math/Vector2.hpp"

#include <vector>

namespace cfd
{

class CellScalarField;
class CellVectorField;
class Mesh;
class ScalarBoundaryConditions;

/// Constant-isotropic finite-volume diffusion operator in two dimensions.
///
/// The operator precomputes fixed face geometry for `-div(Gamma grad(phi))`
/// and returns one integrated outward diffusive-flux balance per cell. Internal
/// faces use an over-relaxed corrected two-point decomposition and contribute
/// equal-and-opposite fluxes to their adjacent cells.
///
/// Dirichlet boundary values prescribe `phi` at the face center. Neumann
/// values prescribe the outward normal derivative `d(phi)/dn`, not physical
/// diffusive flux.
///
/// @note The referenced Mesh is not owned and must outlive this operator.
/// @note Repeated valid calls to `compute_flux_balance()` perform no dynamic
///       allocation; the caller owns the fixed-cardinality output field.
class ScalarDiffusionOperator
{
  public:
    /// Precomputes face coefficients for a fixed Mesh and diffusivity.
    ///
    /// @param mesh Validated mesh whose storage must outlive this operator.
    /// @param diffusivity Constant positive isotropic diffusivity `Gamma`.
    /// @throws std::invalid_argument If `diffusivity` is non-finite or not
    ///         strictly positive.
    /// @throws std::runtime_error If the mesh geometry is numerically
    ///         degenerate or an internal-face supporting line does not
    ///         intersect the owner-neighbor segment.
    ScalarDiffusionOperator(const Mesh &mesh, double diffusivity);

    ScalarDiffusionOperator(const ScalarDiffusionOperator &) = delete;
    ScalarDiffusionOperator &operator=(const ScalarDiffusionOperator &) = delete;

    ScalarDiffusionOperator(ScalarDiffusionOperator &&) noexcept = default;
    ScalarDiffusionOperator &operator=(ScalarDiffusionOperator &&) noexcept = delete;

    ~ScalarDiffusionOperator() = default;

    /// Computes the integrated outward diffusive-flux balance of every cell.
    ///
    /// Each internal face is evaluated once in owner orientation. Its flux is
    /// added to the owner and subtracted from the neighbor. Boundary values are
    /// read on every call, so their type and value may change without rebuilding
    /// geometric coefficients.
    ///
    /// @param field Cell-centered values interpreted at area centroids.
    /// @param boundary_conditions One condition per Mesh boundary group.
    /// @param gradient Cell-centered gradients used by non-orthogonal
    ///        corrections.
    /// @param flux_balance Caller-owned output, overwritten in full.
    /// @throws std::invalid_argument If any cardinality is incompatible with
    ///         the Mesh or if `field` and `flux_balance` are the same object.
    void compute_flux_balance(const CellScalarField &field, const ScalarBoundaryConditions &boundary_conditions,
                              const CellVectorField &gradient, CellScalarField &flux_balance) const;

  private:
    struct FaceData
    {
        Vector2 correction_flux_vector;
        double primary_coefficient{};
        double neighbor_gradient_weight{};
    };

    const Mesh *mesh_;
    double diffusivity_;
    std::vector<FaceData> face_data_;
};

} // namespace cfd
