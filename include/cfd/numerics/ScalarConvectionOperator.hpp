#pragma once

#include <span>

namespace cfd
{

class CellScalarField;
class FaceFluxField;
class Mesh;
class ScalarBoundaryConditions;
class ScalarLinearSystem;

/// First-order upwind finite-volume convection operator for a scalar field.
///
/// The supplied face flux is integrated, signed, and oriented outward from the
/// Mesh owner. Each internal face is evaluated once and contributes equal and
/// opposite balances to its owner and neighbor.
///
/// Dirichlet inflow uses the prescribed boundary value. Neumann inflow uses
/// the first-order closure `phi_b = phi_P + (d(phi)/dn) d_n`, where `d_n` is
/// the owner-centroid-to-face-center distance projected on the outward unit
/// normal. All outflow uses the owner value.
///
/// @note The referenced Mesh is not owned and must outlive this operator.
/// @note Repeated valid calls perform no dynamic allocation.
class ScalarConvectionOperator
{
  public:
    /// Constructs an operator referencing a fixed validated Mesh.
    explicit ScalarConvectionOperator(const Mesh &mesh) noexcept;

    ScalarConvectionOperator(const ScalarConvectionOperator &) = delete;
    ScalarConvectionOperator &operator=(const ScalarConvectionOperator &) = delete;

    ScalarConvectionOperator(ScalarConvectionOperator &&) noexcept = default;
    ScalarConvectionOperator &operator=(ScalarConvectionOperator &&) noexcept = delete;

    ~ScalarConvectionOperator() = default;

    /// Computes one integrated outward convective-flux balance per cell.
    ///
    /// The output is overwritten in full.
    ///
    /// @throws std::invalid_argument If a cardinality is incompatible or
    ///         `field` and `flux_balance` are the same object.
    void compute_flux_balance(const CellScalarField &field, const ScalarBoundaryConditions &boundary_conditions,
                              const FaceFluxField &face_flux, CellScalarField &flux_balance) const;

    /// Adds first-order upwind coefficients to an existing system matrix.
    ///
    /// @throws std::invalid_argument If a cardinality is incompatible or
    ///         `system` does not reference this operator's exact Mesh instance.
    void add_matrix_contributions(const ScalarBoundaryConditions &boundary_conditions, const FaceFluxField &face_flux,
                                  ScalarLinearSystem &system) const;

    /// Adds boundary contributions to `rhs` without clearing it.
    ///
    /// With the assembly convention used here,
    /// `A * phi - b_boundary` equals the convective flux balance.
    ///
    /// @throws std::invalid_argument If a cardinality is incompatible.
    void add_boundary_rhs(const ScalarBoundaryConditions &boundary_conditions, const FaceFluxField &face_flux,
                          std::span<double> rhs) const;

  private:
    const Mesh *mesh_;
};

} // namespace cfd
