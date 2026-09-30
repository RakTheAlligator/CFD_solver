#pragma once

#include <cstdint>
#include <span>
#include <vector>

namespace cfd
{

class CellScalarField;
class FaceFluxField;
class Mesh;
class ScalarBoundaryConditions;
class ScalarLinearSystem;

/// Available finite-volume scalar convection interpolation schemes.
enum class ScalarConvectionScheme : std::uint8_t
{
    /// Robust first-order interpolation from the upwind side of each face.
    FirstOrderUpwind,
    /// Geometry-aware linear interpolation between adjacent cell centers.
    ///
    /// On a uniform Cartesian mesh this is classical centered interpolation.
    Linear,
    /// Selects Linear or FirstOrderUpwind independently on each face from the
    /// local mass flux, interpolation weight, and diffusion conductance.
    Hybrid
};

/// Finite-volume convection operator for a scalar field.
///
/// The supplied face flux is integrated, signed, and oriented outward from the
/// Mesh owner. Each internal face is evaluated once and contributes equal and
/// opposite balances to its owner and neighbor.
///
/// FirstOrderUpwind uses the existing flow-directed boundary treatment. Linear
/// applies boundary conditions independently of flow direction: Dirichlet uses
/// the prescribed face value and Neumann uses the first-order closure
/// `phi_b = phi_P + (d(phi)/dn) d_n`. Hybrid selects between those treatments
/// face by face using caller-supplied diffusion conductances.
///
/// @note The referenced Mesh is not owned and must outlive this operator.
/// @note Repeated valid calls perform no dynamic allocation.
class ScalarConvectionOperator
{
  public:
    /// Constructs the historical FirstOrderUpwind operator for a fixed Mesh.
    ///
    /// This overload performs no allocation or geometry preprocessing.
    explicit ScalarConvectionOperator(const Mesh &mesh) noexcept;

    /// Constructs an operator with an explicitly selected interpolation scheme.
    ///
    /// @throws std::invalid_argument If `scheme` is unsupported.
    /// @throws std::runtime_error If Linear or Hybrid interpolation encounters
    ///         unusable internal-face geometry.
    ScalarConvectionOperator(const Mesh &mesh, ScalarConvectionScheme scheme);

    ScalarConvectionOperator(const ScalarConvectionOperator &) = delete;
    ScalarConvectionOperator &operator=(const ScalarConvectionOperator &) = delete;

    ScalarConvectionOperator(ScalarConvectionOperator &&) noexcept = default;
    ScalarConvectionOperator &operator=(ScalarConvectionOperator &&) noexcept = delete;

    ~ScalarConvectionOperator() = default;

    /// Computes one integrated outward convective-flux balance per cell.
    ///
    /// The output is overwritten in full. For Hybrid,
    /// `face_diffusion_conductances[face]` is the integrated principal
    /// diffusion coefficient for that same face and must be finite and positive.
    ///
    /// @throws std::invalid_argument If a cardinality is incompatible, `field`
    ///         and `flux_balance` are the same object, or Hybrid conductances
    ///         are absent, have incorrect cardinality, are non-finite, or are
    ///         not strictly positive.
    void compute_flux_balance(const CellScalarField &field, const ScalarBoundaryConditions &boundary_conditions,
                              const FaceFluxField &face_flux, CellScalarField &flux_balance,
                              std::span<const double> face_diffusion_conductances = {}) const;

    /// Adds coefficients for the selected scheme to an existing system matrix.
    ///
    /// For Hybrid, `face_diffusion_conductances[face]` is the integrated
    /// principal diffusion coefficient for that same face and must be finite
    /// and positive.
    ///
    /// @throws std::invalid_argument If a cardinality is incompatible, `system`
    ///         does not reference this operator's exact Mesh instance, or Hybrid
    ///         conductances are absent, have incorrect cardinality, are
    ///         non-finite, or are not strictly positive.
    void add_matrix_contributions(const ScalarBoundaryConditions &boundary_conditions, const FaceFluxField &face_flux,
                                  ScalarLinearSystem &system,
                                  std::span<const double> face_diffusion_conductances = {}) const;

    /// Adds boundary contributions to `rhs` without clearing it.
    ///
    /// With the assembly convention used here,
    /// `A * phi - b_boundary` equals the convective flux balance. For Hybrid,
    /// `face_diffusion_conductances[face]` is the integrated principal
    /// diffusion coefficient for that same face and must be finite and positive.
    ///
    /// @throws std::invalid_argument If a cardinality is incompatible or Hybrid
    ///         conductances are absent, have incorrect cardinality, are
    ///         non-finite, or are not strictly positive.
    void add_boundary_rhs(const ScalarBoundaryConditions &boundary_conditions, const FaceFluxField &face_flux,
                          std::span<double> rhs, std::span<const double> face_diffusion_conductances = {}) const;

  private:
    const Mesh *mesh_;
    ScalarConvectionScheme scheme_;
    std::vector<double> internal_face_interpolation_weights_;
};

} // namespace cfd
