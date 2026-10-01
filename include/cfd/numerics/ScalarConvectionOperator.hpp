#pragma once

#include "cfd/mesh/Types.hpp"

#include <cstdint>
#include <span>
#include <vector>

namespace cfd
{

class CellScalarField;
class CellVectorField;
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
    Hybrid,
    /// Uses an implicit FirstOrderUpwind matrix plus an explicit WLS-gradient
    /// reconstruction from the upwind cell.
    LinearUpwind
};

/// Optional limiter applied to a scalar convection reconstruction.
enum class ScalarConvectionLimiter : std::uint8_t
{
    /// Leaves the selected convection scheme unchanged.
    None,
    /// Bounds LinearUpwind face reconstructions with one coefficient per cell.
    BarthJespersen
};

/// Counts faces selected by each branch of the Hybrid convection criterion.
struct HybridConvectionFaceCounts
{
    Index internal_linear_faces{};
    Index internal_upwind_faces{};
    Index boundary_linear_faces{};
    Index boundary_upwind_faces{};
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
/// face by face using caller-supplied diffusion conductances. LinearUpwind
/// retains the FirstOrderUpwind matrix and adds its higher-order reconstruction
/// explicitly to the right-hand side.
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
    /// @throws std::invalid_argument If the scheme/limiter combination is unsupported.
    /// @throws std::runtime_error If Linear or Hybrid interpolation encounters
    ///         unusable internal-face geometry.
    ScalarConvectionOperator(const Mesh &mesh, ScalarConvectionScheme scheme,
                             ScalarConvectionLimiter limiter = ScalarConvectionLimiter::None);

    ScalarConvectionOperator(const ScalarConvectionOperator &) = delete;
    ScalarConvectionOperator &operator=(const ScalarConvectionOperator &) = delete;

    ScalarConvectionOperator(ScalarConvectionOperator &&) noexcept = default;
    ScalarConvectionOperator &operator=(ScalarConvectionOperator &&) noexcept = delete;

    ~ScalarConvectionOperator() = default;

    /// Classifies every face using the Hybrid convection criterion.
    ///
    /// This explicit diagnostic performs one face traversal, allocates no
    /// storage, and does not alter the operator or its inputs.
    ///
    /// @throws std::invalid_argument If this operator is not configured for
    ///         Hybrid, a cardinality is incompatible, or a conductance is
    ///         non-finite or not strictly positive.
    [[nodiscard]]
    HybridConvectionFaceCounts classify_hybrid_faces(const FaceFluxField &face_flux,
                                                     std::span<const double> face_diffusion_conductances) const;

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

    /// Computes the complete LinearUpwind convective-flux balance.
    ///
    /// Internal faces reconstruct from the upwind cell. Boundary inflow retains
    /// the FirstOrderUpwind boundary closure, while boundary outflow reconstructs
    /// from the owner cell. The supplied gradient is treated as explicit data.
    ///
    /// @throws std::invalid_argument If this operator is not configured for
    ///         unlimited LinearUpwind, a cardinality is incompatible, or `field` and
    ///         `flux_balance` are the same object.
    void compute_flux_balance(const CellScalarField &field, const ScalarBoundaryConditions &boundary_conditions,
                              const FaceFluxField &face_flux, const CellVectorField &gradient,
                              CellScalarField &flux_balance) const;

    /// Computes the complete optionally limited LinearUpwind convective-flux balance.
    ///
    /// For BarthJespersen, `limiter_workspace` is overwritten with one
    /// coefficient per mesh cell and must not overlap `field`, `face_flux`,
    /// or `flux_balance` storage.
    ///
    /// @throws std::invalid_argument If this operator is not configured for
    ///         LinearUpwind, a cardinality is incompatible, an input required
    ///         by the limiter is non-finite, or storage ranges overlap.
    void compute_flux_balance(const CellScalarField &field, const ScalarBoundaryConditions &boundary_conditions,
                              const FaceFluxField &face_flux, const CellVectorField &gradient,
                              std::span<double> limiter_workspace, CellScalarField &flux_balance) const;

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

    /// Adds the selected scheme's deferred correction to `rhs`.
    ///
    /// LinearUpwind adds `-F delta_phi` to the owner and `+F delta_phi` to the
    /// neighbor of each internal face. Positive-flux boundary faces add the
    /// corresponding owner reconstruction. Other schemes perform no work.
    ///
    /// @throws std::invalid_argument For LinearUpwind, if a cardinality is
    ///         incompatible, this operator has a limiter, or `rhs` overlaps
    ///         `face_flux` storage.
    void add_deferred_correction_rhs(const CellVectorField &gradient, const FaceFluxField &face_flux,
                                     std::span<double> rhs) const;

    /// Adds the optionally limited LinearUpwind deferred correction to `rhs`.
    ///
    /// BarthJespersen derives one coefficient per cell from `field`, its
    /// boundary conditions, and the raw gradient. `limiter_workspace` is
    /// overwritten and must not overlap `field`, `face_flux`, or `rhs`.
    /// The `rhs` storage must not overlap `field` or `face_flux`. Other
    /// schemes perform no work.
    ///
    /// @throws std::invalid_argument For LinearUpwind, if a cardinality is
    ///         incompatible, an input required by the limiter is non-finite,
    ///         or storage ranges overlap.
    void add_deferred_correction_rhs(const CellScalarField &field, const ScalarBoundaryConditions &boundary_conditions,
                                     const CellVectorField &gradient, const FaceFluxField &face_flux,
                                     std::span<double> limiter_workspace, std::span<double> rhs) const;

  private:
    void compute_barth_jespersen_limiter(const CellScalarField &field,
                                         const ScalarBoundaryConditions &boundary_conditions,
                                         const CellVectorField &gradient, std::span<double> limiter_coefficients) const;

    void add_linear_upwind_deferred_correction(const CellVectorField &gradient, const FaceFluxField &face_flux,
                                               std::span<const double> limiter_coefficients,
                                               std::span<double> rhs) const;

    const Mesh *mesh_;
    ScalarConvectionScheme scheme_;
    ScalarConvectionLimiter limiter_;
    std::vector<double> internal_face_interpolation_weights_;
};

} // namespace cfd
