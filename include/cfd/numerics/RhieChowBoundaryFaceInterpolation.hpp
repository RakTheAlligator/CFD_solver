#pragma once

namespace cfd
{

class CellMomentumPressureResponse;
class CellScalarField;
class CellVectorField;
class CellVelocityField;
class FaceFluxField;
class FacePressureResponseField;
class Mesh;
class PressureCorrectionBoundaryConditions;
class ScalarBoundaryConditions;

/// Momentum-weighted interpolation for FixedPressure boundary faces.
///
/// For each FixedPressure boundary face, this component overwrites the
/// integrated outward provisional mass flux and the integrated mass-flux
/// response to pressure correction. FixedMassFlux and internal-face output
/// entries remain exactly under caller ownership.
///
/// @note The referenced Mesh is not owned and must outlive this object.
/// @note Repeated valid updates are O(Nfaces), perform no dynamic allocation,
///       and copy no field-sized arrays.
class RhieChowBoundaryFaceInterpolation
{
  public:
    /// Constructs boundary interpolation for a fixed Mesh and density.
    ///
    /// @throws std::invalid_argument If `density` is not finite and strictly
    ///         positive.
    RhieChowBoundaryFaceInterpolation(const Mesh &mesh, double density);

    RhieChowBoundaryFaceInterpolation(const RhieChowBoundaryFaceInterpolation &) = delete;
    RhieChowBoundaryFaceInterpolation &operator=(const RhieChowBoundaryFaceInterpolation &) = delete;

    RhieChowBoundaryFaceInterpolation(RhieChowBoundaryFaceInterpolation &&) noexcept = default;
    RhieChowBoundaryFaceInterpolation &operator=(RhieChowBoundaryFaceInterpolation &&) noexcept = delete;

    ~RhieChowBoundaryFaceInterpolation() = default;

    /// Overwrites provisional flux and pressure response on FixedPressure faces.
    ///
    /// The prescribed physical pressure is read from a Dirichlet scalar
    /// condition. For owner `P`, `d_b = x_b - x_P`, `q_b = D_P * S_b`,
    /// `c_b = (S_b . q_b) / (S_b . d_b)`, and
    /// `T_b = q_b - c_b * d_b`. With
    /// `H_P = U_P + D_P * grad(p)_P`, the outputs are
    /// `F_b_star = rho * (H_P . S_b - c_b * (p_b - p_P) - grad(p)_P . T_b)`
    /// and `Dp_b = rho * c_b`.
    ///
    /// FixedMassFlux boundaries and internal entries in both outputs are
    /// neither read nor modified.
    ///
    /// All used inputs and computed results are validated before either output
    /// field is modified.
    ///
    /// @throws std::invalid_argument If an input cardinality is incompatible,
    ///         or a FixedPressure boundary lacks a physical Dirichlet pressure
    ///         condition.
    /// @throws std::runtime_error If a used input or computed quantity is
    ///         non-finite, a used response is not strictly positive, or the
    ///         boundary geometry is unusable.
    void update_fixed_pressure_boundaries(
        const CellVelocityField &velocity, const CellScalarField &pressure, const CellVectorField &pressure_gradient,
        const CellMomentumPressureResponse &momentum_response,
        const ScalarBoundaryConditions &pressure_boundary_conditions,
        const PressureCorrectionBoundaryConditions &pressure_correction_boundary_conditions, FaceFluxField &mass_flux,
        FacePressureResponseField &face_pressure_response) const;

  private:
    const Mesh *mesh_;
    double density_;
};

} // namespace cfd
