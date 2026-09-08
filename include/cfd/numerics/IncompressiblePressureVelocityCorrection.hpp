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

/// Applies a solved pressure correction to incompressible-flow state fields.
///
/// The component updates caller-owned face flux, pressure, and velocity
/// storage independently. It owns no numerical fields and performs no gradient
/// reconstruction or linear solve.
///
/// @note The referenced Mesh is not owned and must outlive this object.
/// @note Repeated valid corrections perform no dynamic allocation and copy no
///       field-sized arrays.
class IncompressiblePressureVelocityCorrection
{
  public:
    /// Constructs a correction component for a fixed Mesh.
    explicit IncompressiblePressureVelocityCorrection(const Mesh &mesh) noexcept;

    IncompressiblePressureVelocityCorrection(const IncompressiblePressureVelocityCorrection &) = delete;
    IncompressiblePressureVelocityCorrection &operator=(const IncompressiblePressureVelocityCorrection &) = delete;

    IncompressiblePressureVelocityCorrection(IncompressiblePressureVelocityCorrection &&) noexcept = default;
    IncompressiblePressureVelocityCorrection &operator=(IncompressiblePressureVelocityCorrection &&) noexcept = delete;

    ~IncompressiblePressureVelocityCorrection() = default;

    /// Corrects integrated owner-oriented mass fluxes in place.
    ///
    /// An internal owner `P` / neighbor `N` face receives
    /// `Dp_f * (p'_P - p'_N)`. A `FixedPressure` boundary receives
    /// `Dp_b * p'_P`, because `p'_b = 0`. `FixedMassFlux` entries remain
    /// exactly unchanged, and their pressure-response values are neither read
    /// nor validated. Pressure relaxation does not apply to face-flux
    /// correction.
    ///
    /// All used inputs and computed corrected fluxes are validated before any
    /// face value is modified.
    ///
    /// @throws std::invalid_argument If a field or boundary-condition
    ///         cardinality is incompatible with the Mesh.
    /// @throws std::runtime_error If a used input or corrected flux is
    ///         non-finite, or a used pressure response is not strictly
    ///         positive.
    void correct_face_mass_flux(const CellScalarField &pressure_correction,
                                const PressureCorrectionBoundaryConditions &boundary_conditions,
                                const FacePressureResponseField &face_pressure_response,
                                FaceFluxField &mass_flux) const;

    /// Applies the relaxed pressure update `p = p + alpha_p * p'` in place.
    ///
    /// The relaxation factor is used only for this pressure update. All input
    /// and computed pressure values are validated before pressure is modified.
    ///
    /// @throws std::invalid_argument If a field cardinality is incompatible,
    ///         the fields alias, or `pressure_relaxation_factor` is not finite
    ///         and in `(0, 1]`.
    /// @throws std::runtime_error If an input or corrected pressure is
    ///         non-finite.
    void correct_pressure(const CellScalarField &pressure_correction, double pressure_relaxation_factor,
                          CellScalarField &pressure) const;

    /// Applies `U = U_star - D_P * grad(p')` component-wise in place.
    ///
    /// Specifically, `u -= d_Pu * grad(p').x` and
    /// `v -= d_Pv * grad(p').y`. Pressure relaxation does not apply to this
    /// velocity correction. The caller supplies the already reconstructed
    /// pressure-correction gradient.
    ///
    /// All used inputs and computed velocities are validated before either
    /// velocity component is modified.
    ///
    /// @throws std::invalid_argument If a field cardinality is incompatible
    ///         with the Mesh.
    /// @throws std::runtime_error If a used input or corrected velocity is
    ///         non-finite, or a momentum response is not strictly positive.
    void correct_velocity(const CellVectorField &pressure_correction_gradient,
                          const CellMomentumPressureResponse &momentum_response, CellVelocityField &velocity) const;

  private:
    const Mesh *mesh_;
};

} // namespace cfd
