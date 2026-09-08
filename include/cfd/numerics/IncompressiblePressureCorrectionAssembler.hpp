#pragma once

namespace cfd
{

class FaceFluxField;
class FacePressureResponseField;
class Mesh;
class ScalarLinearSystem;

/// Assembles incompressible pressure-correction contributions.
///
/// `F_star` is the integrated owner-oriented provisional mass flux, and `Dp_f`
/// is the integrated mass-flux response to the pressure difference. For each
/// internal owner `P` / neighbor `N` face, the correction is
/// `F_new = F_star + Dp_f * (p'_P - p'_N)`.
///
/// @note The referenced Mesh is not owned and must outlive this assembler.
/// @note Repeated valid calls perform no dynamic allocation and copy no
///       field-sized arrays.
class IncompressiblePressureCorrectionAssembler
{
  public:
    /// Constructs a pressure-correction assembler for a fixed Mesh.
    explicit IncompressiblePressureCorrectionAssembler(const Mesh &mesh) noexcept;

    IncompressiblePressureCorrectionAssembler(const IncompressiblePressureCorrectionAssembler &) = delete;
    IncompressiblePressureCorrectionAssembler &operator=(const IncompressiblePressureCorrectionAssembler &) = delete;

    IncompressiblePressureCorrectionAssembler(IncompressiblePressureCorrectionAssembler &&) noexcept = default;
    IncompressiblePressureCorrectionAssembler &operator=(IncompressiblePressureCorrectionAssembler &&) noexcept =
        delete;

    ~IncompressiblePressureCorrectionAssembler() = default;

    /// Adds all internal-face matrix and provisional-flux contributions.
    ///
    /// Each face additively contributes `Dp_f * [[1, -1], [-1, 1]]` to the
    /// local `(P, N)` matrix and `[-F_star, +F_star]` to the local right-hand
    /// side. Existing system entries are never cleared.
    ///
    /// Boundary values in both input face fields are neither read nor
    /// validated. All used values and output-system invariants are validated
    /// before the system is modified.
    ///
    /// @note This internal-face component alone does not form a complete
    ///       solvable pressure-correction system. Boundary pressure-response
    ///       treatment and pressure-reference/gauge fixing remain intentionally
    ///       absent.
    ///
    /// @throws std::invalid_argument If a face-field cardinality is
    ///         incompatible, or `system` does not reference this assembler's
    ///         exact Mesh with matching cardinalities.
    /// @throws std::runtime_error If a used provisional mass flux is not finite,
    ///         or a used face pressure response is not finite or not strictly
    ///         positive.
    void add_internal_face_contributions(const FaceFluxField &provisional_mass_flux,
                                         const FacePressureResponseField &face_pressure_response,
                                         ScalarLinearSystem &system) const;

    /// Adds the provisional boundary mass imbalance to the right-hand side.
    ///
    /// `F_b_star` is the integrated outward, owner-oriented provisional mass
    /// flux. Each boundary face contributes `-F_b_star` to its owner RHS; the
    /// matrix is unchanged. Combined with `add_internal_face_contributions()`,
    /// the provisional-flux increment is
    /// `delta_rhs[P] = -sum_outward(F_f_star)` before later boundary
    /// pressure-response terms are added.
    ///
    /// Internal face values are neither read nor validated. This method does
    /// not define how pressure correction changes boundary fluxes.
    ///
    /// @throws std::invalid_argument If the face-field cardinality is
    ///         incompatible, or `system` does not reference this assembler's
    ///         exact Mesh with matching cardinalities.
    /// @throws std::runtime_error If a used boundary provisional mass flux is
    ///         not finite.
    void add_boundary_provisional_flux_rhs(const FaceFluxField &provisional_mass_flux,
                                           ScalarLinearSystem &system) const;

  private:
    const Mesh *mesh_;
};

} // namespace cfd
