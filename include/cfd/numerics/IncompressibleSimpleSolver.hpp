#pragma once

#include "cfd/field/CellMomentumPressureResponse.hpp"
#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVectorField.hpp"
#include "cfd/field/CellVelocityField.hpp"
#include "cfd/field/FaceFluxField.hpp"
#include "cfd/field/FacePressureResponseField.hpp"
#include "cfd/linear_algebra/EigenBiCGSTABSolver.hpp"
#include "cfd/linear_algebra/EigenConjugateGradientSolver.hpp"
#include "cfd/linear_algebra/ScalarLinearSystem.hpp"
#include "cfd/mesh/Types.hpp"
#include "cfd/numerics/IncompressibleMomentumAssembler.hpp"
#include "cfd/numerics/IncompressiblePressureCorrectionAssembler.hpp"
#include "cfd/numerics/IncompressiblePressureVelocityCorrection.hpp"
#include "cfd/numerics/RhieChowBoundaryFaceInterpolation.hpp"
#include "cfd/numerics/RhieChowInternalFaceInterpolation.hpp"
#include "cfd/numerics/ScalarConvectionOperator.hpp"

#include <functional>
#include <vector>

namespace cfd
{

class Mesh;
class PressureCorrectionBoundaryConditions;
class ScalarBoundaryConditions;

/// Controls for the steady incompressible SIMPLE v1 iteration.
struct IncompressibleSimpleOptions
{
    Index maximum_iterations{500};
    double momentum_relaxation_factor{1.0};
    double pressure_relaxation_factor{0.3};
    /// Relaxes computed integrated owner-oriented Rhie-Chow mass fluxes before pressure correction.
    ///
    /// Internal and `FixedPressure` boundary faces use
    /// `F_provisional = alpha_rc * F_RhieChow + (1 - alpha_rc) * F_previous`,
    /// where `F_previous` is the corrected flux from the preceding iteration.
    /// `FixedMassFlux` boundary values and face pressure-response coefficients
    /// are unchanged. A value of one preserves the unrelaxed path.
    double rhie_chow_flux_relaxation_factor{1.0};
    double velocity_relative_tolerance{1.0e-8};
    /// Tolerance for provisional continuity, which controls outer convergence.
    double continuity_relative_tolerance{1.0e-10};
    BiCGSTABOptions momentum_linear_solver{};
    ConjugateGradientOptions pressure_correction_linear_solver{};
};

/// Linear-solve outcomes and nonlinear diagnostics for one completed SIMPLE iteration.
///
/// The x- and y-velocity equation residuals are normalized assembled-momentum
/// imbalances `sum(|b - A*x|) / (sum(|b|) + sum(|A*x|))`, evaluated with the
/// iteration-start velocity before the corresponding linear solve. They are
/// Fluent-like outer monitoring quantities, not Eigen inner-solver residuals
/// and not a reproduction of another solver's exact residual scaling.
struct SimpleIterationInfo
{
    Index iteration{};
    LinearSolveResult u_solve{};
    LinearSolveResult v_solve{};
    LinearSolveResult pressure_correction_solve{};
    double x_velocity_equation_residual{};
    double y_velocity_equation_residual{};
    double velocity_relative_change{};
    /// Default continuity monitor: provisional continuity before pressure correction.
    double provisional_continuity_relative_residual{};
    double corrected_continuity_relative_residual{};
    double maximum_pressure_correction{};
};

/// Observer invoked after each fully completed SIMPLE outer iteration.
using SimpleIterationCallback = std::function<void(const SimpleIterationInfo &)>;

/// Outcome and final diagnostics of one steady SIMPLE solve.
struct IncompressibleSimpleResult
{
    bool converged{};
    Index iteration_count{};
    double velocity_relative_change{};
    /// Relative continuity residual of the provisional flux before correction.
    double provisional_continuity_relative_residual{};
    /// Relative continuity residual of the corrected final flux.
    double continuity_relative_residual{};
    double maximum_mass_imbalance{};
    double maximum_pressure_correction{};
};

/// Reusable steady incompressible SIMPLE v1 solver.
///
/// This class orchestrates the existing momentum, Rhie-Chow,
/// pressure-correction, field-correction, and linear-solver components. It owns
/// fixed-cardinality workspace and reuses it across outer iterations.
///
/// SIMPLE v1 deliberately requires a momentum relaxation factor of exactly
/// one. The current momentum-weighted interpolation uses the final momentum
/// diagonal and does not yet implement a relaxation-consistent Majumdar
/// treatment. Pressure under-relaxation remains supported.
///
/// Computed provisional face mass fluxes may be relaxed independently after
/// Rhie-Chow interpolation and before pressure-correction assembly. This is
/// provisional face-flux relaxation, not Majumdar momentum under-relaxation.
///
/// @note The referenced Mesh is not owned and must outlive this solver.
/// @note Construction allocates the field and system workspace. Outer
///       iterations do not allocate field-sized temporary arrays; Eigen sparse
///       matrix preparation may allocate as documented by the linear backends.
class IncompressibleSimpleSolver
{
  public:
    /// Constructs a SIMPLE solver and its fixed workspace.
    ///
    /// @throws std::invalid_argument If a physical coefficient or option is
    ///         invalid, momentum relaxation differs from one, or the convection
    ///         scheme is unsupported.
    /// @throws std::runtime_error If an existing numerical component rejects
    ///         the Mesh geometry.
    IncompressibleSimpleSolver(const Mesh &mesh, double density, double dynamic_viscosity,
                               ScalarConvectionScheme convection_scheme, IncompressibleSimpleOptions options = {});

    IncompressibleSimpleSolver(const IncompressibleSimpleSolver &) = delete;
    IncompressibleSimpleSolver &operator=(const IncompressibleSimpleSolver &) = delete;
    IncompressibleSimpleSolver(IncompressibleSimpleSolver &&) = delete;
    IncompressibleSimpleSolver &operator=(IncompressibleSimpleSolver &&) = delete;

    ~IncompressibleSimpleSolver() = default;

    /// Iterates until velocity change and provisional continuity satisfy their
    /// tolerances. Corrected continuity remains available as a diagnostic.
    ///
    /// `FixedPressure` pressure-correction conditions must correspond to
    /// physical pressure Dirichlet conditions. `FixedMassFlux` conditions must
    /// correspond to physical pressure Neumann conditions, and their boundary
    /// mass fluxes remain exactly caller-owned and unchanged.
    ///
    /// For reconstructing `grad(p')` only, `FixedPressure` maps to homogeneous
    /// scalar Dirichlet data and `FixedMassFlux` maps to homogeneous scalar
    /// Neumann data. This reconstruction mapping does not replace the distinct
    /// face-flux correction semantics.
    ///
    /// @return Final diagnostics. Reaching the outer iteration limit returns
    ///         `converged == false`.
    /// @param iteration_callback Optional observer called after convergence is
    ///        evaluated for every successfully completed outer iteration,
    ///        including the final one. An iteration interrupted by an inner
    ///        linear-solver failure is not reported.
    /// @throws std::invalid_argument If a cardinality, initial field value, or
    ///         boundary-condition pairing is invalid.
    /// @throws std::runtime_error If imposed all-`FixedMassFlux` boundary flow
    ///         is globally incompatible, an inner linear solve fails to
    ///         converge, or an existing numerical component rejects its data.
    [[nodiscard]]
    IncompressibleSimpleResult solve(
        const ScalarBoundaryConditions &u_boundary_conditions, const ScalarBoundaryConditions &v_boundary_conditions,
        const ScalarBoundaryConditions &pressure_boundary_conditions,
        const PressureCorrectionBoundaryConditions &pressure_correction_boundary_conditions,
        CellVelocityField &velocity, CellScalarField &pressure, FaceFluxField &mass_flux,
        SimpleIterationCallback iteration_callback = {});

  private:
    const Mesh *mesh_;
    double density_;
    IncompressibleSimpleOptions options_;

    EigenBiCGSTABSolver u_momentum_solver_;
    EigenBiCGSTABSolver v_momentum_solver_;
    EigenConjugateGradientSolver pressure_correction_solver_;

    IncompressibleMomentumAssembler momentum_assembler_;
    RhieChowInternalFaceInterpolation internal_face_interpolation_;
    RhieChowBoundaryFaceInterpolation boundary_face_interpolation_;
    IncompressiblePressureCorrectionAssembler pressure_correction_assembler_;
    IncompressiblePressureVelocityCorrection pressure_velocity_corrector_;

    CellVelocityField previous_velocity_;
    FaceFluxField previous_mass_flux_;
    CellVectorField u_gradient_;
    CellVectorField v_gradient_;
    CellVectorField pressure_gradient_;
    CellVectorField pressure_correction_gradient_;
    CellMomentumPressureResponse momentum_response_;
    FacePressureResponseField face_pressure_response_;
    CellScalarField pressure_correction_;
    CellScalarField mass_imbalance_;
    std::vector<double> momentum_matrix_product_workspace_;

    ScalarLinearSystem u_momentum_system_;
    ScalarLinearSystem v_momentum_system_;
    ScalarLinearSystem pressure_correction_system_;
};

} // namespace cfd
