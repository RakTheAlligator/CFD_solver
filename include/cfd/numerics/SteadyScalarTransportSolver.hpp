#pragma once

#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVectorField.hpp"
#include "cfd/linear_algebra/EigenBiCGSTABSolver.hpp"
#include "cfd/linear_algebra/ScalarLinearSystem.hpp"
#include "cfd/numerics/ScalarTransportAssembler.hpp"

#include <vector>

namespace cfd
{

/// Controls the explicit-correction fixed-point iteration and its inner solve.
struct SteadyScalarTransportOptions
{
    Index maximum_iterations{1000};
    double field_relative_tolerance{1.0e-8};
    double discrete_relative_tolerance{1.0e-8};
    BiCGSTABOptions linear_solver{};
};

/// Diagnostics from the last completed transport iteration.
struct SteadyScalarTransportResult
{
    bool converged{};
    Index iteration_count{};
    /// Inner solve for the RHS frozen before the last field update.
    LinearSolveResult linear_solve;
    double field_relative_change{};
    /// Full equation residual with corrections rebuilt from the updated field.
    double discrete_relative_residual{};
};

/// Steady scalar transport with fixed supplied flux and constant Gamma.
///
/// Owns composed operators, a linear system and persistent work buffers. Mesh
/// is not owned and must outlive this solver. The caller must supply a well-posed
/// problem; no general singularity detector or Neumann gauge is provided.
/// No field/equation relaxation or coupling to SIMPLE is performed.
class SteadyScalarTransportSolver
{
  public:
    /// @throws std::invalid_argument If options, coefficient or scheme/limiter
    ///         selection are invalid.
    /// @throws std::runtime_error If composed operators reject face geometry.
    SteadyScalarTransportSolver(const Mesh &mesh, double diffusion_coefficient,
                                ScalarConvectionScheme scheme = ScalarConvectionScheme::FirstOrderUpwind,
                                ScalarConvectionLimiter limiter = ScalarConvectionLimiter::None,
                                SteadyScalarTransportOptions options = {});

    SteadyScalarTransportSolver(const SteadyScalarTransportSolver &) = delete;
    SteadyScalarTransportSolver &operator=(const SteadyScalarTransportSolver &) = delete;
    SteadyScalarTransportSolver(SteadyScalarTransportSolver &&) = delete;
    SteadyScalarTransportSolver &operator=(SteadyScalarTransportSolver &&) = delete;
    ~SteadyScalarTransportSolver() = default;

    /// Prepares the matrix once, then closes deferred/non-orthogonal corrections.
    ///
    /// field is caller-owned, supplies the initial guess and is updated in place.
    /// Flux and BC are borrowed only for this call and must remain unchanged.
    /// Convergence requires both `||phi_new-phi_old|| / max(||phi_new||,||phi_old||)`
    /// and `||A phi_new-rhs(phi_new)|| / max(||A phi_new||,||rhs(phi_new)||)` to
    /// meet their tolerances. Norms are Euclidean; 0/0 is zero, nonzero/0 is +inf.
    /// The inner Eigen residual is a separate diagnostic, not the outer criterion.
    ///
    /// Reaching the outer limit returns converged=false and retains the latest
    /// field. No explicit allocation is performed in the outer loop; the
    /// existing Eigen backend may allocate its own solve work storage.
    /// @throws std::invalid_argument If cardinalities differ or inputs are non-finite.
    /// @throws std::runtime_error If WLS or the inner solve fails, or diagnostics
    ///         become non-finite. A failed solve may have modified field.
    [[nodiscard]]
    SteadyScalarTransportResult solve(const FaceFluxField &face_flux,
                                      const ScalarBoundaryConditions &boundary_conditions, CellScalarField &field);

  private:
    const Mesh *mesh_;
    SteadyScalarTransportOptions options_;
    ScalarTransportAssembler assembler_;
    ScalarLinearSystem system_;
    EigenBiCGSTABSolver linear_solver_;
    CellVectorField gradient_;
    CellScalarField previous_field_;
    std::vector<double> matrix_product_;
};

} // namespace cfd
