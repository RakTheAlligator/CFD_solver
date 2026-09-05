#pragma once

#include "cfd/linear_algebra/LinearSolveResult.hpp"
#include "cfd/mesh/Types.hpp"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>

#include <span>

namespace cfd
{

class ScalarLinearSystem;

/// Runtime controls for the Eigen BiCGSTAB backend.
struct BiCGSTABOptions
{
    double relative_tolerance{1.0e-10};
    Index maximum_iterations{1000};
};

/// Eigen BiCGSTAB solver for nonsymmetric scalar systems.
///
/// `compute_matrix()` copies the directed finite-volume coefficients into
/// Eigen sparse storage and prepares a diagonal preconditioner. Repeated
/// `solve()` calls reuse both objects and use the supplied solution as the
/// initial guess.
class EigenBiCGSTABSolver
{
  public:
    explicit EigenBiCGSTABSolver(BiCGSTABOptions options = {});

    EigenBiCGSTABSolver(const EigenBiCGSTABSolver &) = delete;
    EigenBiCGSTABSolver &operator=(const EigenBiCGSTABSolver &) = delete;
    EigenBiCGSTABSolver(EigenBiCGSTABSolver &&) = delete;
    EigenBiCGSTABSolver &operator=(EigenBiCGSTABSolver &&) = delete;

    ~EigenBiCGSTABSolver() = default;

    /// Converts and prepares a matrix for subsequent solves.
    ///
    /// Directed off-diagonal coefficients are copied without symmetrization.
    /// Boundary-face off-diagonal storage is ignored.
    ///
    /// @throws std::invalid_argument If a matrix coefficient is non-finite or
    ///         the matrix cardinality is unsupported by Eigen.
    /// @throws std::runtime_error If Eigen cannot prepare the matrix.
    void compute_matrix(const ScalarLinearSystem &system);

    /// Solves the prepared system using `solution` as the initial guess.
    ///
    /// Non-convergence within the configured iteration limit is reported by a
    /// result with `converged == false`.
    ///
    /// @throws std::logic_error If `compute_matrix()` has not succeeded.
    /// @throws std::invalid_argument If cardinalities differ, spans overlap,
    ///         or an input value is non-finite.
    /// @throws std::runtime_error If Eigen reports a numerical or input error.
    [[nodiscard]]
    LinearSolveResult solve(std::span<const double> rhs, std::span<double> solution);

  private:
    using SparseMatrix = Eigen::SparseMatrix<double>;
    using Solver = Eigen::BiCGSTAB<SparseMatrix, Eigen::DiagonalPreconditioner<double>>;

    BiCGSTABOptions options_;
    SparseMatrix matrix_;
    Solver solver_;
    Index matrix_size_{};
    bool matrix_is_prepared_{};
};

} // namespace cfd
