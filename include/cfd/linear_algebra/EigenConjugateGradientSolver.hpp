#pragma once

#include "cfd/mesh/Types.hpp"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>

#include <span>

namespace cfd
{

class ScalarLinearSystem;

/// Runtime controls for the Eigen conjugate-gradient backend.
struct ConjugateGradientOptions
{
    double relative_tolerance{1.0e-10};
    Index maximum_iterations{1000};
};

/// Outcome of one iterative linear solve.
struct LinearSolveResult
{
    bool converged{};
    Index iteration_count{};
    double estimated_relative_error{};
};

/// Eigen conjugate-gradient solver for symmetric positive-definite systems.
///
/// `compute_matrix()` copies the finite-volume coefficients into Eigen sparse
/// storage and prepares a diagonal preconditioner. Repeated `solve()` calls
/// reuse both objects and use the supplied solution as the initial guess.
class EigenConjugateGradientSolver
{
  public:
    explicit EigenConjugateGradientSolver(ConjugateGradientOptions options = {});

    EigenConjugateGradientSolver(const EigenConjugateGradientSolver &) = delete;
    EigenConjugateGradientSolver &operator=(const EigenConjugateGradientSolver &) = delete;
    EigenConjugateGradientSolver(EigenConjugateGradientSolver &&) = delete;
    EigenConjugateGradientSolver &operator=(EigenConjugateGradientSolver &&) = delete;

    ~EigenConjugateGradientSolver() = default;

    /// Converts and prepares a matrix for subsequent solves.
    ///
    /// Directed off-diagonal coefficients are copied without symmetrization.
    /// The caller is responsible for satisfying the symmetric
    /// positive-definite precondition of conjugate gradient.
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
    using Solver =
        Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower | Eigen::Upper, Eigen::DiagonalPreconditioner<double>>;

    ConjugateGradientOptions options_;
    SparseMatrix matrix_;
    Solver solver_;
    Index matrix_size_{};
    bool matrix_is_prepared_{};
};

} // namespace cfd
