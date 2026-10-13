#pragma once

#include "cfd/mesh/Types.hpp"

namespace cfd
{

/// Outcome of one iterative linear solve.
struct LinearSolveResult
{
    bool converged{};
    Index iteration_count{};
    /// Explicit `||b - A*x||_2 / ||b||_2` against the prepared matrix, not
    /// Eigen's recursive estimate. Exactly zero RHS/solution gives zero.
    double estimated_relative_error{};
};

} // namespace cfd
