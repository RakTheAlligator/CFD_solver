#pragma once

#include "cfd/mesh/Types.hpp"

namespace cfd
{

/// Outcome of one iterative linear solve.
struct LinearSolveResult
{
    bool converged{};
    Index iteration_count{};
    double estimated_relative_error{};
};

} // namespace cfd
