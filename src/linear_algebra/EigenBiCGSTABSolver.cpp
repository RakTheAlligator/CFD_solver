#include "cfd/linear_algebra/EigenBiCGSTABSolver.hpp"

#include "cfd/linear_algebra/ScalarLinearSystem.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"

#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfd
{

namespace
{

[[nodiscard]]
bool spans_overlap(const std::span<const double> first, const std::span<double> second) noexcept
{
    if (first.empty() || second.empty())
    {
        return false;
    }

    const double *const first_begin{first.data()};
    const double *const first_end{first.data() + first.size()};
    const double *const second_begin{second.data()};
    const double *const second_end{second.data() + second.size()};
    const std::less<const double *> precedes;
    return precedes(first_begin, second_end) && precedes(second_begin, first_end);
}

void require_finite(const std::span<const double> values, const char *const message)
{
    for (const double value : values)
    {
        if (!std::isfinite(value))
        {
            throw std::invalid_argument(message);
        }
    }
}

[[nodiscard]]
Eigen::Index to_eigen_index(const Index value)
{
    constexpr auto maximum_eigen_index{static_cast<Index>(std::numeric_limits<Eigen::Index>::max())};
    if (value > maximum_eigen_index)
    {
        throw std::invalid_argument("Scalar linear system cardinality exceeds the Eigen index range.");
    }
    return static_cast<Eigen::Index>(value);
}

} // namespace

EigenBiCGSTABSolver::EigenBiCGSTABSolver(const BiCGSTABOptions options) : options_(options)
{
    if (!std::isfinite(options_.relative_tolerance) || !(options_.relative_tolerance > 0.0) ||
        !(options_.relative_tolerance < 1.0))
    {
        throw std::invalid_argument("BiCGSTAB relative tolerance must be finite and in (0, 1).");
    }
    if (options_.maximum_iterations == 0 ||
        options_.maximum_iterations > static_cast<Index>(std::numeric_limits<Eigen::Index>::max()))
    {
        throw std::invalid_argument("BiCGSTAB maximum iterations must fit a positive Eigen index.");
    }

    solver_.setTolerance(options_.relative_tolerance);
    solver_.setMaxIterations(static_cast<Eigen::Index>(options_.maximum_iterations));
}

void EigenBiCGSTABSolver::compute_matrix(const ScalarLinearSystem &system)
{
    matrix_is_prepared_ = false;
    require_finite(system.diagonal(), "Scalar linear system diagonal must contain only finite values.");
    require_finite(system.owner_neighbor_coefficients(),
                   "Scalar linear system owner-neighbor coefficients must contain only finite values.");
    require_finite(system.neighbor_owner_coefficients(),
                   "Scalar linear system neighbor-owner coefficients must contain only finite values.");

    const Eigen::Index cell_count{to_eigen_index(system.cell_count())};
    Index internal_face_count{};
    for (const FaceAdjacency &adjacency : system.mesh().face_adjacencies())
    {
        internal_face_count += static_cast<Index>(!adjacency.is_boundary());
    }

    std::vector<Eigen::Triplet<double>> entries;
    entries.reserve(system.cell_count() + 2 * internal_face_count);
    for (Index cell_id = 0; cell_id < system.cell_count(); ++cell_id)
    {
        entries.emplace_back(to_eigen_index(cell_id), to_eigen_index(cell_id), system.diagonal()[cell_id]);
    }

    const auto face_adjacencies{system.mesh().face_adjacencies()};
    for (Index face_id = 0; face_id < system.face_count(); ++face_id)
    {
        const FaceAdjacency &adjacency{face_adjacencies[face_id]};
        if (adjacency.is_boundary())
        {
            continue;
        }

        entries.emplace_back(to_eigen_index(adjacency.owner), to_eigen_index(adjacency.neighbor),
                             system.owner_neighbor_coefficients()[face_id]);
        entries.emplace_back(to_eigen_index(adjacency.neighbor), to_eigen_index(adjacency.owner),
                             system.neighbor_owner_coefficients()[face_id]);
    }

    matrix_.resize(cell_count, cell_count);
    matrix_.setFromTriplets(entries.begin(), entries.end());
    matrix_.makeCompressed();
    solver_.compute(matrix_);
    if (solver_.info() != Eigen::Success)
    {
        throw std::runtime_error("Eigen BiCGSTAB matrix preparation failed.");
    }

    matrix_size_ = system.cell_count();
    matrix_is_prepared_ = true;
}

LinearSolveResult EigenBiCGSTABSolver::solve(const std::span<const double> rhs, const std::span<double> solution)
{
    if (!matrix_is_prepared_)
    {
        throw std::logic_error("Eigen BiCGSTAB solve requires a prepared matrix.");
    }
    if (rhs.size() != matrix_size_ || solution.size() != matrix_size_)
    {
        throw std::invalid_argument("Eigen BiCGSTAB vectors must match the prepared matrix size.");
    }
    if (spans_overlap(rhs, solution))
    {
        throw std::invalid_argument("Eigen BiCGSTAB RHS and solution spans must not overlap.");
    }
    require_finite(rhs, "Eigen BiCGSTAB RHS must contain only finite values.");
    require_finite(solution, "Eigen BiCGSTAB initial guess must contain only finite values.");

    const Eigen::Index size{to_eigen_index(matrix_size_)};
    const Eigen::Map<const Eigen::VectorXd> rhs_map{rhs.data(), size};
    Eigen::Map<Eigen::VectorXd> solution_map{solution.data(), size};
    solution_map = solver_.solveWithGuess(rhs_map, solution_map);

    const Eigen::ComputationInfo information{solver_.info()};
    const double estimated_relative_error{solver_.error()};
    if ((information != Eigen::Success && information != Eigen::NoConvergence) ||
        !std::isfinite(estimated_relative_error))
    {
        throw std::runtime_error("Eigen BiCGSTAB solve failed numerically.");
    }
    require_finite(solution, "Eigen BiCGSTAB produced a non-finite solution.");

    return {
        information == Eigen::Success,
        static_cast<Index>(solver_.iterations()),
        estimated_relative_error,
    };
}

} // namespace cfd
