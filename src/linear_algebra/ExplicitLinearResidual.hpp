#pragma once

#include <Eigen/SparseCore>

#include <cmath>
#include <limits>
#include <span>
#include <stdexcept>

namespace cfd::detail
{

// Scaled sum of squares: never form the square of an unscaled component.
struct ScaledL2Norm
{
    long double scale{};
    long double squared_sum{1.0L};

    void add(const long double value)
    {
        if (!std::isfinite(value))
        {
            throw std::runtime_error("Explicit linear residual contains a non-finite value.");
        }
        const long double magnitude{std::abs(value)};
        if (magnitude == 0.0L)
        {
            return;
        }
        if (magnitude > scale)
        {
            const long double ratio{scale / magnitude};
            squared_sum = 1.0L + squared_sum * ratio * ratio;
            scale = magnitude;
        }
        else
        {
            const long double ratio{magnitude / scale};
            squared_sum += ratio * ratio;
        }
    }
};

// The workspace is sized during matrix preparation, never during solve().
inline long double explicit_relative_residual(const Eigen::SparseMatrix<double> &matrix,
                                              const std::span<const double> rhs, const std::span<const double> solution,
                                              Eigen::VectorXd &workspace)
{
    const Eigen::Map<const Eigen::VectorXd> solution_map{solution.data(), matrix.cols()};
    if (!solution_map.allFinite())
    {
        throw std::runtime_error("Explicit linear residual requires a finite solution.");
    }
    workspace.noalias() = matrix * solution_map;

    ScaledL2Norm residual_norm;
    ScaledL2Norm rhs_norm;
    for (Eigen::Index index = 0; index < matrix.rows(); ++index)
    {
        const long double rhs_value{rhs[static_cast<std::size_t>(index)]};
        residual_norm.add(rhs_value - static_cast<long double>(workspace[index]));
        rhs_norm.add(rhs_value);
    }
    if (residual_norm.scale == 0.0L)
    {
        return 0.0;
    }
    if (rhs_norm.scale == 0.0L)
    {
        throw std::runtime_error("Explicit linear residual is nonzero for a zero RHS.");
    }

    // Divide norms in mantissa/exponent form; even the individual norms or
    // their scale ratio need not fit in double.
    int residual_exponent{};
    int rhs_exponent{};
    const long double residual_fraction{std::frexp(residual_norm.scale, &residual_exponent)};
    const long double rhs_fraction{std::frexp(rhs_norm.scale, &rhs_exponent)};
    const long double relative_residual{
        std::scalbn((residual_fraction / rhs_fraction) * std::sqrt(residual_norm.squared_sum / rhs_norm.squared_sum),
                    residual_exponent - rhs_exponent)};
    if (!std::isfinite(relative_residual) || relative_residual > std::numeric_limits<double>::max())
    {
        throw std::runtime_error("Explicit relative linear residual is outside the finite double range.");
    }
    return relative_residual;
}

} // namespace cfd::detail
