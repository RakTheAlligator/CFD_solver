#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>

namespace cfd::verification
{

/// Qualitative behavior of three successively refined scalar solutions.
enum class GridConvergenceBehavior : std::uint8_t
{
    Monotonic,
    Oscillatory,
};

/// Three-grid Richardson extrapolation and fine-grid GCI estimates.
struct GridConvergenceIndexResult
{
    GridConvergenceBehavior behavior{};
    double apparent_order{};
    double extrapolated_value{};
    double approximate_relative_error{};
    double extrapolated_relative_error{};
    double fine_grid_convergence_index{};
};

/// Applies the generalized three-grid procedure of Celik et al.
///
/// Grid 1 is fine, grid 2 medium, and grid 3 coarse. Refinement ratios are
/// `r21 = h2 / h1` and `r32 = h3 / h2` and must both exceed one. Unequal
/// ratios are supported through the fixed-point apparent-order equation.
///
/// @return No result when either solution difference is indistinguishable
///         from zero at floating-point scale, or when the apparent-order
///         iteration is numerically ill-conditioned.
[[nodiscard]]
inline std::optional<GridConvergenceIndexResult> compute_grid_convergence_index(const double fine_value,
                                                                                const double medium_value,
                                                                                const double coarse_value,
                                                                                const double r21, const double r32)
{
    if (!std::isfinite(fine_value) || !std::isfinite(medium_value) || !std::isfinite(coarse_value) ||
        !std::isfinite(r21) || !std::isfinite(r32) || !(r21 > 1.0) || !(r32 > 1.0))
    {
        throw std::invalid_argument("GCI inputs must be finite and refinement ratios must exceed one.");
    }

    const double epsilon21{medium_value - fine_value};
    const double epsilon32{coarse_value - medium_value};
    const double solution_scale{std::max({1.0, std::abs(fine_value), std::abs(medium_value), std::abs(coarse_value)})};
    const double difference_tolerance{128.0 * std::numeric_limits<double>::epsilon() * solution_scale};
    if (std::abs(epsilon21) <= difference_tolerance || std::abs(epsilon32) <= difference_tolerance)
    {
        return std::nullopt;
    }

    const double epsilon_ratio{epsilon32 / epsilon21};
    const double sign{epsilon_ratio > 0.0 ? 1.0 : -1.0};
    const double logarithmic_refinement_ratio{std::log(r21)};
    const double logarithmic_error_ratio{std::log(std::abs(epsilon_ratio))};
    double apparent_order{std::abs(logarithmic_error_ratio) / logarithmic_refinement_ratio};

    constexpr std::size_t maximum_iterations{100};
    constexpr double iteration_tolerance{128.0 * std::numeric_limits<double>::epsilon()};
    bool converged{};
    for (std::size_t iteration = 0; iteration < maximum_iterations; ++iteration)
    {
        const double numerator{std::pow(r21, apparent_order) - sign};
        const double denominator{std::pow(r32, apparent_order) - sign};
        const double correction_ratio{numerator / denominator};
        if (!std::isfinite(correction_ratio) || !(correction_ratio > 0.0))
        {
            return std::nullopt;
        }

        const double next_order{std::abs(logarithmic_error_ratio + std::log(correction_ratio)) /
                                logarithmic_refinement_ratio};
        if (!std::isfinite(next_order) || !(next_order > 0.0))
        {
            return std::nullopt;
        }
        if (std::abs(next_order - apparent_order) <=
            iteration_tolerance * std::max({1.0, std::abs(next_order), std::abs(apparent_order)}))
        {
            apparent_order = next_order;
            converged = true;
            break;
        }
        apparent_order = next_order;
    }
    if (!converged)
    {
        return std::nullopt;
    }

    const double powered_refinement_ratio{std::pow(r21, apparent_order)};
    const double extrapolation_denominator{powered_refinement_ratio - 1.0};
    if (!std::isfinite(powered_refinement_ratio) || !(extrapolation_denominator > 0.0) || fine_value == 0.0)
    {
        return std::nullopt;
    }

    const double extrapolated_value{(powered_refinement_ratio * fine_value - medium_value) / extrapolation_denominator};
    if (!std::isfinite(extrapolated_value) || extrapolated_value == 0.0)
    {
        return std::nullopt;
    }

    const double approximate_relative_error{std::abs((fine_value - medium_value) / fine_value)};
    const double extrapolated_relative_error{std::abs((extrapolated_value - fine_value) / extrapolated_value)};
    const double fine_grid_convergence_index{1.25 * approximate_relative_error / extrapolation_denominator};
    if (!std::isfinite(approximate_relative_error) || !std::isfinite(extrapolated_relative_error) ||
        !std::isfinite(fine_grid_convergence_index))
    {
        return std::nullopt;
    }

    return GridConvergenceIndexResult{
        .behavior = epsilon_ratio > 0.0 ? GridConvergenceBehavior::Monotonic : GridConvergenceBehavior::Oscillatory,
        .apparent_order = apparent_order,
        .extrapolated_value = extrapolated_value,
        .approximate_relative_error = approximate_relative_error,
        .extrapolated_relative_error = extrapolated_relative_error,
        .fine_grid_convergence_index = fine_grid_convergence_index,
    };
}

} // namespace cfd::verification
