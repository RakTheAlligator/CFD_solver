#pragma once

#include "cfd/math/Point2.hpp"
#include "cfd/math/Vector2.hpp"

#include <array>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace cfd::verification::kovasznay
{

// Shared analytical reference for the isolated-momentum and coupled-SIMPLE studies.
constexpr double domain_length{1.0};
constexpr double domain_height{1.0};
constexpr double density{1.0};
constexpr double reynolds_number{40.0};
constexpr double dynamic_viscosity{1.0 / reynolds_number};
constexpr double wave_number{2.0 * std::numbers::pi_v<double>};
inline const double lambda{reynolds_number / 2.0 -
                           std::sqrt(reynolds_number * reynolds_number / 4.0 + wave_number * wave_number)};

[[nodiscard]]
inline double analytical_u(const Point2 &point) noexcept
{
    return 1.0 - std::exp(lambda * point.x) * std::cos(wave_number * point.y);
}

[[nodiscard]]
inline double analytical_v(const Point2 &point) noexcept
{
    return lambda / wave_number * std::exp(lambda * point.x) * std::sin(wave_number * point.y);
}

[[nodiscard]]
inline double analytical_pressure(const Point2 &point) noexcept
{
    return 0.5 * (1.0 - std::exp(2.0 * lambda * point.x));
}

[[nodiscard]]
inline Vector2 analytical_u_gradient(const Point2 &point) noexcept
{
    const double exponential{std::exp(lambda * point.x)};
    return {-lambda * exponential * std::cos(wave_number * point.y),
            wave_number * exponential * std::sin(wave_number * point.y)};
}

[[nodiscard]]
inline Vector2 analytical_v_gradient(const Point2 &point) noexcept
{
    const double exponential{std::exp(lambda * point.x)};
    return {lambda * lambda / wave_number * exponential * std::sin(wave_number * point.y),
            lambda * exponential * std::cos(wave_number * point.y)};
}

[[nodiscard]]
inline Vector2 analytical_pressure_gradient(const Point2 &point) noexcept
{
    return {-lambda * std::exp(2.0 * lambda * point.x), 0.0};
}

[[nodiscard]]
inline double analytical_u_laplacian(const Point2 &point) noexcept
{
    return (wave_number * wave_number - lambda * lambda) * std::exp(lambda * point.x) * std::cos(wave_number * point.y);
}

[[nodiscard]]
inline double analytical_v_laplacian(const Point2 &point) noexcept
{
    return lambda * (lambda * lambda / wave_number - wave_number) * std::exp(lambda * point.x) *
           std::sin(wave_number * point.y);
}

inline void verify_analytical_solution()
{
    constexpr std::array sample_points{Point2{0.15, 0.17}, Point2{0.48, 0.39}, Point2{0.83, 0.71}};
    constexpr double identity_tolerance{1.0e-12};
    constexpr double derivative_step{1.0e-6};
    constexpr double derivative_tolerance{1.0e-8};

    for (const Point2 &point : sample_points)
    {
        const double u{analytical_u(point)};
        const double v{analytical_v(point)};
        const Vector2 u_gradient{analytical_u_gradient(point)};
        const Vector2 v_gradient{analytical_v_gradient(point)};
        const Vector2 pressure_gradient{analytical_pressure_gradient(point)};
        const double continuity_residual{u_gradient.x + v_gradient.y};
        const double u_momentum_residual{density * (u * u_gradient.x + v * u_gradient.y) -
                                         dynamic_viscosity * analytical_u_laplacian(point) + pressure_gradient.x};
        const double v_momentum_residual{density * (u * v_gradient.x + v * v_gradient.y) -
                                         dynamic_viscosity * analytical_v_laplacian(point) + pressure_gradient.y};
        const Point2 point_before{point.x - derivative_step, point.y};
        const Point2 point_after{point.x + derivative_step, point.y};
        const double centered_pressure_derivative{
            (analytical_pressure(point_after) - analytical_pressure(point_before)) / (2.0 * derivative_step)};

        if (std::abs(continuity_residual) > identity_tolerance || std::abs(u_momentum_residual) > identity_tolerance ||
            std::abs(v_momentum_residual) > identity_tolerance ||
            std::abs(centered_pressure_derivative - pressure_gradient.x) > derivative_tolerance)
        {
            throw std::runtime_error("The analytical Kovasznay formulas failed their independent identity check.");
        }
    }
}

} // namespace cfd::verification::kovasznay
