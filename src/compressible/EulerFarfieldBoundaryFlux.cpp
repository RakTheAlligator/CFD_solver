#include "cfd/compressible/EulerFarfieldBoundaryFlux.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cfd
{
namespace
{
double checked(const long double value)
{
    if (!std::isfinite(value) || std::abs(value) > std::numeric_limits<double>::max())
        throw std::runtime_error("Euler farfield: non-representable arithmetic.");
    const double stored{static_cast<double>(value)};
    if (value != 0.0L && stored == 0.0)
        throw std::runtime_error("Euler farfield: arithmetic underflows to zero.");
    return stored;
}

struct NormalState
{
    long double rho, un, ut, pressure, sound;
};

NormalState normal_state(const EvaluatedGasState &state, const Vector2 n)
{
    const auto p{state.primitive};
    return {p.density, static_cast<long double>(p.velocity_x) * n.x + static_cast<long double>(p.velocity_r) * n.y,
            -static_cast<long double>(p.velocity_x) * n.y + static_cast<long double>(p.velocity_r) * n.x, p.pressure,
            state.thermodynamics.sound_speed};
}

long double wave_curve(const NormalState &state, const long double log_pressure, const long double gamma)
{
    const long double log_ratio{log_pressure - std::log(state.pressure)};
    if (log_ratio <= 0.0L)
        return 2.0L * state.sound / (gamma - 1.0L) * std::expm1((gamma - 1.0L) / (2.0L * gamma) * log_ratio);
    const long double pressure{std::exp(log_pressure)};
    return (pressure - state.pressure) * std::sqrt(2.0L / ((gamma + 1.0L) * state.rho) /
                                                   (pressure + (gamma - 1.0L) / (gamma + 1.0L) * state.pressure));
}

long double star_pressure(const NormalState &left, const NormalState &right, const long double gamma)
{
    if (right.un - left.un >= 2.0L * (left.sound + right.sound) / (gamma - 1.0L))
        throw std::runtime_error("Euler farfield: vacuum-generating normal Riemann data are unsupported.");
    const auto function = [&](const long double log_pressure) {
        return wave_curve(left, log_pressure, gamma) + wave_curve(right, log_pressure, gamma) + right.un - left.un;
    };
    long double lower{std::log(static_cast<long double>(std::numeric_limits<double>::denorm_min()))};
    long double upper{std::log(std::max(left.pressure, right.pressure))};
    const long double maximum{std::log(static_cast<long double>(std::numeric_limits<double>::max()))};
    const long double lower_value{function(lower)};
    if (lower_value == 0.0L)
        return std::exp(lower);
    if (!(lower_value < 0.0L))
        throw std::runtime_error("Euler farfield: positive star pressure is below the representable range.");
    while (function(upper) < 0.0L)
    {
        if (upper == maximum)
            throw std::runtime_error("Euler farfield: star pressure exceeds the representable range.");
        upper = std::min(maximum, upper + std::log(2.0L));
    }
    if (function(upper) == 0.0L)
        return std::exp(upper);
    // Log-pressure bisection resolves even near-vacuum positive roots without a
    // pressure floor. Stop on representational stagnation, not an empirical epsilon.
    for (int iteration = 0; iteration < 128; ++iteration)
    {
        const double low_pressure{checked(std::exp(lower))}, high_pressure{checked(std::exp(upper))};
        if (std::nextafter(low_pressure, std::numeric_limits<double>::infinity()) >= high_pressure)
            return std::abs(function(lower)) < std::abs(function(upper)) ? low_pressure : high_pressure;
        const long double middle{lower + (upper - lower) / 2.0L};
        if (middle == lower || middle == upper)
            return std::exp(middle);
        const long double value{function(middle)};
        if (!std::isfinite(value))
            throw std::runtime_error("Euler farfield: invalid pressure wave curve.");
        if (value == 0.0L)
            return std::exp(middle);
        if (value > 0.0L)
            upper = middle;
        else
            lower = middle;
    }
    throw std::runtime_error("Euler farfield: pressure root did not resolve.");
}

NormalState star_state(const NormalState &side, const long double pressure, const long double un,
                       const long double gamma)
{
    const long double ratio{pressure / side.pressure}, c{(gamma - 1.0L) / (gamma + 1.0L)};
    const long double rho{side.rho * (ratio > 1.0L ? (ratio + c) / (c * ratio + 1.0L) : std::pow(ratio, 1.0L / gamma))};
    if (!std::isfinite(rho) || !(rho > 0.0L))
        throw std::runtime_error("Euler farfield: invalid star density.");
    return {rho, un, side.ut, pressure, std::sqrt(gamma * pressure / rho)};
}

long double shock_speed(const NormalState &side, const long double pressure, const long double sign,
                        const long double gamma)
{
    return side.un +
           sign * side.sound *
               std::sqrt((gamma + 1.0L) / (2.0L * gamma) * pressure / side.pressure + (gamma - 1.0L) / (2.0L * gamma));
}

NormalState sample(const NormalState &side, const NormalState &star, const long double sign, const long double gamma)
{
    if (star.pressure > side.pressure)
        return sign * shock_speed(side, star.pressure, sign, gamma) <= 0.0L ? side : star;
    if (sign * (side.un + sign * side.sound) <= 0.0L)
        return side;
    if (sign * (star.un + sign * star.sound) >= 0.0L)
        return star;
    // The boundary lies inside a rarefaction: its local characteristic speed is zero.
    const long double sound{2.0L / (gamma + 1.0L) * (side.sound - sign * (gamma - 1.0L) * side.un / 2.0L)};
    const long double ratio{sound / side.sound};
    return {side.rho * std::pow(ratio, 2.0L / (gamma - 1.0L)), -sign * sound, side.ut,
            side.pressure * std::pow(ratio, 2.0L * gamma / (gamma - 1.0L)), sound};
}

EulerFarfieldRegime regime(const double mach)
{
    if (mach <= -1.0)
        return EulerFarfieldRegime::SupersonicInflow;
    if (mach < 0.0)
        return EulerFarfieldRegime::SubsonicInflow;
    if (mach < 1.0)
        return EulerFarfieldRegime::SubsonicOutflow;
    return EulerFarfieldRegime::SupersonicOutflow;
}
} // namespace

EulerFarfieldBoundaryFluxResult compute_euler_farfield_boundary_flux(const ConservativeState &inside,
                                                                     const Vector2 &outward_unit_normal,
                                                                     const EulerBoundaryCondition &condition,
                                                                     const CaloricallyPerfectGas &gas)
{
    const Vector2 n{outward_unit_normal};
    if (condition.type() != EulerBoundaryConditionType::Farfield)
        throw std::invalid_argument("Euler farfield flux requires a Farfield condition.");
    if (!std::isfinite(n.x) || !std::isfinite(n.y) ||
        std::abs(std::hypot(n.x, n.y) - 1.0) > 64.0 * std::numeric_limits<double>::epsilon())
        throw std::invalid_argument("Euler farfield requires a finite unit normal.");
    const auto left{normal_state(gas.evaluate(inside), n)};
    const auto outside{gas.to_conservative(condition.prescribed_state())};
    const auto right{normal_state(gas.evaluate(outside), n)};
    const long double gamma{gas.gamma()};
    const long double pressure{star_pressure(left, right, gamma)};
    const long double un{(left.un + right.un + wave_curve(right, std::log(pressure), gamma) -
                          wave_curve(left, std::log(pressure), gamma)) /
                         2.0L};
    const auto left_star{star_state(left, pressure, un, gamma)}, right_star{star_state(right, pressure, un, gamma)};
    const auto trace{un >= 0.0L ? sample(left, left_star, -1.0L, gamma) : sample(right, right_star, 1.0L, gamma)};
    const auto boundary{gas.to_conservative({checked(trace.rho), checked(trace.un * n.x - trace.ut * n.y),
                                             checked(trace.un * n.y + trace.ut * n.x), checked(trace.pressure)})};
    const auto evaluated{gas.evaluate(boundary)};
    const double mach{checked((static_cast<long double>(evaluated.primitive.velocity_x) * n.x +
                               static_cast<long double>(evaluated.primitive.velocity_r) * n.y) /
                              evaluated.thermodynamics.sound_speed)};
    long double wave{std::max({std::abs(left.un) + left.sound, std::abs(right.un) + right.sound,
                               std::abs(un) + left_star.sound, std::abs(un) + right_star.sound})};
    wave = std::max(wave, std::abs(static_cast<long double>(evaluated.primitive.velocity_x) * n.x +
                                   static_cast<long double>(evaluated.primitive.velocity_r) * n.y) +
                              evaluated.thermodynamics.sound_speed);
    if (pressure > left.pressure)
        wave = std::max(wave, std::abs(shock_speed(left, pressure, -1.0L, gamma)));
    if (pressure > right.pressure)
        wave = std::max(wave, std::abs(shock_speed(right, pressure, 1.0L, gamma)));
    double stored_wave{checked(wave)};
    if (static_cast<long double>(stored_wave) < wave)
        stored_wave = std::nextafter(stored_wave, std::numeric_limits<double>::infinity());
    if (!std::isfinite(stored_wave) || !(stored_wave > 0.0))
        throw std::runtime_error("Euler farfield requires a finite positive wave bound.");
    return {boundary, compute_physical_euler_flux(boundary, n, gas), stored_wave, mach, regime(mach)};
}
} // namespace cfd
