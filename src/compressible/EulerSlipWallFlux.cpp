#include "cfd/compressible/EulerSlipWallFlux.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace cfd
{
namespace
{

double checked_value(const long double value, const char *quantity)
{
    if (!std::isfinite(value) || std::abs(value) > std::numeric_limits<double>::max())
        throw std::runtime_error(std::string{"Euler slip wall: "} + quantity +
                                 " is not representable as finite double.");
    const double result{static_cast<double>(value)};
    if (value != 0.0L && result == 0.0)
        throw std::runtime_error(std::string{"Euler slip wall: "} + quantity + " underflows to zero.");
    return result;
}

double checked_wave_bound(const long double value)
{
    double result{checked_value(value, "wave bound")};
    // Outward rounding keeps a representable bound conservative after narrowing.
    if (static_cast<long double>(result) < value)
        result = std::nextafter(result, std::numeric_limits<double>::infinity());
    if (!std::isfinite(result) || !(result > 0.0))
        throw std::runtime_error("Euler slip wall requires a finite positive wave bound.");
    return result;
}

} // namespace

EulerSlipWallFluxResult compute_euler_slip_wall_flux(const ConservativeState &inside, const Vector2 &unit_normal,
                                                     const CaloricallyPerfectGas &gas)
{
    constexpr double tolerance{64.0 * std::numeric_limits<double>::epsilon()};
    const double length{std::hypot(unit_normal.x, unit_normal.y)};
    if (!std::isfinite(length) || std::abs(length - 1.0) > tolerance)
        throw std::invalid_argument("Euler slip wall requires a finite unit normal.");
    const auto evaluated{gas.evaluate(inside)};
    const auto primitive{evaluated.primitive};
    const long double gamma{gas.gamma()}, rho{primitive.density}, pressure{primitive.pressure};
    const long double sound{evaluated.thermodynamics.sound_speed};
    const long double un{static_cast<long double>(primitive.velocity_x) * unit_normal.x +
                         static_cast<long double>(primitive.velocity_r) * unit_normal.y};
    long double wall_pressure{pressure};
    long double wave{std::abs(un) + sound};
    if (un > 0.0L)
    {
        // f_shock(p_wall)=un. Eliminating the pressure ratio yields
        // c_shock=k*un+hypot(a,k*un), delta_p=rho*un*c_shock, k=(gamma+1)/4.
        // Unlike a Mach-number quadratic, this avoids squaring un/a.
        const long double scaled_velocity{((gamma + 1.0L) / 4.0L) * un};
        const long double relative_shock_speed{scaled_velocity + std::hypot(sound, scaled_velocity)};
        wall_pressure += rho * un * relative_shock_speed;
        const long double sound_bound{std::sqrt(gamma) * (std::sqrt(wall_pressure) / std::sqrt(rho))};
        wave = std::max({wave, relative_shock_speed, sound_bound});
    }
    else if (un < 0.0L)
    {
        const long double drop{((gamma - 1.0L) / 2.0L) * (-un / sound)};
        if (drop >= 1.0L)
            return {{}, 0.0, checked_wave_bound(wave)}; // Exact wall vacuum, not a pressure floor.
        const long double exponent{2.0L * gamma / (gamma - 1.0L)};
        const long double log_ratio{exponent * std::log1p(-drop)};
        const long double ratio{std::exp(log_ratio)};
        // Multiplication preserves p exactly in the vanishing-expansion limit.
        // Use logs if a ratio underflows but multiplication by p could rescue it.
        wall_pressure = ratio > 0.0L ? pressure * ratio : std::exp(std::log(pressure) + log_ratio);
        if (!(wall_pressure > 0.0L))
            throw std::runtime_error("Euler slip wall: non-vacuum pressure underflows to zero.");
    }
    const double stored_pressure{checked_value(wall_pressure, "wall pressure")};
    return {{0.0, checked_value(wall_pressure * unit_normal.x, "axial momentum flux"),
             checked_value(wall_pressure * unit_normal.y, "radial momentum flux"), 0.0},
            stored_pressure,
            checked_wave_bound(wave)};
}

} // namespace cfd
