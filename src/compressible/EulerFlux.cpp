#include "cfd/compressible/EulerFlux.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cfd
{
namespace
{

void validate_normal(const Vector2 &normal)
{
    constexpr double tolerance{64.0 * std::numeric_limits<double>::epsilon()};
    const double length{std::hypot(normal.x, normal.y)};
    if (!std::isfinite(normal.x) || !std::isfinite(normal.y) || !std::isfinite(length) ||
        std::abs(length - 1.0) > tolerance)
        throw std::invalid_argument("Euler flux requires a finite unit normal.");
}

double checked_component(const long double value)
{
    if (!std::isfinite(value) || std::abs(value) > std::numeric_limits<double>::max())
        throw std::runtime_error("Euler flux or wave speed is not representable as finite double.");
    const double result{static_cast<double>(value)};
    if (value != 0.0L && result == 0.0)
        throw std::runtime_error("Euler flux or wave speed underflows to zero.");
    return result;
}

long double normal_velocity(const PrimitiveState &state, const Vector2 &normal)
{
    return static_cast<long double>(state.velocity_x) * normal.x +
           static_cast<long double>(state.velocity_r) * normal.y;
}

EulerFlux physical_flux(const ConservativeState &state, const PrimitiveState &primitive, const Vector2 &normal,
                        const long double un)
{
    const long double pressure{primitive.pressure};
    return {checked_component(static_cast<long double>(state.density) * un),
            checked_component(static_cast<long double>(state.momentum_x) * un + pressure * normal.x),
            checked_component(static_cast<long double>(state.momentum_r) * un + pressure * normal.y),
            checked_component((static_cast<long double>(state.total_energy_density) + pressure) * un)};
}

double rusanov_component(const double flux_left, const double flux_right, const double state_left,
                         const double state_right, const double speed)
{
    // Widen before subtraction/addition: finite operands can have an
    // unrepresentable double difference or sum while the final flux fits.
    return checked_component(0.5L * flux_left + 0.5L * flux_right -
                             0.5L * speed * (static_cast<long double>(state_right) - state_left));
}

} // namespace

EulerFlux compute_physical_euler_flux(const ConservativeState &state, const Vector2 &unit_normal,
                                      const CaloricallyPerfectGas &gas)
{
    validate_normal(unit_normal);
    const auto evaluated{gas.evaluate(state)};
    return physical_flux(state, evaluated.primitive, unit_normal, normal_velocity(evaluated.primitive, unit_normal));
}

RusanovFluxResult compute_rusanov_flux(const ConservativeState &left, const ConservativeState &right,
                                       const Vector2 &unit_normal, const CaloricallyPerfectGas &gas)
{
    validate_normal(unit_normal);
    const auto evaluated_left{gas.evaluate(left)}, evaluated_right{gas.evaluate(right)};
    const long double un_left{normal_velocity(evaluated_left.primitive, unit_normal)},
        un_right{normal_velocity(evaluated_right.primitive, unit_normal)};
    const double speed{checked_component(std::max(std::abs(un_left) + evaluated_left.thermodynamics.sound_speed,
                                                  std::abs(un_right) + evaluated_right.thermodynamics.sound_speed))};
    const auto flux_left{physical_flux(left, evaluated_left.primitive, unit_normal, un_left)},
        flux_right{physical_flux(right, evaluated_right.primitive, unit_normal, un_right)};
    return {{rusanov_component(flux_left.mass, flux_right.mass, left.density, right.density, speed),
             rusanov_component(flux_left.momentum_x, flux_right.momentum_x, left.momentum_x, right.momentum_x, speed),
             rusanov_component(flux_left.momentum_r, flux_right.momentum_r, left.momentum_r, right.momentum_r, speed),
             rusanov_component(flux_left.energy, flux_right.energy, left.total_energy_density,
                               right.total_energy_density, speed)},
            speed};
}

} // namespace cfd
