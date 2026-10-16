#include "cfd/compressible/EulerCharacteristicBoundaryFlux.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace cfd
{
namespace
{
double checked(const long double value)
{
    if (!std::isfinite(value) || std::abs(value) > std::numeric_limits<double>::max())
        throw std::runtime_error("Euler characteristic BC: non-representable arithmetic.");
    const double result{static_cast<double>(value)};
    if (value != 0.0L && result == 0.0)
        throw std::runtime_error("Euler characteristic BC: arithmetic underflows to zero.");
    return result;
}

long double dot(const PrimitiveState &state, const Vector2 normal)
{
    return static_cast<long double>(state.velocity_x) * normal.x +
           static_cast<long double>(state.velocity_r) * normal.y;
}

void require_subsonic_inlet(const long double un, const long double sound)
{
    if (un > 0.0L || un <= -sound)
        throw std::runtime_error(
            "Euler subsonic inlet requires -a < un <= 0; reverse flow or sonic inflow is unsupported.");
}

// Strictly decreasing f(q) on the normal-subsonic interval. The sonic upper
// endpoint follows from mu^2*q^2 = a0^2 - (gamma-1)*q^2/2, not a Mach clip.
long double inlet_speed(const long double invariant, const long double mu, const long double a0_squared,
                        const long double gm)
{
    if (!(mu < 0.0L))
        throw std::runtime_error("Euler inlet direction must point into every face (direction.normal < 0).");
    const auto function = [=](const long double speed) {
        return mu * speed + 2.0L * std::sqrt(a0_squared - 0.5L * gm * speed * speed) / gm - invariant;
    };
    long double lower{}, upper{std::sqrt(a0_squared / (mu * mu + 0.5L * gm))};
    if (!std::isfinite(upper) || !(upper > 0.0L) || !(function(lower) >= 0.0L) || !(function(upper) < 0.0L))
        throw std::runtime_error("Euler inlet has no compatible strictly normal-subsonic root (stagnation/choking).");
    if (function(lower) == 0.0L)
        return lower;
    for (int iteration = 0; iteration < 80; ++iteration)
    {
        const long double middle{lower + 0.5L * (upper - lower)};
        if (middle == lower || middle == upper)
            break;
        if (function(middle) > 0.0L)
            lower = middle;
        else
            upper = middle;
    }
    return lower + 0.5L * (upper - lower);
}

ConservativeState inlet_trace(const EvaluatedGasState &evaluated, const long double un, const Vector2 normal,
                              const EulerBoundaryCondition &condition, const CaloricallyPerfectGas &gas,
                              const double group_area)
{
    require_subsonic_inlet(un, evaluated.thermodynamics.sound_speed);
    const long double gamma{gas.gamma()}, gm{gamma - 1.0L}, gas_constant{gas.gas_constant()};
    const auto direction{condition.direction()};
    const long double mu{static_cast<long double>(direction.x) * normal.x +
                         static_cast<long double>(direction.y) * normal.y};
    const long double a0_squared{gamma * gas_constant * condition.total_temperature()};
    const long double invariant{un + 2.0L * evaluated.thermodynamics.sound_speed / gm};
    const long double speed{inlet_speed(invariant, mu, a0_squared, gm)};
    const long double temperature{(a0_squared - 0.5L * gm * speed * speed) / (gamma * gas_constant)};
    long double pressure{}, density{};
    if (condition.type() == EulerBoundaryConditionType::MassFlowInlet)
    {
        if (!std::isfinite(group_area) || !(group_area > 0.0))
            throw std::invalid_argument("Euler mass inlet requires finite positive group revolution area.");
        if (!(speed > 0.0L))
            throw std::runtime_error("Euler mass inlet cannot impose positive mdot at zero inflow speed.");
        density = static_cast<long double>(condition.inward_mass_flow()) / group_area / (-mu * speed);
        pressure = density * gas_constant * temperature;
    }
    else
    {
        pressure = condition.total_pressure() * std::pow(temperature / condition.total_temperature(), gamma / gm);
        density = pressure / (gas_constant * temperature);
    }
    const auto trace{gas.to_conservative(
        {checked(density), checked(speed * direction.x), checked(speed * direction.y), checked(pressure)})};
    const auto boundary{gas.evaluate(trace)};
    require_subsonic_inlet(dot(boundary.primitive, normal), boundary.thermodynamics.sound_speed);
    return trace;
}

ConservativeState outlet_trace(const EvaluatedGasState &evaluated, const long double un, const Vector2 normal,
                               const EulerBoundaryCondition &condition, const CaloricallyPerfectGas &gas)
{
    if (un < 0.0L)
        throw std::runtime_error(
            "Euler pressure outlet backflow requires unavailable incoming thermal/directional data.");
    const auto primitive{evaluated.primitive};
    const long double gamma{gas.gamma()}, gm{gamma - 1.0L}, pressure{condition.static_pressure()};
    const long double density{primitive.density * std::pow(pressure / primitive.pressure, 1.0L / gamma)};
    const long double sound{std::sqrt(gamma * pressure / density)};
    const long double normal_speed{un + 2.0L * (evaluated.thermodynamics.sound_speed - sound) / gm};
    if (!std::isfinite(normal_speed) || normal_speed < 0.0L || !(normal_speed < sound))
        throw std::runtime_error("Euler pressure outlet trace is backflow or sonic/supersonic; incompatible pressure.");
    const auto trace{
        gas.to_conservative({checked(density), checked(primitive.velocity_x + (normal_speed - un) * normal.x),
                             checked(primitive.velocity_r + (normal_speed - un) * normal.y), checked(pressure)})};
    const auto boundary{gas.evaluate(trace)};
    const long double recovered_un{dot(boundary.primitive, normal)};
    if (recovered_un < 0.0L || !(recovered_un < boundary.thermodynamics.sound_speed))
        throw std::runtime_error("Euler pressure outlet rounded trace is outside the subsonic outflow regime.");
    return trace;
}
} // namespace

EulerCharacteristicBoundaryFluxResult compute_euler_characteristic_boundary_flux(
    const ConservativeState &inside, const Vector2 &outward_unit_normal, const EulerBoundaryCondition &condition,
    const CaloricallyPerfectGas &gas, const double group_revolution_area)
{
    constexpr double normal_tolerance{64.0 * std::numeric_limits<double>::epsilon()};
    if (!std::isfinite(outward_unit_normal.x) || !std::isfinite(outward_unit_normal.y) ||
        std::abs(std::hypot(outward_unit_normal.x, outward_unit_normal.y) - 1.0) > normal_tolerance)
        throw std::invalid_argument("Euler characteristic BC requires a finite unit normal.");
    const auto evaluated{gas.evaluate(inside)};
    const long double un{dot(evaluated.primitive, outward_unit_normal)}, sound{evaluated.thermodynamics.sound_speed};
    ConservativeState trace{inside};
    switch (condition.type())
    {
    case EulerBoundaryConditionType::TotalConditionsInlet:
    case EulerBoundaryConditionType::MassFlowInlet:
        trace = inlet_trace(evaluated, un, outward_unit_normal, condition, gas, group_revolution_area);
        break;
    case EulerBoundaryConditionType::SupersonicInlet: {
        trace = gas.to_conservative(condition.prescribed_state());
        const auto prescribed{gas.evaluate(trace)};
        if (!(un < -sound) ||
            !(dot(prescribed.primitive, outward_unit_normal) < -prescribed.thermodynamics.sound_speed))
            throw std::runtime_error(
                "Euler supersonic inlet requires strictly supersonic normal inflow inside and at trace.");
        break;
    }
    case EulerBoundaryConditionType::PressureOutlet:
        if (un < sound)
            trace = outlet_trace(evaluated, un, outward_unit_normal, condition, gas);
        break;
    case EulerBoundaryConditionType::SupersonicOutlet:
        if (!(un >= sound))
            throw std::runtime_error("Euler supersonic outlet requires outward normal velocity un >= a.");
        break;
    case EulerBoundaryConditionType::Axis:
    case EulerBoundaryConditionType::SlipWall:
    case EulerBoundaryConditionType::Transmissive:
        throw std::invalid_argument("Euler characteristic flux requires an open characteristic BC type.");
    }
    const auto boundary{gas.evaluate(trace)};
    const long double un_boundary{dot(boundary.primitive, outward_unit_normal)};
    return {trace, compute_physical_euler_flux(trace, outward_unit_normal, gas),
            checked(std::max(std::abs(un) + sound, std::abs(un_boundary) + boundary.thermodynamics.sound_speed)),
            checked(un_boundary / boundary.thermodynamics.sound_speed)};
}
} // namespace cfd
