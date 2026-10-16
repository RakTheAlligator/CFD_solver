#include "cfd/compressible/CaloricallyPerfectGas.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace cfd
{
namespace
{

void require_positive(const double value, const char *name)
{
    if (!std::isfinite(value) || !(value > 0.0))
        throw std::runtime_error(std::string{"Perfect gas: "} + name + " must be finite and strictly positive.");
}

double checked_value(const long double value, const char *name)
{
    if (!std::isfinite(value) || std::abs(value) > std::numeric_limits<double>::max())
        throw std::runtime_error(std::string{"Perfect gas: "} + name + " is not representable as finite double.");
    const double result{static_cast<double>(value)};
    if (value != 0.0L && result == 0.0)
        throw std::runtime_error(std::string{"Perfect gas: "} + name + " underflows to zero.");
    return result;
}

void require_finite_velocities(const PrimitiveState &state)
{
    if (!std::isfinite(state.velocity_x) || !std::isfinite(state.velocity_r))
        throw std::runtime_error("Perfect gas: velocities must be finite.");
}

} // namespace

CaloricallyPerfectGas::CaloricallyPerfectGas(const double gamma, const double gas_constant)
    : gamma_(gamma), gas_constant_(gas_constant)
{
    if (!std::isfinite(gamma) || !(gamma > 1.0) || !std::isfinite(gas_constant) || !(gas_constant > 0.0))
        throw std::invalid_argument("Perfect gas requires finite gamma > 1 and gas constant R > 0.");
}

ThermodynamicProperties CaloricallyPerfectGas::properties(const PrimitiveState &state) const
{
    require_positive(state.density, "density");
    require_positive(state.pressure, "pressure");
    require_finite_velocities(state);
    const long double rho{state.density}, pressure{state.pressure}, gamma{gamma_}, gas_constant{gas_constant_};
    const ThermodynamicProperties result{checked_value(pressure / ((gamma - 1.0L) * rho), "specific internal energy"),
                                         checked_value(pressure / (rho * gas_constant), "temperature"),
                                         checked_value(std::sqrt(gamma * pressure / rho), "sound speed")};
    require_positive(result.specific_internal_energy, "specific internal energy");
    require_positive(result.temperature, "temperature");
    require_positive(result.sound_speed, "sound speed");
    return result;
}

EvaluatedGasState CaloricallyPerfectGas::evaluate(const ConservativeState &state) const
{
    require_positive(state.density, "density");
    if (!std::isfinite(state.momentum_x) || !std::isfinite(state.momentum_r) ||
        !std::isfinite(state.total_energy_density))
        throw std::runtime_error("Perfect gas: conservative momenta and total energy density must be finite.");
    const long double rho{state.density}, mx{state.momentum_x}, mr{state.momentum_r};
    const long double ux{mx / rho}, ur{mr / rho};
    const double velocity_x{checked_value(ux, "x velocity")}, velocity_r{checked_value(ur, "radial velocity")};
    // Multiplying momentum by velocity avoids squaring momenta in double.
    const long double kinetic{0.5L * (mx * ux + mr * ur)};
    const long double internal{static_cast<long double>(state.total_energy_density) - kinetic};
    if (!std::isfinite(internal) || !(internal > 0.0L))
        throw std::runtime_error("Perfect gas: internal energy density must be finite and strictly positive.");
    const double pressure{checked_value((static_cast<long double>(gamma_) - 1.0L) * internal, "pressure")};
    const PrimitiveState primitive{state.density, velocity_x, velocity_r, pressure};
    return {primitive, properties(primitive)};
}

ConservativeState CaloricallyPerfectGas::to_conservative(const PrimitiveState &state) const
{
    // Validate all derived thermodynamic properties, even though only p is
    // needed below. No primitive/property cache is kept in the field.
    static_cast<void>(properties(state));
    const long double rho{state.density}, ux{state.velocity_x}, ur{state.velocity_r};
    const long double internal{static_cast<long double>(state.pressure) / (static_cast<long double>(gamma_) - 1.0L)};
    const long double kinetic{0.5L * rho * (ux * ux + ur * ur)};
    const ConservativeState result{state.density, checked_value(rho * ux, "x momentum"),
                                   checked_value(rho * ur, "radial momentum"),
                                   checked_value(internal + kinetic, "total energy density")};
    // Rounding rhoE/momenta to double must not destroy admissibility.
    static_cast<void>(evaluate(result));
    return result;
}

} // namespace cfd
