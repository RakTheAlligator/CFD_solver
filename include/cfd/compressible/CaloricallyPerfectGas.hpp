#pragma once

#include "cfd/compressible/ConservativeField.hpp"

namespace cfd
{

/// Independent primitive variables; velocities are Cartesian x and radial r.
struct PrimitiveState
{
    double density{};    ///< [kg/m^3].
    double velocity_x{}; ///< [m/s].
    double velocity_r{}; ///< [m/s].
    double pressure{};   ///< [Pa].
};

/// Local derived properties, never stored redundantly in ConservativeField.
struct ThermodynamicProperties
{
    double specific_internal_energy{}; ///< e [J/kg], excluding kinetic energy.
    double temperature{};              ///< [K].
    double sound_speed{};              ///< [m/s].
};

/// Conservative-state evaluation returns primitives and properties together.
struct EvaluatedGasState
{
    PrimitiveState primitive;
    ThermodynamicProperties thermodynamics;
};

/// Concrete, temporary non-reacting calorically perfect gas closure.
///
/// Uses constant gamma and specific gas constant R [J/(kg*K)]. Future
/// temperature/composition-dependent thermodynamics is intentionally not modeled.
/// All successful evaluations have finite values and strictly positive rho,
/// p, e, T and a. Invalid states are rejected, never clipped or repaired.
///
/// Products/ratios use long double intermediates (additional range/precision
/// depends on the platform). Results must be representable as double, including
/// nonzero values without underflow to zero. Kinetic/internal energy separation
/// remains ill-conditioned at very high Mach number: double rhoE cannot encode
/// arbitrarily small internal energy alongside large kinetic energy. Conversion
/// rejects a non-positive recovered internal energy, but cannot guarantee small
/// relative pressure error when only a few significant internal-energy bits remain.
/// Subnormal results likewise have reduced relative precision.
/// Valid local operations allocate nothing and return small values, not fields.
class CaloricallyPerfectGas
{
  public:
    /// @throws std::invalid_argument Unless gamma > 1 and R > 0, both finite.
    CaloricallyPerfectGas(double gamma, double gas_constant);

    [[nodiscard]] double gamma() const noexcept
    {
        return gamma_;
    }

    [[nodiscard]] double gas_constant() const noexcept
    {
        return gas_constant_;
    }

    /// Converts U to primitives and computes e, T and a in a single evaluation.
    /// p = (gamma-1)*(rhoE - |rho*u|^2/(2*rho)).
    /// @throws std::runtime_error For invalid/non-representable states or properties.
    [[nodiscard]] EvaluatedGasState evaluate(const ConservativeState &state) const;

    /// Converts primitives to [rho, rho*ux, rho*ur, rhoE], checking that the
    /// stored total energy still yields positive internal energy on recovery.
    /// @throws std::runtime_error For invalid/non-representable states or properties.
    [[nodiscard]] ConservativeState to_conservative(const PrimitiveState &state) const;

    /// e = p/((gamma-1)*rho), T = p/(rho*R), a = sqrt(gamma*p/rho).
    /// @throws std::runtime_error For invalid/non-representable states or properties.
    [[nodiscard]] ThermodynamicProperties properties(const PrimitiveState &state) const;

  private:
    double gamma_;
    double gas_constant_;
};

} // namespace cfd
