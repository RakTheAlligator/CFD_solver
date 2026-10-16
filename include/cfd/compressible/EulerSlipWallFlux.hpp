#pragma once

#include "cfd/compressible/EulerFlux.hpp"

namespace cfd
{

struct EulerSlipWallFluxResult
{
    EulerFlux flux;              ///< Exactly [0, p_wall*nx, p_wall*nr, 0], per unit area.
    double wall_pressure{};      ///< Positive except for a physical wall vacuum (exactly zero).
    double maximum_wave_speed{}; ///< Conservative bound for the normal wall Riemann problem [m/s].
};

/// Exact stationary slip-wall normal Riemann flux for a calorically perfect gas.
///
/// The finite unit normal points out of the gas; un>0 is compression toward
/// the wall. Tangential velocity does not enter the normal Riemann problem.
/// For un<0, p_wall = p*[1+(gamma-1)*un/(2*a)]^(2*gamma/(gamma-1)).
/// If -un >= 2*a/(gamma-1), the exact expansion leaves vacuum at the wall:
/// p_wall=0 and every flux component is zero. No vacuum cell state is created.
/// Compression uses the closed-form Rankine-Hugoniot shock relation, not Rusanov.
///
/// The CFL bound includes |un|+a. In compression it also bounds the shock
/// speed and star-state sound speed using the shock-relative speed and
/// sqrt(gamma*p_wall/rho_inside), since rho_star >= rho_inside.
/// Wider intermediates and log1p avoid unnecessary range/cancellation losses;
/// additional long-double range/precision is platform dependent. No allocation
/// on valid calls, no state clipping, and no reflected ghost-state evaluation.
/// @throws std::invalid_argument If the normal violates the EulerFlux unit-normal contract.
/// @throws std::runtime_error For invalid gas states, non-representable pressure,
///         momentum flux or wave bound, including non-vacuum pressure underflow.
[[nodiscard]] EulerSlipWallFluxResult compute_euler_slip_wall_flux(const ConservativeState &inside,
                                                                   const Vector2 &unit_normal,
                                                                   const CaloricallyPerfectGas &gas);

} // namespace cfd
