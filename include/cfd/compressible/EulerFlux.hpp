#pragma once

#include "cfd/compressible/CaloricallyPerfectGas.hpp"
#include "cfd/math/Vector2.hpp"

namespace cfd
{

/// Normal Euler flux per unit surface, not a conservative state or integrated flux.
struct EulerFlux
{
    double mass{};       ///< [kg/(m^2*s)].
    double momentum_x{}; ///< [Pa].
    double momentum_r{}; ///< [Pa].
    double energy{};     ///< [W/m^2].
};

struct RusanovFluxResult
{
    EulerFlux flux;
    double maximum_wave_speed{}; ///< max(|un_L|+a_L, |un_R|+a_R) [m/s], used in dissipation.
};

/// Fn = [rho*un, rho*ux*un+p*nx, rho*ur*un+p*nr, (rhoE+p)*un].
///
/// The normal is (nx,nr), with nr stored in Vector2::y. It must be finite and
/// unit length within 64*epsilon; it is never
/// silently normalized. Valid calls allocate nothing and leave inputs unchanged.
/// @throws std::invalid_argument If the normal is invalid.
/// @throws std::runtime_error For invalid gas states or non-representable fluxes
///         (including a nonzero component underflowing to zero).
[[nodiscard]] EulerFlux compute_physical_euler_flux(const ConservativeState &state, const Vector2 &unit_normal,
                                                    const CaloricallyPerfectGas &gas);

/// Fhat = (Fn_L+Fn_R)/2 - lambda*(U_R-U_L)/2, per unit surface.
///
/// Same normal/state contract as compute_physical_euler_flux().
/// The normal points from the left state toward the right state.
/// Thermodynamics is evaluated once per side by the supplied model; no gamma/R appears here.
/// No face area, Mesh, boundary condition or time integration is involved.
/// A future assembler multiplies by revolved area and handles zero-area axis
/// faces without dividing by that area to construct a normal.
/// @throws std::invalid_argument If the normal is invalid.
/// @throws std::runtime_error For invalid gas states or non-representable fluxes/speed.
[[nodiscard]] RusanovFluxResult compute_rusanov_flux(const ConservativeState &left, const ConservativeState &right,
                                                     const Vector2 &unit_normal, const CaloricallyPerfectGas &gas);

} // namespace cfd
