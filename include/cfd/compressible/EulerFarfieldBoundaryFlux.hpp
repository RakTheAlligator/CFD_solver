#pragma once

#include "cfd/compressible/EulerBoundaryConditions.hpp"
#include "cfd/compressible/EulerFlux.hpp"

namespace cfd
{

enum class EulerFarfieldRegime : std::uint8_t
{
    SupersonicInflow,
    SubsonicInflow,
    SubsonicOutflow,
    SupersonicOutflow
};

struct EulerFarfieldBoundaryFluxResult
{
    ConservativeState boundary_state; ///< Normal Riemann sample at x_normal/t=0, not a ghost state.
    EulerFlux flux;                   ///< Physical sample flux per unit revolution area.
    double maximum_wave_speed{};      ///< Covers initial/star characteristic and shock speeds, for CFL.
    double normal_mach{};             ///< Sample un/a, outward normal; NOT interior/total Mach.
    EulerFarfieldRegime regime{};     ///< Sample: M<=-1, -1<M<0, 0<=M<1, M>=1 respectively.
};

/// Full-state, steady, monogas farfield using the entropy solution of the local
/// normal perfect-gas Riemann problem (interior left, exterior right). Physical
/// trace flux only; no additional Rusanov dissipation. Internal fluxes unchanged.
/// Handles shocks, sonic fans and incoming entropy/tangent, face by face.
/// At a stationary contact choose the interior-side trace; its flux is unique.
/// Purely outgoing waves give the interior flux independently of exterior data;
/// a strong incoming shock can invalidate a classification based on interior
/// Mach alone. Not an exact multidimensional nonreflecting boundary condition.
/// See docs/AxisymmetricEulerFarfield.md for equations and domain of validity.
/// Same finite unit-normal contract as EulerFlux. No allocation on valid calls.
/// @throws std::invalid_argument For a non-Farfield condition or invalid normal.
/// @throws std::runtime_error For inadmissible states, vacuum-generating data,
///         unrepresentable arithmetic or an unresolved pressure root; never clipped.
[[nodiscard]] EulerFarfieldBoundaryFluxResult compute_euler_farfield_boundary_flux(
    const ConservativeState &inside, const Vector2 &outward_unit_normal, const EulerBoundaryCondition &condition,
    const CaloricallyPerfectGas &gas);

} // namespace cfd
