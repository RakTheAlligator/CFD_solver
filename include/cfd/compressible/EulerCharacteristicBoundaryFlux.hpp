#pragma once

#include "cfd/compressible/EulerBoundaryConditions.hpp"
#include "cfd/compressible/EulerFlux.hpp"

namespace cfd
{

struct EulerCharacteristicBoundaryFluxResult
{
    ConservativeState boundary_state; ///< Reconstructed trace, not a Rusanov ghost state.
    EulerFlux flux;                   ///< Physical trace flux per unit revolution surface.
    double maximum_wave_speed{};      ///< max(|un_i|+a_i, |un_b|+a_b), for CFL.
    double normal_mach{};             ///< Trace un_b/a_b; outward normal convention.
};

/// First-order local normal characteristic closure for the five new open BCs.
/// See docs/CompressibleBoundaryConditions.md for derivation and restrictions.
/// The unit normal has the EulerFlux contract. Existing Axis/SlipWall/
/// Transmissive must use their existing flux routines, not this function.
/// MassFlowInlet requires the finite positive sum of physical revolution face
/// areas over its group, supplied once by the caller; other types ignore it.
/// Its positive inward group mdot is distributed uniformly per physical area.
/// Subsonic data extrapolate J+=un_i+2*a_i/(gamma-1); no second Rusanov
/// dissipation is applied to the resulting trace. PressureOutlet extrapolates
/// at un_i>=a_i, ignoring its prescribed pressure. Backflow, incompatible
/// inlet/outlet roots and inlet choking are rejected, never clipped.
/// Data are steady; no time dependence or hidden state. No allocations on
/// successful calls. This is not an exact multidimensional nonreflecting BC.
/// @throws std::invalid_argument For an invalid normal, type or required group area.
/// @throws std::runtime_error For invalid states, incompatible characteristic
///         regimes or non-representable arithmetic. Caller must propagate BC
///         failures, not conceal them as time-step backtracking.
[[nodiscard]] EulerCharacteristicBoundaryFluxResult compute_euler_characteristic_boundary_flux(
    const ConservativeState &inside, const Vector2 &outward_unit_normal, const EulerBoundaryCondition &condition,
    const CaloricallyPerfectGas &gas, double group_revolution_area = 0.0);

} // namespace cfd
