#pragma once

#include "cfd/math/Vector2.hpp"

namespace cfd
{
class Mesh;
class ScalarBoundaryConditions;
class PressureCorrectionBoundaryConditions;
class CellVelocityField;
class CellScalarField;
} // namespace cfd

namespace cfd::app
{

/// Initializes a uniform bulk cell state from resolved user boundary data.
///
/// Velocity averages only faces with two Dirichlet velocity components and
/// `Ub . Sf < -64 * epsilon * |Ub| * |Sf|`, using face-length weights.
/// Pressure averages all Dirichlet pressure faces with the same weights.
/// Without eligible inflow/pressure faces, the respective internalField
/// fallback is used. Pressure-correction types must match physical pressure
/// Dirichlet (FixedPressure) or Neumann (FixedMassFlux) conditions.
///
/// All inputs are borrowed for this call; only caller-owned u/v/p are written.
/// No face flux is accepted or changed, and no CFD equation is solved.
/// Validation and averaging precede all writes. The calculation allocates
/// nothing and uses O(Nboundaries + Nfaces + Ncells) work.
/// @throws std::invalid_argument For incompatible cardinalities, aliased
///         outputs, non-finite boundary/fallback values or inconsistent BC types.
/// @throws std::runtime_error For unusable face geometry or non-finite averages.
void initialize_from_boundary_conditions(
    const Mesh &mesh, const ScalarBoundaryConditions &u_boundary_conditions,
    const ScalarBoundaryConditions &v_boundary_conditions, const ScalarBoundaryConditions &pressure_boundary_conditions,
    const PressureCorrectionBoundaryConditions &pressure_correction_boundary_conditions, Vector2 fallback_velocity,
    double fallback_pressure, CellVelocityField &velocity, CellScalarField &pressure);

} // namespace cfd::app
