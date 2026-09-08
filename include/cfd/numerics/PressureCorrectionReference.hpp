#pragma once

#include "cfd/mesh/Types.hpp"

namespace cfd
{

class ScalarLinearSystem;

/// Applies the exact zero constraint `p'_R = 0` to a pressure-correction system.
///
/// Both directed couplings on every internal face incident to `reference_cell`
/// are set to zero, preserving matrix symmetry. The existing positive reference
/// diagonal is retained and its right-hand side is set to zero. All unrelated
/// entries remain unchanged.
///
/// This utility performs no allocation and does not decide whether a reference
/// is required; callers should use it only for a system with an unfixed
/// pressure-correction gauge.
///
/// @throws std::invalid_argument If `reference_cell` is outside the system.
/// @throws std::runtime_error If its diagonal is not finite and strictly
///         positive.
void apply_zero_pressure_correction_reference(Index reference_cell, ScalarLinearSystem &system);

} // namespace cfd
