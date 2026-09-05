#pragma once

namespace cfd
{

class CellMomentumPressureResponse;
class Mesh;
class ScalarLinearSystem;

/// Computes cell-centered momentum responses from final momentum diagonals.
///
/// For each cell `P`, the output is overwritten with
/// `d_P,u = A_P / a_P,u` and `d_P,v = A_P / a_P,v`. Here `A_P` is the
/// two-dimensional cell area per unit depth, while `a_P,u` and `a_P,v` are
/// used exactly as stored in the final assembled and equation-under-relaxed
/// momentum systems. The two component responses are intentionally independent.
///
/// Repeated valid calls perform two sequential O(N) passes, allocate no memory,
/// and reuse the caller-owned output storage. All validation completes before
/// either output component is modified.
///
/// @throws std::invalid_argument If an object cardinality is incompatible or a
///         momentum system does not reference the exact supplied Mesh.
/// @throws std::runtime_error If a cell area, momentum diagonal, or computed
///         response is non-finite or not strictly positive.
void compute_momentum_pressure_response(const Mesh &mesh, const ScalarLinearSystem &u_momentum_system,
                                        const ScalarLinearSystem &v_momentum_system,
                                        CellMomentumPressureResponse &response);

} // namespace cfd
