#pragma once

#include "cfd/compressible/CaloricallyPerfectGas.hpp"
#include "cfd/compressible/EulerBoundaryConditions.hpp"
#include "cfd/geometry/AxisymmetricGeometry.hpp"

#include <array>
#include <vector>

namespace cfd
{

struct AxisymmetricEulerOptions
{
    double cfl{0.2};                 ///< Finite, 0 < CFL <= 1; positivity may require a smaller step.
    Index maximum_step_attempts{20}; ///< Includes the first candidate; subsequent attempts halve dt globally.
    Index maximum_steps{1000000};    ///< Total accepted-step limit over the solver lifetime.
};

/// Full-revolution integrals: mass [kg], momenta [kg*m/s], energy [J].
/// For a boundary rate the respective units are per second.
struct EulerIntegralQuantities
{
    double mass{};
    double momentum_x{};
    double momentum_r{};
    double energy{};
};

/// Diagnostics of the last accepted step; counters are cumulative.
struct AxisymmetricEulerResult
{
    double time{};
    double time_step{};
    double cfl{}; ///< Configured CFL upper bound, before final-time clipping/backtracking.
    Index accepted_steps{};
    Index rejected_steps{};
    double minimum_density{};
    double maximum_density{};
    double minimum_pressure{};
    double maximum_pressure{};
    EulerIntegralQuantities integrals;
    EulerIntegralQuantities boundary_outflow_rate;
    double radial_pressure_source_rate{};
    /// I_new - I_old + dt*boundary_outflow - dt*source (radial momentum only).
    EulerIntegralQuantities balance_error;
};

/// First-order, no-swirl axisymmetric Euler / Rusanov / global forward Euler.
///
/// Owns precomputed revolution metrics, BCs and persistent current/candidate,
/// residual and wave-sum buffers. The initial state is deep-copied once; steps
/// alternate buffers without field copies or allocations on the valid path.
/// Mesh must outlive this object and must not be moved from or replaced.
/// Views from state() remain valid only until the next accepted step or move.
/// No I/O, implicit solve, clipping, viscosity or characteristic outlet model.
class AxisymmetricEulerSolver
{
  public:
    /// @throws std::invalid_argument For invalid options, cardinality or BC geometry.
    /// @throws std::runtime_error For invalid geometry or inadmissible initial gas states.
    AxisymmetricEulerSolver(const Mesh &mesh, const ConservativeField &initial_state, CaloricallyPerfectGas gas,
                            EulerBoundaryConditions boundary_conditions, AxisymmetricEulerOptions options = {});

    AxisymmetricEulerSolver(const AxisymmetricEulerSolver &) = delete;
    AxisymmetricEulerSolver &operator=(const AxisymmetricEulerSolver &) = delete;
    AxisymmetricEulerSolver(AxisymmetricEulerSolver &&) noexcept = default;
    AxisymmetricEulerSolver &operator=(AxisymmetricEulerSolver &&) noexcept = delete;
    ~AxisymmetricEulerSolver() = default;

    [[nodiscard]] const ConservativeField &state() const noexcept
    {
        return states_.at(current_);
    }

    [[nodiscard]] const AxisymmetricEulerResult &diagnostics() const noexcept
    {
        return result_;
    }

    /// Accepts one step, clipped to final_time, or does nothing if already there.
    /// dt = CFL*min(V/sum(A*lambda)); zero-area axis faces contribute nothing.
    /// The radial source is explicit; this face CFL alone is not a positivity
    /// theorem. Inadmissible candidates halve dt, with a bounded attempt count.
    /// Only candidate-state failures are retried: flux/BC/assembly errors propagate.
    /// On failure the last accepted state/time are preserved; rejected-step counts
    /// include all tested invalid candidates. Diagnostics describe accepted steps.
    /// @throws std::invalid_argument If final_time is non-finite or earlier than time.
    /// @throws std::runtime_error For non-representable arithmetic, no positive
    ///         time advance, exhausted attempts (with cell/state/dt), or step limit.
    const AxisymmetricEulerResult &step(double final_time);

    /// Repeatedly calls step() up to final_time. No callbacks or hidden output.
    const AxisymmetricEulerResult &advance_to(double final_time);

  private:
    double assemble();
    AxisymmetricEulerResult summarize(const ConservativeField &field) const;

    const Mesh *mesh_;
    AxisymmetricGeometry geometry_;
    CaloricallyPerfectGas gas_;
    EulerBoundaryConditions boundary_conditions_;
    AxisymmetricEulerOptions options_;
    std::array<ConservativeField, 2> states_;
    ConservativeField residual_;
    std::vector<double> wave_sums_;
    std::vector<double> face_areas_;
    Index current_{};
    AxisymmetricEulerResult result_;
    EulerIntegralQuantities boundary_outflow_;
    double radial_source_{};
};

} // namespace cfd
