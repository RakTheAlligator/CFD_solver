#pragma once

#include "cfd/compressible/CaloricallyPerfectGas.hpp"
#include "cfd/compressible/EulerBoundaryConditions.hpp"
#include "cfd/geometry/AxisymmetricGeometry.hpp"

#include <array>
#include <limits>
#include <optional>
#include <span>
#include <vector>

namespace cfd
{

struct AxisymmetricEulerOptions
{
    double cfl{0.2};                 ///< Finite, 0 < CFL <= 1; positivity may require a smaller step.
    Index maximum_step_attempts{20}; ///< Includes the first candidate; subsequent attempts halve dt globally.
    Index maximum_steps{1000000};    ///< Total accepted-step limit over the solver lifetime.
    bool collect_timings{false};     ///< Optional wall-clock instrumentation, including per-boundary-face clocks.
    Index maximum_rejected_steps{std::numeric_limits<Index>::max()}; ///< Cumulative safety limit, never exceeded.
};

struct AxisymmetricEulerTimings
{
    double assembly_seconds{};      ///< Source, face fluxes, BCs and CFL; includes explicit residual preparation.
    double boundary_flux_seconds{}; ///< Subset of assembly_seconds; instrumentation overhead is not subtracted.
    double time_update_seconds{};   ///< Candidates, admissibility, accepted-state summary and commit.
};

/// Full-revolution integrals: mass [kg], momenta [kg*m/s], energy [J].
/// For a boundary rate the respective units are per second.
struct EulerIntegralQuantities
{
    double mass{};
    double momentum_x{};
    double momentum_r{};
    double energy{};
    bool operator==(const EulerIntegralQuantities &) const = default;
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
    bool operator==(const AxisymmetricEulerResult &) const = default;
};

/// Owned accepted state, not a visualization field. Constructed only on demand.
struct AxisymmetricEulerRestartState
{
    ConservativeField conservative;
    AxisymmetricEulerResult diagnostics;
};

/// Current-state integrated residual L1 sums: [kg/s, N, N, W]. No dt factor.
/// Boundary rates use exactly the same numerical flux as the update.
struct EulerResidualEvaluation
{
    std::array<double, 4> absolute_residual{};
    double mass_flow_in{};  ///< Positive inward mass flow over all boundary faces [kg/s].
    double mass_flow_out{}; ///< Positive outward mass flow over all boundary faces [kg/s].
};

/// Current-state CFL proposal, before final-time clipping or positivity retries.
/// The read-only spectral sums belong to the solver; the view is invalidated by
/// its next accepted step, move or destruction. No allocation or field copy.
struct EulerCflEvaluation
{
    double time_step{};
    std::span<const double> cell_wave_sums; ///< sum_f(|S_axi,f| * lambda_f), using the actual flux/BC waves.
};

/// First-order, no-swirl axisymmetric Euler / Rusanov / global forward Euler.
///
/// Owns precomputed revolution metrics, BCs and persistent current/candidate,
/// residual and wave-sum buffers. The initial state is deep-copied once; steps
/// alternate buffers without field copies or allocations on the valid path.
/// Mesh must outlive this object and must not be moved from or replaced.
/// Views from state() remain valid only until the next accepted step or move.
/// No I/O, implicit solve, clipping or viscosity. Open BCs use the local normal
/// characteristic closures documented in CompressibleBoundaryConditions.md.
/// Farfield uses the normal Riemann boundary in AxisymmetricEulerFarfield.md.
class AxisymmetricEulerSolver
{
  public:
    /// @throws std::invalid_argument For invalid options, cardinality or BC geometry.
    /// @throws std::runtime_error For invalid geometry or inadmissible initial gas states.
    AxisymmetricEulerSolver(const Mesh &mesh, const ConservativeField &initial_state, CaloricallyPerfectGas gas,
                            EulerBoundaryConditions boundary_conditions, AxisymmetricEulerOptions options = {});

    /// Restores a validated accepted state; scratch is recomputed, timing counters restart at zero.
    /// Mesh/EOS/BC compatibility across files belongs to the checkpoint reader.
    /// @throws std::invalid_argument For inconsistent times, counters or diagnostics.
    AxisymmetricEulerSolver(const Mesh &mesh, const AxisymmetricEulerRestartState &restart, CaloricallyPerfectGas gas,
                            EulerBoundaryConditions boundary_conditions, AxisymmetricEulerOptions options = {});

    /// Deep copy of conservative storage, only when explicitly requested (e.g. checkpoint).
    [[nodiscard]] AxisymmetricEulerRestartState restart_state() const;

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

    [[nodiscard]] const AxisymmetricGeometry &geometry() const noexcept
    {
        return geometry_;
    }

    [[nodiscard]] const AxisymmetricEulerTimings &timings() const noexcept
    {
        return timings_;
    }

    /// Evaluates sum(abs(R_C)) on the current accepted state using the existing
    /// assembler, never a stale pre-update residual. Does not change state,
    /// time or accepted-step diagnostics. Prepared scratch can serve the next
    /// step; accepting a step invalidates it. No allocation on the valid path.
    /// Same flux/BC/arithmetic failures as step(), propagated without retry.
    [[nodiscard]] EulerResidualEvaluation evaluate_steady_residual();

    /// Explicitly requested diagnostic using the same cached assembly as step().
    /// Does not advance time or change accepted-state diagnostics. No extra work
    /// in ordinary steps; assembly errors propagate exactly as for step().
    [[nodiscard]] EulerCflEvaluation evaluate_cfl();

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
    double prepare_assembly();
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
    std::vector<double> boundary_group_areas_; ///< Physical areas for group mass-flow inlets, built once.
    Index current_{};
    AxisymmetricEulerResult result_;
    EulerIntegralQuantities boundary_outflow_;
    double radial_source_{};
    double mass_flow_in_{};
    double mass_flow_out_{};
    std::optional<double> assembled_time_step_; ///< Valid only for the current accepted state and steady BCs.
    AxisymmetricEulerTimings timings_;
};

} // namespace cfd
