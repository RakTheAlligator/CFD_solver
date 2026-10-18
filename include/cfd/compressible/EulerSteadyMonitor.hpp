#pragma once

#include "cfd/compressible/AxisymmetricEulerSolver.hpp"

#include <array>

namespace cfd
{

struct EulerSteadyOptions
{
    double reference_density{1.0};  ///< Fixed physical rho_ref [kg/m^3], positive finite.
    double reference_velocity{1.0}; ///< Fixed U_ref [m/s], e.g. reservoir sound speed, positive finite.
    double reference_length{1.0};   ///< Fixed L_ref [m], positive finite.
    double residual_tolerance{1.0e-7};
    double field_change_tolerance{1.0e-7};
    double mass_flow_change_tolerance{1.0e-7};
    double mass_flow_imbalance_tolerance{1.0e-6};
    Index check_interval{50};
    Index minimum_steps{100};
    Index maximum_steps{100000}; ///< Absolute accepted-step counter limit, including pre-monitor steps.
    Index consecutive_checks{3};
    double maximum_time{100.0}; ///< Absolute physical time limit [s], later than initial monitor time.
};

enum class EulerSteadyStatus : std::uint8_t
{
    Running,
    Converged,
    StepLimit,
    TimeLimit
};

struct EulerSteadyObservation
{
    AxisymmetricEulerResult temporal; ///< Last accepted-step diagnostics, unchanged by observation.
    EulerResidualEvaluation current;  ///< Residual/fluxes of the current state, not the preceding update.
    std::array<double, 4> normalized_residual{};
    std::array<double, 4> field_change{};
    double mass_flow_relative_imbalance{};
    double mass_flow_change{};
    EulerSteadyStatus status{EulerSteadyStatus::Running};
};

/// Previous observed field and the actual observation history (never synthesized).
struct EulerSteadyRestartState
{
    ConservativeField previous;
    EulerSteadyObservation last;
    Index stable_checks{};
    bool observed{};
};

/// Non-owning observer of one Euler solver, which must outlive the monitor
/// and must not be moved from while the monitor exists.
/// Owns a once-allocated previous-check SoA snapshot, overwritten in place only
/// at observations. No physical update, I/O or allocation per observation.
/// Qref=(rho,rho*U,rho*U,rho*U^2), tref=L/U, V=sum(V_C).
/// Rnorm_k=sum(abs(R_C,k))/(V*Qref_k/tref), independent of dt and cancellation.
/// Field changes are volume-L1 differences / (V*Qref_k); flow changes use
/// rho*U*V/L. Boundary imbalance is |out-in|/(out+in), zero when both vanish.
/// Every equation/change and both flow tests must pass over consecutive checks
/// after minimum_steps. Small dt or a globally balanced transient cannot alone
/// converge. References must be adapted to physical scales; no epsilon floor.
class EulerSteadyMonitor
{
  public:
    /// @throws std::invalid_argument For invalid scales, tolerances or limits.
    explicit EulerSteadyMonitor(AxisymmetricEulerSolver &solver, EulerSteadyOptions options = {});
    /// Restores history; extending time/step budgets is allowed, not convergence tolerances.
    EulerSteadyMonitor(AxisymmetricEulerSolver &solver, EulerSteadyOptions options,
                       const EulerSteadyRestartState &restart);
    [[nodiscard]] EulerSteadyRestartState restart_state() const;
    EulerSteadyMonitor(const EulerSteadyMonitor &) = delete;
    EulerSteadyMonitor &operator=(const EulerSteadyMonitor &) = delete;
    EulerSteadyMonitor(EulerSteadyMonitor &&) = delete;
    EulerSteadyMonitor &operator=(EulerSteadyMonitor &&) = delete;
    ~EulerSteadyMonitor() = default;

    /// True initially, at scheduled accepted-step checks, and at either limit.
    [[nodiscard]] bool needs_observation() const noexcept;
    /// Observes current state; may be called earlier than the schedule, but
    /// never twice at the same accepted step/time. Convergence has priority
    /// when a limit and all convergence conditions are met simultaneously.
    /// @throws std::invalid_argument For a repeated/non-increasing observation.
    /// @throws std::runtime_error For assembly/BC or non-finite diagnostic failures.
    [[nodiscard]] EulerSteadyObservation observe();

  private:
    EulerSteadyMonitor(AxisymmetricEulerSolver &solver, EulerSteadyOptions options,
                       const EulerSteadyRestartState *restart);
    AxisymmetricEulerSolver *solver_;
    EulerSteadyOptions options_;
    ConservativeField previous_;
    std::array<long double, 4> field_scales_{};
    long double volume_{};
    long double reference_time_{};
    long double mass_flow_scale_{};
    double previous_in_{}, previous_out_{}, previous_time_{};
    Index previous_step_{}, stable_checks_{};
    bool observed_{};
    EulerSteadyObservation last_;
};

} // namespace cfd
