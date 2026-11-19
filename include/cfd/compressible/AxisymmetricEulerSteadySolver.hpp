#pragma once

#include "cfd/compressible/AxisymmetricEulerSolver.hpp"

#include <chrono>
#include <cstdint>
#include <stdexcept>
#include <string>

namespace cfd
{

struct EulerLocalSteadyOptions
{
    double cfl{0.05};
    Index maximum_iterations{10000};
    Index maximum_step_attempts{20};
    Index maximum_rejected_attempts{100};
    Index check_interval{20};
    Index minimum_iterations{20};
    Index consecutive_checks{3};
    double maximum_wall_seconds{60.0}; ///< Real elapsed time, not physical or accumulated pseudo-time.
    double reference_density{1.0};     ///< Fixed rho_ref [kg/m^3], strictly positive finite.
    double reference_velocity{1.0};    ///< Fixed U_ref [m/s], strictly positive finite.
    double reference_length{1.0};      ///< Fixed L_ref [m], strictly positive finite.
    double residual_tolerance{1.0e-7}; ///< Both global L1 and maximum cell-normalized residual, every equation.
    double residual_reduction{1.0e-2}; ///< Required reduction from initial maximum cell residual, if above tolerance.
    double field_change_tolerance{1.0e-7};
    double flow_change_tolerance{1.0e-7};
    double balance_tolerance{1.0e-7};       ///< All four stationary boundary/source balances.
    double maximum_residual_growth{1.0e6};  ///< Growth relative to max(initial residual, residual_tolerance).
    Index stagnation_checks{50};            ///< Observations per complete residual-peak window.
    Index stagnation_windows{3};            ///< Consecutive windows without sufficient progress before STAGNATED.
    double minimum_residual_progress{0.01}; ///< Required per-equation peak decrease from its progress reference.
    bool collect_timings{false};
};

enum class EulerLocalSteadyStatus : std::uint8_t
{
    Running,
    Converged,
    MaximumIterations,
    WallLimit,
    Stagnated,
    Diverged,
    NumericalFailure
};

[[nodiscard]] const char *to_string(EulerLocalSteadyStatus status) noexcept;

struct EulerRegionResidual
{
    Index region{}; ///< Caller-provided label; label zero for the default single region.
    double volume{};
    std::array<double, 4> absolute_residual{}; ///< Integrated L1 residuals [kg/s, N, N, W].
    std::array<double, 4> normalized_residual{};
    std::array<double, 4> maximum_cell_residual{};
};

struct EulerLocalSteadyResult
{
    Index iteration{};
    Index observation_iteration{}; ///< Residuals/fluxes refer to this accepted state, not a preceding update.
    Index rejected_attempts{};
    double minimum_local_step{}, maximum_local_step{}; ///< Last accepted alpha*dtau [s], not physical dt.
    double step_scale{1.0};
    double minimum_density{}, minimum_pressure{};
    EulerIntegralQuantities integrals; ///< Integrals of the current state, NOT a temporal conservation ledger.
    std::array<double, 4> absolute_residual{};
    std::array<double, 4> normalized_residual{};
    std::array<double, 4> maximum_cell_residual{};
    std::array<double, 4> field_change{};
    std::array<double, 4> maximum_cell_change{}; ///< max_C |U_C-U_previous_check,C|/Qref, no volume dilution.
    std::array<double, 4> residual_reduction{};  ///< Rmax/Rmax_initial; Rmax/tolerance if initially below tolerance.
    std::array<double, 4> stationary_balance{};
    EulerIntegralQuantities boundary_outflow_rate;
    double radial_pressure_source_rate{};
    double mass_flow_in{}, mass_flow_out{}, mass_flow_relative_imbalance{}, flow_change{};
    EulerLocalSteadyStatus status{EulerLocalSteadyStatus::Running};
};

/// Explicit local pseudo-time / Jacobi iteration of the existing Euler FV residual.
/// Select this type for steady solves; AxisymmetricEulerSolver remains physical/global.
/// No physical time, temporal ledger, or transient checkpoint/restart is exposed.
/// Mesh must outlive the solver and remain unmoved. An initial accepted field is
/// deep-copied; it may come from a transient state, but monitor history is reset.
/// Optional region labels are copied/remapped once; arbitrary labels are allowed.
/// All state/candidate/residual/step/monitor buffers are persistent: no allocations
/// on an unrejected iterate()/observe(). Rejection/terminal messages may allocate.
///
/// Qref=(rho,rho*U,rho*U,rho*U^2), tref=L/U.
/// Cell residual = |R_C,k|/(V_C*Qref_k/tref); its maximum prevents a large farfield
/// from hiding ENGINE. Global/region L1 use sum(V_C) of the respective set.
/// Stationary balances use Qref_k*U*L^2, including the radial pressure source.
/// All four equations, required reduction, field/flow stability and balances must
/// pass over consecutive observations. No residual is inferred from a small step.
class AxisymmetricEulerSteadySolver
{
  public:
    /// Invalid options/labels throw invalid_argument; initial geometry/gas errors
    /// propagate as in the transient constructor. Iteration failures are reported
    /// as NumericalFailure with the last accepted state and failure_message().
    AxisymmetricEulerSteadySolver(const Mesh &mesh, const ConservativeField &initial, CaloricallyPerfectGas gas,
                                  EulerBoundaryConditions conditions, EulerLocalSteadyOptions options = {},
                                  std::span<const Index> cell_regions = {});
    AxisymmetricEulerSteadySolver(const AxisymmetricEulerSteadySolver &) = delete;
    AxisymmetricEulerSteadySolver &operator=(const AxisymmetricEulerSteadySolver &) = delete;
    AxisymmetricEulerSteadySolver(AxisymmetricEulerSteadySolver &&) = delete;
    AxisymmetricEulerSteadySolver &operator=(AxisymmetricEulerSteadySolver &&) = delete;
    ~AxisymmetricEulerSteadySolver() = default;

    [[nodiscard]] const ConservativeField &state() const noexcept
    {
        return engine_.state();
    }

    [[nodiscard]] const AxisymmetricGeometry &geometry() const noexcept
    {
        return engine_.geometry();
    }

    [[nodiscard]] const EulerLocalSteadyResult &diagnostics() const noexcept
    {
        return result_;
    }

    /// Workspace view: zero before the first iteration, effective alpha*dtau
    /// after acceptance; a terminal failed attempt may leave unscaled proposals.
    [[nodiscard]] std::span<const double> local_steps() const noexcept
    {
        return local_steps_;
    }

    [[nodiscard]] std::span<const EulerRegionResidual> regional_residuals() const noexcept
    {
        return regions_;
    }

    [[nodiscard]] const AxisymmetricEulerTimings &timings() const noexcept
    {
        return engine_.timings();
    }

    [[nodiscard]] const std::string &failure_message() const noexcept
    {
        return failure_;
    }

    /// Fresh residual/flux observation, scheduled automatically by iterate().
    /// Repeated calls at the same iteration are idempotent (no fabricated checks).
    const EulerLocalSteadyResult &observe();
    /// One accepted Jacobi iteration, with bounded GLOBAL halvings on rejection.
    /// Candidates pass the EOS and all nonzero-area boundary flux evaluations.
    /// Only candidate EOS failures and typed characteristic incompatibilities retry.
    /// Other BC/assembly failures remain fatal; compatibility does not imply convergence.
    /// A terminal result is a no-op; observe() cannot resume a terminal solve.
    const EulerLocalSteadyResult &iterate();
    const EulerLocalSteadyResult &solve();

  private:
    void observe_current();
    void update_candidate();
    bool candidate_boundaries_compatible(const ConservativeField &candidate, double alpha);
    bool wall_limit();
    void fail(const std::runtime_error &error);
    AxisymmetricEulerSolver engine_; // Shared FV assembly/storage, never stepped in physical time.
    EulerLocalSteadyOptions options_;
    ConservativeField previous_;
    std::vector<double> local_steps_;
    std::vector<Index> boundary_faces_; ///< Nonzero-area boundary faces, collected once.
    std::vector<Index> region_index_;
    std::vector<EulerRegionResidual> regions_;
    std::vector<long double> region_volumes_;
    std::vector<std::array<long double, 4>> region_sums_;
    std::array<long double, 4> field_scales_{};
    std::array<double, 4> initial_residual_{};
    long double volume_{}, reference_time_{}, flow_scale_{};
    std::array<double, 4> window_peaks_{}, reference_window_peaks_{};
    double previous_in_{}, previous_out_{};
    Index stable_checks_{}, window_checks_{}, unproductive_windows_{};
    bool observed_{};
    bool have_reference_window_{};
    std::optional<std::chrono::steady_clock::time_point> start_;
    EulerLocalSteadyResult result_;
    std::string failure_;
};
} // namespace cfd
