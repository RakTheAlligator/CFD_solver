#include "cfd/compressible/EulerSteadyMonitor.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cfd
{
namespace
{
double checked_diagnostic(const long double value)
{
    if (!std::isfinite(value) || value < 0.0L || value > std::numeric_limits<double>::max())
        throw std::runtime_error("Euler steady diagnostic is not finite/nonnegative and representable.");
    return static_cast<double>(value);
}

void validate(const EulerSteadyOptions &options)
{
    for (const double positive :
         {options.reference_density, options.reference_velocity, options.reference_length, options.residual_tolerance,
          options.field_change_tolerance, options.mass_flow_change_tolerance, options.mass_flow_imbalance_tolerance,
          options.maximum_time})
        if (!std::isfinite(positive) || !(positive > 0.0))
            throw std::invalid_argument(
                "Euler steady references, tolerances and maximum time must be finite positive.");
    if (options.check_interval == 0 || options.minimum_steps == 0 || options.maximum_steps < options.minimum_steps ||
        options.consecutive_checks == 0)
        throw std::invalid_argument("Euler steady step/check limits must be positive and ordered.");
}
} // namespace

EulerSteadyMonitor::EulerSteadyMonitor(AxisymmetricEulerSolver &solver, EulerSteadyOptions options)
    : EulerSteadyMonitor(solver, options, nullptr)
{
}

EulerSteadyMonitor::EulerSteadyMonitor(AxisymmetricEulerSolver &solver, EulerSteadyOptions options,
                                       const EulerSteadyRestartState &restart)
    : EulerSteadyMonitor(solver, options, &restart)
{
}

EulerSteadyMonitor::EulerSteadyMonitor(AxisymmetricEulerSolver &solver, EulerSteadyOptions options,
                                       const EulerSteadyRestartState *restart)
    : solver_(&solver), options_(options), previous_(solver.state()), previous_time_(solver.diagnostics().time),
      previous_step_(solver.diagnostics().accepted_steps)
{
    validate(options_);
    last_.temporal = solver.diagnostics();
    if (restart == nullptr && (!(options_.maximum_time > previous_time_) || !(options_.maximum_steps > previous_step_)))
        throw std::invalid_argument("Euler steady limits must lie after the initial monitor state.");
    const long double density{options_.reference_density}, velocity{options_.reference_velocity};
    field_scales_ = {density, density * velocity, density * velocity, density * velocity * velocity};
    reference_time_ = static_cast<long double>(options_.reference_length) / velocity;
    for (const double volume : solver.geometry().cell_volumes())
        volume_ += volume;
    mass_flow_scale_ = density * velocity * volume_ / options_.reference_length;
    for (const auto scale : field_scales_)
        if (!std::isfinite(scale) || !(scale > 0.0L) || !std::isfinite(volume_ * scale / reference_time_) ||
            !(volume_ * scale / reference_time_ > 0.0L))
            throw std::invalid_argument("Euler steady physical normalization is not representable.");
    if (!std::isfinite(mass_flow_scale_) || !(mass_flow_scale_ > 0.0L))
        throw std::invalid_argument("Euler steady mass-flow normalization is not representable.");
    if (restart != nullptr)
    {
        const auto &old{restart->last.temporal};
        if (restart->previous.size() != solver.state().size() || !std::isfinite(old.time) || old.time < 0 ||
            old.time > solver.diagnostics().time || old.accepted_steps > solver.diagnostics().accepted_steps ||
            (old.accepted_steps == solver.diagnostics().accepted_steps && old.time != solver.diagnostics().time) ||
            (old.accepted_steps < solver.diagnostics().accepted_steps && !(old.time < solver.diagnostics().time)) ||
            restart->stable_checks > options.consecutive_checks ||
            (restart->last.status == EulerSteadyStatus::Converged &&
             old.accepted_steps != solver.diagnostics().accepted_steps) ||
            (!restart->observed && restart->stable_checks != 0) || solver.diagnostics().time > options.maximum_time ||
            solver.diagnostics().accepted_steps > options.maximum_steps)
            throw std::invalid_argument("Euler monitor restart has inconsistent history or limits.");
        for (Index cell = 0; cell < previous_.size(); ++cell)
        {
            const auto value{restart->previous.state(cell)};
            if (!std::isfinite(value.density) || !(value.density > 0) || !std::isfinite(value.momentum_x) ||
                !std::isfinite(value.momentum_r) || !std::isfinite(value.total_energy_density))
                throw std::invalid_argument("Euler monitor restart previous field is invalid.");
            previous_.set_state(cell, value);
            if (old.accepted_steps == solver.diagnostics().accepted_steps)
            {
                const auto current{solver.state().state(cell)};
                if (value.density != current.density || value.momentum_x != current.momentum_x ||
                    value.momentum_r != current.momentum_r ||
                    value.total_energy_density != current.total_energy_density)
                    throw std::invalid_argument("Euler monitor restart current observation field is inconsistent.");
            }
        }
        for (const double value : {restart->last.current.mass_flow_in, restart->last.current.mass_flow_out})
            if (!std::isfinite(value) || value < 0)
                throw std::invalid_argument("Euler monitor restart flow history is invalid.");
        for (const auto values :
             {restart->last.current.absolute_residual, restart->last.normalized_residual, restart->last.field_change})
            for (const double value : values)
                if (!std::isfinite(value) || value < 0)
                    throw std::invalid_argument("Euler monitor restart residual/change history is invalid.");
        if (!std::isfinite(restart->last.mass_flow_relative_imbalance) ||
            restart->last.mass_flow_relative_imbalance < 0 || restart->last.mass_flow_relative_imbalance > 1 ||
            !std::isfinite(restart->last.mass_flow_change) || restart->last.mass_flow_change < 0 ||
            static_cast<unsigned>(restart->last.status) > static_cast<unsigned>(EulerSteadyStatus::TimeLimit) ||
            ((restart->last.status == EulerSteadyStatus::Converged) !=
             (restart->stable_checks == options.consecutive_checks)))
            throw std::invalid_argument("Euler monitor restart convergence history is inconsistent.");
        if (!restart->observed && restart->last.status != EulerSteadyStatus::Running)
            throw std::invalid_argument("Euler monitor restart cannot invent a terminal observation.");
        if (restart->stable_checks > 0)
        {
            bool passed{restart->observed && old.accepted_steps >= options.minimum_steps &&
                        restart->last.mass_flow_relative_imbalance <= options.mass_flow_imbalance_tolerance &&
                        restart->last.mass_flow_change <= options.mass_flow_change_tolerance};
            for (Index k = 0; k < 4; ++k)
                passed = passed && restart->last.normalized_residual.at(k) <= options.residual_tolerance &&
                         restart->last.field_change.at(k) <= options.field_change_tolerance;
            if (!passed)
                throw std::invalid_argument("Euler monitor restart stable checks contradict the diagnostics.");
        }
        last_ = restart->last;
        previous_in_ = last_.current.mass_flow_in;
        previous_out_ = last_.current.mass_flow_out;
        previous_time_ = old.time;
        previous_step_ = old.accepted_steps;
        observed_ = restart->observed;
        stable_checks_ = restart->stable_checks;
    }
}

EulerSteadyRestartState EulerSteadyMonitor::restart_state() const
{
    return {ConservativeField{previous_}, last_, stable_checks_, observed_};
}

bool EulerSteadyMonitor::needs_observation() const noexcept
{
    const auto &current{solver_->diagnostics()};
    return !observed_ || current.accepted_steps - previous_step_ >= options_.check_interval ||
           current.accepted_steps >= options_.maximum_steps || current.time >= options_.maximum_time;
}

EulerSteadyObservation EulerSteadyMonitor::observe()
{
    EulerSteadyObservation observation;
    observation.temporal = solver_->diagnostics();
    const auto &temporal{observation.temporal};
    if (observed_ && (temporal.accepted_steps <= previous_step_ || !(temporal.time > previous_time_)))
        throw std::invalid_argument("Euler steady observations require a new accepted step and increasing time.");
    observation.current = solver_->evaluate_steady_residual();
    std::array<long double, 4> changes{};
    for (Index cell = 0; cell < previous_.size(); ++cell)
    {
        const auto old{previous_.state(cell)}, current{solver_->state().state(cell)};
        const long double volume{solver_->geometry().cell_volumes()[cell]};
        changes[0] += volume * std::abs(static_cast<long double>(current.density) - old.density);
        changes[1] += volume * std::abs(static_cast<long double>(current.momentum_x) - old.momentum_x);
        changes[2] += volume * std::abs(static_cast<long double>(current.momentum_r) - old.momentum_r);
        changes[3] +=
            volume * std::abs(static_cast<long double>(current.total_energy_density) - old.total_energy_density);
    }
    bool passed{observed_ && temporal.accepted_steps >= options_.minimum_steps};
    for (Index component = 0; component < changes.size(); ++component)
    {
        observation.normalized_residual.at(component) =
            checked_diagnostic(observation.current.absolute_residual.at(component) * reference_time_ /
                               (volume_ * field_scales_.at(component)));
        observation.field_change.at(component) =
            checked_diagnostic(changes.at(component) / (volume_ * field_scales_.at(component)));
        passed = passed && observation.normalized_residual.at(component) <= options_.residual_tolerance &&
                 observation.field_change.at(component) <= options_.field_change_tolerance;
    }
    const long double in{observation.current.mass_flow_in}, out{observation.current.mass_flow_out};
    observation.mass_flow_relative_imbalance =
        checked_diagnostic(in + out == 0.0L ? 0.0L : std::abs(out - in) / (in + out));
    observation.mass_flow_change =
        checked_diagnostic(std::max(std::abs(in - previous_in_), std::abs(out - previous_out_)) / mass_flow_scale_);
    passed = passed && observation.mass_flow_relative_imbalance <= options_.mass_flow_imbalance_tolerance &&
             observation.mass_flow_change <= options_.mass_flow_change_tolerance;
    if (passed)
    {
        if (stable_checks_ < options_.consecutive_checks)
            ++stable_checks_;
    }
    else
        stable_checks_ = 0;
    if (stable_checks_ >= options_.consecutive_checks)
        observation.status = EulerSteadyStatus::Converged;
    else if (temporal.accepted_steps >= options_.maximum_steps)
        observation.status = EulerSteadyStatus::StepLimit;
    else if (temporal.time >= options_.maximum_time)
        observation.status = EulerSteadyStatus::TimeLimit;
    previous_in_ = observation.current.mass_flow_in;
    previous_out_ = observation.current.mass_flow_out;
    previous_step_ = temporal.accepted_steps;
    previous_time_ = temporal.time;
    observed_ = true;
    for (Index cell = 0; cell < previous_.size(); ++cell)
        previous_.set_state(cell, solver_->state().state(cell));
    last_ = observation;
    return observation;
}
} // namespace cfd
