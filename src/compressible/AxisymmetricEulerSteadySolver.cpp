#include "cfd/compressible/AxisymmetricEulerSteadySolver.hpp"
#include "cfd/compressible/EulerCharacteristicBoundaryFlux.hpp"
#include "cfd/mesh/Mesh.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <utility>

namespace cfd
{
namespace
{
using Clock = std::chrono::steady_clock;

double checked(const long double value)
{
    if (!std::isfinite(value) || std::abs(value) > std::numeric_limits<double>::max())
        throw std::runtime_error("Euler local steady arithmetic is not finite/representable.");
    return static_cast<double>(value);
}

void validate(const EulerLocalSteadyOptions &o)
{
    for (const double value :
         {o.maximum_wall_seconds, o.reference_density, o.reference_velocity, o.reference_length, o.residual_tolerance,
          o.field_change_tolerance, o.flow_change_tolerance, o.balance_tolerance})
        if (!std::isfinite(value) || !(value > 0.0))
            throw std::invalid_argument(
                "Euler local steady scales, tolerances and wall limit must be finite positive.");
    if (!std::isfinite(o.cfl) || !(o.cfl > 0.0 && o.cfl <= 1.0) || !std::isfinite(o.residual_reduction) ||
        !(o.residual_reduction > 0.0 && o.residual_reduction <= 1.0) || !std::isfinite(o.minimum_residual_progress) ||
        !(o.minimum_residual_progress > 0.0 && o.minimum_residual_progress < 1.0) ||
        !std::isfinite(o.maximum_residual_growth) || !(o.maximum_residual_growth > 1.0))
        throw std::invalid_argument("Euler local steady CFL/reduction/progress/growth is invalid.");
    if (o.maximum_iterations == 0 || o.maximum_step_attempts == 0 || o.maximum_rejected_attempts == 0 ||
        o.check_interval == 0 || o.minimum_iterations == 0 || o.minimum_iterations > o.maximum_iterations ||
        o.consecutive_checks == 0 || o.stagnation_checks == 0 || o.stagnation_windows == 0)
        throw std::invalid_argument("Euler local steady iteration/attempt/check limits must be positive and ordered.");
}

std::array<double, 4> components(const ConservativeState &u)
{
    return {u.density, u.momentum_x, u.momentum_r, u.total_energy_density};
}
} // namespace

const char *to_string(const EulerLocalSteadyStatus status) noexcept
{
    switch (status)
    {
    case EulerLocalSteadyStatus::Running:
        return "RUNNING";
    case EulerLocalSteadyStatus::Converged:
        return "CONVERGED";
    case EulerLocalSteadyStatus::MaximumIterations:
        return "MAX_ITERATIONS";
    case EulerLocalSteadyStatus::WallLimit:
        return "WALL_LIMIT";
    case EulerLocalSteadyStatus::Stagnated:
        return "STAGNATED";
    case EulerLocalSteadyStatus::Diverged:
        return "DIVERGED";
    case EulerLocalSteadyStatus::NumericalFailure:
        return "NUMERICAL_FAILURE";
    }
    return "UNKNOWN";
}

AxisymmetricEulerSteadySolver::AxisymmetricEulerSteadySolver(const Mesh &mesh, const ConservativeField &initial,
                                                             CaloricallyPerfectGas gas,
                                                             EulerBoundaryConditions conditions,
                                                             EulerLocalSteadyOptions options,
                                                             const std::span<const Index> cell_regions)
    : engine_(mesh, initial, gas, std::move(conditions),
              {.cfl = options.cfl, .collect_timings = options.collect_timings}),
      options_(options), previous_(initial), local_steps_(mesh.cell_count()), region_index_(mesh.cell_count())
{
    validate(options_);
    for (Index face = 0; face < mesh.face_count(); ++face)
        if (mesh.face_adjacencies()[face].is_boundary() && engine_.face_areas_[face] > 0.0)
            boundary_faces_.push_back(face);
    if (!cell_regions.empty() && cell_regions.size() != mesh.cell_count())
        throw std::invalid_argument("Euler local steady region cardinality must match cells.");
    std::vector<Index> labels;
    if (cell_regions.empty())
        labels.push_back(0);
    else
    {
        labels.assign(cell_regions.begin(), cell_regions.end());
        std::ranges::sort(labels);
        labels.erase(std::unique(labels.begin(), labels.end()), labels.end());
    }
    regions_.resize(labels.size());
    region_sums_.resize(labels.size());
    region_volumes_.resize(labels.size());
    for (Index r = 0; r < labels.size(); ++r)
        regions_[r].region = labels[r];
    for (Index cell = 0; cell < mesh.cell_count(); ++cell)
    {
        const Index label{cell_regions.empty() ? 0 : cell_regions[cell]};
        const auto region{static_cast<Index>(std::lower_bound(labels.begin(), labels.end(), label) - labels.begin())};
        region_index_[cell] = region;
        const double v{geometry().cell_volumes()[cell]};
        region_volumes_[region] += v;
        volume_ += v;
    }
    const long double rho{options_.reference_density}, u{options_.reference_velocity},
        length{options_.reference_length};
    field_scales_ = {rho, rho * u, rho * u, rho * u * u};
    reference_time_ = length / u;
    flow_scale_ = rho * u * length * length;
    for (Index r = 0; r < regions_.size(); ++r)
    {
        regions_[r].volume = checked(region_volumes_[r]);
        for (const auto q : field_scales_)
            if (!std::isfinite(region_volumes_[r] * q / reference_time_) ||
                !(region_volumes_[r] * q / reference_time_ > 0) || !std::isfinite(q * u * length * length) ||
                !(q * u * length * length > 0))
                throw std::invalid_argument("Euler local steady normalization is not representable.");
    }
    if (!std::isfinite(volume_) || !(volume_ > 0) || !std::isfinite(flow_scale_) || !(flow_scale_ > 0))
        throw std::invalid_argument("Euler local steady volume/flow scale is invalid.");
    result_.integrals = engine_.result_.integrals;
    result_.minimum_density = engine_.result_.minimum_density;
    result_.minimum_pressure = engine_.result_.minimum_pressure;
}

bool AxisymmetricEulerSteadySolver::wall_limit()
{
    if (result_.status != EulerLocalSteadyStatus::Running)
        return false;
    if (!start_)
        start_ = Clock::now();
    if (std::chrono::duration<double>(Clock::now() - *start_).count() >= options_.maximum_wall_seconds)
    {
        result_.status = EulerLocalSteadyStatus::WallLimit;
        return true;
    }
    return false;
}

void AxisymmetricEulerSteadySolver::fail(const std::runtime_error &error)
{
    failure_ = error.what();
    result_.status = EulerLocalSteadyStatus::NumericalFailure;
}

const EulerLocalSteadyResult &AxisymmetricEulerSteadySolver::observe()
{
    if (observed_ && result_.observation_iteration == result_.iteration)
        return result_;
    (void)wall_limit();
    const auto terminal_status{result_.status};
    try
    {
        observe_current();
        if (terminal_status != EulerLocalSteadyStatus::Running)
            result_.status = terminal_status; // Refresh diagnostics, never resume a terminal solve.
    }
    catch (const std::runtime_error &error)
    {
        fail(error);
    }
    return result_;
}

void AxisymmetricEulerSteadySolver::observe_current()
{
    const auto current{engine_.evaluate_steady_residual()};
    result_.absolute_residual = current.absolute_residual;
    for (auto &sums : region_sums_)
        sums.fill(0);
    for (auto &r : regions_)
        r.maximum_cell_residual.fill(0);
    std::array<long double, 4> changes{};
    result_.maximum_cell_residual.fill(0);
    result_.maximum_cell_change.fill(0);
    for (Index cell = 0; cell < state().size(); ++cell)
    {
        const auto r{components(engine_.residual_.state(cell))}, now{components(state().state(cell))},
            old{components(previous_.state(cell))};
        const long double v{geometry().cell_volumes()[cell]};
        const auto region{region_index_[cell]};
        for (Index k = 0; k < 4; ++k)
        {
            const long double absolute{std::abs(static_cast<long double>(r.at(k)))};
            const double local{checked(absolute * reference_time_ / (v * field_scales_.at(k)))};
            region_sums_[region].at(k) += absolute;
            regions_[region].maximum_cell_residual.at(k) =
                std::max(regions_[region].maximum_cell_residual.at(k), local);
            result_.maximum_cell_residual.at(k) = std::max(result_.maximum_cell_residual.at(k), local);
            changes.at(k) += v * std::abs(static_cast<long double>(now.at(k)) - old.at(k));
            result_.maximum_cell_change.at(k) =
                std::max(result_.maximum_cell_change.at(k),
                         checked(std::abs(static_cast<long double>(now.at(k)) - old.at(k)) / field_scales_.at(k)));
        }
    }
    bool passed{observed_ && result_.iteration >= options_.minimum_iterations};
    double worst{};
    const auto boundary{engine_.boundary_outflow_};
    const std::array rates{boundary.mass, boundary.momentum_x, boundary.momentum_r, boundary.energy};
    for (Index k = 0; k < 4; ++k)
    {
        const long double q{field_scales_.at(k)};
        result_.normalized_residual.at(k) = checked(current.absolute_residual.at(k) * reference_time_ / (volume_ * q));
        result_.field_change.at(k) = checked(changes.at(k) / (volume_ * q));
        for (Index r = 0; r < regions_.size(); ++r)
        {
            regions_[r].absolute_residual.at(k) = checked(region_sums_[r].at(k));
            regions_[r].normalized_residual.at(k) =
                checked(region_sums_[r].at(k) * reference_time_ / (region_volumes_[r] * q));
            passed = passed && regions_[r].normalized_residual.at(k) <= options_.residual_tolerance;
        }
        const double residual{result_.maximum_cell_residual.at(k)};
        if (!observed_)
            initial_residual_.at(k) = residual;
        const double initial{initial_residual_.at(k)};
        result_.residual_reduction.at(k) =
            checked(static_cast<long double>(residual) / std::max(initial, options_.residual_tolerance));
        const long double source{k == 2 ? engine_.radial_source_ : 0.0};
        result_.stationary_balance.at(k) =
            checked(std::abs(static_cast<long double>(rates.at(k)) - source) /
                    (q * options_.reference_velocity * options_.reference_length * options_.reference_length));
        passed = passed && residual <= options_.residual_tolerance &&
                 result_.normalized_residual.at(k) <= options_.residual_tolerance &&
                 (initial <= options_.residual_tolerance ||
                  result_.residual_reduction.at(k) <= options_.residual_reduction) &&
                 result_.field_change.at(k) <= options_.field_change_tolerance &&
                 result_.maximum_cell_change.at(k) <= options_.field_change_tolerance &&
                 result_.stationary_balance.at(k) <= options_.balance_tolerance;
        worst = std::max(worst, residual);
        window_peaks_.at(k) = std::max(window_peaks_.at(k), residual);
    }
    const double initial_worst{*std::max_element(initial_residual_.begin(), initial_residual_.end())};
    if (worst > static_cast<long double>(options_.maximum_residual_growth) *
                    std::max(initial_worst, options_.residual_tolerance))
        result_.status = EulerLocalSteadyStatus::Diverged;
    result_.mass_flow_in = current.mass_flow_in;
    result_.mass_flow_out = current.mass_flow_out;
    const long double in{current.mass_flow_in}, out{current.mass_flow_out};
    result_.mass_flow_relative_imbalance = checked(in + out == 0 ? 0 : std::abs(out - in) / (in + out));
    result_.flow_change = checked(std::max(std::abs(in - previous_in_), std::abs(out - previous_out_)) / flow_scale_);
    passed = passed && result_.flow_change <= options_.flow_change_tolerance;
    stable_checks_ = passed ? stable_checks_ + 1 : 0;
    bool stagnated{};
    if (observed_ && ++window_checks_ >= options_.stagnation_checks)
    {
        // A rising envelope can be one phase of a damped acoustic oscillation.
        // Require sustained lack of progress; keep the reference until a
        // significant decrease so small improvements can accumulate.
        bool progress{!have_reference_window_};
        for (Index k = 0; k < window_peaks_.size(); ++k)
        {
            auto &reference{reference_window_peaks_.at(k)};
            const double peak{window_peaks_.at(k)};
            if (!have_reference_window_ || reference <= options_.residual_tolerance)
                reference = peak; // Track an equation which becomes unresolved, without calling growth progress.
            else if (peak < (1.0 - options_.minimum_residual_progress) * reference)
            {
                reference = peak;
                progress = true;
            }
        }
        if (progress)
            unproductive_windows_ = 0;
        else if (unproductive_windows_ < options_.stagnation_windows)
            ++unproductive_windows_;
        stagnated = worst > options_.residual_tolerance && unproductive_windows_ >= options_.stagnation_windows;
        window_peaks_.fill(0);
        window_checks_ = 0;
        have_reference_window_ = true;
    }
    if (result_.status == EulerLocalSteadyStatus::Running)
    {
        if (stable_checks_ >= options_.consecutive_checks)
            result_.status = EulerLocalSteadyStatus::Converged;
        else if (result_.iteration >= options_.maximum_iterations)
            result_.status = EulerLocalSteadyStatus::MaximumIterations;
        else if (stagnated)
            result_.status = EulerLocalSteadyStatus::Stagnated;
    }
    result_.boundary_outflow_rate = boundary;
    result_.radial_pressure_source_rate = engine_.radial_source_;
    result_.observation_iteration = result_.iteration;
    previous_in_ = current.mass_flow_in;
    previous_out_ = current.mass_flow_out;
    for (Index cell = 0; cell < state().size(); ++cell)
        previous_.set_state(cell, state().state(cell));
    observed_ = true;
}

bool AxisymmetricEulerSteadySolver::candidate_boundaries_compatible(const ConservativeField &candidate,
                                                                    const double alpha)
{
    for (const Index face : boundary_faces_)
    {
        try
        {
            engine_.validate_boundary_state(candidate, face);
        }
        catch (const EulerCharacteristicCompatibilityError &error)
        {
            std::ostringstream message;
            message.precision(17);
            message << "Euler local steady incompatible candidate BC: face=" << face
                    << ", cell=" << engine_.mesh_->face_adjacencies()[face].owner
                    << ", group=" << engine_.mesh_->face_boundary_ids()[face] << ", alpha=" << alpha << ": "
                    << error.what();
            failure_ = message.str();
            return false;
        }
    }
    return true;
}

void AxisymmetricEulerSteadySolver::update_candidate()
{
    for (Index cell = 0; cell < state().size(); ++cell)
    {
        const double volume{geometry().cell_volumes()[cell]}, wave{engine_.wave_sums_[cell]};
        if (!std::isfinite(volume) || !(volume > 0) || !std::isfinite(wave) || !(wave > 0))
            throw std::runtime_error("Euler local steady requires positive finite cell volume and spectral sum.");
        local_steps_[cell] = checked(static_cast<long double>(options_.cfl) * volume / wave);
        if (!(local_steps_[cell] > 0))
            throw std::runtime_error("Euler local steady pseudo-step underflows to zero.");
    }
    auto &candidate{engine_.states_.at(1 - engine_.current_)};
    double alpha{1.0};
    for (Index attempt = 0; attempt < options_.maximum_step_attempts; ++attempt)
    {
        bool admissible{true};
        for (Index cell = 0; cell < state().size(); ++cell)
        {
            const auto u{components(state().state(cell))}, r{components(engine_.residual_.state(cell))};
            const double dt{checked(static_cast<long double>(alpha) * local_steps_[cell])};
            if (!(dt > 0))
                throw std::runtime_error("Euler local steady effective pseudo-step underflows to zero.");
            const long double scale{static_cast<long double>(dt) / geometry().cell_volumes()[cell]};
            try
            {
                const ConservativeState value{checked(u[0] + scale * r[0]), checked(u[1] + scale * r[1]),
                                              checked(u[2] + scale * r[2]), checked(u[3] + scale * r[3])};
                (void)engine_.gas_.evaluate(value);
                candidate.set_state(cell, value);
            }
            catch (const std::runtime_error &error)
            {
                std::ostringstream message;
                message.precision(17);
                message << "Euler local steady inadmissible candidate: cell=" << cell << ", alpha=" << alpha
                        << ", dtau=" << dt << ": " << error.what();
                failure_ = message.str();
                admissible = false;
                break;
            }
        }
        if (admissible && candidate_boundaries_compatible(candidate, alpha))
        {
            const auto summary{engine_.summarize(candidate)};
            engine_.current_ = 1 - engine_.current_;
            engine_.assembled_time_step_.reset();
            ++result_.iteration;
            result_.step_scale = alpha;
            result_.minimum_local_step = std::numeric_limits<double>::infinity();
            result_.maximum_local_step = 0;
            for (auto &dt : local_steps_)
            {
                dt = checked(static_cast<long double>(alpha) * dt);
                result_.minimum_local_step = std::min(result_.minimum_local_step, dt);
                result_.maximum_local_step = std::max(result_.maximum_local_step, dt);
            }
            result_.integrals = summary.integrals;
            result_.minimum_density = summary.minimum_density;
            result_.minimum_pressure = summary.minimum_pressure;
            failure_.clear();
            return;
        }
        ++result_.rejected_attempts;
        if (result_.rejected_attempts >= options_.maximum_rejected_attempts)
            throw std::runtime_error("Euler local steady rejection limit. " + failure_);
        alpha *= 0.5;
    }
    throw std::runtime_error("Euler local steady exhausted candidate attempts. " + failure_);
}

const EulerLocalSteadyResult &AxisymmetricEulerSteadySolver::iterate()
{
    if (result_.status != EulerLocalSteadyStatus::Running)
        return result_;
    if (wall_limit())
        return observe();
    if (!observed_)
        observe();
    if (result_.status != EulerLocalSteadyStatus::Running)
        return result_;
    try
    {
        (void)engine_.prepare_assembly();
        const auto begin{options_.collect_timings ? Clock::now() : Clock::time_point{}};
        update_candidate();
        if (options_.collect_timings)
            engine_.timings_.time_update_seconds += std::chrono::duration<double>(Clock::now() - begin).count();
    }
    catch (const std::runtime_error &error)
    {
        fail(error);
        return result_;
    }
    if (result_.iteration % options_.check_interval == 0 || result_.iteration >= options_.maximum_iterations)
        observe();
    (void)wall_limit();
    if (result_.status == EulerLocalSteadyStatus::WallLimit)
        observe();
    return result_;
}

const EulerLocalSteadyResult &AxisymmetricEulerSteadySolver::solve()
{
    while (result_.status == EulerLocalSteadyStatus::Running)
        iterate();
    return result_;
}
} // namespace cfd
