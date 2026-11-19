#include "cfd/compressible/AxisymmetricEulerSolver.hpp"

#include "cfd/compressible/EulerCharacteristicBoundaryFlux.hpp"
#include "cfd/compressible/EulerFarfieldBoundaryFlux.hpp"
#include "cfd/compressible/EulerFlux.hpp"
#include "cfd/compressible/EulerSlipWallFlux.hpp"
#include "cfd/mesh/Mesh.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace cfd
{
namespace
{
using Clock = std::chrono::steady_clock;

double finite_value(const long double value)
{
    if (!std::isfinite(value) || std::abs(value) > std::numeric_limits<double>::max())
        throw std::runtime_error("Axisymmetric Euler arithmetic is not representable as finite double.");
    return static_cast<double>(value);
}

void validate_options(const AxisymmetricEulerOptions &options)
{
    if (!std::isfinite(options.cfl) || !(options.cfl > 0.0) || options.cfl > 1.0 ||
        options.maximum_step_attempts == 0 || options.maximum_steps == 0 || options.maximum_rejected_steps == 0)
        throw std::invalid_argument("Euler requires 0 < CFL <= 1 and positive step/attempt limits.");
}

void validate_boundary_geometry(const Mesh &mesh, const EulerBoundaryConditions &conditions)
{
    if (conditions.size() != mesh.boundary_groups().size())
        throw std::invalid_argument("Euler boundary condition count must match Mesh boundary groups.");
    for (Index face = 0; face < mesh.face_count(); ++face)
    {
        if (!mesh.face_adjacencies()[face].is_boundary())
            continue;
        const auto nodes{mesh.faces()[face].node_ids};
        const bool on_axis{mesh.nodes()[nodes[0]].y == 0.0 && mesh.nodes()[nodes[1]].y == 0.0};
        if ((conditions[mesh.face_boundary_ids()[face]] == EulerBoundaryConditionType::Axis) != on_axis)
            throw std::invalid_argument("Euler Axis BC must identify exactly the faces on r=0 (face " +
                                        std::to_string(face) + ").");
    }
}

RusanovFluxResult boundary_flux(const ConservativeState &inside, const Vector2 normal,
                                const EulerBoundaryCondition &condition, const CaloricallyPerfectGas &gas,
                                const double group_area)
{
    switch (condition.type())
    {
    case EulerBoundaryConditionType::SlipWall: {
        const auto wall{compute_euler_slip_wall_flux(inside, normal, gas)};
        return {wall.flux, wall.maximum_wave_speed};
    }
    case EulerBoundaryConditionType::Transmissive:
        return compute_rusanov_flux(inside, inside, normal, gas);
    case EulerBoundaryConditionType::Axis:
        break; // Axis faces must have been skipped before requesting any flux.
    case EulerBoundaryConditionType::Farfield: {
        const auto farfield{compute_euler_farfield_boundary_flux(inside, normal, condition, gas)};
        return {farfield.flux, farfield.maximum_wave_speed};
    }
    case EulerBoundaryConditionType::TotalConditionsInlet:
    case EulerBoundaryConditionType::MassFlowInlet:
    case EulerBoundaryConditionType::SupersonicInlet:
    case EulerBoundaryConditionType::PressureOutlet:
    case EulerBoundaryConditionType::SupersonicOutlet: {
        const auto open{compute_euler_characteristic_boundary_flux(inside, normal, condition, gas, group_area)};
        return {open.flux, open.maximum_wave_speed};
    }
    }
    throw std::runtime_error("Unexpected zero-area Axis flux evaluation.");
}

RusanovFluxResult timed_boundary_flux(const ConservativeState &inside, const Vector2 normal,
                                      const EulerBoundaryCondition &condition, const CaloricallyPerfectGas &gas,
                                      const double group_area, AxisymmetricEulerTimings *timings)
{
    const auto start{timings != nullptr ? Clock::now() : Clock::time_point{}};
    const auto flux{boundary_flux(inside, normal, condition, gas, group_area)};
    if (timings != nullptr)
        timings->boundary_flux_seconds += std::chrono::duration<double>(Clock::now() - start).count();
    return flux;
}

void add_flux(ConservativeField &residual, const Index cell, const EulerFlux &flux, const double signed_area)
{
    auto value{residual.state(cell)};
    value.density =
        finite_value(static_cast<long double>(value.density) + signed_area * static_cast<long double>(flux.mass));
    value.momentum_x = finite_value(static_cast<long double>(value.momentum_x) +
                                    signed_area * static_cast<long double>(flux.momentum_x));
    value.momentum_r = finite_value(static_cast<long double>(value.momentum_r) +
                                    signed_area * static_cast<long double>(flux.momentum_r));
    value.total_energy_density = finite_value(static_cast<long double>(value.total_energy_density) +
                                              signed_area * static_cast<long double>(flux.energy));
    residual.set_state(cell, value);
}

std::string invalid_candidate(const Index cell, const std::array<long double, 4> &value, const double dt,
                              const std::string &reason)
{
    std::ostringstream message;
    message.precision(17);
    message << "Euler inadmissible candidate: cell=" << cell << ", dt=" << dt << ", U=[" << value[0] << ", " << value[1]
            << ", " << value[2] << ", " << value[3] << "]: " << reason;
    return message.str();
}

EulerIntegralQuantities balance_error(const EulerIntegralQuantities &before, const EulerIntegralQuantities &after,
                                      const EulerIntegralQuantities &outflow, const double source, const double dt)
{
    return {
        finite_value(static_cast<long double>(after.mass) - before.mass + static_cast<long double>(dt) * outflow.mass),
        finite_value(static_cast<long double>(after.momentum_x) - before.momentum_x +
                     static_cast<long double>(dt) * outflow.momentum_x),
        finite_value(static_cast<long double>(after.momentum_r) - before.momentum_r +
                     static_cast<long double>(dt) * (static_cast<long double>(outflow.momentum_r) - source)),
        finite_value(static_cast<long double>(after.energy) - before.energy +
                     static_cast<long double>(dt) * outflow.energy)};
}

} // namespace

AxisymmetricEulerSolver::AxisymmetricEulerSolver(const Mesh &mesh, const ConservativeField &initial_state,
                                                 CaloricallyPerfectGas gas, EulerBoundaryConditions boundary_conditions,
                                                 AxisymmetricEulerOptions options)
    : mesh_(&mesh), geometry_(mesh), gas_(gas), boundary_conditions_(std::move(boundary_conditions)), options_(options),
      states_{ConservativeField{initial_state}, ConservativeField{mesh.cell_count()}}, residual_(mesh.cell_count()),
      wave_sums_(mesh.cell_count()), face_areas_(mesh.face_count()),
      boundary_group_areas_(mesh.boundary_groups().size())
{
    validate_options(options_);
    validate_boundary_geometry(mesh, boundary_conditions_);
    if (initial_state.size() != mesh.cell_count())
        throw std::invalid_argument("Euler initial state cardinality must match Mesh cells.");
    for (Index face = 0; face < mesh.face_count(); ++face)
    {
        const auto area{geometry_.face_area_vectors()[face]};
        face_areas_[face] = std::hypot(area.x, area.y);
        if (!std::isfinite(face_areas_[face]))
            throw std::runtime_error("Euler revolved face area is not finite.");
        if (mesh.face_adjacencies()[face].is_boundary())
        {
            const auto group{mesh.face_boundary_ids()[face]};
            if (boundary_conditions_[group] == EulerBoundaryConditionType::MassFlowInlet)
                boundary_group_areas_[group] =
                    finite_value(static_cast<long double>(boundary_group_areas_[group]) + face_areas_[face]);
        }
    }
    for (Index group = 0; group < boundary_conditions_.size(); ++group)
        if (boundary_conditions_[group] == EulerBoundaryConditionType::MassFlowInlet &&
            !(boundary_group_areas_[group] > 0.0))
            throw std::invalid_argument("Euler mass inlet group must have a positive revolution area.");
    result_ = summarize(state());
    result_.cfl = options_.cfl;
}

void AxisymmetricEulerSolver::validate_boundary_state(const ConservativeField &field, const Index face) const
{
    const auto area{mesh_->face_area_vectors()[face]};
    const double length{mesh_->face_lengths()[face]};
    const Vector2 normal{area.x / length, area.y / length};
    const auto group{mesh_->face_boundary_ids()[face]};
    (void)boundary_flux(field.state(mesh_->face_adjacencies()[face].owner), normal,
                        boundary_conditions_.condition(group), gas_, boundary_group_areas_[group]);
}

AxisymmetricEulerResult AxisymmetricEulerSolver::summarize(const ConservativeField &field) const
{
    AxisymmetricEulerResult summary;
    summary.minimum_density = std::numeric_limits<double>::infinity();
    summary.minimum_pressure = std::numeric_limits<double>::infinity();
    long double mass{}, mx{}, mr{}, energy{};
    for (Index cell = 0; cell < mesh_->cell_count(); ++cell)
    {
        const auto value{field.state(cell)};
        const auto primitive{gas_.evaluate(value).primitive};
        summary.minimum_density = std::min(summary.minimum_density, primitive.density);
        summary.maximum_density = std::max(summary.maximum_density, primitive.density);
        summary.minimum_pressure = std::min(summary.minimum_pressure, primitive.pressure);
        summary.maximum_pressure = std::max(summary.maximum_pressure, primitive.pressure);
        const long double volume{geometry_.cell_volumes()[cell]};
        mass += volume * value.density;
        mx += volume * value.momentum_x;
        mr += volume * value.momentum_r;
        energy += volume * value.total_energy_density;
    }
    summary.integrals = {finite_value(mass), finite_value(mx), finite_value(mr), finite_value(energy)};
    return summary;
}

AxisymmetricEulerRestartState AxisymmetricEulerSolver::restart_state() const
{
    return {ConservativeField{state()}, result_};
}

AxisymmetricEulerSolver::AxisymmetricEulerSolver(const Mesh &mesh, const AxisymmetricEulerRestartState &restart,
                                                 CaloricallyPerfectGas gas, EulerBoundaryConditions conditions,
                                                 AxisymmetricEulerOptions options)
    : AxisymmetricEulerSolver(mesh, restart.conservative, gas, std::move(conditions), options)
{
    const auto &d{restart.diagnostics};
    for (const double value :
         {d.time, d.time_step, d.cfl, d.radial_pressure_source_rate, d.boundary_outflow_rate.mass,
          d.boundary_outflow_rate.momentum_x, d.boundary_outflow_rate.momentum_r, d.boundary_outflow_rate.energy,
          d.balance_error.mass, d.balance_error.momentum_x, d.balance_error.momentum_r, d.balance_error.energy})
        if (!std::isfinite(value))
            throw std::invalid_argument("Euler restart diagnostics must be finite.");
    if (d.time < 0 || d.time_step < 0 || d.time_step > d.time || d.cfl != options.cfl ||
        (d.accepted_steps == 0 ? (d.time != 0 || d.time_step != 0) : !(d.time > 0 && d.time_step > 0)) ||
        d.accepted_steps > options.maximum_steps || d.rejected_steps > options.maximum_rejected_steps ||
        d.integrals != result_.integrals || d.minimum_density != result_.minimum_density ||
        d.maximum_density != result_.maximum_density || d.minimum_pressure != result_.minimum_pressure ||
        d.maximum_pressure != result_.maximum_pressure)
        throw std::invalid_argument("Euler restart has inconsistent times, counters or accepted-state diagnostics.");
    if (d.accepted_steps == 0 && (d.boundary_outflow_rate != EulerIntegralQuantities{} ||
                                  d.radial_pressure_source_rate != 0 || d.balance_error != EulerIntegralQuantities{}))
        throw std::invalid_argument("Euler restart initial state cannot contain accepted-step rates or defects.");
    result_ = d;
}

double AxisymmetricEulerSolver::assemble()
{
    std::fill(wave_sums_.begin(), wave_sums_.end(), 0.0);
    long double source{};
    for (Index cell = 0; cell < mesh_->cell_count(); ++cell)
    {
        const double pressure{gas_.evaluate(state().state(cell)).primitive.pressure};
        const double radial{
            finite_value(static_cast<long double>(geometry_.cell_radial_source_coefficients()[cell]) * pressure)};
        residual_.set_state(cell, {0.0, 0.0, radial, 0.0});
        source += radial;
    }
    radial_source_ = finite_value(source);
    long double mass{}, mx{}, mr{}, energy{}, mass_in{}, mass_out{};
    for (Index face = 0; face < mesh_->face_count(); ++face)
    {
        const double area{face_areas_[face]};
        if (area == 0.0)
            continue; // No flux, wave contribution, or axis-area normalization.
        const auto adjacency{mesh_->face_adjacencies()[face]};
        const auto area_2d{mesh_->face_area_vectors()[face]};
        const double length{mesh_->face_lengths()[face]};
        const Vector2 normal{area_2d.x / length, area_2d.y / length};
        const auto left{state().state(adjacency.owner)};
        const auto numerical{adjacency.is_boundary()
                                 ? timed_boundary_flux(left, normal,
                                                       boundary_conditions_.condition(mesh_->face_boundary_ids()[face]),
                                                       gas_, boundary_group_areas_[mesh_->face_boundary_ids()[face]],
                                                       options_.collect_timings ? &timings_ : nullptr)
                                 : compute_rusanov_flux(left, state().state(adjacency.neighbor), normal, gas_)};
        const double wave{finite_value(static_cast<long double>(area) * numerical.maximum_wave_speed)};
        wave_sums_[adjacency.owner] = finite_value(static_cast<long double>(wave_sums_[adjacency.owner]) + wave);
        add_flux(residual_, adjacency.owner, numerical.flux, -area);
        if (!adjacency.is_boundary())
        {
            add_flux(residual_, adjacency.neighbor, numerical.flux, area);
            wave_sums_[adjacency.neighbor] =
                finite_value(static_cast<long double>(wave_sums_[adjacency.neighbor]) + wave);
        }
        else
        {
            const long double boundary_mass{static_cast<long double>(area) * numerical.flux.mass};
            mass += boundary_mass;
            mass_in += std::max(-boundary_mass, 0.0L);
            mass_out += std::max(boundary_mass, 0.0L);
            mx += static_cast<long double>(area) * numerical.flux.momentum_x;
            mr += static_cast<long double>(area) * numerical.flux.momentum_r;
            energy += static_cast<long double>(area) * numerical.flux.energy;
        }
    }
    boundary_outflow_ = {finite_value(mass), finite_value(mx), finite_value(mr), finite_value(energy)};
    mass_flow_in_ = finite_value(mass_in);
    mass_flow_out_ = finite_value(mass_out);
    long double dt{std::numeric_limits<long double>::infinity()};
    for (Index cell = 0; cell < mesh_->cell_count(); ++cell)
    {
        if (!(wave_sums_[cell] > 0.0))
            throw std::runtime_error("Euler cell has no finite positive wave scale.");
        dt = std::min(dt, static_cast<long double>(options_.cfl) * geometry_.cell_volumes()[cell] / wave_sums_[cell]);
    }
    return finite_value(dt);
}

double AxisymmetricEulerSolver::prepare_assembly()
{
    if (!assembled_time_step_)
    {
        const auto start{options_.collect_timings ? Clock::now() : Clock::time_point{}};
        assembled_time_step_ = assemble();
        if (options_.collect_timings)
            timings_.assembly_seconds += std::chrono::duration<double>(Clock::now() - start).count();
    }
    return *assembled_time_step_;
}

EulerCflEvaluation AxisymmetricEulerSolver::evaluate_cfl()
{
    const double time_step{prepare_assembly()};
    return {time_step, wave_sums_};
}

EulerResidualEvaluation AxisymmetricEulerSolver::evaluate_steady_residual()
{
    (void)prepare_assembly();
    std::array<long double, 4> absolute{};
    for (Index cell = 0; cell < state().size(); ++cell)
    {
        const auto value{residual_.state(cell)};
        absolute[0] += std::abs(static_cast<long double>(value.density));
        absolute[1] += std::abs(static_cast<long double>(value.momentum_x));
        absolute[2] += std::abs(static_cast<long double>(value.momentum_r));
        absolute[3] += std::abs(static_cast<long double>(value.total_energy_density));
    }
    return {
        {finite_value(absolute[0]), finite_value(absolute[1]), finite_value(absolute[2]), finite_value(absolute[3])},
        mass_flow_in_,
        mass_flow_out_};
}

const AxisymmetricEulerResult &AxisymmetricEulerSolver::step(const double final_time)
{
    if (!std::isfinite(final_time) || final_time < result_.time)
        throw std::invalid_argument("Euler final time must be finite and not earlier than current time.");
    if (final_time == result_.time)
        return result_;
    if (result_.accepted_steps >= options_.maximum_steps)
        throw std::runtime_error("Euler maximum accepted-step count reached.");
    if (result_.rejected_steps >= options_.maximum_rejected_steps)
        throw std::runtime_error("Euler maximum rejected-step count reached.");
    const double remaining{final_time - result_.time};
    double dt{std::min(prepare_assembly(), remaining)};
    const auto update_start{options_.collect_timings ? Clock::now() : Clock::time_point{}};
    auto &candidate{states_.at(1 - current_)};
    std::string failure;
    for (Index attempt = 0; attempt < options_.maximum_step_attempts; ++attempt)
    {
        if (!(dt > 0.0) || !std::isfinite(dt) || !(result_.time + dt > result_.time))
            throw std::runtime_error("Euler time step cannot advance physical time. " + failure);
        bool admissible{true};
        for (Index cell = 0; cell < mesh_->cell_count(); ++cell)
        {
            const auto old{state().state(cell)}, rate{residual_.state(cell)};
            const long double scale{static_cast<long double>(dt) / geometry_.cell_volumes()[cell]};
            const std::array<long double, 4> attempted{
                old.density + scale * rate.density, old.momentum_x + scale * rate.momentum_x,
                old.momentum_r + scale * rate.momentum_r, old.total_energy_density + scale * rate.total_energy_density};
            ConservativeState value{};
            try
            {
                value = {finite_value(attempted[0]), finite_value(attempted[1]), finite_value(attempted[2]),
                         finite_value(attempted[3])};
                (void)gas_.evaluate(value);
            }
            catch (const std::runtime_error &error)
            {
                failure = invalid_candidate(cell, attempted, dt, error.what());
                admissible = false;
                break;
            }
            candidate.set_state(cell, value);
        }
        if (admissible)
        {
            auto next{summarize(candidate)};
            next.time = dt == remaining ? final_time : result_.time + dt;
            next.time_step = dt;
            next.cfl = options_.cfl;
            next.accepted_steps = result_.accepted_steps + 1;
            next.rejected_steps = result_.rejected_steps;
            next.boundary_outflow_rate = boundary_outflow_;
            next.radial_pressure_source_rate = radial_source_;
            next.balance_error =
                balance_error(result_.integrals, next.integrals, boundary_outflow_, radial_source_, dt);
            current_ = 1 - current_;
            result_ = next;
            assembled_time_step_.reset();
            if (options_.collect_timings)
                timings_.time_update_seconds += std::chrono::duration<double>(Clock::now() - update_start).count();
            return result_;
        }
        ++result_.rejected_steps;
        if (result_.rejected_steps >= options_.maximum_rejected_steps)
            throw std::runtime_error("Euler maximum rejected-step count reached. " + failure);
        dt *= 0.5;
    }
    throw std::runtime_error("Euler exhausted time-step attempts. " + failure);
}

const AxisymmetricEulerResult &AxisymmetricEulerSolver::advance_to(const double final_time)
{
    // step() also validates a target equal to or earlier than the current time.
    step(final_time);
    while (result_.time < final_time)
        step(final_time);
    return result_;
}

} // namespace cfd
