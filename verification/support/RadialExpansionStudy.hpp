#pragma once

#include "EulerTube.hpp"
#include "RadialExpansionExact.hpp"
#include "cfd/compressible/AxisymmetricEulerSolver.hpp"
#include "cfd/mesh/MeshBuilder.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace cfd::verification
{

struct RadialGrid
{
    Index nx{5};
    Index nr{20};
    double length{1.0};
    double radius{1.0};
};

inline MeshBuildResult build_radial_grid(const RadialGrid &grid)
{
    if (!std::isfinite(grid.length) || !(grid.length > 0.0) || !std::isfinite(grid.radius) || !(grid.radius > 0.0))
        throw std::invalid_argument("Radial grid requires finite positive length and radius.");
    auto raw{make_euler_tube(grid.nx, grid.nr)};
    for (auto &node : raw.nodes)
    {
        node.x *= grid.length;
        node.y *= grid.radius;
    }
    return build_mesh(std::move(raw));
}

inline std::array<double, 2> radial_bounds(const RadialGrid &grid, const Index cell)
{
    const Index row{cell / grid.nx};
    return {grid.radius * (static_cast<double>(row) / static_cast<double>(grid.nr)),
            grid.radius * (static_cast<double>(row + 1) / static_cast<double>(grid.nr))};
}

inline double radial_cell_volume(const RadialGrid &grid, const Index cell)
{
    const auto bounds{radial_bounds(grid, cell)};
    return std::numbers::pi * (grid.length / static_cast<double>(grid.nx)) * (bounds[1] - bounds[0]) *
           (bounds[1] + bounds[0]);
}

inline ConservativeField radial_initial_field(const RadialGrid &grid, const RadialExpansionExact &exact)
{
    ConservativeField initial{grid.nx * grid.nr};
    for (Index cell = 0; cell < initial.size(); ++cell)
    {
        const auto bounds{radial_bounds(grid, cell)};
        initial.set_state(cell, exact.cell_average(bounds[0], bounds[1], 0.0));
    }
    return initial;
}

inline EulerBoundaryConditions radial_conditions()
{
    using Type = EulerBoundaryConditionType;
    return {4, {Type::Axis, Type::Transmissive, Type::Transmissive, Type::Transmissive}};
}

struct RadialErrorNorms
{
    double l1{};
    double rms{};
    double linf{};
};

/// Verification-only volume norms of rho, rho*ur and rhoE. No production API.
struct RadialErrorAccumulator
{
    long double volume{};
    std::array<long double, 3> absolute{}, squared{};
    std::array<double, 3> maximum{};

    void add(const double weight, const ConservativeState &numerical, const ConservativeState &reference)
    {
        const std::array errors{numerical.density - reference.density, numerical.momentum_r - reference.momentum_r,
                                numerical.total_energy_density - reference.total_energy_density};
        if (!std::isfinite(weight) || !(weight > 0.0))
            throw std::runtime_error("Radial error norms require finite positive volumes.");
        volume += weight;
        for (Index component = 0; component < errors.size(); ++component)
        {
            const double error{errors.at(component)};
            if (!std::isfinite(error))
                throw std::runtime_error("Radial error norm encountered a non-finite error.");
            absolute.at(component) += weight * std::abs(static_cast<long double>(error));
            squared.at(component) += weight * static_cast<long double>(error) * error;
            maximum.at(component) = std::max(maximum.at(component), std::abs(error));
        }
    }

    [[nodiscard]] std::array<RadialErrorNorms, 3> finish() const
    {
        if (!std::isfinite(volume) || !(volume > 0.0L))
            throw std::runtime_error("Radial error norm has an invalid total volume.");
        std::array<RadialErrorNorms, 3> result{};
        for (Index component = 0; component < result.size(); ++component)
        {
            result.at(component) = {static_cast<double>(absolute.at(component) / volume),
                                    static_cast<double>(std::sqrt(squared.at(component) / volume)),
                                    maximum.at(component)};
            if (!std::isfinite(result.at(component).l1) || !std::isfinite(result.at(component).rms))
                throw std::runtime_error("Radial error norm is not representable.");
        }
        return result;
    }
};

inline std::array<RadialErrorNorms, 3> radial_errors(const RadialGrid &grid, const ConservativeField &state,
                                                     const RadialExpansionExact &exact, const double time)
{
    RadialErrorAccumulator accumulator;
    for (Index cell = 0; cell < state.size(); ++cell)
    {
        const auto bounds{radial_bounds(grid, cell)};
        accumulator.add(radial_cell_volume(grid, cell), state.state(cell),
                        exact.cell_average(bounds[0], bounds[1], time));
    }
    return accumulator.finish();
}

inline std::array<RadialErrorNorms, 3> radial_differences(const RadialGrid &grid, const ConservativeField &first,
                                                          const ConservativeField &second)
{
    if (first.size() != grid.nx * grid.nr || second.size() != first.size())
        throw std::invalid_argument("Radial field comparison has incorrect cardinality.");
    RadialErrorAccumulator accumulator;
    for (Index cell = 0; cell < first.size(); ++cell)
        accumulator.add(radial_cell_volume(grid, cell), first.state(cell), second.state(cell));
    return accumulator.finish();
}

inline double independent_radial_source(const RadialGrid &grid, const ConservativeField &field, const double gamma)
{
    long double source{};
    for (Index cell = 0; cell < field.size(); ++cell)
    {
        const auto u{field.state(cell)};
        const long double kinetic{(static_cast<long double>(u.momentum_x) * u.momentum_x +
                                   static_cast<long double>(u.momentum_r) * u.momentum_r) /
                                  (2.0L * u.density)};
        const long double pressure{(gamma - 1.0L) * (u.total_energy_density - kinetic)};
        const auto bounds{radial_bounds(grid, cell)};
        source +=
            2.0L * std::numbers::pi * (grid.length / static_cast<double>(grid.nx)) * (bounds[1] - bounds[0]) * pressure;
    }
    return static_cast<double>(source);
}

struct RadialRun
{
    ConservativeField field; ///< One final copy, never a copy per step.
    AxisymmetricEulerResult result;
    std::array<RadialErrorNorms, 3> errors;
    std::array<double, 4> global_balance_error{}, maximum_step_balance_error{}; ///< Scale-relative, per component.
    double source_relative_error{};
    double initial_boundary_mass_relative_error{};
    double pressure_point_rms{}, pressure_projected_rms{};
    double axis_density_error{}, outer_density_error{}, maximum_axial_velocity{};
};

inline RadialRun run_radial_expansion(const RadialGrid &grid, const RadialExpansionExact &exact,
                                      const double final_time, const double cfl)
{
    if (!std::isfinite(final_time) || !(final_time > 0.0))
        throw std::invalid_argument("Radial study requires finite positive final time.");
    const auto built{build_radial_grid(grid)};
    const auto initial{radial_initial_field(grid, exact)};
    const CaloricallyPerfectGas gas{exact.gamma(), 1.0};
    AxisymmetricEulerSolver solver{built.mesh, initial, gas, radial_conditions(), {.cfl = cfl}};
    const long double volume{std::numbers::pi * grid.length * grid.radius * grid.radius};
    const auto average{exact.cell_average(0.0, grid.radius, 0.0)};
    const std::array<long double, 4> initial_integrals{volume * average.density, 0.0L, volume * average.momentum_r,
                                                       volume * average.total_energy_density};
    const auto values{exact.values(0.0)};
    const long double momentum_scale{
        volume * values.density *
        (values.rate * grid.radius + std::sqrt(exact.gamma() * values.pressure / values.density))};
    const std::array<long double, 4> scales{initial_integrals[0], momentum_scale, momentum_scale, initial_integrals[3]};
    std::array<long double, 4> increments{};
    std::array<double, 4> step_errors{};
    double source_error{};
    // The first outflow diagnoses Transmissive against the exact boundary trace.
    const double exact_mass_rate{2.0 * std::numbers::pi * grid.radius * grid.length * values.density * values.rate *
                                 grid.radius};
    double initial_mass_error{};
    while (solver.diagnostics().time < final_time)
    {
        const double source{independent_radial_source(grid, solver.state(), exact.gamma())};
        const auto &step{solver.step(final_time)};
        if (step.accepted_steps == 1 && exact_mass_rate != 0.0)
            initial_mass_error = step.boundary_outflow_rate.mass / exact_mass_rate - 1.0;
        source_error = std::max(source_error, std::abs(step.radial_pressure_source_rate / source - 1.0));
        const long double dt{step.time_step};
        const std::array<long double, 4> rates{-step.boundary_outflow_rate.mass, -step.boundary_outflow_rate.momentum_x,
                                               source - static_cast<long double>(step.boundary_outflow_rate.momentum_r),
                                               -step.boundary_outflow_rate.energy};
        const std::array defects{step.balance_error.mass, step.balance_error.momentum_x, step.balance_error.momentum_r,
                                 step.balance_error.energy};
        for (Index component = 0; component < increments.size(); ++component)
        {
            increments.at(component) += dt * rates.at(component);
            step_errors.at(component) = std::max(
                step_errors.at(component), static_cast<double>(std::abs(defects.at(component) / scales.at(component))));
        }
    }
    const auto result{solver.diagnostics()};
    const std::array totals{result.integrals.mass, result.integrals.momentum_x, result.integrals.momentum_r,
                            result.integrals.energy};
    std::array<double, 4> balance{};
    for (Index component = 0; component < balance.size(); ++component)
        balance.at(component) = static_cast<double>(
            (totals.at(component) - initial_integrals.at(component) - increments.at(component)) / scales.at(component));
    long double point_pressure_error{}, projected_pressure_error{};
    const auto final_values{exact.values(final_time)};
    double axis_error{}, outer_error{}, max_ux{};
    for (Index cell = 0; cell < initial.size(); ++cell)
    {
        const auto bounds{radial_bounds(grid, cell)};
        const auto primitive{gas.evaluate(solver.state().state(cell)).primitive};
        const long double weight{radial_cell_volume(grid, cell)};
        const long double point_error{primitive.pressure - final_values.pressure};
        const long double projected_error{primitive.pressure -
                                          exact.projected_pressure(bounds[0], bounds[1], final_time)};
        point_pressure_error += weight * point_error * point_error;
        projected_pressure_error += weight * projected_error * projected_error;
        const double density_error{std::abs(primitive.density - final_values.density)};
        if (cell < grid.nx)
            axis_error = std::max(axis_error, density_error);
        if (cell >= (grid.nr - 1) * grid.nx)
            outer_error = std::max(outer_error, density_error);
        max_ux = std::max(max_ux, std::abs(primitive.velocity_x));
    }
    return {ConservativeField{solver.state()},
            result,
            radial_errors(grid, solver.state(), exact, final_time),
            balance,
            step_errors,
            source_error,
            initial_mass_error,
            static_cast<double>(std::sqrt(point_pressure_error / volume)),
            static_cast<double>(std::sqrt(projected_pressure_error / volume)),
            axis_error,
            outer_error,
            max_ux};
}

struct RadialBoundaryComparison
{
    std::array<RadialErrorNorms, 3> difference, original_error, extended_error;
    Index matched_steps{};
};

/// Same dr/dx and identical accepted time steps in the common domain r<=R.
/// The extended outer boundary at 2R is initially supersonic outflow here;
/// no analytic incoming characteristic enters there. This relocates the
/// approximate boundary, not replaces it by an exact analytical condition.
inline RadialBoundaryComparison compare_radial_outer_boundary(const Index nr)
{
    const RadialExpansionExact exact;
    const RadialGrid grid{.nx = nr / 4, .nr = nr};
    const RadialGrid extended_grid{.nx = grid.nx, .nr = 2 * nr, .radius = 2.0};
    const auto built{build_radial_grid(grid)}, extended_built{build_radial_grid(extended_grid)};
    const CaloricallyPerfectGas gas{exact.gamma(), 1.0};
    AxisymmetricEulerSolver original{
        built.mesh, radial_initial_field(grid, exact), gas, radial_conditions(), {.cfl = 0.1}};
    AxisymmetricEulerSolver extended{
        extended_built.mesh, radial_initial_field(extended_grid, exact), gas, radial_conditions(), {.cfl = 0.1}};
    while (extended.diagnostics().time < 0.2)
    {
        const auto target{extended.step(0.2).time};
        original.step(target);
        if (original.diagnostics().time != target)
            throw std::runtime_error("Boundary control could not use the same accepted time step.");
    }
    RadialErrorAccumulator difference, original_error, extended_error;
    for (Index cell = 0; cell < grid.nx * grid.nr; ++cell)
    {
        const double volume{radial_cell_volume(grid, cell)};
        const auto bounds{radial_bounds(grid, cell)};
        const auto reference{exact.cell_average(bounds[0], bounds[1], 0.2)};
        difference.add(volume, original.state().state(cell), extended.state().state(cell));
        original_error.add(volume, original.state().state(cell), reference);
        extended_error.add(volume, extended.state().state(cell), reference);
    }
    return {difference.finish(), original_error.finish(), extended_error.finish(),
            extended.diagnostics().accepted_steps};
}

} // namespace cfd::verification
