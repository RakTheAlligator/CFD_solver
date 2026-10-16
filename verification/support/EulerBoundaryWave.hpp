#pragma once

#include "EulerTube.hpp"
#include "cfd/compressible/AxisymmetricEulerSolver.hpp"
#include "cfd/mesh/MeshBuilder.hpp"

#include <algorithm>
#include <array>
#include <cmath>

namespace cfd::verification
{

/// Small isentropic right-running pulse in a straight tube. This is initial
/// data, not an asserted exact nonlinear solution at later times.
struct EulerBoundaryWaveParameters
{
    double gamma{1.4};
    double gas_constant{1.0};
    double density{1.0};
    double pressure{1.0};
    double mach{0.3};
    double amplitude{1.0e-4}; ///< Relative pressure perturbation.
    double radius{0.1};
};

inline MeshBuildResult build_boundary_wave_mesh(const Index nx, const Index nr, const double length,
                                                const double radius)
{
    auto raw{make_euler_tube(nx, nr)};
    for (auto &node : raw.nodes)
    {
        node.x *= length;
        node.y *= radius;
    }
    return build_mesh(std::move(raw));
}

inline ConservativeField boundary_wave_initial(const Mesh &mesh, const EulerBoundaryWaveParameters parameters)
{
    const CaloricallyPerfectGas gas{parameters.gamma, parameters.gas_constant};
    const double sound{std::sqrt(parameters.gamma * parameters.pressure / parameters.density)};
    ConservativeField field{mesh.cell_count()};
    for (Index cell = 0; cell < field.size(); ++cell)
    {
        const double coordinate{(mesh.cell_centers()[cell].x - 0.65) / 0.04};
        const double pressure{parameters.pressure * (1.0 + parameters.amplitude * std::exp(-coordinate * coordinate))};
        const double density{parameters.density * std::pow(pressure / parameters.pressure, 1.0 / parameters.gamma)};
        const double local_sound{std::sqrt(parameters.gamma * pressure / density)};
        const double velocity{parameters.mach * sound + 2.0 * (local_sound - sound) / (parameters.gamma - 1.0)};
        field.set_state(cell, gas.to_conservative({density, velocity, 0.0, pressure}));
    }
    return field;
}

inline EulerBoundaryConditions boundary_wave_conditions(const EulerBoundaryWaveParameters parameters,
                                                        const bool transmissive_outlet = false)
{
    const double temperature{parameters.pressure / (parameters.density * parameters.gas_constant)};
    const double ratio{1.0 + 0.5 * (parameters.gamma - 1.0) * parameters.mach * parameters.mach};
    const auto inlet{EulerBoundaryCondition::total_conditions_inlet(
        parameters.pressure * std::pow(ratio, parameters.gamma / (parameters.gamma - 1.0)), temperature * ratio,
        {1.0, 0.0})};
    const auto outlet{transmissive_outlet ? EulerBoundaryCondition{EulerBoundaryConditionType::Transmissive}
                                          : EulerBoundaryCondition::pressure_outlet(parameters.pressure)};
    return {4, {EulerBoundaryConditionType::Axis, EulerBoundaryConditionType::SlipWall, inlet, outlet}};
}

struct EulerBoundaryWaveBalance
{
    EulerIntegralQuantities initial;
    std::array<long double, 4> integrated_boundary{};
    long double integrated_source{};
    double maximum_step_error{};

    void add(const AxisymmetricEulerResult &result)
    {
        const long double dt{result.time_step};
        integrated_boundary[0] += dt * result.boundary_outflow_rate.mass;
        integrated_boundary[1] += dt * result.boundary_outflow_rate.momentum_x;
        integrated_boundary[2] += dt * result.boundary_outflow_rate.momentum_r;
        integrated_boundary[3] += dt * result.boundary_outflow_rate.energy;
        integrated_source += dt * result.radial_pressure_source_rate;
        maximum_step_error = std::max(
            {maximum_step_error, std::abs(result.balance_error.mass), std::abs(result.balance_error.momentum_x),
             std::abs(result.balance_error.momentum_r), std::abs(result.balance_error.energy)});
    }

    [[nodiscard]] std::array<double, 4> defect(const EulerIntegralQuantities final) const
    {
        return {static_cast<double>(final.mass - static_cast<long double>(initial.mass) + integrated_boundary[0]),
                static_cast<double>(final.momentum_x - static_cast<long double>(initial.momentum_x) +
                                    integrated_boundary[1]),
                static_cast<double>(final.momentum_r - static_cast<long double>(initial.momentum_r) +
                                    integrated_boundary[2] - integrated_source),
                static_cast<double>(final.energy - static_cast<long double>(initial.energy) + integrated_boundary[3])};
    }
};

struct EulerBoundaryWaveComparison
{
    double incident_pressure{};            ///< Initial max |delta p| in the common tube.
    double outgoing_pressure{};            ///< Max |p-p_base| on the final extended domain.
    double pressure_rms_difference{};      ///< Short minus extended on the whole short domain.
    double reflected_pressure{};           ///< Max |delta p_short - delta p_extended|.
    double minimum_pressure_difference{};  ///< Sign of the reflected pressure perturbation.
    double reflected_incoming_amplitude{}; ///< max |delta u-delta p/(rho0*a0)|, short minus extended.
    double reflection_ratio{}; ///< Reflected pressure / initial incident pressure, including numerical diffusion.
    double maximum_balance_defect{};
    double maximum_step_balance_defect{};
    AxisymmetricEulerResult result;
};

/// Compare identical dx/dr/initial data at the same physical comparison times;
/// the reference may subcycle to obey its own CFL. Continuum characteristics
/// from the extended boundary at x=2 cannot reach x<=1 during t=0.85/a0;
/// this is not a claim of exactly compact numerical diffusion tails.
/// The initial pulse centre=0.65 and base M=0.3 place the reflected pulse
/// inside the short domain at that time and before it reaches the inlet.
inline EulerBoundaryWaveComparison compare_boundary_wave(const Index nx, const bool transmissive_outlet = false)
{
    constexpr EulerBoundaryWaveParameters parameters{};
    const CaloricallyPerfectGas gas{parameters.gamma, parameters.gas_constant};
    const double sound{std::sqrt(parameters.gamma)};
    const auto short_mesh{build_boundary_wave_mesh(nx, 3, 1.0, parameters.radius)};
    const auto extended_mesh{build_boundary_wave_mesh(2 * nx, 3, 2.0, parameters.radius)};
    const auto initial{boundary_wave_initial(short_mesh.mesh, parameters)};
    AxisymmetricEulerSolver short_solver{short_mesh.mesh, initial, gas,
                                         boundary_wave_conditions(parameters, transmissive_outlet)};
    AxisymmetricEulerSolver extended_solver{extended_mesh.mesh, boundary_wave_initial(extended_mesh.mesh, parameters),
                                            gas, boundary_wave_conditions(parameters)};
    EulerBoundaryWaveBalance balance{.initial = short_solver.diagnostics().integrals};
    EulerBoundaryWaveComparison comparison;
    for (Index cell = 0; cell < initial.size(); ++cell)
        comparison.incident_pressure =
            std::max(comparison.incident_pressure,
                     std::abs(gas.evaluate(initial.state(cell)).primitive.pressure - parameters.pressure));
    const double final_time{0.85 / sound};
    while (short_solver.diagnostics().time < final_time)
    {
        const auto &step{short_solver.step(final_time)};
        balance.add(step);
        extended_solver.advance_to(step.time);
    }
    const AxisymmetricGeometry geometry{short_mesh.mesh};
    long double squared{}, total_volume{};
    for (Index cell = 0; cell < short_mesh.mesh.cell_count(); ++cell)
    {
        const Index i{cell % nx}, j{cell / nx};
        const auto small{gas.evaluate(short_solver.state().state(cell)).primitive};
        const auto large{gas.evaluate(extended_solver.state().state(j * 2 * nx + i)).primitive};
        const double pressure_difference{small.pressure - large.pressure};
        squared += geometry.cell_volumes()[cell] * static_cast<long double>(pressure_difference) * pressure_difference;
        total_volume += geometry.cell_volumes()[cell];
        comparison.reflected_pressure = std::max(comparison.reflected_pressure, std::abs(pressure_difference));
        comparison.minimum_pressure_difference = std::min(comparison.minimum_pressure_difference, pressure_difference);
        comparison.reflected_incoming_amplitude = std::max(
            comparison.reflected_incoming_amplitude,
            std::abs(small.velocity_x - large.velocity_x - pressure_difference / (parameters.density * sound)));
    }
    for (Index cell = 0; cell < extended_solver.state().size(); ++cell)
        comparison.outgoing_pressure = std::max(
            comparison.outgoing_pressure,
            std::abs(gas.evaluate(extended_solver.state().state(cell)).primitive.pressure - parameters.pressure));
    comparison.pressure_rms_difference = static_cast<double>(std::sqrt(squared / total_volume));
    comparison.reflection_ratio = comparison.reflected_pressure / comparison.incident_pressure;
    for (const double component : balance.defect(short_solver.diagnostics().integrals))
        comparison.maximum_balance_defect = std::max(comparison.maximum_balance_defect, std::abs(component));
    comparison.maximum_step_balance_defect = balance.maximum_step_error;
    comparison.result = short_solver.diagnostics();
    return comparison;
}

} // namespace cfd::verification
