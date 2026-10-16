#include "cfd/compressible/EulerSteadyMonitor.hpp"
#include "cfd/io/EulerSteadyCsvWriter.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "support/EulerTube.hpp"
#include "support/TestUtils.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <limits>
#include <numbers>

namespace
{
using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws_with_message;
using Type = cfd::EulerBoundaryConditionType;
using Status = cfd::EulerSteadyStatus;
const cfd::CaloricallyPerfectGas gas{1.4, 1.0};
constexpr double tolerance{3.0e-12};

cfd::EulerBoundaryConditions closed_conditions()
{
    return {4, {Type::Axis, Type::SlipWall, Type::SlipWall, Type::SlipWall}};
}

cfd::EulerSteadyOptions short_options()
{
    return {.check_interval = 1, .minimum_steps = 3, .maximum_steps = 10, .consecutive_checks = 2, .maximum_time = 1.0};
}

void test_uniform_steady_and_observer_does_not_change_evolution()
{
    const auto built{cfd::build_mesh(cfd::verification::make_euler_tube(3, 2))};
    const auto uniform{gas.to_conservative({1.0, 0.3, 0.0, 1.0})};
    const cfd::ConservativeField initial{built.mesh.cell_count(), uniform};
    cfd::AxisymmetricEulerSolver observed{
        built.mesh, initial, gas, cfd::verification::euler_tube_conditions(), {.collect_timings = true}};
    cfd::AxisymmetricEulerSolver control{built.mesh, initial, gas, cfd::verification::euler_tube_conditions()};
    cfd::EulerSteadyMonitor monitor{observed, short_options()};
    require(monitor.needs_observation(), "Initial observation missing.");
    auto observation{monitor.observe()};
    require(observation.status == Status::Running, "Initial observation alone cannot establish convergence.");
    for (const double residual : observation.normalized_residual)
        require(residual < tolerance, "Uniform residual is nonzero.");
    require_throws_with_message<std::invalid_argument>([&] { (void)monitor.observe(); }, "new accepted step",
                                                       "Duplicate checks counted as convergence.");
    while (observation.status == Status::Running)
    {
        observed.step(1.0);
        control.step(1.0);
        const auto before{observed.diagnostics()};
        observation = monitor.observe();
        require(observed.diagnostics().time == before.time &&
                    observed.diagnostics().accepted_steps == before.accepted_steps &&
                    observed.diagnostics().time_step == before.time_step &&
                    observed.diagnostics().integrals.mass == before.integrals.mass,
                "Observation changed accepted-step diagnostics.");
        for (cfd::Index cell = 0; cell < initial.size(); ++cell)
        {
            const auto a{observed.state().state(cell)}, b{control.state().state(cell)};
            require(a.density == b.density && a.momentum_x == b.momentum_x && a.momentum_r == b.momentum_r &&
                        a.total_energy_density == b.total_energy_density,
                    "Observation/cache changed the physical update.");
        }
    }
    require(observation.status == Status::Converged && observation.temporal.accepted_steps == 4,
            "Minimum steps and consecutive checks were not enforced.");
    require(observed.timings().assembly_seconds > 0.0 && observed.timings().boundary_flux_seconds > 0.0 &&
                observed.timings().boundary_flux_seconds <= observed.timings().assembly_seconds &&
                observed.timings().time_update_seconds > 0.0 && control.timings().assembly_seconds == 0.0,
            "Optional timing attribution is invalid or enabled by default.");
}

void test_local_residual_balanced_transient_and_tiny_dt()
{
    const auto built{cfd::build_mesh(cfd::verification::make_euler_tube(2, 1))};
    cfd::ConservativeField initial{2, gas.to_conservative({1.0, 0.0, 0.0, 1.0})};
    initial.set_state(1, gas.to_conservative({2.0, 0.0, 0.0, 1.0}));
    cfd::AxisymmetricEulerSolver solver{built.mesh, initial, gas, closed_conditions(), {.cfl = 1.0e-12}};
    cfd::EulerSteadyMonitor monitor{solver, short_options()};
    const auto first{monitor.observe()};
    // Independent two-cell Rusanov mass rate: F=-sqrt(1.4)/2,
    // area=pi; opposite cell rates cannot cancel in the L1 residual.
    require_near(first.current.absolute_residual[0], std::numbers::pi * std::sqrt(1.4), tolerance,
                 "Wrong integrated mass residual.");
    require_near(first.normalized_residual[0], std::sqrt(1.4), tolerance, "Wrong physical residual normalization.");
    require(first.mass_flow_relative_imbalance == 0.0 && first.normalized_residual[0] > 1.0,
            "A global balance hid local transience.");
    solver.step(1.0);
    const auto second{monitor.observe()};
    require(second.field_change[0] < 1.0e-10 && second.normalized_residual[0] > 1.0 && second.status == Status::Running,
            "A tiny time step falsely established stationarity.");
    const double left{solver.state().density()[0]}, right{solver.state().density()[1]};
    const double expected{std::numbers::pi * std::sqrt(1.4 / left) * (right - left)};
    require_near(second.current.absolute_residual[0], expected, 2.0e-15,
                 "Residual was stale instead of using the current state.");
    require(second.current.absolute_residual[0] < first.current.absolute_residual[0],
            "Current and old residuals were not distinguished.");
}

void test_normalization_under_physical_rescaling()
{
    std::array<double, 4> reference{};
    for (const double length : {1.0, 0.01, 100.0})
    {
        const double density_scale{length == 1.0 ? 1.0 : 7.0}, velocity_scale{length == 0.01 ? 3.0 : 0.4};
        auto raw{cfd::verification::make_euler_tube(2, 1)};
        for (auto &node : raw.nodes)
        {
            node.x *= length;
            node.y *= length;
        }
        const auto built{cfd::build_mesh(std::move(raw))};
        cfd::ConservativeField initial{2,
                                       gas.to_conservative({density_scale, 0.2 * velocity_scale, 0.1 * velocity_scale,
                                                            density_scale * velocity_scale * velocity_scale})};
        initial.set_state(1, gas.to_conservative({2.0 * density_scale, 0.1 * velocity_scale, -0.15 * velocity_scale,
                                                  0.7 * density_scale * velocity_scale * velocity_scale}));
        cfd::AxisymmetricEulerSolver solver{built.mesh, initial, gas, closed_conditions()};
        auto options{short_options()};
        options.reference_density = density_scale;
        options.reference_velocity = velocity_scale;
        options.reference_length = length;
        cfd::EulerSteadyMonitor monitor{solver, options};
        const auto observation{monitor.observe()};
        if (length == 1.0)
            reference = observation.normalized_residual;
        for (cfd::Index component = 0; component < reference.size(); ++component)
        {
            require(reference.at(component) > 1.0e-8, "Scaling test must exercise every equation.");
            require_near(observation.normalized_residual.at(component), reference.at(component), tolerance,
                         "Component residual normalization is not scale-consistent.");
        }
    }
}

void test_limits_options_and_csv()
{
    const auto built{cfd::build_mesh(cfd::verification::make_euler_tube(2, 1))};
    cfd::ConservativeField initial{2, gas.to_conservative({1.0, 0.0, 0.0, 1.0})};
    initial.set_state(1, gas.to_conservative({2.0, 0.0, 0.0, 1.0}));
    for (const bool step_limit : {false, true})
    {
        cfd::AxisymmetricEulerSolver solver{built.mesh, initial, gas, closed_conditions()};
        auto options{short_options()};
        options.minimum_steps = 1;
        options.maximum_steps = step_limit ? 1 : 100;
        options.maximum_time = step_limit ? 1.0 : 1.0e-6;
        cfd::EulerSteadyMonitor monitor{solver, options};
        (void)monitor.observe();
        solver.step(options.maximum_time);
        const auto observation{monitor.observe()};
        require(observation.status == (step_limit ? Status::StepLimit : Status::TimeLimit),
                "Non-converged limit returned the wrong status.");
    }
    cfd::AxisymmetricEulerSolver solver{built.mesh, initial, gas, closed_conditions()};
    for (const double invalid :
         {0.0, -1.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()})
    {
        auto options{short_options()};
        options.reference_velocity = invalid;
        require_throws_with_message<std::invalid_argument>(
            [&] { const cfd::EulerSteadyMonitor monitor{solver, options}; }, "positive", "Invalid reference accepted.");
        options = short_options();
        options.residual_tolerance = invalid;
        require_throws_with_message<std::invalid_argument>(
            [&] { const cfd::EulerSteadyMonitor monitor{solver, options}; }, "positive", "Invalid tolerance accepted.");
    }
    auto options{short_options()};
    options.check_interval = 0;
    require_throws_with_message<std::invalid_argument>([&] { const cfd::EulerSteadyMonitor monitor{solver, options}; },
                                                       "limits", "Invalid check interval accepted.");
    options = short_options();
    options.maximum_steps = 1;
    require_throws_with_message<std::invalid_argument>([&] { const cfd::EulerSteadyMonitor monitor{solver, options}; },
                                                       "limits", "Reversed step limits accepted.");
    cfd::EulerSteadyMonitor monitor{solver, short_options()};
    const auto observation{monitor.observe()};
    const auto path{std::filesystem::temp_directory_path() /
                    ("cfd-euler-steady-monitor-" +
                     std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".csv")};
    {
        cfd::EulerSteadyCsvWriter writer{path};
        writer.write(observation);
        const auto content{cfd::test::read_text_file(path)};
        cfd::test::require_contains(content, "residual_rhoux", "Missing Euler history header.");
        require(std::count(content.begin(), content.end(), '\n') == 2, "History row not flushed.");
    }
    std::filesystem::remove(path);
    require_throws_with_message<std::runtime_error>(
        [] { const cfd::EulerSteadyCsvWriter writer{std::filesystem::temp_directory_path()}; }, "CSV",
        "Output error not propagated.");
}
} // namespace

int main()
{
    int failures{};
    failures += cfd::test::run_test("uniform steady and non-mutating observer",
                                    test_uniform_steady_and_observer_does_not_change_evolution);
    failures += cfd::test::run_test("local transient residual and tiny time step",
                                    test_local_residual_balanced_transient_and_tiny_dt);
    failures += cfd::test::run_test("physical residual scaling", test_normalization_under_physical_rescaling);
    failures += cfd::test::run_test("steady limits invalid options and CSV", test_limits_options_and_csv);
    return cfd::test::finish_tests(failures, "Euler steady monitor");
}
