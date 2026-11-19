#include "cfd/compressible/AxisymmetricEulerSteadySolver.hpp"
#include "cfd/compressible/EulerCharacteristicBoundaryFlux.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "support/EulerNozzle.hpp"
#include "support/EulerTube.hpp"
#include "support/TestUtils.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <thread>
#include <type_traits>

namespace
{
using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws_with_message;
using Type = cfd::EulerBoundaryConditionType;
using Status = cfd::EulerLocalSteadyStatus;
const cfd::CaloricallyPerfectGas gas{1.4, 1.0};
constexpr double tolerance{2.0e-12};
static_assert(!std::is_copy_constructible_v<cfd::AxisymmetricEulerSteadySolver>);
template <class T>
concept PhysicalStep = requires(T &solver) { solver.step(1.0); };
template <class T>
concept PhysicalTime = requires(T result) { result.time; };
static_assert(PhysicalStep<cfd::AxisymmetricEulerSolver>);
static_assert(!PhysicalStep<cfd::AxisymmetricEulerSteadySolver>);
static_assert(!PhysicalTime<cfd::EulerLocalSteadyResult>);

cfd::EulerLocalSteadyOptions options()
{
    cfd::EulerLocalSteadyOptions o;
    o.maximum_iterations = 20;
    o.minimum_iterations = 1;
    o.check_interval = 1;
    o.consecutive_checks = 2;
    return o;
}

void same_field(const cfd::ConservativeField &a, const cfd::ConservativeField &b)
{
    require(a.size() == b.size(), "Cardinality changed.");
    for (cfd::Index cell = 0; cell < a.size(); ++cell)
    {
        const auto u{a.state(cell)}, v{b.state(cell)};
        require_near(u.density, v.density, 0, "Density changed.");
        require_near(u.momentum_x, v.momentum_x, 0, "Axial momentum changed.");
        require_near(u.momentum_r, v.momentum_r, 0, "Radial momentum changed.");
        require_near(u.total_energy_density, v.total_energy_density, 0, "Energy changed.");
    }
}

void test_uniform_stationary_tri_quad_axis_source()
{
    for (const auto type : {cfd::verification::EulerTubeCells::Quad, cfd::verification::EulerTubeCells::Triangle,
                            cfd::verification::EulerTubeCells::Mixed})
    {
        const auto built{cfd::build_mesh(cfd::verification::make_euler_tube(3, 2, type))};
        for (const double ux : {0.0, .3})
        {
            const cfd::ConservativeField initial{built.mesh.cell_count(), gas.to_conservative({1, ux, 0, 1})};
            cfd::AxisymmetricEulerSteadySolver solver{built.mesh, initial, gas,
                                                      cfd::verification::euler_tube_conditions(), options()};
            require(solver.observe().status == Status::Running, "Initial observation cannot fabricate convergence.");
            solver.observe();
            require(solver.iterate().status == Status::Running, "Need two distinct stable observations.");
            require(solver.iterate().status == Status::Converged, "Uniform compatible state must converge.");
            const auto &d{solver.diagnostics()};
            require(d.iteration == 2 && d.observation_iteration == 2 && d.rejected_attempts == 0,
                    "Incorrect steady counters.");
            for (const double r : d.maximum_cell_residual)
                require(r < 1e-10, "Uniform state has a nonzero discrete residual.");
            for (const double b : d.stationary_balance)
                require(b < 1e-10, "Uniform boundary/source balance.");
            require_near(d.radial_pressure_source_rate, 2 * std::numbers::pi, tolerance,
                         "Axisymmetric pressure source.");
            require_near(d.boundary_outflow_rate.momentum_r, 2 * std::numbers::pi, tolerance, "Radial wall flux.");
            for (cfd::Index cell = 0; cell < initial.size(); ++cell)
            {
                const auto u{solver.state().state(cell)}, exact{initial.state(cell)};
                require_near(u.density, exact.density, tolerance, "Uniform density.");
                require_near(u.momentum_x, exact.momentum_x, tolerance, "Uniform axial momentum.");
                require_near(u.momentum_r, 0, tolerance, "Spurious radial flow.");
                require_near(u.total_energy_density, exact.total_energy_density, tolerance, "Uniform energy.");
            }
            require(solver.iterate().iteration == 2, "Terminal solver must not advance.");
        }
    }
}

void test_local_cfl_jacobi_and_internal_conservation()
{
    auto raw{cfd::verification::make_euler_tube(2, 1)};
    for (auto &node : raw.nodes)
        if (node.x == .5)
            node.x = .25;
    const auto built{cfd::build_mesh(std::move(raw))};
    cfd::ConservativeField field{2, gas.to_conservative({1, 0, 0, 1})};
    field.set_state(1, gas.to_conservative({2, 0, 0, 1}));
    const cfd::EulerBoundaryConditions walls{4, {Type::Axis, Type::SlipWall, Type::SlipWall, Type::SlipWall}};
    const std::array<cfd::Index, 2> regions{11, 42};
    auto o{options()};
    cfd::AxisymmetricEulerSteadySolver solver{built.mesh, field, gas, walls, o, regions};
    cfd::AxisymmetricEulerSolver transient{built.mesh, field, gas, walls, {.cfl = o.cfl}};
    const auto wave{transient.evaluate_cfl()};
    const double a{std::sqrt(1.4)};
    const double flux{-.5 * a}; // Independent Rusanov: same p and zero velocity, rho_right-rho_left=1.
    const double rate{-std::numbers::pi * flux};
    const auto initial{solver.observe()};
    require_near(initial.absolute_residual[0], 2 * rate, tolerance, "Integrated mass L1 residual.");
    require_near(wave.cell_wave_sums[0], 2.5 * std::numbers::pi * a, tolerance, "Owner spectral sum.");
    require_near(wave.cell_wave_sums[1], std::numbers::pi * (a + 2.5 * std::sqrt(.7)), tolerance,
                 "Neighbor spectral sum must include both states on the internal face.");
    solver.iterate();
    std::array<double, 2> dt{};
    std::array<double, 2> change{};
    for (cfd::Index cell = 0; cell < 2; ++cell)
    {
        const double v{solver.geometry().cell_volumes()[cell]};
        dt.at(cell) = static_cast<double>(static_cast<long double>(o.cfl) * v / wave.cell_wave_sums[cell]);
        require(solver.local_steps()[cell] == dt.at(cell), "Local step must use the actual assembled spectral sum.");
        const double expected{(cell == 0 ? rate : -rate) * dt.at(cell) / v};
        change.at(cell) = solver.state().density()[cell] - field.density()[cell];
        require_near(change.at(cell), expected, tolerance, "Both Jacobi updates must use the old state.");
        require_near(solver.state().momentum_r()[cell], 0, tolerance, "Pressure source must cancel radial wall flux.");
    }
    require(dt[0] != dt[1], "Fixture must exercise different pseudo-steps.");
    const auto v{solver.geometry().cell_volumes()};
    require_near(v[0] * change[0] / dt[0] + v[1] * change[1] / dt[1], 0, tolerance,
                 "Unweighted internal residual contributions must cancel.");
    require(std::abs(v[0] * change[0] + v[1] * change[1]) > 1e-5,
            "Local pseudo-time is not a globally conservative physical step.");
    require(solver.regional_residuals().size() == 2 && solver.regional_residuals()[0].region == 11 &&
                solver.regional_residuals()[1].region == 42,
            "Region labels must be retained.");
    require_near(initial.normalized_residual[0], initial.absolute_residual[0] / std::numbers::pi, tolerance,
                 "Fixed physical normalization.");
}

void test_rejection_is_global_and_transactional()
{
    const auto built{cfd::build_mesh(cfd::verification::make_euler_tube(1, 1))};
    const cfd::ConservativeField field{1, gas.to_conservative({1, 0, 10, .001})};
    const cfd::EulerBoundaryConditions bc{4, {Type::Axis, Type::Transmissive, Type::Transmissive, Type::Transmissive}};
    auto o{options()};
    o.cfl = 1;
    o.maximum_step_attempts = 1;
    cfd::AxisymmetricEulerSteadySolver rejected{built.mesh, field, gas, bc, o};
    const auto d{rejected.iterate()};
    require(d.status == Status::NumericalFailure && d.iteration == 0 && d.rejected_attempts == 1,
            "Exhausted candidate attempts must be explicit.");
    same_field(rejected.state(), field);
    cfd::test::require_contains(rejected.failure_message(), "cell=0", "Rejection must identify the failing cell.");
    require(rejected.iterate().rejected_attempts == 1, "Terminal failure may not retry silently.");
    o.maximum_step_attempts = 20;
    cfd::AxisymmetricEulerSteadySolver retried{built.mesh, field, gas, bc, o};
    const auto r{retried.iterate()};
    require(r.iteration == 1 && r.rejected_attempts == 1 && r.step_scale == .5, "Expected global halving.");
    require_near(r.minimum_local_step, .5 / (20 + 4 * std::sqrt(.0014)), tolerance, "Incorrect reduced pseudo-step.");
    require(r.minimum_density > 0 && r.minimum_pressure > 0, "Accepted state must remain admissible.");
    cfd::AxisymmetricEulerSolver physical{built.mesh, field, gas, bc, {.cfl = 1}};
    physical.step(1.0);
    same_field(retried.state(), physical.state()); // One cell: the two update formulas coincide exactly.
    o.maximum_rejected_attempts = 1;
    cfd::AxisymmetricEulerSteadySolver capped{built.mesh, field, gas, bc, o};
    require(capped.iterate().status == Status::NumericalFailure, "Rejection budget must stop the solve.");
    same_field(capped.state(), field);

    auto raw{cfd::verification::make_euler_tube(2, 1)};
    for (auto &node : raw.nodes)
        if (node.x == .5)
            node.x = .001; // First candidate is admissible; the wider second cell fails.
    const auto pair{cfd::build_mesh(std::move(raw))};
    const cfd::ConservativeField two{2, field.state(0)};
    o.maximum_step_attempts = 1;
    cfd::AxisymmetricEulerSteadySolver partial{pair.mesh, two, gas, bc, o};
    require(partial.iterate().status == Status::NumericalFailure, "Expected multicell candidate rejection.");
    cfd::test::require_contains(partial.failure_message(), "cell=1", "Must fail after writing an earlier candidate.");
    same_field(partial.state(), two);
}

cfd::ConservativeState analytical_boundary_trial(const double alpha)
{
    // One unit cylinder, gamma=2: initial rho=1, ux=1/2, p=2.
    // TCI P0=2,T0=25/16 gives 3*q^2+7*q-1/4=0; the right wall
    // Rankine-Hugoniot relation gives p_wall=2+(3+sqrt(265))/16.
    const double q{(2 * std::sqrt(13.0) - 7) / 6}, temperature{25.0 / 16 - q * q / 4};
    const double pressure{2 * std::pow(temperature / (25.0 / 16), 2)}, rho{pressure / temperature};
    const double wall_pressure{2 + (3 + std::sqrt(265.0)) / 16};
    // Left/right areas=pi, outer area=2*pi, volume=pi. Wave bounds
    // are 5/2, sqrt(2*p_wall), and 2, respectively; axis area is zero.
    const double dt{alpha / (2.5 + std::sqrt(2 * wall_pressure) + 4)};
    return {1 + dt * rho * q, .5 + dt * (rho * q * q + pressure - wall_pressure), 0,
            2.125 + dt * (2 * pressure + .5 * rho * q * q) * q};
}

void test_candidate_boundary_rejection_and_global_retry()
{
    const auto built{cfd::build_mesh(cfd::verification::make_euler_tube(1, 1))};
    const cfd::CaloricallyPerfectGas model{2, 1};
    const cfd::ConservativeField initial{1, model.to_conservative({1, .5, 0, 2})};
    const auto inlet{cfd::EulerBoundaryCondition::total_conditions_inlet(2, 25.0 / 16, {1, 0})};
    const cfd::EulerBoundaryConditions bc{4, {Type::Axis, Type::SlipWall, inlet, Type::SlipWall}};
    for (const double alpha : {1.0, .5, .25, .125})
    {
        const auto value{model.evaluate(analytical_boundary_trial(alpha))};
        const double invariant{-value.primitive.velocity_x + 2 * value.thermodynamics.sound_speed};
        require(value.primitive.pressure > 0 && value.primitive.density > 0,
                "Boundary-rejection fixture must remain thermodynamically admissible.");
        require((invariant > 5 / std::sqrt(2.0)) == (alpha >= .25),
                "Independent analytical trial does not cross the intended stagnation bound.");
    }
    auto o{options()};
    o.cfl = 1;
    o.maximum_step_attempts = 3;
    cfd::AxisymmetricEulerSteadySolver exhausted{built.mesh, initial, model, bc, o};
    const auto before{exhausted.observe()};
    const auto failure{exhausted.iterate()};
    require(failure.status == Status::NumericalFailure && failure.iteration == 0 && failure.rejected_attempts == 3,
            "Every incompatible candidate must be counted and refused before commitment.");
    same_field(exhausted.state(), initial);
    require(failure.integrals == before.integrals && failure.absolute_residual == before.absolute_residual &&
                failure.maximum_cell_residual == before.maximum_cell_residual &&
                failure.boundary_outflow_rate == before.boundary_outflow_rate &&
                failure.observation_iteration == before.observation_iteration,
            "Rejected BC trials changed the accepted-state diagnostics.");
    cfd::test::require_contains(exhausted.failure_message(), "face=", "Missing rejected face context.");
    cfd::test::require_contains(exhausted.failure_message(), "cell=0", "Missing rejected owner context.");
    cfd::test::require_contains(exhausted.failure_message(), "stagnation bound", "Wrong incompatibility cause.");
    require(exhausted.iterate().rejected_attempts == 3, "An exhausted solve must remain terminal.");

    o.maximum_step_attempts = 4;
    cfd::AxisymmetricEulerSteadySolver retry{built.mesh, initial, model, bc, o};
    const auto accepted{retry.iterate()};
    require(accepted.status == Status::Running && accepted.iteration == 1 && accepted.rejected_attempts == 3 &&
                accepted.step_scale == .125,
            "Compatible fourth candidate must be accepted after three global halvings.");
    const auto expected{analytical_boundary_trial(.125)}, actual{retry.state().state(0)};
    require_near(actual.density, expected.density, tolerance, "Retried density differs from analytical Jacobi trial.");
    require_near(actual.momentum_x, expected.momentum_x, tolerance, "Retried momentum differs from analytical trial.");
    require_near(actual.momentum_r, 0, 0, "Retried trial broke radial equilibrium.");
    require_near(actual.total_energy_density, expected.total_energy_density, tolerance, "Retried energy mismatch.");
    (void)cfd::compute_euler_characteristic_boundary_flux(actual, {-1, 0}, inlet, model);
    require(retry.failure_message().empty(), "Recovered rejection retained a terminal error.");

    o.maximum_rejected_attempts = 2;
    cfd::AxisymmetricEulerSteadySolver capped{built.mesh, initial, model, bc, o};
    require(capped.iterate().status == Status::NumericalFailure && capped.diagnostics().rejected_attempts == 2,
            "BC retries must respect the existing cumulative rejection cap.");
    same_field(capped.state(), initial);
}

void test_initial_boundary_incompatibility_and_fatal_errors()
{
    const auto built{cfd::build_mesh(cfd::verification::make_euler_tube(1, 1))};
    const cfd::CaloricallyPerfectGas model{2, 1};
    const cfd::ConservativeField initial{1, model.to_conservative({1, 0, 0, 2})};
    const auto incompatible{cfd::EulerBoundaryCondition::total_conditions_inlet(2, .5, {1, 0})};
    const auto wrong_direction{cfd::EulerBoundaryCondition::total_conditions_inlet(2, 2, {-1, 0})};
    for (const auto &condition : {incompatible, wrong_direction})
    {
        cfd::AxisymmetricEulerSteadySolver solver{
            built.mesh, initial, model, {4, {Type::Axis, Type::SlipWall, condition, Type::SlipWall}}, options()};
        const auto result{solver.iterate()};
        require(result.status == Status::NumericalFailure && result.iteration == 0 && result.rejected_attempts == 0,
                "An initial-state or configuration failure must not trigger candidate retries.");
        same_field(solver.state(), initial);
        cfd::test::require_contains(solver.failure_message(),
                                    condition.direction().x > 0 ? "stagnation bound" : "direction",
                                    "Initial BC failure lost its cause.");
    }
    const cfd::ConservativeField rest{1, model.to_conservative({1, 0, 0, .5})};
    const auto overflow{
        cfd::EulerBoundaryCondition::total_conditions_inlet(std::numeric_limits<double>::max(), .5, {1, 0})};
    cfd::AxisymmetricEulerSteadySolver fatal{
        built.mesh, rest, model, {4, {Type::Axis, Type::SlipWall, overflow, Type::SlipWall}}, options()};
    require(fatal.iterate().status == Status::NumericalFailure && fatal.diagnostics().rejected_attempts == 0,
            "Unrepresentable BC arithmetic must remain fatal without candidate retries.");
    cfd::test::require_contains(fatal.failure_message(), "non-representable", "Arithmetic error was concealed.");
    same_field(fatal.state(), rest);
}

void test_small_steps_cannot_fake_convergence_and_stop_statuses()
{
    const auto built{cfd::build_mesh(cfd::verification::make_euler_tube(2, 1))};
    cfd::ConservativeField field{2, gas.to_conservative({1, 0, 0, 1})};
    field.set_state(1, gas.to_conservative({2, 0, 0, 1}));
    auto o{options()};
    o.cfl = 1e-200;
    o.maximum_iterations = 2;
    cfd::AxisymmetricEulerSteadySolver limited{built.mesh, field, gas, cfd::verification::euler_tube_conditions(), o};
    require(limited.solve().status == Status::MaximumIterations,
            "Zero field change with finite residual cannot converge.");
    require(limited.diagnostics().field_change[0] == 0 && limited.diagnostics().maximum_cell_residual[0] > .1,
            "Fixture must have negligible updates but a substantial actual residual.");
    o.maximum_iterations = 20;
    o.stagnation_checks = 2;
    cfd::AxisymmetricEulerSteadySolver stagnant{built.mesh, field, gas, cfd::verification::euler_tube_conditions(), o};
    require(stagnant.solve().status == Status::Stagnated, "Residual plateau is not convergence.");
    require(stagnant.diagnostics().iteration == (1 + o.stagnation_windows) * o.stagnation_checks * o.check_interval,
            "A true plateau must stop after the reference window and the configured unproductive windows.");
    o.maximum_wall_seconds = .001;
    cfd::AxisymmetricEulerSteadySolver wall{built.mesh, field, gas, cfd::verification::euler_tube_conditions(), o};
    wall.observe();
    std::this_thread::sleep_for(std::chrono::milliseconds{2});
    require(wall.solve().status == Status::WallLimit, "Real wall budget must be enforced.");
    same_field(wall.state(), field);

    const auto single{cfd::build_mesh(cfd::verification::make_euler_tube(1, 1))};
    const cfd::ConservativeField inward{1, gas.to_conservative({1, 0, -10, .001})};
    const cfd::EulerBoundaryConditions bc{4, {Type::Axis, Type::Transmissive, Type::Transmissive, Type::Transmissive}};
    o = options();
    o.cfl = 1;
    o.maximum_residual_growth = 1.1;
    cfd::AxisymmetricEulerSteadySolver growing{single.mesh, inward, gas, bc, o};
    require(growing.solve().status == Status::Diverged, "Residual growth must have a distinct stop status.");
}

void test_validation_and_unrecoverable_assembly()
{
    const auto built{cfd::build_mesh(cfd::verification::make_euler_tube(1, 1))};
    const cfd::ConservativeField field{1, gas.to_conservative({1, 0, 0, 1})};
    for (const double bad :
         {0.0, -1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
    {
        for (int which = 0; which < 6; ++which)
        {
            auto o{options()};
            switch (which)
            {
            case 0:
                o.cfl = bad;
                break;
            case 1:
                o.reference_density = bad;
                break;
            case 2:
                o.maximum_wall_seconds = bad;
                break;
            case 3:
                o.residual_tolerance = bad;
                break;
            case 4:
                o.reference_length = bad;
                break;
            default:
                o.reference_velocity = bad;
                break;
            }
            cfd::test::require_throws<std::invalid_argument>(
                [&] {
                    const cfd::AxisymmetricEulerSteadySolver invalid{built.mesh, field, gas,
                                                                     cfd::verification::euler_tube_conditions(), o};
                },
                "Invalid steady scalar option accepted.");
        }
    }
    for (int which = 0; which < 8; ++which)
    {
        auto o{options()};
        switch (which)
        {
        case 0:
            o.maximum_iterations = 0;
            break;
        case 1:
            o.check_interval = 0;
            break;
        case 2:
            o.consecutive_checks = 0;
            break;
        case 3:
            o.maximum_step_attempts = 0;
            break;
        case 4:
            o.maximum_rejected_attempts = 0;
            break;
        case 5:
            o.residual_reduction = 1.1;
            break;
        case 6:
            o.stagnation_windows = 0;
            break;
        default:
            o.minimum_residual_progress = 1;
            break;
        }
        cfd::test::require_throws<std::invalid_argument>(
            [&] {
                const cfd::AxisymmetricEulerSteadySolver invalid{built.mesh, field, gas,
                                                                 cfd::verification::euler_tube_conditions(), o};
            },
            "Invalid steady count/reduction accepted.");
    }
    const std::array<cfd::Index, 2> labels{0, 1};
    require_throws_with_message<std::invalid_argument>(
        [&] {
            const cfd::AxisymmetricEulerSteadySolver invalid{
                built.mesh, field, gas, cfd::verification::euler_tube_conditions(), options(), labels};
        },
        "region cardinality", "Region cardinality accepted.");
    const auto extreme{gas.to_conservative({1, 1e110, 0, 1e215})};
    const cfd::ConservativeField huge{1, extreme};
    cfd::AxisymmetricEulerSteadySolver failed{built.mesh, huge, gas, cfd::verification::euler_tube_conditions(),
                                              options()};
    require(failed.iterate().status == Status::NumericalFailure && failed.diagnostics().iteration == 0,
            "Assembly failure must not be retried or hidden.");
    same_field(failed.state(), huge);
    auto tiny{options()};
    tiny.cfl = std::numeric_limits<double>::denorm_min();
    cfd::AxisymmetricEulerSteadySolver underflow{built.mesh, field, gas, cfd::verification::euler_tube_conditions(),
                                                 tiny};
    require(underflow.iterate().status == Status::NumericalFailure, "Unrepresentable local step must be rejected.");
    cfd::test::require_contains(underflow.failure_message(), "underflows", "Local underflow diagnostic.");
    same_field(underflow.state(), field);
}

void test_nozzle_global_and_local_use_the_same_residual()
{
    const cfd::verification::EulerNozzleReference reference{};
    const auto built{cfd::build_mesh(cfd::verification::make_euler_nozzle(32, 4))};
    const auto initial{cfd::verification::nozzle_initial_field(32, 4, reference, false)};
    auto o{options()};
    o.reference_velocity = std::sqrt(1.4);
    cfd::AxisymmetricEulerSteadySolver local{built.mesh, initial, gas, reference.conditions(), o};
    cfd::AxisymmetricEulerSolver global{built.mesh, initial, gas, reference.conditions(), {.cfl = o.cfl}};
    const auto local_residual{local.observe()};
    const auto global_residual{global.evaluate_steady_residual()};
    require(local_residual.absolute_residual == global_residual.absolute_residual,
            "The same nozzle state must produce identical four-equation residuals in both modes.");
    const auto wave{global.evaluate_cfl()};
    local.iterate();
    for (cfd::Index cell = 0; cell < initial.size(); ++cell)
    {
        const double expected{static_cast<double>(static_cast<long double>(o.cfl) *
                                                  local.geometry().cell_volumes()[cell] / wave.cell_wave_sums[cell])};
        require(local.local_steps()[cell] == expected, "Local CFL must use the same assembled face/BC wave sums.");
    }
    const double global_step{wave.time_step};
    global.step(1.0);
    require(global.diagnostics().time_step == global_step && global.diagnostics().rejected_steps == 0 &&
                local.diagnostics().rejected_attempts == 0,
            "Same-state comparison must use the proposed steps without clipping or rejection.");
    require(local_residual.radial_pressure_source_rate == global.diagnostics().radial_pressure_source_rate &&
                local_residual.boundary_outflow_rate == global.diagnostics().boundary_outflow_rate,
            "Nozzle boundary fluxes and radial source must be identical.");
    for (cfd::Index cell = 0; cell < initial.size(); ++cell)
    {
        const auto old{initial.state(cell)}, g{global.state().state(cell)}, l{local.state().state(cell)};
        const double ratio{local.local_steps()[cell] / global_step};
        require_near(l.density, old.density + ratio * (g.density - old.density), 2e-13, "Local/global mass update.");
        require_near(l.momentum_x, old.momentum_x + ratio * (g.momentum_x - old.momentum_x), 2e-13,
                     "Local/global axial update.");
        require_near(l.momentum_r, old.momentum_r + ratio * (g.momentum_r - old.momentum_r), 2e-13,
                     "Local/global radial update.");
        require_near(l.total_energy_density,
                     old.total_energy_density + ratio * (g.total_energy_density - old.total_energy_density), 2e-13,
                     "Local/global energy update.");
    }
}

void test_farfield_volume_cannot_hide_engine_residual()
{
    auto raw{cfd::verification::make_euler_tube(2, 1)};
    for (auto &node : raw.nodes)
        node.x = node.x == .5 ? .01 : node.x * 1e6;
    const auto built{cfd::build_mesh(std::move(raw))};
    cfd::ConservativeField initial{2, gas.to_conservative({1, 0, 0, 1})};
    initial.set_state(1, gas.to_conservative({2, 0, 0, 1}));
    auto o{options()};
    o.cfl = 1e-200;
    o.maximum_iterations = 2;
    o.residual_tolerance = 1e-4;
    const std::array<cfd::Index, 2> labels{0, 1};
    for (const bool regional : {false, true})
    {
        cfd::AxisymmetricEulerSteadySolver solver{
            built.mesh, initial,
            gas,        cfd::verification::euler_tube_conditions(),
            o,          regional ? std::span<const cfd::Index>{labels} : std::span<const cfd::Index>{}};
        const auto &d{solver.solve()};
        require(d.normalized_residual[0] < o.residual_tolerance && d.maximum_cell_residual[0] > o.residual_tolerance,
                "Fixture must expose global-volume dilution of the engine residual.");
        require(d.status == Status::MaximumIterations, "Large farfield volume must not create false convergence.");
        if (regional)
            require(solver.regional_residuals()[0].normalized_residual[0] > o.residual_tolerance,
                    "Engine regional residual must remain visible.");
    }
}

void test_convergence_from_a_nonstationary_supersonic_tube()
{
    const auto built{cfd::build_mesh(cfd::verification::make_euler_tube(8, 1))};
    const cfd::PrimitiveState exact{1, 2, 0, 1};
    cfd::ConservativeField field{built.mesh.cell_count()};
    for (cfd::Index cell = 0; cell < field.size(); ++cell)
        field.set_state(cell, gas.to_conservative(
                                  {1 + .05 * std::sin(std::numbers::pi * built.mesh.cell_centers()[cell].x), 2, 0, 1}));
    const cfd::EulerBoundaryConditions bc{
        4, {Type::Axis, Type::SlipWall, cfd::EulerBoundaryCondition::supersonic_inlet(exact), Type::SupersonicOutlet}};
    auto o{options()};
    o.maximum_iterations = 8000;
    o.check_interval = 20;
    o.stagnation_checks = 50;
    cfd::AxisymmetricEulerSteadySolver solver{built.mesh, field, gas, bc, o};
    const auto &d{solver.solve()};
    require(d.status == Status::Converged, std::string{"Perturbed tube did not converge: "} + cfd::to_string(d.status));
    require(d.iteration > 20 && d.rejected_attempts == 0, "Tube must genuinely reduce a nonzero residual.");
    for (cfd::Index cell = 0; cell < field.size(); ++cell)
    {
        const auto p{gas.evaluate(solver.state().state(cell)).primitive};
        require_near(p.density, exact.density, 1e-7, "Stationary tube density.");
        require_near(p.velocity_x, exact.velocity_x, 1e-7, "Stationary tube velocity.");
        require_near(p.pressure, exact.pressure, 1e-7, "Stationary tube pressure.");
    }
}

void test_damped_nozzle_oscillations_are_not_a_plateau()
{
    const cfd::verification::EulerNozzleReference reference{};
    const auto built{cfd::build_mesh(cfd::verification::make_euler_nozzle(16, 2))};
    const auto initial{cfd::verification::nozzle_initial_field(16, 2, reference, false)};
    auto o{options()};
    o.maximum_iterations = 3000;
    o.check_interval = 20;
    o.reference_velocity = std::sqrt(1.4);
    cfd::AxisymmetricEulerSteadySolver solver{built.mesh, initial, gas, reference.conditions(), o};
    require(solver.solve().status == Status::MaximumIterations,
            "Damped nozzle oscillations must not be classified from an isolated trough or window endpoint.");
    require(solver.diagnostics().maximum_cell_residual[3] < .1, "Nozzle residual must actually decrease.");
}

void test_damped_refined_nozzle_oscillations_are_not_stagnation()
{
    const cfd::verification::EulerNozzleReference reference{};
    const auto built{cfd::build_mesh(cfd::verification::make_euler_nozzle(32, 4))};
    const auto initial{cfd::verification::nozzle_initial_field(32, 4, reference, false)};
    for (const double cfl : {.05, .025})
    {
        cfd::EulerLocalSteadyOptions o;
        o.cfl = cfl;
        o.maximum_iterations = cfl == .05 ? 14000 : 8000;
        o.reference_velocity = std::sqrt(1.4);
        cfd::AxisymmetricEulerSteadySolver solver{built.mesh, initial, gas, reference.conditions(), o};
        solver.observe();
        const double initial_residual{*std::ranges::max_element(solver.diagnostics().maximum_cell_residual)};
        require(solver.solve().status == Status::MaximumIterations,
                "Refined nozzle must survive isolated rising residual envelopes while the oscillations decay.");
        require(*std::ranges::max_element(solver.diagnostics().maximum_cell_residual) < initial_residual / 20,
                "Refined nozzle must genuinely reduce the actual stationary residual.");
        require(solver.diagnostics().rejected_attempts == 0, "Fixture must not rely on positivity retries.");
    }
}
} // namespace

int main()
{
    int failures{};
    failures += cfd::test::run_test("Local steady uniform TRI/QUAD and radial source",
                                    test_uniform_stationary_tri_quad_axis_source);
    failures += cfd::test::run_test("Local CFL, Jacobi and conservative internal residual",
                                    test_local_cfl_jacobi_and_internal_conservation);
    failures +=
        cfd::test::run_test("Local steady global transactional rejection", test_rejection_is_global_and_transactional);
    failures += cfd::test::run_test("Local steady candidate BC rejection and retry",
                                    test_candidate_boundary_rejection_and_global_retry);
    failures += cfd::test::run_test("Local steady initial BC incompatibility and fatal errors",
                                    test_initial_boundary_incompatibility_and_fatal_errors);
    failures += cfd::test::run_test("Local steady residual and stop statuses",
                                    test_small_steps_cannot_fake_convergence_and_stop_statuses);
    failures +=
        cfd::test::run_test("Local steady validation and assembly failure", test_validation_and_unrecoverable_assembly);
    failures += cfd::test::run_test("Local steady protects Engine from farfield dilution",
                                    test_farfield_volume_cannot_hide_engine_residual);
    failures += cfd::test::run_test("Local steady nonuniform tube converges to known solution",
                                    test_convergence_from_a_nonstationary_supersonic_tube);
    failures += cfd::test::run_test("Local steady damped nozzle oscillations",
                                    test_damped_nozzle_oscillations_are_not_a_plateau);
    failures += cfd::test::run_test("Local steady damped refined nozzle oscillations",
                                    test_damped_refined_nozzle_oscillations_are_not_stagnation);
    failures += cfd::test::run_test("Nozzle identical global/local operator and wave sums",
                                    test_nozzle_global_and_local_use_the_same_residual);
    return cfd::test::finish_tests(failures, "axisymmetric Euler local steady");
}
