#include "cfd/compressible/AxisymmetricEulerSolver.hpp"
#include "cfd/compressible/EulerCharacteristicBoundaryFlux.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "support/EulerBoundaryWave.hpp"
#include "support/EulerTube.hpp"
#include "support/TestUtils.hpp"

#include <cmath>
#include <limits>

namespace
{
using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws_with_message;
using BC = cfd::EulerBoundaryCondition;
using Type = cfd::EulerBoundaryConditionType;
const cfd::CaloricallyPerfectGas gas{1.4, 1.0};
constexpr double tolerance{3.0e-12};

double total_temperature(const cfd::PrimitiveState state)
{
    return state.pressure / state.density +
           (state.velocity_x * state.velocity_x + state.velocity_r * state.velocity_r) / 7.0;
}

double total_pressure(const cfd::PrimitiveState state)
{
    return state.pressure * std::pow(total_temperature(state) / (state.pressure / state.density), 3.5);
}

void require_trace_and_flux(const cfd::EulerCharacteristicBoundaryFluxResult &actual,
                            const cfd::PrimitiveState expected, const cfd::Vector2 normal)
{
    const double un{expected.velocity_x * normal.x + expected.velocity_r * normal.y};
    const double energy{expected.pressure / 0.4 +
                        0.5 * expected.density *
                            (expected.velocity_x * expected.velocity_x + expected.velocity_r * expected.velocity_r)};
    require_near(actual.boundary_state.density, expected.density, tolerance, "Trace density mismatch.");
    require_near(actual.boundary_state.momentum_x, expected.density * expected.velocity_x, tolerance,
                 "Trace x momentum.");
    require_near(actual.boundary_state.momentum_r, expected.density * expected.velocity_r, tolerance,
                 "Trace r momentum.");
    require_near(actual.boundary_state.total_energy_density, energy, tolerance, "Trace total energy.");
    require_near(actual.flux.mass, expected.density * un, tolerance, "Actual boundary mass flux.");
    require_near(actual.flux.momentum_x, expected.density * expected.velocity_x * un + expected.pressure * normal.x,
                 tolerance, "Actual boundary x momentum flux.");
    require_near(actual.flux.momentum_r, expected.density * expected.velocity_r * un + expected.pressure * normal.y,
                 tolerance, "Actual boundary r momentum flux.");
    require_near(actual.flux.energy, (energy + expected.pressure) * un, tolerance, "Actual boundary energy flux.");
    const double sound{std::sqrt(1.4 * expected.pressure / expected.density)};
    require_near(actual.normal_mach, un / sound, tolerance, "Normal Mach uses the wrong projection.");
    require(actual.maximum_wave_speed >= (std::abs(un) + sound) * (1.0 - tolerance), "Trace CFL bound too small.");
}

void test_uniform_subsonic_states_and_oblique_normals()
{
    for (const auto normal : {cfd::Vector2{-1.0, 0.0}, cfd::Vector2{0.6, 0.8}, cfd::Vector2{-0.8, 0.6}})
        for (const double mach : {0.05, 0.3, 0.8, 0.999999})
        {
            const double speed{mach * std::sqrt(1.4 * 1.2 / 0.9)};
            const cfd::PrimitiveState expected{0.9, -speed * normal.x, -speed * normal.y, 1.2};
            const auto condition{BC::total_conditions_inlet(total_pressure(expected), total_temperature(expected),
                                                            {-10.0 * normal.x, -10.0 * normal.y})};
            const auto inside{gas.to_conservative(expected)};
            require_trace_and_flux(cfd::compute_euler_characteristic_boundary_flux(inside, normal, condition, gas),
                                   expected, normal);
            require_trace_and_flux(cfd::compute_euler_characteristic_boundary_flux(
                                       inside, {-normal.x, -normal.y}, BC::pressure_outlet(expected.pressure), gas),
                                   expected, {-normal.x, -normal.y});
        }
    // Total Mach > 1, yet normally subsonic: a tangential component must not
    // change how many characteristics enter through this oblique inlet.
    const cfd::PrimitiveState oblique{1.0, 2.0, 0.0, 1.0};
    const cfd::Vector2 normal{-0.2, std::sqrt(0.96)};
    require_trace_and_flux(
        cfd::compute_euler_characteristic_boundary_flux(
            gas.to_conservative(oblique), normal,
            BC::total_conditions_inlet(total_pressure(oblique), total_temperature(oblique), {1.0, 0.0}), gas),
        oblique, normal);
}

void test_nonuniform_characteristic_reconstruction()
{
    // Independent target trace: a=1, un=-1/4, gamma=1.4. Choose the
    // interior sound speed 1.05 so that J+=4.75 is unchanged.
    const auto inside{gas.to_conservative({1.0, -0.5, 0.2, 1.05 * 1.05 / 1.4})};
    const cfd::PrimitiveState target{1.2, -0.25, 0.0, 1.2 / 1.4};
    const auto inlet{cfd::compute_euler_characteristic_boundary_flux(
        inside, {1.0, 0.0}, BC::total_conditions_inlet(total_pressure(target), total_temperature(target), {-1.0, 0.0}),
        gas)};
    require_trace_and_flux(inlet, target, {1.0, 0.0});
    const auto interior{gas.to_conservative({1.0, 0.4, 0.7, 1.0})};
    constexpr double pb{0.9};
    const double rho{std::pow(pb, 1.0 / 1.4)}, ab{std::sqrt(1.4 * pb / rho)};
    const cfd::PrimitiveState outlet{rho, 0.4 + 5.0 * (std::sqrt(1.4) - ab), 0.7, pb};
    require_trace_and_flux(
        cfd::compute_euler_characteristic_boundary_flux(interior, {1.0, 0.0}, BC::pressure_outlet(pb), gas), outlet,
        {1.0, 0.0});
    // Rotated nonuniform outlet: preserve the tangential velocity, not the
    // individual Cartesian components.
    const cfd::Vector2 n{0.6, 0.8};
    const cfd::PrimitiveState rotated{1.0, 0.4 * n.x - 0.7 * n.y, 0.4 * n.y + 0.7 * n.x, 1.0};
    const cfd::PrimitiveState rotated_out{rho, outlet.velocity_x * n.x - 0.7 * n.y, outlet.velocity_x * n.y + 0.7 * n.x,
                                          pb};
    require_trace_and_flux(
        cfd::compute_euler_characteristic_boundary_flux(gas.to_conservative(rotated), n, BC::pressure_outlet(pb), gas),
        rotated_out, n);
    // Supersonic total speed, but subsonic normal outflow: the prescribed
    // pressure must still act. A total-Mach switch would incorrectly ignore it.
    const cfd::Vector2 grazing{0.2, std::sqrt(0.96)};
    const double increment{5.0 * (std::sqrt(1.4) - ab)};
    require_trace_and_flux(cfd::compute_euler_characteristic_boundary_flux(gas.to_conservative({1.0, 2.0, 0.0, 1.0}),
                                                                           grazing, BC::pressure_outlet(pb), gas),
                           {rho, 2.0 + increment * grazing.x, increment * grazing.y, pb}, grazing);
}

void test_supersonic_and_pressure_ignored()
{
    const cfd::PrimitiveState prescribed{1.0, -2.0, 0.7, 1.0};
    const auto inside{gas.to_conservative({0.8, -2.5, 0.2, 0.7})};
    require_trace_and_flux(
        cfd::compute_euler_characteristic_boundary_flux(inside, {1.0, 0.0}, BC::supersonic_inlet(prescribed), gas),
        prescribed, {1.0, 0.0});
    const cfd::PrimitiveState outgoing{1.0, 2.0, 0.7, 1.0};
    for (const double pressure : {1.0e-200, 0.5, 1.0e200})
        require_trace_and_flux(cfd::compute_euler_characteristic_boundary_flux(
                                   gas.to_conservative(outgoing), {1.0, 0.0}, BC::pressure_outlet(pressure), gas),
                               outgoing, {1.0, 0.0});
    require_trace_and_flux(cfd::compute_euler_characteristic_boundary_flux(gas.to_conservative(outgoing), {1.0, 0.0},
                                                                           Type::SupersonicOutlet, gas),
                           outgoing, {1.0, 0.0});
    const cfd::Vector2 n{0.6, 0.8};
    const cfd::PrimitiveState oblique{1.0, -2.0 * n.x - 0.7 * n.y, -2.0 * n.y + 0.7 * n.x, 1.0};
    require_trace_and_flux(cfd::compute_euler_characteristic_boundary_flux(gas.to_conservative(oblique), n,
                                                                           BC::supersonic_inlet(oblique), gas),
                           oblique, n);
    require_throws_with_message<std::runtime_error>(
        [&] {
            (void)cfd::compute_euler_characteristic_boundary_flux(
                gas.to_conservative(outgoing), {-0.2, std::sqrt(0.96)}, BC::supersonic_inlet(outgoing), gas);
        },
        "supersonic inlet", "Total Mach was used instead of normal Mach at supersonic inlet.");
}

void test_mass_flow_unequal_areas_and_normals()
{
    for (const cfd::Index nr : {2U, 5U, 9U})
    {
        auto raw{cfd::verification::make_euler_tube(2, nr)};
        for (auto &node : raw.nodes)
        {
            node.y *= node.y;
            node.x += 0.2 * node.y * node.y;
        }
        const auto built{cfd::build_mesh(std::move(raw))};
        const cfd::AxisymmetricGeometry geometry{built.mesh};
        double area{};
        for (cfd::Index face = 0; face < built.mesh.face_count(); ++face)
            if (built.mesh.face_adjacencies()[face].is_boundary() && built.mesh.face_boundary_ids()[face] == 2)
            {
                const auto vector{geometry.face_area_vectors()[face]};
                area += std::hypot(vector.x, vector.y);
            }
        const auto condition{BC::mass_flow_inlet(0.7, 1.0 + 0.09 / 7.0, {3.0, 0.0})};
        const auto inside{gas.to_conservative({1.0, 0.3, 0.0, 1.0})};
        double integrated{};
        for (cfd::Index face = 0; face < built.mesh.face_count(); ++face)
        {
            if (!built.mesh.face_adjacencies()[face].is_boundary() || built.mesh.face_boundary_ids()[face] != 2)
                continue;
            const auto sf{built.mesh.face_area_vectors()[face]};
            const double length{built.mesh.face_lengths()[face]};
            const auto vector{geometry.face_area_vectors()[face]};
            const double surface{std::hypot(vector.x, vector.y)};
            const auto result{cfd::compute_euler_characteristic_boundary_flux(inside, {sf.x / length, sf.y / length},
                                                                              condition, gas, area)};
            require_near(result.flux.mass, -0.7 / area, tolerance,
                         "Nonuniform group did not get uniform normal mass flux.");
            integrated += surface * result.flux.mass;
        }
        require_near(integrated, -0.7, tolerance,
                     "Group mass inlet used meridional lengths instead of revolution areas.");
        // Also exercise the actual solver's once-built group area and global
        // flux ledger. The right boundary is a wall, so total mass outflow is
        // exactly the inlet contribution, independently of the interior update.
        cfd::AxisymmetricEulerSolver solver{built.mesh,
                                            cfd::ConservativeField{built.mesh.cell_count(), inside},
                                            gas,
                                            {4, {Type::Axis, Type::SlipWall, condition, Type::SlipWall}}};
        const auto before{solver.diagnostics().integrals.mass};
        const auto result{solver.step(1.0e-5)};
        require_near(result.boundary_outflow_rate.mass, -0.7, tolerance, "Solver did not assemble imposed group mdot.");
        require_near(result.integrals.mass - before, 0.7 * result.time_step, tolerance,
                     "Inlet mass was not conserved.");
        require_near(result.balance_error.energy, 0.0, tolerance, "Mass inlet energy ledger mismatch.");
    }
}

void test_uniform_tube_evolution()
{
    const auto built{cfd::build_mesh(cfd::verification::make_euler_tube(6, 3))};
    for (const double mach : {0.05, 0.3, 0.8, 2.0})
    {
        const cfd::PrimitiveState primitive{1.0, mach * std::sqrt(1.4), 0.0, 1.0};
        const auto value{gas.to_conservative(primitive)};
        const auto inlet{
            mach < 1.0 ? BC::total_conditions_inlet(total_pressure(primitive), total_temperature(primitive), {1.0, 0.0})
                       : BC::supersonic_inlet(primitive)};
        cfd::AxisymmetricEulerSolver solver{built.mesh,
                                            cfd::ConservativeField{built.mesh.cell_count(), value},
                                            gas,
                                            {4, {Type::Axis, Type::SlipWall, inlet, BC::pressure_outlet(1.0)}}};
        solver.advance_to(0.1);
        require(solver.diagnostics().rejected_steps == 0, "Uniform characteristic tube rejected steps.");
        for (cfd::Index cell = 0; cell < built.mesh.cell_count(); ++cell)
        {
            const auto actual{solver.state().state(cell)};
            require_near(actual.density, value.density, tolerance, "Uniform density changed.");
            require_near(actual.momentum_x, value.momentum_x, tolerance, "Uniform axial momentum changed.");
            require_near(actual.momentum_r, 0.0, tolerance, "Uniform tube gained radial momentum.");
            require_near(actual.total_energy_density, value.total_energy_density, tolerance,
                         "Uniform total energy changed.");
        }
    }
}

void test_regime_changes_sonic_limits_and_backflow()
{
    const auto outlet{BC::pressure_outlet(1.0)};
    // gamma=2, rho=1,p=1/2 gives exactly representable sound speed=1.
    const cfd::CaloricallyPerfectGas model{2.0, 1.0};
    const auto rest{model.to_conservative({1.0, 0.0, 0.0, 0.5})};
    const auto stagnant{cfd::compute_euler_characteristic_boundary_flux(
        rest, {-1.0, 0.0}, BC::total_conditions_inlet(0.5, 0.5, {1.0, 0.0}), model)};
    require(stagnant.flux.mass == 0.0 && stagnant.flux.energy == 0.0,
            "Stagnant inlet limit gained a mass/energy flux.");
    require_throws_with_message<std::runtime_error>(
        [&] {
            (void)cfd::compute_euler_characteristic_boundary_flux(
                rest, {-1.0, 0.0}, BC::mass_flow_inlet(1.0, 0.5, {1.0, 0.0}), model, 1.0);
        },
        "zero inflow speed", "Positive mdot was imposed at zero characteristic inflow speed.");
    for (const double speed : {std::nextafter(-1.0, -2.0), -1.0})
        require_throws_with_message<std::runtime_error>(
            [&] {
                (void)cfd::compute_euler_characteristic_boundary_flux(
                    model.to_conservative({1.0, speed, 0.0, 0.5}), {1.0, 0.0},
                    BC::total_conditions_inlet(1.0, 1.0, {-1.0, 0.0}), model);
            },
            "sonic inflow", "Normal sonic inlet was silently clipped.");
    for (const double speed : {std::nextafter(1.0, 0.0), 1.0, std::nextafter(1.0, 2.0)})
    {
        const auto result{cfd::compute_euler_characteristic_boundary_flux(model.to_conservative({1.0, speed, 0.0, 0.5}),
                                                                          {1.0, 0.0}, BC::pressure_outlet(0.5), model)};
        require_near(result.flux.mass, speed, tolerance, "Sonic outlet limit discontinuity.");
    }
    for (const double speed : {0.0, 1.0e-14, 0.3, 2.0})
        (void)cfd::compute_euler_characteristic_boundary_flux(gas.to_conservative({1.0, speed, 0.0, 1.0}), {1.0, 0.0},
                                                              outlet, gas);
    const auto reverse{gas.to_conservative({1.0, -0.1, 0.0, 1.0})};
    require_throws_with_message<std::runtime_error>(
        [&] { (void)cfd::compute_euler_characteristic_boundary_flux(reverse, {1.0, 0.0}, outlet, gas); }, "backflow",
        "Outlet reverse flow silently extrapolated.");
    const auto inward{gas.to_conservative({1.0, -0.3, 0.0, 1.0})};
    const auto inlet{BC::total_conditions_inlet(1.1, 1.1, {-1.0, 0.0})};
    require_throws_with_message<std::runtime_error>(
        [&] {
            (void)cfd::compute_euler_characteristic_boundary_flux(reverse, {1.0, 0.0}, Type::SupersonicOutlet, gas);
        },
        "supersonic outlet", "Invalid supersonic outlet accepted.");
    require_throws_with_message<std::runtime_error>(
        [&] {
            (void)cfd::compute_euler_characteristic_boundary_flux(inward, {1.0, 0.0},
                                                                  BC::supersonic_inlet({1.0, -2.0, 0.0, 1.0}), gas);
        },
        "supersonic inlet", "Subsonic interior passed supersonic inlet.");
    require_throws_with_message<std::runtime_error>(
        [&] {
            (void)cfd::compute_euler_characteristic_boundary_flux(gas.to_conservative({1.0, 0.1, 0.0, 1.0}), {1.0, 0.0},
                                                                  inlet, gas);
        },
        "reverse flow", "Inlet changing to outflow was silently accepted.");
    require_throws_with_message<std::runtime_error>(
        [&] {
            (void)cfd::compute_euler_characteristic_boundary_flux(gas.to_conservative({1.0, 0.3, 0.0, 1.0}), {1.0, 0.0},
                                                                  BC::pressure_outlet(1.0e-20), gas);
        },
        "sonic/supersonic", "Low outlet pressure was silently clipped.");
    require_throws_with_message<std::runtime_error>(
        [&] {
            (void)cfd::compute_euler_characteristic_boundary_flux(gas.to_conservative({1.0, 0.3, 0.0, 1.0}), {1.0, 0.0},
                                                                  BC::pressure_outlet(10.0), gas);
        },
        "backflow", "Pressure-induced backflow accepted without inlet data.");
    require_throws_with_message<std::runtime_error>(
        [&] {
            (void)cfd::compute_euler_characteristic_boundary_flux(
                inward, {1.0, 0.0}, BC::total_conditions_inlet(1.0, 0.001, {-1.0, 0.0}), gas);
        },
        "root", "Incompatible total inlet was clipped to a root.");
}

void test_total_inlet_stagnation_existence_boundary()
{
    // gamma=2, R=1, T0=2 give a0=2 and Jmax=4 exactly.
    const cfd::CaloricallyPerfectGas model{2.0, 1.0};
    const auto condition{BC::total_conditions_inlet(2.0, 2.0, {1.0, 0.0})};
    const cfd::Vector2 normal{-1.0, 0.0};
    for (const double speed : {1.0e-8, 4.0 - std::nextafter(4.0, 0.0)})
    {
        const auto inside{model.to_conservative({1.0, speed, 0.0, 2.0})};
        const auto result{cfd::compute_euler_characteristic_boundary_flux(inside, normal, condition, model)};
        const auto trace{model.evaluate(result.boundary_state)};
        require(trace.primitive.velocity_x > 0.0 && result.flux.mass < 0.0,
                "Compatible near-stagnant inlet was replaced by zero or reverse flow.");
        require_near(trace.primitive.velocity_x, speed, 4.0 * std::numeric_limits<double>::epsilon(),
                     "Near-stagnant trace has the wrong velocity scale.");
        require_near(trace.thermodynamics.temperature + 0.25 * trace.primitive.velocity_x * trace.primitive.velocity_x,
                     2.0, tolerance, "Near-stagnant trace does not preserve total temperature.");
        require_near(trace.primitive.pressure * std::pow(2.0 / trace.thermodynamics.temperature, 2.0), 2.0, tolerance,
                     "Near-stagnant trace does not preserve total pressure.");
    }
    const auto rest{model.to_conservative({1.0, 0.0, 0.0, 2.0})};
    const auto stagnant{cfd::compute_euler_characteristic_boundary_flux(rest, normal, condition, model)};
    require(stagnant.boundary_state.momentum_x == 0.0 && stagnant.flux.mass == 0.0 && stagnant.flux.energy == 0.0,
            "The exact stagnation endpoint must have zero velocity, mass flux and energy flux.");
    require_near(stagnant.flux.momentum_x, -2.0, 0.0, "Stagnation pressure force has the wrong sign.");

    const double larger_sound{std::nextafter(2.0, 3.0)};
    const auto above{model.to_conservative({1.0, 0.0, 0.0, 0.5 * larger_sound * larger_sound})};
    require(model.evaluate(above).thermodynamics.sound_speed > 2.0, "Test invariant is not above Jmax.");
    require_throws_with_message<std::runtime_error>(
        [&] { (void)cfd::compute_euler_characteristic_boundary_flux(above, normal, condition, model); }, "root",
        "An invariant above the stagnation bound was accepted.");
}

void test_total_inlet_trace_sonic_existence_boundary()
{
    // With gamma=2 and T0=3/4, qsonic=asonic=Jsonic=1.
    // Interior un=-1/2, a=3/4 is strictly subsonic but has J=1.
    const cfd::CaloricallyPerfectGas model{2.0, 1.0};
    const auto condition{BC::total_conditions_inlet(1.0, 0.75, {1.0, 0.0})};
    const cfd::Vector2 normal{-1.0, 0.0};
    const double delta{std::ldexp(1.0, -20)};
    const auto compatible{model.to_conservative({1.0, 0.5 - delta, 0.0, 0.28125})};
    const auto result{cfd::compute_euler_characteristic_boundary_flux(compatible, normal, condition, model)};
    require(result.normal_mach > -1.0 && result.normal_mach < 0.0,
            "An invariant above Jsonic must admit a strictly subsonic incoming trace.");
    for (const double speed : {0.5, 0.5 + delta})
    {
        const auto inside{model.to_conservative({1.0, speed, 0.0, 0.28125})};
        require_throws_with_message<std::runtime_error>(
            [&] { (void)cfd::compute_euler_characteristic_boundary_flux(inside, normal, condition, model); }, "root",
            "A sonic or missing trace root was accepted for a strictly subsonic interior.");
    }
}

void test_total_inlet_f4_saved_stagnation_incompatibility()
{
    // Binary64 conservative cell 0 from the saved F4 field, not rounded primitives.
    const cfd::ConservativeState inside{0x1.b61192a1ce633p-1, 0x1.ee591d4065b6dp+0, 0x1.15d1a636347e2p-3,
                                        0x1.998db681d37aep+22};
    const cfd::CaloricallyPerfectGas model{1.177, 414.2};
    const auto evaluated{model.evaluate(inside)};
    require(evaluated.primitive.velocity_x > 0.0 &&
                evaluated.primitive.velocity_x < evaluated.thermodynamics.sound_speed,
            "Saved F4 interior must still be normal-subsonic inflow.");
    require_throws_with_message<std::runtime_error>(
        [&] {
            (void)cfd::compute_euler_characteristic_boundary_flux(
                inside, {-1.0, 0.0}, BC::total_conditions_inlet(1185747.2501228245, 3350.2948821164091, {1.0, 0.0}),
                model);
        },
        "root", "Saved incompatible F4 state was silently accepted.");
}

void test_characteristic_compatibility_errors_are_distinct_from_fatal_errors()
{
    const cfd::CaloricallyPerfectGas model{2, 1};
    const auto rest{model.to_conservative({1, 0, 0, 2})};
    try
    {
        (void)cfd::compute_euler_characteristic_boundary_flux(rest, {-1, 0}, BC::total_conditions_inlet(2, .5, {1, 0}),
                                                              model);
        cfd::test::fail("Missing recoverable stagnation incompatibility.");
    }
    catch (const cfd::EulerCharacteristicCompatibilityError &error)
    {
        cfd::test::require_contains(error.what(), "stagnation bound", "Missing stagnation diagnostic.");
        for (const char *text : {"J=", "J0=", "Js=", "f(0)=", "f(qs)="})
            cfd::test::require_contains(error.what(), text, "Missing invariant/bound in root diagnostic.");
    }
    require_throws_with_message<cfd::EulerCharacteristicCompatibilityError>(
        [&] {
            (void)cfd::compute_euler_characteristic_boundary_flux(model.to_conservative({1, .5, 0, .28125}), {-1, 0},
                                                                  BC::total_conditions_inlet(1, .75, {1, 0}), model);
        },
        "sonic bound", "Sonic trace did not produce a distinct compatibility diagnostic.");
    require_throws_with_message<cfd::EulerCharacteristicCompatibilityError>(
        [&] {
            (void)cfd::compute_euler_characteristic_boundary_flux(model.to_conservative({1, -.25, 0, 2}), {-1, 0},
                                                                  BC::total_conditions_inlet(2, 2, {1, 0}), model);
        },
        "interior regime", "Reverse flow did not identify the incompatible interior regime.");
    // A valid stagnant root would require rho=2*DBL_MAX: arithmetic, not a missing root.
    try
    {
        (void)cfd::compute_euler_characteristic_boundary_flux(
            model.to_conservative({1, 0, 0, .5}), {-1, 0},
            BC::total_conditions_inlet(std::numeric_limits<double>::max(), .5, {1, 0}), model);
        cfd::test::fail("Unrepresentable trace density was accepted.");
    }
    catch (const cfd::EulerCharacteristicCompatibilityError &)
    {
        cfd::test::fail("An arithmetic failure was incorrectly classified as recoverable.");
    }
    catch (const std::runtime_error &error)
    {
        cfd::test::require_contains(error.what(), "non-representable", "Missing arithmetic diagnostic.");
    }
}

void test_invalid_data_and_direction()
{
    for (const double invalid :
         {0.0, -1.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()})
    {
        require_throws_with_message<std::invalid_argument>([&] { (void)BC::pressure_outlet(invalid); }, "positive",
                                                           "Invalid pressure accepted.");
        require_throws_with_message<std::invalid_argument>(
            [&] { (void)BC::total_conditions_inlet(1.0, invalid, {1.0, 0.0}); }, "positive", "Invalid T0 accepted.");
        require_throws_with_message<std::invalid_argument>([&] { (void)BC::mass_flow_inlet(invalid, 1.0, {1.0, 0.0}); },
                                                           "positive", "Invalid mdot accepted.");
    }
    for (const auto d : {cfd::Vector2{}, cfd::Vector2{std::numeric_limits<double>::quiet_NaN(), 0.0},
                         cfd::Vector2{0.0, std::numeric_limits<double>::infinity()}})
        require_throws_with_message<std::invalid_argument>([&] { (void)BC::total_conditions_inlet(1.0, 1.0, d); },
                                                           "direction", "Invalid direction accepted.");
    const auto huge{BC::mass_flow_inlet(1.0, 1.0, {1.0e308, 1.0e308})};
    require_near(std::hypot(huge.direction().x, huge.direction().y), 1.0, tolerance,
                 "Direction normalization overflowed.");
    const auto inside{gas.to_conservative({1.0, -0.3, 0.0, 1.0})};
    require_throws_with_message<std::runtime_error>(
        [&] {
            (void)cfd::compute_euler_characteristic_boundary_flux(inside, {1.0, 0.0},
                                                                  BC::mass_flow_inlet(1.0, 1.1, {1.0, 0.0}), gas, 1.0);
        },
        "direction", "Outward inlet direction accepted.");
    require_throws_with_message<std::invalid_argument>(
        [&] {
            (void)cfd::compute_euler_characteristic_boundary_flux(inside, {1.0, 0.0},
                                                                  BC::mass_flow_inlet(1.0, 1.1, {-1.0, 0.0}), gas);
        },
        "area", "Missing group area accepted.");
    require_throws_with_message<std::invalid_argument>(
        [&] {
            (void)cfd::compute_euler_characteristic_boundary_flux(inside, {2.0, 0.0}, BC::pressure_outlet(1.0), gas);
        },
        "unit normal", "Invalid normal accepted.");
    require_throws_with_message<std::invalid_argument>([] { const BC missing{Type::TotalConditionsInlet}; },
                                                       "physical data", "Under-specified BC accepted.");
    require_throws_with_message<std::runtime_error>(
        [&] { (void)cfd::compute_euler_characteristic_boundary_flux({}, {1.0, 0.0}, BC::pressure_outlet(1.0), gas); },
        "Perfect gas", "Invalid interior accepted.");
}

void test_acoustic_pressure_reflection_against_extended_domain()
{
    const auto result{cfd::verification::compare_boundary_wave(60)};
    require(result.result.minimum_density > 0.0 && result.result.minimum_pressure > 0.0 &&
                result.result.rejected_steps == 0,
            "Acoustic pressure boundary lost admissibility.");
    require(result.maximum_balance_defect < 1.0e-12 && result.maximum_step_balance_defect < 1.0e-12,
            "Acoustic boundary violates global conservation.");
    // A fixed pressure boundary reflects pressure with reversed sign in
    // linear acoustics. Diffusion reduces the measured amplitude, but a
    // substantial reflected wave, absent in the extended domain, must exist.
    require(result.reflection_ratio > 0.05 && result.reflection_ratio < 1.5,
            "Fixed-pressure reflection is absent or spuriously amplified.");
    require(result.minimum_pressure_difference < -0.05 * result.incident_pressure,
            "Fixed-pressure reflected pressure must reverse the incident compression sign.");
    require(result.reflected_incoming_amplitude > 0.05 * result.incident_pressure,
            "No incoming acoustic characteristic measured after boundary interaction.");
}
} // namespace

int main()
{
    int failures{};
    failures += cfd::test::run_test("uniform subsonic and oblique characteristic fluxes",
                                    test_uniform_subsonic_states_and_oblique_normals);
    failures +=
        cfd::test::run_test("nonuniform characteristic reconstruction", test_nonuniform_characteristic_reconstruction);
    failures +=
        cfd::test::run_test("supersonic inflow and pressure-independent outflow", test_supersonic_and_pressure_ignored);
    failures += cfd::test::run_test("mass flow on unequal revolution areas and normals",
                                    test_mass_flow_unequal_areas_and_normals);
    failures += cfd::test::run_test("uniform characteristic tube evolution", test_uniform_tube_evolution);
    failures +=
        cfd::test::run_test("regime changes sonic limits and backflow", test_regime_changes_sonic_limits_and_backflow);
    failures += cfd::test::run_test("total inlet stagnation existence boundary",
                                    test_total_inlet_stagnation_existence_boundary);
    failures += cfd::test::run_test("total inlet trace sonic existence boundary",
                                    test_total_inlet_trace_sonic_existence_boundary);
    failures += cfd::test::run_test("saved F4 total inlet stagnation incompatibility",
                                    test_total_inlet_f4_saved_stagnation_incompatibility);
    failures += cfd::test::run_test("invalid characteristic data", test_invalid_data_and_direction);
    failures += cfd::test::run_test("characteristic compatibility versus fatal errors",
                                    test_characteristic_compatibility_errors_are_distinct_from_fatal_errors);
    failures += cfd::test::run_test("acoustic reflection against extended domain",
                                    test_acoustic_pressure_reflection_against_extended_domain);
    return cfd::test::finish_tests(failures, "Euler characteristic boundary");
}
