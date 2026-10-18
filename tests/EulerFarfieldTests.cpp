#include "cfd/compressible/AxisymmetricEulerSolver.hpp"
#include "cfd/compressible/EulerFarfieldBoundaryFlux.hpp"
#include "support/EulerBoundaryWave.hpp"
#include "support/TestUtils.hpp"

#include <array>
#include <cmath>
#include <limits>

namespace
{
using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws;
using BC = cfd::EulerBoundaryCondition;
using Regime = cfd::EulerFarfieldRegime;
const cfd::CaloricallyPerfectGas gas{1.4, 1.0};
constexpr double tolerance{3.0e-12};

auto evaluate(const cfd::PrimitiveState interior, const cfd::PrimitiveState exterior,
              const cfd::Vector2 normal = {1.0, 0.0})
{
    return cfd::compute_euler_farfield_boundary_flux(gas.to_conservative(interior), normal, BC::farfield(exterior),
                                                     gas);
}

void require_flux(const cfd::EulerFarfieldBoundaryFluxResult &actual, const cfd::PrimitiveState expected,
                  const cfd::Vector2 n)
{
    // Independent primitive flux reference, not the production flux evaluator.
    const double un{expected.velocity_x * n.x + expected.velocity_r * n.y};
    const double energy{expected.pressure / 0.4 +
                        0.5 * expected.density *
                            (expected.velocity_x * expected.velocity_x + expected.velocity_r * expected.velocity_r)};
    require_near(actual.flux.mass, expected.density * un, tolerance, "Mass flux.");
    require_near(actual.flux.momentum_x, expected.density * expected.velocity_x * un + expected.pressure * n.x,
                 tolerance, "Axial flux.");
    require_near(actual.flux.momentum_r, expected.density * expected.velocity_r * un + expected.pressure * n.y,
                 tolerance, "Radial flux.");
    require_near(actual.flux.energy, (energy + expected.pressure) * un, tolerance, "Energy flux.");
    const double sound{std::sqrt(1.4 * expected.pressure / expected.density)};
    require_near(actual.normal_mach, un / sound, tolerance, "Trace normal Mach.");
    require(actual.maximum_wave_speed >= std::abs(un) + sound - tolerance, "CFL misses trace speed.");
}

void test_uniform_states_and_normal_not_total_mach()
{
    for (const auto n : {cfd::Vector2{1.0, 0.0}, cfd::Vector2{-1.0, 0.0}, cfd::Vector2{0.6, 0.8}})
        for (const double mach : {-2.0, -0.3, 0.0, 0.3, 2.0})
        {
            const double un{mach * std::sqrt(1.4)};
            const cfd::PrimitiveState p{1.0, un * n.x - 3.0 * n.y, un * n.y + 3.0 * n.x, 1.0};
            const auto result{evaluate(p, p, n)};
            require_flux(result, p, n);
            const auto expected{mach <= -1.0 ? Regime::SupersonicInflow
                                : mach < 0.0 ? Regime::SubsonicInflow
                                : mach < 1.0 ? Regime::SubsonicOutflow
                                             : Regime::SupersonicOutflow};
            if (mach != 0.0)
                require(result.regime == expected, "Regime incorrectly uses total Mach.");
            else
                require(std::abs(result.normal_mach) < tolerance &&
                            (result.regime == Regime::SubsonicInflow || result.regime == Regime::SubsonicOutflow),
                        "Tangential roundoff must remain subsonic, without a sign-clipping threshold.");
        }
    require_flux(evaluate({1.0, 0.0, 0.0, 1.0}, {1.0, 0.0, 0.0, 1.0}), {1.0, 0.0, 0.0, 1.0}, {1.0, 0.0});
    // A root arbitrarily close to log(p*)=0 need not resolve in log space:
    // stop when pressure itself is resolved to adjacent representable values.
    require_flux(evaluate({1.0, 0.3 + 1.0e-15, 0.0, 1.0}, {1.0, 0.3, 0.0, 1.0}), {1.0, 0.3, 0.0, 1.0}, {1.0, 0.0});
}

void test_supersonic_causality_and_incoming_shock()
{
    const cfd::PrimitiveState inside{1.0, 8.0, 0.2, 1.0};
    for (const double pressure : {0.8, 1.0, 1.2})
    {
        const auto result{evaluate(inside, {0.9, 8.0, -0.7, pressure})};
        require(result.regime == Regime::SupersonicOutflow, "Expected all-outgoing branch.");
        require_flux(result, inside, {1.0, 0.0});
        const cfd::PrimitiveState outside{0.9, -8.0, -0.7, pressure};
        const auto incoming{evaluate({1.0, -8.0, 0.2, 1.0}, outside)};
        require(incoming.regime == Regime::SupersonicInflow, "Expected all-incoming branch.");
        require_flux(incoming, outside, {1.0, 0.0});
    }
    const auto strong{evaluate({1.0, 2.0, 0.0, 1.0}, {1.0, 0.0, 0.0, 100.0})};
    require(strong.flux.mass < 0.0, "An incoming shock must not be classified using interior Mach alone.");
}

void test_acoustic_invariants_entropy_tangent_and_backflow()
{
    constexpr double pressure{0.8};
    const double a_left{std::sqrt(1.4)}, a_right{std::sqrt(0.7)};
    const double factor{std::pow(pressure, 1.0 / 7.0)};
    for (const double speed : {-0.2, 0.2})
    {
        // Two independent simple waves meet at p*=0.8, un*=speed.
        const cfd::PrimitiveState left{1.0, speed + 5.0 * a_left * (factor - 1.0), 0.6, 1.0};
        const cfd::PrimitiveState right{2.0, speed - 5.0 * a_right * (factor - 1.0), -0.4, 1.0};
        const auto result{evaluate(left, right)};
        const auto p{gas.evaluate(result.boundary_state).primitive};
        const double source_density{speed > 0.0 ? 1.0 : 2.0};
        require_near(p.density, source_density * std::pow(pressure, 5.0 / 7.0), tolerance, "Entropy side.");
        require_near(p.pressure, pressure, tolerance, "Simple-wave pressure.");
        require_near(p.velocity_x, speed, tolerance, "Normal contact speed.");
        require_near(p.velocity_r, speed > 0.0 ? 0.6 : -0.4, tolerance, "Incoming/outgoing tangent.");
        require_near(p.velocity_x + 5.0 * std::sqrt(1.4 * pressure / std::pow(pressure, 5.0 / 7.0)),
                     left.velocity_x + 5.0 * a_left, tolerance, "Left acoustic invariant.");
        require_near(p.velocity_x - 5.0 * a_right * factor, right.velocity_x - 5.0 * a_right, tolerance,
                     "Right acoustic invariant.");
        const cfd::Vector2 n{0.6, 0.8};
        const cfd::PrimitiveState rotated_left{left.density, left.velocity_x * n.x - left.velocity_r * n.y,
                                               left.velocity_x * n.y + left.velocity_r * n.x, left.pressure};
        const cfd::PrimitiveState rotated_right{right.density, right.velocity_x * n.x - right.velocity_r * n.y,
                                                right.velocity_x * n.y + right.velocity_r * n.x, right.pressure};
        const double tangent{speed > 0.0 ? 0.6 : -0.4};
        require_flux(evaluate(rotated_left, rotated_right, n),
                     {source_density * std::pow(pressure, 5.0 / 7.0), speed * n.x - tangent * n.y,
                      speed * n.y + tangent * n.x, pressure},
                     n);
    }
    // Entropy and tangent differ strongly: switching a stationary contact must
    // not jump pressure flux. At zero speed the density trace may be discontinuous.
    for (const double speed : {-1.0e-9, 0.0, 1.0e-9})
    {
        const auto result{evaluate({1.0, speed, 0.6, 1.0}, {3.0, speed, -0.4, 1.0})};
        const auto p{gas.evaluate(result.boundary_state).primitive};
        require_near(p.density, speed < 0.0 ? 3.0 : 1.0, tolerance, "Backflow entropy source.");
        require_near(result.flux.momentum_x, 1.0, tolerance, "Contact crossing pressure jump.");
        require(std::abs(result.flux.mass) < 4.0e-9 && std::abs(result.flux.energy) < 1.0e-8,
                "Contact flux is not continuous.");
    }
}

void test_shocks_sonic_fans_and_orientation()
{
    // Published Sod normal Riemann sample at xi=0, independently tabulated.
    const auto sod{evaluate({1.0, 0.0, 0.3, 1.0}, {0.125, 0.0, -0.2, 0.1})};
    const auto p{gas.evaluate(sod.boundary_state).primitive};
    require_near(p.pressure, 0.303130178050647, tolerance, "Sod star pressure.");
    require_near(p.density, 0.426319428178495, tolerance, "Sod star density.");
    require_near(p.velocity_x, 0.927452620048950, tolerance, "Sod contact velocity.");
    require_near(p.velocity_r, 0.3, tolerance, "Sod tangent.");
    const auto reversed{evaluate({0.125, 0.0, -0.2, 0.1}, {1.0, 0.0, 0.3, 1.0}, {-1.0, 0.0})};
    require_near(reversed.flux.mass, -sod.flux.mass, tolerance, "Normal reversal mass.");
    require_near(reversed.flux.momentum_x, -sod.flux.momentum_x, tolerance, "Normal reversal momentum.");
    require_near(reversed.flux.momentum_r, -sod.flux.momentum_r, tolerance, "Normal reversal tangent.");
    require_near(reversed.flux.energy, -sod.flux.energy, tolerance, "Normal reversal energy.");
    // Independent symmetric compression: p*=1+rho*w*(k*w+hypot(a,k*w)).
    const double sound{std::sqrt(1.4)}, w{0.4}, k{0.6};
    const auto compressed{evaluate({1.0, w, 0.0, 1.0}, {1.0, -w, 0.0, 1.0})};
    require_near(compressed.flux.momentum_x, 1.0 + w * (k * w + std::hypot(sound, k * w)), tolerance,
                 "Shock compression pressure.");
    // Crossing a sonic rarefaction head: below sonic the boundary is inside the
    // fan, above it the initial state is sampled. Flux must approach the same limit.
    for (const double sign : {-1.0, 1.0})
    {
        const double un{sound * (1.0 + sign * 1.0e-8)};
        const auto result{evaluate({1.0, un, 0.0, 1.0}, {0.6, 2.0 * sound, 0.0, 0.5})};
        require_near(result.flux.mass, sound, 3.0e-8, "Sonic-head mass discontinuity.");
        require_near(result.flux.momentum_x, 2.4, 3.0e-8, "Sonic-head momentum discontinuity.");
        const auto incoming{evaluate({0.6, -2.0 * sound, 0.0, 0.5}, {1.0, -un, 0.0, 1.0})};
        require_near(incoming.flux.mass, -sound, 3.0e-8, "Incoming sonic-head mass discontinuity.");
        require_near(incoming.flux.momentum_x, 2.4, 3.0e-8, "Incoming sonic-head momentum discontinuity.");
    }
    const auto fan{evaluate({1.0, 0.7 * sound, 0.0, 1.0}, {0.1, 2.0 * sound, 0.0, 0.05})};
    require_near(fan.normal_mach, 1.0, tolerance, "Sonic fan not sampled at zero wave speed.");
}

void test_invalid_and_near_vacuum_states()
{
    const auto valid{gas.to_conservative({1.0, 0.0, 0.0, 1.0})};
    for (const double invalid :
         {0.0, -1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
    {
        require_throws<std::invalid_argument>([&] { (void)BC::farfield({invalid, 0.0, 0.0, 1.0}); }, "Invalid rho.");
        require_throws<std::invalid_argument>([&] { (void)BC::farfield({1.0, 0.0, 0.0, invalid}); }, "Invalid p.");
    }
    require_throws<std::invalid_argument>([&] { (void)BC::farfield({1.0, INFINITY, 0.0, 1.0}); }, "Invalid velocity.");
    require_throws<std::invalid_argument>(
        [&] { (void)cfd::compute_euler_farfield_boundary_flux(valid, {2.0, 0.0}, BC::farfield({1, 0, 0, 1}), gas); },
        "Invalid normal.");
    require_throws<std::invalid_argument>(
        [&] { (void)cfd::compute_euler_farfield_boundary_flux(valid, {1, 0}, BC::pressure_outlet(1.0), gas); },
        "Wrong BC.");
    require_throws<std::runtime_error>([&] { (void)evaluate({1, -6, 0, 1}, {1, 6, 0, 1}); }, "Vacuum must fail.");
    const auto rarefied{evaluate({1, -5, 0, 1}, {1, 5, 0, 1})};
    require(gas.evaluate(rarefied.boundary_state).primitive.pressure > 0.0, "Positive expansion was clipped.");
    require_throws<std::runtime_error>(
        [&] { (void)cfd::compute_euler_farfield_boundary_flux({}, {1, 0}, BC::farfield({1, 0, 0, 1}), gas); },
        "Invalid interior.");
    require_throws<std::runtime_error>([&] { (void)evaluate({1, 0, 0, 1}, {1.0e-308, 0, 0, 1.0e308}); },
                                       "Unrepresentable exterior thermal data.");
    for (const double scale : {1.0e-200, 1.0e200})
    {
        const auto result{evaluate({scale, 0.3, 0.0, scale}, {scale, 0.3, 0.0, scale})};
        require_near(result.flux.mass / scale, 0.3, tolerance, "Scaled uniform mass flux.");
        require_near(result.flux.momentum_x / scale, 1.09, tolerance, "Scaled uniform momentum flux.");
    }
    const double smallest_pressure{std::numeric_limits<double>::denorm_min()};
    const auto smallest{evaluate({1.0, 0.0, 0.0, smallest_pressure}, {1.0, 0.0, 0.0, smallest_pressure})};
    require(smallest.flux.momentum_x == smallest_pressure, "Representable pressure endpoint must not be rejected.");
}

void test_solver_uniform_preservation_and_mixed_group()
{
    const auto built{cfd::verification::build_boundary_wave_mesh(8, 4, 1.0, 0.1)};
    const cfd::PrimitiveState uniform{1.0, 0.3, 0.0, 1.0};
    cfd::ConservativeField initial{built.mesh.cell_count(), gas.to_conservative(uniform)};
    const cfd::EulerBoundaryConditions conditions{4,
                                                  {cfd::EulerBoundaryConditionType::Axis,
                                                   cfd::EulerBoundaryConditionType::SlipWall, BC::farfield(uniform),
                                                   BC::farfield(uniform)}};
    cfd::AxisymmetricEulerSolver solver{built.mesh, initial, gas, conditions};
    solver.advance_to(0.1);
    for (cfd::Index cell = 0; cell < initial.size(); ++cell)
    {
        const auto p{gas.evaluate(solver.state().state(cell)).primitive};
        require_near(p.density, 1.0, tolerance, "Uniform density drift.");
        require_near(p.velocity_x, 0.3, tolerance, "Uniform axial drift.");
        require_near(p.velocity_r, 0.0, tolerance, "Axisymmetric source imbalance.");
        require_near(p.pressure, 1.0, tolerance, "Uniform pressure drift.");
    }
    // Four actual faces in one Mesh group, with different adjacent cell states.
    // Choose regimes away from sonic endpoints so their count is not sensitive
    // to the last bit of a sampled sonic fan's recovered Mach number.
    const auto condition{BC::farfield({1.0, -2.0, 0.0, 1.0})};
    constexpr std::array speeds{-8.0, 1.0, 3.0, 8.0};
    std::array<bool, 4> seen{};
    cfd::Index count{};
    for (cfd::Index face = 0; face < built.mesh.face_count(); ++face)
    {
        const auto adjacency{built.mesh.face_adjacencies()[face]};
        if (!adjacency.is_boundary() || built.mesh.face_boundary_ids()[face] != 3)
            continue;
        initial.set_state(adjacency.owner, gas.to_conservative({1.0, speeds.at(count), 0.0, 1.0}));
        ++count;
        const auto sf{built.mesh.face_area_vectors()[face]};
        const double length{built.mesh.face_lengths()[face]};
        const auto result{cfd::compute_euler_farfield_boundary_flux(initial.state(adjacency.owner),
                                                                    {sf.x / length, sf.y / length}, condition, gas)};
        seen.at(static_cast<cfd::Index>(result.regime)) = true;
    }
    require(count == 4, "Missing mixed-regime boundary faces.");
    for (const bool present : seen)
        require(present, "Mixed group did not exercise all four normal regimes.");
    cfd::AxisymmetricEulerSolver mixed{
        built.mesh,
        initial,
        gas,
        {4, {cfd::EulerBoundaryConditionType::Axis, cfd::EulerBoundaryConditionType::SlipWall, condition, condition}}};
    cfd::verification::EulerBoundaryWaveBalance balance{.initial = mixed.diagnostics().integrals};
    balance.add(mixed.step(1.0e-5));
    for (const double defect : balance.defect(mixed.diagnostics().integrals))
        require_near(defect, 0.0, tolerance, "Mixed-group conservative update.");
}
} // namespace

int main()
{
    int failures{};
    failures += cfd::test::run_test("uniform normal/oblique states", test_uniform_states_and_normal_not_total_mach);
    failures += cfd::test::run_test("supersonic causality", test_supersonic_causality_and_incoming_shock);
    failures +=
        cfd::test::run_test("invariants/entropy/backflow", test_acoustic_invariants_entropy_tangent_and_backflow);
    failures += cfd::test::run_test("shocks/sonic transitions/orientation", test_shocks_sonic_fans_and_orientation);
    failures += cfd::test::run_test("invalid/vacuum states", test_invalid_and_near_vacuum_states);
    failures += cfd::test::run_test("uniform solver/mixed group", test_solver_uniform_preservation_and_mixed_group);
    return cfd::test::finish_tests(failures, "Euler farfield");
}
