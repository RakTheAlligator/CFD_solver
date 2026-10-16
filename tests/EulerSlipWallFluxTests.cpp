#include "cfd/compressible/EulerSlipWallFlux.hpp"
#include "support/TestUtils.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace
{
using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws_with_message;
constexpr double tolerance{5.0e-13};
const cfd::CaloricallyPerfectGas gas{1.4, 1.0};

void require_wall_flux(const cfd::EulerSlipWallFluxResult &result, const cfd::Vector2 normal)
{
    require(result.flux.mass == 0.0 && result.flux.energy == 0.0, "Wall mass and energy flux must be exactly zero.");
    require(std::isfinite(result.wall_pressure) && result.wall_pressure >= 0.0, "Invalid wall pressure.");
    require(std::isfinite(result.maximum_wave_speed) && result.maximum_wave_speed > 0.0, "Invalid wall wave bound.");
    require_near(result.flux.momentum_x, result.wall_pressure * normal.x, tolerance, "Axial wall force is not normal.");
    require_near(result.flux.momentum_r, result.wall_pressure * normal.y, tolerance,
                 "Radial wall force is not normal.");
    const double scale{result.wall_pressure > 0.0 ? result.wall_pressure : 1.0};
    require_near((result.flux.momentum_x / scale) * normal.y - (result.flux.momentum_r / scale) * normal.x, 0.0,
                 tolerance, "Wall force has a tangential component.");
}

void test_rest_and_tangential_flow()
{
    for (const double tangent : {0.0, 0.7, -3.0})
    {
        const auto state{gas.to_conservative({2.0, tangent, 0.0, 3.0})};
        const auto result{cfd::compute_euler_slip_wall_flux(state, {0.0, 1.0}, gas)};
        require_wall_flux(result, {0.0, 1.0});
        require_near(result.wall_pressure, 3.0, tolerance, "Tangential motion changed wall pressure.");
        require_near(result.maximum_wave_speed, std::sqrt(1.4 * 3.0 / 2.0), tolerance,
                     "Tangential speed must not enter the normal CFL bound.");
    }
}

void test_independent_negative_rusanov_counterexample()
{
    const double sound{std::sqrt(1.4)}, un{-0.8 * sound};
    const auto inside{gas.to_conservative({1.0, un, 0.0, 1.0})};
    const auto mirror{gas.to_conservative({1.0, -un, 0.0, 1.0})};
    const auto old_flux{cfd::compute_rusanov_flux(inside, mirror, {1.0, 0.0}, gas)};
    require_near(old_flux.flux.momentum_x, -0.12, tolerance, "Counterexample no longer reproduces negative pressure.");
    const auto result{cfd::compute_euler_slip_wall_flux(inside, {1.0, 0.0}, gas)};
    require_wall_flux(result, {1.0, 0.0});
    require_near(result.wall_pressure, 0.29509034655744, tolerance, "Incorrect isentropic receding-wall pressure.");
    require(result.wall_pressure > 0.0 && result.wall_pressure < 1.0,
            "Expansion must reduce pressure without tension.");
}

double independent_shock_pressure(const double rho, const double pressure, const double un, const double gamma)
{
    // Independently solve the monotone shock wave curve f(p_star)=un by
    // bisection, rather than repeating the production closed-form quadratic.
    const auto velocity_drop = [=](const double star) {
        return (star - pressure) *
               std::sqrt(2.0 / ((gamma + 1.0) * rho) / (star + (gamma - 1.0) / (gamma + 1.0) * pressure));
    };
    double lower{pressure}, upper{2.0 * pressure};
    while (velocity_drop(upper) < un)
        upper *= 2.0;
    for (int iteration = 0; iteration < 100; ++iteration)
    {
        const double middle{lower + 0.5 * (upper - lower)};
        if (velocity_drop(middle) > un)
            upper = middle;
        else
            lower = middle;
    }
    return lower + 0.5 * (upper - lower);
}

void test_compression_rankine_hugoniot_and_cfl()
{
    for (const double gamma : {1.01, 1.4, 5.0 / 3.0, 3.0, 10.0})
    {
        const cfd::CaloricallyPerfectGas model{gamma, 1.0};
        constexpr double rho{2.0}, pressure{3.0};
        const double sound{std::sqrt(gamma * pressure / rho)};
        for (const double mach : {1.0e-4, 0.5, 2.0, 20.0})
        {
            const double un{mach * sound};
            const auto result{
                cfd::compute_euler_slip_wall_flux(model.to_conservative({rho, un, 0.0, pressure}), {1.0, 0.0}, model)};
            require_wall_flux(result, {1.0, 0.0});
            const double star{independent_shock_pressure(rho, pressure, un, gamma)};
            require_near(result.wall_pressure, star, tolerance, "Wall pressure violates the shock wave curve.");
            require(result.wall_pressure > pressure, "Compression must increase pressure.");
            const double ratio{star / pressure}, beta{(gamma - 1.0) / (gamma + 1.0)};
            const double star_density{rho * (ratio + beta) / (beta * ratio + 1.0)};
            const double star_sound{std::sqrt(gamma * star / star_density)};
            const double shock_speed{std::abs(
                un - sound * std::sqrt((gamma + 1.0) / (2.0 * gamma) * ratio + (gamma - 1.0) / (2.0 * gamma)))};
            const double required{std::max({un + sound, star_sound, shock_speed})};
            require(result.maximum_wave_speed >= required * (1.0 - tolerance),
                    "CFL bound misses a normal Riemann speed.");
        }
    }
    // An exactly representable RH example: gamma=2, rho=1, p=1/2, un=1.
    const cfd::CaloricallyPerfectGas model{2.0, 1.0};
    const auto result{cfd::compute_euler_slip_wall_flux({1.0, 1.0, 0.0, 1.0}, {1.0, 0.0}, model)};
    require(result.wall_pressure == 2.5, "Exact compression pressure should be 5/2.");
    require(result.maximum_wave_speed >= std::sqrt(5.0), "Compression needs a stronger CFL bound than |un|+a alone.");
}

void test_oblique_wall_and_tangential_invariance()
{
    constexpr cfd::Vector2 normal{0.6, 0.8};
    for (const double un : {-0.5, 0.0, 1.0})
    {
        const auto aligned{
            cfd::compute_euler_slip_wall_flux(gas.to_conservative({1.0, un, 0.0, 1.0}), {1.0, 0.0}, gas)};
        for (const double tangent : {0.0, -2.0, 7.5})
        {
            const auto state{gas.to_conservative(
                {1.0, un * normal.x - tangent * normal.y, un * normal.y + tangent * normal.x, 1.0})};
            const auto result{cfd::compute_euler_slip_wall_flux(state, normal, gas)};
            require_wall_flux(result, normal);
            require_near(result.wall_pressure, aligned.wall_pressure, tolerance,
                         "Oblique/tangential pressure dependence.");
            require_near(result.maximum_wave_speed, aligned.maximum_wave_speed, tolerance,
                         "Oblique/tangential CFL dependence.");
        }
    }
}

void test_vacuum_limit_and_non_vacuum_underflow()
{
    // gamma=2, rho=1, p=1/2 => a=1 and exact vacuum threshold un=-2.
    const cfd::CaloricallyPerfectGas model{2.0, 1.0};
    for (const double un : {-1.0, std::nextafter(-2.0, 0.0), -2.0, std::nextafter(-2.0, -3.0), -3.0})
    {
        const auto result{
            cfd::compute_euler_slip_wall_flux(model.to_conservative({1.0, un, 0.0, 0.5}), {1.0, 0.0}, model)};
        require_wall_flux(result, {1.0, 0.0});
        if (un > -2.0)
        {
            const long double theta{1.0L + static_cast<long double>(un) / 2.0L};
            const double expected{static_cast<double>(0.5L * theta * theta * theta * theta)};
            require(result.wall_pressure > 0.0, "Near-vacuum positive pressure must not be clipped.");
            require_near(result.wall_pressure / expected, 1.0, tolerance,
                         "Near-vacuum pressure lost relative accuracy.");
        }
        else
            require(result.wall_pressure == 0.0 && result.flux.momentum_x == 0.0 && result.flux.momentum_r == 0.0,
                    "Physical wall vacuum must produce exactly zero force.");
        require_near(result.maximum_wave_speed, std::abs(un) + 1.0, tolerance, "Vacuum fan CFL bound is inconsistent.");
    }
    const cfd::CaloricallyPerfectGas nearly_isothermal{1.01, 1.0};
    const double un{-0.99 * 2.0 * std::sqrt(1.01) / (1.01 - 1.0)};
    require_throws_with_message<std::runtime_error>(
        [&] {
            (void)cfd::compute_euler_slip_wall_flux(nearly_isothermal.to_conservative({1.0, un, 0.0, 1.0}), {1.0, 0.0},
                                                    nearly_isothermal);
        },
        "underflows", "Non-vacuum unrepresentable pressure must fail, not be replaced by zero.");
}

void test_invalid_and_extreme_states()
{
    const auto valid{gas.to_conservative({1.0, 0.0, 0.0, 1.0})};
    for (const auto normal :
         {cfd::Vector2{}, cfd::Vector2{1.1, 0.0}, cfd::Vector2{std::numeric_limits<double>::quiet_NaN(), 0.0},
          cfd::Vector2{0.0, std::numeric_limits<double>::infinity()}})
        require_throws_with_message<std::invalid_argument>(
            [&] { (void)cfd::compute_euler_slip_wall_flux(valid, normal, gas); }, "unit normal",
            "Invalid wall normal accepted.");
    for (const auto state : {cfd::ConservativeState{}, cfd::ConservativeState{-1.0, 0.0, 0.0, 1.0},
                             cfd::ConservativeState{1.0, 2.0, 0.0, 1.0},
                             cfd::ConservativeState{1.0, std::numeric_limits<double>::quiet_NaN(), 0.0, 1.0},
                             cfd::ConservativeState{1.0, 0.0, 0.0, std::numeric_limits<double>::infinity()}})
        require_throws_with_message<std::runtime_error>(
            [&] { (void)cfd::compute_euler_slip_wall_flux(state, {1.0, 0.0}, gas); }, "Perfect gas",
            "Invalid gas state accepted by wall flux.");
    const auto large{gas.to_conservative({1.0, 1.5e154, 0.0, 1.0e306})};
    require_throws_with_message<std::runtime_error>(
        [&] { (void)cfd::compute_euler_slip_wall_flux(large, {1.0, 0.0}, gas); }, "wall pressure",
        "Unrepresentable shock pressure must fail explicitly.");
    // Very small/large density and pressure at the same finite sound speed.
    for (const double scale : {1.0e-200, 1.0e200})
    {
        const auto state{gas.to_conservative({scale, 0.3, 0.0, scale})};
        const auto result{cfd::compute_euler_slip_wall_flux(state, {1.0, 0.0}, gas)};
        const auto unit{cfd::compute_euler_slip_wall_flux(gas.to_conservative({1.0, 0.3, 0.0, 1.0}), {1.0, 0.0}, gas)};
        require(result.wall_pressure > 0.0, "Scaled pressure was lost.");
        require_near(result.wall_pressure / scale, unit.wall_pressure, tolerance,
                     "Wall pressure scaling is inconsistent.");
        require_near(result.maximum_wave_speed, unit.maximum_wave_speed, tolerance,
                     "Wave bound depends on density units.");
    }
}

} // namespace

int main()
{
    int failures{};
    failures += cfd::test::run_test("wall rest/tangential flow", test_rest_and_tangential_flow);
    failures +=
        cfd::test::run_test("negative Rusanov wall counterexample", test_independent_negative_rusanov_counterexample);
    failures += cfd::test::run_test("wall shock pressure and CFL", test_compression_rankine_hugoniot_and_cfl);
    failures += cfd::test::run_test("oblique wall invariance", test_oblique_wall_and_tangential_invariance);
    failures += cfd::test::run_test("wall expansion, vacuum and underflow", test_vacuum_limit_and_non_vacuum_underflow);
    failures += cfd::test::run_test("wall invalid/extreme states", test_invalid_and_extreme_states);
    return cfd::test::finish_tests(failures, "Euler slip wall flux");
}
