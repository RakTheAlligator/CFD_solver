#include "cfd/compressible/EulerFlux.hpp"
#include "support/TestUtils.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace
{
using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws;
constexpr double tolerance{1.0e-12};
constexpr cfd::ConservativeState left{2.0, 6.0, 8.0, 35.0};  // ux=3, ur=4, p=10 at gamma=2.
constexpr cfd::ConservativeState right{1.0, -1.0, 2.0, 6.5}; // ux=-1, ur=2, p=4.

void require_flux(const cfd::EulerFlux &actual, const cfd::EulerFlux &expected)
{
    require_near(actual.mass, expected.mass, tolerance, "Incorrect mass flux.");
    require_near(actual.momentum_x, expected.momentum_x, tolerance, "Incorrect x momentum flux.");
    require_near(actual.momentum_r, expected.momentum_r, tolerance, "Incorrect r momentum flux.");
    require_near(actual.energy, expected.energy, tolerance, "Incorrect energy flux.");
}

void test_analytical_physical_fluxes()
{
    const cfd::CaloricallyPerfectGas gas{2.0, 3.0};
    require_flux(cfd::compute_physical_euler_flux(left, {1.0, 0.0}, gas), {6.0, 28.0, 24.0, 135.0});
    require_flux(cfd::compute_physical_euler_flux(left, {0.0, 1.0}, gas), {8.0, 24.0, 42.0, 180.0});
    require_flux(cfd::compute_physical_euler_flux(left, {0.6, 0.8}, gas), {10.0, 36.0, 48.0, 225.0});
    require_flux(cfd::compute_physical_euler_flux(left, {-0.6, -0.8}, gas), {-10.0, -36.0, -48.0, -225.0});
    const cfd::CaloricallyPerfectGas other_gas{1.4, 17.0};
    require_flux(cfd::compute_physical_euler_flux(left, {1.0, 0.0}, other_gas), {6.0, 22.0, 24.0, 117.0});
    require_flux(cfd::compute_physical_euler_flux({2.0, 0.0, 0.0, 10.0}, {0.6, 0.8}, gas), {0.0, 6.0, 8.0, 0.0});
}

void test_uniform_state_consistency_and_wave_speed()
{
    const cfd::CaloricallyPerfectGas gas{2.0, 3.0};
    const auto x{cfd::compute_rusanov_flux(left, left, {1.0, 0.0}, gas)};
    require_flux(x.flux, {6.0, 28.0, 24.0, 135.0});
    require_near(x.maximum_wave_speed, 3.0 + std::sqrt(10.0), tolerance, "Incorrect x signal speed.");
    const auto r{cfd::compute_rusanov_flux(left, left, {0.0, 1.0}, gas)};
    require_flux(r.flux, {8.0, 24.0, 42.0, 180.0});
    require_near(r.maximum_wave_speed, 4.0 + std::sqrt(10.0), tolerance, "Incorrect r signal speed.");
    const auto oblique{cfd::compute_rusanov_flux(left, left, {0.6, 0.8}, gas)};
    require_flux(oblique.flux, {10.0, 36.0, 48.0, 225.0});
    require_near(oblique.maximum_wave_speed, 5.0 + std::sqrt(10.0), tolerance, "Incorrect oblique signal speed.");
    require_flux(oblique.flux, cfd::compute_physical_euler_flux(left, {0.6, 0.8}, gas));
}

void test_different_states_with_independent_reference()
{
    const cfd::CaloricallyPerfectGas gas{2.0, 3.0};
    const double speed{3.0 + std::sqrt(10.0)};
    const auto result{cfd::compute_rusanov_flux(left, right, {1.0, 0.0}, gas)};
    // FnL=[6,28,24,135], FnR=[-1,5,-2,-10.5], deltaU=[-1,-7,-6,-28.5].
    require_flux(result.flux, {2.5 + 0.5 * speed, 16.5 + 3.5 * speed, 11.0 + 3.0 * speed, 62.25 + 14.25 * speed});
    require_near(result.maximum_wave_speed, speed, tolerance, "Left signal speed must dominate.");

    const double oblique_speed{5.0 + std::sqrt(10.0)};
    const auto oblique{cfd::compute_rusanov_flux(left, right, {0.6, 0.8}, gas)};
    require_flux(oblique.flux, {5.5 + 0.5 * oblique_speed, 18.7 + 3.5 * oblique_speed, 26.6 + 3.0 * oblique_speed,
                                117.75 + 14.25 * oblique_speed});
    const cfd::ConservativeState fast_right{1.0, 20.0, 0.0, 204.0};
    const double right_speed{20.0 + std::sqrt(8.0)};
    const auto right_dominant{cfd::compute_rusanov_flux(left, fast_right, {1.0, 0.0}, gas)};
    require_near(right_dominant.maximum_wave_speed, right_speed, tolerance, "Right signal speed must dominate.");
    require_flux(right_dominant.flux, {13.0 + 0.5 * right_speed, 216.0 - 7.0 * right_speed, 12.0 + 4.0 * right_speed,
                                       2147.5 - 84.5 * right_speed});
}

void test_conservative_symmetry()
{
    const cfd::CaloricallyPerfectGas gas{2.0, 3.0};
    for (const cfd::Vector2 normal :
         {cfd::Vector2{1.0, 0.0}, cfd::Vector2{0.0, 1.0}, cfd::Vector2{0.6, 0.8}, cfd::Vector2{-0.6, 0.8}})
    {
        const auto forward{cfd::compute_rusanov_flux(left, right, normal, gas)};
        const auto reverse{cfd::compute_rusanov_flux(right, left, {-normal.x, -normal.y}, gas)};
        require_flux(reverse.flux,
                     {-forward.flux.mass, -forward.flux.momentum_x, -forward.flux.momentum_r, -forward.flux.energy});
        require_near(forward.maximum_wave_speed, reverse.maximum_wave_speed, tolerance,
                     "Signal speed must be orientation independent.");
    }
}

void test_stationary_contact_is_diffusive()
{
    const cfd::CaloricallyPerfectGas gas{2.0, 3.0};
    // Equal p=2 and zero velocity, different density. Rusanov intentionally
    // diffuses the contact: lambda=2, mass flux=-lambda*(4-1)/2=-3.
    const auto result{cfd::compute_rusanov_flux({1.0, 0.0, 0.0, 2.0}, {4.0, 0.0, 0.0, 2.0}, {1.0, 0.0}, gas)};
    require_flux(result.flux, {-3.0, 2.0, 0.0, 0.0});
    require_near(result.maximum_wave_speed, 2.0, tolerance, "Incorrect stationary contact speed.");
}

void test_invalid_normals_and_states()
{
    const cfd::CaloricallyPerfectGas gas{2.0, 3.0};
    for (const cfd::Vector2 normal :
         {cfd::Vector2{0.0, 0.0}, cfd::Vector2{2.0, 0.0}, cfd::Vector2{0.6, 0.6}, cfd::Vector2{1.0e-200, 0.0},
          cfd::Vector2{std::numeric_limits<double>::quiet_NaN(), 1.0},
          cfd::Vector2{0.0, std::numeric_limits<double>::infinity()}})
    {
        require_throws<std::invalid_argument>(
            [&] { static_cast<void>(cfd::compute_physical_euler_flux(left, normal, gas)); },
            "Invalid physical-flux normal accepted.");
        require_throws<std::invalid_argument>(
            [&] { static_cast<void>(cfd::compute_rusanov_flux(left, right, normal, gas)); },
            "Invalid Rusanov normal accepted.");
    }
    for (const cfd::ConservativeState invalid :
         {cfd::ConservativeState{}, cfd::ConservativeState{-1.0, 0.0, 0.0, 2.0},
          cfd::ConservativeState{1.0, 2.0, 0.0, 1.0},
          cfd::ConservativeState{1.0, 0.0, 0.0, std::numeric_limits<double>::quiet_NaN()},
          cfd::ConservativeState{1.0, std::numeric_limits<double>::infinity(), 0.0, 2.0}})
    {
        require_throws<std::runtime_error>(
            [&] { static_cast<void>(cfd::compute_physical_euler_flux(invalid, {1.0, 0.0}, gas)); },
            "Invalid physical state accepted.");
        require_throws<std::runtime_error>(
            [&] { static_cast<void>(cfd::compute_rusanov_flux(invalid, right, {1.0, 0.0}, gas)); },
            "Invalid left state accepted.");
        require_throws<std::runtime_error>(
            [&] { static_cast<void>(cfd::compute_rusanov_flux(left, invalid, {1.0, 0.0}, gas)); },
            "Invalid right state accepted.");
    }
}

void test_flux_ranges_and_avoidable_intermediate_overflow()
{
    const cfd::CaloricallyPerfectGas gas{2.0, 3.0};
    const cfd::ConservativeState huge{1.0e308, 1.0e208, 0.0, 1.0e308};
    const auto flux{cfd::compute_rusanov_flux(huge, huge, {1.0, 0.0}, gas).flux};
    // rhoE+p would overflow double, but its product with un=1e-100 fits.
    require_near(flux.energy / 1.0e208, 2.0, tolerance, "Avoidable enthalpy-sum overflow.");
    require_near(flux.mass / 1.0e208, 1.0, tolerance, "Large mass flux mismatch.");
    const double scale{1.0e-200};
    const cfd::ConservativeState tiny{2.0 * scale, 6.0 * scale, 8.0 * scale, 35.0 * scale};
    const auto tiny_flux{cfd::compute_rusanov_flux(tiny, tiny, {1.0, 0.0}, gas).flux};
    require_flux(
        {tiny_flux.mass / scale, tiny_flux.momentum_x / scale, tiny_flux.momentum_r / scale, tiny_flux.energy / scale},
        {6.0, 28.0, 24.0, 135.0});
    const cfd::ConservativeState overflow{1.0, 1.0e110, 0.0, 1.0e220};
    require_throws<std::runtime_error>(
        [&] { static_cast<void>(cfd::compute_physical_euler_flux(overflow, {1.0, 0.0}, gas)); },
        "Unrepresentable physical flux accepted.");
    require_throws<std::runtime_error>(
        [&] { static_cast<void>(cfd::compute_rusanov_flux(overflow, overflow, {1.0, 0.0}, gas)); },
        "Unrepresentable numerical flux accepted.");
}

void test_inputs_unchanged()
{
    const cfd::CaloricallyPerfectGas gas{2.0, 3.0};
    cfd::ConservativeState l{left}, r{right};
    cfd::Vector2 normal{0.6, 0.8};
    static_cast<void>(cfd::compute_physical_euler_flux(l, normal, gas));
    static_cast<void>(cfd::compute_rusanov_flux(l, r, normal, gas));
    require(l.density == left.density && l.momentum_x == left.momentum_x && l.momentum_r == left.momentum_r &&
                l.total_energy_density == left.total_energy_density,
            "Left state changed.");
    require(r.density == right.density && r.momentum_x == right.momentum_x && r.momentum_r == right.momentum_r &&
                r.total_energy_density == right.total_energy_density,
            "Right state changed.");
    require(normal.x == 0.6 && normal.y == 0.8, "Normal changed.");
}
} // namespace

int main()
{
    int failures{};
    failures += cfd::test::run_test("analytical physical Euler fluxes", test_analytical_physical_fluxes);
    failures += cfd::test::run_test("uniform Rusanov flux and speed", test_uniform_state_consistency_and_wave_speed);
    failures += cfd::test::run_test("unequal Rusanov states", test_different_states_with_independent_reference);
    failures += cfd::test::run_test("Rusanov conservative symmetry", test_conservative_symmetry);
    failures += cfd::test::run_test("diffusive stationary contact", test_stationary_contact_is_diffusive);
    failures += cfd::test::run_test("invalid Euler normals and states", test_invalid_normals_and_states);
    failures +=
        cfd::test::run_test("Euler flux numerical ranges", test_flux_ranges_and_avoidable_intermediate_overflow);
    failures += cfd::test::run_test("Euler inputs unchanged", test_inputs_unchanged);
    return cfd::test::finish_tests(failures, "Euler flux");
}
