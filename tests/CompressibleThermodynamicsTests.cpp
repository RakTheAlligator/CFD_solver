#include "cfd/compressible/CaloricallyPerfectGas.hpp"
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

void test_analytical_rest_and_moving_states()
{
    const cfd::CaloricallyPerfectGas gas{1.4, 2.0};
    const auto moving{gas.to_conservative({2.0, 3.0, 4.0, 10.0})};
    require_near(moving.density, 2.0, tolerance, "Density changed.");
    require_near(moving.momentum_x, 6.0, tolerance, "Incorrect rho*ux.");
    require_near(moving.momentum_r, 8.0, tolerance, "Incorrect rho*ur.");
    require_near(moving.total_energy_density, 50.0, tolerance, "rhoE must be volumetric, not specific energy.");
    // Independent conservative data: 25 J/m^3 internal + 25 J/m^3 kinetic.
    const auto evaluated{gas.evaluate({2.0, 6.0, 8.0, 50.0})};
    require_near(evaluated.primitive.pressure, 10.0, tolerance, "Incorrect pressure.");
    require_near(evaluated.primitive.velocity_x, 3.0, tolerance, "Incorrect ux.");
    require_near(evaluated.primitive.velocity_r, 4.0, tolerance, "Incorrect ur.");
    require_near(evaluated.thermodynamics.specific_internal_energy, 12.5, tolerance, "Incorrect e.");
    require_near(50.0 / 2.0 - evaluated.thermodynamics.specific_internal_energy, 12.5, tolerance,
                 "Incorrect kinetic energy.");
    require_near(evaluated.thermodynamics.temperature, 2.5, tolerance, "Incorrect temperature.");
    require_near(evaluated.thermodynamics.sound_speed, std::sqrt(7.0), tolerance, "Incorrect sound speed.");
    const auto rest{gas.evaluate({2.0, 0.0, 0.0, 25.0})};
    require(rest.primitive.velocity_x == 0.0 && rest.primitive.velocity_r == 0.0, "Rest velocity must be zero.");
    require_near(rest.primitive.pressure, 10.0, tolerance, "Rest pressure is incorrect.");
}

void test_round_trips_and_configurable_gas()
{
    for (const cfd::PrimitiveState input :
         {cfd::PrimitiveState{1.2, -25.0, 12.0, 101325.0}, cfd::PrimitiveState{0.02, 0.0, -0.0, 350.0},
          cfd::PrimitiveState{4.0, -2.0, 1.0, 12.0}})
    {
        const cfd::CaloricallyPerfectGas gas{2.0, 3.0};
        const auto conservative{gas.to_conservative(input)};
        const auto output{gas.evaluate(conservative).primitive};
        require_near(output.density / input.density, 1.0, tolerance, "Round-trip rho mismatch.");
        require_near(output.velocity_x, input.velocity_x, tolerance, "Round-trip ux mismatch.");
        require_near(output.velocity_r, input.velocity_r, tolerance, "Round-trip ur mismatch.");
        require_near(output.pressure / input.pressure, 1.0, tolerance, "Round-trip p mismatch.");
    }
    const cfd::CaloricallyPerfectGas gas{2.0, 3.0};
    require(gas.gamma() == 2.0 && gas.gas_constant() == 3.0, "Gas parameters changed.");
    const auto state{gas.evaluate({4.0, -8.0, 4.0, 22.0})};
    require_near(state.primitive.pressure, 12.0, tolerance, "Configured gamma was ignored.");
    require_near(state.thermodynamics.specific_internal_energy, 3.0, tolerance, "Configured e mismatch.");
    require_near(state.thermodynamics.temperature, 1.0, tolerance, "Configured R was ignored.");
    require_near(state.thermodynamics.sound_speed, std::sqrt(6.0), tolerance, "Configured a mismatch.");
}

void test_physical_scales_without_squared_momentum_underflow_or_overflow()
{
    const cfd::CaloricallyPerfectGas gas{2.0, 2.0};
    for (const double scale : {1.0e-200, 1.0e-100, 1.0, 1.0e100, 1.0e200})
    {
        const cfd::ConservativeState input{2.0 * scale, 6.0 * scale, 8.0 * scale, 35.0 * scale};
        const auto output{gas.evaluate(input)};
        require_near(output.primitive.pressure / scale, 10.0, tolerance, "Scaled pressure mismatch.");
        require_near(output.primitive.velocity_x, 3.0, tolerance, "Scaled ux mismatch.");
        require_near(output.primitive.velocity_r, 4.0, tolerance, "Scaled ur mismatch.");
        require_near(output.thermodynamics.specific_internal_energy, 5.0, tolerance, "Scaled e mismatch.");
        require_near(output.thermodynamics.temperature, 2.5, tolerance, "Scaled T mismatch.");
        require_near(output.thermodynamics.sound_speed, std::sqrt(10.0), tolerance, "Scaled a mismatch.");
        const auto rebuilt{gas.to_conservative(output.primitive)};
        require_near(rebuilt.total_energy_density / scale, 35.0, tolerance, "Scaled rhoE mismatch.");
    }
    // No hidden density/pressure floor: an exactly representable subnormal
    // rest state is admissible when its derived properties remain representable.
    const cfd::CaloricallyPerfectGas subnormal_gas{2.0, 1.0};
    const double smallest{std::numeric_limits<double>::denorm_min()};
    const auto subnormal{subnormal_gas.evaluate({smallest, 0.0, 0.0, smallest})};
    require(subnormal.primitive.density == smallest && subnormal.primitive.pressure == smallest,
            "Representable subnormal density/pressure must not be clipped.");
    require_near(subnormal.thermodynamics.specific_internal_energy, 1.0, tolerance, "Subnormal-state e mismatch.");
    require_near(subnormal.thermodynamics.temperature, 1.0, tolerance, "Subnormal-state T mismatch.");
}

void test_invalid_parameters()
{
    for (const double invalid : {0.0, -1.0, 0.5, 1.0, std::numeric_limits<double>::quiet_NaN(),
                                 std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()})
        require_throws<std::invalid_argument>([&] { const cfd::CaloricallyPerfectGas gas{invalid, 287.0}; },
                                              "Invalid gamma accepted.");
    for (const double invalid : {0.0, -1.0, std::numeric_limits<double>::quiet_NaN(),
                                 std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()})
        require_throws<std::invalid_argument>([&] { const cfd::CaloricallyPerfectGas gas{1.4, invalid}; },
                                              "Invalid R accepted.");
    const cfd::CaloricallyPerfectGas near_one{std::nextafter(1.0, 2.0), 1.0};
    require(near_one.properties({1.0, 0.0, 0.0, 1.0}).specific_internal_energy > 0.0,
            "Representable gamma immediately above one is valid.");
}

void test_invalid_primitives()
{
    const cfd::CaloricallyPerfectGas gas{1.4, 287.0};
    for (const double invalid : {0.0, -1.0, std::numeric_limits<double>::quiet_NaN(),
                                 std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()})
    {
        require_throws<std::runtime_error>([&] { static_cast<void>(gas.to_conservative({invalid, 0.0, 0.0, 1.0})); },
                                           "Invalid primitive rho accepted.");
        require_throws<std::runtime_error>([&] { static_cast<void>(gas.to_conservative({1.0, 0.0, 0.0, invalid})); },
                                           "Invalid primitive pressure accepted.");
    }
    for (const double invalid : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
                                 -std::numeric_limits<double>::infinity()})
    {
        require_throws<std::runtime_error>([&] { static_cast<void>(gas.to_conservative({1.0, invalid, 0.0, 1.0})); },
                                           "Non-finite ux accepted.");
        require_throws<std::runtime_error>([&] { static_cast<void>(gas.to_conservative({1.0, 0.0, invalid, 1.0})); },
                                           "Non-finite ur accepted.");
    }
}

void test_invalid_conservatives()
{
    const cfd::CaloricallyPerfectGas gas{1.4, 287.0};
    for (const double invalid : {0.0, -1.0, std::numeric_limits<double>::quiet_NaN(),
                                 std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()})
    {
        require_throws<std::runtime_error>([&] { static_cast<void>(gas.evaluate({invalid, 0.0, 0.0, 1.0})); },
                                           "Invalid conservative rho accepted.");
        require_throws<std::runtime_error>([&] { static_cast<void>(gas.evaluate({1.0, 0.0, 0.0, invalid})); },
                                           "Invalid rhoE accepted.");
    }
    for (const double invalid : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
                                 -std::numeric_limits<double>::infinity()})
    {
        require_throws<std::runtime_error>([&] { static_cast<void>(gas.evaluate({1.0, invalid, 0.0, 1.0})); },
                                           "Non-finite x momentum accepted.");
        require_throws<std::runtime_error>([&] { static_cast<void>(gas.evaluate({1.0, 0.0, invalid, 1.0})); },
                                           "Non-finite r momentum accepted.");
    }
    for (const double energy : {1.0, 2.0})
        require_throws<std::runtime_error>([&] { static_cast<void>(gas.evaluate({1.0, 2.0, 0.0, energy})); },
                                           "Non-positive internal energy accepted.");
}

void test_nonrepresentable_properties_and_high_mach_cancellation()
{
    const cfd::CaloricallyPerfectGas gas{2.0, 1.0};
    require_throws<std::runtime_error>([&] { static_cast<void>(gas.properties({1.0e308, 0.0, 0.0, 1.0e-308})); },
                                       "Underflowing e accepted.");
    require_throws<std::runtime_error>([&] { static_cast<void>(gas.properties({1.0e-308, 0.0, 0.0, 1.0e308})); },
                                       "Overflowing e accepted.");
    const cfd::CaloricallyPerfectGas huge_R{2.0, 1.0e308};
    require_throws<std::runtime_error>([&] { static_cast<void>(huge_R.properties({1.0e308, 0.0, 0.0, 1.0})); },
                                       "Underflowing T accepted.");
    const cfd::CaloricallyPerfectGas tiny_R{2.0, 1.0e-308};
    require_throws<std::runtime_error>([&] { static_cast<void>(tiny_R.properties({1.0, 0.0, 0.0, 1.0e308})); },
                                       "Overflowing T accepted.");
    require_throws<std::runtime_error>([&] { static_cast<void>(gas.to_conservative({1.0, 1.0e10, 0.0, 1.0})); },
                                       "Lost internal energy must not be repaired silently.");
    require_throws<std::runtime_error>([&] { static_cast<void>(gas.evaluate({1.0, 1.0e10, 0.0, 5.0e19})); },
                                       "Zero internal energy at high Mach accepted.");
    require_throws<std::runtime_error>(
        [&] { static_cast<void>(gas.to_conservative({1.0e200, 1.0e100, 0.0, 1.0e200})); },
        "Overflowing conservative energy accepted.");
}

void test_inputs_unchanged()
{
    const cfd::CaloricallyPerfectGas gas{2.0, 3.0};
    const cfd::PrimitiveState primitive{4.0, -2.0, 1.0, 12.0};
    const cfd::ConservativeState conservative{4.0, -8.0, 4.0, 22.0};
    static_cast<void>(gas.to_conservative(primitive));
    static_cast<void>(gas.evaluate(conservative));
    require(primitive.density == 4.0 && primitive.velocity_x == -2.0 && primitive.velocity_r == 1.0 &&
                primitive.pressure == 12.0,
            "Primitive input changed.");
    require(conservative.density == 4.0 && conservative.momentum_x == -8.0 && conservative.momentum_r == 4.0 &&
                conservative.total_energy_density == 22.0,
            "Conservative input changed.");
}
} // namespace

int main()
{
    int failures{};
    failures += cfd::test::run_test("analytical perfect-gas states", test_analytical_rest_and_moving_states);
    failures += cfd::test::run_test("perfect-gas round trips and parameters", test_round_trips_and_configurable_gas);
    failures += cfd::test::run_test("perfect-gas physical scales",
                                    test_physical_scales_without_squared_momentum_underflow_or_overflow);
    failures += cfd::test::run_test("invalid gas parameters", test_invalid_parameters);
    failures += cfd::test::run_test("invalid primitive states", test_invalid_primitives);
    failures += cfd::test::run_test("invalid conservative states", test_invalid_conservatives);
    failures += cfd::test::run_test("property range and cancellation",
                                    test_nonrepresentable_properties_and_high_mach_cancellation);
    failures += cfd::test::run_test("thermodynamic inputs unchanged", test_inputs_unchanged);
    return cfd::test::finish_tests(failures, "compressible thermodynamics");
}
