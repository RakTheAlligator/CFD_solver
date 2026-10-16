#include "support/EulerNozzle.hpp"
#include "support/TestUtils.hpp"

#include <cmath>
#include <numbers>

namespace
{
using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws_with_message;
using Reference = cfd::verification::EulerNozzleReference;
constexpr double tolerance{4.0e-12};

void test_area_mach_independent_identities()
{
    // For gamma=1.4 the area relation is exactly ((5+M^2)/6)^3/M.
    for (const double m : {0.1, 0.5, 0.9, 1.1, 2.0, 4.0})
    {
        const double ratio{std::pow((5.0 + m * m) / 6.0, 3.0) / m};
        require_near(cfd::verification::nozzle_mach_from_area(ratio, m > 1.0), m, tolerance,
                     "Incorrect area-Mach branch.");
    }
    require(cfd::verification::nozzle_mach_from_area(1.0, false) == 1.0 &&
                cfd::verification::nozzle_mach_from_area(1.0, true) == 1.0,
            "Sonic throat is not exact.");
    require_throws_with_message<std::invalid_argument>(
        [] { (void)cfd::verification::nozzle_mach_from_area(0.99, false); }, "A/Astar", "Subcritical area accepted.");
    for (const bool choked : {false, true})
    {
        const Reference reference{choked};
        double previous{};
        for (int i = 0; i <= 10; ++i)
        {
            const double x{static_cast<double>(i) / 10.0};
            const auto primitive{reference.primitive(x)};
            const double flow{primitive.density * primitive.velocity_x * reference.area(x)};
            require_near(flow, reference.mass_flow(), tolerance, "Independent section continuity is not constant.");
            const double temperature{primitive.pressure / primitive.density};
            require_near(temperature + primitive.velocity_x * primitive.velocity_x / 7.0, 1.0, tolerance,
                         "Reference total temperature varies.");
            require_near(primitive.pressure / std::pow(primitive.density, 1.4), 1.0, tolerance,
                         "Reference entropy varies.");
            if (!choked)
                require(reference.mach(x) <= 0.5 + tolerance, "Subsonic reference choked unexpectedly.");
            else if (i > 0)
                require(reference.mach(x) > previous, "Choked reference did not accelerate monotonically.");
            previous = reference.mach(x);
        }
        if (choked)
        {
            require_near(reference.mass_flow(), std::numbers::pi * 0.01 * std::sqrt(1.4) * std::pow(5.0 / 6.0, 3.0),
                         tolerance, "Independent critical mass-flow formula mismatch.");
            require(reference.mach(0.5) == 1.0 && reference.mach(1.0) > 1.0, "Wrong choked branches.");
        }
    }
}

void test_geometry_projection_initialization_and_bc()
{
    const cfd::CaloricallyPerfectGas gas{1.4, 1.0};
    double previous_geometry_error{1.0};
    for (const cfd::Index nx : {24U, 48U, 96U})
    {
        const auto built{cfd::build_mesh(cfd::verification::make_euler_nozzle(nx, nx / 8))};
        const auto &mesh{built.mesh};
        const cfd::AxisymmetricGeometry geometry{mesh};
        long double volume{};
        for (const double cell_volume : geometry.cell_volumes())
            volume += cell_volume;
        const double geometry_error{std::abs(static_cast<double>(volume) - 1.5 * std::numbers::pi * 0.01)};
        require(geometry_error < previous_geometry_error, "Polygonal geometry did not approach the analytic nozzle.");
        previous_geometry_error = geometry_error;
        for (const bool choked : {false, true})
        {
            const Reference reference{choked};
            const auto initial{cfd::verification::nozzle_initial_field(nx, nx / 8, reference, false)};
            for (cfd::Index cell = 0; cell < initial.size(); ++cell)
            {
                const auto state{gas.evaluate(initial.state(cell))};
                require(state.primitive.density > 0.0 && state.primitive.pressure > 0.0,
                        "Non-admissible projected initial data.");
            }
            cfd::AxisymmetricEulerSolver solver{mesh, initial, gas, reference.conditions()};
            require(solver.evaluate_steady_residual().mass_flow_in > 0.0, "Nozzle inlet orientation is wrong.");
            solver.step(1.0e-5);
            require(std::abs(solver.diagnostics().balance_error.mass) < tolerance &&
                        std::abs(solver.diagnostics().balance_error.energy) < tolerance,
                    "Nozzle axisymmetric update is not conservative.");
        }
    }
    const auto uniform{cfd::verification::nozzle_initial_field(24, 3, {}, true)};
    require(uniform.density()[0] == uniform.density()[23], "Uniform initialization is not uniform.");
    require_throws_with_message<std::invalid_argument>(
        [] { (void)cfd::verification::nozzle_initial_field(24, 3, {true}, true); }, "startup",
        "Incompatible uniform choked startup accepted.");
}

void test_projection_against_independent_midpoint_quadrature()
{
    const Reference reference;
    constexpr cfd::Index nx{24}, samples{4096}, column{6};
    const auto projected{cfd::verification::nozzle_initial_field(nx, 3, reference, false)};
    const double xl{static_cast<double>(column) / static_cast<double>(nx)},
        xr{static_cast<double>(column + 1) / static_cast<double>(nx)};
    const double rl{reference.radius(xl)}, rr{reference.radius(xr)};
    long double weight_sum{}, density{}, momentum{}, energy{};
    for (cfd::Index sample = 0; sample < samples; ++sample)
    {
        const double fraction{(static_cast<double>(sample) + 0.5) / static_cast<double>(samples)};
        const double radius{rl + fraction * (rr - rl)};
        const auto primitive{reference.primitive(xl + fraction * (xr - xl))};
        const long double weight{std::numbers::pi * radius * radius};
        weight_sum += weight;
        density += weight * primitive.density;
        momentum += weight * primitive.density * primitive.velocity_x;
        energy += weight *
                  (primitive.pressure / 0.4 + 0.5L * primitive.density * primitive.velocity_x * primitive.velocity_x);
    }
    const auto actual{projected.state(column)};
    require_near(actual.density, static_cast<double>(density / weight_sum), 2.0e-9,
                 "Projected density disagrees with independent quadrature.");
    require_near(actual.momentum_x, static_cast<double>(momentum / weight_sum), 2.0e-9,
                 "Projected momentum disagrees with independent quadrature.");
    require_near(actual.total_energy_density, static_cast<double>(energy / weight_sum), 2.0e-9,
                 "Projected energy disagrees with independent quadrature.");
    require(std::abs(actual.density - reference.primitive(0.5 * (xl + xr)).density) > 1.0e-7,
            "Conservative projection was replaced by centre sampling.");
}
} // namespace

int main()
{
    int failures{};
    failures +=
        cfd::test::run_test("independent area-Mach critical flow and branches", test_area_mach_independent_identities);
    failures += cfd::test::run_test("nozzle geometry initial projection and boundaries",
                                    test_geometry_projection_initialization_and_bc);
    failures += cfd::test::run_test("conservative projection independent quadrature",
                                    test_projection_against_independent_midpoint_quadrature);
    return cfd::test::finish_tests(failures, "Euler nozzle");
}
