#include "cfd/compressible/CaloricallyPerfectGas.hpp"
#include "support/RadialExpansionExact.hpp"
#include "support/RadialExpansionStudy.hpp"
#include "support/TestUtils.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <numbers>

namespace
{
using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws;

void test_exact_moments_and_volume_projection()
{
    const cfd::verification::RadialExpansionExact exact;
    const cfd::CaloricallyPerfectGas gas{1.4, 1.0};
    // Gauss integration of r*U(r), degree <=3, independent of the moment formulas.
    const std::array abscissae{-std::sqrt(3.0 / 5.0), 0.0, std::sqrt(3.0 / 5.0)};
    constexpr std::array weights{5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
    for (const auto bounds : {std::array{0.0, 0.2}, std::array{0.2, 0.7}, std::array{1.0, 1.000001}})
    {
        const double inner{bounds[0]}, outer{bounds[1]}, middle{inner + (outer - inner) / 2.0};
        for (const double time : {0.0, 0.2})
        {
            long double weight{}, first{}, second{}, radial_momentum{}, energy{};
            for (cfd::Index point = 0; point < weights.size(); ++point)
            {
                const double r{middle + (outer - inner) / 2.0 * abscissae.at(point)};
                const long double w{weights.at(point) * static_cast<long double>(r)};
                const auto state{exact.point_state(r, time)};
                weight += w;
                first += w * r;
                second += w * r * r;
                radial_momentum += w * state.momentum_r;
                energy += w * state.total_energy_density;
            }
            const auto moments{exact.moments(inner, outer)};
            const auto average{exact.cell_average(inner, outer, time)};
            require_near(moments.mean_radius, static_cast<double>(first / weight), 2.0e-15, "Volume mean radius.");
            require_near(moments.mean_squared_radius, static_cast<double>(second / weight), 2.0e-15,
                         "Volume mean radius squared.");
            require_near(average.momentum_r, static_cast<double>(radial_momentum / weight), 2.0e-15,
                         "Projected momentum.");
            require_near(average.total_energy_density, static_cast<double>(energy / weight), 2.0e-15,
                         "Projected energy.");
            require_near(gas.evaluate(average).primitive.pressure, exact.projected_pressure(inner, outer, time),
                         2.0e-15, "EOS of projected state differs from variance formula.");
        }
    }
    const auto axis{exact.moments(0.0, 1.0)};
    require_near(axis.mean_radius, 2.0 / 3.0, 1.0e-15, "Axis volume centroid is not the meridional midpoint.");
    require(axis.mean_squared_radius == 0.5, "Axis quadratic moment.");
    require_near(axis.radius_variance, 1.0 / 18.0, 1.0e-15, "Axis radius variance.");
    const double outer{std::nextafter(1.0, 2.0)};
    const auto thin{exact.moments(1.0, outer)};
    const double h{outer - 1.0};
    require(thin.radius_variance > 0.0, "Thin annulus variance was lost through cancellation.");
    require_near(thin.radius_variance / (h * h / 12.0), 1.0, 1.0e-15, "Thin annulus variance is unstable.");
}

void test_exact_dynamics_and_invalid_parameters()
{
    const cfd::verification::RadialExpansionExact exact;
    const auto initial{exact.values(0.0)}, later{exact.values(1.0)};
    require(later.density == 0.25 && later.rate == 0.5, "Expansion time scaling.");
    require_near(later.pressure, std::pow(2.0, -2.8), 1.0e-15, "Isentropic pressure scaling.");
    require(initial.density > later.density && initial.pressure > later.pressure, "Expansion must change rho and p.");
    const cfd::verification::RadialExpansionExact rest{{.alpha0 = 0.0}};
    require(rest.cell_average(0.0, 1.0, 1.0).momentum_r == 0.0, "Rest limit is not regular.");
    for (const auto parameters :
         {cfd::verification::RadialExpansionParameters{.density = 0.0},
          cfd::verification::RadialExpansionParameters{.pressure = -1.0},
          cfd::verification::RadialExpansionParameters{.gamma = 1.0},
          cfd::verification::RadialExpansionParameters{.alpha0 = -1.0},
          cfd::verification::RadialExpansionParameters{.alpha0 = std::numeric_limits<double>::infinity()},
          cfd::verification::RadialExpansionParameters{.density = std::numeric_limits<double>::quiet_NaN()}})
        require_throws<std::invalid_argument>([&] { const cfd::verification::RadialExpansionExact bad{parameters}; },
                                              "Invalid analytical parameters accepted.");
    for (const double time : {-1.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()})
        require_throws<std::invalid_argument>([&] { (void)exact.values(time); }, "Invalid reference time accepted.");
    for (const auto bounds : {std::array{-1.0, 1.0}, std::array{1.0, 1.0}, std::array{2.0, 1.0}})
        require_throws<std::invalid_argument>([&] { (void)exact.cell_average(bounds[0], bounds[1], 0.0); },
                                              "Invalid cell interval accepted.");
    require_throws<std::invalid_argument>([&] { (void)exact.point_state(-1.0, 0.0); }, "Negative radius accepted.");
}

void test_volume_weighted_norms()
{
    cfd::verification::RadialErrorAccumulator errors;
    errors.add(1.0, {2.0, 0.0, -1.0, 3.0}, {});
    errors.add(3.0, {-5.0, 0.0, 2.0, -4.0}, {});
    const auto norms{errors.finish()};
    require_near(norms[0].l1, 4.25, 1.0e-15, "Volume weighted L1.");
    require_near(norms[0].rms, std::sqrt(79.0 / 4.0), 1.0e-15, "Volume weighted RMS.");
    require(norms[0].linf == 5.0, "Signed Linf error.");
    require_throws<std::runtime_error>([&] { errors.add(0.0, {}, {}); }, "Zero volume accepted.");
    require_throws<std::runtime_error>([&] { errors.add(1.0, {std::numeric_limits<double>::quiet_NaN()}, {}); },
                                       "NaN error accepted.");
}

void test_moving_expansion_source_geometry_and_boundary_flux()
{
    const cfd::verification::RadialGrid grid{.nx = 2, .nr = 4};
    const auto built{cfd::verification::build_radial_grid(grid)};
    const cfd::AxisymmetricGeometry geometry{built.mesh};
    const cfd::verification::RadialExpansionExact exact;
    const auto initial{cfd::verification::radial_initial_field(grid, exact)};
    const cfd::CaloricallyPerfectGas gas{1.4, 1.0};
    for (cfd::Index cell = 0; cell < initial.size(); ++cell)
    {
        require_near(geometry.cell_volumes()[cell], cfd::verification::radial_cell_volume(grid, cell), 1.0e-14,
                     "Revolved cell volume must use radial weighting.");
        require_near(geometry.cell_radial_source_coefficients()[cell], 2.0 * std::numbers::pi / 8.0, 1.0e-14,
                     "G must equal 2*pi*dx*dr.");
    }
    cfd::AxisymmetricEulerSolver solver{built.mesh, initial, gas, cfd::verification::radial_conditions()};
    const auto before{solver.diagnostics().integrals};
    const double expected_source{cfd::verification::independent_radial_source(grid, initial, 1.4)};
    const auto result{solver.step(1.0e-5)};
    // This identity uses measured integral change and an independently formed
    // G*p_EOS. Missing/reversed production source cannot pass by reporting a
    // correspondingly modified diagnostic source.
    const double recovered_source{(result.integrals.momentum_r - before.momentum_r) / result.time_step +
                                  result.boundary_outflow_rate.momentum_r};
    require_near(recovered_source, expected_source, 2.0e-11, "Moving radial momentum source sign/omission.");
    require_near(result.radial_pressure_source_rate, expected_source, 2.0e-15, "Reported radial source.");
    const auto outer{initial.state(initial.size() - 1)};
    const double velocity{outer.momentum_r / outer.density};
    const double pressure{0.4 * (outer.total_energy_density - 0.5 * outer.momentum_r * velocity)};
    const double area{2.0 * std::numbers::pi};
    require_near(result.boundary_outflow_rate.mass, area * outer.momentum_r, 2.0e-14,
                 "Transmissive boundary mass flux.");
    require_near(result.boundary_outflow_rate.energy, area * (outer.total_energy_density + pressure) * velocity,
                 2.0e-14, "Transmissive boundary energy flux.");
    require_near(result.boundary_outflow_rate.momentum_r, area * (outer.momentum_r * velocity + pressure), 2.0e-14,
                 "Transmissive boundary radial momentum flux.");
    require(result.integrals.mass < before.mass && result.integrals.energy < before.energy &&
                result.integrals.momentum_r < before.momentum_r,
            "Expansion must evolve nontrivially and discharge mass/energy.");
    require_near(result.integrals.mass - before.mass + result.time_step * result.boundary_outflow_rate.mass, 0.0,
                 2.0e-14, "Integrated mass conservation with outflow.");
    require_near(result.integrals.energy - before.energy + result.time_step * result.boundary_outflow_rate.energy, 0.0,
                 2.0e-14, "Integrated energy conservation with outflow.");
}

void require_run_valid(const cfd::verification::RadialRun &run)
{
    require(run.result.accepted_steps > 1 && run.result.rejected_steps == 0,
            "Smooth expansion should accept its CFL steps.");
    require(run.result.minimum_density > 0.0 && run.result.minimum_pressure > 0.0, "Expansion lost positivity.");
    require(run.maximum_axial_velocity < 1.0e-12, "Spurious axial flow.");
    require(run.source_relative_error < 1.0e-13, "Geometric source coefficient mismatch.");
    for (const double defect : run.global_balance_error)
        require(std::isfinite(defect) && std::abs(defect) < 1.0e-12, "Global boundary/source balance failure.");
    for (const double defect : run.maximum_step_balance_error)
        require(std::isfinite(defect) && defect < 1.0e-12, "Step balance failure.");
    for (const auto &norm : run.errors)
        require(std::isfinite(norm.l1) && std::isfinite(norm.rms) && std::isfinite(norm.linf),
                "Non-finite expansion error.");
}

void test_strength_rest_refinement_and_geometric_scaling()
{
    const cfd::verification::RadialGrid grid{.nx = 3, .nr = 12, .radius = 2.0};
    for (const double rate : {0.0, 0.05, 2.0})
    {
        const cfd::verification::RadialExpansionExact exact{{.alpha0 = rate}};
        const auto run{cfd::verification::run_radial_expansion(grid, exact, 0.2, 0.2)};
        require_run_valid(run);
        if (rate == 0.0)
            for (const auto &norm : run.errors)
                require(norm.linf < 1.0e-13, "Rest pressure/source cancellation failed.");
        else
            require(run.result.integrals.mass < 4.0 * std::numbers::pi,
                    "Expansion failed to lose mass through outer boundary.");
    }
    const cfd::verification::RadialExpansionExact exact;
    const auto coarse{cfd::verification::run_radial_expansion({.nx = 4, .nr = 16, .radius = 2.0}, exact, 0.2, 0.1)};
    const auto fine{cfd::verification::run_radial_expansion({.nx = 8, .nr = 32, .radius = 2.0}, exact, 0.2, 0.1)};
    require_run_valid(coarse);
    require_run_valid(fine);
    for (cfd::Index component = 0; component < coarse.errors.size(); ++component)
        require(fine.errors.at(component).rms < coarse.errors.at(component).rms,
                "Radial refinement did not reduce RMS error.");
    const auto unit{cfd::verification::run_radial_expansion(grid, exact, 0.2, 0.1)};
    for (const double scale : {0.01, 100.0})
    {
        const cfd::verification::RadialExpansionExact scaled_exact{{.alpha0 = 1.0 / scale}};
        const auto scaled{cfd::verification::run_radial_expansion(
            {.nx = 3, .nr = 12, .length = scale, .radius = 2.0 * scale}, scaled_exact, 0.2 * scale, 0.1)};
        require_run_valid(scaled);
        const auto differences{cfd::verification::radial_differences(grid, scaled.field, unit.field)};
        for (const auto &norm : differences)
            require(norm.linf < 2.0e-12, "Geometry/time scaling changed the dimensionless solution.");
    }
}
} // namespace

int main()
{
    int failures{};
    failures += cfd::test::run_test("radial exact volume projection", test_exact_moments_and_volume_projection);
    failures += cfd::test::run_test("radial exact dynamics and validation", test_exact_dynamics_and_invalid_parameters);
    failures += cfd::test::run_test("radial volume error norms", test_volume_weighted_norms);
    failures += cfd::test::run_test("moving radial source and boundary flux",
                                    test_moving_expansion_source_geometry_and_boundary_flux);
    failures += cfd::test::run_test("radial rest/strength/refinement/geometric scales",
                                    test_strength_rest_refinement_and_geometric_scaling);
    return cfd::test::finish_tests(failures, "radial expansion");
}
