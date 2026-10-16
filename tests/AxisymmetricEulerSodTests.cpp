#include "cfd/compressible/AxisymmetricEulerSolver.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "support/EulerTube.hpp"
#include "support/SodRiemann.hpp"
#include "support/TestUtils.hpp"
#include "support/VerificationStatistics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>

namespace
{
using cfd::test::require;
using cfd::test::require_near;

void test_independent_riemann_reference()
{
    const cfd::verification::SodRiemann exact;
    // Independent tabulated classical Sod star values, not computed by EulerFlux.
    const auto left_star{exact.sample(0.5, 0.2)};
    const auto right_star{exact.sample(0.75, 0.2)};
    require_near(left_star.pressure, 0.303130178050647, 1.0e-13, "Sod star pressure.");
    require_near(left_star.velocity, 0.927452620048950, 1.0e-13, "Sod contact speed.");
    require_near(left_star.density, 0.426319428178495, 1.0e-13, "Sod left star density.");
    require_near(right_star.density, 0.265573711705307, 1.0e-13, "Sod right star density.");
    require_near(exact.shock_speed(), 1.75215573203018, 1.0e-13, "Sod shock speed.");
    require(exact.sample(0.0, 0.2).density == 1.0 && exact.sample(1.0, 0.2).density == 0.125,
            "Sod reference boundaries must be untouched at t=0.2.");
    const auto fan{exact.sample(0.35, 0.2)};
    require_near(fan.pressure / std::pow(fan.density, 1.4), 1.0, 1.0e-13, "Rarefaction must be isentropic.");
}

std::array<double, 3> solve_sod(const cfd::Index nx)
{
    const auto built{cfd::build_mesh(cfd::verification::make_euler_tube(nx, 2))};
    const cfd::CaloricallyPerfectGas gas{1.4, 1.0};
    const cfd::verification::SodRiemann exact;
    cfd::ConservativeField initial{built.mesh.cell_count()};
    for (cfd::Index cell = 0; cell < initial.size(); ++cell)
    {
        const auto primitive{exact.sample(built.mesh.cell_centers()[cell].x, 0.0)};
        initial.set_state(cell, gas.to_conservative({primitive.density, primitive.velocity, 0.0, primitive.pressure}));
    }
    cfd::AxisymmetricEulerSolver solver{built.mesh, initial, gas, cfd::verification::euler_tube_conditions()};
    double maximum_balance{};
    while (solver.diagnostics().time < 0.15)
    {
        const auto &step{solver.step(0.15)};
        maximum_balance =
            std::max({maximum_balance, std::abs(step.balance_error.mass), std::abs(step.balance_error.momentum_x),
                      std::abs(step.balance_error.momentum_r), std::abs(step.balance_error.energy)});
    }
    require(maximum_balance < 1.0e-13, "Sod integral balance has excessive roundoff defect.");
    std::array<cfd::verification::ErrorAccumulator, 3> errors{};
    for (cfd::Index cell = 0; cell < initial.size(); ++cell)
    {
        const auto expected{exact.sample(built.mesh.cell_centers()[cell].x, 0.15)};
        const auto actual{gas.evaluate(solver.state().state(cell)).primitive};
        const double area{built.mesh.cell_areas()[cell]};
        errors[0].add(area, actual.density - expected.density);
        errors[1].add(area, actual.velocity_x - expected.velocity);
        errors[2].add(area, actual.pressure - expected.pressure);
        require(std::abs(actual.velocity_r) < 1.0e-13, "Sod generated radial velocity.");
    }
    const std::array<double, 3> rms{errors[0].finish().area_weighted_rms_error,
                                    errors[1].finish().area_weighted_rms_error,
                                    errors[2].finish().area_weighted_rms_error};
    std::cout << "Sod N=" << nx << " steps=" << solver.diagnostics().accepted_steps << " RMS(rho,u,p)=" << rms[0] << ','
              << rms[1] << ',' << rms[2] << " maximum balance defect=" << maximum_balance << '\n';
    return rms;
}

void test_sod_two_resolutions()
{
    const auto coarse{solve_sod(40)}, fine{solve_sod(80)};
    for (cfd::Index component = 0; component < coarse.size(); ++component)
    {
        require(fine.at(component) < coarse.at(component), "Sod RMS error must decrease under refinement.");
        // Broad quantitative guard for first-order Rusanov at t=0.15, not a
        // smooth second-order expectation across shocks/contact discontinuities.
        require(fine.at(component) < 0.12, "Sod solution differs substantially from the exact Riemann solution.");
    }
}
} // namespace

int main()
{
    int failures{};
    failures += cfd::test::run_test("independent exact Sod Riemann solution", test_independent_riemann_reference);
    failures += cfd::test::run_test("axisymmetric Sod refinement and conservation", test_sod_two_resolutions);
    return cfd::test::finish_tests(failures, "axisymmetric Euler Sod");
}
