#include "support/EulerRpaExploratory.hpp"
#include "support/EulerRpaInitialization.hpp"
#include "support/TestUtils.hpp"

#include "cfd/geometry/AxisymmetricGeometry.hpp"

#include <array>
#include <bit>
#include <cstdint>
#include <limits>

namespace
{
using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws;
using Model = cfd::verification::EulerRpaModel;

std::array<std::uint64_t, 4> bits(const cfd::ConservativeState &state)
{
    return {std::bit_cast<std::uint64_t>(state.density), std::bit_cast<std::uint64_t>(state.momentum_x),
            std::bit_cast<std::uint64_t>(state.momentum_r), std::bit_cast<std::uint64_t>(state.total_energy_density)};
}

cfd::MeshBuildResult nozzle(const bool named_wall = true)
{
    auto raw{cfd::verification::make_euler_tube(4, 2)};
    if (named_wall)
        raw.boundary_groups[1].name = "inner_wall"; // Deliberately not boundary ID 2.
    for (auto &point : raw.nodes)
        point.y *= .1 + .2 * std::abs(point.x - .5);
    return cfd::build_mesh(std::move(raw));
}

void test_quasi1d_geometry_totals_branches_and_exterior()
{
    const auto built{nozzle()};
    const auto &mesh{built.mesh};
    const cfd::AxisymmetricGeometry geometry{mesh};
    const std::array<cfd::Index, 8> regions{2, 1, 1, 1, 3, 4, 1, 5};
    for (const double gamma : {1.4, Model::gamma})
    {
        const cfd::CaloricallyPerfectGas gas{gamma, Model::gas_constant};
        const double p0{Model::total_pressure()}, t0{Model::total_temperature()};
        cfd::ConservativeField field{mesh.cell_count()};
        for (cfd::Index cell = 0; cell < field.size(); ++cell)
            field.set_state(cell, gas.to_conservative({1 + .1 * static_cast<double>(cell), 2, -3, 10000}));
        const cfd::ConservativeField original{field};
        const auto result{cfd::verification::initialize_rpa_quasi1d(mesh, regions, gas, p0, t0, field)};
        require_near(result.throat_x, .5, 0, "Minimum wall radius at known throat plane.");
        require_near(result.throat_radius, .1, 0, "Known throat radius.");
        require_near(result.critical_area, .01 * std::numbers::pi, 2e-15, "Disk throat area.");
        double section_area{};
        for (cfd::Index face = 0; face < mesh.face_count(); ++face)
        {
            const auto nodes{mesh.faces()[face].node_ids};
            if (mesh.nodes()[nodes[0]].x == .5 && mesh.nodes()[nodes[1]].x == .5)
            {
                const auto surface{geometry.face_area_vectors()[face]};
                section_area += std::hypot(surface.x, surface.y);
            }
        }
        require_near(result.critical_area, section_area, 2e-15, "Actual revolved throat section.");
        // Independent sonic state: critical flow = rho_star * sound_speed_star * disk area.
        const double tstar{2 * t0 / (gamma + 1)}, pstar{p0 * std::pow(tstar / t0, gamma / (gamma - 1))};
        const double expected_flow{(pstar / (gas.gas_constant() * tstar)) *
                                   std::sqrt(gamma * gas.gas_constant() * tstar) * (.01 * std::numbers::pi)};
        require_near(result.critical_mass_flow, expected_flow, 3e-14, "Critical flow from independent sonic state.");
        bool subsonic{}, supersonic{};
        for (cfd::Index cell = 0; cell < field.size(); ++cell)
        {
            if (regions.at(cell) != 1)
            {
                require(bits(field.state(cell)) == bits(original.state(cell)), "Exterior conservative bits changed.");
                continue;
            }
            const auto value{gas.evaluate(field.state(cell))};
            const auto &p{value.primitive};
            const double t{value.thermodynamics.temperature}, mach{p.velocity_x / value.thermodynamics.sound_speed};
            const double x{mesh.cell_centers()[cell].x};
            require(p.density > 0 && p.pressure > 0 && t > 0 && p.velocity_r == 0, "Admissible axial state.");
            require(std::isfinite(field.state(cell).total_energy_density), "Finite conservative energy.");
            if (x < .5)
            {
                require(mach > 0 && mach < 1, "Upstream branch must be subsonic.");
                subsonic = true;
            }
            else
            {
                require(mach > 1, "Downstream branch must be supersonic.");
                supersonic = true;
            }
            const double theta{1 + (gamma - 1) * mach * mach / 2};
            require_near(t + p.velocity_x * p.velocity_x / (2 * gamma * gas.gas_constant() / (gamma - 1)), t0, 3e-14,
                         "Total enthalpy / temperature.");
            require_near(p.pressure * std::pow(theta, gamma / (gamma - 1)), p0, 3e-14, "Total pressure.");
            const double radius{.1 + .2 * std::abs(x - .5)}, area{std::numbers::pi * radius * radius};
            const double area_ratio{std::pow(2 * theta / (gamma + 1), (gamma + 1) / (2 * (gamma - 1))) / mach};
            require_near(area / (.01 * std::numbers::pi), area_ratio, 3e-14, "Area-Mach equation, both branches.");
            require_near(field.state(cell).momentum_x * area, expected_flow, 3e-14, "Analytical axial mass flow.");
        }
        require(subsonic && supersonic, "Both branches must be exercised.");
    }
}

void test_uniform_historical_state_is_bit_exact()
{
    auto raw{cfd::verification::make_euler_tube(4, 1)};
    for (auto &point : raw.nodes)
    {
        point.x *= 2;
        point.y *= .02;
    }
    const auto built{cfd::build_mesh(std::move(raw))};
    const cfd::RpaNozzleDomain domain{.inner_wall = {{0, .1}, {.5, .05}, {1, .1}},
                                      .outer_wall = {{0, .12}, {1, .12}},
                                      .curve_junctions = {},
                                      .throat = {.5, .05},
                                      .downstream_x = 2,
                                      .farfield_radius = 1};
    const cfd::RpaRegionLayout layout{domain, {.band_thickness_factor = .5, .plume_outlet_radius_factor = 6}};
    const auto ambient{Model::ambient(300, 0, 0)};
    const cfd::CaloricallyPerfectGas gas{Model::gamma, Model::gas_constant};
    cfd::ConservativeField field{4};
    std::array<double, 4> regions{};
    require(cfd::verification::initialize_rpa(built.mesh, layout, ambient, field, regions) == 2,
            "Historical partition must be retained.");
    for (cfd::Index cell = 0; cell < field.size(); ++cell)
    {
        const auto expected{gas.to_conservative(cell < 2 ? Model::hot() : ambient)};
        require(bits(field.state(cell)) == bits(expected), "Historical uniform conservative state changed.");
        require(regions.at(cell) == (cell < 2 ? 1 : 3), "Historical region label changed.");
    }
}

void test_invalid_initialization_inputs()
{
    const auto built{nozzle()}, unnamed{nozzle(false)};
    const cfd::CaloricallyPerfectGas gas{Model::gamma, Model::gas_constant};
    const std::array<cfd::Index, 8> regions{1, 1, 1, 1, 2, 3, 4, 5};
    cfd::ConservativeField field{8, gas.to_conservative(Model::ambient(300, 0, 0))};
    for (const double invalid :
         {0.0, -1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
    {
        require_throws<std::invalid_argument>(
            [&] { (void)cfd::verification::initialize_rpa_quasi1d(built.mesh, regions, gas, invalid, 300, field); },
            "Invalid total pressure.");
        require_throws<std::invalid_argument>(
            [&] { (void)cfd::verification::initialize_rpa_quasi1d(built.mesh, regions, gas, 1e6, invalid, field); },
            "Invalid total temperature.");
    }
    require_throws<std::invalid_argument>(
        [&] {
            (void)cfd::verification::initialize_rpa_quasi1d(built.mesh, std::span{regions}.first(7), gas, 1e6, 300,
                                                            field);
        },
        "Invalid region cardinality.");
    cfd::ConservativeField short_field{7};
    require_throws<std::invalid_argument>(
        [&] { (void)cfd::verification::initialize_rpa_quasi1d(built.mesh, regions, gas, 1e6, 300, short_field); },
        "Invalid field cardinality.");
    require_throws<std::invalid_argument>(
        [&] { (void)cfd::verification::initialize_rpa_quasi1d(unnamed.mesh, regions, gas, 1e6, 300, field); },
        "Missing wall name.");
    const std::array<cfd::Index, 8> exterior{2, 2, 3, 3, 4, 4, 5, 5};
    require_throws<std::invalid_argument>(
        [&] { (void)cfd::verification::initialize_rpa_quasi1d(built.mesh, exterior, gas, 1e6, 300, field); },
        "No Engine cells.");
}
} // namespace

int main()
{
    int failures{};
    failures += cfd::test::run_test("quasi1d geometry, totals, branches and exterior",
                                    test_quasi1d_geometry_totals_branches_and_exterior);
    failures += cfd::test::run_test("uniform historical conservative bits", test_uniform_historical_state_is_bit_exact);
    failures += cfd::test::run_test("invalid quasi1d data", test_invalid_initialization_inputs);
    return cfd::test::finish_tests(failures, "RPA initialization");
}
