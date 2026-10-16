#include "cfd/compressible/AxisymmetricEulerSolver.hpp"
#include "cfd/io/VtkWriter.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "support/EulerTube.hpp"
#include "support/SodRiemann.hpp"
#include "support/VerificationStatistics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <span>
#include <stdexcept>
#include <string>

namespace
{
constexpr std::array snapshot_times{0.0, 0.05, 0.1, 0.15, 0.2};
constexpr std::array<cfd::Index, 2> grid_levels{100, 200};
constexpr double balance_tolerance{2.0e-12};

struct OutputFields
{
    explicit OutputFields(const cfd::Index count)
        : pressure(count), temperature(count), ux(count), ur(count), speed(count), mach(count)
    {
    }

    cfd::CellScalarField pressure, temperature, ux, ur, speed, mach;
};

void write_snapshot(const cfd::Mesh &mesh, const cfd::ConservativeField &state, const cfd::CaloricallyPerfectGas &gas,
                    OutputFields &fields, const std::filesystem::path &path)
{
    for (cfd::Index cell = 0; cell < mesh.cell_count(); ++cell)
    {
        const auto evaluated{gas.evaluate(state.state(cell))};
        fields.pressure[cell] = evaluated.primitive.pressure;
        fields.temperature[cell] = evaluated.thermodynamics.temperature;
        fields.ux[cell] = evaluated.primitive.velocity_x;
        fields.ur[cell] = evaluated.primitive.velocity_r;
        fields.speed[cell] = std::hypot(fields.ux[cell], fields.ur[cell]);
        fields.mach[cell] = fields.speed[cell] / evaluated.thermodynamics.sound_speed;
    }
    const std::array data{cfd::VtkCellScalarData{"rho", state.density()},
                          cfd::VtkCellScalarData{"p", fields.pressure.values()},
                          cfd::VtkCellScalarData{"T", fields.temperature.values()},
                          cfd::VtkCellScalarData{"ux", fields.ux.values()},
                          cfd::VtkCellScalarData{"ur", fields.ur.values()},
                          cfd::VtkCellScalarData{"velocity_magnitude", fields.speed.values()},
                          cfd::VtkCellScalarData{"Mach", fields.mach.values()},
                          cfd::VtkCellScalarData{"rhoE", state.total_energy_density()}};
    cfd::write_vtu(mesh, path, {.scalars = data});
}

void write_series(const std::filesystem::path &directory)
{
    std::ofstream output{directory / "sod.pvd"};
    if (!output)
        throw std::runtime_error("Cannot open Sod PVD series.");
    output << std::setprecision(17)
           << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
              "  <Collection>\n";
    for (cfd::Index snapshot = 0; snapshot < snapshot_times.size(); ++snapshot)
        output << "    <DataSet timestep=\"" << snapshot_times.at(snapshot) << "\" group=\"\" part=\"0\" file=\"sod_"
               << snapshot << ".vtu\"/>\n";
    output << "  </Collection>\n</VTKFile>\n";
    output.close();
    if (!output)
        throw std::runtime_error("Cannot write Sod PVD series.");
}

struct SodErrors
{
    std::array<cfd::verification::ErrorStatistics, 3> fields;
};

SodErrors solve_level(const cfd::Index nx, const std::filesystem::path &root)
{
    const auto built{cfd::build_mesh(cfd::verification::make_euler_tube(nx, 4))};
    const auto &mesh{built.mesh};
    const cfd::CaloricallyPerfectGas gas{1.4, 1.0};
    const cfd::verification::SodRiemann exact;
    if (0.5 - std::sqrt(1.4) * snapshot_times.back() <= 0.0 || 0.5 + exact.shock_speed() * snapshot_times.back() >= 1.0)
        throw std::runtime_error("Sod verification must stop before boundary interaction.");
    cfd::ConservativeField initial{mesh.cell_count()};
    for (cfd::Index cell = 0; cell < mesh.cell_count(); ++cell)
    {
        const auto primitive{exact.sample(mesh.cell_centers()[cell].x, 0.0)};
        initial.set_state(cell, gas.to_conservative({primitive.density, primitive.velocity, 0.0, primitive.pressure}));
    }
    cfd::AxisymmetricEulerSolver solver{mesh, initial, gas, cfd::verification::euler_tube_conditions()};
    OutputFields output_fields{mesh.cell_count()};
    const auto directory{root / ("N" + std::to_string(nx))};
    std::filesystem::create_directories(directory);
    const auto initial_integrals{solver.diagnostics().integrals};
    // Independent time integration of measured boundary rates and radial source.
    std::array<long double, 4> integrated_rates{};
    double maximum_step_defect{}, maximum_radial_velocity{};
    for (cfd::Index snapshot = 0; snapshot < snapshot_times.size(); ++snapshot)
    {
        while (solver.diagnostics().time < snapshot_times.at(snapshot))
        {
            const auto &result{solver.step(snapshot_times.at(snapshot))};
            const long double dt{result.time_step};
            integrated_rates[0] -= dt * result.boundary_outflow_rate.mass;
            integrated_rates[1] -= dt * result.boundary_outflow_rate.momentum_x;
            integrated_rates[2] += dt * (static_cast<long double>(result.radial_pressure_source_rate) -
                                         result.boundary_outflow_rate.momentum_r);
            integrated_rates[3] -= dt * result.boundary_outflow_rate.energy;
            maximum_step_defect = std::max(
                {maximum_step_defect, std::abs(result.balance_error.mass), std::abs(result.balance_error.momentum_x),
                 std::abs(result.balance_error.momentum_r), std::abs(result.balance_error.energy)});
        }
        write_snapshot(mesh, solver.state(), gas, output_fields,
                       directory / ("sod_" + std::to_string(snapshot) + ".vtu"));
    }
    write_series(directory);
    std::array<cfd::verification::ErrorAccumulator, 3> errors{};
    for (cfd::Index cell = 0; cell < mesh.cell_count(); ++cell)
    {
        const auto expected{exact.sample(mesh.cell_centers()[cell].x, snapshot_times.back())};
        const auto actual{gas.evaluate(solver.state().state(cell)).primitive};
        const double area{mesh.cell_areas()[cell]};
        errors[0].add(area, actual.density - expected.density);
        errors[1].add(area, actual.velocity_x - expected.velocity);
        errors[2].add(area, actual.pressure - expected.pressure);
        maximum_radial_velocity = std::max(maximum_radial_velocity, std::abs(actual.velocity_r));
    }
    const auto result{solver.diagnostics()};
    const std::array<long double, 4> global_defects{
        result.integrals.mass - static_cast<long double>(initial_integrals.mass) - integrated_rates[0],
        result.integrals.momentum_x - static_cast<long double>(initial_integrals.momentum_x) - integrated_rates[1],
        result.integrals.momentum_r - static_cast<long double>(initial_integrals.momentum_r) - integrated_rates[2],
        result.integrals.energy - static_cast<long double>(initial_integrals.energy) - integrated_rates[3]};
    if (maximum_step_defect > balance_tolerance || maximum_radial_velocity > balance_tolerance)
        throw std::runtime_error("Sod per-step conservation or radial-invariance verification failed.");
    for (const auto defect : global_defects)
        if (!std::isfinite(defect) || std::abs(defect) > balance_tolerance)
            throw std::runtime_error("Sod accumulated boundary/source balance verification failed.");
    const SodErrors statistics{{errors[0].finish(), errors[1].finish(), errors[2].finish()}};
    std::cout << "N=" << nx << " cells=" << mesh.cell_count() << " time=" << result.time << " CFL=" << result.cfl
              << " steps=" << result.accepted_steps << " rejected=" << result.rejected_steps << " rho=["
              << result.minimum_density << ',' << result.maximum_density << "] p=[" << result.minimum_pressure << ','
              << result.maximum_pressure << "]\n";
    constexpr std::array names{"rho", "ux", "p"};
    for (cfd::Index component = 0; component < names.size(); ++component)
        std::cout << "  " << names.at(component) << " RMS=" << statistics.fields.at(component).area_weighted_rms_error
                  << " Linf=" << statistics.fields.at(component).linf_error << '\n';
    std::cout << "  max per-step balance defect=" << maximum_step_defect << " max|ur|=" << maximum_radial_velocity
              << "\n  global defects (mass,mx,mr,energy)=" << global_defects[0] << ',' << global_defects[1] << ','
              << global_defects[2] << ',' << global_defects[3] << "\n  ParaView: " << (directory / "sod.pvd") << '\n';
    return statistics;
}
} // namespace

int main(int argc, char **argv)
{
    try
    {
        if (argc > 2)
            throw std::invalid_argument("Usage: cfd_axisymmetric_euler_sod [output_directory]");
        const std::span arguments{argv, static_cast<cfd::Index>(argc)};
        const std::filesystem::path output{argc == 2 ? arguments[1] : "euler_sod_results"};
        std::cout << std::scientific << std::setprecision(9);
        const auto coarse{solve_level(grid_levels[0], output)};
        const auto fine{solve_level(grid_levels[1], output)};
        for (cfd::Index component = 0; component < coarse.fields.size(); ++component)
        {
            const double coarse_error{coarse.fields.at(component).area_weighted_rms_error};
            const double fine_error{fine.fields.at(component).area_weighted_rms_error};
            if (!(fine_error < coarse_error) || fine_error > 0.1)
                throw std::runtime_error("Sod quantitative RMS/refinement guard failed.");
            std::cout << "Observed RMS order component " << component << ": " << std::log2(coarse_error / fine_error)
                      << '\n';
        }
        std::cout << "Sod verification passed (no smooth-order requirement at discontinuities).\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "Sod verification failed: " << error.what() << '\n';
        return 1;
    }
}
