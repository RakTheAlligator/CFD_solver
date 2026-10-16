#include "cfd/compressible/EulerCharacteristicBoundaryFlux.hpp"
#include "cfd/io/VtkWriter.hpp"
#include "support/EulerBoundaryWave.hpp"

#include <array>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <span>
#include <string>

namespace
{
constexpr std::array dimensionless_times{0.0, 0.1, 0.25, 0.45, 0.65, 0.85};

struct OutputFields
{
    explicit OutputFields(const cfd::Index count)
        : pressure(count), temperature(count), ux(count), ur(count), speed(count), mach(count)
    {
    }

    cfd::CellScalarField pressure, temperature, ux, ur, speed, mach;
};

void write_snapshot(const cfd::Mesh &mesh, const cfd::AxisymmetricEulerSolver &solver,
                    const cfd::CaloricallyPerfectGas &gas, const cfd::EulerBoundaryConditions &conditions,
                    const std::filesystem::path &path, OutputFields &output, std::ostream &boundary_csv)
{
    for (cfd::Index cell = 0; cell < mesh.cell_count(); ++cell)
    {
        const auto evaluated{gas.evaluate(solver.state().state(cell))};
        output.pressure[cell] = evaluated.primitive.pressure;
        output.temperature[cell] = evaluated.thermodynamics.temperature;
        output.ux[cell] = evaluated.primitive.velocity_x;
        output.ur[cell] = evaluated.primitive.velocity_r;
        output.speed[cell] = std::hypot(output.ux[cell], output.ur[cell]);
        output.mach[cell] = output.speed[cell] / evaluated.thermodynamics.sound_speed;
    }
    const std::array fields{cfd::VtkCellScalarData{"rho", solver.state().density()},
                            cfd::VtkCellScalarData{"p", output.pressure.values()},
                            cfd::VtkCellScalarData{"T", output.temperature.values()},
                            cfd::VtkCellScalarData{"ux", output.ux.values()},
                            cfd::VtkCellScalarData{"ur", output.ur.values()},
                            cfd::VtkCellScalarData{"velocity_magnitude", output.speed.values()},
                            cfd::VtkCellScalarData{"Mach", output.mach.values()},
                            cfd::VtkCellScalarData{"rhoE", solver.state().total_energy_density()}};
    cfd::write_vtu(mesh, path, {.scalars = fields});
    const cfd::AxisymmetricGeometry geometry{mesh};
    for (cfd::Index face = 0; face < mesh.face_count(); ++face)
    {
        if (!mesh.face_adjacencies()[face].is_boundary())
            continue;
        const auto group{mesh.face_boundary_ids()[face]};
        if (group != 2 && group != 3)
            continue;
        const auto sf{mesh.face_area_vectors()[face]};
        const double length{mesh.face_lengths()[face]};
        const auto result{cfd::compute_euler_characteristic_boundary_flux(
            solver.state().state(mesh.face_adjacencies()[face].owner), {sf.x / length, sf.y / length},
            conditions.condition(group), gas)};
        const auto area{geometry.face_area_vectors()[face]};
        boundary_csv << solver.diagnostics().time << ',' << face << ',' << mesh.boundary_groups()[group].name << ','
                     << static_cast<int>(conditions[group]) << ',' << result.normal_mach << ','
                     << result.flux.mass * std::hypot(area.x, area.y) << '\n';
    }
}

void write_demo(const std::filesystem::path &directory)
{
    // Equivalent already-hot gas frozen at the RPA nozzle inlet, not a
    // reacting or geometrically representative nozzle simulation.
    constexpr double gamma{1.177}, gas_constant{414.2}, pressure{1170000.0}, temperature{3343.5658};
    const cfd::verification::EulerBoundaryWaveParameters parameters{.gamma = gamma,
                                                                    .gas_constant = gas_constant,
                                                                    .density = pressure / (gas_constant * temperature),
                                                                    .pressure = pressure,
                                                                    .mach = 0.1508};
    const cfd::CaloricallyPerfectGas gas{gamma, gas_constant};
    const double sound{std::sqrt(gamma * gas_constant * temperature)};
    const auto built{cfd::verification::build_boundary_wave_mesh(120, 3, 1.0, parameters.radius)};
    const auto conditions{cfd::verification::boundary_wave_conditions(parameters)};
    cfd::AxisymmetricEulerSolver solver{built.mesh, cfd::verification::boundary_wave_initial(built.mesh, parameters),
                                        gas, conditions};
    cfd::verification::EulerBoundaryWaveBalance balance{.initial = solver.diagnostics().integrals};
    OutputFields fields{built.mesh.cell_count()};
    std::filesystem::create_directories(directory);
    std::ofstream pvd{directory / "boundary_demo.pvd"}, csv{directory / "boundary_faces.csv"};
    if (!pvd || !csv)
        throw std::runtime_error("Cannot open Euler boundary demonstration outputs.");
    pvd << std::setprecision(17)
        << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" "
           "byte_order=\"LittleEndian\">\n<Collection>\n";
    csv << std::setprecision(17) << "time,face,group,boundary_type,normal_mach,integrated_outward_mass_flux\n";
    double maximum_defect{};
    for (cfd::Index index = 0; index < dimensionless_times.size(); ++index)
    {
        const double time{dimensionless_times.at(index) / sound};
        while (solver.diagnostics().time < time)
            balance.add(solver.step(time));
        const auto name{"boundary_" + std::to_string(index) + ".vtu"};
        write_snapshot(built.mesh, solver, gas, conditions, directory / name, fields, csv);
        pvd << "<DataSet timestep=\"" << time << "\" group=\"\" part=\"0\" file=\"" << name << "\"/>\n";
    }
    pvd << "</Collection>\n</VTKFile>\n";
    pvd.close();
    csv.close();
    if (!pvd || !csv)
        throw std::runtime_error("Cannot finish Euler boundary demonstration outputs.");
    for (const double component : balance.defect(solver.diagnostics().integrals))
        maximum_defect = std::max(maximum_defect, std::abs(component));
    if (!std::isfinite(maximum_defect) || maximum_defect > 1.0e-8 || solver.diagnostics().rejected_steps != 0)
        throw std::runtime_error("Hot-gas demonstration conservation or positivity verification failed.");
    std::cout << "HOT GAS gamma=" << gamma << " R=" << gas_constant << " rho=" << parameters.density
              << " steps=" << solver.diagnostics().accepted_steps << " rejects=" << solver.diagnostics().rejected_steps
              << " rho_min=" << solver.diagnostics().minimum_density
              << " p_min=" << solver.diagnostics().minimum_pressure << " global_balance_max=" << maximum_defect
              << " max_step_balance=" << balance.maximum_step_error << "\nPVD: " << directory / "boundary_demo.pvd"
              << '\n';
}
} // namespace

int main(int argc, char **argv)
{
    try
    {
        if (argc > 2)
            throw std::invalid_argument("Usage: cfd_axisymmetric_euler_boundary_demo [output_directory]");
        const std::span<char *> arguments{argv, static_cast<std::size_t>(argc)};
        const std::filesystem::path directory{argc == 2 ? arguments[1] : "output/euler_boundary"};
        std::cout << std::scientific << std::setprecision(10);
        for (const cfd::Index nx : {60U, 120U})
            for (const bool transmissive : {false, true})
            {
                const auto result{cfd::verification::compare_boundary_wave(nx, transmissive)};
                if (result.result.rejected_steps != 0 || result.maximum_balance_defect > 1.0e-12 ||
                    !std::isfinite(result.reflection_ratio) || !(result.incident_pressure > 0.0))
                    throw std::runtime_error("Acoustic boundary comparison conservation/admissibility check failed.");
                std::cout << "WAVE Nx=" << nx << " outlet=" << (transmissive ? "Transmissive" : "PressureOutlet")
                          << " steps=" << result.result.accepted_steps << " rejects=" << result.result.rejected_steps
                          << " incident_p=" << result.incident_pressure << " outgoing_p=" << result.outgoing_pressure
                          << " short_extended_RMS=" << result.pressure_rms_difference
                          << " reflected_p=" << result.reflected_pressure
                          << " min_delta_p=" << result.minimum_pressure_difference
                          << " incoming_amplitude=" << result.reflected_incoming_amplitude
                          << " reflection_ratio=" << result.reflection_ratio
                          << " balance_max=" << result.maximum_balance_defect << '\n';
            }
        write_demo(directory);
        std::cout << "Euler boundary verification passed. Fixed pressure is not a nonreflecting outlet.\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
