#include "cfd/compressible/EulerFarfieldBoundaryFlux.hpp"
#include "cfd/io/VtkWriter.hpp"
#include "support/EulerBoundaryWave.hpp"

#include <array>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <span>
#include <string>

namespace
{
using Clock = std::chrono::steady_clock;

struct OutputFields
{
    explicit OutputFields(const cfd::Index size)
        : p(size), temperature(size), ux(size), ur(size), speed(size), mach(size)
    {
    }

    cfd::CellScalarField p, temperature, ux, ur, speed, mach;
};

void snapshot(const cfd::Mesh &mesh, const cfd::AxisymmetricEulerSolver &solver, const cfd::CaloricallyPerfectGas &gas,
              OutputFields &output, const std::filesystem::path &path)
{
    for (cfd::Index cell = 0; cell < mesh.cell_count(); ++cell)
    {
        const auto state{gas.evaluate(solver.state().state(cell))};
        output.p[cell] = state.primitive.pressure;
        output.temperature[cell] = state.thermodynamics.temperature;
        output.ux[cell] = state.primitive.velocity_x;
        output.ur[cell] = state.primitive.velocity_r;
        output.speed[cell] = std::hypot(output.ux[cell], output.ur[cell]);
        output.mach[cell] = output.speed[cell] / state.thermodynamics.sound_speed;
    }
    const std::array fields{cfd::VtkCellScalarData{"rho", solver.state().density()},
                            cfd::VtkCellScalarData{"p", output.p.values()},
                            cfd::VtkCellScalarData{"T", output.temperature.values()},
                            cfd::VtkCellScalarData{"ux", output.ux.values()},
                            cfd::VtkCellScalarData{"ur", output.ur.values()},
                            cfd::VtkCellScalarData{"velocity_magnitude", output.speed.values()},
                            cfd::VtkCellScalarData{"Mach", output.mach.values()},
                            cfd::VtkCellScalarData{"rhoE", solver.state().total_energy_density()}};
    cfd::write_vtu(mesh, path, {.scalars = fields});
}

struct BoundaryCounts
{
    std::array<cfd::Index, 4> regimes{};
    double mass_flow{};
    cfd::Index evaluations{};
    double evaluation_seconds{}; ///< Observation-only farfield calls, excluding CSV; includes clock overhead.
};

BoundaryCounts boundary_diagnostics(const cfd::Mesh &mesh, const cfd::AxisymmetricEulerSolver &solver,
                                    const cfd::CaloricallyPerfectGas &gas,
                                    const cfd::EulerBoundaryConditions &conditions, std::ostream &csv)
{
    BoundaryCounts right;
    for (cfd::Index face = 0; face < mesh.face_count(); ++face)
    {
        const auto adjacency{mesh.face_adjacencies()[face]};
        if (!adjacency.is_boundary())
            continue;
        const auto group{mesh.face_boundary_ids()[face]};
        if (conditions[group] != cfd::EulerBoundaryConditionType::Farfield)
            continue;
        const auto sf{mesh.face_area_vectors()[face]};
        const auto area{solver.geometry().face_area_vectors()[face]};
        const auto start{Clock::now()};
        const auto result{cfd::compute_euler_farfield_boundary_flux(
            solver.state().state(adjacency.owner), {sf.x / mesh.face_lengths()[face], sf.y / mesh.face_lengths()[face]},
            conditions.condition(group), gas)};
        right.evaluation_seconds += std::chrono::duration<double>(Clock::now() - start).count();
        ++right.evaluations;
        const double flux{result.flux.mass * std::hypot(area.x, area.y)};
        csv << solver.diagnostics().time << ',' << face << ',' << mesh.boundary_groups()[group].name << ','
            << static_cast<int>(result.regime) << ',' << result.normal_mach << ',' << flux << '\n';
        if (group == 3)
        {
            ++right.regimes.at(static_cast<cfd::Index>(result.regime));
            right.mass_flow += flux;
        }
    }
    return right;
}

void write_backflow_demo(const std::filesystem::path &directory)
{
    if (std::filesystem::exists(directory) &&
        (!std::filesystem::is_directory(directory) || !std::filesystem::is_empty(directory)))
        throw std::invalid_argument("Farfield output directory must be empty; existing results are never removed.");
    // Explicit normalized verification gas, NOT an RPA atmosphere. The pulse
    // temporarily reverses a weak baseline inflow through the right farfield.
    constexpr cfd::verification::EulerBoundaryWaveParameters parameters{.mach = -1.0e-4, .amplitude = 1.0e-3};
    const cfd::CaloricallyPerfectGas gas{parameters.gamma, parameters.gas_constant};
    const double sound{std::sqrt(parameters.gamma)};
    const auto built{cfd::verification::build_boundary_wave_mesh(120, 3, 1.0, parameters.radius)};
    const cfd::PrimitiveState exterior{1.0, parameters.mach * sound, 0.0, 1.0};
    const auto farfield{cfd::EulerBoundaryCondition::farfield(exterior)};
    const cfd::EulerBoundaryConditions conditions{
        4, {cfd::EulerBoundaryConditionType::Axis, cfd::EulerBoundaryConditionType::SlipWall, farfield, farfield}};
    cfd::AxisymmetricEulerSolver solver{built.mesh,
                                        cfd::verification::boundary_wave_initial(built.mesh, parameters),
                                        gas,
                                        conditions,
                                        {.collect_timings = true}};
    cfd::verification::EulerBoundaryWaveBalance balance{.initial = solver.diagnostics().integrals};
    OutputFields output{built.mesh.cell_count()};
    std::filesystem::create_directories(directory);
    std::ofstream pvd{directory / "farfield.pvd"}, csv{directory / "farfield_faces.csv"};
    if (!pvd || !csv)
        throw std::runtime_error("Cannot open farfield demonstration outputs.");
    pvd << std::setprecision(17)
        << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" "
           "byte_order=\"LittleEndian\">\n<Collection>\n";
    csv << std::setprecision(17) << "time,face,group,trace_regime,trace_normal_mach,integrated_outward_mass_flux\n";
    double diagnostic_seconds{}, export_seconds{}, farfield_sample_seconds{};
    cfd::Index farfield_sample_calls{};
    bool observed_inflow{}, observed_outflow{};
    const auto diagnose = [&] {
        const auto start{Clock::now()};
        const auto counts{boundary_diagnostics(built.mesh, solver, gas, conditions, csv)};
        farfield_sample_seconds += counts.evaluation_seconds;
        farfield_sample_calls += counts.evaluations;
        diagnostic_seconds += std::chrono::duration<double>(Clock::now() - start).count();
        observed_inflow = observed_inflow || counts.regimes.at(1) > 0;
        observed_outflow = observed_outflow || counts.regimes.at(2) > 0;
    };
    const auto start{Clock::now()};
    constexpr std::array times{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.65, 0.85};
    for (cfd::Index index = 0; index < times.size(); ++index)
    {
        const double time{times.at(index) / sound};
        while (solver.diagnostics().time < time)
        {
            balance.add(solver.step(time));
            if (solver.diagnostics().accepted_steps % 10 == 0)
                diagnose();
        }
        diagnose();
        const auto export_start{Clock::now()};
        const auto name{"farfield_" + std::to_string(index) + ".vtu"};
        snapshot(built.mesh, solver, gas, output, directory / name);
        pvd << "<DataSet timestep=\"" << solver.diagnostics().time << "\" group=\"\" part=\"0\" file=\"" << name
            << "\"/>\n";
        export_seconds += std::chrono::duration<double>(Clock::now() - export_start).count();
    }
    pvd << "</Collection>\n</VTKFile>\n";
    pvd.close();
    csv.close();
    if (!pvd || !csv)
        throw std::runtime_error("Cannot finish farfield outputs.");
    double defect{};
    for (const double component : balance.defect(solver.diagnostics().integrals))
        defect = std::max(defect, std::abs(component));
    if (!observed_inflow || !observed_outflow || !std::isfinite(defect) || defect > 1.0e-12 ||
        solver.diagnostics().rejected_steps != 0)
        throw std::runtime_error("Farfield backflow, conservation or positivity verification failed.");
    std::cout << "DEMO cells=" << built.mesh.cell_count() << " steps=" << solver.diagnostics().accepted_steps
              << " rejects=" << solver.diagnostics().rejected_steps << " inflow=" << observed_inflow
              << " outflow=" << observed_outflow << " rho_min=" << solver.diagnostics().minimum_density
              << " p_min=" << solver.diagnostics().minimum_pressure << " balance=" << defect
              << " step_balance=" << balance.maximum_step_error
              << " elapsed_s=" << std::chrono::duration<double>(Clock::now() - start).count()
              << " assembly_s=" << solver.timings().assembly_seconds
              << " BC_subset_s=" << solver.timings().boundary_flux_seconds
              << " update_s=" << solver.timings().time_update_seconds << " diagnostic_s=" << diagnostic_seconds
              << " export_s=" << export_seconds << " farfield_sample_calls=" << farfield_sample_calls
              << " farfield_sample_s=" << farfield_sample_seconds << "\nPVD: " << directory / "farfield.pvd" << '\n';
}
} // namespace

int main(int argc, char **argv)
{
    try
    {
        if (argc > 2)
            throw std::invalid_argument("Usage: cfd_axisymmetric_euler_farfield [output_directory]");
        const std::span<char *> arguments{argv, static_cast<cfd::Index>(argc)};
        const std::filesystem::path output{argc == 2 ? arguments[1] : "/tmp/cfd-f3-farfield"};
        std::cout << std::scientific << std::setprecision(10);
        for (const cfd::Index nx : {60U, 120U})
        {
            double pressure_reflection{};
            for (const int mode : {0, 1, 2})
            {
                const auto result{cfd::verification::compare_boundary_wave(nx, mode == 1, mode == 2)};
                if (mode == 0)
                    pressure_reflection = result.reflection_ratio;
                if (!std::isfinite(result.reflection_ratio) || !(result.incident_pressure > 0.0) ||
                    result.result.rejected_steps != 0 || result.maximum_balance_defect > 1.0e-12 ||
                    (mode == 2 && !(result.reflection_ratio < pressure_reflection)))
                    throw std::runtime_error("Farfield acoustic reflection/conservation verification failed.");
                std::cout << "ACOUSTIC nx=" << nx << " BC="
                          << (mode == 0   ? "PressureOutlet"
                              : mode == 1 ? "Transmissive"
                                          : "Farfield")
                          << " incident=" << result.incident_pressure << " reflected=" << result.reflected_pressure
                          << " reflection_ratio=" << result.reflection_ratio
                          << " incoming_amplitude=" << result.reflected_incoming_amplitude
                          << " short_extended_RMS=" << result.pressure_rms_difference
                          << " balance=" << result.maximum_balance_defect << '\n';
            }
        }
        write_backflow_demo(output);
        std::cout << "Farfield verification passed; no RPA simulation or exact 2D nonreflection claimed.\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
