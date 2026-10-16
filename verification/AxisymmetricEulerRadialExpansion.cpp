#include "cfd/io/VtkWriter.hpp"
#include "support/RadialExpansionStudy.hpp"

#include <array>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <span>
#include <string>

namespace
{
namespace verification = cfd::verification;
constexpr std::array names{"rho", "rho_ur", "rhoE"};
constexpr std::array snapshot_times{0.0, 0.05, 0.1, 0.2};
constexpr double conservation_tolerance{1.0e-12};

void validate_run(const verification::RadialRun &run)
{
    if (run.result.time != 0.2 || run.result.accepted_steps == 0 || !(run.result.minimum_density > 0.0) ||
        !(run.result.minimum_pressure > 0.0) || !std::isfinite(run.source_relative_error) ||
        run.source_relative_error > conservation_tolerance || run.maximum_axial_velocity > conservation_tolerance)
        throw std::runtime_error("Radial time/admissibility/source/axial-invariance check failed.");
    for (cfd::Index component = 0; component < run.global_balance_error.size(); ++component)
        if (!std::isfinite(run.global_balance_error.at(component)) ||
            std::abs(run.global_balance_error.at(component)) > conservation_tolerance ||
            !std::isfinite(run.maximum_step_balance_error.at(component)) ||
            run.maximum_step_balance_error.at(component) > conservation_tolerance)
            throw std::runtime_error("Independent radial boundary/source conservation check failed.");
}

void spatial_campaign(const double radius)
{
    const verification::RadialExpansionExact exact;
    std::array<verification::RadialErrorNorms, 3> previous{};
    std::cout << "SPATIAL R=" << radius << " ALL cells, volume-weighted, CFL=0.1\n";
    for (const cfd::Index nr : {20U, 40U, 80U})
    {
        const auto run{verification::run_radial_expansion({.nx = nr / 4, .nr = nr, .radius = radius}, exact, 0.2, 0.1)};
        validate_run(run);
        std::cout << "Nr=" << nr << " Nx=" << nr / 4 << " steps=" << run.result.accepted_steps
                  << " rejected=" << run.result.rejected_steps << " rho_min=" << run.result.minimum_density
                  << " p_min=" << run.result.minimum_pressure << '\n';
        for (cfd::Index component = 0; component < names.size(); ++component)
        {
            const auto error{run.errors.at(component)};
            std::cout << "  " << names.at(component) << " L1=" << error.l1 << " RMS=" << error.rms
                      << " Linf=" << error.linf;
            if (nr > 20)
            {
                const auto coarse{previous.at(component)};
                if (!(error.l1 < coarse.l1) || !(error.rms < coarse.rms) || !(error.linf < coarse.linf))
                    throw std::runtime_error("Radial spatial L1/RMS/Linf error did not strictly decrease.");
                std::cout << " orders L1/RMS/Linf=" << std::log2(coarse.l1 / error.l1) << '/'
                          << std::log2(coarse.rms / error.rms) << '/' << std::log2(coarse.linf / error.linf);
            }
            std::cout << '\n';
        }
        previous = run.errors;
        std::cout << "  initial_boundary_mass_relative_error=" << run.initial_boundary_mass_relative_error
                  << " p_point_RMS=" << run.pressure_point_rms << " p_projected_RMS=" << run.pressure_projected_rms
                  << "\n  axis_density_error=" << run.axis_density_error
                  << " outer_density_error=" << run.outer_density_error << " max|ux|=" << run.maximum_axial_velocity
                  << " source_relative_error=" << run.source_relative_error
                  << "\n  scale-relative global defects (M,Px,Pr,E): ";
        for (const auto defect : run.global_balance_error)
            std::cout << defect << ' ';
        std::cout << "\n  maximum scale-relative step defects: ";
        for (const auto defect : run.maximum_step_balance_error)
            std::cout << defect << ' ';
        std::cout << '\n';
    }
}

struct RadialOutput
{
    explicit RadialOutput(const cfd::Index count)
        : pressure(count), ux(count), ur(count), speed(count), mach(count), density_error(count), momentum_error(count),
          energy_error(count), projected_pressure(count), point_pressure(count)
    {
    }

    cfd::CellScalarField pressure, ux, ur, speed, mach, density_error, momentum_error, energy_error, projected_pressure,
        point_pressure;
};

void write_snapshot(const std::filesystem::path &path, const cfd::Mesh &mesh, const verification::RadialGrid &grid,
                    const cfd::ConservativeField &field, const double time, RadialOutput &output)
{
    const verification::RadialExpansionExact exact;
    const cfd::CaloricallyPerfectGas gas{exact.gamma(), 1.0};
    for (cfd::Index cell = 0; cell < field.size(); ++cell)
    {
        const auto state{field.state(cell)};
        const auto evaluated{gas.evaluate(state)};
        const auto bounds{verification::radial_bounds(grid, cell)};
        const auto reference{exact.cell_average(bounds[0], bounds[1], time)};
        output.pressure[cell] = evaluated.primitive.pressure;
        output.ux[cell] = evaluated.primitive.velocity_x;
        output.ur[cell] = evaluated.primitive.velocity_r;
        output.speed[cell] = std::hypot(output.ux[cell], output.ur[cell]);
        output.mach[cell] = output.speed[cell] / evaluated.thermodynamics.sound_speed;
        output.density_error[cell] = state.density - reference.density;
        output.momentum_error[cell] = state.momentum_r - reference.momentum_r;
        output.energy_error[cell] = state.total_energy_density - reference.total_energy_density;
        output.projected_pressure[cell] = exact.projected_pressure(bounds[0], bounds[1], time);
        output.point_pressure[cell] = exact.values(time).pressure;
    }
    const std::array data{cfd::VtkCellScalarData{"rho", field.density()},
                          cfd::VtkCellScalarData{"p", output.pressure.values()},
                          cfd::VtkCellScalarData{"ux", output.ux.values()},
                          cfd::VtkCellScalarData{"ur", output.ur.values()},
                          cfd::VtkCellScalarData{"velocity_magnitude", output.speed.values()},
                          cfd::VtkCellScalarData{"Mach", output.mach.values()},
                          cfd::VtkCellScalarData{"rhoE", field.total_energy_density()},
                          cfd::VtkCellScalarData{"rho_error", output.density_error.values()},
                          cfd::VtkCellScalarData{"rho_ur_error", output.momentum_error.values()},
                          cfd::VtkCellScalarData{"rhoE_error", output.energy_error.values()},
                          cfd::VtkCellScalarData{"p_exact_projected", output.projected_pressure.values()},
                          cfd::VtkCellScalarData{"p_exact_point", output.point_pressure.values()}};
    cfd::write_vtu(mesh, path, {.scalars = data});
}

void write_visualization(const std::filesystem::path &directory)
{
    const verification::RadialExpansionExact exact;
    const verification::RadialGrid grid{.nx = 20, .nr = 80, .radius = 2.0};
    const auto built{verification::build_radial_grid(grid)};
    cfd::AxisymmetricEulerSolver solver{built.mesh,
                                        verification::radial_initial_field(grid, exact),
                                        cfd::CaloricallyPerfectGas{exact.gamma(), 1.0},
                                        verification::radial_conditions(),
                                        {.cfl = 0.1}};
    RadialOutput output{built.mesh.cell_count()};
    std::filesystem::create_directories(directory);
    std::ofstream pvd{directory / "radial_expansion.pvd"};
    if (!pvd)
        throw std::runtime_error("Cannot open radial expansion PVD.");
    pvd << std::setprecision(17)
        << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" "
           "byte_order=\"LittleEndian\">\n<Collection>\n";
    for (cfd::Index index = 0; index < snapshot_times.size(); ++index)
    {
        const double time{snapshot_times.at(index)};
        solver.advance_to(time);
        const auto name{"radial_" + std::to_string(index) + ".vtu"};
        write_snapshot(directory / name, built.mesh, grid, solver.state(), time, output);
        pvd << "<DataSet timestep=\"" << time << "\" group=\"\" part=\"0\" file=\"" << name << "\"/>\n";
    }
    pvd << "</Collection>\n</VTKFile>\n";
    pvd.close();
    if (!pvd)
        throw std::runtime_error("Cannot write radial expansion PVD.");
    std::cout << "ParaView (separate solve with output-time clipping): " << directory / "radial_expansion.pvd" << '\n';
}

} // namespace

int main(int argc, char **argv)
{
    try
    {
        if (argc > 2)
            throw std::invalid_argument("Usage: cfd_axisymmetric_euler_radial_expansion [output_directory]");
        const std::span arguments{argv, static_cast<cfd::Index>(argc)};
        const std::filesystem::path output{argc == 2 ? arguments[1] : "radial_expansion_results"};
        const cfd::verification::RadialExpansionExact exact;
        std::cout << std::scientific << std::setprecision(9);
        spatial_campaign(1.0);
        for (const cfd::Index nr : {40U, 80U})
        {
            const auto comparison{cfd::verification::compare_radial_outer_boundary(nr)};
            std::cout << "BOUNDARY CONTROL Nr=" << nr << " matched_steps=" << comparison.matched_steps << '\n';
            for (cfd::Index component = 0; component < comparison.difference.size(); ++component)
                std::cout << "  " << names.at(component)
                          << " difference RMS=" << comparison.difference.at(component).rms
                          << " original RMS=" << comparison.original_error.at(component).rms
                          << " extended RMS=" << comparison.extended_error.at(component).rms << '\n';
        }
        spatial_campaign(2.0);
        const cfd::verification::RadialGrid time_grid{.nx = 20, .nr = 80, .radius = 2.0};
        const std::array temporal{cfd::verification::run_radial_expansion(time_grid, exact, 0.2, 0.4),
                                  cfd::verification::run_radial_expansion(time_grid, exact, 0.2, 0.2),
                                  cfd::verification::run_radial_expansion(time_grid, exact, 0.2, 0.1),
                                  cfd::verification::run_radial_expansion(time_grid, exact, 0.2, 0.025)};
        constexpr std::array cfls{0.4, 0.2, 0.1, 0.025};
        std::array<cfd::verification::RadialErrorNorms, 3> previous_time_error{};
        for (cfd::Index i = 0; i < temporal.size(); ++i)
        {
            validate_run(temporal.at(i));
            std::cout << "TIME cfl=" << cfls.at(i) << " steps=" << temporal.at(i).result.accepted_steps
                      << " rejects=" << temporal.at(i).result.rejected_steps << '\n';
            const auto delta{
                cfd::verification::radial_differences(time_grid, temporal.at(i).field, temporal.back().field)};
            for (cfd::Index k = 0; k < delta.size(); ++k)
            {
                std::cout << "  " << names.at(k) << " analytic_RMS=" << temporal.at(i).errors.at(k).rms
                          << " time_reference_L1/RMS/Linf=" << delta.at(k).l1 << '/' << delta.at(k).rms << '/'
                          << delta.at(k).linf << '\n';
                if (i > 0 && i + 1 < temporal.size() && !(delta.at(k).rms < previous_time_error.at(k).rms))
                    throw std::runtime_error("Temporal RMS differences to the fine-step reference did not decrease.");
            }
            previous_time_error = delta;
        }
        const auto first{cfd::verification::radial_differences(time_grid, temporal[0].field, temporal[1].field)};
        const auto second{cfd::verification::radial_differences(time_grid, temporal[1].field, temporal[2].field)};
        for (cfd::Index k = 0; k < first.size(); ++k)
            std::cout << names.at(k) << " TIME SELF ORDER=" << std::log2(first.at(k).rms / second.at(k).rms) << '\n';
        write_visualization(output);
        std::cout << "Radial expansion verification passed. No production BC or solver changes.\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
