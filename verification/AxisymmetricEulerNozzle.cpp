#include "cfd/compressible/EulerCharacteristicBoundaryFlux.hpp"
#include "cfd/compressible/EulerSteadyMonitor.hpp"
#include "cfd/io/EulerSteadyCsvWriter.hpp"
#include "cfd/io/VtkWriter.hpp"
#include "support/EulerBoundaryWave.hpp"
#include "support/EulerNozzle.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace
{
using Clock = std::chrono::steady_clock;
using Reference = cfd::verification::EulerNozzleReference;
constexpr std::array<cfd::Index, 3> levels{24, 48, 96};
constexpr std::array output_times{0.1, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0}; ///< In acoustic crossing times.

double seconds(const Clock::time_point start)
{
    return std::chrono::duration<double>(Clock::now() - start).count();
}

struct OutputFields
{
    explicit OutputFields(const cfd::Index count)
        : pressure(count), temperature(count), ux(count), ur(count), speed(count), mach(count), total_pressure(count),
          total_temperature(count)
    {
    }

    cfd::CellScalarField pressure, temperature, ux, ur, speed, mach, total_pressure, total_temperature;
};

void snapshot(const std::filesystem::path &path, const cfd::Mesh &mesh, const cfd::ConservativeField &state,
              const cfd::CaloricallyPerfectGas &gas, OutputFields &fields)
{
    for (cfd::Index cell = 0; cell < state.size(); ++cell)
    {
        const auto evaluated{gas.evaluate(state.state(cell))};
        const auto primitive{evaluated.primitive};
        fields.pressure[cell] = primitive.pressure;
        fields.temperature[cell] = evaluated.thermodynamics.temperature;
        fields.ux[cell] = primitive.velocity_x;
        fields.ur[cell] = primitive.velocity_r;
        fields.speed[cell] = std::hypot(primitive.velocity_x, primitive.velocity_r);
        fields.mach[cell] = fields.speed[cell] / evaluated.thermodynamics.sound_speed;
        const double ratio{1.0 + 0.5 * (gas.gamma() - 1.0) * fields.mach[cell] * fields.mach[cell]};
        fields.total_temperature[cell] = fields.temperature[cell] * ratio;
        fields.total_pressure[cell] = primitive.pressure * std::pow(ratio, gas.gamma() / (gas.gamma() - 1.0));
    }
    const std::array data{cfd::VtkCellScalarData{"rho", state.density()},
                          cfd::VtkCellScalarData{"p", fields.pressure.values()},
                          cfd::VtkCellScalarData{"T", fields.temperature.values()},
                          cfd::VtkCellScalarData{"ux", fields.ux.values()},
                          cfd::VtkCellScalarData{"ur", fields.ur.values()},
                          cfd::VtkCellScalarData{"velocity_magnitude", fields.speed.values()},
                          cfd::VtkCellScalarData{"Mach", fields.mach.values()},
                          cfd::VtkCellScalarData{"rhoE", state.total_energy_density()},
                          cfd::VtkCellScalarData{"pressure_total", fields.total_pressure.values()},
                          cfd::VtkCellScalarData{"temperature_total", fields.total_temperature.values()}};
    cfd::write_vtu(mesh, path, {.scalars = data});
}

struct Section
{
    double area{}, mach{}, pressure{}, temperature{}, density{}, velocity{}, mass_flow{};
};

std::vector<Section> sections(const cfd::Mesh &mesh, const cfd::ConservativeField &field,
                              const cfd::AxisymmetricGeometry &geometry, const cfd::EulerBoundaryConditions &conditions,
                              const cfd::CaloricallyPerfectGas &gas, const cfd::Index nx)
{
    std::vector<Section> profile(nx + 1);
    for (cfd::Index face = 0; face < mesh.face_count(); ++face)
    {
        const auto nodes{mesh.faces()[face].node_ids};
        const cfd::Index column{nodes[0] % (nx + 1)};
        if (nodes[1] % (nx + 1) != column)
            continue; // Only the explicit axial mesh planes are integrated.
        const auto adjacency{mesh.face_adjacencies()[face]};
        const auto sf{mesh.face_area_vectors()[face]};
        const double length{mesh.face_lengths()[face]};
        const cfd::Vector2 normal{sf.x / length, sf.y / length};
        auto left{gas.evaluate(field.state(adjacency.owner)).primitive}, right{left};
        cfd::EulerFlux flux;
        if (adjacency.is_boundary())
        {
            const auto numerical{cfd::compute_euler_characteristic_boundary_flux(
                field.state(adjacency.owner), normal, conditions.condition(mesh.face_boundary_ids()[face]), gas)};
            flux = numerical.flux;
            left = gas.evaluate(numerical.boundary_state).primitive;
            right = left;
        }
        else
        {
            right = gas.evaluate(field.state(adjacency.neighbor)).primitive;
            flux = cfd::compute_rusanov_flux(field.state(adjacency.owner), field.state(adjacency.neighbor), normal, gas)
                       .flux;
        }
        const double area{std::abs(geometry.face_area_vectors()[face].x)};
        auto &section{profile[column]};
        section.area += area;
        section.pressure += area * 0.5 * (left.pressure + right.pressure);
        section.density += area * 0.5 * (left.density + right.density);
        section.velocity += area * 0.5 * (left.velocity_x + right.velocity_x);
        section.temperature += area * 0.5 *
                               (left.pressure / (gas.gas_constant() * left.density) +
                                right.pressure / (gas.gas_constant() * right.density));
        section.mach += area * 0.5 *
                        (std::hypot(left.velocity_x, left.velocity_r) / gas.properties(left).sound_speed +
                         std::hypot(right.velocity_x, right.velocity_r) / gas.properties(right).sound_speed);
        section.mass_flow += (normal.x > 0.0 ? area : -area) * flux.mass;
    }
    for (auto &section : profile)
    {
        if (!(section.area > 0.0))
            throw std::runtime_error("Nozzle axial section has no physical surface.");
        section.pressure /= section.area;
        section.density /= section.area;
        section.velocity /= section.area;
        section.temperature /= section.area;
        section.mach /= section.area;
    }
    return profile;
}

void write_profiles(const std::filesystem::path &directory, const std::vector<Section> &profile,
                    const Reference reference)
{
    std::ofstream output{directory / "profiles.csv"}, exact{directory / "quasi1d.csv"};
    if (!output || !exact)
        throw std::runtime_error("Cannot open nozzle profiles.");
    output << std::setprecision(17)
           << "x,area,Mach_CFD,Mach_quasi1D,p_CFD,p_quasi1D,T_CFD,T_quasi1D,rho_CFD,rho_quasi1D,ux_CFD,ux_quasi1D,mdot_"
              "CFD,mdot_quasi1D\n";
    exact << std::setprecision(17) << "x,area,Mach,p,T,rho,ux,mdot\n";
    for (cfd::Index i = 0; i < profile.size(); ++i)
    {
        const double x{static_cast<double>(i) * reference.length / static_cast<double>(profile.size() - 1)};
        const auto primitive{reference.primitive(x)};
        const auto &section{profile[i]};
        output << x << ',' << section.area << ',' << section.mach << ',' << reference.mach(x) << ',' << section.pressure
               << ',' << primitive.pressure << ',' << section.temperature << ','
               << primitive.pressure / (reference.gas_constant * primitive.density) << ',' << section.density << ','
               << primitive.density << ',' << section.velocity << ',' << primitive.velocity_x << ','
               << section.mass_flow << ',' << reference.mass_flow() << '\n';
        exact << x << ',' << reference.area(x) << ',' << reference.mach(x) << ',' << primitive.pressure << ','
              << primitive.pressure / (reference.gas_constant * primitive.density) << ',' << primitive.density << ','
              << primitive.velocity_x << ',' << reference.mass_flow() << '\n';
    }
    output.close();
    exact.close();
    if (!output || !exact)
        throw std::runtime_error("Cannot finish nozzle profiles.");
}

bool run_level(const cfd::Index nx, const Reference reference, const bool uniform, const std::filesystem::path &root)
{
    const cfd::Index nr{nx / 8};
    const auto directory{root / (reference.choked ? "choked" : "subsonic") / ("N" + std::to_string(nx))};
    std::filesystem::create_directories(directory);
    auto start{Clock::now()};
    auto raw{cfd::verification::make_euler_nozzle(nx, nr)};
    const double mesh_generation{seconds(start)};
    start = Clock::now();
    const auto built{cfd::build_mesh(std::move(raw))};
    const double mesh_build{seconds(start)};
    start = Clock::now();
    const auto initial{cfd::verification::nozzle_initial_field(nx, nr, reference, uniform)};
    const double initialization{seconds(start)};
    const cfd::CaloricallyPerfectGas gas{reference.gamma, reference.gas_constant};
    const auto conditions{reference.conditions()};
    start = Clock::now();
    cfd::AxisymmetricEulerSolver solver{built.mesh, initial, gas, conditions, {.collect_timings = true}};
    const double construction{seconds(start)};
    const double sound{std::sqrt(reference.gamma * reference.gas_constant * reference.total_temperature)};
    cfd::EulerSteadyOptions options{.reference_density = reference.total_pressure /
                                                         (reference.gas_constant * reference.total_temperature),
                                    .reference_velocity = sound,
                                    .reference_length = reference.length,
                                    .check_interval = 50,
                                    .minimum_steps = 100,
                                    .maximum_steps = 200000,
                                    .maximum_time = 80.0 * reference.length / sound};
    cfd::EulerSteadyMonitor monitor{solver, options};
    cfd::EulerSteadyCsvWriter history{directory / "convergence.csv"};
    OutputFields fields{built.mesh.cell_count()};
    std::ofstream pvd{directory / "nozzle.pvd"};
    if (!pvd)
        throw std::runtime_error("Cannot open nozzle PVD.");
    pvd << std::setprecision(17)
        << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" "
           "byte_order=\"LittleEndian\">\n<Collection>\n";
    cfd::Index snapshot_count{}, next_output{};
    double previous_snapshot_time{-1.0}, export_time{}, monitor_time{}, stepping_time{};
    cfd::verification::EulerBoundaryWaveBalance balance{.initial = solver.diagnostics().integrals};
    const auto write = [&] {
        const double time{solver.diagnostics().time};
        if (time == previous_snapshot_time)
            return;
        const auto begin{Clock::now()};
        const auto name{"nozzle_" + std::to_string(snapshot_count++) + ".vtu"};
        snapshot(directory / name, built.mesh, solver.state(), gas, fields);
        pvd << "<DataSet timestep=\"" << time << "\" group=\"\" part=\"0\" file=\"" << name << "\"/>\n";
        previous_snapshot_time = time;
        export_time += seconds(begin);
    };
    write();
    auto observation{monitor.observe()};
    history.write(observation);
    const auto campaign_start{Clock::now()};
    bool near_steady_written{};
    while (observation.status == cfd::EulerSteadyStatus::Running)
    {
        start = Clock::now();
        const auto &step{solver.step(options.maximum_time)};
        stepping_time += seconds(start);
        balance.add(step);
        if (next_output < output_times.size() && step.time * sound / reference.length >= output_times.at(next_output))
        {
            write();
            ++next_output;
        }
        if (!monitor.needs_observation())
            continue;
        start = Clock::now();
        observation = monitor.observe();
        monitor_time += seconds(start);
        history.write(observation);
        const double maximum_residual{
            *std::max_element(observation.normalized_residual.begin(), observation.normalized_residual.end())};
        if (maximum_residual < 1.0e-5 && !near_steady_written)
        {
            write();
            near_steady_written = true;
        }
        if (step.accepted_steps % 2000 == 0)
            std::cout << "PROGRESS " << (reference.choked ? "choked" : "subsonic") << " N=" << nx
                      << " step=" << step.accepted_steps << " t=" << step.time << " Rmax=" << maximum_residual
                      << " imbalance=" << observation.mass_flow_relative_imbalance << '\n'
                      << std::flush;
    }
    const double elapsed{seconds(campaign_start)};
    write();
    pvd << "</Collection>\n</VTKFile>\n";
    pvd.close();
    if (!pvd)
        throw std::runtime_error("Cannot finish nozzle PVD.");
    const auto profile{sections(built.mesh, solver.state(), solver.geometry(), conditions, gas, nx)};
    write_profiles(directory, profile, reference);
    std::array<long double, 5> squared{};
    double minimum_flow{profile.front().mass_flow}, maximum_flow{minimum_flow};
    for (cfd::Index i = 0; i < profile.size(); ++i)
    {
        const double x{static_cast<double>(i) / static_cast<double>(nx)};
        const auto exact{reference.primitive(x)};
        const auto section{profile[i]};
        const std::array error{section.mach - reference.mach(x),
                               section.pressure / reference.total_pressure - exact.pressure / reference.total_pressure,
                               section.temperature / reference.total_temperature -
                                   exact.pressure /
                                       (exact.density * reference.gas_constant * reference.total_temperature),
                               section.density / options.reference_density - exact.density / options.reference_density,
                               (section.velocity - exact.velocity_x) / sound};
        const long double weight{i == 0 || i == nx ? 0.5L : 1.0L};
        for (cfd::Index k = 0; k < squared.size(); ++k)
            squared.at(k) += weight * error.at(k) * error.at(k) / static_cast<long double>(nx);
        minimum_flow = std::min(minimum_flow, section.mass_flow);
        maximum_flow = std::max(maximum_flow, section.mass_flow);
    }
    // Summing the assembled cell mass residuals between two planes gives
    // their numerical flux difference (walls are impermeable). Its magnitude
    // is bounded by sum(abs(R_mass)); only summation roundoff is added here.
    const double flux_roundoff{128.0 * std::numeric_limits<double>::epsilon() * static_cast<double>(nx + 1) *
                               std::max(std::abs(minimum_flow), std::abs(maximum_flow))};
    if (!std::isfinite(maximum_flow - minimum_flow) ||
        maximum_flow - minimum_flow > observation.current.absolute_residual[0] + flux_roundoff)
        throw std::runtime_error("Nozzle section fluxes violate the assembled continuity residual bound.");
    std::cout << "RESULT " << (reference.choked ? "choked" : "subsonic") << " N=" << nx << " nr=" << nr
              << " cells=" << built.mesh.cell_count() << " status=" << static_cast<int>(observation.status)
              << " steps=" << observation.temporal.accepted_steps << " rejects=" << observation.temporal.rejected_steps
              << " time=" << observation.temporal.time << " CFL=" << observation.temporal.cfl << " residuals=";
    for (const double residual : observation.normalized_residual)
        std::cout << residual << ',';
    std::cout << " in=" << observation.current.mass_flow_in << " out=" << observation.current.mass_flow_out
              << " imbalance=" << observation.mass_flow_relative_imbalance
              << " max_step_balance=" << balance.maximum_step_error << " throat_M=" << profile[nx / 2].mach
              << " exit_M=" << profile.back().mach << " mdot_ref=" << reference.mass_flow()
              << " mdot_rel_error=" << (observation.current.mass_flow_out / reference.mass_flow() - 1.0)
              << " sectional_flow_spread=" << (maximum_flow - minimum_flow) / reference.mass_flow()
              << " RMS_M_p_T_rho_u=";
    for (const auto value : squared)
        std::cout << std::sqrt(value) << ',';
    std::cout << " elapsed=" << elapsed << " stepping=" << stepping_time << " monitoring=" << monitor_time
              << " export=" << export_time << " mesh_generation=" << mesh_generation << " mesh_build=" << mesh_build
              << " initialization=" << initialization << " construction=" << construction << '\n'
              << std::flush;
    std::cout << "TIMING assembly_including_BC=" << solver.timings().assembly_seconds
              << " boundary_subset=" << solver.timings().boundary_flux_seconds
              << " update=" << solver.timings().time_update_seconds << "\n"
              << std::flush;
    std::cout << "GLOBAL_BALANCE mass_mx_mr_energy=";
    for (const double defect : balance.defect(solver.diagnostics().integrals))
    {
        if (!std::isfinite(defect))
            throw std::runtime_error("Nozzle global balance is non-finite.");
        std::cout << defect << ',';
    }
    std::cout << '\n' << std::flush;
    return observation.status == cfd::EulerSteadyStatus::Converged;
}
} // namespace

int main(int argc, char **argv)
{
    try
    {
        const std::span<char *> arguments{argv, static_cast<std::size_t>(argc)};
        if (argc < 2 || argc > 5)
            throw std::invalid_argument("Usage: cfd_axisymmetric_euler_nozzle output_directory [subsonic|choked|all] "
                                        "[quasi1d|uniform] [single_nx]");
        const std::filesystem::path directory{arguments[1]};
        const std::string mode{argc >= 3 ? arguments[2] : "all"}, initialization{argc >= 4 ? arguments[3] : "quasi1d"};
        if ((mode != "all" && mode != "subsonic" && mode != "choked") ||
            (initialization != "uniform" && initialization != "quasi1d"))
            throw std::invalid_argument("Unknown nozzle case or initialization.");
        if (initialization == "uniform" && mode != "subsonic")
            throw std::invalid_argument(
                "Uniform startup supports only the subsonic case; select quasi1d for choked/all.");
        cfd::Index single{};
        if (argc == 5)
        {
            const std::string_view number{arguments[4]};
            const auto parsed{std::from_chars(number.begin(), number.end(), single)};
            if (parsed.ec != std::errc{} || parsed.ptr != number.end() ||
                std::find(levels.begin(), levels.end(), single) == levels.end())
                throw std::invalid_argument("single_nx must be 24, 48 or 96.");
        }
        std::cout << std::scientific << std::setprecision(10);
        for (const bool choked : {false, true})
        {
            if ((mode == "subsonic" && choked) || (mode == "choked" && !choked))
                continue;
            for (const cfd::Index nx : levels)
            {
                if (single != 0 && single != nx)
                    continue;
                if (!run_level(nx, {choked}, initialization == "uniform", directory))
                    throw std::runtime_error(
                        "Nozzle did not satisfy stationarity: inspect its history and final profiles.");
            }
        }
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
