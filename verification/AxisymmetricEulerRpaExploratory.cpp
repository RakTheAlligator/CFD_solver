// Exploratory RPA nozzle + plume Euler calculation, NOT a validated space-plume model.
#include "cfd/compressible/AxisymmetricEulerSolver.hpp"
#include "cfd/compressible/EulerBoundaryConditions.hpp"
#include "cfd/compressible/EulerCharacteristicBoundaryFlux.hpp"
#include "cfd/compressible/EulerFarfieldBoundaryFlux.hpp"
#include "cfd/compressible/EulerSteadyMonitor.hpp"
#include "cfd/input/OpenFOAMCaseReader.hpp"
#include "cfd/io/EulerCheckpoint.hpp"
#include "cfd/io/EulerSteadyCsvWriter.hpp"
#include "cfd/io/VtkWriter.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/meshing/GmshMesher.hpp"
#include "cfd/meshing/RpaNozzleGeometry.hpp"
#include "cfd/meshing/RpaRegionLayout.hpp"
#include "support/EulerBoundaryWave.hpp"
#include "support/EulerRpaExploratory.hpp"
#include "support/EulerRpaTemporal.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace fs = std::filesystem;

namespace
{

using Model = cfd::verification::EulerRpaModel;
using Clock = std::chrono::steady_clock;
constexpr double gamma_hot{Model::gamma}, gas_constant_hot{Model::gas_constant};
constexpr double p_hot{Model::pressure}, t_hot{Model::temperature}, mach_hot{Model::mach};
constexpr double ambient_pressure{Model::ambient_pressure};

using BCType = cfd::EulerBoundaryConditionType;

struct Options
{
    fs::path case_path;
    fs::path output_path;
    double ambient_temperature{}; // Explicit scenario assumption: NOT measured RPA data.
    double ambient_ux{};          // Explicit scenario assumption.
    double cfl{};
    unsigned long steps{};
    unsigned long save_every{};
    double ambient_ur{}, maximum_time{1.0}, maximum_wall_seconds{600.0};
    unsigned long diagnostic_every{25}, maximum_rejected{20};
    bool physical_exports{};
    std::optional<double> comparison_time{};
    std::optional<fs::path> restart{};
    unsigned long checkpoint_every{}, stop_after{};
};

[[nodiscard]] double real_number(const char *input, const char *name)
{
    const std::string token(input);
    std::size_t consumed{};
    const double value = std::stod(token, &consumed);
    if (consumed != token.size() || !std::isfinite(value))
        throw std::invalid_argument(std::string(name) + " must be a finite number.");
    return value;
}

[[nodiscard]] unsigned long positive_integer(const char *input, const char *name)
{
    const std::string_view token{input};
    unsigned long result{};
    const auto parsed{std::from_chars(token.data(), token.data() + token.size(), result)};
    if (token.empty() || parsed.ec != std::errc{} || parsed.ptr != token.data() + token.size() || result == 0)
        throw std::invalid_argument(std::string(name) + " must be a positive integer.");
    return result;
}

[[nodiscard]] Options parse_options(int argc, char **argv)
{
    const int full_count{argc};
    const std::span<char *> all_args{argv, static_cast<cfd::Index>(full_count)};
    for (int k = 8; k < full_count; ++k)
        if (std::string_view{all_args[k]} == "--restart" || std::string_view{all_args[k]} == "--checkpoint-every" ||
            std::string_view{all_args[k]} == "--stop-after")
        {
            argc = k;
            break;
        }
    if (argc < 8 || argc > 15)
    {
        throw std::invalid_argument("Usage: cfd_axisymmetric_euler_rpa_exploratory CASE_DIR OUTPUT_DIR "
                                    "AMBIENT_T_K AMBIENT_UX_M_PER_S CFL ACCEPTED_STEPS SAVE_EVERY "
                                    "[AMBIENT_UR [MAX_TIME [DIAG_EVERY [MAX_REJECTED [MAX_WALL_SECONDS "
                                    "[--physical-exports [COMPARISON_TIME]]]]]]] "
                                    "[--restart FILE] [--checkpoint-every STEPS] [--stop-after SEGMENT_STEPS]");
    }
    const std::span<char *> args{argv, static_cast<cfd::Index>(argc)};
    Options o{args[1],
              args[2],
              real_number(args[3], "AMBIENT_T_K"),
              real_number(args[4], "AMBIENT_UX_M_PER_S"),
              real_number(args[5], "CFL"),
              positive_integer(args[6], "ACCEPTED_STEPS"),
              positive_integer(args[7], "SAVE_EVERY")};
    if (argc > 8)
        o.ambient_ur = real_number(args[8], "AMBIENT_UR");
    if (argc > 9)
        o.maximum_time = real_number(args[9], "MAX_TIME");
    if (argc > 10)
        o.diagnostic_every = positive_integer(args[10], "DIAG_EVERY");
    if (argc > 11)
        o.maximum_rejected = positive_integer(args[11], "MAX_REJECTED");
    if (argc > 12)
        o.maximum_wall_seconds = real_number(args[12], "MAX_WALL_SECONDS");
    if (argc > 13)
    {
        if (std::string_view{args[13]} != "--physical-exports")
            throw std::invalid_argument("The last optional argument must be --physical-exports.");
        o.physical_exports = true;
        if (argc > 14)
            o.comparison_time = real_number(args[14], "COMPARISON_TIME");
        (void)cfd::verification::RpaPhysicalOutputs{o.maximum_time, o.comparison_time};
    }
    for (int k = argc; k < full_count; k += 2)
    {
        if (k + 1 >= full_count)
            throw std::invalid_argument("Missing restart/checkpoint option value.");
        const std::string_view flag{all_args[k]};
        if (flag == "--restart" && !o.restart)
            o.restart = fs::path{all_args[k + 1]};
        else if (flag == "--checkpoint-every" && o.checkpoint_every == 0)
            o.checkpoint_every = positive_integer(all_args[k + 1], "checkpoint interval");
        else if (flag == "--stop-after" && o.stop_after == 0)
            o.stop_after = positive_integer(all_args[k + 1], "segment steps");
        else
            throw std::invalid_argument("Unknown or duplicate restart/checkpoint option.");
    }
    if (!(o.maximum_time > 0) || !(o.maximum_wall_seconds > 0) || o.steps > 1000000 ||
        (!o.physical_exports && o.steps / o.save_every > 1000))
        throw std::invalid_argument("Positive time limits, at most 1e6 steps and at most 1002 snapshots required.");
    if (!(o.ambient_temperature > 0.0) || !(o.cfl > 0.0 && o.cfl <= 1.0))
        throw std::invalid_argument("AMBIENT_T_K > 0 and 0 < CFL <= 1 are required.");
    if (!fs::exists(o.case_path / "system" / "meshDict"))
        throw std::runtime_error("Missing RPA system/meshDict; case directory does not exist.");
    cfd::verification::require_external_empty_output(o.case_path, o.output_path);
    return o;
}

struct OutputBuffers
{
    explicit OutputBuffers(cfd::Index count)
        : p(count), T(count), ux(count), ur(count), speed(count), mach(count), log10p(count), log10rho(count),
          region(count)
    {
    }

    std::vector<double> p, T, ux, ur, speed, mach, log10p, log10rho, region;
};

void write_snapshot(const fs::path &path, const cfd::Mesh &mesh, const cfd::ConservativeField &field,
                    const cfd::CaloricallyPerfectGas &gas, OutputBuffers &b)
{
    for (cfd::Index cell = 0; cell < mesh.cell_count(); ++cell)
    {
        const auto evaluated = gas.evaluate(field.state(cell));
        b.p[cell] = evaluated.primitive.pressure;
        b.T[cell] = evaluated.thermodynamics.temperature;
        b.ux[cell] = evaluated.primitive.velocity_x;
        b.ur[cell] = evaluated.primitive.velocity_r;
        b.speed[cell] = std::hypot(b.ux[cell], b.ur[cell]);
        b.mach[cell] = b.speed[cell] / evaluated.thermodynamics.sound_speed;
        b.log10p[cell] = std::log10(b.p[cell]);
        b.log10rho[cell] = std::log10(evaluated.primitive.density);
    }
    const std::array scalars{
        cfd::VtkCellScalarData{"rho", field.density()},
        cfd::VtkCellScalarData{"p", b.p},
        cfd::VtkCellScalarData{"T", b.T},
        cfd::VtkCellScalarData{"ux", b.ux},
        cfd::VtkCellScalarData{"ur", b.ur},
        cfd::VtkCellScalarData{"velocity_magnitude", b.speed},
        cfd::VtkCellScalarData{"Mach", b.mach},
        cfd::VtkCellScalarData{"rhoE", field.total_energy_density()},
        cfd::VtkCellScalarData{"log10_p_Pa", b.log10p},
        cfd::VtkCellScalarData{"log10_rho", b.log10rho},
        cfd::VtkCellScalarData{"rpa_region", b.region},
    };
    const std::array vectors{cfd::VtkCellVectorComponentData{"U", b.ux, b.ur}};
    cfd::write_vtu(mesh, path, {.scalars = scalars, .component_vectors = vectors});
}

struct Snapshot
{
    double time;
    std::string filename;
};

void write_pvd(const fs::path &path, const std::vector<Snapshot> &snapshots)
{
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("Cannot write PVD: " + path.string());
    out << std::setprecision(17)
        << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" "
           "byte_order=\"LittleEndian\">\n<Collection>\n";
    for (const auto &s : snapshots)
        out << "<DataSet timestep=\"" << s.time << "\" group=\"\" part=\"0\" file=\"" << s.filename << "\"/>\n";
    out << "</Collection>\n</VTKFile>\n";
    out.flush();
    if (!out)
        throw std::runtime_error("PVD write failed.");
}

struct GroupStats
{
    std::array<cfd::Index, 4> regimes{};
    double mass{}, area{}, pressure{}, min_mach{std::numeric_limits<double>::infinity()},
        max_mach{-std::numeric_limits<double>::infinity()};
    double min_pressure{std::numeric_limits<double>::infinity()}, max_pressure{};
};

double write_boundary_diagnostics(const cfd::Mesh &mesh, const cfd::AxisymmetricEulerSolver &solver,
                                  const cfd::EulerBoundaryConditions &bcs, const cfd::CaloricallyPerfectGas &gas,
                                  std::ostream &csv)
{
    std::array<GroupStats, 8> groups{};
    double farfield_seconds{};
    for (cfd::Index face = 0; face < mesh.face_count(); ++face)
    {
        const auto adjacency{mesh.face_adjacencies()[face]};
        if (!adjacency.is_boundary())
            continue;
        const auto id{mesh.face_boundary_ids()[face]};
        const auto &bc{bcs.condition(id)};
        if (bc.type() != BCType::Farfield && bc.type() != BCType::TotalConditionsInlet)
            continue;
        const auto sf{mesh.face_area_vectors()[face]}, sa{solver.geometry().face_area_vectors()[face]};
        const double length{mesh.face_lengths()[face]}, area{std::hypot(sa.x, sa.y)};
        const cfd::Vector2 n{sf.x / length, sf.y / length};
        cfd::ConservativeState trace;
        cfd::EulerFlux flux;
        double mach{};
        cfd::Index regime{};
        try
        {
            if (bc.type() == BCType::Farfield)
            {
                const auto start{Clock::now()};
                const auto r{
                    cfd::compute_euler_farfield_boundary_flux(solver.state().state(adjacency.owner), n, bc, gas)};
                farfield_seconds += std::chrono::duration<double>(Clock::now() - start).count();
                trace = r.boundary_state;
                flux = r.flux;
                mach = r.normal_mach;
                regime = static_cast<cfd::Index>(r.regime);
            }
            else
            {
                const auto r{
                    cfd::compute_euler_characteristic_boundary_flux(solver.state().state(adjacency.owner), n, bc, gas)};
                trace = r.boundary_state;
                flux = r.flux;
                mach = r.normal_mach;
                regime = mach <= -1 ? 0 : mach < 0 ? 1 : mach < 1 ? 2 : 3;
            }
        }
        catch (const std::exception &error)
        {
            throw std::runtime_error("Boundary face=" + std::to_string(face) +
                                     " owner=" + std::to_string(adjacency.owner) +
                                     " group=" + mesh.boundary_groups()[id].name + ": " + error.what());
        }
        auto &g{groups.at(id)};
        ++g.regimes.at(regime);
        g.mass += area * flux.mass;
        g.area += area;
        const double pressure{gas.evaluate(trace).primitive.pressure};
        g.pressure += area * pressure;
        g.min_pressure = std::min(g.min_pressure, pressure);
        g.max_pressure = std::max(g.max_pressure, pressure);
        g.min_mach = std::min(g.min_mach, mach);
        g.max_mach = std::max(g.max_mach, mach);
    }
    for (cfd::Index id = 0; id < groups.size(); ++id)
    {
        const auto &g{groups.at(id)};
        if (g.area == 0)
            continue;
        csv << solver.diagnostics().accepted_steps << ',' << solver.diagnostics().time << ','
            << mesh.boundary_groups()[id].name << ',' << g.mass;
        for (const auto count : g.regimes)
            csv << ',' << count;
        csv << ',' << g.min_mach << ',' << g.max_mach << ',' << g.pressure / g.area << ',' << g.min_pressure << ','
            << g.max_pressure << '\n';
    }
    csv.flush();
    if (!csv)
        throw std::runtime_error("Cannot write boundary diagnostics.");
    return farfield_seconds;
}

int run(const Options &o)
{
    const cfd::CaloricallyPerfectGas gas{gamma_hot, gas_constant_hot};
    const auto hot{Model::hot()}, ambient{Model::ambient(o.ambient_temperature, o.ambient_ux, o.ambient_ur)};
    // Output ownership is acquired before any file is written; existing runs are never cleaned.
    cfd::verification::require_external_empty_output(o.case_path, o.output_path);
    fs::create_directories(o.output_path);
    std::ofstream metadata{o.output_path / "assumptions.txt"};
    metadata << std::setprecision(17) << "EXPLORATORY; NOT A VALIDATED RAREFIED SPACE PLUME\n"
             << "case=" << fs::absolute(o.case_path) << "\ngamma=" << gamma_hot << "\nR=" << gas_constant_hot
             << "\np_inlet_static=" << p_hot << "\nT_inlet_static=" << t_hot << "\nMach_inlet=" << mach_hot
             << "\nP0_derived=" << Model::total_pressure() << "\nT0_derived=" << Model::total_temperature()
             << "\np_ambient_nominal=" << ambient_pressure << "\nT_ambient_ASSUMED=" << o.ambient_temperature
             << "\nux_ambient_ASSUMED=" << o.ambient_ux << "\nur_ambient_ASSUMED=" << o.ambient_ur
             << "\nambient_same_hot_gas_EOS_ASSUMED=true\ninitial_engine_hot_exterior_ambient=true\nCFL=" << o.cfl
             << "\nstep_cap=" << o.steps << "\nmaximum_time=" << o.maximum_time
             << "\nmaximum_wall_seconds=" << o.maximum_wall_seconds << "\nmaximum_rejected=" << o.maximum_rejected
             << "\nsave_every=" << o.save_every << "\ndiag_every=" << o.diagnostic_every
             << "\nrestart=" << (o.restart ? "true" : "false") << "\ncheckpoint_every=" << o.checkpoint_every << '\n';
    if (o.restart)
        metadata << "restart_from=" << fs::absolute(*o.restart) << '\n';
    metadata << "physical_exports=" << o.physical_exports << '\n';
    metadata.flush();
    if (!metadata)
        throw std::runtime_error("Cannot write assumptions.");
    const auto elapsed = [](Clock::time_point begin) {
        return std::chrono::duration<double>(Clock::now() - begin).count();
    };
    const auto input{cfd::input::read_mesh_dict(o.case_path / "system" / "meshDict")};
    const auto *rpa{std::get_if<cfd::RpaNozzleGeometry>(&input.geometry)};
    if (rpa == nullptr)
        throw std::invalid_argument("meshDict must specify rpaNozzle.");
    auto start{Clock::now()};
    auto raw{cfd::generate_mesh(input.geometry, input.generation_options, input.automatic_meshing,
                                input.backward_facing_step_meshing, input.rpa_nozzle_meshing)};
    const double generation_seconds{elapsed(start)};
    start = Clock::now();
    auto built{cfd::build_mesh(std::move(raw))};
    const double mesh_seconds{elapsed(start)};
    const auto &mesh{built.mesh};
    cfd::verification::require_connected_rpa(mesh);
    const auto domain{cfd::read_rpa_nozzle_domain(*rpa)};
    const cfd::RpaRegionLayout layout{domain, input.rpa_nozzle_meshing.regions};
    const auto bcs{cfd::verification::rpa_euler_conditions(mesh.boundary_groups(), ambient)};
    start = Clock::now();
    cfd::ConservativeField initial{mesh.cell_count()};
    OutputBuffers buffers{mesh.cell_count()};
    const auto hot_cells{cfd::verification::initialize_rpa(mesh, layout, ambient, initial, buffers.region)};
    const double initialization_seconds{elapsed(start)};
    const double length{domain.inner_wall.back().x - domain.inner_wall.front().x};
    const double sound{std::sqrt(gamma_hot * gas_constant_hot * t_hot)};
    const double acoustic_time{cfd::verification::rpa_nozzle_acoustic_time(length)};
    const cfd::AxisymmetricEulerOptions solver_options{
        .cfl = o.cfl, .maximum_steps = o.steps, .collect_timings = true, .maximum_rejected_steps = o.maximum_rejected};
    cfd::EulerSteadyOptions steady;
    steady.reference_density = hot.density;
    steady.reference_velocity = sound;
    steady.reference_length = length;
    steady.check_interval = o.diagnostic_every;
    steady.minimum_steps = std::min<unsigned long>(100, o.steps);
    steady.maximum_steps = o.steps;
    steady.maximum_time = o.maximum_time;
    std::ostringstream context;
    context << std::setprecision(17) << "RPA frozen equivalent gas; exploratory ambient T/u; hot uniform Engine; "
            << p_hot << ' ' << t_hot << ' ' << mach_hot << ' ' << o.ambient_temperature << ' ' << o.ambient_ux << ' '
            << o.ambient_ur << '\n';
    for (const double region : buffers.region)
        context << region << ' ';
    const auto physical_context{context.str()};
    const cfd::EulerCheckpointConfiguration checkpoint_configuration{mesh,           gas,    bcs,
                                                                     solver_options, steady, physical_context};
    std::optional<cfd::EulerCheckpoint> restored;
    double checkpoint_read_seconds{}, checkpoint_write_seconds{};
    if (o.restart)
    {
        const auto begin{Clock::now()};
        restored.emplace(cfd::read_euler_checkpoint(*o.restart, checkpoint_configuration));
        checkpoint_read_seconds = elapsed(begin);
    }
    const auto *const restart_state{restored.has_value() ? &restored.value() : nullptr};
    start = Clock::now();
    cfd::AxisymmetricEulerSolver solver =
        restart_state != nullptr ? cfd::AxisymmetricEulerSolver{mesh, restart_state->solver, gas, bcs, solver_options}
                                 : cfd::AxisymmetricEulerSolver{mesh, initial, gas, bcs, solver_options};
    const double construction_seconds{elapsed(start)};
    const auto segment_initial{solver.diagnostics()};
    metadata << "segment_initial_time=" << segment_initial.time
             << "\nsegment_initial_step=" << segment_initial.accepted_steps << '\n';
    const std::array sections{
        cfd::verification::RpaFaceSection{mesh, solver.geometry(), buffers.region, domain.inner_wall.front().x,
                                          domain.inner_wall.front().y},
        cfd::verification::RpaFaceSection{mesh, solver.geometry(), buffers.region, domain.throat.x, domain.throat.y},
        cfd::verification::RpaFaceSection{mesh, solver.geometry(), buffers.region, domain.inner_wall.back().x,
                                          domain.inner_wall.back().y}};
    constexpr std::array<std::string_view, 3> section_names{"inlet", "throat", "exit"};
    std::optional<cfd::verification::RpaPhysicalOutputs> physical_outputs;
    cfd::verification::RpaPhysicalOutputs *physical_schedule{};
    if (o.physical_exports)
    {
        physical_outputs.emplace(o.maximum_time, o.comparison_time);
        physical_schedule = &physical_outputs.value();
        if (restart_state != nullptr)
            physical_schedule->resume_after(segment_initial.time);
        for (const double time : physical_schedule->times())
            if (time > segment_initial.time)
                metadata << "requested_output_time=" << time << '\n';
    }
    cfd::EulerSteadyMonitor monitor = restart_state != nullptr
                                          ? cfd::EulerSteadyMonitor{solver, steady, restart_state->monitor}
                                          : cfd::EulerSteadyMonitor{solver, steady};
    cfd::EulerSteadyCsvWriter convergence{o.output_path / "convergence.csv"};
    std::ofstream history{o.output_path / "steps.csv"}, boundary{o.output_path / "boundary_groups.csv"},
        cells{o.output_path / "cell_diagnostics.csv"}, section_csv{o.output_path / "sections.csv"};
    history << std::setprecision(17)
            << "step,time,dt,CFL_upper_bound,minimum_density,minimum_pressure,maximum_density,maximum_pressure,"
               "mass_outflow_rate,rejected_steps,mass,momentum_x,momentum_r,energy,balance_error_mass,balance_error_mx,"
               "balance_error_mr,balance_error_energy,boundary_mx,boundary_mr,boundary_energy,radial_source\n";
    boundary << std::setprecision(17)
             << "step,time,group,outward_mass_flow,supersonic_in,subsonic_in,subsonic_out,supersonic_out,min_trace_"
                "mach,max_trace_mach,area_mean_trace_pressure,min_trace_pressure,max_trace_pressure\n";
    cells << std::setprecision(17)
          << "step,time,minimum_temperature,maximum_mach,maximum_temperature,engine_maximum_mach,time_acoustic,"
             "dt,CFL_upper_bound,rejected_steps,minimum_pressure,maximum_pressure,minimum_density,maximum_density,"
             "mass_flow_in,mass_flow_out,residual_max,elapsed_wall_s,engine_mass\n";
    section_csv << std::setprecision(17)
                << "step,time,time_acoustic,section,faces,area,mean_cell_Mach,mean_cell_axial_Mach,"
                   "minimum_cell_Mach,maximum_cell_Mach,mean_cell_pressure,mean_cell_temperature,numerical_mass_flow\n";
    if (!history || !boundary || !cells || !section_csv)
        throw std::runtime_error("Cannot open RPA histories.");
    metadata << "cells=" << mesh.cell_count() << "\nfaces=" << mesh.face_count() << "\nhot_cells=" << hot_cells
             << "\nnozzle_length=" << length << "\ndomain_length=" << domain.downstream_x - domain.inner_wall.front().x
             << "\ninlet_x=" << domain.inner_wall.front().x << "\ninlet_r=" << domain.inner_wall.front().y
             << "\nfarfield_radius=" << domain.farfield_radius << "\ndownstream_x=" << domain.downstream_x
             << "\nthroat_x=" << domain.throat.x << "\nthroat_r=" << domain.throat.y
             << "\nexit_x=" << domain.inner_wall.back().x << "\nexit_r=" << domain.inner_wall.back().y
             << "\nacoustic_time_nozzle=" << acoustic_time
             << "\nconvective_time_nozzle_inlet=" << length / hot.velocity_x << "\nminimum_shape_quality="
             << *std::min_element(mesh.cell_qualities().begin(), mesh.cell_qualities().end())
             << "\ngeneration_s=" << generation_seconds << "\nmesh_build_s=" << mesh_seconds
             << "\ninitialization_s=" << initialization_seconds << "\nsolver_construction_s=" << construction_seconds
             << '\n';
    metadata.flush();
    if (!metadata)
        throw std::runtime_error("Cannot finish assumptions.");
    std::cout << "MESH cells=" << mesh.cell_count() << " faces=" << mesh.face_count() << " hot_cells=" << hot_cells
              << '\n'
              << std::flush;
    std::vector<Snapshot> snapshots;
    snapshots.reserve(o.physical_exports ? 8 : 2 + o.steps / o.save_every);
    double export_seconds{}, monitoring_seconds{}, farfield_sample_seconds{}, section_seconds{};
    const auto save = [&] {
        const auto begin{Clock::now()};
        const auto name{"rpa_" + std::to_string(solver.diagnostics().accepted_steps) + ".vtu"};
        write_snapshot(o.output_path / name, mesh, solver.state(), gas, buffers);
        snapshots.push_back({solver.diagnostics().time, name});
        write_pvd(o.output_path / "rpa.pvd", snapshots);
        export_seconds += elapsed(begin);
    };
    cfd::verification::EulerBoundaryWaveBalance balance{.initial = solver.diagnostics().integrals};
    std::array<long double, 4> balance_scales{};
    cfd::EulerSteadyObservation observation;
    if (restart_state != nullptr)
    {
        balance.initial = restart_state->history.initial;
        balance.integrated_boundary = restart_state->history.integrated_boundary;
        balance.integrated_source = restart_state->history.integrated_source;
        balance.maximum_step_error = restart_state->history.maximum_step_error;
        balance_scales = restart_state->history.balance_scales;
        observation = restart_state->monitor.last;
        if (physical_schedule != nullptr &&
            std::ranges::equal(physical_schedule->times(), restart_state->history.output_times) &&
            physical_schedule->next_index() != restart_state->history.next_output)
            throw std::runtime_error("RPA checkpoint output cursor mismatch.");
    }
    std::optional<std::pair<cfd::Index, cfd::Index>> last_checkpoint;
    const auto checkpoint = [&] {
        const auto step{solver.diagnostics().accepted_steps};
        const auto identity{std::pair{step, solver.diagnostics().rejected_steps}};
        if (last_checkpoint == identity)
            return;
        cfd::EulerCheckpointHistory history_state{balance.initial,
                                                  balance.integrated_boundary,
                                                  balance_scales,
                                                  balance.integrated_source,
                                                  balance.maximum_step_error,
                                                  {},
                                                  0};
        if (physical_schedule != nullptr)
        {
            history_state.output_times.assign(physical_schedule->times().begin(), physical_schedule->times().end());
            history_state.next_output = physical_schedule->next_index();
            // A failed visualization must not make the accepted state
            // uncheckpointable. A resumed segment exports this state itself.
            while (history_state.next_output < history_state.output_times.size() &&
                   history_state.output_times[history_state.next_output] <= solver.diagnostics().time)
                ++history_state.next_output;
        }
        const cfd::EulerCheckpoint state{solver.restart_state(), monitor.restart_state(), std::move(history_state)};
        const auto begin{Clock::now()};
        const auto suffix{identity.second == 0 ? std::string{} : "_r" + std::to_string(identity.second)};
        cfd::write_euler_checkpoint(o.output_path / ("checkpoint_" + std::to_string(step) + suffix + ".chk"),
                                    checkpoint_configuration, state);
        checkpoint_write_seconds += elapsed(begin);
        last_checkpoint = identity;
    };
    std::string status{"STEP_LIMIT"}, phase{"initial diagnostics"}, failure;
    const auto observe = [&](const double wall_seconds, const bool refresh = true) {
        const auto begin{Clock::now()};
        if (refresh)
            observation = monitor.observe();
        convergence.write(observation);
        double tmin{std::numeric_limits<double>::infinity()}, tmax{}, mmax{}, engine_mmax{};
        long double engine_mass{};
        for (cfd::Index cell = 0; cell < mesh.cell_count(); ++cell)
        {
            const auto e{gas.evaluate(solver.state().state(cell))};
            tmin = std::min(tmin, e.thermodynamics.temperature);
            tmax = std::max(tmax, e.thermodynamics.temperature);
            const double mach{std::hypot(e.primitive.velocity_x, e.primitive.velocity_r) /
                              e.thermodynamics.sound_speed};
            mmax = std::max(mmax, mach);
            if (buffers.region[cell] == static_cast<double>(cfd::RpaRegion::Engine))
            {
                engine_mmax = std::max(engine_mmax, mach);
                engine_mass += static_cast<long double>(solver.geometry().cell_volumes()[cell]) * e.primitive.density;
            }
        }
        const auto &d{solver.diagnostics()};
        const double residual_max{*std::ranges::max_element(observation.normalized_residual)};
        cells << d.accepted_steps << ',' << d.time << ',' << tmin << ',' << mmax << ',' << tmax << ',' << engine_mmax
              << ',' << d.time / acoustic_time << ',' << d.time_step << ',' << d.cfl << ',' << d.rejected_steps << ','
              << d.minimum_pressure << ',' << d.maximum_pressure << ',' << d.minimum_density << ',' << d.maximum_density
              << ',' << observation.current.mass_flow_in << ',' << observation.current.mass_flow_out << ','
              << residual_max << ',' << wall_seconds << ',' << static_cast<double>(engine_mass) << '\n';
        cells.flush();
        if (!cells)
            throw std::runtime_error("Cannot write cell diagnostics.");
        farfield_sample_seconds += write_boundary_diagnostics(mesh, solver, bcs, gas, boundary);
        const auto section_start{Clock::now()};
        double throat_mach{};
        for (cfd::Index k = 0; k < sections.size(); ++k)
        {
            const auto sample{sections.at(k).sample(mesh, solver.state(), gas, bcs)};
            section_csv << d.accepted_steps << ',' << d.time << ',' << d.time / acoustic_time << ','
                        << section_names.at(k) << ',' << sections.at(k).face_count() << ',' << sample.area << ','
                        << sample.mach << ',' << sample.axial_mach << ',' << sample.minimum_mach << ','
                        << sample.maximum_mach << ',' << sample.pressure << ',' << sample.temperature << ','
                        << sample.mass_flow << '\n';
            if (k == 1)
                throat_mach = sample.mach;
        }
        section_csv.flush();
        if (!section_csv)
            throw std::runtime_error("Cannot write RPA section diagnostics.");
        section_seconds += elapsed(section_start);
        std::cout << "OBS step=" << d.accepted_steps << " t=" << d.time << " t/t_a=" << d.time / acoustic_time
                  << " dt=" << d.time_step << " Mach_max=" << mmax << " Mach_throat=" << throat_mach
                  << " residual_max=" << residual_max << " wall_s=" << wall_seconds << " rejects=" << d.rejected_steps
                  << " Engine_mass=" << static_cast<double>(engine_mass) << '\n'
                  << std::flush;
        monitoring_seconds += elapsed(begin);
    };
    save(); // Preserve the admissible initial state even if BC preflight fails.
    const auto loop_start{Clock::now()};
    try
    {
        if (restart_state == nullptr)
            observe(elapsed(loop_start));
        else
        {
            (void)solver.evaluate_steady_residual(); // Preflight without inventing a monitor observation.
            if (restart_state->monitor.observed &&
                observation.temporal.accepted_steps == segment_initial.accepted_steps)
                observe(elapsed(loop_start), false);
        }
        if (restart_state != nullptr && observation.status == cfd::EulerSteadyStatus::Converged)
            status = "CONVERGED";
        while (solver.diagnostics().accepted_steps < o.steps && solver.diagnostics().time < o.maximum_time)
        {
            if (status == "CONVERGED")
                break;
            if (o.stop_after > 0 &&
                solver.diagnostics().accepted_steps - segment_initial.accepted_steps >= o.stop_after)
            {
                status = "SEGMENT_STOP";
                break;
            }
            if (elapsed(loop_start) >= o.maximum_wall_seconds)
            {
                status = "WALL_LIMIT";
                break;
            }
            phase = "Euler update";
            const auto &d{solver.step(physical_schedule != nullptr ? physical_schedule->next_time() : o.maximum_time)};
            balance.add(d);
            const std::array rates{d.boundary_outflow_rate.mass, d.boundary_outflow_rate.momentum_x,
                                   d.boundary_outflow_rate.momentum_r, d.boundary_outflow_rate.energy};
            for (cfd::Index k = 0; k < 4; ++k)
                balance_scales.at(k) += static_cast<long double>(d.time_step) * std::abs(rates.at(k));
            balance_scales.at(2) += static_cast<long double>(d.time_step) * std::abs(d.radial_pressure_source_rate);
            history << d.accepted_steps << ',' << d.time << ',' << d.time_step << ',' << d.cfl << ','
                    << d.minimum_density << ',' << d.minimum_pressure << ',' << d.maximum_density << ','
                    << d.maximum_pressure << ',' << d.boundary_outflow_rate.mass << ',' << d.rejected_steps << ','
                    << d.integrals.mass << ',' << d.integrals.momentum_x << ',' << d.integrals.momentum_r << ','
                    << d.integrals.energy << ',' << d.balance_error.mass << ',' << d.balance_error.momentum_x << ','
                    << d.balance_error.momentum_r << ',' << d.balance_error.energy << ','
                    << d.boundary_outflow_rate.momentum_x << ',' << d.boundary_outflow_rate.momentum_r << ','
                    << d.boundary_outflow_rate.energy << ',' << d.radial_pressure_source_rate << '\n';
            history.flush();
            if (!history)
                throw std::runtime_error("Cannot write accepted-step history.");
            phase = "current-state diagnostics";
            const bool physical_export_due{physical_schedule != nullptr && physical_schedule->due(d.time)};
            if (monitor.needs_observation() || physical_export_due)
                observe(elapsed(loop_start));
            phase = "accepted-state export";
            if (physical_export_due || (physical_schedule == nullptr && d.accepted_steps % o.save_every == 0))
            {
                save();
                if (physical_schedule != nullptr && physical_export_due)
                    physical_schedule->exported(d.time);
            }
            if (observation.status == cfd::EulerSteadyStatus::Converged)
            {
                status = "CONVERGED";
                break;
            }
            if (d.rejected_steps >= o.maximum_rejected)
            {
                status = "REJECT_LIMIT";
                break;
            }
            if (o.checkpoint_every > 0 && d.accepted_steps % o.checkpoint_every == 0)
            {
                phase = "checkpoint publication";
                checkpoint();
            }
        }
        if (status == "STEP_LIMIT" && solver.diagnostics().time >= o.maximum_time)
            status = o.physical_exports ? "TARGET_REACHED" : "TIME_LIMIT";
        phase = "final accepted-state diagnostics";
        if (status != "SEGMENT_STOP" && status != "WALL_LIMIT" &&
            solver.diagnostics().accepted_steps > observation.temporal.accepted_steps)
            observe(elapsed(loop_start));
        if (status != "WALL_LIMIT" && status != "REJECT_LIMIT" &&
            observation.status == cfd::EulerSteadyStatus::Converged)
            status = "CONVERGED";
    }
    catch (const std::exception &error)
    {
        status = solver.diagnostics().rejected_steps >= o.maximum_rejected ? "REJECT_LIMIT" : "FAILED";
        failure = phase + ": " + error.what();
    }
    if (o.checkpoint_every > 0 || o.stop_after > 0 || o.restart)
        checkpoint();
    if (snapshots.back().time != solver.diagnostics().time)
        save(); // Only the last accepted state, never a rejected candidate.
    const auto defect{balance.defect(solver.diagnostics().integrals)};
    const auto first{balance.initial}, last{solver.diagnostics().integrals};
    const std::array initial_values{first.mass, first.momentum_x, first.momentum_r, first.energy};
    const std::array final_values{last.mass, last.momentum_x, last.momentum_r, last.energy};
    std::array<double, 4> relative_defects{};
    bool conserved{true};
    for (cfd::Index k = 0; k < 4; ++k)
    {
        const long double scale{balance_scales.at(k) + std::abs(initial_values.at(k)) + std::abs(final_values.at(k))};
        const double relative{scale > 0 ? static_cast<double>(std::abs(defect.at(k)) / scale)
                                        : (defect.at(k) == 0 ? 0.0 : std::numeric_limits<double>::infinity())};
        conserved = conserved && std::isfinite(relative) &&
                    relative <= 256.0 * std::numeric_limits<double>::epsilon() *
                                    (static_cast<double>(solver.diagnostics().accepted_steps) + 1.0);
        relative_defects.at(k) = relative;
    }
    if (!conserved)
    {
        status = "FAILED_CONSERVATION";
        failure += " Four-equation conservation gate failed.";
    }
    std::ofstream result{o.output_path / "status.txt"};
    result << std::setprecision(17) << "status=" << status << "\nreason=" << failure
           << "\nsteps=" << solver.diagnostics().accepted_steps << "\ntime_seconds=" << solver.diagnostics().time
           << "\nrejected_steps=" << solver.diagnostics().rejected_steps
           << "\nsteady_status=" << static_cast<int>(observation.status)
           << "\nlast_observation_step=" << observation.temporal.accepted_steps << "\ntarget_time=" << o.maximum_time
           << "\ntarget_reached=" << (solver.diagnostics().time >= o.maximum_time)
           << "\nacoustic_crossings_nozzle=" << solver.diagnostics().time * sound / length
           << "\nconvective_crossings_nozzle_inlet=" << solver.diagnostics().time * hot.velocity_x / length << '\n';
    for (cfd::Index k = 0; k < 4; ++k)
        result << "balance_absolute_" << k << '=' << defect.at(k) << "\nbalance_relative_" << k << '='
               << relative_defects.at(k) << '\n';
    const auto timing{solver.timings()};
    result << "conserved=" << conserved << "\nloop_s=" << elapsed(loop_start)
           << "\nassembly_s=" << timing.assembly_seconds << "\nall_BC_subset_s=" << timing.boundary_flux_seconds
           << "\nupdate_s=" << timing.time_update_seconds << "\nmonitoring_including_sample_BC_s=" << monitoring_seconds
           << "\nfarfield_diagnostic_subset_s=" << farfield_sample_seconds
           << "\nsection_diagnostics_subset_s=" << section_seconds << "\nexport_s=" << export_seconds
           << "\ncheckpoint_write_s=" << checkpoint_write_seconds << "\ncheckpoint_read_s=" << checkpoint_read_seconds
           << '\n';
    result.close();
    if (!result)
        throw std::runtime_error("Cannot write final status.");
    std::cout << "STATUS " << status << " steps=" << solver.diagnostics().accepted_steps
              << " time=" << solver.diagnostics().time << " conserved=" << conserved << " reason=" << failure
              << "\nPVD " << o.output_path / "rpa.pvd" << '\n';
    return status == "FAILED" || status == "REJECT_LIMIT" || !conserved ? 2 : 0;
}

} // namespace

int main(int argc, char **argv)
{
    std::optional<Options> options;
    try
    {
        options.emplace(parse_options(argc, argv));
        return run(*options);
    }
    catch (const std::exception &e)
    {
        if (options && fs::exists(options->output_path / "assumptions.txt") &&
            !fs::exists(options->output_path / "status.txt"))
        {
            std::ofstream status{options->output_path / "status.txt"};
            status << "status=SETUP_OR_OUTPUT_FAILURE\nreason=" << e.what() << '\n';
            status.close();
            if (!status)
                std::cerr << "Cannot write failure status.\n";
        }
        std::cerr << "RPA exploratory driver: " << e.what() << '\n';
        return 1;
    }
}
