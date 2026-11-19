#include "cfd/compressible/AxisymmetricEulerSteadySolver.hpp"
#include "cfd/input/OpenFOAMCaseReader.hpp"
#include "cfd/io/VtkWriter.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "support/EulerNozzle.hpp"
#include "support/EulerRpaExploratory.hpp"
#include "support/EulerRpaInitialization.hpp"
#include "support/EulerRpaTemporal.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <span>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

namespace
{
namespace fs = std::filesystem;
using Clock = std::chrono::steady_clock;
using Model = cfd::verification::EulerRpaModel;

struct Problem
{
    cfd::MeshBuildResult built;
    cfd::ConservativeField initial;
    cfd::CaloricallyPerfectGas gas;
    cfd::EulerBoundaryConditions conditions;
    std::vector<cfd::Index> regions;
    std::array<double, 3> section_x, section_radius;
    double rho_ref{}, velocity_ref{}, length_ref{};
};

Problem make_problem(const std::string_view name, const fs::path &rpa_case, const std::string_view initialization)
{
    if (name == "rpa")
    {
        const auto input{cfd::input::read_mesh_dict(rpa_case / "system/meshDict")};
        const auto domain{cfd::read_rpa_nozzle_domain(std::get<cfd::RpaNozzleGeometry>(input.geometry))};
        auto built{cfd::build_mesh(cfd::generate_mesh(input.geometry, input.generation_options, input.automatic_meshing,
                                                      input.backward_facing_step_meshing, input.rpa_nozzle_meshing))};
        cfd::verification::require_connected_rpa(built.mesh);
        const cfd::RpaRegionLayout layout{domain, input.rpa_nozzle_meshing.regions};
        const auto ambient{Model::ambient(300, 0, 0)}; // Same explicitly exploratory exterior as previous RPA runs.
        auto conditions{cfd::verification::rpa_euler_conditions(built.mesh.boundary_groups(), ambient)};
        cfd::ConservativeField field{built.mesh.cell_count()};
        std::vector<double> labels(field.size());
        (void)cfd::verification::initialize_rpa(built.mesh, layout, ambient, field, labels);
        std::vector<cfd::Index> regions(labels.size());
        std::ranges::transform(labels, regions.begin(), [](const double r) { return static_cast<cfd::Index>(r); });
        if (initialization == "quasi1d")
            (void)cfd::verification::initialize_rpa_quasi1d(
                built.mesh, regions, cfd::CaloricallyPerfectGas{Model::gamma, Model::gas_constant},
                Model::total_pressure(), Model::total_temperature(), field);
        return {std::move(built),
                std::move(field),
                {Model::gamma, Model::gas_constant},
                std::move(conditions),
                std::move(regions),
                {domain.inner_wall.front().x, domain.throat.x, domain.inner_wall.back().x},
                {domain.inner_wall.front().y, domain.throat.y, domain.inner_wall.back().y},
                Model::hot().density,
                std::sqrt(Model::gamma * Model::gas_constant * Model::temperature),
                domain.inner_wall.back().x - domain.inner_wall.front().x};
    }
    if (name == "nozzle")
    {
        constexpr cfd::Index nx{16}, nr{2};
        const cfd::verification::EulerNozzleReference reference{}; // Smooth, subsonic, same quasi-1D reference.
        auto built{cfd::build_mesh(cfd::verification::make_euler_nozzle(nx, nr))};
        auto initial{cfd::verification::nozzle_initial_field(nx, nr, reference, false)};
        std::vector<cfd::Index> regions(initial.size(), 1);
        return {std::move(built),
                std::move(initial),
                {reference.gamma, reference.gas_constant},
                reference.conditions(),
                std::move(regions),
                {0, .5, 1},
                {reference.radius(0), reference.radius(.5), reference.radius(1)},
                1,
                std::sqrt(reference.gamma),
                1};
    }
    if (name != "tube")
        throw std::invalid_argument("CASE must be tube, nozzle or rpa.");
    auto raw{cfd::verification::make_euler_tube(24, 2, cfd::verification::EulerTubeCells::Mixed)};
    for (auto &node : raw.nodes)
        node.x *= node.x; // Graded spacing: local/global share exactly this same mesh.
    auto built{cfd::build_mesh(std::move(raw))};
    const cfd::CaloricallyPerfectGas gas{1.4, 1};
    cfd::ConservativeField initial{built.mesh.cell_count()};
    for (cfd::Index cell = 0; cell < initial.size(); ++cell)
    {
        const double x{built.mesh.cell_centers()[cell].x};
        initial.set_state(cell, gas.to_conservative({1 + .05 * std::sin(std::numbers::pi * x), 2, 0, 1}));
    }
    using Type = cfd::EulerBoundaryConditionType;
    cfd::EulerBoundaryConditions conditions{4,
                                            {Type::Axis, Type::SlipWall,
                                             cfd::EulerBoundaryCondition::supersonic_inlet({1, 2, 0, 1}),
                                             Type::SupersonicOutlet}};
    std::vector<cfd::Index> regions(initial.size(), 1);
    // x=.25 is a mesh plane after squaring x=.5. It is a diagnostic section, not a throat.
    return {std::move(built), std::move(initial), gas, std::move(conditions), std::move(regions),
            {0, .25, 1},      {1, 1, 1},          1,   std::sqrt(1.4),        1};
}

double real(const char *token)
{
    std::size_t end{};
    const std::string value{token};
    const double result{std::stod(value, &end)};
    if (end != value.size() || !std::isfinite(result) || !(result > 0))
        throw std::invalid_argument("Expected a positive finite real.");
    return result;
}

cfd::Index count(const char *token)
{
    const std::string_view value{token};
    cfd::Index result{};
    const auto parsed{std::from_chars(value.data(), value.data() + value.size(), result)};
    if (parsed.ec != std::errc{} || parsed.ptr != value.data() + value.size() || result == 0)
        throw std::invalid_argument("Expected a positive iteration count.");
    return result;
}

class Output
{
  public:
    Output(const Problem &problem, const fs::path &directory, const std::string &mode)
        : problem_(&problem), directory_(directory), mode_(mode), geometry_(problem.built.mesh),
          regions_(problem.regions.begin(), problem.regions.end()), p_(problem.initial.size()), t_(p_.size()),
          ux_(p_.size()), ur_(p_.size()), mach_(p_.size()),
          sections_{cfd::verification::RpaFaceSection{problem.built.mesh, geometry_, regions_, problem.section_x[0],
                                                      problem.section_radius[0]},
                    cfd::verification::RpaFaceSection{problem.built.mesh, geometry_, regions_, problem.section_x[1],
                                                      problem.section_radius[1]},
                    cfd::verification::RpaFaceSection{problem.built.mesh, geometry_, regions_, problem.section_x[2],
                                                      problem.section_radius[2]}}
    {
        fs::create_directories(directory_);
        history_.open(directory_ / "iterations.csv");
        regional_.open(directory_ / "regions.csv");
        manifest_.open(directory_ / "snapshots.csv");
        for (auto *stream : {&history_, &regional_, &manifest_})
        {
            stream->exceptions(std::ios::badbit | std::ios::failbit);
            *stream << std::setprecision(17);
        }
        history_
            << "iteration,step_seconds,loop_cpu_seconds,residual_mass,residual_mx,residual_mr,residual_energy,max_cell_"
               "mass,max_cell_mx,"
               "max_cell_mr,max_cell_energy,balance_mass,balance_mx,balance_mr,balance_energy,rho_min,p_min,rejected,"
               "dtau_min,dtau_max,alpha,Qin,Qmid,Qout,Mmid,Mout,pmid,pout,Tout,source_r\n";
        regional_ << "iteration,region,volume,Rmass,Rmx,Rmr,Renergy,max_mass,max_mx,max_mr,max_energy\n";
        manifest_ << "iteration,file,kind\n"; // Not a physical-time PVD collection.
    }

    void observe(const cfd::Index iteration, const cfd::ConservativeField &state, const cfd::EulerLocalSteadyResult &d,
                 const std::span<const cfd::EulerRegionResidual> regions, const double step_seconds,
                 const double cpu_seconds)
    {
        if (last_observed_ && *last_observed_ == iteration)
            return;
        history_ << iteration << ',' << step_seconds << ',' << cpu_seconds;
        for (const auto values : {d.normalized_residual, d.maximum_cell_residual, d.stationary_balance})
            for (const double value : values)
                history_ << ',' << value;
        history_ << ',' << d.minimum_density << ',' << d.minimum_pressure << ',' << d.rejected_attempts << ','
                 << d.minimum_local_step << ',' << d.maximum_local_step << ',' << d.step_scale;
        std::array<cfd::verification::RpaSectionSample, 3> samples;
        for (cfd::Index k = 0; k < sections_.size(); ++k)
            samples.at(k) = sections_.at(k).sample(problem_->built.mesh, state, problem_->gas, problem_->conditions);
        history_ << ',' << samples[0].mass_flow << ',' << samples[1].mass_flow << ',' << samples[2].mass_flow << ','
                 << samples[1].mach << ',' << samples[2].mach << ',' << samples[1].pressure << ','
                 << samples[2].pressure << ',' << samples[2].temperature << ',' << d.radial_pressure_source_rate
                 << '\n';
        for (const auto &region : regions)
        {
            regional_ << iteration << ',' << region.region << ',' << region.volume;
            for (const auto values : {region.normalized_residual, region.maximum_cell_residual})
                for (const double value : values)
                    regional_ << ',' << value;
            regional_ << '\n';
        }
        last_observed_ = iteration;
    }

    void snapshot(const cfd::Index iteration, const cfd::ConservativeField &state)
    {
        if (last_export_ && *last_export_ == iteration)
            return;
        for (cfd::Index cell = 0; cell < state.size(); ++cell)
        {
            const auto value{problem_->gas.evaluate(state.state(cell))};
            p_[cell] = value.primitive.pressure;
            t_[cell] = value.thermodynamics.temperature;
            ux_[cell] = value.primitive.velocity_x;
            ur_[cell] = value.primitive.velocity_r;
            mach_[cell] = std::hypot(ux_[cell], ur_[cell]) / value.thermodynamics.sound_speed;
        }
        const auto filename{mode_ + "_iteration_" + std::to_string(iteration) + ".vtu"};
        const std::array scalars{cfd::VtkCellScalarData{"rho", state.density()},
                                 cfd::VtkCellScalarData{"rho_ux", state.momentum_x()},
                                 cfd::VtkCellScalarData{"rho_ur", state.momentum_r()},
                                 cfd::VtkCellScalarData{"rhoE", state.total_energy_density()},
                                 cfd::VtkCellScalarData{"p", p_},
                                 cfd::VtkCellScalarData{"T", t_},
                                 cfd::VtkCellScalarData{"ux", ux_},
                                 cfd::VtkCellScalarData{"ur", ur_},
                                 cfd::VtkCellScalarData{"Mach", mach_},
                                 cfd::VtkCellScalarData{"region", regions_}};
        cfd::write_vtu(problem_->built.mesh, directory_ / filename, {.scalars = scalars});
        manifest_ << iteration << ',' << filename << ',' << mode_ << "_iteration_not_physical_time\n";
        last_export_ = iteration;
    }

  private:
    const Problem *problem_;
    fs::path directory_;
    std::string mode_;
    cfd::AxisymmetricGeometry geometry_;
    std::vector<double> regions_, p_, t_, ux_, ur_, mach_;
    std::array<cfd::verification::RpaFaceSection, 3> sections_;
    std::ofstream history_, regional_, manifest_;
    std::optional<cfd::Index> last_export_;
    std::optional<cfd::Index> last_observed_;
};

int run(const Problem &problem, const std::string &mode, const fs::path &directory,
        cfd::EulerLocalSteadyOptions options, const std::string_view initialization = {})
{
    Output output{problem, directory, mode};
    output.snapshot(0, problem.initial);
    const auto start{Clock::now()};
    const auto cpu_start{std::clock()};
    double step_seconds{}, diagnostic_seconds{}, export_seconds{};
    const auto elapsed = [&] { return std::chrono::duration<double>(Clock::now() - start).count(); };
    const auto export_every{std::max(options.check_interval, options.maximum_iterations / 10)};
    cfd::Index iteration{};
    std::string status;
    std::string failure;
    cfd::EulerLocalSteadyResult final;
    const auto sample = [&](const cfd::Index i, const cfd::ConservativeField &state,
                            const cfd::EulerLocalSteadyResult &d,
                            const std::span<const cfd::EulerRegionResidual> regions) {
        const auto begin{Clock::now()};
        output.observe(i, state, d, regions, step_seconds,
                       static_cast<double>(std::clock() - cpu_start) / CLOCKS_PER_SEC);
        diagnostic_seconds += std::chrono::duration<double>(Clock::now() - begin).count();
    };
    const auto save = [&](const cfd::Index i, const cfd::ConservativeField &state) {
        const auto begin{Clock::now()};
        output.snapshot(i, state);
        export_seconds += std::chrono::duration<double>(Clock::now() - begin).count();
    };
    if (mode == "steady")
    {
        cfd::AxisymmetricEulerSteadySolver solver{problem.built.mesh, problem.initial, problem.gas,
                                                  problem.conditions, options,         problem.regions};
        solver.observe();
        sample(0, solver.state(), solver.diagnostics(), solver.regional_residuals());
        while (solver.diagnostics().status == cfd::EulerLocalSteadyStatus::Running)
        {
            const auto begin{Clock::now()};
            solver.iterate();
            step_seconds += std::chrono::duration<double>(Clock::now() - begin).count();
            iteration = solver.diagnostics().iteration;
            if (solver.diagnostics().observation_iteration == iteration)
                sample(iteration, solver.state(), solver.diagnostics(), solver.regional_residuals());
            if (iteration != 0 && iteration % export_every == 0)
                save(iteration, solver.state());
        }
        save(iteration, solver.state());
        final = solver.diagnostics();
        status = cfd::to_string(final.status);
        failure = solver.failure_message();
    }
    else
    {
        cfd::AxisymmetricEulerSolver solver{problem.built.mesh,
                                            problem.initial,
                                            problem.gas,
                                            problem.conditions,
                                            {.cfl = options.cfl, .maximum_steps = options.maximum_iterations}};
        // Comparison-only observer: reuse exactly the steady residual normalization
        // through the same assembler. Its setup is timed as diagnostics, not stepping.
        std::optional<cfd::Index> last_inspected;
        const auto inspect = [&] {
            if (last_inspected && *last_inspected == iteration)
                return;
            const auto begin{Clock::now()};
            cfd::AxisymmetricEulerSteadySolver observer{problem.built.mesh, solver.state(), problem.gas,
                                                        problem.conditions, options,        problem.regions};
            observer.observe();
            final = observer.diagnostics();
            final.observation_iteration = iteration;
            final.rejected_attempts = solver.diagnostics().rejected_steps;
            diagnostic_seconds += std::chrono::duration<double>(Clock::now() - begin).count();
            sample(iteration, solver.state(), final, observer.regional_residuals());
            last_inspected = iteration;
        };
        inspect();
        status = "MAX_ITERATIONS"; // A bounded physical reference, never a claim of steady convergence.
        try
        {
            while (iteration < options.maximum_iterations)
            {
                if (elapsed() >= options.maximum_wall_seconds)
                {
                    status = "WALL_LIMIT";
                    break;
                }
                const auto begin{Clock::now()};
                solver.step(std::numeric_limits<double>::max());
                step_seconds += std::chrono::duration<double>(Clock::now() - begin).count();
                iteration = solver.diagnostics().accepted_steps;
                if (iteration % options.check_interval == 0)
                    inspect();
                if (iteration % export_every == 0)
                    save(iteration, solver.state());
            }
        }
        catch (const std::runtime_error &error)
        {
            status = "NUMERICAL_FAILURE";
            failure = error.what();
        }
        save(iteration, solver.state());
        if (status != "NUMERICAL_FAILURE")
            inspect();
    }
    std::ofstream summary{directory / "summary.txt"};
    summary.exceptions(std::ios::badbit | std::ios::failbit);
    summary << std::setprecision(17) << "mode=" << mode << "\nstatus=" << status << "\nfailure=" << failure
            << "\ncells=" << problem.initial.size() << "\nfaces=" << problem.built.mesh.face_count()
            << "\niterations=" << iteration << "\nrejected=" << final.rejected_attempts
            << "\nobservation_iteration=" << final.observation_iteration << "\nstep_seconds=" << step_seconds
            << "\ndiagnostic_seconds=" << diagnostic_seconds << "\nexport_seconds=" << export_seconds
            << "\nloop_wall_seconds=" << elapsed()
            << "\nloop_cpu_seconds=" << static_cast<double>(std::clock() - cpu_start) / CLOCKS_PER_SEC
            << "\ncfl=" << options.cfl << "\nrho_ref=" << options.reference_density
            << "\nvelocity_ref=" << options.reference_velocity << "\nlength_ref=" << options.reference_length << '\n';
    if (!initialization.empty())
        summary << "initialization=" << initialization << '\n';
    std::cout << mode << ' ' << status << " iterations=" << iteration << " step_seconds=" << step_seconds << '\n';
    if (status == "NUMERICAL_FAILURE")
        std::cerr << failure << '\n';
    return status == "NUMERICAL_FAILURE" || status == "DIVERGED" ? 2 : 0;
}
} // namespace

int main(int argc, char **argv)
{
    try
    {
        if (argc < 8 || argc > 10)
            throw std::invalid_argument("Usage: cfd_axisymmetric_euler_steady MODE CASE NEW_OUTPUT CFL MAX_ITER "
                                        "CHECK_INTERVAL WALL_SECONDS [RPA_CASE [uniform|quasi1d]]\n"
                                        "MODE=steady|transient, CASE=tube|nozzle|rpa");
        const std::span<char *> args{argv, static_cast<cfd::Index>(argc)};
        const std::string mode{args[1]}, name{args[2]};
        if (mode != "steady" && mode != "transient")
            throw std::invalid_argument("MODE must be steady or transient.");
        if ((name == "rpa") != (argc >= 9))
            throw std::invalid_argument("RPA requires exactly one explicit case directory.");
        const std::string_view initialization{argc == 10 ? args[9] : "uniform"};
        if (initialization != "uniform" && initialization != "quasi1d")
            throw std::invalid_argument("RPA initialization must be uniform or quasi1d.");
        const fs::path output{args[3]}, input{argc >= 9 ? fs::path{args[8]} : fs::path{"."}};
        if (name == "rpa")
            cfd::verification::require_external_empty_output(input, output);
        else if (fs::exists(output) && (!fs::is_directory(output) || !fs::is_empty(output)))
            throw std::invalid_argument("Output directory must be empty; existing results are never overwritten.");
        cfd::EulerLocalSteadyOptions options;
        options.cfl = real(args[4]);
        options.maximum_iterations = count(args[5]);
        options.minimum_iterations = std::min<cfd::Index>(20, options.maximum_iterations);
        options.check_interval = count(args[6]);
        options.maximum_wall_seconds = real(args[7]);
        if (options.maximum_iterations > 100000 || options.maximum_wall_seconds > 120)
            throw std::invalid_argument("Verification guard: at most 100000 iterations and 120 s per invocation.");
        const auto problem{make_problem(name, input, initialization)};
        options.reference_density = problem.rho_ref;
        options.reference_velocity = problem.velocity_ref;
        options.reference_length = problem.length_ref;
        return run(problem, mode, output, options, name == "rpa" ? initialization : std::string_view{});
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
