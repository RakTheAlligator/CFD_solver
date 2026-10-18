#include "cfd/io/EulerCheckpoint.hpp"
#include "support/EulerBoundaryWave.hpp"
#include "support/TestUtils.hpp"

#include <bit>
#include <chrono>
#include <fstream>
#include <sstream>

namespace
{
using cfd::test::require;
using cfd::test::require_throws;
using cfd::test::require_throws_with_message;

struct Directory
{
    std::filesystem::path path{
        std::filesystem::temp_directory_path() /
        ("cfd-euler-restart-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()))};

    Directory()
    {
        std::filesystem::create_directory(path);
    }

    ~Directory()
    {
        std::error_code error;
        std::filesystem::remove_all(path, error);
    }

    Directory(const Directory &) = delete;
    Directory &operator=(const Directory &) = delete;
    Directory(Directory &&) = delete;
    Directory &operator=(Directory &&) = delete;
};

cfd::EulerSteadyOptions monitor_options()
{
    return {.residual_tolerance = 1e-14,
            .field_change_tolerance = 1e-14,
            .check_interval = 7,
            .minimum_steps = 100,
            .maximum_steps = 1000,
            .maximum_time = 1};
}

void advance(cfd::AxisymmetricEulerSolver &solver, cfd::EulerSteadyMonitor &monitor,
             cfd::EulerCheckpointHistory &history, cfd::Index count)
{
    for (cfd::Index step = 0; step < count; ++step)
    {
        const double target{
            history.next_output < history.output_times.size() ? history.output_times[history.next_output] : 1.0};
        const auto &d{solver.step(target)};
        const std::array rates{d.boundary_outflow_rate.mass, d.boundary_outflow_rate.momentum_x,
                               d.boundary_outflow_rate.momentum_r, d.boundary_outflow_rate.energy};
        for (cfd::Index k = 0; k < 4; ++k)
        {
            history.integrated_boundary.at(k) += static_cast<long double>(d.time_step) * rates.at(k);
            history.balance_scales.at(k) += static_cast<long double>(d.time_step) * std::abs(rates.at(k));
        }
        history.integrated_source += static_cast<long double>(d.time_step) * d.radial_pressure_source_rate;
        history.balance_scales.at(2) += static_cast<long double>(d.time_step) * std::abs(d.radial_pressure_source_rate);
        if (monitor.needs_observation())
            (void)monitor.observe();
        if (history.next_output < history.output_times.size() && d.time == target)
            ++history.next_output;
    }
}

void equal_fields(const cfd::ConservativeField &a, const cfd::ConservativeField &b)
{
    require(a.size() == b.size(), "Equal field cardinalities.");
    for (cfd::Index c = 0; c < a.size(); ++c)
    {
        const auto x{a.state(c)}, y{b.state(c)};
        const std::array left{x.density, x.momentum_x, x.momentum_r, x.total_energy_density};
        const std::array right{y.density, y.momentum_x, y.momentum_r, y.total_energy_density};
        for (cfd::Index k = 0; k < 4; ++k)
            require(std::bit_cast<std::uint64_t>(left.at(k)) == std::bit_cast<std::uint64_t>(right.at(k)),
                    "Conservative values must match bit for bit, including signed zero.");
    }
}

void write_checked_payload(const std::filesystem::path &path, const std::string &payload)
{
    std::uint64_t hash{14695981039346656037ULL};
    for (const unsigned char byte : payload)
    {
        hash ^= byte;
        hash *= 1099511628211ULL;
    }
    std::ofstream out{path};
    out << payload << "FNV1A64 " << hash << '\n';
}

void test_continuous_and_checkpointed_are_identical()
{
    Directory directory;
    const auto built{cfd::verification::build_boundary_wave_mesh(6, 3, 1, .1)};
    const cfd::CaloricallyPerfectGas gas{1.4, 1};
    const auto initial{cfd::verification::boundary_wave_initial(built.mesh, {})};
    const auto bc{cfd::verification::boundary_wave_conditions({}, false, true)};
    const auto options{monitor_options()};
    const cfd::AxisymmetricEulerOptions solver_options{.cfl = .2, .maximum_steps = 1000};
    const cfd::EulerCheckpointConfiguration config{built.mesh,     gas,     bc,
                                                   solver_options, options, "synthetic normal pulse"};
    cfd::AxisymmetricEulerSolver continuous{built.mesh, initial, gas, bc, solver_options};
    cfd::AxisymmetricEulerSolver interrupted{built.mesh, initial, gas, bc, solver_options};
    cfd::EulerSteadyMonitor a{continuous, options}, b{interrupted, options};
    (void)a.observe();
    (void)b.observe();
    cfd::EulerCheckpointHistory ha{.initial = continuous.diagnostics().integrals, .output_times = {.001, .02, .04}};
    auto hb{ha};
    advance(continuous, a, ha, 100);
    advance(interrupted, b, hb, 50);
    require(b.restart_state().last.temporal.accepted_steps == 49,
            "Checkpoint must preserve a non-current monitor observation.");
    const cfd::EulerCheckpoint checkpoint{interrupted.restart_state(), b.restart_state(), hb};
    const auto path{directory.path / "state with spaces.chk"};
    cfd::write_euler_checkpoint(path, config, checkpoint);
    const auto loaded{cfd::read_euler_checkpoint(path, config)};
    equal_fields(loaded.solver.conservative, interrupted.state());
    require(loaded.history.integrated_boundary == hb.integrated_boundary &&
                loaded.history.integrated_source == hb.integrated_source,
            "Extended-precision cumulative ledger must round-trip exactly.");
    cfd::AxisymmetricEulerSolver resumed{built.mesh, loaded.solver, gas, bc, solver_options};
    cfd::EulerSteadyMonitor c{resumed, options, loaded.monitor};
    auto hc{loaded.history};
    advance(resumed, c, hc, 50);
    equal_fields(continuous.state(), resumed.state());
    require(continuous.diagnostics() == resumed.diagnostics(), "All accepted-step diagnostics must be identical.");
    require(ha.integrated_boundary == hc.integrated_boundary && ha.integrated_source == hc.integrated_source &&
                ha.balance_scales == hc.balance_scales && ha.next_output == hc.next_output,
            "Ledger and schedule must be identical.");
    const auto ma{a.restart_state()}, mc{c.restart_state()};
    equal_fields(ma.previous, mc.previous);
    require(ma.last.temporal == mc.last.temporal && ma.last.normalized_residual == mc.last.normalized_residual &&
                ma.last.field_change == mc.last.field_change && ma.stable_checks == mc.stable_checks &&
                ma.last.status == mc.last.status,
            "Actual monitor history, residuals and status must be identical.");
    require(continuous.evaluate_steady_residual().absolute_residual ==
                resumed.evaluate_steady_residual().absolute_residual,
            "Recomputed current residual must match.");
    for (cfd::Index cell = 0; cell < initial.size(); ++cell)
        require(gas.evaluate(continuous.state().state(cell)).primitive.pressure ==
                    gas.evaluate(resumed.state().state(cell)).primitive.pressure,
                "Primitive pressure must match.");
    require_throws<std::runtime_error>([&] { cfd::write_euler_checkpoint(path, config, checkpoint); },
                                       "Never overwrite a checkpoint.");
    const auto bytes{cfd::test::read_text_file(path)};
    auto unknown_version{bytes.substr(0, bytes.rfind("FNV1A64 "))};
    unknown_version.replace(0, unknown_version.find('\n'), "CFD_EULER_CHECKPOINT 99");
    write_checked_payload(directory.path / "version.chk", unknown_version);
    require_throws_with_message<std::runtime_error>(
        [&] { (void)cfd::read_euler_checkpoint(directory.path / "version.chk", config); }, "unknown version",
        "Version rejection must not just be a checksum failure.");
    auto bad_field{bytes.substr(0, bytes.rfind("FNV1A64 "))};
    std::istringstream positions{bad_field};
    std::string magic;
    unsigned version{};
    std::size_t metadata_size{};
    positions >> magic >> version >> metadata_size;
    positions.get();
    positions.seekg(static_cast<std::streamoff>(metadata_size), std::ios::cur);
    std::string line;
    for (int row = 0; row < 5; ++row)
        std::getline(positions, line);
    const auto count_position{static_cast<std::size_t>(positions.tellg())};
    std::getline(positions, line);
    bad_field.replace(count_position, line.size(), "1");
    write_checked_payload(directory.path / "count.chk", bad_field);
    require_throws_with_message<std::runtime_error>(
        [&] { (void)cfd::read_euler_checkpoint(directory.path / "count.chk", config); }, "cardinality",
        "Field cardinality must be checked before reads.");
    for (const auto &bad :
         {bytes.substr(0, bytes.size() / 2), bytes + "unexpected", std::string{"CFD_EULER_CHECKPOINT 99\n"}})
    {
        const auto corrupt{directory.path / "corrupt.chk"};
        {
            std::ofstream stream{corrupt};
            stream << bad;
        }
        require_throws<std::runtime_error>([&] { (void)cfd::read_euler_checkpoint(corrupt, config); },
                                           "Truncation/version/trailing content.");
    }
    auto corrupt{bytes};
    corrupt.at(corrupt.size() / 2) ^= 1;
    const auto corrupt_path{directory.path / "bitflip.chk"};
    {
        std::ofstream stream{corrupt_path};
        stream << corrupt;
    }
    require_throws_with_message<std::runtime_error>([&] { (void)cfd::read_euler_checkpoint(corrupt_path, config); },
                                                    "checksum", "Bit corruption.");
    const auto other_mesh{cfd::verification::build_boundary_wave_mesh(6, 3, 1.01, .1)};
    auto incompatible{config};
    incompatible.solver.cfl = .1;
    require_throws<std::runtime_error>([&] { (void)cfd::read_euler_checkpoint(path, incompatible); }, "Wrong CFL.");
    incompatible = config; // Reference wrappers retain borrowed ownership without rebinding the underlying objects.
    incompatible.monitor.residual_tolerance *= 2;
    require_throws<std::runtime_error>([&] { (void)cfd::read_euler_checkpoint(path, incompatible); },
                                       "Wrong monitor tolerance.");
    require_throws<std::runtime_error>(
        [&] {
            (void)cfd::read_euler_checkpoint(
                path, {other_mesh.mesh, gas, bc, solver_options, options, "synthetic normal pulse"});
        },
        "Wrong mesh with same cell count.");
    auto reordered_raw{cfd::verification::make_euler_tube(6, 3)};
    for (auto &node : reordered_raw.nodes)
        node.y *= .1;
    const std::span<cfd::Index> nodes{reordered_raw.cell_nodes};
    std::ranges::swap_ranges(nodes.first(4), nodes.subspan(4, 4));
    const auto reordered{cfd::build_mesh(std::move(reordered_raw))};
    require_throws<std::runtime_error>(
        [&] {
            (void)cfd::read_euler_checkpoint(
                path, {reordered.mesh, gas, bc, solver_options, options, "synthetic normal pulse"});
        },
        "Changed cell ordering.");
    require_throws<std::runtime_error>(
        [&] {
            (void)cfd::read_euler_checkpoint(
                path, {built.mesh, gas, bc, solver_options, options, "different physical assumptions"});
        },
        "Physical assumptions mismatch.");
    require_throws<std::runtime_error>(
        [&] {
            (void)cfd::read_euler_checkpoint(
                path, {built.mesh, {1.41, 1}, bc, solver_options, options, "synthetic normal pulse"});
        },
        "Wrong EOS.");
    const auto other_bc{cfd::verification::euler_tube_conditions()};
    require_throws<std::runtime_error>(
        [&] {
            (void)cfd::read_euler_checkpoint(
                path, {built.mesh, gas, other_bc, solver_options, options, "synthetic normal pulse"});
        },
        "Wrong boundary conditions.");
    require_throws<std::filesystem::filesystem_error>(
        [&] { (void)cfd::read_euler_checkpoint(directory.path / "missing", config); }, "Missing path.");
    auto invalid{loaded};
    invalid.solver.conservative.density()[0] = -1;
    require_throws<std::runtime_error>(
        [&] { cfd::write_euler_checkpoint(directory.path / "inadmissible.chk", config, invalid); },
        "Invalid state cannot be published.");
    invalid.solver.conservative.density()[0] = std::numeric_limits<double>::quiet_NaN();
    require_throws<std::runtime_error>(
        [&] { cfd::write_euler_checkpoint(directory.path / "nonfinite.chk", config, invalid); },
        "Nonfinite checkpoint cannot be published.");
    auto false_convergence{loaded};
    false_convergence.monitor.stable_checks = options.consecutive_checks;
    false_convergence.monitor.last.status = cfd::EulerSteadyStatus::Converged;
    require_throws<std::invalid_argument>(
        [&] { cfd::write_euler_checkpoint(directory.path / "false-converged.chk", config, false_convergence); },
        "Convergence history must agree with the retained residuals.");
    require_throws<std::filesystem::filesystem_error>(
        [&] { cfd::write_euler_checkpoint(directory.path / "missing-parent" / "state.chk", config, loaded); },
        "Missing output parent.");
}

void test_terminal_and_initial_monitor_states()
{
    const auto built{cfd::verification::build_boundary_wave_mesh(2, 1, 1, .1)};
    const cfd::CaloricallyPerfectGas gas{1.4, 1};
    const auto bc{cfd::verification::euler_tube_conditions()};
    const cfd::ConservativeField field{2, gas.to_conservative({1, 0, 0, 1})};
    const cfd::AxisymmetricEulerOptions opt{.maximum_steps = 1};
    cfd::AxisymmetricEulerSolver solver{built.mesh, field, gas, bc, opt};
    cfd::EulerSteadyOptions mon{.minimum_steps = 1, .maximum_steps = 1, .maximum_time = 1};
    cfd::EulerSteadyMonitor monitor{solver, mon};
    const auto virgin{monitor.restart_state()};
    cfd::EulerSteadyMonitor restored_initial{solver, mon, virgin};
    require(restored_initial.observe().status == cfd::EulerSteadyStatus::Running,
            "No invented initial convergence checks.");
    (void)monitor.observe();
    solver.step(1);
    (void)monitor.observe();
    cfd::AxisymmetricEulerSolver terminal{built.mesh, solver.restart_state(), gas, bc, opt};
    cfd::EulerSteadyMonitor restored{terminal, mon, monitor.restart_state()};
    require(restored.restart_state().last.status == cfd::EulerSteadyStatus::StepLimit, "Terminal status is preserved.");
    require_throws<std::runtime_error>([&] { terminal.step(1); },
                                       "A terminal step budget cannot be bypassed by restart.");
}

void test_rejections_tiny_states_and_converged_history()
{
    Directory directory;
    const auto built{cfd::verification::build_boundary_wave_mesh(1, 1, 1, 1)};
    const cfd::CaloricallyPerfectGas gas{1.4, 1};
    using Type = cfd::EulerBoundaryConditionType;
    const cfd::EulerBoundaryConditions bc{4, {Type::Axis, Type::Transmissive, Type::Transmissive, Type::Transmissive}};
    const cfd::AxisymmetricEulerOptions opt{.cfl = 1, .maximum_steps = 1000};
    const auto mon{monitor_options()};
    cfd::AxisymmetricEulerSolver solver{built.mesh, cfd::ConservativeField{1, gas.to_conservative({1, 0, 10, .001})},
                                        gas, bc, opt};
    cfd::EulerSteadyMonitor monitor{solver, mon};
    (void)monitor.observe();
    cfd::EulerCheckpointHistory h{.initial = solver.diagnostics().integrals};
    advance(solver, monitor, h, 1);
    require(solver.diagnostics().rejected_steps == 1, "Exercise an actual rejected candidate before checkpoint.");
    const cfd::EulerCheckpointConfiguration config{built.mesh, gas, bc, opt, mon, "radial retry"};
    cfd::write_euler_checkpoint(directory.path / "retry.chk", config,
                                {solver.restart_state(), monitor.restart_state(), h});
    const auto checkpoint{cfd::read_euler_checkpoint(directory.path / "retry.chk", config)};
    cfd::AxisymmetricEulerSolver resumed{built.mesh, checkpoint.solver, gas, bc, opt};
    require(resumed.diagnostics().rejected_steps == 1, "Rejection counts must not reset on restart.");
    solver.step(1);
    resumed.step(1);
    equal_fields(solver.state(), resumed.state());
    require(solver.diagnostics() == resumed.diagnostics(), "Retry continuation must be exact.");
    const cfd::AxisymmetricEulerOptions quiet_opt{};
    const cfd::ConservativeField tiny{1, {1e-200, -0.0, 0.0, 1e-200}};
    cfd::AxisymmetricEulerSolver small{built.mesh, tiny, gas, bc, quiet_opt};
    cfd::EulerSteadyMonitor small_monitor{small, mon};
    const cfd::EulerCheckpointConfiguration small_config{built.mesh, gas, bc, quiet_opt, mon, "tiny exact values"};
    cfd::write_euler_checkpoint(
        directory.path / "tiny.chk", small_config,
        {small.restart_state(), small_monitor.restart_state(), {.initial = small.diagnostics().integrals}});
    equal_fields(cfd::read_euler_checkpoint(directory.path / "tiny.chk", small_config).solver.conservative, tiny);
    cfd::EulerSteadyOptions steady{
        .check_interval = 1, .minimum_steps = 1, .maximum_steps = 10, .consecutive_checks = 1, .maximum_time = 1};
    cfd::AxisymmetricEulerSolver rest{built.mesh, cfd::ConservativeField{1, gas.to_conservative({1, 0, 0, 1})}, gas, bc,
                                      quiet_opt};
    cfd::EulerSteadyMonitor observer{rest, steady};
    (void)observer.observe();
    cfd::EulerCheckpointHistory hs{.initial = rest.diagnostics().integrals};
    advance(rest, observer, hs, 1);
    require(observer.restart_state().last.status == cfd::EulerSteadyStatus::Converged,
            "Synthetic resting state must converge.");
    const cfd::EulerCheckpointConfiguration steady_config{built.mesh, gas, bc, quiet_opt, steady, "steady terminal"};
    cfd::write_euler_checkpoint(directory.path / "steady.chk", steady_config,
                                {rest.restart_state(), observer.restart_state(), hs});
    const auto loaded{cfd::read_euler_checkpoint(directory.path / "steady.chk", steady_config)};
    cfd::AxisymmetricEulerSolver terminal{built.mesh, loaded.solver, gas, bc, quiet_opt};
    cfd::EulerSteadyMonitor terminal_monitor{terminal, steady, loaded.monitor};
    require(terminal_monitor.restart_state().last.status == cfd::EulerSteadyStatus::Converged &&
                terminal_monitor.restart_state().stable_checks == 1,
            "Restart preserves a genuinely converged terminal history.");
}
} // namespace

int main()
{
    int failures{};
    failures += cfd::test::run_test("exact 100 vs 50+checkpoint+50, integrity and compatibility",
                                    test_continuous_and_checkpointed_are_identical);
    failures += cfd::test::run_test("initial and terminal restart history", test_terminal_and_initial_monitor_states);
    failures += cfd::test::run_test("rejections, tiny values and converged history",
                                    test_rejections_tiny_states_and_converged_history);
    return cfd::test::finish_tests(failures, "Euler checkpoint");
}
