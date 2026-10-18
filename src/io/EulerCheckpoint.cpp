#include "cfd/io/EulerCheckpoint.hpp"

#include "cfd/mesh/Mesh.hpp"

#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <type_traits>

#ifdef __linux__
#include <cerrno>
#include <fcntl.h>
#endif

namespace cfd
{
namespace
{
template <typename... T> void put(std::ostream &out, const T &...values)
{
    ((out << values << ' '), ...);
    out << '\n';
}

template <typename T> void get(std::istream &in, T &value)
{
    std::string token;
    if (!(in >> token))
        throw std::runtime_error("Euler checkpoint: truncated payload.");
    const std::string_view digits{token};
    const auto result{std::from_chars(digits.begin(), digits.end(), value)};
    if (result.ec != std::errc{} || result.ptr != digits.end())
        throw std::runtime_error("Euler checkpoint: invalid scalar.");
    if constexpr (std::is_floating_point_v<T>)
        if (!std::isfinite(value))
            throw std::runtime_error("Euler checkpoint: nonfinite scalar.");
    if constexpr (std::is_same_v<T, long double>)
    {
        std::ostringstream roundtrip;
        roundtrip.imbue(std::locale::classic());
        roundtrip << std::setprecision(std::numeric_limits<long double>::max_digits10) << value;
        if (roundtrip.str() != token)
            throw std::runtime_error(
                "Euler checkpoint: extended-precision ledger is not exactly round-trippable on this platform.");
    }
}

void put_integral(std::ostream &out, const EulerIntegralQuantities &q)
{
    put(out, q.mass, q.momentum_x, q.momentum_r, q.energy);
}

void get_integral(std::istream &in, EulerIntegralQuantities &q)
{
    get(in, q.mass);
    get(in, q.momentum_x);
    get(in, q.momentum_r);
    get(in, q.energy);
}

void put_result(std::ostream &out, const AxisymmetricEulerResult &d)
{
    put(out, d.time, d.time_step, d.cfl, d.accepted_steps, d.rejected_steps, d.minimum_density, d.maximum_density,
        d.minimum_pressure, d.maximum_pressure);
    put_integral(out, d.integrals);
    put_integral(out, d.boundary_outflow_rate);
    put(out, d.radial_pressure_source_rate);
    put_integral(out, d.balance_error);
}

void get_result(std::istream &in, AxisymmetricEulerResult &d)
{
    get(in, d.time);
    get(in, d.time_step);
    get(in, d.cfl);
    get(in, d.accepted_steps);
    get(in, d.rejected_steps);
    get(in, d.minimum_density);
    get(in, d.maximum_density);
    get(in, d.minimum_pressure);
    get(in, d.maximum_pressure);
    get_integral(in, d.integrals);
    get_integral(in, d.boundary_outflow_rate);
    get(in, d.radial_pressure_source_rate);
    get_integral(in, d.balance_error);
}

void put_field(std::ostream &out, const ConservativeField &field)
{
    put(out, field.size());
    for (Index c = 0; c < field.size(); ++c)
    {
        const auto v{field.state(c)};
        put(out, v.density, v.momentum_x, v.momentum_r, v.total_energy_density);
    }
}

void get_field(std::istream &in, ConservativeField &field, const CaloricallyPerfectGas &gas)
{
    Index count{};
    get(in, count);
    if (count != field.size())
        throw std::runtime_error("Euler checkpoint: field cardinality mismatch.");
    for (Index c = 0; c < count; ++c)
    {
        ConservativeState v;
        get(in, v.density);
        get(in, v.momentum_x);
        get(in, v.momentum_r);
        get(in, v.total_energy_density);
        (void)gas.evaluate(v);
        field.set_state(c, v);
    }
}

std::string configuration_text(const EulerCheckpointConfiguration &c)
{
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
    put(out, c.gas.gamma(), c.gas.gas_constant(), c.solver.cfl, c.solver.maximum_step_attempts,
        c.solver.maximum_rejected_steps);
    const auto &m{c.monitor};
    put(out, m.reference_density, m.reference_velocity, m.reference_length, m.residual_tolerance,
        m.field_change_tolerance, m.mass_flow_change_tolerance, m.mass_flow_imbalance_tolerance, m.check_interval,
        m.minimum_steps, m.consecutive_checks);
    out << std::quoted(std::string{c.physical_context}) << '\n';
    const auto &mesh{c.mesh.get()};
    put(out, mesh.node_count(), mesh.cell_count(), mesh.face_count(), mesh.boundary_groups().size());
    for (const auto &n : mesh.nodes())
        put(out, n.x, n.y);
    for (const auto v : mesh.cell_types())
        put(out, static_cast<int>(v));
    for (const auto v : mesh.cell_node_offsets())
        put(out, v);
    for (const auto v : mesh.cell_nodes())
        put(out, v);
    for (const auto v : mesh.cell_faces())
        put(out, v);
    for (Index cell = 0; cell < mesh.cell_count(); ++cell)
        put(out, mesh.cell_areas()[cell], mesh.cell_centers()[cell].x, mesh.cell_centers()[cell].y);
    for (Index f = 0; f < mesh.face_count(); ++f)
    {
        const auto a{mesh.face_adjacencies()[f]};
        const auto nodes{mesh.faces()[f].node_ids};
        put(out, nodes[0], nodes[1], a.owner, a.neighbor, mesh.face_boundary_ids()[f], mesh.face_lengths()[f],
            mesh.face_centers()[f].x, mesh.face_centers()[f].y, mesh.face_area_vectors()[f].x,
            mesh.face_area_vectors()[f].y);
    }
    if (c.boundaries.get().size() != mesh.boundary_groups().size())
        throw std::runtime_error("Euler checkpoint: BC cardinality.");
    for (const auto &group : mesh.boundary_groups())
    {
        const auto &b{c.boundaries.get().condition(group.id)};
        const auto p{b.prescribed_state()};
        out << group.id << ' ' << std::quoted(group.name) << '\n';
        put(out, static_cast<int>(b.type()), b.total_pressure(), b.total_temperature(), b.inward_mass_flow(),
            b.direction().x, b.direction().y, p.density, p.velocity_x, p.velocity_r, p.pressure);
    }
    return out.str();
}

// FNV-1a-64 detects accidental byte corruption, not malicious modification.
// Compatibility itself is an exact full metadata comparison, not this hash.
std::uint64_t checksum(const std::string_view bytes)
{
    std::uint64_t value{14695981039346656037ULL};
    for (const unsigned char byte : bytes)
    {
        value ^= byte;
        value *= 1099511628211ULL;
    }
    return value;
}

EulerCheckpoint decode(std::istream &in, const EulerCheckpointConfiguration &c, const std::string &metadata)
{
    std::string magic;
    unsigned version{};
    Index length{};
    if (!(in >> magic) || magic != "CFD_EULER_CHECKPOINT")
        throw std::runtime_error("Euler checkpoint: invalid magic.");
    get(in, version);
    if (version != 1)
        throw std::runtime_error("Euler checkpoint: unknown version.");
    get(in, length);
    if (length != metadata.size() || in.get() != '\n')
        throw std::runtime_error("Euler checkpoint: incompatible configuration.");
    std::string saved(length, '\0');
    in.read(saved.data(), static_cast<std::streamsize>(length));
    if (!in || saved != metadata)
        throw std::runtime_error("Euler checkpoint: incompatible mesh/EOS/BC/CFL/model/monitor configuration.");
    EulerCheckpoint checkpoint{{ConservativeField{c.mesh.get().cell_count()}, {}},
                               {ConservativeField{c.mesh.get().cell_count()}, {}, 0, false},
                               {}};
    get_result(in, checkpoint.solver.diagnostics);
    get_field(in, checkpoint.solver.conservative, c.gas);
    auto &m{checkpoint.monitor};
    unsigned observed{}, status{};
    get(in, observed);
    if (observed > 1)
        throw std::runtime_error("Euler checkpoint: invalid observed flag.");
    m.observed = observed != 0;
    get(in, m.stable_checks);
    get_field(in, m.previous, c.gas);
    get_result(in, m.last.temporal);
    for (auto &v : m.last.current.absolute_residual)
        get(in, v);
    get(in, m.last.current.mass_flow_in);
    get(in, m.last.current.mass_flow_out);
    for (auto &v : m.last.normalized_residual)
        get(in, v);
    for (auto &v : m.last.field_change)
        get(in, v);
    get(in, m.last.mass_flow_relative_imbalance);
    get(in, m.last.mass_flow_change);
    get(in, status);
    if (status > static_cast<unsigned>(EulerSteadyStatus::TimeLimit))
        throw std::runtime_error("Euler checkpoint: invalid monitor status.");
    m.last.status = static_cast<EulerSteadyStatus>(status);
    auto &h{checkpoint.history};
    get_integral(in, h.initial);
    for (auto &v : h.integrated_boundary)
        get(in, v);
    for (auto &v : h.balance_scales)
    {
        get(in, v);
        if (v < 0)
            throw std::runtime_error("Euler checkpoint: invalid balance scale.");
    }
    get(in, h.integrated_source);
    get(in, h.maximum_step_error);
    Index count{};
    get(in, count);
    get(in, h.next_output);
    if (count > 256 || h.next_output > count)
        throw std::runtime_error("Euler checkpoint: invalid output schedule.");
    double previous{};
    for (Index k = 0; k < count; ++k)
    {
        double t{};
        get(in, t);
        if (!(t > previous) ||
            (k < h.next_output ? t > checkpoint.solver.diagnostics.time : t <= checkpoint.solver.diagnostics.time))
            throw std::runtime_error("Euler checkpoint: inconsistent physical output schedule.");
        h.output_times.push_back(t);
        previous = t;
    }
    in >> std::ws;
    if (!in.eof())
        throw std::runtime_error("Euler checkpoint: unexpected trailing content.");
    // Restore through the classes that own the invariants; never poke private state.
    AxisymmetricEulerSolver solver{c.mesh.get(), checkpoint.solver, c.gas, c.boundaries.get(), c.solver};
    const AxisymmetricEulerRestartState previous_state{ConservativeField{m.previous}, m.last.temporal};
    const AxisymmetricEulerSolver previous_solver{c.mesh.get(), previous_state, c.gas, c.boundaries.get(), c.solver};
    if (m.last.temporal.rejected_steps > checkpoint.solver.diagnostics.rejected_steps)
        throw std::runtime_error("Euler checkpoint: monitor rejection history exceeds current count.");
    const EulerSteadyMonitor monitor{solver, c.monitor, m};
    if (h.maximum_step_error < 0 || h.initial.mass <= 0 || h.initial.energy <= 0 || h.integrated_source < 0)
        throw std::runtime_error("Euler checkpoint: invalid conservation history.");
    // Preserve the ledger exactly, but do not invent a universal relative
    // conservation gate here: a zero net momentum scale can hide cancellation
    // of finite face forces. Physical verification belongs to the driver/tests.
    return checkpoint;
}

std::string payload(const EulerCheckpointConfiguration &c, const EulerCheckpoint &p)
{
    const auto metadata{configuration_text(c)};
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << std::setprecision(std::numeric_limits<long double>::max_digits10);
    out << "CFD_EULER_CHECKPOINT 1\n" << metadata.size() << '\n' << metadata;
    put_result(out, p.solver.diagnostics);
    put_field(out, p.solver.conservative);
    const auto &m{p.monitor};
    put(out, m.observed ? 1 : 0, m.stable_checks);
    put_field(out, m.previous);
    put_result(out, m.last.temporal);
    for (const auto v : m.last.current.absolute_residual)
        put(out, v);
    put(out, m.last.current.mass_flow_in, m.last.current.mass_flow_out);
    for (const auto v : m.last.normalized_residual)
        put(out, v);
    for (const auto v : m.last.field_change)
        put(out, v);
    put(out, m.last.mass_flow_relative_imbalance, m.last.mass_flow_change, static_cast<unsigned>(m.last.status));
    const auto &h{p.history};
    put_integral(out, h.initial);
    for (const auto v : h.integrated_boundary)
        put(out, v);
    for (const auto v : h.balance_scales)
        put(out, v);
    put(out, h.integrated_source, h.maximum_step_error, h.output_times.size(), h.next_output);
    for (const auto t : h.output_times)
        put(out, t);
    auto bytes{out.str()};
    std::istringstream verify{bytes};
    verify.imbue(std::locale::classic());
    (void)decode(verify, c, metadata); // Invalid snapshots are never published.
    return bytes;
}

struct TemporaryFile
{
    std::filesystem::path directory;

    explicit TemporaryFile(const std::filesystem::path &path) : directory(path.string() + ".pending")
    {
        if (!std::filesystem::create_directory(directory))
            throw std::runtime_error("Euler checkpoint: pending publication already exists.");
    }

    ~TemporaryFile()
    {
        std::error_code error;
        std::filesystem::remove_all(directory, error);
    }

    TemporaryFile(const TemporaryFile &) = delete;
    TemporaryFile &operator=(const TemporaryFile &) = delete;
    TemporaryFile(TemporaryFile &&) = delete;
    TemporaryFile &operator=(TemporaryFile &&) = delete;
};
} // namespace

void write_euler_checkpoint(const std::filesystem::path &path, const EulerCheckpointConfiguration &c,
                            const EulerCheckpoint &p)
{
    if (path.filename().empty() || path.filename() == "." || path.filename() == ".." ||
        std::filesystem::exists(std::filesystem::symlink_status(path)))
        throw std::runtime_error("Euler checkpoint: destination must be a new file.");
    const auto bytes{payload(c, p)};
    const TemporaryFile temporary{path};
    const auto name{temporary.directory / "payload"};
    std::ofstream out{name, std::ios::binary};
    out.imbue(std::locale::classic());
    out << bytes << "FNV1A64 " << checksum(bytes) << '\n';
    out.close();
    if (!out)
        throw std::runtime_error("Euler checkpoint: write/close failed.");
#ifdef __linux__
    if (::renameat2(AT_FDCWD, name.c_str(), AT_FDCWD, path.c_str(), RENAME_NOREPLACE) != 0)
        throw std::filesystem::filesystem_error("Euler checkpoint: atomic no-replace publication failed", name, path,
                                                std::error_code{errno, std::generic_category()});
#else
    // Atomic no-replace publication where the Linux rename primitive is unavailable.
    std::filesystem::create_hard_link(name, path);
#endif
}

EulerCheckpoint read_euler_checkpoint(const std::filesystem::path &path, const EulerCheckpointConfiguration &c)
{
    const auto metadata{configuration_text(c)};
    const auto size{std::filesystem::file_size(path)};
    if (size > static_cast<long double>(metadata.size()) + static_cast<long double>(c.mesh.get().cell_count()) * 4096 +
                   1048576)
        throw std::runtime_error("Euler checkpoint: implausible file size.");
    std::ifstream file{path, std::ios::binary};
    const std::string bytes{std::istreambuf_iterator<char>{file}, std::istreambuf_iterator<char>{}};
    if (file.bad() || bytes.size() != size)
        throw std::runtime_error("Euler checkpoint: read failed.");
    const auto footer{bytes.rfind("FNV1A64 ")};
    if (footer == std::string::npos)
        throw std::runtime_error("Euler checkpoint: missing integrity footer.");
    std::istringstream tail{bytes.substr(footer + 8)};
    std::uint64_t expected{};
    get(tail, expected);
    tail >> std::ws;
    if (!tail.eof() || expected != checksum(std::string_view{bytes}.substr(0, footer)))
        throw std::runtime_error("Euler checkpoint: checksum mismatch or trailing content.");
    std::istringstream in{bytes.substr(0, footer)};
    in.imbue(std::locale::classic());
    return decode(in, c, metadata);
}
} // namespace cfd
