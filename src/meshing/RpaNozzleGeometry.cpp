#include "cfd/meshing/RpaNozzleGeometry.hpp"

#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstring>
#include <fcntl.h>
#include <iomanip>
#include <limits>
#include <spawn.h>
#include <sstream>
#include <stdexcept>
#include <sys/wait.h>
#include <unistd.h>

namespace cfd
{
namespace
{

std::string import_cad(const RpaNozzleGeometry &geometry)
{
    if (!std::filesystem::is_regular_file(geometry.file))
    {
        throw std::runtime_error("RPA DXF file is absent: " + geometry.file.string());
    }
    std::string interpreter{CFD_RPA_PYTHON_EXECUTABLE};
    std::string adapter{CFD_RPA_DXF_ADAPTER};
    std::string file{std::filesystem::absolute(geometry.file).string()};
    std::ostringstream tolerance_stream;
    tolerance_stream << std::setprecision(std::numeric_limits<double>::max_digits10) << geometry.geometric_tolerance;
    std::string tolerance{tolerance_stream.str()};
    std::string units{geometry.dxf_units};
    std::array<char *, 6> arguments{interpreter.data(), adapter.data(), file.data(),
                                    tolerance.data(),   units.data(),   nullptr};
    std::array<int, 2> descriptors{};
    if (::pipe2(descriptors.data(), O_CLOEXEC) != 0)
    {
        throw std::runtime_error("Unable to create RPA import pipe: " + std::string{std::strerror(errno)});
    }
    posix_spawn_file_actions_t actions;
    int error{::posix_spawn_file_actions_init(&actions)};
    const bool initialized{error == 0};
    if (error == 0)
    {
        error = ::posix_spawn_file_actions_addclose(&actions, descriptors[0]);
    }
    if (error == 0)
    {
        error = ::posix_spawn_file_actions_adddup2(&actions, descriptors[1], STDOUT_FILENO);
    }
    if (error == 0 && descriptors[1] != STDOUT_FILENO)
    {
        error = ::posix_spawn_file_actions_addclose(&actions, descriptors[1]);
    }
    pid_t child{};
    if (error == 0)
    {
        error = ::posix_spawnp(&child, interpreter.c_str(), &actions, nullptr, arguments.data(), environ);
    }
    if (initialized)
    {
        static_cast<void>(::posix_spawn_file_actions_destroy(&actions));
    }
    static_cast<void>(::close(descriptors[1]));
    if (error != 0)
    {
        static_cast<void>(::close(descriptors[0]));
        throw std::runtime_error("Unable to launch RPA CAD adapter: " + std::string{std::strerror(error)});
    }
    std::string output;
    std::array<char, 4096> buffer{};
    bool read_failed{};
    try
    {
        while (true)
        {
            const ssize_t count{::read(descriptors[0], buffer.data(), buffer.size())};
            if (count == 0)
            {
                break;
            }
            if (count < 0)
            {
                if (errno == EINTR)
                {
                    continue;
                }
                read_failed = true;
                break;
            }
            output.append(buffer.data(), static_cast<std::size_t>(count));
        }
    }
    catch (...)
    {
        static_cast<void>(::close(descriptors[0]));
        int status{};
        while (::waitpid(child, &status, 0) < 0 && errno == EINTR)
        {
        }
        throw;
    }
    static_cast<void>(::close(descriptors[0]));
    int status{};
    pid_t waited{};
    waited = ::waitpid(child, &status, 0);
    while (waited < 0 && errno == EINTR)
    {
        waited = ::waitpid(child, &status, 0);
    }
    if (waited < 0 || read_failed || !WIFEXITED(status) || WEXITSTATUS(status) != 0)
    {
        throw std::runtime_error("RPA CAD adapter failed for " + file + ": " + output);
    }
    return output;
}

double positive_factor(double value, const char *name)
{
    if (!std::isfinite(value) || !(value > 0.0))
        throw std::invalid_argument(std::string{"RPA "} + name + " must be finite and positive.");
    return value;
}

void trim_upstream(std::vector<Point2> &points, const double x)
{
    const Point2 first{interpolate_rpa_profile(points, x)};
    const auto after{std::upper_bound(points.begin(), points.end(), x,
                                      [](double value, const Point2 &point) { return value < point.x; })};
    points.erase(points.begin(), after);
    points.insert(points.begin(), first);
}

} // namespace

Point2 interpolate_rpa_profile(const std::span<const Point2> points, double x)
{
    if (points.empty())
        throw std::invalid_argument("RPA profile must not be empty.");
    const auto upper{std::upper_bound(points.begin(), points.end(), x,
                                      [](double value, const Point2 &point) { return value < point.x; })};
    if (upper == points.begin() || upper == points.end())
        return upper == points.begin() ? points.front() : points.back();
    const auto &left{*(upper - 1)};
    return {x, std::lerp(left.y, upper->y, (x - left.x) / (upper->x - left.x))};
}

RpaNozzleDomain read_rpa_nozzle_domain(const RpaNozzleGeometry &geometry)
{
    if (!std::isfinite(geometry.geometric_tolerance) || !(geometry.geometric_tolerance > 0.0))
    {
        throw std::invalid_argument("RPA geometricTolerance must be finite and positive.");
    }
    if (geometry.dxf_units != "auto" && geometry.dxf_units != "mm" && geometry.dxf_units != "cm" &&
        geometry.dxf_units != "m" && geometry.dxf_units != "in")
    {
        throw std::invalid_argument("RPA dxfUnits must be auto, mm, cm, m or in.");
    }
    std::istringstream input{import_cad(geometry)};
    std::string magic;
    std::size_t count{};
    if (!(input >> magic >> count) || magic != "CFD_RPA_PROFILE_V1" || count < 3 || count > 1000000)
    {
        throw std::runtime_error("Invalid RPA CAD adapter profile header.");
    }
    RpaNozzleDomain domain;
    domain.geometric_tolerance = geometry.geometric_tolerance;
    domain.inner_wall.resize(count);
    double previous_x{-std::numeric_limits<double>::infinity()};
    for (Point2 &point : domain.inner_wall)
    {
        if (!(input >> point.x >> point.y) || !std::isfinite(point.x) || !std::isfinite(point.y) || !(point.y > 0.0) ||
            !(point.x > previous_x))
        {
            throw std::runtime_error("Invalid RPA profile: finite, positive-radius, continuous r(x) required.");
        }
        previous_x = point.x;
    }
    std::size_t junction_count{};
    if (!(input >> magic >> junction_count) || magic != "CURVE_JUNCTIONS" || junction_count < 2 ||
        junction_count > count)
    {
        throw std::runtime_error("Invalid RPA CAD curve junctions.");
    }
    for (std::size_t junction = 0; junction < junction_count; ++junction)
    {
        std::size_t index{};
        if (!(input >> index) || index >= count ||
            (junction > 0 && domain.inner_wall[index].x <= domain.curve_junctions.back()) ||
            (junction == 0 && index != 0) || (junction == junction_count - 1 && index != count - 1))
        {
            throw std::runtime_error("Invalid RPA CAD curve junction index.");
        }
        domain.curve_junctions.push_back(domain.inner_wall[index].x);
    }
    if (input >> magic)
    {
        throw std::runtime_error("Unexpected data after RPA CAD profile.");
    }
    const auto throat{std::min_element(domain.inner_wall.begin(), domain.inner_wall.end(),
                                       [](const Point2 &a, const Point2 &b) { return a.y < b.y; })};
    if (throat == domain.inner_wall.begin() || throat == domain.inner_wall.end() - 1)
    {
        throw std::runtime_error("RPA profile must contain an interior throat.");
    }
    domain.throat = *throat;
    const Point2 exit{domain.inner_wall.back()};
    const double available{exit.x - domain.inner_wall.front().x};
    const double exterior_length{geometry.upstream_factor
                                     ? positive_factor(*geometry.upstream_factor, "box.upstreamFactor") * exit.y
                                     : available};
    if (!std::isfinite(exterior_length) || exterior_length > available || !(exterior_length > 0.0))
        throw std::invalid_argument("RPA box.upstreamFactor exceeds the available engine profile.");
    const double thickness{positive_factor(geometry.wall_thickness_factor, "engine.wallThicknessFactor") * exit.y};
    domain.outer_wall = domain.inner_wall;
    trim_upstream(domain.outer_wall, exit.x - exterior_length);
    for (auto &point : domain.outer_wall)
    {
        point.y += thickness;
        if (!std::isfinite(point.y) || !(point.y > interpolate_rpa_profile(domain.inner_wall, point.x).y))
            throw std::invalid_argument("RPA fictitious shell thickness is not representable.");
    }
    domain.farfield_radius = positive_factor(geometry.outer_radius_factor, "box.outerRadiusFactor") * exit.y;
    for (const auto &point : domain.outer_wall)
        if (!std::isfinite(domain.farfield_radius) || !(point.y < domain.farfield_radius))
            throw std::invalid_argument("RPA box.outerRadiusFactor must enclose the entire exterior wall.");
    domain.downstream_x = exit.x + positive_factor(geometry.downstream_factor, "box.downstreamFactor") * exit.y;
    if (!std::isfinite(domain.downstream_x) || !(domain.downstream_x > exit.x))
        throw std::invalid_argument("RPA downstream extent is not representable.");
    return domain;
}

} // namespace cfd
