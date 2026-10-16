#pragma once

#include "EulerTube.hpp"
#include "cfd/compressible/AxisymmetricEulerSolver.hpp"
#include "cfd/mesh/MeshBuilder.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <numbers>

namespace cfd::verification
{

/// Independent quasi-1D isentropic reference; no CFD flux/assembler calls.
inline long double nozzle_log_area_ratio(const long double mach, const long double gamma)
{
    return -std::log(mach) + (gamma + 1.0L) / (2.0L * (gamma - 1.0L)) *
                                 (std::log(2.0L / (gamma + 1.0L)) + std::log1p(0.5L * (gamma - 1.0L) * mach * mach));
}

inline double nozzle_mach_from_area(const double ratio, const bool supersonic, const double gamma = 1.4)
{
    if (!std::isfinite(ratio) || ratio < 1.0 || !std::isfinite(gamma) || !(gamma > 1.0))
        throw std::invalid_argument("Nozzle area/Mach requires finite A/Astar >= 1 and gamma > 1.");
    if (ratio == 1.0)
        return 1.0;
    const long double target{std::log(static_cast<long double>(ratio))};
    long double lower{supersonic ? 1.0L : 0.0L}, upper{supersonic ? 2.0L : 1.0L};
    if (supersonic)
        while (nozzle_log_area_ratio(upper, gamma) < target)
        {
            upper *= 2.0L;
            if (!std::isfinite(upper * upper))
                throw std::runtime_error("Nozzle supersonic reference could not bracket Mach.");
        }
    for (int iteration = 0; iteration < 120; ++iteration)
    {
        const long double middle{lower + 0.5L * (upper - lower)};
        if (middle == lower || middle == upper)
            break;
        const bool above{nozzle_log_area_ratio(middle, gamma) > target};
        if (above == supersonic)
            upper = middle;
        else
            lower = middle;
    }
    const double result{static_cast<double>(lower + 0.5L * (upper - lower))};
    if (!std::isfinite(result) || !(result > 0.0))
        throw std::runtime_error("Nozzle reference Mach is not representable.");
    return result;
}

struct EulerNozzleReference
{
    bool choked{};
    static constexpr double length{1.0}, throat_radius{0.1}, gamma{1.4}, gas_constant{1.0};
    static constexpr double total_pressure{1.0}, total_temperature{1.0};

    [[nodiscard]] static double area(const double x)
    {
        if (!std::isfinite(x) || x < 0.0 || x > length)
            throw std::invalid_argument("Nozzle reference coordinate must lie in [0,L].");
        const double cosine{std::cos((2.0 * x / length - 1.0) * std::numbers::pi)};
        const double ratio{x <= length / 2.0 ? 1.75 - 0.75 * cosine : 1.25 - 0.25 * cosine};
        return std::numbers::pi * throat_radius * throat_radius * ratio;
    }

    [[nodiscard]] static double radius(const double x)
    {
        return std::sqrt(area(x) / std::numbers::pi);
    }

    [[nodiscard]] double critical_area() const
    {
        const double throat{area(length / 2.0)};
        return choked ? throat : throat / static_cast<double>(std::exp(nozzle_log_area_ratio(0.5L, gamma)));
    }

    [[nodiscard]] double mach(const double x) const
    {
        return nozzle_mach_from_area(area(x) / critical_area(), choked && x > length / 2.0, gamma);
    }

    [[nodiscard]] PrimitiveState primitive(const double x) const
    {
        const double m{mach(x)}, theta{1.0 + 0.5 * (gamma - 1.0) * m * m};
        const double temperature{total_temperature / theta};
        const double pressure{total_pressure * std::pow(theta, -gamma / (gamma - 1.0))};
        return {pressure / (gas_constant * temperature), m * std::sqrt(gamma * gas_constant * temperature), 0.0,
                pressure};
    }

    [[nodiscard]] double mass_flow() const
    {
        return critical_area() * total_pressure / std::sqrt(total_temperature) * std::sqrt(gamma / gas_constant) *
               std::pow(2.0 / (gamma + 1.0), (gamma + 1.0) / (2.0 * (gamma - 1.0)));
    }

    [[nodiscard]] EulerBoundaryConditions conditions() const
    {
        using Type = EulerBoundaryConditionType;
        const auto inlet{EulerBoundaryCondition::total_conditions_inlet(total_pressure, total_temperature, {1.0, 0.0})};
        const auto outlet{choked ? EulerBoundaryCondition{Type::SupersonicOutlet}
                                 : EulerBoundaryCondition::pressure_outlet(primitive(length).pressure)};
        return {4, {Type::Axis, Type::SlipWall, inlet, outlet}};
    }
};

inline RawMeshData make_euler_nozzle(const Index nx, const Index nr)
{
    if (nx < 4 || nx % 2 != 0 || nr == 0)
        throw std::invalid_argument("Euler nozzle requires even nx >= 4 and nr > 0.");
    auto raw{make_euler_tube(nx, nr)};
    raw.boundary_groups = {{0, "axis"}, {1, "wall"}, {2, "inlet"}, {3, "outlet"}};
    for (auto &node : raw.nodes)
    {
        node.x *= EulerNozzleReference::length;
        node.y *= EulerNozzleReference::radius(node.x);
    }
    return raw;
}

/// Four-point Gauss integration of conservative quasi-1D initial data using
/// the revolution area of the actual straight-sided trapezoid at x. This is
/// a convergent numerical projection, not an exact analytic cell average.
inline ConservativeField nozzle_initial_field(const Index nx, const Index nr, const EulerNozzleReference reference,
                                              const bool uniform)
{
    const CaloricallyPerfectGas gas{reference.gamma, reference.gas_constant};
    if (uniform && reference.choked)
        throw std::invalid_argument(
            "Uniform choked startup cannot satisfy both subsonic inlet and supersonic outlet; use quasi1d.");
    if (uniform)
        return ConservativeField{nx * nr, gas.to_conservative(reference.primitive(0.0))};
    constexpr std::array points{-0.8611363115940525752, -0.3399810435848562648, 0.3399810435848562648,
                                0.8611363115940525752};
    constexpr std::array weights{0.3478548451374538574, 0.6521451548625461426, 0.6521451548625461426,
                                 0.3478548451374538574};
    ConservativeField field{nx * nr};
    for (Index i = 0; i < nx; ++i)
    {
        const double xl{static_cast<double>(i) * reference.length / static_cast<double>(nx)};
        const double xr{static_cast<double>(i + 1) * reference.length / static_cast<double>(nx)};
        const double rl{reference.radius(xl)}, rr{reference.radius(xr)};
        long double volume{}, rho{}, momentum{}, energy{};
        for (Index point = 0; point < points.size(); ++point)
        {
            const double fraction{0.5 * (1.0 + points.at(point))};
            const auto primitive{reference.primitive(xl + fraction * (xr - xl))};
            const long double radius{rl + fraction * (rr - rl)};
            // A radial-fraction ring multiplies every integral by the same
            // constant; it cancels from its conservative cell average.
            const long double weight{weights.at(point) * radius * radius};
            volume += weight;
            rho += weight * primitive.density;
            momentum += weight * primitive.density * primitive.velocity_x;
            energy += weight * (primitive.pressure / (reference.gamma - 1.0) +
                                0.5L * primitive.density * primitive.velocity_x * primitive.velocity_x);
        }
        const ConservativeState value{static_cast<double>(rho / volume), static_cast<double>(momentum / volume), 0.0,
                                      static_cast<double>(energy / volume)};
        (void)gas.evaluate(value);
        for (Index j = 0; j < nr; ++j)
            field.set_state(j * nx + i, value);
    }
    return field;
}

} // namespace cfd::verification
