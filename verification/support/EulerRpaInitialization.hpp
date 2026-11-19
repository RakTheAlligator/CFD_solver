#pragma once

#include "EulerNozzle.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/meshing/RpaRegionLayout.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <map>
#include <numbers>
#include <span>
#include <stdexcept>

namespace cfd::verification
{
struct RpaQuasi1dInitialization
{
    double throat_x{}, throat_radius{}, critical_area{}, critical_mass_flow{};
};

/// Replace only ENGINE (region labels from RpaRegion) in an already initialized field.
/// The saved straight-sided inner_wall defines R(x); samples are point values at
/// planar cell centroids, with zero radial velocity, NOT exact 2D cell averages.
/// P0/T0 and gas determine the critical flow, subsonic before the minimum radius,
/// supersonic after it. Exterior conservative components are never written.
/// This is case preparation only: no Euler assembly, BC change or solver advance.
inline RpaQuasi1dInitialization initialize_rpa_quasi1d(const Mesh &mesh, const std::span<const Index> regions,
                                                       const CaloricallyPerfectGas &gas, const double total_pressure,
                                                       const double total_temperature, ConservativeField &field)
{
    if (field.size() != mesh.cell_count() || regions.size() != mesh.cell_count())
        throw std::invalid_argument("RPA quasi1d initialization cardinality mismatch.");
    if (!std::isfinite(total_pressure) || !(total_pressure > 0) || !std::isfinite(total_temperature) ||
        !(total_temperature > 0))
        throw std::invalid_argument("RPA quasi1d requires finite positive total pressure and temperature.");
    const auto groups{mesh.boundary_groups()};
    const auto inner_wall{std::ranges::find_if(groups, [](const auto &group) { return group.name == "inner_wall"; })};
    if (inner_wall == groups.end())
        throw std::invalid_argument("RPA quasi1d requires the inner_wall boundary.");
    std::map<double, double> wall;
    for (Index face = 0; face < mesh.face_count(); ++face)
        if (mesh.face_adjacencies()[face].is_boundary() && mesh.face_boundary_ids()[face] == inner_wall->id)
            for (const auto node : mesh.faces()[face].node_ids)
            {
                const auto point{mesh.nodes()[node]};
                if (!std::isfinite(point.x) || !std::isfinite(point.y) || !(point.y > 0))
                    throw std::invalid_argument("RPA quasi1d requires finite wall coordinates and positive radii.");
                const auto [position, inserted]{wall.emplace(point.x, point.y)};
                if (!inserted && position->second != point.y)
                    throw std::invalid_argument("RPA quasi1d inner_wall must be single-valued in x.");
            }
    if (wall.size() < 3)
        throw std::invalid_argument("RPA quasi1d requires a complete inner_wall profile.");
    const auto throat{
        std::min_element(wall.begin(), wall.end(), [](const auto &a, const auto &b) { return a.second < b.second; })};
    if (throat == wall.begin() || std::next(throat) == wall.end())
        throw std::invalid_argument("RPA quasi1d requires an interior throat.");
    RpaQuasi1dInitialization result;
    result.throat_x = throat->first;
    result.throat_radius = throat->second;
    result.critical_area = std::numbers::pi * result.throat_radius * result.throat_radius;
    const double gamma = gas.gamma(), R = gas.gas_constant(), g = gamma - 1;
    result.critical_mass_flow = result.critical_area * total_pressure / std::sqrt(total_temperature) *
                                std::sqrt(gamma / R) * std::pow(2 / (gamma + 1), (gamma + 1) / (2 * g));
    if (!std::isfinite(result.critical_area) || !(result.critical_area > 0) ||
        !std::isfinite(result.critical_mass_flow) || !(result.critical_mass_flow > 0))
        throw std::invalid_argument("RPA quasi1d critical area/flow is not positive and representable.");
    Index engine_cells{};
    for (Index cell = 0; cell < mesh.cell_count(); ++cell)
    {
        if (regions[cell] != static_cast<Index>(RpaRegion::Engine))
            continue;
        ++engine_cells;
        const double x{mesh.cell_centers()[cell].x};
        const auto right{wall.lower_bound(x)};
        if (right == wall.end() || x < wall.begin()->first)
            throw std::invalid_argument("RPA quasi1d Engine cell lies outside the inner_wall profile.");
        double radius{right->second};
        if (right->first != x)
        {
            const auto left{std::prev(right)};
            const double fraction{(x - left->first) / (right->first - left->first)};
            radius = std::lerp(left->second, right->second, fraction);
        }
        const double area = std::numbers::pi * radius * radius;
        const double mach = nozzle_mach_from_area(area / result.critical_area, x > result.throat_x, gamma);
        const double theta = 1 + .5 * g * mach * mach, temperature = total_temperature / theta,
                     pressure = total_pressure * std::pow(theta, -gamma / g), density = pressure / (R * temperature);
        field.set_state(cell, gas.to_conservative({density, mach * std::sqrt(gamma * R * temperature), 0, pressure}));
    }
    if (engine_cells == 0)
        throw std::invalid_argument("RPA quasi1d requires Engine cells.");
    return result;
}
} // namespace cfd::verification
