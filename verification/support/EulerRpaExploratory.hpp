#pragma once

#include "cfd/compressible/EulerBoundaryConditions.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/meshing/RpaRegionLayout.hpp"

#include <array>
#include <cmath>
#include <filesystem>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace cfd::verification
{
/// Case-specific constants, not a second thermodynamic model. Ambient T/u are
/// explicitly exploratory; only its pressure is nominal mission input.
struct EulerRpaModel
{
    static constexpr double gamma{1.177}, gas_constant{414.2};
    static constexpr double pressure{1.1700e6}, temperature{3343.5658}, mach{0.1508};
    static constexpr double ambient_pressure{1.0e-4};

    [[nodiscard]] static double total_temperature()
    {
        return temperature * (1.0 + 0.5 * (gamma - 1.0) * mach * mach);
    }

    [[nodiscard]] static double total_pressure()
    {
        return pressure * std::pow(total_temperature() / temperature, gamma / (gamma - 1.0));
    }

    [[nodiscard]] static PrimitiveState hot()
    {
        return {pressure / (gas_constant * temperature), mach * std::sqrt(gamma * gas_constant * temperature), 0.0,
                pressure};
    }

    [[nodiscard]] static PrimitiveState ambient(const double t, const double ux, const double ur)
    {
        if (!std::isfinite(t) || !(t > 0.0) || !std::isfinite(ux) || !std::isfinite(ur))
            throw std::invalid_argument("RPA ambient requires finite T>0 and finite velocities.");
        const PrimitiveState state{ambient_pressure / (gas_constant * t), ux, ur, ambient_pressure};
        (void)CaloricallyPerfectGas{gamma, gas_constant}.to_conservative(state);
        return state;
    }
};

inline void require_external_empty_output(const std::filesystem::path &input, const std::filesystem::path &output)
{
    const auto source{std::filesystem::weakly_canonical(input)}, target{std::filesystem::weakly_canonical(output)};
    auto position{target.begin()};
    bool contained{true};
    for (const auto &part : source)
    {
        if (position == target.end() || *position != part)
        {
            contained = false;
            break;
        }
        ++position;
    }
    if (contained)
        throw std::invalid_argument("RPA output must be outside the input case, including through symlinks.");
    if (std::filesystem::exists(target) &&
        (!std::filesystem::is_directory(target) || !std::filesystem::is_empty(target)))
        throw std::invalid_argument("RPA output directory must be empty; existing results are never removed.");
}

inline EulerBoundaryConditions rpa_euler_conditions(const std::span<const BoundaryGroup> groups,
                                                    const PrimitiveState ambient)
{
    constexpr std::array<std::string_view, 8> expected{"inlet",           "axis",
                                                       "inner_wall",      "outer_wall",
                                                       "nozzle_lip",      "farfield_upstream",
                                                       "farfield_radial", "farfield_downstream"};
    if (groups.size() != expected.size())
        throw std::invalid_argument("RPA Euler requires exactly eight boundary groups.");
    std::array<bool, 8> seen{};
    std::vector<EulerBoundaryCondition> conditions;
    conditions.reserve(groups.size());
    for (const auto &group : groups)
    {
        Index match{};
        while (match < expected.size() && group.name != expected.at(match))
            ++match;
        if (match == expected.size() || seen.at(match) || group.id != conditions.size())
            throw std::invalid_argument("Unknown, repeated or misindexed RPA boundary group: " + group.name);
        seen.at(match) = true;
        if (match == 0)
            conditions.push_back(EulerBoundaryCondition::total_conditions_inlet(
                EulerRpaModel::total_pressure(), EulerRpaModel::total_temperature(), {1.0, 0.0}));
        else if (match == 1)
            conditions.emplace_back(EulerBoundaryConditionType::Axis);
        else if (match < 5)
            conditions.emplace_back(EulerBoundaryConditionType::SlipWall);
        else
            conditions.push_back(EulerBoundaryCondition::farfield(ambient));
    }
    return {groups.size(), std::move(conditions)};
}

inline void require_connected_rpa(const Mesh &mesh)
{
    if (mesh.cell_count() == 0)
        throw std::invalid_argument("RPA mesh must contain fluid cells.");
    std::vector<bool> reached(mesh.cell_count());
    std::vector<Index> stack;
    stack.reserve(mesh.cell_count());
    stack.push_back(0);
    reached[0] = true;
    Index count{};
    while (!stack.empty())
    {
        const auto cell{stack.back()};
        stack.pop_back();
        ++count;
        for (Index entry = mesh.cell_node_offsets()[cell]; entry < mesh.cell_node_offsets()[cell + 1]; ++entry)
        {
            const auto face{mesh.face_adjacencies()[mesh.cell_faces()[entry]]};
            if (face.is_boundary())
                continue;
            const auto next{face.owner == cell ? face.neighbor : face.owner};
            if (!reached[next])
            {
                reached[next] = true;
                stack.push_back(next);
            }
        }
    }
    if (count != mesh.cell_count())
        throw std::invalid_argument("RPA Euler requires a connected fluid domain.");
}

inline Index initialize_rpa(const Mesh &mesh, const RpaRegionLayout &layout, const PrimitiveState ambient,
                            ConservativeField &field, const std::span<double> regions)
{
    if (field.size() != mesh.cell_count() || regions.size() != mesh.cell_count())
        throw std::invalid_argument("RPA initialization cardinality mismatch.");
    const CaloricallyPerfectGas gas{EulerRpaModel::gamma, EulerRpaModel::gas_constant};
    const auto hot{gas.to_conservative(EulerRpaModel::hot())}, outside{gas.to_conservative(ambient)};
    Index hot_cells{};
    for (Index cell = 0; cell < field.size(); ++cell)
    {
        const auto region{layout.region(mesh.cell_centers()[cell])};
        regions[cell] = static_cast<double>(region);
        field.set_state(cell, region == RpaRegion::Engine ? hot : outside);
        hot_cells += region == RpaRegion::Engine ? 1 : 0;
    }
    if (hot_cells == 0 || hot_cells == field.size())
        throw std::runtime_error("RPA Engine/exterior partition must contain both states.");
    return hot_cells;
}
} // namespace cfd::verification
