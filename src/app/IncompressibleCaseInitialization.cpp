#include "app/IncompressibleCaseInitialization.hpp"
#include "app/BoundaryBasedInitialization.hpp"

#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVectorField.hpp"
#include "cfd/field/CellVelocityField.hpp"
#include "cfd/field/FaceFluxField.hpp"
#include "cfd/field/PressureCorrectionBoundaryConditions.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/input/OpenFOAMCaseReader.hpp"
#include "cfd/mesh/Boundary.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/numerics/CellFieldTransfer.hpp"
#include "cfd/numerics/IncompressibleSimpleSolver.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cfd::app
{

[[nodiscard]]
cfd::PressureCorrectionBoundaryConditions make_pressure_correction_boundary_conditions(
    const cfd::ScalarBoundaryConditions &pressure_boundary_conditions)
{
    std::vector<cfd::PressureCorrectionBoundaryConditionType> conditions;
    conditions.reserve(pressure_boundary_conditions.size());

    for (cfd::BoundaryId boundary_id = 0; boundary_id < pressure_boundary_conditions.size(); ++boundary_id)
    {
        switch (pressure_boundary_conditions[boundary_id].type)
        {
        case cfd::ScalarBoundaryConditionType::Dirichlet:
            conditions.push_back(cfd::PressureCorrectionBoundaryConditionType::FixedPressure);
            break;

        case cfd::ScalarBoundaryConditionType::Neumann:
            conditions.push_back(cfd::PressureCorrectionBoundaryConditionType::FixedMassFlux);
            break;
        }
    }

    return {pressure_boundary_conditions.size(), std::move(conditions)};
}

void initialize_fixed_mass_flux_boundaries(
    const cfd::Mesh &mesh, const double density, const cfd::ScalarBoundaryConditions &u_boundary_conditions,
    const cfd::ScalarBoundaryConditions &v_boundary_conditions,
    const cfd::PressureCorrectionBoundaryConditions &pressure_correction_boundary_conditions,
    cfd::FaceFluxField &mass_flux)
{
    const auto face_adjacencies{mesh.face_adjacencies()};
    const auto face_boundary_ids{mesh.face_boundary_ids()};
    const auto face_area_vectors{mesh.face_area_vectors()};
    const auto boundary_groups{mesh.boundary_groups()};

    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (!face_adjacencies[face_id].is_boundary())
        {
            continue;
        }

        const cfd::BoundaryId boundary_id{face_boundary_ids[face_id]};
        switch (pressure_correction_boundary_conditions[boundary_id])
        {
        case cfd::PressureCorrectionBoundaryConditionType::FixedMassFlux: {
            const cfd::ScalarBoundaryCondition &u_condition{u_boundary_conditions[boundary_id]};
            const cfd::ScalarBoundaryCondition &v_condition{v_boundary_conditions[boundary_id]};
            if (u_condition.type != cfd::ScalarBoundaryConditionType::Dirichlet ||
                v_condition.type != cfd::ScalarBoundaryConditionType::Dirichlet)
            {
                throw std::invalid_argument("FixedMassFlux boundary '" + boundary_groups[boundary_id].name +
                                            "' requires fixedValue conditions for both u and v.");
            }

            const cfd::Vector2 &area_vector{face_area_vectors[face_id]};
            const double boundary_mass_flux{density *
                                            (u_condition.value * area_vector.x + v_condition.value * area_vector.y)};
            if (!std::isfinite(boundary_mass_flux))
            {
                throw std::runtime_error("Computed mass flux is non-finite on boundary '" +
                                         boundary_groups[boundary_id].name + "'.");
            }
            mass_flux[face_id] = boundary_mass_flux;
            break;
        }

        case cfd::PressureCorrectionBoundaryConditionType::FixedPressure:
            break;
        }
    }
}

std::vector<GridLevelPlan> make_grid_sequencing_plan(Index final_cell_count, double final_mesh_size,
                                                     const input::InitializationInput &initialization)
{
    switch (initialization.type)
    {
    case input::InitializationType::Zero:
        if (initialization.target_coarse_cell_count.has_value())
        {
            throw std::invalid_argument("targetCoarseCellCount is incompatible with initialization type zero.");
        }
        return {};
    case input::InitializationType::CoarseMesh:
        break;
    default:
        throw std::invalid_argument("Unsupported case initialization type.");
    }
    if (final_cell_count == 0 || !std::isfinite(final_mesh_size) || !(final_mesh_size > 0.0))
    {
        throw std::invalid_argument("Grid sequencing requires a nonempty final mesh and positive finite mesh size.");
    }
    constexpr Index level_count{3};
    const auto target{initialization.target_coarse_cell_count};
    if (target.has_value() && (*target == 0 || *target >= final_cell_count))
    {
        throw std::invalid_argument("targetCoarseCellCount must be positive and smaller than the final cell count.");
    }
    const double first_size{target.has_value() ? final_mesh_size * std::sqrt(static_cast<double>(final_cell_count) /
                                                                             static_cast<double>(*target))
                                               : final_mesh_size * std::pow(2.0, static_cast<double>(level_count - 1))};
    if (!std::isfinite(first_size) || !(first_size > final_mesh_size))
    {
        throw std::invalid_argument("Grid-sequencing mesh sizes must be finite and strictly decreasing.");
    }
    const double ratio{std::pow(first_size / final_mesh_size, 1.0 / static_cast<double>(level_count - 1))};
    std::vector<GridLevelPlan> plan;
    plan.reserve(level_count);
    for (Index level = 0; level < level_count; ++level)
    {
        const bool final_level{level + 1 == level_count};
        const double size{final_level ? final_mesh_size : first_size / std::pow(ratio, static_cast<double>(level))};
        if (!std::isfinite(size) || !(size > 0.0) || (!plan.empty() && !(plan.back().mesh_size > size)))
        {
            throw std::invalid_argument("Grid-sequencing mesh sizes must be finite and strictly decreasing.");
        }
        const long double scale{static_cast<long double>(size) / final_mesh_size};
        const long double estimate{static_cast<long double>(final_cell_count) / (scale * scale)};
        // Clamp before conversion so the final Index range is never narrowed through floating point.
        Index count{estimate >= static_cast<long double>(final_cell_count)
                        ? final_cell_count
                        : static_cast<Index>(std::max(1.0L, std::floor(estimate)))};
        if (level == 0 && target.has_value())
        {
            count = *target;
        }
        if (final_level)
        {
            count = final_cell_count;
        }
        plan.push_back({count, size});
    }
    return plan;
}

namespace
{

struct GridLevelState
{
    MeshBuildResult built;
    ScalarBoundaryConditions u_boundary;
    ScalarBoundaryConditions v_boundary;
    ScalarBoundaryConditions pressure_boundary;
    PressureCorrectionBoundaryConditions correction_boundary;
    CellVelocityField velocity;
    CellScalarField pressure;

    GridLevelState(MeshBuildResult build, const input::ScalarFieldInput &u_input,
                   const input::ScalarFieldInput &v_input, const input::ScalarFieldInput &pressure_input)
        : built(std::move(build)), u_boundary(input::resolve_boundary_conditions(built.mesh, u_input)),
          v_boundary(input::resolve_boundary_conditions(built.mesh, v_input)),
          pressure_boundary(input::resolve_boundary_conditions(built.mesh, pressure_input)),
          correction_boundary(make_pressure_correction_boundary_conditions(pressure_boundary)),
          velocity(built.mesh.cell_count(), Vector2{u_input.internal_value, v_input.internal_value}),
          pressure(built.mesh.cell_count(), pressure_input.internal_value)
    {
    }
};

void hierarchy_failure(GridSequencingResult &result, const input::InitializationInput &initialization,
                       const std::string &reason)
{
    if (initialization.target_coarse_cell_count.has_value())
    {
        throw std::invalid_argument("Explicit targetCoarseCellCount cannot form a usable increasing grid hierarchy: " +
                                    reason);
    }
    result.automatic_fallback = true;
    result.fallback_reason = reason;
}

void solve_auxiliary_level(GridLevelState &state, GridLevelReport &report, Index level, double density,
                           double dynamic_viscosity, ScalarConvectionScheme scheme,
                           const IncompressibleSimpleOptions &options)
{
    FaceFluxField flux{state.built.mesh.face_count()};
    initialize_fixed_mass_flux_boundaries(state.built.mesh, density, state.u_boundary, state.v_boundary,
                                          state.correction_boundary, flux);
    IncompressibleSimpleSolver solver{state.built.mesh, density, dynamic_viscosity, scheme, options};
    const auto solved{solver.solve(state.u_boundary, state.v_boundary, state.pressure_boundary,
                                   state.correction_boundary, state.velocity, state.pressure, flux)};
    if (!solved.converged)
    {
        throw std::runtime_error("Grid-sequencing SIMPLE level " + std::to_string(level) + " did not converge within " +
                                 std::to_string(solved.iteration_count) + " iterations.");
    }
    report.iteration_count.emplace(solved.iteration_count);
    report.solve_seconds = solved.timings.total_seconds;
}

} // namespace

GridSequencingResult initialize_with_grid_sequencing(const Mesh &final_mesh, const input::MeshInput &mesh_input,
                                                     const input::InitializationInput &initialization,
                                                     const input::ScalarFieldInput &u_input,
                                                     const input::ScalarFieldInput &v_input,
                                                     const input::ScalarFieldInput &pressure_input, double density,
                                                     double dynamic_viscosity, ScalarConvectionScheme scheme,
                                                     const IncompressibleSimpleOptions &simple_options,
                                                     CellVelocityField &velocity, CellScalarField &pressure)
{
    GridSequencingResult result;
    const auto plan{
        make_grid_sequencing_plan(final_mesh.cell_count(), mesh_input.generation_options.mesh_size, initialization)};
    if (plan.empty())
    {
        return result;
    }
    result.levels.reserve(plan.size());
    for (const auto &level : plan)
    {
        result.levels.push_back({.plan = level});
    }
    result.levels.back().actual_cell_count = final_mesh.cell_count();
    if (final_mesh.cell_count() < plan.size())
    {
        hierarchy_failure(result, initialization, "Too few final cells for a strictly increasing hierarchy.");
        return result;
    }
    if (velocity.size() != final_mesh.cell_count() || pressure.size() != final_mesh.cell_count())
    {
        throw std::invalid_argument("Grid-sequencing output cardinalities must match the final Mesh.");
    }

    // Moving only owning pointers keeps Mesh addresses stable while transfer borrows them.
    std::unique_ptr<GridLevelState> source;
    for (Index level = 0; level < plan.size(); ++level)
    {
        const bool final_level{level + 1 == plan.size()};
        std::unique_ptr<GridLevelState> target;
        if (!final_level)
        {
            std::optional<MeshBuildResult> built;
            try
            {
                MeshGenerationOptions options{mesh_input.generation_options};
                options.mesh_size = plan.at(level).mesh_size;
                built.emplace(build_mesh(generate_mesh(mesh_input.geometry, options, mesh_input.automatic_meshing,
                                                       mesh_input.backward_facing_step_meshing)));
            }
            catch (const std::runtime_error &error)
            {
                hierarchy_failure(result, initialization,
                                  "Meshing level " + std::to_string(level) + ": " + error.what());
                return result;
            }
            const Index actual{built->mesh.cell_count()};
            result.levels.at(level).actual_cell_count = actual;
            if (actual >= final_mesh.cell_count() || (source && actual <= source->built.mesh.cell_count()))
            {
                hierarchy_failure(result, initialization,
                                  "Actual cell counts are not strictly increasing below the final Mesh.");
                return result;
            }
            target = std::make_unique<GridLevelState>(std::move(*built), u_input, v_input, pressure_input);
            if (level == 0)
            {
                initialize_from_boundary_conditions(target->built.mesh, target->u_boundary, target->v_boundary,
                                                    target->pressure_boundary, target->correction_boundary,
                                                    {u_input.internal_value, v_input.internal_value},
                                                    pressure_input.internal_value, target->velocity, target->pressure);
            }
        }
        if (source)
        {
            const Mesh &target_mesh{final_level ? final_mesh : target->built.mesh};
            CellVelocityField &target_velocity{final_level ? velocity : target->velocity};
            CellScalarField &target_pressure{final_level ? pressure : target->pressure};
            const auto transfer_start{std::chrono::steady_clock::now()};
            {
                std::optional<CellFieldTransfer> mapping;
                try
                {
                    // Coverage is validated before touching target fields, including the final ones.
                    mapping.emplace(source->built.mesh, target_mesh);
                }
                catch (const std::runtime_error &error)
                {
                    hierarchy_failure(result, initialization,
                                      "Transfer to level " + std::to_string(level) + ": " + error.what());
                    return result;
                }
                CellVectorField workspace{source->built.mesh.cell_count()};
                mapping->apply_linear_reconstruction(source->velocity.u(), source->u_boundary, workspace,
                                                     target_velocity.u());
                mapping->apply_linear_reconstruction(source->velocity.v(), source->v_boundary, workspace,
                                                     target_velocity.v());
                mapping->apply_linear_reconstruction(source->pressure, source->pressure_boundary, workspace,
                                                     target_pressure);
            }
            result.levels.at(level - 1).transfer_seconds =
                std::chrono::duration<double>(std::chrono::steady_clock::now() - transfer_start).count();
        }
        // The old source, its fields and BCs are freed before solving the new level.
        source = std::move(target);
        if (!final_level)
        {
            solve_auxiliary_level(*source, result.levels.at(level), level, density, dynamic_viscosity, scheme,
                                  simple_options);
        }
    }
    result.used_grid_sequencing = true;
    return result;
}

} // namespace cfd::app
