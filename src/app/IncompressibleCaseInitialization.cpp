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

CoarseMeshPlan make_coarse_mesh_plan(Index final_cell_count, double final_mesh_size,
                                     const input::InitializationInput &initialization)
{
    if (final_cell_count == 0 || !std::isfinite(final_mesh_size) || !(final_mesh_size > 0.0))
    {
        throw std::invalid_argument(
            "Coarse initialization requires a nonempty final mesh and positive finite mesh size.");
    }
    const Index target{initialization.target_coarse_cell_count.value_or(std::max(Index{1}, final_cell_count / 4))};
    if (target == 0 || (initialization.target_coarse_cell_count.has_value() && target >= final_cell_count))
    {
        throw std::invalid_argument("targetCoarseCellCount must be positive and smaller than the final cell count.");
    }
    const double mesh_size{final_mesh_size *
                           std::sqrt(static_cast<double>(final_cell_count) / static_cast<double>(target))};
    if (!std::isfinite(mesh_size) || !(mesh_size > 0.0))
    {
        throw std::invalid_argument("Computed coarse mesh size must be finite and positive.");
    }
    return {.target_cell_count = target, .mesh_size = mesh_size};
}

CoarseInitializationResult initialize_from_coarse_mesh(const Mesh &final_mesh, const input::MeshInput &mesh_input,
                                                       const input::InitializationInput &initialization,
                                                       const input::ScalarFieldInput &u_input,
                                                       const input::ScalarFieldInput &v_input,
                                                       const input::ScalarFieldInput &pressure_input, double density,
                                                       double dynamic_viscosity, ScalarConvectionScheme scheme,
                                                       const IncompressibleSimpleOptions &simple_options,
                                                       CellVelocityField &velocity, CellScalarField &pressure)
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
    CoarseInitializationResult result;
    result.plan =
        make_coarse_mesh_plan(final_mesh.cell_count(), mesh_input.generation_options.mesh_size, initialization);
    if (result.plan.target_cell_count >= final_mesh.cell_count())
    {
        result.automatic_fallback = true;
        result.fallback_reason = "The final mesh has only one cell.";
        return result;
    }
    MeshGenerationOptions coarse_options{mesh_input.generation_options};
    coarse_options.mesh_size = result.plan.mesh_size;
    std::optional<MeshBuildResult> built;
    try
    {
        built.emplace(build_mesh(generate_mesh(mesh_input.geometry, coarse_options, mesh_input.automatic_meshing,
                                               mesh_input.backward_facing_step_meshing)));
    }
    catch (const std::runtime_error &error)
    {
        if (initialization.target_coarse_cell_count.has_value())
        {
            throw std::invalid_argument("Explicit targetCoarseCellCount could not produce a usable coarse mesh: " +
                                        std::string{error.what()});
        }
        result.automatic_fallback = true;
        result.fallback_reason = error.what();
        return result;
    }
    const Mesh &coarse_mesh{built->mesh};
    result.actual_cell_count = coarse_mesh.cell_count();
    if (coarse_mesh.cell_count() >= final_mesh.cell_count())
    {
        if (initialization.target_coarse_cell_count.has_value())
        {
            throw std::invalid_argument("Explicit targetCoarseCellCount did not produce a smaller coarse mesh.");
        }
        result.automatic_fallback = true;
        result.fallback_reason = "The generated coarse mesh is not smaller than the final mesh.";
        return result;
    }
    const ScalarBoundaryConditions u_boundary{input::resolve_boundary_conditions(coarse_mesh, u_input)};
    const ScalarBoundaryConditions v_boundary{input::resolve_boundary_conditions(coarse_mesh, v_input)};
    const ScalarBoundaryConditions pressure_boundary{input::resolve_boundary_conditions(coarse_mesh, pressure_input)};
    const PressureCorrectionBoundaryConditions correction_boundary{
        make_pressure_correction_boundary_conditions(pressure_boundary)};
    CellVelocityField coarse_velocity{coarse_mesh.cell_count(),
                                      Vector2{u_input.internal_value, v_input.internal_value}};
    CellScalarField coarse_pressure{coarse_mesh.cell_count(), pressure_input.internal_value};
    initialize_from_boundary_conditions(coarse_mesh, u_boundary, v_boundary, pressure_boundary, correction_boundary,
                                        {u_input.internal_value, v_input.internal_value}, pressure_input.internal_value,
                                        coarse_velocity, coarse_pressure);
    {
        FaceFluxField coarse_flux{coarse_mesh.face_count()};
        initialize_fixed_mass_flux_boundaries(coarse_mesh, density, u_boundary, v_boundary, correction_boundary,
                                              coarse_flux);
        IncompressibleSimpleSolver coarse_solver{coarse_mesh, density, dynamic_viscosity, scheme, simple_options};
        const IncompressibleSimpleResult coarse_result{coarse_solver.solve(u_boundary, v_boundary, pressure_boundary,
                                                                           correction_boundary, coarse_velocity,
                                                                           coarse_pressure, coarse_flux)};
        if (!coarse_result.converged)
        {
            throw std::runtime_error("Coarse SIMPLE did not converge within " +
                                     std::to_string(coarse_result.iteration_count) + " iterations.");
        }
        result.iteration_count = coarse_result.iteration_count;
        result.solve_seconds = coarse_result.timings.total_seconds;
    }
    const auto transfer_start{std::chrono::steady_clock::now()};
    const CellFieldTransfer mapping{coarse_mesh, final_mesh};
    CellVectorField gradient_workspace{coarse_mesh.cell_count()};
    mapping.apply_linear_reconstruction(coarse_velocity.u(), u_boundary, gradient_workspace, velocity.u());
    mapping.apply_linear_reconstruction(coarse_velocity.v(), v_boundary, gradient_workspace, velocity.v());
    mapping.apply_linear_reconstruction(coarse_pressure, pressure_boundary, gradient_workspace, pressure);
    result.transfer_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - transfer_start).count();
    result.used_coarse_mesh = true;
    return result;
}

} // namespace cfd::app
