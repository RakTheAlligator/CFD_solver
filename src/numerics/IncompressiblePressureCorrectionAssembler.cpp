#include "cfd/numerics/IncompressiblePressureCorrectionAssembler.hpp"

#include "cfd/field/FaceFluxField.hpp"
#include "cfd/field/FacePressureResponseField.hpp"
#include "cfd/field/PressureCorrectionBoundaryConditions.hpp"
#include "cfd/linear_algebra/ScalarLinearSystem.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/Types.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace cfd
{
namespace
{

void validate_system(const Mesh &mesh, const ScalarLinearSystem &system)
{
    if (&system.mesh() != &mesh || system.cell_count() != mesh.cell_count() ||
        system.diagonal().size() != mesh.cell_count() || system.rhs().size() != mesh.cell_count() ||
        system.face_count() != mesh.face_count() || system.owner_neighbor_coefficients().size() != mesh.face_count() ||
        system.neighbor_owner_coefficients().size() != mesh.face_count())
    {
        throw std::invalid_argument("Pressure-correction system must reference the assembler Mesh with matching "
                                    "cardinalities.");
    }
}

void validate_flux_cardinality(const Mesh &mesh, const FaceFluxField &provisional_mass_flux)
{
    if (provisional_mass_flux.size() != mesh.face_count())
    {
        throw std::invalid_argument("Provisional mass-flux size must match the pressure-correction Mesh face count.");
    }
}

void validate_response_cardinality(const Mesh &mesh, const FacePressureResponseField &face_pressure_response)
{
    if (face_pressure_response.size() != mesh.face_count())
    {
        throw std::invalid_argument("Face pressure-response size must match the pressure-correction Mesh face count.");
    }
}

void validate_boundary_condition_cardinality(const Mesh &mesh,
                                             const PressureCorrectionBoundaryConditions &boundary_conditions)
{
    if (boundary_conditions.size() != mesh.boundary_groups().size())
    {
        throw std::invalid_argument(
            "Pressure-correction boundary condition count must match the pressure-correction Mesh boundary count.");
    }
}

[[noreturn]]
void throw_invalid_internal_face_value(const Index face_id, const std::string &reason)
{
    throw std::runtime_error("Pressure-correction assembly rejected face " + std::to_string(face_id) + ": " + reason);
}

[[noreturn]]
void throw_invalid_boundary_face_value(const Index face_id, const std::string &reason)
{
    throw std::runtime_error("Pressure-correction boundary assembly rejected face " + std::to_string(face_id) + ": " +
                             reason);
}

void validate_internal_face_values(const Mesh &mesh, const FaceFluxField &provisional_mass_flux,
                                   const FacePressureResponseField &face_pressure_response)
{
    const auto face_adjacencies{mesh.face_adjacencies()};
    const auto flux_values{provisional_mass_flux.values()};
    const auto response_values{face_pressure_response.values()};
    for (Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (face_adjacencies[face_id].is_boundary())
        {
            continue;
        }

        if (!std::isfinite(flux_values[face_id]))
        {
            throw_invalid_internal_face_value(face_id, "the provisional mass flux must be finite.");
        }
        if (!std::isfinite(response_values[face_id]) || !(response_values[face_id] > 0.0))
        {
            throw_invalid_internal_face_value(face_id,
                                              "the face pressure response must be finite and strictly positive.");
        }
    }
}

void validate_boundary_flux_values(const Mesh &mesh, const FaceFluxField &provisional_mass_flux)
{
    const auto face_adjacencies{mesh.face_adjacencies()};
    const auto flux_values{provisional_mass_flux.values()};
    for (Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (!face_adjacencies[face_id].is_boundary())
        {
            continue;
        }

        if (!std::isfinite(flux_values[face_id]))
        {
            throw_invalid_boundary_face_value(face_id, "the provisional mass flux must be finite.");
        }
    }
}

void validate_boundary_response_values(const Mesh &mesh,
                                       const PressureCorrectionBoundaryConditions &boundary_conditions,
                                       const FacePressureResponseField &face_pressure_response)
{
    const auto face_adjacencies{mesh.face_adjacencies()};
    const auto face_boundary_ids{mesh.face_boundary_ids()};
    const auto response_values{face_pressure_response.values()};
    for (Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (!face_adjacencies[face_id].is_boundary())
        {
            continue;
        }

        switch (boundary_conditions[face_boundary_ids[face_id]])
        {
        case PressureCorrectionBoundaryConditionType::FixedMassFlux:
            break;

        case PressureCorrectionBoundaryConditionType::FixedPressure:
            if (!std::isfinite(response_values[face_id]) || !(response_values[face_id] > 0.0))
            {
                throw_invalid_boundary_face_value(face_id,
                                                  "the face pressure response must be finite and strictly positive.");
            }
            break;
        }
    }
}

void add_internal_face_contributions_unchecked(const Mesh &mesh, const FaceFluxField &provisional_mass_flux,
                                               const FacePressureResponseField &face_pressure_response,
                                               ScalarLinearSystem &system) noexcept
{
    const auto face_adjacencies{mesh.face_adjacencies()};
    const auto flux_values{provisional_mass_flux.values()};
    const auto response_values{face_pressure_response.values()};
    auto diagonal{system.diagonal()};
    auto owner_neighbor_coefficients{system.owner_neighbor_coefficients()};
    auto neighbor_owner_coefficients{system.neighbor_owner_coefficients()};
    auto rhs{system.rhs()};
    for (Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        const FaceAdjacency &adjacency{face_adjacencies[face_id]};
        if (adjacency.is_boundary())
        {
            continue;
        }

        const double provisional_flux{flux_values[face_id]};
        const double pressure_response{response_values[face_id]};
        diagonal[adjacency.owner] += pressure_response;
        owner_neighbor_coefficients[face_id] -= pressure_response;
        diagonal[adjacency.neighbor] += pressure_response;
        neighbor_owner_coefficients[face_id] -= pressure_response;

        rhs[adjacency.owner] -= provisional_flux;
        rhs[adjacency.neighbor] += provisional_flux;
    }
}

void add_boundary_provisional_flux_rhs_unchecked(const Mesh &mesh, const FaceFluxField &provisional_mass_flux,
                                                 ScalarLinearSystem &system) noexcept
{
    const auto face_adjacencies{mesh.face_adjacencies()};
    const auto flux_values{provisional_mass_flux.values()};
    auto rhs{system.rhs()};
    for (Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        const FaceAdjacency &adjacency{face_adjacencies[face_id]};
        if (adjacency.is_boundary())
        {
            rhs[adjacency.owner] -= flux_values[face_id];
        }
    }
}

void add_boundary_pressure_response_unchecked(const Mesh &mesh,
                                              const PressureCorrectionBoundaryConditions &boundary_conditions,
                                              const FacePressureResponseField &face_pressure_response,
                                              ScalarLinearSystem &system) noexcept
{
    const auto face_adjacencies{mesh.face_adjacencies()};
    const auto face_boundary_ids{mesh.face_boundary_ids()};
    const auto response_values{face_pressure_response.values()};
    auto diagonal{system.diagonal()};
    for (Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        const FaceAdjacency &adjacency{face_adjacencies[face_id]};
        if (!adjacency.is_boundary())
        {
            continue;
        }

        switch (boundary_conditions[face_boundary_ids[face_id]])
        {
        case PressureCorrectionBoundaryConditionType::FixedMassFlux:
            break;

        case PressureCorrectionBoundaryConditionType::FixedPressure:
            diagonal[adjacency.owner] += response_values[face_id];
            break;
        }
    }
}

} // namespace

IncompressiblePressureCorrectionAssembler::IncompressiblePressureCorrectionAssembler(const Mesh &mesh) noexcept
    : mesh_(&mesh)
{
}

void IncompressiblePressureCorrectionAssembler::assemble(
    const FaceFluxField &provisional_mass_flux, const PressureCorrectionBoundaryConditions &boundary_conditions,
    const FacePressureResponseField &face_pressure_response, ScalarLinearSystem &system) const
{
    validate_flux_cardinality(*mesh_, provisional_mass_flux);
    validate_response_cardinality(*mesh_, face_pressure_response);
    validate_boundary_condition_cardinality(*mesh_, boundary_conditions);
    validate_system(*mesh_, system);
    validate_internal_face_values(*mesh_, provisional_mass_flux, face_pressure_response);
    validate_boundary_flux_values(*mesh_, provisional_mass_flux);
    validate_boundary_response_values(*mesh_, boundary_conditions, face_pressure_response);

    system.clear();
    add_internal_face_contributions_unchecked(*mesh_, provisional_mass_flux, face_pressure_response, system);
    add_boundary_provisional_flux_rhs_unchecked(*mesh_, provisional_mass_flux, system);
    add_boundary_pressure_response_unchecked(*mesh_, boundary_conditions, face_pressure_response, system);
}

void IncompressiblePressureCorrectionAssembler::add_internal_face_contributions(
    const FaceFluxField &provisional_mass_flux, const FacePressureResponseField &face_pressure_response,
    ScalarLinearSystem &system) const
{
    validate_flux_cardinality(*mesh_, provisional_mass_flux);
    validate_response_cardinality(*mesh_, face_pressure_response);
    validate_system(*mesh_, system);
    validate_internal_face_values(*mesh_, provisional_mass_flux, face_pressure_response);
    add_internal_face_contributions_unchecked(*mesh_, provisional_mass_flux, face_pressure_response, system);
}

void IncompressiblePressureCorrectionAssembler::add_boundary_provisional_flux_rhs(
    const FaceFluxField &provisional_mass_flux, ScalarLinearSystem &system) const
{
    validate_flux_cardinality(*mesh_, provisional_mass_flux);
    validate_system(*mesh_, system);
    validate_boundary_flux_values(*mesh_, provisional_mass_flux);
    add_boundary_provisional_flux_rhs_unchecked(*mesh_, provisional_mass_flux, system);
}

void IncompressiblePressureCorrectionAssembler::add_boundary_pressure_response(
    const PressureCorrectionBoundaryConditions &boundary_conditions,
    const FacePressureResponseField &face_pressure_response, ScalarLinearSystem &system) const
{
    validate_boundary_condition_cardinality(*mesh_, boundary_conditions);
    validate_response_cardinality(*mesh_, face_pressure_response);
    validate_system(*mesh_, system);
    validate_boundary_response_values(*mesh_, boundary_conditions, face_pressure_response);
    add_boundary_pressure_response_unchecked(*mesh_, boundary_conditions, face_pressure_response, system);
}

} // namespace cfd
