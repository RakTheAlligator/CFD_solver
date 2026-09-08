#include "cfd/numerics/IncompressiblePressureVelocityCorrection.hpp"

#include "cfd/field/CellMomentumPressureResponse.hpp"
#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVectorField.hpp"
#include "cfd/field/CellVelocityField.hpp"
#include "cfd/field/FaceFluxField.hpp"
#include "cfd/field/FacePressureResponseField.hpp"
#include "cfd/field/PressureCorrectionBoundaryConditions.hpp"
#include "cfd/math/Vector2.hpp"
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

[[noreturn]]
void throw_invalid_face_correction(const Index face_id, const std::string &reason)
{
    throw std::runtime_error("Pressure-velocity correction rejected face " + std::to_string(face_id) + ": " + reason);
}

[[noreturn]]
void throw_invalid_cell_correction(const Index cell_id, const std::string &reason)
{
    throw std::runtime_error("Pressure-velocity correction rejected cell " + std::to_string(cell_id) + ": " + reason);
}

void validate_positive_finite_response(const double response, const Index face_id)
{
    if (!std::isfinite(response) || !(response > 0.0))
    {
        throw_invalid_face_correction(face_id, "the face pressure response must be finite and strictly positive.");
    }
}

} // namespace

IncompressiblePressureVelocityCorrection::IncompressiblePressureVelocityCorrection(const Mesh &mesh) noexcept
    : mesh_(&mesh)
{
}

void IncompressiblePressureVelocityCorrection::correct_face_mass_flux(
    const CellScalarField &pressure_correction, const PressureCorrectionBoundaryConditions &boundary_conditions,
    const FacePressureResponseField &face_pressure_response, FaceFluxField &mass_flux) const
{
    const Index cell_count{mesh_->cell_count()};
    const Index face_count{mesh_->face_count()};
    if (pressure_correction.size() != cell_count)
    {
        throw std::invalid_argument("Pressure-correction size must match the correction Mesh cell count.");
    }
    if (boundary_conditions.size() != mesh_->boundary_groups().size())
    {
        throw std::invalid_argument(
            "Pressure-correction boundary condition count must match the correction Mesh boundary count.");
    }
    if (face_pressure_response.size() != face_count)
    {
        throw std::invalid_argument("Face pressure-response size must match the correction Mesh face count.");
    }
    if (mass_flux.size() != face_count)
    {
        throw std::invalid_argument("Mass-flux size must match the correction Mesh face count.");
    }

    const auto pressure_correction_values{pressure_correction.values()};
    const auto response_values{face_pressure_response.values()};
    const auto mass_flux_values{mass_flux.values()};
    const auto face_adjacencies{mesh_->face_adjacencies()};
    const auto face_boundary_ids{mesh_->face_boundary_ids()};

    const auto evaluate_internal_face = [&](const Index face_id) {
        const FaceAdjacency &adjacency{face_adjacencies[face_id]};
        const double owner_pressure_correction{pressure_correction_values[adjacency.owner]};
        const double neighbor_pressure_correction{pressure_correction_values[adjacency.neighbor]};
        const double pressure_response{response_values[face_id]};
        const double current_mass_flux{mass_flux_values[face_id]};
        if (!std::isfinite(owner_pressure_correction) || !std::isfinite(neighbor_pressure_correction))
        {
            throw_invalid_face_correction(face_id, "the used cell pressure corrections must be finite.");
        }
        validate_positive_finite_response(pressure_response, face_id);
        if (!std::isfinite(current_mass_flux))
        {
            throw_invalid_face_correction(face_id, "the current mass flux must be finite.");
        }

        const double corrected_mass_flux{
            current_mass_flux + pressure_response * (owner_pressure_correction - neighbor_pressure_correction)};
        if (!std::isfinite(corrected_mass_flux))
        {
            throw_invalid_face_correction(face_id, "the corrected mass flux must be finite.");
        }
        return corrected_mass_flux;
    };

    const auto evaluate_fixed_pressure_face = [&](const Index face_id) {
        const Index owner_id{face_adjacencies[face_id].owner};
        const double owner_pressure_correction{pressure_correction_values[owner_id]};
        const double pressure_response{response_values[face_id]};
        const double current_mass_flux{mass_flux_values[face_id]};
        if (!std::isfinite(owner_pressure_correction))
        {
            throw_invalid_face_correction(face_id, "the owner pressure correction must be finite.");
        }
        validate_positive_finite_response(pressure_response, face_id);
        if (!std::isfinite(current_mass_flux))
        {
            throw_invalid_face_correction(face_id, "the current mass flux must be finite.");
        }

        const double corrected_mass_flux{current_mass_flux + pressure_response * owner_pressure_correction};
        if (!std::isfinite(corrected_mass_flux))
        {
            throw_invalid_face_correction(face_id, "the corrected mass flux must be finite.");
        }
        return corrected_mass_flux;
    };

    // Validate every corrected face before modifying caller-owned flux storage.
    for (Index face_id = 0; face_id < face_count; ++face_id)
    {
        if (!face_adjacencies[face_id].is_boundary())
        {
            static_cast<void>(evaluate_internal_face(face_id));
            continue;
        }

        switch (boundary_conditions[face_boundary_ids[face_id]])
        {
        case PressureCorrectionBoundaryConditionType::FixedMassFlux:
            break;

        case PressureCorrectionBoundaryConditionType::FixedPressure:
            static_cast<void>(evaluate_fixed_pressure_face(face_id));
            break;
        }
    }

    for (Index face_id = 0; face_id < face_count; ++face_id)
    {
        if (!face_adjacencies[face_id].is_boundary())
        {
            mass_flux_values[face_id] = evaluate_internal_face(face_id);
            continue;
        }

        switch (boundary_conditions[face_boundary_ids[face_id]])
        {
        case PressureCorrectionBoundaryConditionType::FixedMassFlux:
            break;

        case PressureCorrectionBoundaryConditionType::FixedPressure:
            mass_flux_values[face_id] = evaluate_fixed_pressure_face(face_id);
            break;
        }
    }
}

void IncompressiblePressureVelocityCorrection::correct_pressure(const CellScalarField &pressure_correction,
                                                                const double pressure_relaxation_factor,
                                                                CellScalarField &pressure) const
{
    const Index cell_count{mesh_->cell_count()};
    if (pressure_correction.size() != cell_count || pressure.size() != cell_count)
    {
        throw std::invalid_argument("Pressure fields must match the correction Mesh cell count.");
    }
    if (&pressure_correction == &pressure)
    {
        throw std::invalid_argument("Pressure and pressure-correction fields must not alias.");
    }
    if (!std::isfinite(pressure_relaxation_factor) || !(pressure_relaxation_factor > 0.0) ||
        !(pressure_relaxation_factor <= 1.0))
    {
        throw std::invalid_argument("Pressure relaxation factor must be finite and in (0, 1].");
    }

    const auto pressure_correction_values{pressure_correction.values()};
    const auto pressure_values{pressure.values()};
    const auto evaluate_cell = [&](const Index cell_id) {
        const double old_pressure{pressure_values[cell_id]};
        const double correction{pressure_correction_values[cell_id]};
        if (!std::isfinite(old_pressure) || !std::isfinite(correction))
        {
            throw_invalid_cell_correction(cell_id, "pressure and pressure correction must be finite.");
        }

        const double corrected_pressure{old_pressure + pressure_relaxation_factor * correction};
        if (!std::isfinite(corrected_pressure))
        {
            throw_invalid_cell_correction(cell_id, "the corrected pressure must be finite.");
        }
        return corrected_pressure;
    };

    for (Index cell_id = 0; cell_id < cell_count; ++cell_id)
    {
        static_cast<void>(evaluate_cell(cell_id));
    }
    for (Index cell_id = 0; cell_id < cell_count; ++cell_id)
    {
        pressure_values[cell_id] = evaluate_cell(cell_id);
    }
}

void IncompressiblePressureVelocityCorrection::correct_velocity(const CellVectorField &pressure_correction_gradient,
                                                                const CellMomentumPressureResponse &momentum_response,
                                                                CellVelocityField &velocity) const
{
    const Index cell_count{mesh_->cell_count()};
    if (pressure_correction_gradient.size() != cell_count)
    {
        throw std::invalid_argument("Pressure-correction gradient size must match the correction Mesh cell count.");
    }
    if (momentum_response.size() != cell_count || momentum_response.u().size() != cell_count ||
        momentum_response.v().size() != cell_count)
    {
        throw std::invalid_argument("Momentum-response size must match the correction Mesh cell count.");
    }
    if (velocity.size() != cell_count || velocity.u().size() != cell_count || velocity.v().size() != cell_count)
    {
        throw std::invalid_argument("Velocity size must match the correction Mesh cell count.");
    }

    const auto gradient_values{pressure_correction_gradient.values()};
    const auto u_response_values{momentum_response.u().values()};
    const auto v_response_values{momentum_response.v().values()};
    const auto u_values{velocity.u().values()};
    const auto v_values{velocity.v().values()};
    const auto evaluate_cell = [&](const Index cell_id) {
        const Vector2 &gradient{gradient_values[cell_id]};
        const double u_response{u_response_values[cell_id]};
        const double v_response{v_response_values[cell_id]};
        const double old_u{u_values[cell_id]};
        const double old_v{v_values[cell_id]};
        if (!std::isfinite(gradient.x) || !std::isfinite(gradient.y))
        {
            throw_invalid_cell_correction(cell_id, "pressure-correction gradient components must be finite.");
        }
        if (!std::isfinite(u_response) || !(u_response > 0.0) || !std::isfinite(v_response) || !(v_response > 0.0))
        {
            throw_invalid_cell_correction(cell_id, "momentum responses must be finite and strictly positive.");
        }
        if (!std::isfinite(old_u) || !std::isfinite(old_v))
        {
            throw_invalid_cell_correction(cell_id, "velocity components must be finite.");
        }

        const Vector2 corrected_velocity{
            old_u - u_response * gradient.x,
            old_v - v_response * gradient.y,
        };
        if (!std::isfinite(corrected_velocity.x) || !std::isfinite(corrected_velocity.y))
        {
            throw_invalid_cell_correction(cell_id, "corrected velocity components must be finite.");
        }
        return corrected_velocity;
    };

    for (Index cell_id = 0; cell_id < cell_count; ++cell_id)
    {
        static_cast<void>(evaluate_cell(cell_id));
    }
    for (Index cell_id = 0; cell_id < cell_count; ++cell_id)
    {
        const Vector2 corrected_velocity{evaluate_cell(cell_id)};
        u_values[cell_id] = corrected_velocity.x;
        v_values[cell_id] = corrected_velocity.y;
    }
}

} // namespace cfd
