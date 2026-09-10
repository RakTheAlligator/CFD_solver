#include "cfd/numerics/RhieChowBoundaryFaceInterpolation.hpp"

#include "cfd/field/CellMomentumPressureResponse.hpp"
#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVectorField.hpp"
#include "cfd/field/CellVelocityField.hpp"
#include "cfd/field/FaceFluxField.hpp"
#include "cfd/field/FacePressureResponseField.hpp"
#include "cfd/field/PressureCorrectionBoundaryConditions.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/math/Point2.hpp"
#include "cfd/math/Vector2.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/Types.hpp"
#include "cfd/numerics/FiniteVolumeFaceGeometry.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace cfd
{
namespace
{

struct BoundaryFaceResult
{
    double mass_flux{};
    double pressure_response{};
};

[[noreturn]]
void throw_invalid_boundary_cell_value(const Index cell_id, const std::string &reason)
{
    throw std::runtime_error("Boundary Rhie-Chow interpolation rejected cell " + std::to_string(cell_id) + ": " +
                             reason);
}

[[noreturn]]
void throw_invalid_boundary_face_result(const Index face_id, const std::string &reason)
{
    throw std::runtime_error("Boundary Rhie-Chow interpolation rejected face " + std::to_string(face_id) + ": " +
                             reason);
}

[[nodiscard]]
double dot(const Vector2 &first, const Vector2 &second) noexcept
{
    return first.x * second.x + first.y * second.y;
}

} // namespace

RhieChowBoundaryFaceInterpolation::RhieChowBoundaryFaceInterpolation(const Mesh &mesh, const double density)
    : mesh_(&mesh), density_(density)
{
    if (!std::isfinite(density_) || !(density_ > 0.0))
    {
        throw std::invalid_argument("Boundary Rhie-Chow density must be finite and strictly positive.");
    }
}

void RhieChowBoundaryFaceInterpolation::update_fixed_pressure_boundaries(
    const CellVelocityField &velocity, const CellScalarField &pressure, const CellVectorField &pressure_gradient,
    const CellMomentumPressureResponse &momentum_response, const ScalarBoundaryConditions &pressure_boundary_conditions,
    const PressureCorrectionBoundaryConditions &pressure_correction_boundary_conditions, FaceFluxField &mass_flux,
    FacePressureResponseField &face_pressure_response) const
{
    const Index cell_count{mesh_->cell_count()};
    const Index face_count{mesh_->face_count()};
    const Index boundary_count{mesh_->boundary_groups().size()};

    if (velocity.size() != cell_count || velocity.u().size() != cell_count || velocity.v().size() != cell_count)
    {
        throw std::invalid_argument("Boundary Rhie-Chow velocity size must match the mesh cell count.");
    }
    if (pressure.size() != cell_count)
    {
        throw std::invalid_argument("Boundary Rhie-Chow pressure size must match the mesh cell count.");
    }
    if (pressure_gradient.size() != cell_count)
    {
        throw std::invalid_argument("Boundary Rhie-Chow pressure-gradient size must match the mesh cell count.");
    }
    if (momentum_response.size() != cell_count || momentum_response.u().size() != cell_count ||
        momentum_response.v().size() != cell_count)
    {
        throw std::invalid_argument("Boundary Rhie-Chow momentum-response size must match the mesh cell count.");
    }
    if (pressure_boundary_conditions.size() != boundary_count)
    {
        throw std::invalid_argument("Pressure boundary-condition count must match the mesh boundary count.");
    }
    if (pressure_correction_boundary_conditions.size() != boundary_count)
    {
        throw std::invalid_argument("Pressure-correction boundary-condition count must match the mesh boundary count.");
    }
    if (mass_flux.size() != face_count)
    {
        throw std::invalid_argument("Boundary Rhie-Chow mass-flux size must match the mesh face count.");
    }
    if (face_pressure_response.size() != face_count)
    {
        throw std::invalid_argument("Boundary Rhie-Chow face pressure-response size must match the mesh face count.");
    }

    const auto u_values{velocity.u().values()};
    const auto v_values{velocity.v().values()};
    const auto pressure_values{pressure.values()};
    const auto gradient_values{pressure_gradient.values()};
    const auto u_response_values{momentum_response.u().values()};
    const auto v_response_values{momentum_response.v().values()};
    const auto face_adjacencies{mesh_->face_adjacencies()};
    const auto face_boundary_ids{mesh_->face_boundary_ids()};
    const auto cell_centers{mesh_->cell_centers()};
    const auto face_centers{mesh_->face_centers()};
    const auto face_lengths{mesh_->face_lengths()};
    const auto area_vectors{mesh_->face_area_vectors()};
    const auto evaluate_fixed_pressure_face = [&](const Index face_id) {
        const FaceAdjacency &adjacency{face_adjacencies[face_id]};
        const Index owner_id{adjacency.owner};
        const BoundaryId boundary_id{face_boundary_ids[face_id]};
        const ScalarBoundaryCondition &pressure_condition{pressure_boundary_conditions[boundary_id]};
        if (pressure_condition.type != ScalarBoundaryConditionType::Dirichlet)
        {
            throw std::invalid_argument("FixedPressure pressure correction requires a physical Dirichlet pressure "
                                        "condition on face " +
                                        std::to_string(face_id) + ".");
        }

        const double owner_u{u_values[owner_id]};
        const double owner_v{v_values[owner_id]};
        const double owner_pressure{pressure_values[owner_id]};
        const Vector2 &owner_gradient{gradient_values[owner_id]};
        const double owner_u_response{u_response_values[owner_id]};
        const double owner_v_response{v_response_values[owner_id]};
        const double boundary_pressure{pressure_condition.value};
        if (!std::isfinite(owner_u) || !std::isfinite(owner_v))
        {
            throw_invalid_boundary_cell_value(owner_id, "velocity components must be finite.");
        }
        if (!std::isfinite(owner_pressure))
        {
            throw_invalid_boundary_cell_value(owner_id, "pressure must be finite.");
        }
        if (!std::isfinite(owner_gradient.x) || !std::isfinite(owner_gradient.y))
        {
            throw_invalid_boundary_cell_value(owner_id, "pressure-gradient components must be finite.");
        }
        if (!std::isfinite(owner_u_response) || !(owner_u_response > 0.0) || !std::isfinite(owner_v_response) ||
            !(owner_v_response > 0.0))
        {
            throw_invalid_boundary_cell_value(owner_id, "momentum responses must be finite and strictly positive.");
        }
        if (!std::isfinite(boundary_pressure))
        {
            throw_invalid_boundary_face_result(face_id, "prescribed pressure must be finite.");
        }

        const Point2 &owner_center{cell_centers[owner_id]};
        const Point2 &face_center{face_centers[face_id]};
        const Vector2 &area_vector{area_vectors[face_id]};
        const double face_length{face_lengths[face_id]};
        const FaceProjectionGeometry geometry{
            compute_face_projection_geometry(face_id, owner_center, face_center, area_vector, face_length)};
        const Vector2 &displacement{geometry.displacement};
        const double area_dot_displacement{geometry.area_dot_displacement};

        const Vector2 response_area_vector{
            owner_u_response * area_vector.x,
            owner_v_response * area_vector.y,
        };
        const double pressure_coefficient{dot(area_vector, response_area_vector) / area_dot_displacement};
        const Vector2 tangential_response{
            response_area_vector.x - pressure_coefficient * displacement.x,
            response_area_vector.y - pressure_coefficient * displacement.y,
        };
        const double integrated_pressure_response{density_ * pressure_coefficient};
        if (!std::isfinite(response_area_vector.x) || !std::isfinite(response_area_vector.y) ||
            !std::isfinite(pressure_coefficient) || !(pressure_coefficient > 0.0) ||
            !std::isfinite(tangential_response.x) || !std::isfinite(tangential_response.y) ||
            !std::isfinite(integrated_pressure_response) || !(integrated_pressure_response > 0.0))
        {
            throw_invalid_boundary_face_result(face_id, "the pressure-response decomposition is unusable.");
        }

        const Vector2 pressure_free_velocity{
            owner_u + owner_u_response * owner_gradient.x,
            owner_v + owner_v_response * owner_gradient.y,
        };
        const double pressure_difference{boundary_pressure - owner_pressure};
        const double flux_without_density{dot(pressure_free_velocity, area_vector) -
                                          pressure_coefficient * pressure_difference -
                                          dot(owner_gradient, tangential_response)};
        const double integrated_mass_flux{density_ * flux_without_density};
        if (!std::isfinite(pressure_free_velocity.x) || !std::isfinite(pressure_free_velocity.y) ||
            !std::isfinite(pressure_difference) || !std::isfinite(flux_without_density) ||
            !std::isfinite(integrated_mass_flux))
        {
            throw_invalid_boundary_face_result(face_id, "the interpolated mass flux is non-finite.");
        }

        return BoundaryFaceResult{integrated_mass_flux, integrated_pressure_response};
    };

    // Validate every FixedPressure result before either caller-owned output changes.
    for (Index face_id = 0; face_id < face_count; ++face_id)
    {
        if (!face_adjacencies[face_id].is_boundary())
        {
            continue;
        }

        switch (pressure_correction_boundary_conditions[face_boundary_ids[face_id]])
        {
        case PressureCorrectionBoundaryConditionType::FixedMassFlux:
            break;

        case PressureCorrectionBoundaryConditionType::FixedPressure:
            static_cast<void>(evaluate_fixed_pressure_face(face_id));
            break;
        }
    }

    auto mass_flux_values{mass_flux.values()};
    auto face_response_values{face_pressure_response.values()};
    for (Index face_id = 0; face_id < face_count; ++face_id)
    {
        if (!face_adjacencies[face_id].is_boundary())
        {
            continue;
        }

        switch (pressure_correction_boundary_conditions[face_boundary_ids[face_id]])
        {
        case PressureCorrectionBoundaryConditionType::FixedMassFlux:
            break;

        case PressureCorrectionBoundaryConditionType::FixedPressure: {
            const BoundaryFaceResult result{evaluate_fixed_pressure_face(face_id)};
            mass_flux_values[face_id] = result.mass_flux;
            face_response_values[face_id] = result.pressure_response;
            break;
        }
        }
    }
}

} // namespace cfd
