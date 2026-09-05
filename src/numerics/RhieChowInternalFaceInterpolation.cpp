#include "cfd/numerics/RhieChowInternalFaceInterpolation.hpp"

#include "cfd/field/CellMomentumPressureResponse.hpp"
#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVectorField.hpp"
#include "cfd/field/CellVelocityField.hpp"
#include "cfd/field/FaceFluxField.hpp"
#include "cfd/field/FacePressureResponseField.hpp"
#include "cfd/math/Point2.hpp"
#include "cfd/math/Vector2.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/Types.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace cfd
{
namespace
{

constexpr double projection_safety_factor{64.0};

struct InternalFaceResult
{
    double mass_flux{};
    double pressure_response{};
};

[[noreturn]]
void throw_unusable_face_geometry(const Index face_id, const std::string &reason)
{
    throw std::runtime_error("Rhie-Chow interpolation rejected face " + std::to_string(face_id) + ": " + reason);
}

[[noreturn]]
void throw_invalid_cell_value(const Index cell_id, const std::string &reason)
{
    throw std::runtime_error("Rhie-Chow interpolation rejected cell " + std::to_string(cell_id) + ": " + reason);
}

[[noreturn]]
void throw_invalid_face_result(const Index face_id, const std::string &reason)
{
    throw std::runtime_error("Rhie-Chow interpolation rejected face " + std::to_string(face_id) + ": " + reason);
}

[[nodiscard]]
double dot(const Vector2 &first, const Vector2 &second) noexcept
{
    return first.x * second.x + first.y * second.y;
}

} // namespace

RhieChowInternalFaceInterpolation::RhieChowInternalFaceInterpolation(const Mesh &mesh, const double density)
    : mesh_(&mesh), density_(density)
{
    if (!std::isfinite(density_) || !(density_ > 0.0))
    {
        throw std::invalid_argument("Rhie-Chow density must be finite and strictly positive.");
    }

    internal_face_geometry_.resize(mesh_->face_count());
    const auto face_adjacencies{mesh_->face_adjacencies()};
    const auto cell_centers{mesh_->cell_centers()};
    const auto face_centers{mesh_->face_centers()};
    const auto face_lengths{mesh_->face_lengths()};
    const auto face_area_vectors{mesh_->face_area_vectors()};
    constexpr double relative_projection_tolerance{projection_safety_factor * std::numeric_limits<double>::epsilon()};

    for (Index face_id = 0; face_id < mesh_->face_count(); ++face_id)
    {
        const FaceAdjacency &adjacency{face_adjacencies[face_id]};
        if (adjacency.is_boundary())
        {
            continue;
        }

        const Point2 &owner_center{cell_centers[adjacency.owner]};
        const Point2 &neighbor_center{cell_centers[adjacency.neighbor]};
        const Point2 &face_center{face_centers[face_id]};
        const Vector2 &area_vector{face_area_vectors[face_id]};
        const Vector2 center_displacement{
            neighbor_center.x - owner_center.x,
            neighbor_center.y - owner_center.y,
        };
        const Vector2 owner_to_face{
            face_center.x - owner_center.x,
            face_center.y - owner_center.y,
        };
        const double face_length{face_lengths[face_id]};
        const double center_distance{std::hypot(center_displacement.x, center_displacement.y)};
        const double projection_scale{face_length * center_distance};
        const double area_dot_displacement{dot(area_vector, center_displacement)};
        const double minimum_projection{relative_projection_tolerance * projection_scale};

        if (!std::isfinite(projection_scale) || !(projection_scale > 0.0) || !std::isfinite(area_dot_displacement) ||
            !(area_dot_displacement > minimum_projection))
        {
            throw_unusable_face_geometry(face_id, "Sf dot d is not a usable positive projection.");
        }

        const double owner_to_face_projection{dot(area_vector, owner_to_face)};
        const double owner_to_face_distance{std::hypot(owner_to_face.x, owner_to_face.y)};
        const double intersection_tolerance{relative_projection_tolerance * face_length *
                                            (center_distance + owner_to_face_distance)};
        if (!std::isfinite(owner_to_face_projection) || !std::isfinite(intersection_tolerance) ||
            owner_to_face_projection < -intersection_tolerance ||
            owner_to_face_projection > area_dot_displacement + intersection_tolerance)
        {
            throw_unusable_face_geometry(face_id,
                                         "the face-line intersection lies outside the owner-neighbor segment.");
        }

        const double interpolation_weight{owner_to_face_projection / area_dot_displacement};
        const double inverse_area_dot_displacement{1.0 / area_dot_displacement};
        if (!std::isfinite(interpolation_weight) || !std::isfinite(inverse_area_dot_displacement))
        {
            throw_unusable_face_geometry(face_id, "the cached interpolation geometry is non-finite.");
        }
        internal_face_geometry_[face_id] = {interpolation_weight, inverse_area_dot_displacement};
    }
}

void RhieChowInternalFaceInterpolation::update_internal_faces(const CellVelocityField &velocity,
                                                              const CellScalarField &pressure,
                                                              const CellVectorField &pressure_gradient,
                                                              const CellMomentumPressureResponse &momentum_response,
                                                              FaceFluxField &mass_flux,
                                                              FacePressureResponseField &face_pressure_response) const
{
    const Index cell_count{mesh_->cell_count()};
    const Index face_count{mesh_->face_count()};

    if (velocity.size() != cell_count || velocity.u().size() != cell_count || velocity.v().size() != cell_count)
    {
        throw std::invalid_argument("Rhie-Chow velocity size must match the mesh cell count.");
    }
    if (pressure.size() != cell_count)
    {
        throw std::invalid_argument("Rhie-Chow pressure size must match the mesh cell count.");
    }
    if (pressure_gradient.size() != cell_count)
    {
        throw std::invalid_argument("Rhie-Chow pressure-gradient size must match the mesh cell count.");
    }
    if (momentum_response.size() != cell_count || momentum_response.u().size() != cell_count ||
        momentum_response.v().size() != cell_count)
    {
        throw std::invalid_argument("Rhie-Chow momentum-response size must match the mesh cell count.");
    }
    if (mass_flux.size() != face_count)
    {
        throw std::invalid_argument("Rhie-Chow mass-flux size must match the mesh face count.");
    }
    if (face_pressure_response.size() != face_count)
    {
        throw std::invalid_argument("Rhie-Chow face pressure-response size must match the mesh face count.");
    }

    const auto u_values{velocity.u().values()};
    const auto v_values{velocity.v().values()};
    const auto pressure_values{pressure.values()};
    const auto gradient_values{pressure_gradient.values()};
    const auto u_response_values{momentum_response.u().values()};
    const auto v_response_values{momentum_response.v().values()};

    for (Index cell_id = 0; cell_id < cell_count; ++cell_id)
    {
        if (!std::isfinite(u_values[cell_id]) || !std::isfinite(v_values[cell_id]))
        {
            throw_invalid_cell_value(cell_id, "velocity components must be finite.");
        }
        if (!std::isfinite(pressure_values[cell_id]))
        {
            throw_invalid_cell_value(cell_id, "pressure must be finite.");
        }
        if (!std::isfinite(gradient_values[cell_id].x) || !std::isfinite(gradient_values[cell_id].y))
        {
            throw_invalid_cell_value(cell_id, "pressure-gradient components must be finite.");
        }
        if (!std::isfinite(u_response_values[cell_id]) || !(u_response_values[cell_id] > 0.0) ||
            !std::isfinite(v_response_values[cell_id]) || !(v_response_values[cell_id] > 0.0))
        {
            throw_invalid_cell_value(cell_id, "momentum responses must be finite and strictly positive.");
        }
    }

    const auto face_adjacencies{mesh_->face_adjacencies()};
    const auto cell_centers{mesh_->cell_centers()};
    const auto area_vectors{mesh_->face_area_vectors()};

    const auto evaluate_face = [&](const Index face_id) {
        const FaceAdjacency &adjacency{face_adjacencies[face_id]};
        const Index owner_id{adjacency.owner};
        const Index neighbor_id{adjacency.neighbor};
        const double lambda{internal_face_geometry_[face_id].interpolation_weight};
        const double owner_weight{1.0 - lambda};
        const Vector2 &area_vector{area_vectors[face_id]};
        const Point2 &owner_center{cell_centers[owner_id]};
        const Point2 &neighbor_center{cell_centers[neighbor_id]};
        const Vector2 center_displacement{
            neighbor_center.x - owner_center.x,
            neighbor_center.y - owner_center.y,
        };

        const double face_u_response{owner_weight * u_response_values[owner_id] +
                                     lambda * u_response_values[neighbor_id]};
        const double face_v_response{owner_weight * v_response_values[owner_id] +
                                     lambda * v_response_values[neighbor_id]};
        if (!std::isfinite(face_u_response) || !std::isfinite(face_v_response))
        {
            throw_invalid_face_result(face_id, "interpolated momentum responses are non-finite.");
        }

        const Vector2 response_area_vector{
            face_u_response * area_vector.x,
            face_v_response * area_vector.y,
        };
        const double pressure_coefficient{dot(area_vector, response_area_vector) *
                                          internal_face_geometry_[face_id].inverse_area_dot_displacement};
        const double integrated_pressure_response{density_ * pressure_coefficient};
        const Vector2 tangential_response{
            response_area_vector.x - pressure_coefficient * center_displacement.x,
            response_area_vector.y - pressure_coefficient * center_displacement.y,
        };
        if (!std::isfinite(response_area_vector.x) || !std::isfinite(response_area_vector.y) ||
            !std::isfinite(pressure_coefficient) || !(pressure_coefficient > 0.0) ||
            !std::isfinite(integrated_pressure_response) || !(integrated_pressure_response > 0.0) ||
            !std::isfinite(tangential_response.x) || !std::isfinite(tangential_response.y))
        {
            throw_invalid_face_result(face_id, "the pressure-response decomposition is unusable.");
        }

        const Vector2 &owner_gradient{gradient_values[owner_id]};
        const Vector2 &neighbor_gradient{gradient_values[neighbor_id]};
        const Vector2 face_gradient{
            owner_weight * owner_gradient.x + lambda * neighbor_gradient.x,
            owner_weight * owner_gradient.y + lambda * neighbor_gradient.y,
        };
        const Vector2 owner_pressure_free_velocity{
            u_values[owner_id] + u_response_values[owner_id] * owner_gradient.x,
            v_values[owner_id] + v_response_values[owner_id] * owner_gradient.y,
        };
        const Vector2 neighbor_pressure_free_velocity{
            u_values[neighbor_id] + u_response_values[neighbor_id] * neighbor_gradient.x,
            v_values[neighbor_id] + v_response_values[neighbor_id] * neighbor_gradient.y,
        };
        const Vector2 face_pressure_free_velocity{
            owner_weight * owner_pressure_free_velocity.x + lambda * neighbor_pressure_free_velocity.x,
            owner_weight * owner_pressure_free_velocity.y + lambda * neighbor_pressure_free_velocity.y,
        };
        const double pressure_jump{pressure_values[neighbor_id] - pressure_values[owner_id]};
        const double flux_without_density{dot(face_pressure_free_velocity, area_vector) -
                                          pressure_coefficient * pressure_jump -
                                          dot(face_gradient, tangential_response)};
        const double integrated_mass_flux{density_ * flux_without_density};
        if (!std::isfinite(face_gradient.x) || !std::isfinite(face_gradient.y) ||
            !std::isfinite(owner_pressure_free_velocity.x) || !std::isfinite(owner_pressure_free_velocity.y) ||
            !std::isfinite(neighbor_pressure_free_velocity.x) || !std::isfinite(neighbor_pressure_free_velocity.y) ||
            !std::isfinite(face_pressure_free_velocity.x) || !std::isfinite(face_pressure_free_velocity.y) ||
            !std::isfinite(pressure_jump) || !std::isfinite(flux_without_density) ||
            !std::isfinite(integrated_mass_flux))
        {
            throw_invalid_face_result(face_id, "the interpolated mass flux is non-finite.");
        }

        return InternalFaceResult{integrated_mass_flux, integrated_pressure_response};
    };

    // Validate every internal result before either caller-owned output changes.
    for (Index face_id = 0; face_id < face_count; ++face_id)
    {
        if (!face_adjacencies[face_id].is_boundary())
        {
            static_cast<void>(evaluate_face(face_id));
        }
    }

    auto mass_flux_values{mass_flux.values()};
    auto face_response_values{face_pressure_response.values()};
    for (Index face_id = 0; face_id < face_count; ++face_id)
    {
        if (face_adjacencies[face_id].is_boundary())
        {
            continue;
        }

        const InternalFaceResult result{evaluate_face(face_id)};
        mass_flux_values[face_id] = result.mass_flux;
        face_response_values[face_id] = result.pressure_response;
    }
}

} // namespace cfd
