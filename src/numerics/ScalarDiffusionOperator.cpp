#include "cfd/numerics/ScalarDiffusionOperator.hpp"

#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVectorField.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/math/Point2.hpp"
#include "cfd/math/Vector2.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/Types.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace cfd
{

namespace
{

constexpr double projection_safety_factor{64.0};

[[noreturn]]
void throw_unusable_face_geometry(const Index face_id, const std::string &reason)
{
    throw std::runtime_error("Scalar diffusion operator rejected face " + std::to_string(face_id) + ": " + reason);
}

[[nodiscard]]
double dot(const Vector2 &first, const Vector2 &second) noexcept
{
    return first.x * second.x + first.y * second.y;
}

} // namespace

ScalarDiffusionOperator::ScalarDiffusionOperator(const Mesh &mesh, const double diffusivity)
    : mesh_(&mesh), diffusivity_(diffusivity)
{
    if (!std::isfinite(diffusivity_) || !(diffusivity_ > 0.0))
    {
        throw std::invalid_argument("Scalar diffusion diffusivity must be finite and strictly positive.");
    }

    const Index face_count{mesh_->face_count()};
    face_data_.resize(face_count);

    const auto face_adjacencies{mesh_->face_adjacencies()};
    const auto cell_centers{mesh_->cell_centers()};
    const auto face_centers{mesh_->face_centers()};
    const auto face_lengths{mesh_->face_lengths()};
    const auto face_area_vectors{mesh_->face_area_vectors()};

    constexpr double relative_projection_tolerance{projection_safety_factor * std::numeric_limits<double>::epsilon()};

    for (Index face_id = 0; face_id < face_count; ++face_id)
    {
        const FaceAdjacency &adjacency{face_adjacencies[face_id]};
        const Point2 &owner_center{cell_centers[adjacency.owner]};
        const Point2 &face_center{face_centers[face_id]};
        const Vector2 &area_vector{face_area_vectors[face_id]};
        const double face_length{face_lengths[face_id]};

        Vector2 center_displacement;
        if (adjacency.is_boundary())
        {
            center_displacement = {
                face_center.x - owner_center.x,
                face_center.y - owner_center.y,
            };
        }
        else
        {
            const Point2 &neighbor_center{cell_centers[adjacency.neighbor]};
            center_displacement = {
                neighbor_center.x - owner_center.x,
                neighbor_center.y - owner_center.y,
            };
        }

        const double center_distance{std::hypot(center_displacement.x, center_displacement.y)};
        const double projection_scale{face_length * center_distance};
        const double area_dot_displacement{dot(area_vector, center_displacement)};
        const double minimum_projection{relative_projection_tolerance * projection_scale};

        if (!std::isfinite(projection_scale) || !(projection_scale > 0.0) || !std::isfinite(area_dot_displacement) ||
            !(area_dot_displacement > minimum_projection))
        {
            throw_unusable_face_geometry(face_id, "Sf dot d is not a usable positive projection.");
        }

        const double area_vector_squared{dot(area_vector, area_vector)};
        const double beta{area_vector_squared / area_dot_displacement};
        const Vector2 correction{
            area_vector.x - beta * center_displacement.x,
            area_vector.y - beta * center_displacement.y,
        };

        FaceData &data{face_data_[face_id]};
        data.primary_coefficient = diffusivity_ * beta;
        data.correction_flux_vector = {
            -diffusivity_ * correction.x,
            -diffusivity_ * correction.y,
        };

        if (!std::isfinite(data.primary_coefficient) || !(data.primary_coefficient > 0.0) ||
            !std::isfinite(data.correction_flux_vector.x) || !std::isfinite(data.correction_flux_vector.y))
        {
            throw_unusable_face_geometry(face_id, "precomputed coefficients are non-finite.");
        }

        if (adjacency.is_boundary())
        {
            continue;
        }

        const Vector2 owner_to_face{
            face_center.x - owner_center.x,
            face_center.y - owner_center.y,
        };
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

        data.neighbor_gradient_weight = owner_to_face_projection / area_dot_displacement;
        if (!std::isfinite(data.neighbor_gradient_weight))
        {
            throw_unusable_face_geometry(face_id, "the face-gradient interpolation weight is non-finite.");
        }
    }
}

void ScalarDiffusionOperator::compute_flux_balance(const CellScalarField &field,
                                                   const ScalarBoundaryConditions &boundary_conditions,
                                                   const CellVectorField &gradient, CellScalarField &flux_balance) const
{
    const Index cell_count{mesh_->cell_count()};

    if (field.size() != cell_count)
    {
        throw std::invalid_argument("Scalar diffusion field size must match the mesh cell count.");
    }
    if (gradient.size() != cell_count)
    {
        throw std::invalid_argument("Scalar diffusion gradient size must match the mesh cell count.");
    }
    if (flux_balance.size() != cell_count)
    {
        throw std::invalid_argument("Scalar diffusion output size must match the mesh cell count.");
    }
    if (boundary_conditions.size() != mesh_->boundary_groups().size())
    {
        throw std::invalid_argument("Scalar diffusion boundary-condition count must match the mesh boundary count.");
    }
    if (&field == &flux_balance)
    {
        throw std::invalid_argument("Scalar diffusion output must not alias the input scalar field.");
    }

    const auto face_adjacencies{mesh_->face_adjacencies()};
    const auto face_boundary_ids{mesh_->face_boundary_ids()};
    const auto face_lengths{mesh_->face_lengths()};
    const auto field_values{field.values()};
    const auto gradient_values{gradient.values()};
    auto balance_values{flux_balance.values()};

    std::fill(balance_values.begin(), balance_values.end(), 0.0);

    for (Index face_id = 0; face_id < mesh_->face_count(); ++face_id)
    {
        const FaceAdjacency &adjacency{face_adjacencies[face_id]};
        const Index owner_id{adjacency.owner};
        const FaceData &data{face_data_[face_id]};

        double flux{};
        if (!adjacency.is_boundary())
        {
            const Index neighbor_id{adjacency.neighbor};
            const Vector2 &owner_gradient{gradient_values[owner_id]};
            const Vector2 &neighbor_gradient{gradient_values[neighbor_id]};
            const double lambda{data.neighbor_gradient_weight};
            const Vector2 face_gradient{
                owner_gradient.x + lambda * (neighbor_gradient.x - owner_gradient.x),
                owner_gradient.y + lambda * (neighbor_gradient.y - owner_gradient.y),
            };

            flux = data.primary_coefficient * (field_values[owner_id] - field_values[neighbor_id]) +
                   dot(face_gradient, data.correction_flux_vector);

            balance_values[owner_id] += flux;
            balance_values[neighbor_id] -= flux;
            continue;
        }

        const ScalarBoundaryCondition &condition{boundary_conditions[face_boundary_ids[face_id]]};
        switch (condition.type)
        {
        case ScalarBoundaryConditionType::Dirichlet:
            flux = data.primary_coefficient * (field_values[owner_id] - condition.value) +
                   dot(gradient_values[owner_id], data.correction_flux_vector);
            break;

        case ScalarBoundaryConditionType::Neumann:
            flux = -diffusivity_ * condition.value * face_lengths[face_id];
            break;
        }

        balance_values[owner_id] += flux;
    }
}

} // namespace cfd
