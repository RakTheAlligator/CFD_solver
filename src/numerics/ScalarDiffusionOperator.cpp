#include "cfd/numerics/ScalarDiffusionOperator.hpp"

#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVectorField.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/linear_algebra/ScalarLinearSystem.hpp"
#include "cfd/math/Point2.hpp"
#include "cfd/math/Vector2.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/Types.hpp"
#include "cfd/numerics/FiniteVolumeFaceGeometry.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace cfd
{

namespace
{

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

void validate_boundary_condition_count(const Mesh &mesh, const ScalarBoundaryConditions &boundary_conditions)
{
    if (boundary_conditions.size() != mesh.boundary_groups().size())
    {
        throw std::invalid_argument("Scalar diffusion boundary-condition count must match the mesh boundary count.");
    }
}

void validate_rhs_size(const Mesh &mesh, const std::span<double> rhs)
{
    if (rhs.size() != mesh.cell_count())
    {
        throw std::invalid_argument("Scalar diffusion RHS size must match the mesh cell count.");
    }
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

    for (Index face_id = 0; face_id < face_count; ++face_id)
    {
        const FaceAdjacency &adjacency{face_adjacencies[face_id]};
        const Point2 &owner_center{cell_centers[adjacency.owner]};
        const Point2 &face_center{face_centers[face_id]};
        const Vector2 &area_vector{face_area_vectors[face_id]};
        const double face_length{face_lengths[face_id]};

        Vector2 center_displacement;
        double area_dot_displacement{};
        double neighbor_gradient_weight{};
        if (adjacency.is_boundary())
        {
            const FaceProjectionGeometry geometry{
                compute_face_projection_geometry(face_id, owner_center, face_center, area_vector, face_length)};
            center_displacement = geometry.displacement;
            area_dot_displacement = geometry.area_dot_displacement;
        }
        else
        {
            const Point2 &neighbor_center{cell_centers[adjacency.neighbor]};
            const InternalFaceInterpolationGeometry geometry{compute_internal_face_interpolation_geometry(
                face_id, owner_center, neighbor_center, face_center, area_vector, face_length)};
            center_displacement = geometry.owner_to_neighbor;
            area_dot_displacement = geometry.area_dot_owner_to_neighbor;
            neighbor_gradient_weight = geometry.neighbor_weight;
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

        data.neighbor_gradient_weight = neighbor_gradient_weight;
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

void ScalarDiffusionOperator::add_matrix_contributions(const ScalarBoundaryConditions &boundary_conditions,
                                                       ScalarLinearSystem &system) const
{
    validate_boundary_condition_count(*mesh_, boundary_conditions);
    if (&system.mesh() != mesh_ || system.cell_count() != mesh_->cell_count() ||
        system.face_count() != mesh_->face_count())
    {
        throw std::invalid_argument("Scalar diffusion system must reference the operator Mesh.");
    }

    auto diagonal{system.diagonal()};
    auto owner_neighbor_coefficients{system.owner_neighbor_coefficients()};
    auto neighbor_owner_coefficients{system.neighbor_owner_coefficients()};
    const auto face_adjacencies{mesh_->face_adjacencies()};
    const auto face_boundary_ids{mesh_->face_boundary_ids()};

    for (Index face_id = 0; face_id < mesh_->face_count(); ++face_id)
    {
        const FaceAdjacency &adjacency{face_adjacencies[face_id]};
        const double coefficient{face_data_[face_id].primary_coefficient};
        if (!adjacency.is_boundary())
        {
            diagonal[adjacency.owner] += coefficient;
            diagonal[adjacency.neighbor] += coefficient;
            owner_neighbor_coefficients[face_id] -= coefficient;
            neighbor_owner_coefficients[face_id] -= coefficient;
            continue;
        }

        if (boundary_conditions[face_boundary_ids[face_id]].type == ScalarBoundaryConditionType::Dirichlet)
        {
            diagonal[adjacency.owner] += coefficient;
        }
    }
}

void ScalarDiffusionOperator::add_boundary_rhs(const ScalarBoundaryConditions &boundary_conditions,
                                               const std::span<double> rhs) const
{
    validate_boundary_condition_count(*mesh_, boundary_conditions);
    validate_rhs_size(*mesh_, rhs);

    const auto face_adjacencies{mesh_->face_adjacencies()};
    const auto face_boundary_ids{mesh_->face_boundary_ids()};
    const auto face_lengths{mesh_->face_lengths()};
    for (Index face_id = 0; face_id < mesh_->face_count(); ++face_id)
    {
        const FaceAdjacency &adjacency{face_adjacencies[face_id]};
        if (!adjacency.is_boundary())
        {
            continue;
        }

        const ScalarBoundaryCondition &condition{boundary_conditions[face_boundary_ids[face_id]]};
        switch (condition.type)
        {
        case ScalarBoundaryConditionType::Dirichlet:
            rhs[adjacency.owner] += face_data_[face_id].primary_coefficient * condition.value;
            break;

        case ScalarBoundaryConditionType::Neumann:
            rhs[adjacency.owner] += diffusivity_ * condition.value * face_lengths[face_id];
            break;
        }
    }
}

void ScalarDiffusionOperator::add_non_orthogonal_rhs(const ScalarBoundaryConditions &boundary_conditions,
                                                     const CellVectorField &gradient, const std::span<double> rhs) const
{
    validate_boundary_condition_count(*mesh_, boundary_conditions);
    if (gradient.size() != mesh_->cell_count())
    {
        throw std::invalid_argument("Scalar diffusion gradient size must match the mesh cell count.");
    }
    validate_rhs_size(*mesh_, rhs);

    const auto face_adjacencies{mesh_->face_adjacencies()};
    const auto face_boundary_ids{mesh_->face_boundary_ids()};
    const auto gradients{gradient.values()};
    for (Index face_id = 0; face_id < mesh_->face_count(); ++face_id)
    {
        const FaceAdjacency &adjacency{face_adjacencies[face_id]};
        const FaceData &data{face_data_[face_id]};
        if (!adjacency.is_boundary())
        {
            const Vector2 &owner_gradient{gradients[adjacency.owner]};
            const Vector2 &neighbor_gradient{gradients[adjacency.neighbor]};
            const double lambda{data.neighbor_gradient_weight};
            const Vector2 face_gradient{
                owner_gradient.x + lambda * (neighbor_gradient.x - owner_gradient.x),
                owner_gradient.y + lambda * (neighbor_gradient.y - owner_gradient.y),
            };
            const double correction_flux{dot(face_gradient, data.correction_flux_vector)};
            rhs[adjacency.owner] -= correction_flux;
            rhs[adjacency.neighbor] += correction_flux;
            continue;
        }

        if (boundary_conditions[face_boundary_ids[face_id]].type == ScalarBoundaryConditionType::Dirichlet)
        {
            const double correction_flux{dot(gradients[adjacency.owner], data.correction_flux_vector)};
            rhs[adjacency.owner] -= correction_flux;
        }
    }
}

} // namespace cfd
