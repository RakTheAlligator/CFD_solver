#include "cfd/numerics/ScalarConvectionOperator.hpp"

#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/FaceFluxField.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/linear_algebra/ScalarLinearSystem.hpp"
#include "cfd/math/Point2.hpp"
#include "cfd/math/Vector2.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/Types.hpp"

#include <algorithm>
#include <stdexcept>

namespace cfd
{
namespace
{

void validate_boundary_condition_count(const Mesh &mesh, const ScalarBoundaryConditions &boundary_conditions)
{
    if (boundary_conditions.size() != mesh.boundary_groups().size())
    {
        throw std::invalid_argument("Scalar convection boundary-condition count must match the mesh boundary count.");
    }
}

void validate_face_flux_count(const Mesh &mesh, const FaceFluxField &face_flux)
{
    if (face_flux.size() != mesh.face_count())
    {
        throw std::invalid_argument("Scalar convection face-flux size must match the mesh face count.");
    }
}

void validate_rhs_size(const Mesh &mesh, const std::span<double> rhs)
{
    if (rhs.size() != mesh.cell_count())
    {
        throw std::invalid_argument("Scalar convection RHS size must match the mesh cell count.");
    }
}

[[nodiscard]]
double boundary_normal_distance(const Point2 &cell_center, const Point2 &face_center, const Vector2 &area_vector,
                                const double face_length) noexcept
{
    return ((face_center.x - cell_center.x) * area_vector.x + (face_center.y - cell_center.y) * area_vector.y) /
           face_length;
}

} // namespace

ScalarConvectionOperator::ScalarConvectionOperator(const Mesh &mesh) noexcept : mesh_(&mesh)
{
}

void ScalarConvectionOperator::compute_flux_balance(const CellScalarField &field,
                                                    const ScalarBoundaryConditions &boundary_conditions,
                                                    const FaceFluxField &face_flux, CellScalarField &flux_balance) const
{
    const Index cell_count{mesh_->cell_count()};

    if (field.size() != cell_count)
    {
        throw std::invalid_argument("Scalar convection field size must match the mesh cell count.");
    }
    if (flux_balance.size() != cell_count)
    {
        throw std::invalid_argument("Scalar convection output size must match the mesh cell count.");
    }
    validate_face_flux_count(*mesh_, face_flux);
    validate_boundary_condition_count(*mesh_, boundary_conditions);
    if (&field == &flux_balance)
    {
        throw std::invalid_argument("Scalar convection output must not alias the input scalar field.");
    }

    const auto face_adjacencies{mesh_->face_adjacencies()};
    const auto face_boundary_ids{mesh_->face_boundary_ids()};
    const auto cell_centers{mesh_->cell_centers()};
    const auto face_centers{mesh_->face_centers()};
    const auto face_lengths{mesh_->face_lengths()};
    const auto face_area_vectors{mesh_->face_area_vectors()};
    const auto field_values{field.values()};
    const auto flux_values{face_flux.values()};
    auto balance_values{flux_balance.values()};

    std::fill(balance_values.begin(), balance_values.end(), 0.0);

    for (Index face_id = 0; face_id < mesh_->face_count(); ++face_id)
    {
        const FaceAdjacency &adjacency{face_adjacencies[face_id]};
        const Index owner_id{adjacency.owner};
        const double carrier_flux{flux_values[face_id]};

        if (!adjacency.is_boundary())
        {
            const Index neighbor_id{adjacency.neighbor};
            const double upwind_value{carrier_flux >= 0.0 ? field_values[owner_id] : field_values[neighbor_id]};
            const double convective_flux{carrier_flux * upwind_value};
            balance_values[owner_id] += convective_flux;
            balance_values[neighbor_id] -= convective_flux;
            continue;
        }

        double face_value{field_values[owner_id]};
        if (carrier_flux < 0.0)
        {
            const ScalarBoundaryCondition &condition{boundary_conditions[face_boundary_ids[face_id]]};
            switch (condition.type)
            {
            case ScalarBoundaryConditionType::Dirichlet:
                face_value = condition.value;
                break;

            case ScalarBoundaryConditionType::Neumann:
                face_value +=
                    condition.value * boundary_normal_distance(cell_centers[owner_id], face_centers[face_id],
                                                               face_area_vectors[face_id], face_lengths[face_id]);
                break;
            }
        }

        balance_values[owner_id] += carrier_flux * face_value;
    }
}

void ScalarConvectionOperator::add_matrix_contributions(const ScalarBoundaryConditions &boundary_conditions,
                                                        const FaceFluxField &face_flux,
                                                        ScalarLinearSystem &system) const
{
    validate_boundary_condition_count(*mesh_, boundary_conditions);
    validate_face_flux_count(*mesh_, face_flux);
    if (&system.mesh() != mesh_ || system.cell_count() != mesh_->cell_count() ||
        system.face_count() != mesh_->face_count())
    {
        throw std::invalid_argument("Scalar convection system must reference the operator Mesh.");
    }

    auto diagonal{system.diagonal()};
    auto owner_neighbor_coefficients{system.owner_neighbor_coefficients()};
    auto neighbor_owner_coefficients{system.neighbor_owner_coefficients()};
    const auto face_adjacencies{mesh_->face_adjacencies()};
    const auto face_boundary_ids{mesh_->face_boundary_ids()};
    const auto flux_values{face_flux.values()};

    for (Index face_id = 0; face_id < mesh_->face_count(); ++face_id)
    {
        const FaceAdjacency &adjacency{face_adjacencies[face_id]};
        const double carrier_flux{flux_values[face_id]};

        if (!adjacency.is_boundary())
        {
            const double positive_flux{std::max(carrier_flux, 0.0)};
            const double negative_flux{std::min(carrier_flux, 0.0)};
            diagonal[adjacency.owner] += positive_flux;
            owner_neighbor_coefficients[face_id] += negative_flux;
            diagonal[adjacency.neighbor] -= negative_flux;
            neighbor_owner_coefficients[face_id] -= positive_flux;
            continue;
        }

        const ScalarBoundaryCondition &condition{boundary_conditions[face_boundary_ids[face_id]]};
        if (carrier_flux >= 0.0 || condition.type == ScalarBoundaryConditionType::Neumann)
        {
            diagonal[adjacency.owner] += carrier_flux;
        }
    }
}

void ScalarConvectionOperator::add_boundary_rhs(const ScalarBoundaryConditions &boundary_conditions,
                                                const FaceFluxField &face_flux, const std::span<double> rhs) const
{
    validate_boundary_condition_count(*mesh_, boundary_conditions);
    validate_face_flux_count(*mesh_, face_flux);
    validate_rhs_size(*mesh_, rhs);

    const auto face_adjacencies{mesh_->face_adjacencies()};
    const auto face_boundary_ids{mesh_->face_boundary_ids()};
    const auto cell_centers{mesh_->cell_centers()};
    const auto face_centers{mesh_->face_centers()};
    const auto face_lengths{mesh_->face_lengths()};
    const auto face_area_vectors{mesh_->face_area_vectors()};
    const auto flux_values{face_flux.values()};

    for (Index face_id = 0; face_id < mesh_->face_count(); ++face_id)
    {
        const FaceAdjacency &adjacency{face_adjacencies[face_id]};
        const double carrier_flux{flux_values[face_id]};
        if (!adjacency.is_boundary() || carrier_flux >= 0.0)
        {
            continue;
        }

        const ScalarBoundaryCondition &condition{boundary_conditions[face_boundary_ids[face_id]]};
        switch (condition.type)
        {
        case ScalarBoundaryConditionType::Dirichlet:
            rhs[adjacency.owner] -= carrier_flux * condition.value;
            break;

        case ScalarBoundaryConditionType::Neumann:
            rhs[adjacency.owner] -= carrier_flux * condition.value *
                                    boundary_normal_distance(cell_centers[adjacency.owner], face_centers[face_id],
                                                             face_area_vectors[face_id], face_lengths[face_id]);
            break;
        }
    }
}

} // namespace cfd
