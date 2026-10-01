#include "cfd/numerics/ScalarConvectionOperator.hpp"

#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVectorField.hpp"
#include "cfd/field/FaceFluxField.hpp"
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
#include <functional>
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

void validate_gradient_size(const Mesh &mesh, const CellVectorField &gradient)
{
    if (gradient.size() != mesh.cell_count())
    {
        throw std::invalid_argument("Scalar convection gradient size must match the mesh cell count.");
    }
}

[[nodiscard]]
bool spans_overlap(const std::span<const double> first, const std::span<const double> second) noexcept
{
    if (first.empty() || second.empty())
    {
        return false;
    }

    const double *const first_end{first.data() + first.size()};
    const double *const second_end{second.data() + second.size()};
    const std::less<const double *> pointer_less;
    return pointer_less(first.data(), second_end) && pointer_less(second.data(), first_end);
}

void validate_barth_jespersen_inputs(const Mesh &mesh, const CellScalarField &field, const CellVectorField &gradient,
                                     const FaceFluxField &face_flux, const std::span<double> limiter_workspace)
{
    if (field.size() != mesh.cell_count())
    {
        throw std::invalid_argument("Scalar convection field size must match the mesh cell count.");
    }
    if (limiter_workspace.size() != mesh.cell_count())
    {
        throw std::invalid_argument("Barth-Jespersen limiter workspace size must match the mesh cell count.");
    }
    for (const double value : field.values())
    {
        if (!std::isfinite(value))
        {
            throw std::invalid_argument("Barth-Jespersen scalar field values must be finite.");
        }
    }
    for (const Vector2 &value : gradient.values())
    {
        if (!std::isfinite(value.x) || !std::isfinite(value.y))
        {
            throw std::invalid_argument("Barth-Jespersen gradient values must be finite.");
        }
    }
    for (const double value : face_flux.values())
    {
        if (!std::isfinite(value))
        {
            throw std::invalid_argument("Barth-Jespersen face flux values must be finite.");
        }
    }
}

void validate_hybrid_conductance_count(const Mesh &mesh, const ScalarConvectionScheme scheme,
                                       const std::span<const double> face_diffusion_conductances)
{
    if (scheme == ScalarConvectionScheme::Hybrid && face_diffusion_conductances.size() != mesh.face_count())
    {
        throw std::invalid_argument("Hybrid convection conductance count must match the mesh face count.");
    }
}

[[nodiscard]]
double validated_hybrid_conductance(const std::span<const double> face_diffusion_conductances, const Index face_id)
{
    const double conductance{face_diffusion_conductances[face_id]};
    if (!std::isfinite(conductance) || !(conductance > 0.0))
    {
        throw std::invalid_argument("Hybrid convection conductances must be finite and strictly positive.");
    }
    return conductance;
}

[[nodiscard]]
bool hybrid_uses_linear_internal(const double carrier_flux, const double neighbor_weight,
                                 const double diffusion_conductance) noexcept
{
    if (carrier_flux >= 0.0)
    {
        return carrier_flux * neighbor_weight <= diffusion_conductance;
    }
    return -carrier_flux * (1.0 - neighbor_weight) <= diffusion_conductance;
}

[[nodiscard]]
bool hybrid_uses_linear_boundary(const double carrier_flux, const double diffusion_conductance) noexcept
{
    return carrier_flux <= diffusion_conductance;
}

[[nodiscard]]
double boundary_normal_distance(const Point2 &cell_center, const Point2 &face_center, const Vector2 &area_vector,
                                const double face_length) noexcept
{
    return ((face_center.x - cell_center.x) * area_vector.x + (face_center.y - cell_center.y) * area_vector.y) /
           face_length;
}

[[nodiscard]]
double linear_upwind_correction(const Point2 &upwind_center, const Point2 &face_center, const Vector2 &upwind_gradient,
                                const double limiter_coefficient = 1.0) noexcept
{
    return limiter_coefficient * (upwind_gradient.x * (face_center.x - upwind_center.x) +
                                  upwind_gradient.y * (face_center.y - upwind_center.y));
}
} // namespace

ScalarConvectionOperator::ScalarConvectionOperator(const Mesh &mesh) noexcept
    : mesh_(&mesh), scheme_(ScalarConvectionScheme::FirstOrderUpwind), limiter_(ScalarConvectionLimiter::None)
{
}

ScalarConvectionOperator::ScalarConvectionOperator(const Mesh &mesh, const ScalarConvectionScheme scheme,
                                                   const ScalarConvectionLimiter limiter)
    : mesh_(&mesh), scheme_(scheme), limiter_(limiter)
{
    switch (scheme_)
    {
    case ScalarConvectionScheme::FirstOrderUpwind:
    case ScalarConvectionScheme::LinearUpwind:
    case ScalarConvectionScheme::Linear:
    case ScalarConvectionScheme::Hybrid:
        break;

    default:
        throw std::invalid_argument("Scalar convection scheme is unsupported.");
    }

    switch (limiter_)
    {
    case ScalarConvectionLimiter::None:
        break;

    case ScalarConvectionLimiter::BarthJespersen:
        if (scheme_ != ScalarConvectionScheme::LinearUpwind)
        {
            throw std::invalid_argument("Barth-Jespersen limiting is supported only with LinearUpwind convection.");
        }
        break;

    default:
        throw std::invalid_argument("Scalar convection limiter is unsupported.");
    }

    if (scheme_ == ScalarConvectionScheme::FirstOrderUpwind || scheme_ == ScalarConvectionScheme::LinearUpwind)
    {
        return;
    }

    internal_face_interpolation_weights_.resize(mesh_->face_count());
    const auto face_adjacencies{mesh_->face_adjacencies()};
    const auto cell_centers{mesh_->cell_centers()};
    const auto face_centers{mesh_->face_centers()};
    const auto face_lengths{mesh_->face_lengths()};
    const auto face_area_vectors{mesh_->face_area_vectors()};
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
        const double face_length{face_lengths[face_id]};
        const InternalFaceInterpolationGeometry geometry{compute_internal_face_interpolation_geometry(
            face_id, owner_center, neighbor_center, face_center, area_vector, face_length)};
        internal_face_interpolation_weights_[face_id] = geometry.neighbor_weight;
    }
}

HybridConvectionFaceCounts ScalarConvectionOperator::classify_hybrid_faces(
    const FaceFluxField &face_flux, const std::span<const double> face_diffusion_conductances) const
{
    if (scheme_ != ScalarConvectionScheme::Hybrid)
    {
        throw std::invalid_argument("Hybrid face classification requires a Hybrid convection operator.");
    }
    validate_face_flux_count(*mesh_, face_flux);
    validate_hybrid_conductance_count(*mesh_, scheme_, face_diffusion_conductances);

    HybridConvectionFaceCounts counts;
    const auto face_adjacencies{mesh_->face_adjacencies()};
    const auto flux_values{face_flux.values()};
    for (Index face_id = 0; face_id < mesh_->face_count(); ++face_id)
    {
        const FaceAdjacency &adjacency{face_adjacencies[face_id]};
        const double carrier_flux{flux_values[face_id]};
        const double diffusion_conductance{validated_hybrid_conductance(face_diffusion_conductances, face_id)};
        if (adjacency.is_boundary())
        {
            if (hybrid_uses_linear_boundary(carrier_flux, diffusion_conductance))
            {
                ++counts.boundary_linear_faces;
            }
            else
            {
                ++counts.boundary_upwind_faces;
            }
            continue;
        }

        const double neighbor_weight{internal_face_interpolation_weights_[face_id]};
        if (hybrid_uses_linear_internal(carrier_flux, neighbor_weight, diffusion_conductance))
        {
            ++counts.internal_linear_faces;
        }
        else
        {
            ++counts.internal_upwind_faces;
        }
    }
    return counts;
}

void ScalarConvectionOperator::compute_flux_balance(const CellScalarField &field,
                                                    const ScalarBoundaryConditions &boundary_conditions,
                                                    const FaceFluxField &face_flux, CellScalarField &flux_balance,
                                                    const std::span<const double> face_diffusion_conductances) const
{
    const Index cell_count{mesh_->cell_count()};

    if (scheme_ == ScalarConvectionScheme::LinearUpwind)
    {
        throw std::invalid_argument("LinearUpwind convection flux balance requires a cell gradient.");
    }

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
    validate_hybrid_conductance_count(*mesh_, scheme_, face_diffusion_conductances);
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

    if (scheme_ == ScalarConvectionScheme::Linear)
    {
        for (Index face_id = 0; face_id < mesh_->face_count(); ++face_id)
        {
            const FaceAdjacency &adjacency{face_adjacencies[face_id]};
            const Index owner_id{adjacency.owner};
            const double carrier_flux{flux_values[face_id]};

            if (!adjacency.is_boundary())
            {
                const Index neighbor_id{adjacency.neighbor};
                const double neighbor_weight{internal_face_interpolation_weights_[face_id]};
                const double face_value{(1.0 - neighbor_weight) * field_values[owner_id] +
                                        neighbor_weight * field_values[neighbor_id]};
                const double convective_flux{carrier_flux * face_value};
                balance_values[owner_id] += convective_flux;
                balance_values[neighbor_id] -= convective_flux;
                continue;
            }

            const ScalarBoundaryCondition &condition{boundary_conditions[face_boundary_ids[face_id]]};
            double face_value{field_values[owner_id]};
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
            balance_values[owner_id] += carrier_flux * face_value;
        }
        return;
    }

    if (scheme_ == ScalarConvectionScheme::Hybrid)
    {
        for (Index face_id = 0; face_id < mesh_->face_count(); ++face_id)
        {
            const FaceAdjacency &adjacency{face_adjacencies[face_id]};
            const Index owner_id{adjacency.owner};
            const double carrier_flux{flux_values[face_id]};
            const double diffusion_conductance{validated_hybrid_conductance(face_diffusion_conductances, face_id)};

            if (!adjacency.is_boundary())
            {
                const Index neighbor_id{adjacency.neighbor};
                const double neighbor_weight{internal_face_interpolation_weights_[face_id]};
                double face_value{};
                if (hybrid_uses_linear_internal(carrier_flux, neighbor_weight, diffusion_conductance))
                {
                    face_value =
                        (1.0 - neighbor_weight) * field_values[owner_id] + neighbor_weight * field_values[neighbor_id];
                }
                else
                {
                    face_value = carrier_flux >= 0.0 ? field_values[owner_id] : field_values[neighbor_id];
                }
                const double convective_flux{carrier_flux * face_value};
                balance_values[owner_id] += convective_flux;
                balance_values[neighbor_id] -= convective_flux;
                continue;
            }

            double face_value{field_values[owner_id]};
            if (hybrid_uses_linear_boundary(carrier_flux, diffusion_conductance))
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
        return;
    }

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

void ScalarConvectionOperator::compute_barth_jespersen_limiter(const CellScalarField &field,
                                                               const ScalarBoundaryConditions &boundary_conditions,
                                                               const CellVectorField &gradient,
                                                               const std::span<double> limiter_coefficients) const
{
    const auto cell_offsets{mesh_->cell_node_offsets()};
    const auto cell_faces{mesh_->cell_faces()};
    const auto face_adjacencies{mesh_->face_adjacencies()};
    const auto face_boundary_ids{mesh_->face_boundary_ids()};
    const auto cell_centers{mesh_->cell_centers()};
    const auto face_centers{mesh_->face_centers()};
    const auto face_lengths{mesh_->face_lengths()};
    const auto face_area_vectors{mesh_->face_area_vectors()};
    const auto field_values{field.values()};
    const auto gradient_values{gradient.values()};

    for (Index cell_id = 0; cell_id < mesh_->cell_count(); ++cell_id)
    {
        const double cell_value{field_values[cell_id]};
        double minimum_value{cell_value};
        double maximum_value{cell_value};
        const Index face_begin{cell_offsets[cell_id]};
        const Index face_end{cell_offsets[cell_id + 1]};

        for (Index position = face_begin; position < face_end; ++position)
        {
            const Index face_id{cell_faces[position]};
            const FaceAdjacency &adjacency{face_adjacencies[face_id]};
            double stencil_value{};
            if (!adjacency.is_boundary())
            {
                const Index other_cell_id{adjacency.owner == cell_id ? adjacency.neighbor : adjacency.owner};
                stencil_value = field_values[other_cell_id];
            }
            else
            {
                const ScalarBoundaryCondition &condition{boundary_conditions[face_boundary_ids[face_id]]};
                switch (condition.type)
                {
                case ScalarBoundaryConditionType::Dirichlet:
                    stencil_value = condition.value;
                    break;

                case ScalarBoundaryConditionType::Neumann:
                    stencil_value = cell_value + condition.value * boundary_normal_distance(cell_centers[cell_id],
                                                                                            face_centers[face_id],
                                                                                            face_area_vectors[face_id],
                                                                                            face_lengths[face_id]);
                    break;
                }
            }
            minimum_value = std::min(minimum_value, stencil_value);
            maximum_value = std::max(maximum_value, stencil_value);
        }

        double coefficient{1.0};
        for (Index position = face_begin; position < face_end; ++position)
        {
            const Index face_id{cell_faces[position]};
            const double correction{
                linear_upwind_correction(cell_centers[cell_id], face_centers[face_id], gradient_values[cell_id])};
            if (correction > 0.0)
            {
                coefficient = std::min(coefficient, (maximum_value - cell_value) / correction);
            }
            else if (correction < 0.0)
            {
                coefficient = std::min(coefficient, (minimum_value - cell_value) / correction);
            }
        }
        limiter_coefficients[cell_id] = std::clamp(coefficient, 0.0, 1.0);
    }
}

void ScalarConvectionOperator::compute_flux_balance(const CellScalarField &field,
                                                    const ScalarBoundaryConditions &boundary_conditions,
                                                    const FaceFluxField &face_flux, const CellVectorField &gradient,
                                                    CellScalarField &flux_balance) const
{
    if (scheme_ != ScalarConvectionScheme::LinearUpwind)
    {
        throw std::invalid_argument("Gradient-aware convection flux balance requires a LinearUpwind operator.");
    }

    if (limiter_ != ScalarConvectionLimiter::None)
    {
        throw std::invalid_argument("Limited LinearUpwind flux balance requires a limiter workspace.");
    }

    const Index cell_count{mesh_->cell_count()};
    if (field.size() != cell_count)
    {
        throw std::invalid_argument("Scalar convection field size must match the mesh cell count.");
    }
    if (flux_balance.size() != cell_count)
    {
        throw std::invalid_argument("Scalar convection output size must match the mesh cell count.");
    }
    validate_gradient_size(*mesh_, gradient);
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
    const auto gradient_values{gradient.values()};
    const auto flux_values{face_flux.values()};
    auto balance_values{flux_balance.values()};

    std::fill(balance_values.begin(), balance_values.end(), 0.0);

    for (Index face_id = 0; face_id < mesh_->face_count(); ++face_id)
    {
        const FaceAdjacency &adjacency{face_adjacencies[face_id]};
        const Index owner_id{adjacency.owner};
        const double carrier_flux{flux_values[face_id]};
        if (carrier_flux == 0.0)
        {
            continue;
        }

        if (!adjacency.is_boundary())
        {
            const Index neighbor_id{adjacency.neighbor};
            const Index upwind_id{carrier_flux >= 0.0 ? owner_id : neighbor_id};
            const double face_value{field_values[upwind_id] + linear_upwind_correction(cell_centers[upwind_id],
                                                                                       face_centers[face_id],
                                                                                       gradient_values[upwind_id])};
            const double convective_flux{carrier_flux * face_value};
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
        else
        {
            face_value +=
                linear_upwind_correction(cell_centers[owner_id], face_centers[face_id], gradient_values[owner_id]);
        }
        balance_values[owner_id] += carrier_flux * face_value;
    }
}

void ScalarConvectionOperator::compute_flux_balance(const CellScalarField &field,
                                                    const ScalarBoundaryConditions &boundary_conditions,
                                                    const FaceFluxField &face_flux, const CellVectorField &gradient,
                                                    const std::span<double> limiter_workspace,
                                                    CellScalarField &flux_balance) const
{
    if (scheme_ != ScalarConvectionScheme::LinearUpwind)
    {
        throw std::invalid_argument("Gradient-aware convection flux balance requires a LinearUpwind operator.");
    }

    const Index cell_count{mesh_->cell_count()};
    if (field.size() != cell_count)
    {
        throw std::invalid_argument("Scalar convection field size must match the mesh cell count.");
    }
    if (flux_balance.size() != cell_count)
    {
        throw std::invalid_argument("Scalar convection output size must match the mesh cell count.");
    }
    validate_gradient_size(*mesh_, gradient);
    validate_face_flux_count(*mesh_, face_flux);
    validate_boundary_condition_count(*mesh_, boundary_conditions);
    if (&field == &flux_balance)
    {
        throw std::invalid_argument("Scalar convection output must not alias the input scalar field.");
    }
    if (spans_overlap(limiter_workspace, field.values()) || spans_overlap(limiter_workspace, face_flux.values()) ||
        spans_overlap(limiter_workspace, flux_balance.values()))
    {
        throw std::invalid_argument("Barth-Jespersen limiter workspace must not overlap input or output storage.");
    }

    if (limiter_ == ScalarConvectionLimiter::BarthJespersen)
    {
        validate_barth_jespersen_inputs(*mesh_, field, gradient, face_flux, limiter_workspace);
        compute_barth_jespersen_limiter(field, boundary_conditions, gradient, limiter_workspace);
    }

    const auto face_adjacencies{mesh_->face_adjacencies()};
    const auto face_boundary_ids{mesh_->face_boundary_ids()};
    const auto cell_centers{mesh_->cell_centers()};
    const auto face_centers{mesh_->face_centers()};
    const auto face_lengths{mesh_->face_lengths()};
    const auto face_area_vectors{mesh_->face_area_vectors()};
    const auto field_values{field.values()};
    const auto gradient_values{gradient.values()};
    const auto flux_values{face_flux.values()};
    auto balance_values{flux_balance.values()};

    std::fill(balance_values.begin(), balance_values.end(), 0.0);

    for (Index face_id = 0; face_id < mesh_->face_count(); ++face_id)
    {
        const FaceAdjacency &adjacency{face_adjacencies[face_id]};
        const Index owner_id{adjacency.owner};
        const double carrier_flux{flux_values[face_id]};
        if (carrier_flux == 0.0)
        {
            continue;
        }

        if (!adjacency.is_boundary())
        {
            const Index neighbor_id{adjacency.neighbor};
            const Index upwind_id{carrier_flux >= 0.0 ? owner_id : neighbor_id};
            const double limiter_coefficient{
                limiter_ == ScalarConvectionLimiter::BarthJespersen ? limiter_workspace[upwind_id] : 1.0};
            const double face_value{field_values[upwind_id] +
                                    linear_upwind_correction(cell_centers[upwind_id], face_centers[face_id],
                                                             gradient_values[upwind_id], limiter_coefficient)};
            const double convective_flux{carrier_flux * face_value};
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
        else
        {
            const double limiter_coefficient{
                limiter_ == ScalarConvectionLimiter::BarthJespersen ? limiter_workspace[owner_id] : 1.0};
            face_value += linear_upwind_correction(cell_centers[owner_id], face_centers[face_id],
                                                   gradient_values[owner_id], limiter_coefficient);
        }
        balance_values[owner_id] += carrier_flux * face_value;
    }
}

void ScalarConvectionOperator::add_matrix_contributions(const ScalarBoundaryConditions &boundary_conditions,
                                                        const FaceFluxField &face_flux, ScalarLinearSystem &system,
                                                        const std::span<const double> face_diffusion_conductances) const
{
    validate_boundary_condition_count(*mesh_, boundary_conditions);
    validate_face_flux_count(*mesh_, face_flux);
    validate_hybrid_conductance_count(*mesh_, scheme_, face_diffusion_conductances);
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

    if (scheme_ == ScalarConvectionScheme::Linear)
    {
        for (Index face_id = 0; face_id < mesh_->face_count(); ++face_id)
        {
            const FaceAdjacency &adjacency{face_adjacencies[face_id]};
            const double carrier_flux{flux_values[face_id]};

            if (!adjacency.is_boundary())
            {
                const double neighbor_weight{internal_face_interpolation_weights_[face_id]};
                const double owner_weight{1.0 - neighbor_weight};
                diagonal[adjacency.owner] += carrier_flux * owner_weight;
                owner_neighbor_coefficients[face_id] += carrier_flux * neighbor_weight;
                neighbor_owner_coefficients[face_id] -= carrier_flux * owner_weight;
                diagonal[adjacency.neighbor] -= carrier_flux * neighbor_weight;
                continue;
            }

            const ScalarBoundaryCondition &condition{boundary_conditions[face_boundary_ids[face_id]]};
            if (condition.type == ScalarBoundaryConditionType::Neumann)
            {
                diagonal[adjacency.owner] += carrier_flux;
            }
        }
        return;
    }

    if (scheme_ == ScalarConvectionScheme::Hybrid)
    {
        for (Index face_id = 0; face_id < mesh_->face_count(); ++face_id)
        {
            const FaceAdjacency &adjacency{face_adjacencies[face_id]};
            const double carrier_flux{flux_values[face_id]};
            const double diffusion_conductance{validated_hybrid_conductance(face_diffusion_conductances, face_id)};

            if (!adjacency.is_boundary())
            {
                const double neighbor_weight{internal_face_interpolation_weights_[face_id]};
                if (hybrid_uses_linear_internal(carrier_flux, neighbor_weight, diffusion_conductance))
                {
                    const double owner_weight{1.0 - neighbor_weight};
                    diagonal[adjacency.owner] += carrier_flux * owner_weight;
                    owner_neighbor_coefficients[face_id] += carrier_flux * neighbor_weight;
                    neighbor_owner_coefficients[face_id] -= carrier_flux * owner_weight;
                    diagonal[adjacency.neighbor] -= carrier_flux * neighbor_weight;
                }
                else
                {
                    const double positive_flux{std::max(carrier_flux, 0.0)};
                    const double negative_flux{std::min(carrier_flux, 0.0)};
                    diagonal[adjacency.owner] += positive_flux;
                    owner_neighbor_coefficients[face_id] += negative_flux;
                    diagonal[adjacency.neighbor] -= negative_flux;
                    neighbor_owner_coefficients[face_id] -= positive_flux;
                }
                continue;
            }

            const ScalarBoundaryCondition &condition{boundary_conditions[face_boundary_ids[face_id]]};
            if (hybrid_uses_linear_boundary(carrier_flux, diffusion_conductance))
            {
                if (condition.type == ScalarBoundaryConditionType::Neumann)
                {
                    diagonal[adjacency.owner] += carrier_flux;
                }
            }
            else
            {
                diagonal[adjacency.owner] += carrier_flux;
            }
        }
        return;
    }

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
                                                const FaceFluxField &face_flux, const std::span<double> rhs,
                                                const std::span<const double> face_diffusion_conductances) const
{
    validate_boundary_condition_count(*mesh_, boundary_conditions);
    validate_face_flux_count(*mesh_, face_flux);
    validate_rhs_size(*mesh_, rhs);
    validate_hybrid_conductance_count(*mesh_, scheme_, face_diffusion_conductances);

    const auto face_adjacencies{mesh_->face_adjacencies()};
    const auto face_boundary_ids{mesh_->face_boundary_ids()};
    const auto cell_centers{mesh_->cell_centers()};
    const auto face_centers{mesh_->face_centers()};
    const auto face_lengths{mesh_->face_lengths()};
    const auto face_area_vectors{mesh_->face_area_vectors()};
    const auto flux_values{face_flux.values()};

    if (scheme_ == ScalarConvectionScheme::Linear)
    {
        for (Index face_id = 0; face_id < mesh_->face_count(); ++face_id)
        {
            const FaceAdjacency &adjacency{face_adjacencies[face_id]};
            if (!adjacency.is_boundary())
            {
                continue;
            }

            const double carrier_flux{flux_values[face_id]};
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
        return;
    }

    if (scheme_ == ScalarConvectionScheme::Hybrid)
    {
        for (Index face_id = 0; face_id < mesh_->face_count(); ++face_id)
        {
            const FaceAdjacency &adjacency{face_adjacencies[face_id]};
            const double diffusion_conductance{validated_hybrid_conductance(face_diffusion_conductances, face_id)};
            if (!adjacency.is_boundary())
            {
                continue;
            }

            const double carrier_flux{flux_values[face_id]};
            if (!hybrid_uses_linear_boundary(carrier_flux, diffusion_conductance))
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
        return;
    }

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

void ScalarConvectionOperator::add_deferred_correction_rhs(const CellVectorField &gradient,
                                                           const FaceFluxField &face_flux,
                                                           const std::span<double> rhs) const
{
    if (scheme_ != ScalarConvectionScheme::LinearUpwind)
    {
        return;
    }
    if (limiter_ != ScalarConvectionLimiter::None)
    {
        throw std::invalid_argument("Limited LinearUpwind assembly requires field data and a limiter workspace.");
    }

    validate_gradient_size(*mesh_, gradient);
    validate_face_flux_count(*mesh_, face_flux);
    validate_rhs_size(*mesh_, rhs);
    if (spans_overlap(rhs, face_flux.values()))
    {
        throw std::invalid_argument("Scalar convection RHS must not overlap face-flux storage.");
    }
    add_linear_upwind_deferred_correction(gradient, face_flux, {}, rhs);
}

void ScalarConvectionOperator::add_linear_upwind_deferred_correction(const CellVectorField &gradient,
                                                                     const FaceFluxField &face_flux,
                                                                     const std::span<const double> limiter_coefficients,
                                                                     const std::span<double> rhs) const
{
    const auto face_adjacencies{mesh_->face_adjacencies()};
    const auto cell_centers{mesh_->cell_centers()};
    const auto face_centers{mesh_->face_centers()};
    const auto gradient_values{gradient.values()};
    const auto flux_values{face_flux.values()};
    for (Index face_id = 0; face_id < mesh_->face_count(); ++face_id)
    {
        const FaceAdjacency &adjacency{face_adjacencies[face_id]};
        const double carrier_flux{flux_values[face_id]};
        if (carrier_flux == 0.0)
        {
            continue;
        }
        if (!adjacency.is_boundary())
        {
            const Index upwind_id{carrier_flux >= 0.0 ? adjacency.owner : adjacency.neighbor};
            const double limiter_coefficient{limiter_coefficients.empty() ? 1.0 : limiter_coefficients[upwind_id]};
            const double correction{linear_upwind_correction(cell_centers[upwind_id], face_centers[face_id],
                                                             gradient_values[upwind_id], limiter_coefficient)};
            const double correction_flux{carrier_flux * correction};
            rhs[adjacency.owner] -= correction_flux;
            rhs[adjacency.neighbor] += correction_flux;
            continue;
        }
        if (carrier_flux > 0.0)
        {
            const double limiter_coefficient{limiter_coefficients.empty() ? 1.0
                                                                          : limiter_coefficients[adjacency.owner]};
            const double correction{linear_upwind_correction(cell_centers[adjacency.owner], face_centers[face_id],
                                                             gradient_values[adjacency.owner], limiter_coefficient)};
            rhs[adjacency.owner] -= carrier_flux * correction;
        }
    }
}

void ScalarConvectionOperator::add_deferred_correction_rhs(
    const CellScalarField &field, const ScalarBoundaryConditions &boundary_conditions, const CellVectorField &gradient,
    const FaceFluxField &face_flux, const std::span<double> limiter_workspace, const std::span<double> rhs) const
{
    if (scheme_ != ScalarConvectionScheme::LinearUpwind)
    {
        return;
    }

    validate_gradient_size(*mesh_, gradient);
    validate_face_flux_count(*mesh_, face_flux);
    validate_rhs_size(*mesh_, rhs);
    if (spans_overlap(rhs, field.values()) || spans_overlap(rhs, face_flux.values()))
    {
        throw std::invalid_argument("Scalar convection RHS must not overlap scalar-field or face-flux storage.");
    }
    if (spans_overlap(limiter_workspace, field.values()) || spans_overlap(limiter_workspace, face_flux.values()) ||
        spans_overlap(limiter_workspace, rhs))
    {
        throw std::invalid_argument("Barth-Jespersen limiter workspace must not overlap input or RHS storage.");
    }

    if (limiter_ == ScalarConvectionLimiter::None)
    {
        add_linear_upwind_deferred_correction(gradient, face_flux, {}, rhs);
        return;
    }

    validate_boundary_condition_count(*mesh_, boundary_conditions);
    validate_barth_jespersen_inputs(*mesh_, field, gradient, face_flux, limiter_workspace);
    compute_barth_jespersen_limiter(field, boundary_conditions, gradient, limiter_workspace);
    add_linear_upwind_deferred_correction(gradient, face_flux, limiter_workspace, rhs);
}
} // namespace cfd
