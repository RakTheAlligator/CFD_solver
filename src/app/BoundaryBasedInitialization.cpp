#include "app/BoundaryBasedInitialization.hpp"

#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVelocityField.hpp"
#include "cfd/field/PressureCorrectionBoundaryConditions.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/mesh/Mesh.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace cfd::app
{
namespace
{

void require_finite_condition(const ScalarBoundaryCondition &condition)
{
    switch (condition.type)
    {
    case ScalarBoundaryConditionType::Dirichlet:
    case ScalarBoundaryConditionType::Neumann:
        break;
    default:
        throw std::invalid_argument("Boundary initialization received an unsupported scalar condition type.");
    }
    if (!std::isfinite(condition.value))
    {
        throw std::invalid_argument("Boundary initialization requires finite boundary values.");
    }
}

} // namespace

void initialize_from_boundary_conditions(
    const Mesh &mesh, const ScalarBoundaryConditions &u_boundary_conditions,
    const ScalarBoundaryConditions &v_boundary_conditions, const ScalarBoundaryConditions &pressure_boundary_conditions,
    const PressureCorrectionBoundaryConditions &pressure_correction_boundary_conditions, Vector2 fallback_velocity,
    double fallback_pressure, CellVelocityField &velocity, CellScalarField &pressure)
{
    const Index boundary_count{mesh.boundary_groups().size()};
    if (u_boundary_conditions.size() != boundary_count || v_boundary_conditions.size() != boundary_count ||
        pressure_boundary_conditions.size() != boundary_count ||
        pressure_correction_boundary_conditions.size() != boundary_count || velocity.size() != mesh.cell_count() ||
        pressure.size() != mesh.cell_count())
    {
        throw std::invalid_argument("Boundary initialization cardinalities must match the Mesh.");
    }
    if (&pressure == &velocity.u() || &pressure == &velocity.v())
    {
        throw std::invalid_argument("Boundary initialization requires distinct velocity and pressure outputs.");
    }
    if (!std::isfinite(fallback_velocity.x) || !std::isfinite(fallback_velocity.y) || !std::isfinite(fallback_pressure))
    {
        throw std::invalid_argument("Boundary initialization requires finite internalField fallback values.");
    }
    for (BoundaryId boundary = 0; boundary < boundary_count; ++boundary)
    {
        require_finite_condition(u_boundary_conditions[boundary]);
        require_finite_condition(v_boundary_conditions[boundary]);
        require_finite_condition(pressure_boundary_conditions[boundary]);
        const auto pressure_type{pressure_boundary_conditions[boundary].type};
        switch (pressure_correction_boundary_conditions[boundary])
        {
        case PressureCorrectionBoundaryConditionType::FixedPressure:
            if (pressure_type != ScalarBoundaryConditionType::Dirichlet)
            {
                throw std::invalid_argument("Boundary initialization requires Dirichlet pressure for FixedPressure.");
            }
            break;
        case PressureCorrectionBoundaryConditionType::FixedMassFlux:
            if (pressure_type != ScalarBoundaryConditionType::Neumann)
            {
                throw std::invalid_argument("Boundary initialization requires Neumann pressure for FixedMassFlux.");
            }
            break;
        default:
            throw std::invalid_argument("Boundary initialization received an unsupported pressure-correction type.");
        }
    }

    // Wider accumulators avoid overflow/cancellation in the weighted sums;
    // the inflow threshold uses double epsilon, matching the stored geometry/data.
    long double inflow_area{};
    long double weighted_u{};
    long double weighted_v{};
    long double pressure_area{};
    long double weighted_pressure{};
    constexpr long double relative_flux_tolerance{64.0L * std::numeric_limits<double>::epsilon()};
    for (Index face = 0; face < mesh.face_count(); ++face)
    {
        if (!mesh.face_adjacencies()[face].is_boundary())
        {
            continue;
        }
        const double length{mesh.face_lengths()[face]};
        const Vector2 &area_vector{mesh.face_area_vectors()[face]};
        if (!std::isfinite(length) || !(length > 0.0) || !std::isfinite(area_vector.x) || !std::isfinite(area_vector.y))
        {
            throw std::runtime_error("Boundary initialization requires finite positive face geometry.");
        }
        const BoundaryId boundary{mesh.face_boundary_ids()[face]};
        const auto &u{u_boundary_conditions[boundary]};
        const auto &v{v_boundary_conditions[boundary]};
        if (u.type == ScalarBoundaryConditionType::Dirichlet && v.type == ScalarBoundaryConditionType::Dirichlet)
        {
            const long double q{static_cast<long double>(u.value) * area_vector.x +
                                static_cast<long double>(v.value) * area_vector.y};
            const long double scale{std::hypot(static_cast<long double>(u.value), static_cast<long double>(v.value)) *
                                    length};
            if (!std::isfinite(q) || !std::isfinite(scale))
            {
                throw std::runtime_error("Boundary initialization computed a non-finite boundary flux scale.");
            }
            if (q < -relative_flux_tolerance * scale)
            {
                inflow_area += length;
                weighted_u += static_cast<long double>(length) * u.value;
                weighted_v += static_cast<long double>(length) * v.value;
            }
        }
        const auto &p{pressure_boundary_conditions[boundary]};
        if (p.type == ScalarBoundaryConditionType::Dirichlet)
        {
            pressure_area += length;
            weighted_pressure += static_cast<long double>(length) * p.value;
        }
    }
    if (!std::isfinite(inflow_area) || !std::isfinite(pressure_area))
    {
        throw std::runtime_error("Boundary initialization accumulated a non-finite face area.");
    }
    Vector2 bulk{fallback_velocity};
    if (inflow_area > 0.0L)
    {
        bulk = {static_cast<double>(weighted_u / inflow_area), static_cast<double>(weighted_v / inflow_area)};
    }
    const double initial_pressure{pressure_area > 0.0L ? static_cast<double>(weighted_pressure / pressure_area)
                                                       : fallback_pressure};
    if (!std::isfinite(bulk.x) || !std::isfinite(bulk.y) || !std::isfinite(initial_pressure))
    {
        throw std::runtime_error("Boundary initialization computed a non-finite bulk state.");
    }
    for (Index cell = 0; cell < mesh.cell_count(); ++cell)
    {
        velocity.u()[cell] = bulk.x;
        velocity.v()[cell] = bulk.y;
        pressure[cell] = initial_pressure;
    }
}

} // namespace cfd::app
