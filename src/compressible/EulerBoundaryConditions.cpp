#include "cfd/compressible/EulerBoundaryConditions.hpp"

#include <algorithm>
#include <cmath>
#include <string>

namespace cfd
{
namespace
{
void require_positive(const double value, const char *name)
{
    if (!std::isfinite(value) || !(value > 0.0))
        throw std::invalid_argument(std::string{"Euler BC: "} + name + " must be finite and strictly positive.");
}

Vector2 normalized_direction(const Vector2 direction)
{
    const double scale{std::max(std::abs(direction.x), std::abs(direction.y))};
    if (!std::isfinite(direction.x) || !std::isfinite(direction.y) || !(scale > 0.0))
        throw std::invalid_argument("Euler BC requires a finite nonzero direction.");
    const Vector2 scaled{direction.x / scale, direction.y / scale};
    const double norm{std::hypot(scaled.x, scaled.y)};
    return {scaled.x / norm, scaled.y / norm};
}
} // namespace

EulerBoundaryCondition::EulerBoundaryCondition(const EulerBoundaryConditionType type) : type_(type)
{
    switch (type)
    {
    case EulerBoundaryConditionType::Axis:
    case EulerBoundaryConditionType::SlipWall:
    case EulerBoundaryConditionType::Transmissive:
    case EulerBoundaryConditionType::SupersonicOutlet:
        return;
    case EulerBoundaryConditionType::TotalConditionsInlet:
    case EulerBoundaryConditionType::MassFlowInlet:
    case EulerBoundaryConditionType::SupersonicInlet:
    case EulerBoundaryConditionType::PressureOutlet:
    case EulerBoundaryConditionType::Farfield:
        throw std::invalid_argument("Euler boundary condition type requires physical data from its factory.");
    }
    throw std::invalid_argument("Unsupported Euler boundary condition type.");
}

EulerBoundaryCondition::EulerBoundaryCondition(const EulerBoundaryConditionType type, const double pressure,
                                               const double temperature, const double mass_flow,
                                               const Vector2 direction, const PrimitiveState state)
    : type_(type), pressure_(pressure), temperature_(temperature), mass_flow_(mass_flow), direction_(direction),
      state_(state)
{
}

EulerBoundaryCondition EulerBoundaryCondition::total_conditions_inlet(const double total_pressure,
                                                                      const double total_temperature,
                                                                      const Vector2 direction)
{
    require_positive(total_pressure, "total pressure");
    require_positive(total_temperature, "total temperature");
    return {EulerBoundaryConditionType::TotalConditionsInlet,
            total_pressure,
            total_temperature,
            0.0,
            normalized_direction(direction),
            {}};
}

EulerBoundaryCondition EulerBoundaryCondition::mass_flow_inlet(const double inward_mass_flow,
                                                               const double total_temperature, const Vector2 direction)
{
    require_positive(inward_mass_flow, "inward group mass flow");
    require_positive(total_temperature, "total temperature");
    return {EulerBoundaryConditionType::MassFlowInlet, 0.0, total_temperature, inward_mass_flow,
            normalized_direction(direction),           {}};
}

EulerBoundaryCondition EulerBoundaryCondition::supersonic_inlet(const PrimitiveState state)
{
    require_positive(state.density, "density");
    require_positive(state.pressure, "static pressure");
    if (!std::isfinite(state.velocity_x) || !std::isfinite(state.velocity_r))
        throw std::invalid_argument("Euler BC velocities must be finite.");
    return {EulerBoundaryConditionType::SupersonicInlet, 0.0, 0.0, 0.0, {}, state};
}

EulerBoundaryCondition EulerBoundaryCondition::pressure_outlet(const double pressure)
{
    require_positive(pressure, "static pressure");
    return {EulerBoundaryConditionType::PressureOutlet, pressure, 0.0, 0.0, {}, {}};
}

EulerBoundaryCondition EulerBoundaryCondition::farfield(const PrimitiveState exterior)
{
    require_positive(exterior.density, "farfield density");
    require_positive(exterior.pressure, "farfield static pressure");
    if (!std::isfinite(exterior.velocity_x) || !std::isfinite(exterior.velocity_r))
        throw std::invalid_argument("Euler BC velocities must be finite.");
    return {EulerBoundaryConditionType::Farfield, 0.0, 0.0, 0.0, {}, exterior};
}
} // namespace cfd
