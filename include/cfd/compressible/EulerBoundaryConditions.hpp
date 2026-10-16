#pragma once

#include "cfd/compressible/CaloricallyPerfectGas.hpp"
#include "cfd/math/Vector2.hpp"
#include "cfd/mesh/Boundary.hpp"

#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cfd
{

enum class EulerBoundaryConditionType : std::uint8_t
{
    Axis,                 ///< Zero integrated flux; both face endpoints must lie on r=0.
    SlipWall,             ///< Exact normal wall Riemann flux; free tangential velocity, zero mass/energy flux.
    Transmissive,         ///< Interior extrapolation, not a general characteristic BC.
    TotalConditionsInlet, ///< Subsonic normal inflow: P0, T0, direction; one outgoing invariant.
    MassFlowInlet,        ///< Subsonic normal inflow: inward group mdot, T0, direction.
    SupersonicInlet,      ///< Full static primitive state, strictly supersonic normal inflow.
    PressureOutlet,       ///< Static absolute pressure at subsonic outflow; ignored at un >= a.
    SupersonicOutlet      ///< Interior extrapolation, requiring un >= a.
};

/// Immutable physical boundary data. Factories validate finite positive inputs
/// and normalize directions once. Gas/regime/face compatibility is checked by
/// the flux evaluator. No chemical composition or time-dependent state is stored.
class EulerBoundaryCondition
{
  public:
    /// Legacy data-free conditions and SupersonicOutlet only.
    /// @throws std::invalid_argument For a type requiring physical data or an invalid enum.
    EulerBoundaryCondition(EulerBoundaryConditionType type);

    /// P0 [Pa absolute], T0 [K], meridional flow direction (normalized internally).
    [[nodiscard]] static EulerBoundaryCondition total_conditions_inlet(double total_pressure, double total_temperature,
                                                                       Vector2 direction);
    /// Positive inward group mdot [kg/s], T0 [K], direction. Uniform normal
    /// mass flux per physical area; no independently prescribed pressure.
    [[nodiscard]] static EulerBoundaryCondition mass_flow_inlet(double inward_mass_flow, double total_temperature,
                                                                Vector2 direction);
    /// Complete static primitive state; EOS and normal regime checked at evaluation.
    [[nodiscard]] static EulerBoundaryCondition supersonic_inlet(PrimitiveState state);
    /// Positive static absolute pressure [Pa], never a gauge pressure.
    [[nodiscard]] static EulerBoundaryCondition pressure_outlet(double pressure);

    [[nodiscard]] EulerBoundaryConditionType type() const noexcept
    {
        return type_;
    }

    [[nodiscard]] double total_pressure() const noexcept
    {
        return pressure_;
    }

    [[nodiscard]] double static_pressure() const noexcept
    {
        return pressure_;
    }

    [[nodiscard]] double total_temperature() const noexcept
    {
        return temperature_;
    }

    [[nodiscard]] double inward_mass_flow() const noexcept
    {
        return mass_flow_;
    }

    [[nodiscard]] Vector2 direction() const noexcept
    {
        return direction_;
    }

    [[nodiscard]] PrimitiveState prescribed_state() const noexcept
    {
        return state_;
    }

  private:
    EulerBoundaryCondition(EulerBoundaryConditionType type, double pressure, double temperature, double mass_flow,
                           Vector2 direction, PrimitiveState state);
    EulerBoundaryConditionType type_;
    double pressure_{};
    double temperature_{};
    double mass_flow_{};
    Vector2 direction_{};
    PrimitiveState state_{};
};

/// Owning fixed collection indexed by BoundaryId. Every group is explicit.
/// No Mesh is owned or referenced; the solver checks geometric compatibility.
class EulerBoundaryConditions
{
  public:
    /// @throws std::invalid_argument For missing conditions or unsupported types.
    EulerBoundaryConditions(Index boundary_count, std::vector<EulerBoundaryCondition> conditions)
        : conditions_(std::move(conditions))
    {
        if (conditions_.size() != boundary_count)
            throw std::invalid_argument("Euler boundary condition count must match boundary count.");
    }

    EulerBoundaryConditions(const EulerBoundaryConditions &) = default;
    EulerBoundaryConditions &operator=(const EulerBoundaryConditions &) = delete;
    EulerBoundaryConditions(EulerBoundaryConditions &&) noexcept = default;
    EulerBoundaryConditions &operator=(EulerBoundaryConditions &&) noexcept = delete;
    ~EulerBoundaryConditions() = default;

    [[nodiscard]] Index size() const noexcept
    {
        return conditions_.size();
    }

    /// @pre boundary < size().
    [[nodiscard]] EulerBoundaryConditionType operator[](BoundaryId boundary) const noexcept
    {
        return conditions_[boundary].type();
    }

    /// @pre boundary < size(). Read-only physical data; no allocation.
    [[nodiscard]] const EulerBoundaryCondition &condition(BoundaryId boundary) const noexcept
    {
        return conditions_[boundary];
    }

  private:
    std::vector<EulerBoundaryCondition> conditions_;
};

} // namespace cfd
