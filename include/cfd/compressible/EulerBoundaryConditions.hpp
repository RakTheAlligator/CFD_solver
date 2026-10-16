#pragma once

#include "cfd/mesh/Boundary.hpp"

#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cfd
{

enum class EulerBoundaryConditionType : std::uint8_t
{
    Axis,        ///< Zero integrated flux; both face endpoints must lie on r=0.
    SlipWall,    ///< Exact normal wall Riemann flux; free tangential velocity, zero mass/energy flux.
    Transmissive ///< Interior extrapolation, not a general characteristic BC.
};

/// Owning fixed collection indexed by BoundaryId. Every group is explicit.
/// No Mesh is owned or referenced; the solver checks geometric compatibility.
class EulerBoundaryConditions
{
  public:
    /// @throws std::invalid_argument For missing conditions or unsupported types.
    EulerBoundaryConditions(Index boundary_count, std::vector<EulerBoundaryConditionType> conditions)
        : conditions_(std::move(conditions))
    {
        if (conditions_.size() != boundary_count)
            throw std::invalid_argument("Euler boundary condition count must match boundary count.");
        for (const auto condition : conditions_)
        {
            switch (condition)
            {
            case EulerBoundaryConditionType::Axis:
            case EulerBoundaryConditionType::SlipWall:
            case EulerBoundaryConditionType::Transmissive:
                break;
            default:
                throw std::invalid_argument("Unsupported Euler boundary condition type.");
            }
        }
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
        return conditions_[boundary];
    }

  private:
    std::vector<EulerBoundaryConditionType> conditions_;
};

} // namespace cfd
