#pragma once

#include "cfd/mesh/Boundary.hpp"
#include "cfd/mesh/Types.hpp"

#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cfd
{

/// Algorithmic pressure-correction condition applied to a boundary group.
enum class PressureCorrectionBoundaryConditionType : std::uint8_t
{
    /// The boundary mass flux is imposed and is not pressure-corrected:
    /// `F'_b = 0`.
    FixedMassFlux,

    /// The physical pressure is prescribed, hence `p'_b = 0`.
    FixedPressure
};

/// Fixed collection of pressure-correction conditions indexed by BoundaryId.
///
/// PressureCorrectionBoundaryConditions owns its conditions but neither owns
/// nor references a Mesh. The Mesh is used by the caller only to determine
/// `boundary_count` during construction.
class PressureCorrectionBoundaryConditions
{
  public:
    /// Constructs and validates one condition per boundary group.
    ///
    /// @throws std::invalid_argument If the number of conditions differs from
    ///         `boundary_count`, or a condition type is unsupported.
    PressureCorrectionBoundaryConditions(Index boundary_count,
                                         std::vector<PressureCorrectionBoundaryConditionType> conditions)
        : conditions_(std::move(conditions))
    {
        if (conditions_.size() != boundary_count)
        {
            throw std::invalid_argument("Pressure-correction boundary condition count must match boundary count.");
        }

        for (const PressureCorrectionBoundaryConditionType condition : conditions_)
        {
            switch (condition)
            {
            case PressureCorrectionBoundaryConditionType::FixedMassFlux:
            case PressureCorrectionBoundaryConditionType::FixedPressure:
                break;

            default:
                throw std::invalid_argument("Pressure-correction boundary condition type is unsupported.");
            }
        }
    }

    PressureCorrectionBoundaryConditions(const PressureCorrectionBoundaryConditions &) = default;
    PressureCorrectionBoundaryConditions &operator=(const PressureCorrectionBoundaryConditions &) = delete;

    PressureCorrectionBoundaryConditions(PressureCorrectionBoundaryConditions &&) noexcept = default;
    PressureCorrectionBoundaryConditions &operator=(PressureCorrectionBoundaryConditions &&) noexcept = delete;

    ~PressureCorrectionBoundaryConditions() = default;

    /// Returns the number of boundary conditions.
    [[nodiscard]]
    Index size() const noexcept
    {
        return conditions_.size();
    }

    /// Returns the condition associated with `boundary_id`.
    ///
    /// @pre `boundary_id < size()`.
    [[nodiscard]]
    PressureCorrectionBoundaryConditionType operator[](BoundaryId boundary_id) const noexcept
    {
        return conditions_[boundary_id];
    }

  private:
    std::vector<PressureCorrectionBoundaryConditionType> conditions_;
};

} // namespace cfd
