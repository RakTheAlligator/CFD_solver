#pragma once

#include "cfd/mesh/Boundary.hpp"
#include "cfd/mesh/Types.hpp"

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cfd
{

/// Mathematical type of a scalar boundary condition.
enum class ScalarBoundaryConditionType : std::uint8_t
{
    /// Prescribes the scalar value: `phi = value`.
    Dirichlet,

    /// Prescribes its outward normal derivative: `d(phi)/dn = value`.
    Neumann
};

/// Constant mathematical condition applied to a scalar field boundary.
///
/// For Dirichlet conditions, `value` is the prescribed value of `phi`. For
/// Neumann conditions, it is the prescribed derivative `d(phi)/dn`, where `n`
/// is the outward boundary normal.
struct ScalarBoundaryCondition
{
    ScalarBoundaryCondition(ScalarBoundaryConditionType condition_type, double condition_value) noexcept
        : type(condition_type), value(condition_value)
    {
    }

    ScalarBoundaryConditionType type;
    double value;
};

/// Fixed collection of scalar conditions indexed directly by BoundaryId.
///
/// ScalarBoundaryConditions owns its conditions but neither owns nor references
/// a Mesh. The Mesh is used by the caller only to determine `boundary_count`
/// during construction.
class ScalarBoundaryConditions
{
  public:
    /// Constructs and validates one condition per boundary group.
    ///
    /// @throws std::invalid_argument If the number of conditions differs from
    ///         `boundary_count`, a condition type is unsupported, or a value is
    ///         non-finite.
    ScalarBoundaryConditions(Index boundary_count, std::vector<ScalarBoundaryCondition> conditions)
        : conditions_(std::move(conditions))
    {
        if (conditions_.size() != boundary_count)
        {
            throw std::invalid_argument("Scalar boundary condition count must match boundary count.");
        }

        for (const ScalarBoundaryCondition &condition : conditions_)
        {
            switch (condition.type)
            {
            case ScalarBoundaryConditionType::Dirichlet:
            case ScalarBoundaryConditionType::Neumann:
                break;

            default:
                throw std::invalid_argument("Scalar boundary condition type is unsupported.");
            }

            if (!std::isfinite(condition.value))
            {
                throw std::invalid_argument("Scalar boundary condition value must be finite.");
            }
        }
    }

    ScalarBoundaryConditions(const ScalarBoundaryConditions &) = default;
    ScalarBoundaryConditions &operator=(const ScalarBoundaryConditions &) = delete;

    ScalarBoundaryConditions(ScalarBoundaryConditions &&) noexcept = default;
    ScalarBoundaryConditions &operator=(ScalarBoundaryConditions &&) noexcept = delete;

    ~ScalarBoundaryConditions() = default;

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
    const ScalarBoundaryCondition &operator[](BoundaryId boundary_id) const noexcept
    {
        return conditions_[boundary_id];
    }

  private:
    std::vector<ScalarBoundaryCondition> conditions_;
};

} // namespace cfd
