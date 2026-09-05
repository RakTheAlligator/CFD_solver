#pragma once

#include "cfd/field/CellScalarField.hpp"
#include "cfd/math/Vector2.hpp"
#include "cfd/mesh/Types.hpp"

namespace cfd
{

/// Two-dimensional cell-centered velocity field with component-wise storage.
///
/// The representation is Structure of Arrays: all `u` values are contiguous
/// in one CellScalarField and all `v` values are contiguous in another. This
/// matches the component-wise assembly and solution expected for the momentum
/// equations, allowing existing scalar operators to consume either component
/// directly without conversion or temporary fields. It is not a claim that
/// SoA is universally faster than an interleaved representation.
///
/// CellVelocityField owns exactly two dynamic component arrays but neither owns
/// nor references a Mesh. For nonzero cardinality, construction performs one
/// allocation through each component field (two total); accessors do not
/// allocate.
///
/// @note Copy construction is an intentional deep O(N) copy that allocates and
///       copies both component arrays. Numerical hot paths should pass velocity
///       fields by reference.
class CellVelocityField
{
  public:
    /// Constructs uniformly initialized `u` and `v` component fields.
    explicit CellVelocityField(Index cell_count, Vector2 initial_value = {})
        : u_(cell_count, initial_value.x), v_(cell_count, initial_value.y)
    {
    }

    CellVelocityField(const CellVelocityField &) = default;
    CellVelocityField &operator=(const CellVelocityField &) = delete;

    CellVelocityField(CellVelocityField &&) noexcept = default;
    CellVelocityField &operator=(CellVelocityField &&) noexcept = delete;

    ~CellVelocityField() = default;

    /// Returns the common number of cell-centered component values.
    [[nodiscard]]
    Index size() const noexcept
    {
        return u_.size();
    }

    /// Returns the mutable x-velocity component field.
    [[nodiscard]]
    CellScalarField &u() noexcept
    {
        return u_;
    }

    /// Returns the read-only x-velocity component field.
    [[nodiscard]]
    const CellScalarField &u() const noexcept
    {
        return u_;
    }

    /// Returns the mutable y-velocity component field.
    [[nodiscard]]
    CellScalarField &v() noexcept
    {
        return v_;
    }

    /// Returns the read-only y-velocity component field.
    [[nodiscard]]
    const CellScalarField &v() const noexcept
    {
        return v_;
    }

  private:
    CellScalarField u_;
    CellScalarField v_;
};

} // namespace cfd
