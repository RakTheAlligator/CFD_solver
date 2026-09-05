#pragma once

#include "cfd/field/CellScalarField.hpp"
#include "cfd/mesh/Types.hpp"

namespace cfd
{

/// Cell-centered diagonal momentum response to a pressure gradient.
///
/// For each cell `P`, this field stores the two independent coefficients
/// `d_P,u = A_P / a_P,u` and `d_P,v = A_P / a_P,v`, where `A_P` is the
/// two-dimensional cell area per unit depth and each `a_P` is the corresponding
/// final assembled and equation-under-relaxed momentum diagonal.
///
/// These values are pressure-response coefficients, not velocity components.
/// Their Structure-of-Arrays representation preserves independent contiguous
/// storage for the segregated Cartesian momentum equations and is intended for
/// future pressure-velocity coupling and momentum interpolation.
///
/// CellMomentumPressureResponse owns both scalar arrays but neither owns nor
/// references a Mesh. For nonzero cardinality, construction performs one
/// allocation through each component field; accessors do not allocate.
///
/// @note Copy construction is a deep O(N) copy that allocates and copies both
///       response arrays. Numerical hot paths should pass this field by
///       reference.
class CellMomentumPressureResponse
{
  public:
    /// Constructs zero-initialized `u` and `v` pressure-response arrays.
    explicit CellMomentumPressureResponse(Index cell_count) : u_response_(cell_count), v_response_(cell_count)
    {
    }

    CellMomentumPressureResponse(const CellMomentumPressureResponse &) = default;
    CellMomentumPressureResponse &operator=(const CellMomentumPressureResponse &) = delete;

    CellMomentumPressureResponse(CellMomentumPressureResponse &&) noexcept = default;
    CellMomentumPressureResponse &operator=(CellMomentumPressureResponse &&) noexcept = delete;

    ~CellMomentumPressureResponse() = default;

    /// Returns the common number of cell-centered response values.
    [[nodiscard]]
    Index size() const noexcept
    {
        return u_response_.size();
    }

    /// Returns the mutable response associated with the u-momentum diagonal.
    [[nodiscard]]
    CellScalarField &u() noexcept
    {
        return u_response_;
    }

    /// Returns the read-only response associated with the u-momentum diagonal.
    [[nodiscard]]
    const CellScalarField &u() const noexcept
    {
        return u_response_;
    }

    /// Returns the mutable response associated with the v-momentum diagonal.
    [[nodiscard]]
    CellScalarField &v() noexcept
    {
        return v_response_;
    }

    /// Returns the read-only response associated with the v-momentum diagonal.
    [[nodiscard]]
    const CellScalarField &v() const noexcept
    {
        return v_response_;
    }

  private:
    CellScalarField u_response_;
    CellScalarField v_response_;
};

} // namespace cfd
