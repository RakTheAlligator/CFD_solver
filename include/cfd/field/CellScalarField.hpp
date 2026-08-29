#pragma once

#include "cfd/mesh/Types.hpp"

#include <span>
#include <vector>

namespace cfd
{

/// Scalar numerical field storing one contiguous value per mesh cell.
///
/// CellScalarField owns its values but neither owns nor references a Mesh. The
/// caller is responsible for constructing it with the appropriate cell count.
///
/// @note Spans returned by this class do not own their data and must not outlive
///       the CellScalarField storage from which they were obtained.
class CellScalarField
{
  public:
    /// Constructs a field with `cell_count` values initialized uniformly.
    explicit CellScalarField(Index cell_count, double initial_value = 0.0) : values_(cell_count, initial_value)
    {
    }

    CellScalarField(const CellScalarField &) = default;
    CellScalarField &operator=(const CellScalarField &) = delete;

    CellScalarField(CellScalarField &&) noexcept = default;
    CellScalarField &operator=(CellScalarField &&) noexcept = delete;

    ~CellScalarField() = default;

    /// Returns the number of cell values in the field.
    [[nodiscard]]
    Index size() const noexcept
    {
        return values_.size();
    }

    /// Returns the mutable value associated with `cell_id`.
    ///
    /// @pre `cell_id < size()`.
    [[nodiscard]]
    double &operator[](Index cell_id) noexcept
    {
        return values_[cell_id];
    }

    /// Returns the value associated with `cell_id`.
    ///
    /// @pre `cell_id < size()`.
    [[nodiscard]]
    const double &operator[](Index cell_id) const noexcept
    {
        return values_[cell_id];
    }

    /// Returns a non-owning mutable view of all cell values.
    [[nodiscard]]
    std::span<double> values() noexcept
    {
        return values_;
    }

    /// Returns a non-owning read-only view of all cell values.
    [[nodiscard]]
    std::span<const double> values() const noexcept
    {
        return values_;
    }

  private:
    std::vector<double> values_;
};

} // namespace cfd
