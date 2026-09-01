#pragma once

#include "cfd/math/Vector2.hpp"
#include "cfd/mesh/Types.hpp"

#include <span>
#include <vector>

namespace cfd
{

/// Vector numerical field storing one contiguous value per mesh cell.
///
/// CellVectorField owns its fixed-cardinality storage but neither owns nor
/// references a Mesh. The caller is responsible for constructing it with the
/// appropriate cell count.
///
/// @note Spans returned by this class do not own their data and must not outlive
///       the CellVectorField storage from which they were obtained.
class CellVectorField
{
  public:
    /// Constructs a field with `cell_count` values initialized uniformly.
    explicit CellVectorField(Index cell_count, Vector2 initial_value = {}) : values_(cell_count, initial_value)
    {
    }

    CellVectorField(const CellVectorField &) = default;
    CellVectorField &operator=(const CellVectorField &) = delete;

    CellVectorField(CellVectorField &&) noexcept = default;
    CellVectorField &operator=(CellVectorField &&) noexcept = delete;

    ~CellVectorField() = default;

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
    Vector2 &operator[](Index cell_id) noexcept
    {
        return values_[cell_id];
    }

    /// Returns the value associated with `cell_id`.
    ///
    /// @pre `cell_id < size()`.
    [[nodiscard]]
    const Vector2 &operator[](Index cell_id) const noexcept
    {
        return values_[cell_id];
    }

    /// Returns a non-owning mutable view of all cell values.
    [[nodiscard]]
    std::span<Vector2> values() noexcept
    {
        return values_;
    }

    /// Returns a non-owning read-only view of all cell values.
    [[nodiscard]]
    std::span<const Vector2> values() const noexcept
    {
        return values_;
    }

  private:
    std::vector<Vector2> values_;
};

} // namespace cfd
