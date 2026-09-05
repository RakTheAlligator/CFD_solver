#pragma once

#include "cfd/mesh/Types.hpp"

#include <span>
#include <vector>

namespace cfd
{

/// Integrated mass-flux response to a pressure difference per mesh face.
///
/// On an internal face, a value `Dp_f` maps the owner-neighbor pressure
/// difference to the owner-oriented contribution `-Dp_f * (p_N - p_P)`.
/// Boundary values are not defined by the internal-face Rhie-Chow component
/// and remain under caller ownership.
///
/// FacePressureResponseField owns its contiguous fixed-cardinality storage but
/// neither owns nor references a Mesh. The caller supplies the appropriate
/// face count.
///
/// @note Spans returned by this class do not own their data and must not outlive
///       the FacePressureResponseField storage from which they were obtained.
class FacePressureResponseField
{
  public:
    /// Constructs a field with `face_count` values initialized uniformly.
    explicit FacePressureResponseField(Index face_count, double initial_value = 0.0)
        : values_(face_count, initial_value)
    {
    }

    FacePressureResponseField(const FacePressureResponseField &) = default;
    FacePressureResponseField &operator=(const FacePressureResponseField &) = delete;

    FacePressureResponseField(FacePressureResponseField &&) noexcept = default;
    FacePressureResponseField &operator=(FacePressureResponseField &&) noexcept = delete;

    ~FacePressureResponseField() = default;

    /// Returns the number of face pressure-response values.
    [[nodiscard]]
    Index size() const noexcept
    {
        return values_.size();
    }

    /// Returns the mutable pressure response associated with `face_id`.
    ///
    /// @pre `face_id < size()`.
    [[nodiscard]]
    double &operator[](Index face_id) noexcept
    {
        return values_[face_id];
    }

    /// Returns the pressure response associated with `face_id`.
    ///
    /// @pre `face_id < size()`.
    [[nodiscard]]
    const double &operator[](Index face_id) const noexcept
    {
        return values_[face_id];
    }

    /// Returns a non-owning mutable view of all face pressure responses.
    [[nodiscard]]
    std::span<double> values() noexcept
    {
        return values_;
    }

    /// Returns a non-owning read-only view of all face pressure responses.
    [[nodiscard]]
    std::span<const double> values() const noexcept
    {
        return values_;
    }

  private:
    std::vector<double> values_;
};

} // namespace cfd
