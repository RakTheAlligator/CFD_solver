#pragma once

#include "cfd/mesh/Types.hpp"

#include <span>
#include <vector>

namespace cfd
{

/// Integrated signed convective carrier flux stored per global mesh face.
///
/// Flux orientation follows the Mesh owner orientation. On an internal face,
/// a positive value denotes transport from owner to neighbor and a negative
/// value denotes transport from neighbor to owner. On a boundary face,
/// positive denotes domain outflow and negative denotes domain inflow.
///
/// Units depend on the transport formulation, for example volumetric or mass
/// flux. The convection operator requires only that they are used consistently.
/// FaceFluxField owns its contiguous values but neither owns nor references a
/// Mesh; the caller supplies the appropriate fixed face count.
///
/// @note Spans returned by this class do not own their data and must not outlive
///       the FaceFluxField storage from which they were obtained.
class FaceFluxField
{
  public:
    /// Constructs a field with `face_count` fluxes initialized uniformly.
    explicit FaceFluxField(Index face_count, double initial_value = 0.0) : values_(face_count, initial_value)
    {
    }

    FaceFluxField(const FaceFluxField &) = default;
    FaceFluxField &operator=(const FaceFluxField &) = delete;

    FaceFluxField(FaceFluxField &&) noexcept = default;
    FaceFluxField &operator=(FaceFluxField &&) noexcept = delete;

    ~FaceFluxField() = default;

    /// Returns the number of face fluxes.
    [[nodiscard]]
    Index size() const noexcept
    {
        return values_.size();
    }

    /// Returns the mutable flux associated with `face_id`.
    ///
    /// @pre `face_id < size()`.
    [[nodiscard]]
    double &operator[](Index face_id) noexcept
    {
        return values_[face_id];
    }

    /// Returns the flux associated with `face_id`.
    ///
    /// @pre `face_id < size()`.
    [[nodiscard]]
    const double &operator[](Index face_id) const noexcept
    {
        return values_[face_id];
    }

    /// Returns a non-owning mutable view of all face fluxes.
    [[nodiscard]]
    std::span<double> values() noexcept
    {
        return values_;
    }

    /// Returns a non-owning read-only view of all face fluxes.
    [[nodiscard]]
    std::span<const double> values() const noexcept
    {
        return values_;
    }

  private:
    std::vector<double> values_;
};

} // namespace cfd
