#pragma once

#include "cfd/meshing/RpaNozzleGeometry.hpp"

#include <cstdint>

namespace cfd
{
/// Five fluid regions; IDs are also used by the VTU rpa_region array.
enum class RpaRegion : std::uint8_t
{
    Engine = 1,
    Band,
    Plume,
    Transition,
    Farfield
};

/// Geometric interfaces only: no cell counts or resolution policies.
/// Non-owning: domain must outlive the layout.
class RpaRegionLayout
{
  public:
    explicit RpaRegionLayout(const RpaNozzleDomain &domain, const RpaRegionLayoutOptions &options = {});

    [[nodiscard]] double band_thickness() const noexcept
    {
        return band_thickness_;
    }

    [[nodiscard]] double plume_inlet_radius() const noexcept
    {
        return plume_inlet_radius_;
    }

    [[nodiscard]] double plume_outlet_radius() const noexcept
    {
        return plume_outlet_radius_;
    }

    [[nodiscard]] double plume_radius(double x) const noexcept;
    /// Rejects points in the fictitious solid or outside the fluid domain.
    [[nodiscard]] RpaRegion region(const Point2 &point) const;

  private:
    const RpaNozzleDomain *domain_;
    double band_thickness_{};
    double plume_inlet_radius_{};
    double plume_outlet_radius_{};
};
} // namespace cfd
