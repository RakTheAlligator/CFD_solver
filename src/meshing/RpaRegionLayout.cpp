#include "cfd/meshing/RpaRegionLayout.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cfd
{
RpaRegionLayout::RpaRegionLayout(const RpaNozzleDomain &domain, const RpaRegionLayoutOptions &options)
    : domain_(&domain)
{
    if (domain.inner_wall.size() < 3 || domain.outer_wall.size() < 2)
        throw std::invalid_argument("RPA layout requires a validated engine/shell profile.");
    band_thickness_ = options.band_thickness_factor * domain.inner_wall.back().y;
    plume_inlet_radius_ = domain.outer_wall.back().y + band_thickness_;
    plume_outlet_radius_ = options.plume_outlet_radius_factor * domain.inner_wall.back().y;
    if (!std::isfinite(band_thickness_) || !(band_thickness_ > 0.0) || !std::isfinite(plume_outlet_radius_) ||
        !(plume_inlet_radius_ < plume_outlet_radius_) || !(plume_outlet_radius_ < domain.farfield_radius))
        throw std::invalid_argument(
            "RPA regions require positive band thickness and plume inlet < outlet < box radius.");
    for (const auto &p : domain.outer_wall)
        if (!(p.y + band_thickness_ < plume_outlet_radius_))
            throw std::invalid_argument(
                "RPA BAND intersects FARFIELD; increase plumeOutletRadiusFactor or reduce band thickness.");
}

double RpaRegionLayout::plume_radius(double x) const noexcept
{
    const double t{std::clamp(
        (x - domain_->inner_wall.back().x) / (domain_->downstream_x - domain_->inner_wall.back().x), 0.0, 1.0)};
    return std::lerp(plume_inlet_radius_, plume_outlet_radius_, t);
}

RpaRegion RpaRegionLayout::region(const Point2 &p) const
{
    const auto &d{*domain_};
    const double inner{interpolate_rpa_profile(d.inner_wall, p.x).y};
    const bool upstream{p.x < d.inner_wall.back().x};
    if (p.x < d.inner_wall.front().x || p.x > d.downstream_x || p.y < 0.0 || p.y > d.farfield_radius ||
        (p.x < d.outer_wall.front().x && p.y > inner) ||
        (upstream && p.x >= d.outer_wall.front().x && p.y > inner &&
         p.y < interpolate_rpa_profile(d.outer_wall, p.x).y))
        throw std::runtime_error("RPA point lies outside fluid domain or inside fictitious shell.");
    if (upstream && p.y <= inner)
        return RpaRegion::Engine;
    if (upstream && p.y <= interpolate_rpa_profile(d.outer_wall, p.x).y + band_thickness())
        return RpaRegion::Band;
    if (p.y >= plume_outlet_radius())
        return RpaRegion::Farfield;
    if (!upstream && p.y <= plume_radius(p.x))
        return RpaRegion::Plume;
    return RpaRegion::Transition;
}
} // namespace cfd
