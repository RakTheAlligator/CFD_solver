#pragma once

#include "cfd/math/Point2.hpp"

#include <filesystem>
#include <optional>
#include <span>
#include <string>
#include <vector>

namespace cfd
{

/// RPA meridian CAD input; all lengths other than DXF coordinates are metres.
/// The reader resolves a relative file path against the case directory.
struct RpaNozzleGeometry
{
    std::filesystem::path file;
    double geometric_tolerance{1.0e-6};
    /// Overrides absent DXF units: auto, mm, cm, m or in. Conflicts are rejected.
    std::string dxf_units{"auto"};
    /// Radial offset of the fictitious wall, not a mechanical wall design.
    double wall_thickness_factor{0.02};
    /// Absent means the complete available engine length.
    std::optional<double> upstream_factor{};
    double downstream_factor{80.0};
    double outer_radius_factor{65.0};
};

/// Geometry of fluid-region interfaces, independent of cell counts.
struct RpaRegionLayoutOptions
{
    double band_thickness_factor{0.5};
    double plume_outlet_radius_factor{16.0};
};

/// Validated, flattened meridian and fictitious exterior wall in Cartesian (x,r).
/// Interior/exterior fluid communicate downstream of exit; there are no solids
/// or axisymmetric finite-volume metrics in this representation.
struct RpaNozzleDomain
{
    std::vector<Point2> inner_wall;
    std::vector<Point2> outer_wall;
    /// Source-entity junctions remain mandatory geometric vertices in Gmsh.
    std::vector<double> curve_junctions;
    Point2 throat;
    double downstream_x{};
    double farfield_radius{};
    double geometric_tolerance{1.0e-6};
};

/// Linear interpolation of the validated, increasing-x flattened CAD profile.
[[nodiscard]] Point2 interpolate_rpa_profile(std::span<const Point2> profile, double x);

/// Imports CAD using the configured Python/ezdxf adapter, validates r(x), and
/// applies the optional domain extents. No meshing is performed here.
/// @throws std::invalid_argument If domain options are inconsistent.
/// @throws std::runtime_error If CAD import, units or contour validation fails.
[[nodiscard]]
RpaNozzleDomain read_rpa_nozzle_domain(const RpaNozzleGeometry &geometry);

} // namespace cfd
