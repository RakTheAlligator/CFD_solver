#pragma once

#include "cfd/math/Point2.hpp"

#include <filesystem>
#include <optional>
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
    std::optional<double> inlet_position{};
    std::optional<double> outer_upstream_length{};
    /// Radial offset of the fictitious wall, not a mechanical wall design.
    std::optional<double> wall_thickness{};
    std::optional<double> downstream_length{};
    std::optional<double> farfield_radius{};
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
};

/// Imports CAD using the configured Python/ezdxf adapter, validates r(x), and
/// applies the optional domain extents. No meshing is performed here.
/// @throws std::invalid_argument If domain options are inconsistent.
/// @throws std::runtime_error If CAD import, units or contour validation fails.
[[nodiscard]]
RpaNozzleDomain read_rpa_nozzle_domain(const RpaNozzleGeometry &geometry);

} // namespace cfd
