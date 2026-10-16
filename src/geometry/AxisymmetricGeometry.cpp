#include "cfd/geometry/AxisymmetricGeometry.hpp"

#include "cfd/mesh/Mesh.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>

namespace cfd
{
namespace
{

[[noreturn]]
void throw_geometry_error(const std::string &message)
{
    throw std::runtime_error("Axisymmetric geometry: " + message);
}

void validate_cell_identities(const Mesh &mesh, const std::span<const Vector2> face_vectors,
                              const std::span<const double> radial_sources)
{
    constexpr double tolerance{128.0 * std::numeric_limits<double>::epsilon()};
    for (Index cell = 0; cell < mesh.cell_count(); ++cell)
    {
        const Index begin{mesh.cell_node_offsets()[cell]}, end{mesh.cell_node_offsets()[cell + 1]};
        Vector2 scale{0.0, radial_sources[cell]};
        for (Index local = begin; local < end; ++local)
        {
            const Vector2 &area{face_vectors[mesh.cell_faces()[local]]};
            scale.x = std::max(scale.x, std::abs(area.x));
            scale.y = std::max(scale.y, std::abs(area.y));
        }

        // Normalize before summation to avoid overflow in absolute face sums.
        // The radial source provides a strictly positive radial scale; a zero
        // axial scale means every axial component is exactly zero.
        Vector2 balance{}, absolute_sum{};
        for (Index local = begin; local < end; ++local)
        {
            const Index face{mesh.cell_faces()[local]};
            const double sign{mesh.face_adjacencies()[face].owner == cell ? 1.0 : -1.0};
            const Vector2 &area{face_vectors[face]};
            if (scale.x > 0.0)
            {
                const double axial{area.x / scale.x};
                balance.x += sign * axial;
                absolute_sum.x += std::abs(axial);
            }
            const double radial{area.y / scale.y};
            balance.y += sign * radial;
            absolute_sum.y += std::abs(radial);
        }
        const double source{radial_sources[cell] / scale.y};
        if (std::abs(balance.x) > tolerance * absolute_sum.x ||
            std::abs(balance.y - source) > tolerance * (absolute_sum.y + source))
        {
            throw_geometry_error("cell " + std::to_string(cell) + " violates the integrated face/source identities.");
        }
    }
}

} // namespace

AxisymmetricGeometry::AxisymmetricGeometry(const Mesh &mesh)
    : mesh_(&mesh), cell_volumes_(mesh.cell_count()), face_area_vectors_(mesh.face_count()),
      cell_radial_source_coefficients_(mesh.cell_count())
{
    constexpr double two_pi{2.0 * std::numbers::pi};
    for (Index node = 0; node < mesh.node_count(); ++node)
    {
        const auto &point{mesh.nodes()[node]};
        if (!std::isfinite(point.x) || !std::isfinite(point.y) || point.y < 0.0)
        {
            throw_geometry_error("node " + std::to_string(node) + " requires finite coordinates and radius r >= 0.");
        }
    }
    for (Index cell = 0; cell < mesh.cell_count(); ++cell)
    {
        const auto &center{mesh.cell_centers()[cell]};
        const double area{mesh.cell_areas()[cell]};
        const double source{two_pi * area};
        const double volume{source * center.y};
        if (!std::isfinite(center.x) || !std::isfinite(center.y) || center.y < 0.0 || !std::isfinite(source) ||
            !(source > 0.0) || !std::isfinite(volume) || !(volume > 0.0))
        {
            throw_geometry_error("cell " + std::to_string(cell) + " has invalid radius, volume or radial source.");
        }
        cell_volumes_[cell] = volume;
        cell_radial_source_coefficients_[cell] = source;
    }
    for (Index face = 0; face < mesh.face_count(); ++face)
    {
        const auto &center{mesh.face_centers()[face]};
        const auto &planar_area{mesh.face_area_vectors()[face]};
        const double circumference{two_pi * center.y};
        const Vector2 area{circumference * planar_area.x, circumference * planar_area.y};
        const double magnitude{std::hypot(area.x, area.y)};
        const auto &nodes{mesh.faces()[face].node_ids};
        const bool on_axis{mesh.nodes()[nodes[0]].y == 0.0 && mesh.nodes()[nodes[1]].y == 0.0};
        if (!std::isfinite(center.x) || !std::isfinite(center.y) || center.y < 0.0 || !std::isfinite(magnitude) ||
            (magnitude == 0.0 && !on_axis) || (on_axis && magnitude != 0.0))
        {
            throw_geometry_error("face " + std::to_string(face) + " has invalid revolved area (zero only on r=0).");
        }
        face_area_vectors_[face] = area;
    }
    validate_cell_identities(mesh, face_area_vectors_, cell_radial_source_coefficients_);
}

} // namespace cfd
