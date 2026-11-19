#include "RpaMeshDiscretization.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace cfd::detail
{
namespace
{
constexpr double maximum_growth{1.2};
constexpr Index minimum_engine_rows{8}, minimum_band_rows{4};

/// Miter intersections of adjacent parallel offsets. End caps retain x.
/// Distances are normal to the flattened CAD; physical walls are unchanged.
Point2 wall_offset_point(const std::vector<Point2> &profile, double x, double distance)
{
    const auto right{std::upper_bound(profile.begin(), profile.end(), x,
                                      [](double value, const Point2 &point) { return value < point.x; })};
    const auto end{right == profile.begin() ? std::next(right) : right == profile.end() ? std::prev(right) : right};
    const auto index{static_cast<Index>(end - profile.begin())};
    const auto normal = [&](Index a, Index b) {
        const double dx{profile[b].x - profile[a].x}, dy{profile[b].y - profile[a].y};
        const double length{std::hypot(dx, dy)};
        return Point2{-dy / length, dx / length};
    };
    const auto offset = [&](Index vertex) {
        const auto &point{profile[vertex]};
        if (vertex == 0 || vertex + 1 == profile.size())
        {
            const auto n{vertex == 0 ? normal(0, 1) : normal(vertex - 1, vertex)};
            return Point2{point.x, point.y + distance / n.y};
        }
        const auto left{normal(vertex - 1, vertex)}, next{normal(vertex, vertex + 1)};
        const double denominator{1.0 + left.x * next.x + left.y * next.y};
        if (!(denominator > 0.0))
            throw std::invalid_argument("No usable normal at RPA wall cusp.");
        return Point2{point.x + distance * (left.x + next.x) / denominator,
                      point.y + distance * (left.y + next.y) / denominator};
    };
    const auto a{offset(index - 1)}, b{offset(index)};
    const double fraction{std::clamp((x - profile[index - 1].x) / (profile[index].x - profile[index - 1].x), 0.0, 1.0)};
    return {std::lerp(a.x, b.x, fraction), std::lerp(a.y, b.y, fraction)};
}

double layer_thickness(const WallLayerOptions &layers, double maximum_growth)
{
    if (!layers.active())
        return 0.0;
    if (layers.layer_count == 0 || layers.layer_count >= static_cast<Index>(std::numeric_limits<int>::max()) ||
        !std::isfinite(layers.growth_rate) || layers.growth_rate < 1.0 || layers.growth_rate > maximum_growth ||
        (!layers.first_cell_height && !layers.thickness))
        throw std::invalid_argument(
            "Enabled wall layers need a height/thickness, positive count and growth in [1,maximumGrowthRate].");
    const double weight{layers.growth_rate == 1.0
                            ? static_cast<double>(layers.layer_count)
                            : std::expm1(static_cast<double>(layers.layer_count) * std::log(layers.growth_rate)) /
                                  std::expm1(std::log(layers.growth_rate))};
    const double total{layers.first_cell_height ? *layers.first_cell_height * weight : *layers.thickness};
    if (!std::isfinite(weight) || !std::isfinite(total) || !(total > 0.0) ||
        (layers.first_cell_height &&
         (!std::isfinite(*layers.first_cell_height) || !(*layers.first_cell_height > 0.0))) ||
        (layers.thickness && (!std::isfinite(*layers.thickness) || !(*layers.thickness > 0.0) ||
                              std::abs(*layers.thickness - total) > 1e-10 * total)))
        throw std::invalid_argument("Wall-layer height/count/growth/thickness are invalid or inconsistent.");
    return total;
}

Index cell_count(double length, double size)
{
    const double n{std::ceil(length / size)};
    if (!std::isfinite(n) || n < 1.0 || n >= static_cast<double>(std::numeric_limits<int>::max() - 1))
        throw std::invalid_argument("RPA discretization exceeds Gmsh's supported index range.");
    return static_cast<Index>(n);
}

std::vector<double> fractions(Index count, double q)
{
    std::vector<double> result(count + 1);
    double weight{1.0}, total{};
    for (Index i = 0; i < count; ++i)
    {
        total += weight;
        weight *= q;
    }
    weight = 1.0;
    for (Index i = 1; i < count; ++i)
    {
        result[i] = result[i - 1] + weight / total;
        weight *= q;
    }
    result.back() = 1.0;
    return result;
}

// Capped geometric spacing, redistributed as a whole to avoid a tiny last cell.
std::vector<double> growing_positions(double length, double first, double cap, double growth)
{
    std::vector<double> result{0.0};
    double step{first};
    while (result.back() < length)
    {
        const double next{result.back() + step};
        if (!std::isfinite(next) || !(next > result.back()) ||
            result.size() >= static_cast<Index>(std::numeric_limits<int>::max() - 1))
            throw std::invalid_argument("RPA mapped spacing exceeds the supported range.");
        result.push_back(next);
        step = std::min(cap, step * growth);
    }
    const double total{result.back()};
    for (auto &x : result)
        x *= length / total;
    result.back() = length;
    return result;
}

void finite_farfield(RpaMeshDiscretization &mesh, const RpaNozzleDomain &domain, const RpaRegionLayout &layout,
                     double interface_size, double cap, double growth)
{
    auto &upper{mesh.farfield};
    const double right{layout.farfield_begin_x()}, low{layout.transition_outer_radius(right)};
    const auto heights{growing_positions(domain.farfield_radius - low, interface_size, cap, growth)};
    upper.radial_cells = heights.size() - 1;
    upper.x.push_back(domain.outer_wall.front().x);
    // Preserve the envelope's corners, independently of the inner block columns.
    for (const double end : {domain.inner_wall.back().x, layout.plume_end_x(), right})
    {
        const double start{upper.x.back()};
        const Index count{cell_count(end - start, interface_size)};
        for (Index i = 1; i <= count; ++i)
            upper.x.push_back(std::lerp(start, end, static_cast<double>(i) / static_cast<double>(count)));
    }
    for (double x : upper.x)
        for (double h : heights)
            upper.nodes.push_back(
                {x, std::lerp(layout.transition_outer_radius(x), domain.farfield_radius, h / heights.back())});

    auto &downstream{mesh.farfield_downstream};
    mesh.farfield_split_row = cell_count(low, interface_size);
    downstream.radial_cells = mesh.farfield_split_row + upper.radial_cells;
    const auto distances{growing_positions(domain.downstream_x - right, interface_size, cap, growth)};
    for (double dx : distances)
    {
        downstream.x.push_back(right + dx);
        for (Index j = 0; j <= mesh.farfield_split_row; ++j)
            downstream.nodes.push_back(
                {right + dx, low * static_cast<double>(j) / static_cast<double>(mesh.farfield_split_row)});
        for (Index j = 1; j <= upper.radial_cells; ++j)
            downstream.nodes.push_back({right + dx, upper.at(upper.x.size() - 1, j).y});
    }
}

struct RadialDistribution
{
    Index cells, core_cells;
    std::vector<double> core, distances;
    bool inner;
};

std::vector<double> multiblock_plume_positions(double length, double near_length, double first, double cap,
                                               double growth, std::vector<Index> &breaks)
{
    near_length = std::min(length, near_length);
    // A tail too short for two terminal cells is folded into the near block;
    // do not create a tiny sliver surface just to keep three named blocks.
    if (length - near_length < 2 * cap)
        near_length = length;
    const Index near_cells{std::max(Index{2}, cell_count(near_length, first))};
    std::vector<double> positions;
    positions.reserve(near_cells + 1);
    for (Index i = 0; i <= near_cells; ++i)
        positions.push_back(near_length * static_cast<double>(i) / static_cast<double>(near_cells));
    breaks = {0, near_cells};
    if (near_length < length)
    {
        // Redistribute the complete tail, never leave a tiny remainder cell.
        const double next{std::min(cap, (near_length / static_cast<double>(near_cells)) * growth)};
        const auto tail{growing_positions(length - near_length, next, cap, growth)};
        Index graded{};
        for (double step = next; step < cap && graded + 1 < tail.size(); step = std::min(cap, step * growth))
            ++graded;
        for (Index i = 1; i < tail.size(); ++i)
            positions.push_back(near_length + tail[i]);
        if (graded > 0 && graded + 1 < tail.size())
            breaks.push_back(near_cells + graded);
        breaks.push_back(positions.size() - 1);
    }
    positions.back() = length;
    return positions;
}

RadialDistribution radial_distribution(double extent, double base, const WallLayerOptions &layers, bool inner,
                                       std::optional<Index> requested = {})
{
    const double thickness{layer_thickness(layers, maximum_growth)};
    const Index count{
        requested.value_or(std::max({inner ? minimum_engine_rows : minimum_band_rows, cell_count(extent, base),
                                     layers.active() ? layers.layer_count + 2 : Index{0}}))};
    if (count < 2 || count >= static_cast<Index>(std::numeric_limits<int>::max() - 1) ||
        (layers.active() && count <= layers.layer_count + 1))
        throw std::invalid_argument("RPA radial cell counts must fit the supported range and prescribed layers.");
    const Index core_count{count - (layers.active() ? layers.layer_count : 0)};
    RadialDistribution result{
        count,
        core_count,
        fractions(core_count,
                  layers.active() ? 1.0 : std::pow(inner ? 0.25 : 4.0, 1.0 / static_cast<double>(core_count - 1))),
        {},
        inner};
    if (layers.active())
    {
        result.distances = fractions(layers.layer_count, layers.growth_rate);
        for (auto &d : result.distances)
            d *= thickness;
    }
    return result;
}

void column(const std::vector<Point2> &wall, double x, double band_height, const RadialDistribution &radial,
            std::span<Point2> output)
{
    const double radius{interpolate_rpa_profile(wall, x).y};
    const double thickness{radial.distances.empty() ? 0.0 : radial.distances.back()};
    const Point2 edge{thickness > 0.0 ? wall_offset_point(wall, x, radial.inner ? -thickness : thickness)
                                      : Point2{x, radius}};
    if (!std::isfinite(edge.x) || !std::isfinite(edge.y) ||
        (radial.inner ? !(edge.y > 0.0) : !(edge.y < radius + band_height)))
        throw std::invalid_argument("RPA prescribed wall layers do not fit inside ENGINE/BAND.");
    if (radial.inner)
    {
        for (Index j = 0; j <= radial.core_cells; ++j)
            output[j] = {x, radial.core[j] * edge.y};
        output[radial.core_cells] = edge;
        for (Index j = 1; j < radial.distances.size(); ++j)
            output[radial.core_cells + j] =
                wall_offset_point(wall, x, -radial.distances[radial.distances.size() - 1 - j]);
    }
    else
    {
        const Index layers{radial.cells - radial.core_cells};
        for (Index j = 0; j < radial.distances.size(); ++j)
            output[j] = wall_offset_point(wall, x, radial.distances[j]);
        for (Index j = 0; j <= radial.core_cells; ++j)
            output[layers + j] = {x, std::lerp(edge.y, radius + band_height, radial.core[j])};
        output[layers] = edge;
    }
}

// March with the arithmetic mean of ALL radial heights, then redistribute
// the interval to its mandatory endpoint. Geometry/quality never refine dx.
RpaStructuredGrid wall_grid(const RpaNozzleDomain &domain, const std::vector<Point2> &wall,
                            const RadialDistribution &radial, double band_height, double quality, bool even,
                            double axial_spacing_factor, bool preserve_vertices)
{
    RpaStructuredGrid grid{{}, radial.cells, {}};
    std::vector<double> stations{wall.front().x, wall.back().x};
    for (double x : domain.curve_junctions)
        if (x > stations.front() && x < wall.back().x)
            stations.push_back(x);
    if (preserve_vertices)
        for (const auto &point : wall)
            stations.push_back(point.x);
    if (domain.throat.x > wall.front().x && domain.throat.x < wall.back().x)
        stations.push_back(domain.throat.x);
    std::ranges::sort(stations);
    stations.erase(std::unique(stations.begin(), stations.end()), stations.end());
    std::vector<Point2> points(radial.cells + 1);
    std::vector<double> candidate;
    grid.x.push_back(wall.front().x);
    for (Index segment = 1; segment < stations.size(); ++segment)
    {
        const double a{stations[segment - 1]}, b{stations[segment]};
        candidate.assign(1, a);
        while (candidate.back() < b)
        {
            const double x{candidate.back()};
            column(wall, x, band_height, radial, points);
            const double step{axial_spacing_factor * (points.back().y - points.front().y) /
                              static_cast<double>(radial.cells)};
            const double next{x + step};
            if (!std::isfinite(step) || !(step > 0.0) || !std::isfinite(next) || !(next > x))
                throw std::runtime_error("RPA unusable mean radial height for axial spacing.");
            if (candidate.size() >= static_cast<Index>(std::numeric_limits<int>::max() - 1))
                throw std::invalid_argument("RPA column count exceeds Gmsh's supported index range.");
            candidate.push_back(next);
        }
        // Round to the nearer endpoint instead of leaving a tiny final cell.
        // The first interval is retained even if a CAD segment is very short.
        if (candidate.size() > 2 && b - candidate[candidate.size() - 2] < candidate.back() - b)
            candidate.pop_back();
        const Index nominal{candidate.size() - 1};
        const Index count{nominal + (even && nominal % 2 != 0 ? 1 : 0)};
        // Full-QUAD needs even counts per CAD curve. Resample the same march,
        // distributing the parity adjustment over the complete CAD interval.
        for (Index i = 1; i < count; ++i)
        {
            const double index{static_cast<double>(i) * static_cast<double>(nominal) / static_cast<double>(count)};
            const Index k{std::min(static_cast<Index>(index), nominal - 1)};
            const double x{std::lerp(candidate[k], candidate[k + 1], index - static_cast<double>(k))};
            grid.x.push_back(std::lerp(a, b, (x - a) / (candidate.back() - a)));
        }
        grid.x.push_back(b);
    }
    if (grid.x.size() > static_cast<Index>(std::numeric_limits<int>::max()) ||
        grid.x.size() > std::numeric_limits<Index>::max() / (radial.cells + 1))
        throw std::invalid_argument("RPA grid size overflow.");
    grid.nodes.resize(grid.x.size() * (radial.cells + 1));
    for (Index i = 0; i < grid.x.size(); ++i)
        column(wall, grid.x[i], band_height, radial,
               std::span{grid.nodes}.subspan(i * (radial.cells + 1), radial.cells + 1));
    // Validation only: never alter the requested spacing to repair a cell.
    for (Index i = 1; i < grid.x.size(); ++i)
        for (Index j = 0; j < radial.cells; ++j)
        {
            const std::array<Point2, 4> quad{grid.at(i - 1, j), grid.at(i, j), grid.at(i, j + 1),
                                             grid.at(i - 1, j + 1)};
            if (!(rpa_angular_quality(quad) >= quality))
            {
                std::ostringstream message;
                message << std::setprecision(17) << "RPA mean-height column violates minCellQuality at x=["
                        << grid.x[i - 1] << ", " << grid.x[i] << "], radial cell " << j
                        << "; spacing and wall layers were not modified.";
                throw std::runtime_error(message.str());
            }
        }
    return grid;
}

} // namespace

double rpa_wall_chord_error(const std::span<const Point2> wall, Point2 a, Point2 b)
{
    if (a.x > b.x)
        std::swap(a, b);
    if (!(a.x < b.x))
        throw std::invalid_argument("RPA wall chord needs strictly increasing x.");
    double error{std::max(std::abs(a.y - interpolate_rpa_profile(wall, a.x).y),
                          std::abs(b.y - interpolate_rpa_profile(wall, b.x).y))};
    auto point{std::upper_bound(wall.begin(), wall.end(), a.x, [](double x, const Point2 &p) { return x < p.x; })};
    for (; point != wall.end() && point->x < b.x; ++point)
        error = std::max(error, std::abs(point->y - std::lerp(a.y, b.y, (point->x - a.x) / (b.x - a.x))));
    return error;
}

double rpa_angular_quality(const std::span<const Point2> vertices)
{
    if (vertices.size() != 3 && vertices.size() != 4)
        return 0.0;
    double minimum{1.0};
    for (Index i = 0; i < vertices.size(); ++i)
    {
        const auto &a{vertices[i]}, &b{vertices[(i + 1) % vertices.size()]}, &c{vertices[(i + 2) % vertices.size()]};
        const double first{std::hypot(b.x - a.x, b.y - a.y)}, second{std::hypot(c.x - b.x, c.y - b.y)};
        if (!std::isfinite(first) || !std::isfinite(second) || !(first > 0.0) || !(second > 0.0))
            return 0.0;
        const double q{((b.x - a.x) / first) * ((c.y - b.y) / second) - ((b.y - a.y) / first) * ((c.x - b.x) / second)};
        if (!std::isfinite(q) || !(q > 0.0))
            return 0.0;
        minimum = std::min(minimum, q);
    }
    return minimum;
}

RpaMeshDiscretization discretize_rpa(const RpaNozzleDomain &d, const RpaRegionLayout &layout,
                                     const MeshGenerationOptions &generation, const RpaNozzleMeshingOptions &options)
{
    const double base{generation.mesh_size};
    if (!std::isfinite(base) || !(base > 0.0) || !std::isfinite(options.min_cell_quality) ||
        !(options.min_cell_quality > 0.0) || options.min_cell_quality > 1.0)
        throw std::invalid_argument("RPA baseCellSize must be positive; minCellQuality must be in (0,1].");
    const bool even{generation.cell_type == CellType::Quadrilateral && !generation.hybrid};
    const bool finite_plume{options.regions.plume_length_factor.has_value()};
    if ((options.plume_near_length_factor && (!finite_plume || !std::isfinite(*options.plume_near_length_factor) ||
                                              !(*options.plume_near_length_factor > 0.0))) ||
        !std::isfinite(options.plume_axial_spacing_factor) || !(options.plume_axial_spacing_factor > 0.0) ||
        (!options.plume_near_length_factor &&
         (options.plume_axial_spacing_factor != 1.0 || options.plume_radial_growth_rate)) ||
        (options.plume_radial_growth_rate &&
         (!std::isfinite(*options.plume_radial_growth_rate) || *options.plume_radial_growth_rate < 1.0 ||
          *options.plume_radial_growth_rate > options.maximum_growth_rate)))
        throw std::invalid_argument("RPA multiblock PLUME requires finite plumeNearLengthFactor > 0, positive "
                                    "plumeAxialSpacingFactor and plumeRadialGrowthRate in [1,maximumGrowthRate].");
    if ((finite_plume && (!generation.hybrid || !(layout.plume_end_x() < d.downstream_x))) ||
        !std::isfinite(options.engine_axial_spacing_factor) || !(options.engine_axial_spacing_factor > 0) ||
        !std::isfinite(options.maximum_growth_rate) || !(options.maximum_growth_rate > 1.0) ||
        options.plume_refinement_radius_factor.has_value() != options.plume_refinement_cell_size.has_value())
        throw std::invalid_argument("RPA finite plume requires HYBRID and room downstream; invalid spacing options.");
    for (const auto value : {options.farfield_cell_size, options.farfield_interface_cell_size,
                             options.plume_refinement_radius_factor, options.plume_refinement_cell_size})
        if (value && (!std::isfinite(*value) || !(*value > 0)))
            throw std::invalid_argument("RPA exterior size/refinement values must be finite and positive.");
    if (!finite_plume && (options.farfield_cell_size || options.farfield_interface_cell_size ||
                          options.plume_refinement_radius_factor || options.maximum_growth_rate != maximum_growth))
        throw std::invalid_argument("RPA exterior refinement options require a finite plumeLengthFactor.");
    const double run{layout.plume_end_x() - d.inner_wall.back().x};
    const double radius{d.inner_wall.back().y};
    if (options.plume_refinement_radius_factor)
    {
        const double envelope{*options.plume_refinement_radius_factor * radius};
        if (!std::isfinite(envelope) || !(envelope > 0.0) || envelope > layout.plume_outlet_radius())
            throw std::invalid_argument(
                "RPA near-plume refinement radius must be positive and no larger than its outlet radius.");
    }
    const double downstream_factor{run / radius};
    const double requested_factor{layout.plume_outlet_radius() / radius};
    const double outer_factor{d.farfield_radius / radius};
    const double quality{options.min_cell_quality};
    // Invert the outlet-angle sine bound. q == 1 needs an infinite rise for
    // positive run; the factored radicand avoids cancellation near that limit.
    const double minimum_factor{quality == 1.0
                                    ? std::numeric_limits<double>::infinity()
                                    : layout.plume_inlet_radius() / radius +
                                          downstream_factor * (quality / std::sqrt((1.0 - quality) * (1.0 + quality)))};
    if (!finite_plume &&
        (!std::isfinite(minimum_factor) || requested_factor < minimum_factor || !(minimum_factor < outer_factor)))
    {
        std::ostringstream message;
        message << std::setprecision(17)
                << "RPA TRANSITION outlet angle requires plumeOutletRadiusFactor >= minimum=" << minimum_factor
                << "; requested=" << requested_factor
                << ", downstreamFactor=" << (d.downstream_x - d.inner_wall.back().x) / radius
                << ", minCellQuality=" << quality << ", outerRadiusFactor=" << outer_factor << '.';
        if (finite_plume)
            message << " plumeLengthFactor=" << downstream_factor << '.';
        if (!(minimum_factor < outer_factor))
            message << " The minimum plume radius must be strictly below the outer box radius.";
        throw std::runtime_error(message.str());
    }
    RpaMeshDiscretization result;
    result.engine = wall_grid(
        d, d.inner_wall,
        radial_distribution(d.inner_wall.back().y, base, options.inner_wall_layers, true, options.engine_radial_cells),
        0.0, options.min_cell_quality, false, options.engine_axial_spacing_factor, options.preserve_wall_vertices);
    result.band = wall_grid(
        d, d.outer_wall,
        radial_distribution(layout.band_thickness(), base, options.outer_wall_layers, false, options.band_radial_cells),
        layout.band_thickness(), options.min_cell_quality, even, 1.0, options.preserve_wall_vertices);
    // Lip is a physical wall, not a fluid interface. A single interval suffices:
    // it is not on the unstructured transition boundary, even in pure QUAD.
    result.lip_y = {d.inner_wall.back().y, d.outer_wall.back().y};
    const Index rows{result.engine.radial_cells + 1 + result.band.radial_cells};
    std::vector<double> left;
    for (Index j = 0; j <= result.engine.radial_cells; ++j)
        left.push_back(result.engine.at(result.engine.x.size() - 1, j).y);
    left.push_back(result.lip_y.back());
    for (Index j = 1; j <= result.band.radial_cells; ++j)
        left.push_back(result.band.at(result.band.x.size() - 1, j).y);
    const double q{
        std::min(maximum_growth, std::pow(std::max(1.0, layout.plume_outlet_radius() / layout.plume_inlet_radius()),
                                          1.0 / static_cast<double>(rows - 1)))};
    const double terminal_growth{options.plume_radial_growth_rate.value_or(
        options.plume_near_length_factor ? std::min(q, options.maximum_growth_rate) : q)};
    auto right{fractions(rows, terminal_growth)};
    for (auto &v : right)
        v *= layout.plume_outlet_radius();
    const double inlet_step{result.engine.x.back() - result.engine.x[result.engine.x.size() - 2]};
    const double length{run};
    const double axial_target_size{
        finite_plume ? std::sqrt(inlet_step) * std::sqrt(layout.plume_outlet_radius() / static_cast<double>(rows))
                     : std::max(inlet_step, layout.plume_outlet_radius() / static_cast<double>(rows))};
    Index columns{std::max(Index{2}, cell_count(length, axial_target_size))};
    if (even && columns % 2 != 0)
        ++columns;
    const double axial_q{
        std::min(maximum_growth, std::pow(std::max(1.0, layout.plume_outlet_radius() / layout.plume_inlet_radius()),
                                          1.0 / static_cast<double>(columns - 1)))};
    auto axial{fractions(columns, axial_q)};
    result.plume_block_columns = {0, columns};
    if (options.plume_near_length_factor)
    {
        const double first{inlet_step * options.plume_axial_spacing_factor};
        const double cap{std::max(first, options.plume_axial_spacing_factor * layout.plume_outlet_radius() /
                                             static_cast<double>(rows))};
        const double near_length{*options.plume_near_length_factor * radius};
        if (!std::isfinite(first) || !(first > 0) || !std::isfinite(cap) || !std::isfinite(near_length))
            throw std::invalid_argument("RPA multiblock PLUME lengths/sizes are not representable.");
        axial = multiblock_plume_positions(length, near_length, first, cap, options.maximum_growth_rate,
                                           result.plume_block_columns);
        columns = axial.size() - 1;
        for (auto &position : axial)
            position /= length;
    }
    result.plume.radial_cells = rows;
    result.plume.nodes.reserve((columns + 1) * (rows + 1));
    for (double t : axial)
    {
        const double x{std::lerp(d.inner_wall.back().x, layout.plume_end_x(), t)};
        result.plume.x.push_back(x);
        for (Index j = 0; j <= rows; ++j)
            result.plume.nodes.push_back({x, std::lerp(left[j], right[j], t)});
    }
    if (finite_plume)
    {
        result.transition_size = options.farfield_cell_size.value_or(base);
        const double interface_size{options.farfield_interface_cell_size.value_or(base)};
        if (interface_size > result.transition_size)
            throw std::invalid_argument("RPA farfieldInterfaceCellSize must not exceed farfieldCellSize.");
        finite_farfield(result, d, layout, interface_size, result.transition_size, options.maximum_growth_rate);
        return result;
    }
    const double thickness{d.farfield_radius - layout.plume_outlet_radius()};
    const double first{std::max(base, right.back() - right[right.size() - 2])};
    Index outer_rows{2};
    double series{1.0 + maximum_growth}, weight{maximum_growth * maximum_growth};
    while (series + weight <= thickness / first)
    {
        series += weight;
        weight *= maximum_growth;
        ++outer_rows;
        if (outer_rows >= static_cast<Index>(std::numeric_limits<int>::max() - 1))
            throw std::invalid_argument("RPA farfield row overflow.");
    }
    const auto radial{fractions(outer_rows, maximum_growth)};
    const double minimum_height{thickness * (radial[1] - radial[0])};
    const double maximum_height{thickness * (radial.back() - radial[radial.size() - 2])};
    const double axial_target{std::sqrt(minimum_height) * std::sqrt(maximum_height)};
    Index outer_columns{cell_count(d.downstream_x - d.outer_wall.front().x, axial_target)};
    if (even && outer_columns % 2 != 0)
        ++outer_columns;
    result.farfield.radial_cells = outer_rows;
    for (Index i = 0; i <= outer_columns; ++i)
    {
        const double x{std::lerp(d.outer_wall.front().x, d.downstream_x,
                                 static_cast<double>(i) / static_cast<double>(outer_columns))};
        result.farfield.x.push_back(x);
        for (double t : radial)
            result.farfield.nodes.push_back({x, std::lerp(layout.plume_outlet_radius(), d.farfield_radius, t)});
    }
    result.transition_size = thickness / static_cast<double>(outer_rows);
    result.transition_upstream_cells = cell_count(
        layout.plume_outlet_radius() - result.band.at(0, result.band.radial_cells).y, result.transition_size);
    if (even && result.transition_upstream_cells % 2 != 0)
        ++result.transition_upstream_cells;
    return result;
}

void require_rpa_finite_mapped_topology(const RawMeshData &raw, const RpaMeshDiscretization &grids,
                                        const RpaNozzleDomain &domain, const RpaRegionLayout &layout)
{
    Index cell{};
    const double tolerance{64 * std::numeric_limits<double>::epsilon() *
                           std::max(domain.downstream_x, domain.farfield_radius)};
    for (const auto *grid : {&grids.engine, &grids.band, &grids.plume, &grids.farfield, &grids.farfield_downstream})
        for (Index i = 1; i < grid->x.size(); ++i)
            for (Index j = 0; j < grid->radial_cells; ++j, ++cell)
            {
                if (cell >= raw.cell_types.size() || raw.cell_types[cell] != CellType::Quadrilateral ||
                    raw.cell_node_offsets[cell + 1] - raw.cell_node_offsets[cell] != 4)
                    throw std::runtime_error("RPA mapped block lost its logical QUAD topology.");
                const std::array expected{grid->at(i - 1, j), grid->at(i, j), grid->at(i, j + 1),
                                          grid->at(i - 1, j + 1)};
                Point2 center{};
                for (Index k = 0; k < expected.size(); ++k)
                {
                    const auto &point{raw.nodes[raw.cell_nodes[raw.cell_node_offsets[cell] + k]]};
                    if (std::hypot(point.x - expected.at(k).x, point.y - expected.at(k).y) > tolerance)
                        throw std::runtime_error("RPA mapped block no longer matches its logical grid coordinates.");
                    center.x += point.x / 4;
                    center.y += point.y / 4;
                }
                if (layout.region(center) == RpaRegion::Transition)
                    throw std::runtime_error("RPA QUAD incorrectly lies in TRANSITION.");
            }
    if (cell == raw.cell_types.size())
        throw std::runtime_error("RPA finite plume requires a nonempty TRI TRANSITION.");
    for (; cell < raw.cell_types.size(); ++cell)
    {
        if (raw.cell_types[cell] != CellType::Triangle)
            throw std::runtime_error("RPA TRANSITION must contain TRI only.");
        Point2 center{};
        for (Index k = raw.cell_node_offsets[cell]; k < raw.cell_node_offsets[cell + 1]; ++k)
        {
            const auto &point{raw.nodes[raw.cell_nodes[k]]};
            center.x += point.x / 3;
            center.y += point.y / 3;
        }
        if (layout.region(center) != RpaRegion::Transition)
        {
            std::ostringstream message;
            message << "RPA TRI " << cell << " outside TRANSITION at (" << center.x << ',' << center.y << ").";
            throw std::runtime_error(message.str());
        }
    }
}

void require_rpa_mesh_quality(const RawMeshData &mesh, const RpaNozzleMeshingOptions &options)
{
    for (Index cell = 0; cell < mesh.cell_types.size(); ++cell)
    {
        std::array<Point2, 4> vertices{};
        const Index begin{mesh.cell_node_offsets[cell]}, count{mesh.cell_node_offsets[cell + 1] - begin};
        for (Index i = 0; i < count; ++i)
            vertices.at(i) = mesh.nodes[mesh.cell_nodes[begin + i]];
        const double quality{rpa_angular_quality(std::span{vertices}.first(count))};
        if (!(quality >= options.min_cell_quality))
        {
            std::ostringstream message;
            message << std::setprecision(17) << "RPA generated cell " << cell << " violates minCellQuality: " << quality
                    << " < " << options.min_cell_quality << "; vertices:";
            for (const auto &p : std::span{vertices}.first(count))
                message << " (" << p.x << ',' << p.y << ')';
            message << ". No threshold was relaxed.";
            throw std::runtime_error(message.str());
        }
    }
}
} // namespace cfd::detail
