#include "RpaMeshDiscretization.hpp"
#include "cfd/input/OpenFOAMCaseReader.hpp"
#include "cfd/io/VtkWriter.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/mesh/MeshStatistics.hpp"
#include "cfd/meshing/RpaMeshDiagnostics.hpp"
#include "support/TestUtils.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <numbers>
#include <numeric>
#include <set>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace
{
using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws;

class TemporaryDirectory
{
  public:
    TemporaryDirectory()
    {
        std::string pattern{(std::filesystem::temp_directory_path() / "cfd-rpa-tests-XXXXXX").string()};
        const char *created{::mkdtemp(pattern.data())};
        if (!created)
        {
            throw std::runtime_error("Cannot create RPA test directory.");
        }
        path_ = created;
    }

    ~TemporaryDirectory()
    {
        std::error_code error;
        std::filesystem::remove_all(path_, error);
    }

    TemporaryDirectory(const TemporaryDirectory &) = delete;
    TemporaryDirectory &operator=(const TemporaryDirectory &) = delete;
    TemporaryDirectory(TemporaryDirectory &&) = delete;
    TemporaryDirectory &operator=(TemporaryDirectory &&) = delete;

    const std::filesystem::path &path() const noexcept
    {
        return path_;
    }

  private:
    std::filesystem::path path_;
};

double radius_at(const std::vector<cfd::Point2> &wall, const double x)
{
    const auto right{
        std::upper_bound(wall.begin(), wall.end(), x, [](double value, const cfd::Point2 &p) { return value < p.x; })};
    if (right == wall.begin())
    {
        return wall.front().y;
    }
    if (right == wall.end())
    {
        return wall.back().y;
    }
    const auto &left{*(right - 1)};
    return left.y + (right->y - left.y) * (x - left.x) / (right->x - left.x);
}

void validate_mesh(const cfd::Mesh &mesh, const cfd::RpaNozzleDomain &domain, const cfd::CellType type,
                   const bool hybrid = false)
{
    require(mesh.boundary_groups().size() == 8, "Missing RPA physical boundaries.");
    constexpr std::array<std::string_view, 8> names{"inlet",           "axis",
                                                    "inner_wall",      "outer_wall",
                                                    "nozzle_lip",      "farfield_upstream",
                                                    "farfield_radial", "farfield_downstream"};
    std::array<cfd::Index, 8> counts{};
    const auto exit{domain.inner_wall.back()};
    cfd::Vector2 boundary_sum{};
    for (cfd::Index face = 0; face < mesh.face_count(); ++face)
    {
        if (!mesh.face_adjacencies()[face].is_boundary())
        {
            continue;
        }
        const auto id{mesh.face_boundary_ids()[face]};
        require(id < names.size() && mesh.boundary_groups()[id].name == names.at(id),
                "Incorrect RPA boundary classification.");
        ++counts.at(id);
        const auto centre{mesh.face_centers()[face]};
        boundary_sum.x += mesh.face_area_vectors()[face].x;
        boundary_sum.y += mesh.face_area_vectors()[face].y;
        if (id == 0)
        {
            require_near(centre.x, domain.inner_wall.front().x, 1e-12, "Inlet misplaced.");
        }
        if (id == 1)
        {
            require_near(centre.y, 0.0, 1e-12, "Axis misplaced.");
        }
        if (id == 4)
        {
            require_near(centre.x, exit.x, 1e-12, "Lip misplaced.");
        }
        if (id == 5)
        {
            require_near(centre.x, domain.outer_wall.front().x, 1e-12, "Upstream farfield misplaced.");
        }
        if (id == 6)
        {
            require_near(centre.y, domain.farfield_radius, 1e-12, "Radial farfield misplaced.");
        }
        if (id == 7)
        {
            require_near(centre.x, domain.downstream_x, 1e-12, "Downstream farfield misplaced.");
        }
        if (std::abs(centre.x - exit.x) < 1e-12 && centre.y > 0 && centre.y < exit.y)
        {
            require(false, "Geometric nozzle exit was incorrectly made a CFD boundary.");
        }
    }
    for (const auto count : counts)
    {
        require(count > 0, "Empty RPA boundary group.");
    }
    require_near(boundary_sum.x, 0.0, 1e-11, "Fluid boundary is not closed in x.");
    require_near(boundary_sum.y, 0.0, 1e-11, "Fluid boundary is not closed in r.");
    for (cfd::Index cell = 0; cell < mesh.cell_count(); ++cell)
    {
        require((hybrid || mesh.cell_types()[cell] == type) && mesh.cell_areas()[cell] > 0,
                "Unexpected topology or invalid area.");
        const auto point{mesh.cell_centers()[cell]};
        require(point.y >= 0 && point.y <= domain.farfield_radius, "Cell outside positive-radius domain.");
        if (point.x >= domain.outer_wall.front().x && point.x <= exit.x)
        {
            const double inner{radius_at(domain.inner_wall, point.x)};
            const double outer{radius_at(domain.outer_wall, point.x)};
            require(point.y <= inner || point.y >= outer,
                    "Solid wall contains fluid-cell centre " + std::to_string(cell) + " at (" +
                        std::to_string(point.x) + "," + std::to_string(point.y) + "), inner=" + std::to_string(inner) +
                        ", outer=" + std::to_string(outer));
        }
    }
    std::vector<bool> visited(mesh.cell_count(), false);
    std::vector<cfd::Index> queue{0};
    visited[0] = true;
    for (cfd::Index head = 0; head < queue.size(); ++head)
    {
        const auto cell{queue[head]};
        for (cfd::Index offset = mesh.cell_node_offsets()[cell]; offset < mesh.cell_node_offsets()[cell + 1]; ++offset)
        {
            const auto &adjacency{mesh.face_adjacencies()[mesh.cell_faces()[offset]]};
            if (adjacency.is_boundary())
            {
                continue;
            }
            const auto next{adjacency.owner == cell ? adjacency.neighbor : adjacency.owner};
            if (!visited[next])
            {
                visited[next] = true;
                queue.push_back(next);
            }
        }
    }
    require(queue.size() == mesh.cell_count(), "Interior and exterior fluid are disconnected.");
    // Independent integral of the meridian polygon: upstream inner fluid,
    // then farfield rectangle minus the radial fictitious wall strip.
    const double exterior_start{domain.outer_wall.front().x};
    double expected_area{domain.farfield_radius * (domain.downstream_x - exterior_start) -
                         (domain.outer_wall.back().y - exit.y) * (exit.x - exterior_start)};
    for (cfd::Index index = 1; index < domain.inner_wall.size(); ++index)
    {
        const auto &left{domain.inner_wall[index - 1]};
        if (left.x >= exterior_start)
        {
            break;
        }
        const double right_x{std::min(domain.inner_wall[index].x, exterior_start)};
        expected_area += 0.5 * (left.y + radius_at(domain.inner_wall, right_x)) * (right_x - left.x);
    }
    const auto statistics{cfd::compute_mesh_statistics(mesh)};
    // First-order boundary chords add a discretization error independent of
    // the CAD flattening tolerance; do not demand exact curved-domain area.
    require_near(statistics.total_cell_area, expected_area, 1e-4 * std::max(1.0, expected_area),
                 "Fluid area disagrees with the meridian domain.");
}

cfd::RpaNozzleGeometry geometry(const std::filesystem::path &file = CFD_RPA_REFERENCE_DXF)
{
    return {.file = file, .downstream_factor = 20.0, .outer_radius_factor = 10.0};
}

cfd::input::MeshInput dictionary(const std::filesystem::path &directory, const std::string &body,
                                 const std::string &mesh = "cellType hybrid; baseCellSize 0.04;",
                                 const std::string &regions = "")
{
    std::filesystem::create_directories(directory / "system");
    const auto path{directory / "system/meshDict"};
    {
        std::ofstream out{path};
        out << "FoamFile {format ascii; class dictionary; object meshDict;}\ngeometry {type rpaNozzle; " << body
            << "}\n"
            << regions << "\nmesh {" << mesh << "}\n";
        require(static_cast<bool>(out), "Could not write meshDict.");
    }
    return cfd::input::read_mesh_dict(path);
}

void test_reader_and_migration()
{
    const TemporaryDirectory dir;
    const auto parsed{dictionary(dir.path(),
                                 "box { outerRadiusFactor 10; downstreamFactor 20; upstreamFactor auto; } engine { "
                                 "wallThicknessFactor 0.02; file \"constant/geometry/engine.dxf\"; }")};
    const auto &g{std::get<cfd::RpaNozzleGeometry>(parsed.geometry)};
    require(g.file == dir.path() / "constant/geometry/engine.dxf" && !g.upstream_factor && g.downstream_factor == 20 &&
                g.outer_radius_factor == 10,
            "Nested geometry parsing failed.");
    require(parsed.rpa_nozzle_meshing.regions.plume_outlet_radius_factor == 16 &&
                parsed.rpa_nozzle_meshing.min_cell_quality == .05,
            "RPA defaults changed.");
    const auto defaults{dictionary(dir.path(), "engine {file test.dxf;}")};
    const auto &dg{std::get<cfd::RpaNozzleGeometry>(defaults.geometry)};
    require(dg.downstream_factor == 80 && dg.outer_radius_factor == 65, "Extended box defaults changed.");
    const auto custom{dictionary(dir.path(), "engine {file test.dxf;}",
                                 "minCellQuality 0.1; baseCellSize 0.02; cellType quadrilateral; outerWallLayers true; "
                                 "outerWallFirstCellHeight 0.001; outerWallLayerCount 4; outerWallGrowthRate 1.1;",
                                 "regions {plumeOutletRadiusFactor 4; bandThicknessFactor 0.25;}")};
    require(custom.rpa_nozzle_meshing.regions.band_thickness_factor == .25 &&
                custom.rpa_nozzle_meshing.outer_wall_layers.active(),
            "Regions/layers overrides failed.");
    for (const auto *old :
         {"size", "automaticMeshing", "throatRefinementFactor", "convergentRefinementFactor",
          "divergentRefinementFactor", "exitRefinementFactor", "nearPlumeLength", "plumeCorridorLength",
          "plumeOpeningSlope", "chamberRefinementFactor", "lipRefinementFactor", "plumeRefinementFactor",
          "farPlumeRefinementFactor", "farfieldRefinementFactor", "resolutionLevel", "resolutionPolicy"})
        require_throws<std::runtime_error>(
            [&] {
                static_cast<void>(dictionary(dir.path(), "engine {file a.dxf;}",
                                             "cellType hybrid; baseCellSize 0.04; " + std::string{old} + " 0.5;"));
            },
            "Removed RPA option silently accepted.");
    for (const auto *bad :
         {"file a.dxf;", "engine {file a.dxf; file b.dxf;}", "engine {file a.dxf; wallThicknessFactor 0;}",
          "engine {file a.dxf;} box {downstreamFactor -1;}", "engine {file a.dxf;} engine {file b.dxf;}"})
        require_throws<std::runtime_error>([&] { static_cast<void>(dictionary(dir.path(), bad)); },
                                           "Invalid geometry accepted.");
    for (const auto *bad :
         {"cellType hybrid; baseCellSize 0;", "cellType hybrid; baseCellSize 0.04; minCellQuality 1.1;",
          "cellType hybrid; baseCellSize 0.04; baseCellSize 0.02;",
          "cellType hybrid; baseCellSize 0.04; innerWallFirstCellHeight .001;"})
        require_throws<std::runtime_error>(
            [&] { static_cast<void>(dictionary(dir.path(), "engine {file a.dxf;}", bad)); },
            "Invalid mesh options accepted.");
}

void test_geometry_and_layout()
{
    auto g{geometry()};
    g.outer_radius_factor = 65;
    const auto d{cfd::read_rpa_nozzle_domain(g)};
    const auto exit{d.inner_wall.back()};
    require_near(d.downstream_x - exit.x, 20 * exit.y, 1e-12, "Wrong downstream scaling.");
    require_near(d.farfield_radius, 65 * exit.y, 1e-12, "Wrong radial box scaling.");
    require_near(d.outer_wall.front().x, d.inner_wall.front().x, 1e-12,
                 "Auto upstream extent did not use the complete engine.");
    for (const auto &p : d.outer_wall)
        require_near(p.y - radius_at(d.inner_wall, p.x), .02 * exit.y, 1e-12, "Shell offset changed.");
    auto alternate{g};
    alternate.wall_thickness_factor = .04;
    alternate.upstream_factor = 2;
    const auto other{cfd::read_rpa_nozzle_domain(alternate)};
    require(other.inner_wall.size() == d.inner_wall.size() && other.curve_junctions == d.curve_junctions,
            "CAD profile/junctions changed.");
    for (cfd::Index i = 0; i < d.inner_wall.size(); ++i)
        require(other.inner_wall[i].x == d.inner_wall[i].x && other.inner_wall[i].y == d.inner_wall[i].y,
                "DXF coordinates changed.");
    const cfd::RpaRegionLayout l{d};
    require_near(l.band_thickness(), .5 * exit.y, 1e-12, "Band factor failed.");
    require_near(l.plume_inlet_radius(), 1.52 * exit.y, 1e-12, "Plume inlet radius failed.");
    require_near(l.plume_outlet_radius(), 16 * exit.y, 1e-12, "Default plume outlet factor is not 16.");
    require_near(l.plume_radius(std::midpoint(exit.x, d.downstream_x)),
                 .5 * (l.plume_inlet_radius() + l.plume_outlet_radius()), 1e-12, "Nonlinear plume envelope.");
    for (double value : {0., -1., std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
    {
        auto bad{g};
        bad.downstream_factor = value;
        require_throws<std::invalid_argument>([&] { static_cast<void>(cfd::read_rpa_nozzle_domain(bad)); },
                                              "Invalid box factor accepted.");
    }
    auto bad{g};
    bad.upstream_factor = 1e6;
    require_throws<std::invalid_argument>([&] { static_cast<void>(cfd::read_rpa_nozzle_domain(bad)); },
                                          "Unavailable upstream contour accepted.");
    bad = g;
    bad.outer_radius_factor = .5;
    require_throws<std::invalid_argument>([&] { static_cast<void>(cfd::read_rpa_nozzle_domain(bad)); },
                                          "Box intersecting shell accepted.");
    require_throws<std::invalid_argument>(
        [&] { cfd::RpaRegionLayout invalid{d, {.band_thickness_factor = .5, .plume_outlet_radius_factor = 1.}}; },
        "Invalid plume envelope accepted.");
    require_throws<std::runtime_error>(
        [&] { static_cast<void>(l.region({exit.x - .01, radius_at(d.inner_wall, exit.x - .01) + .01 * exit.y})); },
        "Solid point accepted.");
}

void test_discretization_and_scale()
{
    cfd::RpaNozzleDomain d{.inner_wall = {{0, 1}, {1, 1}, {2, .25}, {3, 1}},
                           .outer_wall = {{0, 1.02}, {1, 1.02}, {2, .27}, {3, 1.02}},
                           .curve_junctions = {0, 1, 2, 3},
                           .throat = {2, .25},
                           .downstream_x = 23,
                           .farfield_radius = 10,
                           .geometric_tolerance = 1e-5};
    const cfd::RpaNozzleMeshingOptions o{.regions = {.plume_outlet_radius_factor = 6}};
    const cfd::RpaRegionLayout l{d, o.regions};
    const auto m{cfd::detail::discretize_rpa(
        d, l, {.mesh_size = .1, .cell_type = cfd::CellType::Quadrilateral, .hybrid = true}, o)};
    require(m.engine.radial_cells == 10 && m.band.radial_cells == 5, "Constant radial count formula changed.");
    require(m.plume.radial_cells == m.engine.radial_cells + 1 + m.band.radial_cells, "Plume is not ENGINE+LIP+BAND.");
    for (const auto *grid : {&m.engine, &m.band, &m.plume, &m.farfield})
        for (cfd::Index i = 0; i < grid->x.size(); ++i)
            for (cfd::Index j = 0; j <= grid->radial_cells; ++j)
                require(grid->at(i, j).x == grid->x[i], "Core column is not vertical.");
    for (double x : d.curve_junctions)
        require(std::binary_search(m.engine.x.begin(), m.engine.x.end(), x), "CAD junction lost.");
    double chamber{}, throat{1e100};
    for (cfd::Index i = 1; i < m.engine.x.size(); ++i)
    {
        double x = std::midpoint(m.engine.x[i - 1], m.engine.x[i]);
        double dx = m.engine.x[i] - m.engine.x[i - 1];
        if (x < .8)
            chamber = std::max(chamber, dx);
        if (x > 1.8 && x < 2.2)
            throat = std::min(throat, dx);
    }
    require(throat < .6 * chamber, "Axial width did not follow the smaller throat core height.");
    auto scaled{d};
    for (auto &p : scaled.inner_wall)
    {
        p.x *= 10;
        p.y *= 10;
    }
    for (auto &p : scaled.outer_wall)
    {
        p.x *= 10;
        p.y *= 10;
    }
    for (auto &x : scaled.curve_junctions)
        x *= 10;
    scaled.throat = {20, 2.5};
    scaled.downstream_x *= 10;
    scaled.farfield_radius *= 10;
    scaled.geometric_tolerance *= 10;
    const auto sm{
        cfd::detail::discretize_rpa(scaled, cfd::RpaRegionLayout{scaled, o.regions},
                                    {.mesh_size = 1., .cell_type = cfd::CellType::Quadrilateral, .hybrid = true}, o)};
    require(sm.engine.nodes.size() == m.engine.nodes.size(), "Scaling changed dimensionless resolution.");
    for (cfd::Index i = 0; i < m.engine.nodes.size(); ++i)
    {
        require_near(sm.engine.nodes[i].x, 10 * m.engine.nodes[i].x, 1e-10, "Scaled x mismatch.");
        require_near(sm.engine.nodes[i].y, 10 * m.engine.nodes[i].y, 1e-10, "Scaled radius mismatch.");
    }
    const std::array<cfd::Point2, 4> anisotropic{{{0, 0}, {1000, 0}, {1000, 1e-6}, {0, 1e-6}}};
    require_near(cfd::detail::rpa_angular_quality(anisotropic), 1., 1e-14,
                 "Angular quality penalizes orthogonal anisotropy.");
    const std::array<cfd::Point2, 4> folded{{{0, 0}, {1, 0}, {.2, -.1}, {0, 1}}};
    require(cfd::detail::rpa_angular_quality(folded) == 0, "Folded QUAD accepted.");
    const double angle{std::numbers::pi / 180.0};
    const std::array<cfd::Point2, 3> acute{{{0, 0}, {1, 0}, {1, std::tan(angle)}}};
    require_near(cfd::detail::rpa_angular_quality(acute), std::sin(angle), 1e-14,
                 "A 1/89/90-degree triangle was reported as excellent.");
    cfd::RawMeshData acute_raw;
    acute_raw.nodes.assign(acute.begin(), acute.end());
    acute_raw.cell_types = {cfd::CellType::Triangle};
    acute_raw.cell_node_offsets = {0, 3};
    acute_raw.cell_nodes = {0, 1, 2};
    require_throws<std::runtime_error>(
        [&] { cfd::detail::require_rpa_mesh_quality(acute_raw, {.min_cell_quality = .1}); },
        "A very acute but positively oriented TRI passed the angular threshold.");
    cfd::detail::require_rpa_mesh_quality(acute_raw, {.min_cell_quality = .01});
    auto invalid{o};
    invalid.min_cell_quality = 1.1;
    require_throws<std::invalid_argument>(
        [&] { static_cast<void>(cfd::detail::discretize_rpa(d, l, {.mesh_size = .1}, invalid)); },
        "Invalid quality option accepted.");
    invalid.min_cell_quality = .1;
    invalid.regions.plume_outlet_radius_factor = 3;
    require_throws<std::runtime_error>(
        [&] {
            static_cast<void>(cfd::detail::discretize_rpa(d, cfd::RpaRegionLayout{d, invalid.regions},
                                                          {.mesh_size = .1, .cell_type = cfd::CellType::Quadrilateral},
                                                          invalid));
        },
        "Impossible QUAD outlet-angle threshold accepted.");
}

void test_mean_height_columns_and_chord_diagnostics()
{
    const cfd::RpaNozzleDomain flat{.inner_wall = {{0, 1}, {1, 1}, {2, .25}, {3, 1}},
                                    .outer_wall = {{0, 1.02}, {1, 1.02}, {2, .27}, {3, 1.02}},
                                    .curve_junctions = {0, 1, 2, 3},
                                    .throat = {2, .25},
                                    .downstream_x = 23,
                                    .farfield_radius = 10,
                                    .geometric_tolerance = 1e-6};
    const cfd::MeshGenerationOptions gen{.mesh_size = .1, .cell_type = cfd::CellType::Quadrilateral, .hybrid = true};
    const cfd::RpaNozzleMeshingOptions options{.regions = {.plume_outlet_radius_factor = 6}};
    const auto baseline{cfd::detail::discretize_rpa(flat, cfd::RpaRegionLayout{flat, options.regions}, gen, options)};
    auto curved{flat};
    curved.inner_wall.insert(curved.inner_wall.begin() + 1,
                             {{.5, 1}, {.5005, .9995}, {.501, .998}, {.5015, .9995}, {.502, 1}});
    curved.outer_wall.clear();
    for (const auto &point : curved.inner_wall)
        curved.outer_wall.push_back({point.x, point.y + .02});
    const auto diagnostic{
        cfd::detail::discretize_rpa(curved, cfd::RpaRegionLayout{curved, options.regions}, gen, options)};
    require(diagnostic.band.x.size() == baseline.band.x.size(),
            "A local contour error still triggers extra BAND columns.");
    for (cfd::Index i = 0; i < baseline.band.x.size(); ++i)
        require_near(diagnostic.band.x[i], baseline.band.x[i], 1e-12,
                     "Constant-height BAND no longer follows the arithmetic mean.");
    for (const auto *grid : {&diagnostic.engine, &diagnostic.band})
    {
        const auto &wall{grid == &diagnostic.engine ? curved.inner_wall : curved.outer_wall};
        double maximum{};
        for (cfd::Index i = 1; i < grid->x.size(); ++i)
            maximum = std::max(
                maximum, cfd::detail::rpa_wall_chord_error(wall, {grid->x[i - 1], radius_at(wall, grid->x[i - 1])},
                                                           {grid->x[i], radius_at(wall, grid->x[i])}));
        require(maximum > curved.geometric_tolerance, "The deliberate geometric-tolerance incompatibility was hidden.");
    }
    auto tighter{curved};
    tighter.geometric_tolerance = 1e-12;
    const auto same{cfd::detail::discretize_rpa(tighter, cfd::RpaRegionLayout{tighter, options.regions}, gen, options)};
    require(same.engine.x == diagnostic.engine.x && same.band.x == diagnostic.band.x,
            "Geometric tolerance still determines column widths.");
    for (const auto *grid : {&baseline.engine, &baseline.band})
        for (cfd::Index i = 1; i < grid->x.size(); ++i)
        {
            const auto &wall{grid == &baseline.engine ? flat.inner_wall : flat.outer_wall};
            require(cfd::detail::rpa_wall_chord_error(wall, {grid->x[i - 1], radius_at(wall, grid->x[i - 1])},
                                                      {grid->x[i], radius_at(wall, grid->x[i])}) <=
                        flat.geometric_tolerance,
                    "Straight CAD segment exceeds the chord tolerance.");
        }

    const std::array<cfd::Point2, 3> known_curve{{{0, 1}, {.5, 1.1}, {1, 1}}};
    require_near(cfd::detail::rpa_wall_chord_error(known_curve, {0, 1}, {1, 1}), .1, 1e-14,
                 "Chord-error diagnostic misses a known profile deviation.");
    require_near(cfd::detail::rpa_wall_chord_error(known_curve, {1, 1}, {0, 1}), .1, 1e-14,
                 "Chord-error diagnostic depends on face orientation.");
    const std::array<cfd::Point2, 2> line{{{0, 1}, {1, 1}}};
    require_near(cfd::detail::rpa_wall_chord_error(line, {0, 1.01}, {1, 1.01}), .01, 1e-14,
                 "Chord-error diagnostic ignores displaced endpoints.");

    auto straight{flat};
    straight.inner_wall = {{0, 1}, {.5, 1}, {1.001, 1}};
    straight.outer_wall = {{0, 1.02}, {1.001, 1.02}};
    straight.curve_junctions = {0, 1.001};
    straight.throat = {0, 1};
    const auto redistributed{
        cfd::detail::discretize_rpa(straight, cfd::RpaRegionLayout{straight, options.regions}, gen, options)};
    for (const auto *grid : {&redistributed.engine, &redistributed.band})
    {
        require(grid->x.size() == 11, "Arithmetic mean did not select the nearest interval count.");
        for (cfd::Index i = 1; i < grid->x.size(); ++i)
            require_near(grid->x[i] - grid->x[i - 1], .1001, 1e-12, "Endpoint redistribution left a tiny last cell.");
    }
    auto layers{options};
    layers.inner_wall_layers = {.enabled = true, .first_cell_height = .0025, .layer_count = 3, .growth_rate = 1.1};
    layers.outer_wall_layers = {.enabled = true, .first_cell_height = .003, .layer_count = 3, .growth_rate = 1.1};
    const auto layered{
        cfd::detail::discretize_rpa(straight, cfd::RpaRegionLayout{straight, layers.regions}, gen, layers)};
    require(layered.engine.x == redistributed.engine.x && layered.band.x == redistributed.band.x,
            "Physical-layer heights replaced the full arithmetic radial mean.");

    const auto d{cfd::read_rpa_nozzle_domain(geometry())};
    const cfd::RpaNozzleMeshingOptions old_plume{.regions = {.plume_outlet_radius_factor = 3}};
    const auto plan{cfd::detail::discretize_rpa(
        d, cfd::RpaRegionLayout{d, old_plume.regions},
        {.mesh_size = .01, .cell_type = cfd::CellType::Quadrilateral, .hybrid = true}, old_plume)};
    const auto &f{plan.farfield};
    require(f.radial_cells == 12, "Existing FARFIELD radial distribution changed.");
    const double dr_min{f.at(0, 1).y - f.at(0, 0).y};
    const double dr_max{f.at(0, f.radial_cells).y - f.at(0, f.radial_cells - 1).y};
    const auto expected{static_cast<cfd::Index>(std::ceil((f.x.back() - f.x.front()) / std::sqrt(dr_min * dr_max)))};
    require(f.x.size() - 1 == expected, "FARFIELD axial count still follows the PLUME outlet step.");
    for (cfd::Index i = 2; i < f.x.size(); ++i)
        require_near(f.x[i] - f.x[i - 1], f.x[1] - f.x[0], 1e-12, "FARFIELD x spacing is not uniform.");
}

void test_transition_minimum_outlet_factor()
{
    cfd::RpaNozzleDomain domain{.inner_wall = {{0, 1}, {1, 1}, {2, 1}},
                                .outer_wall = {{0, 1.02}, {2, 1.02}},
                                .curve_junctions = {0, 2},
                                .throat = {0, 1},
                                .downstream_x = 22,
                                .farfield_radius = 10};
    const cfd::MeshGenerationOptions generation{.mesh_size = .1, .cell_type = cfd::CellType::Quadrilateral};
    cfd::RpaNozzleMeshingOptions options{.regions = {.plume_outlet_radius_factor = 3}};
    static_cast<void>(
        cfd::detail::discretize_rpa(domain, cfd::RpaRegionLayout{domain, options.regions}, generation, options));
    require(options.regions.plume_outlet_radius_factor == 3, "Admissible outlet factor was changed.");

    const auto reject = [&](double downstream, double quality, double outer, double expected_minimum) {
        domain.downstream_x = 2 + downstream;
        domain.farfield_radius = outer;
        options.min_cell_quality = quality;
        try
        {
            static_cast<void>(cfd::detail::discretize_rpa(domain, cfd::RpaRegionLayout{domain, options.regions},
                                                          generation, options));
        }
        catch (const std::runtime_error &error)
        {
            const std::string message{error.what()};
            for (const auto *key : {"requested=3", "minimum=", "downstreamFactor=", "minCellQuality="})
                require(message.find(key) != std::string::npos, "Missing outlet-factor diagnostic value.");
            const double minimum{std::stod(message.substr(message.find("minimum=") + 8))};
            if (std::isinf(expected_minimum))
                require(minimum == expected_minimum, "q=1 must require an infinite outlet factor.");
            else
                require_near(minimum, expected_minimum, expected_minimum * 1e-14, "Incorrect minimum outlet factor.");
            if (expected_minimum >= outer)
                require(message.find("outer box radius") != std::string::npos,
                        "Unavailable outer radius was not diagnosed.");
            require(options.regions.plume_outlet_radius_factor == 3, "Rejected outlet factor was changed.");
            return;
        }
        require(false, "Incompatible TRANSITION outlet factor was accepted.");
    };
    reject(80, .05, 65, 1.52 + 80 * .05 / std::sqrt(1 - .05 * .05));
    reject(20, .1, 10, 1.52 + 20 * .1 / std::sqrt(1 - .1 * .1));
    reject(80, .05, 5, 1.52 + 80 * .05 / std::sqrt(1 - .05 * .05));
    reject(20, 1, 10, std::numeric_limits<double>::infinity());
    const double near_one{std::nextafter(1., 0.)};
    reject(20, near_one, 10, 1.52 + 20 * near_one / std::sqrt((1 - near_one) * (1 + near_one)));
    // A positive, subnormal q must not create a zero denominator or NaN bound.
    options.min_cell_quality = std::numeric_limits<double>::denorm_min();
    static_cast<void>(
        cfd::detail::discretize_rpa(domain, cfd::RpaRegionLayout{domain, options.regions}, generation, options));
}

void validate_generated_mesh(const cfd::Mesh &m, const cfd::RpaNozzleDomain &d, const cfd::MeshGenerationOptions &gen,
                             const cfd::RpaNozzleMeshingOptions &o)
{
    validate_mesh(m, d, gen.cell_type, gen.hybrid);
    const auto plan{cfd::detail::discretize_rpa(d, cfd::RpaRegionLayout{d, o.regions}, gen, o)};
    std::vector<cfd::Point2> sorted_nodes(m.nodes().begin(), m.nodes().end());
    std::ranges::sort(sorted_nodes, {}, &cfd::Point2::x);
    // A local TRANSITION repair must not move any prescribed grid node.
    // Gmsh can reproject curved CAD edges within the import tolerance.
    for (const auto *grid : {&plan.engine, &plan.band, &plan.plume, &plan.farfield})
        for (const auto &point : grid->nodes)
        {
            const double tolerance{grid == &plan.engine || grid == &plan.band ? d.geometric_tolerance : 1e-10};
            const auto begin{std::ranges::lower_bound(sorted_nodes, point.x - tolerance, {}, &cfd::Point2::x)};
            const auto end{std::ranges::upper_bound(sorted_nodes, point.x + tolerance, {}, &cfd::Point2::x)};
            require(std::any_of(begin, end, [&](const auto &p) { return std::abs(p.y - point.y) <= tolerance; }),
                    "TRANSITION repair moved a prescribed structured-grid node.");
        }
    std::vector<double> regions(m.cell_count());
    std::ostringstream report;
    cfd::write_rpa_mesh_diagnostics(report, m, d, gen, o, regions);
    require(report.str().find("maximum vertical chord error") != std::string::npos &&
                report.str().find("EXCEEDED") != std::string::npos &&
                report.str().find("Topological validity does not certify") != std::string::npos,
            "Diagnostic mesh was silently certified against the geometric tolerance.");
    std::array<std::set<cfd::Index>, 5> nodes;
    std::array<cfd::Index, 5> counts{};
    for (cfd::Index cell = 0; cell < m.cell_count(); ++cell)
    {
        auto r = static_cast<cfd::Index>(regions[cell]) - 1;
        ++counts.at(r);
        if (gen.hybrid)
            require((m.cell_types()[cell] == cfd::CellType::Triangle) == (r == 3),
                    "HYBRID element escaped its region.");
        for (auto j = m.cell_node_offsets()[cell]; j < m.cell_node_offsets()[cell + 1]; ++j)
            nodes.at(r).insert(m.cell_nodes()[j]);
    }
    for (auto count : counts)
        require(count > 0, "Empty fluid region.");
    std::set<cfd::Index> left, joined;
    std::vector<double> right;
    for (auto i : nodes.at(2))
    {
        const auto &p{m.nodes()[i]};
        if (std::abs(p.x - d.inner_wall.back().x) < 1e-10)
            left.insert(i);
        if (std::abs(p.x - d.downstream_x) < 1e-10)
            right.push_back(p.y);
    }
    for (auto r : {0U, 1U})
        for (auto i : nodes.at(r))
            if (std::abs(m.nodes()[i].x - d.inner_wall.back().x) < 1e-10)
                joined.insert(i);
    for (cfd::Index face = 0; face < m.face_count(); ++face)
        if (m.face_boundary_ids()[face] == 4)
            for (auto i : m.faces()[face].node_ids)
                joined.insert(i);
    require(left == joined && left.size() == right.size(), "Plume did not reuse the ENGINE/LIP/BAND node union.");
    std::vector<double> ys;
    ys.reserve(left.size());
    for (auto i : left)
        ys.push_back(m.nodes()[i].y);
    std::ranges::sort(ys);
    std::ranges::sort(right);
    for (auto i : nodes.at(2))
    {
        const auto &p{m.nodes()[i]};
        double t = (p.x - d.inner_wall.back().x) / (d.downstream_x - d.inner_wall.back().x);
        bool on_row{};
        for (cfd::Index j = 0; j < ys.size(); ++j)
            on_row = on_row || std::abs(p.y - std::lerp(ys[j], right[j], t)) < 1e-10;
        require(on_row, "Plume row is not straight.");
    }
    for (cfd::Index j = 2; j < right.size(); ++j)
        require(right[j] - right[j - 1] > right[j - 1] - right[j - 2], "Outlet law is not growing.");

    std::vector<bool> used(m.node_count());
    for (const auto node : m.cell_nodes())
        used[node] = true;
    std::set<std::pair<double, double>> coordinates;
    for (cfd::Index node = 0; node < m.node_count(); ++node)
        if (used[node])
            require(coordinates.emplace(m.nodes()[node].x, m.nodes()[node].y).second,
                    "Distinct used nodes share the same coordinates.");
    const std::array<const cfd::detail::RpaStructuredGrid *, 4> grids{&plan.engine, &plan.band, &plan.plume,
                                                                      &plan.farfield};
    const std::array<cfd::Index, 4> region_indices{0, 1, 2, 4};
    for (cfd::Index i = 0; i < grids.size(); ++i)
        require(counts.at(region_indices.at(i)) ==
                    (grids.at(i)->x.size() - 1) * grids.at(i)->radial_cells *
                        (gen.cell_type == cfd::CellType::Triangle && !gen.hybrid ? 2 : 1),
                "Structured region cell/radial cardinality mismatch.");
    for (cfd::Index cell = 0; cell < m.cell_count(); ++cell)
    {
        std::array<cfd::Point2, 4> corners{};
        const cfd::Index begin{m.cell_node_offsets()[cell]}, count{m.cell_node_offsets()[cell + 1] - begin};
        for (cfd::Index j = 0; j < count; ++j)
            corners.at(j) = m.nodes()[m.cell_nodes()[begin + j]];
        require(cfd::detail::rpa_angular_quality(std::span{corners}.first(count)) >= o.min_cell_quality,
                "Generated cell violates the minimum corner-sine contract.");
    }
}

struct ParametricCase
{
    const char *name;
    double size{.020};
    double wall{.02};
    const char *upstream{"auto"};
    double downstream{20};
    double outer{10};
    double band{.5};
    double plume{6};
    const char *type{"quadrilateral"};
    double quality{.1};
    unsigned layers{}; // 1: inner, 2: outer, 3: both; zero: automatic distribution.
    bool integration{};
};

std::vector<ParametricCase> parametric_cases()
{
    // One-factor changes isolate failures; two combined cases cover interactions.
    // TRI/HYBRID use an explicit 0.01 minimum corner sine, not the old maximum.
    return {
        {.name = "reference", .size = .005, .integration = true},
        {.name = "resolution_010", .size = .010, .integration = true},
        {.name = "resolution_020"},
        {.name = "wall_001", .wall = .01},
        {.name = "wall_004", .wall = .04, .integration = true},
        {.name = "upstream_2", .upstream = "2"},
        {.name = "upstream_4", .upstream = "4"},
        {.name = "downstream_12", .downstream = 12},
        {.name = "downstream_30", .downstream = 30},
        {.name = "outer_8", .outer = 8},
        {.name = "outer_15", .outer = 15},
        {.name = "band_025", .band = .25},
        {.name = "band_075", .band = .75},
        {.name = "plume_3", .plume = 3, .quality = .05, .integration = true},
        {.name = "plume_4", .plume = 4, .quality = .05},
        {.name = "box_combined", .upstream = "4", .downstream = 30, .outer = 15, .quality = .05, .integration = true},
        {.name = "regions_combined", .band = .75, .plume = 4, .quality = .05, .integration = true},
        {.name = "triangle", .type = "triangle", .quality = .01, .integration = true},
        {.name = "hybrid", .type = "hybrid", .quality = .01, .integration = true},
        // At the reference resolution these first layers exceed the wall-
        // chord errors while fitting the curved-wall offset geometry.
        // The incompatible coarse/thin combination has its own negative test.
        {.name = "inner_layers", .size = .005, .layers = 1},
        {.name = "outer_layers", .size = .005, .layers = 2},
        {.name = "both_layers", .size = .005, .layers = 3, .integration = true}};
}

cfd::input::MeshInput parametric_input(const std::filesystem::path &directory, const ParametricCase &test)
{
    if (std::string_view{test.name} == "reference")
    {
        const auto original{std::filesystem::path{CFD_RPA_REFERENCE_DXF}.parent_path().parent_path().parent_path() /
                            "system/meshDict"};
        auto contents{cfd::test::read_text_file(original)};
        const std::string relative{"constant/geometry/engine.dxf"};
        const auto position{contents.find(relative)};
        require(position != std::string::npos, "Reference DXF entry was not found.");
        contents.replace(position, relative.size(), CFD_RPA_REFERENCE_DXF);
        std::filesystem::create_directories(directory / "system");
        const auto copy{directory / "system/meshDict"};
        {
            std::ofstream output{copy};
            output << contents;
            require(static_cast<bool>(output), "Could not create isolated reference meshDict.");
        }
        return cfd::input::read_mesh_dict(copy);
    }
    std::ostringstream body, mesh, regions;
    body << "engine { file \"" << CFD_RPA_REFERENCE_DXF << "\"; wallThicknessFactor " << test.wall
         << "; } box { upstreamFactor " << test.upstream << "; downstreamFactor " << test.downstream
         << "; outerRadiusFactor " << test.outer << "; }";
    mesh << "cellType " << test.type << "; baseCellSize " << test.size << "; minCellQuality " << test.quality << ';';
    if ((test.layers & 1U) != 0)
        mesh << " innerWallLayers true; innerWallFirstCellHeight 0.0001; innerWallLayerCount 3;"
                " innerWallGrowthRate 1.1; innerWallThickness 0.000331;";
    if ((test.layers & 2U) != 0)
        mesh << " outerWallLayers true; outerWallFirstCellHeight 0.0008; outerWallLayerCount 3;"
                " outerWallGrowthRate 1.1; outerWallThickness 0.002648;";
    regions << "regions { bandThicknessFactor " << test.band << "; plumeOutletRadiusFactor " << test.plume << "; }";
    return dictionary(directory, body.str(), mesh.str(), regions.str());
}

void validate_parametric_plan(const cfd::input::MeshInput &input, const cfd::RpaNozzleDomain &d,
                              const ParametricCase &test)
{
    const auto &g{std::get<cfd::RpaNozzleGeometry>(input.geometry)};
    const auto &options{input.rpa_nozzle_meshing};
    const cfd::RpaRegionLayout layout{d, options.regions};
    const auto plan{cfd::detail::discretize_rpa(d, layout, input.generation_options, options)};
    const double radius{d.inner_wall.back().y}, exit_x{d.inner_wall.back().x};
    require_near(input.generation_options.mesh_size, test.size, 1e-14, "Resolution was not read from meshDict.");
    require_near(options.min_cell_quality, test.quality, 1e-14, "User quality threshold changed.");
    require_near(d.downstream_x - exit_x, test.downstream * radius, 1e-12, "Downstream box factor mismatch.");
    require_near(d.farfield_radius, test.outer * radius, 1e-12, "Radial box factor mismatch.");
    const double upstream_x{g.upstream_factor ? exit_x - *g.upstream_factor * radius : d.inner_wall.front().x};
    require_near(d.outer_wall.front().x, upstream_x, 1e-12, "Upstream box factor mismatch.");
    require(g.upstream_factor.has_value() == (std::string_view{test.upstream} != "auto"), "Auto upstream changed.");
    if (g.upstream_factor)
        require_near(*g.upstream_factor, std::stod(test.upstream), 1e-14, "Explicit upstream factor changed.");
    for (const auto &p : d.outer_wall)
        require_near(p.y - radius_at(d.inner_wall, p.x), test.wall * radius, 1e-12, "Shell thickness mismatch.");
    require_near(plan.lip_y.back() - plan.lip_y.front(), test.wall * radius, 1e-12, "Lip thickness mismatch.");
    require_near(layout.band_thickness(), test.band * radius, 1e-12, "Band factor mismatch.");
    require_near(layout.plume_inlet_radius(), (1 + test.wall + test.band) * radius, 1e-12,
                 "Plume inlet did not follow ENGINE/LIP/BAND.");
    require_near(layout.plume_outlet_radius(), test.plume * radius, 1e-12, "Plume outlet factor mismatch.");
    require(plan.engine.radial_cells == std::max(cfd::Index{8}, static_cast<cfd::Index>(std::ceil(radius / test.size))),
            "Unexpected ENGINE row count for these compatible three-layer prescriptions.");
    require(plan.band.radial_cells ==
                std::max({cfd::Index{4}, static_cast<cfd::Index>(std::ceil(layout.band_thickness() / test.size)),
                          options.outer_wall_layers.active() ? cfd::Index{5} : cfd::Index{0}}),
            "Unexpected BAND row count.");
    const std::array<const cfd::detail::RpaStructuredGrid *, 4> grids{&plan.engine, &plan.band, &plan.plume,
                                                                      &plan.farfield};
    for (const auto *grid : grids)
    {
        require(grid->nodes.size() == grid->x.size() * (grid->radial_cells + 1), "Nonconstant radial cardinality.");
        for (cfd::Index i = 1; i < grid->x.size(); ++i)
        {
            require(grid->x[i] > grid->x[i - 1], "Axial columns overlap.");
            for (cfd::Index j = 0; j < grid->radial_cells; ++j)
            {
                const std::array<cfd::Point2, 4> corners{grid->at(i - 1, j), grid->at(i, j), grid->at(i, j + 1),
                                                         grid->at(i - 1, j + 1)};
                require(cfd::detail::rpa_angular_quality(corners) > 0, "Structured cell has an invalid Jacobian.");
            }
        }
    }
    for (const auto *grid : {&plan.engine, &plan.band})
        for (const double x : d.curve_junctions)
            if (x >= grid->x.front() && x <= grid->x.back())
                require(std::binary_search(grid->x.begin(), grid->x.end(), x), "Mandatory CAD junction lost.");
    require(plan.plume.radial_cells == plan.engine.radial_cells + 1 + plan.band.radial_cells,
            "Plume left cardinality does not include the lip.");
    for (cfd::Index i = 0; i < plan.plume.x.size(); ++i)
        for (cfd::Index j = 0; j <= plan.plume.radial_cells; ++j)
        {
            const double t{(plan.plume.x[i] - exit_x) / (d.downstream_x - exit_x)};
            require_near(plan.plume.at(i, j).y,
                         std::lerp(plan.plume.at(0, j).y, plan.plume.at(plan.plume.x.size() - 1, j).y, t), 1e-12,
                         "Plume row continuity changed.");
        }
    // Independent normal-height checks in the straight chamber, where the
    // prescribed normal is exactly vertical. All curved nodes are also checked
    // against the prescribed grid during full-mesh validation.
    for (const bool inner : {true, false})
    {
        const auto &layers{inner ? options.inner_wall_layers : options.outer_wall_layers};
        if (!layers.active())
            continue;
        const auto &grid{inner ? plan.engine : plan.band};
        const cfd::Index wall{inner ? grid.radial_cells : 0};
        const cfd::Index first{inner ? wall - 1 : 1}, last{inner ? wall - layers.layer_count : layers.layer_count};
        require_near(std::abs(grid.at(0, first).y - grid.at(0, wall).y), *layers.first_cell_height, 1e-12,
                     "Prescribed first normal height changed.");
        require_near(std::abs(grid.at(0, last).y - grid.at(0, wall).y), *layers.thickness, 1e-12,
                     "Prescribed total normal thickness changed.");
    }
}

void test_default_outlet_generation()
{
    const TemporaryDirectory directory;
    for (double size : {.020, .005})
        for (const auto *type : {"quadrilateral", "hybrid"})
        {
            const auto input{
                dictionary(directory.path() / type, "engine { file \"" + std::string{CFD_RPA_REFERENCE_DXF} + "\"; }",
                           "cellType " + std::string{type} + "; baseCellSize " + std::to_string(size) + ";")};
            const auto &geometry{std::get<cfd::RpaNozzleGeometry>(input.geometry)};
            require(geometry.downstream_factor == 80 && geometry.outer_radius_factor == 65 &&
                        input.rpa_nozzle_meshing.regions.plume_outlet_radius_factor == 16 &&
                        input.rpa_nozzle_meshing.min_cell_quality == .05,
                    "Default-generation regression is not using the unmodified defaults.");
            const auto domain{cfd::read_rpa_nozzle_domain(geometry)};
            const auto built{
                cfd::build_mesh(cfd::generate_mesh(input.geometry, input.generation_options, input.automatic_meshing,
                                                   input.backward_facing_step_meshing, input.rpa_nozzle_meshing))};
            validate_generated_mesh(built.mesh, domain, input.generation_options, input.rpa_nozzle_meshing);
            std::cout << "[MESH] defaults_80_65_16 " << type << " size=" << size << " cells=" << built.mesh.cell_count()
                      << '\n';
        }
}

int run_parametric_cases(bool generate, bool complete)
{
    int failures{};
    const TemporaryDirectory directory;
    for (const auto &test : parametric_cases())
    {
        if (generate && !complete && !test.integration)
            continue;
        failures += cfd::test::run_test(test.name, [&] {
            const auto input{parametric_input(directory.path() / test.name, test)};
            const auto domain{cfd::read_rpa_nozzle_domain(std::get<cfd::RpaNozzleGeometry>(input.geometry))};
            validate_parametric_plan(input, domain, test);
            if (!generate)
                return;
            const auto built{
                cfd::build_mesh(cfd::generate_mesh(input.geometry, input.generation_options, input.automatic_meshing,
                                                   input.backward_facing_step_meshing, input.rpa_nozzle_meshing))};
            const auto &mesh{built.mesh};
            validate_generated_mesh(mesh, domain, input.generation_options, input.rpa_nozzle_meshing);
            const auto statistics{cfd::compute_mesh_statistics(mesh)};
            const double mean{std::accumulate(mesh.cell_qualities().begin(), mesh.cell_qualities().end(), 0.0) /
                              static_cast<double>(mesh.cell_count())};
            std::cout << "[MESH] " << test.name << " cells=" << mesh.cell_count()
                      << " quality_min=" << *std::ranges::min_element(mesh.cell_qualities()) << " quality_mean=" << mean
                      << " nonorth_max=" << statistics.internal_face_non_orthogonality_degrees.maximum << '\n';
            if (std::string_view{test.name} == "reference")
            {
                const auto path{directory.path() / "reference.vtu"};
                cfd::write_vtu(mesh, path);
                require(std::filesystem::file_size(path) > 0, "Empty reference VTU.");
            }
        });
    }
    return failures;
}

void test_invalid_parametric_inputs()
{
    const TemporaryDirectory directory;
    const std::string engine{"engine { file \"" + std::string{CFD_RPA_REFERENCE_DXF} + "\"; }"};
    for (const auto *mesh : {"cellType quadrilateral; baseCellSize 0;", "cellType quadrilateral; baseCellSize -1;",
                             "cellType quadrilateral; baseCellSize .02; minCellQuality 0;",
                             "cellType quadrilateral; baseCellSize .02; minCellQuality -0.1;",
                             "cellType quadrilateral; baseCellSize .02; minCellQuality 1.1;"})
        require_throws<std::runtime_error>([&] { static_cast<void>(dictionary(directory.path(), engine, mesh)); },
                                           "A: invalid meshDict configuration accepted.");
    for (const auto *factor : {"wallThicknessFactor 0;", "wallThicknessFactor -1;"})
        require_throws<std::runtime_error>(
            [&] {
                static_cast<void>(dictionary(directory.path(), "engine { file test.dxf; " + std::string{factor} + "}"));
            },
            "A: nonpositive geometry factor accepted.");
    for (const auto &test : std::vector<ParametricCase>{{.name = "small_plume", .plume = 1},
                                                        {.name = "small_box", .outer = 4},
                                                        {.name = "large_band", .band = 8},
                                                        {.name = "unavailable_upstream", .upstream = "6"}})
        require_throws<std::invalid_argument>(
            [&] {
                const auto input{parametric_input(directory.path() / test.name, test)};
                const auto domain{cfd::read_rpa_nozzle_domain(std::get<cfd::RpaNozzleGeometry>(input.geometry))};
                const cfd::RpaRegionLayout layout{domain, input.rpa_nozzle_meshing.regions};
                static_cast<void>(layout);
            },
            "B: impossible geometry accepted.");
    const auto input{parametric_input(directory.path() / "layers", {.name = "layers"})};
    const auto domain{cfd::read_rpa_nozzle_domain(std::get<cfd::RpaNozzleGeometry>(input.geometry))};
    auto options{input.rpa_nozzle_meshing};
    options.inner_wall_layers = {.enabled = true, .first_cell_height = 1, .layer_count = 3, .growth_rate = 1.1};
    require_throws<std::invalid_argument>(
        [&] {
            static_cast<void>(cfd::detail::discretize_rpa(domain, cfd::RpaRegionLayout{domain, options.regions},
                                                          input.generation_options, options));
        },
        "B: impossible wall layers accepted.");
    options = input.rpa_nozzle_meshing;
    options.regions.plume_outlet_radius_factor = 3;
    cfd::test::require_throws_with_message<std::runtime_error>(
        [&] {
            static_cast<void>(cfd::detail::discretize_rpa(domain, cfd::RpaRegionLayout{domain, options.regions},
                                                          input.generation_options, options));
        },
        "outlet angle", "C: incompatible angular threshold accepted.");
    const cfd::RpaNozzleGeometry unused{.file = "not-opened.dxf"};
    cfd::test::require_throws_with_message<std::invalid_argument>(
        [&] {
            static_cast<void>(
                cfd::generate_mesh(unused, {.mesh_size = .02, .cell_type = cfd::CellType::Triangle, .hybrid = true}));
        },
        "HYBRID requires", "A: contradictory RPA topology options accepted.");
}

void test_external_size_grid_format()
{
    const TemporaryDirectory directory;
    const auto path{directory.path() / "sizes with spaces.txt"};
    {
        std::ofstream output{path};
        output << "2 2\n0 2 0 4\n1 3\n5 7\n";
    }
    const auto field{cfd::read_mesh_size_field(path)};
    require_near(field({1, 2}), 4, 1e-14, "Bilinear size-grid interpolation changed.");
    require_near(field({-1, -1}), 1, 1e-14, "Size-grid constant extension changed.");
    require_near(field({3, 5}), 7, 1e-14, "Size-grid upper extension changed.");
    for (const auto *invalid : {"CFD_SIZE_GRID_V1 2 2 0 2 0 4 1 3 5 7", "2 2 0 2 0 4 1 3 -5 7", "2 2 0 2 0 4 1 3 5",
                                "2 2 0 2 0 4 1 3 5 7 extra"})
    {
        {
            std::ofstream output{path};
            output << invalid;
        }
        require_throws<std::invalid_argument>([&] { static_cast<void>(cfd::read_mesh_size_field(path)); },
                                              "Invalid external size-grid format accepted.");
    }
}

void test_acute_transition_rejects_strict_triangle_threshold()
{
    const TemporaryDirectory directory;
    const auto input{
        parametric_input(directory.path(), {.name = "acute_hybrid", .plume = 3, .type = "hybrid", .quality = .1})};
    cfd::test::require_throws_with_message<std::runtime_error>(
        [&] {
            static_cast<void>(cfd::generate_mesh(input.geometry, input.generation_options, input.automatic_meshing,
                                                 input.backward_facing_step_meshing, input.rpa_nozzle_meshing));
        },
        "outlet angle", "C: acute TRANSITION was falsely accepted as a high-quality TRI mesh.");
    std::cout << "[EXPECTED REJECTION C] acute Hybrid TRANSITION with minimum corner sine 0.1\n";
}

void test_thin_layers_are_not_certified_inside_the_cad_shell()
{
    const TemporaryDirectory directory;
    // These offsets fit the region heights, but the 0.020 m discretization's
    // wall-chord error exceeds the first-layer scale. Preserve this failure as
    // a negative case instead of relaxing the physical-domain invariant.
    const auto input{dictionary(
        directory.path(),
        "engine { file \"" + std::string{CFD_RPA_REFERENCE_DXF} +
            "\"; } box { downstreamFactor 20; outerRadiusFactor 10; }",
        "cellType quadrilateral; baseCellSize .020; minCellQuality .1; "
        "innerWallLayers true; innerWallFirstCellHeight .0001; innerWallLayerCount 3; innerWallGrowthRate 1.1; "
        "outerWallLayers true; outerWallFirstCellHeight .00015; outerWallLayerCount 3; outerWallGrowthRate 1.1;",
        "regions {bandThicknessFactor .5; plumeOutletRadiusFactor 6;}")};
    const auto built{
        cfd::build_mesh(cfd::generate_mesh(input.geometry, input.generation_options, input.automatic_meshing,
                                           input.backward_facing_step_meshing, input.rpa_nozzle_meshing))};
    const auto domain{cfd::read_rpa_nozzle_domain(std::get<cfd::RpaNozzleGeometry>(input.geometry))};
    std::ostringstream output;
    cfd::test::require_throws_with_message<std::runtime_error>(
        [&] {
            cfd::write_rpa_mesh_diagnostics(output, built.mesh, domain, input.generation_options,
                                            input.rpa_nozzle_meshing);
        },
        "outside fluid domain or inside fictitious shell",
        "C: the preview diagnostics certified a thin-layer mesh with cell centres in the CAD shell.");
    std::cout << "[EXPECTED REJECTION C] thin layers versus coarse wall chords, rejected by preview diagnostics\n";
}

void test_profiles_layers_and_failures()
{
    for (const auto *file : {"conical.dxf", "cubic.dxf"})
    {
        const auto g{geometry(std::filesystem::path{CFD_RPA_TEST_FIXTURES} / file)};
        const auto d{cfd::read_rpa_nozzle_domain(g)};
        const auto built{cfd::build_mesh(
            cfd::generate_mesh(g, {.mesh_size = .15, .cell_type = cfd::CellType::Quadrilateral, .hybrid = true},
                               {.regions = {.plume_outlet_radius_factor = 6}, .min_cell_quality = .01}))};
        validate_mesh(built.mesh, d, cfd::CellType::Quadrilateral, true);
    }
    const auto g{geometry(std::filesystem::path{CFD_RPA_TEST_FIXTURES} / "conical.dxf")};
    const auto d{cfd::read_rpa_nozzle_domain(g)};
    cfd::RpaNozzleMeshingOptions o{.regions = {.plume_outlet_radius_factor = 6}};
    o.inner_wall_layers = {.enabled = true, .first_cell_height = .0025, .layer_count = 3, .growth_rate = 1.1};
    o.outer_wall_layers = {.enabled = true, .first_cell_height = .003, .layer_count = 4, .growth_rate = 1.1};
    const auto built{
        cfd::build_mesh(cfd::generate_mesh(g, {.mesh_size = .15, .cell_type = cfd::CellType::Quadrilateral}, o))};
    validate_mesh(built.mesh, d, cfd::CellType::Quadrilateral);
    const auto &m{built.mesh};
    for (cfd::Index face = 0; face < m.face_count(); ++face)
        if (m.face_boundary_ids()[face] == 2 || m.face_boundary_ids()[face] == 3)
        {
            auto cell = m.face_adjacencies()[face].owner;
            const auto &c{m.face_centers()[face]};
            const auto &n{m.face_area_vectors()[face]};
            double height{};
            for (auto j = m.cell_node_offsets()[cell]; j < m.cell_node_offsets()[cell + 1]; ++j)
            {
                const auto &p{m.nodes()[m.cell_nodes()[j]]};
                height = std::max(height, ((c.x - p.x) * n.x + (c.y - p.y) * n.y) / m.face_lengths()[face]);
            }
            require_near(height, m.face_boundary_ids()[face] == 2 ? .0025 : .003, 1e-8,
                         "Prescribed first normal height changed.");
        }
    auto bad{o};
    bad.inner_wall_layers.first_cell_height = 10;
    require_throws<std::invalid_argument>(
        [&] {
            static_cast<void>(
                cfd::generate_mesh(g, {.mesh_size = .15, .cell_type = cfd::CellType::Quadrilateral}, bad));
        },
        "Impossible layers altered silently.");
    bad = {};
    bad.regions = o.regions;
    bad.target_cell_budget = 1;
    require_throws<std::runtime_error>(
        [&] {
            static_cast<void>(
                cfd::generate_mesh(g, {.mesh_size = .15, .cell_type = cfd::CellType::Quadrilateral}, bad));
        },
        "Impossible budget accepted.");
    cfd::AutomaticMeshingOptions legacy;
    legacy.wall_refinement_factor = .25;
    require_throws<std::invalid_argument>(
        [&] { static_cast<void>(cfd::generate_mesh(cfd::GeometryInput{g}, {.mesh_size = .15}, legacy, {})); },
        "Removed nondefault automatic policy ignored by RPA dispatch.");
    require_throws<std::invalid_argument>(
        [&] {
            static_cast<void>(cfd::generate_mesh(g, {.mesh_size = .15, .cell_type = cfd::CellType::Quadrilateral},
                                                 {.regions = o.regions}, [](const cfd::Point2 &) { return 0.; }));
        },
        "Invalid external field accepted.");
}
} // namespace

int main(const int argc, char *argv[])
try
{
    int failures{};
    const std::span<char *> arguments{argv, static_cast<std::size_t>(argc)};
    const std::string_view mode{arguments.size() == 2 ? arguments[1] : ""};
    if (arguments.size() > 2 || (!mode.empty() && mode != "--integration" && mode != "--parametric"))
    {
        std::cerr << "Usage: cfd_rpa_nozzle_tests [--integration|--parametric]\n";
        return 1;
    }
    if (!mode.empty())
    {
        failures += cfd::test::run_test("RPA default 80/65/16 generation", test_default_outlet_generation);
        failures += run_parametric_cases(true, mode == "--parametric");
        failures += cfd::test::run_test("strict TRI threshold rejects acute outlet",
                                        test_acute_transition_rejects_strict_triangle_threshold);
        if (mode == "--parametric")
        {
            failures += cfd::test::run_test("thin layers cannot bypass CAD-domain diagnostic",
                                            test_thin_layers_are_not_certified_inside_the_cad_shell);
            failures += cfd::test::run_test("RPA profiles/layers/budget/field", test_profiles_layers_and_failures);
        }
        return cfd::test::finish_tests(failures, "RPA integration/parametric");
    }
    failures += cfd::test::run_test("RPA nested reader/defaults/migration", test_reader_and_migration);
    failures += cfd::test::run_test("RPA DXF geometry/factors/layout", test_geometry_and_layout);
    failures += cfd::test::run_test("RPA pure radial/axial/quality/scale", test_discretization_and_scale);
    failures += cfd::test::run_test("RPA TRANSITION minimum outlet factor", test_transition_minimum_outlet_factor);
    failures += cfd::test::run_test("RPA arithmetic mean/chord diagnostic/FARFIELD invariance",
                                    test_mean_height_columns_and_chord_diagnostics);
    failures += cfd::test::run_test("RPA invalid configurations A/B/C", test_invalid_parametric_inputs);
    failures += cfd::test::run_test("RPA external-size format", test_external_size_grid_format);
    failures += run_parametric_cases(false, true);
    return cfd::test::finish_tests(failures, "RPA nozzle");
}
catch (const std::exception &error)
{
    std::cerr << "[UNEXPECTED ERROR D] RPA test setup: " << error.what() << '\n';
    return 1;
}
