#include "cfd/input/OpenFOAMCaseReader.hpp"
#include "cfd/io/VtkWriter.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/mesh/MeshStatistics.hpp"
#include "cfd/meshing/GmshMesher.hpp"
#include "cfd/meshing/RawMeshValidation.hpp"
#include "cfd/meshing/RpaNozzleGeometry.hpp"
#include "support/TestUtils.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
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

cfd::RpaNozzleGeometry reference_geometry()
{
    return {.file = CFD_RPA_REFERENCE_DXF};
}

cfd::input::MeshInput read_dict(const std::filesystem::path &directory, const std::string_view geometry,
                                const std::string_view meshing = "")
{
    std::filesystem::create_directories(directory / "system");
    const auto file{directory / "system/meshDict"};
    {
        std::ofstream output{file};
        output << "FoamFile { format ascii; class dictionary; object meshDict; }\ngeometry { " << geometry
               << " }\nmesh { cellType quadrilateral; size 0.08; " << meshing << " }\n";
        require(static_cast<bool>(output), "Unable to write RPA test dictionary.");
    }
    return cfd::input::read_mesh_dict(file);
}

void test_parser_defaults_overrides_paths_and_rejections()
{
    const TemporaryDirectory directory;
    const auto input{read_dict(directory.path(), "file \"constant/geometry/engine.dxf\"; type rpaNozzle;")};
    const auto &geometry{std::get<cfd::RpaNozzleGeometry>(input.geometry)};
    require(geometry.file == directory.path() / "constant/geometry/engine.dxf", "DXF path was not case-relative.");
    require(!geometry.inlet_position && !geometry.outer_upstream_length && !geometry.wall_thickness &&
                !geometry.downstream_length && !geometry.farfield_radius,
            "Geometry defaults are not automatic.");
    require(input.automatic_meshing.enabled, "Automatic RPA meshing was not selected by default.");
    require_near(input.rpa_nozzle_meshing.throat_refinement_factor, 0.25, 0.0, "Invalid default throat policy.");
    const auto overridden{read_dict(
        directory.path(),
        "type rpaNozzle; file \"constant/geometry/engine.dxf\"; dxfUnits mm; geometricTolerance 2e-6; "
        "inletPosition 0.01; outerUpstreamLength 0.3; wallThickness 0.004; downstreamLength 2; farfieldRadius 1;",
        "automaticMeshing false; maximumGrowthRate 1.3; wallRefinementFactor 0.4; wallRefinementLayers 5; "
        "convergentRefinementFactor 0.8; throatRefinementFactor 0.2; divergentRefinementFactor 0.9; "
        "exitRefinementFactor 0.3; plumeRefinementFactor 0.6; throatRefinementLength 0.08; "
        "exitRefinementLength 0.1; nearPlumeLength 0.4;")};
    const auto &overridden_geometry{std::get<cfd::RpaNozzleGeometry>(overridden.geometry)};
    require(overridden_geometry.inlet_position == 0.01 && overridden_geometry.outer_upstream_length == 0.3 &&
                overridden_geometry.wall_thickness == 0.004 && overridden_geometry.downstream_length == 2.0 &&
                overridden_geometry.farfield_radius == 1.0 && overridden_geometry.dxf_units == "mm" &&
                overridden_geometry.geometric_tolerance == 2.0e-6,
            "Geometry overrides were lost.");
    const auto &policy{overridden.rpa_nozzle_meshing};
    require(!overridden.automatic_meshing.enabled && overridden.automatic_meshing.maximum_growth_rate == 1.3 &&
                overridden.automatic_meshing.wall_refinement_factor == 0.4 &&
                overridden.automatic_meshing.wall_refinement_layers == 5 &&
                policy.convergent_refinement_factor == 0.8 && policy.throat_refinement_factor == 0.2 &&
                policy.divergent_refinement_factor == 0.9 && policy.exit_refinement_factor == 0.3 &&
                policy.plume_refinement_factor == 0.6 && policy.throat_refinement_length == 0.08 &&
                policy.exit_refinement_length == 0.1 && policy.near_plume_length == 0.4,
            "Mesh overrides were lost.");
    for (const std::string_view invalid :
         {"type rpaNozzle;", "type rpaNozzle; file x; wallThickness 0;", "type rpaNozzle; file x; downstreamLength -1;",
          "type rpaNozzle; file x; dxfUnits ambiguous;", "type rpaNozzle; file x; inletPosition 0; inletPosition 0.1;",
          "type rpaNozzle; file x; length 1;", "type rectangle; length 1; height 1; file x;",
          "type rpaNozzle; file x; geometricTolerance nan;"})
    {
        require_throws<std::runtime_error>([&]() { static_cast<void>(read_dict(directory.path(), invalid)); },
                                           "Invalid RPA dictionary was accepted.");
    }
    for (const std::string_view invalid :
         {"throatRefinementFactor 0;", "exitRefinementFactor 2;",
          "plumeRefinementFactor 0.5; plumeRefinementFactor 0.2;", "nearPlumeLength -1;", "unknown 1;"})
    {
        require_throws<std::runtime_error>(
            [&]() { static_cast<void>(read_dict(directory.path(), "type rpaNozzle; file x;", invalid)); },
            "Invalid RPA mesh policy was accepted.");
    }
    std::filesystem::create_directories(directory.path() / "constant/geometry");
    const auto quoted_file{directory.path() / "constant/geometry/engine ' profile.dxf"};
    std::filesystem::copy_file(CFD_RPA_REFERENCE_DXF, quoted_file);
    const auto quoted{read_dict(directory.path(), "type rpaNozzle; file \"constant/geometry/engine ' profile.dxf\";")};
    require(cfd::read_rpa_nozzle_domain(std::get<cfd::RpaNozzleGeometry>(quoted.geometry)).inner_wall.size() > 3,
            "CAD invocation did not handle spaces/apostrophe in a case-relative path.");
}

void test_real_domain_and_invalid_options()
{
    const auto geometry{reference_geometry()};
    const auto domain{cfd::read_rpa_nozzle_domain(geometry)};
    const auto exit{domain.inner_wall.back()};
    const double diameter{2.0 * exit.y};
    require_near(2.0 * domain.throat.y, 0.02858, 3e-6, "Incorrect DXF throat diameter.");
    require_near(diameter, 0.28582, 5e-6, "Incorrect DXF exit diameter.");
    require_near((exit.y / domain.throat.y) * (exit.y / domain.throat.y), 100.0, 1e-4, "Incorrect DXF area ratio.");
    require_near(exit.x - domain.throat.x, 0.38458, 3e-6, "Incorrect divergent length.");
    require_near(domain.outer_wall.front().x, exit.x - 1.5 * diameter, 1e-12, "Incorrect exterior upstream default.");
    require_near(domain.outer_wall.back().y - exit.y, 0.01 * diameter, 1e-12, "Wall thickness is not radial.");
    require_near(domain.downstream_x - exit.x, 10.0 * diameter, 1e-12, "Incorrect downstream default.");
    require_near(domain.farfield_radius, 5.0 * diameter, 1e-12, "Incorrect radial farfield default.");
    auto absent{geometry};
    absent.file = geometry.file.parent_path() / "absent.dxf";
    require_throws<std::runtime_error>([&]() { static_cast<void>(cfd::read_rpa_nozzle_domain(absent)); },
                                       "Missing DXF accepted.");
    for (const double value :
         {0.0, -1.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()})
    {
        auto bad{geometry};
        bad.wall_thickness = value;
        require_throws<std::invalid_argument>([&]() { static_cast<void>(cfd::read_rpa_nozzle_domain(bad)); },
                                              "Invalid thickness accepted.");
    }
    for (const auto &bad : {cfd::RpaNozzleGeometry{.file = geometry.file, .inlet_position = domain.throat.x},
                            cfd::RpaNozzleGeometry{.file = geometry.file, .outer_upstream_length = 2.0},
                            cfd::RpaNozzleGeometry{.file = geometry.file, .farfield_radius = exit.y}})
    {
        require_throws<std::invalid_argument>([&]() { static_cast<void>(cfd::read_rpa_nozzle_domain(bad)); },
                                              "Impossible domain extent accepted.");
    }
}

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

void validate_mesh(const cfd::Mesh &mesh, const cfd::RpaNozzleDomain &domain, const cfd::CellType type)
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
        require(mesh.cell_types()[cell] == type && mesh.cell_areas()[cell] > 0, "Unexpected topology or invalid area.");
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
    require_near(statistics.total_cell_area, expected_area, 1e-4, "Fluid area disagrees with the meridian domain.");
}

double wall_size(const cfd::Mesh &mesh, const double x, const double width)
{
    double sum{};
    cfd::Index count{};
    for (cfd::Index face = 0; face < mesh.face_count(); ++face)
    {
        const auto id{mesh.face_boundary_ids()[face]};
        if ((id == 2 || id == 3) && std::abs(mesh.face_centers()[face].x - x) < width)
        {
            sum += mesh.face_lengths()[face];
            ++count;
        }
    }
    require(count > 0, "Size diagnostic has no samples.");
    return sum / static_cast<double>(count);
}

void test_real_mesh_types_refinement_and_export()
{
    const auto geometry{reference_geometry()};
    const auto domain{cfd::read_rpa_nozzle_domain(geometry)};
    const TemporaryDirectory output;
    for (const auto type : {cfd::CellType::Triangle, cfd::CellType::Quadrilateral})
    {
        const cfd::MeshGenerationOptions options{0.08, type};
        auto base{cfd::build_mesh(cfd::generate_mesh(geometry, options, {.enabled = false}))};
        auto refined{cfd::build_mesh(cfd::generate_mesh(geometry, options))};
        cfd::AutomaticMeshingOptions wall_policy;
        wall_policy.wall_refinement_factor = 0.25;
        cfd::RpaNozzleMeshingOptions policy;
        // Choose targets below the lip/junction geometry-driven size cap so
        // this exercises a real user-requested refinement, not an inactive cap.
        policy.throat_refinement_factor = 0.03125;
        policy.exit_refinement_factor = 0.03125;
        auto overridden{cfd::build_mesh(cfd::generate_mesh(geometry, options, wall_policy, policy))};
        std::cout << "Validating RPA base type=" << static_cast<int>(type) << " cells=" << base.mesh.cell_count()
                  << '\n';
        validate_mesh(base.mesh, domain, type);
        validate_mesh(refined.mesh, domain, type);
        validate_mesh(overridden.mesh, domain, type);
        for (const double x : {0.1, domain.throat.x, domain.inner_wall.back().x})
        {
            const double basic{wall_size(base.mesh, x, 0.06)};
            const double automatic{wall_size(refined.mesh, x, 0.06)};
            const double override_size{wall_size(overridden.mesh, x, 0.06)};
            std::cout << "RPA type=" << static_cast<int>(type) << " x=" << x << " mean boundary size: " << basic
                      << " -> " << automatic << " -> " << override_size << '\n';
            require(automatic < 0.9 * basic, "Automatic meshing did not measurably refine wall/throat/exit.");
            require(override_size < 0.9 * automatic, "User refinement override did not measurably reduce sizes.");
        }
        const auto stats{cfd::compute_mesh_statistics(refined.mesh)};
        std::cout << "RPA cells=" << refined.mesh.cell_count() << " area=" << stats.total_cell_area
                  << " neighbor ratio mean/max=" << stats.internal_face_neighbor_cell_size_ratios.mean << '/'
                  << stats.internal_face_neighbor_cell_size_ratios.maximum << '\n';
        const auto file{output.path() / "mesh.vtu"};
        cfd::write_vtu(refined.mesh, file);
        require(std::filesystem::file_size(file) > 0, "RPA VTU export is empty.");
    }
    cfd::RpaNozzleMeshingOptions bad;
    bad.near_plume_length = 1000;
    require_throws<std::invalid_argument>(
        [&]() { static_cast<void>(cfd::generate_mesh(geometry, {0.08, cfd::CellType::Triangle}, {}, bad)); },
        "Impossible plume refinement extent accepted.");
    require_throws<std::invalid_argument>(
        [&]() {
            static_cast<void>(
                cfd::generate_mesh(cfd::GeometryInput{geometry}, {0.08, cfd::CellType::Triangle}, {}, {}, bad));
        },
        "Geometry dispatch silently ignored the RPA policy.");
}

void test_synthetic_meshes()
{
    for (const std::string_view file : {"conical.dxf", "cubic.dxf"})
    {
        const cfd::RpaNozzleGeometry geometry{.file = std::filesystem::path{CFD_RPA_TEST_FIXTURES} / file};
        const auto domain{cfd::read_rpa_nozzle_domain(geometry)};
        for (const auto type : {cfd::CellType::Triangle, cfd::CellType::Quadrilateral})
        {
            auto built{cfd::build_mesh(cfd::generate_mesh(geometry, {0.3, type}, {.enabled = false}))};
            validate_mesh(built.mesh, domain, type);
        }
    }
}
} // namespace

int main()
{
    int failures{};
    failures +=
        cfd::test::run_test("RPA parser/defaults/overrides/paths", test_parser_defaults_overrides_paths_and_rejections);
    failures += cfd::test::run_test("RPA domain/defaults/rejections", test_real_domain_and_invalid_options);
    failures +=
        cfd::test::run_test("RPA triangles/quads/refinement/export", test_real_mesh_types_refinement_and_export);
    failures += cfd::test::run_test("RPA synthetic conical/cubic meshes", test_synthetic_meshes);
    return cfd::test::finish_tests(failures, "RPA nozzle");
}
