#include "support/EulerRpaExploratory.hpp"
#include "support/EulerRpaTemporal.hpp"

#include "cfd/input/OpenFOAMCaseReader.hpp"
#include "cfd/io/MeshReport.hpp"
#include "cfd/io/VtkWriter.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/mesh/MeshStatistics.hpp"
#include "cfd/meshing/RpaMeshDiagnostics.hpp"

#include <array>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <span>
#include <tuple>
#include <variant>

// Geometry-only preflight. All geometric/resolution settings come from meshDict.
int main(int argc, char **argv)
{
    try
    {
        if (argc != 3)
            throw std::invalid_argument("Usage: cfd_rpa_hybrid_mesh CASE_DIR NEW_OUTPUT_DIR");
        const std::span<char *> arguments{argv, static_cast<cfd::Index>(argc)};
        const std::filesystem::path input{arguments[1]}, output{arguments[2]};
        cfd::verification::require_external_empty_output(input, output);
        const auto settings{cfd::input::read_mesh_dict(input / "system" / "meshDict")};
        const auto *geometry{std::get_if<cfd::RpaNozzleGeometry>(&settings.geometry)};
        if (geometry == nullptr)
            throw std::invalid_argument("Expected rpaNozzle meshDict.");
        using Clock = std::chrono::steady_clock;
        auto begin{Clock::now()};
        auto raw{cfd::generate_mesh(settings.geometry, settings.generation_options, settings.automatic_meshing,
                                    settings.backward_facing_step_meshing, settings.rpa_nozzle_meshing)};
        const double generation{std::chrono::duration<double>(Clock::now() - begin).count()};
        begin = Clock::now();
        const auto built{cfd::build_mesh(std::move(raw))};
        const double build{std::chrono::duration<double>(Clock::now() - begin).count()};
        const auto &mesh{built.mesh};
        cfd::verification::require_connected_rpa(mesh);
        const cfd::AxisymmetricGeometry metrics{mesh};
        const auto domain{cfd::read_rpa_nozzle_domain(*geometry)};
        const cfd::RpaRegionLayout layout{domain, settings.rpa_nozzle_meshing.regions};
        const auto conditions{cfd::verification::rpa_euler_conditions(
            mesh.boundary_groups(), cfd::verification::EulerRpaModel::ambient(300, 0, 0))};
        const auto statistics{cfd::compute_mesh_statistics(mesh)};
        if (settings.generation_options.hybrid &&
            (statistics.triangle_cell_count == 0 || statistics.quadrilateral_cell_count == 0))
            throw std::runtime_error("Requested HYBRID must contain both TRI and QUAD.");
        std::filesystem::create_directories(output);
        std::ofstream report{output / "mesh-report.txt"};
        report.exceptions(std::ios::failbit | std::ios::badbit);
        report << std::setprecision(17) << "generation_s=" << generation << "\nmesh_build_s=" << build
               << "\nDe=" << 2 * domain.inner_wall.back().y << "\nexit_x=" << domain.inner_wall.back().x
               << "\nx_upstream=" << domain.outer_wall.front().x << "\nx_downstream=" << domain.downstream_x
               << "\nr_max=" << domain.farfield_radius << "\nplume_inlet_radius=" << layout.plume_inlet_radius()
               << "\nplume_outlet_radius=" << layout.plume_outlet_radius() << "\nplume_end_x=" << layout.plume_end_x()
               << '\n';
        cfd::write_mesh_report(report, mesh, statistics, built.timings);
        std::vector<double> regions(mesh.cell_count()), sizes(mesh.cell_count()), blocks(mesh.cell_count());
        cfd::write_rpa_mesh_diagnostics(report, mesh, domain, settings.generation_options, settings.rpa_nozzle_meshing,
                                        regions, blocks);
        for (const auto &[name, x, radius] : std::array<std::tuple<const char *, double, double>, 3>{
                 {{"inlet", domain.inner_wall.front().x, domain.inner_wall.front().y},
                  {"throat", domain.throat.x, domain.throat.y},
                  {"exit", domain.inner_wall.back().x, domain.inner_wall.back().y}}})
        {
            const cfd::verification::RpaFaceSection section{mesh, metrics, regions, x, radius};
            report << name << "_section_faces=" << section.face_count() << '\n';
            if (settings.rpa_nozzle_meshing.engine_radial_cells &&
                section.face_count() != *settings.rpa_nozzle_meshing.engine_radial_cells)
                throw std::runtime_error("Organized section count differs from engineRadialCells.");
        }
        std::array<cfd::Index, 8> counts{};
        for (cfd::Index face = 0; face < mesh.face_count(); ++face)
            if (mesh.face_adjacencies()[face].is_boundary())
                ++counts.at(mesh.face_boundary_ids()[face]);
        for (const auto &group : mesh.boundary_groups())
        {
            if (counts.at(group.id) == 0)
                throw std::runtime_error("Empty physical boundary group.");
            report << "boundary " << group.id << ' ' << group.name << " faces=" << counts.at(group.id)
                   << " EulerBoundaryConditionType=" << static_cast<int>(conditions[group.id]) << '\n';
        }
        for (cfd::Index cell = 0; cell < mesh.cell_count(); ++cell)
            sizes[cell] = std::sqrt(mesh.cell_areas()[cell]);
        const std::array scalars{cfd::VtkCellScalarData{"rpa_region", regions},
                                 cfd::VtkCellScalarData{"cell_size_m", sizes},
                                 cfd::VtkCellScalarData{"axisymmetric_volume", metrics.cell_volumes()},
                                 cfd::VtkCellScalarData{"rpa_block", blocks}};
        cfd::write_vtu(mesh, output / "mesh.vtu", {.scalars = scalars});
        report << "PASS: connected TRI/QUAD mesh, axisymmetric identities and three exact sections.\n";
        std::cout << "PASS cells=" << mesh.cell_count() << " TRI=" << statistics.triangle_cell_count
                  << " QUAD=" << statistics.quadrilateral_cell_count << " generation_s=" << generation
                  << " build_s=" << build << " output=" << output << '\n';
    }
    catch (const std::exception &error)
    {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
