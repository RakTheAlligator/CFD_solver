#include "cfd/input/OpenFOAMCaseReader.hpp"
#include "cfd/io/MeshReport.hpp"
#include "cfd/io/VtkWriter.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/mesh/MeshStatistics.hpp"
#include "cfd/meshing/GmshMesher.hpp"
#include "cfd/meshing/RawMeshData.hpp"
#include "cfd/meshing/RpaMeshDiagnostics.hpp"

#include <array>
#include <exception>
#include <filesystem>
#include <iostream>
#include <optional>
#include <span>
#include <utility>
#include <vector>

int main(const int argc, char *argv[])
{
    const std::span<char *> arguments{argv, static_cast<std::size_t>(argc)};
    if (arguments.size() != 2)
    {
        std::cerr << "Usage: cfd_mesh_preview <case-directory>\n";
        return 1;
    }

    try
    {
        const std::filesystem::path case_directory{arguments[1]};
        const cfd::input::MeshInput input{cfd::input::read_mesh_dict(case_directory / "system" / "meshDict")};
        cfd::RawMeshData raw_mesh{cfd::generate_mesh(input.geometry, input.generation_options, input.automatic_meshing,
                                                     input.backward_facing_step_meshing, input.rpa_nozzle_meshing)};

        std::cout << "[Mesh generation]\n"
                  << "  Nodes             : " << raw_mesh.nodes.size() << '\n'
                  << "  Cells             : " << raw_mesh.cell_types.size() << '\n';

        cfd::MeshBuildResult build_result{cfd::build_mesh(std::move(raw_mesh))};
        const cfd::Mesh &mesh{build_result.mesh};
        const cfd::MeshStatistics statistics{cfd::compute_mesh_statistics(mesh)};
        cfd::write_mesh_report(std::cout, mesh, statistics, build_result.timings);

        const std::filesystem::path output_directory{case_directory / "results"};
        const std::filesystem::path output_file{output_directory / "mesh.vtu"};
        std::filesystem::create_directories(output_directory);
        std::vector<double> regions;
        if (const auto *geometry = std::get_if<cfd::RpaNozzleGeometry>(&input.geometry))
        {
            regions.resize(mesh.cell_count());
            cfd::write_rpa_mesh_diagnostics(std::cout, mesh, cfd::read_rpa_nozzle_domain(*geometry),
                                            input.generation_options, input.rpa_nozzle_meshing, regions);
        }
        const std::array<cfd::VtkCellScalarData, 1> scalars{cfd::VtkCellScalarData{"rpa_region", regions}};
        cfd::write_vtu(mesh, output_file,
                       cfd::VtkCellData{.scalars = std::span{scalars}.first(regions.empty() ? 0 : 1)});

        std::cout << "\n[Output]\n"
                  << "  Mesh              : " << output_file.string() << '\n';
    }
    catch (const std::exception &error)
    {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }

    return 0;
}
