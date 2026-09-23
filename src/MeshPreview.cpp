#include "cfd/input/OpenFOAMCaseReader.hpp"
#include "cfd/io/MeshReport.hpp"
#include "cfd/io/VtkWriter.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/mesh/MeshStatistics.hpp"
#include "cfd/meshing/GmshMesher.hpp"
#include "cfd/meshing/RawMeshData.hpp"

#include <exception>
#include <filesystem>
#include <iostream>
#include <span>
#include <utility>

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
                                                     input.backward_facing_step_meshing)};

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
        cfd::write_vtu(mesh, output_file);

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
