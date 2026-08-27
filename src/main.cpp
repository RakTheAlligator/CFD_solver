#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/mesh/MeshStatistics.hpp"

#include "cfd/meshing/GmshMesher.hpp"

#include "cfd/io/MeshReport.hpp"
#include "cfd/io/VtkWriter.hpp"

#include <chrono>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string_view>
#include <utility>

namespace
{

[[nodiscard]]
std::string_view cell_type_name(const cfd::CellType cell_type)
{
    switch (cell_type)
    {
    case cfd::CellType::Triangle:
        return "triangles";

    case cfd::CellType::Quadrilateral:
        return "quadrilaterals";

    default:
        return "unknown";
    }
}

} // namespace

int main()
{
    try
    {
        const cfd::RectangleGeometry geometry{
            .length = 5.0,
            .height = 1.0,
        };

        const cfd::MeshGenerationOptions options{
            .mesh_size = 0.2,
            .cell_type = cfd::CellType::Quadrilateral,
        };

        std::cout << "============================================================\n"
                  << " CFD Solver\n"
                  << "============================================================\n";

        const auto mesh_generation_start{std::chrono::steady_clock::now()};

        cfd::RawMeshData raw_mesh{cfd::generate_mesh(geometry, options)};

        const auto mesh_generation_end{std::chrono::steady_clock::now()};

        const auto mesh_generation_duration{
            std::chrono::duration<double, std::milli>(mesh_generation_end - mesh_generation_start)};

        std::cout << "\n[Mesh generation]\n";

        std::cout << std::fixed << std::setprecision(3) << "  Domain            : rectangle " << geometry.length
                  << " x " << geometry.height << " m\n"
                  << "  Target mesh size  : " << options.mesh_size << " m\n";

        std::cout << "  Cell type         : " << cell_type_name(options.cell_type) << '\n'
                  << "  Nodes             : " << raw_mesh.nodes.size() << '\n'
                  << "  Cells             : " << raw_mesh.cell_types.size() << '\n';

        std::cout << std::fixed << std::setprecision(2) << "  Time              : " << mesh_generation_duration.count()
                  << " ms\n";

        cfd::MeshBuildResult mesh_build_result{cfd::build_mesh(std::move(raw_mesh))};

        const cfd::MeshStatistics mesh_statistics{cfd::compute_mesh_statistics(mesh_build_result.mesh)};

        cfd::write_mesh_report(std::cout, mesh_build_result.mesh, mesh_statistics, mesh_build_result.timings);

        const std::filesystem::path output_directory{"results"};
        const std::filesystem::path mesh_output_file{output_directory / "mesh.vtu"};

        std::filesystem::create_directories(output_directory);

        cfd::write_vtu(mesh_build_result.mesh, mesh_output_file);

        std::cout << "\n[Summary]\n"
                  << "  Mesh              : " << mesh_build_result.mesh.node_count() << " nodes | "
                  << mesh_build_result.mesh.cell_count() << " cells | " << mesh_build_result.mesh.face_count()
                  << " faces\n";

        std::cout << "\n[Output]\n"
                  << "  Mesh              : " << mesh_output_file.string() << '\n';

        std::cout << "\n============================================================\n"
                  << " Mesh preprocessing complete\n"
                  << "============================================================\n";
    }
    catch (const std::exception &error)
    {
        std::cerr << "\n[Error]\n"
                  << "  " << error.what() << '\n';

        return 1;
    }

    return 0;
}