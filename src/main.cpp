#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/meshing/GmshMesher.hpp"

#include "cfd/io/VtkWriter.hpp"

#include <chrono>
#include <exception>
#include <iomanip>
#include <iostream>
#include <string_view>
#include <utility>
#include <filesystem>

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
            .cell_type = cfd::CellType::Triangle,
        };

        std::cout << "============================================================\n"
                  << " CFD Solver\n"
                  << "============================================================\n";

        const auto generation_start{std::chrono::steady_clock::now()};

        cfd::RawMeshData raw_mesh{cfd::generate_mesh(geometry, options)};

        const auto generation_end{std::chrono::steady_clock::now()};

        const auto generation_elapsed{std::chrono::duration<double, std::milli>(generation_end - generation_start)};

        std::cout << "\n[Mesh generation]\n";

        std::cout << std::fixed << std::setprecision(3) << "  Domain            : rectangle " << geometry.length
                  << " x " << geometry.height << " m\n"
                  << "  Target mesh size  : " << options.mesh_size << " m\n";

        std::cout << "  Cell type         : " << cell_type_name(options.cell_type) << '\n'
                  << "  Nodes             : " << raw_mesh.nodes.size() << '\n'
                  << "  Cells             : " << raw_mesh.cell_types.size() << '\n';

        std::cout << std::fixed << std::setprecision(2) << "  Time              : " << generation_elapsed.count()
                  << " ms\n";

        cfd::Mesh mesh{cfd::build_mesh(std::move(raw_mesh))};

        const std::filesystem::path output_directory{"results"};
        const std::filesystem::path mesh_output_file{output_directory / "mesh.vtu"};

        std::filesystem::create_directories(output_directory);

        cfd::write_vtu(mesh, mesh_output_file);

        std::cout << "\n[Summary]\n"
                << "  Mesh              : " << mesh.node_count() << " nodes | " << mesh.cell_count() << " cells | "
                << mesh.face_count() << " faces\n";

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