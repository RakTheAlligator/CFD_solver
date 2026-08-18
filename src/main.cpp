#include "cfd/GmshMesher.hpp"
#include "cfd/MeshBuilder.hpp"

#include <exception>
#include <iostream>
#include <utility>

int main()
{
    try {
        const cfd::RectangleGeometry geometry{
            .length = 5.0,
            .height = 1.0
        };

        const cfd::MeshGenerationOptions options{
            .mesh_size = 0.2,
            .cell_type = cfd::CellType::Triangle
        };

        cfd::RawMeshData raw_mesh =
            cfd::generate_mesh(geometry, options);

        cfd::Mesh mesh =
            cfd::build_mesh(std::move(raw_mesh));

        std::cout
            << "Number of nodes: "
            << mesh.node_count() << '\n';

        std::cout
            << "Number of cells: "
            << mesh.cell_count() << '\n';

        std::cout
            << "Number of faces: "
            << mesh.face_count() << '\n';
    }
    catch (const std::exception& error) {
        std::cerr
            << "Error: "
            << error.what()
            << '\n';

        return 1;
    }

    return 0;
}