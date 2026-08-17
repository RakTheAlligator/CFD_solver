#include "cfd/GmshMesher.hpp"
#include "cfd/RawMeshValidation.hpp"

#include <exception>
#include <iostream>

int main()
{
    try {
        const cfd::RectangleGeometry geometry{      //  Define the rectangle geometry
            .length = 5.0,
            .height = 1.0
        };

        const cfd::MeshGenerationOptions options{   //  Set mesh generation options
            .mesh_size = 0.2,
            .cell_type = cfd::CellType::Triangle
        };

        const cfd::RawMeshData raw_mesh =           //  Generate the mesh using Gmsh
            cfd::generate_mesh(geometry, options);

        std::cout
            << "Number of nodes: "
            << raw_mesh.nodes.size()
            << '\n';

        std::cout
            << "Number of cells: "
            << raw_mesh.cell_types.size()
            << '\n';

        std::cout
            << "Number of boundary edges: "
            << raw_mesh.boundary_edges.size()
            << '\n';

        cfd::validate_raw_mesh(raw_mesh);  //  Validate the generated mesh
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