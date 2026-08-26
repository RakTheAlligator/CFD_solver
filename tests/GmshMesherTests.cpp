#include "cfd/meshing/GmshMesher.hpp"
#include "cfd/meshing/RawMeshValidation.hpp"

#include "support/TestUtils.hpp"

#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>

namespace
{

using cfd::test::require;
using cfd::test::require_throws_with_message;

void require_invalid_meshing_input(const cfd::RectangleGeometry &geometry, const cfd::MeshGenerationOptions &options,
                                   const std::string_view expected_message, const std::string &failure_message)
{
    require_throws_with_message<std::invalid_argument>(
        [&geometry, &options]() { static_cast<void>(cfd::generate_mesh(geometry, options)); }, expected_message,
        failure_message);
}

void test_rejects_non_finite_rectangle_length()
{
    const cfd::RectangleGeometry geometry{
        .length = std::numeric_limits<double>::quiet_NaN(),
        .height = 1.0,
    };

    const cfd::MeshGenerationOptions options{
        .mesh_size = 0.2,
        .cell_type = cfd::CellType::Triangle,
    };

    require_invalid_meshing_input(geometry, options, "Rectangle length must be finite and positive.",
                                  "Gmsh mesher accepted a non-finite rectangle length.");
}

void test_rejects_infinite_rectangle_height()
{
    const cfd::RectangleGeometry geometry{
        .length = 5.0,
        .height = std::numeric_limits<double>::infinity(),
    };

    const cfd::MeshGenerationOptions options{
        .mesh_size = 0.2,
        .cell_type = cfd::CellType::Triangle,
    };

    require_invalid_meshing_input(geometry, options, "Rectangle height must be finite and positive.",
                                  "Gmsh mesher accepted an infinite rectangle height.");
}

void test_rejects_non_finite_mesh_size()
{
    const cfd::RectangleGeometry geometry{
        .length = 5.0,
        .height = 1.0,
    };

    const cfd::MeshGenerationOptions options{
        .mesh_size = std::numeric_limits<double>::quiet_NaN(),
        .cell_type = cfd::CellType::Triangle,
    };

    require_invalid_meshing_input(geometry, options, "Mesh size must be finite and positive.",
                                  "Gmsh mesher accepted a non-finite mesh size.");
}

void test_rejects_non_positive_dimensions()
{
    const cfd::MeshGenerationOptions options{
        .mesh_size = 0.2,
        .cell_type = cfd::CellType::Triangle,
    };

    require_invalid_meshing_input({.length = 0.0, .height = 1.0}, options,
                                  "Rectangle length must be finite and positive.",
                                  "Gmsh mesher accepted a zero rectangle length.");

    require_invalid_meshing_input({.length = 5.0, .height = -1.0}, options,
                                  "Rectangle height must be finite and positive.",
                                  "Gmsh mesher accepted a negative rectangle height.");
}

void test_rejects_non_positive_mesh_size()
{
    const cfd::RectangleGeometry geometry{
        .length = 5.0,
        .height = 1.0,
    };

    require_invalid_meshing_input(geometry,
                                  {
                                      .mesh_size = 0.0,
                                      .cell_type = cfd::CellType::Triangle,
                                  },
                                  "Mesh size must be finite and positive.", "Gmsh mesher accepted a zero mesh size.");
}

void test_generates_triangular_rectangle()
{
    const cfd::RectangleGeometry geometry{
        .length = 1.0,
        .height = 1.0,
    };

    const cfd::MeshGenerationOptions options{
        .mesh_size = 0.2,
        .cell_type = cfd::CellType::Triangle,
    };

    const cfd::RawMeshData raw_mesh{cfd::generate_mesh(geometry, options)};

    cfd::validate_raw_mesh(raw_mesh);

    require(!raw_mesh.cell_types.empty(), "Triangular Gmsh mesh contains no cells.");

    for (cfd::Index cell_id = 0; cell_id < raw_mesh.cell_types.size(); ++cell_id)
    {
        require(raw_mesh.cell_types[cell_id] == cfd::CellType::Triangle,
                "Triangular Gmsh mesh contains a non-triangular cell.");

        const cfd::Index node_count{raw_mesh.cell_node_offsets[cell_id + 1] - raw_mesh.cell_node_offsets[cell_id]};

        require(node_count == 3, "Triangular Gmsh cell does not contain three nodes.");
    }
}

void test_generates_quadrilateral_rectangle()
{
    const cfd::RectangleGeometry geometry{
        .length = 1.0,
        .height = 1.0,
    };

    const cfd::MeshGenerationOptions options{
        .mesh_size = 0.2,
        .cell_type = cfd::CellType::Quadrilateral,
    };

    const cfd::RawMeshData raw_mesh{cfd::generate_mesh(geometry, options)};

    cfd::validate_raw_mesh(raw_mesh);

    require(!raw_mesh.cell_types.empty(), "Quadrilateral Gmsh mesh contains no cells.");

    for (cfd::Index cell_id = 0; cell_id < raw_mesh.cell_types.size(); ++cell_id)
    {
        require(raw_mesh.cell_types[cell_id] == cfd::CellType::Quadrilateral,
                "Quadrilateral Gmsh mesh contains a non-quadrilateral cell.");

        const cfd::Index node_count{raw_mesh.cell_node_offsets[cell_id + 1] - raw_mesh.cell_node_offsets[cell_id]};

        require(node_count == 4, "Quadrilateral Gmsh cell does not contain four nodes.");
    }
}

} // namespace

int main()
{
    int failure_count{};

    failure_count +=
        cfd::test::run_test("reject non-finite rectangle length", test_rejects_non_finite_rectangle_length);

    failure_count += cfd::test::run_test("reject infinite rectangle height", test_rejects_infinite_rectangle_height);

    failure_count += cfd::test::run_test("reject non-finite mesh size", test_rejects_non_finite_mesh_size);

    failure_count += cfd::test::run_test("reject non-positive dimensions", test_rejects_non_positive_dimensions);

    failure_count += cfd::test::run_test("reject non-positive mesh size", test_rejects_non_positive_mesh_size);

    failure_count += cfd::test::run_test("generate triangular rectangle", test_generates_triangular_rectangle);

    failure_count += cfd::test::run_test("generate quadrilateral rectangle", test_generates_quadrilateral_rectangle);

    return cfd::test::finish_tests(failure_count, "Gmsh mesher");
}