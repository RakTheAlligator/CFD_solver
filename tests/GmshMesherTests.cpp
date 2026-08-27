#include "cfd/mesh/Boundary.hpp"
#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/mesh/Types.hpp"
#include "cfd/meshing/GmshMesher.hpp"
#include "cfd/meshing/RawMeshData.hpp"
#include "cfd/meshing/RawMeshValidation.hpp"
#include "cfd/meshing/RectangleGeometry.hpp"

#include "support/TestUtils.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace
{

using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws_with_message;

[[nodiscard]]
cfd::BoundaryId find_boundary_id(const cfd::Mesh &mesh, const std::string_view name)
{
    for (const cfd::BoundaryGroup &group : mesh.boundary_groups())
    {
        if (std::string_view{group.name} == name)
        {
            return group.id;
        }
    }

    return cfd::invalid_boundary_id;
}

[[nodiscard]]
double compute_boundary_length(const cfd::Mesh &mesh, const cfd::BoundaryId boundary_id)
{
    double total_length{};

    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (mesh.face_boundary_ids()[face_id] == boundary_id)
        {
            total_length += mesh.face_lengths()[face_id];
        }
    }

    return total_length;
}

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

    require_invalid_meshing_input(geometry,
                                  {
                                      .mesh_size = -0.2,
                                      .cell_type = cfd::CellType::Triangle,
                                  },
                                  "Mesh size must be finite and positive.",
                                  "Gmsh mesher accepted a negative mesh size.");
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

void test_builds_quadrilateral_rectangle_end_to_end()
{
    constexpr double length{5.0};
    constexpr double height{1.0};
    constexpr double mesh_size{0.2};
    constexpr double tolerance{1.0e-10};

    const cfd::RectangleGeometry geometry{
        .length = length,
        .height = height,
    };

    const cfd::MeshGenerationOptions options{
        .mesh_size = mesh_size,
        .cell_type = cfd::CellType::Quadrilateral,
    };

    cfd::RawMeshData raw_mesh{cfd::generate_mesh(geometry, options)};
    cfd::MeshBuildResult build_result{cfd::build_mesh(std::move(raw_mesh))};

    const cfd::Mesh &mesh{build_result.mesh};

    require(mesh.node_count() > 0, "Quadrilateral end-to-end mesh contains no nodes.");
    require(mesh.cell_count() > 0, "Quadrilateral end-to-end mesh contains no cells.");
    require(mesh.face_count() > 0, "Quadrilateral end-to-end mesh contains no faces.");

    require(mesh.cell_faces().size() == 4 * mesh.cell_count(),
            "Quadrilateral mesh has an incorrect number of cell-face incidences.");

    double total_area{};

    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        require(mesh.cell_types()[cell_id] == cfd::CellType::Quadrilateral,
                "Quadrilateral end-to-end mesh contains a non-quadrilateral cell.");

        const double area{mesh.cell_areas()[cell_id]};

        require(std::isfinite(area) && area > 0.0, "Quadrilateral end-to-end mesh contains an invalid cell area.");

        total_area += area;

        const double quality{mesh.cell_qualities()[cell_id]};

        require(std::isfinite(quality) && quality > 0.0,
                "Quadrilateral end-to-end mesh contains an invalid cell quality.");
        require(quality <= 1.0 + tolerance, "Quadrilateral end-to-end mesh contains a cell quality greater than 1.");
    }

    require_near(total_area, length * height, tolerance, "Quadrilateral end-to-end mesh has an incorrect total area.");

    cfd::Index internal_face_count{};
    cfd::Index boundary_face_count{};

    for (const cfd::FaceAdjacency &adjacency : mesh.face_adjacencies())
    {
        if (adjacency.is_boundary())
        {
            ++boundary_face_count;
        }
        else
        {
            ++internal_face_count;
        }
    }

    // Each quadrilateral contributes four cell-face incidences. Internal faces
    // are counted by two cells and boundary faces by one.
    require(4 * mesh.cell_count() == 2 * internal_face_count + boundary_face_count,
            "Quadrilateral end-to-end mesh violates the face-incidence identity.");

    const cfd::BoundaryId inlet_id{find_boundary_id(mesh, "inlet")};
    const cfd::BoundaryId wall_id{find_boundary_id(mesh, "wall")};
    const cfd::BoundaryId outlet_id{find_boundary_id(mesh, "outlet")};

    require(inlet_id != cfd::invalid_boundary_id, "Quadrilateral mesh does not contain an inlet boundary.");
    require(wall_id != cfd::invalid_boundary_id, "Quadrilateral mesh does not contain a wall boundary.");
    require(outlet_id != cfd::invalid_boundary_id, "Quadrilateral mesh does not contain an outlet boundary.");

    require_near(compute_boundary_length(mesh, inlet_id), height, tolerance,
                 "Quadrilateral inlet boundary length is incorrect.");
    require_near(compute_boundary_length(mesh, outlet_id), height, tolerance,
                 "Quadrilateral outlet boundary length is incorrect.");
    require_near(compute_boundary_length(mesh, wall_id), 2.0 * length, tolerance,
                 "Quadrilateral wall boundary length is incorrect.");
}

void test_rectangle_boundary_groups()
{
    constexpr double length{5.0};
    constexpr double height{1.0};
    constexpr double mesh_size{0.2};
    constexpr double tolerance{1.0e-10};

    const cfd::RectangleGeometry geometry{
        .length = length,
        .height = height,
    };

    const cfd::MeshGenerationOptions options{
        .mesh_size = mesh_size,
        .cell_type = cfd::CellType::Triangle,
    };

    cfd::RawMeshData raw_mesh{cfd::generate_mesh(geometry, options)};
    cfd::MeshBuildResult build_result{cfd::build_mesh(std::move(raw_mesh))};

    const cfd::Mesh &mesh{build_result.mesh};

    const cfd::BoundaryId inlet_id{find_boundary_id(mesh, "inlet")};
    const cfd::BoundaryId wall_id{find_boundary_id(mesh, "wall")};
    const cfd::BoundaryId outlet_id{find_boundary_id(mesh, "outlet")};

    require(inlet_id != cfd::invalid_boundary_id, "Rectangle mesh does not contain an inlet boundary.");
    require(wall_id != cfd::invalid_boundary_id, "Rectangle mesh does not contain a wall boundary.");
    require(outlet_id != cfd::invalid_boundary_id, "Rectangle mesh does not contain an outlet boundary.");

    require_near(compute_boundary_length(mesh, inlet_id), height, tolerance, "Inlet boundary length is incorrect.");
    require_near(compute_boundary_length(mesh, outlet_id), height, tolerance, "Outlet boundary length is incorrect.");
    require_near(compute_boundary_length(mesh, wall_id), 2.0 * length, tolerance, "Wall boundary length is incorrect.");
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
    failure_count +=
        cfd::test::run_test("build quadrilateral rectangle end-to-end", test_builds_quadrilateral_rectangle_end_to_end);
    failure_count += cfd::test::run_test("rectangle boundary groups", test_rectangle_boundary_groups);

    return cfd::test::finish_tests(failure_count, "Gmsh mesher");
}