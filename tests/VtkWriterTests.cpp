#include "cfd/io/VtkWriter.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"

#include "support/MeshFixtures.hpp"
#include "support/TestUtils.hpp"

#include <filesystem>
#include <string>
#include <utility>

namespace
{

using cfd::test::make_single_triangle_raw_mesh;
using cfd::test::read_text_file;
using cfd::test::require;
using cfd::test::require_contains;

void test_single_triangle_vtu_export()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};
    cfd::MeshBuildResult build_result{cfd::build_mesh(std::move(raw_mesh))};

    const cfd::Mesh &mesh{build_result.mesh};

    const std::filesystem::path file_path{std::filesystem::temp_directory_path() / "cfd_single_triangle_test.vtu"};

    std::filesystem::remove(file_path);

    cfd::write_vtu(mesh, file_path);

    require(std::filesystem::exists(file_path), "VTU writer did not create the output file.");

    require(std::filesystem::file_size(file_path) > 0, "VTU writer created an empty output file.");

    const std::string file_content{read_text_file(file_path)};

    require_contains(file_content, "<VTKFile type=\"UnstructuredGrid\"",
                     "VTU output does not declare an UnstructuredGrid.");

    require_contains(file_content, "<Piece NumberOfPoints=\"3\" NumberOfCells=\"1\">",
                     "VTU output contains incorrect mesh dimensions.");

    require_contains(file_content, "Name=\"connectivity\"", "VTU output does not contain cell connectivity.");

    require_contains(file_content, "Name=\"offsets\"", "VTU output does not contain cell offsets.");

    require_contains(file_content, "Name=\"types\"", "VTU output does not contain VTK cell types.");

    require_contains(file_content, "Name=\"cell_id\"", "VTU output does not contain cell IDs.");

    require_contains(file_content, "Name=\"cell_area\"", "VTU output does not contain cell areas.");

    require_contains(file_content, "Name=\"cell_quality\"", "VTU output does not contain cell qualities.");

    std::filesystem::remove(file_path);
}

} // namespace

int main()
{
    int failure_count{};

    failure_count += cfd::test::run_test("single triangle VTU export", test_single_triangle_vtu_export);

    return cfd::test::finish_tests(failure_count, "vtk writer");
}