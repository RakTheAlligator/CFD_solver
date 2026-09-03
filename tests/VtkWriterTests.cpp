#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVectorField.hpp"
#include "cfd/io/VtkWriter.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/meshing/RawMeshData.hpp"

#include "support/MeshFixtures.hpp"
#include "support/TestUtils.hpp"

#include <array>
#include <filesystem>
#include <fstream>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>

namespace
{

using cfd::test::make_single_quadrilateral_raw_mesh;
using cfd::test::make_single_triangle_raw_mesh;
using cfd::test::make_two_triangle_raw_mesh;
using cfd::test::read_text_file;
using cfd::test::require;
using cfd::test::require_contains;
using cfd::test::require_throws;

class TemporaryDirectory
{
  public:
    TemporaryDirectory()
    {
        constexpr std::size_t maximum_attempt_count{64};

        const std::filesystem::path temporary_root{std::filesystem::temp_directory_path()};
        std::random_device random_source;
        std::uniform_int_distribution<unsigned long long> token_distribution;
        std::error_code last_error;

        for (std::size_t attempt = 0; attempt < maximum_attempt_count; ++attempt)
        {
            const std::filesystem::path candidate{
                temporary_root / ("cfd_vtk_writer_tests_" + std::to_string(token_distribution(random_source)))};
            std::error_code creation_error;

            // create_directory() is the atomic uniqueness check; the random
            // token only makes collisions unlikely before this check.
            if (std::filesystem::create_directory(candidate, creation_error))
            {
                path_ = candidate;
                return;
            }

            last_error = creation_error;
        }

        std::string message{"Unable to create a unique temporary directory for VtkWriter tests."};
        if (last_error)
        {
            message += " Last filesystem error: " + last_error.message();
        }
        throw std::runtime_error(message);
    }

    ~TemporaryDirectory() noexcept
    {
        std::error_code cleanup_error;
        std::filesystem::remove_all(path_, cleanup_error);
    }

    TemporaryDirectory(const TemporaryDirectory &) = delete;
    TemporaryDirectory &operator=(const TemporaryDirectory &) = delete;
    TemporaryDirectory(TemporaryDirectory &&) = delete;
    TemporaryDirectory &operator=(TemporaryDirectory &&) = delete;

    [[nodiscard]]
    const std::filesystem::path &path() const noexcept
    {
        return path_;
    }

  private:
    std::filesystem::path path_;
};

[[nodiscard]]
std::string export_mesh_and_read(cfd::RawMeshData raw_mesh, const std::string_view file_name,
                                 const cfd::VtkCellData *const cell_data = nullptr)
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(std::move(raw_mesh))};
    const cfd::Mesh &mesh{build_result.mesh};

    const TemporaryDirectory temporary_directory;
    const std::filesystem::path file_path{temporary_directory.path() / file_name};

    // Prepopulate the destination to verify that write_vtu() replaces existing
    // file content rather than appending to it.
    {
        std::ofstream existing_file{file_path};
        existing_file << "stale test content";
    }

    if (cell_data == nullptr)
    {
        cfd::write_vtu(mesh, file_path);
    }
    else
    {
        cfd::write_vtu(mesh, file_path, *cell_data);
    }

    require(std::filesystem::exists(file_path), "VTU writer did not create the output file.");
    require(std::filesystem::file_size(file_path) > 0, "VTU writer created an empty output file.");

    std::string file_content{read_text_file(file_path)};

    require(file_content.find("stale test content") == std::string::npos,
            "VTU writer did not replace existing file content.");

    return file_content;
}

void require_cell_data_rejected(const cfd::VtkCellData &cell_data, const std::string_view file_name,
                                const std::string &message)
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_triangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const TemporaryDirectory temporary_directory;
    const std::filesystem::path file_path{temporary_directory.path() / file_name};

    require_throws<std::invalid_argument>(
        [&mesh, &file_path, &cell_data]() { cfd::write_vtu(mesh, file_path, cell_data); }, message);

    require(!std::filesystem::exists(file_path), "Invalid VTU cell data created an output file.");
}

void require_common_vtu_content(const std::string &file_content)
{
    require_contains(file_content, "<VTKFile type=\"UnstructuredGrid\"",
                     "VTU output does not declare an UnstructuredGrid.");

    require_contains(file_content, "Name=\"connectivity\"", "VTU output does not contain cell connectivity.");
    require_contains(file_content, "Name=\"offsets\"", "VTU output does not contain cell offsets.");
    require_contains(file_content, "Name=\"types\"", "VTU output does not contain VTK cell types.");

    require_contains(file_content, "Name=\"cell_id\"", "VTU output does not contain cell IDs.");
    require_contains(file_content, "Name=\"cell_area\"", "VTU output does not contain cell areas.");
    require_contains(file_content, "Name=\"cell_quality\"", "VTU output does not contain cell qualities.");
}

void test_single_triangle_vtu_export()
{
    const std::string file_content{
        export_mesh_and_read(make_single_triangle_raw_mesh(), "cfd_single_triangle_test.vtu")};

    require_common_vtu_content(file_content);

    require_contains(file_content, "<Piece NumberOfPoints=\"3\" NumberOfCells=\"1\">",
                     "Triangle VTU output contains incorrect mesh dimensions.");

    require_contains(file_content, "Name=\"connectivity\" format=\"ascii\">\n          0 1 2\n",
                     "Triangle VTU output contains incorrect cell connectivity.");

    // Mesh offsets are {0, 3}; VTU stores only the exclusive end offset.
    require_contains(file_content, "Name=\"offsets\" format=\"ascii\">\n          3\n",
                     "Triangle VTU output contains an incorrect cell offset.");

    require_contains(file_content, "Name=\"types\" format=\"ascii\">\n          5\n",
                     "Triangle VTU output does not contain VTK_TRIANGLE type 5.");
}

void test_single_quadrilateral_vtu_export()
{
    const std::string file_content{
        export_mesh_and_read(make_single_quadrilateral_raw_mesh(), "cfd_single_quadrilateral_test.vtu")};

    require_common_vtu_content(file_content);

    require_contains(file_content, "<Piece NumberOfPoints=\"4\" NumberOfCells=\"1\">",
                     "Quadrilateral VTU output contains incorrect mesh dimensions.");

    require_contains(file_content, "Name=\"connectivity\" format=\"ascii\">\n          0 1 2 3\n",
                     "Quadrilateral VTU output contains incorrect cell connectivity.");

    // Mesh offsets are {0, 4}; VTU stores only the exclusive end offset.
    require_contains(file_content, "Name=\"offsets\" format=\"ascii\">\n          4\n",
                     "Quadrilateral VTU output contains an incorrect cell offset.");

    require_contains(file_content, "Name=\"types\" format=\"ascii\">\n          9\n",
                     "Quadrilateral VTU output does not contain VTK_QUAD type 9.");
}

void test_scalar_cell_field_vtu_export()
{
    cfd::CellScalarField phi{2};
    phi[0] = 1.25;
    phi[1] = -2.5;

    const std::array scalar_fields{
        cfd::VtkCellScalarData{"phi", phi.values()},
    };
    const cfd::VtkCellData cell_data{.scalars = scalar_fields};

    const std::string file_content{
        export_mesh_and_read(make_two_triangle_raw_mesh(), "cfd_scalar_cell_field_test.vtu", &cell_data)};

    require_contains(file_content, "<CellData>\n", "Scalar field was not written under CellData.");
    require_contains(file_content, "Name=\"phi\" format=\"ascii\">\n          1.25 -2.5\n",
                     "VTU output contains incorrect scalar field values.");
}

void test_vector_cell_field_vtu_export()
{
    cfd::CellVectorField gradient{2};
    gradient[0] = {1.0, 2.0};
    gradient[1] = {-3.5, 4.25};

    const std::array vector_fields{
        cfd::VtkCellVectorData{"grad_phi", gradient.values()},
    };
    const cfd::VtkCellData cell_data{.vectors = vector_fields};

    const std::string file_content{
        export_mesh_and_read(make_two_triangle_raw_mesh(), "cfd_vector_cell_field_test.vtu", &cell_data)};

    require_contains(file_content,
                     "Name=\"grad_phi\" NumberOfComponents=\"3\" format=\"ascii\">\n"
                     "          1 2 0\n"
                     "          -3.5 4.25 0\n",
                     "VTU output contains incorrect three-component Vector2 field values.");
}

void test_component_vector_cell_field_vtu_export()
{
    cfd::CellScalarField u{2};
    u[0] = 1.0;
    u[1] = -2.0;

    cfd::CellScalarField v{2};
    v[0] = 3.5;
    v[1] = 4.0;

    const std::array component_vectors{
        cfd::VtkCellVectorComponentData{"velocity", u.values(), v.values()},
    };
    const cfd::VtkCellData cell_data{.component_vectors = component_vectors};

    const std::string file_content{
        export_mesh_and_read(make_two_triangle_raw_mesh(), "cfd_component_vector_cell_field_test.vtu", &cell_data)};

    require_contains(file_content,
                     "Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n"
                     "          1 3.5 0\n"
                     "          -2 4 0\n",
                     "VTU output contains incorrect separate-component vector values.");
}

void test_multiple_user_cell_fields_coexist()
{
    cfd::CellScalarField phi{2};
    phi[0] = 4.0;
    phi[1] = 5.0;

    cfd::CellVectorField gradient{2};
    gradient[0] = {2.0, 3.0};
    gradient[1] = {2.0, 3.0};

    const std::array scalar_fields{
        cfd::VtkCellScalarData{"phi", phi.values()},
    };
    const std::array vector_fields{
        cfd::VtkCellVectorData{"grad_phi", gradient.values()},
    };
    const cfd::VtkCellData cell_data{
        .scalars = scalar_fields,
        .vectors = vector_fields,
    };

    const std::string file_content{
        export_mesh_and_read(make_two_triangle_raw_mesh(), "cfd_multiple_cell_fields_test.vtu", &cell_data)};

    require_contains(file_content, "Name=\"phi\"", "VTU output omitted a scalar field from a multi-field export.");
    require_contains(file_content, "Name=\"grad_phi\"", "VTU output omitted a vector field from a multi-field export.");
}

void test_rejects_scalar_field_cardinality_mismatch()
{
    const cfd::CellScalarField scalar_field{1};
    const std::array scalar_fields{
        cfd::VtkCellScalarData{"phi", scalar_field.values()},
    };
    const cfd::VtkCellData cell_data{.scalars = scalar_fields};

    require_cell_data_rejected(cell_data, "cfd_invalid_scalar_cardinality_test.vtu",
                               "VTU writer accepted an incorrect scalar field cardinality.");
}

void test_rejects_vector_field_cardinality_mismatch()
{
    const cfd::CellVectorField vector_field{1};
    const std::array vector_fields{
        cfd::VtkCellVectorData{"gradient", vector_field.values()},
    };
    const cfd::VtkCellData cell_data{.vectors = vector_fields};

    require_cell_data_rejected(cell_data, "cfd_invalid_vector_cardinality_test.vtu",
                               "VTU writer accepted an incorrect Vector2 field cardinality.");
}

void test_rejects_component_vector_cardinality_mismatch()
{
    const cfd::CellScalarField one_value{1};
    const cfd::CellScalarField two_values{2};

    const std::array wrong_x_components{
        cfd::VtkCellVectorComponentData{"velocity", one_value.values(), two_values.values()},
    };
    const cfd::VtkCellData wrong_x_data{.component_vectors = wrong_x_components};

    require_cell_data_rejected(wrong_x_data, "cfd_invalid_x_component_cardinality_test.vtu",
                               "VTU writer accepted an incorrect x-component cardinality.");

    const std::array wrong_y_components{
        cfd::VtkCellVectorComponentData{"velocity", two_values.values(), one_value.values()},
    };
    const cfd::VtkCellData wrong_y_data{.component_vectors = wrong_y_components};

    require_cell_data_rejected(wrong_y_data, "cfd_invalid_y_component_cardinality_test.vtu",
                               "VTU writer accepted an incorrect y-component cardinality.");
}

void test_rejects_empty_cell_field_name()
{
    const cfd::CellScalarField scalar_field{2};
    const std::array scalar_fields{
        cfd::VtkCellScalarData{"", scalar_field.values()},
    };
    const cfd::VtkCellData cell_data{.scalars = scalar_fields};

    require_cell_data_rejected(cell_data, "cfd_empty_field_name_test.vtu",
                               "VTU writer accepted an empty CellData field name.");
}

void test_rejects_duplicate_cell_field_names()
{
    const cfd::CellScalarField scalar_field{2};
    const cfd::CellVectorField vector_field{2};

    const std::array scalar_fields{
        cfd::VtkCellScalarData{"solution", scalar_field.values()},
    };
    const std::array vector_fields{
        cfd::VtkCellVectorData{"solution", vector_field.values()},
    };
    const cfd::VtkCellData cell_data{
        .scalars = scalar_fields,
        .vectors = vector_fields,
    };

    require_cell_data_rejected(cell_data, "cfd_duplicate_field_name_test.vtu",
                               "VTU writer accepted duplicate CellData field names.");
}

void test_rejects_intrinsic_cell_field_name_collision()
{
    const cfd::CellScalarField scalar_field{2};
    const std::array scalar_fields{
        cfd::VtkCellScalarData{"cell_area", scalar_field.values()},
    };
    const cfd::VtkCellData cell_data{.scalars = scalar_fields};

    require_cell_data_rejected(cell_data, "cfd_intrinsic_field_name_collision_test.vtu",
                               "VTU writer accepted a name colliding with intrinsic CellData.");
}

void test_escapes_xml_special_characters_in_field_name()
{
    const cfd::CellScalarField scalar_field{2};
    const std::array scalar_fields{
        cfd::VtkCellScalarData{"phi<&\"'>", scalar_field.values()},
    };
    const cfd::VtkCellData cell_data{.scalars = scalar_fields};

    const std::string file_content{
        export_mesh_and_read(make_two_triangle_raw_mesh(), "cfd_xml_field_name_test.vtu", &cell_data)};

    require_contains(file_content, "Name=\"phi&lt;&amp;&quot;&apos;&gt;\"",
                     "VTU writer did not escape XML-special characters in a field name.");
}

} // namespace

int main()
{
    int failure_count{};

    failure_count += cfd::test::run_test("single triangle VTU export", test_single_triangle_vtu_export);
    failure_count += cfd::test::run_test("single quadrilateral VTU export", test_single_quadrilateral_vtu_export);
    failure_count += cfd::test::run_test("scalar cell field VTU export", test_scalar_cell_field_vtu_export);
    failure_count += cfd::test::run_test("vector cell field VTU export", test_vector_cell_field_vtu_export);
    failure_count +=
        cfd::test::run_test("component-vector cell field VTU export", test_component_vector_cell_field_vtu_export);
    failure_count +=
        cfd::test::run_test("multiple user cell fields VTU export", test_multiple_user_cell_fields_coexist);
    failure_count +=
        cfd::test::run_test("reject scalar field cardinality mismatch", test_rejects_scalar_field_cardinality_mismatch);
    failure_count +=
        cfd::test::run_test("reject vector field cardinality mismatch", test_rejects_vector_field_cardinality_mismatch);
    failure_count += cfd::test::run_test("reject component-vector cardinality mismatch",
                                         test_rejects_component_vector_cardinality_mismatch);
    failure_count += cfd::test::run_test("reject empty cell field name", test_rejects_empty_cell_field_name);
    failure_count += cfd::test::run_test("reject duplicate cell field names", test_rejects_duplicate_cell_field_names);
    failure_count += cfd::test::run_test("reject intrinsic cell field name collision",
                                         test_rejects_intrinsic_cell_field_name_collision);
    failure_count +=
        cfd::test::run_test("escape XML-special field name", test_escapes_xml_special_characters_in_field_name);

    return cfd::test::finish_tests(failure_count, "vtk writer");
}
