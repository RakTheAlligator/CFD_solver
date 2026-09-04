#include "cfd/input/OpenFOAMCaseReader.hpp"
#include "cfd/mesh/Boundary.hpp"
#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/meshing/RawMeshData.hpp"

#include "support/TestUtils.hpp"

#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>

namespace
{

using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws;
using cfd::test::require_throws_with_message;

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
            const std::filesystem::path candidate{temporary_root / ("cfd_openfoam_case_reader_tests_" +
                                                                    std::to_string(token_distribution(random_source)))};
            std::error_code creation_error;

            if (std::filesystem::create_directory(candidate, creation_error))
            {
                path_ = candidate;
                return;
            }

            last_error = creation_error;
        }

        std::string message{"Unable to create a unique temporary directory for OpenFOAMCaseReader tests."};
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

void write_text_file(const std::filesystem::path &file_path, const std::string_view content)
{
    std::ofstream output{file_path};
    output << content;

    if (!output)
    {
        throw std::runtime_error("Unable to write test input file: " + file_path.string());
    }
}

[[nodiscard]]
cfd::input::MeshInput read_mesh_input(const std::string_view content)
{
    const TemporaryDirectory temporary_directory;
    const std::filesystem::path file_path{temporary_directory.path() / "meshDict"};
    write_text_file(file_path, content);
    return cfd::input::read_mesh_dict(file_path);
}

[[nodiscard]]
cfd::input::ScalarFieldInput read_scalar_input(const std::string_view content)
{
    const TemporaryDirectory temporary_directory;
    const std::filesystem::path file_path{temporary_directory.path() / "c"};
    write_text_file(file_path, content);
    return cfd::input::read_scalar_field(file_path);
}

[[nodiscard]]
cfd::Mesh make_named_boundary_mesh()
{
    cfd::RawMeshData raw_mesh;

    raw_mesh.nodes = {
        {0.0, 0.0},
        {1.0, 0.0},
        {1.0, 1.0},
        {0.0, 1.0},
    };
    raw_mesh.cell_types = {cfd::CellType::Quadrilateral};
    raw_mesh.cell_nodes = {0, 1, 2, 3};
    raw_mesh.cell_node_offsets = {0, 4};

    raw_mesh.boundary_groups = {
        {0, "wall"},
        {1, "outlet"},
        {2, "inlet"},
    };
    raw_mesh.boundary_edges = {
        {{0, 1}, 0},
        {{1, 2}, 1},
        {{2, 3}, 0},
        {{3, 0}, 2},
    };

    return cfd::build_mesh(std::move(raw_mesh)).mesh;
}

[[nodiscard]]
const cfd::input::NamedScalarBoundaryCondition &find_named_condition(const cfd::input::ScalarFieldInput &field_input,
                                                                     const std::string_view name)
{
    for (const cfd::input::NamedScalarBoundaryCondition &condition : field_input.boundary_conditions)
    {
        if (condition.boundary_name == name)
        {
            return condition;
        }
    }

    cfd::test::fail("Parsed scalar field does not contain boundary '" + std::string(name) + "'.");
}

[[nodiscard]]
std::string make_mesh_dictionary(const std::string_view header_entries, const std::string_view geometry_entries,
                                 const std::string_view mesh_entries)
{
    std::ostringstream dictionary;
    dictionary << "FoamFile\n{\n"
               << header_entries << "\n}\ngeometry\n{\n"
               << geometry_entries << "\n}\nmesh\n{\n"
               << mesh_entries << "\n}\n";
    return dictionary.str();
}

[[nodiscard]]
std::string make_scalar_dictionary(const std::string_view header_entries, const std::string_view boundary_entries)
{
    std::ostringstream dictionary;
    dictionary << "FoamFile\n{\n"
               << header_entries << "\n}\ndimensions [0 0 0 0 0 0 0];\ninternalField uniform 0;\nboundaryField\n{\n"
               << boundary_entries << "\n}\n";
    return dictionary.str();
}

void test_reads_valid_quadrilateral_mesh_dict()
{
    constexpr std::string_view content{R"(
FoamFile
{
    version 2.0;
    format "ascii";
    class "dictionary";
    object "meshDict";
}

geometry
{
    type rectangle;
    length +5.0e0;
    height 1e0;
}

mesh
{
    cellType quadrilateral;
    size 2e-1;
}
)"};

    const cfd::input::MeshInput input{read_mesh_input(content)};

    require_near(input.geometry.length, 5.0, 0.0, "meshDict rectangle length was parsed incorrectly.");
    require_near(input.geometry.height, 1.0, 0.0, "meshDict rectangle height was parsed incorrectly.");
    require_near(input.generation_options.mesh_size, 0.2, 0.0, "meshDict mesh size was parsed incorrectly.");
    require(input.generation_options.cell_type == cfd::CellType::Quadrilateral,
            "meshDict quadrilateral cell type was parsed incorrectly.");
}

void test_reads_triangle_mesh_dict()
{
    constexpr std::string_view content{R"(
FoamFile
{
    version 2.0;
    format ascii;
    class dictionary;
    object meshDict;
}

geometry
{
    type rectangle;
    length 2;
    height 1;
}
mesh
{
    cellType triangle;
    size 0.25;
}
)"};

    const cfd::input::MeshInput input{read_mesh_input(content)};

    require(input.generation_options.cell_type == cfd::CellType::Triangle,
            "meshDict triangle cell type was parsed incorrectly.");
}

void test_accepts_reordered_foam_file_entries()
{
    const std::string content{make_mesh_dictionary(
        R"(object meshDict;
class dictionary;
version 2.0;
format ascii;)",
        R"(type rectangle;
length 5;
height 1;)",
        R"(cellType quadrilateral;
size 0.2;)")};

    const cfd::input::MeshInput input{read_mesh_input(content)};

    require(input.generation_options.cell_type == cfd::CellType::Quadrilateral,
            "Reordered FoamFile entries changed the parsed mesh input.");
}

void test_accepts_omitted_optional_version()
{
    const std::string content{make_mesh_dictionary(
        R"(format ascii;
class dictionary;
object meshDict;)",
        R"(type rectangle;
length 5;
height 1;)",
        R"(cellType triangle;
size 0.2;)")};

    const cfd::input::MeshInput input{read_mesh_input(content)};

    require(input.generation_options.cell_type == cfd::CellType::Triangle,
            "Omitting the optional FoamFile version changed the parsed mesh input.");
}

void test_accepts_optional_location()
{
    const std::string content{make_scalar_dictionary(
        R"(location "0";
object c;
format ascii;
class volScalarField;
version 2.0;)",
        R"(wall
{
    type zeroGradient;
})")};

    const cfd::input::ScalarFieldInput input{read_scalar_input(content)};

    require(input.object_name == "c", "Optional FoamFile location changed the parsed object name.");
}

void test_accepts_reordered_geometry_entries()
{
    const std::string content{make_mesh_dictionary(
        R"(version 2.0;
format ascii;
class dictionary;
object meshDict;)",
        R"(height 1;
type rectangle;
length 5;)",
        R"(cellType triangle;
size 0.2;)")};

    const cfd::input::MeshInput input{read_mesh_input(content)};

    require_near(input.geometry.length, 5.0, 0.0, "Reordered geometry length was parsed incorrectly.");
    require_near(input.geometry.height, 1.0, 0.0, "Reordered geometry height was parsed incorrectly.");
}

void test_accepts_reordered_mesh_entries()
{
    const std::string content{make_mesh_dictionary(
        R"(version 2.0;
format ascii;
class dictionary;
object meshDict;)",
        R"(type rectangle;
length 5;
height 1;)",
        R"(size 0.125;
cellType quadrilateral;)")};

    const cfd::input::MeshInput input{read_mesh_input(content)};

    require_near(input.generation_options.mesh_size, 0.125, 0.0, "Reordered mesh size was parsed incorrectly.");
    require(input.generation_options.cell_type == cfd::CellType::Quadrilateral,
            "Reordered mesh cellType was parsed incorrectly.");
}

void test_accepts_reordered_fixed_value_entries()
{
    const std::string content{make_scalar_dictionary(
        R"(version 2.0;
format ascii;
class volScalarField;
object c;)",
        R"(inlet
{
    value uniform 1.5;
    type fixedValue;
})")};

    const cfd::input::ScalarFieldInput input{read_scalar_input(content)};
    const cfd::ScalarBoundaryCondition &condition{input.boundary_conditions.front().condition};

    require(condition.type == cfd::ScalarBoundaryConditionType::Dirichlet && condition.value == 1.5,
            "Reordered fixedValue entries were mapped incorrectly.");
}

void test_accepts_reordered_fixed_gradient_entries()
{
    const std::string content{make_scalar_dictionary(
        R"(version 2.0;
format ascii;
class volScalarField;
object c;)",
        R"(outlet
{
    gradient uniform -0.75;
    type fixedGradient;
})")};

    const cfd::input::ScalarFieldInput input{read_scalar_input(content)};
    const cfd::ScalarBoundaryCondition &condition{input.boundary_conditions.front().condition};

    require(condition.type == cfd::ScalarBoundaryConditionType::Neumann && condition.value == -0.75,
            "Reordered fixedGradient entries were mapped incorrectly.");
}

void test_rejects_duplicate_supported_key()
{
    const std::string content{make_mesh_dictionary(
        R"(version 2.0;
format ascii;
format ascii;
class dictionary;
object meshDict;)",
        R"(type rectangle;
length 5;
height 1;)",
        R"(cellType triangle;
size 0.2;)")};

    require_throws<std::runtime_error>([&content]() { static_cast<void>(read_mesh_input(content)); },
                                       "meshDict reader accepted a duplicate supported key.");
}

void test_rejects_missing_required_key()
{
    const std::string content{make_mesh_dictionary(
        R"(version 2.0;
class dictionary;
object meshDict;)",
        R"(type rectangle;
length 5;
height 1;)",
        R"(cellType triangle;
size 0.2;)")};

    require_throws<std::runtime_error>([&content]() { static_cast<void>(read_mesh_input(content)); },
                                       "meshDict reader accepted a missing required key.");
}

void test_rejects_unknown_key()
{
    const std::string content{make_mesh_dictionary(
        R"(version 2.0;
format ascii;
class dictionary;
object meshDict;)",
        R"(type rectangle;
length 5;
height 1;)",
        R"(cellType triangle;
algorithm automatic;
size 0.2;)")};

    require_throws<std::runtime_error>([&content]() { static_cast<void>(read_mesh_input(content)); },
                                       "meshDict reader accepted an unknown key.");
}

void test_reads_scalar_field_and_maps_supported_conditions()
{
    constexpr std::string_view content{R"(
FoamFile
{
    version 2.0;
    format ascii;
    class volScalarField;
    object c;
}
dimensions [0 0 0 0 0 0 0];
internalField uniform -1.25e+2;
boundaryField
{
    inlet
    {
        type fixedValue;
        value uniform 1;
    }
    wall
    {
        type zeroGradient;
    }
    outlet
    {
        type fixedGradient;
        gradient uniform -2.5e-1;
    }
}
)"};

    const cfd::input::ScalarFieldInput input{read_scalar_input(content)};

    require(input.object_name == "c", "Scalar field object name was not retained.");
    require_near(input.internal_value, -125.0, 0.0, "Scalar internalField was parsed incorrectly.");
    require(input.boundary_conditions.size() == 3, "Scalar field has an incorrect boundary-condition count.");

    const cfd::ScalarBoundaryCondition &inlet{find_named_condition(input, "inlet").condition};
    const cfd::ScalarBoundaryCondition &wall{find_named_condition(input, "wall").condition};
    const cfd::ScalarBoundaryCondition &outlet{find_named_condition(input, "outlet").condition};

    require(inlet.type == cfd::ScalarBoundaryConditionType::Dirichlet && inlet.value == 1.0,
            "fixedValue was not mapped to the expected Dirichlet condition.");
    require(wall.type == cfd::ScalarBoundaryConditionType::Neumann && wall.value == 0.0,
            "zeroGradient was not mapped to the expected zero Neumann condition.");
    require(outlet.type == cfd::ScalarBoundaryConditionType::Neumann && outlet.value == -0.25,
            "fixedGradient was not mapped to the expected Neumann condition.");
}

void test_accepts_line_and_block_comments()
{
    constexpr std::string_view content{R"(
// Header comment
FoamFile
{
    version 2.0; /* inline block comment */
    format ascii;
    class volScalarField;
    object c;
}
dimensions [0 0 0 0 0 0 0];
internalField /* distribution */ uniform 0;
boundaryField
{
    wall
    {
        // Homogeneous outward-normal derivative.
        type zeroGradient;
    }
}
)"};

    const cfd::input::ScalarFieldInput input{read_scalar_input(content)};

    require(input.boundary_conditions.size() == 1, "Comments changed the parsed boundary-condition count.");
    require(input.boundary_conditions.front().condition.type == cfd::ScalarBoundaryConditionType::Neumann,
            "Comments changed the parsed boundary-condition type.");
}

void test_resolves_boundaries_independently_of_input_order()
{
    cfd::Mesh mesh{make_named_boundary_mesh()};
    const cfd::input::ScalarFieldInput field_input{
        .object_name = "c",
        .internal_value = 0.0,
        .boundary_conditions =
            {
                {"outlet", {cfd::ScalarBoundaryConditionType::Dirichlet, 0.0}},
                {"inlet", {cfd::ScalarBoundaryConditionType::Dirichlet, 1.0}},
                {"wall", {cfd::ScalarBoundaryConditionType::Neumann, 0.0}},
            },
    };

    const cfd::ScalarBoundaryConditions conditions{cfd::input::resolve_boundary_conditions(mesh, field_input)};

    require(conditions[0].type == cfd::ScalarBoundaryConditionType::Neumann,
            "Named wall condition was not assigned to the wall BoundaryId.");
    require(conditions[1].type == cfd::ScalarBoundaryConditionType::Dirichlet && conditions[1].value == 0.0,
            "Named outlet condition was not assigned to the outlet BoundaryId.");
    require(conditions[2].type == cfd::ScalarBoundaryConditionType::Dirichlet && conditions[2].value == 1.0,
            "Named inlet condition was not assigned to the inlet BoundaryId.");
}

void test_rejects_unknown_boundary_name()
{
    cfd::Mesh mesh{make_named_boundary_mesh()};
    const cfd::input::ScalarFieldInput field_input{
        .object_name = "c",
        .boundary_conditions =
            {
                {"wall", {cfd::ScalarBoundaryConditionType::Neumann, 0.0}},
                {"outlet", {cfd::ScalarBoundaryConditionType::Dirichlet, 0.0}},
                {"farfield", {cfd::ScalarBoundaryConditionType::Dirichlet, 1.0}},
            },
    };

    require_throws<std::invalid_argument>(
        [&mesh, &field_input]() { static_cast<void>(cfd::input::resolve_boundary_conditions(mesh, field_input)); },
        "Boundary resolution accepted an unknown boundary name.");
}

void test_rejects_duplicate_boundary_name()
{
    cfd::Mesh mesh{make_named_boundary_mesh()};
    const cfd::input::ScalarFieldInput field_input{
        .object_name = "c",
        .boundary_conditions =
            {
                {"wall", {cfd::ScalarBoundaryConditionType::Neumann, 0.0}},
                {"outlet", {cfd::ScalarBoundaryConditionType::Dirichlet, 0.0}},
                {"inlet", {cfd::ScalarBoundaryConditionType::Dirichlet, 1.0}},
                {"wall", {cfd::ScalarBoundaryConditionType::Dirichlet, 2.0}},
            },
    };

    require_throws<std::invalid_argument>(
        [&mesh, &field_input]() { static_cast<void>(cfd::input::resolve_boundary_conditions(mesh, field_input)); },
        "Boundary resolution accepted a duplicate boundary name.");
}

void test_rejects_missing_boundary_condition()
{
    cfd::Mesh mesh{make_named_boundary_mesh()};
    const cfd::input::ScalarFieldInput field_input{
        .object_name = "c",
        .boundary_conditions =
            {
                {"wall", {cfd::ScalarBoundaryConditionType::Neumann, 0.0}},
                {"inlet", {cfd::ScalarBoundaryConditionType::Dirichlet, 1.0}},
            },
    };

    require_throws<std::invalid_argument>(
        [&mesh, &field_input]() { static_cast<void>(cfd::input::resolve_boundary_conditions(mesh, field_input)); },
        "Boundary resolution accepted a missing boundary condition.");
}

void test_rejects_malformed_syntax_with_file_and_line()
{
    constexpr std::string_view content{R"(FoamFile
{
    version 2.0;
    format ascii;
    class volScalarField;
    object c;
}
dimensions [0 0 0 0 0 0 0];
internalField uniform 0
boundaryField {}
)"};

    const TemporaryDirectory temporary_directory;
    const std::filesystem::path file_path{temporary_directory.path() / "malformed_c"};
    write_text_file(file_path, content);

    require_throws_with_message<std::runtime_error>(
        [&file_path]() { static_cast<void>(cfd::input::read_scalar_field(file_path)); },
        file_path.string() + ":10:", "Malformed scalar syntax did not report its file and line.");
}

void test_rejects_unsupported_boundary_condition_type()
{
    constexpr std::string_view content{R"(
FoamFile
{
    version 2.0;
    format ascii;
    class volScalarField;
    object c;
}
dimensions [0 0 0 0 0 0 0];
internalField uniform 0;
boundaryField
{
    wall
    {
        type mixed;
    }
}
)"};

    require_throws<std::runtime_error>([content]() { static_cast<void>(read_scalar_input(content)); },
                                       "Scalar reader accepted an unsupported boundary-condition type.");
}

void test_rejects_nonuniform_scalar_field()
{
    constexpr std::string_view content{R"(
FoamFile
{
    version 2.0;
    format ascii;
    class volScalarField;
    object c;
}
dimensions [0 0 0 0 0 0 0];
internalField nonuniform 0;
boundaryField {}
)"};

    require_throws_with_message<std::runtime_error>([content]() { static_cast<void>(read_scalar_input(content)); },
                                                    "Nonuniform",
                                                    "Scalar reader accepted a nonuniform internal field.");
}

void test_rejects_non_positive_mesh_input()
{
    constexpr std::string_view content{R"(
FoamFile
{
    version 2.0;
    format ascii;
    class dictionary;
    object meshDict;
}
geometry
{
    type rectangle;
    length 0;
    height 1;
}
mesh
{
    cellType triangle;
    size 0.2;
}
)"};

    require_throws<std::runtime_error>([content]() { static_cast<void>(read_mesh_input(content)); },
                                       "meshDict reader accepted a non-positive rectangle dimension.");
}

} // namespace

int main()
{
    int failure_count{};

    failure_count += cfd::test::run_test("read valid quadrilateral meshDict", test_reads_valid_quadrilateral_mesh_dict);
    failure_count += cfd::test::run_test("read triangle meshDict", test_reads_triangle_mesh_dict);
    failure_count += cfd::test::run_test("accept reordered FoamFile entries", test_accepts_reordered_foam_file_entries);
    failure_count +=
        cfd::test::run_test("accept omitted optional FoamFile version", test_accepts_omitted_optional_version);
    failure_count += cfd::test::run_test("accept optional FoamFile location", test_accepts_optional_location);
    failure_count += cfd::test::run_test("accept reordered geometry entries", test_accepts_reordered_geometry_entries);
    failure_count += cfd::test::run_test("accept reordered mesh entries", test_accepts_reordered_mesh_entries);
    failure_count +=
        cfd::test::run_test("accept reordered fixedValue entries", test_accepts_reordered_fixed_value_entries);
    failure_count +=
        cfd::test::run_test("accept reordered fixedGradient entries", test_accepts_reordered_fixed_gradient_entries);
    failure_count += cfd::test::run_test("reject duplicate supported key", test_rejects_duplicate_supported_key);
    failure_count += cfd::test::run_test("reject missing required key", test_rejects_missing_required_key);
    failure_count += cfd::test::run_test("reject unknown key", test_rejects_unknown_key);
    failure_count += cfd::test::run_test("read scalar field and map supported conditions",
                                         test_reads_scalar_field_and_maps_supported_conditions);
    failure_count += cfd::test::run_test("accept line and block comments", test_accepts_line_and_block_comments);
    failure_count += cfd::test::run_test("resolve boundary names independently of input order",
                                         test_resolves_boundaries_independently_of_input_order);
    failure_count += cfd::test::run_test("reject unknown boundary name", test_rejects_unknown_boundary_name);
    failure_count += cfd::test::run_test("reject duplicate boundary name", test_rejects_duplicate_boundary_name);
    failure_count += cfd::test::run_test("reject missing boundary condition", test_rejects_missing_boundary_condition);
    failure_count +=
        cfd::test::run_test("report malformed syntax file and line", test_rejects_malformed_syntax_with_file_and_line);
    failure_count += cfd::test::run_test("reject unsupported scalar boundary condition",
                                         test_rejects_unsupported_boundary_condition_type);
    failure_count += cfd::test::run_test("reject nonuniform scalar field", test_rejects_nonuniform_scalar_field);
    failure_count += cfd::test::run_test("reject non-positive mesh input", test_rejects_non_positive_mesh_input);

    return cfd::test::finish_tests(failure_count, "OpenFOAMCaseReader");
}
