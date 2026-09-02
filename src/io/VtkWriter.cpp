#include "cfd/io/VtkWriter.hpp"

#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/Node.hpp"
#include "cfd/mesh/Types.hpp"

#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>

namespace cfd
{

namespace
{

constexpr std::string_view cell_id_name{"cell_id"};
constexpr std::string_view cell_area_name{"cell_area"};
constexpr std::string_view cell_quality_name{"cell_quality"};

void validate_field_name(const std::string_view name, std::unordered_set<std::string_view> &field_names)
{
    if (name.empty())
    {
        throw std::invalid_argument("VTU CellData field names must not be empty.");
    }

    if (!field_names.insert(name).second)
    {
        throw std::invalid_argument("VTU CellData field name must be unique: " + std::string{name});
    }
}

void validate_cell_data(const Mesh &mesh, const VtkCellData &cell_data)
{
    std::unordered_set<std::string_view> field_names{
        cell_id_name,
        cell_area_name,
        cell_quality_name,
    };

    const Index cell_count{mesh.cell_count()};

    for (const VtkCellScalarData &field : cell_data.scalars)
    {
        validate_field_name(field.name, field_names);

        if (field.values.size() != cell_count)
        {
            throw std::invalid_argument("VTU scalar CellData field has incorrect cardinality: " +
                                        std::string{field.name});
        }
    }

    for (const VtkCellVectorData &field : cell_data.vectors)
    {
        validate_field_name(field.name, field_names);

        if (field.values.size() != cell_count)
        {
            throw std::invalid_argument("VTU vector CellData field has incorrect cardinality: " +
                                        std::string{field.name});
        }
    }

    for (const VtkCellVectorComponentData &field : cell_data.component_vectors)
    {
        validate_field_name(field.name, field_names);

        if (field.x_values.size() != cell_count)
        {
            throw std::invalid_argument("VTU component-vector x field has incorrect cardinality: " +
                                        std::string{field.name});
        }

        if (field.y_values.size() != cell_count)
        {
            throw std::invalid_argument("VTU component-vector y field has incorrect cardinality: " +
                                        std::string{field.name});
        }
    }
}

void write_xml_attribute_value(std::ofstream &output, const std::string_view value)
{
    for (const char character : value)
    {
        switch (character)
        {
        case '&':
            output << "&amp;";
            break;

        case '<':
            output << "&lt;";
            break;

        case '>':
            output << "&gt;";
            break;

        case '"':
            output << "&quot;";
            break;

        case '\'':
            output << "&apos;";
            break;

        default:
            output << character;
            break;
        }
    }
}

[[nodiscard]]
std::uint8_t vtk_cell_type(const CellType cell_type)
{
    // VTK legacy cell-type identifiers used by the XML UnstructuredGrid
    // format: 5 = VTK_TRIANGLE, 9 = VTK_QUAD.
    switch (cell_type)
    {
    case CellType::Triangle:
        return 5;

    case CellType::Quadrilateral:
        return 9;

    default:
        throw std::runtime_error("Unsupported cell type for VTU export.");
    }
}

[[nodiscard]]
std::int64_t vtk_index(const Index index)
{
    constexpr auto vtk_index_max{static_cast<std::uintmax_t>(std::numeric_limits<std::int64_t>::max())};

    // Mesh Index is unsigned, whereas this VTU writer stores connectivity and
    // offsets as signed Int64. Check representability before conversion.
    if (static_cast<std::uintmax_t>(index) > vtk_index_max)
    {
        throw std::runtime_error("Mesh index is too large for VTU Int64 storage.");
    }

    return static_cast<std::int64_t>(index);
}

void write_points(std::ofstream &output, const Mesh &mesh)
{
    output << "      <Points>\n";
    output << "        <DataArray type=\"Float64\" "
              "NumberOfComponents=\"3\" format=\"ascii\">\n";

    for (const Node &node : mesh.nodes())
    {
        // VTK point coordinates are three-dimensional. Embed the solver's 2D
        // mesh in the z = 0 plane.
        output << "          " << node.x << ' ' << node.y << " 0\n";
    }

    output << "        </DataArray>\n";
    output << "      </Points>\n";
}

void write_cell_connectivity(std::ofstream &output, const Mesh &mesh)
{
    output << "        <DataArray type=\"Int64\" "
              "Name=\"connectivity\" format=\"ascii\">\n";

    const auto cell_nodes{mesh.cell_nodes()};
    const auto cell_node_offsets{mesh.cell_node_offsets()};

    for (Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        const Index cell_node_begin_offset{cell_node_offsets[cell_id]};
        const Index cell_node_end_offset{cell_node_offsets[cell_id + 1]};

        output << "          ";

        for (Index cell_node_position = cell_node_begin_offset; cell_node_position < cell_node_end_offset;
             ++cell_node_position)
        {
            output << vtk_index(cell_nodes[cell_node_position]);

            if (cell_node_position + 1 < cell_node_end_offset)
            {
                output << ' ';
            }
        }

        output << '\n';
    }

    output << "        </DataArray>\n";
}

void write_cell_offsets(std::ofstream &output, const Mesh &mesh)
{
    output << "        <DataArray type=\"Int64\" "
              "Name=\"offsets\" format=\"ascii\">\n";
    output << "          ";

    const auto cell_node_offsets{mesh.cell_node_offsets()};

    // Mesh connectivity uses CSR-like offsets including the initial zero:
    //
    //   {0, 3, 6, 9}
    //
    // VTK instead stores only the exclusive end offset of each cell:
    //
    //   {3, 6, 9}
    for (Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        output << vtk_index(cell_node_offsets[cell_id + 1]);

        if (cell_id + 1 < mesh.cell_count())
        {
            output << ' ';
        }
    }

    output << '\n';
    output << "        </DataArray>\n";
}

void write_cell_types(std::ofstream &output, const Mesh &mesh)
{
    output << "        <DataArray type=\"UInt8\" "
              "Name=\"types\" format=\"ascii\">\n";
    output << "          ";

    const auto cell_types{mesh.cell_types()};

    for (Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        // std::uint8_t is typically an alias of unsigned char. Cast before
        // streaming so the numeric VTK identifier is written, not a character.
        output << static_cast<unsigned int>(vtk_cell_type(cell_types[cell_id]));

        if (cell_id + 1 < mesh.cell_count())
        {
            output << ' ';
        }
    }

    output << '\n';
    output << "        </DataArray>\n";
}

void write_cells(std::ofstream &output, const Mesh &mesh)
{
    output << "      <Cells>\n";

    write_cell_connectivity(output, mesh);
    write_cell_offsets(output, mesh);
    write_cell_types(output, mesh);

    output << "      </Cells>\n";
}

void write_cell_ids(std::ofstream &output, const Mesh &mesh)
{
    output << "        <DataArray type=\"Int64\" "
              "Name=\"cell_id\" format=\"ascii\">\n";
    output << "          ";

    for (Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        output << vtk_index(cell_id);

        if (cell_id + 1 < mesh.cell_count())
        {
            output << ' ';
        }
    }

    output << '\n';
    output << "        </DataArray>\n";
}

void write_cell_areas(std::ofstream &output, const Mesh &mesh)
{
    output << "        <DataArray type=\"Float64\" "
              "Name=\"cell_area\" format=\"ascii\">\n";
    output << "          ";

    const auto cell_areas{mesh.cell_areas()};

    for (Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        output << cell_areas[cell_id];

        if (cell_id + 1 < mesh.cell_count())
        {
            output << ' ';
        }
    }

    output << '\n';
    output << "        </DataArray>\n";
}

void write_cell_qualities(std::ofstream &output, const Mesh &mesh)
{
    output << "        <DataArray type=\"Float64\" "
              "Name=\"cell_quality\" format=\"ascii\">\n";
    output << "          ";

    const auto cell_qualities{mesh.cell_qualities()};

    for (Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        output << cell_qualities[cell_id];

        if (cell_id + 1 < mesh.cell_count())
        {
            output << ' ';
        }
    }

    output << '\n';
    output << "        </DataArray>\n";
}

void write_cell_scalar_data(std::ofstream &output, const VtkCellScalarData &field)
{
    output << "        <DataArray type=\"Float64\" Name=\"";
    write_xml_attribute_value(output, field.name);
    output << "\" format=\"ascii\">\n";
    output << "          ";

    for (Index cell_id = 0; cell_id < field.values.size(); ++cell_id)
    {
        output << field.values[cell_id];

        if (cell_id + 1 < field.values.size())
        {
            output << ' ';
        }
    }

    output << '\n';
    output << "        </DataArray>\n";
}

void write_cell_vector_data(std::ofstream &output, const VtkCellVectorData &field)
{
    output << "        <DataArray type=\"Float64\" Name=\"";
    write_xml_attribute_value(output, field.name);
    output << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";

    for (const Vector2 &value : field.values)
    {
        output << "          " << value.x << ' ' << value.y << " 0\n";
    }

    output << "        </DataArray>\n";
}

void write_cell_vector_component_data(std::ofstream &output, const VtkCellVectorComponentData &field)
{
    output << "        <DataArray type=\"Float64\" Name=\"";
    write_xml_attribute_value(output, field.name);
    output << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";

    for (Index cell_id = 0; cell_id < field.x_values.size(); ++cell_id)
    {
        output << "          " << field.x_values[cell_id] << ' ' << field.y_values[cell_id] << " 0\n";
    }

    output << "        </DataArray>\n";
}

void write_cell_data(std::ofstream &output, const Mesh &mesh, const VtkCellData &cell_data)
{
    output << "      <CellData>\n";

    write_cell_ids(output, mesh);
    write_cell_areas(output, mesh);
    write_cell_qualities(output, mesh);

    for (const VtkCellScalarData &field : cell_data.scalars)
    {
        write_cell_scalar_data(output, field);
    }

    for (const VtkCellVectorData &field : cell_data.vectors)
    {
        write_cell_vector_data(output, field);
    }

    for (const VtkCellVectorComponentData &field : cell_data.component_vectors)
    {
        write_cell_vector_component_data(output, field);
    }

    output << "      </CellData>\n";
}

} // namespace

void write_vtu(const Mesh &mesh, const std::filesystem::path &file_path)
{
    write_vtu(mesh, file_path, VtkCellData{});
}

void write_vtu(const Mesh &mesh, const std::filesystem::path &file_path, const VtkCellData &cell_data)
{
    validate_cell_data(mesh, cell_data);

    std::ofstream output{file_path};

    if (!output)
    {
        throw std::runtime_error("Unable to open VTU output file: " + file_path.string());
    }

    // max_digits10 preserves enough decimal digits for a written double to
    // round-trip back to the same binary value.
    output << std::setprecision(std::numeric_limits<double>::max_digits10);

    // ASCII keeps this diagnostic export simple and inspectable. If VTU output
    // becomes a significant cost for large or transient simulations, binary or
    // appended VTK data should be evaluated from measurements.
    output << "<?xml version=\"1.0\"?>\n";
    output << "<VTKFile type=\"UnstructuredGrid\" "
              "version=\"0.1\" byte_order=\"LittleEndian\">\n";
    output << "  <UnstructuredGrid>\n";
    output << "    <Piece NumberOfPoints=\"" << mesh.node_count() << "\" NumberOfCells=\"" << mesh.cell_count()
           << "\">\n";

    write_points(output, mesh);
    write_cells(output, mesh);
    write_cell_data(output, mesh, cell_data);

    output << "    </Piece>\n";
    output << "  </UnstructuredGrid>\n";
    output << "</VTKFile>\n";

    output.close();

    if (!output)
    {
        throw std::runtime_error("Error while writing VTU output file: " + file_path.string());
    }
}

} // namespace cfd
