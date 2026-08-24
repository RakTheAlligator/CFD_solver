#include "cfd/io/VtkWriter.hpp"

#include "cfd/mesh/Mesh.hpp"

#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>

namespace cfd
{
namespace
{

[[nodiscard]]
std::uint8_t vtk_cell_type(const CellType cell_type)
{
    // VTK cell-type identifiers:
    // 5 = VTK_TRIANGLE
    // 9 = VTK_QUAD
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
        // VTK coordinates are three-dimensional.
        // The solver is 2D, so z = 0.
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
    const auto offsets{mesh.cell_node_offsets()};

    for (Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        const Index begin{offsets[cell_id]};
        const Index end{offsets[cell_id + 1]};

        output << "          ";

        for (Index local_index = begin; local_index < end; ++local_index)
        {
            output << vtk_index(cell_nodes[local_index]);

            if (local_index + 1 < end)
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

    const auto offsets{mesh.cell_node_offsets()};

    // Our compressed connectivity begins with offset 0:
    //
    //   {0, 3, 6, 9}
    //
    // VTK expects the end offset of each cell:
    //
    //   {3, 6, 9}
    for (Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        output << vtk_index(offsets[cell_id + 1]);

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
        // std::uint8_t is generally an alias of unsigned char.
        // Cast to unsigned int so operator<< writes a number rather than
        // a character.
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

void write_cell_data(std::ofstream &output, const Mesh &mesh)
{
    output << "      <CellData>\n";

    write_cell_ids(output, mesh);
    write_cell_areas(output, mesh);
    write_cell_qualities(output, mesh);

    output << "      </CellData>\n";
}

} // namespace

void write_vtu(const Mesh &mesh, const std::filesystem::path &file_path)
{
    std::ofstream output{file_path};

    if (!output)
    {
        throw std::runtime_error("Unable to open VTU output file: " + file_path.string());
    }

    // Keep enough digits to reconstruct stored double values accurately.
    // This is intentionally different from the concise terminal formatting.
    output << std::setprecision(std::numeric_limits<double>::max_digits10);

    output << "<?xml version=\"1.0\"?>\n";
    output << "<VTKFile type=\"UnstructuredGrid\" "
              "version=\"0.1\" byte_order=\"LittleEndian\">\n";
    output << "  <UnstructuredGrid>\n";
    output << "    <Piece NumberOfPoints=\"" << mesh.node_count() << "\" NumberOfCells=\"" << mesh.cell_count()
           << "\">\n";

    write_points(output, mesh);
    write_cells(output, mesh);
    write_cell_data(output, mesh);

    output << "    </Piece>\n";
    output << "  </UnstructuredGrid>\n";
    output << "</VTKFile>\n";

    if (!output)
    {
        throw std::runtime_error("Error while writing VTU output file: " + file_path.string());
    }
}

} // namespace cfd