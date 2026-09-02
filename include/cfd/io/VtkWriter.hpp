#pragma once

#include "cfd/math/Vector2.hpp"

#include <filesystem>
#include <span>
#include <string_view>

namespace cfd
{

class Mesh;

/// Non-owning description of one cell-centered scalar array for VTU export.
///
/// The name and values must remain valid for the duration of the synchronous
/// `write_vtu()` call.
struct VtkCellScalarData
{
    std::string_view name;
    std::span<const double> values;
};

/// Non-owning description of one cell-centered Vector2 array for VTU export.
///
/// Each value is serialized as the three-component VTK vector `(x, y, 0)`.
/// The name and values must remain valid for the duration of `write_vtu()`.
struct VtkCellVectorData
{
    std::string_view name;
    std::span<const Vector2> values;
};

/// Non-owning description of a cell-centered vector stored as scalar components.
///
/// Corresponding x and y values are serialized directly as `(x, y, 0)`, without
/// constructing an interleaved vector field. The name and both component views
/// must remain valid for the duration of `write_vtu()`.
struct VtkCellVectorComponentData
{
    std::string_view name;
    std::span<const double> x_values;
    std::span<const double> y_values;
};

/// Non-owning collection of additional cell-centered arrays for VTU export.
///
/// The descriptor spans and every name/value view referenced by them must remain
/// valid for the duration of the synchronous `write_vtu()` call.
struct VtkCellData
{
    std::span<const VtkCellScalarData> scalars{};
    std::span<const VtkCellVectorData> vectors{};
    std::span<const VtkCellVectorComponentData> component_vectors{};
};

/// Writes a mesh to a VTK XML UnstructuredGrid (.vtu) file in ASCII format.
///
/// Existing content at `file_path` is replaced. The export includes mesh
/// connectivity together with cell IDs, areas, and quality values.
///
/// @param mesh Validated mesh to export.
/// @param file_path Destination VTU file.
/// @throws std::runtime_error If the file cannot be opened or written, or if
/// mesh data cannot be represented by the supported VTU format.
void write_vtu(const Mesh &mesh, const std::filesystem::path &file_path);

/// Writes a mesh and additional cell-centered arrays to an ASCII VTU file.
///
/// Existing intrinsic cell IDs, areas, and quality values are preserved.
/// Additional scalar arrays are written directly, while both supported 2D
/// vector representations are serialized as three VTK components `(x, y, 0)`.
/// This function does not take ownership of or copy descriptor values.
///
/// @param mesh Validated mesh to export.
/// @param file_path Destination VTU file.
/// @param cell_data Non-owning descriptors valid throughout this call.
/// @throws std::invalid_argument If a field has the wrong cardinality or an
///         empty name, or if its name conflicts with another user field or an
///         intrinsic CellData name (`cell_id`, `cell_area`, `cell_quality`).
/// @throws std::runtime_error If the file cannot be opened or written, or if
///         mesh data cannot be represented by the supported VTU format.
void write_vtu(const Mesh &mesh, const std::filesystem::path &file_path, const VtkCellData &cell_data);

} // namespace cfd
