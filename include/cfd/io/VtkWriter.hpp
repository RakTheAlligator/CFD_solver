#pragma once

#include <filesystem>

namespace cfd
{

class Mesh;

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

} // namespace cfd