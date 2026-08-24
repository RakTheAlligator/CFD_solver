#pragma once

#include <filesystem>

namespace cfd
{

class Mesh;

void write_vtu(const Mesh &mesh, const std::filesystem::path &file_path);

} // namespace cfd
