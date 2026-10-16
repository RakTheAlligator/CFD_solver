#pragma once

#include "cfd/math/Point2.hpp"

#include <filesystem>
#include <functional>

namespace cfd
{
/// Spatial target length in metres, finite and strictly positive at every
/// sampled point. Used during generation only; not retained by Mesh.
using MeshSizeField = std::function<double(const Point2 &)>;

/// Loads a bilinear Cartesian size grid: nx nr, xmin xmax rmin rmax, then
/// nx*nr positive lengths (x varies fastest). Outside the grid, edge values
/// are extended constantly. No CFD variables or units conversion are involved.
[[nodiscard]] MeshSizeField read_mesh_size_field(const std::filesystem::path &file);
} // namespace cfd
