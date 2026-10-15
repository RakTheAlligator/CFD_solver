#include "cfd/meshing/MeshSizeField.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfd
{
MeshSizeField read_mesh_size_field(const std::filesystem::path &file)
{
    std::ifstream input{file};
    std::size_t nx{}, nr{};
    double xmin{}, xmax{}, rmin{}, rmax{};
    if (!(input >> nx >> nr >> xmin >> xmax >> rmin >> rmax) || nx < 2 || nr < 2 || nx > 10000000 / nr ||
        !std::isfinite(xmin) || !std::isfinite(xmax) || !std::isfinite(rmin) || !std::isfinite(rmax) ||
        !(xmax > xmin) || !(rmax > rmin) || !std::isfinite(xmax - xmin) || !std::isfinite(rmax - rmin))
    {
        throw std::invalid_argument("Invalid external mesh-size grid header: " + file.string());
    }
    std::vector<double> sizes(nx * nr);
    for (double &value : sizes)
    {
        if (!(input >> value) || !std::isfinite(value) || !(value > 0.0))
        {
            throw std::invalid_argument("External mesh-size grid needs nx*nr finite positive lengths.");
        }
    }
    std::string trailing;
    if (input >> trailing)
    {
        throw std::invalid_argument("Unexpected data after external mesh-size grid.");
    }
    return [nx, nr, xmin, xmax, rmin, rmax, sizes = std::move(sizes)](const Point2 &point) {
        if (!std::isfinite(point.x) || !std::isfinite(point.y))
            throw std::invalid_argument("Mesh-size grid requires finite sampling coordinates.");
        const double x{std::clamp((point.x - xmin) / (xmax - xmin), 0.0, 1.0) * static_cast<double>(nx - 1)};
        const double r{std::clamp((point.y - rmin) / (rmax - rmin), 0.0, 1.0) * static_cast<double>(nr - 1)};
        const auto i{std::min(static_cast<std::size_t>(x), nx - 2)};
        const auto j{std::min(static_cast<std::size_t>(r), nr - 2)};
        const double a{std::lerp(sizes[j * nx + i], sizes[j * nx + i + 1], x - static_cast<double>(i))};
        const double b{std::lerp(sizes[(j + 1) * nx + i], sizes[(j + 1) * nx + i + 1], x - static_cast<double>(i))};
        return std::lerp(a, b, r - static_cast<double>(j));
    };
}
} // namespace cfd
