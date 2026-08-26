#include "cfd/mesh/MeshStatistics.hpp"

#include "cfd/mesh/Mesh.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace cfd
{

namespace
{

struct ScalarAccumulator
{
    double minimum{std::numeric_limits<double>::infinity()};
    double maximum{-std::numeric_limits<double>::infinity()};
    double sum{};
    Index count{};
};

void add_value(ScalarAccumulator &accumulator, const double value) noexcept
{
    accumulator.minimum = std::min(accumulator.minimum, value);
    accumulator.maximum = std::max(accumulator.maximum, value);
    accumulator.sum += value;
    ++accumulator.count;
}

[[nodiscard]]
ScalarStats finalize(const ScalarAccumulator &accumulator) noexcept
{
    if (accumulator.count == 0)
    {
        return {};
    }

    return {
        .minimum = accumulator.minimum,
        .maximum = accumulator.maximum,
        .mean = accumulator.sum / static_cast<double>(accumulator.count),
    };
}

} // namespace

MeshStatistics compute_mesh_statistics(const Mesh &mesh)
{
    MeshStatistics statistics;

    ScalarAccumulator cell_area_statistics;
    ScalarAccumulator cell_size_statistics;
    ScalarAccumulator face_length_statistics;
    ScalarAccumulator cell_quality_statistics;

    for (const FaceAdjacency &adjacency : mesh.face_adjacencies())
    {
        if (adjacency.is_boundary())
        {
            ++statistics.boundary_face_count;
        }
        else
        {
            ++statistics.internal_face_count;
        }
    }

    for (const double area : mesh.cell_areas())
    {
        add_value(cell_area_statistics, area);
        add_value(cell_size_statistics, std::sqrt(area));

        statistics.total_cell_area += area;
    }

    for (const double length : mesh.face_lengths())
    {
        add_value(face_length_statistics, length);
    }

    double worst_quality{std::numeric_limits<double>::infinity()};

    for (Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        const double quality{mesh.cell_qualities()[cell_id]};

        add_value(cell_quality_statistics, quality);

        if (quality < worst_quality)
        {
            worst_quality = quality;
            statistics.worst_quality_cell = cell_id;
        }
    }

    statistics.cell_areas = finalize(cell_area_statistics);
    statistics.cell_sizes = finalize(cell_size_statistics);
    statistics.face_lengths = finalize(face_length_statistics);
    statistics.cell_quality = finalize(cell_quality_statistics);

    return statistics;
}

} // namespace cfd