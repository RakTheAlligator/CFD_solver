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
    Index value_count{};
};

void add_value(ScalarAccumulator &accumulator, const double value) noexcept
{
    accumulator.minimum = std::min(accumulator.minimum, value);
    accumulator.maximum = std::max(accumulator.maximum, value);
    accumulator.sum += value;
    ++accumulator.value_count;
}

[[nodiscard]]
ScalarStatistics finalize_statistics(const ScalarAccumulator &accumulator) noexcept
{
    if (accumulator.value_count == 0)
    {
        return {};
    }

    return {
        .minimum = accumulator.minimum,
        .maximum = accumulator.maximum,
        .mean = accumulator.sum / static_cast<double>(accumulator.value_count),
    };
}

} // namespace

MeshStatistics compute_mesh_statistics(const Mesh &mesh)
{
    MeshStatistics statistics;

    ScalarAccumulator cell_area_accumulator;
    ScalarAccumulator cell_size_accumulator;
    ScalarAccumulator face_length_accumulator;
    ScalarAccumulator cell_quality_accumulator;

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
        add_value(cell_area_accumulator, area);
        add_value(cell_size_accumulator, std::sqrt(area));

        statistics.total_cell_area += area;
    }

    for (const double length : mesh.face_lengths())
    {
        add_value(face_length_accumulator, length);
    }

    double worst_quality{std::numeric_limits<double>::infinity()};

    for (Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        const double quality{mesh.cell_qualities()[cell_id]};

        add_value(cell_quality_accumulator, quality);

        if (quality < worst_quality)
        {
            worst_quality = quality;
            statistics.worst_quality_cell_id = cell_id;
        }
    }

    statistics.cell_areas = finalize_statistics(cell_area_accumulator);
    statistics.cell_sizes = finalize_statistics(cell_size_accumulator);
    statistics.face_lengths = finalize_statistics(face_length_accumulator);
    statistics.cell_quality = finalize_statistics(cell_quality_accumulator);

    return statistics;
}

} // namespace cfd