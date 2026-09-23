#include "cfd/mesh/MeshStatistics.hpp"

#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/Types.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

namespace cfd
{

namespace
{

// Running accumulator for one-pass descriptive statistics. Derived quantities
// do not need to be stored solely for reporting their minimum, maximum, and mean.
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
    ScalarAccumulator non_orthogonality_accumulator;
    ScalarAccumulator neighbor_cell_size_ratio_accumulator;

    double maximum_non_orthogonality{-std::numeric_limits<double>::infinity()};
    double maximum_neighbor_cell_size_ratio{-std::numeric_limits<double>::infinity()};

    for (Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        const FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};
        if (adjacency.is_boundary())
        {
            ++statistics.boundary_face_count;
            continue;
        }

        ++statistics.internal_face_count;

        const Point2 &owner_center{mesh.cell_centers()[adjacency.owner]};
        const Point2 &neighbor_center{mesh.cell_centers()[adjacency.neighbor]};
        const double displacement_x{neighbor_center.x - owner_center.x};
        const double displacement_y{neighbor_center.y - owner_center.y};
        const Vector2 &area_vector{mesh.face_area_vectors()[face_id]};
        const double cosine{std::clamp((displacement_x * area_vector.x + displacement_y * area_vector.y) /
                                           (std::hypot(displacement_x, displacement_y) * mesh.face_lengths()[face_id]),
                                       -1.0, 1.0)};
        const double non_orthogonality{std::acos(cosine) * 180.0 / std::numbers::pi_v<double>};
        add_value(non_orthogonality_accumulator, non_orthogonality);
        if (non_orthogonality > maximum_non_orthogonality)
        {
            maximum_non_orthogonality = non_orthogonality;
            statistics.maximum_non_orthogonality_face_id = face_id;
        }

        const double owner_size{std::sqrt(mesh.cell_areas()[adjacency.owner])};
        const double neighbor_size{std::sqrt(mesh.cell_areas()[adjacency.neighbor])};
        const double neighbor_cell_size_ratio{std::max(owner_size / neighbor_size, neighbor_size / owner_size)};
        add_value(neighbor_cell_size_ratio_accumulator, neighbor_cell_size_ratio);
        if (neighbor_cell_size_ratio > maximum_neighbor_cell_size_ratio)
        {
            maximum_neighbor_cell_size_ratio = neighbor_cell_size_ratio;
            statistics.maximum_neighbor_cell_size_ratio_face_id = face_id;
        }
    }

    for (const double area : mesh.cell_areas())
    {
        add_value(cell_area_accumulator, area);

        // sqrt(area) provides a topology-independent characteristic linear
        // cell size; it is not intended to represent an actual edge length.
        add_value(cell_size_accumulator, std::sqrt(area));

        statistics.total_cell_area += area;
    }

    for (const double length : mesh.face_lengths())
    {
        add_value(face_length_accumulator, length);
    }

    double worst_quality{std::numeric_limits<double>::infinity()};

    // A strict comparison keeps the first cell ID when several cells share the
    // same minimum quality, making the reported result deterministic.
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
    statistics.internal_face_non_orthogonality_degrees = finalize_statistics(non_orthogonality_accumulator);
    statistics.internal_face_neighbor_cell_size_ratios = finalize_statistics(neighbor_cell_size_ratio_accumulator);

    return statistics;
}

} // namespace cfd
