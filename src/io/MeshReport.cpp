#include "cfd/io/MeshReport.hpp"

#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/mesh/MeshStatistics.hpp"
#include "cfd/mesh/Types.hpp"

#include <iomanip>
#include <ostream>
#include <sstream>
#include <string_view>

namespace cfd
{

namespace
{

void write_statistics_row(std::ostream &output, const std::string_view quantity, const std::string_view unit,
                          const ScalarStatistics &statistics)
{
    output << "  " << std::left << std::setw(22) << quantity << std::setw(8) << unit << std::right << std::defaultfloat
           << std::setprecision(5) << std::setw(12) << statistics.minimum << std::setw(12) << statistics.mean
           << std::setw(12) << statistics.maximum << '\n';
}

} // namespace

void write_mesh_report(std::ostream &output, const Mesh &mesh, const MeshStatistics &statistics,
                       const MeshBuildTimings &timings)
{
    // Format into a local stream so manipulators such as fixed, left, and
    // setprecision do not alter the formatting state of the caller's stream.
    std::ostringstream report;

    report << "\n[Mesh topology]\n"
           << "  Faces             : " << mesh.face_count() << '\n'
           << "    Internal         : " << statistics.internal_face_count << '\n'
           << "    Boundary         : " << statistics.boundary_face_count << '\n'
           << std::fixed << std::setprecision(2) << "  Time              : " << timings.topology.count() << " ms\n";

    report << "\n[Mesh geometry]\n"
           << std::fixed << std::setprecision(4) << "  Total area        : " << statistics.total_cell_area << " m^2\n"
           << std::setprecision(2) << "  Time              : " << timings.geometry.count() << " ms\n";

    report << "\n  " << std::left << std::setw(22) << "Quantity" << std::setw(8) << "Unit" << std::right
           << std::setw(12) << "min" << std::setw(12) << "mean" << std::setw(12) << "max" << '\n';

    report << "  ------------------------------------------------------------------\n";

    write_statistics_row(report, "Cell area", "m^2", statistics.cell_areas);
    write_statistics_row(report, "Cell size", "m", statistics.cell_sizes);
    write_statistics_row(report, "Face length", "m", statistics.face_lengths);
    write_statistics_row(report, "Internal non-orth.", "deg", statistics.internal_face_non_orthogonality_degrees);
    write_statistics_row(report, "Neighbor size ratio", "-", statistics.internal_face_neighbor_cell_size_ratios);

    if (statistics.worst_quality_cell_id != invalid_index)
    {
        write_statistics_row(report, "Cell quality", "-", statistics.cell_quality);
        report << "\n  Worst cell        : cell " << statistics.worst_quality_cell_id << '\n';
    }
    if (statistics.maximum_non_orthogonality_face_id != invalid_index)
    {
        report << "  Max non-orth. face : face " << statistics.maximum_non_orthogonality_face_id << '\n';
    }
    if (statistics.maximum_neighbor_cell_size_ratio_face_id != invalid_index)
    {
        const Index face_id{statistics.maximum_neighbor_cell_size_ratio_face_id};
        const FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};
        const Point2 &center{mesh.face_centers()[face_id]};
        report << "  Max size-ratio face: face " << face_id << '\n'
               << "    Owner / neighbor : cell " << adjacency.owner << " / cell " << adjacency.neighbor << '\n'
               << "    Ratio            : " << std::setprecision(5)
               << statistics.internal_face_neighbor_cell_size_ratios.maximum << '\n'
               << "    Face center      : (" << center.x << ", " << center.y << ") m\n";
    }

    output << report.str();
}

} // namespace cfd
