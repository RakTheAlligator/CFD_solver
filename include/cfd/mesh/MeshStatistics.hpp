#pragma once

#include "cfd/mesh/Types.hpp"

namespace cfd
{

class Mesh;

struct ScalarStats
{
    double minimum{};
    double maximum{};
    double mean{};
};

struct MeshStatistics
{
    Index internal_face_count{};
    Index boundary_face_count{};

    double total_cell_area{};

    ScalarStats cell_areas;
    ScalarStats cell_sizes;
    ScalarStats face_lengths;
    ScalarStats triangle_quality;

    Index worst_quality_cell{invalid_index};
};

[[nodiscard]]
MeshStatistics compute_mesh_statistics(const Mesh &mesh);

} // namespace cfd