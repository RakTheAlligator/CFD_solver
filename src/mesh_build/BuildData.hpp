#pragma once

#include "cfd/mesh/Boundary.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Types.hpp"
#include "cfd/mesh/Vector2.hpp"

#include <vector>

namespace cfd::detail
{
struct TopologyBuildData
{
    std::vector<Face> faces;

    // Same compressed layout as RawMeshData::cell_nodes.
    std::vector<Index> cell_faces;

    std::vector<FaceAdjacency> face_adjacencies;
    std::vector<BoundaryId> face_boundary_ids;
};

struct GeometryBuildData
{
    std::vector<double> cell_areas;
    std::vector<Vector2> cell_centers;
    std::vector<double> cell_qualities;

    std::vector<Vector2> face_centers;
    std::vector<double> face_lengths;
    std::vector<Vector2> face_area_vectors;
};
} // namespace cfd::detail