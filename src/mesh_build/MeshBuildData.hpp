#pragma once

#include "cfd/math/Point2.hpp"
#include "cfd/math/Vector2.hpp"
#include "cfd/mesh/Boundary.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Types.hpp"

#include <vector>

namespace cfd::detail
{

// Temporary topology storage produced before validation and transfer into Mesh.
struct TopologyBuildData
{
    std::vector<Face> faces;

    // Reuse RawMeshData::cell_node_offsets: in a 2D polygon, each cell has one
    // face per node, so cell_faces has the same flattened layout as cell_nodes.
    std::vector<Index> cell_faces;

    std::vector<FaceAdjacency> face_adjacencies;
    std::vector<BoundaryId> face_boundary_ids;
};

// Temporary geometry storage indexed directly by internal cell and face IDs.
struct GeometryBuildData
{
    std::vector<double> cell_areas;
    std::vector<Point2> cell_centers;
    std::vector<double> cell_qualities;

    std::vector<Point2> face_centers;
    std::vector<double> face_lengths;
    std::vector<Vector2> face_area_vectors;
};

} // namespace cfd::detail
