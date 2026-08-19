#pragma once

#include "cfd/mesh/Boundary.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Types.hpp"

#include <vector>

namespace cfd::detail {

struct TopologyBuildData {
    std::vector<Face> faces;

    // Same compressed layout as RawMeshData::cell_nodes.
    std::vector<Index> cell_faces;

    std::vector<FaceAdjacency> face_adjacencies;
    std::vector<BoundaryId> face_boundary_ids;
};

struct TopologyStats {
    Index internal_face_count{};
    Index boundary_face_count{};
};

} // namespace cfd::detail