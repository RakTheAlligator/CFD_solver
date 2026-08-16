#pragma once

#include "cfd/Cell.hpp"
#include "cfd/Node.hpp"
#include "cfd/Types.hpp"

#include <array>
#include <string>
#include <vector>

namespace cfd {

using BoundaryId = Index;

struct BoundaryGroup {
    BoundaryId id{};
    std::string name;
};

struct BoundaryEdge {
    std::array<Index, 2> node_ids{};
    BoundaryId boundary_id{};
};

struct RawMeshData {
    // Coordinates
    std::vector<Node> nodes;

    // Cell -> node connectivity
    std::vector<CellType> cell_types;
    std::vector<Index> cell_nodes;
    std::vector<Index> cell_node_offsets;

    // Boundary information supplied by Gmsh
    std::vector<BoundaryGroup> boundary_groups;
    std::vector<BoundaryEdge> boundary_edges;
};

} // namespace cfd