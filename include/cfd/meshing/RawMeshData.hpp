#pragma once

#include "cfd/mesh/Boundary.hpp"
#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Node.hpp"
#include "cfd/mesh/Types.hpp"

#include <array>
#include <vector>

namespace cfd
{

struct BoundaryEdge
{
    std::array<Index, 2> node_ids{};
    BoundaryId boundary_id{};
};

struct RawMeshData
{
    std::vector<Node> nodes;

    std::vector<CellType> cell_types;
    std::vector<Index> cell_nodes;
    std::vector<Index> cell_node_offsets;

    std::vector<BoundaryGroup> boundary_groups;
    std::vector<BoundaryEdge> boundary_edges;
};

} // namespace cfd