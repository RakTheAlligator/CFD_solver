#pragma once

#include "cfd/Cell.hpp"
#include "cfd/Face.hpp"
#include "cfd/Node.hpp"
#include "cfd/Types.hpp"

#include <vector>

namespace cfd {

class Mesh {
public:
    Mesh() = default;

private:
    std::vector<Node> nodes_;

    // Cell connectivity will be implemented progressively.
    std::vector<CellType> cell_types_;
    std::vector<Index> cell_nodes_;
    std::vector<Index> cell_node_offsets_;

    // Face connectivity will be implemented after the first manual mesh.
    std::vector<Face> faces_;
};

} // namespace cfd
