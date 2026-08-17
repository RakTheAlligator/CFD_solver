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

struct RawMeshData {                                  // Structure containing all raw mesh data
                                                      // imported from Gmsh and converted
                                                      // to our own CFD indexing.

    // ---------- NODES ----------

    std::vector<Node> nodes;                          // Array containing all mesh nodes.
                                                      //
                                                      // The index inside this vector is
                                                      // our internal CFD node identifier.
                                                      //
                                                      // Example:
                                                      // nodes[12] is CFD node number 12.
                                                      //
                                                      // Each Node currently stores
                                                      // its x and y coordinates.


    // ---------- CELLS ----------

    std::vector<CellType> cell_types;                 // Type of each cell.
                                                      //
                                                      // Example:
                                                      // cell_types[5] == CellType::Triangle
                                                      //
                                                      // There is exactly one entry per cell.
                                                      //
                                                      // This will later allow us to support
                                                      // both triangles and quadrilaterals.


    std::vector<Index> cell_nodes;                    // Flat array containing the node indices
                                                      // of all cells.
                                                      //
                                                      // Example:
                                                      //
                                                      // cell 0 -> nodes {0, 1, 4}
                                                      // cell 1 -> nodes {1, 5, 4}
                                                      //
                                                      // Then cell_nodes contains:
                                                      //
                                                      // {0, 1, 4, 1, 5, 4}
                                                      //
                                                      // We use one flat vector instead of
                                                      // one std::vector per cell to avoid
                                                      // many small dynamic allocations.


    std::vector<Index> cell_node_offsets;             // Array indicating where the node list
                                                      // of each cell starts inside cell_nodes.
                                                      //
                                                      // For the previous example:
                                                      //
                                                      // cell_node_offsets = {0, 3, 6}
                                                      //
                                                      // Cell 0 uses:
                                                      // cell_nodes[0 ... 3[
                                                      //
                                                      // Cell 1 uses:
                                                      // cell_nodes[3 ... 6[
                                                      //
                                                      // Therefore:
                                                      // cell_node_offsets.back()
                                                      // must always equal
                                                      // cell_nodes.size().


    // ---------- BOUNDARIES ----------

    std::vector<BoundaryGroup> boundary_groups;       // List of the physical boundary groups
                                                      // defined for the domain.
                                                      //
                                                      // Example:
                                                      //
                                                      // {id = 0, name = "inlet"}
                                                      // {id = 1, name = "wall"}
                                                      // {id = 2, name = "outlet"}
                                                      //
                                                      // This describes the categories
                                                      // of boundaries, not the individual
                                                      // boundary edges.


    std::vector<BoundaryEdge> boundary_edges;         // List of all mesh edges located
                                                      // on a physical boundary.
                                                      //
                                                      // Each BoundaryEdge stores:
                                                      //
                                                      // - the two node indices of the edge;
                                                      // - the boundary group identifier.
                                                      //
                                                      // Example:
                                                      //
                                                      // node_ids = {4, 17}
                                                      // boundary_id = 0
                                                      //
                                                      // means that the edge connecting
                                                      // CFD nodes 4 and 17 belongs
                                                      // to the "inlet" boundary group.
};

} // namespace cfd