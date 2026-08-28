#pragma once

#include "cfd/mesh/Boundary.hpp"
#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Node.hpp"
#include "cfd/mesh/Types.hpp"

#include <array>
#include <vector>

namespace cfd
{

/// Boundary edge extracted from the external mesh representation.
///
/// Node IDs use the internal zero-based numbering stored in `RawMeshData::nodes`.
/// `boundary_id` is a compact zero-based index into
/// `RawMeshData::boundary_groups`.
struct BoundaryEdge
{
    std::array<Index, 2> node_ids{};
    BoundaryId boundary_id{};
};

/// Mesh data extracted from the mesher before internal topology construction.
///
/// RawMeshData contains only the information imported from the external mesh:
/// nodes, cell-to-node connectivity, cell types, and physical boundary data.
/// Faces, cell-to-face connectivity, adjacency, and geometric quantities are
/// constructed later by the preprocessing pipeline.
///
/// Cell-to-node connectivity uses a flattened, CSR-like representation. Node
/// IDs for cell `c` occupy
/// `[cell_node_offsets[c], cell_node_offsets[c + 1])`.
///
/// All node and boundary IDs stored here use compact solver-internal
/// zero-based numbering. Node IDs index `nodes`, and boundary IDs index
/// `boundary_groups`; external Gmsh tags are converted before this
/// representation is returned.
///
/// @note RawMeshData is not assumed to be valid merely because it has been
///       constructed. `build_mesh()` validates its contents before building
///       the final Mesh representation.
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