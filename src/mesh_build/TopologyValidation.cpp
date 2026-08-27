#include "mesh_build/TopologyValidation.hpp"

#include "cfd/meshing/RawMeshData.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfd::detail
{

namespace
{

struct CellPair
{
    Index cell_0_id{};
    Index cell_1_id{};

    bool operator==(const CellPair &other) const noexcept = default;

    bool operator<(const CellPair &other) const noexcept
    {
        if (cell_0_id != other.cell_0_id)
        {
            return cell_0_id < other.cell_0_id;
        }

        return cell_1_id < other.cell_1_id;
    }
};

struct FaceCounts
{
    Index internal_face_count{};
    Index boundary_face_count{};
};

[[noreturn]]
void throw_topology_validation_error(const std::string &message)
{
    throw std::runtime_error("Topology validation failed: " + message);
}

// Canonicalize the cell pair so that adjacency identity is independent of
// which cell was assigned as owner during topology construction.
[[nodiscard]]
CellPair make_cell_pair(const Index cell_0_id, const Index cell_1_id) noexcept
{
    if (cell_0_id < cell_1_id)
    {
        return {cell_0_id, cell_1_id};
    }

    return {cell_1_id, cell_0_id};
}

void validate_topology_storage(const RawMeshData &raw_mesh, const TopologyBuildData &topology)
{
    if (topology.faces.empty())
    {
        throw_topology_validation_error("no faces were constructed.");
    }

    if (topology.cell_faces.size() != raw_mesh.cell_nodes.size())
    {
        throw_topology_validation_error("cell_faces size must equal cell_nodes size.");
    }

    if (topology.face_adjacencies.size() != topology.faces.size())
    {
        throw_topology_validation_error("face_adjacencies size must equal faces size.");
    }

    if (topology.face_boundary_ids.size() != topology.faces.size())
    {
        throw_topology_validation_error("face_boundary_ids size must equal faces size.");
    }
}

[[nodiscard]]
bool cell_references_face(const RawMeshData &raw_mesh, const TopologyBuildData &topology, const Index cell_id,
                          const Index face_id) noexcept
{
    const Index cell_face_begin_offset{raw_mesh.cell_node_offsets[cell_id]};
    const Index cell_face_end_offset{raw_mesh.cell_node_offsets[cell_id + 1]};

    for (Index cell_face_position = cell_face_begin_offset; cell_face_position < cell_face_end_offset;
         ++cell_face_position)
    {
        if (topology.cell_faces[cell_face_position] == face_id)
        {
            return true;
        }
    }

    return false;
}

void validate_cell_face_connectivity(const RawMeshData &raw_mesh, const TopologyBuildData &topology)
{
    for (Index cell_id = 0; cell_id < raw_mesh.cell_types.size(); ++cell_id)
    {
        const Index cell_face_begin_offset{raw_mesh.cell_node_offsets[cell_id]};
        const Index cell_face_end_offset{raw_mesh.cell_node_offsets[cell_id + 1]};
        const Index local_face_count{cell_face_end_offset - cell_face_begin_offset};

        for (Index local_face_index = 0; local_face_index < local_face_count; ++local_face_index)
        {
            const Index cell_face_position{cell_face_begin_offset + local_face_index};
            const Index face_id{topology.cell_faces[cell_face_position]};

            if (face_id == invalid_index || face_id >= topology.faces.size())
            {
                throw_topology_validation_error("cell " + std::to_string(cell_id) + " references an invalid face.");
            }

            // Cells currently have only three or four local faces. A bounded
            // O(n^2) scan detects duplicates without allocating temporary
            // storage for every cell.
            for (Index previous_cell_face_position = cell_face_begin_offset;
                 previous_cell_face_position < cell_face_position; ++previous_cell_face_position)
            {
                if (topology.cell_faces[previous_cell_face_position] == face_id)
                {
                    throw_topology_validation_error("cell " + std::to_string(cell_id) +
                                                    " references the same face more than once.");
                }
            }

            const FaceAdjacency &adjacency{topology.face_adjacencies[face_id]};

            if (adjacency.owner != cell_id && adjacency.neighbor != cell_id)
            {
                throw_topology_validation_error("cell-to-face and face-to-cell connectivities are inconsistent.");
            }

            const Index next_local_node_index{(local_face_index + 1) % local_face_count};
            const Index expected_node_0_id{raw_mesh.cell_nodes[cell_face_begin_offset + local_face_index]};
            const Index expected_node_1_id{raw_mesh.cell_nodes[cell_face_begin_offset + next_local_node_index]};

            const Face &face{topology.faces[face_id]};

            // Face identity is orientation-independent: either stored node
            // ordering is valid as long as it matches the local cell edge.
            const bool has_expected_nodes{
                (face.node_ids[0] == expected_node_0_id && face.node_ids[1] == expected_node_1_id) ||
                (face.node_ids[0] == expected_node_1_id && face.node_ids[1] == expected_node_0_id)};

            if (!has_expected_nodes)
            {
                throw_topology_validation_error("cell " + std::to_string(cell_id) +
                                                " references a face with incorrect nodes.");
            }
        }
    }
}

[[nodiscard]]
FaceCounts validate_faces(const RawMeshData &raw_mesh, const TopologyBuildData &topology)
{
    FaceCounts face_counts;

    for (Index face_id = 0; face_id < topology.faces.size(); ++face_id)
    {
        const Face &face{topology.faces[face_id]};

        const Index node_0_id{face.node_ids[0]};
        const Index node_1_id{face.node_ids[1]};

        if (node_0_id >= raw_mesh.nodes.size() || node_1_id >= raw_mesh.nodes.size())
        {
            throw_topology_validation_error("face " + std::to_string(face_id) +
                                            " references a node outside the nodes array.");
        }

        if (node_0_id == node_1_id)
        {
            throw_topology_validation_error("face " + std::to_string(face_id) + " references the same node twice.");
        }

        const FaceAdjacency &adjacency{topology.face_adjacencies[face_id]};

        if (adjacency.owner >= raw_mesh.cell_types.size())
        {
            throw_topology_validation_error("face " + std::to_string(face_id) + " has an invalid owner cell.");
        }

        if (!cell_references_face(raw_mesh, topology, adjacency.owner, face_id))
        {
            throw_topology_validation_error("face " + std::to_string(face_id) +
                                            " is not referenced by its owner cell.");
        }

        const BoundaryId boundary_id{topology.face_boundary_ids[face_id]};

        // Boundary status is encoded topologically by the absence of a
        // neighbor. Boundary assignment must agree with that topology.
        if (adjacency.neighbor == invalid_index)
        {
            ++face_counts.boundary_face_count;

            if (boundary_id == invalid_boundary_id)
            {
                throw_topology_validation_error("external face " + std::to_string(face_id) +
                                                " has no boundary assignment.");
            }

            continue;
        }

        ++face_counts.internal_face_count;

        if (adjacency.neighbor >= raw_mesh.cell_types.size())
        {
            throw_topology_validation_error("face " + std::to_string(face_id) + " has an invalid neighbor cell.");
        }

        if (adjacency.owner == adjacency.neighbor)
        {
            throw_topology_validation_error("face " + std::to_string(face_id) +
                                            " has identical owner and neighbor cells.");
        }

        if (boundary_id != invalid_boundary_id)
        {
            throw_topology_validation_error("internal face " + std::to_string(face_id) + " has a boundary assignment.");
        }

        if (!cell_references_face(raw_mesh, topology, adjacency.neighbor, face_id))
        {
            throw_topology_validation_error("face " + std::to_string(face_id) +
                                            " is not referenced by its neighbor cell.");
        }
    }

    if (face_counts.boundary_face_count != raw_mesh.boundary_edges.size())
    {
        throw_topology_validation_error("number of external faces does not match the number of boundary edges.");
    }

    return face_counts;
}

void validate_unique_cell_neighbors(const TopologyBuildData &topology)
{
    std::vector<CellPair> neighboring_cell_pairs;
    neighboring_cell_pairs.reserve(topology.faces.size());

    for (const FaceAdjacency &adjacency : topology.face_adjacencies)
    {
        if (adjacency.neighbor == invalid_index)
        {
            continue;
        }

        neighboring_cell_pairs.push_back(make_cell_pair(adjacency.owner, adjacency.neighbor));
    }

    // This O(F log F) validation runs only during preprocessing. Sorting
    // canonical cell pairs keeps the final Mesh free of additional lookup
    // structures while detecting cells that incorrectly share multiple faces.
    std::sort(neighboring_cell_pairs.begin(), neighboring_cell_pairs.end());

    const auto duplicate_pair_iterator{
        std::adjacent_find(neighboring_cell_pairs.begin(), neighboring_cell_pairs.end())};

    if (duplicate_pair_iterator != neighboring_cell_pairs.end())
    {
        throw_topology_validation_error("cells " + std::to_string(duplicate_pair_iterator->cell_0_id) + " and " +
                                        std::to_string(duplicate_pair_iterator->cell_1_id) +
                                        " share more than one face.");
    }
}

} // namespace

void validate_topology(const RawMeshData &raw_mesh, const TopologyBuildData &topology)
{
    validate_topology_storage(raw_mesh, topology);
    validate_cell_face_connectivity(raw_mesh, topology);

    const FaceCounts face_counts{validate_faces(raw_mesh, topology)};

    validate_unique_cell_neighbors(topology);

    // Count cell-face incidences independently of the constructed face count.
    // Every internal face contributes twice and every boundary face once.
    if (raw_mesh.cell_nodes.size() != 2 * face_counts.internal_face_count + face_counts.boundary_face_count)
    {
        throw_topology_validation_error("local face count is inconsistent with internal and boundary face counts.");
    }
}

} // namespace cfd::detail