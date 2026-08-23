#include "mesh_build/TopologyValidation.hpp"

#include "cfd/meshing/RawMeshData.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfd::detail {

namespace {

[[noreturn]]
void throw_topology_validation_error(
    const std::string& message)
{
    throw std::runtime_error(
        "Topology validation failed: " + message);
}

struct CellPair {
    Index cell_0{};
    Index cell_1{};

    bool operator==(const CellPair& other) const noexcept = default;

    bool operator<(const CellPair& other) const noexcept
    {
        if (cell_0 != other.cell_0)
        {
            return cell_0 < other.cell_0;
        }

        return cell_1 < other.cell_1;
    }
};

[[nodiscard]]
CellPair make_cell_pair(
    const Index cell_0,
    const Index cell_1) noexcept
{
    if (cell_0 < cell_1)
    {
        return { cell_0,
                 cell_1 };
    }

    return { cell_1,
             cell_0 };
}

void validate_topology_storage(
    const RawMeshData& raw_mesh,
    const TopologyBuildData& topology)
{
    if (topology.faces.empty())
    {
        throw_topology_validation_error(
            "no faces were constructed.");
    }
    if (topology.cell_faces.size()!= raw_mesh.cell_nodes.size())
    {
        throw_topology_validation_error(
            "cell_faces size must equal cell_nodes size.");
    }
    if (topology.face_adjacencies.size()!= topology.faces.size())
    {
        throw_topology_validation_error(
            "face_adjacencies size must equal faces size.");
    }
    if (topology.face_boundary_ids.size() != topology.faces.size())
    {
        throw_topology_validation_error(
            "face_boundary_ids size must equal faces size.");
    }
}

[[nodiscard]]
bool cell_references_face(
    const RawMeshData& raw_mesh,
    const TopologyBuildData& topology,
    const Index cell_id,
    const Index face_id) noexcept
{
    const Index begin { raw_mesh.cell_node_offsets[cell_id] };
    const Index end { raw_mesh.cell_node_offsets[cell_id + 1] };

    for (Index position = begin; position < end; ++position)
    {
        if (topology.cell_faces[position] == face_id)
        {
            return true;
        }
    }

    return false;
}

void validate_cell_face_connectivity(
    const RawMeshData& raw_mesh,
    const TopologyBuildData& topology)
{
    for (Index cell_id = 0; cell_id < raw_mesh.cell_types.size(); ++cell_id)
    {

        const Index begin { raw_mesh.cell_node_offsets[cell_id] };
        const Index end { raw_mesh.cell_node_offsets[cell_id + 1] };
        const Index local_face_count { end - begin };

        for (Index local_face = 0; local_face < local_face_count; ++local_face)
        {
            const Index position { begin + local_face };
            const Index face_id { topology.cell_faces[position] };

            if (face_id == invalid_index || face_id >= topology.faces.size())
            {
                throw_topology_validation_error(
                    "cell " + std::to_string(cell_id) + " references an invalid face.");
            }

            for (Index previous_position = begin; previous_position < position; ++previous_position)
            {
                if (topology.cell_faces[previous_position] == face_id)
                {
                    throw_topology_validation_error(
                        "cell " + std::to_string(cell_id) + " references the same face more than once.");
                }
            }

            const FaceAdjacency& adjacency { topology.face_adjacencies[face_id] };

            if (adjacency.owner != cell_id && adjacency.neighbor != cell_id)
            {
                throw_topology_validation_error(
                    "cell-to-face and face-to-cell connectivities are inconsistent.");
            }

            const Index next_local_node { (local_face + 1) % local_face_count };
            const Index expected_node_0 { raw_mesh.cell_nodes[begin + local_face] };
            const Index expected_node_1 { raw_mesh.cell_nodes[begin + next_local_node] };
            const Face& face { topology.faces[face_id] };

            const auto matches {
                (face.node_ids[0] == expected_node_0 && face.node_ids[1] == expected_node_1)
                ||
                (face.node_ids[0] == expected_node_1 && face.node_ids[1] == expected_node_0) };

            if (!matches)
            {
                throw_topology_validation_error(
                    "cell " + std::to_string(cell_id) + " references a face with incorrect nodes.");
            }
        }
    }
}

[[nodiscard]]
TopologyStats validate_faces(
    const RawMeshData& raw_mesh,
    const TopologyBuildData& topology)
{
    TopologyStats stats;

    for (Index face_id = 0; face_id < topology.faces.size(); ++face_id)
    {
        const Face& face { topology.faces[face_id] };
        const Index node_0 { face.node_ids[0] };
        const Index node_1 { face.node_ids[1] };

        if (node_0 >= raw_mesh.nodes.size() || node_1 >= raw_mesh.nodes.size())
        {
            throw_topology_validation_error(
                "face " + std::to_string(face_id) + " references a node outside the nodes array.");
        }
        if (node_0 == node_1)
        {
            throw_topology_validation_error(
                "face " + std::to_string(face_id) + " references the same node twice.");
        }

        const FaceAdjacency& adjacency {topology.face_adjacencies[face_id]};

        if (adjacency.owner >= raw_mesh.cell_types.size())
        {
            throw_topology_validation_error(
                "face " + std::to_string(face_id) + " has an invalid owner cell.");
        }
        if (!cell_references_face(raw_mesh, topology, adjacency.owner, face_id))
        {
            throw_topology_validation_error(
                "face " + std::to_string(face_id) + " is not referenced by its owner cell.");
        }

        const BoundaryId boundary_id { topology.face_boundary_ids[face_id] };

        if (adjacency.neighbor == invalid_index)
        {
            ++stats.boundary_face_count;

            if (boundary_id == invalid_boundary_id)
            {
                throw_topology_validation_error(
                    "external face " + std::to_string(face_id) + " has no boundary assignment.");
            }

            continue;
        }

        ++stats.internal_face_count;

        if (adjacency.neighbor >= raw_mesh.cell_types.size())
        {
            throw_topology_validation_error(
                "face " + std::to_string(face_id) + " has an invalid neighbor cell.");
        }
        if (adjacency.owner == adjacency.neighbor)
        {
            throw_topology_validation_error(
                "face " + std::to_string(face_id) + " has identical owner and neighbor cells.");
        }
        if (boundary_id != invalid_boundary_id)
        {
            throw_topology_validation_error(
                "internal face " + std::to_string(face_id) + " has a boundary assignment.");
        }
        if (!cell_references_face(raw_mesh, topology, adjacency.neighbor, face_id))
        {
            throw_topology_validation_error(
                "face "+ std::to_string(face_id) + " is not referenced by its neighbor cell.");
        }
    }

    if (stats.boundary_face_count != raw_mesh.boundary_edges.size())
    {
        throw_topology_validation_error(
            "number of external faces does not match the number of boundary edges.");
    }

    return stats;
}

void validate_unique_cell_neighbors(
    const TopologyBuildData& topology)
{
    std::vector<CellPair> neighboring_cell_pairs;

    neighboring_cell_pairs.reserve(topology.faces.size());

    for (const FaceAdjacency& adjacency : topology.face_adjacencies)
    {
        if (adjacency.neighbor == invalid_index)
        {
            continue;
        }

        neighboring_cell_pairs.push_back(make_cell_pair(adjacency.owner, adjacency.neighbor));
    }

    std::sort(neighboring_cell_pairs.begin(), neighboring_cell_pairs.end());

    const auto duplicate { std::adjacent_find(neighboring_cell_pairs.begin(), neighboring_cell_pairs.end()) };

    if (duplicate != neighboring_cell_pairs.end())
    {
        throw_topology_validation_error(
            "cells " + std::to_string(duplicate->cell_0) + " and " + std::to_string(duplicate->cell_1)+ " share more than one face.");
    }
}
} // namespace

TopologyStats validate_topology(
    const RawMeshData& raw_mesh,
    const TopologyBuildData& topology)
{
    validate_topology_storage(raw_mesh, topology);

    validate_cell_face_connectivity(raw_mesh, topology);

    const TopologyStats stats { validate_faces(raw_mesh, topology) };

    validate_unique_cell_neighbors(topology);

    if (raw_mesh.cell_nodes.size() != 2 * stats.internal_face_count + stats.boundary_face_count)
    {
        throw_topology_validation_error(
            "local face count is inconsistent with internal and boundary face counts.");
    }

    return stats;
}
} // namespace cfd::detail