#include "mesh_build/MeshBuildData.hpp"

#include "cfd/mesh/Types.hpp"
#include "cfd/meshing/RawMeshData.hpp"

#include <functional>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace cfd::detail
{

namespace
{

[[noreturn]]
void throw_topology_build_error(const std::string &message)
{
    throw std::runtime_error("Topology construction failed: " + message);
}

struct FaceKey
{
    Index node_0_id{};
    Index node_1_id{};

    bool operator==(const FaceKey &other) const noexcept = default;
};

[[nodiscard]]
FaceKey make_face_key(const Index node_0_id, const Index node_1_id) noexcept
{
    if (node_0_id < node_1_id)
    {
        return {node_0_id, node_1_id};
    }

    return {node_1_id, node_0_id};
}

struct FaceKeyHash
{
    [[nodiscard]]
    std::size_t operator()(const FaceKey &face_key) const noexcept
    {
        const std::size_t hash_0{std::hash<Index>{}(face_key.node_0_id)};
        const std::size_t hash_1{std::hash<Index>{}(face_key.node_1_id)};

        return hash_0 ^ (hash_1 + 0x9e3779b9U + (hash_0 << 6U) + (hash_0 >> 2U));
    }
};

[[nodiscard]]
Index estimate_face_count(const RawMeshData &raw_mesh) noexcept
{
    return raw_mesh.cell_nodes.size() / 2 + raw_mesh.boundary_edges.size();
}

} // namespace

TopologyBuildData build_topology(const RawMeshData &raw_mesh)
{
    TopologyBuildData topology;

    const Index estimated_face_count{estimate_face_count(raw_mesh)};

    topology.faces.reserve(estimated_face_count);
    topology.face_adjacencies.reserve(estimated_face_count);
    topology.face_boundary_ids.reserve(estimated_face_count);

    topology.cell_faces.resize(raw_mesh.cell_nodes.size(), invalid_index);

    std::unordered_map<FaceKey, Index, FaceKeyHash> face_id_by_key;
    face_id_by_key.reserve(estimated_face_count);

    // -------------------------------------------------------------------------
    // Build unique faces and cell <-> face connectivity.
    // -------------------------------------------------------------------------

    for (Index cell_id = 0; cell_id < raw_mesh.cell_types.size(); ++cell_id)
    {
        const Index cell_node_begin_offset{raw_mesh.cell_node_offsets[cell_id]};
        const Index cell_node_end_offset{raw_mesh.cell_node_offsets[cell_id + 1]};
        const Index local_face_count{cell_node_end_offset - cell_node_begin_offset};

        for (Index local_face_index = 0; local_face_index < local_face_count; ++local_face_index)
        {
            const Index next_local_node_index{(local_face_index + 1) % local_face_count};
            const Index node_0_id{raw_mesh.cell_nodes[cell_node_begin_offset + local_face_index]};
            const Index node_1_id{raw_mesh.cell_nodes[cell_node_begin_offset + next_local_node_index]};

            const FaceKey face_key{make_face_key(node_0_id, node_1_id)};

            const Index candidate_face_id{topology.faces.size()};

            const auto [iterator, is_inserted]{face_id_by_key.try_emplace(face_key, candidate_face_id)};

            const Index face_id{iterator->second};

            if (is_inserted)
            {
                topology.faces.push_back(Face{{node_0_id, node_1_id}});

                topology.face_adjacencies.push_back(FaceAdjacency{.owner = cell_id, .neighbor = invalid_index});

                topology.face_boundary_ids.push_back(invalid_boundary_id);
            }
            else
            {
                FaceAdjacency &adjacency = topology.face_adjacencies[face_id];

                if (adjacency.owner == cell_id || adjacency.neighbor == cell_id)
                {
                    throw_topology_build_error("cell " + std::to_string(cell_id) +
                                               " references the same face more than once.");
                }

                if (adjacency.neighbor != invalid_index)
                {
                    throw_topology_build_error("face " + std::to_string(face_id) + " belongs to more than two cells.");
                }

                adjacency.neighbor = cell_id;
            }

            const Index cell_face_position{cell_node_begin_offset + local_face_index};

            if (topology.cell_faces[cell_face_position] != invalid_index)
            {
                throw_topology_build_error("cell-face connectivity position was assigned more than once.");
            }

            topology.cell_faces[cell_face_position] = face_id;
        }
    }

    if (face_id_by_key.size() != topology.faces.size())
    {
        throw_topology_build_error("internal face-indexing inconsistency.");
    }

    // -------------------------------------------------------------------------
    // Attach physical boundary information.
    // -------------------------------------------------------------------------

    for (const BoundaryEdge &boundary_edge : raw_mesh.boundary_edges)
    {
        const FaceKey face_key{make_face_key(boundary_edge.node_ids[0], boundary_edge.node_ids[1])};
        const auto iterator{face_id_by_key.find(face_key)};

        if (iterator == face_id_by_key.end())
        {
            throw_topology_build_error("a boundary edge does not correspond to any cell face.");
        }

        const Index face_id{iterator->second};
        const FaceAdjacency &adjacency{topology.face_adjacencies[face_id]};

        if (adjacency.neighbor != invalid_index)
        {
            throw_topology_build_error("face " + std::to_string(face_id) + " is internal but is marked as a boundary.");
        }

        if (topology.face_boundary_ids[face_id] != invalid_boundary_id)
        {
            throw_topology_build_error("face " + std::to_string(face_id) +
                                       " received more than one boundary assignment.");
        }

        topology.face_boundary_ids[face_id] = boundary_edge.boundary_id;
    }

    return topology;
}

} // namespace cfd::detail
