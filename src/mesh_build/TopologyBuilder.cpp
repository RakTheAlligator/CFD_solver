#include "mesh_build/TopologyBuilder.hpp"

#include "cfd/mesh/Types.hpp"
#include "cfd/meshing/RawMeshData.hpp"

#include <functional>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace cfd::detail {

namespace {

[[noreturn]]
void throw_topology_build_error(
    const std::string& message)
{
    throw std::runtime_error(
        "Topology construction failed: " + message);
}
struct EdgeKey {
    Index node_0{};
    Index node_1{};

    bool operator==(const EdgeKey& other) const noexcept = default;
};

[[nodiscard]]
EdgeKey make_edge_key(
    const Index node_0,
    const Index node_1) noexcept
{
    if (node_0 < node_1)
    {
        return {node_0, node_1};
    }

    return {node_1, node_0};
}

struct EdgeKeyHash {
    [[nodiscard]]
    std::size_t operator()(const EdgeKey& edge) const noexcept
    {
        const std::size_t hash_0 { std::hash<Index>{}(edge.node_0) };
        const std::size_t hash_1 { std::hash<Index>{}(edge.node_1) };

        return hash_0
            ^ (hash_1
               + 0x9e3779b9U
               + (hash_0 << 6U)
               + (hash_0 >> 2U));
    }
};

[[nodiscard]]
Index estimate_face_count(
    const RawMeshData& raw_mesh) noexcept
{
    return raw_mesh.cell_nodes.size() / 2
        + raw_mesh.boundary_edges.size();
}

} // namespace

TopologyBuildData build_topology(
    const RawMeshData& raw_mesh)
{
    TopologyBuildData topology;

    const Index estimated_face_count { estimate_face_count(raw_mesh) };

    topology.faces.reserve(estimated_face_count);
    topology.face_adjacencies.reserve(estimated_face_count);
    topology.face_boundary_ids.reserve(estimated_face_count);

    topology.cell_faces.resize(raw_mesh.cell_nodes.size(), invalid_index);

    std::unordered_map<EdgeKey, Index, EdgeKeyHash> face_ids;
    face_ids.reserve(estimated_face_count);

    // -------------------------------------------------------------------------
    // Build unique faces and cell <-> face connectivity.
    // -------------------------------------------------------------------------

    for (Index cell_id = 0; cell_id < raw_mesh.cell_types.size(); ++cell_id)
    {
        const Index begin { raw_mesh.cell_node_offsets[cell_id] };
        const Index end { raw_mesh.cell_node_offsets[cell_id + 1] };
        const Index node_count { end - begin };

        for (Index local_face = 0; local_face < node_count; ++local_face)
        {
            const Index next_local_node { (local_face + 1) % node_count };
            const Index node_0 { raw_mesh.cell_nodes[begin + local_face] };
            const Index node_1 { raw_mesh.cell_nodes[begin + next_local_node] };

            const EdgeKey edge_key { make_edge_key(node_0, node_1) };

            const Index candidate_face_id { topology.faces.size() };

            const auto [iterator, inserted] { face_ids.try_emplace(edge_key, candidate_face_id) };

            const Index face_id { iterator->second };

            if (inserted)
            {
                topology.faces.push_back(Face{{node_0, node_1}});

                topology.face_adjacencies.push_back(
                    FaceAdjacency{
                        .owner = cell_id,
                        .neighbor = invalid_index
                    });

                topology.face_boundary_ids.push_back(
                    invalid_boundary_id);
            }
            else
            {
                FaceAdjacency& adjacency = topology.face_adjacencies[face_id];

                if (adjacency.owner == cell_id || adjacency.neighbor == cell_id)
                {
                    throw_topology_build_error(
                        "cell " + std::to_string(cell_id) + " references the same face more than once.");
                }
                if (adjacency.neighbor != invalid_index)
                {
                    throw_topology_build_error(
                        "face "+ std::to_string(face_id) + " belongs to more than two cells.");
                }

                adjacency.neighbor = cell_id;
            }

            const Index cell_face_position { begin + local_face };

            if (topology.cell_faces[cell_face_position] != invalid_index)
            {
                throw_topology_build_error(
                    "cell-face connectivity position was assigned more than once.");
            }

            topology.cell_faces[cell_face_position] = face_id;
        }
    }

    if (face_ids.size() != topology.faces.size())
    {
        throw_topology_build_error(
            "internal face-indexing inconsistency.");
    }

    // -------------------------------------------------------------------------
    // Attach physical boundary information.
    // -------------------------------------------------------------------------

    for (const BoundaryEdge& boundary_edge : raw_mesh.boundary_edges)
    {
        const EdgeKey edge_key { make_edge_key(boundary_edge.node_ids[0], boundary_edge.node_ids[1]) };
        const auto iterator { face_ids.find(edge_key) };

        if (iterator == face_ids.end())
        {
            throw_topology_build_error(
                "a boundary edge does not correspond to any cell face.");
        }

        const Index face_id { iterator->second };
        const FaceAdjacency& adjacency { topology.face_adjacencies[face_id] };

        if (adjacency.neighbor != invalid_index)
        {
            throw_topology_build_error(
                "face " + std::to_string(face_id) + " is internal but is marked as a boundary.");
        }
        if (topology.face_boundary_ids[face_id] != invalid_boundary_id)
        {
            throw_topology_build_error(
                "face " + std::to_string(face_id) + " received more than one boundary assignment.");
        }

        topology.face_boundary_ids[face_id] = boundary_edge.boundary_id;
    }

    return topology;
}
} // namespace cfd::detail