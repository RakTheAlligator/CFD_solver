#include "cfd/MeshBuilder.hpp"

#include "cfd/RawMeshData.hpp"
#include "cfd/RawMeshValidation.hpp"

#include <algorithm>
#include <chrono>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace cfd {

namespace {

// -----------------------------------------------------------------------------
// Error reporting
// -----------------------------------------------------------------------------

[[noreturn]]
void throw_topology_build_error(const std::string& message)
{
    throw std::runtime_error(
        "Topology construction failed: " + message);
}

[[noreturn]]
void throw_topology_validation_error(const std::string& message)
{
    throw std::runtime_error(
        "Topology validation failed: " + message);
}

// -----------------------------------------------------------------------------
// Canonical edge representation
// -----------------------------------------------------------------------------

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
    if (node_0 < node_1) {
        return {node_0, node_1};
    }

    return {node_1, node_0};
}

struct EdgeKeyHash {
    [[nodiscard]]
    std::size_t operator()(const EdgeKey& edge) const noexcept
    {
        const std::size_t hash_0 =
            std::hash<Index>{}(edge.node_0);

        const std::size_t hash_1 =
            std::hash<Index>{}(edge.node_1);

        // Combine the two node hashes.
        return hash_0
            ^ (hash_1
               + 0x9e3779b9U
               + (hash_0 << 6U)
               + (hash_0 >> 2U));
    }
};

// -----------------------------------------------------------------------------
// Cell-pair representation
//
// Used during validation to verify that two cells never share more than
// one face.
// -----------------------------------------------------------------------------

struct CellPair {
    Index cell_0{};
    Index cell_1{};

    bool operator==(const CellPair& other) const noexcept = default;

    bool operator<(const CellPair& other) const noexcept
    {
        if (cell_0 != other.cell_0) {
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
    if (cell_0 < cell_1) {
        return {cell_0, cell_1};
    }

    return {cell_1, cell_0};
}

// -----------------------------------------------------------------------------
// Temporary topology
//
// This object exists only while Mesh is being built.
// If construction or validation throws, all its vectors are automatically
// destroyed.
// -----------------------------------------------------------------------------

struct TopologyBuildData {
    std::vector<Face> faces;

    // Same compressed layout as RawMeshData::cell_nodes.
    // cell_node_offsets therefore also delimit cell_faces.
    std::vector<Index> cell_faces;

    std::vector<FaceAdjacency> face_adjacencies;
    std::vector<BoundaryId> face_boundary_ids;
};

// -----------------------------------------------------------------------------
// Topology statistics
//
// Produced as a by-product of validation.
// They are useful both for consistency checks and concise diagnostics.
// -----------------------------------------------------------------------------

struct TopologyStats {
    Index internal_face_count{};
    Index boundary_face_count{};
};

// -----------------------------------------------------------------------------
// Reservation estimate
//
// For a valid 2D manifold mesh, most local cell edges are shared by two cells.
// This is only a capacity estimate. Correctness never depends on it.
// -----------------------------------------------------------------------------

[[nodiscard]]
Index estimate_face_count(const RawMeshData& raw_mesh) noexcept
{
    return raw_mesh.cell_nodes.size() / 2
        + raw_mesh.boundary_edges.size();
}

// -----------------------------------------------------------------------------
// Topology construction
// -----------------------------------------------------------------------------

[[nodiscard]]
TopologyBuildData build_topology(const RawMeshData& raw_mesh)
{
    TopologyBuildData topology;

    const Index estimated_face_count =
        estimate_face_count(raw_mesh);

    topology.faces.reserve(estimated_face_count);
    topology.face_adjacencies.reserve(estimated_face_count);
    topology.face_boundary_ids.reserve(estimated_face_count);

    // In our 2D polygonal mesh:
    //
    // triangle -> 3 nodes -> 3 local faces
    // quad     -> 4 nodes -> 4 local faces
    //
    // Therefore cell_faces uses exactly the same compressed layout
    // as cell_nodes.
    topology.cell_faces.resize(
        raw_mesh.cell_nodes.size(),
        invalid_index);

    std::unordered_map<EdgeKey, Index, EdgeKeyHash> face_ids;
    face_ids.reserve(estimated_face_count);

    // -------------------------------------------------------------------------
    // Build unique faces and cell <-> face connectivity.
    // -------------------------------------------------------------------------

    for (Index cell_id = 0;
         cell_id < raw_mesh.cell_types.size();
         ++cell_id) {

        const Index begin =
            raw_mesh.cell_node_offsets[cell_id];

        const Index end =
            raw_mesh.cell_node_offsets[cell_id + 1];

        const Index node_count =
            end - begin;

        for (Index local_face = 0;
             local_face < node_count;
             ++local_face) {

            const Index next_local_node =
                (local_face + 1) % node_count;

            const Index node_0 =
                raw_mesh.cell_nodes[
                    begin + local_face];

            const Index node_1 =
                raw_mesh.cell_nodes[
                    begin + next_local_node];

            const EdgeKey edge_key =
                make_edge_key(node_0, node_1);

            // This value will only be used if the EdgeKey is new.
            const Index candidate_face_id =
                topology.faces.size();

            const auto [iterator, inserted] =
                face_ids.try_emplace(
                    edge_key,
                    candidate_face_id);

            const Index face_id =
                iterator->second;

            if (inserted) {
                // First occurrence of this edge.
                //
                // EdgeKey is canonical only for identifying an
                // undirected edge.
                //
                // Face itself keeps the node order encountered from
                // its owner cell, which preserves orientation information.
                topology.faces.push_back(
                    Face{{node_0, node_1}});

                topology.face_adjacencies.push_back(
                    FaceAdjacency{
                        .owner = cell_id,
                        .neighbor = invalid_index
                    });

                topology.face_boundary_ids.push_back(
                    invalid_boundary_id);
            } else {
                // The face already exists.
                FaceAdjacency& adjacency =
                    topology.face_adjacencies[face_id];

                if (adjacency.owner == cell_id
                    || adjacency.neighbor == cell_id) {

                    throw_topology_build_error(
                        "cell "
                        + std::to_string(cell_id)
                        + " references the same face more than once.");
                }

                if (adjacency.neighbor != invalid_index) {
                    throw_topology_build_error(
                        "face "
                        + std::to_string(face_id)
                        + " belongs to more than two cells.");
                }

                adjacency.neighbor =
                    cell_id;
            }

            const Index cell_face_position =
                begin + local_face;

            if (topology.cell_faces[cell_face_position]
                != invalid_index) {

                throw_topology_build_error(
                    "cell-face connectivity position was assigned "
                    "more than once.");
            }

            topology.cell_faces[cell_face_position] =
                face_id;
        }
    }

    // A Face is created exactly when a new EdgeKey is inserted.
    if (face_ids.size() != topology.faces.size()) {
        throw_topology_build_error(
            "internal face-indexing inconsistency.");
    }

    // -------------------------------------------------------------------------
    // Attach physical boundary information.
    // -------------------------------------------------------------------------

    for (const BoundaryEdge& boundary_edge :
         raw_mesh.boundary_edges) {

        const EdgeKey edge_key =
            make_edge_key(
                boundary_edge.node_ids[0],
                boundary_edge.node_ids[1]);

        const auto iterator =
            face_ids.find(edge_key);

        if (iterator == face_ids.end()) {
            throw_topology_build_error(
                "a boundary edge does not correspond "
                "to any cell face.");
        }

        const Index face_id =
            iterator->second;

        const FaceAdjacency& adjacency =
            topology.face_adjacencies[face_id];

        if (adjacency.neighbor != invalid_index) {
            throw_topology_build_error(
                "face "
                + std::to_string(face_id)
                + " is internal but is marked as a boundary.");
        }

        if (topology.face_boundary_ids[face_id]
            != invalid_boundary_id) {

            throw_topology_build_error(
                "face "
                + std::to_string(face_id)
                + " received more than one boundary assignment.");
        }

        topology.face_boundary_ids[face_id] =
            boundary_edge.boundary_id;
    }

    return topology;
}

// -----------------------------------------------------------------------------
// Topology validation helpers
// -----------------------------------------------------------------------------

void validate_topology_storage(
    const RawMeshData& raw_mesh,
    const TopologyBuildData& topology)
{
    if (topology.faces.empty()) {
        throw_topology_validation_error(
            "no faces were constructed.");
    }

    if (topology.cell_faces.size()
        != raw_mesh.cell_nodes.size()) {

        throw_topology_validation_error(
            "cell_faces size must equal cell_nodes size.");
    }

    if (topology.face_adjacencies.size()
        != topology.faces.size()) {

        throw_topology_validation_error(
            "face_adjacencies size must equal faces size.");
    }

    if (topology.face_boundary_ids.size()
        != topology.faces.size()) {

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
    const Index begin =
        raw_mesh.cell_node_offsets[cell_id];

    const Index end =
        raw_mesh.cell_node_offsets[cell_id + 1];

    for (Index position = begin;
         position < end;
         ++position) {

        if (topology.cell_faces[position]
            == face_id) {

            return true;
        }
    }

    return false;
}

void validate_cell_face_connectivity(
    const RawMeshData& raw_mesh,
    const TopologyBuildData& topology)
{
    for (Index cell_id = 0;
         cell_id < raw_mesh.cell_types.size();
         ++cell_id) {

        const Index begin =
            raw_mesh.cell_node_offsets[cell_id];

        const Index end =
            raw_mesh.cell_node_offsets[cell_id + 1];

        const Index local_face_count =
            end - begin;

        for (Index local_face = 0;
             local_face < local_face_count;
             ++local_face) {

            const Index position =
                begin + local_face;

            const Index face_id =
                topology.cell_faces[position];

            if (face_id == invalid_index
                || face_id >= topology.faces.size()) {

                throw_topology_validation_error(
                    "cell "
                    + std::to_string(cell_id)
                    + " references an invalid face.");
            }

            // A cell must not reference the same face twice.
            //
            // Cells contain only 3 or 4 faces, so this small local
            // O(n^2) check is simpler and cheaper than allocating
            // a set for every cell.
            for (Index previous_position = begin;
                 previous_position < position;
                 ++previous_position) {

                if (topology.cell_faces[previous_position]
                    == face_id) {

                    throw_topology_validation_error(
                        "cell "
                        + std::to_string(cell_id)
                        + " references the same face more than once.");
                }
            }

            // cell -> face must agree with face -> cell.
            const FaceAdjacency& adjacency =
                topology.face_adjacencies[face_id];

            if (adjacency.owner != cell_id
                && adjacency.neighbor != cell_id) {

                throw_topology_validation_error(
                    "cell-to-face and face-to-cell "
                    "connectivities are inconsistent.");
            }

            // Verify that cell_faces[position] actually corresponds
            // to the expected local edge of this cell.
            const Index next_local_node =
                (local_face + 1)
                % local_face_count;

            const Index expected_node_0 =
                raw_mesh.cell_nodes[
                    begin + local_face];

            const Index expected_node_1 =
                raw_mesh.cell_nodes[
                    begin + next_local_node];

            const EdgeKey expected_edge =
                make_edge_key(
                    expected_node_0,
                    expected_node_1);

            const Face& face =
                topology.faces[face_id];

            const EdgeKey actual_edge =
                make_edge_key(
                    face.node_ids[0],
                    face.node_ids[1]);

            if (!(actual_edge == expected_edge)) {
                throw_topology_validation_error(
                    "cell "
                    + std::to_string(cell_id)
                    + " references a face with incorrect nodes.");
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

    for (Index face_id = 0;
         face_id < topology.faces.size();
         ++face_id) {

        const Face& face =
            topology.faces[face_id];

        const Index node_0 =
            face.node_ids[0];

        const Index node_1 =
            face.node_ids[1];

        if (node_0 >= raw_mesh.nodes.size()
            || node_1 >= raw_mesh.nodes.size()) {

            throw_topology_validation_error(
                "face "
                + std::to_string(face_id)
                + " references a node outside the nodes array.");
        }

        if (node_0 == node_1) {
            throw_topology_validation_error(
                "face "
                + std::to_string(face_id)
                + " references the same node twice.");
        }

        const FaceAdjacency& adjacency =
            topology.face_adjacencies[face_id];

        if (adjacency.owner
            >= raw_mesh.cell_types.size()) {

            throw_topology_validation_error(
                "face "
                + std::to_string(face_id)
                + " has an invalid owner cell.");
        }

        if (!cell_references_face(
                raw_mesh,
                topology,
                adjacency.owner,
                face_id)) {

            throw_topology_validation_error(
                "face "
                + std::to_string(face_id)
                + " is not referenced by its owner cell.");
        }

        const BoundaryId boundary_id =
            topology.face_boundary_ids[face_id];

        if (adjacency.neighbor == invalid_index) {
            // Boundary face.
            ++stats.boundary_face_count;

            if (boundary_id == invalid_boundary_id) {
                throw_topology_validation_error(
                    "external face "
                    + std::to_string(face_id)
                    + " has no boundary assignment.");
            }

            continue;
        }

        // Internal face.
        ++stats.internal_face_count;

        if (adjacency.neighbor
            >= raw_mesh.cell_types.size()) {

            throw_topology_validation_error(
                "face "
                + std::to_string(face_id)
                + " has an invalid neighbor cell.");
        }

        if (adjacency.owner
            == adjacency.neighbor) {

            throw_topology_validation_error(
                "face "
                + std::to_string(face_id)
                + " has identical owner and neighbor cells.");
        }

        if (boundary_id
            != invalid_boundary_id) {

            throw_topology_validation_error(
                "internal face "
                + std::to_string(face_id)
                + " has a boundary assignment.");
        }

        if (!cell_references_face(
                raw_mesh,
                topology,
                adjacency.neighbor,
                face_id)) {

            throw_topology_validation_error(
                "face "
                + std::to_string(face_id)
                + " is not referenced by its neighbor cell.");
        }
    }

    // Every raw boundary edge must correspond to exactly one
    // constructed external face.
    if (stats.boundary_face_count
        != raw_mesh.boundary_edges.size()) {

        throw_topology_validation_error(
            "number of external faces does not match "
            "the number of boundary edges.");
    }

    return stats;
}

void validate_unique_cell_neighbors(
    const TopologyBuildData& topology)
{
    std::vector<CellPair> neighboring_cell_pairs;

    neighboring_cell_pairs.reserve(
        topology.faces.size());

    for (const FaceAdjacency& adjacency :
         topology.face_adjacencies) {

        if (adjacency.neighbor
            == invalid_index) {

            continue;
        }

        neighboring_cell_pairs.push_back(
            make_cell_pair(
                adjacency.owner,
                adjacency.neighbor));
    }

    std::sort(
        neighboring_cell_pairs.begin(),
        neighboring_cell_pairs.end());

    const auto duplicate =
        std::adjacent_find(
            neighboring_cell_pairs.begin(),
            neighboring_cell_pairs.end());

    if (duplicate
        != neighboring_cell_pairs.end()) {

        throw_topology_validation_error(
            "cells "
            + std::to_string(duplicate->cell_0)
            + " and "
            + std::to_string(duplicate->cell_1)
            + " share more than one face.");
    }
}

[[nodiscard]]
TopologyStats validate_topology(
    const RawMeshData& raw_mesh,
    const TopologyBuildData& topology)
{
    validate_topology_storage(
        raw_mesh,
        topology);

    validate_cell_face_connectivity(
        raw_mesh,
        topology);

    const TopologyStats stats =
        validate_faces(
            raw_mesh,
            topology);

    validate_unique_cell_neighbors(
        topology);

    // Each internal face appears in two cells.
    // Each boundary face appears in one cell.
    //
    // In our 2D representation, cell_nodes.size() is also
    // the total number of local cell faces.
    if (raw_mesh.cell_nodes.size()
        != 2 * stats.internal_face_count
         + stats.boundary_face_count) {

        throw_topology_validation_error(
            "local face count is inconsistent with "
            "internal and boundary face counts.");
    }

    return stats;
}

} // namespace

// -----------------------------------------------------------------------------
// Public mesh construction
// -----------------------------------------------------------------------------

Mesh build_mesh(RawMeshData&& raw_mesh)
{
    // Raw data must remain intact until all construction and validation
    // steps have succeeded.
    validate_raw_mesh(raw_mesh);

    const auto topology_start =
        std::chrono::steady_clock::now();

    auto topology =
        build_topology(raw_mesh);

    const TopologyStats stats =
        validate_topology(
            raw_mesh,
            topology);

    const auto topology_end =
        std::chrono::steady_clock::now();

    const auto topology_elapsed =
        std::chrono::duration<double, std::milli>(
            topology_end - topology_start);

    std::cout
        << "[CFD] Topology: "
        << topology.faces.size()
        << " faces ("
        << stats.internal_face_count
        << " internal, "
        << stats.boundary_face_count
        << " boundary) ["
        << topology_elapsed.count()
        << " ms]\n";

    // -------------------------------------------------------------------------
    // All topology invariants are satisfied from this point onward.
    // Transfer ownership of persistent buffers into the final Mesh.
    // -------------------------------------------------------------------------

    Mesh mesh;

    mesh.nodes_ =
        std::move(raw_mesh.nodes);

    mesh.cell_types_ =
        std::move(raw_mesh.cell_types);

    mesh.cell_nodes_ =
        std::move(raw_mesh.cell_nodes);

    mesh.cell_node_offsets_ =
        std::move(raw_mesh.cell_node_offsets);

    mesh.boundary_groups_ =
        std::move(raw_mesh.boundary_groups);

    mesh.faces_ =
        std::move(topology.faces);

    mesh.cell_faces_ =
        std::move(topology.cell_faces);

    mesh.face_adjacencies_ =
        std::move(topology.face_adjacencies);

    mesh.face_boundary_ids_ =
        std::move(topology.face_boundary_ids);

    // raw_mesh.boundary_edges is intentionally not transferred.
    //
    // BoundaryEdge is an import-time representation. Once physical
    // boundary IDs have been attached to the final faces, these raw
    // edges are no longer needed.

    return mesh;
}

} // namespace cfd