#include "cfd/meshing/RawMeshValidation.hpp"

#include "cfd/meshing/RawMeshData.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

namespace cfd
{

namespace
{

[[noreturn]]
void throw_raw_mesh_validation_error(const std::string &message)
{
    throw std::runtime_error("Raw mesh validation failed: " + message);
}

void validate_nodes(const RawMeshData &raw_mesh)
{
    if (raw_mesh.nodes.empty())
    {
        throw_raw_mesh_validation_error("nodes array is empty.");
    }

    for (Index node_id = 0; node_id < raw_mesh.nodes.size(); ++node_id)
    {
        const Node &node{raw_mesh.nodes[node_id]};

        if (!std::isfinite(node.x) || !std::isfinite(node.y))
        {
            throw_raw_mesh_validation_error("node " + std::to_string(node_id) + " contains a non-finite coordinate.");
        }
    }
}

void validate_cell_storage(const RawMeshData &raw_mesh)
{
    if (raw_mesh.cell_types.empty())
    {
        throw_raw_mesh_validation_error("cell_types array is empty.");
    }

    if (raw_mesh.cell_nodes.empty())
    {
        throw_raw_mesh_validation_error("cell_nodes array is empty.");
    }

    if (raw_mesh.cell_node_offsets.empty())
    {
        throw_raw_mesh_validation_error("cell_node_offsets array is empty.");
    }

    if (raw_mesh.cell_node_offsets.size() != raw_mesh.cell_types.size() + 1)
    {
        throw_raw_mesh_validation_error("cell_node_offsets size must equal the number of cells + 1.");
    }

    if (raw_mesh.cell_node_offsets.front() != 0)
    {
        throw_raw_mesh_validation_error("cell_node_offsets must start at 0.");
    }

    for (std::size_t offset_index = 1; offset_index < raw_mesh.cell_node_offsets.size(); ++offset_index)
    {
        const Index previous_offset{raw_mesh.cell_node_offsets[offset_index - 1]};
        const Index current_offset{raw_mesh.cell_node_offsets[offset_index]};

        if (current_offset <= previous_offset)
        {
            throw_raw_mesh_validation_error("cell_node_offsets must be strictly increasing.");
        }

        if (current_offset > raw_mesh.cell_nodes.size())
        {
            throw_raw_mesh_validation_error("cell_node_offsets contains an offset outside cell_nodes.");
        }
    }

    if (raw_mesh.cell_node_offsets.back() != raw_mesh.cell_nodes.size())
    {
        throw_raw_mesh_validation_error("last cell_node_offset must equal cell_nodes.size().");
    }
}

[[nodiscard]]
Index expected_cell_node_count(const CellType cell_type, const Index cell_id)
{
    switch (cell_type)
    {
    case CellType::Triangle:
        return 3;

    case CellType::Quadrilateral:
        return 4;
    }

    throw_raw_mesh_validation_error("cell " + std::to_string(cell_id) + " has an unsupported CellType.");
}

void validate_cells(const RawMeshData &raw_mesh)
{
    for (Index cell_id = 0; cell_id < raw_mesh.cell_types.size(); ++cell_id)
    {
        const Index cell_node_begin_offset{raw_mesh.cell_node_offsets[cell_id]};
        const Index cell_node_end_offset{raw_mesh.cell_node_offsets[cell_id + 1]};
        const Index node_count{cell_node_end_offset - cell_node_begin_offset};
        const Index required_node_count{expected_cell_node_count(raw_mesh.cell_types[cell_id], cell_id)};

        if (node_count != required_node_count)
        {
            throw_raw_mesh_validation_error("cell " + std::to_string(cell_id) + " has " + std::to_string(node_count) +
                                            " nodes, but its type requires " + std::to_string(required_node_count) +
                                            ".");
        }

        for (Index cell_node_position = cell_node_begin_offset; cell_node_position < cell_node_end_offset;
             ++cell_node_position)
        {
            const Index node_id{raw_mesh.cell_nodes[cell_node_position]};

            if (node_id >= raw_mesh.nodes.size())
            {
                throw_raw_mesh_validation_error("cell " + std::to_string(cell_id) + " references node " +
                                                std::to_string(node_id) + ", which is outside the nodes array.");
            }

            // Cells currently contain only three or four nodes. A bounded O(n^2)
            // local scan avoids allocating a temporary set for every cell.
            for (Index previous_cell_node_position = cell_node_begin_offset;
                 previous_cell_node_position < cell_node_position; ++previous_cell_node_position)
            {
                if (raw_mesh.cell_nodes[previous_cell_node_position] == node_id)
                {
                    throw_raw_mesh_validation_error("cell " + std::to_string(cell_id) +
                                                    " contains duplicate node indices.");
                }
            }
        }
    }
}

using BoundaryIdSet = std::unordered_set<BoundaryId>;

[[nodiscard]]
BoundaryIdSet validate_boundary_groups(const RawMeshData &raw_mesh)
{
    if (raw_mesh.boundary_groups.empty())
    {
        throw_raw_mesh_validation_error("boundary_groups array is empty.");
    }

    BoundaryIdSet boundary_ids;
    boundary_ids.reserve(raw_mesh.boundary_groups.size());

    // Non-owning views avoid copying group names. Their lifetime is safe here
    // because the strings remain owned and unmodified by raw_mesh throughout
    // this validation.
    std::unordered_set<std::string_view> boundary_names;
    boundary_names.reserve(raw_mesh.boundary_groups.size());

    for (Index boundary_group_index = 0; boundary_group_index < raw_mesh.boundary_groups.size(); ++boundary_group_index)
    {
        const BoundaryGroup &group{raw_mesh.boundary_groups[boundary_group_index]};
        if (group.id == invalid_boundary_id)
        {
            throw_raw_mesh_validation_error("boundary group \"" + group.name +
                                            "\" uses the reserved invalid boundary ID.");
        }

        if (group.name.empty())
        {
            throw_raw_mesh_validation_error("a boundary group has an empty name.");
        }

        const bool is_id_inserted{boundary_ids.insert(group.id).second};

        if (!is_id_inserted)
        {
            throw_raw_mesh_validation_error("boundary group ID " + std::to_string(group.id) + " is duplicated.");
        }
        if (group.id != boundary_group_index)
        {
            throw_raw_mesh_validation_error("boundary group ID " + std::to_string(group.id) +
                                            " does not match its zero-based position " +
                                            std::to_string(boundary_group_index) + ".");
        }

        const bool is_name_inserted{boundary_names.insert(group.name).second};

        if (!is_name_inserted)
        {
            throw_raw_mesh_validation_error("boundary group name \"" + group.name + "\" is duplicated.");
        }
    }

    return boundary_ids;
}

void validate_boundary_edges(const RawMeshData &raw_mesh, const BoundaryIdSet &valid_boundary_ids)
{
    if (raw_mesh.boundary_edges.empty())
    {
        throw_raw_mesh_validation_error("boundary_edges array is empty.");
    }

    std::unordered_set<BoundaryId> used_boundary_ids;
    used_boundary_ids.reserve(valid_boundary_ids.size());

    // Boundary edges are undirected. Canonicalizing each node pair before
    // sorting makes duplicate detection independent of the stored orientation.
    std::vector<std::array<Index, 2>> boundary_edge_keys;
    boundary_edge_keys.reserve(raw_mesh.boundary_edges.size());

    for (std::size_t boundary_edge_index = 0; boundary_edge_index < raw_mesh.boundary_edges.size();
         ++boundary_edge_index)
    {
        const BoundaryEdge &boundary_edge{raw_mesh.boundary_edges[boundary_edge_index]};
        const Index node_0_id{boundary_edge.node_ids[0]};
        const Index node_1_id{boundary_edge.node_ids[1]};

        if (node_0_id >= raw_mesh.nodes.size() || node_1_id >= raw_mesh.nodes.size())
        {
            throw_raw_mesh_validation_error("boundary edge " + std::to_string(boundary_edge_index) +
                                            " references a node outside the nodes array.");
        }

        if (node_0_id == node_1_id)
        {
            throw_raw_mesh_validation_error("boundary edge " + std::to_string(boundary_edge_index) +
                                            " references the same node twice.");
        }

        if (!valid_boundary_ids.contains(boundary_edge.boundary_id))
        {
            throw_raw_mesh_validation_error("boundary edge " + std::to_string(boundary_edge_index) +
                                            " references unknown boundary group ID " +
                                            std::to_string(boundary_edge.boundary_id) + ".");
        }

        used_boundary_ids.insert(boundary_edge.boundary_id);

        if (node_0_id < node_1_id)
        {
            boundary_edge_keys.push_back({node_0_id, node_1_id});
        }
        else
        {
            boundary_edge_keys.push_back({node_1_id, node_0_id});
        }
    }

    std::sort(boundary_edge_keys.begin(), boundary_edge_keys.end());

    const auto duplicate_edge_iterator{std::adjacent_find(boundary_edge_keys.begin(), boundary_edge_keys.end())};

    if (duplicate_edge_iterator != boundary_edge_keys.end())
    {
        throw_raw_mesh_validation_error("the same boundary edge appears more than once.");
    }

    for (const BoundaryGroup &group : raw_mesh.boundary_groups)
    {
        if (!used_boundary_ids.contains(group.id))
        {
            throw_raw_mesh_validation_error("boundary group \"" + group.name + "\" contains no boundary edge.");
        }
    }
}

} // namespace

void validate_raw_mesh(const RawMeshData &raw_mesh)
{
    // Validate storage-level invariants before traversing cell connectivity.
    validate_nodes(raw_mesh);
    validate_cell_storage(raw_mesh);
    validate_cells(raw_mesh);

    // Boundary edges may reference only IDs from validated boundary groups.
    const BoundaryIdSet valid_boundary_ids{validate_boundary_groups(raw_mesh)};

    validate_boundary_edges(raw_mesh, valid_boundary_ids);
}

} // namespace cfd