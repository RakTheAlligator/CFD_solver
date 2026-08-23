#include "cfd/meshing/RawMeshValidation.hpp"

#include "cfd/meshing/RawMeshData.hpp"

#include <algorithm>
#include <array>
#include <cmath>
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
void throw_validation_error(const std::string &message)
{
    throw std::runtime_error("Raw mesh validation failed: " + message);
}

void validate_nodes(const RawMeshData &raw_mesh)
{
    if (raw_mesh.nodes.empty())
    {
        throw_validation_error("nodes array is empty.");
    }

    for (Index node_id = 0; node_id < raw_mesh.nodes.size(); ++node_id)
    {
        const Node &node = raw_mesh.nodes[node_id];

        if (!std::isfinite(node.x) || !std::isfinite(node.y))
        {
            throw_validation_error("node " + std::to_string(node_id) + " contains a non-finite coordinate.");
        }
    }
}

void validate_cell_storage(const RawMeshData &raw_mesh)
{
    if (raw_mesh.cell_types.empty())
    {
        throw_validation_error("cell_types array is empty.");
    }
    if (raw_mesh.cell_nodes.empty())
    {
        throw_validation_error("cell_nodes array is empty.");
    }
    if (raw_mesh.cell_node_offsets.empty())
    {
        throw_validation_error("cell_node_offsets array is empty.");
    }
    if (raw_mesh.cell_node_offsets.size() != raw_mesh.cell_types.size() + 1)
    {
        throw_validation_error("cell_node_offsets size must equal the number of cells + 1.");
    }
    if (raw_mesh.cell_node_offsets.front() != 0)
    {
        throw_validation_error("cell_node_offsets must start at 0.");
    }

    for (std::size_t i = 1; i < raw_mesh.cell_node_offsets.size(); ++i)
    {
        const Index previous{raw_mesh.cell_node_offsets[i - 1]};
        const Index current{raw_mesh.cell_node_offsets[i]};

        if (current <= previous)
        {
            throw_validation_error("cell_node_offsets must be strictly increasing.");
        }
        if (current > raw_mesh.cell_nodes.size())
        {
            throw_validation_error("cell_node_offsets contains an offset outside cell_nodes.");
        }
    }

    if (raw_mesh.cell_node_offsets.back() != raw_mesh.cell_nodes.size())
    {
        throw_validation_error("last cell_node_offset must equal cell_nodes.size().");
    }
}

[[nodiscard]]
Index expected_node_count(const CellType cell_type, const Index cell_id)
{
    switch (cell_type)
    {
    case CellType::Triangle:
        return 3;

    case CellType::Quadrilateral:
        return 4;
    }

    throw_validation_error("cell " + std::to_string(cell_id) + " has an unsupported CellType.");
}

void validate_cells(const RawMeshData &raw_mesh)
{
    for (Index cell_id = 0; cell_id < raw_mesh.cell_types.size(); ++cell_id)
    {
        const Index begin{raw_mesh.cell_node_offsets[cell_id]};
        const Index end{raw_mesh.cell_node_offsets[cell_id + 1]};

        const Index node_count{end - begin};
        const Index expected_count{expected_node_count(raw_mesh.cell_types[cell_id], cell_id)};

        if (node_count != expected_count)
        {
            throw_validation_error("cell " + std::to_string(cell_id) + " has " + std::to_string(node_count) +
                                   " nodes, but its type requires " + std::to_string(expected_count) + ".");
        }

        for (Index position = begin; position < end; ++position)
        {
            const Index node_id{raw_mesh.cell_nodes[position]};

            if (node_id >= raw_mesh.nodes.size())
            {
                throw_validation_error("cell " + std::to_string(cell_id) + " references node " +
                                       std::to_string(node_id) + ", which is outside the nodes array.");
            }

            // Cells contain only 3 or 4 nodes, so a small
            // O(n^2) local search is simpler and cheaper than
            // allocating a set for every cell.
            for (Index previous_position = begin; previous_position < position; ++previous_position)
            {
                if (raw_mesh.cell_nodes[previous_position] == node_id)
                {
                    throw_validation_error("cell " + std::to_string(cell_id) + " contains duplicate node indices.");
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
        throw_validation_error("boundary_groups array is empty.");
    }

    BoundaryIdSet boundary_ids;
    boundary_ids.reserve(raw_mesh.boundary_groups.size());

    std::unordered_set<std::string_view> boundary_names;
    boundary_names.reserve(raw_mesh.boundary_groups.size());

    for (const BoundaryGroup &group : raw_mesh.boundary_groups)
    {
        if (group.id == invalid_boundary_id)
        {
            throw_validation_error("boundary group \"" + group.name + "\" uses the reserved invalid boundary ID.");
        }
        if (group.name.empty())
        {
            throw_validation_error("a boundary group has an empty name.");
        }

        const bool id_inserted{boundary_ids.insert(group.id).second};

        if (!id_inserted)
        {
            throw_validation_error("boundary group ID " + std::to_string(group.id) + " is duplicated.");
        }

        const bool name_inserted{boundary_names.insert(group.name).second};

        if (!name_inserted)
        {
            throw_validation_error("boundary group name \"" + group.name + "\" is duplicated.");
        }
    }

    return boundary_ids;
}

void validate_boundary_edges(const RawMeshData &raw_mesh, const BoundaryIdSet &valid_boundary_ids)
{
    if (raw_mesh.boundary_edges.empty())
    {
        throw_validation_error("boundary_edges array is empty.");
    }

    std::unordered_set<BoundaryId> used_boundary_ids;
    used_boundary_ids.reserve(valid_boundary_ids.size());

    // Canonical edge representation:
    //
    // {3, 8} and {8, 3} represent the same undirected edge.
    std::vector<std::array<Index, 2>> boundary_edge_keys;
    boundary_edge_keys.reserve(raw_mesh.boundary_edges.size());

    for (std::size_t edge_index = 0; edge_index < raw_mesh.boundary_edges.size(); ++edge_index)
    {

        const BoundaryEdge &edge{raw_mesh.boundary_edges[edge_index]};

        const Index node_a{edge.node_ids[0]};
        const Index node_b{edge.node_ids[1]};

        if (node_a >= raw_mesh.nodes.size() || node_b >= raw_mesh.nodes.size())
        {
            throw_validation_error("boundary edge " + std::to_string(edge_index) +
                                   " references a node outside the nodes array.");
        }
        if (node_a == node_b)
        {
            throw_validation_error("boundary edge " + std::to_string(edge_index) + " references the same node twice.");
        }
        if (!valid_boundary_ids.contains(edge.boundary_id))
        {
            throw_validation_error("boundary edge " + std::to_string(edge_index) +
                                   " references unknown boundary group ID " + std::to_string(edge.boundary_id) + ".");
        }

        used_boundary_ids.insert(edge.boundary_id);

        if (node_a < node_b)
        {
            boundary_edge_keys.push_back({node_a, node_b});
        }
        else
        {
            boundary_edge_keys.push_back({node_b, node_a});
        }
    }

    std::sort(boundary_edge_keys.begin(), boundary_edge_keys.end());

    const auto duplicate_edge{std::adjacent_find(boundary_edge_keys.begin(), boundary_edge_keys.end())};

    if (duplicate_edge != boundary_edge_keys.end())
    {
        throw_validation_error("the same boundary edge appears more than once.");
    }

    for (const BoundaryGroup &group : raw_mesh.boundary_groups)
    {
        if (!used_boundary_ids.contains(group.id))
        {
            throw_validation_error("boundary group \"" + group.name + "\" contains no boundary edge.");
        }
    }
}
} // namespace

void validate_raw_mesh(const RawMeshData &raw_mesh)
{
    validate_nodes(raw_mesh);
    validate_cell_storage(raw_mesh);
    validate_cells(raw_mesh);

    const BoundaryIdSet boundary_ids{validate_boundary_groups(raw_mesh)};

    validate_boundary_edges(raw_mesh, boundary_ids);
}
} // namespace cfd