#pragma once

#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Types.hpp"
#include "cfd/meshing/RawMeshData.hpp"

#include <cmath>

namespace cfd::test
{

[[nodiscard]]
inline RawMeshData make_single_triangle_raw_mesh()
{
    constexpr BoundaryId boundary_id{0};

    RawMeshData raw_mesh;

    raw_mesh.nodes = {
        {0.0, 0.0},
        {1.0, 0.0},
        {0.0, 1.0},
    };

    raw_mesh.cell_types = {
        CellType::Triangle,
    };

    raw_mesh.cell_nodes = {
        0,
        1,
        2,
    };

    raw_mesh.cell_node_offsets = {
        0,
        3,
    };

    raw_mesh.boundary_groups = {
        {boundary_id, "wall"},
    };

    raw_mesh.boundary_edges = {
        {{0, 1}, boundary_id},
        {{1, 2}, boundary_id},
        {{2, 0}, boundary_id},
    };

    return raw_mesh;
}

[[nodiscard]]
inline RawMeshData make_equilateral_triangle_raw_mesh(const double origin_x, const double origin_y, const double side)
{
    constexpr BoundaryId boundary_id{0};

    const double height{std::sqrt(3.0) * side / 2.0};

    RawMeshData raw_mesh;

    raw_mesh.nodes = {
        {origin_x, origin_y},
        {origin_x + side, origin_y},
        {origin_x + side / 2.0, origin_y + height},
    };

    raw_mesh.cell_types = {
        CellType::Triangle,
    };

    raw_mesh.cell_nodes = {
        0,
        1,
        2,
    };

    raw_mesh.cell_node_offsets = {
        0,
        3,
    };

    raw_mesh.boundary_groups = {
        {boundary_id, "wall"},
    };

    raw_mesh.boundary_edges = {
        {{0, 1}, boundary_id},
        {{1, 2}, boundary_id},
        {{2, 0}, boundary_id},
    };

    return raw_mesh;
}

[[nodiscard]]
inline RawMeshData make_non_manifold_raw_mesh()
{
    constexpr BoundaryId boundary_id{0};

    RawMeshData raw_mesh;

    raw_mesh.nodes = {
        {0.0, 0.0}, {1.0, 0.0}, {0.5, 1.0}, {0.5, -1.0}, {0.5, 2.0},
    };

    raw_mesh.cell_types = {
        CellType::Triangle,
        CellType::Triangle,
        CellType::Triangle,
    };

    raw_mesh.cell_nodes = {
        0, 1, 2, 1, 0, 3, 0, 1, 4,
    };

    raw_mesh.cell_node_offsets = {
        0,
        3,
        6,
        9,
    };

    raw_mesh.boundary_groups = {
        {boundary_id, "wall"},
    };

    raw_mesh.boundary_edges = {
        {{1, 2}, boundary_id}, {{2, 0}, boundary_id}, {{0, 3}, boundary_id},
        {{3, 1}, boundary_id}, {{1, 4}, boundary_id}, {{4, 0}, boundary_id},
    };

    return raw_mesh;
}

} // namespace cfd::test