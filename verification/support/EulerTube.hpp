#pragma once

#include "cfd/compressible/EulerBoundaryConditions.hpp"
#include "cfd/meshing/RawMeshData.hpp"

#include <cstdint>
#include <stdexcept>

namespace cfd::verification
{

enum class EulerTubeCells : std::uint8_t
{
    Quad,
    Triangle,
    Mixed
};

/// Deterministic straight unit-length/unit-radius tube, independent of Gmsh.
/// Boundary IDs: axis=0, outer wall=1, left=2, right=3.
inline RawMeshData make_euler_tube(const Index nx, const Index nr, const EulerTubeCells type = EulerTubeCells::Quad)
{
    if (nx == 0 || nr == 0)
        throw std::invalid_argument("Euler verification tube requires positive grid sizes.");
    RawMeshData raw;
    raw.boundary_groups = {{0, "axis"}, {1, "wall"}, {2, "left"}, {3, "right"}};
    for (Index j = 0; j <= nr; ++j)
        for (Index i = 0; i <= nx; ++i)
            raw.nodes.push_back(
                {static_cast<double>(i) / static_cast<double>(nx), static_cast<double>(j) / static_cast<double>(nr)});
    raw.cell_node_offsets.push_back(0);
    for (Index j = 0; j < nr; ++j)
    {
        for (Index i = 0; i < nx; ++i)
        {
            const Index a{j * (nx + 1) + i}, b{a + 1}, d{a + nx + 1}, c{d + 1};
            if (type == EulerTubeCells::Triangle || (type == EulerTubeCells::Mixed && (i + j) % 2 == 0))
            {
                raw.cell_types.push_back(CellType::Triangle);
                raw.cell_nodes.insert(raw.cell_nodes.end(), {a, b, c});
                raw.cell_node_offsets.push_back(raw.cell_nodes.size());
                raw.cell_types.push_back(CellType::Triangle);
                raw.cell_nodes.insert(raw.cell_nodes.end(), {a, c, d});
            }
            else
            {
                raw.cell_types.push_back(CellType::Quadrilateral);
                raw.cell_nodes.insert(raw.cell_nodes.end(), {a, b, c, d});
            }
            raw.cell_node_offsets.push_back(raw.cell_nodes.size());
        }
    }
    for (Index i = 0; i < nx; ++i)
    {
        raw.boundary_edges.push_back({{i, i + 1}, 0});
        raw.boundary_edges.push_back({{nr * (nx + 1) + i, nr * (nx + 1) + i + 1}, 1});
    }
    for (Index j = 0; j < nr; ++j)
    {
        raw.boundary_edges.push_back({{j * (nx + 1), (j + 1) * (nx + 1)}, 2});
        raw.boundary_edges.push_back({{j * (nx + 1) + nx, (j + 1) * (nx + 1) + nx}, 3});
    }
    return raw;
}

inline EulerBoundaryConditions euler_tube_conditions()
{
    using Type = EulerBoundaryConditionType;
    return {4, {Type::Axis, Type::SlipWall, Type::Transmissive, Type::Transmissive}};
}

} // namespace cfd::verification
