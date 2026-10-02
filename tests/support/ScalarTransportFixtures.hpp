#pragma once

#include "cfd/field/FaceFluxField.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/meshing/RawMeshData.hpp"

#include <utility>
#include <vector>

namespace cfd::test
{

/// Two QUAD cells with a separate group for each boundary face, so affine
/// Dirichlet data can also be prescribed on the sheared variant.
[[nodiscard]]
inline RawMeshData make_transport_raw_mesh(const bool sheared = false)
{
    const double shear{sheared ? 0.35 : 0.0};
    RawMeshData raw;
    raw.nodes = {{0.0, 0.0}, {0.5, 0.0}, {1.0, 0.0}, {shear, 1.0}, {0.5 + shear, 1.0}, {1.0 + shear, 1.0}};
    raw.cell_types = {CellType::Quadrilateral, CellType::Quadrilateral};
    raw.cell_nodes = {0, 1, 4, 3, 1, 2, 5, 4};
    raw.cell_node_offsets = {0, 4, 8};
    raw.boundary_groups = {{0, "left"},         {1, "right"},     {2, "bottom_left"},
                           {3, "bottom_right"}, {4, "top_right"}, {5, "top_left"}};
    raw.boundary_edges = {{{3, 0}, 0}, {{2, 5}, 1}, {{0, 1}, 2}, {{1, 2}, 3}, {{5, 4}, 4}, {{4, 3}, 5}};
    return raw;
}

[[nodiscard]]
inline double transport_affine_value(const Point2 &point) noexcept
{
    return 2.0 * point.x - 3.0 * point.y + 1.5;
}

[[nodiscard]]
inline ScalarBoundaryConditions transport_affine_conditions(const Mesh &mesh)
{
    std::vector<ScalarBoundaryCondition> conditions(mesh.boundary_groups().size(),
                                                    {ScalarBoundaryConditionType::Dirichlet, 0.0});
    for (Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (mesh.face_adjacencies()[face_id].is_boundary())
        {
            conditions[mesh.face_boundary_ids()[face_id]].value = transport_affine_value(mesh.face_centers()[face_id]);
        }
    }
    return {conditions.size(), std::move(conditions)};
}

/// Unit interval end values, insulating horizontal walls; groups come from
/// make_transport_raw_mesh(), not from an assumed face ordering.
[[nodiscard]]
inline ScalarBoundaryConditions transport_channel_conditions()
{
    return {6,
            {{ScalarBoundaryConditionType::Dirichlet, 0.0},
             {ScalarBoundaryConditionType::Dirichlet, 1.0},
             {ScalarBoundaryConditionType::Neumann, 0.0},
             {ScalarBoundaryConditionType::Neumann, 0.0},
             {ScalarBoundaryConditionType::Neumann, 0.0},
             {ScalarBoundaryConditionType::Neumann, 0.0}}};
}

[[nodiscard]]
inline FaceFluxField transport_flux(const Mesh &mesh, const double speed = 1.0)
{
    FaceFluxField flux{mesh.face_count()};
    for (Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        flux[face_id] = speed * mesh.face_area_vectors()[face_id].x;
    }
    return flux;
}

} // namespace cfd::test
