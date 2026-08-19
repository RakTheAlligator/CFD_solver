#include "cfd/mesh/MeshBuilder.hpp"

#include "cfd/meshing/RawMeshData.hpp"
#include "cfd/meshing/RawMeshValidation.hpp"

#include "mesh_build/TopologyBuilder.hpp"
#include "mesh_build/TopologyValidation.hpp"

#include <chrono>
#include <iostream>
#include <utility>

namespace cfd {

Mesh build_mesh(RawMeshData&& raw_mesh)
{
    validate_raw_mesh(raw_mesh);

    const auto topology_start =
        std::chrono::steady_clock::now();

    auto topology =
        detail::build_topology(raw_mesh);

    const detail::TopologyStats stats =
        detail::validate_topology(
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

    return mesh;
}

} // namespace cfd