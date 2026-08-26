#include "cfd/mesh/MeshBuilder.hpp"

#include "cfd/meshing/RawMeshData.hpp"
#include "cfd/meshing/RawMeshValidation.hpp"

#include "mesh_build/GeometryBuilder.hpp"
#include "mesh_build/GeometryValidation.hpp"
#include "mesh_build/TopologyBuilder.hpp"
#include "mesh_build/TopologyValidation.hpp"

#include <chrono>
#include <utility>

namespace cfd
{

MeshBuildResult build_mesh(RawMeshData &&raw_mesh)
{
    validate_raw_mesh(raw_mesh);

    // -------------------------------------------------------------------------
    // Topology
    // -------------------------------------------------------------------------

    const auto topology_start{std::chrono::steady_clock::now()};

    auto topology{detail::build_topology(raw_mesh)};
    detail::validate_topology(raw_mesh, topology);

    const auto topology_end{std::chrono::steady_clock::now()};

    // -------------------------------------------------------------------------
    // Geometry
    // -------------------------------------------------------------------------

    const auto geometry_start{std::chrono::steady_clock::now()};

    auto geometry{detail::build_geometry(raw_mesh, topology)};
    detail::validate_geometry(raw_mesh, topology, geometry);

    const auto geometry_end{std::chrono::steady_clock::now()};

    // -------------------------------------------------------------------------
    // Final mesh
    // -------------------------------------------------------------------------

    Mesh mesh;

    mesh.nodes_ = std::move(raw_mesh.nodes);

    mesh.cell_types_ = std::move(raw_mesh.cell_types);
    mesh.cell_nodes_ = std::move(raw_mesh.cell_nodes);
    mesh.cell_node_offsets_ = std::move(raw_mesh.cell_node_offsets);

    mesh.boundary_groups_ = std::move(raw_mesh.boundary_groups);

    mesh.faces_ = std::move(topology.faces);
    mesh.cell_faces_ = std::move(topology.cell_faces);
    mesh.face_adjacencies_ = std::move(topology.face_adjacencies);
    mesh.face_boundary_ids_ = std::move(topology.face_boundary_ids);

    mesh.cell_areas_ = std::move(geometry.cell_areas);
    mesh.cell_centers_ = std::move(geometry.cell_centers);
    mesh.cell_qualities_ = std::move(geometry.cell_qualities);
    mesh.face_centers_ = std::move(geometry.face_centers);
    mesh.face_lengths_ = std::move(geometry.face_lengths);
    mesh.face_area_vectors_ = std::move(geometry.face_area_vectors);

    MeshBuildTimings timings{
        .topology = topology_end - topology_start,
        .geometry = geometry_end - geometry_start,
    };

    return {
        .mesh = std::move(mesh),
        .timings = timings,
    };
}

} // namespace cfd