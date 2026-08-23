#include "cfd/mesh/MeshBuilder.hpp"

#include "cfd/meshing/RawMeshData.hpp"
#include "cfd/meshing/RawMeshValidation.hpp"

#include "mesh_build/GeometryBuilder.hpp"
#include "mesh_build/GeometryValidation.hpp"
#include "mesh_build/TopologyBuilder.hpp"
#include "mesh_build/TopologyValidation.hpp"

#include <chrono>
#include <iostream>
#include <utility>

namespace cfd {

Mesh build_mesh(RawMeshData&& raw_mesh)
{
    validate_raw_mesh(raw_mesh);

    // -------------------------------------------------------------------------
    // Topology
    // -------------------------------------------------------------------------

    const auto topology_start { std::chrono::steady_clock::now() };

    auto topology { detail::build_topology(raw_mesh) };

    const detail::TopologyStats topology_stats { detail::validate_topology(raw_mesh, topology) };

    const auto topology_end { std::chrono::steady_clock::now() };

    const auto topology_elapsed { std::chrono::duration<double, std::milli>(topology_end - topology_start) };

    std::cout
        << "[CFD] Topology: "
        << topology.faces.size()
        << " faces ("
        << topology_stats.internal_face_count
        << " internal, "
        << topology_stats.boundary_face_count
        << " boundary) ["
        << topology_elapsed.count()
        << " ms]\n";

    // -------------------------------------------------------------------------
    // Geometry
    // -------------------------------------------------------------------------

    const auto geometry_start { std::chrono::steady_clock::now() };

    auto geometry { detail::build_geometry(raw_mesh, topology) };

    const detail::GeometryStats geometry_stats { detail::validate_geometry(raw_mesh, topology, geometry) };

    const auto geometry_end { std::chrono::steady_clock::now() };

    const auto geometry_elapsed { std::chrono::duration<double, std::milli>(geometry_end - geometry_start) };

    std::cout
        << "[CFD] Geometry: total area="
        << geometry_stats.total_cell_area
        << " ["
        << geometry_elapsed.count()
        << " ms]\n";

    std::cout
        << "      Cell area:   min="
        << geometry_stats.cell_areas.minimum
        << ", mean="
        << geometry_stats.cell_areas.mean
        << ", max="
        << geometry_stats.cell_areas.maximum
        << '\n';

    std::cout
        << "      Cell size:   min="
        << geometry_stats.cell_sizes.minimum
        << ", mean="
        << geometry_stats.cell_sizes.mean
        << ", max="
        << geometry_stats.cell_sizes.maximum
        << '\n';

    std::cout
        << "      Face length: min="
        << geometry_stats.face_lengths.minimum
        << ", mean="
        << geometry_stats.face_lengths.mean
        << ", max="
        << geometry_stats.face_lengths.maximum
        << '\n';

    if (geometry_stats.worst_quality_cell != invalid_index)
    {
        std::cout
            << "      Triangle q:  min="
            << geometry_stats.triangle_quality.minimum
            << ", mean="
            << geometry_stats.triangle_quality.mean
            << ", max="
            << geometry_stats.triangle_quality.maximum
            << ", worst cell="
            << geometry_stats.worst_quality_cell
            << '\n';
    }

    // -------------------------------------------------------------------------
    // Final Mesh
    // -------------------------------------------------------------------------

    Mesh mesh;

    mesh.nodes_ = std::move(raw_mesh.nodes);

    mesh.cell_types_ = std::move(raw_mesh.cell_types);
    mesh.cell_nodes_ = std::move(raw_mesh.cell_nodes);
    mesh.cell_node_offsets_ = std::move(raw_mesh.cell_node_offsets);

    mesh.boundary_groups_ = std::move(raw_mesh.boundary_groups);

    // Topology
    mesh.faces_ = std::move(topology.faces);
    mesh.cell_faces_ = std::move(topology.cell_faces);
    mesh.face_adjacencies_ = std::move(topology.face_adjacencies);
    mesh.face_boundary_ids_ = std::move(topology.face_boundary_ids);

    // Geometry
    mesh.cell_areas_ = std::move(geometry.cell_areas);
    mesh.cell_centers_ = std::move(geometry.cell_centers);

    mesh.face_centers_ = std::move(geometry.face_centers);
    mesh.face_lengths_ = std::move(geometry.face_lengths);
    mesh.face_area_vectors_ = std::move(geometry.face_area_vectors);
    return mesh;
}

} // namespace cfd