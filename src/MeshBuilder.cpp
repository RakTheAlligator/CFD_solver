#include "cfd/mesh/MeshBuilder.hpp"

#include "cfd/meshing/RawMeshData.hpp"
#include "cfd/meshing/RawMeshValidation.hpp"

#include "mesh_build/GeometryBuilder.hpp"
#include "mesh_build/GeometryValidation.hpp"
#include "mesh_build/TopologyBuilder.hpp"
#include "mesh_build/TopologyValidation.hpp"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <string_view>
#include <utility>

namespace cfd
{

namespace
{

void print_statistics_row(const std::string_view quantity, const std::string_view unit,
                          const detail::ScalarStats &statistics)
{
    std::cout << "  " << std::left << std::setw(22) << quantity << std::setw(8) << unit << std::right
              << std::defaultfloat << std::setprecision(5) << std::setw(12) << statistics.minimum << std::setw(12)
              << statistics.mean << std::setw(12) << statistics.maximum << '\n';
}

} // namespace

Mesh build_mesh(RawMeshData &&raw_mesh)
{
    validate_raw_mesh(raw_mesh);

    // -------------------------------------------------------------------------
    // Topology
    // -------------------------------------------------------------------------

    const auto topology_start{std::chrono::steady_clock::now()};

    auto topology{detail::build_topology(raw_mesh)};
    const detail::TopologyStats topology_stats{detail::validate_topology(raw_mesh, topology)};

    const auto topology_end{std::chrono::steady_clock::now()};

    const auto topology_elapsed{std::chrono::duration<double, std::milli>(topology_end - topology_start)};

    std::cout << "\n[Mesh topology]\n"
              << "  Faces             : " << topology.faces.size() << '\n'
              << "    Internal         : " << topology_stats.internal_face_count << '\n'
              << "    Boundary         : " << topology_stats.boundary_face_count << '\n'
              << std::fixed << std::setprecision(2) << "  Time              : " << topology_elapsed.count() << " ms\n";

    // -------------------------------------------------------------------------
    // Geometry
    // -------------------------------------------------------------------------

    const auto geometry_start{std::chrono::steady_clock::now()};

    auto geometry{detail::build_geometry(raw_mesh, topology)};
    const detail::GeometryStats geometry_stats{detail::validate_geometry(raw_mesh, topology, geometry)};

    const auto geometry_end{std::chrono::steady_clock::now()};

    const auto geometry_elapsed{std::chrono::duration<double, std::milli>(geometry_end - geometry_start)};

    std::cout << "\n[Mesh geometry]\n"
              << std::fixed << std::setprecision(4) << "  Total area        : " << geometry_stats.total_cell_area
              << " m^2\n"
              << std::setprecision(2) << "  Time              : " << geometry_elapsed.count() << " ms\n";

    std::cout << "\n  " << std::left << std::setw(22) << "Quantity" << std::setw(8) << "Unit" << std::right
              << std::setw(12) << "min" << std::setw(12) << "mean" << std::setw(12) << "max" << '\n';

    std::cout << "  ------------------------------------------------------------------\n";

    print_statistics_row("Cell area", "m^2", geometry_stats.cell_areas);
    print_statistics_row("Cell size", "m", geometry_stats.cell_sizes);
    print_statistics_row("Face length", "m", geometry_stats.face_lengths);

    if (geometry_stats.worst_quality_cell != invalid_index)
    {
        print_statistics_row("Triangle quality", "-", geometry_stats.triangle_quality);

        std::cout << "\n  Worst triangle    : cell " << geometry_stats.worst_quality_cell << '\n';
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