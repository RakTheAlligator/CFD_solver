#include "cfd/GmshMesher.hpp"

#include <gmsh.h>

#include <cstddef>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace cfd {

namespace {

constexpr int gmsh_line_2_nodes = 1;
constexpr int gmsh_triangle_3_nodes = 2;

constexpr BoundaryId inlet_boundary_id{0};
constexpr BoundaryId wall_boundary_id{1};
constexpr BoundaryId outlet_boundary_id{2};


class GmshSession {
public:
    GmshSession()
    {
        gmsh::initialize();
    }

    ~GmshSession()
    {
        gmsh::finalize();
    }

    GmshSession(const GmshSession&) = delete;
    GmshSession& operator=(const GmshSession&) = delete;
};


struct RectangleGmshTags {
    int surface{};
    int inlet_physical_group{};
    int wall_physical_group{};
    int outlet_physical_group{};
};


void validate(
    const RectangleGeometry& geometry,
    const MeshGenerationOptions& options)
{
    if (geometry.length <= 0.0) {
        throw std::invalid_argument(
            "Rectangle length must be positive.");
    }

    if (geometry.height <= 0.0) {
        throw std::invalid_argument(
            "Rectangle height must be positive.");
    }

    if (options.mesh_size <= 0.0) {
        throw std::invalid_argument(
            "Mesh size must be positive.");
    }

    if (options.cell_type != CellType::Triangle) {
        throw std::invalid_argument(
            "Only triangular meshes are implemented for now.");
    }
}


RectangleGmshTags create_rectangle(
    const RectangleGeometry& geometry,
    double mesh_size)
{
    const int p0 =
        gmsh::model::geo::addPoint(
            0.0, 0.0, 0.0, mesh_size);

    const int p1 =
        gmsh::model::geo::addPoint(
            geometry.length, 0.0, 0.0, mesh_size);

    const int p2 =
        gmsh::model::geo::addPoint(
            geometry.length, geometry.height, 0.0, mesh_size);

    const int p3 =
        gmsh::model::geo::addPoint(
            0.0, geometry.height, 0.0, mesh_size);


    const int bottom =
        gmsh::model::geo::addLine(p0, p1);

    const int right =
        gmsh::model::geo::addLine(p1, p2);

    const int top =
        gmsh::model::geo::addLine(p2, p3);

    const int left =
        gmsh::model::geo::addLine(p3, p0);


    const int contour =
        gmsh::model::geo::addCurveLoop(
            {bottom, right, top, left});

    const int surface =
        gmsh::model::geo::addPlaneSurface({contour});


    // Makes the entities created in the geometry kernel
    // available to the rest of the Gmsh model.
    gmsh::model::geo::synchronize();


    const int inlet_group =
        gmsh::model::addPhysicalGroup(
            1, {left});

    gmsh::model::setPhysicalName(
        1, inlet_group, "inlet");


    const int wall_group =
        gmsh::model::addPhysicalGroup(
            1, {bottom, top});

    gmsh::model::setPhysicalName(
        1, wall_group, "wall");


    const int outlet_group =
        gmsh::model::addPhysicalGroup(
            1, {right});

    gmsh::model::setPhysicalName(
        1, outlet_group, "outlet");


    const int fluid_group =
        gmsh::model::addPhysicalGroup(
            2, {surface});

    gmsh::model::setPhysicalName(
        2, fluid_group, "fluid");


    return {
        surface,
        inlet_group,
        wall_group,
        outlet_group
    };
}


Index local_node_index(
    const std::unordered_map<std::size_t, Index>& gmsh_to_local,
    std::size_t gmsh_node_tag)
{
    const auto iterator = gmsh_to_local.find(gmsh_node_tag);

    if (iterator == gmsh_to_local.end()) {
        throw std::runtime_error(
            "Unknown Gmsh node tag.");
    }

    return iterator->second;
}


std::unordered_map<std::size_t, Index>
extract_nodes(RawMeshData& raw_mesh)
{
    std::vector<std::size_t> node_tags;
    std::vector<double> coordinates;
    std::vector<double> parametric_coordinates;

    gmsh::model::mesh::getNodes(
        node_tags,
        coordinates,
        parametric_coordinates,
        -1,
        -1,
        false,
        false
    );

    if (coordinates.size() != 3 * node_tags.size()) {
        throw std::runtime_error(
            "Invalid node coordinate data returned by Gmsh.");
    }


    raw_mesh.nodes.reserve(node_tags.size());

    std::unordered_map<std::size_t, Index> gmsh_to_local;
    gmsh_to_local.reserve(node_tags.size());


    for (Index i_node = 0; i_node < node_tags.size(); ++i_node) {

        const std::size_t gmsh_tag = node_tags[i_node];

        gmsh_to_local.emplace(gmsh_tag, i_node);

        raw_mesh.nodes.push_back({
            coordinates[3 * i_node],
            coordinates[3 * i_node + 1]
        });
    }


    return gmsh_to_local;
}


void extract_triangles(
    RawMeshData& raw_mesh,
    const std::unordered_map<std::size_t, Index>& gmsh_to_local,
    int surface_tag)
{
    std::vector<int> element_types;

    std::vector<std::vector<std::size_t>> element_tags;
    std::vector<std::vector<std::size_t>> element_node_tags;

    gmsh::model::mesh::getElements(
        element_types,
        element_tags,
        element_node_tags,
        2,
        surface_tag
    );


    raw_mesh.cell_node_offsets.clear();
    raw_mesh.cell_node_offsets.push_back(0);


    for (std::size_t i_type = 0;
         i_type < element_types.size();
         ++i_type) {

        if (element_types[i_type] != gmsh_triangle_3_nodes) {
            throw std::runtime_error(
                "Gmsh generated a non-triangular 2D element.");
        }


        const auto& nodes = element_node_tags[i_type];

        if (nodes.size() % 3 != 0) {
            throw std::runtime_error(
                "Invalid triangle connectivity returned by Gmsh.");
        }


        for (std::size_t i = 0; i < nodes.size(); i += 3) {

            raw_mesh.cell_types.push_back(CellType::Triangle);

            raw_mesh.cell_nodes.push_back(
                local_node_index(gmsh_to_local, nodes[i]));

            raw_mesh.cell_nodes.push_back(
                local_node_index(gmsh_to_local, nodes[i + 1]));

            raw_mesh.cell_nodes.push_back(
                local_node_index(gmsh_to_local, nodes[i + 2]));

            raw_mesh.cell_node_offsets.push_back(
                raw_mesh.cell_nodes.size());
        }
    }
}


void extract_boundary_edges(
    RawMeshData& raw_mesh,
    const std::unordered_map<std::size_t, Index>& gmsh_to_local,
    int physical_group_tag,
    BoundaryId boundary_id)
{
    std::vector<int> curve_tags;

    gmsh::model::getEntitiesForPhysicalGroup(
        1,
        physical_group_tag,
        curve_tags
    );


    for (const int curve_tag : curve_tags) {

        std::vector<int> element_types;

        std::vector<std::vector<std::size_t>> element_tags;
        std::vector<std::vector<std::size_t>> element_node_tags;

        gmsh::model::mesh::getElements(
            element_types,
            element_tags,
            element_node_tags,
            1,
            curve_tag
        );


        for (std::size_t i_type = 0;
             i_type < element_types.size();
             ++i_type) {

            if (element_types[i_type] != gmsh_line_2_nodes) {
                throw std::runtime_error(
                    "Gmsh generated a non-linear boundary element.");
            }


            const auto& nodes = element_node_tags[i_type];

            if (nodes.size() % 2 != 0) {
                throw std::runtime_error(
                    "Invalid boundary connectivity returned by Gmsh.");
            }


            for (std::size_t i = 0; i < nodes.size(); i += 2) {

                raw_mesh.boundary_edges.push_back({
                    {
                        local_node_index(gmsh_to_local, nodes[i]),
                        local_node_index(gmsh_to_local, nodes[i + 1])
                    },
                    boundary_id
                });
            }
        }
    }
}

} // namespace


RawMeshData generate_mesh(
    const RectangleGeometry& geometry,
    const MeshGenerationOptions& options)
{
    validate(geometry, options);

    GmshSession gmsh_session;

    gmsh::model::add("rectangle");


    const RectangleGmshTags tags =
        create_rectangle(
            geometry,
            options.mesh_size);


    gmsh::model::mesh::generate(2);


    RawMeshData raw_mesh;

    raw_mesh.boundary_groups = {
        {inlet_boundary_id,  "inlet"},
        {wall_boundary_id,   "wall"},
        {outlet_boundary_id, "outlet"}
    };


    const auto gmsh_to_local =
        extract_nodes(raw_mesh);


    extract_triangles(
        raw_mesh,
        gmsh_to_local,
        tags.surface);


    extract_boundary_edges(
        raw_mesh,
        gmsh_to_local,
        tags.inlet_physical_group,
        inlet_boundary_id);

    extract_boundary_edges(
        raw_mesh,
        gmsh_to_local,
        tags.wall_physical_group,
        wall_boundary_id);

    extract_boundary_edges(
        raw_mesh,
        gmsh_to_local,
        tags.outlet_physical_group,
        outlet_boundary_id);


    return raw_mesh;
}

} // namespace cfd