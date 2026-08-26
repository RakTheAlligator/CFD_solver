#include "cfd/meshing/GmshMesher.hpp"

#include <gmsh.h>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace cfd
{

namespace
{

constexpr int gmsh_line_2_nodes = 1;
constexpr int gmsh_triangle_3_nodes = 2;
constexpr int gmsh_quadrilateral_4_nodes = 3;

constexpr BoundaryId inlet_boundary_id{0};
constexpr BoundaryId wall_boundary_id{1};
constexpr BoundaryId outlet_boundary_id{2};

class GmshSession
{
  public:
    GmshSession()
    {
        gmsh::initialize();
        gmsh::option::setNumber("General.Terminal", 0);
    }

    ~GmshSession()
    {
        gmsh::finalize();
    }

    GmshSession(const GmshSession &) = delete;
    GmshSession &operator=(const GmshSession &) = delete;

    GmshSession(GmshSession &&) = delete;
    GmshSession &operator=(GmshSession &&) = delete;
};

struct RectangleGmshTags
{
    int surface{};
    int inlet_physical_group{};
    int wall_physical_group{};
    int outlet_physical_group{};
};

struct CellExtractionInfo
{
    int gmsh_element_type{};
    Index node_count{};
    CellType cell_type{};
};

[[nodiscard]]
CellExtractionInfo cell_extraction_info(const CellType cell_type)
{
    switch (cell_type)
    {
    case CellType::Triangle:
        return {
            .gmsh_element_type = gmsh_triangle_3_nodes,
            .node_count = 3,
            .cell_type = CellType::Triangle,
        };

    case CellType::Quadrilateral:
        return {
            .gmsh_element_type = gmsh_quadrilateral_4_nodes,
            .node_count = 4,
            .cell_type = CellType::Quadrilateral,
        };
    }

    throw std::invalid_argument("Unsupported cell type.");
}

void validate(const RectangleGeometry &geometry, const MeshGenerationOptions &options)
{
    if (!std::isfinite(geometry.length) || !(geometry.length > 0.0))
    {
        throw std::invalid_argument("Rectangle length must be finite and positive.");
    }

    if (!std::isfinite(geometry.height) || !(geometry.height > 0.0))
    {
        throw std::invalid_argument("Rectangle height must be finite and positive.");
    }

    if (!std::isfinite(options.mesh_size) || !(options.mesh_size > 0.0))
    {
        throw std::invalid_argument("Mesh size must be finite and positive.");
    }

    static_cast<void>(cell_extraction_info(options.cell_type));
}

RectangleGmshTags create_rectangle(const RectangleGeometry &geometry, const double mesh_size)
{
    const int p0{gmsh::model::geo::addPoint(0.0, 0.0, 0.0, mesh_size)};

    const int p1{gmsh::model::geo::addPoint(geometry.length, 0.0, 0.0, mesh_size)};

    const int p2{gmsh::model::geo::addPoint(geometry.length, geometry.height, 0.0, mesh_size)};

    const int p3{gmsh::model::geo::addPoint(0.0, geometry.height, 0.0, mesh_size)};

    const int bottom{gmsh::model::geo::addLine(p0, p1)};

    const int right{gmsh::model::geo::addLine(p1, p2)};

    const int top{gmsh::model::geo::addLine(p2, p3)};

    const int left{gmsh::model::geo::addLine(p3, p0)};

    const int contour{gmsh::model::geo::addCurveLoop({bottom, right, top, left})};

    const int surface{gmsh::model::geo::addPlaneSurface({contour})};

    // Make the entities created in the geometry kernel
    // available to the rest of the Gmsh model.
    gmsh::model::geo::synchronize();

    const int inlet_group{gmsh::model::addPhysicalGroup(1, {left})};

    gmsh::model::setPhysicalName(1, inlet_group, "inlet");

    const int wall_group{gmsh::model::addPhysicalGroup(1, {bottom, top})};

    gmsh::model::setPhysicalName(1, wall_group, "wall");

    const int outlet_group{gmsh::model::addPhysicalGroup(1, {right})};

    gmsh::model::setPhysicalName(1, outlet_group, "outlet");

    const int fluid_group{gmsh::model::addPhysicalGroup(2, {surface})};

    gmsh::model::setPhysicalName(2, fluid_group, "fluid");

    return {
        surface,
        inlet_group,
        wall_group,
        outlet_group,
    };
}

void configure_surface_mesh(const CellType cell_type, const int surface_tag)
{
    switch (cell_type)
    {
    case CellType::Triangle:
        return;

    case CellType::Quadrilateral:
        // Blossom recombination:
        // generate a triangular background mesh, then recombine
        // pairs of triangles into quadrilateral elements.
        gmsh::option::setNumber("Mesh.RecombinationAlgorithm", 1);

        gmsh::model::mesh::setRecombine(2, surface_tag);

        return;
    }

    throw std::invalid_argument("Unsupported cell type.");
}

[[nodiscard]]
Index local_node_index(const std::unordered_map<std::size_t, Index> &gmsh_to_local, const std::size_t gmsh_node_tag)
{
    const auto iterator{gmsh_to_local.find(gmsh_node_tag)};

    if (iterator == gmsh_to_local.end())
    {
        throw std::runtime_error("Unknown Gmsh node tag.");
    }

    return iterator->second;
}

[[nodiscard]]
std::unordered_map<std::size_t, Index> extract_nodes(RawMeshData &raw_mesh)
{
    std::vector<std::size_t> node_tags;
    std::vector<double> coordinates;
    std::vector<double> parametric_coordinates;

    gmsh::model::mesh::getNodes(node_tags, coordinates, parametric_coordinates, -1, -1, false, false);

    if (coordinates.size() != 3 * node_tags.size())
    {
        throw std::runtime_error("Invalid node coordinate data returned by Gmsh.");
    }

    raw_mesh.nodes.reserve(node_tags.size());

    std::unordered_map<std::size_t, Index> gmsh_to_local;
    gmsh_to_local.reserve(node_tags.size());

    for (Index node_id = 0; node_id < node_tags.size(); ++node_id)
    {
        const std::size_t gmsh_tag{node_tags[node_id]};

        gmsh_to_local.emplace(gmsh_tag, node_id);

        raw_mesh.nodes.push_back({
            coordinates[3 * node_id],
            coordinates[3 * node_id + 1],
        });
    }

    return gmsh_to_local;
}

void extract_cells(RawMeshData &raw_mesh, const std::unordered_map<std::size_t, Index> &gmsh_to_local,
                   const int surface_tag, const CellType requested_cell_type)
{
    std::vector<int> element_types;
    std::vector<std::vector<std::size_t>> element_tags;
    std::vector<std::vector<std::size_t>> element_node_tags;

    gmsh::model::mesh::getElements(element_types, element_tags, element_node_tags, 2, surface_tag);

    const CellExtractionInfo extraction{cell_extraction_info(requested_cell_type)};

    Index total_cell_count{};
    Index total_connectivity_size{};

    for (std::size_t type_index = 0; type_index < element_types.size(); ++type_index)
    {
        if (element_types[type_index] != extraction.gmsh_element_type)
        {
            throw std::runtime_error("Gmsh generated a 2D element type that does not "
                                     "match the requested cell type.");
        }

        const auto &nodes{element_node_tags[type_index]};

        if (nodes.size() % extraction.node_count != 0)
        {
            throw std::runtime_error("Invalid cell connectivity returned by Gmsh.");
        }

        total_cell_count += nodes.size() / extraction.node_count;

        total_connectivity_size += nodes.size();
    }

    raw_mesh.cell_types.reserve(total_cell_count);
    raw_mesh.cell_nodes.reserve(total_connectivity_size);
    raw_mesh.cell_node_offsets.reserve(total_cell_count + 1);

    raw_mesh.cell_node_offsets.push_back(0);

    for (std::size_t type_index = 0; type_index < element_types.size(); ++type_index)
    {
        const auto &nodes{element_node_tags[type_index]};

        for (std::size_t first_node = 0; first_node < nodes.size(); first_node += extraction.node_count)
        {
            raw_mesh.cell_types.push_back(extraction.cell_type);

            for (Index local_node = 0; local_node < extraction.node_count; ++local_node)
            {
                raw_mesh.cell_nodes.push_back(local_node_index(gmsh_to_local, nodes[first_node + local_node]));
            }

            raw_mesh.cell_node_offsets.push_back(raw_mesh.cell_nodes.size());
        }
    }
}

void extract_boundary_edges(RawMeshData &raw_mesh, const std::unordered_map<std::size_t, Index> &gmsh_to_local,
                            const int physical_group_tag, const BoundaryId boundary_id)
{
    std::vector<int> curve_tags;

    gmsh::model::getEntitiesForPhysicalGroup(1, physical_group_tag, curve_tags);

    for (const int curve_tag : curve_tags)
    {
        std::vector<int> element_types;
        std::vector<std::vector<std::size_t>> element_tags;
        std::vector<std::vector<std::size_t>> element_node_tags;

        gmsh::model::mesh::getElements(element_types, element_tags, element_node_tags, 1, curve_tag);

        for (std::size_t type_index = 0; type_index < element_types.size(); ++type_index)
        {
            if (element_types[type_index] != gmsh_line_2_nodes)
            {
                throw std::runtime_error("Gmsh generated a non-linear boundary element.");
            }

            const auto &nodes{element_node_tags[type_index]};

            if (nodes.size() % 2 != 0)
            {
                throw std::runtime_error("Invalid boundary connectivity returned by Gmsh.");
            }

            for (std::size_t node_index = 0; node_index < nodes.size(); node_index += 2)
            {
                raw_mesh.boundary_edges.push_back({
                    {
                        local_node_index(gmsh_to_local, nodes[node_index]),
                        local_node_index(gmsh_to_local, nodes[node_index + 1]),
                    },
                    boundary_id,
                });
            }
        }
    }
}

} // namespace

RawMeshData generate_mesh(const RectangleGeometry &geometry, const MeshGenerationOptions &options)
{
    validate(geometry, options);

    GmshSession gmsh_session;

    gmsh::model::add("rectangle");

    const RectangleGmshTags tags{create_rectangle(geometry, options.mesh_size)};

    configure_surface_mesh(options.cell_type, tags.surface);

    gmsh::model::mesh::generate(2);

    RawMeshData raw_mesh;

    raw_mesh.boundary_groups = {
        {inlet_boundary_id, "inlet"},
        {wall_boundary_id, "wall"},
        {outlet_boundary_id, "outlet"},
    };

    const auto gmsh_to_local{extract_nodes(raw_mesh)};

    extract_cells(raw_mesh, gmsh_to_local, tags.surface, options.cell_type);

    extract_boundary_edges(raw_mesh, gmsh_to_local, tags.inlet_physical_group, inlet_boundary_id);

    extract_boundary_edges(raw_mesh, gmsh_to_local, tags.wall_physical_group, wall_boundary_id);

    extract_boundary_edges(raw_mesh, gmsh_to_local, tags.outlet_physical_group, outlet_boundary_id);

    return raw_mesh;
}

} // namespace cfd