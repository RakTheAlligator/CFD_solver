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

constexpr int gmsh_line_2_element_type = 1;
constexpr int gmsh_triangle_3_element_type = 2;
constexpr int gmsh_quadrilateral_4_element_type = 3;

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

struct RectangleModelTags
{
    int surface_tag{};
    int inlet_group_tag{};
    int wall_group_tag{};
    int outlet_group_tag{};
};

struct CellExtractionSpec
{
    int gmsh_element_type{};
    Index node_count{};
    CellType cell_type{};
};

[[nodiscard]]
CellExtractionSpec cell_extraction_spec(const CellType cell_type)
{
    switch (cell_type)
    {
    case CellType::Triangle:
        return {
            .gmsh_element_type = gmsh_triangle_3_element_type,
            .node_count = 3,
            .cell_type = CellType::Triangle,
        };

    case CellType::Quadrilateral:
        return {
            .gmsh_element_type = gmsh_quadrilateral_4_element_type,
            .node_count = 4,
            .cell_type = CellType::Quadrilateral,
        };
    }

    throw std::invalid_argument("Unsupported cell type.");
}

void validate_mesh_generation_inputs(const RectangleGeometry &geometry, const MeshGenerationOptions &options)
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

    static_cast<void>(cell_extraction_spec(options.cell_type));
}

RectangleModelTags create_rectangle(const RectangleGeometry &geometry, const double mesh_size)
{
    const int bottom_left_point_tag{gmsh::model::geo::addPoint(0.0, 0.0, 0.0, mesh_size)};

    const int bottom_right_point_tag{gmsh::model::geo::addPoint(geometry.length, 0.0, 0.0, mesh_size)};

    const int top_right_point_tag{gmsh::model::geo::addPoint(geometry.length, geometry.height, 0.0, mesh_size)};

    const int top_left_point_tag{gmsh::model::geo::addPoint(0.0, geometry.height, 0.0, mesh_size)};

    const int bottom_curve_tag{gmsh::model::geo::addLine(bottom_left_point_tag, bottom_right_point_tag)};

    const int right_curve_tag{gmsh::model::geo::addLine(bottom_right_point_tag, top_right_point_tag)};

    const int top_curve_tag{gmsh::model::geo::addLine(top_right_point_tag, top_left_point_tag)};

    const int left_curve_tag{gmsh::model::geo::addLine(top_left_point_tag, bottom_left_point_tag)};

    const int curve_loop_tag{gmsh::model::geo::addCurveLoop({
        bottom_curve_tag,
        right_curve_tag,
        top_curve_tag,
        left_curve_tag,
    })};

    const int surface_tag{gmsh::model::geo::addPlaneSurface({curve_loop_tag})};

    // Make the entities created in the geometry kernel
    // available to the rest of the Gmsh model.
    gmsh::model::geo::synchronize();

    const int inlet_group_tag{gmsh::model::addPhysicalGroup(1, {left_curve_tag})};

    gmsh::model::setPhysicalName(1, inlet_group_tag, "inlet");

    const int wall_group_tag{gmsh::model::addPhysicalGroup(1, {
                                                                  bottom_curve_tag,
                                                                  top_curve_tag,
                                                              })};

    gmsh::model::setPhysicalName(1, wall_group_tag, "wall");

    const int outlet_group_tag{gmsh::model::addPhysicalGroup(1, {right_curve_tag})};

    gmsh::model::setPhysicalName(1, outlet_group_tag, "outlet");

    const int fluid_group_tag{gmsh::model::addPhysicalGroup(2, {surface_tag})};

    gmsh::model::setPhysicalName(2, fluid_group_tag, "fluid");

    return {
        surface_tag,
        inlet_group_tag,
        wall_group_tag,
        outlet_group_tag,
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
Index node_id_from_gmsh_tag(const std::unordered_map<std::size_t, Index> &node_id_by_gmsh_tag,
                            const std::size_t gmsh_node_tag)
{
    const auto iterator{node_id_by_gmsh_tag.find(gmsh_node_tag)};

    if (iterator == node_id_by_gmsh_tag.end())
    {
        throw std::runtime_error("Unknown Gmsh node tag.");
    }

    return iterator->second;
}

[[nodiscard]]
std::unordered_map<std::size_t, Index> extract_nodes(RawMeshData &raw_mesh)
{
    std::vector<std::size_t> gmsh_node_tags;
    std::vector<double> coordinates;
    std::vector<double> parametric_coordinates;

    gmsh::model::mesh::getNodes(gmsh_node_tags, coordinates, parametric_coordinates, -1, -1, false, false);

    if (coordinates.size() != 3 * gmsh_node_tags.size())
    {
        throw std::runtime_error("Invalid node coordinate data returned by Gmsh.");
    }

    raw_mesh.nodes.reserve(gmsh_node_tags.size());

    std::unordered_map<std::size_t, Index> node_id_by_gmsh_tag;

    node_id_by_gmsh_tag.reserve(gmsh_node_tags.size());

    for (Index node_id = 0; node_id < gmsh_node_tags.size(); ++node_id)
    {
        const std::size_t gmsh_node_tag{gmsh_node_tags[node_id]};

        node_id_by_gmsh_tag.emplace(gmsh_node_tag, node_id);

        raw_mesh.nodes.push_back({
            coordinates[3 * node_id],
            coordinates[3 * node_id + 1],
        });
    }

    return node_id_by_gmsh_tag;
}

void extract_cells(RawMeshData &raw_mesh, const std::unordered_map<std::size_t, Index> &node_id_by_gmsh_tag,
                   const int surface_tag, const CellType requested_cell_type)
{
    std::vector<int> gmsh_element_types;
    std::vector<std::vector<std::size_t>> gmsh_element_tags;
    std::vector<std::vector<std::size_t>> gmsh_element_node_tags;

    gmsh::model::mesh::getElements(gmsh_element_types, gmsh_element_tags, gmsh_element_node_tags, 2, surface_tag);

    const CellExtractionSpec extraction_spec{cell_extraction_spec(requested_cell_type)};

    Index total_cell_count{};
    Index total_cell_node_count{};

    for (std::size_t element_type_index = 0; element_type_index < gmsh_element_types.size(); ++element_type_index)
    {
        if (gmsh_element_types[element_type_index] != extraction_spec.gmsh_element_type)
        {
            throw std::runtime_error("Gmsh generated a 2D element type that does not "
                                     "match the requested cell type.");
        }

        const auto &gmsh_node_tags{gmsh_element_node_tags[element_type_index]};

        if (gmsh_node_tags.size() % extraction_spec.node_count != 0)
        {
            throw std::runtime_error("Invalid cell connectivity returned by Gmsh.");
        }

        total_cell_count += gmsh_node_tags.size() / extraction_spec.node_count;

        total_cell_node_count += gmsh_node_tags.size();
    }

    raw_mesh.cell_types.reserve(total_cell_count);

    raw_mesh.cell_nodes.reserve(total_cell_node_count);

    raw_mesh.cell_node_offsets.reserve(total_cell_count + 1);

    raw_mesh.cell_node_offsets.push_back(0);

    for (std::size_t element_type_index = 0; element_type_index < gmsh_element_types.size(); ++element_type_index)
    {
        const auto &gmsh_node_tags{gmsh_element_node_tags[element_type_index]};

        for (std::size_t first_node_position = 0; first_node_position < gmsh_node_tags.size();
             first_node_position += extraction_spec.node_count)
        {
            raw_mesh.cell_types.push_back(extraction_spec.cell_type);

            for (Index local_node_index = 0; local_node_index < extraction_spec.node_count; ++local_node_index)
            {
                raw_mesh.cell_nodes.push_back(
                    node_id_from_gmsh_tag(node_id_by_gmsh_tag, gmsh_node_tags[first_node_position + local_node_index]));
            }

            raw_mesh.cell_node_offsets.push_back(raw_mesh.cell_nodes.size());
        }
    }
}

void extract_boundary_edges(RawMeshData &raw_mesh, const std::unordered_map<std::size_t, Index> &node_id_by_gmsh_tag,
                            const int boundary_group_tag, const BoundaryId boundary_id)
{
    std::vector<int> curve_tags;

    gmsh::model::getEntitiesForPhysicalGroup(1, boundary_group_tag, curve_tags);

    for (const int curve_tag : curve_tags)
    {
        std::vector<int> gmsh_element_types;
        std::vector<std::vector<std::size_t>> gmsh_element_tags;
        std::vector<std::vector<std::size_t>> gmsh_element_node_tags;

        gmsh::model::mesh::getElements(gmsh_element_types, gmsh_element_tags, gmsh_element_node_tags, 1, curve_tag);

        for (std::size_t element_type_index = 0; element_type_index < gmsh_element_types.size(); ++element_type_index)
        {
            if (gmsh_element_types[element_type_index] != gmsh_line_2_element_type)
            {
                throw std::runtime_error("Gmsh generated a non-linear boundary element.");
            }

            const auto &gmsh_node_tags{gmsh_element_node_tags[element_type_index]};

            if (gmsh_node_tags.size() % 2 != 0)
            {
                throw std::runtime_error("Invalid boundary connectivity returned by Gmsh.");
            }

            for (std::size_t connectivity_position = 0; connectivity_position < gmsh_node_tags.size();
                 connectivity_position += 2)
            {
                raw_mesh.boundary_edges.push_back({
                    {
                        node_id_from_gmsh_tag(node_id_by_gmsh_tag, gmsh_node_tags[connectivity_position]),
                        node_id_from_gmsh_tag(node_id_by_gmsh_tag, gmsh_node_tags[connectivity_position + 1]),
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
    validate_mesh_generation_inputs(geometry, options);

    GmshSession gmsh_session;

    gmsh::model::add("rectangle");

    const RectangleModelTags tags{create_rectangle(geometry, options.mesh_size)};

    configure_surface_mesh(options.cell_type, tags.surface_tag);

    gmsh::model::mesh::generate(2);

    RawMeshData raw_mesh;

    raw_mesh.boundary_groups = {
        {inlet_boundary_id, "inlet"},
        {wall_boundary_id, "wall"},
        {outlet_boundary_id, "outlet"},
    };

    const auto node_id_by_gmsh_tag{extract_nodes(raw_mesh)};

    extract_cells(raw_mesh, node_id_by_gmsh_tag, tags.surface_tag, options.cell_type);

    extract_boundary_edges(raw_mesh, node_id_by_gmsh_tag, tags.inlet_group_tag, inlet_boundary_id);

    extract_boundary_edges(raw_mesh, node_id_by_gmsh_tag, tags.wall_group_tag, wall_boundary_id);

    extract_boundary_edges(raw_mesh, node_id_by_gmsh_tag, tags.outlet_group_tag, outlet_boundary_id);

    return raw_mesh;
}

} // namespace cfd