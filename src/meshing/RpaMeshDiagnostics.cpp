#include "cfd/meshing/RpaMeshDiagnostics.hpp"
#include "RpaMeshDiscretization.hpp"

#include "cfd/mesh/Mesh.hpp"
#include "cfd/meshing/RpaRegionLayout.hpp"

#include <algorithm>
#include <array>
#include <iomanip>
#include <limits>
#include <ostream>
#include <sstream>
#include <vector>

namespace cfd
{
void write_rpa_mesh_diagnostics(std::ostream &output, const Mesh &mesh, const RpaNozzleDomain &domain,
                                const MeshGenerationOptions &generation, const RpaNozzleMeshingOptions &options,
                                const std::span<double> region_ids, const std::span<double> block_ids)
{
    if ((!region_ids.empty() && region_ids.size() != mesh.cell_count()) ||
        (!block_ids.empty() && block_ids.size() != mesh.cell_count()))
        throw std::invalid_argument("RPA region output cardinality mismatch.");
    const RpaRegionLayout layout{domain, options.regions};
    std::vector<double> plume_interfaces;
    if (options.regions.plume_length_factor)
    {
        const auto grids{detail::discretize_rpa(domain, layout, generation, options)};
        output << "[RPA mapped blocks: columns / rows / cells; all logical coordinates verified during generation]\n";
        const std::array<const char *, 5> block_names{"engine", "band", "plume", "farfield_upper",
                                                      "farfield_downstream"};
        Index block{};
        for (const auto *grid : {&grids.engine, &grids.band, &grids.plume, &grids.farfield, &grids.farfield_downstream})
            output << "  " << block_names.at(block++) << " : " << grid->x.size() - 1 << " / " << grid->radial_cells
                   << " / " << (grid->x.size() - 1) * grid->radial_cells << '\n';
        for (Index i = 1; i < grids.plume_block_columns.size(); ++i)
        {
            const Index begin{grids.plume_block_columns[i - 1]}, end{grids.plume_block_columns[i]};
            plume_interfaces.push_back(grids.plume.x[end]);
            output << "  plume_block " << i << " : x=[" << grids.plume.x[begin] << ',' << grids.plume.x[end]
                   << "] columns=" << end - begin << " rows=" << grids.plume.radial_cells
                   << " cells=" << (end - begin) * grids.plume.radial_cells << '\n';
        }
        output << "transition_outer_radius_upstream=" << layout.transition_outer_radius(domain.outer_wall.front().x)
               << "\ntransition_outer_radius_end=" << layout.transition_outer_radius(layout.plume_end_x())
               << "\nfarfield_begin_x=" << layout.farfield_begin_x() << '\n';
    }
    constexpr std::array<const char *, 5> names{"engine", "band", "plume", "transition", "farfield"};
    std::array<Index, 5> counts{}, triangles{}, quads{};
    std::array<double, 5> minimum{};
    minimum.fill(std::numeric_limits<double>::infinity());
    auto angular_minimum{minimum};
    std::array<double, 5> angular_sum{};
    for (Index cell = 0; cell < mesh.cell_count(); ++cell)
    {
        const Index region{static_cast<Index>(layout.region(mesh.cell_centers()[cell])) - 1};
        ++counts.at(region);
        if (mesh.cell_types()[cell] == CellType::Triangle)
            ++triangles.at(region);
        else
            ++quads.at(region);
        minimum.at(region) = std::min(minimum.at(region), mesh.cell_qualities()[cell]);
        std::array<Point2, 4> vertices{};
        const Index begin{mesh.cell_node_offsets()[cell]}, count{mesh.cell_node_offsets()[cell + 1] - begin};
        for (Index i = 0; i < count; ++i)
            vertices.at(i) = mesh.nodes()[mesh.cell_nodes()[begin + i]];
        const double angular{detail::rpa_angular_quality(std::span{vertices}.first(count))};
        angular_minimum.at(region) = std::min(angular_minimum.at(region), angular);
        angular_sum.at(region) += angular;
        if (!region_ids.empty())
            region_ids[cell] = static_cast<double>(region + 1);
        if (!block_ids.empty())
        {
            double id{static_cast<double>(region + 1)};
            if (options.regions.plume_length_factor)
            {
                if (region == 2)
                    id = 3 + static_cast<double>(std::upper_bound(plume_interfaces.begin(), plume_interfaces.end(),
                                                                  mesh.cell_centers()[cell].x) -
                                                 plume_interfaces.begin());
                else if (region == 3)
                    id = 8;
                else if (region == 4)
                    id = mesh.cell_centers()[cell].x < layout.farfield_begin_x() ? 6 : 7;
            }
            block_ids[cell] = id;
        }
    }
    std::array<double, 2> wall_error{};
    for (Index face = 0; face < mesh.face_count(); ++face)
    {
        const auto boundary{mesh.face_boundary_ids()[face]};
        if (boundary == 2 || boundary == 3)
        {
            const auto &nodes{mesh.faces()[face].node_ids};
            const Index region{boundary - 2};
            wall_error.at(region) =
                std::max(wall_error.at(region),
                         detail::rpa_wall_chord_error(boundary == 2 ? domain.inner_wall : domain.outer_wall,
                                                      mesh.nodes()[nodes[0]], mesh.nodes()[nodes[1]]));
        }
    }
    std::ostringstream report;
    report << "\n[RPA baseCellSize: " << generation.mesh_size << " m; angular threshold: " << options.min_cell_quality
           << "]\n[RPA logical regions: cells / TRI / QUAD / min Mesh shape quality / angular min / angular mean]\n";
    for (Index region = 0; region < names.size(); ++region)
        report << "  " << region + 1 << ' ' << names.at(region) << " : " << counts.at(region) << " / "
               << triangles.at(region) << " / " << quads.at(region) << " / "
               << (counts.at(region) ? minimum.at(region) : 0.0) << " / "
               << (counts.at(region) ? angular_minimum.at(region) : 0.0) << " / "
               << (counts.at(region) ? angular_sum.at(region) / static_cast<double>(counts.at(region)) : 0.0) << '\n';
    report << "[RPA wall-face approximation: maximum vertical chord error]\n" << std::setprecision(10);
    for (Index region = 0; region < wall_error.size(); ++region)
        report << "  " << names.at(region) << " : " << wall_error.at(region) << " m; tolerance "
               << domain.geometric_tolerance << " m; "
               << (wall_error.at(region) <= domain.geometric_tolerance ? "WITHIN TOLERANCE" : "EXCEEDED") << '\n';
    report << "  CAD import tolerance is unchanged. ";
    if (options.preserve_wall_vertices)
        report << "All flattened wall vertices are preserved; this does not claim exact curved-CAD walls.\n";
    else
        report << "Topological validity does not certify wall-face geometric accuracy; mean-height spacing is "
                  "not refined to enforce this tolerance.\n";
    output << report.str();
}
} // namespace cfd
