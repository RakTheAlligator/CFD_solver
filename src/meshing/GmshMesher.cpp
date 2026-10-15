#include "cfd/meshing/GmshMesher.hpp"
#include "RpaMeshDiscretization.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <gmsh.h>
#include <iterator>
#include <limits>
#include <map>
#include <numeric>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

namespace cfd
{

namespace
{

// Gmsh element type codes used by the API for the supported first-order
// one- and two-dimensional elements.
constexpr int gmsh_line_2_element_type = 1;
constexpr int gmsh_triangle_3_element_type = 2;
constexpr int gmsh_quadrilateral_4_element_type = 3;

// Solver boundary IDs are deliberately independent of the physical-group tags
// assigned by Gmsh.
constexpr BoundaryId inlet_boundary_id{0};
constexpr BoundaryId wall_boundary_id{1};
constexpr BoundaryId outlet_boundary_id{2};

class GmshSession
{
  public:
    GmshSession()
    {
        gmsh::initialize();
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

struct PhysicalBoundaryTags
{
    BoundaryGroup group;
    int physical_tag{};
};

struct ModelTags
{
    std::vector<int> surface_tags;
    std::vector<PhysicalBoundaryTags> boundaries;
    bool transfinite_quadrilateral_mesh{};
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

void validate_mesh_generation_options(const MeshGenerationOptions &options)
{
    if (!std::isfinite(options.mesh_size) || !(options.mesh_size > 0.0))
    {
        throw std::invalid_argument("Mesh size must be finite and positive.");
    }

    // Validate the enum value through the same mapping used during extraction.
    static_cast<void>(cell_extraction_spec(options.cell_type));
    if (options.hybrid && options.cell_type != CellType::Quadrilateral)
        throw std::invalid_argument("HYBRID requires cell_type == Quadrilateral.");
}

void validate_automatic_meshing_options(const AutomaticMeshingOptions &options)
{
    if (!std::isfinite(options.maximum_growth_rate) || !(options.maximum_growth_rate > 1.0))
    {
        throw std::invalid_argument("Maximum mesh growth rate must be finite and greater than one.");
    }
    if (!std::isfinite(options.wall_refinement_factor) || !(options.wall_refinement_factor > 0.0) ||
        options.wall_refinement_factor > 1.0)
    {
        throw std::invalid_argument("Wall refinement factor must be finite and in (0, 1].");
    }
    if (options.wall_refinement_layers == 0 ||
        options.wall_refinement_layers >= static_cast<Index>(std::numeric_limits<int>::max()))
    {
        throw std::invalid_argument("Wall refinement layer count must be positive and representable by Gmsh.");
    }
}

void validate_geometry(const RectangleGeometry &geometry)
{
    if (!std::isfinite(geometry.length) || !(geometry.length > 0.0))
    {
        throw std::invalid_argument("Rectangle length must be finite and positive.");
    }

    if (!std::isfinite(geometry.height) || !(geometry.height > 0.0))
    {
        throw std::invalid_argument("Rectangle height must be finite and positive.");
    }
}

void validate_geometry(const BackwardFacingStepGeometry &geometry)
{
    if (!std::isfinite(geometry.upstream_length) || !(geometry.upstream_length > 0.0))
    {
        throw std::invalid_argument("Backward-facing-step upstream length must be finite and positive.");
    }
    if (!std::isfinite(geometry.downstream_length) || !(geometry.downstream_length > 0.0))
    {
        throw std::invalid_argument("Backward-facing-step downstream length must be finite and positive.");
    }
    if (!std::isfinite(geometry.channel_height) || !(geometry.channel_height > 0.0))
    {
        throw std::invalid_argument("Backward-facing-step channel height must be finite and positive.");
    }
    if (!std::isfinite(geometry.step_height) || !(geometry.step_height > 0.0) ||
        !(geometry.step_height < geometry.channel_height))
    {
        throw std::invalid_argument(
            "Backward-facing-step step height must be finite, positive, and less than the channel height.");
    }
}

void validate_backward_facing_step_meshing_options(const BackwardFacingStepMeshingOptions &options)
{
    if (!std::isfinite(options.step_refinement_factor) || !(options.step_refinement_factor > 0.0) ||
        options.step_refinement_factor > 1.0)
    {
        throw std::invalid_argument("Step refinement factor must be finite and in (0, 1].");
    }
    if (options.step_refinement_layers == 0 ||
        options.step_refinement_layers >= static_cast<Index>(std::numeric_limits<int>::max()))
    {
        throw std::invalid_argument("Step refinement layer count must be positive and representable by Gmsh.");
    }
}

[[nodiscard]]
ModelTags add_physical_groups(const std::vector<int> &surface_tags, const std::vector<int> &inlet_curve_tags,
                              const std::vector<int> &wall_curve_tags, const std::vector<int> &outlet_curve_tags)
{
    // Entities created in the geometry kernel must be synchronized before they
    // can be referenced by the model API, including physical groups.
    gmsh::model::geo::synchronize();

    const int inlet_group_tag{gmsh::model::addPhysicalGroup(1, inlet_curve_tags)};
    gmsh::model::setPhysicalName(1, inlet_group_tag, "inlet");

    const int wall_group_tag{gmsh::model::addPhysicalGroup(1, wall_curve_tags)};
    gmsh::model::setPhysicalName(1, wall_group_tag, "wall");

    const int outlet_group_tag{gmsh::model::addPhysicalGroup(1, outlet_curve_tags)};
    gmsh::model::setPhysicalName(1, outlet_group_tag, "outlet");

    const int fluid_group_tag{gmsh::model::addPhysicalGroup(2, surface_tags)};
    gmsh::model::setPhysicalName(2, fluid_group_tag, "fluid");

    return {surface_tags,
            {{{inlet_boundary_id, "inlet"}, inlet_group_tag},
             {{wall_boundary_id, "wall"}, wall_group_tag},
             {{outlet_boundary_id, "outlet"}, outlet_group_tag}},
            false};
}

ModelTags create_rectangle(const RectangleGeometry &geometry, const double mesh_size)
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

    return add_physical_groups({surface_tag}, {left_curve_tag}, {bottom_curve_tag, top_curve_tag}, {right_curve_tag});
}

ModelTags create_backward_facing_step(const BackwardFacingStepGeometry &geometry, const double mesh_size)
{
    const double outlet_x{geometry.upstream_length + geometry.downstream_length};
    const int inlet_bottom_point_tag{gmsh::model::geo::addPoint(0.0, geometry.step_height, 0.0, mesh_size)};
    const int step_corner_point_tag{
        gmsh::model::geo::addPoint(geometry.upstream_length, geometry.step_height, 0.0, mesh_size)};
    const int step_bottom_point_tag{gmsh::model::geo::addPoint(geometry.upstream_length, 0.0, 0.0, mesh_size)};
    const int outlet_bottom_point_tag{gmsh::model::geo::addPoint(outlet_x, 0.0, 0.0, mesh_size)};
    const int outlet_top_point_tag{gmsh::model::geo::addPoint(outlet_x, geometry.channel_height, 0.0, mesh_size)};
    const int inlet_top_point_tag{gmsh::model::geo::addPoint(0.0, geometry.channel_height, 0.0, mesh_size)};

    const int upstream_bottom_curve_tag{gmsh::model::geo::addLine(inlet_bottom_point_tag, step_corner_point_tag)};
    const int step_curve_tag{gmsh::model::geo::addLine(step_corner_point_tag, step_bottom_point_tag)};
    const int downstream_bottom_curve_tag{gmsh::model::geo::addLine(step_bottom_point_tag, outlet_bottom_point_tag)};
    const int outlet_curve_tag{gmsh::model::geo::addLine(outlet_bottom_point_tag, outlet_top_point_tag)};
    const int top_curve_tag{gmsh::model::geo::addLine(outlet_top_point_tag, inlet_top_point_tag)};
    const int inlet_curve_tag{gmsh::model::geo::addLine(inlet_top_point_tag, inlet_bottom_point_tag)};

    const int curve_loop_tag{
        gmsh::model::geo::addCurveLoop({upstream_bottom_curve_tag, step_curve_tag, downstream_bottom_curve_tag,
                                        outlet_curve_tag, top_curve_tag, inlet_curve_tag})};
    const int surface_tag{gmsh::model::geo::addPlaneSurface({curve_loop_tag})};

    return add_physical_groups({surface_tag}, {inlet_curve_tag},
                               {upstream_bottom_curve_tag, step_curve_tag, downstream_bottom_curve_tag, top_curve_tag},
                               {outlet_curve_tag});
}

struct AxisSegment
{
    double begin{};
    double end{};
    int cell_count{};
    double progression{1.0};
};

[[nodiscard]]
int nominal_cell_count(const double length, const double mesh_size)
{
    const double count{std::ceil(length / mesh_size)};
    if (!std::isfinite(count) || count >= static_cast<double>(std::numeric_limits<int>::max()))
    {
        throw std::invalid_argument("Requested mesh resolution exceeds Gmsh's supported transfinite node count.");
    }
    return std::max(1, static_cast<int>(count));
}

struct TransitionSpecification
{
    int cell_count{};
    double growth_rate{};
    double nominal_width{};
};

[[nodiscard]]
TransitionSpecification make_transition_specification(const double mesh_size, const double refinement_factor,
                                                      const Index nominal_layers, const double maximum_growth_rate)
{
    const double minimum_layer_count{std::ceil(std::log(1.0 / refinement_factor) / std::log(maximum_growth_rate))};
    if (!std::isfinite(minimum_layer_count) ||
        minimum_layer_count >= static_cast<double>(std::numeric_limits<int>::max()))
    {
        throw std::invalid_argument("Requested mesh transition exceeds Gmsh's supported transfinite node count.");
    }

    const int cell_count{std::max(static_cast<int>(nominal_layers), static_cast<int>(minimum_layer_count))};
    const double growth_rate{std::pow(1.0 / refinement_factor, 1.0 / static_cast<double>(cell_count))};
    const double nominal_width{growth_rate == 1.0
                                   ? static_cast<double>(cell_count) * refinement_factor * mesh_size
                                   : refinement_factor * mesh_size *
                                         std::expm1(static_cast<double>(cell_count) * std::log(growth_rate)) /
                                         std::expm1(std::log(growth_rate))};
    return {cell_count, growth_rate, nominal_width};
}

void append_uniform_segment(std::vector<AxisSegment> &segments, const double begin, const double end,
                            const double mesh_size)
{
    if (end > begin)
    {
        segments.push_back({begin, end, nominal_cell_count(end - begin, mesh_size), 1.0});
    }
}

[[nodiscard]]
double geometric_series_sum(double first, double ratio, int term_count);

[[nodiscard]]
double relative_transition_width(const TransitionSpecification &transition)
{
    // N refined cells are followed by a bulk cell. Since r^N = 1/factor,
    // the first refined cell is factor times the compatible bulk size and the
    // final transition-to-bulk ratio is r.
    const double first_relative_size{std::pow(transition.growth_rate, -static_cast<double>(transition.cell_count))};
    return geometric_series_sum(first_relative_size, transition.growth_rate, transition.cell_count);
}

[[nodiscard]]
int compatible_uniform_cell_count(const double interval_length, const double mesh_size,
                                  const double transition_relative_width)
{
    const double count{std::ceil(interval_length / mesh_size - transition_relative_width)};
    if (!std::isfinite(count) || count >= static_cast<double>(std::numeric_limits<int>::max()))
    {
        throw std::invalid_argument("Requested mesh resolution exceeds Gmsh's supported transfinite node count.");
    }
    return std::max(0, static_cast<int>(count));
}

void append_one_sided_refined_interval(std::vector<AxisSegment> &segments, const double begin, const double end,
                                       const bool refine_at_begin, const double mesh_size,
                                       const double refinement_factor, const Index nominal_layers,
                                       const double maximum_growth_rate)
{
    if (refinement_factor == 1.0)
    {
        append_uniform_segment(segments, begin, end, mesh_size);
        return;
    }

    const TransitionSpecification transition{
        make_transition_specification(mesh_size, refinement_factor, nominal_layers, maximum_growth_rate)};
    const double interval_length{end - begin};
    const double transition_weight{relative_transition_width(transition)};
    const int uniform_cell_count{compatible_uniform_cell_count(interval_length, mesh_size, transition_weight)};
    const double bulk_cell_size{interval_length / (transition_weight + static_cast<double>(uniform_cell_count))};
    const double width{bulk_cell_size * transition_weight};
    if (refine_at_begin)
    {
        segments.push_back({begin, begin + width, transition.cell_count, transition.growth_rate});
        if (uniform_cell_count > 0)
        {
            segments.push_back({begin + width, end, uniform_cell_count, 1.0});
        }
    }
    else
    {
        if (uniform_cell_count > 0)
        {
            segments.push_back({begin, end - width, uniform_cell_count, 1.0});
        }
        segments.push_back({end - width, end, transition.cell_count, 1.0 / transition.growth_rate});
    }
}

struct StreamwiseLayout
{
    std::vector<AxisSegment> segments;
    double downstream_wall_relaxation_end{};
};

[[nodiscard]]
double geometric_series_sum(const double first, const double ratio, const int term_count)
{
    if (ratio == 1.0)
    {
        return first * static_cast<double>(term_count);
    }
    return first * std::expm1(static_cast<double>(term_count) * std::log(ratio)) / std::expm1(std::log(ratio));
}

[[nodiscard]]
StreamwiseLayout make_streamwise_layout(const BackwardFacingStepGeometry &geometry,
                                        const MeshGenerationOptions &options,
                                        const AutomaticMeshingOptions &automatic_options,
                                        const BackwardFacingStepMeshingOptions &step_options)
{
    StreamwiseLayout layout;
    const double refinement_factor{
        std::min(automatic_options.wall_refinement_factor, step_options.step_refinement_factor)};
    append_one_sided_refined_interval(layout.segments, 0.0, geometry.upstream_length, false, options.mesh_size,
                                      refinement_factor, step_options.step_refinement_layers,
                                      automatic_options.maximum_growth_rate);

    const double outlet_x{geometry.upstream_length + geometry.downstream_length};
    const AxisSegment &upstream_step_segment{layout.segments.back()};
    const double upstream_first_cell_size{
        (upstream_step_segment.end - upstream_step_segment.begin) /
        geometric_series_sum(1.0, upstream_step_segment.progression, upstream_step_segment.cell_count)};
    const double step_adjacent_cell_size{
        upstream_first_cell_size *
        std::pow(upstream_step_segment.progression, static_cast<double>(upstream_step_segment.cell_count - 1))};
    const double downstream_refinement_factor{std::min(refinement_factor, step_adjacent_cell_size / options.mesh_size)};
    std::vector<AxisSegment> downstream_segments;
    append_one_sided_refined_interval(downstream_segments, geometry.upstream_length, outlet_x, true, options.mesh_size,
                                      downstream_refinement_factor, step_options.step_refinement_layers,
                                      automatic_options.maximum_growth_rate);

    if (refinement_factor == 1.0)
    {
        layout.segments.insert(layout.segments.end(), downstream_segments.begin(), downstream_segments.end());
        layout.downstream_wall_relaxation_end = geometry.upstream_length;
        return layout;
    }

    if (automatic_options.wall_refinement_factor == 1.0)
    {
        layout.segments.insert(layout.segments.end(), downstream_segments.begin(), downstream_segments.end());
        layout.downstream_wall_relaxation_end = geometry.upstream_length;
        return layout;
    }

    const TransitionSpecification wall_transition{
        make_transition_specification(options.mesh_size, automatic_options.wall_refinement_factor,
                                      automatic_options.wall_refinement_layers, automatic_options.maximum_growth_rate)};
    const AxisSegment &downstream_transition{downstream_segments.front()};
    const double transition_width{downstream_transition.end - downstream_transition.begin};
    const double first_cell_size{transition_width / geometric_series_sum(1.0, downstream_transition.progression,
                                                                         downstream_transition.cell_count)};
    const double requested_relaxation_width{std::min(wall_transition.nominal_width, transition_width)};
    // Unit progression uses the uniform-spacing limit of the geometric series.
    const double layer_count_value{std::ceil(
        downstream_transition.progression == 1.0
            ? requested_relaxation_width / first_cell_size
            : std::log1p(requested_relaxation_width * (downstream_transition.progression - 1.0) / first_cell_size) /
                  std::log(downstream_transition.progression))};
    if (!std::isfinite(layer_count_value))
    {
        throw std::invalid_argument("Requested mesh wall relaxation cell count must be finite.");
    }
    const int relaxation_cell_count{
        static_cast<int>(std::clamp(layer_count_value, 1.0, static_cast<double>(downstream_transition.cell_count)))};
    const double relaxation_width{
        geometric_series_sum(first_cell_size, downstream_transition.progression, relaxation_cell_count)};
    layout.downstream_wall_relaxation_end = geometry.upstream_length + relaxation_width;
    layout.segments.push_back({geometry.upstream_length, layout.downstream_wall_relaxation_end, relaxation_cell_count,
                               downstream_transition.progression});

    const int remaining_transition_cells{downstream_transition.cell_count - relaxation_cell_count};
    if (remaining_transition_cells > 0)
    {
        layout.segments.push_back({layout.downstream_wall_relaxation_end, downstream_transition.end,
                                   remaining_transition_cells, downstream_transition.progression});
    }
    layout.segments.insert(layout.segments.end(), std::next(downstream_segments.begin()), downstream_segments.end());
    return layout;
}

void append_two_sided_refined_interval(std::vector<AxisSegment> &segments, const double begin, const double end,
                                       const double mesh_size, const double refinement_factor,
                                       const Index nominal_layers, const double maximum_growth_rate)
{
    if (refinement_factor == 1.0)
    {
        append_uniform_segment(segments, begin, end, mesh_size);
        return;
    }

    const TransitionSpecification transition{
        make_transition_specification(mesh_size, refinement_factor, nominal_layers, maximum_growth_rate)};
    const double interval_length{end - begin};
    const double transition_weight{relative_transition_width(transition)};
    const int uniform_cell_count{compatible_uniform_cell_count(interval_length, mesh_size, 2.0 * transition_weight)};
    const double bulk_cell_size{interval_length / (2.0 * transition_weight + static_cast<double>(uniform_cell_count))};
    const double width{bulk_cell_size * transition_weight};
    segments.push_back({begin, begin + width, transition.cell_count, transition.growth_rate});
    if (uniform_cell_count > 0)
    {
        segments.push_back({begin + width, end - width, uniform_cell_count, 1.0});
    }
    segments.push_back({end - width, end, transition.cell_count, 1.0 / transition.growth_rate});
}

[[nodiscard]]
std::vector<AxisSegment> make_wall_normal_segments(const BackwardFacingStepGeometry &geometry,
                                                   const MeshGenerationOptions &options,
                                                   const AutomaticMeshingOptions &automatic_options)
{
    std::vector<AxisSegment> segments;
    append_one_sided_refined_interval(segments, 0.0, geometry.step_height, true, options.mesh_size,
                                      automatic_options.wall_refinement_factor,
                                      automatic_options.wall_refinement_layers, automatic_options.maximum_growth_rate);
    append_two_sided_refined_interval(segments, geometry.step_height, geometry.channel_height, options.mesh_size,
                                      automatic_options.wall_refinement_factor,
                                      automatic_options.wall_refinement_layers, automatic_options.maximum_growth_rate);
    return segments;
}

[[nodiscard]]
std::vector<AxisSegment> make_relaxed_wall_normal_segments(const BackwardFacingStepGeometry &geometry,
                                                           const AutomaticMeshingOptions &automatic_options,
                                                           const std::vector<AxisSegment> &near_step_segments)
{
    std::vector<AxisSegment> relaxed_segments{near_step_segments};
    Index first_upper_segment{};
    while (first_upper_segment < near_step_segments.size() &&
           near_step_segments[first_upper_segment].end <= geometry.step_height)
    {
        ++first_upper_segment;
    }

    int upper_cell_count{};
    for (Index segment_index = first_upper_segment; segment_index < near_step_segments.size(); ++segment_index)
    {
        upper_cell_count += near_step_segments[segment_index].cell_count;
    }
    const double progression{
        automatic_options.wall_refinement_factor == 1.0
            ? 1.0
            : std::pow(automatic_options.wall_refinement_factor, 1.0 / static_cast<double>(upper_cell_count))};
    const double total_weight{geometric_series_sum(1.0, progression, upper_cell_count)};
    const double scale{(geometry.channel_height - geometry.step_height) / total_weight};

    double segment_begin{geometry.step_height};
    int completed_cells{};
    for (Index segment_index = first_upper_segment; segment_index < relaxed_segments.size(); ++segment_index)
    {
        AxisSegment &segment{relaxed_segments[segment_index]};
        const double first_weight{std::pow(progression, static_cast<double>(completed_cells))};
        const double segment_width{scale * geometric_series_sum(first_weight, progression, segment.cell_count)};
        segment.begin = segment_begin;
        segment.end =
            segment_index + 1 == relaxed_segments.size() ? geometry.channel_height : segment_begin + segment_width;
        segment.progression = progression;
        segment_begin = segment.end;
        completed_cells += segment.cell_count;
    }
    return relaxed_segments;
}

[[nodiscard]]
std::vector<double> axis_coordinates(const std::vector<AxisSegment> &segments)
{
    std::vector<double> coordinates;
    coordinates.reserve(segments.size() + 1);
    coordinates.push_back(segments.front().begin);
    for (const AxisSegment &segment : segments)
    {
        coordinates.push_back(segment.end);
    }
    return coordinates;
}

ModelTags create_structured_backward_facing_step(const BackwardFacingStepGeometry &geometry,
                                                 const MeshGenerationOptions &options,
                                                 const AutomaticMeshingOptions &automatic_options,
                                                 const BackwardFacingStepMeshingOptions &step_options)
{
    const StreamwiseLayout streamwise_layout{
        make_streamwise_layout(geometry, options, automatic_options, step_options)};
    const std::vector<AxisSegment> &x_segments{streamwise_layout.segments};
    const std::vector<AxisSegment> near_step_y_segments{
        make_wall_normal_segments(geometry, options, automatic_options)};
    const std::vector<AxisSegment> relaxed_y_segments{
        make_relaxed_wall_normal_segments(geometry, automatic_options, near_step_y_segments)};
    const std::vector<double> x_coordinates{axis_coordinates(x_segments)};
    const Index x_point_count{x_coordinates.size()};
    const Index y_point_count{near_step_y_segments.size() + 1};
    Index first_upper_segment{};
    while (first_upper_segment < near_step_y_segments.size() &&
           near_step_y_segments[first_upper_segment].end <= geometry.step_height)
    {
        ++first_upper_segment;
    }

    std::vector<std::vector<AxisSegment>> vertical_profiles;
    vertical_profiles.reserve(x_point_count);
    for (const double x_coordinate : x_coordinates)
    {
        if (x_coordinate <= geometry.upstream_length ||
            streamwise_layout.downstream_wall_relaxation_end == geometry.upstream_length)
        {
            vertical_profiles.push_back(near_step_y_segments);
            continue;
        }
        if (x_coordinate >= streamwise_layout.downstream_wall_relaxation_end)
        {
            vertical_profiles.push_back(relaxed_y_segments);
            continue;
        }

        const double fraction{(x_coordinate - geometry.upstream_length) /
                              (streamwise_layout.downstream_wall_relaxation_end - geometry.upstream_length)};
        std::vector<AxisSegment> profile{near_step_y_segments};
        for (Index segment_index = 0; segment_index < profile.size(); ++segment_index)
        {
            profile[segment_index].begin =
                std::lerp(near_step_y_segments[segment_index].begin, relaxed_y_segments[segment_index].begin, fraction);
            profile[segment_index].end =
                std::lerp(near_step_y_segments[segment_index].end, relaxed_y_segments[segment_index].end, fraction);
            profile[segment_index].progression =
                std::exp(std::lerp(std::log(near_step_y_segments[segment_index].progression),
                                   std::log(relaxed_y_segments[segment_index].progression), fraction));
        }
        vertical_profiles.push_back(std::move(profile));
    }

    std::vector<int> point_tags(x_point_count * y_point_count);
    std::vector<int> horizontal_curve_tags(x_segments.size() * y_point_count);
    std::vector<int> vertical_curve_tags(x_point_count * near_step_y_segments.size());

    const auto point_tag = [&](const Index x_index, const Index y_index) {
        int &tag{point_tags[x_index * y_point_count + y_index]};
        if (tag == 0)
        {
            const double y_coordinate{y_index == 0 ? vertical_profiles[x_index].front().begin
                                                   : vertical_profiles[x_index][y_index - 1].end};
            tag = gmsh::model::geo::addPoint(x_coordinates[x_index], y_coordinate, 0.0, options.mesh_size);
        }
        return tag;
    };

    const auto horizontal_curve_tag = [&](const Index x_index, const Index y_index) {
        int &tag{horizontal_curve_tags[x_index * y_point_count + y_index]};
        if (tag == 0)
        {
            tag = gmsh::model::geo::addLine(point_tag(x_index, y_index), point_tag(x_index + 1, y_index));
            const AxisSegment &segment{x_segments[x_index]};
            gmsh::model::geo::mesh::setTransfiniteCurve(tag, segment.cell_count + 1, "Progression",
                                                        segment.progression);
        }
        return tag;
    };

    const auto vertical_curve_tag = [&](const Index x_index, const Index y_index) {
        int &tag{vertical_curve_tags[x_index * near_step_y_segments.size() + y_index]};
        if (tag == 0)
        {
            tag = gmsh::model::geo::addLine(point_tag(x_index, y_index), point_tag(x_index, y_index + 1));
            const AxisSegment &segment{vertical_profiles[x_index][y_index]};
            gmsh::model::geo::mesh::setTransfiniteCurve(tag, segment.cell_count + 1, "Progression",
                                                        segment.progression);
        }
        return tag;
    };

    std::vector<int> surface_tags;
    for (Index x_index = 0; x_index < x_segments.size(); ++x_index)
    {
        for (Index y_index = 0; y_index < near_step_y_segments.size(); ++y_index)
        {
            if (x_segments[x_index].end <= geometry.upstream_length && y_index < first_upper_segment)
            {
                continue;
            }

            const int bottom_curve_tag{horizontal_curve_tag(x_index, y_index)};
            const int right_curve_tag{vertical_curve_tag(x_index + 1, y_index)};
            const int top_curve_tag{horizontal_curve_tag(x_index, y_index + 1)};
            const int left_curve_tag{vertical_curve_tag(x_index, y_index)};
            const int curve_loop_tag{
                gmsh::model::geo::addCurveLoop({bottom_curve_tag, right_curve_tag, -top_curve_tag, -left_curve_tag})};
            const int surface_tag{gmsh::model::geo::addPlaneSurface({curve_loop_tag})};
            gmsh::model::geo::mesh::setTransfiniteSurface(surface_tag, "Left",
                                                          {point_tag(x_index, y_index), point_tag(x_index + 1, y_index),
                                                           point_tag(x_index + 1, y_index + 1),
                                                           point_tag(x_index, y_index + 1)});
            gmsh::model::geo::mesh::setRecombine(2, surface_tag);
            surface_tags.push_back(surface_tag);
        }
    }

    std::vector<int> inlet_curve_tags;
    std::vector<int> wall_curve_tags;
    std::vector<int> outlet_curve_tags;
    const double outlet_x{geometry.upstream_length + geometry.downstream_length};
    for (Index x_index = 0; x_index < x_segments.size(); ++x_index)
    {
        for (Index y_index = 0; y_index < y_point_count; ++y_index)
        {
            const int curve_tag{horizontal_curve_tags[x_index * y_point_count + y_index]};
            if (curve_tag == 0)
            {
                continue;
            }
            if (y_index + 1 == y_point_count ||
                (y_index == 0 && x_segments[x_index].begin >= geometry.upstream_length) ||
                (y_index == first_upper_segment && x_segments[x_index].end <= geometry.upstream_length))
            {
                wall_curve_tags.push_back(curve_tag);
            }
        }
    }
    for (Index x_index = 0; x_index < x_point_count; ++x_index)
    {
        for (Index y_index = 0; y_index < near_step_y_segments.size(); ++y_index)
        {
            const int curve_tag{vertical_curve_tags[x_index * near_step_y_segments.size() + y_index]};
            if (curve_tag == 0)
            {
                continue;
            }
            if (x_coordinates[x_index] == 0.0)
            {
                inlet_curve_tags.push_back(curve_tag);
            }
            else if (x_coordinates[x_index] == outlet_x)
            {
                outlet_curve_tags.push_back(curve_tag);
            }
            else if (x_coordinates[x_index] == geometry.upstream_length && y_index < first_upper_segment)
            {
                wall_curve_tags.push_back(curve_tag);
            }
        }
    }

    ModelTags tags{add_physical_groups(surface_tags, inlet_curve_tags, wall_curve_tags, outlet_curve_tags)};
    tags.transfinite_quadrilateral_mesh = true;
    return tags;
}

void configure_surface_mesh(const CellType cell_type, const std::vector<int> &surface_tags)
{
    switch (cell_type)
    {
    case CellType::Triangle:
        return;

    case CellType::Quadrilateral:
        // Blossom recombination generates a triangular background mesh and
        // recombines pairs of triangles into quadrilateral elements.
        gmsh::option::setNumber("Mesh.RecombinationAlgorithm", 1);
        for (const int surface_tag : surface_tags)
        {
            gmsh::model::mesh::setRecombine(2, surface_tag);
        }

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

    // Gmsh returns Cartesian coordinates as x/y/z triples even for a 2D model.
    if (coordinates.size() != 3 * gmsh_node_tags.size())
    {
        throw std::runtime_error("Invalid node coordinate data returned by Gmsh.");
    }

    raw_mesh.nodes.reserve(gmsh_node_tags.size());

    // Gmsh tags are external identifiers and are not assumed to be contiguous
    // or zero-based. Convert them once to compact IDs suitable for direct
    // indexing throughout the solver.
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
                   const std::vector<int> &surface_tags, const CellType requested_cell_type, const bool hybrid)
{
    const CellExtractionSpec extraction_spec{cell_extraction_spec(requested_cell_type)};
    if (extraction_spec.node_count == 0)
    {
        throw std::logic_error("Cell extraction specification has no nodes.");
    }

    // Validate every returned element block first and determine exact storage
    // requirements before filling the flattened connectivity arrays.
    Index total_cell_count{};
    Index total_cell_node_count{};
    std::vector<std::pair<CellExtractionSpec, std::vector<std::size_t>>> cell_connectivity_blocks;

    for (const int surface_tag : surface_tags)
    {
        std::vector<int> gmsh_element_types;
        std::vector<std::vector<std::size_t>> gmsh_element_tags;
        std::vector<std::vector<std::size_t>> gmsh_element_node_tags;

        // Query only the fluid surfaces explicitly created by the selected
        // geometry. Element tags are not retained because compact solver cell
        // IDs are assigned below.
        gmsh::model::mesh::getElements(gmsh_element_types, gmsh_element_tags, gmsh_element_node_tags, 2, surface_tag);
        for (std::size_t element_type_index = 0; element_type_index < gmsh_element_types.size(); ++element_type_index)
        {
            const int type{gmsh_element_types[element_type_index]};
            if ((!hybrid && type != extraction_spec.gmsh_element_type) ||
                (hybrid && type != gmsh_triangle_3_element_type && type != gmsh_quadrilateral_4_element_type))
            {
                throw std::runtime_error("Gmsh generated a 2D element type that does not "
                                         "match the requested cell type.");
            }

            std::vector<std::size_t> &gmsh_node_tags{gmsh_element_node_tags[element_type_index]};
            const auto spec{cell_extraction_spec(type == gmsh_triangle_3_element_type ? CellType::Triangle
                                                                                      : CellType::Quadrilateral)};

            if (gmsh_node_tags.size() % spec.node_count != 0)
            {
                throw std::runtime_error("Invalid cell connectivity returned by Gmsh.");
            }

            total_cell_count += gmsh_node_tags.size() / spec.node_count;
            total_cell_node_count += gmsh_node_tags.size();
            cell_connectivity_blocks.emplace_back(spec, std::move(gmsh_node_tags));
        }
    }

    raw_mesh.cell_types.reserve(total_cell_count);
    raw_mesh.cell_nodes.reserve(total_cell_node_count);
    raw_mesh.cell_node_offsets.reserve(total_cell_count + 1);

    raw_mesh.cell_node_offsets.push_back(0);

    for (const auto &[spec, gmsh_node_tags] : cell_connectivity_blocks)
    {
        for (std::size_t first_node_position = 0; first_node_position < gmsh_node_tags.size();
             first_node_position += spec.node_count)
        {
            raw_mesh.cell_types.push_back(spec.cell_type);

            for (Index local_node_index = 0; local_node_index < spec.node_count; ++local_node_index)
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

    // A physical boundary group contains geometric curves; boundary
    // connectivity is carried by the one-dimensional mesh elements belonging
    // to those curves.
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

[[nodiscard]]
RawMeshData generate_and_extract_mesh(const ModelTags &tags, const MeshGenerationOptions &options)
{
    if (!tags.transfinite_quadrilateral_mesh && !options.hybrid)
    {
        configure_surface_mesh(options.cell_type, tags.surface_tags);
    }

    // Generate the two-dimensional mesh. Gmsh also creates the lower-dimensional
    // boundary mesh needed for subsequent boundary-edge extraction.
    gmsh::model::mesh::generate(2);

    RawMeshData raw_mesh;

    // Expose stable solver boundary IDs rather than the physical-group tags
    // assigned internally by Gmsh.
    for (const PhysicalBoundaryTags &boundary : tags.boundaries)
    {
        raw_mesh.boundary_groups.push_back(boundary.group);
    }

    const auto node_id_by_gmsh_tag{extract_nodes(raw_mesh)};

    extract_cells(raw_mesh, node_id_by_gmsh_tag, tags.surface_tags, options.cell_type, options.hybrid);

    for (const PhysicalBoundaryTags &boundary : tags.boundaries)
    {
        extract_boundary_edges(raw_mesh, node_id_by_gmsh_tag, boundary.physical_tag, boundary.group.id);
    }

    return raw_mesh;
}

// RPA backend: all geometry/discretization decisions have already been made.
// Explicit grids share the same Gmsh boundary nodes; only TRANSITION is meshed
// by Gmsh. MeshOnlyEmpty prevents remeshing those supplied structured cells.
struct RpaCadPoint
{
    int tag;
    std::size_t node;
    Point2 position;
};

struct RpaCadCurve
{
    int tag;
    std::vector<std::size_t> nodes;
};

struct RpaNodeGroup
{
    int dimension, entity;
    std::vector<std::size_t> tags;
    std::vector<double> xyz;
};

struct RpaElementGroup
{
    int entity, type;
    std::vector<std::size_t> nodes;
};

class RpaGmshBuilder
{
  public:
    RpaCadPoint point(Point2 p)
    {
        const int tag{gmsh::model::geo::addPoint(p.x, p.y, 0.0)};
        const auto node{next_node_++};
        nodes_.push_back({0, tag, {node}, {p.x, p.y, 0.0}});
        return {tag, node, p};
    }

    RpaCadCurve curve(const RpaCadPoint &a, const RpaCadPoint &b, const std::span<const Point2> mesh,
                      const std::span<const Point2> cad = {})
    {
        std::vector<int> points{a.tag};
        for (const auto &p : cad)
            if (p.x > a.position.x && p.x < b.position.x)
                points.push_back(point(p).tag);
        points.push_back(b.tag);
        const int tag{points.size() == 2 ? gmsh::model::geo::addLine(a.tag, b.tag)
                                         : gmsh::model::geo::addPolyline(points)};
        RpaCadCurve result{tag, {a.node}};
        RpaNodeGroup group{1, tag, {}, {}};
        for (Index i = 1; i + 1 < mesh.size(); ++i)
        {
            result.nodes.push_back(next_node_++);
            group.tags.push_back(result.nodes.back());
            group.xyz.insert(group.xyz.end(), {mesh[i].x, mesh[i].y, 0.0});
        }
        result.nodes.push_back(b.node);
        nodes_.push_back(std::move(group));
        RpaElementGroup edges{tag, gmsh_line_2_element_type, {}};
        for (Index i = 1; i < result.nodes.size(); ++i)
            edges.nodes.insert(edges.nodes.end(), {result.nodes[i - 1], result.nodes[i]});
        elements_.push_back(std::move(edges));
        return result;
    }

    // Preserve source-entity junctions as separate CAD curves and mesh vertices.
    std::vector<RpaCadCurve> wall(const detail::RpaStructuredGrid &grid, Index row, const std::vector<Point2> &profile,
                                  const RpaNozzleDomain &domain, RpaCadPoint first, const RpaCadPoint &last,
                                  bool physical)
    {
        std::vector<RpaCadCurve> result;
        Index begin{};
        for (Index end = 1; end < grid.x.size(); ++end)
            if (end + 1 == grid.x.size() ||
                std::binary_search(domain.curve_junctions.begin(), domain.curve_junctions.end(), grid.x[end]))
            {
                const auto next{end + 1 == grid.x.size() ? last : point(grid.at(end, row))};
                std::vector<Point2> coordinates;
                for (Index i = begin; i <= end; ++i)
                    coordinates.push_back(grid.at(i, row));
                std::vector<Point2> cad;
                for (const auto &p : profile)
                    if (p.x > grid.x[begin] && p.x < grid.x[end])
                        cad.push_back({p.x, physical ? p.y
                                                     : p.y + (grid.at(begin, row).y -
                                                              interpolate_rpa_profile(profile, grid.x[begin]).y)});
                result.push_back(curve(first, next, coordinates, cad));
                begin = end;
                first = next;
            }
        return result;
    }

    RpaCadCurve edge(const detail::RpaStructuredGrid &grid, bool axial, Index at, const RpaCadPoint &a,
                     const RpaCadPoint &b)
    {
        std::vector<Point2> coordinates;
        const Index count{axial ? grid.x.size() : grid.radial_cells + 1};
        coordinates.reserve(count);
        for (Index i = 0; i < count; ++i)
            coordinates.push_back(axial ? grid.at(i, at) : grid.at(at, i));
        return curve(a, b, coordinates);
    }

    int surface(const std::vector<int> &loop)
    {
        const int tag{gmsh::model::geo::addPlaneSurface({gmsh::model::geo::addCurveLoop(loop)})};
        surfaces.push_back(tag);
        return tag;
    }

    void grid(const detail::RpaStructuredGrid &grid, const std::vector<RpaCadCurve> &bottom,
              const std::vector<RpaCadCurve> &left, const std::vector<RpaCadCurve> &right,
              const std::vector<RpaCadCurve> &top, const std::vector<int> &corners, bool triangles)
    {
        std::vector<int> loop;
        loop.reserve(bottom.size() + right.size() + top.size() + left.size());
        for (const auto &c : bottom)
            loop.push_back(c.tag);
        for (const auto &c : right)
            loop.push_back(c.tag);
        for (auto it = top.rbegin(); it != top.rend(); ++it)
            loop.push_back(-it->tag);
        for (auto it = left.rbegin(); it != left.rend(); ++it)
            loop.push_back(-it->tag);
        const int tag{surface(loop)};
        gmsh::model::geo::mesh::setTransfiniteSurface(tag, "Left", corners);
        const auto chain = [](const std::vector<RpaCadCurve> &curves) {
            std::vector<std::size_t> result;
            for (const auto &c : curves)
                result.insert(result.end(), c.nodes.begin() + (result.empty() ? 0 : 1), c.nodes.end());
            return result;
        };
        const auto down{chain(bottom)}, up{chain(top)}, start{chain(left)}, end{chain(right)};
        if (down.size() != grid.x.size() || up.size() != grid.x.size() || start.size() != grid.radial_cells + 1 ||
            start.size() != end.size())
            throw std::logic_error("RPA shared structured boundary cardinality mismatch.");
        std::vector<std::size_t> ids(grid.nodes.size());
        const Index stride{grid.radial_cells + 1};
        const auto assign = [&](Index i, std::size_t node) {
            if (ids[i] != 0 && ids[i] != node)
                throw std::logic_error("RPA corner is not a shared Gmsh node.");
            ids[i] = node;
        };
        for (Index i = 0; i < grid.x.size(); ++i)
        {
            assign(i * stride, down[i]);
            assign(i * stride + grid.radial_cells, up[i]);
        }
        for (Index j = 0; j < stride; ++j)
        {
            assign(j, start[j]);
            assign((grid.x.size() - 1) * stride + j, end[j]);
        }
        RpaNodeGroup group{2, tag, {}, {}};
        for (Index i = 0; i < grid.nodes.size(); ++i)
            if (ids[i] == 0)
            {
                ids[i] = next_node_++;
                group.tags.push_back(ids[i]);
                group.xyz.insert(group.xyz.end(), {grid.nodes[i].x, grid.nodes[i].y, 0.0});
            }
        nodes_.push_back(std::move(group));
        RpaElementGroup cells{tag, triangles ? gmsh_triangle_3_element_type : gmsh_quadrilateral_4_element_type, {}};
        for (Index i = 1; i < grid.x.size(); ++i)
            for (Index j = 0; j < grid.radial_cells; ++j)
            {
                const auto a{ids[(i - 1) * stride + j]}, b{ids[i * stride + j]}, c{ids[i * stride + j + 1]},
                    d{ids[(i - 1) * stride + j + 1]};
                if (triangles)
                    cells.nodes.insert(cells.nodes.end(), {a, b, c, a, c, d});
                else
                    cells.nodes.insert(cells.nodes.end(), {a, b, c, d});
            }
        elements_.push_back(std::move(cells));
    }

    void install()
    {
        gmsh::model::geo::synchronize();
        // Let Gmsh finish its native 1D stage so subsequent 2D generation does
        // not invalidate the supplied grids. Only the counts come from Gmsh;
        // the discretization owns every prescribed boundary coordinate.
        for (const auto &group : elements_)
            if (group.type == gmsh_line_2_element_type)
                gmsh::model::mesh::setTransfiniteCurve(group.entity, static_cast<int>(group.nodes.size() / 2 + 1));
        gmsh::model::mesh::generate(1);
        std::vector<std::size_t> actual(next_node_);
        std::size_t next{};
        for (const auto &group : nodes_)
            if (group.dimension < 2)
            {
                std::vector<std::size_t> tags;
                std::vector<double> xyz, parameters;
                gmsh::model::mesh::getNodes(tags, xyz, parameters, group.dimension, group.entity);
                if (tags.size() != group.tags.size())
                    throw std::runtime_error("RPA Gmsh boundary node count differs from the prescribed grid.");
                std::vector<Index> order(tags.size());
                std::iota(order.begin(), order.end(), 0);
                if (group.dimension == 1)
                    std::ranges::sort(order, [&](Index a, Index b) { return parameters[a] < parameters[b]; });
                for (Index i = 0; i < order.size(); ++i)
                {
                    const auto tag{tags[order[i]]};
                    actual[group.tags[i]] = tag;
                    next = std::max(next, tag);
                    const std::vector<double> coordinate{group.xyz[3 * i], group.xyz[3 * i + 1], group.xyz[3 * i + 2]};
                    std::vector<double> parameter;
                    if (group.dimension == 1)
                        gmsh::model::getParametrization(1, group.entity, coordinate, parameter);
                    gmsh::model::mesh::setNode(tag, coordinate, parameter);
                }
            }
        for (const auto &group : nodes_)
            if (group.dimension == 2 && !group.tags.empty())
            {
                std::vector<std::size_t> tags;
                tags.reserve(group.tags.size());
                for (const auto tag : group.tags)
                {
                    actual[tag] = ++next;
                    tags.push_back(next);
                }
                gmsh::model::mesh::addNodes(2, group.entity, tags, group.xyz);
            }
        for (const auto &group : elements_)
            if (group.type != gmsh_line_2_element_type)
            {
                std::vector<std::size_t> nodes;
                nodes.reserve(group.nodes.size());
                for (const auto tag : group.nodes)
                    nodes.push_back(actual[tag]);
                gmsh::model::mesh::addElementsByType(group.entity, group.type, {}, nodes);
            }
        gmsh::option::setNumber("Mesh.MeshOnlyEmpty", 1);
    }

    std::vector<int> surfaces;

  private:
    std::size_t next_node_{1};
    std::vector<RpaNodeGroup> nodes_;
    std::vector<RpaElementGroup> elements_;
};

// Coarsening FARFIELD can leave a poor full-QUAD recombination at the acute
// TRANSITION outlet. Relocate only free nodes near rejected cells, with an
// aspect-sensitive objective: angular smoothing alone can collapse edges.
// All prescribed grids and physical boundary nodes remain fixed.
void improve_rpa_transition_quality(RawMeshData &raw, Index first_transition, double minimum_quality)
{
    const auto cell_quality = [&](Index cell) {
        std::array<Point2, 4> corners{};
        for (Index j = 0; j < corners.size(); ++j)
            corners.at(j) = raw.nodes[raw.cell_nodes[raw.cell_node_offsets[cell] + j]];
        return detail::rpa_angular_quality(corners);
    };
    std::vector<bool> active(raw.nodes.size());
    bool needed{};
    for (Index cell = first_transition; cell < raw.cell_types.size(); ++cell)
        if (cell_quality(cell) < minimum_quality)
        {
            needed = true;
            for (Index j = raw.cell_node_offsets[cell]; j < raw.cell_node_offsets[cell + 1]; ++j)
                active[raw.cell_nodes[j]] = true;
        }
    if (!needed)
        return;

    // Include a small neighborhood so adjacent free nodes can redistribute
    // angles rather than moving only a subdivision's three-valent centre.
    for (Index ring = 0; ring < 3; ++ring)
    {
        const auto previous{active};
        for (Index cell = first_transition; cell < raw.cell_types.size(); ++cell)
        {
            bool touches{};
            for (Index j = raw.cell_node_offsets[cell]; j < raw.cell_node_offsets[cell + 1]; ++j)
                touches = touches || previous[raw.cell_nodes[j]];
            if (touches)
                for (Index j = raw.cell_node_offsets[cell]; j < raw.cell_node_offsets[cell + 1]; ++j)
                    active[raw.cell_nodes[j]] = true;
        }
    }
    std::vector<bool> pinned(raw.nodes.size());
    std::vector<Index> offsets(raw.nodes.size() + 1);
    for (Index cell = 0; cell < raw.cell_types.size(); ++cell)
        for (Index j = raw.cell_node_offsets[cell]; j < raw.cell_node_offsets[cell + 1]; ++j)
        {
            const Index node{raw.cell_nodes[j]};
            ++offsets[node + 1];
            if (cell < first_transition)
                pinned[node] = true;
        }
    for (const auto &edge : raw.boundary_edges)
        for (const auto node : edge.node_ids)
            pinned[node] = true;
    std::partial_sum(offsets.begin(), offsets.end(), offsets.begin());
    auto next{offsets};
    std::vector<Index> incident(raw.cell_nodes.size());
    for (Index cell = 0; cell < raw.cell_types.size(); ++cell)
        for (Index j = raw.cell_node_offsets[cell]; j < raw.cell_node_offsets[cell + 1]; ++j)
            incident[next[raw.cell_nodes[j]]++] = cell;
    std::vector<Index> movable;
    for (Index node = 0; node < raw.nodes.size(); ++node)
        if (active[node] && !pinned[node])
            movable.push_back(node);

    const auto energy = [&](Index node) {
        double result{};
        for (Index k = offsets[node]; k < offsets[node + 1]; ++k)
        {
            const Index begin{raw.cell_node_offsets[incident[k]]};
            for (Index j = 0; j < 4; ++j)
            {
                const auto &a{raw.nodes[raw.cell_nodes[begin + j]]}, &b{raw.nodes[raw.cell_nodes[begin + (j + 1) % 4]]},
                    &c{raw.nodes[raw.cell_nodes[begin + (j + 3) % 4]]};
                const double ax{b.x - a.x}, ay{b.y - a.y}, bx{c.x - a.x}, by{c.y - a.y};
                const double first{std::hypot(ax, ay)}, second{std::hypot(bx, by)};
                if (!(first > 0.0) || !(second > 0.0))
                    return std::numeric_limits<double>::infinity();
                const double ratio{std::min(first, second) / std::max(first, second)};
                const double sine{(ax / first) * (by / second) - (ay / first) * (bx / second)};
                const double shape{2.0 * ratio * sine / (1.0 + ratio * ratio)};
                if (!std::isfinite(shape) || !(shape > 0.0))
                    return std::numeric_limits<double>::infinity();
                result += 1.0 / (shape * shape);
            }
        }
        return result;
    };
    const double diagonal{std::sqrt(0.5)};
    const std::array<Point2, 8> directions{{{1, 0},
                                            {-1, 0},
                                            {0, 1},
                                            {0, -1},
                                            {diagonal, diagonal},
                                            {-diagonal, diagonal},
                                            {diagonal, -diagonal},
                                            {-diagonal, -diagonal}}};
    // Bounded local coordinate search. The final hard quality check remains
    // authoritative; neither geometry nor the requested threshold is relaxed.
    for (Index pass = 0; pass < 300; ++pass)
    {
        const Index step_reductions{pass / 30};
        const double fraction{0.2 * std::pow(0.5, static_cast<double>(step_reductions))};
        for (const auto node : movable)
        {
            const auto start{raw.nodes[node]};
            double step{std::numeric_limits<double>::infinity()};
            for (Index k = offsets[node]; k < offsets[node + 1]; ++k)
                for (Index j = raw.cell_node_offsets[incident[k]]; j < raw.cell_node_offsets[incident[k] + 1]; ++j)
                    if (raw.cell_nodes[j] != node)
                    {
                        const auto &p{raw.nodes[raw.cell_nodes[j]]};
                        step = std::min(step, std::hypot(p.x - start.x, p.y - start.y));
                    }
            step *= fraction;
            auto best{start};
            double best_energy{energy(node)};
            for (const auto &direction : directions)
            {
                raw.nodes[node] = {start.x + step * direction.x, start.y + step * direction.y};
                const double candidate{energy(node)};
                if (candidate < best_energy)
                {
                    best = raw.nodes[node];
                    best_energy = candidate;
                }
            }
            raw.nodes[node] = best;
        }
    }
}

ModelTags create_rpa_nozzle(const RpaNozzleDomain &domain, const detail::RpaMeshDiscretization &mesh,
                            const MeshGenerationOptions &options)
{
    RpaGmshBuilder b;
    gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
    gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);
    gmsh::option::setNumber("Mesh.MeshSizeMax", mesh.transition_size);
    const auto &e{mesh.engine}, &band{mesh.band}, &p{mesh.plume}, &f{mesh.farfield};
    const auto e00{b.point(e.at(0, 0))}, e10{b.point(e.at(e.x.size() - 1, 0))}, e01{b.point(e.at(0, e.radial_cells))},
        e11{b.point(e.at(e.x.size() - 1, e.radial_cells))};
    const auto eb{b.edge(e, true, 0, e00, e10)}, el{b.edge(e, false, 0, e00, e01)},
        er{b.edge(e, false, e.x.size() - 1, e10, e11)};
    const auto ew{b.wall(e, e.radial_cells, domain.inner_wall, domain, e01, e11, true)};
    const auto b00{b.point(band.at(0, 0))}, b10{b.point(band.at(band.x.size() - 1, 0))},
        b01{b.point(band.at(0, band.radial_cells))}, b11{b.point(band.at(band.x.size() - 1, band.radial_cells))};
    const auto bw{b.wall(band, 0, domain.outer_wall, domain, b00, b10, true)},
        bt{b.wall(band, band.radial_cells, domain.outer_wall, domain, b01, b11, false)};
    const auto bl{b.edge(band, false, 0, b00, b01)}, br{b.edge(band, false, band.x.size() - 1, b10, b11)};
    const std::array<Point2, 2> lip_points{e11.position, b10.position};
    const auto lip{b.curve(e11, b10, lip_points)};
    const auto p10{b.point(p.at(p.x.size() - 1, 0))}, p11{b.point(p.at(p.x.size() - 1, p.radial_cells))};
    const auto pb{b.edge(p, true, 0, e10, p10)}, pt{b.edge(p, true, p.radial_cells, b11, p11)},
        pr{b.edge(p, false, p.x.size() - 1, p10, p11)};
    const auto f00{b.point(f.at(0, 0))}, f01{b.point(f.at(0, f.radial_cells))},
        f11{b.point(f.at(f.x.size() - 1, f.radial_cells))};
    const auto fb{b.edge(f, true, 0, f00, p11)}, ft{b.edge(f, true, f.radial_cells, f01, f11)},
        fl{b.edge(f, false, 0, f00, f01)}, fr{b.edge(f, false, f.x.size() - 1, p11, f11)};
    std::vector<Point2> connection;
    for (Index j = 0; j <= mesh.transition_upstream_cells; ++j)
        connection.push_back(
            {b01.position.x, std::lerp(b01.position.y, f00.position.y,
                                       static_cast<double>(j) / static_cast<double>(mesh.transition_upstream_cells))});
    const auto upstream{b.curve(b01, f00, connection)};
    const bool triangles{options.cell_type == CellType::Triangle && !options.hybrid};
    b.grid(e, {eb}, {el}, {er}, ew, {e00.tag, e10.tag, e11.tag, e01.tag}, triangles);
    b.grid(band, bw, {bl}, {br}, bt, {b00.tag, b10.tag, b11.tag, b01.tag}, triangles);
    b.grid(p, {pb}, {er, lip, br}, {pr}, {pt}, {e10.tag, p10.tag, p11.tag, b11.tag}, triangles);
    b.grid(f, {fb}, {fl}, {fr}, {ft}, {f00.tag, p11.tag, f11.tag, f01.tag}, triangles);
    std::vector<int> loop{pt.tag, -fb.tag, -upstream.tag};
    for (const auto &c : bt)
        loop.push_back(c.tag);
    const int transition{b.surface(loop)};
    if (options.cell_type == CellType::Quadrilateral && !options.hybrid)
    {
        gmsh::model::geo::mesh::setRecombine(2, transition);
        gmsh::option::setNumber("Mesh.Algorithm", 5);
        gmsh::option::setNumber("Mesh.RecombinationAlgorithm", 3);
    }
    b.install();
    std::array<std::vector<int>, 8> boundaries{
        {{el.tag}, {eb.tag, pb.tag}, {}, {}, {lip.tag}, {bl.tag, upstream.tag, fl.tag}, {ft.tag}, {pr.tag, fr.tag}}};
    for (const auto &c : ew)
        boundaries.at(2).push_back(c.tag);
    for (const auto &c : bw)
        boundaries.at(3).push_back(c.tag);
    constexpr std::array<const char *, 8> names{"inlet",           "axis",
                                                "inner_wall",      "outer_wall",
                                                "nozzle_lip",      "farfield_upstream",
                                                "farfield_radial", "farfield_downstream"};
    ModelTags result{b.surfaces, {}, true};
    for (Index i = 0; i < boundaries.size(); ++i)
    {
        const int tag{gmsh::model::addPhysicalGroup(1, boundaries.at(i))};
        gmsh::model::setPhysicalName(1, tag, names.at(i));
        result.boundaries.push_back({{i, names.at(i)}, tag});
    }
    const int fluid{gmsh::model::addPhysicalGroup(2, result.surface_tags)};
    gmsh::model::setPhysicalName(2, fluid, "fluid");
    return result;
}

void require_growth_rate(const std::vector<double> &spacings, const double maximum_growth_rate,
                         const double ratio_tolerance, const char *const line_direction, const Index line_index)
{
    for (Index index = 1; index < spacings.size(); ++index)
    {
        const double ratio{std::max(spacings[index] / spacings[index - 1], spacings[index - 1] / spacings[index])};
        if (ratio > maximum_growth_rate + ratio_tolerance)
        {
            throw std::runtime_error("Generated structured mesh spacing ratio " + std::to_string(ratio) +
                                     " exceeds maximumGrowthRate " + std::to_string(maximum_growth_rate) + " on " +
                                     line_direction + " line " + std::to_string(line_index) + ".");
        }
    }
}

struct StructuredNodeColumn
{
    double x{};
    Index first_node{};
    Index node_count{};
};

void validate_structured_region_growth(const std::vector<Node> &nodes,
                                       const std::vector<StructuredNodeColumn> &all_columns, const double minimum_x,
                                       const double minimum_y, const double maximum_y,
                                       const double coordinate_tolerance, const double maximum_growth_rate,
                                       const double ratio_tolerance)
{
    std::vector<std::span<const Node>> columns;
    columns.reserve(all_columns.size());
    for (const StructuredNodeColumn &source_column : all_columns)
    {
        if (source_column.x < minimum_x - coordinate_tolerance)
        {
            continue;
        }

        Index first_node{source_column.first_node};
        const Index source_end{source_column.first_node + source_column.node_count};
        while (first_node < source_end && nodes[first_node].y < minimum_y - coordinate_tolerance)
        {
            ++first_node;
        }
        Index region_end{first_node};
        while (region_end < source_end && nodes[region_end].y <= maximum_y + coordinate_tolerance)
        {
            ++region_end;
        }

        const std::span<const Node> column{std::span<const Node>{nodes}.subspan(first_node, region_end - first_node)};
        if (column.size() < 2 || (!columns.empty() && column.size() != columns.front().size()))
        {
            throw std::runtime_error("Generated structured mesh has inconsistent logical columns.");
        }

        std::vector<double> vertical_spacings;
        vertical_spacings.reserve(column.size() - 1);
        for (Index node_index = 1; node_index < column.size(); ++node_index)
        {
            vertical_spacings.push_back(column[node_index].y - column[node_index - 1].y);
        }
        require_growth_rate(vertical_spacings, maximum_growth_rate, ratio_tolerance, "vertical", columns.size());
        columns.push_back(column);
    }

    for (Index row_index = 0; row_index < columns.front().size(); ++row_index)
    {
        std::vector<double> horizontal_spacings;
        horizontal_spacings.reserve(columns.size() - 1);
        for (Index column_index = 1; column_index < columns.size(); ++column_index)
        {
            const double delta_x{columns[column_index][row_index].x - columns[column_index - 1][row_index].x};
            const double delta_y{columns[column_index][row_index].y - columns[column_index - 1][row_index].y};
            horizontal_spacings.push_back(std::hypot(delta_x, delta_y));
        }
        require_growth_rate(horizontal_spacings, maximum_growth_rate, ratio_tolerance, "horizontal", row_index);
    }
}

void validate_structured_backward_facing_step_growth(const RawMeshData &raw_mesh,
                                                     const BackwardFacingStepGeometry &geometry,
                                                     const double maximum_growth_rate)
{
    const double coordinate_scale{
        std::max({1.0, geometry.upstream_length + geometry.downstream_length, geometry.channel_height})};
    const double coordinate_tolerance{std::sqrt(std::numeric_limits<double>::epsilon()) * coordinate_scale};
    const double ratio_tolerance{std::sqrt(std::numeric_limits<double>::epsilon()) * maximum_growth_rate};

    // Sort and group the generated nodes once. Subsequent upper/lower-region
    // validation traverses these structured columns without rescanning the
    // complete node array for every x-coordinate.
    std::vector<Node> sorted_nodes{raw_mesh.nodes};
    std::ranges::sort(sorted_nodes, [](const Node &left, const Node &right) {
        return left.x < right.x || (left.x == right.x && left.y < right.y);
    });

    std::vector<StructuredNodeColumn> columns;
    for (Index first_node = 0; first_node < sorted_nodes.size();)
    {
        Index column_end{first_node + 1};
        while (column_end < sorted_nodes.size() &&
               std::abs(sorted_nodes[column_end].x - sorted_nodes[first_node].x) <= coordinate_tolerance)
        {
            ++column_end;
        }
        std::ranges::sort(std::span<Node>{sorted_nodes}.subspan(first_node, column_end - first_node), {}, &Node::y);
        columns.push_back({sorted_nodes[first_node].x, first_node, column_end - first_node});
        first_node = column_end;
    }

    validate_structured_region_growth(sorted_nodes, columns, 0.0, geometry.step_height, geometry.channel_height,
                                      coordinate_tolerance, maximum_growth_rate, ratio_tolerance);
    validate_structured_region_growth(sorted_nodes, columns, geometry.upstream_length, 0.0, geometry.step_height,
                                      coordinate_tolerance, maximum_growth_rate, ratio_tolerance);
}

} // namespace

RawMeshData generate_mesh(const RectangleGeometry &geometry, const MeshGenerationOptions &options)
{
    validate_geometry(geometry);
    validate_mesh_generation_options(options);
    if (options.hybrid)
        throw std::invalid_argument("HYBRID meshing is supported only for RPA, not rectangle.");

    GmshSession gmsh_session;
    gmsh::option::setNumber("General.Terminal", 0);
    gmsh::model::add("rectangle");

    return generate_and_extract_mesh(create_rectangle(geometry, options.mesh_size), options);
}

RawMeshData generate_mesh(const BackwardFacingStepGeometry &geometry, const MeshGenerationOptions &options)
{
    return generate_mesh(geometry, options, AutomaticMeshingOptions{}, BackwardFacingStepMeshingOptions{});
}

RawMeshData generate_mesh(const BackwardFacingStepGeometry &geometry, const MeshGenerationOptions &options,
                          const AutomaticMeshingOptions &automatic_options,
                          const BackwardFacingStepMeshingOptions &step_options)
{
    validate_geometry(geometry);
    validate_mesh_generation_options(options);
    validate_automatic_meshing_options(automatic_options);
    validate_backward_facing_step_meshing_options(step_options);
    if (options.hybrid)
        throw std::invalid_argument("HYBRID meshing is supported only for RPA, not backward-facing step.");

    GmshSession gmsh_session;
    gmsh::option::setNumber("General.Terminal", 0);
    gmsh::model::add("backward_facing_step");

    if (options.cell_type == CellType::Quadrilateral && automatic_options.enabled)
    {
        RawMeshData raw_mesh{generate_and_extract_mesh(
            create_structured_backward_facing_step(geometry, options, automatic_options, step_options), options)};
        validate_structured_backward_facing_step_growth(raw_mesh, geometry, automatic_options.maximum_growth_rate);
        return raw_mesh;
    }

    return generate_and_extract_mesh(create_backward_facing_step(geometry, options.mesh_size), options);
}

RawMeshData generate_mesh(const GeometryInput &geometry, const MeshGenerationOptions &options)
{
    return generate_mesh(geometry, options, AutomaticMeshingOptions{}, BackwardFacingStepMeshingOptions{});
}

RawMeshData generate_mesh(const RpaNozzleGeometry &geometry, const MeshGenerationOptions &options,
                          const RpaNozzleMeshingOptions &rpa_options, const MeshSizeField &size_field)
{
    validate_mesh_generation_options(options);
    const auto domain{read_rpa_nozzle_domain(geometry)};
    const RpaRegionLayout layout{domain, rpa_options.regions};
    const MeshSizeField file_field{rpa_options.external_size_field_file
                                       ? read_mesh_size_field(*rpa_options.external_size_field_file)
                                       : MeshSizeField{}};
    MeshSizeField combined;
    if (file_field || size_field)
        combined = [&](const Point2 &p) {
            const double a{file_field ? file_field(p) : std::numeric_limits<double>::max()};
            const double b{size_field ? size_field(p) : std::numeric_limits<double>::max()};
            if (!std::isfinite(a) || !std::isfinite(b) || !(a > 0.0) || !(b > 0.0))
                throw std::invalid_argument("External RPA size field must be finite and positive.");
            return std::min(a, b);
        };
    const auto discretization{detail::discretize_rpa(domain, layout, options, rpa_options)};
    GmshSession session;
    gmsh::option::setNumber("General.Terminal", 0);
    gmsh::model::add("rpa_nozzle");
    try
    {
        const auto model{create_rpa_nozzle(domain, discretization, options)};
        std::atomic<bool> invalid_field{};
        gmsh::model::mesh::setSizeCallback([&](int, int, double x, double y, double, double size) {
            try
            {
                double target{discretization.transition_size};
                for (const auto *grid : {&discretization.band, &discretization.plume})
                {
                    const auto it{std::upper_bound(grid->x.begin(), grid->x.end(), x)};
                    const Index i{std::clamp(static_cast<Index>(it - grid->x.begin()), Index{1}, grid->x.size() - 1)};
                    const auto &a{grid->at(i - 1, grid->radial_cells)}, &b{grid->at(i, grid->radial_cells)};
                    const double length{std::hypot(b.x - a.x, b.y - a.y)};
                    const double ux{(b.x - a.x) / length}, uy{(b.y - a.y) / length};
                    const double along{std::clamp((x - a.x) * ux + (y - a.y) * uy, 0.0, length)};
                    const double distance{std::hypot(x - a.x - along * ux, y - a.y - along * uy)};
                    target = std::min(target, length + 0.2 * distance);
                }
                return combined ? std::min(target, combined({x, y})) : target;
            }
            catch (...)
            {
                invalid_field.store(true);
                return size;
            }
        });
        RawMeshData raw;
        try
        {
            raw = generate_and_extract_mesh(model, options);
        }
        catch (...)
        {
            gmsh::model::mesh::removeSizeCallback();
            throw;
        }
        gmsh::model::mesh::removeSizeCallback();
        if (invalid_field.load())
            throw std::invalid_argument("External RPA size field failed during Gmsh generation.");
        if (options.cell_type == CellType::Quadrilateral && !options.hybrid)
        {
            // Extraction preserves surface order. Verify each of the four
            // supplied surfaces before treating the remaining cells as free
            // TRANSITION: a changed backend order must fail, not move grids.
            if (model.surface_tags.size() != 5)
                throw std::logic_error("RPA expected four structured surfaces followed by TRANSITION.");
            Index first_transition{};
            Index surface{};
            for (const auto *grid :
                 {&discretization.engine, &discretization.band, &discretization.plume, &discretization.farfield})
            {
                const Index expected{(grid->x.size() - 1) * grid->radial_cells};
                std::vector<std::size_t> elements, nodes;
                gmsh::model::mesh::getElementsByType(gmsh_quadrilateral_4_element_type, elements, nodes,
                                                     model.surface_tags.at(surface++));
                if (elements.size() != expected || nodes.size() / 4 != expected || nodes.size() % 4 != 0 ||
                    first_transition > raw.cell_types.size() || expected > raw.cell_types.size() - first_transition)
                    throw std::logic_error("RPA structured surface counts do not match the extracted cell prefix.");
                first_transition += expected;
            }
            if (first_transition >= raw.cell_types.size() ||
                !std::all_of(raw.cell_types.begin(),
                             raw.cell_types.begin() + static_cast<std::ptrdiff_t>(first_transition),
                             [](CellType type) { return type == CellType::Quadrilateral; }))
                throw std::logic_error("RPA expected a QUAD structured prefix followed by a nonempty TRANSITION.");
            improve_rpa_transition_quality(raw, first_transition, rpa_options.min_cell_quality);
        }
        detail::require_rpa_mesh_quality(raw, rpa_options);
        if (rpa_options.target_cell_budget && raw.cell_types.size() > *rpa_options.target_cell_budget)
            throw std::runtime_error(
                "RPA targetCellBudget exceeded; geometric and quality constraints were not relaxed.");
        return raw;
    }
    catch (const std::string &error)
    {
        throw std::runtime_error("RPA Gmsh generation failed: " + error);
    }
}

RawMeshData generate_mesh(const GeometryInput &geometry, const MeshGenerationOptions &options,
                          const AutomaticMeshingOptions &automatic_options,
                          const BackwardFacingStepMeshingOptions &step_options,
                          const RpaNozzleMeshingOptions &rpa_options, const MeshSizeField &size_field)
{
    return std::visit(
        [&options, &automatic_options, &step_options, &rpa_options, &size_field](const auto &selected_geometry) {
            using Geometry = std::remove_cvref_t<decltype(selected_geometry)>;
            if constexpr (std::is_same_v<Geometry, BackwardFacingStepGeometry>)
            {
                if (size_field)
                    throw std::invalid_argument("External spatial size fields currently require rpaNozzle geometry.");
                return generate_mesh(selected_geometry, options, automatic_options, step_options);
            }
            else if constexpr (std::is_same_v<Geometry, RpaNozzleGeometry>)
            {
                const AutomaticMeshingOptions defaults;
                if (automatic_options.enabled != defaults.enabled ||
                    automatic_options.maximum_growth_rate != defaults.maximum_growth_rate ||
                    automatic_options.wall_refinement_factor != defaults.wall_refinement_factor ||
                    automatic_options.wall_refinement_layers != defaults.wall_refinement_layers)
                    throw std::invalid_argument(
                        "RPA automaticMeshing is removed; use baseCellSize and optional wall layers.");
                return generate_mesh(selected_geometry, options, rpa_options, size_field);
            }
            else
            {
                if (size_field)
                    throw std::invalid_argument("External spatial size fields currently require rpaNozzle geometry.");
                return generate_mesh(selected_geometry, options);
            }
        },
        geometry);
}

} // namespace cfd
