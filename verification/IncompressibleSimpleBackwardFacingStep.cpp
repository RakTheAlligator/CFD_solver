#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVelocityField.hpp"
#include "cfd/field/FaceFluxField.hpp"
#include "cfd/field/PressureCorrectionBoundaryConditions.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/math/Point2.hpp"
#include "cfd/math/Vector2.hpp"
#include "cfd/mesh/Boundary.hpp"
#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/mesh/MeshStatistics.hpp"
#include "cfd/mesh/Types.hpp"
#include "cfd/meshing/BackwardFacingStepGeometry.hpp"
#include "cfd/meshing/GmshMesher.hpp"
#include "cfd/meshing/RawMeshData.hpp"
#include "cfd/numerics/IncompressibleSimpleSolver.hpp"
#include "cfd/numerics/ScalarConvectionOperator.hpp"

#include "support/BackwardFacingStepReattachment.hpp"
#include "support/GridConvergenceIndex.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace
{

constexpr double downstream_height{1.0};
constexpr double inlet_height{0.52 / 1.01};
constexpr double step_height{downstream_height - inlet_height};
constexpr double upstream_length{20.0 / 1.01};
constexpr double downstream_length{25.0};
constexpr double step_x{upstream_length};
constexpr double reynolds_number{100.0};
constexpr double mean_inlet_velocity{1.0};
constexpr double density{1.0};
constexpr double kinematic_viscosity{2.0 * inlet_height * mean_inlet_velocity / reynolds_number};
constexpr double dynamic_viscosity{density * kinematic_viscosity};
constexpr double inlet_flow_relative_tolerance{1.0e-12};

constexpr cfd::BoundaryId wall_boundary_id{0};
constexpr cfd::BoundaryId outlet_boundary_id{1};

constexpr std::array<double, 14> armaly_stations{
    0.00, 2.55, 3.06, 3.57, 4.18, 4.80, 5.41, 6.12, 7.76, 12.04, 16.33, 20.92, 30.31, 44.90,
};

struct GridDefinition
{
    std::string_view name;
    double nominal_size;
    double wall_refinement_factor;
    double step_refinement_factor;
};

constexpr std::array orthogonal_grids{
    GridDefinition{"coarse", 0.10, 1.0, 1.0},
    GridDefinition{"medium", 0.05, 1.0, 1.0},
    GridDefinition{"fine", 0.025, 1.0, 1.0},
};

constexpr GridDefinition refined_grid{"fine_D", 0.025, 0.5, 0.5};

struct BoundaryData
{
    cfd::ScalarBoundaryConditions u;
    cfd::ScalarBoundaryConditions v;
    cfd::ScalarBoundaryConditions pressure;
    cfd::PressureCorrectionBoundaryConditions pressure_correction;
};

struct LevelResult
{
    std::string name;
    cfd::Index cell_count{};
    cfd::Index face_count{};
    double effective_grid_size{};
    double maximum_non_orthogonality{};
    double maximum_neighbor_size_ratio{};
    cfd::IncompressibleSimpleResult simple;
    double reattachment_length_over_step{};
    double inlet_volumetric_flow{};
    double inlet_flow_relative_error{};
    double outlet_volumetric_flow{};
    double maximum_velocity_magnitude{};
    std::filesystem::path profile_file;
};

void split_inlet_boundary_groups(cfd::RawMeshData &raw_mesh)
{
    std::vector<std::string> original_names;
    original_names.reserve(raw_mesh.boundary_groups.size());
    for (const cfd::BoundaryGroup &group : raw_mesh.boundary_groups)
    {
        original_names.push_back(group.name);
    }

    std::vector<cfd::BoundaryGroup> groups{
        {wall_boundary_id, "wall"},
        {outlet_boundary_id, "outlet"},
    };
    cfd::Index inlet_face_number{};
    for (cfd::BoundaryEdge &edge : raw_mesh.boundary_edges)
    {
        const std::string &name{original_names.at(edge.boundary_id)};
        if (name == "wall")
        {
            edge.boundary_id = wall_boundary_id;
        }
        else if (name == "outlet")
        {
            edge.boundary_id = outlet_boundary_id;
        }
        else if (name == "inlet")
        {
            const cfd::BoundaryId boundary_id{groups.size()};
            groups.push_back({boundary_id, "inlet_" + std::to_string(inlet_face_number)});
            edge.boundary_id = boundary_id;
            ++inlet_face_number;
        }
        else
        {
            throw std::runtime_error("Armaly mesh contains an unexpected boundary group.");
        }
    }
    if (inlet_face_number == 0)
    {
        throw std::runtime_error("Armaly mesh contains no inlet face.");
    }
    raw_mesh.boundary_groups = std::move(groups);
}

[[nodiscard]]
cfd::MeshBuildResult make_mesh(const GridDefinition &grid)
{
    constexpr cfd::BackwardFacingStepGeometry geometry{
        .upstream_length = upstream_length,
        .downstream_length = downstream_length,
        .channel_height = downstream_height,
        .step_height = step_height,
    };
    const cfd::MeshGenerationOptions generation{
        .mesh_size = grid.nominal_size,
        .cell_type = cfd::CellType::Quadrilateral,
    };
    const cfd::AutomaticMeshingOptions automatic{
        .enabled = true,
        .maximum_growth_rate = 1.2,
        .wall_refinement_factor = grid.wall_refinement_factor,
        .wall_refinement_layers = 4,
    };
    const cfd::BackwardFacingStepMeshingOptions step{
        .step_refinement_factor = grid.step_refinement_factor,
        .step_refinement_layers = 6,
    };
    cfd::RawMeshData raw_mesh{cfd::generate_mesh(geometry, generation, automatic, step)};
    split_inlet_boundary_groups(raw_mesh);
    return cfd::build_mesh(std::move(raw_mesh));
}

[[nodiscard]]
BoundaryData make_boundary_data(const cfd::Mesh &mesh)
{
    const cfd::Index boundary_count{mesh.boundary_groups().size()};
    std::vector<cfd::ScalarBoundaryCondition> u_conditions;
    std::vector<cfd::ScalarBoundaryCondition> v_conditions;
    std::vector<cfd::ScalarBoundaryCondition> pressure_conditions;
    u_conditions.reserve(boundary_count);
    v_conditions.reserve(boundary_count);
    pressure_conditions.reserve(boundary_count);
    for (cfd::Index boundary_id = 0; boundary_id < boundary_count; ++boundary_id)
    {
        u_conditions.emplace_back(cfd::ScalarBoundaryConditionType::Dirichlet, 0.0);
        v_conditions.emplace_back(cfd::ScalarBoundaryConditionType::Dirichlet, 0.0);
        pressure_conditions.emplace_back(cfd::ScalarBoundaryConditionType::Neumann, 0.0);
    }
    std::vector<cfd::PressureCorrectionBoundaryConditionType> pressure_correction_conditions(
        boundary_count, cfd::PressureCorrectionBoundaryConditionType::FixedMassFlux);
    std::vector<bool> assigned(boundary_count, false);

    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (!mesh.face_adjacencies()[face_id].is_boundary())
        {
            continue;
        }
        const cfd::BoundaryId boundary_id{mesh.face_boundary_ids()[face_id]};
        const std::string &name{mesh.boundary_groups()[boundary_id].name};
        if (name.starts_with("inlet_"))
        {
            if (assigned[boundary_id])
            {
                throw std::runtime_error("An Armaly inlet BoundaryId is shared by multiple faces.");
            }
            const cfd::Face &face{mesh.faces()[face_id]};
            const cfd::Point2 &first_node{mesh.nodes()[face.node_ids[0]]};
            const cfd::Point2 &second_node{mesh.nodes()[face.node_ids[1]]};
            const double value{cfd::verification::parabolic_inlet_face_average(first_node.y, second_node.y, step_height,
                                                                               inlet_height, mean_inlet_velocity)};
            u_conditions[boundary_id] = {cfd::ScalarBoundaryConditionType::Dirichlet, value};
            v_conditions[boundary_id] = {cfd::ScalarBoundaryConditionType::Dirichlet, 0.0};
            pressure_conditions[boundary_id] = {cfd::ScalarBoundaryConditionType::Neumann, 0.0};
        }
        else if (name == "wall")
        {
            u_conditions[boundary_id] = {cfd::ScalarBoundaryConditionType::Dirichlet, 0.0};
            v_conditions[boundary_id] = {cfd::ScalarBoundaryConditionType::Dirichlet, 0.0};
            pressure_conditions[boundary_id] = {cfd::ScalarBoundaryConditionType::Neumann, 0.0};
        }
        else if (name == "outlet")
        {
            u_conditions[boundary_id] = {cfd::ScalarBoundaryConditionType::Neumann, 0.0};
            v_conditions[boundary_id] = {cfd::ScalarBoundaryConditionType::Neumann, 0.0};
            pressure_conditions[boundary_id] = {cfd::ScalarBoundaryConditionType::Dirichlet, 0.0};
            pressure_correction_conditions[boundary_id] = cfd::PressureCorrectionBoundaryConditionType::FixedPressure;
        }
        else
        {
            throw std::runtime_error("Armaly mesh contains an unexpected solver boundary group.");
        }
        assigned[boundary_id] = true;
    }
    if (std::find(assigned.begin(), assigned.end(), false) != assigned.end())
    {
        throw std::runtime_error("An Armaly boundary group has no face.");
    }

    return {
        cfd::ScalarBoundaryConditions{boundary_count, std::move(u_conditions)},
        cfd::ScalarBoundaryConditions{boundary_count, std::move(v_conditions)},
        cfd::ScalarBoundaryConditions{boundary_count, std::move(pressure_conditions)},
        cfd::PressureCorrectionBoundaryConditions{boundary_count, std::move(pressure_correction_conditions)},
    };
}

[[nodiscard]]
cfd::FaceFluxField make_initial_mass_flux(const cfd::Mesh &mesh, const BoundaryData &boundary)
{
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};
        if (!adjacency.is_boundary())
        {
            continue;
        }
        const cfd::BoundaryId boundary_id{mesh.face_boundary_ids()[face_id]};
        if (boundary.pressure_correction[boundary_id] == cfd::PressureCorrectionBoundaryConditionType::FixedPressure)
        {
            continue;
        }
        const cfd::Vector2 &area_vector{mesh.face_area_vectors()[face_id]};
        mass_flux[face_id] =
            density * (boundary.u[boundary_id].value * area_vector.x + boundary.v[boundary_id].value * area_vector.y);
    }
    return mass_flux;
}

[[nodiscard]]
double inlet_volumetric_flow(const cfd::Mesh &mesh, const BoundaryData &boundary)
{
    double flow{};
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (!mesh.face_adjacencies()[face_id].is_boundary())
        {
            continue;
        }
        const cfd::BoundaryId boundary_id{mesh.face_boundary_ids()[face_id]};
        if (!mesh.boundary_groups()[boundary_id].name.starts_with("inlet_"))
        {
            continue;
        }
        const cfd::Vector2 &area_vector{mesh.face_area_vectors()[face_id]};
        flow -= boundary.u[boundary_id].value * area_vector.x + boundary.v[boundary_id].value * area_vector.y;
    }
    return flow;
}

[[nodiscard]]
cfd::IncompressibleSimpleOptions simple_options(const double tolerance)
{
    return {
        .maximum_iterations = 10'000,
        .momentum_relaxation_factor = 1.0,
        .pressure_relaxation_factor = 0.10,
        .rhie_chow_flux_relaxation_factor = 0.10,
        .velocity_relative_tolerance = tolerance,
        .rhie_chow_flux_relative_tolerance = tolerance,
        .continuity_relative_tolerance = tolerance,
        .momentum_linear_solver = {.relative_tolerance = 1.0e-10, .maximum_iterations = 5000},
        .pressure_correction_linear_solver = {.relative_tolerance = 1.0e-10, .maximum_iterations = 5000},
    };
}

[[nodiscard]]
double outlet_volumetric_flow(const cfd::Mesh &mesh, const cfd::FaceFluxField &mass_flux)
{
    double outlet_mass_flow{};
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (mesh.face_adjacencies()[face_id].is_boundary() && mesh.face_boundary_ids()[face_id] == outlet_boundary_id)
        {
            outlet_mass_flow += mass_flux[face_id];
        }
    }
    return outlet_mass_flow / density;
}

[[nodiscard]]
double maximum_velocity_magnitude(const cfd::CellVelocityField &velocity)
{
    double maximum{};
    for (cfd::Index cell_id = 0; cell_id < velocity.size(); ++cell_id)
    {
        maximum = std::max(maximum, std::hypot(velocity.u()[cell_id], velocity.v()[cell_id]));
    }
    return maximum;
}

[[nodiscard]]
double reattachment_length_over_step(const cfd::Mesh &mesh, const cfd::CellVelocityField &velocity)
{
    constexpr double coordinate_tolerance{256.0 * std::numeric_limits<double>::epsilon()};
    std::vector<cfd::verification::NearWallVelocitySample> samples;
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};
        if (!adjacency.is_boundary() || mesh.face_boundary_ids()[face_id] != wall_boundary_id)
        {
            continue;
        }
        const cfd::Point2 &face_center{mesh.face_centers()[face_id]};
        if (std::abs(face_center.y) <= coordinate_tolerance && face_center.x > step_x)
        {
            samples.push_back({mesh.cell_centers()[adjacency.owner].x, velocity.u()[adjacency.owner]});
        }
    }
    std::sort(samples.begin(), samples.end(), [](const auto &left, const auto &right) { return left.x < right.x; });
    const double reattachment_x{cfd::verification::estimate_primary_reattachment_x(samples)};
    return (reattachment_x - step_x) / step_height;
}

using ProfileSamples = std::vector<std::pair<double, double>>;

[[nodiscard]]
double interpolate_in_y(const ProfileSamples &samples, const double y)
{
    const auto upper{
        std::lower_bound(samples.begin(), samples.end(), y,
                         [](const auto &sample, const double coordinate) { return sample.first < coordinate; })};
    if (upper != samples.end() && upper->first == y)
    {
        return upper->second;
    }
    if (upper == samples.begin() || upper == samples.end())
    {
        throw std::runtime_error("Armaly profile interpolation lacks a wall-normal bracket.");
    }
    const auto lower{std::prev(upper)};
    const double weight{(y - lower->first) / (upper->first - lower->first)};
    return lower->second + weight * (upper->second - lower->second);
}

void write_armaly_profiles(const cfd::Mesh &mesh, const cfd::CellVelocityField &velocity,
                           const std::filesystem::path &output_file)
{
    std::map<double, ProfileSamples> columns;
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        const cfd::Point2 &center{mesh.cell_centers()[cell_id]};
        if (center.x > step_x)
        {
            columns[center.x].push_back({center.y, velocity.u()[cell_id]});
        }
    }
    if (columns.size() < 2)
    {
        throw std::runtime_error("Armaly profile extraction requires at least two downstream cell columns.");
    }
    for (auto &[x, samples] : columns)
    {
        static_cast<void>(x);
        std::sort(samples.begin(), samples.end());
    }

    std::ofstream output{output_file};
    if (!output)
    {
        throw std::runtime_error("Cannot create Armaly profile CSV: " + output_file.string());
    }
    output << "x_over_S,y_over_H,u_over_V\n" << std::setprecision(17);

    for (const double station : armaly_stations)
    {
        const double target_x{step_x + station * step_height};
        auto upper_column{columns.lower_bound(target_x)};
        auto lower_column{upper_column};
        if (upper_column == columns.begin())
        {
            lower_column = upper_column;
            ++upper_column;
        }
        else if (upper_column == columns.end())
        {
            upper_column = std::prev(columns.end());
            lower_column = std::prev(upper_column);
        }
        else if (upper_column->first == target_x)
        {
            lower_column = upper_column;
        }
        else
        {
            lower_column = std::prev(upper_column);
        }

        const double lower_x{lower_column->first};
        const double upper_x{upper_column->first};
        const ProfileSamples &lower_samples{lower_column->second};
        const ProfileSamples &upper_samples{upper_column->second};
        const double minimum_common_y{std::max(lower_samples.front().first, upper_samples.front().first)};
        const double maximum_common_y{std::min(lower_samples.back().first, upper_samples.back().first)};
        for (const auto &[y, lower_u] : lower_samples)
        {
            if (y < minimum_common_y || y > maximum_common_y)
            {
                continue;
            }
            double interpolated_u{lower_u};
            if (lower_column != upper_column)
            {
                const double upper_u{interpolate_in_y(upper_samples, y)};
                const double weight{(target_x - lower_x) / (upper_x - lower_x)};
                interpolated_u = lower_u + weight * (upper_u - lower_u);
            }
            output << station << ',' << y / downstream_height << ',' << interpolated_u / mean_inlet_velocity << '\n';
        }
    }
}

[[nodiscard]]
LevelResult run_level(const GridDefinition &grid, const double tolerance,
                      const std::optional<std::filesystem::path> &profile_file)
{
    cfd::MeshBuildResult build_result{make_mesh(grid)};
    const cfd::Mesh &mesh{build_result.mesh};
    const BoundaryData boundary{make_boundary_data(mesh)};
    const cfd::MeshStatistics statistics{cfd::compute_mesh_statistics(mesh)};
    double total_area{};
    for (const double area : mesh.cell_areas())
    {
        total_area += area;
    }
    const double effective_grid_size{std::sqrt(total_area / static_cast<double>(mesh.cell_count()))};
    const double inlet_flow{inlet_volumetric_flow(mesh, boundary)};
    const double expected_inlet_flow{mean_inlet_velocity * inlet_height};
    const double inlet_flow_relative_error{std::abs(inlet_flow - expected_inlet_flow) / expected_inlet_flow};
    if (inlet_flow_relative_error > inlet_flow_relative_tolerance)
    {
        throw std::runtime_error("Discrete Armaly inlet flow does not match the prescribed analytical flow.");
    }
    std::cout << "Grid setup " << grid.name << ": cells=" << mesh.cell_count() << " faces=" << mesh.face_count()
              << std::scientific << std::setprecision(8) << " h_eff=" << effective_grid_size << " Qin=" << inlet_flow
              << " Qin_relative_error=" << inlet_flow_relative_error << '\n';
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    cfd::FaceFluxField mass_flux{make_initial_mass_flux(mesh, boundary)};
    const cfd::IncompressibleSimpleOptions options{simple_options(tolerance)};
    cfd::IncompressibleSimpleSolver solver{mesh, density, dynamic_viscosity,
                                           cfd::ScalarConvectionScheme::FirstOrderUpwind, options};
    cfd::Index last_completed_iteration{};
    cfd::IncompressibleSimpleResult simple;
    try
    {
        simple =
            solver.solve(boundary.u, boundary.v, boundary.pressure, boundary.pressure_correction, velocity, pressure,
                         mass_flux, [&last_completed_iteration](const cfd::SimpleIterationInfo &information) {
                             last_completed_iteration = information.iteration;
                         });
    }
    catch (const std::exception &error)
    {
        throw std::runtime_error("Armaly SIMPLE failed on grid " + std::string{grid.name} + " (" +
                                 std::to_string(mesh.cell_count()) + " cells) after " +
                                 std::to_string(last_completed_iteration) + " completed iterations: " + error.what());
    }
    if (!simple.converged)
    {
        throw std::runtime_error("Armaly SIMPLE did not converge on grid " + std::string{grid.name} + " after " +
                                 std::to_string(simple.iteration_count) + " iterations.");
    }

    if (profile_file)
    {
        write_armaly_profiles(mesh, velocity, *profile_file);
    }

    return {
        .name = std::string{grid.name},
        .cell_count = mesh.cell_count(),
        .face_count = mesh.face_count(),
        .effective_grid_size = effective_grid_size,
        .maximum_non_orthogonality = statistics.internal_face_non_orthogonality_degrees.maximum,
        .maximum_neighbor_size_ratio = statistics.internal_face_neighbor_cell_size_ratios.maximum,
        .simple = simple,
        .reattachment_length_over_step = reattachment_length_over_step(mesh, velocity),
        .inlet_volumetric_flow = inlet_flow,
        .inlet_flow_relative_error = inlet_flow_relative_error,
        .outlet_volumetric_flow = outlet_volumetric_flow(mesh, mass_flux),
        .maximum_velocity_magnitude = maximum_velocity_magnitude(velocity),
        .profile_file = profile_file.value_or(std::filesystem::path{}),
    };
}

[[nodiscard]]
double relative_difference(const double reference, const double candidate)
{
    if (reference == 0.0)
    {
        return candidate == 0.0 ? 0.0 : std::numeric_limits<double>::infinity();
    }
    return std::abs(candidate - reference) / std::abs(reference);
}

[[nodiscard]]
std::string_view convergence_behavior_name(const cfd::verification::GridConvergenceBehavior behavior) noexcept
{
    switch (behavior)
    {
    case cfd::verification::GridConvergenceBehavior::Monotonic:
        return "monotonic";
    case cfd::verification::GridConvergenceBehavior::Oscillatory:
        return "oscillatory";
    }
    return "unknown";
}

void print_level(const LevelResult &level)
{
    std::cout << std::left << std::setw(9) << level.name << std::right << std::setw(8) << level.cell_count
              << std::setw(9) << level.face_count << std::scientific << std::setprecision(5) << std::setw(13)
              << level.effective_grid_size << std::fixed << std::setprecision(3) << std::setw(11)
              << level.maximum_non_orthogonality << std::setw(11) << level.maximum_neighbor_size_ratio << std::setw(9)
              << level.simple.iteration_count << std::setprecision(6) << std::setw(12)
              << level.reattachment_length_over_step << std::scientific << std::setprecision(5) << std::setw(13)
              << level.inlet_volumetric_flow << std::setw(13) << level.inlet_flow_relative_error << std::setw(13)
              << level.outlet_volumetric_flow << std::setw(13) << level.maximum_velocity_magnitude << '\n';
}

} // namespace

int main()
{
    try
    {
        const std::filesystem::path output_directory{"output/verification/backward_facing_step"};
        std::filesystem::create_directories(output_directory);

        std::cout << "Backward-facing-step verification — Armaly Re=100\n\n"
                  << "Geometry:\n"
                  << "  H=" << downstream_height << " h_inlet=" << inlet_height << " S=" << step_height
                  << " H/h=" << downstream_height / inlet_height << "\n"
                  << "  upstream=" << upstream_length << " downstream=" << downstream_length << " x_step=" << step_x
                  << "\n"
                  << "  rho=" << density << " mu=" << dynamic_viscosity << " V=" << mean_inlet_velocity << "\n\n";

        std::cout << "Grid and solution:\n"
                  << "name        cells    faces        h_eff  nonorth°  sizeRatio    iter       xR/S         Qin"
                     "      QinErr         Qout         Umax\n";

        std::vector<LevelResult> levels;
        levels.reserve(orthogonal_grids.size());
        for (const GridDefinition &grid : orthogonal_grids)
        {
            const std::filesystem::path profile_file{output_directory /
                                                     ("armaly_profiles_" + std::string{grid.name} + ".csv")};
            levels.push_back(run_level(grid, 1.0e-8, profile_file));
            print_level(levels.back());
        }

        const LevelResult &coarse{levels[0]};
        const LevelResult &medium{levels[1]};
        const LevelResult &fine{levels[2]};
        const double r21{medium.effective_grid_size / fine.effective_grid_size};
        const double r32{coarse.effective_grid_size / medium.effective_grid_size};
        const std::optional<cfd::verification::GridConvergenceIndexResult> gci{
            cfd::verification::compute_grid_convergence_index(fine.reattachment_length_over_step,
                                                              medium.reattachment_length_over_step,
                                                              coarse.reattachment_length_over_step, r21, r32)};

        std::cout << "\nGrid convergence:\n"
                  << std::scientific << std::setprecision(8) << "  r21=" << r21 << " r32=" << r32 << '\n';
        if (gci)
        {
            std::cout << "  behavior=" << convergence_behavior_name(gci->behavior) << " p=" << gci->apparent_order
                      << " Richardson=" << gci->extrapolated_value << " ea21=" << gci->approximate_relative_error
                      << " eext21=" << gci->extrapolated_relative_error
                      << " GCI_fine=" << gci->fine_grid_convergence_index << '\n';
        }
        else
        {
            std::cout << "  GCI unavailable: nearly identical differences or ill-conditioned apparent order.\n";
        }

        const LevelResult medium_tight{run_level(orthogonal_grids[1], 1.0e-10, std::nullopt)};
        std::cout << "\nIteration sensitivity (medium):\n"
                  << "  iterations 1e-8/1e-10=" << medium.simple.iteration_count << '/'
                  << medium_tight.simple.iteration_count << '\n'
                  << "  xR/S 1e-8/1e-10=" << medium.reattachment_length_over_step << '/'
                  << medium_tight.reattachment_length_over_step << " relative_difference="
                  << relative_difference(medium_tight.reattachment_length_over_step,
                                         medium.reattachment_length_over_step)
                  << '\n'
                  << "  Qin relative_difference="
                  << relative_difference(medium_tight.inlet_volumetric_flow, medium.inlet_volumetric_flow) << '\n'
                  << "  Qout relative_difference="
                  << relative_difference(medium_tight.outlet_volumetric_flow, medium.outlet_volumetric_flow) << '\n'
                  << "  Umax relative_difference="
                  << relative_difference(medium_tight.maximum_velocity_magnitude, medium.maximum_velocity_magnitude)
                  << '\n';

        const std::filesystem::path refined_profile_file{output_directory / "armaly_profiles_fine_D.csv"};
        const LevelResult refined{run_level(refined_grid, 1.0e-8, refined_profile_file)};
        std::cout << "\nRefined automatic mesh D:\n";
        print_level(fine);
        print_level(refined);
        std::cout << "  relative xR/S difference="
                  << relative_difference(fine.reattachment_length_over_step, refined.reattachment_length_over_step)
                  << " Qout difference="
                  << relative_difference(fine.outlet_volumetric_flow, refined.outlet_volumetric_flow)
                  << " Umax difference="
                  << relative_difference(fine.maximum_velocity_magnitude, refined.maximum_velocity_magnitude) << '\n';

        std::cout << "\nCSV profiles:\n";
        for (const LevelResult &level : levels)
        {
            std::cout << "  " << level.profile_file.string() << '\n';
        }
        std::cout << "  " << refined.profile_file.string() << '\n'
                  << "  x/S=0 uses one-sided linear extrapolation from the first two downstream cell columns.\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "Backward-facing-step verification failed: " << error.what() << '\n';
        return 1;
    }
}
