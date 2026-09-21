#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVelocityField.hpp"
#include "cfd/field/FaceFluxField.hpp"
#include "cfd/field/PressureCorrectionBoundaryConditions.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/math/Point2.hpp"
#include "cfd/math/Vector2.hpp"
#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/meshing/RawMeshData.hpp"
#include "cfd/numerics/IncompressibleSimpleSolver.hpp"

#include "support/KovasznayAnalytical.hpp"
#include "support/VerificationStatistics.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace
{

using cfd::verification::ErrorAccumulator;
using cfd::verification::ErrorStatistics;
using cfd::verification::observed_order;
using cfd::verification::kovasznay::analytical_pressure;
using cfd::verification::kovasznay::analytical_pressure_gradient;
using cfd::verification::kovasznay::analytical_u;
using cfd::verification::kovasznay::analytical_v;
using cfd::verification::kovasznay::density;
using cfd::verification::kovasznay::domain_height;
using cfd::verification::kovasznay::domain_length;
using cfd::verification::kovasznay::dynamic_viscosity;
using cfd::verification::kovasznay::lambda;
using cfd::verification::kovasznay::reynolds_number;
using cfd::verification::kovasznay::verify_analytical_solution;
using cfd::verification::kovasznay::wave_number;

constexpr std::array<cfd::Index, 4> grid_levels{8, 16, 32, 64};
constexpr double minimum_final_rms_order{1.5};
constexpr double boundary_flux_relative_tolerance{1024.0 * std::numeric_limits<double>::epsilon()};

[[nodiscard]]
cfd::Index node_id(const cfd::Index i, const cfd::Index j, const cfd::Index n) noexcept
{
    return j * (n + 1) + i;
}

void append_boundary_edge(cfd::RawMeshData &raw_mesh, const cfd::Index first_node, const cfd::Index second_node,
                          const std::string_view side)
{
    const cfd::BoundaryId boundary_id{raw_mesh.boundary_groups.size()};
    raw_mesh.boundary_groups.push_back({boundary_id, std::string{side} + "_" + std::to_string(boundary_id)});
    raw_mesh.boundary_edges.push_back({{first_node, second_node}, boundary_id});
}

[[nodiscard]]
cfd::RawMeshData make_raw_mesh(const cfd::Index n)
{
    const double delta_x{domain_length / static_cast<double>(n)};
    const double delta_y{domain_height / static_cast<double>(n)};
    cfd::RawMeshData raw_mesh;
    raw_mesh.nodes.reserve((n + 1) * (n + 1));
    raw_mesh.cell_types.reserve(n * n);
    raw_mesh.cell_nodes.reserve(4 * n * n);
    raw_mesh.cell_node_offsets.reserve(n * n + 1);
    raw_mesh.cell_node_offsets.push_back(0);
    for (cfd::Index j = 0; j <= n; ++j)
    {
        for (cfd::Index i = 0; i <= n; ++i)
        {
            raw_mesh.nodes.push_back({static_cast<double>(i) * delta_x, static_cast<double>(j) * delta_y});
        }
    }
    for (cfd::Index j = 0; j < n; ++j)
    {
        for (cfd::Index i = 0; i < n; ++i)
        {
            raw_mesh.cell_types.push_back(cfd::CellType::Quadrilateral);
            raw_mesh.cell_nodes.insert(raw_mesh.cell_nodes.end(), {node_id(i, j, n), node_id(i + 1, j, n),
                                                                   node_id(i + 1, j + 1, n), node_id(i, j + 1, n)});
            raw_mesh.cell_node_offsets.push_back(raw_mesh.cell_nodes.size());
        }
    }
    raw_mesh.boundary_groups.reserve(4 * n);
    raw_mesh.boundary_edges.reserve(4 * n);
    for (cfd::Index j = 0; j < n; ++j)
    {
        append_boundary_edge(raw_mesh, node_id(0, j, n), node_id(0, j + 1, n), "left");
        append_boundary_edge(raw_mesh, node_id(n, j, n), node_id(n, j + 1, n), "right");
    }
    for (cfd::Index i = 0; i < n; ++i)
    {
        append_boundary_edge(raw_mesh, node_id(i, 0, n), node_id(i + 1, 0, n), "bottom");
        append_boundary_edge(raw_mesh, node_id(i, n, n), node_id(i + 1, n, n), "top");
    }
    return raw_mesh;
}

struct BoundaryData
{
    cfd::ScalarBoundaryConditions u;
    cfd::ScalarBoundaryConditions v;
    cfd::ScalarBoundaryConditions pressure;
    cfd::PressureCorrectionBoundaryConditions pressure_correction;
};

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
    std::vector<bool> assigned(boundary_count, false);
    for (cfd::Index boundary_id = 0; boundary_id < boundary_count; ++boundary_id)
    {
        u_conditions.emplace_back(cfd::ScalarBoundaryConditionType::Dirichlet, 0.0);
        v_conditions.emplace_back(cfd::ScalarBoundaryConditionType::Dirichlet, 0.0);
        pressure_conditions.emplace_back(cfd::ScalarBoundaryConditionType::Neumann, 0.0);
    }
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (!mesh.face_adjacencies()[face_id].is_boundary())
        {
            continue;
        }
        const cfd::BoundaryId boundary_id{mesh.face_boundary_ids()[face_id]};
        if (assigned[boundary_id])
        {
            throw std::runtime_error("A Kovasznay BoundaryId is shared by multiple faces.");
        }
        assigned[boundary_id] = true;
        const cfd::Point2 &point{mesh.face_centers()[face_id]};
        const cfd::Vector2 &area_vector{mesh.face_area_vectors()[face_id]};
        const cfd::Vector2 pressure_gradient{analytical_pressure_gradient(point)};
        const double inverse_face_length{1.0 / mesh.face_lengths()[face_id]};
        u_conditions[boundary_id] = {cfd::ScalarBoundaryConditionType::Dirichlet, analytical_u(point)};
        v_conditions[boundary_id] = {cfd::ScalarBoundaryConditionType::Dirichlet, analytical_v(point)};
        pressure_conditions[boundary_id] = {
            cfd::ScalarBoundaryConditionType::Neumann,
            (pressure_gradient.x * area_vector.x + pressure_gradient.y * area_vector.y) * inverse_face_length};
    }
    if (std::find(assigned.begin(), assigned.end(), false) != assigned.end())
    {
        throw std::runtime_error("A Kovasznay boundary group has no face.");
    }
    return {cfd::ScalarBoundaryConditions{boundary_count, std::move(u_conditions)},
            cfd::ScalarBoundaryConditions{boundary_count, std::move(v_conditions)},
            cfd::ScalarBoundaryConditions{boundary_count, std::move(pressure_conditions)},
            cfd::PressureCorrectionBoundaryConditions{
                boundary_count, std::vector<cfd::PressureCorrectionBoundaryConditionType>{
                                    boundary_count, cfd::PressureCorrectionBoundaryConditionType::FixedMassFlux}}};
}

[[nodiscard]]
double vertical_flux_antiderivative(const double x, const double y) noexcept
{
    return y - std::exp(lambda * x) * std::sin(wave_number * y) / wave_number;
}

[[nodiscard]]
double horizontal_flux_antiderivative(const double x, const double y) noexcept
{
    return std::exp(lambda * x) * std::sin(wave_number * y) / wave_number;
}

[[nodiscard]]
double analytical_boundary_flux(const cfd::Mesh &mesh, const cfd::Index face_id)
{
    const cfd::Face &face{mesh.faces()[face_id]};
    const cfd::Point2 &first{mesh.nodes()[face.node_ids[0]]};
    const cfd::Point2 &second{mesh.nodes()[face.node_ids[1]]};
    const cfd::Vector2 &area_vector{mesh.face_area_vectors()[face_id]};
    if (first.x == second.x)
    {
        return density * area_vector.x / mesh.face_lengths()[face_id] *
               (vertical_flux_antiderivative(first.x, std::max(first.y, second.y)) -
                vertical_flux_antiderivative(first.x, std::min(first.y, second.y)));
    }
    if (first.y == second.y)
    {
        return density * area_vector.y / mesh.face_lengths()[face_id] *
               (horizontal_flux_antiderivative(std::max(first.x, second.x), first.y) -
                horizontal_flux_antiderivative(std::min(first.x, second.x), first.y));
    }
    throw std::runtime_error("The Kovasznay boundary flux requires axis-aligned faces.");
}

[[nodiscard]]
cfd::FaceFluxField make_initial_mass_flux(const cfd::Mesh &mesh)
{
    cfd::FaceFluxField flux{mesh.face_count()};
    long double boundary_flux_sum{};
    long double boundary_flux_magnitude{};
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (!mesh.face_adjacencies()[face_id].is_boundary())
        {
            continue;
        }
        flux[face_id] = analytical_boundary_flux(mesh, face_id);
        boundary_flux_sum += flux[face_id];
        boundary_flux_magnitude += std::abs(flux[face_id]);
    }
    if (std::abs(boundary_flux_sum) > 1.0e-12L * boundary_flux_magnitude)
    {
        throw std::runtime_error("The analytical Kovasznay boundary flux is not globally compatible.");
    }
    return flux;
}

[[nodiscard]]
cfd::IncompressibleSimpleOptions simple_options()
{
    return {
        .maximum_iterations = 5000,
        .momentum_relaxation_factor = 1.0,
        .pressure_relaxation_factor = 0.1,
        .rhie_chow_flux_relaxation_factor = 0.1,
        .velocity_relative_tolerance = 1.0e-10,
        .rhie_chow_flux_relative_tolerance = 1.0e-10,
        .continuity_relative_tolerance = 1.0e-10,
        .momentum_linear_solver = {.relative_tolerance = 1.0e-12, .maximum_iterations = 5000},
        .pressure_correction_linear_solver = {.relative_tolerance = 1.0e-12, .maximum_iterations = 5000},
    };
}

struct LevelResult
{
    cfd::Index n{};
    cfd::Index cell_count{};
    cfd::IncompressibleSimpleResult simple;
    ErrorStatistics u;
    ErrorStatistics v;
    ErrorStatistics pressure;
    double boundary_flux_relative_error{};
    double boundary_flux_relative_imbalance{};
};

[[nodiscard]]
LevelResult run_level(const cfd::Index n, const cfd::IncompressibleSimpleOptions &options)
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_raw_mesh(n))};
    const cfd::Mesh &mesh{build_result.mesh};
    const BoundaryData boundary{make_boundary_data(mesh)};
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    cfd::FaceFluxField flux{make_initial_mass_flux(mesh)};
    cfd::IncompressibleSimpleSolver solver{mesh, density, dynamic_viscosity, cfd::ScalarConvectionScheme::Linear,
                                           options};
    const cfd::IncompressibleSimpleResult simple{solver.solve(boundary.u, boundary.v, boundary.pressure,
                                                              boundary.pressure_correction, velocity, pressure, flux)};

    double total_area{};
    double pressure_difference_integral{};
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        const double area{mesh.cell_areas()[cell_id]};
        total_area += area;
        pressure_difference_integral += area * (pressure[cell_id] - analytical_pressure(mesh.cell_centers()[cell_id]));
    }
    // Remove the area-weighted mean difference: all-Neumann pressure has an arbitrary gauge.
    const double pressure_offset{pressure_difference_integral / total_area};
    ErrorAccumulator u_errors;
    ErrorAccumulator v_errors;
    ErrorAccumulator pressure_errors;
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        const cfd::Point2 &point{mesh.cell_centers()[cell_id]};
        const double area{mesh.cell_areas()[cell_id]};
        u_errors.add(area, velocity.u()[cell_id] - analytical_u(point));
        v_errors.add(area, velocity.v()[cell_id] - analytical_v(point));
        pressure_errors.add(area, pressure[cell_id] - analytical_pressure(point) - pressure_offset);
    }

    double boundary_flux_error{};
    double boundary_flux_sum{};
    double boundary_flux_magnitude{};
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (!mesh.face_adjacencies()[face_id].is_boundary())
        {
            continue;
        }
        const double exact_flux{analytical_boundary_flux(mesh, face_id)};
        boundary_flux_error += std::abs(flux[face_id] - exact_flux);
        boundary_flux_sum += flux[face_id];
        boundary_flux_magnitude += std::abs(exact_flux);
    }

    return {.n = n,
            .cell_count = mesh.cell_count(),
            .simple = simple,
            .u = u_errors.finish(),
            .v = v_errors.finish(),
            .pressure = pressure_errors.finish(),
            .boundary_flux_relative_error = boundary_flux_error / boundary_flux_magnitude,
            .boundary_flux_relative_imbalance = std::abs(boundary_flux_sum) / boundary_flux_magnitude};
}

void print_error_row(const LevelResult &level, const LevelResult *const previous, const std::string_view field,
                     const ErrorStatistics &errors, const ErrorStatistics *const previous_errors)
{
    std::optional<double> rms_order;
    std::optional<double> linf_order;
    if (previous != nullptr && previous_errors != nullptr)
    {
        rms_order = observed_order(previous_errors->area_weighted_rms_error, errors.area_weighted_rms_error,
                                   domain_length / static_cast<double>(previous->n),
                                   domain_length / static_cast<double>(level.n));
        linf_order = observed_order(previous_errors->linf_error, errors.linf_error,
                                    domain_length / static_cast<double>(previous->n),
                                    domain_length / static_cast<double>(level.n));
    }
    std::cout << std::setw(3) << level.n << ' ' << std::setw(5) << level.cell_count << ' ' << std::setw(2) << field
              << ' ' << std::scientific << std::setprecision(6) << std::setw(13) << errors.area_weighted_rms_error
              << ' ' << std::setw(13) << errors.linf_error << ' ';
    if (previous == nullptr)
    {
        std::cout << "   -       -";
    }
    else
    {
        if (rms_order)
        {
            std::cout << std::fixed << std::setprecision(3) << std::setw(6) << *rms_order;
        }
        else
        {
            std::cout << "   n/a";
        }
        if (linf_order)
        {
            std::cout << std::fixed << std::setprecision(3) << std::setw(8) << *linf_order;
        }
        else
        {
            std::cout << "     n/a";
        }
    }
    std::cout << '\n';
}

void require_strict_decrease(const double coarse_error, const double fine_error, const std::string_view quantity)
{
    if (!(fine_error < coarse_error))
    {
        throw std::runtime_error("Kovasznay " + std::string{quantity} + " did not decrease under refinement.");
    }
}

void require_final_rms_order(const ErrorStatistics &coarse, const ErrorStatistics &fine, const std::string_view field)
{
    const std::optional<double> order{
        observed_order(coarse.area_weighted_rms_error, fine.area_weighted_rms_error,
                       domain_length / static_cast<double>(grid_levels[grid_levels.size() - 2]),
                       domain_length / static_cast<double>(grid_levels.back()))};
    if (!order || !(*order > minimum_final_rms_order))
    {
        throw std::runtime_error("Kovasznay final " + std::string{field} + " RMS order did not exceed 1.5.");
    }
}

void verify_convergence(const std::vector<LevelResult> &levels)
{
    if (levels.size() != grid_levels.size())
    {
        throw std::runtime_error("Kovasznay verification did not produce every configured grid level.");
    }

    for (const LevelResult &level : levels)
    {
        if (level.boundary_flux_relative_error > boundary_flux_relative_tolerance ||
            level.boundary_flux_relative_imbalance > boundary_flux_relative_tolerance)
        {
            throw std::runtime_error("Kovasznay boundary mass-flux verification exceeded its tolerance.");
        }
    }

    for (std::size_t level_index = 1; level_index < levels.size(); ++level_index)
    {
        const LevelResult &coarse{levels[level_index - 1]};
        const LevelResult &fine{levels[level_index]};
        require_strict_decrease(coarse.u.area_weighted_rms_error, fine.u.area_weighted_rms_error, "u RMS error");
        require_strict_decrease(coarse.u.linf_error, fine.u.linf_error, "u Linf error");
        require_strict_decrease(coarse.v.area_weighted_rms_error, fine.v.area_weighted_rms_error, "v RMS error");
        require_strict_decrease(coarse.v.linf_error, fine.v.linf_error, "v Linf error");
        require_strict_decrease(coarse.pressure.area_weighted_rms_error, fine.pressure.area_weighted_rms_error,
                                "pressure RMS error");
        require_strict_decrease(coarse.pressure.linf_error, fine.pressure.linf_error, "pressure Linf error");
    }

    const LevelResult &coarse{levels[levels.size() - 2]};
    const LevelResult &fine{levels.back()};
    require_final_rms_order(coarse.u, fine.u, "u");
    require_final_rms_order(coarse.v, fine.v, "v");
    require_final_rms_order(coarse.pressure, fine.pressure, "pressure");
}

} // namespace

int main()
{
    try
    {
        const auto campaign_start{std::chrono::steady_clock::now()};
        verify_analytical_solution();
        const cfd::IncompressibleSimpleOptions options{simple_options()};
        std::cout << "Kovasznay SIMPLE: Re=" << reynolds_number
                  << ", Linear, alpha_u=" << options.momentum_relaxation_factor
                  << ", alpha_p=" << options.pressure_relaxation_factor
                  << ", alpha_rc=" << options.rhie_chow_flux_relaxation_factor << '\n';
        std::cout << "N cells field       RMS error    Linf error  order RMS  order Linf\n";
        std::vector<LevelResult> levels;
        levels.reserve(grid_levels.size());
        for (const cfd::Index n : grid_levels)
        {
            levels.push_back(run_level(n, options));
            const LevelResult &level{levels.back()};
            const LevelResult *const previous{levels.size() > 1 ? &levels[levels.size() - 2] : nullptr};
            print_error_row(level, previous, "u", level.u, previous != nullptr ? &previous->u : nullptr);
            print_error_row(level, previous, "v", level.v, previous != nullptr ? &previous->v : nullptr);
            print_error_row(level, previous, "p", level.pressure, previous != nullptr ? &previous->pressure : nullptr);
            std::cout << "  SIMPLE=" << (level.simple.converged ? "yes" : "no")
                      << " iterations=" << level.simple.iteration_count << std::scientific << std::setprecision(3)
                      << " r_U=" << level.simple.velocity_relative_change
                      << " r_RC=" << level.simple.rhie_chow_flux_relative_residual
                      << " continuity*=" << level.simple.provisional_continuity_relative_residual
                      << " continuity=" << level.simple.continuity_relative_residual
                      << " boundary_flux_error=" << level.boundary_flux_relative_error
                      << " boundary_imbalance=" << level.boundary_flux_relative_imbalance << '\n';
            if (!level.simple.converged)
            {
                throw std::runtime_error("Kovasznay SIMPLE did not converge at N=" + std::to_string(n));
            }
        }
        verify_convergence(levels);
        std::cout << "Acceptance checks: passed\n";
        std::cout << "Campaign seconds="
                  << std::chrono::duration<double>(std::chrono::steady_clock::now() - campaign_start).count() << '\n';
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "Kovasznay SIMPLE verification failed: " << error.what() << '\n';
        return 1;
    }
}
