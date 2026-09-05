#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVectorField.hpp"
#include "cfd/field/FaceFluxField.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/linear_algebra/EigenBiCGSTABSolver.hpp"
#include "cfd/linear_algebra/ScalarLinearSystem.hpp"
#include "cfd/math/Point2.hpp"
#include "cfd/mesh/Boundary.hpp"
#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/mesh/Types.hpp"
#include "cfd/meshing/RawMeshData.hpp"
#include "cfd/numerics/LeastSquaresGradient.hpp"
#include "cfd/numerics/ScalarConvectionOperator.hpp"
#include "cfd/numerics/ScalarDiffusionOperator.hpp"

#include "support/VerificationStatistics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <span>
#include <sstream>
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

constexpr double domain_length{2.0};
constexpr double domain_height{1.0};
constexpr double velocity_x{1.0};
constexpr double velocity_y{0.5};
constexpr double diffusivity{0.1};
constexpr double algebraic_residual_tolerance{1.0e-11};
constexpr double full_residual_tolerance{1.0e-8};
constexpr double flux_imbalance_tolerance{1.0e-12};

struct GridLevel
{
    double delta;
    cfd::Index nx;
    cfd::Index ny;
};

constexpr std::array grid_levels{
    GridLevel{0.20, 10, 5},   GridLevel{0.10, 20, 10},    GridLevel{0.05, 40, 20},
    GridLevel{0.025, 80, 40}, GridLevel{0.0125, 160, 80},
};

struct LevelResult
{
    double delta{};
    double characteristic_mesh_length{};
    cfd::Index cell_count{};
    cfd::Index internal_face_count{};
    ErrorStatistics solution_errors;
    std::optional<double> rms_order;
    std::optional<double> linf_order;
    cfd::Index iteration_count{};
    double estimated_relative_error{};
    double normalized_algebraic_residual{};
    double normalized_full_residual{};
    double maximum_flux_imbalance{};
};

// Manufactured solution and derivatives:
//   phi = sin(x + 2y) + 0.2 x^2 + 0.3 xy + 0.1 y^2 + 1
//   grad(phi) = (cos(x + 2y) + 0.4x + 0.3y,
//                2 cos(x + 2y) + 0.3x + 0.2y)
//   laplacian(phi) = -5 sin(x + 2y) + 0.6
// For u=(1,0.5), Gamma=0.1, and div(u)=0, the conservative source is
//   s = 2 cos(x + 2y) + 0.55x + 0.4y + Gamma(5 sin(x + 2y) - 0.6).
[[nodiscard]]
double analytical_phi(const cfd::Point2 &point) noexcept
{
    return std::sin(point.x + 2.0 * point.y) + 0.2 * point.x * point.x + 0.3 * point.x * point.y +
           0.1 * point.y * point.y + 1.0;
}

[[nodiscard]]
double source_density(const cfd::Point2 &point) noexcept
{
    return 2.0 * std::cos(point.x + 2.0 * point.y) + 0.55 * point.x + 0.4 * point.y +
           diffusivity * (5.0 * std::sin(point.x + 2.0 * point.y) - 0.6);
}

[[nodiscard]]
cfd::Index node_id(const cfd::Index i, const cfd::Index j, const cfd::Index nx) noexcept
{
    return j * (nx + 1) + i;
}

void append_boundary_edge(cfd::RawMeshData &raw_mesh, const cfd::Index first_node, const cfd::Index second_node,
                          const std::string_view side)
{
    const cfd::BoundaryId boundary_id{raw_mesh.boundary_groups.size()};
    raw_mesh.boundary_groups.push_back({boundary_id, std::string{side} + "_" + std::to_string(boundary_id)});
    raw_mesh.boundary_edges.push_back({{first_node, second_node}, boundary_id});
}

[[nodiscard]]
cfd::RawMeshData make_raw_mesh(const GridLevel &level)
{
    constexpr double extent_tolerance{64.0 * std::numeric_limits<double>::epsilon()};
    if (std::abs(static_cast<double>(level.nx) * level.delta - domain_length) > extent_tolerance * domain_length ||
        std::abs(static_cast<double>(level.ny) * level.delta - domain_height) > extent_tolerance * domain_height)
    {
        throw std::runtime_error("Convection-diffusion grid dimensions do not match the domain.");
    }

    const cfd::Index node_count{(level.nx + 1) * (level.ny + 1)};
    const cfd::Index cell_count{level.nx * level.ny};
    cfd::RawMeshData raw_mesh;
    raw_mesh.nodes.reserve(node_count);
    raw_mesh.cell_types.reserve(cell_count);
    raw_mesh.cell_nodes.reserve(4 * cell_count);
    raw_mesh.cell_node_offsets.reserve(cell_count + 1);
    raw_mesh.cell_node_offsets.push_back(0);

    for (cfd::Index j = 0; j <= level.ny; ++j)
    {
        for (cfd::Index i = 0; i <= level.nx; ++i)
        {
            raw_mesh.nodes.push_back({static_cast<double>(i) * level.delta, static_cast<double>(j) * level.delta});
        }
    }

    for (cfd::Index j = 0; j < level.ny; ++j)
    {
        for (cfd::Index i = 0; i < level.nx; ++i)
        {
            const cfd::Index lower_left{node_id(i, j, level.nx)};
            const cfd::Index lower_right{node_id(i + 1, j, level.nx)};
            const cfd::Index upper_right{node_id(i + 1, j + 1, level.nx)};
            const cfd::Index upper_left{node_id(i, j + 1, level.nx)};
            raw_mesh.cell_types.push_back(cfd::CellType::Quadrilateral);
            raw_mesh.cell_nodes.insert(raw_mesh.cell_nodes.end(), {lower_left, lower_right, upper_right, upper_left});
            raw_mesh.cell_node_offsets.push_back(raw_mesh.cell_nodes.size());
        }
    }

    raw_mesh.boundary_groups.reserve(2 * (level.nx + level.ny));
    raw_mesh.boundary_edges.reserve(2 * (level.nx + level.ny));
    for (cfd::Index j = 0; j < level.ny; ++j)
    {
        append_boundary_edge(raw_mesh, node_id(0, j, level.nx), node_id(0, j + 1, level.nx), "left");
        append_boundary_edge(raw_mesh, node_id(level.nx, j, level.nx), node_id(level.nx, j + 1, level.nx), "right");
    }
    for (cfd::Index i = 0; i < level.nx; ++i)
    {
        append_boundary_edge(raw_mesh, node_id(i, 0, level.nx), node_id(i + 1, 0, level.nx), "bottom");
        append_boundary_edge(raw_mesh, node_id(i, level.ny, level.nx), node_id(i + 1, level.ny, level.nx), "top");
    }
    return raw_mesh;
}

[[nodiscard]]
cfd::ScalarBoundaryConditions make_boundary_conditions(const cfd::Mesh &mesh)
{
    std::vector<cfd::ScalarBoundaryCondition> conditions;
    conditions.reserve(mesh.boundary_groups().size());
    for (cfd::Index boundary_id = 0; boundary_id < mesh.boundary_groups().size(); ++boundary_id)
    {
        conditions.emplace_back(cfd::ScalarBoundaryConditionType::Dirichlet, 0.0);
    }
    std::vector<bool> assigned(mesh.boundary_groups().size(), false);

    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (!mesh.face_adjacencies()[face_id].is_boundary())
        {
            continue;
        }
        const cfd::BoundaryId boundary_id{mesh.face_boundary_ids()[face_id]};
        if (assigned[boundary_id])
        {
            throw std::runtime_error("A convection-diffusion BoundaryId is shared by multiple faces.");
        }
        assigned[boundary_id] = true;
        conditions[boundary_id] = {
            cfd::ScalarBoundaryConditionType::Dirichlet,
            analytical_phi(mesh.face_centers()[face_id]),
        };
    }

    if (std::find(assigned.begin(), assigned.end(), false) != assigned.end())
    {
        throw std::runtime_error("A convection-diffusion boundary group has no face.");
    }
    return {mesh.boundary_groups().size(), std::move(conditions)};
}

[[nodiscard]]
double euclidean_norm(const std::span<const double> values) noexcept
{
    double squared_norm{};
    for (const double value : values)
    {
        squared_norm += value * value;
    }
    return std::sqrt(squared_norm);
}

[[nodiscard]]
double normalized_algebraic_residual(const cfd::ScalarLinearSystem &system, const cfd::CellScalarField &solution,
                                     const std::span<const double> rhs, const std::span<double> matrix_product)
{
    system.apply_matrix(solution.values(), matrix_product);
    double squared_residual_norm{};
    for (cfd::Index cell_id = 0; cell_id < system.cell_count(); ++cell_id)
    {
        const double residual{matrix_product[cell_id] - rhs[cell_id]};
        squared_residual_norm += residual * residual;
    }
    const double rhs_norm{euclidean_norm(rhs)};
    if (!(rhs_norm > 0.0) || !std::isfinite(rhs_norm))
    {
        throw std::runtime_error("Convection-diffusion algebraic RHS norm is invalid.");
    }
    return std::sqrt(squared_residual_norm) / rhs_norm;
}

[[nodiscard]]
double normalized_full_residual(const cfd::Mesh &mesh, const cfd::ScalarConvectionOperator &convection,
                                const cfd::ScalarDiffusionOperator &diffusion,
                                const cfd::ScalarBoundaryConditions &boundary_conditions,
                                const cfd::FaceFluxField &face_flux, const cfd::CellScalarField &solution,
                                const cfd::CellVectorField &gradient, const std::span<const double> source_integrals,
                                cfd::CellScalarField &convection_balance, cfd::CellScalarField &diffusion_balance)
{
    convection.compute_flux_balance(solution, boundary_conditions, face_flux, convection_balance);
    diffusion.compute_flux_balance(solution, boundary_conditions, gradient, diffusion_balance);

    double source_squared_norm{};
    double residual_squared_norm{};
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        const double source{source_integrals[cell_id]};
        const double residual{convection_balance[cell_id] + diffusion_balance[cell_id] - source};
        source_squared_norm += source * source;
        residual_squared_norm += residual * residual;
    }
    if (!(source_squared_norm > 0.0) || !std::isfinite(source_squared_norm))
    {
        throw std::runtime_error("Convection-diffusion source norm is invalid.");
    }
    return std::sqrt(residual_squared_norm / source_squared_norm);
}

[[nodiscard]]
std::pair<cfd::FaceFluxField, double> make_face_flux(const cfd::Mesh &mesh)
{
    cfd::FaceFluxField face_flux{mesh.face_count()};
    std::vector<double> cell_flux_imbalance(mesh.cell_count());
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};
        const auto &area_vector{mesh.face_area_vectors()[face_id]};
        const double flux{velocity_x * area_vector.x + velocity_y * area_vector.y};
        face_flux[face_id] = flux;
        cell_flux_imbalance[adjacency.owner] += flux;
        if (!adjacency.is_boundary())
        {
            cell_flux_imbalance[adjacency.neighbor] -= flux;
        }
    }

    double maximum_imbalance{};
    for (const double imbalance : cell_flux_imbalance)
    {
        maximum_imbalance = std::max(maximum_imbalance, std::abs(imbalance));
    }
    return {std::move(face_flux), maximum_imbalance};
}

[[nodiscard]]
cfd::Index count_internal_faces(const cfd::Mesh &mesh) noexcept
{
    cfd::Index internal_face_count{};
    for (const cfd::FaceAdjacency &adjacency : mesh.face_adjacencies())
    {
        internal_face_count += static_cast<cfd::Index>(!adjacency.is_boundary());
    }
    return internal_face_count;
}

[[nodiscard]]
LevelResult run_level(const GridLevel &level)
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_raw_mesh(level))};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarBoundaryConditions boundary_conditions{make_boundary_conditions(mesh)};
    auto [face_flux, maximum_flux_imbalance]{make_face_flux(mesh)};
    const cfd::ScalarDiffusionOperator diffusion{mesh, diffusivity};
    const cfd::ScalarConvectionOperator convection{mesh};
    cfd::ScalarLinearSystem system{mesh};

    system.clear();
    diffusion.add_matrix_contributions(boundary_conditions, system);
    convection.add_matrix_contributions(boundary_conditions, face_flux, system);

    std::vector<double> source_integrals(mesh.cell_count());
    double total_area{};
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        const double cell_area{mesh.cell_areas()[cell_id]};
        source_integrals[cell_id] = cell_area * source_density(mesh.cell_centers()[cell_id]);
        system.rhs()[cell_id] = source_integrals[cell_id];
        total_area += cell_area;
    }
    diffusion.add_boundary_rhs(boundary_conditions, system.rhs());
    convection.add_boundary_rhs(boundary_conditions, face_flux, system.rhs());

    cfd::EigenBiCGSTABSolver solver{{1.0e-14, 10000}};
    solver.compute_matrix(system);
    cfd::CellScalarField numerical_solution{mesh.cell_count()};
    const cfd::LinearSolveResult solve_result{solver.solve(system.rhs(), numerical_solution.values())};
    if (!solve_result.converged)
    {
        throw std::runtime_error("BiCGSTAB did not converge for the convection-diffusion system.");
    }

    std::vector<double> matrix_product(mesh.cell_count());
    const double algebraic_residual{
        normalized_algebraic_residual(system, numerical_solution, system.rhs(), matrix_product)};

    cfd::CellVectorField numerical_gradient{mesh.cell_count()};
    cfd::compute_least_squares_gradient(mesh, numerical_solution, boundary_conditions, numerical_gradient);
    cfd::CellScalarField convection_balance{mesh.cell_count()};
    cfd::CellScalarField diffusion_balance{mesh.cell_count()};
    const double full_residual{normalized_full_residual(mesh, convection, diffusion, boundary_conditions, face_flux,
                                                        numerical_solution, numerical_gradient, source_integrals,
                                                        convection_balance, diffusion_balance)};

    if (algebraic_residual > algebraic_residual_tolerance || full_residual > full_residual_tolerance ||
        maximum_flux_imbalance > flux_imbalance_tolerance)
    {
        std::ostringstream message;
        message << "Convection-diffusion consistency diagnostic exceeded its tolerance at delta=" << level.delta
                << "; algebraic_residual=" << algebraic_residual << ", full_residual=" << full_residual
                << ", maximum_flux_imbalance=" << maximum_flux_imbalance << '.';
        throw std::runtime_error(message.str());
    }

    ErrorAccumulator error_accumulator;
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        error_accumulator.add(mesh.cell_areas()[cell_id],
                              std::abs(numerical_solution[cell_id] - analytical_phi(mesh.cell_centers()[cell_id])));
    }

    return {
        .delta = level.delta,
        .characteristic_mesh_length = std::sqrt(total_area / static_cast<double>(mesh.cell_count())),
        .cell_count = mesh.cell_count(),
        .internal_face_count = count_internal_faces(mesh),
        .solution_errors = error_accumulator.finish(),
        .rms_order = std::nullopt,
        .linf_order = std::nullopt,
        .iteration_count = solve_result.iteration_count,
        .estimated_relative_error = solve_result.estimated_relative_error,
        .normalized_algebraic_residual = algebraic_residual,
        .normalized_full_residual = full_residual,
        .maximum_flux_imbalance = maximum_flux_imbalance,
    };
}

void compute_orders(std::vector<LevelResult> &results)
{
    for (std::size_t level_index = 1; level_index < results.size(); ++level_index)
    {
        const LevelResult &coarse{results[level_index - 1]};
        LevelResult &fine{results[level_index]};
        fine.rms_order =
            observed_order(coarse.solution_errors.area_weighted_rms_error, fine.solution_errors.area_weighted_rms_error,
                           coarse.characteristic_mesh_length, fine.characteristic_mesh_length);
        fine.linf_order = observed_order(coarse.solution_errors.linf_error, fine.solution_errors.linf_error,
                                         coarse.characteristic_mesh_length, fine.characteristic_mesh_length);
    }
}

void require_expected_convergence(const std::vector<LevelResult> &results)
{
    for (std::size_t level_index = 1; level_index < results.size(); ++level_index)
    {
        const LevelResult &coarse{results[level_index - 1]};
        const LevelResult &fine{results[level_index]};
        if (!(fine.solution_errors.area_weighted_rms_error < coarse.solution_errors.area_weighted_rms_error) ||
            !(fine.solution_errors.linf_error < coarse.solution_errors.linf_error))
        {
            throw std::runtime_error("Convection-diffusion solution errors did not decrease monotonically.");
        }
    }

    const LevelResult &finest{results.back()};
    if (!finest.rms_order.has_value() || !finest.linf_order.has_value() || !(*finest.rms_order > 0.5) ||
        !(*finest.rms_order < 1.5) || !(*finest.linf_order > 0.5) || !(*finest.linf_order < 1.5))
    {
        throw std::runtime_error("Convection-diffusion asymptotic order is not consistent with first-order upwind.");
    }
}

void print_order(const std::optional<double> order, const bool first_level)
{
    if (order.has_value())
    {
        std::cout << std::fixed << std::setprecision(3) << std::setw(9) << *order;
        return;
    }
    std::cout << std::setw(9) << (first_level ? "-" : "n/a");
}

void print_results(const std::vector<LevelResult> &results)
{
    std::cout << "\nScalar convection-diffusion solution convergence - Cartesian QUAD - all Dirichlet\n"
              << "div(u phi) - div(Gamma grad(phi)) = s, u=(1,0.5), Gamma=0.1\n\n"
              << std::left << std::setw(10) << "delta" << std::setw(12) << "h_char" << std::setw(9) << "cells"
              << std::setw(11) << "int_faces" << std::setw(17) << "RMS(phi)" << std::setw(9) << "p_RMS" << std::setw(17)
              << "Linf(phi)" << std::setw(9) << "p_Linf" << std::setw(12) << "iterations" << std::setw(15)
              << "Eigen_error" << std::setw(15) << "alg_res" << std::setw(15) << "full_res" << std::setw(15)
              << "max_divF" << '\n'
              << std::string(185, '-') << '\n';

    for (std::size_t level_index = 0; level_index < results.size(); ++level_index)
    {
        const LevelResult &result{results[level_index]};
        std::cout << std::fixed << std::setprecision(4) << std::setw(10) << result.delta << std::scientific
                  << std::setprecision(4) << std::setw(12) << result.characteristic_mesh_length << std::fixed
                  << std::setw(9) << result.cell_count << std::setw(11) << result.internal_face_count << std::scientific
                  << std::setprecision(7) << std::setw(17) << result.solution_errors.area_weighted_rms_error;
        print_order(result.rms_order, level_index == 0);
        std::cout << std::scientific << std::setprecision(7) << std::setw(17) << result.solution_errors.linf_error;
        print_order(result.linf_order, level_index == 0);
        std::cout << std::fixed << std::setw(12) << result.iteration_count << std::scientific << std::setprecision(4)
                  << std::setw(15) << result.estimated_relative_error << std::setw(15)
                  << result.normalized_algebraic_residual << std::setw(15) << result.normalized_full_residual
                  << std::setw(15) << result.maximum_flux_imbalance << '\n';
    }

    std::cout << "\nScope: production diffusion and first-order upwind convection assemble consistently; the directed "
                 "nonsymmetric system is solved by BiCGSTAB; exact constant face flux is discretely conservative; "
                 "Cartesian solution error approaches first order because of upwind convection.\n";
}

} // namespace

int main()
{
    try
    {
        std::vector<LevelResult> results;
        results.reserve(grid_levels.size());
        for (const GridLevel &level : grid_levels)
        {
            results.push_back(run_level(level));
        }
        compute_orders(results);
        require_expected_convergence(results);
        print_results(results);
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "Scalar convection-diffusion solution verification failed: " << error.what() << '\n';
        return 1;
    }
}
