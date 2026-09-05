#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVectorField.hpp"
#include "cfd/field/CellVelocityField.hpp"
#include "cfd/field/FaceFluxField.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/linear_algebra/EigenBiCGSTABSolver.hpp"
#include "cfd/linear_algebra/ScalarLinearSystem.hpp"
#include "cfd/math/Point2.hpp"
#include "cfd/math/Vector2.hpp"
#include "cfd/mesh/Boundary.hpp"
#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/mesh/Types.hpp"
#include "cfd/meshing/RawMeshData.hpp"
#include "cfd/numerics/IncompressibleMomentumAssembler.hpp"
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
#include <numbers>
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

constexpr double domain_length{1.0};
constexpr double domain_height{1.0};
constexpr double density{1.0};
constexpr double reynolds_number{40.0};
constexpr double dynamic_viscosity{1.0 / reynolds_number};
constexpr double wave_number{2.0 * std::numbers::pi_v<double>};
const double kovasznay_lambda{reynolds_number / 2.0 -
                              std::sqrt(reynolds_number * reynolds_number / 4.0 + wave_number * wave_number)};

constexpr double algebraic_residual_tolerance{1.0e-10};
constexpr double full_residual_tolerance{1.0e-8};
constexpr double mass_imbalance_tolerance{1.0e-12};
constexpr double fixed_point_tolerance{1.0e-10};

struct GridLevel
{
    cfd::Index cells_per_direction;
};

constexpr std::array grid_levels{
    GridLevel{8}, GridLevel{16}, GridLevel{32}, GridLevel{64}, GridLevel{128},
};

struct SchemeStudy
{
    cfd::ScalarConvectionScheme scheme;
    std::string_view name;
    double minimum_final_order;
    double maximum_final_order;
};

constexpr std::array scheme_studies{
    SchemeStudy{cfd::ScalarConvectionScheme::FirstOrderUpwind, "FirstOrderUpwind", 0.5, 1.5},
    SchemeStudy{cfd::ScalarConvectionScheme::Linear, "Linear", 1.5, 2.5},
};

struct ComponentResult
{
    ErrorStatistics solution_errors;
    std::optional<double> rms_order;
    std::optional<double> linf_order;
    cfd::Index iteration_count{};
    double normalized_algebraic_residual{};
    double normalized_full_residual{};
};

struct LevelResult
{
    double characteristic_mesh_length{};
    cfd::Index cell_count{};
    ComponentResult u;
    ComponentResult v;
    double maximum_cell_mass_imbalance{};
};

struct FixedPointResult
{
    cfd::Index cell_count{};
    double relaxation_factor{};
    double maximum_u_difference{};
    double maximum_v_difference{};
};

[[nodiscard]]
double analytical_u(const cfd::Point2 &point) noexcept
{
    return 1.0 - std::exp(kovasznay_lambda * point.x) * std::cos(wave_number * point.y);
}

[[nodiscard]]
double analytical_v(const cfd::Point2 &point) noexcept
{
    return kovasznay_lambda / wave_number * std::exp(kovasznay_lambda * point.x) * std::sin(wave_number * point.y);
}

[[nodiscard]]
double analytical_pressure(const cfd::Point2 &point) noexcept
{
    return 0.5 * (1.0 - std::exp(2.0 * kovasznay_lambda * point.x));
}

[[nodiscard]]
cfd::Vector2 analytical_u_gradient(const cfd::Point2 &point) noexcept
{
    const double exponential{std::exp(kovasznay_lambda * point.x)};
    return {
        -kovasznay_lambda * exponential * std::cos(wave_number * point.y),
        wave_number * exponential * std::sin(wave_number * point.y),
    };
}

[[nodiscard]]
cfd::Vector2 analytical_v_gradient(const cfd::Point2 &point) noexcept
{
    const double exponential{std::exp(kovasznay_lambda * point.x)};
    return {
        kovasznay_lambda * kovasznay_lambda / wave_number * exponential * std::sin(wave_number * point.y),
        kovasznay_lambda * exponential * std::cos(wave_number * point.y),
    };
}

[[nodiscard]]
cfd::Vector2 analytical_pressure_gradient(const cfd::Point2 &point) noexcept
{
    return {
        -kovasznay_lambda * std::exp(2.0 * kovasznay_lambda * point.x),
        0.0,
    };
}

[[nodiscard]]
double analytical_u_laplacian(const cfd::Point2 &point) noexcept
{
    return (wave_number * wave_number - kovasznay_lambda * kovasznay_lambda) * std::exp(kovasznay_lambda * point.x) *
           std::cos(wave_number * point.y);
}

[[nodiscard]]
double analytical_v_laplacian(const cfd::Point2 &point) noexcept
{
    return kovasznay_lambda * (kovasznay_lambda * kovasznay_lambda / wave_number - wave_number) *
           std::exp(kovasznay_lambda * point.x) * std::sin(wave_number * point.y);
}

void verify_analytical_solution()
{
    constexpr std::array sample_points{
        cfd::Point2{0.15, 0.17},
        cfd::Point2{0.48, 0.39},
        cfd::Point2{0.83, 0.71},
    };
    constexpr double identity_tolerance{1.0e-12};
    constexpr double derivative_step{1.0e-6};
    constexpr double derivative_tolerance{1.0e-8};

    for (const cfd::Point2 &point : sample_points)
    {
        const double u{analytical_u(point)};
        const double v{analytical_v(point)};
        const cfd::Vector2 u_gradient{analytical_u_gradient(point)};
        const cfd::Vector2 v_gradient{analytical_v_gradient(point)};
        const cfd::Vector2 pressure_gradient{analytical_pressure_gradient(point)};
        const double continuity_residual{u_gradient.x + v_gradient.y};
        const double u_momentum_residual{density * (u * u_gradient.x + v * u_gradient.y) -
                                         dynamic_viscosity * analytical_u_laplacian(point) + pressure_gradient.x};
        const double v_momentum_residual{density * (u * v_gradient.x + v * v_gradient.y) -
                                         dynamic_viscosity * analytical_v_laplacian(point) + pressure_gradient.y};
        const cfd::Point2 point_before{point.x - derivative_step, point.y};
        const cfd::Point2 point_after{point.x + derivative_step, point.y};
        const double centered_pressure_derivative{
            (analytical_pressure(point_after) - analytical_pressure(point_before)) / (2.0 * derivative_step)};

        if (std::abs(continuity_residual) > identity_tolerance || std::abs(u_momentum_residual) > identity_tolerance ||
            std::abs(v_momentum_residual) > identity_tolerance ||
            std::abs(centered_pressure_derivative - pressure_gradient.x) > derivative_tolerance)
        {
            throw std::runtime_error("The analytical Kovasznay formulas failed their independent identity check.");
        }
    }
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
    const cfd::Index nx{level.cells_per_direction};
    const cfd::Index ny{level.cells_per_direction};
    const double delta_x{domain_length / static_cast<double>(nx)};
    const double delta_y{domain_height / static_cast<double>(ny)};
    const cfd::Index node_count{(nx + 1) * (ny + 1)};
    const cfd::Index cell_count{nx * ny};
    cfd::RawMeshData raw_mesh;
    raw_mesh.nodes.reserve(node_count);
    raw_mesh.cell_types.reserve(cell_count);
    raw_mesh.cell_nodes.reserve(4 * cell_count);
    raw_mesh.cell_node_offsets.reserve(cell_count + 1);
    raw_mesh.cell_node_offsets.push_back(0);

    for (cfd::Index j = 0; j <= ny; ++j)
    {
        for (cfd::Index i = 0; i <= nx; ++i)
        {
            raw_mesh.nodes.push_back({static_cast<double>(i) * delta_x, static_cast<double>(j) * delta_y});
        }
    }

    for (cfd::Index j = 0; j < ny; ++j)
    {
        for (cfd::Index i = 0; i < nx; ++i)
        {
            const cfd::Index lower_left{node_id(i, j, nx)};
            const cfd::Index lower_right{node_id(i + 1, j, nx)};
            const cfd::Index upper_right{node_id(i + 1, j + 1, nx)};
            const cfd::Index upper_left{node_id(i, j + 1, nx)};
            raw_mesh.cell_types.push_back(cfd::CellType::Quadrilateral);
            raw_mesh.cell_nodes.insert(raw_mesh.cell_nodes.end(), {lower_left, lower_right, upper_right, upper_left});
            raw_mesh.cell_node_offsets.push_back(raw_mesh.cell_nodes.size());
        }
    }

    raw_mesh.boundary_groups.reserve(2 * (nx + ny));
    raw_mesh.boundary_edges.reserve(2 * (nx + ny));
    for (cfd::Index j = 0; j < ny; ++j)
    {
        append_boundary_edge(raw_mesh, node_id(0, j, nx), node_id(0, j + 1, nx), "left");
        append_boundary_edge(raw_mesh, node_id(nx, j, nx), node_id(nx, j + 1, nx), "right");
    }
    for (cfd::Index i = 0; i < nx; ++i)
    {
        append_boundary_edge(raw_mesh, node_id(i, 0, nx), node_id(i + 1, 0, nx), "bottom");
        append_boundary_edge(raw_mesh, node_id(i, ny, nx), node_id(i + 1, ny, nx), "top");
    }
    return raw_mesh;
}

[[nodiscard]]
std::pair<cfd::ScalarBoundaryConditions, cfd::ScalarBoundaryConditions> make_velocity_boundary_conditions(
    const cfd::Mesh &mesh)
{
    const cfd::Index boundary_count{mesh.boundary_groups().size()};
    std::vector<cfd::ScalarBoundaryCondition> u_conditions;
    std::vector<cfd::ScalarBoundaryCondition> v_conditions;
    u_conditions.reserve(boundary_count);
    v_conditions.reserve(boundary_count);
    for (cfd::Index boundary_id = 0; boundary_id < boundary_count; ++boundary_id)
    {
        u_conditions.emplace_back(cfd::ScalarBoundaryConditionType::Dirichlet, 0.0);
        v_conditions.emplace_back(cfd::ScalarBoundaryConditionType::Dirichlet, 0.0);
    }
    std::vector<bool> assigned(boundary_count, false);

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
        const cfd::Point2 &face_center{mesh.face_centers()[face_id]};
        u_conditions[boundary_id] = {
            cfd::ScalarBoundaryConditionType::Dirichlet,
            analytical_u(face_center),
        };
        v_conditions[boundary_id] = {
            cfd::ScalarBoundaryConditionType::Dirichlet,
            analytical_v(face_center),
        };
    }

    if (std::find(assigned.begin(), assigned.end(), false) != assigned.end())
    {
        throw std::runtime_error("A Kovasznay boundary group has no face.");
    }
    return {
        cfd::ScalarBoundaryConditions{boundary_count, std::move(u_conditions)},
        cfd::ScalarBoundaryConditions{boundary_count, std::move(v_conditions)},
    };
}

[[nodiscard]]
double vertical_flux_antiderivative(const double x, const double y) noexcept
{
    return y - std::exp(kovasznay_lambda * x) * std::sin(wave_number * y) / wave_number;
}

[[nodiscard]]
double horizontal_flux_antiderivative(const double x, const double y) noexcept
{
    return std::exp(kovasznay_lambda * x) * std::sin(wave_number * y) / wave_number;
}

[[nodiscard]]
std::pair<cfd::FaceFluxField, double> make_integrated_mass_flux(const cfd::Mesh &mesh)
{
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    std::vector<double> cell_mass_imbalance(mesh.cell_count());
    const auto nodes{mesh.nodes()};
    const auto faces{mesh.faces()};
    const auto face_area_vectors{mesh.face_area_vectors()};
    const auto face_lengths{mesh.face_lengths()};
    constexpr double orientation_tolerance{256.0 * std::numeric_limits<double>::epsilon()};

    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        const cfd::Face &face{faces[face_id]};
        const cfd::Point2 &first{nodes[face.node_ids[0]]};
        const cfd::Point2 &second{nodes[face.node_ids[1]]};
        const cfd::Vector2 &area_vector{face_area_vectors[face_id]};
        const double coordinate_scale{
            std::max({1.0, std::abs(first.x), std::abs(first.y), std::abs(second.x), std::abs(second.y)})};
        const double coordinate_tolerance{orientation_tolerance * coordinate_scale};
        const double delta_x{second.x - first.x};
        const double delta_y{second.y - first.y};
        double flux{};

        if (std::abs(delta_x) <= coordinate_tolerance && std::abs(delta_y) > coordinate_tolerance)
        {
            const double y_min{std::min(first.y, second.y)};
            const double y_max{std::max(first.y, second.y)};
            const double owner_normal_sign{area_vector.x / face_lengths[face_id]};
            if (std::abs(std::abs(owner_normal_sign) - 1.0) > orientation_tolerance)
            {
                throw std::runtime_error("A vertical Kovasznay face has an invalid owner normal.");
            }
            flux = density * owner_normal_sign *
                   (vertical_flux_antiderivative(first.x, y_max) - vertical_flux_antiderivative(first.x, y_min));
        }
        else if (std::abs(delta_y) <= coordinate_tolerance && std::abs(delta_x) > coordinate_tolerance)
        {
            const double x_min{std::min(first.x, second.x)};
            const double x_max{std::max(first.x, second.x)};
            const double owner_normal_sign{area_vector.y / face_lengths[face_id]};
            if (std::abs(std::abs(owner_normal_sign) - 1.0) > orientation_tolerance)
            {
                throw std::runtime_error("A horizontal Kovasznay face has an invalid owner normal.");
            }
            flux = density * owner_normal_sign *
                   (horizontal_flux_antiderivative(x_max, first.y) - horizontal_flux_antiderivative(x_min, first.y));
        }
        else
        {
            throw std::runtime_error("The Kovasznay mass-flux integration requires axis-aligned faces.");
        }

        mass_flux[face_id] = flux;
        const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};
        cell_mass_imbalance[adjacency.owner] += flux;
        if (!adjacency.is_boundary())
        {
            cell_mass_imbalance[adjacency.neighbor] -= flux;
        }
    }

    double maximum_cell_mass_imbalance{};
    for (const double imbalance : cell_mass_imbalance)
    {
        maximum_cell_mass_imbalance = std::max(maximum_cell_mass_imbalance, std::abs(imbalance));
    }
    return {std::move(mass_flux), maximum_cell_mass_imbalance};
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
                                     const std::span<double> matrix_product)
{
    system.apply_matrix(solution.values(), matrix_product);
    double residual_squared_norm{};
    for (cfd::Index cell_id = 0; cell_id < system.cell_count(); ++cell_id)
    {
        const double residual{matrix_product[cell_id] - system.rhs()[cell_id]};
        residual_squared_norm += residual * residual;
    }
    const double rhs_norm{euclidean_norm(system.rhs())};
    if (!(rhs_norm > 0.0) || !std::isfinite(rhs_norm))
    {
        throw std::runtime_error("A Kovasznay momentum algebraic RHS norm is invalid.");
    }
    return std::sqrt(residual_squared_norm) / rhs_norm;
}

[[nodiscard]]
double normalized_full_residual(const cfd::Mesh &mesh, const cfd::ScalarConvectionOperator &convection,
                                const cfd::ScalarDiffusionOperator &diffusion,
                                const cfd::ScalarBoundaryConditions &boundary_conditions,
                                const cfd::FaceFluxField &mass_flux, const cfd::CellScalarField &solution,
                                const cfd::CellVectorField &gradient, const std::span<const double> pressure_source,
                                cfd::CellScalarField &convection_balance, cfd::CellScalarField &diffusion_balance)
{
    convection.compute_flux_balance(solution, boundary_conditions, mass_flux, convection_balance);
    diffusion.compute_flux_balance(solution, boundary_conditions, gradient, diffusion_balance);

    double convection_squared_norm{};
    double diffusion_squared_norm{};
    double pressure_source_squared_norm{};
    double residual_squared_norm{};
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        const double convection_value{convection_balance[cell_id]};
        const double diffusion_value{diffusion_balance[cell_id]};
        const double source_value{pressure_source[cell_id]};
        const double residual{convection_value + diffusion_value - source_value};
        convection_squared_norm += convection_value * convection_value;
        diffusion_squared_norm += diffusion_value * diffusion_value;
        pressure_source_squared_norm += source_value * source_value;
        residual_squared_norm += residual * residual;
    }
    const double normalization{std::sqrt(convection_squared_norm) + std::sqrt(diffusion_squared_norm) +
                               std::sqrt(pressure_source_squared_norm)};
    // This remains a useful scale for v momentum, whose pressure source is
    // identically zero in the Kovasznay solution.
    if (!(normalization > 0.0) || !std::isfinite(normalization))
    {
        throw std::runtime_error("A Kovasznay full momentum residual normalization is invalid.");
    }
    return std::sqrt(residual_squared_norm) / normalization;
}

[[nodiscard]]
LevelResult run_level(const GridLevel &level, const SchemeStudy &study)
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_raw_mesh(level))};
    const cfd::Mesh &mesh{build_result.mesh};
    auto [u_boundary_conditions, v_boundary_conditions]{make_velocity_boundary_conditions(mesh)};
    auto [mass_flux, maximum_cell_mass_imbalance]{make_integrated_mass_flux(mesh)};

    cfd::CellVelocityField previous_velocity{mesh.cell_count()};
    cfd::CellVectorField u_gradient{mesh.cell_count()};
    cfd::CellVectorField v_gradient{mesh.cell_count()};
    cfd::CellVectorField pressure_gradient{mesh.cell_count()};
    std::vector<double> u_pressure_source(mesh.cell_count());
    std::vector<double> v_pressure_source(mesh.cell_count());
    double total_area{};
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        const cfd::Point2 &cell_center{mesh.cell_centers()[cell_id]};
        u_gradient[cell_id] = analytical_u_gradient(cell_center);
        v_gradient[cell_id] = analytical_v_gradient(cell_center);
        pressure_gradient[cell_id] = analytical_pressure_gradient(cell_center);
        const double cell_area{mesh.cell_areas()[cell_id]};
        u_pressure_source[cell_id] = -cell_area * pressure_gradient[cell_id].x;
        v_pressure_source[cell_id] = -cell_area * pressure_gradient[cell_id].y;
        total_area += cell_area;
    }

    const cfd::IncompressibleMomentumAssembler assembler{mesh, dynamic_viscosity, study.scheme};
    cfd::ScalarLinearSystem u_system{mesh};
    cfd::ScalarLinearSystem v_system{mesh};
    assembler.assemble(previous_velocity, u_gradient, v_gradient, pressure_gradient, u_boundary_conditions,
                       v_boundary_conditions, mass_flux, 1.0, u_system, v_system);

    cfd::EigenBiCGSTABSolver u_solver{{1.0e-13, 10000}};
    cfd::EigenBiCGSTABSolver v_solver{{1.0e-13, 10000}};
    u_solver.compute_matrix(u_system);
    v_solver.compute_matrix(v_system);
    cfd::CellScalarField numerical_u{mesh.cell_count()};
    cfd::CellScalarField numerical_v{mesh.cell_count()};
    const cfd::LinearSolveResult u_solve_result{u_solver.solve(u_system.rhs(), numerical_u.values())};
    const cfd::LinearSolveResult v_solve_result{v_solver.solve(v_system.rhs(), numerical_v.values())};
    if (!u_solve_result.converged || !v_solve_result.converged)
    {
        throw std::runtime_error(std::string{study.name} + " Kovasznay BiCGSTAB solve did not converge.");
    }

    std::vector<double> u_matrix_product(mesh.cell_count());
    std::vector<double> v_matrix_product(mesh.cell_count());
    const double u_algebraic_residual{normalized_algebraic_residual(u_system, numerical_u, u_matrix_product)};
    const double v_algebraic_residual{normalized_algebraic_residual(v_system, numerical_v, v_matrix_product)};

    cfd::CellVectorField numerical_u_gradient{mesh.cell_count()};
    cfd::CellVectorField numerical_v_gradient{mesh.cell_count()};
    cfd::compute_least_squares_gradient(mesh, numerical_u, u_boundary_conditions, numerical_u_gradient);
    cfd::compute_least_squares_gradient(mesh, numerical_v, v_boundary_conditions, numerical_v_gradient);
    const cfd::ScalarConvectionOperator convection{mesh, study.scheme};
    const cfd::ScalarDiffusionOperator diffusion{mesh, dynamic_viscosity};
    cfd::CellScalarField u_convection_balance{mesh.cell_count()};
    cfd::CellScalarField u_diffusion_balance{mesh.cell_count()};
    cfd::CellScalarField v_convection_balance{mesh.cell_count()};
    cfd::CellScalarField v_diffusion_balance{mesh.cell_count()};
    const double u_full_residual{normalized_full_residual(mesh, convection, diffusion, u_boundary_conditions, mass_flux,
                                                          numerical_u, numerical_u_gradient, u_pressure_source,
                                                          u_convection_balance, u_diffusion_balance)};
    const double v_full_residual{normalized_full_residual(mesh, convection, diffusion, v_boundary_conditions, mass_flux,
                                                          numerical_v, numerical_v_gradient, v_pressure_source,
                                                          v_convection_balance, v_diffusion_balance)};

    if (u_algebraic_residual > algebraic_residual_tolerance || v_algebraic_residual > algebraic_residual_tolerance ||
        u_full_residual > full_residual_tolerance || v_full_residual > full_residual_tolerance ||
        maximum_cell_mass_imbalance > mass_imbalance_tolerance)
    {
        std::ostringstream message;
        message << study.name << " Kovasznay diagnostic exceeded its tolerance at N=" << level.cells_per_direction
                << "; u_alg=" << u_algebraic_residual << ", v_alg=" << v_algebraic_residual
                << ", u_fv=" << u_full_residual << ", v_fv=" << v_full_residual
                << ", max_divF=" << maximum_cell_mass_imbalance << '.';
        throw std::runtime_error(message.str());
    }

    ErrorAccumulator u_error_accumulator;
    ErrorAccumulator v_error_accumulator;
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        const cfd::Point2 &cell_center{mesh.cell_centers()[cell_id]};
        const double cell_area{mesh.cell_areas()[cell_id]};
        u_error_accumulator.add(cell_area, std::abs(numerical_u[cell_id] - analytical_u(cell_center)));
        v_error_accumulator.add(cell_area, std::abs(numerical_v[cell_id] - analytical_v(cell_center)));
    }

    return {
        .characteristic_mesh_length = std::sqrt(total_area / static_cast<double>(mesh.cell_count())),
        .cell_count = mesh.cell_count(),
        .u =
            {
                .solution_errors = u_error_accumulator.finish(),
                .rms_order = std::nullopt,
                .linf_order = std::nullopt,
                .iteration_count = u_solve_result.iteration_count,
                .normalized_algebraic_residual = u_algebraic_residual,
                .normalized_full_residual = u_full_residual,
            },
        .v =
            {
                .solution_errors = v_error_accumulator.finish(),
                .rms_order = std::nullopt,
                .linf_order = std::nullopt,
                .iteration_count = v_solve_result.iteration_count,
                .normalized_algebraic_residual = v_algebraic_residual,
                .normalized_full_residual = v_full_residual,
            },
        .maximum_cell_mass_imbalance = maximum_cell_mass_imbalance,
    };
}

[[nodiscard]]
FixedPointResult run_relaxation_fixed_point_check()
{
    constexpr GridLevel representative_level{32};
    constexpr double relaxation_factor{0.6};
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_raw_mesh(representative_level))};
    const cfd::Mesh &mesh{build_result.mesh};
    auto [u_boundary_conditions, v_boundary_conditions]{make_velocity_boundary_conditions(mesh)};
    auto [mass_flux, maximum_cell_mass_imbalance]{make_integrated_mass_flux(mesh)};
    if (maximum_cell_mass_imbalance > mass_imbalance_tolerance)
    {
        throw std::runtime_error("The relaxation fixture mass flux is not discretely conservative.");
    }

    cfd::CellVectorField u_gradient{mesh.cell_count()};
    cfd::CellVectorField v_gradient{mesh.cell_count()};
    cfd::CellVectorField pressure_gradient{mesh.cell_count()};
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        const cfd::Point2 &cell_center{mesh.cell_centers()[cell_id]};
        u_gradient[cell_id] = analytical_u_gradient(cell_center);
        v_gradient[cell_id] = analytical_v_gradient(cell_center);
        pressure_gradient[cell_id] = analytical_pressure_gradient(cell_center);
    }

    const cfd::IncompressibleMomentumAssembler assembler{mesh, dynamic_viscosity, cfd::ScalarConvectionScheme::Linear};
    const cfd::CellVelocityField initial_previous_velocity{mesh.cell_count()};
    cfd::ScalarLinearSystem unrelaxed_u_system{mesh};
    cfd::ScalarLinearSystem unrelaxed_v_system{mesh};
    assembler.assemble(initial_previous_velocity, u_gradient, v_gradient, pressure_gradient, u_boundary_conditions,
                       v_boundary_conditions, mass_flux, 1.0, unrelaxed_u_system, unrelaxed_v_system);

    cfd::EigenBiCGSTABSolver unrelaxed_u_solver{{1.0e-13, 10000}};
    cfd::EigenBiCGSTABSolver unrelaxed_v_solver{{1.0e-13, 10000}};
    unrelaxed_u_solver.compute_matrix(unrelaxed_u_system);
    unrelaxed_v_solver.compute_matrix(unrelaxed_v_system);
    cfd::CellScalarField unrelaxed_u{mesh.cell_count()};
    cfd::CellScalarField unrelaxed_v{mesh.cell_count()};
    const cfd::LinearSolveResult unrelaxed_u_result{
        unrelaxed_u_solver.solve(unrelaxed_u_system.rhs(), unrelaxed_u.values())};
    const cfd::LinearSolveResult unrelaxed_v_result{
        unrelaxed_v_solver.solve(unrelaxed_v_system.rhs(), unrelaxed_v.values())};
    if (!unrelaxed_u_result.converged || !unrelaxed_v_result.converged)
    {
        throw std::runtime_error("The unrelaxed Kovasznay fixed-point solve did not converge.");
    }

    cfd::CellVelocityField discrete_previous_velocity{mesh.cell_count()};
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        discrete_previous_velocity.u()[cell_id] = unrelaxed_u[cell_id];
        discrete_previous_velocity.v()[cell_id] = unrelaxed_v[cell_id];
    }
    cfd::ScalarLinearSystem relaxed_u_system{mesh};
    cfd::ScalarLinearSystem relaxed_v_system{mesh};
    assembler.assemble(discrete_previous_velocity, u_gradient, v_gradient, pressure_gradient, u_boundary_conditions,
                       v_boundary_conditions, mass_flux, relaxation_factor, relaxed_u_system, relaxed_v_system);

    cfd::EigenBiCGSTABSolver relaxed_u_solver{{1.0e-13, 10000}};
    cfd::EigenBiCGSTABSolver relaxed_v_solver{{1.0e-13, 10000}};
    relaxed_u_solver.compute_matrix(relaxed_u_system);
    relaxed_v_solver.compute_matrix(relaxed_v_system);
    cfd::CellScalarField relaxed_u{mesh.cell_count()};
    cfd::CellScalarField relaxed_v{mesh.cell_count()};
    const cfd::LinearSolveResult relaxed_u_result{relaxed_u_solver.solve(relaxed_u_system.rhs(), relaxed_u.values())};
    const cfd::LinearSolveResult relaxed_v_result{relaxed_v_solver.solve(relaxed_v_system.rhs(), relaxed_v.values())};
    if (!relaxed_u_result.converged || !relaxed_v_result.converged)
    {
        throw std::runtime_error("The relaxed Kovasznay fixed-point solve did not converge.");
    }

    double maximum_u_difference{};
    double maximum_v_difference{};
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        maximum_u_difference = std::max(maximum_u_difference, std::abs(relaxed_u[cell_id] - unrelaxed_u[cell_id]));
        maximum_v_difference = std::max(maximum_v_difference, std::abs(relaxed_v[cell_id] - unrelaxed_v[cell_id]));
    }
    if (maximum_u_difference > fixed_point_tolerance || maximum_v_difference > fixed_point_tolerance)
    {
        throw std::runtime_error("Equation under-relaxation changed the discrete Kovasznay fixed point.");
    }

    return {
        .cell_count = mesh.cell_count(),
        .relaxation_factor = relaxation_factor,
        .maximum_u_difference = maximum_u_difference,
        .maximum_v_difference = maximum_v_difference,
    };
}

void compute_orders(std::vector<LevelResult> &results)
{
    for (std::size_t level_index = 1; level_index < results.size(); ++level_index)
    {
        const LevelResult &coarse{results[level_index - 1]};
        LevelResult &fine{results[level_index]};
        fine.u.rms_order = observed_order(coarse.u.solution_errors.area_weighted_rms_error,
                                          fine.u.solution_errors.area_weighted_rms_error,
                                          coarse.characteristic_mesh_length, fine.characteristic_mesh_length);
        fine.u.linf_order = observed_order(coarse.u.solution_errors.linf_error, fine.u.solution_errors.linf_error,
                                           coarse.characteristic_mesh_length, fine.characteristic_mesh_length);
        fine.v.rms_order = observed_order(coarse.v.solution_errors.area_weighted_rms_error,
                                          fine.v.solution_errors.area_weighted_rms_error,
                                          coarse.characteristic_mesh_length, fine.characteristic_mesh_length);
        fine.v.linf_order = observed_order(coarse.v.solution_errors.linf_error, fine.v.solution_errors.linf_error,
                                           coarse.characteristic_mesh_length, fine.characteristic_mesh_length);
    }
}

void require_expected_convergence(const SchemeStudy &study, const std::vector<LevelResult> &results)
{
    for (std::size_t level_index = 1; level_index < results.size(); ++level_index)
    {
        const LevelResult &coarse{results[level_index - 1]};
        const LevelResult &fine{results[level_index]};
        if (!(fine.u.solution_errors.area_weighted_rms_error < coarse.u.solution_errors.area_weighted_rms_error) ||
            !(fine.u.solution_errors.linf_error < coarse.u.solution_errors.linf_error) ||
            !(fine.v.solution_errors.area_weighted_rms_error < coarse.v.solution_errors.area_weighted_rms_error) ||
            !(fine.v.solution_errors.linf_error < coarse.v.solution_errors.linf_error))
        {
            throw std::runtime_error(std::string{study.name} +
                                     " Kovasznay velocity errors did not decrease monotonically.");
        }
    }

    const LevelResult &finest{results.back()};
    const auto order_is_expected = [&study](const std::optional<double> order) {
        return order.has_value() && *order > study.minimum_final_order && *order < study.maximum_final_order;
    };
    if (!order_is_expected(finest.u.rms_order) || !order_is_expected(finest.u.linf_order) ||
        !order_is_expected(finest.v.rms_order) || !order_is_expected(finest.v.linf_order))
    {
        throw std::runtime_error(std::string{study.name} +
                                 " Kovasznay asymptotic velocity order is outside its expected range.");
    }
}

void print_order(const std::optional<double> order, const bool first_level)
{
    if (order.has_value())
    {
        std::cout << std::fixed << std::setprecision(3) << std::setw(8) << *order;
        return;
    }
    std::cout << std::setw(8) << (first_level ? "-" : "n/a");
}

void print_results(const SchemeStudy &study, const std::vector<LevelResult> &results)
{
    std::cout << "\nIncompressible momentum solution convergence - Kovasznay - " << study.name << " - Cartesian QUAD\n"
              << "rho=1, Re=40, mu=1/Re, exact integrated face mass flux, exact velocity Dirichlet data\n\n"
              << std::left << std::setw(11) << "h" << std::setw(9) << "cells" << std::setw(15) << "RMS(u)"
              << std::setw(8) << "p" << std::setw(15) << "Linf(u)" << std::setw(8) << "p" << std::setw(15) << "RMS(v)"
              << std::setw(8) << "p" << std::setw(15) << "Linf(v)" << std::setw(8) << "p" << std::setw(8) << "u_it"
              << std::setw(8) << "v_it" << std::setw(13) << "u_alg" << std::setw(13) << "v_alg" << std::setw(13)
              << "u_FV" << std::setw(13) << "v_FV" << std::setw(13) << "max_divF" << '\n'
              << std::string(203, '-') << '\n';

    for (std::size_t level_index = 0; level_index < results.size(); ++level_index)
    {
        const LevelResult &result{results[level_index]};
        std::cout << std::scientific << std::setprecision(3) << std::setw(11) << result.characteristic_mesh_length
                  << std::fixed << std::setw(9) << result.cell_count << std::scientific << std::setprecision(6)
                  << std::setw(15) << result.u.solution_errors.area_weighted_rms_error;
        print_order(result.u.rms_order, level_index == 0);
        std::cout << std::scientific << std::setprecision(6) << std::setw(15) << result.u.solution_errors.linf_error;
        print_order(result.u.linf_order, level_index == 0);
        std::cout << std::scientific << std::setprecision(6) << std::setw(15)
                  << result.v.solution_errors.area_weighted_rms_error;
        print_order(result.v.rms_order, level_index == 0);
        std::cout << std::scientific << std::setprecision(6) << std::setw(15) << result.v.solution_errors.linf_error;
        print_order(result.v.linf_order, level_index == 0);
        std::cout << std::fixed << std::setw(8) << result.u.iteration_count << std::setw(8) << result.v.iteration_count
                  << std::scientific << std::setprecision(3) << std::setw(13) << result.u.normalized_algebraic_residual
                  << std::setw(13) << result.v.normalized_algebraic_residual << std::setw(13)
                  << result.u.normalized_full_residual << std::setw(13) << result.v.normalized_full_residual
                  << std::setw(13) << result.maximum_cell_mass_imbalance << '\n';
    }

    std::cout << "\nInterpretation: on this Cartesian Kovasznay study, " << study.name
              << " reaches final observed orders within [" << study.minimum_final_order << ", "
              << study.maximum_final_order << "] for RMS and Linf errors of both velocity components" << ".\n";
}

} // namespace

int main()
{
    try
    {
        verify_analytical_solution();
        for (const SchemeStudy &study : scheme_studies)
        {
            std::vector<LevelResult> results;
            results.reserve(grid_levels.size());
            for (const GridLevel &level : grid_levels)
            {
                results.push_back(run_level(level, study));
            }
            compute_orders(results);
            require_expected_convergence(study, results);
            print_results(study, results);
        }
        const FixedPointResult fixed_point{run_relaxation_fixed_point_check()};
        std::cout << "\nEquation under-relaxation fixed-point check - Linear - " << fixed_point.cell_count
                  << " cells - alpha=" << fixed_point.relaxation_factor << "\n"
                  << "max|u_alpha-u_1|=" << std::scientific << fixed_point.maximum_u_difference
                  << ", max|v_alpha-v_1|=" << fixed_point.maximum_v_difference << '\n';
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "Incompressible momentum solution verification failed: " << error.what() << '\n';
        return 1;
    }
}
