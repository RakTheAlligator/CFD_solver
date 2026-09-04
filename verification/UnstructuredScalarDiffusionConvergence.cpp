#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVectorField.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/io/VtkWriter.hpp"
#include "cfd/linear_algebra/EigenConjugateGradientSolver.hpp"
#include "cfd/linear_algebra/ScalarLinearSystem.hpp"
#include "cfd/math/Point2.hpp"
#include "cfd/math/Vector2.hpp"
#include "cfd/mesh/Boundary.hpp"
#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/mesh/Types.hpp"
#include "cfd/meshing/GmshMesher.hpp"
#include "cfd/meshing/RawMeshData.hpp"
#include "cfd/meshing/RectangleGeometry.hpp"
#include "cfd/numerics/LeastSquaresGradient.hpp"
#include "cfd/numerics/ScalarDiffusionOperator.hpp"

#include "support/VerificationStatistics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
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

constexpr double domain_length{2.0};
constexpr double domain_height{1.0};
constexpr double diffusivity{1.0};
constexpr double grading_strength{0.5};
constexpr double full_residual_tolerance{1.0e-10};
constexpr cfd::Index maximum_correction_count{200};
constexpr double radians_to_degrees{180.0 / std::numbers::pi};

struct GridLevel
{
    double target_size;
    cfd::Index nx;
    cfd::Index ny;
};

constexpr std::array grid_levels{
    GridLevel{0.20, 10, 5},   GridLevel{0.10, 20, 10},    GridLevel{0.05, 40, 20},
    GridLevel{0.025, 80, 40}, GridLevel{0.0125, 160, 80},
};

enum class MeshFamilyKind : std::uint8_t
{
    GmshTriangle,
    GmshQuadrilateral,
    GradedOrthogonalQuadrilateral
};

struct MeshFamily
{
    MeshFamilyKind kind;
    std::string_view name;
    std::string_view file_prefix;
};

constexpr std::array mesh_families{
    MeshFamily{MeshFamilyKind::GmshTriangle, "unstructured Gmsh TRI", "gmsh_tri"},
    MeshFamily{MeshFamilyKind::GmshQuadrilateral, "recombined Gmsh QUAD", "gmsh_quad"},
    MeshFamily{MeshFamilyKind::GradedOrthogonalQuadrilateral, "graded orthogonal QUAD", "graded_quad"},
};

struct GeometryDiagnostics
{
    double minimum_normalized_projection{std::numeric_limits<double>::infinity()};
    double mean_normalized_projection{};
    double maximum_normalized_projection{-std::numeric_limits<double>::infinity()};
    double maximum_non_orthogonality_degrees{};
    double minimum_lambda{std::numeric_limits<double>::infinity()};
    double maximum_lambda{-std::numeric_limits<double>::infinity()};
    double maximum_correction_ratio{};
    double minimum_cell_area{std::numeric_limits<double>::infinity()};
    double maximum_cell_area{};
    double cell_area_ratio{};
};

struct GradingDiagnostics
{
    double minimum_cell_width{};
    double maximum_cell_width{};
    double width_ratio{};
};

struct ExactGradientDiagnostic
{
    ErrorStatistics solution_errors;
    std::optional<double> rms_order;
    std::optional<double> linf_order;
    cfd::Index correction_count{};
    cfd::Index total_cg_iterations{};
    double normalized_full_residual{};
    double normalized_algebraic_residual{};
};

struct LevelResult
{
    double target_size{};
    double characteristic_mesh_length{};
    cfd::Index cell_count{};
    cfd::Index internal_face_count{};
    cfd::Index matrix_nonzero_count{};
    ErrorStatistics solution_errors;
    std::optional<double> rms_order;
    std::optional<double> linf_order;
    cfd::Index correction_count{};
    cfd::Index total_cg_iterations{};
    double normalized_full_residual{};
    double normalized_algebraic_residual{};
    GeometryDiagnostics geometry;
    std::optional<GradingDiagnostics> grading;
    std::optional<ExactGradientDiagnostic> exact_gradient_diagnostic;
};

[[nodiscard]]
double analytical_phi(const cfd::Point2 &point) noexcept
{
    return std::sin(point.x + 2.0 * point.y) + 0.2 * point.x * point.x + 0.3 * point.x * point.y +
           0.1 * point.y * point.y + 1.0;
}

[[nodiscard]]
cfd::Vector2 analytical_gradient(const cfd::Point2 &point) noexcept
{
    const double cosine{std::cos(point.x + 2.0 * point.y)};
    return {
        cosine + 0.4 * point.x + 0.3 * point.y,
        2.0 * cosine + 0.3 * point.x + 0.2 * point.y,
    };
}

[[nodiscard]]
double source_density(const cfd::Point2 &point) noexcept
{
    return 5.0 * std::sin(point.x + 2.0 * point.y) - 0.6;
}

[[nodiscard]]
double dot(const cfd::Vector2 &first, const cfd::Vector2 &second) noexcept
{
    return first.x * second.x + first.y * second.y;
}

[[nodiscard]]
cfd::Index node_id(const cfd::Index i, const cfd::Index j, const cfd::Index nx) noexcept
{
    return j * (nx + 1) + i;
}

[[nodiscard]]
double graded_x(const double xi) noexcept
{
    const double normalized_coordinate{xi / domain_length};
    return domain_length * ((1.0 - grading_strength) * normalized_coordinate +
                            grading_strength * normalized_coordinate * normalized_coordinate);
}

void append_boundary_edge(cfd::RawMeshData &raw_mesh, const cfd::Index first_node, const cfd::Index second_node,
                          const std::string_view side)
{
    const cfd::BoundaryId boundary_id{raw_mesh.boundary_groups.size()};
    raw_mesh.boundary_groups.push_back({boundary_id, std::string{side} + "_" + std::to_string(boundary_id)});
    raw_mesh.boundary_edges.push_back({{first_node, second_node}, boundary_id});
}

[[nodiscard]]
cfd::RawMeshData make_graded_raw_mesh(const GridLevel &level)
{
    cfd::RawMeshData raw_mesh;
    raw_mesh.nodes.reserve((level.nx + 1) * (level.ny + 1));
    raw_mesh.cell_types.reserve(level.nx * level.ny);
    raw_mesh.cell_nodes.reserve(4 * level.nx * level.ny);
    raw_mesh.cell_node_offsets.reserve(level.nx * level.ny + 1);
    raw_mesh.cell_node_offsets.push_back(0);

    for (cfd::Index j = 0; j <= level.ny; ++j)
    {
        for (cfd::Index i = 0; i <= level.nx; ++i)
        {
            raw_mesh.nodes.push_back({
                graded_x(static_cast<double>(i) * level.target_size),
                static_cast<double>(j) * level.target_size,
            });
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

void split_gmsh_boundary_groups_by_edge(cfd::RawMeshData &raw_mesh)
{
    std::vector<cfd::BoundaryGroup> boundary_groups;
    boundary_groups.reserve(raw_mesh.boundary_edges.size());
    for (cfd::Index edge_id = 0; edge_id < raw_mesh.boundary_edges.size(); ++edge_id)
    {
        cfd::BoundaryEdge &edge{raw_mesh.boundary_edges[edge_id]};
        const cfd::BoundaryId original_boundary_id{edge.boundary_id};
        const cfd::BoundaryId boundary_id{boundary_groups.size()};
        boundary_groups.push_back(
            {boundary_id, raw_mesh.boundary_groups[original_boundary_id].name + "_" + std::to_string(edge_id)});
        edge.boundary_id = boundary_id;
    }
    raw_mesh.boundary_groups = std::move(boundary_groups);
}

[[nodiscard]]
cfd::RawMeshData make_raw_mesh(const MeshFamily &family, const GridLevel &level)
{
    if (family.kind == MeshFamilyKind::GradedOrthogonalQuadrilateral)
    {
        return make_graded_raw_mesh(level);
    }

    const cfd::CellType cell_type{family.kind == MeshFamilyKind::GmshTriangle ? cfd::CellType::Triangle
                                                                              : cfd::CellType::Quadrilateral};
    cfd::RawMeshData raw_mesh{cfd::generate_mesh({.length = domain_length, .height = domain_height},
                                                 {.mesh_size = level.target_size, .cell_type = cell_type})};
    // Gmsh's three physical groups are split locally because the manufactured
    // Dirichlet value is face-dependent and ScalarBoundaryConditions is constant per group.
    split_gmsh_boundary_groups_by_edge(raw_mesh);
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
            throw std::runtime_error("A solution-verification BoundaryId is shared by multiple faces.");
        }
        assigned[boundary_id] = true;
        conditions[boundary_id] = {
            cfd::ScalarBoundaryConditionType::Dirichlet,
            analytical_phi(mesh.face_centers()[face_id]),
        };
    }

    if (std::find(assigned.begin(), assigned.end(), false) != assigned.end())
    {
        throw std::runtime_error("A solution-verification boundary group has no face.");
    }
    return {mesh.boundary_groups().size(), std::move(conditions)};
}

[[nodiscard]]
GeometryDiagnostics compute_geometry_diagnostics(const cfd::Mesh &mesh)
{
    GeometryDiagnostics diagnostics;
    double normalized_projection_sum{};
    cfd::Index internal_face_count{};

    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};
        if (adjacency.is_boundary())
        {
            continue;
        }
        ++internal_face_count;

        const cfd::Point2 &owner_center{mesh.cell_centers()[adjacency.owner]};
        const cfd::Point2 &neighbor_center{mesh.cell_centers()[adjacency.neighbor]};
        const cfd::Point2 &face_center{mesh.face_centers()[face_id]};
        const cfd::Vector2 &area_vector{mesh.face_area_vectors()[face_id]};
        const cfd::Vector2 displacement{
            neighbor_center.x - owner_center.x,
            neighbor_center.y - owner_center.y,
        };
        const cfd::Vector2 owner_to_face{
            face_center.x - owner_center.x,
            face_center.y - owner_center.y,
        };
        const double area_dot_displacement{dot(area_vector, displacement)};
        const double area_norm{mesh.face_lengths()[face_id]};
        const double displacement_norm{std::hypot(displacement.x, displacement.y)};
        const double normalized_projection{area_dot_displacement / (area_norm * displacement_norm)};
        const double lambda{dot(area_vector, owner_to_face) / area_dot_displacement};
        const double beta{dot(area_vector, area_vector) / area_dot_displacement};
        const cfd::Vector2 correction{
            area_vector.x - beta * displacement.x,
            area_vector.y - beta * displacement.y,
        };
        const double correction_ratio{std::hypot(correction.x, correction.y) / area_norm};
        const double bounded_projection{std::clamp(normalized_projection, -1.0, 1.0)};

        diagnostics.minimum_normalized_projection =
            std::min(diagnostics.minimum_normalized_projection, normalized_projection);
        diagnostics.maximum_normalized_projection =
            std::max(diagnostics.maximum_normalized_projection, normalized_projection);
        normalized_projection_sum += normalized_projection;
        diagnostics.maximum_non_orthogonality_degrees =
            std::max(diagnostics.maximum_non_orthogonality_degrees, std::acos(bounded_projection) * radians_to_degrees);
        diagnostics.minimum_lambda = std::min(diagnostics.minimum_lambda, lambda);
        diagnostics.maximum_lambda = std::max(diagnostics.maximum_lambda, lambda);
        diagnostics.maximum_correction_ratio = std::max(diagnostics.maximum_correction_ratio, correction_ratio);
    }

    for (const double area : mesh.cell_areas())
    {
        diagnostics.minimum_cell_area = std::min(diagnostics.minimum_cell_area, area);
        diagnostics.maximum_cell_area = std::max(diagnostics.maximum_cell_area, area);
    }
    if (internal_face_count == 0 || !(diagnostics.minimum_cell_area > 0.0))
    {
        throw std::runtime_error("Solution verification cannot compute geometry diagnostics.");
    }
    diagnostics.mean_normalized_projection = normalized_projection_sum / static_cast<double>(internal_face_count);
    diagnostics.cell_area_ratio = diagnostics.maximum_cell_area / diagnostics.minimum_cell_area;

    if (!std::isfinite(diagnostics.minimum_normalized_projection) ||
        !std::isfinite(diagnostics.mean_normalized_projection) ||
        !std::isfinite(diagnostics.maximum_normalized_projection) ||
        !std::isfinite(diagnostics.maximum_non_orthogonality_degrees) || !std::isfinite(diagnostics.minimum_lambda) ||
        !std::isfinite(diagnostics.maximum_lambda) || !std::isfinite(diagnostics.maximum_correction_ratio) ||
        !std::isfinite(diagnostics.maximum_cell_area) || !std::isfinite(diagnostics.cell_area_ratio))
    {
        throw std::runtime_error("Solution verification produced non-finite geometry diagnostics.");
    }
    return diagnostics;
}

[[nodiscard]]
GradingDiagnostics compute_grading_diagnostics(const GridLevel &level)
{
    double minimum_width{std::numeric_limits<double>::infinity()};
    double maximum_width{};
    for (cfd::Index i = 0; i < level.nx; ++i)
    {
        const double left{graded_x(static_cast<double>(i) * level.target_size)};
        const double right{graded_x(static_cast<double>(i + 1) * level.target_size)};
        const double width{right - left};
        minimum_width = std::min(minimum_width, width);
        maximum_width = std::max(maximum_width, width);
    }
    return {
        .minimum_cell_width = minimum_width,
        .maximum_cell_width = maximum_width,
        .width_ratio = maximum_width / minimum_width,
    };
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
        const double residual{rhs[cell_id] - matrix_product[cell_id]};
        squared_residual_norm += residual * residual;
    }
    const double rhs_norm{euclidean_norm(rhs)};
    if (!(rhs_norm > 0.0))
    {
        throw std::runtime_error("Solution verification algebraic RHS has zero norm.");
    }
    return std::sqrt(squared_residual_norm) / rhs_norm;
}

[[nodiscard]]
double normalized_full_residual(const cfd::Mesh &mesh, const cfd::ScalarDiffusionOperator &diffusion,
                                const cfd::ScalarBoundaryConditions &boundary_conditions,
                                const cfd::CellScalarField &solution, const cfd::CellVectorField &gradient,
                                const std::span<const double> source_integrals, const double source_norm,
                                cfd::CellScalarField &flux_balance)
{
    diffusion.compute_flux_balance(solution, boundary_conditions, gradient, flux_balance);
    double squared_residual_norm{};
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        const double residual{source_integrals[cell_id] - flux_balance[cell_id]};
        squared_residual_norm += residual * residual;
    }
    return std::sqrt(squared_residual_norm) / source_norm;
}

[[nodiscard]]
ExactGradientDiagnostic run_exact_gradient_diagnostic(const cfd::Mesh &mesh,
                                                      const cfd::ScalarDiffusionOperator &diffusion,
                                                      const cfd::ScalarBoundaryConditions &boundary_conditions,
                                                      const std::span<const double> source_integrals,
                                                      const double source_norm)
{
    cfd::ScalarLinearSystem system{mesh};
    diffusion.add_matrix_contributions(boundary_conditions, system);
    diffusion.add_boundary_rhs(boundary_conditions, system.rhs());
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        system.rhs()[cell_id] += source_integrals[cell_id];
    }

    cfd::CellVectorField exact_gradient{mesh.cell_count()};
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        exact_gradient[cell_id] = analytical_gradient(mesh.cell_centers()[cell_id]);
    }
    diffusion.add_non_orthogonal_rhs(boundary_conditions, exact_gradient, system.rhs());

    cfd::EigenConjugateGradientSolver solver{{1.0e-14, 10000}};
    solver.compute_matrix(system);
    cfd::CellScalarField numerical_solution{mesh.cell_count()};
    const cfd::LinearSolveResult solve_result{solver.solve(system.rhs(), numerical_solution.values())};
    if (!solve_result.converged)
    {
        throw std::runtime_error("Exact-gradient diagnostic conjugate-gradient solve did not converge.");
    }

    std::vector<double> matrix_product(mesh.cell_count());
    const double algebraic_residual{
        normalized_algebraic_residual(system, numerical_solution, system.rhs(), matrix_product)};
    cfd::CellScalarField flux_balance{mesh.cell_count()};
    const double full_residual{normalized_full_residual(mesh, diffusion, boundary_conditions, numerical_solution,
                                                        exact_gradient, source_integrals, source_norm, flux_balance)};
    if (full_residual > full_residual_tolerance)
    {
        throw std::runtime_error("Exact-gradient diagnostic did not reach the full discrete residual tolerance.");
    }

    ErrorAccumulator error_accumulator;
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        const double error{std::abs(numerical_solution[cell_id] - analytical_phi(mesh.cell_centers()[cell_id]))};
        error_accumulator.add(mesh.cell_areas()[cell_id], error);
    }
    return {
        .solution_errors = error_accumulator.finish(),
        .rms_order = std::nullopt,
        .linf_order = std::nullopt,
        .correction_count = 1,
        .total_cg_iterations = solve_result.iteration_count,
        .normalized_full_residual = full_residual,
        .normalized_algebraic_residual = algebraic_residual,
    };
}

void write_verification_vtu(const cfd::Mesh &mesh, const cfd::CellScalarField &exact_solution,
                            const cfd::CellScalarField &numerical_solution, const cfd::CellScalarField &solution_error,
                            const cfd::CellScalarField &absolute_solution_error,
                            const cfd::CellVectorField &numerical_gradient, const cfd::CellVectorField &exact_gradient,
                            const cfd::CellVectorField &gradient_error,
                            const cfd::CellScalarField &full_residual_density, const std::filesystem::path &output_path)
{
    const std::array scalar_fields{
        cfd::VtkCellScalarData{"phi_exact", exact_solution.values()},
        cfd::VtkCellScalarData{"phi_numerical", numerical_solution.values()},
        cfd::VtkCellScalarData{"phi_error", solution_error.values()},
        cfd::VtkCellScalarData{"abs_phi_error", absolute_solution_error.values()},
        cfd::VtkCellScalarData{"full_residual_density", full_residual_density.values()},
    };
    const std::array vector_fields{
        cfd::VtkCellVectorData{"gradient_numerical", numerical_gradient.values()},
        cfd::VtkCellVectorData{"gradient_exact", exact_gradient.values()},
        cfd::VtkCellVectorData{"gradient_error", gradient_error.values()},
    };
    cfd::write_vtu(mesh, output_path, cfd::VtkCellData{.scalars = scalar_fields, .vectors = vector_fields});
}

[[nodiscard]]
LevelResult run_level(const MeshFamily &family, const GridLevel &level, const std::size_t level_index,
                      const bool write_vtu, const std::filesystem::path &output_directory,
                      std::vector<std::filesystem::path> &written_paths)
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_raw_mesh(family, level))};
    const cfd::Mesh &mesh{build_result.mesh};
    const GeometryDiagnostics geometry{compute_geometry_diagnostics(mesh)};
    constexpr double orthogonality_tolerance{1024.0 * std::numeric_limits<double>::epsilon()};
    if (family.kind == MeshFamilyKind::GradedOrthogonalQuadrilateral &&
        geometry.maximum_correction_ratio > orthogonality_tolerance)
    {
        throw std::runtime_error("Graded quadrilateral fixture is not orthogonal to roundoff.");
    }
    const cfd::ScalarBoundaryConditions boundary_conditions{make_boundary_conditions(mesh)};
    const cfd::ScalarDiffusionOperator diffusion{mesh, diffusivity};
    cfd::ScalarLinearSystem system{mesh};
    diffusion.add_matrix_contributions(boundary_conditions, system);
    diffusion.add_boundary_rhs(boundary_conditions, system.rhs());

    std::vector<double> source_integrals(mesh.cell_count());
    double total_area{};
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        source_integrals[cell_id] = mesh.cell_areas()[cell_id] * source_density(mesh.cell_centers()[cell_id]);
        system.rhs()[cell_id] += source_integrals[cell_id];
        total_area += mesh.cell_areas()[cell_id];
    }
    if (!std::isfinite(total_area) || std::abs(total_area - domain_length * domain_height) > 1.0e-10)
    {
        throw std::runtime_error("Solution-verification mesh has an invalid total area.");
    }
    const std::vector<double> fixed_rhs{system.rhs().begin(), system.rhs().end()};
    const double source_norm{euclidean_norm(source_integrals)};
    if (!(source_norm > 0.0) || !std::isfinite(source_norm))
    {
        throw std::runtime_error("Solution-verification source norm is invalid.");
    }

    cfd::EigenConjugateGradientSolver solver{{1.0e-14, 10000}};
    solver.compute_matrix(system);
    cfd::CellScalarField numerical_solution{mesh.cell_count()};
    cfd::CellVectorField numerical_gradient{mesh.cell_count()};
    cfd::CellScalarField flux_balance{mesh.cell_count()};
    std::vector<double> matrix_product(mesh.cell_count());
    cfd::Index correction_count{};
    cfd::Index total_cg_iterations{};
    double final_full_residual{std::numeric_limits<double>::infinity()};
    double final_algebraic_residual{std::numeric_limits<double>::infinity()};

    const bool requires_correction{family.kind != MeshFamilyKind::GradedOrthogonalQuadrilateral};
    const cfd::Index solve_limit{requires_correction ? maximum_correction_count : cfd::Index{1}};
    for (cfd::Index correction_index = 0; correction_index < solve_limit; ++correction_index)
    {
        if (requires_correction)
        {
            cfd::compute_least_squares_gradient(mesh, numerical_solution, boundary_conditions, numerical_gradient);
        }
        std::copy(fixed_rhs.begin(), fixed_rhs.end(), system.rhs().begin());
        if (requires_correction)
        {
            diffusion.add_non_orthogonal_rhs(boundary_conditions, numerical_gradient, system.rhs());
        }

        const cfd::LinearSolveResult solve_result{solver.solve(system.rhs(), numerical_solution.values())};
        if (!solve_result.converged)
        {
            throw std::runtime_error("Inner conjugate-gradient solve did not converge.");
        }
        ++correction_count;
        total_cg_iterations += solve_result.iteration_count;
        final_algebraic_residual =
            normalized_algebraic_residual(system, numerical_solution, system.rhs(), matrix_product);

        cfd::compute_least_squares_gradient(mesh, numerical_solution, boundary_conditions, numerical_gradient);
        final_full_residual = normalized_full_residual(mesh, diffusion, boundary_conditions, numerical_solution,
                                                       numerical_gradient, source_integrals, source_norm, flux_balance);
        if (final_full_residual <= full_residual_tolerance)
        {
            break;
        }
    }
    if (final_full_residual > full_residual_tolerance)
    {
        std::ostringstream message;
        message << "Non-orthogonal correction failed for " << family.name << " at target h=" << level.target_size
                << "; residual=" << std::scientific << final_full_residual << '.';
        throw std::runtime_error(message.str());
    }

    cfd::Index internal_face_count{};
    for (const cfd::FaceAdjacency &adjacency : mesh.face_adjacencies())
    {
        internal_face_count += static_cast<cfd::Index>(!adjacency.is_boundary());
    }

    cfd::CellScalarField exact_solution{mesh.cell_count()};
    cfd::CellScalarField solution_error{mesh.cell_count()};
    cfd::CellScalarField absolute_solution_error{mesh.cell_count()};
    cfd::CellVectorField exact_gradient{mesh.cell_count()};
    cfd::CellVectorField gradient_error{mesh.cell_count()};
    cfd::CellScalarField full_residual_density{mesh.cell_count()};
    ErrorAccumulator solution_error_accumulator;
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        exact_solution[cell_id] = analytical_phi(mesh.cell_centers()[cell_id]);
        solution_error[cell_id] = numerical_solution[cell_id] - exact_solution[cell_id];
        absolute_solution_error[cell_id] = std::abs(solution_error[cell_id]);
        solution_error_accumulator.add(mesh.cell_areas()[cell_id], absolute_solution_error[cell_id]);

        exact_gradient[cell_id] = analytical_gradient(mesh.cell_centers()[cell_id]);
        gradient_error[cell_id] = {
            numerical_gradient[cell_id].x - exact_gradient[cell_id].x,
            numerical_gradient[cell_id].y - exact_gradient[cell_id].y,
        };
        full_residual_density[cell_id] =
            (source_integrals[cell_id] - flux_balance[cell_id]) / mesh.cell_areas()[cell_id];
    }

    if (write_vtu)
    {
        const std::filesystem::path output_path{
            output_directory / (std::string{family.file_prefix} + "_level_" + std::to_string(level_index) + ".vtu")};
        write_verification_vtu(mesh, exact_solution, numerical_solution, solution_error, absolute_solution_error,
                               numerical_gradient, exact_gradient, gradient_error, full_residual_density, output_path);
        written_paths.push_back(output_path);
    }

    const std::optional<GradingDiagnostics> grading{
        family.kind == MeshFamilyKind::GradedOrthogonalQuadrilateral
            ? std::optional<GradingDiagnostics>{compute_grading_diagnostics(level)}
            : std::nullopt};
    const std::optional<ExactGradientDiagnostic> exact_gradient_diagnostic{
        family.kind == MeshFamilyKind::GradedOrthogonalQuadrilateral
            ? std::nullopt
            : std::optional<ExactGradientDiagnostic>{
                  run_exact_gradient_diagnostic(mesh, diffusion, boundary_conditions, source_integrals, source_norm)}};
    return {
        .target_size = level.target_size,
        .characteristic_mesh_length = std::sqrt(total_area / static_cast<double>(mesh.cell_count())),
        .cell_count = mesh.cell_count(),
        .internal_face_count = internal_face_count,
        .matrix_nonzero_count = mesh.cell_count() + 2 * internal_face_count,
        .solution_errors = solution_error_accumulator.finish(),
        .rms_order = std::nullopt,
        .linf_order = std::nullopt,
        .correction_count = correction_count,
        .total_cg_iterations = total_cg_iterations,
        .normalized_full_residual = final_full_residual,
        .normalized_algebraic_residual = final_algebraic_residual,
        .geometry = geometry,
        .grading = grading,
        .exact_gradient_diagnostic = exact_gradient_diagnostic,
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

        if (coarse.exact_gradient_diagnostic.has_value() && fine.exact_gradient_diagnostic.has_value())
        {
            fine.exact_gradient_diagnostic->rms_order =
                observed_order(coarse.exact_gradient_diagnostic->solution_errors.area_weighted_rms_error,
                               fine.exact_gradient_diagnostic->solution_errors.area_weighted_rms_error,
                               coarse.characteristic_mesh_length, fine.characteristic_mesh_length);
            fine.exact_gradient_diagnostic->linf_order =
                observed_order(coarse.exact_gradient_diagnostic->solution_errors.linf_error,
                               fine.exact_gradient_diagnostic->solution_errors.linf_error,
                               coarse.characteristic_mesh_length, fine.characteristic_mesh_length);
        }
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

void print_solution_table(const MeshFamily &family, const std::vector<LevelResult> &results)
{
    std::cout << "\nComplete scalar diffusion solution - " << family.name << " - all Dirichlet\n\n"
              << std::left << std::setw(10) << "target_h" << std::setw(12) << "h_char" << std::setw(9) << "cells"
              << std::setw(11) << "int_faces" << std::setw(10) << "nnz" << std::setw(17) << "RMS(phi)" << std::setw(9)
              << "p_RMS" << std::setw(17) << "Linf(phi)" << std::setw(9) << "p_Linf" << std::setw(8) << "corr"
              << std::setw(10) << "CG_iter" << std::setw(15) << "full_res" << std::setw(15) << "alg_res" << '\n'
              << std::string(161, '-') << '\n';

    for (std::size_t level_index = 0; level_index < results.size(); ++level_index)
    {
        const LevelResult &result{results[level_index]};
        std::cout << std::fixed << std::setprecision(4) << std::setw(10) << result.target_size << std::scientific
                  << std::setprecision(4) << std::setw(12) << result.characteristic_mesh_length << std::fixed
                  << std::setw(9) << result.cell_count << std::setw(11) << result.internal_face_count << std::setw(10)
                  << result.matrix_nonzero_count << std::scientific << std::setprecision(7) << std::setw(17)
                  << result.solution_errors.area_weighted_rms_error;
        print_order(result.rms_order, level_index == 0);
        std::cout << std::scientific << std::setprecision(7) << std::setw(17) << result.solution_errors.linf_error;
        print_order(result.linf_order, level_index == 0);
        std::cout << std::fixed << std::setw(8) << result.correction_count << std::setw(10)
                  << result.total_cg_iterations << std::scientific << std::setprecision(4) << std::setw(15)
                  << result.normalized_full_residual << std::setw(15) << result.normalized_algebraic_residual << '\n';
    }
}

void print_geometry_table(const MeshFamily &family, const std::vector<LevelResult> &results)
{
    std::cout << "\nGeometry diagnostics - " << family.name << "\n\n"
              << std::left << std::setw(10) << "target_h" << std::setw(13) << "min_cos" << std::setw(13) << "mean_cos"
              << std::setw(13) << "max_cos" << std::setw(15) << "max_angle" << std::setw(13) << "min_lambda"
              << std::setw(13) << "max_lambda" << std::setw(14) << "max_|T|/|S|" << std::setw(14) << "min_area"
              << std::setw(14) << "max_area" << std::setw(12) << "area_ratio" << '\n'
              << std::string(140, '-') << '\n';
    for (const LevelResult &result : results)
    {
        const GeometryDiagnostics &geometry{result.geometry};
        std::cout << std::fixed << std::setprecision(4) << std::setw(10) << result.target_size << std::scientific
                  << std::setprecision(5) << std::setw(13) << geometry.minimum_normalized_projection << std::setw(13)
                  << geometry.mean_normalized_projection << std::setw(13) << geometry.maximum_normalized_projection
                  << std::fixed << std::setprecision(5) << std::setw(15) << geometry.maximum_non_orthogonality_degrees
                  << std::scientific << std::setprecision(5) << std::setw(13) << geometry.minimum_lambda
                  << std::setw(13) << geometry.maximum_lambda << std::setw(14) << geometry.maximum_correction_ratio
                  << std::setw(14) << geometry.minimum_cell_area << std::setw(14) << geometry.maximum_cell_area
                  << std::fixed << std::setprecision(4) << std::setw(12) << geometry.cell_area_ratio << '\n';
    }
}

void print_exact_gradient_diagnostic_table(const MeshFamily &family, const std::vector<LevelResult> &results)
{
    std::cout << "\nExact-cell-gradient correction diagnostic - " << family.name << " - all Dirichlet\n\n"
              << std::left << std::setw(10) << "target_h" << std::setw(12) << "h_char" << std::setw(9) << "cells"
              << std::setw(17) << "RMS(phi)" << std::setw(9) << "p_RMS" << std::setw(17) << "Linf(phi)" << std::setw(9)
              << "p_Linf" << std::setw(8) << "corr" << std::setw(10) << "CG_iter" << std::setw(15) << "full_res"
              << std::setw(15) << "alg_res" << '\n'
              << std::string(131, '-') << '\n';
    for (std::size_t level_index = 0; level_index < results.size(); ++level_index)
    {
        const LevelResult &result{results[level_index]};
        if (!result.exact_gradient_diagnostic.has_value())
        {
            throw std::runtime_error("Gmsh result has no exact-gradient diagnostic.");
        }
        const ExactGradientDiagnostic &diagnostic{*result.exact_gradient_diagnostic};
        std::cout << std::fixed << std::setprecision(4) << std::setw(10) << result.target_size << std::scientific
                  << std::setprecision(4) << std::setw(12) << result.characteristic_mesh_length << std::fixed
                  << std::setw(9) << result.cell_count << std::scientific << std::setprecision(7) << std::setw(17)
                  << diagnostic.solution_errors.area_weighted_rms_error;
        print_order(diagnostic.rms_order, level_index == 0);
        std::cout << std::scientific << std::setprecision(7) << std::setw(17) << diagnostic.solution_errors.linf_error;
        print_order(diagnostic.linf_order, level_index == 0);
        std::cout << std::fixed << std::setw(8) << diagnostic.correction_count << std::setw(10)
                  << diagnostic.total_cg_iterations << std::scientific << std::setprecision(4) << std::setw(15)
                  << diagnostic.normalized_full_residual << std::setw(15) << diagnostic.normalized_algebraic_residual
                  << '\n';
    }
}

void print_grading_table(const MeshFamily &family, const std::vector<LevelResult> &results)
{
    std::cout << "\nGrading diagnostics - " << family.name << "\n\n"
              << std::left << std::setw(10) << "target_h" << std::setw(16) << "min_width" << std::setw(16)
              << "max_width" << std::setw(14) << "width_ratio" << std::setw(16) << "max_|T|/|S|" << '\n'
              << std::string(72, '-') << '\n';
    for (const LevelResult &result : results)
    {
        if (!result.grading.has_value())
        {
            throw std::runtime_error("Graded mesh result has no grading diagnostics.");
        }
        std::cout << std::fixed << std::setprecision(4) << std::setw(10) << result.target_size << std::scientific
                  << std::setprecision(7) << std::setw(16) << result.grading->minimum_cell_width << std::setw(16)
                  << result.grading->maximum_cell_width << std::fixed << std::setprecision(5) << std::setw(14)
                  << result.grading->width_ratio << std::scientific << std::setprecision(7) << std::setw(16)
                  << result.geometry.maximum_correction_ratio << '\n';
    }
}

[[nodiscard]]
bool parse_write_vtu(const std::span<char *> arguments)
{
    if (arguments.size() == 1)
    {
        return false;
    }
    if (arguments.size() == 2 && std::string_view{arguments[1]} == "--write-vtu")
    {
        return true;
    }
    throw std::invalid_argument(
        "Unsupported arguments.\nUsage: cfd_unstructured_scalar_diffusion_convergence [--write-vtu]");
}

} // namespace

int main(const int argument_count, char **arguments)
{
    try
    {
        const std::span argument_view{arguments, static_cast<std::size_t>(argument_count)};
        const bool write_vtu{parse_write_vtu(argument_view)};
        const std::filesystem::path output_directory{"output/verification/unstructured_scalar_diffusion"};
        if (write_vtu)
        {
            std::filesystem::create_directories(output_directory);
        }

        std::vector<std::filesystem::path> written_paths;
        for (const MeshFamily &family : mesh_families)
        {
            std::vector<LevelResult> results;
            results.reserve(grid_levels.size());
            std::size_t level_index{};
            for (const GridLevel &level : grid_levels)
            {
                results.push_back(run_level(family, level, level_index, write_vtu, output_directory, written_paths));
                ++level_index;
            }
            compute_orders(results);
            print_solution_table(family, results);
            if (family.kind == MeshFamilyKind::GradedOrthogonalQuadrilateral)
            {
                print_grading_table(family, results);
            }
            else
            {
                print_geometry_table(family, results);
                print_exact_gradient_diagnostic_table(family, results);
            }
        }

        if (write_vtu)
        {
            std::cout << "\nVTU files written:\n";
            for (const std::filesystem::path &path : written_paths)
            {
                std::cout << "  " << path.string() << '\n';
            }
        }
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "Unstructured scalar diffusion verification failed: " << error.what() << '\n';
        return 1;
    }
}
