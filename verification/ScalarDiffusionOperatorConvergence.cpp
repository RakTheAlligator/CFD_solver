#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVectorField.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/io/VtkWriter.hpp"
#include "cfd/math/Point2.hpp"
#include "cfd/math/Vector2.hpp"
#include "cfd/mesh/Boundary.hpp"
#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/mesh/Types.hpp"
#include "cfd/meshing/RawMeshData.hpp"
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
constexpr double shear_factor{0.35};
constexpr double radians_to_degrees{180.0 / std::numbers::pi};

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

enum class MeshFamilyKind : std::uint8_t
{
    CartesianQuadrilateral,
    ShearedQuadrilateral,
    StructuredTriangle
};

struct MeshFamily
{
    MeshFamilyKind kind;
    std::string_view name;
    std::string_view file_prefix;
};

constexpr std::array mesh_families{
    MeshFamily{MeshFamilyKind::CartesianQuadrilateral, "orthogonal Cartesian QUAD", "cartesian_quad"},
    MeshFamily{MeshFamilyKind::ShearedQuadrilateral, "affine/sheared QUAD", "sheared_quad"},
    MeshFamily{MeshFamilyKind::StructuredTriangle, "fixed-diagonal structured TRI", "structured_tri"},
};

enum class BoundaryStudyKind : std::uint8_t
{
    AllDirichlet,
    Mixed
};

struct BoundaryStudy
{
    BoundaryStudyKind kind;
    std::string_view name;
    std::string_view file_prefix;
};

constexpr std::array boundary_studies{
    BoundaryStudy{BoundaryStudyKind::AllDirichlet, "all Dirichlet", "all_dirichlet"},
    BoundaryStudy{BoundaryStudyKind::Mixed, "left/bottom Dirichlet; right/top Neumann", "mixed"},
};

struct GeometryDiagnostics
{
    double minimum_normalized_projection{std::numeric_limits<double>::infinity()};
    double maximum_normalized_projection{};
    double maximum_non_orthogonality_degrees{};
    double minimum_lambda{std::numeric_limits<double>::infinity()};
    double maximum_lambda{-std::numeric_limits<double>::infinity()};
    double maximum_correction_ratio{};
};

struct ErrorCategoryResult
{
    ErrorStatistics errors;
    std::optional<double> rms_order;
    std::optional<double> linf_order;
};

struct MethodResult
{
    ErrorCategoryResult all_cells;
    ErrorCategoryResult boundary_cells;
    ErrorCategoryResult interior_cells;
};

struct BoundaryScalingResult
{
    cfd::Index cell_count{};
    double total_area{};
    double area_fraction{};
    std::optional<double> area_fraction_order;
};

struct LevelResult
{
    double delta{};
    double characteristic_mesh_length{};
    cfd::Index cell_count{};
    MethodResult exact_gradient;
    MethodResult wls_gradient;
    BoundaryScalingResult boundary_scaling;
    GeometryDiagnostics geometry;
};

[[nodiscard]]
double analytical_phi(const cfd::Point2 &point) noexcept
{
    const double x{point.x};
    const double y{point.y};
    return 1.0 + 0.3 * x - 0.2 * y + 0.4 * x * x + 0.25 * x * y + 0.35 * y * y + 0.1 * x * x * x + 0.08 * x * x * y -
           0.06 * x * y * y + 0.07 * y * y * y;
}

[[nodiscard]]
cfd::Vector2 analytical_gradient(const cfd::Point2 &point) noexcept
{
    const double x{point.x};
    const double y{point.y};
    return {
        0.3 + 0.8 * x + 0.25 * y + 0.3 * x * x + 0.16 * x * y - 0.06 * y * y,
        -0.2 + 0.25 * x + 0.7 * y + 0.08 * x * x - 0.12 * x * y + 0.21 * y * y,
    };
}

[[nodiscard]]
double source_density(const cfd::Point2 &point) noexcept
{
    return -1.5 - 0.48 * point.x - 0.58 * point.y;
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
cfd::Point2 mapped_node(const MeshFamilyKind kind, const double xi, const double eta) noexcept
{
    if (kind == MeshFamilyKind::ShearedQuadrilateral)
    {
        return {xi + shear_factor * eta, eta};
    }
    return {xi, eta};
}

void append_quadrilateral(cfd::RawMeshData &raw_mesh, const cfd::Index lower_left, const cfd::Index lower_right,
                          const cfd::Index upper_right, const cfd::Index upper_left)
{
    raw_mesh.cell_types.push_back(cfd::CellType::Quadrilateral);
    raw_mesh.cell_nodes.insert(raw_mesh.cell_nodes.end(), {lower_left, lower_right, upper_right, upper_left});
    raw_mesh.cell_node_offsets.push_back(raw_mesh.cell_nodes.size());
}

void append_triangles(cfd::RawMeshData &raw_mesh, const cfd::Index lower_left, const cfd::Index lower_right,
                      const cfd::Index upper_right, const cfd::Index upper_left)
{
    raw_mesh.cell_types.push_back(cfd::CellType::Triangle);
    raw_mesh.cell_nodes.insert(raw_mesh.cell_nodes.end(), {lower_left, lower_right, upper_right});
    raw_mesh.cell_node_offsets.push_back(raw_mesh.cell_nodes.size());

    raw_mesh.cell_types.push_back(cfd::CellType::Triangle);
    raw_mesh.cell_nodes.insert(raw_mesh.cell_nodes.end(), {lower_left, upper_right, upper_left});
    raw_mesh.cell_node_offsets.push_back(raw_mesh.cell_nodes.size());
}

void append_boundary_edge(cfd::RawMeshData &raw_mesh, const cfd::Index first_node, const cfd::Index second_node,
                          const std::string_view side)
{
    const cfd::BoundaryId boundary_id{raw_mesh.boundary_groups.size()};
    raw_mesh.boundary_groups.push_back({boundary_id, std::string{side} + "_" + std::to_string(boundary_id)});
    raw_mesh.boundary_edges.push_back({{first_node, second_node}, boundary_id});
}

[[nodiscard]]
cfd::RawMeshData make_raw_mesh(const MeshFamily &family, const GridLevel &level)
{
    constexpr double extent_tolerance{64.0 * std::numeric_limits<double>::epsilon()};
    if (std::abs(static_cast<double>(level.nx) * level.delta - domain_length) > extent_tolerance * domain_length ||
        std::abs(static_cast<double>(level.ny) * level.delta - domain_height) > extent_tolerance * domain_height)
    {
        throw std::runtime_error("Diffusion verification grid dimensions do not match the logical domain.");
    }

    const bool triangles{family.kind == MeshFamilyKind::StructuredTriangle};
    const cfd::Index node_count{(level.nx + 1) * (level.ny + 1)};
    const cfd::Index square_count{level.nx * level.ny};
    const cfd::Index cell_count{triangles ? 2 * square_count : square_count};
    const cfd::Index nodes_per_cell{triangles ? cfd::Index{3} : cfd::Index{4}};

    cfd::RawMeshData raw_mesh;
    raw_mesh.nodes.reserve(node_count);
    raw_mesh.cell_types.reserve(cell_count);
    raw_mesh.cell_nodes.reserve(cell_count * nodes_per_cell);
    raw_mesh.cell_node_offsets.reserve(cell_count + 1);
    raw_mesh.cell_node_offsets.push_back(0);

    for (cfd::Index j = 0; j <= level.ny; ++j)
    {
        for (cfd::Index i = 0; i <= level.nx; ++i)
        {
            raw_mesh.nodes.push_back(
                mapped_node(family.kind, static_cast<double>(i) * level.delta, static_cast<double>(j) * level.delta));
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

            if (triangles)
            {
                append_triangles(raw_mesh, lower_left, lower_right, upper_right, upper_left);
            }
            else
            {
                append_quadrilateral(raw_mesh, lower_left, lower_right, upper_right, upper_left);
            }
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
bool is_dirichlet_boundary(const BoundaryStudyKind study, const std::string_view boundary_name)
{
    if (study == BoundaryStudyKind::AllDirichlet)
    {
        return true;
    }
    return boundary_name.starts_with("left_") || boundary_name.starts_with("bottom_");
}

[[nodiscard]]
cfd::ScalarBoundaryConditions make_boundary_conditions(const cfd::Mesh &mesh, const BoundaryStudyKind study)
{
    std::vector<cfd::ScalarBoundaryCondition> conditions;
    conditions.reserve(mesh.boundary_groups().size());
    for (cfd::Index boundary_id = 0; boundary_id < mesh.boundary_groups().size(); ++boundary_id)
    {
        conditions.emplace_back(cfd::ScalarBoundaryConditionType::Neumann, 0.0);
    }
    std::vector<bool> assigned(mesh.boundary_groups().size(), false);

    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};
        if (!adjacency.is_boundary())
        {
            continue;
        }

        const cfd::BoundaryId boundary_id{mesh.face_boundary_ids()[face_id]};
        if (assigned[boundary_id])
        {
            throw std::runtime_error("A diffusion verification BoundaryId is shared by multiple faces.");
        }
        assigned[boundary_id] = true;

        const cfd::Point2 &face_center{mesh.face_centers()[face_id]};
        const std::string_view boundary_name{mesh.boundary_groups()[boundary_id].name};
        if (is_dirichlet_boundary(study, boundary_name))
        {
            conditions[boundary_id] = {
                cfd::ScalarBoundaryConditionType::Dirichlet,
                analytical_phi(face_center),
            };
            continue;
        }

        if (!boundary_name.starts_with("right_") && !boundary_name.starts_with("top_"))
        {
            throw std::runtime_error("Mixed diffusion verification encountered an unknown boundary side.");
        }
        const cfd::Vector2 exact_gradient{analytical_gradient(face_center)};
        const cfd::Vector2 &area_vector{mesh.face_area_vectors()[face_id]};
        const double inverse_face_length{1.0 / mesh.face_lengths()[face_id]};
        const cfd::Vector2 unit_normal{
            area_vector.x * inverse_face_length,
            area_vector.y * inverse_face_length,
        };
        conditions[boundary_id] = {
            cfd::ScalarBoundaryConditionType::Neumann,
            dot(exact_gradient, unit_normal),
        };
    }

    if (std::find(assigned.begin(), assigned.end(), false) != assigned.end())
    {
        throw std::runtime_error("A diffusion verification boundary group has no face.");
    }
    return {mesh.boundary_groups().size(), std::move(conditions)};
}

[[nodiscard]]
std::vector<bool> classify_boundary_cells(const cfd::Mesh &mesh)
{
    std::vector<bool> boundary_cells(mesh.cell_count(), false);
    for (const cfd::FaceAdjacency &adjacency : mesh.face_adjacencies())
    {
        if (adjacency.is_boundary())
        {
            boundary_cells[adjacency.owner] = true;
        }
    }
    return boundary_cells;
}

[[nodiscard]]
GeometryDiagnostics compute_geometry_diagnostics(const cfd::Mesh &mesh)
{
    GeometryDiagnostics diagnostics;
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
        diagnostics.maximum_non_orthogonality_degrees =
            std::max(diagnostics.maximum_non_orthogonality_degrees, std::acos(bounded_projection) * radians_to_degrees);
        diagnostics.minimum_lambda = std::min(diagnostics.minimum_lambda, lambda);
        diagnostics.maximum_lambda = std::max(diagnostics.maximum_lambda, lambda);
        diagnostics.maximum_correction_ratio = std::max(diagnostics.maximum_correction_ratio, correction_ratio);
    }

    if (internal_face_count == 0 || !std::isfinite(diagnostics.minimum_normalized_projection) ||
        !std::isfinite(diagnostics.maximum_normalized_projection) ||
        !std::isfinite(diagnostics.maximum_non_orthogonality_degrees) || !std::isfinite(diagnostics.minimum_lambda) ||
        !std::isfinite(diagnostics.maximum_lambda) || !std::isfinite(diagnostics.maximum_correction_ratio))
    {
        throw std::runtime_error("Diffusion verification produced invalid geometry diagnostics.");
    }
    return diagnostics;
}

void write_verification_vtu(const cfd::Mesh &mesh, const cfd::CellScalarField &phi,
                            const cfd::CellVectorField &exact_gradient, const cfd::CellVectorField &wls_gradient,
                            const cfd::CellVectorField &gradient_error,
                            const cfd::CellScalarField &gradient_error_magnitude,
                            const cfd::CellScalarField &exact_integrated_residual,
                            const cfd::CellScalarField &wls_integrated_residual,
                            const cfd::CellScalarField &exact_residual_density,
                            const cfd::CellScalarField &wls_residual_density, const std::filesystem::path &output_path)
{
    const std::array scalar_fields{
        cfd::VtkCellScalarData{"phi_exact", phi.values()},
        cfd::VtkCellScalarData{"gradient_error_magnitude", gradient_error_magnitude.values()},
        cfd::VtkCellScalarData{"diffusion_residual_density_exact_gradient", exact_residual_density.values()},
        cfd::VtkCellScalarData{"diffusion_residual_density_wls_gradient", wls_residual_density.values()},
        cfd::VtkCellScalarData{"diffusion_integrated_residual_exact_gradient", exact_integrated_residual.values()},
        cfd::VtkCellScalarData{"diffusion_integrated_residual_wls_gradient", wls_integrated_residual.values()},
    };
    const std::array vector_fields{
        cfd::VtkCellVectorData{"gradient_exact", exact_gradient.values()},
        cfd::VtkCellVectorData{"gradient_wls", wls_gradient.values()},
        cfd::VtkCellVectorData{"gradient_error", gradient_error.values()},
    };
    cfd::write_vtu(mesh, output_path, cfd::VtkCellData{.scalars = scalar_fields, .vectors = vector_fields});
}

[[nodiscard]]
LevelResult run_level(const MeshFamily &family, const BoundaryStudy &boundary_study, const GridLevel &level,
                      const std::size_t level_index, const bool write_vtu,
                      const std::filesystem::path &output_directory, std::vector<std::filesystem::path> &written_paths)
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_raw_mesh(family, level))};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarDiffusionOperator diffusion{mesh, diffusivity};
    const cfd::ScalarBoundaryConditions boundary_conditions{make_boundary_conditions(mesh, boundary_study.kind)};

    cfd::CellScalarField phi{mesh.cell_count()};
    cfd::CellVectorField exact_gradient{mesh.cell_count()};
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        phi[cell_id] = analytical_phi(mesh.cell_centers()[cell_id]);
        exact_gradient[cell_id] = analytical_gradient(mesh.cell_centers()[cell_id]);
    }

    cfd::CellVectorField wls_gradient{mesh.cell_count()};
    cfd::compute_least_squares_gradient(mesh, phi, boundary_conditions, wls_gradient);

    cfd::CellScalarField exact_flux_balance{mesh.cell_count()};
    cfd::CellScalarField wls_flux_balance{mesh.cell_count()};
    diffusion.compute_flux_balance(phi, boundary_conditions, exact_gradient, exact_flux_balance);
    diffusion.compute_flux_balance(phi, boundary_conditions, wls_gradient, wls_flux_balance);

    cfd::CellVectorField gradient_error{mesh.cell_count()};
    cfd::CellScalarField gradient_error_magnitude{mesh.cell_count()};
    cfd::CellScalarField exact_integrated_residual{mesh.cell_count()};
    cfd::CellScalarField wls_integrated_residual{mesh.cell_count()};
    cfd::CellScalarField exact_residual_density{mesh.cell_count()};
    cfd::CellScalarField wls_residual_density{mesh.cell_count()};
    const std::vector<bool> boundary_cells{classify_boundary_cells(mesh)};
    ErrorAccumulator exact_all_errors;
    ErrorAccumulator exact_boundary_errors;
    ErrorAccumulator exact_interior_errors;
    ErrorAccumulator wls_all_errors;
    ErrorAccumulator wls_boundary_errors;
    ErrorAccumulator wls_interior_errors;

    double total_area{};
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        const double cell_area{mesh.cell_areas()[cell_id]};
        const double source_integral{cell_area * source_density(mesh.cell_centers()[cell_id])};
        exact_integrated_residual[cell_id] = exact_flux_balance[cell_id] - source_integral;
        wls_integrated_residual[cell_id] = wls_flux_balance[cell_id] - source_integral;
        exact_residual_density[cell_id] = exact_integrated_residual[cell_id] / cell_area;
        wls_residual_density[cell_id] = wls_integrated_residual[cell_id] / cell_area;

        gradient_error[cell_id] = {
            wls_gradient[cell_id].x - exact_gradient[cell_id].x,
            wls_gradient[cell_id].y - exact_gradient[cell_id].y,
        };
        gradient_error_magnitude[cell_id] = std::hypot(gradient_error[cell_id].x, gradient_error[cell_id].y);

        if (!std::isfinite(exact_residual_density[cell_id]) || !std::isfinite(wls_residual_density[cell_id]))
        {
            throw std::runtime_error("Diffusion verification produced a non-finite residual density.");
        }
        const auto accumulate_error = [cell_area, is_boundary = boundary_cells[cell_id]](
                                          ErrorAccumulator &all_errors, ErrorAccumulator &boundary_errors,
                                          ErrorAccumulator &interior_errors, const double error) {
            all_errors.add(cell_area, error);
            (is_boundary ? boundary_errors : interior_errors).add(cell_area, error);
        };
        accumulate_error(exact_all_errors, exact_boundary_errors, exact_interior_errors,
                         std::abs(exact_residual_density[cell_id]));
        accumulate_error(wls_all_errors, wls_boundary_errors, wls_interior_errors,
                         std::abs(wls_residual_density[cell_id]));
        total_area += cell_area;
    }

    constexpr double expected_domain_area{domain_length * domain_height};
    const double area_summation_tolerance{64.0 * std::numeric_limits<double>::epsilon() *
                                          static_cast<double>(mesh.cell_count()) * expected_domain_area};
    if (std::abs(total_area - expected_domain_area) > area_summation_tolerance)
    {
        throw std::runtime_error("Diffusion verification mesh does not preserve the expected domain area.");
    }

    const ErrorStatistics exact_all_statistics{exact_all_errors.finish()};
    const ErrorStatistics exact_boundary_statistics{exact_boundary_errors.finish()};
    const ErrorStatistics exact_interior_statistics{exact_interior_errors.finish()};
    const ErrorStatistics wls_all_statistics{wls_all_errors.finish()};
    const ErrorStatistics wls_boundary_statistics{wls_boundary_errors.finish()};
    const ErrorStatistics wls_interior_statistics{wls_interior_errors.finish()};
    if (exact_boundary_statistics.cell_count == 0 || exact_interior_statistics.cell_count == 0 ||
        exact_boundary_statistics.cell_count + exact_interior_statistics.cell_count != mesh.cell_count())
    {
        throw std::runtime_error("Diffusion verification produced an invalid boundary-cell partition.");
    }
    const double boundary_area_fraction{exact_boundary_statistics.total_area / exact_all_statistics.total_area};
    if (!std::isfinite(boundary_area_fraction) || !(boundary_area_fraction > 0.0) || !(boundary_area_fraction < 1.0))
    {
        throw std::runtime_error("Diffusion verification produced an invalid boundary-cell area fraction.");
    }

    if (write_vtu)
    {
        const std::filesystem::path output_path{output_directory / (std::string{family.file_prefix} + "_" +
                                                                    std::string{boundary_study.file_prefix} +
                                                                    "_level_" + std::to_string(level_index) + ".vtu")};
        write_verification_vtu(mesh, phi, exact_gradient, wls_gradient, gradient_error, gradient_error_magnitude,
                               exact_integrated_residual, wls_integrated_residual, exact_residual_density,
                               wls_residual_density, output_path);
        written_paths.push_back(output_path);
    }

    return {
        .delta = level.delta,
        .characteristic_mesh_length = std::sqrt(total_area / static_cast<double>(mesh.cell_count())),
        .cell_count = mesh.cell_count(),
        .exact_gradient =
            {
                .all_cells = {.errors = exact_all_statistics, .rms_order = std::nullopt, .linf_order = std::nullopt},
                .boundary_cells = {.errors = exact_boundary_statistics,
                                   .rms_order = std::nullopt,
                                   .linf_order = std::nullopt},
                .interior_cells = {.errors = exact_interior_statistics,
                                   .rms_order = std::nullopt,
                                   .linf_order = std::nullopt},
            },
        .wls_gradient =
            {
                .all_cells = {.errors = wls_all_statistics, .rms_order = std::nullopt, .linf_order = std::nullopt},
                .boundary_cells = {.errors = wls_boundary_statistics,
                                   .rms_order = std::nullopt,
                                   .linf_order = std::nullopt},
                .interior_cells = {.errors = wls_interior_statistics,
                                   .rms_order = std::nullopt,
                                   .linf_order = std::nullopt},
            },
        .boundary_scaling =
            {
                .cell_count = exact_boundary_statistics.cell_count,
                .total_area = exact_boundary_statistics.total_area,
                .area_fraction = boundary_area_fraction,
                .area_fraction_order = std::nullopt,
            },
        .geometry = compute_geometry_diagnostics(mesh),
    };
}

void compute_orders(std::vector<LevelResult> &results)
{
    for (std::size_t level_index = 1; level_index < results.size(); ++level_index)
    {
        const LevelResult &coarse{results[level_index - 1]};
        LevelResult &fine{results[level_index]};

        const auto update_category = [&coarse, &fine](ErrorCategoryResult &fine_category,
                                                      const ErrorCategoryResult &coarse_category) {
            fine_category.rms_order = observed_order(
                coarse_category.errors.area_weighted_rms_error, fine_category.errors.area_weighted_rms_error,
                coarse.characteristic_mesh_length, fine.characteristic_mesh_length);
            fine_category.linf_order =
                observed_order(coarse_category.errors.linf_error, fine_category.errors.linf_error,
                               coarse.characteristic_mesh_length, fine.characteristic_mesh_length);
        };
        const auto update_method = [&update_category](MethodResult &fine_method, const MethodResult &coarse_method) {
            update_category(fine_method.all_cells, coarse_method.all_cells);
            update_category(fine_method.boundary_cells, coarse_method.boundary_cells);
            update_category(fine_method.interior_cells, coarse_method.interior_cells);
        };
        update_method(fine.exact_gradient, coarse.exact_gradient);
        update_method(fine.wls_gradient, coarse.wls_gradient);
        fine.boundary_scaling.area_fraction_order =
            observed_order(coarse.boundary_scaling.area_fraction, fine.boundary_scaling.area_fraction,
                           coarse.characteristic_mesh_length, fine.characteristic_mesh_length);
    }
}

void print_order(const std::optional<double> order, const bool first_level)
{
    if (order.has_value())
    {
        std::cout << std::fixed << std::setprecision(3) << std::setw(10) << *order;
        return;
    }
    std::cout << std::setw(10) << (first_level ? "-" : "n/a");
}

void print_convergence_table(const MeshFamily &family, const BoundaryStudy &boundary_study,
                             const std::vector<LevelResult> &results, const bool use_wls_gradient)
{
    std::cout << "\nScalar diffusion residual-density convergence - " << family.name << " - " << boundary_study.name
              << " - " << (use_wls_gradient ? "production WLS gradient" : "exact cell gradient") << "\n\n"
              << std::left << std::setw(11) << "delta" << std::setw(13) << "h_char" << std::setw(10) << "cells"
              << std::setw(18) << "RMS(tau)" << std::setw(10) << "p_RMS" << std::setw(18) << "Linf(tau)"
              << std::setw(10) << "p_Linf" << '\n'
              << std::string(90, '-') << '\n';

    for (std::size_t level_index = 0; level_index < results.size(); ++level_index)
    {
        const LevelResult &result{results[level_index]};
        const MethodResult &method{use_wls_gradient ? result.wls_gradient : result.exact_gradient};
        const ErrorCategoryResult &all_cells{method.all_cells};
        std::cout << std::fixed << std::setprecision(4) << std::setw(11) << result.delta << std::scientific
                  << std::setprecision(5) << std::setw(13) << result.characteristic_mesh_length << std::fixed
                  << std::setw(10) << result.cell_count << std::scientific << std::setprecision(8) << std::setw(18)
                  << all_cells.errors.area_weighted_rms_error;
        print_order(all_cells.rms_order, level_index == 0);
        std::cout << std::scientific << std::setprecision(8) << std::setw(18) << all_cells.errors.linf_error;
        print_order(all_cells.linf_order, level_index == 0);
        std::cout << '\n';
    }
}

void print_category_convergence_table(const MeshFamily &family, const BoundaryStudy &boundary_study,
                                      const std::vector<LevelResult> &results, const bool use_wls_gradient,
                                      const bool use_boundary_cells)
{
    const std::string_view category_name{use_boundary_cells ? "boundary cells" : "interior cells"};
    std::cout << "\nScalar diffusion " << category_name << " residual-density convergence - " << family.name << " - "
              << boundary_study.name << " - " << (use_wls_gradient ? "production WLS gradient" : "exact cell gradient")
              << "\n\n"
              << std::left << std::setw(11) << "delta" << std::setw(13) << "h_char" << std::setw(10) << "cells"
              << std::setw(18) << "RMS(tau)" << std::setw(10) << "p_RMS" << std::setw(18) << "Linf(tau)"
              << std::setw(10) << "p_Linf" << '\n'
              << std::string(90, '-') << '\n';

    for (std::size_t level_index = 0; level_index < results.size(); ++level_index)
    {
        const LevelResult &result{results[level_index]};
        const MethodResult &method{use_wls_gradient ? result.wls_gradient : result.exact_gradient};
        const ErrorCategoryResult &category{use_boundary_cells ? method.boundary_cells : method.interior_cells};
        std::cout << std::fixed << std::setprecision(4) << std::setw(11) << result.delta << std::scientific
                  << std::setprecision(5) << std::setw(13) << result.characteristic_mesh_length << std::fixed
                  << std::setw(10) << category.errors.cell_count << std::scientific << std::setprecision(8)
                  << std::setw(18) << category.errors.area_weighted_rms_error;
        print_order(category.rms_order, level_index == 0);
        std::cout << std::scientific << std::setprecision(8) << std::setw(18) << category.errors.linf_error;
        print_order(category.linf_order, level_index == 0);
        std::cout << '\n';
    }
}

void print_geometry_table(const MeshFamily &family, const std::vector<LevelResult> &results)
{
    std::cout << "\nScalar diffusion geometry - " << family.name << "\n\n"
              << std::left << std::setw(11) << "delta" << std::setw(10) << "cells" << std::setw(15) << "min_cos"
              << std::setw(15) << "max_cos" << std::setw(17) << "max_angle_deg" << std::setw(15) << "min_lambda"
              << std::setw(15) << "max_lambda" << std::setw(15) << "max_|T|/|S|" << '\n'
              << std::string(113, '-') << '\n';
    for (const LevelResult &result : results)
    {
        const GeometryDiagnostics &geometry{result.geometry};
        std::cout << std::fixed << std::setprecision(4) << std::setw(11) << result.delta << std::setw(10)
                  << result.cell_count << std::scientific << std::setprecision(7) << std::setw(15)
                  << geometry.minimum_normalized_projection << std::setw(15) << geometry.maximum_normalized_projection
                  << std::fixed << std::setprecision(6) << std::setw(17) << geometry.maximum_non_orthogonality_degrees
                  << std::scientific << std::setprecision(7) << std::setw(15) << geometry.minimum_lambda
                  << std::setw(15) << geometry.maximum_lambda << std::setw(15) << geometry.maximum_correction_ratio
                  << '\n';
    }
}

void print_boundary_scaling_table(const MeshFamily &family, const std::vector<LevelResult> &results)
{
    std::cout << "\nScalar diffusion boundary-cell scaling - " << family.name << "\n\n"
              << std::left << std::setw(11) << "delta" << std::setw(10) << "cells" << std::setw(16) << "boundary_cells"
              << std::setw(20) << "boundary_area" << std::setw(20) << "boundary_fraction" << std::setw(10)
              << "p_fraction" << '\n'
              << std::string(87, '-') << '\n';
    for (std::size_t level_index = 0; level_index < results.size(); ++level_index)
    {
        const LevelResult &result{results[level_index]};
        const BoundaryScalingResult &scaling{result.boundary_scaling};
        std::cout << std::fixed << std::setprecision(4) << std::setw(11) << result.delta << std::setw(10)
                  << result.cell_count << std::setw(16) << scaling.cell_count << std::scientific << std::setprecision(8)
                  << std::setw(20) << scaling.total_area << std::setw(20) << scaling.area_fraction;
        print_order(scaling.area_fraction_order, level_index == 0);
        std::cout << '\n';
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
        "Unsupported arguments.\nUsage: cfd_scalar_diffusion_operator_convergence [--write-vtu]");
}

} // namespace

int main(const int argument_count, char **arguments)
{
    try
    {
        const std::span argument_view{arguments, static_cast<std::size_t>(argument_count)};
        const bool write_vtu{parse_write_vtu(argument_view)};
        const std::filesystem::path output_directory{"output/verification/scalar_diffusion_operator"};
        if (write_vtu)
        {
            std::filesystem::create_directories(output_directory);
        }

        std::vector<std::filesystem::path> written_paths;
        for (const MeshFamily &family : mesh_families)
        {
            bool geometry_printed{};
            for (const BoundaryStudy &boundary_study : boundary_studies)
            {
                std::vector<LevelResult> results;
                results.reserve(grid_levels.size());
                std::size_t level_index{};
                for (const GridLevel &level : grid_levels)
                {
                    results.push_back(run_level(family, boundary_study, level, level_index, write_vtu, output_directory,
                                                written_paths));
                    ++level_index;
                }
                compute_orders(results);

                if (!geometry_printed)
                {
                    print_geometry_table(family, results);
                    print_boundary_scaling_table(family, results);
                    geometry_printed = true;
                }
                print_convergence_table(family, boundary_study, results, false);
                print_category_convergence_table(family, boundary_study, results, false, true);
                print_category_convergence_table(family, boundary_study, results, false, false);
                print_convergence_table(family, boundary_study, results, true);
                print_category_convergence_table(family, boundary_study, results, true, true);
                print_category_convergence_table(family, boundary_study, results, true, false);
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
        std::cerr << "Scalar diffusion verification failed: " << error.what() << '\n';
        return 1;
    }
}
