#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVectorField.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/io/VtkWriter.hpp"
#include "cfd/math/Point2.hpp"
#include "cfd/math/Vector2.hpp"
#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/meshing/RawMeshData.hpp"
#include "cfd/numerics/LeastSquaresGradient.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace
{

constexpr double domain_length{2.0};
constexpr double domain_height{1.0};

constexpr cfd::BoundaryId left_boundary_id{0};
constexpr cfd::BoundaryId right_boundary_id{1};
constexpr cfd::BoundaryId bottom_boundary_id{2};
constexpr cfd::BoundaryId top_boundary_id{3};

struct GridLevel
{
    double delta;
    cfd::Index nx;
    cfd::Index ny;
};

constexpr std::array grid_levels{
    GridLevel{0.20, 10, 5},
    GridLevel{0.10, 20, 10},
    GridLevel{0.05, 40, 20},
    GridLevel{0.025, 80, 40},
};

struct MeshFamily
{
    cfd::CellType cell_type;
    std::string_view name;
    std::string_view file_prefix;
};

constexpr std::array mesh_families{
    MeshFamily{cfd::CellType::Quadrilateral, "quadrilaterals", "quad"},
    MeshFamily{cfd::CellType::Triangle, "triangles", "triangle"},
};

struct ErrorStatistics
{
    cfd::Index cell_count{};
    double area_weighted_rms_error{};
    double linf_error{};
};

struct ErrorAccumulator
{
    cfd::Index cell_count{};
    double total_area{};
    double area_weighted_squared_error{};
    double linf_error{};

    void add(const double cell_area, const double error_magnitude) noexcept
    {
        ++cell_count;
        total_area += cell_area;
        area_weighted_squared_error += cell_area * error_magnitude * error_magnitude;
        linf_error = std::max(linf_error, error_magnitude);
    }

    [[nodiscard]]
    ErrorStatistics finish() const
    {
        if (cell_count == 0 || !std::isfinite(total_area) || !(total_area > 0.0) ||
            !std::isfinite(area_weighted_squared_error) || area_weighted_squared_error < 0.0 ||
            !std::isfinite(linf_error))
        {
            throw std::runtime_error("Structured verification produced invalid error statistics.");
        }

        return {
            .cell_count = cell_count,
            .area_weighted_rms_error = std::sqrt(area_weighted_squared_error / total_area),
            .linf_error = linf_error,
        };
    }
};

struct CategoryResult
{
    ErrorStatistics errors;
    std::optional<double> rms_order;
    std::optional<double> linf_order;
};

struct TaylorResult
{
    ErrorStatistics leading;
    ErrorStatistics remainder;
    std::optional<double> actual_rms_order;
    std::optional<double> leading_rms_order;
    std::optional<double> remainder_rms_order;
};

struct LevelResult
{
    double delta{};
    cfd::Index cell_count{};
    double characteristic_mesh_length{};
    CategoryResult all_cells;
    CategoryResult boundary_cells;
    CategoryResult interior_cells;
    TaylorResult interior_taylor;
    double maximum_opposite_pair_angular_defect{};
    double maximum_opposite_pair_distance_imbalance{};
    double minimum_cell_quality{};
    double maximum_cell_quality{};
};

struct Hessian2
{
    double m00{};
    double m01{};
    double m11{};
};

struct SymmetricSystem2
{
    double m00{};
    double m01{};
    double m11{};
    cfd::Vector2 right_hand_side{};
};

struct OppositePairMetrics
{
    bool valid{};
    double angular_defect{};
    double distance_imbalance{};
};

struct CellStencilDiagnostic
{
    bool is_boundary_cell{};
    bool leading_metric_valid{};
    cfd::Vector2 leading_predicted_error{};
    OppositePairMetrics opposite_pair_metrics;
};

struct DiagnosticFields
{
    explicit DiagnosticFields(const cfd::Index cell_count)
        : leading_error(cell_count), leading_error_magnitude(cell_count), leading_remainder(cell_count),
          is_boundary_cell(cell_count), leading_metric_valid(cell_count)
    {
    }

    cfd::CellVectorField leading_error;
    cfd::CellScalarField leading_error_magnitude;
    cfd::CellScalarField leading_remainder;
    cfd::CellScalarField is_boundary_cell;
    cfd::CellScalarField leading_metric_valid;
};

[[nodiscard]]
cfd::Index node_id(const cfd::Index i, const cfd::Index j, const cfd::Index nx) noexcept
{
    return j * (nx + 1) + i;
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
    // Both triangles are counter-clockwise and share the LL-to-UR diagonal.
    raw_mesh.cell_types.push_back(cfd::CellType::Triangle);
    raw_mesh.cell_nodes.insert(raw_mesh.cell_nodes.end(), {lower_left, lower_right, upper_right});
    raw_mesh.cell_node_offsets.push_back(raw_mesh.cell_nodes.size());

    raw_mesh.cell_types.push_back(cfd::CellType::Triangle);
    raw_mesh.cell_nodes.insert(raw_mesh.cell_nodes.end(), {lower_left, upper_right, upper_left});
    raw_mesh.cell_node_offsets.push_back(raw_mesh.cell_nodes.size());
}

[[nodiscard]]
cfd::RawMeshData make_structured_raw_mesh(const GridLevel &level, const cfd::CellType cell_type)
{
    constexpr double extent_tolerance{64.0 * std::numeric_limits<double>::epsilon()};
    if (std::abs(static_cast<double>(level.nx) * level.delta - domain_length) > extent_tolerance * domain_length ||
        std::abs(static_cast<double>(level.ny) * level.delta - domain_height) > extent_tolerance * domain_height)
    {
        throw std::runtime_error("Structured grid dimensions do not match the analytical domain.");
    }

    cfd::RawMeshData raw_mesh;

    const cfd::Index node_count{(level.nx + 1) * (level.ny + 1)};
    const cfd::Index square_count{level.nx * level.ny};
    const cfd::Index cell_count{cell_type == cfd::CellType::Quadrilateral ? square_count : 2 * square_count};
    const cfd::Index nodes_per_cell{cell_type == cfd::CellType::Quadrilateral ? cfd::Index{4} : cfd::Index{3}};

    raw_mesh.nodes.reserve(node_count);
    raw_mesh.cell_types.reserve(cell_count);
    raw_mesh.cell_nodes.reserve(cell_count * nodes_per_cell);
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

            if (cell_type == cfd::CellType::Quadrilateral)
            {
                append_quadrilateral(raw_mesh, lower_left, lower_right, upper_right, upper_left);
            }
            else
            {
                append_triangles(raw_mesh, lower_left, lower_right, upper_right, upper_left);
            }
        }
    }

    raw_mesh.boundary_groups = {
        {left_boundary_id, "left"},
        {right_boundary_id, "right"},
        {bottom_boundary_id, "bottom"},
        {top_boundary_id, "top"},
    };
    raw_mesh.boundary_edges.reserve(2 * (level.nx + level.ny));

    for (cfd::Index j = 0; j < level.ny; ++j)
    {
        raw_mesh.boundary_edges.push_back({{node_id(0, j, level.nx), node_id(0, j + 1, level.nx)}, left_boundary_id});
        raw_mesh.boundary_edges.push_back(
            {{node_id(level.nx, j, level.nx), node_id(level.nx, j + 1, level.nx)}, right_boundary_id});
    }

    for (cfd::Index i = 0; i < level.nx; ++i)
    {
        raw_mesh.boundary_edges.push_back({{node_id(i, 0, level.nx), node_id(i + 1, 0, level.nx)}, bottom_boundary_id});
        raw_mesh.boundary_edges.push_back(
            {{node_id(i, level.ny, level.nx), node_id(i + 1, level.ny, level.nx)}, top_boundary_id});
    }

    return raw_mesh;
}

[[nodiscard]]
double analytical_phi(const cfd::Point2 &point) noexcept
{
    return std::sin(point.x) + point.y * point.y;
}

[[nodiscard]]
cfd::Vector2 analytical_gradient(const cfd::Point2 &point) noexcept
{
    return {std::cos(point.x), 2.0 * point.y};
}

[[nodiscard]]
Hessian2 analytical_hessian(const cfd::Point2 &point) noexcept
{
    return {
        .m00 = -std::sin(point.x),
        .m01 = 0.0,
        .m11 = 2.0,
    };
}

[[nodiscard]]
double dot(const cfd::Vector2 &first, const cfd::Vector2 &second) noexcept
{
    return first.x * second.x + first.y * second.y;
}

[[nodiscard]]
double magnitude(const cfd::Vector2 &vector) noexcept
{
    return std::hypot(vector.x, vector.y);
}

[[nodiscard]]
cfd::Vector2 apply(const Hessian2 &matrix, const cfd::Vector2 &vector) noexcept
{
    return {
        matrix.m00 * vector.x + matrix.m01 * vector.y,
        matrix.m01 * vector.x + matrix.m11 * vector.y,
    };
}

void add_leading_observation(SymmetricSystem2 &system, const cfd::Vector2 &direction,
                             const double leading_defect) noexcept
{
    system.m00 += direction.x * direction.x;
    system.m01 += direction.x * direction.y;
    system.m11 += direction.y * direction.y;
    system.right_hand_side.x += direction.x * leading_defect;
    system.right_hand_side.y += direction.y * leading_defect;
}

[[nodiscard]]
cfd::Vector2 solve(const SymmetricSystem2 &system)
{
    constexpr double relative_singularity_tolerance{64.0 * std::numeric_limits<double>::epsilon()};

    const double determinant{system.m00 * system.m11 - system.m01 * system.m01};
    const double trace{system.m00 + system.m11};
    if (!(determinant > relative_singularity_tolerance * trace * trace))
    {
        throw std::runtime_error("Structured leading-Taylor diagnostic encountered a singular stencil.");
    }

    return {
        (system.right_hand_side.x * system.m11 - system.right_hand_side.y * system.m01) / determinant,
        (system.m00 * system.right_hand_side.y - system.m01 * system.right_hand_side.x) / determinant,
    };
}

struct PairMetrics
{
    double angular_defect{};
    double distance_imbalance{};
};

[[nodiscard]]
PairMetrics compute_pair_metrics(const cfd::Vector2 &first, const cfd::Vector2 &second)
{
    const double first_length{magnitude(first)};
    const double second_length{magnitude(second)};
    if (!(first_length > 0.0) || !(second_length > 0.0))
    {
        throw std::runtime_error("Structured symmetry diagnostic encountered a zero displacement.");
    }

    const double cosine{std::clamp(dot(first, second) / (first_length * second_length), -1.0, 1.0)};
    return {
        .angular_defect = 0.5 * (1.0 + cosine),
        .distance_imbalance = std::abs(first_length - second_length) / (first_length + second_length),
    };
}

[[nodiscard]]
OppositePairMetrics compute_opposite_pair_metrics(const std::array<cfd::Vector2, 4> &displacements)
{
    constexpr std::array<std::array<std::size_t, 4>, 3> partitions{
        std::array<std::size_t, 4>{0, 1, 2, 3},
        std::array<std::size_t, 4>{0, 2, 1, 3},
        std::array<std::size_t, 4>{0, 3, 1, 2},
    };

    double best_angular_defect{std::numeric_limits<double>::infinity()};
    double selected_distance_imbalance{};

    for (const std::array<std::size_t, 4> &partition : partitions)
    {
        const PairMetrics first_pair{
            compute_pair_metrics(displacements.at(partition.at(0)), displacements.at(partition.at(1)))};
        const PairMetrics second_pair{
            compute_pair_metrics(displacements.at(partition.at(2)), displacements.at(partition.at(3)))};
        const double angular_defect{std::max(first_pair.angular_defect, second_pair.angular_defect)};

        if (angular_defect < best_angular_defect)
        {
            best_angular_defect = angular_defect;
            selected_distance_imbalance = std::max(first_pair.distance_imbalance, second_pair.distance_imbalance);
        }
    }

    return {
        .valid = true,
        .angular_defect = best_angular_defect,
        .distance_imbalance = selected_distance_imbalance,
    };
}

[[nodiscard]]
CellStencilDiagnostic compute_cell_stencil_diagnostic(const cfd::Mesh &mesh, const cfd::Index cell_id)
{
    const auto cell_offsets{mesh.cell_node_offsets()};
    const auto cell_faces{mesh.cell_faces()};
    const auto face_adjacencies{mesh.face_adjacencies()};
    const auto cell_centers{mesh.cell_centers()};

    const cfd::Point2 &cell_center{cell_centers[cell_id]};
    const Hessian2 hessian{analytical_hessian(cell_center)};
    std::array<cfd::Vector2, 4> displacements{};
    std::size_t displacement_count{};
    bool is_boundary_cell{};

    const cfd::Index begin_offset{cell_offsets[cell_id]};
    const cfd::Index end_offset{cell_offsets[cell_id + 1]};

    for (cfd::Index position = begin_offset; position < end_offset; ++position)
    {
        const cfd::FaceAdjacency &adjacency{face_adjacencies[cell_faces[position]]};
        if (adjacency.is_boundary())
        {
            is_boundary_cell = true;
            continue;
        }

        const cfd::Index other_cell_id{adjacency.owner == cell_id ? adjacency.neighbor : adjacency.owner};
        if (displacement_count >= displacements.size())
        {
            throw std::runtime_error("Structured diagnostic encountered more than four internal neighbors.");
        }

        displacements.at(displacement_count) = {
            cell_centers[other_cell_id].x - cell_center.x,
            cell_centers[other_cell_id].y - cell_center.y,
        };
        ++displacement_count;
    }

    if (is_boundary_cell)
    {
        return {
            .is_boundary_cell = true,
            .leading_metric_valid = false,
            .leading_predicted_error = {},
            .opposite_pair_metrics = {},
        };
    }

    SymmetricSystem2 system;
    for (std::size_t displacement_index = 0; displacement_index < displacement_count; ++displacement_index)
    {
        const cfd::Vector2 &displacement{displacements.at(displacement_index)};
        const double distance{magnitude(displacement)};
        if (!std::isfinite(distance) || !(distance > 0.0))
        {
            throw std::runtime_error("Structured diagnostic encountered coincident cell centers.");
        }

        const cfd::Vector2 direction{displacement.x / distance, displacement.y / distance};
        const double leading_defect{0.5 * distance * dot(direction, apply(hessian, direction))};
        add_leading_observation(system, direction, leading_defect);
    }

    OppositePairMetrics opposite_pair_metrics;
    if (mesh.cell_types()[cell_id] == cfd::CellType::Quadrilateral)
    {
        if (displacement_count != displacements.size())
        {
            throw std::runtime_error("An interior structured quadrilateral does not have four internal neighbors.");
        }
        opposite_pair_metrics = compute_opposite_pair_metrics(displacements);
    }

    return {
        .is_boundary_cell = false,
        .leading_metric_valid = true,
        .leading_predicted_error = solve(system),
        .opposite_pair_metrics = opposite_pair_metrics,
    };
}

[[nodiscard]]
double outward_normal_derivative(const std::string_view boundary_name)
{
    if (boundary_name == "left")
    {
        return -1.0;
    }
    if (boundary_name == "right")
    {
        return std::cos(domain_length);
    }
    if (boundary_name == "bottom")
    {
        return 0.0;
    }
    if (boundary_name == "top")
    {
        return 2.0;
    }

    throw std::runtime_error("Unsupported structured verification boundary: " + std::string{boundary_name});
}

[[nodiscard]]
cfd::ScalarBoundaryConditions make_boundary_conditions(const cfd::Mesh &mesh)
{
    std::vector<cfd::ScalarBoundaryCondition> conditions;
    conditions.reserve(mesh.boundary_groups().size());

    // Phi is supplied analytically; this reconstructs its gradient and does
    // not solve a PDE, so pure Neumann data introduce no nullspace here.
    for (const cfd::BoundaryGroup &group : mesh.boundary_groups())
    {
        if (group.id != conditions.size())
        {
            throw std::runtime_error("Structured mesh boundary groups are not indexed contiguously.");
        }
        conditions.emplace_back(cfd::ScalarBoundaryConditionType::Neumann, outward_normal_derivative(group.name));
    }

    return {mesh.boundary_groups().size(), std::move(conditions)};
}

void write_verification_vtu(const cfd::Mesh &mesh, const cfd::CellScalarField &phi,
                            const cfd::CellVectorField &numerical_gradient, const cfd::CellVectorField &exact_gradient,
                            const cfd::CellVectorField &gradient_error, const cfd::CellScalarField &error_magnitude,
                            const DiagnosticFields &diagnostics, const std::filesystem::path &output_path)
{
    const std::array scalar_fields{
        cfd::VtkCellScalarData{"phi", phi.values()},
        cfd::VtkCellScalarData{"grad_phi_error_magnitude", error_magnitude.values()},
        cfd::VtkCellScalarData{"grad_phi_error_leading_magnitude", diagnostics.leading_error_magnitude.values()},
        cfd::VtkCellScalarData{"grad_phi_error_leading_remainder", diagnostics.leading_remainder.values()},
        cfd::VtkCellScalarData{"is_boundary_cell", diagnostics.is_boundary_cell.values()},
        cfd::VtkCellScalarData{"leading_metric_valid", diagnostics.leading_metric_valid.values()},
    };
    const std::array vector_fields{
        cfd::VtkCellVectorData{"grad_phi", numerical_gradient.values()},
        cfd::VtkCellVectorData{"grad_phi_exact", exact_gradient.values()},
        cfd::VtkCellVectorData{"grad_phi_error", gradient_error.values()},
        cfd::VtkCellVectorData{"grad_phi_error_leading", diagnostics.leading_error.values()},
    };
    const cfd::VtkCellData cell_data{
        .scalars = scalar_fields,
        .vectors = vector_fields,
    };

    cfd::write_vtu(mesh, output_path, cell_data);
}

[[nodiscard]]
LevelResult run_level(const MeshFamily &family, const GridLevel &level, const std::size_t level_index,
                      const bool write_vtu, const std::filesystem::path &output_directory,
                      std::vector<std::filesystem::path> &written_paths)
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_structured_raw_mesh(level, family.cell_type))};
    const cfd::Mesh &mesh{build_result.mesh};

    const cfd::Index square_count{level.nx * level.ny};
    const cfd::Index expected_cell_count{family.cell_type == cfd::CellType::Quadrilateral ? square_count
                                                                                          : 2 * square_count};
    if (mesh.cell_count() != expected_cell_count)
    {
        throw std::runtime_error("Structured mesh cell count differs from its Cartesian construction.");
    }

    cfd::CellScalarField phi{mesh.cell_count()};
    cfd::CellVectorField exact_gradient{mesh.cell_count()};
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        phi[cell_id] = analytical_phi(mesh.cell_centers()[cell_id]);
        exact_gradient[cell_id] = analytical_gradient(mesh.cell_centers()[cell_id]);
    }

    const cfd::ScalarBoundaryConditions boundary_conditions{make_boundary_conditions(mesh)};
    cfd::CellVectorField numerical_gradient{mesh.cell_count()};
    cfd::compute_least_squares_gradient(mesh, phi, boundary_conditions, numerical_gradient);

    cfd::CellVectorField gradient_error{mesh.cell_count()};
    cfd::CellScalarField error_magnitude{mesh.cell_count()};
    DiagnosticFields diagnostics{mesh.cell_count()};

    ErrorAccumulator all_errors;
    ErrorAccumulator boundary_errors;
    ErrorAccumulator interior_errors;
    ErrorAccumulator leading_errors;
    ErrorAccumulator leading_remainders;
    double maximum_angular_defect{};
    double maximum_distance_imbalance{};

    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        const double cell_area{mesh.cell_areas()[cell_id]};
        gradient_error[cell_id] = {
            numerical_gradient[cell_id].x - exact_gradient[cell_id].x,
            numerical_gradient[cell_id].y - exact_gradient[cell_id].y,
        };
        error_magnitude[cell_id] = magnitude(gradient_error[cell_id]);
        if (!std::isfinite(cell_area) || !(cell_area > 0.0) || !std::isfinite(error_magnitude[cell_id]))
        {
            throw std::runtime_error("Structured verification encountered invalid cell error data.");
        }

        all_errors.add(cell_area, error_magnitude[cell_id]);
        const CellStencilDiagnostic cell_diagnostic{compute_cell_stencil_diagnostic(mesh, cell_id)};
        diagnostics.is_boundary_cell[cell_id] = cell_diagnostic.is_boundary_cell ? 1.0 : 0.0;

        if (cell_diagnostic.is_boundary_cell)
        {
            boundary_errors.add(cell_area, error_magnitude[cell_id]);
            continue;
        }

        interior_errors.add(cell_area, error_magnitude[cell_id]);
        diagnostics.leading_metric_valid[cell_id] = 1.0;
        diagnostics.leading_error[cell_id] = cell_diagnostic.leading_predicted_error;
        diagnostics.leading_error_magnitude[cell_id] = magnitude(cell_diagnostic.leading_predicted_error);
        const cfd::Vector2 remainder{
            gradient_error[cell_id].x - cell_diagnostic.leading_predicted_error.x,
            gradient_error[cell_id].y - cell_diagnostic.leading_predicted_error.y,
        };
        diagnostics.leading_remainder[cell_id] = magnitude(remainder);
        leading_errors.add(cell_area, diagnostics.leading_error_magnitude[cell_id]);
        leading_remainders.add(cell_area, diagnostics.leading_remainder[cell_id]);

        if (cell_diagnostic.opposite_pair_metrics.valid)
        {
            maximum_angular_defect =
                std::max(maximum_angular_defect, cell_diagnostic.opposite_pair_metrics.angular_defect);
            maximum_distance_imbalance =
                std::max(maximum_distance_imbalance, cell_diagnostic.opposite_pair_metrics.distance_imbalance);
        }
    }

    const ErrorStatistics all_statistics{all_errors.finish()};
    const ErrorStatistics boundary_statistics{boundary_errors.finish()};
    const ErrorStatistics interior_statistics{interior_errors.finish()};
    const ErrorStatistics leading_statistics{leading_errors.finish()};
    const ErrorStatistics remainder_statistics{leading_remainders.finish()};
    const double characteristic_mesh_length{std::sqrt(all_errors.total_area / static_cast<double>(mesh.cell_count()))};

    const auto [minimum_quality,
                maximum_quality]{std::minmax_element(mesh.cell_qualities().begin(), mesh.cell_qualities().end())};
    if (minimum_quality == mesh.cell_qualities().end() || !std::isfinite(*minimum_quality) ||
        !std::isfinite(*maximum_quality))
    {
        throw std::runtime_error("Structured verification encountered invalid mesh quality data.");
    }

    if (write_vtu)
    {
        const std::filesystem::path output_path{
            output_directory / (std::string{family.file_prefix} + "_level_" + std::to_string(level_index) + ".vtu")};
        write_verification_vtu(mesh, phi, numerical_gradient, exact_gradient, gradient_error, error_magnitude,
                               diagnostics, output_path);
        written_paths.push_back(output_path);
    }

    return {
        .delta = level.delta,
        .cell_count = mesh.cell_count(),
        .characteristic_mesh_length = characteristic_mesh_length,
        .all_cells = {.errors = all_statistics, .rms_order = std::nullopt, .linf_order = std::nullopt},
        .boundary_cells = {.errors = boundary_statistics, .rms_order = std::nullopt, .linf_order = std::nullopt},
        .interior_cells = {.errors = interior_statistics, .rms_order = std::nullopt, .linf_order = std::nullopt},
        .interior_taylor =
            {
                .leading = leading_statistics,
                .remainder = remainder_statistics,
                .actual_rms_order = std::nullopt,
                .leading_rms_order = std::nullopt,
                .remainder_rms_order = std::nullopt,
            },
        .maximum_opposite_pair_angular_defect = maximum_angular_defect,
        .maximum_opposite_pair_distance_imbalance = maximum_distance_imbalance,
        .minimum_cell_quality = *minimum_quality,
        .maximum_cell_quality = *maximum_quality,
    };
}

[[nodiscard]]
std::optional<double> observed_order(const double coarse_error, const double fine_error, const double coarse_delta,
                                     const double fine_delta)
{
    if (!std::isfinite(coarse_error) || !std::isfinite(fine_error) || coarse_error < 0.0 || fine_error < 0.0 ||
        !std::isfinite(coarse_delta) || !std::isfinite(fine_delta) || !(coarse_delta > 0.0) || !(fine_delta > 0.0))
    {
        throw std::runtime_error("Cannot form a structured convergence order from invalid inputs.");
    }
    if (coarse_error == 0.0 || fine_error == 0.0)
    {
        return std::nullopt;
    }

    const double order{std::log(coarse_error / fine_error) / std::log(coarse_delta / fine_delta)};
    if (!std::isfinite(order))
    {
        throw std::runtime_error("Structured convergence order is non-finite.");
    }
    return order;
}

void compute_orders(std::vector<LevelResult> &results)
{
    for (std::size_t level_index = 1; level_index < results.size(); ++level_index)
    {
        LevelResult &fine{results[level_index]};
        const LevelResult &coarse{results[level_index - 1]};

        const auto update_category = [&coarse, &fine](CategoryResult &fine_category,
                                                      const CategoryResult &coarse_category) {
            fine_category.rms_order =
                observed_order(coarse_category.errors.area_weighted_rms_error,
                               fine_category.errors.area_weighted_rms_error, coarse.delta, fine.delta);
            fine_category.linf_order = observed_order(coarse_category.errors.linf_error,
                                                      fine_category.errors.linf_error, coarse.delta, fine.delta);
        };

        update_category(fine.all_cells, coarse.all_cells);
        update_category(fine.boundary_cells, coarse.boundary_cells);
        update_category(fine.interior_cells, coarse.interior_cells);

        fine.interior_taylor.actual_rms_order = fine.interior_cells.rms_order;
        fine.interior_taylor.leading_rms_order =
            observed_order(coarse.interior_taylor.leading.area_weighted_rms_error,
                           fine.interior_taylor.leading.area_weighted_rms_error, coarse.delta, fine.delta);
        fine.interior_taylor.remainder_rms_order =
            observed_order(coarse.interior_taylor.remainder.area_weighted_rms_error,
                           fine.interior_taylor.remainder.area_weighted_rms_error, coarse.delta, fine.delta);
    }
}

void print_optional(const std::optional<double> value, const std::size_t width)
{
    if (value.has_value())
    {
        std::cout << std::fixed << std::setprecision(3) << std::setw(static_cast<int>(width)) << *value;
        return;
    }
    std::cout << std::setw(static_cast<int>(width)) << "-";
}

void print_overall_table(const MeshFamily &family, const std::vector<LevelResult> &results)
{
    std::cout << "\nStructured WLS overall convergence - " << family.name << "\n\n"
              << std::left << std::setw(10) << "delta" << std::setw(10) << "cells" << std::setw(16) << "h_char"
              << std::setw(16) << "RMS_all" << std::setw(10) << "p_RMS" << std::setw(16) << "Linf_all" << std::setw(10)
              << "p_Linf" << '\n'
              << std::string(88, '-') << '\n';

    for (const LevelResult &result : results)
    {
        std::cout << std::fixed << std::setprecision(3) << std::setw(10) << result.delta << std::setw(10)
                  << result.cell_count << std::scientific << std::setprecision(6) << std::setw(16)
                  << result.characteristic_mesh_length << std::setw(16)
                  << result.all_cells.errors.area_weighted_rms_error;
        print_optional(result.all_cells.rms_order, 10);
        std::cout << std::scientific << std::setprecision(6) << std::setw(16) << result.all_cells.errors.linf_error;
        print_optional(result.all_cells.linf_order, 10);
        std::cout << '\n';
    }
}

void print_category_table(const MeshFamily &family, const std::vector<LevelResult> &results)
{
    std::cout << "\nStructured WLS boundary/interior convergence - " << family.name << "\n\n"
              << std::left << std::setw(16) << "family" << std::setw(10) << "delta" << std::setw(12) << "category"
              << std::setw(10) << "cells" << std::setw(16) << "RMS" << std::setw(10) << "p_RMS" << std::setw(16)
              << "Linf" << std::setw(10) << "p_Linf" << '\n'
              << std::string(100, '-') << '\n';

    for (const LevelResult &result : results)
    {
        const std::array categories{
            std::pair{"boundary", &result.boundary_cells},
            std::pair{"interior", &result.interior_cells},
        };
        for (const auto &[category_name, category] : categories)
        {
            std::cout << std::setw(16) << family.name << std::fixed << std::setprecision(3) << std::setw(10)
                      << result.delta << std::setw(12) << category_name << std::setw(10) << category->errors.cell_count
                      << std::scientific << std::setprecision(6) << std::setw(16)
                      << category->errors.area_weighted_rms_error;
            print_optional(category->rms_order, 10);
            std::cout << std::scientific << std::setprecision(6) << std::setw(16) << category->errors.linf_error;
            print_optional(category->linf_order, 10);
            std::cout << '\n';
        }
    }
}

void print_taylor_table(const MeshFamily &family, const std::vector<LevelResult> &results)
{
    std::cout << "\nInterior leading-Taylor diagnostic - " << family.name << "\n\n"
              << std::left << std::setw(16) << "family" << std::setw(10) << "delta" << std::setw(16) << "RMS_actual"
              << std::setw(10) << "p_actual" << std::setw(16) << "RMS_leading" << std::setw(10) << "p_leading"
              << std::setw(16) << "RMS_remainder" << std::setw(12) << "p_remainder" << std::setw(16) << "Linf_actual"
              << std::setw(16) << "Linf_leading" << std::setw(16) << "Linf_remainder" << '\n'
              << std::string(154, '-') << '\n';

    for (const LevelResult &result : results)
    {
        std::cout << std::setw(16) << family.name << std::fixed << std::setprecision(3) << std::setw(10) << result.delta
                  << std::scientific << std::setprecision(6) << std::setw(16)
                  << result.interior_cells.errors.area_weighted_rms_error;
        print_optional(result.interior_taylor.actual_rms_order, 10);
        std::cout << std::scientific << std::setprecision(6) << std::setw(16)
                  << result.interior_taylor.leading.area_weighted_rms_error;
        print_optional(result.interior_taylor.leading_rms_order, 10);
        std::cout << std::scientific << std::setprecision(6) << std::setw(16)
                  << result.interior_taylor.remainder.area_weighted_rms_error;
        print_optional(result.interior_taylor.remainder_rms_order, 12);
        std::cout << std::scientific << std::setprecision(6) << std::setw(16) << result.interior_cells.errors.linf_error
                  << std::setw(16) << result.interior_taylor.leading.linf_error << std::setw(16)
                  << result.interior_taylor.remainder.linf_error << '\n';
    }
}

void print_geometry_table(const MeshFamily &family, const std::vector<LevelResult> &results)
{
    std::cout << "\nStructured geometry diagnostics - " << family.name << "\n\n"
              << std::left << std::setw(16) << "family" << std::setw(10) << "delta" << std::setw(20)
              << "min_cell_quality" << std::setw(20) << "max_cell_quality";
    if (family.cell_type == cfd::CellType::Quadrilateral)
    {
        std::cout << std::setw(24) << "max_angular_defect" << std::setw(24) << "max_distance_imbalance";
    }
    std::cout << '\n' << std::string(family.cell_type == cfd::CellType::Quadrilateral ? 114 : 66, '-') << '\n';

    for (const LevelResult &result : results)
    {
        std::cout << std::setw(16) << family.name << std::fixed << std::setprecision(3) << std::setw(10) << result.delta
                  << std::scientific << std::setprecision(12) << std::setw(20) << result.minimum_cell_quality
                  << std::setw(20) << result.maximum_cell_quality;
        if (family.cell_type == cfd::CellType::Quadrilateral)
        {
            std::cout << std::setw(24) << result.maximum_opposite_pair_angular_defect << std::setw(24)
                      << result.maximum_opposite_pair_distance_imbalance;
        }
        std::cout << '\n';
    }
}

[[nodiscard]]
bool parse_write_vtu_flag(const std::span<char *> arguments)
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
        "Unsupported arguments.\nUsage: cfd_structured_least_squares_gradient_convergence [--write-vtu]");
}

} // namespace

int main(const int argument_count, char **arguments)
{
    try
    {
        const std::span argument_view{arguments, static_cast<std::size_t>(argument_count)};
        const bool write_vtu{parse_write_vtu_flag(argument_view)};
        const std::filesystem::path output_directory{"output/verification/structured_least_squares_gradient"};
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
            print_overall_table(family, results);
            print_category_table(family, results);
            print_taylor_table(family, results);
            print_geometry_table(family, results);
        }

        if (write_vtu)
        {
            std::cout << "\nVTU output:\n";
            for (const std::filesystem::path &path : written_paths)
            {
                std::cout << "  " << path.string() << '\n';
            }
        }
    }
    catch (const std::exception &error)
    {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }

    return 0;
}
