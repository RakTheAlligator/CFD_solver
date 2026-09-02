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
#include <cstdint>
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

// Holding this dimensionless distortion fixed isolates refinement from a
// changing mesh-shape perturbation.
constexpr double transition_alpha{0.25};

// This is far above structured-grid roundoff and far below the deliberate
// O(alpha) distance imbalance, so it separates the two numerical regimes.
constexpr double stencil_transition_tolerance{1.0e-12};

constexpr double geometry_classification_tolerance{4096.0 * std::numeric_limits<double>::epsilon() * domain_length};
constexpr double cell_quality_tolerance{4096.0 * std::numeric_limits<double>::epsilon()};

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
    GridLevel{0.20, 10, 5},   GridLevel{0.10, 20, 10},    GridLevel{0.05, 40, 20},
    GridLevel{0.025, 80, 40}, GridLevel{0.0125, 160, 80},
};

struct MeshVariant
{
    std::string_view name;
    std::string_view file_prefix;
    bool shift_internal_column;
};

constexpr std::array mesh_variants{
    MeshVariant{"regular", "regular", false},
    MeshVariant{"controlled_transition", "transition", true},
};

struct ErrorStatistics
{
    cfd::Index cell_count{};
    double total_area{};
    double area_fraction{};
    double area_weighted_squared_error{};
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
    ErrorStatistics finish(const double domain_area) const
    {
        if (!std::isfinite(total_area) || total_area < 0.0 || !std::isfinite(area_weighted_squared_error) ||
            area_weighted_squared_error < 0.0 || !std::isfinite(linf_error) || !std::isfinite(domain_area) ||
            !(domain_area > 0.0))
        {
            throw std::runtime_error("Controlled-transition verification produced invalid error statistics.");
        }

        if (cell_count == 0)
        {
            return {};
        }

        if (!(total_area > 0.0))
        {
            throw std::runtime_error("A non-empty controlled-transition category has zero area.");
        }

        return {
            .cell_count = cell_count,
            .total_area = total_area,
            .area_fraction = total_area / domain_area,
            .area_weighted_squared_error = area_weighted_squared_error,
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
    ErrorStatistics actual;
    ErrorStatistics leading;
    ErrorStatistics remainder;
    std::optional<double> actual_rms_order;
    std::optional<double> leading_rms_order;
    std::optional<double> remainder_rms_order;
};

struct MetricStatistics
{
    cfd::Index cell_count{};
    double mean_angular_defect{};
    double maximum_angular_defect{};
    double mean_distance_imbalance{};
    double maximum_distance_imbalance{};
};

struct MetricAccumulator
{
    cfd::Index cell_count{};
    double angular_defect_sum{};
    double maximum_angular_defect{};
    double distance_imbalance_sum{};
    double maximum_distance_imbalance{};

    void add(const double angular_defect, const double distance_imbalance) noexcept
    {
        ++cell_count;
        angular_defect_sum += angular_defect;
        maximum_angular_defect = std::max(maximum_angular_defect, angular_defect);
        distance_imbalance_sum += distance_imbalance;
        maximum_distance_imbalance = std::max(maximum_distance_imbalance, distance_imbalance);
    }

    [[nodiscard]]
    MetricStatistics finish() const
    {
        if (cell_count == 0)
        {
            return {};
        }

        const double count{static_cast<double>(cell_count)};
        return {
            .cell_count = cell_count,
            .mean_angular_defect = angular_defect_sum / count,
            .maximum_angular_defect = maximum_angular_defect,
            .mean_distance_imbalance = distance_imbalance_sum / count,
            .maximum_distance_imbalance = maximum_distance_imbalance,
        };
    }
};

struct LevelResult
{
    double delta{};
    cfd::Index cell_count{};
    CategoryResult all_cells;
    CategoryResult boundary_cells;
    CategoryResult regular_interior_cells;
    CategoryResult transition_interior_cells;
    CategoryResult shape_distorted_cells;
    CategoryResult high_quality_transition_cells;
    ErrorStatistics transition_stencil_cells;
    TaylorResult regular_interior_taylor;
    TaylorResult transition_interior_taylor;
    MetricStatistics regular_interior_metrics;
    MetricStatistics transition_interior_metrics;
    double boundary_squared_error_fraction{};
    double regular_interior_squared_error_fraction{};
    double transition_interior_squared_error_fraction{};
    double transition_area_over_delta{};
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

struct CellClassification
{
    std::vector<std::uint8_t> is_boundary_cell;
    std::vector<std::uint8_t> is_shape_distorted_cell;
    std::vector<std::uint8_t> is_transition_stencil_cell;
};

struct DiagnosticFields
{
    explicit DiagnosticFields(const cfd::Index cell_count)
        : leading_error(cell_count), leading_error_magnitude(cell_count), leading_remainder(cell_count),
          is_boundary_cell(cell_count), is_shape_distorted_cell(cell_count), is_transition_stencil_cell(cell_count),
          opposite_pair_angular_defect(cell_count), opposite_pair_distance_imbalance(cell_count),
          opposite_pair_metrics_valid(cell_count)
    {
    }

    cfd::CellVectorField leading_error;
    cfd::CellScalarField leading_error_magnitude;
    cfd::CellScalarField leading_remainder;
    cfd::CellScalarField is_boundary_cell;
    cfd::CellScalarField is_shape_distorted_cell;
    cfd::CellScalarField is_transition_stencil_cell;
    cfd::CellScalarField opposite_pair_angular_defect;
    cfd::CellScalarField opposite_pair_distance_imbalance;
    cfd::CellScalarField opposite_pair_metrics_valid;
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

[[nodiscard]]
cfd::RawMeshData make_raw_mesh(const GridLevel &level, const bool shift_internal_column)
{
    constexpr double extent_tolerance{64.0 * std::numeric_limits<double>::epsilon()};
    if (level.nx % 2 != 0 ||
        std::abs(static_cast<double>(level.nx) * level.delta - domain_length) > extent_tolerance * domain_length ||
        std::abs(static_cast<double>(level.ny) * level.delta - domain_height) > extent_tolerance * domain_height)
    {
        throw std::runtime_error("Controlled-transition grid dimensions do not match the analytical domain.");
    }

    cfd::RawMeshData raw_mesh;
    const cfd::Index node_count{(level.nx + 1) * (level.ny + 1)};
    const cfd::Index cell_count{level.nx * level.ny};
    raw_mesh.nodes.reserve(node_count);
    raw_mesh.cell_types.reserve(cell_count);
    raw_mesh.cell_nodes.reserve(4 * cell_count);
    raw_mesh.cell_node_offsets.reserve(cell_count + 1);
    raw_mesh.cell_node_offsets.push_back(0);

    const cfd::Index transition_node_index{level.nx / 2};
    for (cfd::Index j = 0; j <= level.ny; ++j)
    {
        for (cfd::Index i = 0; i <= level.nx; ++i)
        {
            double x{static_cast<double>(i) * level.delta};
            if (shift_internal_column && i == transition_node_index)
            {
                x += transition_alpha * level.delta;
            }
            raw_mesh.nodes.push_back({x, static_cast<double>(j) * level.delta});
        }
    }

    for (cfd::Index j = 0; j < level.ny; ++j)
    {
        for (cfd::Index i = 0; i < level.nx; ++i)
        {
            append_quadrilateral(raw_mesh, node_id(i, j, level.nx), node_id(i + 1, j, level.nx),
                                 node_id(i + 1, j + 1, level.nx), node_id(i, j + 1, level.nx));
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
        throw std::runtime_error("Controlled-transition leading-Taylor diagnostic encountered a singular stencil.");
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
    if (!std::isfinite(first_length) || !std::isfinite(second_length) || !(first_length > 0.0) ||
        !(second_length > 0.0))
    {
        throw std::runtime_error("Controlled-transition diagnostic encountered an invalid displacement.");
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
bool cell_is_boundary(const cfd::Mesh &mesh, const cfd::Index cell_id)
{
    const auto cell_offsets{mesh.cell_node_offsets()};
    const auto cell_faces{mesh.cell_faces()};
    const auto face_adjacencies{mesh.face_adjacencies()};
    for (cfd::Index position = cell_offsets[cell_id]; position < cell_offsets[cell_id + 1]; ++position)
    {
        if (face_adjacencies[cell_faces[position]].is_boundary())
        {
            return true;
        }
    }
    return false;
}

[[nodiscard]]
bool cell_is_shape_distorted(const cfd::Mesh &mesh, const cfd::Index cell_id, const double delta)
{
    const auto cell_offsets{mesh.cell_node_offsets()};
    const auto cell_nodes{mesh.cell_nodes()};
    const auto nodes{mesh.nodes()};
    double minimum_x{std::numeric_limits<double>::infinity()};
    double maximum_x{-std::numeric_limits<double>::infinity()};
    double minimum_y{std::numeric_limits<double>::infinity()};
    double maximum_y{-std::numeric_limits<double>::infinity()};

    for (cfd::Index position = cell_offsets[cell_id]; position < cell_offsets[cell_id + 1]; ++position)
    {
        const cfd::Node &node{nodes[cell_nodes[position]]};
        minimum_x = std::min(minimum_x, node.x);
        maximum_x = std::max(maximum_x, node.x);
        minimum_y = std::min(minimum_y, node.y);
        maximum_y = std::max(maximum_y, node.y);
    }

    const double width{maximum_x - minimum_x};
    const double height{maximum_y - minimum_y};
    if (!std::isfinite(width) || !std::isfinite(height) || !(width > 0.0) || !(height > 0.0) ||
        std::abs(height - delta) > geometry_classification_tolerance)
    {
        throw std::runtime_error("Controlled-transition cell geometry differs from the intended Cartesian rows.");
    }

    return std::abs(width - delta) > geometry_classification_tolerance;
}

[[nodiscard]]
CellClassification classify_cells(const cfd::Mesh &mesh, const GridLevel &level, const bool shifted)
{
    CellClassification classification{
        .is_boundary_cell = std::vector<std::uint8_t>(mesh.cell_count()),
        .is_shape_distorted_cell = std::vector<std::uint8_t>(mesh.cell_count()),
        .is_transition_stencil_cell = std::vector<std::uint8_t>(mesh.cell_count()),
    };

    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        classification.is_boundary_cell[cell_id] = cell_is_boundary(mesh, cell_id) ? 1U : 0U;
        classification.is_shape_distorted_cell[cell_id] = cell_is_shape_distorted(mesh, cell_id, level.delta) ? 1U : 0U;
    }

    const auto cell_offsets{mesh.cell_node_offsets()};
    const auto cell_faces{mesh.cell_faces()};
    const auto face_adjacencies{mesh.face_adjacencies()};
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        bool is_transition{classification.is_shape_distorted_cell[cell_id] != 0U};
        for (cfd::Index position = cell_offsets[cell_id]; !is_transition && position < cell_offsets[cell_id + 1];
             ++position)
        {
            const cfd::FaceAdjacency &adjacency{face_adjacencies[cell_faces[position]]};
            if (adjacency.is_boundary())
            {
                continue;
            }

            const cfd::Index other_cell_id{adjacency.owner == cell_id ? adjacency.neighbor : adjacency.owner};
            is_transition = classification.is_shape_distorted_cell[other_cell_id] != 0U;
        }
        classification.is_transition_stencil_cell[cell_id] = is_transition ? 1U : 0U;
    }

    const auto count_flagged = [](const std::vector<std::uint8_t> &flags) {
        return static_cast<cfd::Index>(std::ranges::count(flags, std::uint8_t{1}));
    };
    const cfd::Index expected_shape_count{shifted ? 2 * level.ny : 0};
    const cfd::Index expected_transition_count{shifted ? 4 * level.ny : 0};
    if (count_flagged(classification.is_shape_distorted_cell) != expected_shape_count ||
        count_flagged(classification.is_transition_stencil_cell) != expected_transition_count)
    {
        throw std::runtime_error("Geometry-based transition classification does not match the intended band.");
    }

    return classification;
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

    for (cfd::Index position = cell_offsets[cell_id]; position < cell_offsets[cell_id + 1]; ++position)
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
            throw std::runtime_error("Controlled-transition diagnostic encountered more than four neighbors.");
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
    if (displacement_count != displacements.size())
    {
        throw std::runtime_error("An interior controlled quadrilateral does not have four internal neighbors.");
    }

    SymmetricSystem2 system;
    for (const cfd::Vector2 &displacement : displacements)
    {
        const double distance{magnitude(displacement)};
        if (!std::isfinite(distance) || !(distance > 0.0))
        {
            throw std::runtime_error("Controlled-transition diagnostic encountered coincident cell centers.");
        }

        const cfd::Vector2 direction{displacement.x / distance, displacement.y / distance};
        const double leading_defect{0.5 * distance * dot(direction, apply(hessian, direction))};
        add_leading_observation(system, direction, leading_defect);
    }

    return {
        .is_boundary_cell = false,
        .leading_metric_valid = true,
        .leading_predicted_error = solve(system),
        .opposite_pair_metrics = compute_opposite_pair_metrics(displacements),
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
    throw std::runtime_error("Unsupported controlled-transition boundary: " + std::string{boundary_name});
}

[[nodiscard]]
cfd::ScalarBoundaryConditions make_boundary_conditions(const cfd::Mesh &mesh)
{
    std::vector<cfd::ScalarBoundaryCondition> conditions;
    conditions.reserve(mesh.boundary_groups().size());

    // Phi is supplied analytically, so pure Neumann data do not introduce a
    // nullspace in this gradient-reconstruction experiment.
    for (const cfd::BoundaryGroup &group : mesh.boundary_groups())
    {
        if (group.id != conditions.size())
        {
            throw std::runtime_error("Controlled mesh boundary groups are not indexed contiguously.");
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
        cfd::VtkCellScalarData{"is_shape_distorted_cell", diagnostics.is_shape_distorted_cell.values()},
        cfd::VtkCellScalarData{"is_transition_stencil_cell", diagnostics.is_transition_stencil_cell.values()},
        cfd::VtkCellScalarData{"opposite_pair_angular_defect", diagnostics.opposite_pair_angular_defect.values()},
        cfd::VtkCellScalarData{"opposite_pair_distance_imbalance",
                               diagnostics.opposite_pair_distance_imbalance.values()},
        cfd::VtkCellScalarData{"opposite_pair_metrics_valid", diagnostics.opposite_pair_metrics_valid.values()},
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
double squared_error_fraction(const ErrorStatistics &category, const ErrorStatistics &all)
{
    if (!(all.area_weighted_squared_error > 0.0))
    {
        throw std::runtime_error("Controlled-transition total squared error is not positive.");
    }
    return category.area_weighted_squared_error / all.area_weighted_squared_error;
}

[[nodiscard]]
LevelResult run_level(const MeshVariant &variant, const GridLevel &level, const std::size_t level_index,
                      const bool write_vtu, const std::filesystem::path &output_directory,
                      std::vector<std::filesystem::path> &written_paths)
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_raw_mesh(level, variant.shift_internal_column))};
    const cfd::Mesh &mesh{build_result.mesh};
    if (mesh.cell_count() != level.nx * level.ny)
    {
        throw std::runtime_error("Controlled mesh cell count differs from its Cartesian construction.");
    }

    const CellClassification classification{classify_cells(mesh, level, variant.shift_internal_column)};
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
    ErrorAccumulator regular_interior_errors;
    ErrorAccumulator transition_interior_errors;
    ErrorAccumulator shape_distorted_errors;
    ErrorAccumulator high_quality_transition_errors;
    ErrorAccumulator transition_stencil_errors;
    ErrorAccumulator regular_leading_errors;
    ErrorAccumulator regular_remainder_errors;
    ErrorAccumulator transition_leading_errors;
    ErrorAccumulator transition_remainder_errors;
    MetricAccumulator regular_metrics;
    MetricAccumulator transition_metrics;

    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        const double cell_area{mesh.cell_areas()[cell_id]};
        const double cell_quality{mesh.cell_qualities()[cell_id]};
        gradient_error[cell_id] = {
            numerical_gradient[cell_id].x - exact_gradient[cell_id].x,
            numerical_gradient[cell_id].y - exact_gradient[cell_id].y,
        };
        error_magnitude[cell_id] = magnitude(gradient_error[cell_id]);
        if (!std::isfinite(cell_area) || !(cell_area > 0.0) || !std::isfinite(cell_quality) ||
            !std::isfinite(error_magnitude[cell_id]))
        {
            throw std::runtime_error("Controlled-transition verification encountered invalid cell data.");
        }

        const bool is_boundary{classification.is_boundary_cell[cell_id] != 0U};
        const bool is_shape_distorted{classification.is_shape_distorted_cell[cell_id] != 0U};
        const bool is_transition{classification.is_transition_stencil_cell[cell_id] != 0U};
        diagnostics.is_boundary_cell[cell_id] = is_boundary ? 1.0 : 0.0;
        diagnostics.is_shape_distorted_cell[cell_id] = is_shape_distorted ? 1.0 : 0.0;
        diagnostics.is_transition_stencil_cell[cell_id] = is_transition ? 1.0 : 0.0;

        all_errors.add(cell_area, error_magnitude[cell_id]);
        if (is_shape_distorted)
        {
            shape_distorted_errors.add(cell_area, error_magnitude[cell_id]);
            if (!(cell_quality < 1.0 - cell_quality_tolerance))
            {
                throw std::runtime_error("A shape-distorted cell was not identified by the quality metric.");
            }
        }
        if (is_transition)
        {
            transition_stencil_errors.add(cell_area, error_magnitude[cell_id]);
            if (std::abs(cell_quality - 1.0) <= cell_quality_tolerance)
            {
                high_quality_transition_errors.add(cell_area, error_magnitude[cell_id]);
            }
        }

        const CellStencilDiagnostic cell_diagnostic{compute_cell_stencil_diagnostic(mesh, cell_id)};
        if (cell_diagnostic.is_boundary_cell != is_boundary)
        {
            throw std::runtime_error("Controlled-transition boundary classifications disagree.");
        }
        if (is_boundary)
        {
            boundary_errors.add(cell_area, error_magnitude[cell_id]);
            continue;
        }

        if (!cell_diagnostic.leading_metric_valid || !cell_diagnostic.opposite_pair_metrics.valid)
        {
            throw std::runtime_error("An interior controlled-transition diagnostic is unexpectedly invalid.");
        }
        const bool metric_detects_transition{cell_diagnostic.opposite_pair_metrics.distance_imbalance >
                                             stencil_transition_tolerance};
        if (metric_detects_transition != is_transition)
        {
            throw std::runtime_error("Topology/shape transition classification disagrees with stencil geometry.");
        }

        diagnostics.leading_error[cell_id] = cell_diagnostic.leading_predicted_error;
        diagnostics.leading_error_magnitude[cell_id] = magnitude(cell_diagnostic.leading_predicted_error);
        const cfd::Vector2 remainder{
            gradient_error[cell_id].x - cell_diagnostic.leading_predicted_error.x,
            gradient_error[cell_id].y - cell_diagnostic.leading_predicted_error.y,
        };
        diagnostics.leading_remainder[cell_id] = magnitude(remainder);
        diagnostics.opposite_pair_angular_defect[cell_id] = cell_diagnostic.opposite_pair_metrics.angular_defect;
        diagnostics.opposite_pair_distance_imbalance[cell_id] =
            cell_diagnostic.opposite_pair_metrics.distance_imbalance;
        diagnostics.opposite_pair_metrics_valid[cell_id] = 1.0;

        if (is_transition)
        {
            transition_interior_errors.add(cell_area, error_magnitude[cell_id]);
            transition_leading_errors.add(cell_area, diagnostics.leading_error_magnitude[cell_id]);
            transition_remainder_errors.add(cell_area, diagnostics.leading_remainder[cell_id]);
            transition_metrics.add(cell_diagnostic.opposite_pair_metrics.angular_defect,
                                   cell_diagnostic.opposite_pair_metrics.distance_imbalance);
        }
        else
        {
            regular_interior_errors.add(cell_area, error_magnitude[cell_id]);
            regular_leading_errors.add(cell_area, diagnostics.leading_error_magnitude[cell_id]);
            regular_remainder_errors.add(cell_area, diagnostics.leading_remainder[cell_id]);
            regular_metrics.add(cell_diagnostic.opposite_pair_metrics.angular_defect,
                                cell_diagnostic.opposite_pair_metrics.distance_imbalance);
        }
    }

    const double measured_domain_area{all_errors.total_area};
    if (std::abs(measured_domain_area - domain_length * domain_height) > geometry_classification_tolerance)
    {
        throw std::runtime_error("Controlled-transition mesh does not preserve the rectangular domain area.");
    }

    const ErrorStatistics all_statistics{all_errors.finish(measured_domain_area)};
    const ErrorStatistics boundary_statistics{boundary_errors.finish(measured_domain_area)};
    const ErrorStatistics regular_interior_statistics{regular_interior_errors.finish(measured_domain_area)};
    const ErrorStatistics transition_interior_statistics{transition_interior_errors.finish(measured_domain_area)};
    const ErrorStatistics shape_distorted_statistics{shape_distorted_errors.finish(measured_domain_area)};
    const ErrorStatistics high_quality_transition_statistics{
        high_quality_transition_errors.finish(measured_domain_area)};
    const ErrorStatistics transition_stencil_statistics{transition_stencil_errors.finish(measured_domain_area)};
    const ErrorStatistics regular_leading_statistics{regular_leading_errors.finish(measured_domain_area)};
    const ErrorStatistics regular_remainder_statistics{regular_remainder_errors.finish(measured_domain_area)};
    const ErrorStatistics transition_leading_statistics{transition_leading_errors.finish(measured_domain_area)};
    const ErrorStatistics transition_remainder_statistics{transition_remainder_errors.finish(measured_domain_area)};

    const double partition_area{boundary_statistics.total_area + regular_interior_statistics.total_area +
                                transition_interior_statistics.total_area};
    if (std::abs(partition_area - measured_domain_area) > geometry_classification_tolerance)
    {
        throw std::runtime_error("Controlled-transition error categories do not partition the mesh.");
    }
    if (transition_stencil_statistics.cell_count !=
            shape_distorted_statistics.cell_count + high_quality_transition_statistics.cell_count ||
        std::abs(transition_stencil_statistics.total_area - shape_distorted_statistics.total_area -
                 high_quality_transition_statistics.total_area) > geometry_classification_tolerance)
    {
        throw std::runtime_error("Transition cells do not partition into distorted and high-quality cells.");
    }

    const auto [minimum_quality,
                maximum_quality]{std::minmax_element(mesh.cell_qualities().begin(), mesh.cell_qualities().end())};
    if (minimum_quality == mesh.cell_qualities().end() || !std::isfinite(*minimum_quality) ||
        !std::isfinite(*maximum_quality))
    {
        throw std::runtime_error("Controlled-transition verification encountered invalid cell quality.");
    }

    if (write_vtu)
    {
        const std::filesystem::path output_path{
            output_directory / (std::string{variant.file_prefix} + "_level_" + std::to_string(level_index) + ".vtu")};
        write_verification_vtu(mesh, phi, numerical_gradient, exact_gradient, gradient_error, error_magnitude,
                               diagnostics, output_path);
        written_paths.push_back(output_path);
    }

    return {
        .delta = level.delta,
        .cell_count = mesh.cell_count(),
        .all_cells = {.errors = all_statistics, .rms_order = std::nullopt, .linf_order = std::nullopt},
        .boundary_cells = {.errors = boundary_statistics, .rms_order = std::nullopt, .linf_order = std::nullopt},
        .regular_interior_cells = {.errors = regular_interior_statistics,
                                   .rms_order = std::nullopt,
                                   .linf_order = std::nullopt},
        .transition_interior_cells = {.errors = transition_interior_statistics,
                                      .rms_order = std::nullopt,
                                      .linf_order = std::nullopt},
        .shape_distorted_cells = {.errors = shape_distorted_statistics,
                                  .rms_order = std::nullopt,
                                  .linf_order = std::nullopt},
        .high_quality_transition_cells = {.errors = high_quality_transition_statistics,
                                          .rms_order = std::nullopt,
                                          .linf_order = std::nullopt},
        .transition_stencil_cells = transition_stencil_statistics,
        .regular_interior_taylor =
            {
                .actual = regular_interior_statistics,
                .leading = regular_leading_statistics,
                .remainder = regular_remainder_statistics,
                .actual_rms_order = std::nullopt,
                .leading_rms_order = std::nullopt,
                .remainder_rms_order = std::nullopt,
            },
        .transition_interior_taylor =
            {
                .actual = transition_interior_statistics,
                .leading = transition_leading_statistics,
                .remainder = transition_remainder_statistics,
                .actual_rms_order = std::nullopt,
                .leading_rms_order = std::nullopt,
                .remainder_rms_order = std::nullopt,
            },
        .regular_interior_metrics = regular_metrics.finish(),
        .transition_interior_metrics = transition_metrics.finish(),
        .boundary_squared_error_fraction = squared_error_fraction(boundary_statistics, all_statistics),
        .regular_interior_squared_error_fraction = squared_error_fraction(regular_interior_statistics, all_statistics),
        .transition_interior_squared_error_fraction =
            squared_error_fraction(transition_interior_statistics, all_statistics),
        .transition_area_over_delta = transition_stencil_statistics.total_area / level.delta,
        .minimum_cell_quality = *minimum_quality,
        .maximum_cell_quality = *maximum_quality,
    };
}

[[nodiscard]]
std::optional<double> observed_order(const ErrorStatistics &coarse, const ErrorStatistics &fine,
                                     const double coarse_delta, const double fine_delta, const bool use_linf)
{
    if (coarse.cell_count == 0 || fine.cell_count == 0)
    {
        return std::nullopt;
    }
    const double coarse_error{use_linf ? coarse.linf_error : coarse.area_weighted_rms_error};
    const double fine_error{use_linf ? fine.linf_error : fine.area_weighted_rms_error};
    if (!std::isfinite(coarse_error) || !std::isfinite(fine_error) || coarse_error < 0.0 || fine_error < 0.0 ||
        coarse_error == 0.0 || fine_error == 0.0)
    {
        return std::nullopt;
    }
    const double order{std::log(coarse_error / fine_error) / std::log(coarse_delta / fine_delta)};
    if (!std::isfinite(order))
    {
        throw std::runtime_error("Controlled-transition convergence order is non-finite.");
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
                observed_order(coarse_category.errors, fine_category.errors, coarse.delta, fine.delta, false);
            fine_category.linf_order =
                observed_order(coarse_category.errors, fine_category.errors, coarse.delta, fine.delta, true);
        };
        update_category(fine.all_cells, coarse.all_cells);
        update_category(fine.boundary_cells, coarse.boundary_cells);
        update_category(fine.regular_interior_cells, coarse.regular_interior_cells);
        update_category(fine.transition_interior_cells, coarse.transition_interior_cells);
        update_category(fine.shape_distorted_cells, coarse.shape_distorted_cells);
        update_category(fine.high_quality_transition_cells, coarse.high_quality_transition_cells);

        const auto update_taylor = [&coarse, &fine](TaylorResult &fine_taylor, const TaylorResult &coarse_taylor) {
            fine_taylor.actual_rms_order =
                observed_order(coarse_taylor.actual, fine_taylor.actual, coarse.delta, fine.delta, false);
            fine_taylor.leading_rms_order =
                observed_order(coarse_taylor.leading, fine_taylor.leading, coarse.delta, fine.delta, false);
            fine_taylor.remainder_rms_order =
                observed_order(coarse_taylor.remainder, fine_taylor.remainder, coarse.delta, fine.delta, false);
        };
        update_taylor(fine.regular_interior_taylor, coarse.regular_interior_taylor);
        update_taylor(fine.transition_interior_taylor, coarse.transition_interior_taylor);
    }
}

void print_order(const std::optional<double> order, const bool first_level, const int width)
{
    if (order.has_value())
    {
        std::cout << std::fixed << std::setprecision(3) << std::setw(width) << *order;
        return;
    }
    std::cout << std::setw(width) << (first_level ? "-" : "n/a");
}

void print_error_value(const ErrorStatistics &statistics, const double value, const int width)
{
    if (statistics.cell_count == 0)
    {
        std::cout << std::setw(width) << "n/a";
        return;
    }
    std::cout << std::scientific << std::setprecision(6) << std::setw(width) << value;
}

void print_overall_table(const MeshVariant &variant, const std::vector<LevelResult> &results)
{
    std::cout << "\nControlled quadrilateral WLS overall convergence - " << variant.name << "\n\n"
              << std::left << std::setw(12) << "delta" << std::setw(12) << "cells" << std::setw(18) << "RMS_all"
              << std::setw(10) << "p_RMS" << std::setw(18) << "Linf_all" << std::setw(10) << "p_Linf" << '\n'
              << std::string(80, '-') << '\n';

    for (std::size_t level_index = 0; level_index < results.size(); ++level_index)
    {
        const LevelResult &result{results[level_index]};
        std::cout << std::fixed << std::setprecision(4) << std::setw(12) << result.delta << std::setw(12)
                  << result.cell_count << std::scientific << std::setprecision(6) << std::setw(18)
                  << result.all_cells.errors.area_weighted_rms_error;
        print_order(result.all_cells.rms_order, level_index == 0, 10);
        std::cout << std::scientific << std::setprecision(6) << std::setw(18) << result.all_cells.errors.linf_error;
        print_order(result.all_cells.linf_order, level_index == 0, 10);
        std::cout << '\n';
    }
}

void print_category_table(const MeshVariant &variant, const std::vector<LevelResult> &results)
{
    std::cout << "\nControlled quadrilateral error categories - " << variant.name << "\n\n"
              << std::left << std::setw(10) << "delta" << std::setw(26) << "category" << std::setw(10) << "cells"
              << std::setw(16) << "area" << std::setw(14) << "area_frac" << std::setw(16) << "RMS" << std::setw(9)
              << "p_RMS" << std::setw(16) << "Linf" << std::setw(9) << "p_Linf" << std::setw(14) << "S/S_total" << '\n'
              << std::string(150, '-') << '\n';

    for (std::size_t level_index = 0; level_index < results.size(); ++level_index)
    {
        const LevelResult &result{results[level_index]};
        const std::array categories{
            std::pair{"all", &result.all_cells},
            std::pair{"boundary", &result.boundary_cells},
            std::pair{"regular_interior", &result.regular_interior_cells},
            std::pair{"transition_interior", &result.transition_interior_cells},
            std::pair{"shape_distorted", &result.shape_distorted_cells},
            std::pair{"high_quality_transition", &result.high_quality_transition_cells},
        };

        for (const auto &[name, category] : categories)
        {
            const ErrorStatistics &errors{category->errors};
            std::cout << std::fixed << std::setprecision(4) << std::setw(10) << result.delta << std::setw(26) << name
                      << std::setw(10) << errors.cell_count;
            print_error_value(errors, errors.total_area, 16);
            print_error_value(errors, errors.area_fraction, 14);
            print_error_value(errors, errors.area_weighted_rms_error, 16);
            print_order(category->rms_order, level_index == 0, 9);
            print_error_value(errors, errors.linf_error, 16);
            print_order(category->linf_order, level_index == 0, 9);
            if (errors.cell_count == 0)
            {
                std::cout << std::setw(14) << "n/a";
            }
            else
            {
                std::cout << std::scientific << std::setprecision(6) << std::setw(14)
                          << squared_error_fraction(errors, result.all_cells.errors);
            }
            std::cout << '\n';
        }
    }
}

void print_transition_scaling_table(const std::vector<LevelResult> &results)
{
    std::cout << "\nControlled transition-band population scaling\n\n"
              << std::left << std::setw(12) << "delta" << std::setw(18) << "transition_cells" << std::setw(20)
              << "transition_area" << std::setw(22) << "transition_area_frac" << std::setw(20) << "transition_area/d"
              << '\n'
              << std::string(92, '-') << '\n';
    for (const LevelResult &result : results)
    {
        std::cout << std::fixed << std::setprecision(4) << std::setw(12) << result.delta << std::setw(18)
                  << result.transition_stencil_cells.cell_count << std::scientific << std::setprecision(8)
                  << std::setw(20) << result.transition_stencil_cells.total_area << std::setw(22)
                  << result.transition_stencil_cells.area_fraction << std::setw(20) << result.transition_area_over_delta
                  << '\n';
    }
}

void print_squared_error_contributions(const MeshVariant &variant, const std::vector<LevelResult> &results)
{
    std::cout << "\nSquared-error contribution fractions - " << variant.name << "\n\n"
              << std::left << std::setw(12) << "delta" << std::setw(22) << "boundary" << std::setw(22)
              << "regular_interior" << std::setw(22) << "transition_interior" << '\n'
              << std::string(78, '-') << '\n';
    for (const LevelResult &result : results)
    {
        std::cout << std::fixed << std::setprecision(4) << std::setw(12) << result.delta << std::scientific
                  << std::setprecision(8) << std::setw(22) << result.boundary_squared_error_fraction << std::setw(22)
                  << result.regular_interior_squared_error_fraction << std::setw(22)
                  << result.transition_interior_squared_error_fraction << '\n';
    }
}

void print_geometry_metrics(const MeshVariant &variant, const std::vector<LevelResult> &results)
{
    std::cout << "\nOpposite-pair geometry metrics - " << variant.name << "\n\n"
              << std::left << std::setw(10) << "delta" << std::setw(24) << "category" << std::setw(10) << "cells"
              << std::setw(18) << "mean_ang" << std::setw(18) << "max_ang" << std::setw(18) << "mean_dist"
              << std::setw(18) << "max_dist" << '\n'
              << std::string(116, '-') << '\n';
    for (const LevelResult &result : results)
    {
        const std::array categories{
            std::pair{"regular_interior", &result.regular_interior_metrics},
            std::pair{"transition_interior", &result.transition_interior_metrics},
        };
        for (const auto &[name, metrics] : categories)
        {
            std::cout << std::fixed << std::setprecision(4) << std::setw(10) << result.delta << std::setw(24) << name
                      << std::setw(10) << metrics->cell_count;
            if (metrics->cell_count == 0)
            {
                std::cout << std::setw(18) << "n/a" << std::setw(18) << "n/a" << std::setw(18) << "n/a" << std::setw(18)
                          << "n/a";
            }
            else
            {
                std::cout << std::scientific << std::setprecision(8) << std::setw(18) << metrics->mean_angular_defect
                          << std::setw(18) << metrics->maximum_angular_defect << std::setw(18)
                          << metrics->mean_distance_imbalance << std::setw(18) << metrics->maximum_distance_imbalance;
            }
            std::cout << '\n';
        }
    }
}

void print_quality_table(const MeshVariant &variant, const std::vector<LevelResult> &results)
{
    std::cout << "\nCell-quality range - " << variant.name << "\n\n"
              << std::left << std::setw(12) << "delta" << std::setw(24) << "minimum_quality" << std::setw(24)
              << "maximum_quality" << '\n'
              << std::string(60, '-') << '\n';
    for (const LevelResult &result : results)
    {
        std::cout << std::fixed << std::setprecision(4) << std::setw(12) << result.delta << std::scientific
                  << std::setprecision(12) << std::setw(24) << result.minimum_cell_quality << std::setw(24)
                  << result.maximum_cell_quality << '\n';
    }
}

void print_taylor_table(const MeshVariant &variant, const std::vector<LevelResult> &results)
{
    std::cout << "\nInterior leading-Taylor diagnostic - " << variant.name << "\n\n"
              << std::left << std::setw(10) << "delta" << std::setw(24) << "category" << std::setw(10) << "cells"
              << std::setw(17) << "RMS_actual" << std::setw(10) << "p_actual" << std::setw(17) << "RMS_leading"
              << std::setw(10) << "p_leading" << std::setw(17) << "RMS_remainder" << std::setw(12) << "p_remainder"
              << '\n'
              << std::string(128, '-') << '\n';

    for (std::size_t level_index = 0; level_index < results.size(); ++level_index)
    {
        const LevelResult &result{results[level_index]};
        const std::array categories{
            std::pair{"regular_interior", &result.regular_interior_taylor},
            std::pair{"transition_interior", &result.transition_interior_taylor},
        };
        for (const auto &[name, taylor] : categories)
        {
            std::cout << std::fixed << std::setprecision(4) << std::setw(10) << result.delta << std::setw(24) << name
                      << std::setw(10) << taylor->actual.cell_count;
            print_error_value(taylor->actual, taylor->actual.area_weighted_rms_error, 17);
            print_order(taylor->actual_rms_order, level_index == 0, 10);
            print_error_value(taylor->leading, taylor->leading.area_weighted_rms_error, 17);
            print_order(taylor->leading_rms_order, level_index == 0, 10);
            print_error_value(taylor->remainder, taylor->remainder.area_weighted_rms_error, 17);
            print_order(taylor->remainder_rms_order, level_index == 0, 12);
            std::cout << '\n';
        }
    }
}

void print_normalized_scaling(const MeshVariant &variant, const std::vector<LevelResult> &results)
{
    std::cout << "\nNormalized local-error scaling - " << variant.name << "\n\n"
              << std::left << std::setw(12) << "delta" << std::setw(24) << "regular_RMS/d^2" << std::setw(24)
              << "transition_RMS/d" << std::setw(24) << "transition_lead/d" << '\n'
              << std::string(84, '-') << '\n';
    for (const LevelResult &result : results)
    {
        std::cout << std::fixed << std::setprecision(4) << std::setw(12) << result.delta;
        print_error_value(result.regular_interior_taylor.actual,
                          result.regular_interior_taylor.actual.area_weighted_rms_error / (result.delta * result.delta),
                          24);
        print_error_value(result.transition_interior_taylor.actual,
                          result.transition_interior_taylor.actual.area_weighted_rms_error / result.delta, 24);
        print_error_value(result.transition_interior_taylor.leading,
                          result.transition_interior_taylor.leading.area_weighted_rms_error / result.delta, 24);
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
        "Unsupported arguments.\nUsage: cfd_controlled_quad_transition_convergence [--write-vtu]");
}

} // namespace

int main(const int argument_count, char **arguments)
{
    try
    {
        const std::span argument_view{arguments, static_cast<std::size_t>(argument_count)};
        const bool write_vtu{parse_write_vtu_flag(argument_view)};
        const std::filesystem::path output_directory{"output/verification/controlled_quad_transition"};
        if (write_vtu)
        {
            std::filesystem::create_directories(output_directory);
        }

        std::vector<std::filesystem::path> written_paths;
        for (const MeshVariant &variant : mesh_variants)
        {
            std::vector<LevelResult> results;
            results.reserve(grid_levels.size());
            for (std::size_t level_index = 0; level_index < grid_levels.size(); ++level_index)
            {
                results.push_back(run_level(variant, grid_levels.at(level_index), level_index, write_vtu,
                                            output_directory, written_paths));
            }

            compute_orders(results);
            print_overall_table(variant, results);
            print_category_table(variant, results);
            print_squared_error_contributions(variant, results);
            print_geometry_metrics(variant, results);
            print_quality_table(variant, results);
            print_taylor_table(variant, results);
            print_normalized_scaling(variant, results);
            if (variant.shift_internal_column)
            {
                print_transition_scaling_table(results);
            }
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
