#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVectorField.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/io/VtkWriter.hpp"
#include "cfd/math/Point2.hpp"
#include "cfd/math/Vector2.hpp"
#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/meshing/GmshMesher.hpp"
#include "cfd/meshing/RawMeshData.hpp"
#include "cfd/meshing/RectangleGeometry.hpp"
#include "cfd/numerics/LeastSquaresGradient.hpp"

#include "support/GradientVerification.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
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

constexpr double boundary_coordinate_tolerance{1.0e-12};

using cfd::verification::ErrorStatistics;
using ErrorStatisticsAccumulator = cfd::verification::ErrorAccumulator;
using cfd::verification::analytical_gradient;
using cfd::verification::analytical_hessian;
using cfd::verification::analytical_phi;
using cfd::verification::apply;
using cfd::verification::bottom_boundary_id;
using cfd::verification::compute_opposite_pair_metrics;
using cfd::verification::domain_height;
using cfd::verification::domain_length;
using cfd::verification::dot;
using cfd::verification::Hessian2;
using cfd::verification::left_boundary_id;
using cfd::verification::magnitude;
using cfd::verification::make_manufactured_neumann_boundary_conditions;
using cfd::verification::observed_order;
using cfd::verification::OppositePairMetrics;
using cfd::verification::right_boundary_id;
using cfd::verification::top_boundary_id;

constexpr std::array target_mesh_sizes{0.20, 0.10, 0.05, 0.025};

struct MeshFamily
{
    cfd::CellType cell_type;
    std::string_view table_name;
    std::string_view file_prefix;
};

constexpr std::array mesh_families{
    MeshFamily{cfd::CellType::Triangle, "triangles", "triangle"},
    MeshFamily{cfd::CellType::Quadrilateral, "quadrilaterals", "quad"},
};

struct PearsonAccumulator
{
    std::size_t sample_count{};
    double mean_x{};
    double mean_y{};
    double sum_squared_x{};
    double sum_squared_y{};
    double sum_cross{};

    void add(const double x_value, const double y_value) noexcept
    {
        ++sample_count;
        const double sample_count_as_double{static_cast<double>(sample_count)};
        const double delta_x{x_value - mean_x};
        const double delta_y{y_value - mean_y};

        mean_x += delta_x / sample_count_as_double;
        mean_y += delta_y / sample_count_as_double;
        sum_squared_x += delta_x * (x_value - mean_x);
        sum_squared_y += delta_y * (y_value - mean_y);
        sum_cross += delta_x * (y_value - mean_y);
    }

    [[nodiscard]]
    std::optional<double> correlation() const noexcept
    {
        if (sample_count < 2 || !(sum_squared_x > 0.0) || !(sum_squared_y > 0.0))
        {
            return std::nullopt;
        }

        const double denominator{std::sqrt(sum_squared_x * sum_squared_y)};
        if (!std::isfinite(denominator) || !(denominator > 0.0))
        {
            return std::nullopt;
        }

        const double result{sum_cross / denominator};
        if (!std::isfinite(result))
        {
            return std::nullopt;
        }

        return std::clamp(result, -1.0, 1.0);
    }
};

struct QuartileStatistics
{
    ErrorStatistics errors;
    double minimum_distance_imbalance{};
    double maximum_distance_imbalance{};
};

struct QuadLevelDiagnostics
{
    ErrorStatistics boundary_cells;
    ErrorStatistics interior_cells;
    ErrorStatistics eligible_cells;
    double mean_angular_defect{};
    double maximum_angular_defect{};
    double mean_distance_imbalance{};
    double maximum_distance_imbalance{};
    std::optional<double> log_error_angular_correlation;
    std::optional<double> log_error_distance_correlation;
    std::optional<double> log_error_cell_quality_correlation;
    std::array<QuartileStatistics, 4> distance_imbalance_quartiles;
};

struct LevelResult
{
    double target_mesh_size{};
    cfd::Index cell_count{};
    double characteristic_mesh_length{};
    double l2_error{};
    double linf_error{};
    double max_prediction_mismatch{};
    std::optional<double> l2_order;
    std::optional<double> linf_order;
    std::optional<QuadLevelDiagnostics> quad_diagnostics;
};

struct LocalWlsAssembly
{
    double m00{};
    double m01{};
    double m11{};
    cfd::Vector2 exact_defect_forcing{};
    cfd::Vector2 leading_defect_forcing{};
};

struct LocalWlsDiagnostics
{
    cfd::Vector2 predicted_error{};
    cfd::Vector2 leading_predicted_error{};
    double wls_det_over_trace_squared{};
    bool is_boundary_cell{};
    OppositePairMetrics opposite_pair_metrics;
};

struct DiagnosticFields
{
    explicit DiagnosticFields(const cfd::Index cell_count)
        : predicted_error(cell_count), predicted_error_magnitude(cell_count), prediction_mismatch(cell_count),
          leading_error(cell_count), leading_error_magnitude(cell_count), leading_error_difference(cell_count),
          wls_det_over_trace_squared(cell_count), is_boundary_cell(cell_count), opposite_pair_metrics_valid(cell_count),
          opposite_pair_angular_defect(cell_count), opposite_pair_distance_imbalance(cell_count)
    {
    }

    cfd::CellVectorField predicted_error;
    cfd::CellScalarField predicted_error_magnitude;
    cfd::CellScalarField prediction_mismatch;
    cfd::CellVectorField leading_error;
    cfd::CellScalarField leading_error_magnitude;
    cfd::CellScalarField leading_error_difference;
    cfd::CellScalarField wls_det_over_trace_squared;
    cfd::CellScalarField is_boundary_cell;
    cfd::CellScalarField opposite_pair_metrics_valid;
    cfd::CellScalarField opposite_pair_angular_defect;
    cfd::CellScalarField opposite_pair_distance_imbalance;
};

struct EligibleQuadSample
{
    double distance_imbalance{};
    double cell_area{};
    double error_magnitude{};
};

[[nodiscard]]
bool matches_coordinate(const double first, const double second, const double target) noexcept
{
    return std::abs(first - target) <= boundary_coordinate_tolerance &&
           std::abs(second - target) <= boundary_coordinate_tolerance;
}

void split_rectangle_boundary_groups(cfd::RawMeshData &raw_mesh)
{
    // The production rectangle mesher combines top and bottom as `wall`.
    // This analytical case needs different exact Neumann data on those sides,
    // so the verification study relabels only its generated raw boundary edges.
    raw_mesh.boundary_groups = {
        {left_boundary_id, "left"},
        {right_boundary_id, "right"},
        {bottom_boundary_id, "bottom"},
        {top_boundary_id, "top"},
    };

    std::array<cfd::Index, 4> boundary_edge_counts{};

    for (cfd::BoundaryEdge &edge : raw_mesh.boundary_edges)
    {
        const cfd::Point2 &first{raw_mesh.nodes[edge.node_ids[0]]};
        const cfd::Point2 &second{raw_mesh.nodes[edge.node_ids[1]]};

        if (matches_coordinate(first.x, second.x, 0.0))
        {
            edge.boundary_id = left_boundary_id;
        }
        else if (matches_coordinate(first.x, second.x, domain_length))
        {
            edge.boundary_id = right_boundary_id;
        }
        else if (matches_coordinate(first.y, second.y, 0.0))
        {
            edge.boundary_id = bottom_boundary_id;
        }
        else if (matches_coordinate(first.y, second.y, domain_height))
        {
            edge.boundary_id = top_boundary_id;
        }
        else
        {
            throw std::runtime_error("Unable to classify a generated rectangle boundary edge.");
        }

        ++boundary_edge_counts.at(edge.boundary_id);
    }

    if (std::any_of(boundary_edge_counts.begin(), boundary_edge_counts.end(),
                    [](const cfd::Index count) { return count == 0; }))
    {
        throw std::runtime_error("Generated rectangle does not contain all four verification boundaries.");
    }
}

void add_observation(LocalWlsAssembly &assembly, const cfd::Vector2 &direction, const double exact_defect,
                     const double leading_defect) noexcept
{
    assembly.m00 += direction.x * direction.x;
    assembly.m01 += direction.x * direction.y;
    assembly.m11 += direction.y * direction.y;
    assembly.exact_defect_forcing.x += direction.x * exact_defect;
    assembly.exact_defect_forcing.y += direction.y * exact_defect;
    assembly.leading_defect_forcing.x += direction.x * leading_defect;
    assembly.leading_defect_forcing.y += direction.y * leading_defect;
}

[[nodiscard]]
cfd::Vector2 solve(const LocalWlsAssembly &assembly, const cfd::Vector2 &right_hand_side)
{
    constexpr double relative_singularity_tolerance{64.0 * std::numeric_limits<double>::epsilon()};

    const double determinant{assembly.m00 * assembly.m11 - assembly.m01 * assembly.m01};
    const double trace{assembly.m00 + assembly.m11};
    const double determinant_tolerance{relative_singularity_tolerance * trace * trace};

    if (!(determinant > determinant_tolerance))
    {
        throw std::runtime_error("Diagnostic WLS assembly produced a singular stencil.");
    }

    return {
        (right_hand_side.x * assembly.m11 - right_hand_side.y * assembly.m01) / determinant,
        (assembly.m00 * right_hand_side.y - assembly.m01 * right_hand_side.x) / determinant,
    };
}

[[nodiscard]]
LocalWlsDiagnostics compute_local_wls_diagnostics(const cfd::Mesh &mesh, const cfd::CellScalarField &phi,
                                                  const cfd::ScalarBoundaryConditions &boundary_conditions,
                                                  const cfd::Index cell_id, const cfd::Vector2 &exact_gradient_at_cell)
{
    const auto cell_offsets{mesh.cell_node_offsets()};
    const auto cell_faces{mesh.cell_faces()};
    const auto face_adjacencies{mesh.face_adjacencies()};
    const auto face_boundary_ids{mesh.face_boundary_ids()};
    const auto cell_centers{mesh.cell_centers()};
    const auto face_centers{mesh.face_centers()};
    const auto face_lengths{mesh.face_lengths()};
    const auto face_area_vectors{mesh.face_area_vectors()};
    const auto field_values{phi.values()};

    const cfd::Point2 &cell_center{cell_centers[cell_id]};
    const Hessian2 hessian{analytical_hessian(cell_center)};
    LocalWlsAssembly assembly;
    std::array<cfd::Vector2, 4> internal_displacements{};
    std::size_t internal_neighbor_count{};
    bool is_boundary_cell{};

    const cfd::Index cell_face_begin_offset{cell_offsets[cell_id]};
    const cfd::Index cell_face_end_offset{cell_offsets[cell_id + 1]};

    for (cfd::Index cell_face_position = cell_face_begin_offset; cell_face_position < cell_face_end_offset;
         ++cell_face_position)
    {
        const cfd::Index face_id{cell_faces[cell_face_position]};
        const cfd::FaceAdjacency &adjacency{face_adjacencies[face_id]};

        if (!adjacency.is_boundary())
        {
            const cfd::Index other_cell_id{adjacency.owner == cell_id ? adjacency.neighbor : adjacency.owner};
            const cfd::Vector2 displacement{
                cell_centers[other_cell_id].x - cell_center.x,
                cell_centers[other_cell_id].y - cell_center.y,
            };
            const double distance{magnitude(displacement)};
            if (!std::isfinite(distance) || !(distance > 0.0))
            {
                throw std::runtime_error("WLS diagnostics encountered coincident cell centers.");
            }

            const cfd::Vector2 direction{displacement.x / distance, displacement.y / distance};
            const double normalized_difference{(field_values[other_cell_id] - field_values[cell_id]) / distance};
            const double exact_defect{normalized_difference - dot(direction, exact_gradient_at_cell)};
            const double leading_defect{0.5 * distance * dot(direction, apply(hessian, direction))};
            add_observation(assembly, direction, exact_defect, leading_defect);

            if (internal_neighbor_count >= internal_displacements.size())
            {
                throw std::runtime_error("WLS diagnostics encountered more than four internal neighbors.");
            }
            internal_displacements.at(internal_neighbor_count) = displacement;
            ++internal_neighbor_count;
            continue;
        }

        is_boundary_cell = true;
        const cfd::ScalarBoundaryCondition &condition{boundary_conditions[face_boundary_ids[face_id]]};
        if (condition.type != cfd::ScalarBoundaryConditionType::Neumann)
        {
            throw std::runtime_error("This verification diagnostic expects exact Neumann boundary observations.");
        }

        const double face_length{face_lengths[face_id]};
        if (!std::isfinite(face_length) || !(face_length > 0.0))
        {
            throw std::runtime_error("WLS diagnostics encountered an invalid boundary-face length.");
        }

        const cfd::Vector2 &area_vector{face_area_vectors[face_id]};
        const cfd::Vector2 outward_normal{area_vector.x / face_length, area_vector.y / face_length};
        const double exact_defect{condition.value - dot(outward_normal, exact_gradient_at_cell)};
        const cfd::Vector2 boundary_displacement{
            face_centers[face_id].x - cell_center.x,
            face_centers[face_id].y - cell_center.y,
        };
        const double leading_defect{dot(outward_normal, apply(hessian, boundary_displacement))};
        add_observation(assembly, outward_normal, exact_defect, leading_defect);
    }

    const double determinant{assembly.m00 * assembly.m11 - assembly.m01 * assembly.m01};
    const double trace{assembly.m00 + assembly.m11};
    const double trace_squared{trace * trace};
    if (!std::isfinite(determinant) || !std::isfinite(trace_squared) || !(trace > 0.0) || !(trace_squared > 0.0))
    {
        throw std::runtime_error("WLS diagnostics produced an invalid determinant/trace indicator.");
    }

    OppositePairMetrics opposite_pair_metrics;
    if (mesh.cell_types()[cell_id] == cfd::CellType::Quadrilateral && !is_boundary_cell &&
        internal_neighbor_count == internal_displacements.size())
    {
        opposite_pair_metrics = compute_opposite_pair_metrics(internal_displacements);
    }

    return {
        .predicted_error = solve(assembly, assembly.exact_defect_forcing),
        .leading_predicted_error = solve(assembly, assembly.leading_defect_forcing),
        .wls_det_over_trace_squared = determinant / trace_squared,
        .is_boundary_cell = is_boundary_cell,
        .opposite_pair_metrics = opposite_pair_metrics,
    };
}

[[nodiscard]]
std::array<QuartileStatistics, 4> make_distance_imbalance_quartiles(std::vector<EligibleQuadSample> samples)
{
    std::ranges::sort(samples, {}, &EligibleQuadSample::distance_imbalance);

    std::array<QuartileStatistics, 4> quartiles{};
    for (std::size_t quartile_index = 0; quartile_index < quartiles.size(); ++quartile_index)
    {
        const std::size_t begin_index{quartile_index * samples.size() / quartiles.size()};
        const std::size_t end_index{(quartile_index + 1) * samples.size() / quartiles.size()};
        if (begin_index == end_index)
        {
            continue;
        }

        ErrorStatisticsAccumulator error_accumulator;
        for (std::size_t sample_index = begin_index; sample_index < end_index; ++sample_index)
        {
            const EligibleQuadSample &sample{samples.at(sample_index)};
            error_accumulator.add(sample.cell_area, sample.error_magnitude);
        }

        quartiles.at(quartile_index) = {
            .errors = error_accumulator.finish(),
            .minimum_distance_imbalance = samples.at(begin_index).distance_imbalance,
            .maximum_distance_imbalance = samples.at(end_index - 1).distance_imbalance,
        };
    }

    return quartiles;
}

void write_verification_vtu(const cfd::Mesh &mesh, const cfd::CellScalarField &phi,
                            const cfd::CellVectorField &numerical_gradient, const cfd::CellVectorField &exact_gradient,
                            const cfd::CellVectorField &gradient_error, const cfd::CellScalarField &error_magnitude,
                            const DiagnosticFields &diagnostics, const std::filesystem::path &output_path)
{
    const std::array scalar_fields{
        cfd::VtkCellScalarData{"phi", phi.values()},
        cfd::VtkCellScalarData{"grad_phi_error_magnitude", error_magnitude.values()},
        cfd::VtkCellScalarData{"grad_phi_error_predicted_magnitude", diagnostics.predicted_error_magnitude.values()},
        cfd::VtkCellScalarData{"grad_phi_error_prediction_mismatch", diagnostics.prediction_mismatch.values()},
        cfd::VtkCellScalarData{"grad_phi_error_leading_magnitude", diagnostics.leading_error_magnitude.values()},
        cfd::VtkCellScalarData{"grad_phi_error_leading_difference", diagnostics.leading_error_difference.values()},
        cfd::VtkCellScalarData{"wls_det_over_trace_squared", diagnostics.wls_det_over_trace_squared.values()},
        cfd::VtkCellScalarData{"is_boundary_cell", diagnostics.is_boundary_cell.values()},
        cfd::VtkCellScalarData{"opposite_pair_metrics_valid", diagnostics.opposite_pair_metrics_valid.values()},
        cfd::VtkCellScalarData{"opposite_pair_angular_defect", diagnostics.opposite_pair_angular_defect.values()},
        cfd::VtkCellScalarData{"opposite_pair_distance_imbalance",
                               diagnostics.opposite_pair_distance_imbalance.values()},
    };
    const std::array vector_fields{
        cfd::VtkCellVectorData{"grad_phi", numerical_gradient.values()},
        cfd::VtkCellVectorData{"grad_phi_exact", exact_gradient.values()},
        cfd::VtkCellVectorData{"grad_phi_error", gradient_error.values()},
        cfd::VtkCellVectorData{"grad_phi_error_predicted", diagnostics.predicted_error.values()},
        cfd::VtkCellVectorData{"grad_phi_error_leading", diagnostics.leading_error.values()},
    };
    const cfd::VtkCellData cell_data{
        .scalars = scalar_fields,
        .vectors = vector_fields,
    };

    cfd::write_vtu(mesh, output_path, cell_data);
}

[[nodiscard]]
LevelResult run_level(const MeshFamily &family, const double target_mesh_size, const std::size_t level_index,
                      const bool write_vtu, const std::filesystem::path &output_directory,
                      std::vector<std::filesystem::path> &written_paths)
{
    const cfd::RectangleGeometry geometry{
        .length = domain_length,
        .height = domain_height,
    };
    const cfd::MeshGenerationOptions options{
        .mesh_size = target_mesh_size,
        .cell_type = family.cell_type,
    };

    cfd::RawMeshData raw_mesh{cfd::generate_mesh(geometry, options)};
    split_rectangle_boundary_groups(raw_mesh);

    cfd::MeshBuildResult build_result{cfd::build_mesh(std::move(raw_mesh))};
    const cfd::Mesh &mesh{build_result.mesh};

    if (mesh.cell_count() == 0)
    {
        throw std::runtime_error("Gradient verification requires a non-empty mesh.");
    }

    cfd::CellScalarField phi{mesh.cell_count()};
    cfd::CellVectorField exact_gradient{mesh.cell_count()};

    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        phi[cell_id] = analytical_phi(mesh.cell_centers()[cell_id]);
        exact_gradient[cell_id] = analytical_gradient(mesh.cell_centers()[cell_id]);
    }

    const cfd::ScalarBoundaryConditions boundary_conditions{make_manufactured_neumann_boundary_conditions(mesh)};
    cfd::CellVectorField numerical_gradient{mesh.cell_count()};
    cfd::compute_least_squares_gradient(mesh, phi, boundary_conditions, numerical_gradient);

    cfd::CellVectorField gradient_error{mesh.cell_count()};
    cfd::CellScalarField error_magnitude{mesh.cell_count()};
    DiagnosticFields diagnostics{mesh.cell_count()};

    double total_cell_area{};
    double area_weighted_squared_error{};
    double linf_error{};
    double max_prediction_mismatch{};
    double max_predicted_error_magnitude{};

    ErrorStatisticsAccumulator boundary_error_accumulator;
    ErrorStatisticsAccumulator interior_error_accumulator;
    ErrorStatisticsAccumulator eligible_error_accumulator;
    PearsonAccumulator angular_correlation;
    PearsonAccumulator distance_correlation;
    PearsonAccumulator cell_quality_correlation;
    std::vector<EligibleQuadSample> eligible_samples;
    double angular_defect_sum{};
    double maximum_angular_defect{};
    double distance_imbalance_sum{};
    double maximum_distance_imbalance{};

    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        const double cell_area{mesh.cell_areas()[cell_id]};
        if (!std::isfinite(cell_area) || !(cell_area > 0.0))
        {
            throw std::runtime_error("Gradient verification encountered an invalid cell area.");
        }

        gradient_error[cell_id] = {
            numerical_gradient[cell_id].x - exact_gradient[cell_id].x,
            numerical_gradient[cell_id].y - exact_gradient[cell_id].y,
        };

        const double squared_error{gradient_error[cell_id].x * gradient_error[cell_id].x +
                                   gradient_error[cell_id].y * gradient_error[cell_id].y};
        if (!std::isfinite(squared_error) || squared_error < 0.0)
        {
            throw std::runtime_error("Gradient verification produced a non-finite error.");
        }

        error_magnitude[cell_id] = std::sqrt(squared_error);
        const LocalWlsDiagnostics local_diagnostics{
            compute_local_wls_diagnostics(mesh, phi, boundary_conditions, cell_id, exact_gradient[cell_id])};
        diagnostics.predicted_error[cell_id] = local_diagnostics.predicted_error;
        diagnostics.predicted_error_magnitude[cell_id] = magnitude(local_diagnostics.predicted_error);
        diagnostics.leading_error[cell_id] = local_diagnostics.leading_predicted_error;
        diagnostics.leading_error_magnitude[cell_id] = magnitude(local_diagnostics.leading_predicted_error);
        diagnostics.wls_det_over_trace_squared[cell_id] = local_diagnostics.wls_det_over_trace_squared;
        diagnostics.is_boundary_cell[cell_id] = local_diagnostics.is_boundary_cell ? 1.0 : 0.0;
        diagnostics.opposite_pair_metrics_valid[cell_id] = local_diagnostics.opposite_pair_metrics.valid ? 1.0 : 0.0;
        diagnostics.opposite_pair_angular_defect[cell_id] = local_diagnostics.opposite_pair_metrics.angular_defect;
        diagnostics.opposite_pair_distance_imbalance[cell_id] =
            local_diagnostics.opposite_pair_metrics.distance_imbalance;

        const cfd::Vector2 prediction_difference{
            local_diagnostics.predicted_error.x - gradient_error[cell_id].x,
            local_diagnostics.predicted_error.y - gradient_error[cell_id].y,
        };
        diagnostics.prediction_mismatch[cell_id] = magnitude(prediction_difference);

        const cfd::Vector2 leading_difference{
            gradient_error[cell_id].x - local_diagnostics.leading_predicted_error.x,
            gradient_error[cell_id].y - local_diagnostics.leading_predicted_error.y,
        };
        diagnostics.leading_error_difference[cell_id] = magnitude(leading_difference);

        max_prediction_mismatch = std::max(max_prediction_mismatch, diagnostics.prediction_mismatch[cell_id]);
        max_predicted_error_magnitude =
            std::max(max_predicted_error_magnitude, diagnostics.predicted_error_magnitude[cell_id]);
        total_cell_area += cell_area;

        // Normalization by total_cell_area below makes this an area-weighted,
        // domain-normalized (RMS) L2 error measure.
        area_weighted_squared_error += cell_area * squared_error;
        linf_error = std::max(linf_error, error_magnitude[cell_id]);

        if (family.cell_type == cfd::CellType::Quadrilateral)
        {
            if (local_diagnostics.is_boundary_cell)
            {
                boundary_error_accumulator.add(cell_area, error_magnitude[cell_id]);
            }
            else
            {
                interior_error_accumulator.add(cell_area, error_magnitude[cell_id]);
            }

            if (local_diagnostics.opposite_pair_metrics.valid)
            {
                const double angular_defect{local_diagnostics.opposite_pair_metrics.angular_defect};
                const double distance_imbalance{local_diagnostics.opposite_pair_metrics.distance_imbalance};
                eligible_error_accumulator.add(cell_area, error_magnitude[cell_id]);
                angular_defect_sum += angular_defect;
                maximum_angular_defect = std::max(maximum_angular_defect, angular_defect);
                distance_imbalance_sum += distance_imbalance;
                maximum_distance_imbalance = std::max(maximum_distance_imbalance, distance_imbalance);
                eligible_samples.push_back({
                    .distance_imbalance = distance_imbalance,
                    .cell_area = cell_area,
                    .error_magnitude = error_magnitude[cell_id],
                });

                if (error_magnitude[cell_id] > 0.0)
                {
                    const double logarithmic_error{std::log(error_magnitude[cell_id])};
                    angular_correlation.add(logarithmic_error, angular_defect);
                    distance_correlation.add(logarithmic_error, distance_imbalance);
                    cell_quality_correlation.add(logarithmic_error, mesh.cell_qualities()[cell_id]);
                }
            }
        }
    }

    if (!std::isfinite(total_cell_area) || !(total_cell_area > 0.0))
    {
        throw std::runtime_error("Gradient verification produced an invalid total cell area.");
    }

    if (!std::isfinite(area_weighted_squared_error) || area_weighted_squared_error < 0.0 || !std::isfinite(linf_error))
    {
        throw std::runtime_error("Gradient verification produced an invalid error norm.");
    }

    const double l2_error{std::sqrt(area_weighted_squared_error / total_cell_area)};

    // Observed order uses this measured geometric length because a requested
    // Gmsh target size is not the exact spacing of the generated cells.
    const double characteristic_mesh_length{std::sqrt(total_cell_area / static_cast<double>(mesh.cell_count()))};

    if (!std::isfinite(l2_error) || !std::isfinite(characteristic_mesh_length) || !(characteristic_mesh_length > 0.0))
    {
        throw std::runtime_error("Gradient verification produced an invalid L2 error or mesh length.");
    }

    const double identity_scale{std::max({1.0, linf_error, max_predicted_error_magnitude})};
    const double identity_tolerance{1024.0 * std::numeric_limits<double>::epsilon() * identity_scale};
    if (!std::isfinite(max_prediction_mismatch) || max_prediction_mismatch > identity_tolerance)
    {
        std::ostringstream message;
        message << "Exact WLS error identity mismatch exceeds roundoff tolerance: mismatch=" << max_prediction_mismatch
                << ", tolerance=" << identity_tolerance << '.';
        throw std::runtime_error(message.str());
    }

    std::optional<QuadLevelDiagnostics> quad_diagnostics;
    if (family.cell_type == cfd::CellType::Quadrilateral)
    {
        const ErrorStatistics eligible_statistics{eligible_error_accumulator.finish()};
        const double eligible_count{static_cast<double>(eligible_statistics.cell_count)};

        quad_diagnostics = QuadLevelDiagnostics{
            .boundary_cells = boundary_error_accumulator.finish(),
            .interior_cells = interior_error_accumulator.finish(),
            .eligible_cells = eligible_statistics,
            .mean_angular_defect = eligible_count > 0.0 ? angular_defect_sum / eligible_count : 0.0,
            .maximum_angular_defect = maximum_angular_defect,
            .mean_distance_imbalance = eligible_count > 0.0 ? distance_imbalance_sum / eligible_count : 0.0,
            .maximum_distance_imbalance = maximum_distance_imbalance,
            .log_error_angular_correlation = angular_correlation.correlation(),
            .log_error_distance_correlation = distance_correlation.correlation(),
            .log_error_cell_quality_correlation = cell_quality_correlation.correlation(),
            .distance_imbalance_quartiles = make_distance_imbalance_quartiles(std::move(eligible_samples)),
        };
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
        .target_mesh_size = target_mesh_size,
        .cell_count = mesh.cell_count(),
        .characteristic_mesh_length = characteristic_mesh_length,
        .l2_error = l2_error,
        .linf_error = linf_error,
        .max_prediction_mismatch = max_prediction_mismatch,
        .l2_order = std::nullopt,
        .linf_order = std::nullopt,
        .quad_diagnostics = quad_diagnostics,
    };
}

void compute_observed_orders(std::vector<LevelResult> &results)
{
    for (std::size_t level_index = 1; level_index < results.size(); ++level_index)
    {
        const LevelResult &coarse{results[level_index - 1]};
        LevelResult &fine{results[level_index]};

        fine.l2_order = observed_order(coarse.l2_error, fine.l2_error, coarse.characteristic_mesh_length,
                                       fine.characteristic_mesh_length);
        fine.linf_order = observed_order(coarse.linf_error, fine.linf_error, coarse.characteristic_mesh_length,
                                         fine.characteristic_mesh_length);
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

void print_convergence_table(const MeshFamily &family, const std::vector<LevelResult> &results)
{
    std::cout << "\nLeast-squares gradient convergence - " << family.table_name << "\n\n"
              << std::left << std::setw(12) << "target_h" << std::setw(12) << "cells" << std::setw(16) << "h_char"
              << std::setw(16) << "L2_error" << std::setw(16) << "Linf_error" << std::setw(22) << "max_pred_mismatch"
              << std::setw(10) << "p_L2" << std::setw(10) << "p_Linf" << '\n'
              << std::string(114, '-') << '\n';

    for (std::size_t level_index = 0; level_index < results.size(); ++level_index)
    {
        const LevelResult &result{results[level_index]};

        std::cout << std::fixed << std::setprecision(3) << std::setw(12) << result.target_mesh_size << std::setw(12)
                  << result.cell_count << std::scientific << std::setprecision(6) << std::setw(16)
                  << result.characteristic_mesh_length << std::setw(16) << result.l2_error << std::setw(16)
                  << result.linf_error << std::setw(22) << result.max_prediction_mismatch;

        print_order(result.l2_order, level_index == 0);
        print_order(result.linf_order, level_index == 0);
        std::cout << '\n';
    }
}

void print_optional_correlation(const std::optional<double> correlation)
{
    if (correlation.has_value())
    {
        std::cout << std::fixed << std::setprecision(4) << std::setw(18) << *correlation;
        return;
    }

    std::cout << std::setw(18) << "n/a";
}

void print_quadrilateral_diagnostics(const std::vector<LevelResult> &results)
{
    std::cout << "\nQuadrilateral error categories\n\n"
              << std::left << std::setw(12) << "target_h" << std::setw(14) << "category" << std::setw(12) << "cells"
              << std::setw(18) << "RMS_error" << std::setw(18) << "Linf_error" << '\n'
              << std::string(74, '-') << '\n';

    for (const LevelResult &result : results)
    {
        if (!result.quad_diagnostics.has_value())
        {
            throw std::runtime_error("Missing quadrilateral diagnostics.");
        }

        const QuadLevelDiagnostics &diagnostics{*result.quad_diagnostics};
        const std::array categories{
            std::pair{"all", ErrorStatistics{.cell_count = result.cell_count,
                                             .area_weighted_rms_error = result.l2_error,
                                             .linf_error = result.linf_error}},
            std::pair{"boundary", diagnostics.boundary_cells},
            std::pair{"interior", diagnostics.interior_cells},
            std::pair{"eligible", diagnostics.eligible_cells},
        };

        for (const auto &[category_name, statistics] : categories)
        {
            std::cout << std::fixed << std::setprecision(3) << std::setw(12) << result.target_mesh_size << std::setw(14)
                      << category_name << std::setw(12) << statistics.cell_count << std::scientific
                      << std::setprecision(6) << std::setw(18) << statistics.area_weighted_rms_error << std::setw(18)
                      << statistics.linf_error << '\n';
        }
    }

    std::cout << "\nEligible interior-quadrilateral stencil asymmetry\n\n"
              << std::left << std::setw(12) << "target_h" << std::setw(12) << "cells" << std::setw(16) << "mean_ang"
              << std::setw(16) << "max_ang" << std::setw(16) << "mean_dist" << std::setw(16) << "max_dist"
              << std::setw(18) << "corr_logE_ang" << std::setw(18) << "corr_logE_dist" << std::setw(18)
              << "corr_logE_quality" << '\n'
              << std::string(142, '-') << '\n';

    for (const LevelResult &result : results)
    {
        const QuadLevelDiagnostics &diagnostics{*result.quad_diagnostics};
        std::cout << std::fixed << std::setprecision(3) << std::setw(12) << result.target_mesh_size << std::setw(12)
                  << diagnostics.eligible_cells.cell_count << std::scientific << std::setprecision(6) << std::setw(16)
                  << diagnostics.mean_angular_defect << std::setw(16) << diagnostics.maximum_angular_defect
                  << std::setw(16) << diagnostics.mean_distance_imbalance << std::setw(16)
                  << diagnostics.maximum_distance_imbalance;
        print_optional_correlation(diagnostics.log_error_angular_correlation);
        print_optional_correlation(diagnostics.log_error_distance_correlation);
        print_optional_correlation(diagnostics.log_error_cell_quality_correlation);
        std::cout << '\n';
    }

    std::cout << "\nDistance-imbalance rank quartiles for eligible quadrilaterals\n\n"
              << std::left << std::setw(12) << "target_h" << std::setw(10) << "quartile" << std::setw(12) << "cells"
              << std::setw(16) << "min_dist" << std::setw(16) << "max_dist" << std::setw(18) << "RMS_error"
              << std::setw(18) << "Linf_error" << '\n'
              << std::string(102, '-') << '\n';

    for (const LevelResult &result : results)
    {
        const QuadLevelDiagnostics &diagnostics{*result.quad_diagnostics};
        for (std::size_t quartile_index = 0; quartile_index < diagnostics.distance_imbalance_quartiles.size();
             ++quartile_index)
        {
            const QuartileStatistics &quartile{diagnostics.distance_imbalance_quartiles.at(quartile_index)};
            std::cout << std::fixed << std::setprecision(3) << std::setw(12) << result.target_mesh_size << std::setw(10)
                      << quartile_index + 1 << std::setw(12) << quartile.errors.cell_count << std::scientific
                      << std::setprecision(6) << std::setw(16) << quartile.minimum_distance_imbalance << std::setw(16)
                      << quartile.maximum_distance_imbalance << std::setw(18) << quartile.errors.area_weighted_rms_error
                      << std::setw(18) << quartile.errors.linf_error << '\n';
        }
    }
}

void print_controlled_pair_cancellation()
{
    constexpr cfd::Point2 representative_point{1.0, 0.5};
    constexpr double inverse_square_root_two{0.70710678118654752440};
    constexpr cfd::Vector2 direction{inverse_square_root_two, inverse_square_root_two};
    constexpr double h{0.1};
    constexpr std::array distance_ratios{1.0, 0.9, 0.75, 0.5};

    const Hessian2 hessian{analytical_hessian(representative_point)};
    const double plus_defect{0.5 * h * dot(direction, apply(hessian, direction))};

    std::cout << "Controlled opposite-pair leading-defect cancellation\n\n"
              << "point=(1, 0.5), direction=(1/sqrt(2), 1/sqrt(2)), h=0.1\n\n"
              << std::left << std::setw(14) << "k_over_h" << std::setw(32) << "leading_pair_forcing_magnitude" << '\n'
              << std::string(46, '-') << '\n';

    for (const double distance_ratio : distance_ratios)
    {
        const double k{distance_ratio * h};
        const cfd::Vector2 opposite_direction{-direction.x, -direction.y};
        const double minus_defect{0.5 * k * dot(opposite_direction, apply(hessian, opposite_direction))};
        const cfd::Vector2 pair_forcing{
            direction.x * plus_defect + opposite_direction.x * minus_defect,
            direction.y * plus_defect + opposite_direction.y * minus_defect,
        };

        std::cout << std::fixed << std::setprecision(2) << std::setw(14) << distance_ratio << std::scientific
                  << std::setprecision(12) << std::setw(32) << magnitude(pair_forcing) << '\n';
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

    throw std::invalid_argument("Unsupported arguments.\nUsage: cfd_least_squares_gradient_convergence [--write-vtu]");
}

} // namespace

int main(const int argument_count, char **arguments)
{
    try
    {
        const std::span argument_view{arguments, static_cast<std::size_t>(argument_count)};
        const bool write_vtu{parse_write_vtu_flag(argument_view)};
        const std::filesystem::path output_directory{"output/verification/least_squares_gradient"};

        if (write_vtu)
        {
            std::filesystem::create_directories(output_directory);
        }

        std::vector<std::filesystem::path> written_paths;
        print_controlled_pair_cancellation();

        for (const MeshFamily &family : mesh_families)
        {
            std::vector<LevelResult> results;
            results.reserve(target_mesh_sizes.size());

            std::size_t level_index{};
            for (const double target_mesh_size : target_mesh_sizes)
            {
                results.push_back(
                    run_level(family, target_mesh_size, level_index, write_vtu, output_directory, written_paths));
                ++level_index;
            }

            compute_observed_orders(results);
            print_convergence_table(family, results);
            if (family.cell_type == cfd::CellType::Quadrilateral)
            {
                print_quadrilateral_diagnostics(results);
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
