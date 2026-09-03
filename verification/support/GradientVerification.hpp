#pragma once

#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/math/Point2.hpp"
#include "cfd/math/Vector2.hpp"
#include "cfd/mesh/Boundary.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/Types.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace cfd::verification
{

inline constexpr double domain_length{2.0};
inline constexpr double domain_height{1.0};

inline constexpr BoundaryId left_boundary_id{0};
inline constexpr BoundaryId right_boundary_id{1};
inline constexpr BoundaryId bottom_boundary_id{2};
inline constexpr BoundaryId top_boundary_id{3};

/// Hessian of the manufactured scalar field used by the gradient studies.
struct Hessian2
{
    double m00{};
    double m01{};
    double m11{};
};

/// Returns phi = sin(x) + y^2 for the shared manufactured case.
[[nodiscard]]
inline double analytical_phi(const Point2 &point) noexcept
{
    return std::sin(point.x) + point.y * point.y;
}

/// Returns the exact gradient (cos(x), 2y) of the manufactured field.
[[nodiscard]]
inline Vector2 analytical_gradient(const Point2 &point) noexcept
{
    return {std::cos(point.x), 2.0 * point.y};
}

/// Returns the exact Hessian diag(-sin(x), 2) of the manufactured field.
[[nodiscard]]
inline Hessian2 analytical_hessian(const Point2 &point) noexcept
{
    return {
        .m00 = -std::sin(point.x),
        .m01 = 0.0,
        .m11 = 2.0,
    };
}

[[nodiscard]]
inline double dot(const Vector2 &first, const Vector2 &second) noexcept
{
    return first.x * second.x + first.y * second.y;
}

[[nodiscard]]
inline double magnitude(const Vector2 &vector) noexcept
{
    return std::hypot(vector.x, vector.y);
}

[[nodiscard]]
inline Vector2 apply(const Hessian2 &matrix, const Vector2 &vector) noexcept
{
    return {
        matrix.m00 * vector.x + matrix.m01 * vector.y,
        matrix.m01 * vector.x + matrix.m11 * vector.y,
    };
}

/// Area-weighted error statistics for a set of cells.
///
/// The RMS value is the domain-normalized measure
/// `sqrt(sum(A * |error|^2) / sum(A))` over the accumulated cells.
struct ErrorStatistics
{
    Index cell_count{};
    double total_area{};
    double area_weighted_squared_error{};
    double area_weighted_rms_error{};
    double linf_error{};
};

/// Accumulates cellwise errors without changing their summation order.
struct ErrorAccumulator
{
    Index cell_count{};
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
        if (!std::isfinite(total_area) || total_area < 0.0 || !std::isfinite(area_weighted_squared_error) ||
            area_weighted_squared_error < 0.0 || !std::isfinite(linf_error) || linf_error < 0.0)
        {
            throw std::runtime_error("Gradient verification produced invalid error statistics.");
        }

        if (cell_count == 0)
        {
            return {};
        }

        if (!(total_area > 0.0))
        {
            throw std::runtime_error("Gradient verification produced invalid error statistics.");
        }

        return {
            .cell_count = cell_count,
            .total_area = total_area,
            .area_weighted_squared_error = area_weighted_squared_error,
            .area_weighted_rms_error = std::sqrt(area_weighted_squared_error / total_area),
            .linf_error = linf_error,
        };
    }
};

/// Computes p = log(E_coarse/E_fine) / log(h_coarse/h_fine).
///
/// The caller defines the meaning of the refinement lengths. Exact zero error
/// makes the logarithmic ratio undefined and produces `std::nullopt`.
[[nodiscard]]
inline std::optional<double> observed_order(const double coarse_error, const double fine_error,
                                            const double coarse_length, const double fine_length)
{
    if (!std::isfinite(coarse_error) || !std::isfinite(fine_error) || coarse_error < 0.0 || fine_error < 0.0)
    {
        throw std::runtime_error("Cannot compute convergence order from invalid error values.");
    }
    if (!std::isfinite(coarse_length) || !std::isfinite(fine_length) || !(coarse_length > 0.0) || !(fine_length > 0.0))
    {
        throw std::runtime_error("Cannot compute convergence order from invalid refinement lengths.");
    }
    if (coarse_error == 0.0 || fine_error == 0.0)
    {
        return std::nullopt;
    }

    const double logarithmic_length_ratio{std::log(coarse_length / fine_length)};
    if (!std::isfinite(logarithmic_length_ratio) || logarithmic_length_ratio == 0.0)
    {
        throw std::runtime_error("Cannot compute convergence order from an invalid refinement-length ratio.");
    }

    const double order{std::log(coarse_error / fine_error) / logarithmic_length_ratio};
    if (!std::isfinite(order))
    {
        throw std::runtime_error("Computed convergence order is non-finite.");
    }
    return order;
}

[[nodiscard]]
inline double manufactured_outward_normal_derivative(const std::string_view boundary_name)
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

    throw std::runtime_error("Unsupported boundary in manufactured gradient verification: " +
                             std::string{boundary_name});
}

/// Builds exact Neumann data for the shared manufactured gradient case.
[[nodiscard]]
inline ScalarBoundaryConditions make_manufactured_neumann_boundary_conditions(const Mesh &mesh)
{
    std::vector<ScalarBoundaryCondition> conditions;
    conditions.reserve(mesh.boundary_groups().size());

    // Phi is supplied analytically; these studies reconstruct its gradient
    // rather than solving a PDE, so pure Neumann data introduce no nullspace.
    for (const BoundaryGroup &group : mesh.boundary_groups())
    {
        if (group.id != conditions.size())
        {
            throw std::runtime_error("Manufactured-gradient boundary groups are not indexed contiguously.");
        }
        conditions.emplace_back(ScalarBoundaryConditionType::Neumann,
                                manufactured_outward_normal_derivative(group.name));
    }

    return {mesh.boundary_groups().size(), std::move(conditions)};
}

/// Symmetry metrics for the best partition of four neighbors into two pairs.
struct OppositePairMetrics
{
    bool valid{};
    double angular_defect{};
    double distance_imbalance{};
};

namespace detail
{

struct PairMetrics
{
    double angular_defect{};
    double distance_imbalance{};
};

[[nodiscard]]
inline PairMetrics compute_pair_metrics(const Vector2 &first, const Vector2 &second)
{
    const double first_length{magnitude(first)};
    const double second_length{magnitude(second)};
    if (!std::isfinite(first_length) || !std::isfinite(second_length) || !(first_length > 0.0) ||
        !(second_length > 0.0))
    {
        throw std::runtime_error("Opposite-pair diagnostics encountered an invalid neighbor displacement.");
    }

    const double cosine{std::clamp(dot(first, second) / (first_length * second_length), -1.0, 1.0)};
    return {
        .angular_defect = 0.5 * (1.0 + cosine),
        .distance_imbalance = std::abs(first_length - second_length) / (first_length + second_length),
    };
}

} // namespace detail

/// Selects the partition minimizing the worse angular defect of its two pairs.
/// Both returned metrics are maxima over the selected partition.
[[nodiscard]]
inline OppositePairMetrics compute_opposite_pair_metrics(const std::array<Vector2, 4> &displacements)
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
        const detail::PairMetrics first_pair{
            detail::compute_pair_metrics(displacements.at(partition.at(0)), displacements.at(partition.at(1)))};
        const detail::PairMetrics second_pair{
            detail::compute_pair_metrics(displacements.at(partition.at(2)), displacements.at(partition.at(3)))};
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

} // namespace cfd::verification
