#pragma once

#include "cfd/mesh/Types.hpp"

#include <algorithm>
#include <cmath>
#include <optional>
#include <stdexcept>

namespace cfd::verification
{

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
            throw std::runtime_error("Verification produced invalid error statistics.");
        }

        if (cell_count == 0)
        {
            return {};
        }

        if (!(total_area > 0.0))
        {
            throw std::runtime_error("Verification produced invalid error statistics.");
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

} // namespace cfd::verification
