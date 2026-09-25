#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>

namespace cfd::verification
{

/// Returns the analytical face average of a developed parabolic inlet profile.
///
/// The face is bounded by the two physical wall-normal coordinates `first_y`
/// and `second_y`. The profile spans `step_height <= y <= step_height +
/// inlet_height` and has the specified cross-sectional mean velocity.
[[nodiscard]]
inline double parabolic_inlet_face_average(const double first_y, const double second_y, const double step_height,
                                           const double inlet_height, const double mean_velocity)
{
    if (!std::isfinite(first_y) || !std::isfinite(second_y) || !std::isfinite(step_height) ||
        !std::isfinite(inlet_height) || !std::isfinite(mean_velocity) || !(inlet_height > 0.0))
    {
        throw std::invalid_argument("Parabolic inlet face-average inputs must be finite with positive inlet height.");
    }

    const auto [lower_y, upper_y]{std::minmax(first_y, second_y)};
    if (!(upper_y > lower_y))
    {
        throw std::invalid_argument("Parabolic inlet face average requires a positive face length.");
    }
    const double lower_eta{(lower_y - step_height) / inlet_height};
    const double upper_eta{(upper_y - step_height) / inlet_height};
    const auto antiderivative = [](const double eta) { return 0.5 * eta * eta - eta * eta * eta / 3.0; };
    return 6.0 * mean_velocity * (antiderivative(upper_eta) - antiderivative(lower_eta)) / (upper_eta - lower_eta);
}

/// Cell-centered streamwise velocity sampled beside the downstream lower wall.
struct NearWallVelocitySample
{
    double x{};
    double u{};
};

/// Estimates the primary reattachment position from ordered near-wall samples.
///
/// The first negative-to-positive zero crossing downstream of a negative sample
/// is linearly interpolated. The positive sign must persist for the following
/// sample when one is available, rejecting an isolated sign reversal.
[[nodiscard]]
inline double estimate_primary_reattachment_x(const std::span<const NearWallVelocitySample> samples)
{
    if (samples.size() < 2)
    {
        throw std::runtime_error("Reattachment estimation requires at least two near-wall samples.");
    }

    bool found_negative{};
    for (std::size_t sample_id = 0; sample_id < samples.size(); ++sample_id)
    {
        const NearWallVelocitySample &sample{samples[sample_id]};
        if (!std::isfinite(sample.x) || !std::isfinite(sample.u) ||
            (sample_id > 0 && !(sample.x > samples[sample_id - 1].x)))
        {
            throw std::runtime_error("Reattachment samples must be finite and strictly ordered in x.");
        }
        if (sample.u < 0.0)
        {
            found_negative = true;
            continue;
        }
        if (!found_negative || sample_id == 0)
        {
            continue;
        }

        const NearWallVelocitySample &negative_sample{samples[sample_id - 1]};
        if (!(negative_sample.u < 0.0))
        {
            continue;
        }
        if (sample_id + 1 < samples.size() && !(samples[sample_id + 1].u > 0.0))
        {
            continue;
        }

        if (sample.u == 0.0)
        {
            return sample.x;
        }
        const double interpolation_fraction{-negative_sample.u / (sample.u - negative_sample.u)};
        return negative_sample.x + interpolation_fraction * (sample.x - negative_sample.x);
    }

    if (!found_negative)
    {
        throw std::runtime_error("No negative near-wall velocity was found in the primary recirculation region.");
    }
    throw std::runtime_error("No durable primary reattachment was found before the outlet.");
}

} // namespace cfd::verification
