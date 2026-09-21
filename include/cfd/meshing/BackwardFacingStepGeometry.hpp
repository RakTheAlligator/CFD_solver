#pragma once

namespace cfd
{

/// Dimensions of an axis-aligned two-dimensional backward-facing step.
///
/// All dimensions are expressed in metres. The upstream channel occupies
/// `[0, upstream_length] x [step_height, channel_height]`; the downstream
/// channel extends from `upstream_length` to the total streamwise length.
struct BackwardFacingStepGeometry
{
    double upstream_length{};
    double downstream_length{};
    double channel_height{};
    double step_height{};
};

} // namespace cfd
