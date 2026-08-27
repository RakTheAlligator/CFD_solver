#pragma once

namespace cfd
{

/// Dimensions of an axis-aligned rectangular two-dimensional domain.
///
/// Length and height are expressed in metres.
struct RectangleGeometry
{
    double length{};
    double height{};
};

} // namespace cfd