#pragma once

namespace cfd
{

/// Two-dimensional Cartesian vector.
///
/// The physical meaning and units of the components depend on the quantity
/// represented by the vector.
struct Vector2
{
    double x{};
    double y{};
};

} // namespace cfd