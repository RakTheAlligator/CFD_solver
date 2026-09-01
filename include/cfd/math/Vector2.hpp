#pragma once

namespace cfd
{

/// Two-dimensional mathematical vector.
///
/// The physical meaning and units of the components depend on the quantity
/// represented by the vector.
struct Vector2
{
    double x{};
    double y{};
};

} // namespace cfd
