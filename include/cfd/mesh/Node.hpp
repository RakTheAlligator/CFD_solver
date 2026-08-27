#pragma once

namespace cfd
{

/// Position of a mesh node in two-dimensional Cartesian space.
///
/// Coordinates are expressed in metres.
struct Node
{
    double x{};
    double y{};
};

} // namespace cfd