#pragma once

#include "cfd/math/Point2.hpp"

namespace cfd
{

/// Mesh node represented by its position in two-dimensional Cartesian space.
///
/// Coordinates are expressed in metres. A node's identity is its zero-based
/// index in mesh node storage.
using Node = Point2;

} // namespace cfd
