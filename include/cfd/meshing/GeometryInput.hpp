#pragma once

#include "cfd/meshing/BackwardFacingStepGeometry.hpp"
#include "cfd/meshing/RectangleGeometry.hpp"
#include "cfd/meshing/RpaNozzleGeometry.hpp"

#include <variant>

namespace cfd
{

/// Closed set of geometric inputs supported by the built-in mesher.
using GeometryInput = std::variant<RectangleGeometry, BackwardFacingStepGeometry, RpaNozzleGeometry>;

} // namespace cfd
