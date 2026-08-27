#pragma once

#include <cstdint>

namespace cfd
{

/// Supported two-dimensional cell topologies.
///
/// The fixed-width underlying type keeps per-cell topology storage compact.
enum class CellType : std::uint8_t
{
    Triangle,
    Quadrilateral
};

} // namespace cfd