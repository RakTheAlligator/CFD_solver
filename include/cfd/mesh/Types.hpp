#pragma once

#include <cstddef>
#include <limits>

namespace cfd
{

/// Common unsigned integer type for mesh indices and size-related quantities.
using Index = std::size_t;

/// Sentinel representing the absence of a valid Index value.
inline constexpr Index invalid_index{std::numeric_limits<Index>::max()};

} // namespace cfd