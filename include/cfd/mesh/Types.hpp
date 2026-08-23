#pragma once

#include <cstddef>
#include <limits>

namespace cfd {

using Index = std::size_t;

inline constexpr Index invalid_index { std::numeric_limits<Index>::max() };

} // namespace cfd