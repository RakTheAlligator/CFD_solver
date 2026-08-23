#pragma once

#include "cfd/mesh/Types.hpp"

#include <string>

namespace cfd
{

using BoundaryId = Index;

inline constexpr BoundaryId invalid_boundary_id = invalid_index;

struct BoundaryGroup
{
    BoundaryId id{};
    std::string name;
};

} // namespace cfd