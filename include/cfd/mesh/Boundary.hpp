#pragma once

#include "cfd/mesh/Types.hpp"

#include <string>

namespace cfd
{

/// Internal identifier of a logical mesh boundary group.
///
/// Boundary IDs belong to the solver mesh representation and are independent
/// of tags used by external mesh generators.
using BoundaryId = Index;

/// Sentinel indicating that no boundary group is associated with a mesh entity.
inline constexpr BoundaryId invalid_boundary_id{invalid_index};

/// Logical boundary group attached to mesh boundary entities.
///
/// A boundary group provides only geometric/topological identity. Numerical
/// boundary conditions are defined separately from the mesh representation.
struct BoundaryGroup
{
    BoundaryId id{};
    std::string name;
};

} // namespace cfd