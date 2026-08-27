#pragma once

#include "cfd/mesh/Types.hpp"

#include <array>

namespace cfd
{

/// Topological face of a two-dimensional mesh.
///
/// A face is an edge connecting two internal mesh node IDs. The node ordering
/// must not be used to infer the orientation of face-based geometric vectors.
struct Face
{
    std::array<Index, 2> node_ids{};
};

/// Cell adjacency associated with a mesh face.
///
/// In a validated mesh, `owner` identifies one adjacent cell. Internal faces
/// also have a valid `neighbor`, whereas boundary faces use `invalid_index`.
///
/// The owner/neighbor distinction is purely topological and does not imply a
/// flow direction or an upwind/downwind relationship.
struct FaceAdjacency
{
    Index owner{invalid_index};
    Index neighbor{invalid_index};

    [[nodiscard]]
    bool is_boundary() const noexcept
    {
        return neighbor == invalid_index;
    }
};

} // namespace cfd