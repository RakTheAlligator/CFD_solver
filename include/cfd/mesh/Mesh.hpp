#pragma once

#include "cfd/mesh/Boundary.hpp"
#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Node.hpp"
#include "cfd/mesh/Types.hpp"

#include <span>
#include <vector>

namespace cfd {

struct RawMeshData;

class Mesh {
public:
    [[nodiscard]]
    std::span<const Node> nodes() const noexcept
    {
        return nodes_;
    }

    [[nodiscard]]
    std::span<const CellType> cell_types() const noexcept
    {
        return cell_types_;
    }

    [[nodiscard]]
    std::span<const Index> cell_nodes() const noexcept
    {
        return cell_nodes_;
    }

    [[nodiscard]]
    std::span<const Index> cell_node_offsets() const noexcept
    {
        return cell_node_offsets_;
    }

    [[nodiscard]]
    std::span<const Face> faces() const noexcept
    {
        return faces_;
    }

    [[nodiscard]]
    std::span<const Index> cell_faces() const noexcept
    {
        return cell_faces_;
    }

    [[nodiscard]]
    std::span<const FaceAdjacency> face_adjacencies() const noexcept
    {
        return face_adjacencies_;
    }

    [[nodiscard]]
    std::span<const BoundaryId> face_boundary_ids() const noexcept
    {
        return face_boundary_ids_;
    }

    [[nodiscard]]
    std::span<const BoundaryGroup> boundary_groups() const noexcept
    {
        return boundary_groups_;
    }

    [[nodiscard]]
    Index node_count() const noexcept
    {
        return nodes_.size();
    }

    [[nodiscard]]
    Index cell_count() const noexcept
    {
        return cell_types_.size();
    }

    [[nodiscard]]
    Index face_count() const noexcept
    {
        return faces_.size();
    }

private:
    Mesh() = default;

    // Imported topology.
    std::vector<Node> nodes_;

    std::vector<CellType> cell_types_;
    std::vector<Index> cell_nodes_;
    std::vector<Index> cell_node_offsets_;

    std::vector<BoundaryGroup> boundary_groups_;

    // Constructed topology.
    std::vector<Face> faces_;
    std::vector<Index> cell_faces_;
    std::vector<FaceAdjacency> face_adjacencies_;
    std::vector<BoundaryId> face_boundary_ids_;

    friend Mesh build_mesh(RawMeshData&& raw_mesh);
};

} // namespace cfd