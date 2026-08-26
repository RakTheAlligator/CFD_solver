#pragma once

#include "cfd/mesh/Boundary.hpp"
#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Node.hpp"
#include "cfd/mesh/Types.hpp"
#include "cfd/mesh/Vector2.hpp"

#include <span>
#include <vector>

namespace cfd
{

struct RawMeshData;
struct MeshBuildResult;

class Mesh
{
  public:
    Mesh(const Mesh &) = delete;
    Mesh &operator=(const Mesh &) = delete;

    Mesh(Mesh &&) noexcept = default;
    Mesh &operator=(Mesh &&) noexcept = default;

    ~Mesh() = default;

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

    [[nodiscard]]
    std::span<const double> cell_areas() const noexcept
    {
        return cell_areas_;
    }

    [[nodiscard]]
    std::span<const Vector2> cell_centers() const noexcept
    {
        return cell_centers_;
    }

    [[nodiscard]]
    std::span<const Vector2> face_centers() const noexcept
    {
        return face_centers_;
    }

    [[nodiscard]]
    std::span<const double> face_lengths() const noexcept
    {
        return face_lengths_;
    }

    [[nodiscard]]
    std::span<const Vector2> face_area_vectors() const noexcept
    {
        return face_area_vectors_;
    }
    [[nodiscard]]
    std::span<const double> cell_qualities() const noexcept
    {
        return cell_qualities_;
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

    // Constructed geometry
    std::vector<double> cell_areas_;
    std::vector<Vector2> cell_centers_;
    std::vector<double> cell_qualities_;

    std::vector<Vector2> face_centers_;
    std::vector<double> face_lengths_;
    std::vector<Vector2> face_area_vectors_;

    friend MeshBuildResult build_mesh(RawMeshData &&raw_mesh);
};

} // namespace cfd