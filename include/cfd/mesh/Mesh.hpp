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

/// Validated two-dimensional mesh representation.
///
/// Mesh owns the topology and geometry required by the numerical solver and
/// exposes them through non-owning, read-only spans. Internal node, cell, and
/// face IDs are zero-based and index directly into their corresponding arrays.
///
/// Cell-to-node connectivity uses a flattened, CSR-like layout defined by
/// `cell_nodes()` and `cell_node_offsets()`. Cell-to-face connectivity reuses
/// the same offsets because every two-dimensional polygonal cell has one face
/// per node.
///
/// @invariant `cell_node_offsets().size() == cell_count() + 1`.
/// @invariant `cell_faces().size() == cell_nodes().size()`.
/// @invariant Per-cell geometry arrays contain `cell_count()` entries.
/// @invariant Per-face topology and geometry arrays contain `face_count()` entries.
///
/// @note Spans returned by this class do not own their data and must not outlive
///       the Mesh storage from which they were obtained.
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

    /// Returns the flattened cell-to-node connectivity.
    ///
    /// Node IDs for cell `c` occupy
    /// `[cell_node_offsets()[c], cell_node_offsets()[c + 1])`.
    [[nodiscard]]
    std::span<const Index> cell_nodes() const noexcept
    {
        return cell_nodes_;
    }

    /// Returns the offsets delimiting each cell in `cell_nodes()` and
    /// `cell_faces()`.
    ///
    /// The first offset is zero and the final offset equals the size of both
    /// flattened connectivity arrays.
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

    /// Returns the flattened cell-to-face connectivity.
    ///
    /// This array uses `cell_node_offsets()` as its offsets. For each cell,
    /// local face ordering follows local node ordering.
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

    /// Returns the boundary-group ID associated with each face.
    ///
    /// Internal faces contain `invalid_boundary_id`; validated boundary faces
    /// contain a valid ID referring to `boundary_groups()`.
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

    /// Returns cell areas in square metres.
    [[nodiscard]]
    std::span<const double> cell_areas() const noexcept
    {
        return cell_areas_;
    }

    /// Returns cell-centroid coordinates in metres.
    [[nodiscard]]
    std::span<const Vector2> cell_centers() const noexcept
    {
        return cell_centers_;
    }

    /// Returns face-center coordinates in metres.
    [[nodiscard]]
    std::span<const Vector2> face_centers() const noexcept
    {
        return face_centers_;
    }

    /// Returns face lengths in metres.
    [[nodiscard]]
    std::span<const double> face_lengths() const noexcept
    {
        return face_lengths_;
    }

    /// Returns oriented face-area vectors.
    ///
    /// In two dimensions, each vector has magnitude equal to the corresponding
    /// face length and is oriented outward from the owner cell. For an internal
    /// face, it therefore points from owner to neighbor.
    [[nodiscard]]
    std::span<const Vector2> face_area_vectors() const noexcept
    {
        return face_area_vectors_;
    }

    /// Returns the dimensionless quality metric of each cell.
    ///
    /// Valid values lie in `(0, 1]`, with larger values representing better
    /// cell quality and `1` corresponding to the ideal shape for the metric.
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

    // Constructed geometry.
    std::vector<double> cell_areas_;
    std::vector<Vector2> cell_centers_;
    std::vector<double> cell_qualities_;

    std::vector<Vector2> face_centers_;
    std::vector<double> face_lengths_;
    std::vector<Vector2> face_area_vectors_;

    friend MeshBuildResult build_mesh(RawMeshData &&raw_mesh);
};

} // namespace cfd