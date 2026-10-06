#pragma once

#include "cfd/mesh/Types.hpp"

#include <Eigen/SparseCore>

#include <vector>

namespace cfd
{

class Mesh;
struct FaceAdjacency;
class ScalarLinearSystem;

/// Immutable compressed-column pattern and finite-volume coefficient mapping.
///
/// Includes every diagonal and both directions of every internal face, even
/// when their coefficients are zero. Validated Mesh topology guarantees unique
/// cell pairs. Numerical matrices and preconditioners are not shared.
///
/// @note The Mesh is not owned; it must outlive this object and must not be
///       moved from or replaced while the pattern is in use.
class EigenSparseMatrixPattern
{
  public:
    using SparseMatrix = Eigen::SparseMatrix<double>;

    /// Builds the pattern once from validated Mesh topology.
    /// @throws std::invalid_argument If dimensions or entries exceed Eigen's index ranges.
    explicit EigenSparseMatrixPattern(const Mesh &mesh);

    EigenSparseMatrixPattern(const EigenSparseMatrixPattern &) = delete;
    EigenSparseMatrixPattern &operator=(const EigenSparseMatrixPattern &) = delete;
    EigenSparseMatrixPattern(EigenSparseMatrixPattern &&) = delete;
    EigenSparseMatrixPattern &operator=(EigenSparseMatrixPattern &&) = delete;
    ~EigenSparseMatrixPattern() = default;

    /// Creates an independent numerical matrix, retaining explicit zero entries.
    /// @throws std::invalid_argument If the system references another Mesh,
    ///         its topology storage changed, cardinalities differ, or coefficients are non-finite.
    [[nodiscard]]
    SparseMatrix create_matrix(const ScalarLinearSystem &system) const;

  private:
    using StorageIndex = SparseMatrix::StorageIndex;

    struct InternalFaceSlots
    {
        Index face_id{};
        StorageIndex owner_neighbor{};
        StorageIndex neighbor_owner{};
    };

    [[nodiscard]]
    bool matches(const Mesh &mesh) const noexcept;
    void validate_system(const ScalarLinearSystem &system) const;
    void initialize_matrix(SparseMatrix &matrix) const;
    void update_values(const ScalarLinearSystem &system, SparseMatrix &matrix) const noexcept;

    const Mesh *mesh_;
    const FaceAdjacency *face_adjacencies_;
    Index cell_count_;
    Index face_count_;
    std::vector<StorageIndex> column_offsets_;
    std::vector<StorageIndex> row_indices_;
    std::vector<StorageIndex> diagonal_slots_;
    std::vector<InternalFaceSlots> internal_face_slots_;

    friend class EigenBiCGSTABSolver;
    friend class EigenConjugateGradientSolver;
};

} // namespace cfd
