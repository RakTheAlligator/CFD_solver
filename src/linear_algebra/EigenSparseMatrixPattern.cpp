#include "cfd/linear_algebra/EigenSparseMatrixPattern.hpp"

#include "cfd/linear_algebra/ScalarLinearSystem.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>

namespace cfd
{
namespace
{

void require_finite(const std::span<const double> values, const char *message)
{
    for (const double value : values)
    {
        if (!std::isfinite(value))
        {
            throw std::invalid_argument(message);
        }
    }
}

} // namespace

EigenSparseMatrixPattern::EigenSparseMatrixPattern(const Mesh &mesh)
    : mesh_(&mesh), face_adjacencies_(mesh.face_adjacencies().data()), cell_count_(mesh.cell_count()),
      face_count_(mesh.face_count())
{
    if (!std::in_range<Eigen::Index>(cell_count_))
    {
        throw std::invalid_argument("Scalar linear system cardinality exceeds the Eigen index range.");
    }
    constexpr Index maximum_entries{static_cast<Index>(std::numeric_limits<StorageIndex>::max())};
    static_assert(std::in_range<Index>(std::numeric_limits<StorageIndex>::max()));
    if (cell_count_ > maximum_entries)
    {
        throw std::invalid_argument("Scalar linear system cardinality exceeds the Eigen sparse storage index range.");
    }
    const auto adjacencies{mesh.face_adjacencies()};
    const Index internal_count{
        static_cast<Index>(std::count_if(adjacencies.begin(), adjacencies.end(),
                                         [](const FaceAdjacency &adjacency) { return !adjacency.is_boundary(); }))};
    if (internal_count > (maximum_entries - cell_count_) / 2)
    {
        throw std::invalid_argument("Scalar linear system entry count exceeds the Eigen sparse storage index range.");
    }

    using Triplet = Eigen::Triplet<double, StorageIndex>;
    std::vector<Triplet> entries;
    entries.reserve(cell_count_ + 2 * internal_count);
    for (Index cell = 0; cell < cell_count_; ++cell)
    {
        const auto index{static_cast<StorageIndex>(cell)};
        entries.emplace_back(index, index, 0.0);
    }
    for (const FaceAdjacency &adjacency : adjacencies)
    {
        if (!adjacency.is_boundary())
        {
            entries.emplace_back(static_cast<StorageIndex>(adjacency.owner),
                                 static_cast<StorageIndex>(adjacency.neighbor), 0.0);
            entries.emplace_back(static_cast<StorageIndex>(adjacency.neighbor),
                                 static_cast<StorageIndex>(adjacency.owner), 0.0);
        }
    }
    SparseMatrix pattern{static_cast<Eigen::Index>(cell_count_), static_cast<Eigen::Index>(cell_count_)};
    pattern.setFromTriplets(entries.begin(), entries.end());
    const std::span<const StorageIndex> offsets{pattern.outerIndexPtr(), cell_count_ + 1};
    const std::span<const StorageIndex> rows{pattern.innerIndexPtr(), static_cast<Index>(pattern.nonZeros())};
    column_offsets_.assign(offsets.begin(), offsets.end());
    row_indices_.assign(rows.begin(), rows.end());

    const auto find_slot = [this](const Index row, const Index column) {
        const auto begin{row_indices_.begin() + column_offsets_[column]};
        const auto end{row_indices_.begin() + column_offsets_[column + 1]};
        const auto position{std::lower_bound(begin, end, static_cast<StorageIndex>(row))};
        return static_cast<StorageIndex>(position - row_indices_.begin());
    };
    diagonal_slots_.resize(cell_count_);
    for (Index cell = 0; cell < cell_count_; ++cell)
    {
        diagonal_slots_[cell] = find_slot(cell, cell);
    }
    internal_face_slots_.reserve(internal_count);
    for (Index face = 0; face < face_count_; ++face)
    {
        const FaceAdjacency &adjacency{adjacencies[face]};
        if (!adjacency.is_boundary())
        {
            internal_face_slots_.push_back(
                {face, find_slot(adjacency.owner, adjacency.neighbor), find_slot(adjacency.neighbor, adjacency.owner)});
        }
    }
}

bool EigenSparseMatrixPattern::matches(const Mesh &mesh) const noexcept
{
    // Mesh exposes read-only topology. Also detect replacement/move of its storage
    // without hashing the connectivity on every preparation.
    return &mesh == mesh_ && mesh.cell_count() == cell_count_ && mesh.face_count() == face_count_ &&
           mesh.face_adjacencies().data() == face_adjacencies_;
}

void EigenSparseMatrixPattern::validate_system(const ScalarLinearSystem &system) const
{
    if (!matches(system.mesh()) || system.cell_count() != cell_count_ || system.face_count() != face_count_)
    {
        throw std::invalid_argument("Eigen sparse pattern requires the same unchanged Mesh and system cardinalities.");
    }
    require_finite(system.diagonal(), "Scalar linear system diagonal must contain only finite values.");
    require_finite(system.owner_neighbor_coefficients(),
                   "Scalar linear system owner-neighbor coefficients must contain only finite values.");
    require_finite(system.neighbor_owner_coefficients(),
                   "Scalar linear system neighbor-owner coefficients must contain only finite values.");
}

void EigenSparseMatrixPattern::initialize_matrix(SparseMatrix &matrix) const
{
    const auto size{static_cast<Eigen::Index>(cell_count_)};
    SparseMatrix initialized{size, size};
    initialized.reserve(static_cast<Eigen::Index>(row_indices_.size()));
    for (Index column = 0; column < cell_count_; ++column)
    {
        initialized.startVec(static_cast<Eigen::Index>(column));
        for (StorageIndex slot = column_offsets_[column]; slot < column_offsets_[column + 1]; ++slot)
        {
            initialized.insertBackByOuterInner(static_cast<Eigen::Index>(column), row_indices_[slot]) = 0.0;
        }
    }
    initialized.finalize();
    matrix.swap(initialized);
}

void EigenSparseMatrixPattern::update_values(const ScalarLinearSystem &system, SparseMatrix &matrix) const noexcept
{
    const std::span<double> values{matrix.valuePtr(), row_indices_.size()};
    for (Index cell = 0; cell < cell_count_; ++cell)
    {
        values[diagonal_slots_[cell]] = system.diagonal()[cell];
    }
    for (const InternalFaceSlots &slots : internal_face_slots_)
    {
        values[slots.owner_neighbor] = system.owner_neighbor_coefficients()[slots.face_id];
        values[slots.neighbor_owner] = system.neighbor_owner_coefficients()[slots.face_id];
    }
}

EigenSparseMatrixPattern::SparseMatrix EigenSparseMatrixPattern::create_matrix(const ScalarLinearSystem &system) const
{
    validate_system(system);
    SparseMatrix matrix;
    initialize_matrix(matrix);
    update_values(system, matrix);
    return matrix;
}

} // namespace cfd
