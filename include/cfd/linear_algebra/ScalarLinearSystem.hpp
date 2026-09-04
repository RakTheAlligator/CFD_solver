#pragma once

#include "cfd/mesh/Types.hpp"

#include <span>
#include <vector>

namespace cfd
{

class Mesh;

/// Fixed-cardinality finite-volume scalar linear system.
///
/// The matrix stores one diagonal value per cell and two directed
/// off-diagonal values per global face. Boundary-face off-diagonal entries are
/// unused. Matrix-vector products obtain owner/neighbor addressing directly
/// from the referenced Mesh.
///
/// @note The Mesh is not owned and must outlive this system.
/// @note Spans returned by this class do not own their data and must not outlive
///       the ScalarLinearSystem storage from which they were obtained.
class ScalarLinearSystem
{
  public:
    /// Constructs a zero-initialized system matching a fixed Mesh.
    explicit ScalarLinearSystem(const Mesh &mesh);

    ScalarLinearSystem(const ScalarLinearSystem &) = delete;
    ScalarLinearSystem &operator=(const ScalarLinearSystem &) = delete;

    ScalarLinearSystem(ScalarLinearSystem &&) noexcept = default;
    ScalarLinearSystem &operator=(ScalarLinearSystem &&) noexcept = delete;

    ~ScalarLinearSystem() = default;

    [[nodiscard]]
    Index cell_count() const noexcept;

    [[nodiscard]]
    Index face_count() const noexcept;

    [[nodiscard]]
    const Mesh &mesh() const noexcept;

    [[nodiscard]]
    std::span<double> diagonal() noexcept;

    [[nodiscard]]
    std::span<const double> diagonal() const noexcept;

    /// Returns A(owner, neighbor), indexed by global face ID.
    [[nodiscard]]
    std::span<double> owner_neighbor_coefficients() noexcept;

    /// Returns A(owner, neighbor), indexed by global face ID.
    [[nodiscard]]
    std::span<const double> owner_neighbor_coefficients() const noexcept;

    /// Returns A(neighbor, owner), indexed by global face ID.
    [[nodiscard]]
    std::span<double> neighbor_owner_coefficients() noexcept;

    /// Returns A(neighbor, owner), indexed by global face ID.
    [[nodiscard]]
    std::span<const double> neighbor_owner_coefficients() const noexcept;

    [[nodiscard]]
    std::span<double> rhs() noexcept;

    [[nodiscard]]
    std::span<const double> rhs() const noexcept;

    /// Clears all matrix coefficients without changing the right-hand side.
    void clear_matrix() noexcept;

    /// Clears the right-hand side without changing matrix coefficients.
    void clear_rhs() noexcept;

    /// Clears all matrix coefficients and the right-hand side.
    void clear() noexcept;

    /// Computes `output = A * input` without allocation.
    ///
    /// @throws std::invalid_argument If either span cardinality differs from
    ///         `cell_count()` or the spans overlap.
    void apply_matrix(std::span<const double> input, std::span<double> output) const;

  private:
    const Mesh *mesh_;
    std::vector<double> diagonal_;
    std::vector<double> owner_neighbor_coefficients_;
    std::vector<double> neighbor_owner_coefficients_;
    std::vector<double> rhs_;
};

} // namespace cfd
