#pragma once

#include "cfd/math/Vector2.hpp"

#include <span>
#include <vector>

namespace cfd
{

class Mesh;

/// Precomputed finite-volume metrics for a full revolution about the x axis.
///
/// Mesh coordinates are interpreted as (x,r), with r = y >= 0. For the
/// discretized polygons and straight faces, the exact integrals are
/// V_C = 2*pi*A_C*r_centroid, S_f = 2*pi*r_mid*S_f,2D and G_C = 2*pi*A_C.
/// These metrics do not describe the exact revolution of curved CAD walls.
/// Construction checks sum(sign*S_f,x) = 0 and sum(sign*S_f,r) = G_C for
/// every cell, using owner/neighbor signs and scale-relative roundoff bounds.
/// The bound is 128*epsilon times the sum of absolute component contributions
/// (including G_C radially), with no unit-sized tolerance floor.
///
/// @note The Mesh is a non-owning dependency: it must outlive this object and must not be
///       moved from or replaced. No nodes or connectivity are copied.
///       Metric arrays are owned, contiguous and use the original Mesh IDs.
///       Read-only access allocates nothing; returned spans must not outlive
///       their storage, including when this object is moved from.
class AxisymmetricGeometry
{
  public:
    /// Builds and validates all metrics once, with linear work in mesh size.
    /// @throws std::runtime_error If coordinates are non-finite, a radius is
    ///         negative, a volume or source coefficient is not finite positive,
    ///         a face metric is invalid, or a cell integral identity fails.
    ///         Zero face area is allowed only when both endpoints lie on r=0.
    explicit AxisymmetricGeometry(const Mesh &mesh);

    AxisymmetricGeometry(const AxisymmetricGeometry &) = delete;
    AxisymmetricGeometry &operator=(const AxisymmetricGeometry &) = delete;
    AxisymmetricGeometry(AxisymmetricGeometry &&) noexcept = default;
    AxisymmetricGeometry &operator=(AxisymmetricGeometry &&) noexcept = delete;
    ~AxisymmetricGeometry() = default;

    [[nodiscard]]
    const Mesh &mesh() const noexcept
    {
        return *mesh_;
    }

    /// Full-revolution cell volumes in cubic metres, finite and positive.
    [[nodiscard]]
    std::span<const double> cell_volumes() const noexcept
    {
        return cell_volumes_;
    }

    /// Integrated face-area vectors in square metres, outward from the owner.
    /// The neighbor uses the opposite vector. An axis face has a zero vector;
    /// callers must not divide by its magnitude to construct a normal.
    [[nodiscard]]
    std::span<const Vector2> face_area_vectors() const noexcept
    {
        return face_area_vectors_;
    }

    /// G_C = 2*pi*A_C in square metres, for a future radial pressure source.
    /// This class supplies geometry only; no physical source is assembled.
    [[nodiscard]]
    std::span<const double> cell_radial_source_coefficients() const noexcept
    {
        return cell_radial_source_coefficients_;
    }

  private:
    const Mesh *mesh_;
    std::vector<double> cell_volumes_;
    std::vector<Vector2> face_area_vectors_;
    std::vector<double> cell_radial_source_coefficients_;
};

} // namespace cfd
