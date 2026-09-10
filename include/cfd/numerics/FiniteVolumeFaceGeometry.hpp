#pragma once

#include "cfd/math/Vector2.hpp"
#include "cfd/mesh/Types.hpp"

namespace cfd
{

struct Point2;

/// Owner-oriented projection geometry between two points around a face.
///
/// `displacement` points from `reference_point` to `sample_point`, while
/// `area_dot_displacement` is `S_f . displacement` for the owner-oriented face
/// area vector. The value is strictly positive for a usable projection.
struct FaceProjectionGeometry
{
    Vector2 displacement;
    double area_dot_displacement{};
};

/// Geometry for projected interpolation across an internal finite-volume face.
///
/// `neighbor_weight` is the owner-to-neighbor line-intersection fraction
///
/// `S_f . (x_f - x_P) / (S_f . (x_N - x_P))`.
struct InternalFaceInterpolationGeometry
{
    Vector2 owner_to_neighbor;
    double area_dot_owner_to_neighbor{};
    double neighbor_weight{};
};

/// Computes and validates an owner-oriented face projection.
///
/// @param owner_oriented_area_vector Owner-oriented face area vector `S_f`.
/// @param face_length Magnitude `|S_f|` of `owner_oriented_area_vector` for the
///        same face; must be finite and strictly positive.
/// @throws std::runtime_error If the displacement or its positive projection
///         onto the owner-oriented face area vector is numerically unusable.
[[nodiscard]]
FaceProjectionGeometry compute_face_projection_geometry(Index face_id, const Point2 &reference_point,
                                                        const Point2 &sample_point,
                                                        const Vector2 &owner_oriented_area_vector, double face_length);

/// Computes and validates projected interpolation geometry for an internal face.
///
/// The intersection of the owner-neighbor line with the face plane must lie on
/// the owner-neighbor segment, within a roundoff-scaled tolerance.
///
/// @param owner_oriented_area_vector Owner-oriented face area vector `S_f`.
/// @param face_length Magnitude `|S_f|` of `owner_oriented_area_vector` for the
///        same face; must be finite and strictly positive.
/// @throws std::runtime_error If the projected geometry is numerically unusable
///         or the face-plane intersection lies outside the cell-center segment.
[[nodiscard]]
InternalFaceInterpolationGeometry compute_internal_face_interpolation_geometry(
    Index face_id, const Point2 &owner_center, const Point2 &neighbor_center, const Point2 &face_center,
    const Vector2 &owner_oriented_area_vector, double face_length);

} // namespace cfd
