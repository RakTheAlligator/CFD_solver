#include "cfd/numerics/FiniteVolumeFaceGeometry.hpp"

#include "cfd/math/Point2.hpp"
#include "cfd/math/Vector2.hpp"
#include "cfd/mesh/Types.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace cfd
{
namespace
{

constexpr double relative_projection_tolerance{64.0 * std::numeric_limits<double>::epsilon()};

struct ValidatedFaceProjection
{
    FaceProjectionGeometry geometry;
    double displacement_norm{};
};

[[noreturn]]
void throw_unusable_face_geometry(const Index face_id, const std::string &reason)
{
    throw std::runtime_error("Finite-volume geometry rejected face " + std::to_string(face_id) + ": " + reason);
}

[[nodiscard]]
double dot(const Vector2 &first, const Vector2 &second) noexcept
{
    return first.x * second.x + first.y * second.y;
}

[[nodiscard]]
ValidatedFaceProjection compute_validated_face_projection(const Index face_id, const Point2 &reference_point,
                                                          const Point2 &sample_point,
                                                          const Vector2 &owner_oriented_area_vector,
                                                          const double face_length)
{
    const Vector2 displacement{
        sample_point.x - reference_point.x,
        sample_point.y - reference_point.y,
    };
    const double displacement_norm{std::hypot(displacement.x, displacement.y)};
    const double projection_scale{face_length * displacement_norm};
    const double area_dot_displacement{dot(owner_oriented_area_vector, displacement)};
    const double minimum_projection{relative_projection_tolerance * projection_scale};

    if (!std::isfinite(projection_scale) || !(projection_scale > 0.0) || !std::isfinite(area_dot_displacement) ||
        !(area_dot_displacement > minimum_projection))
    {
        throw_unusable_face_geometry(face_id, "Sf dot d is not a usable positive projection.");
    }

    return {{displacement, area_dot_displacement}, displacement_norm};
}

} // namespace

FaceProjectionGeometry compute_face_projection_geometry(const Index face_id, const Point2 &reference_point,
                                                        const Point2 &sample_point,
                                                        const Vector2 &owner_oriented_area_vector,
                                                        const double face_length)
{
    return compute_validated_face_projection(face_id, reference_point, sample_point, owner_oriented_area_vector,
                                             face_length)
        .geometry;
}

InternalFaceInterpolationGeometry compute_internal_face_interpolation_geometry(
    const Index face_id, const Point2 &owner_center, const Point2 &neighbor_center, const Point2 &face_center,
    const Vector2 &owner_oriented_area_vector, const double face_length)
{
    const ValidatedFaceProjection center_projection{compute_validated_face_projection(
        face_id, owner_center, neighbor_center, owner_oriented_area_vector, face_length)};
    const Vector2 owner_to_face{
        face_center.x - owner_center.x,
        face_center.y - owner_center.y,
    };
    const double owner_to_face_projection{dot(owner_oriented_area_vector, owner_to_face)};
    const double owner_to_face_distance{std::hypot(owner_to_face.x, owner_to_face.y)};
    const double intersection_tolerance{relative_projection_tolerance * face_length *
                                        (center_projection.displacement_norm + owner_to_face_distance)};

    if (!std::isfinite(owner_to_face_projection) || !std::isfinite(intersection_tolerance) ||
        owner_to_face_projection < -intersection_tolerance ||
        owner_to_face_projection > center_projection.geometry.area_dot_displacement + intersection_tolerance)
    {
        throw_unusable_face_geometry(face_id, "the face-line intersection lies outside the owner-neighbor segment.");
    }

    const double neighbor_weight{owner_to_face_projection / center_projection.geometry.area_dot_displacement};
    if (!std::isfinite(neighbor_weight))
    {
        throw_unusable_face_geometry(face_id, "the face interpolation weight is non-finite.");
    }

    return {
        center_projection.geometry.displacement,
        center_projection.geometry.area_dot_displacement,
        neighbor_weight,
    };
}

} // namespace cfd
