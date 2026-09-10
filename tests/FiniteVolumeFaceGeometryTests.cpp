#include "cfd/numerics/FiniteVolumeFaceGeometry.hpp"

#include "cfd/math/Point2.hpp"
#include "cfd/math/Vector2.hpp"

#include "support/TestUtils.hpp"

#include <cmath>
#include <stdexcept>
#include <type_traits>

namespace
{

using cfd::test::require_near;
using cfd::test::require_throws;

static_assert(std::is_trivially_copyable_v<cfd::FaceProjectionGeometry>);
static_assert(std::is_trivially_copyable_v<cfd::InternalFaceInterpolationGeometry>);

void test_boundary_projection_geometry()
{
    const cfd::Point2 owner_center{0.0, 0.0};
    const cfd::Point2 face_center{0.75, 0.5};
    const cfd::Vector2 owner_oriented_area_vector{2.0, 1.0};
    const cfd::FaceProjectionGeometry geometry{cfd::compute_face_projection_geometry(
        4, owner_center, face_center, owner_oriented_area_vector, std::sqrt(5.0))};

    require_near(geometry.displacement.x, 0.75, 0.0, "Boundary projection x displacement is incorrect.");
    require_near(geometry.displacement.y, 0.5, 0.0, "Boundary projection y displacement is incorrect.");
    require_near(geometry.area_dot_displacement, 2.0, 0.0,
                 "Boundary owner-oriented projected distance numerator is incorrect.");
}

void test_orthogonal_internal_face_geometry()
{
    const cfd::InternalFaceInterpolationGeometry geometry{
        cfd::compute_internal_face_interpolation_geometry(7, {0.0, 0.0}, {2.0, 0.0}, {0.5, 0.25}, {1.0, 0.0}, 1.0)};

    require_near(geometry.owner_to_neighbor.x, 2.0, 0.0, "Orthogonal owner-neighbor x displacement is incorrect.");
    require_near(geometry.owner_to_neighbor.y, 0.0, 0.0, "Orthogonal owner-neighbor y displacement is incorrect.");
    require_near(geometry.area_dot_owner_to_neighbor, 2.0, 0.0, "Orthogonal owner-neighbor projection is incorrect.");
    require_near(geometry.neighbor_weight, 0.25, 0.0, "Orthogonal neighbor interpolation weight is incorrect.");
}

void test_skewed_internal_face_geometry()
{
    const cfd::InternalFaceInterpolationGeometry geometry{cfd::compute_internal_face_interpolation_geometry(
        9, {0.0, 0.0}, {2.0, 1.0}, {0.8, -0.2}, {1.0, 1.0}, std::sqrt(2.0))};

    require_near(geometry.owner_to_neighbor.x, 2.0, 0.0, "Skewed owner-neighbor x displacement is incorrect.");
    require_near(geometry.owner_to_neighbor.y, 1.0, 0.0, "Skewed owner-neighbor y displacement is incorrect.");
    require_near(geometry.area_dot_owner_to_neighbor, 3.0, 0.0, "Skewed owner-neighbor projection is incorrect.");
    require_near(geometry.neighbor_weight, 0.2, 1.0e-15, "Skewed neighbor interpolation weight is incorrect.");
}

void test_rejects_unusable_projection_and_intersection()
{
    require_throws<std::runtime_error>(
        []() { static_cast<void>(cfd::compute_face_projection_geometry(0, {0.0, 0.0}, {0.0, 0.0}, {1.0, 0.0}, 1.0)); },
        "Finite-volume geometry accepted zero displacement.");

    require_throws<std::runtime_error>(
        []() { static_cast<void>(cfd::compute_face_projection_geometry(1, {0.0, 0.0}, {0.0, 1.0}, {1.0, 0.0}, 1.0)); },
        "Finite-volume geometry accepted a tangential displacement.");

    require_throws<std::runtime_error>(
        []() { static_cast<void>(cfd::compute_face_projection_geometry(2, {0.0, 0.0}, {1.0, 0.0}, {-1.0, 0.0}, 1.0)); },
        "Finite-volume geometry accepted an area vector opposed to the displacement.");

    require_throws<std::runtime_error>(
        []() {
            static_cast<void>(cfd::compute_internal_face_interpolation_geometry(3, {0.0, 0.0}, {1.0, 0.0}, {2.0, 0.0},
                                                                                {1.0, 0.0}, 1.0));
        },
        "Finite-volume geometry accepted an intersection outside the owner-neighbor segment.");
}

} // namespace

int main()
{
    int failure_count{};

    failure_count += cfd::test::run_test("finite-volume boundary projection", test_boundary_projection_geometry);
    failure_count +=
        cfd::test::run_test("finite-volume orthogonal internal geometry", test_orthogonal_internal_face_geometry);
    failure_count += cfd::test::run_test("finite-volume skewed internal geometry", test_skewed_internal_face_geometry);
    failure_count += cfd::test::run_test("finite-volume face geometry validation",
                                         test_rejects_unusable_projection_and_intersection);

    return cfd::test::finish_tests(failure_count, "finite-volume face geometry");
}
