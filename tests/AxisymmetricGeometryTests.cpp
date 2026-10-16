#include "cfd/geometry/AxisymmetricGeometry.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/meshing/RawMeshData.hpp"
#include "support/MeshFixtures.hpp"
#include "support/TestUtils.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace
{

using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws_with_message;
constexpr double pi{std::numbers::pi};
constexpr double tolerance{1.0e-12};

static_assert(!std::is_copy_constructible_v<cfd::AxisymmetricGeometry>);
static_assert(!std::is_copy_assignable_v<cfd::AxisymmetricGeometry>);
static_assert(std::is_nothrow_move_constructible_v<cfd::AxisymmetricGeometry>);
static_assert(!std::is_move_assignable_v<cfd::AxisymmetricGeometry>);
static_assert(std::is_same_v<decltype(std::declval<const cfd::AxisymmetricGeometry &>().cell_volumes()),
                             std::span<const double>>);
static_assert(std::is_same_v<decltype(std::declval<const cfd::AxisymmetricGeometry &>().face_area_vectors()),
                             std::span<const cfd::Vector2>>);
static_assert(
    std::is_same_v<decltype(std::declval<const cfd::AxisymmetricGeometry &>().cell_radial_source_coefficients()),
                   std::span<const double>>);

cfd::RawMeshData rectangle(const double length, const double inner_radius, const double outer_radius)
{
    auto raw{cfd::test::make_single_quadrilateral_raw_mesh()};
    raw.nodes = {{0.0, inner_radius}, {length, inner_radius}, {length, outer_radius}, {0.0, outer_radius}};
    return raw;
}

void require_face_vector(const cfd::Mesh &mesh, const cfd::AxisymmetricGeometry &geometry, const cfd::Index first,
                         const cfd::Index second, const cfd::Vector2 expected)
{
    for (cfd::Index face = 0; face < mesh.face_count(); ++face)
    {
        const auto &nodes{mesh.faces()[face].node_ids};
        if ((nodes[0] == first && nodes[1] == second) || (nodes[0] == second && nodes[1] == first))
        {
            require_near(geometry.face_area_vectors()[face].x, expected.x, tolerance, "Incorrect axial face area.");
            require_near(geometry.face_area_vectors()[face].y, expected.y, tolerance, "Incorrect radial face area.");
            return;
        }
    }
    cfd::test::fail("Expected face was not found.");
}

void require_cell_identities(const cfd::Mesh &mesh, const cfd::AxisymmetricGeometry &geometry,
                             const std::span<const double> expected_sources)
{
    require(expected_sources.size() == mesh.cell_count(), "Invalid test source cardinality.");
    std::vector<cfd::Vector2> sums(mesh.cell_count());
    std::vector<double> scales(expected_sources.begin(), expected_sources.end());
    // Accumulate by face adjacency, independently of the production cell-face traversal.
    for (cfd::Index face = 0; face < mesh.face_count(); ++face)
    {
        const auto &adjacency{mesh.face_adjacencies()[face]};
        const auto &area{geometry.face_area_vectors()[face]};
        sums[adjacency.owner].x += area.x;
        sums[adjacency.owner].y += area.y;
        scales[adjacency.owner] += std::abs(area.x) + std::abs(area.y);
        if (!adjacency.is_boundary())
        {
            sums[adjacency.neighbor].x -= area.x;
            sums[adjacency.neighbor].y -= area.y;
            scales[adjacency.neighbor] += std::abs(area.x) + std::abs(area.y);
        }
    }
    for (cfd::Index cell = 0; cell < mesh.cell_count(); ++cell)
    {
        require_near(sums[cell].x / scales[cell], 0.0, tolerance, "Axial cell closure failed.");
        require_near(sums[cell].y / scales[cell], expected_sources[cell] / scales[cell], tolerance,
                     "Radial cell/source identity failed.");
        require_near(geometry.cell_radial_source_coefficients()[cell] / expected_sources[cell], 1.0, tolerance,
                     "Incorrect integrated radial source coefficient.");
    }
}

void test_cylinder_volume_and_face_integrals()
{
    const auto built{cfd::build_mesh(rectangle(2.0, 0.0, 3.0))};
    const auto &mesh{built.mesh};
    const cfd::AxisymmetricGeometry geometry{mesh};
    require(&geometry.mesh() == &mesh, "Geometry must reference the original Mesh.");
    require(geometry.cell_volumes().size() == mesh.cell_count() &&
                geometry.face_area_vectors().size() == mesh.face_count() &&
                geometry.cell_radial_source_coefficients().size() == mesh.cell_count(),
            "Metric cardinalities do not match Mesh.");
    // A solid cylinder: pi * radius^2 * length.
    require_near(geometry.cell_volumes()[0], 18.0 * pi, tolerance, "Incorrect cylinder volume.");
    require_face_vector(mesh, geometry, 0, 1, {0.0, 0.0});
    require_face_vector(mesh, geometry, 1, 2, {9.0 * pi, 0.0});
    require_face_vector(mesh, geometry, 2, 3, {0.0, 12.0 * pi});
    require_face_vector(mesh, geometry, 3, 0, {-9.0 * pi, 0.0});
    require_cell_identities(mesh, geometry, std::array{12.0 * pi});
}

void test_annular_cylinder()
{
    const auto built{cfd::build_mesh(rectangle(2.0, 1.0, 3.0))};
    const cfd::AxisymmetricGeometry geometry{built.mesh};
    // Difference of cylinders with radii 3 and 1 and length 2.
    require_near(geometry.cell_volumes()[0], 16.0 * pi, tolerance, "Incorrect annular volume.");
    require_face_vector(built.mesh, geometry, 0, 1, {0.0, -4.0 * pi});
    require_face_vector(built.mesh, geometry, 1, 2, {8.0 * pi, 0.0});
    require_face_vector(built.mesh, geometry, 2, 3, {0.0, 12.0 * pi});
    require_cell_identities(built.mesh, geometry, std::array{8.0 * pi});
}

void test_triangle_is_a_cone()
{
    auto raw{cfd::test::make_single_triangle_raw_mesh()};
    raw.nodes = {{0.0, 0.0}, {2.0, 0.0}, {0.0, 3.0}};
    const auto built{cfd::build_mesh(std::move(raw))};
    const cfd::AxisymmetricGeometry geometry{built.mesh};
    // A cone: pi * radius^2 * height / 3, not a midpoint approximation in x.
    require_near(geometry.cell_volumes()[0], 6.0 * pi, tolerance, "Incorrect cone volume.");
    require_face_vector(built.mesh, geometry, 1, 2, {9.0 * pi, 6.0 * pi});
    require_face_vector(built.mesh, geometry, 2, 0, {-9.0 * pi, 0.0});
    require_cell_identities(built.mesh, geometry, std::array{6.0 * pi});
}

void test_inclined_quadrilateral()
{
    auto raw{cfd::test::make_single_quadrilateral_raw_mesh()};
    raw.nodes = {{0.0, 1.0}, {2.0, 2.0}, {2.0, 4.0}, {0.0, 3.0}};
    const auto built{cfd::build_mesh(std::move(raw))};
    const cfd::AxisymmetricGeometry geometry{built.mesh};
    // Washers: pi * integral_0^2 [(3+x/2)^2 - (1+x/2)^2] dx = 20*pi.
    require_near(geometry.cell_volumes()[0], 20.0 * pi, tolerance, "Incorrect inclined annular volume.");
    require_face_vector(built.mesh, geometry, 0, 1, {3.0 * pi, -6.0 * pi});
    require_face_vector(built.mesh, geometry, 1, 2, {12.0 * pi, 0.0});
    require_face_vector(built.mesh, geometry, 2, 3, {-7.0 * pi, 14.0 * pi});
    require_face_vector(built.mesh, geometry, 3, 0, {-8.0 * pi, 0.0});
    require_cell_identities(built.mesh, geometry, std::array{8.0 * pi});
}

void test_axis_face_is_exactly_zero()
{
    auto raw{rectangle(1.0, -0.0, 1.0)};
    const auto built{cfd::build_mesh(std::move(raw))};
    const cfd::AxisymmetricGeometry geometry{built.mesh};
    cfd::Index axis_faces{};
    for (cfd::Index face = 0; face < built.mesh.face_count(); ++face)
    {
        if (built.mesh.face_centers()[face].y == 0.0)
        {
            const auto &area{geometry.face_area_vectors()[face]};
            require(area.x == 0.0 && area.y == 0.0, "Axis face must have exactly zero revolved area.");
            ++axis_faces;
        }
        else
        {
            const auto &area{geometry.face_area_vectors()[face]};
            require(std::hypot(area.x, area.y) > 0.0, "Non-axis face must have positive revolved area.");
        }
    }
    require(axis_faces == 1, "Expected exactly one axis face.");
}

void test_trapezoidal_quad_is_a_conical_frustum()
{
    auto raw{cfd::test::make_single_quadrilateral_raw_mesh()};
    raw.nodes = {{0.0, 0.0}, {2.0, 0.0}, {2.0, 3.0}, {0.0, 1.0}};
    const auto built{cfd::build_mesh(std::move(raw))};
    const cfd::AxisymmetricGeometry geometry{built.mesh};
    // Frustum: pi * height * (R^2 + R*r + r^2) / 3.
    // The polygon centroid is not the arithmetic mean of its vertices.
    require_near(geometry.cell_volumes()[0], 26.0 * pi / 3.0, tolerance, "Incorrect conical-frustum volume.");
    require_face_vector(built.mesh, geometry, 2, 3, {-8.0 * pi, 8.0 * pi});
    require_cell_identities(built.mesh, geometry, std::array{8.0 * pi});
}

void test_negative_radius_is_rejected()
{
    auto raw{cfd::test::make_single_quadrilateral_raw_mesh()};
    raw.nodes[0].y = -1.0e-12;
    const auto built{cfd::build_mesh(std::move(raw))};
    require(built.mesh.cell_centers()[0].y > 0.0, "Regression requires a positive centroid despite a negative node.");
    require_throws_with_message<std::runtime_error>(
        [&] { const cfd::AxisymmetricGeometry rejected{built.mesh}; }, "node 0",
        "A negative nodal radius must be rejected explicitly, not clamped.");
}

cfd::RawMeshData mixed_mesh()
{
    cfd::RawMeshData raw;
    raw.nodes = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}, {2.0, 0.0}, {2.0, 1.0}};
    raw.cell_types = {cfd::CellType::Quadrilateral, cfd::CellType::Triangle, cfd::CellType::Triangle};
    raw.cell_nodes = {0, 1, 2, 3, 1, 4, 5, 1, 5, 2};
    raw.cell_node_offsets = {0, 4, 7, 10};
    raw.boundary_groups = {{0, "wall"}};
    raw.boundary_edges = {{{0, 1}, 0}, {{1, 4}, 0}, {{4, 5}, 0}, {{5, 2}, 0}, {{2, 3}, 0}, {{3, 0}, 0}};
    return raw;
}

void test_mixed_mesh_and_internal_face_orientation()
{
    const auto built{cfd::build_mesh(mixed_mesh())};
    const auto &mesh{built.mesh};
    const cfd::AxisymmetricGeometry geometry{mesh};
    require_near(geometry.cell_volumes()[0], pi, tolerance, "Incorrect QUAD cylinder volume.");
    require_near(geometry.cell_volumes()[1], pi / 3.0, tolerance, "Incorrect TRI cone volume.");
    require_near(geometry.cell_volumes()[2], 2.0 * pi / 3.0, tolerance, "Incorrect cylinder-minus-cone volume.");
    require_cell_identities(mesh, geometry, std::array{2.0 * pi, pi, pi});

    cfd::Index internal_faces{};
    double boundary_axial{}, boundary_radial{};
    for (cfd::Index face = 0; face < mesh.face_count(); ++face)
    {
        const auto &adjacency{mesh.face_adjacencies()[face]};
        const auto &area{geometry.face_area_vectors()[face]};
        if (adjacency.is_boundary())
        {
            boundary_axial += area.x;
            boundary_radial += area.y;
            continue;
        }
        ++internal_faces;
        const auto &owner{mesh.cell_centers()[adjacency.owner]};
        const auto &neighbor{mesh.cell_centers()[adjacency.neighbor]};
        require(area.x * (neighbor.x - owner.x) + area.y * (neighbor.y - owner.y) > 0.0,
                "Internal revolved area must point from owner toward neighbor.");
        const double owner_flux{2.0 * area.x + 3.0 * area.y};
        const double neighbor_flux{-2.0 * area.x - 3.0 * area.y};
        require(owner_flux + neighbor_flux == 0.0, "An internal face must cancel exactly between its cells.");
    }
    require(internal_faces == 2, "Expected QUAD/TRI and TRI/TRI interfaces.");
    require_near(boundary_axial, 0.0, tolerance, "Global axial closure failed.");
    require_near(boundary_radial, 4.0 * pi, tolerance, "Internal faces did not cancel globally.");
}

void test_clockwise_cells_and_move_construction()
{
    const auto built{cfd::build_mesh(cfd::test::make_single_clockwise_quadrilateral_raw_mesh())};
    cfd::AxisymmetricGeometry original{built.mesh};
    const cfd::AxisymmetricGeometry geometry{std::move(original)};
    require(&geometry.mesh() == &built.mesh, "Move construction lost Mesh identity.");
    require_near(geometry.cell_volumes()[0], pi, tolerance, "Clockwise winding changed the revolved volume.");
    require_face_vector(built.mesh, geometry, 0, 1, {0.0, 0.0});
    require_face_vector(built.mesh, geometry, 2, 3, {0.0, 2.0 * pi});
    require_cell_identities(built.mesh, geometry, std::array{2.0 * pi});
}

void test_small_and_large_geometric_scales()
{
    for (const double scale : {1.0e-50, 1.0e-6, 1.0, 1.0e6, 1.0e50})
    {
        const auto built{cfd::build_mesh(rectangle(scale, 0.0, scale))};
        const cfd::AxisymmetricGeometry geometry{built.mesh};
        // Normalize the comparison explicitly: TestUtils has a unit-sized floor.
        require_near(geometry.cell_volumes()[0] / (pi * scale * scale * scale), 1.0, tolerance,
                     "Scaled cylinder volume failed.");
        require_cell_identities(built.mesh, geometry, std::array{2.0 * pi * scale * scale});
    }
    const auto built{cfd::build_mesh(rectangle(2.0, 1.0e8, 1.0e8 + 1.0))};
    const cfd::AxisymmetricGeometry geometry{built.mesh};
    // Expand R_outer^2 - R_inner^2 analytically to avoid subtracting huge squares.
    require_near(geometry.cell_volumes()[0] / (2.0 * pi * (2.0e8 + 1.0)), 1.0, tolerance,
                 "Large radial offset changed the annular volume.");
    require_cell_identities(built.mesh, geometry, std::array{4.0 * pi});
}

void test_overflowing_revolved_volume_is_rejected()
{
    const auto built{cfd::build_mesh(rectangle(1.0e100, 1.0e110, 1.0e110 + 1.0e100))};
    require(std::isfinite(built.mesh.cell_areas()[0]), "The planar mesh must be valid before revolution.");
    require_throws_with_message<std::runtime_error>([&] { const cfd::AxisymmetricGeometry rejected{built.mesh}; },
                                                    "cell 0", "A non-representable revolved volume must be rejected.");
}

void test_non_finite_coordinates_are_rejected_by_mesh_pipeline()
{
    for (const double invalid : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
                                 -std::numeric_limits<double>::infinity()})
    {
        auto raw{cfd::test::make_single_quadrilateral_raw_mesh()};
        raw.nodes[0].y = invalid;
        require_throws_with_message<std::runtime_error>(
            [&] {
                const auto built{cfd::build_mesh(std::move(raw))};
                const cfd::AxisymmetricGeometry rejected{built.mesh};
            },
            "non-finite coordinate", "Non-finite coordinates must not produce axisymmetric metrics.");
    }
}

} // namespace

int main()
{
    int failures{};
    failures += cfd::test::run_test("cylinder volume and face integrals", test_cylinder_volume_and_face_integrals);
    failures += cfd::test::run_test("annular cylinder", test_annular_cylinder);
    failures += cfd::test::run_test("triangle cone integral", test_triangle_is_a_cone);
    failures += cfd::test::run_test("inclined quadrilateral integral", test_inclined_quadrilateral);
    failures += cfd::test::run_test("trapezoidal frustum integral", test_trapezoidal_quad_is_a_conical_frustum);
    failures += cfd::test::run_test("zero-area axis face", test_axis_face_is_exactly_zero);
    failures += cfd::test::run_test("negative radius rejection", test_negative_radius_is_rejected);
    failures +=
        cfd::test::run_test("mixed mesh identities and orientation", test_mixed_mesh_and_internal_face_orientation);
    failures +=
        cfd::test::run_test("clockwise cells and move construction", test_clockwise_cells_and_move_construction);
    failures += cfd::test::run_test("geometric scaling", test_small_and_large_geometric_scales);
    failures += cfd::test::run_test("volume overflow rejection", test_overflowing_revolved_volume_is_rejected);
    failures += cfd::test::run_test("non-finite coordinate rejection",
                                    test_non_finite_coordinates_are_rejected_by_mesh_pipeline);
    return cfd::test::finish_tests(failures, "axisymmetric geometry");
}
