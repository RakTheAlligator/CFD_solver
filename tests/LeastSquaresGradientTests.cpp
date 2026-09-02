#include "cfd/numerics/LeastSquaresGradient.hpp"

#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVectorField.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/math/Point2.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/meshing/RawMeshData.hpp"

#include "support/MeshFixtures.hpp"
#include "support/TestUtils.hpp"

#include <stdexcept>
#include <utility>
#include <vector>

namespace
{

using cfd::test::make_single_quadrilateral_raw_mesh;
using cfd::test::make_single_triangle_raw_mesh;
using cfd::test::make_two_triangle_raw_mesh;
using cfd::test::require_near;
using cfd::test::require_throws;
using cfd::test::test_tolerance;

constexpr cfd::BoundaryId left_boundary_id{0};
constexpr cfd::BoundaryId right_boundary_id{1};
constexpr cfd::BoundaryId bottom_boundary_id{2};
constexpr cfd::BoundaryId top_boundary_id{3};

[[nodiscard]]
cfd::RawMeshData make_four_boundary_quadrilateral_raw_mesh()
{
    cfd::RawMeshData raw_mesh{make_single_quadrilateral_raw_mesh()};

    raw_mesh.boundary_groups = {
        {left_boundary_id, "left"},
        {right_boundary_id, "right"},
        {bottom_boundary_id, "bottom"},
        {top_boundary_id, "top"},
    };

    raw_mesh.boundary_edges = {
        {{3, 0}, left_boundary_id},
        {{1, 2}, right_boundary_id},
        {{0, 1}, bottom_boundary_id},
        {{2, 3}, top_boundary_id},
    };

    return raw_mesh;
}

[[nodiscard]]
cfd::RawMeshData make_skewed_two_triangle_raw_mesh()
{
    cfd::RawMeshData raw_mesh{make_two_triangle_raw_mesh()};

    raw_mesh.nodes = {
        {0.0, 0.0},
        {2.0, 0.2},
        {1.6, 1.4},
        {-0.2, 1.0},
    };

    raw_mesh.boundary_groups = {
        {left_boundary_id, "left"},
        {right_boundary_id, "right"},
        {bottom_boundary_id, "bottom"},
        {top_boundary_id, "top"},
    };

    raw_mesh.boundary_edges = {
        {{3, 0}, left_boundary_id},
        {{1, 2}, right_boundary_id},
        {{0, 1}, bottom_boundary_id},
        {{2, 3}, top_boundary_id},
    };

    return raw_mesh;
}

[[nodiscard]]
cfd::RawMeshData make_inverse_distance_weighting_triangle_raw_mesh()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};

    raw_mesh.nodes = {
        {-1.0, 0.0},
        {1.0, 0.0},
        {0.0, 3.0},
    };

    raw_mesh.boundary_groups = {
        {left_boundary_id, "left"},
        {right_boundary_id, "right"},
        {bottom_boundary_id, "bottom"},
    };

    raw_mesh.boundary_edges = {
        {{0, 1}, bottom_boundary_id},
        {{1, 2}, right_boundary_id},
        {{2, 0}, left_boundary_id},
    };

    return raw_mesh;
}

[[nodiscard]]
cfd::RawMeshData make_two_by_one_rectangle_raw_mesh()
{
    cfd::RawMeshData raw_mesh{make_four_boundary_quadrilateral_raw_mesh()};

    raw_mesh.nodes = {
        {0.0, 0.0},
        {2.0, 0.0},
        {2.0, 1.0},
        {0.0, 1.0},
    };

    return raw_mesh;
}

[[nodiscard]]
double linear_value(const cfd::Point2 &point) noexcept
{
    return 2.0 * point.x + 3.0 * point.y + 4.0;
}

[[nodiscard]]
cfd::ScalarBoundaryConditions make_constant_mixed_boundary_conditions(const double value)
{
    std::vector<cfd::ScalarBoundaryCondition> conditions{
        {cfd::ScalarBoundaryConditionType::Dirichlet, value},
        {cfd::ScalarBoundaryConditionType::Dirichlet, value},
        {cfd::ScalarBoundaryConditionType::Neumann, 0.0},
        {cfd::ScalarBoundaryConditionType::Neumann, 0.0},
    };

    return {4, std::move(conditions)};
}

[[nodiscard]]
cfd::ScalarBoundaryConditions make_linear_mixed_boundary_conditions()
{
    std::vector<cfd::ScalarBoundaryCondition> conditions{
        {cfd::ScalarBoundaryConditionType::Dirichlet, linear_value({0.0, 0.5})},
        {cfd::ScalarBoundaryConditionType::Dirichlet, linear_value({1.0, 0.5})},
        {cfd::ScalarBoundaryConditionType::Neumann, -3.0},
        {cfd::ScalarBoundaryConditionType::Neumann, 3.0},
    };

    return {4, std::move(conditions)};
}

[[nodiscard]]
cfd::ScalarBoundaryConditions make_skewed_dirichlet_boundary_conditions()
{
    std::vector<cfd::ScalarBoundaryCondition> conditions{
        {cfd::ScalarBoundaryConditionType::Dirichlet, linear_value({-0.1, 0.5})},
        {cfd::ScalarBoundaryConditionType::Dirichlet, linear_value({1.8, 0.8})},
        {cfd::ScalarBoundaryConditionType::Dirichlet, linear_value({1.0, 0.1})},
        {cfd::ScalarBoundaryConditionType::Dirichlet, linear_value({0.7, 1.2})},
    };

    return {4, std::move(conditions)};
}

[[nodiscard]]
cfd::ScalarBoundaryConditions make_two_by_one_rectangle_boundary_conditions()
{
    std::vector<cfd::ScalarBoundaryCondition> conditions{
        {cfd::ScalarBoundaryConditionType::Dirichlet, linear_value({0.0, 0.5})},
        {cfd::ScalarBoundaryConditionType::Dirichlet, linear_value({2.0, 0.5})},
        {cfd::ScalarBoundaryConditionType::Neumann, -3.0},
        {cfd::ScalarBoundaryConditionType::Neumann, 3.0},
    };

    return {4, std::move(conditions)};
}

void test_constant_field_on_quadrilateral()
{
    constexpr double constant_value{5.25};

    cfd::MeshBuildResult build_result{cfd::build_mesh(make_four_boundary_quadrilateral_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};

    const cfd::CellScalarField field{mesh.cell_count(), constant_value};
    const cfd::ScalarBoundaryConditions boundary_conditions{make_constant_mixed_boundary_conditions(constant_value)};
    cfd::CellVectorField gradient{mesh.cell_count(), {11.0, -7.0}};

    cfd::compute_least_squares_gradient(mesh, field, boundary_conditions, gradient);

    require_near(gradient[0].x, 0.0, test_tolerance, "Constant field produced a non-zero x-gradient.");
    require_near(gradient[0].y, 0.0, test_tolerance, "Constant field produced a non-zero y-gradient.");
}

void test_mixed_boundaries_reconstruct_linear_field_on_quadrilateral()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_four_boundary_quadrilateral_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};

    cfd::CellScalarField field{mesh.cell_count()};
    field[0] = linear_value(mesh.cell_centers()[0]);

    const cfd::ScalarBoundaryConditions boundary_conditions{make_linear_mixed_boundary_conditions()};
    cfd::CellVectorField gradient{mesh.cell_count()};

    cfd::compute_least_squares_gradient(mesh, field, boundary_conditions, gradient);

    require_near(gradient[0].x, 2.0, test_tolerance,
                 "Mixed-boundary quadrilateral reconstruction produced an incorrect x-gradient.");
    require_near(gradient[0].y, 3.0, test_tolerance,
                 "Mixed-boundary quadrilateral reconstruction produced an incorrect y-gradient.");
}

void test_linear_exactness_on_skewed_triangles()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_skewed_two_triangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};

    cfd::CellScalarField field{mesh.cell_count()};

    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        field[cell_id] = linear_value(mesh.cell_centers()[cell_id]);
    }

    const cfd::ScalarBoundaryConditions boundary_conditions{make_skewed_dirichlet_boundary_conditions()};
    cfd::CellVectorField gradient{mesh.cell_count()};

    cfd::compute_least_squares_gradient(mesh, field, boundary_conditions, gradient);

    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        require_near(gradient[cell_id].x, 2.0, test_tolerance,
                     "Skewed-triangle reconstruction produced an incorrect x-gradient.");
        require_near(gradient[cell_id].y, 3.0, test_tolerance,
                     "Skewed-triangle reconstruction produced an incorrect y-gradient.");
    }
}

void test_inverse_distance_weighting_on_triangle()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_inverse_distance_weighting_triangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};

    const cfd::CellScalarField field{mesh.cell_count(), 0.0};

    // These deliberately inconsistent observations distinguish inverse-distance
    // weighting (gy = 1/2) from unweighted least squares (gy = 1/3).
    std::vector<cfd::ScalarBoundaryCondition> conditions{
        {cfd::ScalarBoundaryConditionType::Dirichlet, 0.0},
        {cfd::ScalarBoundaryConditionType::Dirichlet, 1.0},
        {cfd::ScalarBoundaryConditionType::Dirichlet, 0.0},
    };
    const cfd::ScalarBoundaryConditions boundary_conditions{3, std::move(conditions)};
    cfd::CellVectorField gradient{mesh.cell_count()};

    cfd::compute_least_squares_gradient(mesh, field, boundary_conditions, gradient);

    require_near(gradient[0].x, 1.0, test_tolerance,
                 "Triangle reconstruction did not use inverse-distance weighting for the x-gradient.");
    require_near(gradient[0].y, 0.5, test_tolerance,
                 "Triangle reconstruction did not use inverse-distance weighting for the y-gradient.");
}

void test_neumann_uses_unit_normal_on_two_by_one_rectangle()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_by_one_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};

    cfd::CellScalarField field{mesh.cell_count()};
    field[0] = linear_value(mesh.cell_centers()[0]);

    const cfd::ScalarBoundaryConditions boundary_conditions{make_two_by_one_rectangle_boundary_conditions()};
    cfd::CellVectorField gradient{mesh.cell_count()};

    cfd::compute_least_squares_gradient(mesh, field, boundary_conditions, gradient);

    require_near(gradient[0].x, 2.0, test_tolerance,
                 "Two-by-one rectangle reconstruction produced an incorrect x-gradient.");
    require_near(gradient[0].y, 3.0, test_tolerance, "Neumann reconstruction did not normalize the face-area vector.");
}

void test_rejects_cardinality_mismatches()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_four_boundary_quadrilateral_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};

    const cfd::CellScalarField field{mesh.cell_count()};
    const cfd::ScalarBoundaryConditions boundary_conditions{make_constant_mixed_boundary_conditions(0.0)};
    cfd::CellVectorField gradient{mesh.cell_count()};

    require_throws<std::invalid_argument>(
        [&mesh, &boundary_conditions, &gradient]() {
            const cfd::CellScalarField wrong_field{mesh.cell_count() + 1};
            cfd::compute_least_squares_gradient(mesh, wrong_field, boundary_conditions, gradient);
        },
        "Least-squares gradient accepted a scalar field with incorrect cardinality.");

    require_throws<std::invalid_argument>(
        [&mesh, &field, &boundary_conditions]() {
            cfd::CellVectorField wrong_gradient{mesh.cell_count() + 1};
            cfd::compute_least_squares_gradient(mesh, field, boundary_conditions, wrong_gradient);
        },
        "Least-squares gradient accepted an output field with incorrect cardinality.");

    require_throws<std::invalid_argument>(
        [&mesh, &field, &gradient]() {
            std::vector<cfd::ScalarBoundaryCondition> conditions{
                {cfd::ScalarBoundaryConditionType::Dirichlet, 0.0},
            };
            const cfd::ScalarBoundaryConditions wrong_boundary_conditions{1, std::move(conditions)};

            cfd::compute_least_squares_gradient(mesh, field, wrong_boundary_conditions, gradient);
        },
        "Least-squares gradient accepted boundary conditions with incorrect cardinality.");
}

} // namespace

int main()
{
    int failure_count{};

    failure_count +=
        cfd::test::run_test("least-squares gradient of a constant field", test_constant_field_on_quadrilateral);
    failure_count += cfd::test::run_test("least-squares gradient with mixed boundary conditions",
                                         test_mixed_boundaries_reconstruct_linear_field_on_quadrilateral);
    failure_count +=
        cfd::test::run_test("least-squares gradient on skewed triangles", test_linear_exactness_on_skewed_triangles);
    failure_count += cfd::test::run_test("least-squares gradient inverse-distance weighting",
                                         test_inverse_distance_weighting_on_triangle);
    failure_count += cfd::test::run_test("least-squares gradient Neumann unit normal",
                                         test_neumann_uses_unit_normal_on_two_by_one_rectangle);
    failure_count +=
        cfd::test::run_test("least-squares gradient cardinality validation", test_rejects_cardinality_mismatches);

    return cfd::test::finish_tests(failure_count, "least-squares gradient");
}
