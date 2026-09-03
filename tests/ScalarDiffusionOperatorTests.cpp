#include "cfd/numerics/ScalarDiffusionOperator.hpp"

#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVectorField.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/math/Point2.hpp"
#include "cfd/math/Vector2.hpp"
#include "cfd/mesh/Boundary.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/meshing/RawMeshData.hpp"

#include "support/MeshFixtures.hpp"
#include "support/TestUtils.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace
{

using cfd::test::make_single_quadrilateral_raw_mesh;
using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws;
using cfd::test::test_tolerance;

static_assert(!std::is_copy_constructible_v<cfd::ScalarDiffusionOperator>);
static_assert(!std::is_copy_assignable_v<cfd::ScalarDiffusionOperator>);
static_assert(std::is_nothrow_move_constructible_v<cfd::ScalarDiffusionOperator>);
static_assert(!std::is_move_assignable_v<cfd::ScalarDiffusionOperator>);

[[nodiscard]]
cfd::RawMeshData make_two_cell_rectangle_raw_mesh()
{
    constexpr cfd::BoundaryId wall_boundary_id{0};

    cfd::RawMeshData raw_mesh;
    raw_mesh.nodes = {
        {0.0, 0.0}, {1.0, 0.0}, {2.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}, {2.0, 1.0},
    };
    raw_mesh.cell_types = {
        cfd::CellType::Quadrilateral,
        cfd::CellType::Quadrilateral,
    };
    raw_mesh.cell_nodes = {
        0, 1, 4, 3, 1, 2, 5, 4,
    };
    raw_mesh.cell_node_offsets = {0, 4, 8};
    raw_mesh.boundary_groups = {{wall_boundary_id, "wall"}};
    raw_mesh.boundary_edges = {
        {{0, 1}, wall_boundary_id}, {{1, 2}, wall_boundary_id}, {{2, 5}, wall_boundary_id},
        {{5, 4}, wall_boundary_id}, {{4, 3}, wall_boundary_id}, {{3, 0}, wall_boundary_id},
    };
    return raw_mesh;
}

[[nodiscard]]
cfd::RawMeshData make_two_cell_sheared_raw_mesh()
{
    cfd::RawMeshData raw_mesh{make_two_cell_rectangle_raw_mesh()};
    raw_mesh.nodes = {
        {0.0, 0.0}, {1.0, 0.0}, {2.0, 0.0}, {0.4, 1.0}, {1.4, 1.0}, {2.4, 1.0},
    };

    raw_mesh.boundary_groups.clear();
    raw_mesh.boundary_edges.clear();

    const auto add_boundary = [&raw_mesh](const cfd::Index first_node, const cfd::Index second_node, const char *name) {
        const cfd::BoundaryId boundary_id{raw_mesh.boundary_groups.size()};
        raw_mesh.boundary_groups.push_back({boundary_id, name});
        raw_mesh.boundary_edges.push_back({{first_node, second_node}, boundary_id});
    };

    add_boundary(0, 1, "bottom_0");
    add_boundary(1, 2, "bottom_1");
    add_boundary(2, 5, "right");
    add_boundary(5, 4, "top_1");
    add_boundary(4, 3, "top_0");
    add_boundary(3, 0, "left");
    return raw_mesh;
}

[[nodiscard]]
cfd::RawMeshData make_non_conjunctional_two_cell_raw_mesh()
{
    cfd::RawMeshData raw_mesh;
    raw_mesh.nodes = {
        {0.0, 0.0},
        {1.0, 0.0},
        {2.5, 1.0},
        {2.5, -1.0},
    };
    raw_mesh.cell_types = {
        cfd::CellType::Triangle,
        cfd::CellType::Triangle,
    };
    raw_mesh.cell_nodes = {
        0, 1, 2, 1, 0, 3,
    };
    raw_mesh.cell_node_offsets = {0, 3, 6};
    raw_mesh.boundary_groups = {
        {0, "upper_right"},
        {1, "upper_left"},
        {2, "lower_left"},
        {3, "lower_right"},
    };
    raw_mesh.boundary_edges = {
        {{1, 2}, 0},
        {{2, 0}, 1},
        {{0, 3}, 2},
        {{3, 1}, 3},
    };
    return raw_mesh;
}

[[nodiscard]]
cfd::RawMeshData make_four_boundary_rectangle_raw_mesh(const double width = 1.0)
{
    constexpr cfd::BoundaryId left_boundary_id{0};
    constexpr cfd::BoundaryId right_boundary_id{1};
    constexpr cfd::BoundaryId bottom_boundary_id{2};
    constexpr cfd::BoundaryId top_boundary_id{3};

    cfd::RawMeshData raw_mesh{make_single_quadrilateral_raw_mesh()};
    raw_mesh.nodes[1].x = width;
    raw_mesh.nodes[2].x = width;
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
cfd::ScalarBoundaryConditions make_uniform_neumann_conditions(const cfd::Index boundary_count,
                                                              const double derivative = 0.0)
{
    std::vector<cfd::ScalarBoundaryCondition> conditions;
    conditions.reserve(boundary_count);
    for (cfd::Index boundary_id = 0; boundary_id < boundary_count; ++boundary_id)
    {
        conditions.emplace_back(cfd::ScalarBoundaryConditionType::Neumann, derivative);
    }
    return {boundary_count, std::move(conditions)};
}

[[nodiscard]]
double linear_value(const cfd::Point2 &point) noexcept
{
    return 2.0 * point.x - 3.0 * point.y + 1.5;
}

[[nodiscard]]
cfd::ScalarBoundaryConditions make_exact_linear_dirichlet_conditions(const cfd::Mesh &mesh)
{
    std::vector<cfd::ScalarBoundaryCondition> conditions;
    conditions.reserve(mesh.boundary_groups().size());
    for (cfd::Index boundary_id = 0; boundary_id < mesh.boundary_groups().size(); ++boundary_id)
    {
        conditions.emplace_back(cfd::ScalarBoundaryConditionType::Neumann, 0.0);
    }

    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (!mesh.face_adjacencies()[face_id].is_boundary())
        {
            continue;
        }
        const cfd::BoundaryId boundary_id{mesh.face_boundary_ids()[face_id]};
        conditions[boundary_id] = {
            cfd::ScalarBoundaryConditionType::Dirichlet,
            linear_value(mesh.face_centers()[face_id]),
        };
    }

    return {mesh.boundary_groups().size(), std::move(conditions)};
}

void test_constructor_rejects_invalid_diffusivity()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_four_boundary_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};

    const auto require_rejected = [&mesh](const double diffusivity) {
        require_throws<std::invalid_argument>(
            [&mesh, diffusivity]() { const cfd::ScalarDiffusionOperator diffusion{mesh, diffusivity}; },
            "Scalar diffusion accepted an invalid diffusivity.");
    };

    require_rejected(0.0);
    require_rejected(-1.0);
    require_rejected(std::numeric_limits<double>::quiet_NaN());
    require_rejected(std::numeric_limits<double>::infinity());
    require_rejected(-std::numeric_limits<double>::infinity());
}

void test_accepts_non_conjunctional_internal_face()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_non_conjunctional_two_cell_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarDiffusionOperator diffusion{mesh, 1.0};

    // The centroid line is x = 7/6 and intersects the supporting line y = 0
    // beyond the shared edge endpoint x = 1.
    const double supporting_line_intersection_x{0.5 * (mesh.cell_centers()[0].x + mesh.cell_centers()[1].x)};
    require(supporting_line_intersection_x > 1.0,
            "The non-conjunctional test fixture does not intersect beyond the shared face.");

    cfd::CellScalarField field{mesh.cell_count()};
    field[0] = 1.2;
    field[1] = -0.7;
    cfd::CellVectorField gradient{mesh.cell_count()};
    gradient[0] = {0.4, -1.1};
    gradient[1] = {-2.0, 0.3};
    const cfd::ScalarBoundaryConditions homogeneous_neumann{make_uniform_neumann_conditions(4)};
    cfd::CellScalarField flux_balance{mesh.cell_count()};

    diffusion.compute_flux_balance(field, homogeneous_neumann, gradient, flux_balance);

    require(std::isfinite(flux_balance[0]) && std::isfinite(flux_balance[1]),
            "A non-conjunctional face produced a non-finite flux.");
    require(std::abs(flux_balance[0]) > test_tolerance,
            "The non-conjunctional conservation test produced a trivial internal flux.");
    require_near(flux_balance[0] + flux_balance[1], 0.0, test_tolerance,
                 "A non-conjunctional internal face is not exactly conservative.");

    const cfd::CellVectorField exact_gradient{mesh.cell_count(), {2.0, -3.0}};
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        field[cell_id] = linear_value(mesh.cell_centers()[cell_id]);
    }
    const cfd::ScalarBoundaryConditions exact_dirichlet{make_exact_linear_dirichlet_conditions(mesh)};

    diffusion.compute_flux_balance(field, exact_dirichlet, exact_gradient, flux_balance);

    for (const double balance : flux_balance.values())
    {
        require_near(balance, 0.0, test_tolerance,
                     "Corrected diffusion is not linearly exact across a non-conjunctional face.");
    }
}

void test_constant_field_with_homogeneous_neumann_boundaries()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarDiffusionOperator diffusion{mesh, 1.75};
    const cfd::CellScalarField field{mesh.cell_count(), 4.25};
    const cfd::CellVectorField gradient{mesh.cell_count(), {0.0, 0.0}};
    const cfd::ScalarBoundaryConditions boundary_conditions{make_uniform_neumann_conditions(1)};
    cfd::CellScalarField flux_balance{mesh.cell_count(), 99.0};

    diffusion.compute_flux_balance(field, boundary_conditions, gradient, flux_balance);

    for (const double balance : flux_balance.values())
    {
        require_near(balance, 0.0, test_tolerance,
                     "A constant field with homogeneous Neumann data produced diffusion flux.");
    }
}

void test_orthogonal_internal_face_flux()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarDiffusionOperator diffusion{mesh, 2.0};
    cfd::CellScalarField field{mesh.cell_count()};
    field[0] = 5.0;
    field[1] = 2.0;
    const cfd::CellVectorField gradient{mesh.cell_count(), {7.0, -11.0}};
    const cfd::ScalarBoundaryConditions boundary_conditions{make_uniform_neumann_conditions(1)};
    cfd::CellScalarField flux_balance{mesh.cell_count()};

    diffusion.compute_flux_balance(field, boundary_conditions, gradient, flux_balance);

    require_near(flux_balance[0], 6.0, test_tolerance, "Orthogonal owner flux has an incorrect value.");
    require_near(flux_balance[1], -6.0, test_tolerance, "Orthogonal neighbor flux has an incorrect value.");
}

void test_internal_face_conservation_on_non_orthogonal_mesh()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_sheared_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarDiffusionOperator diffusion{mesh, 0.7};
    cfd::CellScalarField field{mesh.cell_count()};
    field[0] = -1.2;
    field[1] = 3.4;
    cfd::CellVectorField gradient{mesh.cell_count()};
    gradient[0] = {0.8, -2.1};
    gradient[1] = {-1.3, 4.2};
    const cfd::ScalarBoundaryConditions boundary_conditions{
        make_uniform_neumann_conditions(mesh.boundary_groups().size())};
    cfd::CellScalarField flux_balance{mesh.cell_count()};

    diffusion.compute_flux_balance(field, boundary_conditions, gradient, flux_balance);

    require(std::abs(flux_balance[0]) > test_tolerance,
            "Non-orthogonal conservation test produced a trivial internal flux.");
    require_near(flux_balance[0] + flux_balance[1], 0.0, test_tolerance,
                 "Internal-face contributions are not exactly conservative.");
}

void test_linear_exactness_on_sheared_quadrilaterals()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_sheared_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarDiffusionOperator diffusion{mesh, 1.3};
    cfd::CellScalarField field{mesh.cell_count()};
    const cfd::CellVectorField gradient{mesh.cell_count(), {2.0, -3.0}};
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        field[cell_id] = linear_value(mesh.cell_centers()[cell_id]);
    }
    const cfd::ScalarBoundaryConditions boundary_conditions{make_exact_linear_dirichlet_conditions(mesh)};
    cfd::CellScalarField flux_balance{mesh.cell_count()};

    diffusion.compute_flux_balance(field, boundary_conditions, gradient, flux_balance);

    for (const double balance : flux_balance.values())
    {
        require_near(balance, 0.0, test_tolerance,
                     "Corrected diffusion is not linearly exact on a sheared quadrilateral mesh.");
    }
}

void test_dirichlet_flux_sign_and_magnitude()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_four_boundary_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarDiffusionOperator diffusion{mesh, 2.0};
    const cfd::CellScalarField field{mesh.cell_count(), 0.5};
    const cfd::CellVectorField gradient{mesh.cell_count(), {1.0, 0.0}};
    std::vector<cfd::ScalarBoundaryCondition> conditions{
        {cfd::ScalarBoundaryConditionType::Neumann, 0.0},
        {cfd::ScalarBoundaryConditionType::Dirichlet, 1.0},
        {cfd::ScalarBoundaryConditionType::Neumann, 0.0},
        {cfd::ScalarBoundaryConditionType::Neumann, 0.0},
    };
    const cfd::ScalarBoundaryConditions boundary_conditions{4, std::move(conditions)};
    cfd::CellScalarField flux_balance{mesh.cell_count()};

    diffusion.compute_flux_balance(field, boundary_conditions, gradient, flux_balance);

    require_near(flux_balance[0], -2.0, test_tolerance,
                 "Positive outward Dirichlet gradient produced the wrong diffusion-flux sign or magnitude.");
}

void test_neumann_flux_sign_magnitude_and_zero()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_four_boundary_rectangle_raw_mesh(2.0))};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarDiffusionOperator diffusion{mesh, 3.0};
    const cfd::CellScalarField field{mesh.cell_count(), 8.0};
    const cfd::CellVectorField gradient{mesh.cell_count(), {-7.0, 5.0}};
    std::vector<cfd::ScalarBoundaryCondition> conditions{
        {cfd::ScalarBoundaryConditionType::Neumann, 0.0},
        {cfd::ScalarBoundaryConditionType::Neumann, 0.0},
        {cfd::ScalarBoundaryConditionType::Neumann, 0.0},
        {cfd::ScalarBoundaryConditionType::Neumann, 2.0},
    };
    const cfd::ScalarBoundaryConditions boundary_conditions{4, std::move(conditions)};
    cfd::CellScalarField flux_balance{mesh.cell_count()};

    diffusion.compute_flux_balance(field, boundary_conditions, gradient, flux_balance);
    require_near(flux_balance[0], -12.0, test_tolerance,
                 "Neumann derivative was not converted to outward diffusion flux with face length.");

    const cfd::ScalarBoundaryConditions zero_conditions{make_uniform_neumann_conditions(4)};
    diffusion.compute_flux_balance(field, zero_conditions, gradient, flux_balance);
    require_near(flux_balance[0], 0.0, test_tolerance, "Homogeneous Neumann data did not produce exactly zero flux.");
}

void test_rejects_cardinality_mismatches_and_aliasing()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarDiffusionOperator diffusion{mesh, 1.0};
    cfd::CellScalarField field{mesh.cell_count()};
    const cfd::CellVectorField gradient{mesh.cell_count()};
    const cfd::ScalarBoundaryConditions boundary_conditions{make_uniform_neumann_conditions(1)};
    cfd::CellScalarField flux_balance{mesh.cell_count()};

    require_throws<std::invalid_argument>(
        [&diffusion, &boundary_conditions, &gradient, &flux_balance]() {
            const cfd::CellScalarField wrong_field{3};
            diffusion.compute_flux_balance(wrong_field, boundary_conditions, gradient, flux_balance);
        },
        "Scalar diffusion accepted an input field with incorrect cardinality.");

    require_throws<std::invalid_argument>(
        [&diffusion, &field, &boundary_conditions, &flux_balance]() {
            const cfd::CellVectorField wrong_gradient{3};
            diffusion.compute_flux_balance(field, boundary_conditions, wrong_gradient, flux_balance);
        },
        "Scalar diffusion accepted a gradient with incorrect cardinality.");

    require_throws<std::invalid_argument>(
        [&diffusion, &field, &boundary_conditions, &gradient]() {
            cfd::CellScalarField wrong_output{3};
            diffusion.compute_flux_balance(field, boundary_conditions, gradient, wrong_output);
        },
        "Scalar diffusion accepted an output field with incorrect cardinality.");

    require_throws<std::invalid_argument>(
        [&diffusion, &field, &gradient, &flux_balance]() {
            const cfd::ScalarBoundaryConditions wrong_conditions{make_uniform_neumann_conditions(2)};
            diffusion.compute_flux_balance(field, wrong_conditions, gradient, flux_balance);
        },
        "Scalar diffusion accepted boundary conditions with incorrect cardinality.");

    require_throws<std::invalid_argument>(
        [&diffusion, &field, &boundary_conditions, &gradient]() {
            diffusion.compute_flux_balance(field, boundary_conditions, gradient, field);
        },
        "Scalar diffusion accepted an output alias of its input field.");
}

void test_does_not_mutate_inputs()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_sheared_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::ScalarDiffusionOperator diffusion{mesh, 1.0};
    cfd::CellScalarField field{mesh.cell_count()};
    cfd::CellVectorField gradient{mesh.cell_count()};
    field[0] = 1.25;
    field[1] = -0.75;
    gradient[0] = {2.0, -1.0};
    gradient[1] = {-3.0, 4.0};
    const cfd::ScalarBoundaryConditions boundary_conditions{make_exact_linear_dirichlet_conditions(mesh)};
    cfd::CellScalarField flux_balance{mesh.cell_count()};

    diffusion.compute_flux_balance(field, boundary_conditions, gradient, flux_balance);

    require_near(field[0], 1.25, 0.0, "Scalar diffusion modified its first input value.");
    require_near(field[1], -0.75, 0.0, "Scalar diffusion modified its second input value.");
    require_near(gradient[0].x, 2.0, 0.0, "Scalar diffusion modified its first x-gradient.");
    require_near(gradient[0].y, -1.0, 0.0, "Scalar diffusion modified its first y-gradient.");
    require_near(gradient[1].x, -3.0, 0.0, "Scalar diffusion modified its second x-gradient.");
    require_near(gradient[1].y, 4.0, 0.0, "Scalar diffusion modified its second y-gradient.");
}

} // namespace

int main()
{
    int failure_count{};

    failure_count +=
        cfd::test::run_test("scalar diffusion constructor validation", test_constructor_rejects_invalid_diffusivity);
    failure_count += cfd::test::run_test("scalar diffusion non-conjunctional face acceptance",
                                         test_accepts_non_conjunctional_internal_face);
    failure_count += cfd::test::run_test("scalar diffusion of a constant field",
                                         test_constant_field_with_homogeneous_neumann_boundaries);
    failure_count +=
        cfd::test::run_test("scalar diffusion orthogonal internal face", test_orthogonal_internal_face_flux);
    failure_count += cfd::test::run_test("scalar diffusion internal-face conservation",
                                         test_internal_face_conservation_on_non_orthogonal_mesh);
    failure_count += cfd::test::run_test("scalar diffusion linear exactness on sheared quads",
                                         test_linear_exactness_on_sheared_quadrilaterals);
    failure_count += cfd::test::run_test("scalar diffusion Dirichlet sign", test_dirichlet_flux_sign_and_magnitude);
    failure_count += cfd::test::run_test("scalar diffusion Neumann sign", test_neumann_flux_sign_magnitude_and_zero);
    failure_count += cfd::test::run_test("scalar diffusion cardinality and alias validation",
                                         test_rejects_cardinality_mismatches_and_aliasing);
    failure_count += cfd::test::run_test("scalar diffusion input immutability", test_does_not_mutate_inputs);

    return cfd::test::finish_tests(failure_count, "scalar diffusion operator");
}
