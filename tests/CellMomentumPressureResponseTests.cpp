#include "cfd/field/CellMomentumPressureResponse.hpp"

#include "cfd/linear_algebra/ScalarLinearSystem.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/numerics/MomentumPressureResponse.hpp"

#include "support/MeshFixtures.hpp"
#include "support/TestUtils.hpp"

#include <array>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace
{

using cfd::test::make_single_quadrilateral_raw_mesh;
using cfd::test::make_two_triangle_raw_mesh;
using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws;

static_assert(!std::is_default_constructible_v<cfd::CellMomentumPressureResponse>);
static_assert(std::is_copy_constructible_v<cfd::CellMomentumPressureResponse>);
static_assert(!std::is_copy_assignable_v<cfd::CellMomentumPressureResponse>);
static_assert(std::is_nothrow_move_constructible_v<cfd::CellMomentumPressureResponse>);
static_assert(!std::is_move_assignable_v<cfd::CellMomentumPressureResponse>);

static_assert(
    std::is_same_v<decltype(std::declval<cfd::CellMomentumPressureResponse &>().u()), cfd::CellScalarField &>);
static_assert(std::is_same_v<decltype(std::declval<const cfd::CellMomentumPressureResponse &>().u()),
                             const cfd::CellScalarField &>);
static_assert(
    std::is_same_v<decltype(std::declval<cfd::CellMomentumPressureResponse &>().v()), cfd::CellScalarField &>);
static_assert(std::is_same_v<decltype(std::declval<const cfd::CellMomentumPressureResponse &>().v()),
                             const cfd::CellScalarField &>);

void require_response_unchanged(const cfd::CellMomentumPressureResponse &response, const double expected_u,
                                const double expected_v, const char *const context)
{
    for (cfd::Index cell_id = 0; cell_id < response.size(); ++cell_id)
    {
        require_near(response.u()[cell_id], expected_u, 0.0, context);
        require_near(response.v()[cell_id], expected_v, 0.0, context);
    }
}

void seed_system(cfd::ScalarLinearSystem &system, const double diagonal_value, const double state_seed)
{
    for (double &diagonal : system.diagonal())
    {
        diagonal = diagonal_value;
    }
    for (double &coefficient : system.owner_neighbor_coefficients())
    {
        coefficient = state_seed;
    }
    for (double &coefficient : system.neighbor_owner_coefficients())
    {
        coefficient = state_seed + 1.0;
    }
    for (double &rhs_value : system.rhs())
    {
        rhs_value = state_seed + 2.0;
    }
}

void test_field_construction_and_component_access()
{
    cfd::CellMomentumPressureResponse response{4};

    require(response.size() == 4, "Momentum pressure response has an incorrect cardinality.");
    require(response.u().size() == response.size(), "The u-response cardinality differs from the field.");
    require(response.v().size() == response.size(), "The v-response cardinality differs from the field.");
    for (cfd::Index cell_id = 0; cell_id < response.size(); ++cell_id)
    {
        require(response.u()[cell_id] == 0.0, "A u-response value was not initialized to zero.");
        require(response.v()[cell_id] == 0.0, "A v-response value was not initialized to zero.");
    }

    response.u()[1] = 2.5;
    response.v()[2] = -3.5;
    require(response.u()[1] == 2.5 && response.u()[2] == 0.0, "Mutable u-response access changed unrelated storage.");
    require(response.v()[1] == 0.0 && response.v()[2] == -3.5, "Mutable v-response access changed unrelated storage.");

    const cfd::CellMomentumPressureResponse &const_response{response};
    require(const_response.u()[1] == 2.5, "Const u-response access returned an incorrect value.");
    require(const_response.v()[2] == -3.5, "Const v-response access returned an incorrect value.");
}

void test_copy_construction_is_deep()
{
    cfd::CellMomentumPressureResponse original{3};
    original.u()[0] = 1.0;
    original.u()[1] = 2.0;
    original.v()[1] = -2.0;
    original.v()[2] = -3.0;

    cfd::CellMomentumPressureResponse copy{original};
    copy.u()[0] = 8.0;
    copy.v()[2] = -9.0;

    require(copy.u()[0] == 8.0 && copy.u()[1] == 2.0, "The copied u-response contains incorrect values.");
    require(copy.v()[1] == -2.0 && copy.v()[2] == -9.0, "The copied v-response contains incorrect values.");
    require(original.u()[0] == 1.0, "Modifying a copied u-response changed the original.");
    require(original.v()[2] == -3.0, "Modifying a copied v-response changed the original.");
}

void test_computes_one_cell_response_from_independent_diagonals()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_quadrilateral_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    cfd::ScalarLinearSystem u_system{mesh};
    cfd::ScalarLinearSystem v_system{mesh};
    constexpr double u_diagonal{2.0};
    constexpr double v_diagonal{5.0};
    u_system.diagonal()[0] = u_diagonal;
    v_system.diagonal()[0] = v_diagonal;
    cfd::CellMomentumPressureResponse response{mesh.cell_count()};

    cfd::compute_momentum_pressure_response(mesh, u_system, v_system, response);

    const double area{mesh.cell_areas()[0]};
    require_near(response.u()[0], area / u_diagonal, 0.0, "The one-cell u-momentum pressure response is incorrect.");
    require_near(response.v()[0], area / v_diagonal, 0.0, "The one-cell v-momentum pressure response is incorrect.");
    require(response.u()[0] != response.v()[0], "Different component diagonals produced equal responses.");
}

void test_multiple_cells_use_current_diagonals_and_overwrite_output()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_triangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    cfd::ScalarLinearSystem u_system{mesh};
    cfd::ScalarLinearSystem v_system{mesh};
    cfd::CellMomentumPressureResponse response{mesh.cell_count()};
    constexpr std::array first_u_diagonals{2.0, 4.0};
    constexpr std::array first_v_diagonals{5.0, 10.0};
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        u_system.diagonal()[cell_id] = first_u_diagonals.at(cell_id);
        v_system.diagonal()[cell_id] = first_v_diagonals.at(cell_id);
    }
    cfd::compute_momentum_pressure_response(mesh, u_system, v_system, response);
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        require_near(response.u()[cell_id], mesh.cell_areas()[cell_id] / first_u_diagonals.at(cell_id), 0.0,
                     "A multi-cell u response is incorrect.");
        require_near(response.v()[cell_id], mesh.cell_areas()[cell_id] / first_v_diagonals.at(cell_id), 0.0,
                     "A multi-cell v response is incorrect.");
    }

    constexpr std::array second_u_diagonals{8.0, 16.0};
    constexpr std::array second_v_diagonals{20.0, 40.0};
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        u_system.diagonal()[cell_id] = second_u_diagonals.at(cell_id);
        v_system.diagonal()[cell_id] = second_v_diagonals.at(cell_id);
        response.u()[cell_id] = 91.0;
        response.v()[cell_id] = -93.0;
    }
    cfd::compute_momentum_pressure_response(mesh, u_system, v_system, response);
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        require_near(response.u()[cell_id], mesh.cell_areas()[cell_id] / second_u_diagonals.at(cell_id), 0.0,
                     "Repeated computation did not use the current u diagonal.");
        require_near(response.v()[cell_id], mesh.cell_areas()[cell_id] / second_v_diagonals.at(cell_id), 0.0,
                     "Repeated computation did not use the current v diagonal.");
    }
}

void test_valid_computation_does_not_modify_inputs()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_quadrilateral_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    cfd::ScalarLinearSystem u_system{mesh};
    cfd::ScalarLinearSystem v_system{mesh};
    seed_system(u_system, 2.0, 11.0);
    seed_system(v_system, 4.0, -17.0);
    const double cell_area_before{mesh.cell_areas()[0]};
    const cfd::Index node_count_before{mesh.node_count()};
    const cfd::Index face_count_before{mesh.face_count()};
    cfd::CellMomentumPressureResponse response{mesh.cell_count()};

    cfd::compute_momentum_pressure_response(mesh, u_system, v_system, response);

    require_near(mesh.cell_areas()[0], cell_area_before, 0.0, "Response computation modified the Mesh area.");
    require(mesh.node_count() == node_count_before && mesh.face_count() == face_count_before,
            "Response computation modified Mesh cardinalities.");
    require_near(u_system.diagonal()[0], 2.0, 0.0, "Response computation modified the u diagonal.");
    require_near(v_system.diagonal()[0], 4.0, 0.0, "Response computation modified the v diagonal.");
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        require_near(u_system.owner_neighbor_coefficients()[face_id], 11.0, 0.0,
                     "Response computation modified a u owner-neighbor coefficient.");
        require_near(u_system.neighbor_owner_coefficients()[face_id], 12.0, 0.0,
                     "Response computation modified a u neighbor-owner coefficient.");
        require_near(v_system.owner_neighbor_coefficients()[face_id], -17.0, 0.0,
                     "Response computation modified a v owner-neighbor coefficient.");
        require_near(v_system.neighbor_owner_coefficients()[face_id], -16.0, 0.0,
                     "Response computation modified a v neighbor-owner coefficient.");
    }
    require_near(u_system.rhs()[0], 13.0, 0.0, "Response computation modified the u-system RHS.");
    require_near(v_system.rhs()[0], -15.0, 0.0, "Response computation modified the v-system RHS.");
}

void test_rejects_api_mismatches_before_mutating_output()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_quadrilateral_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    cfd::ScalarLinearSystem u_system{mesh};
    cfd::ScalarLinearSystem v_system{mesh};
    u_system.diagonal()[0] = 2.0;
    v_system.diagonal()[0] = 4.0;

    cfd::CellMomentumPressureResponse wrong_size_response{mesh.cell_count() + 1};
    for (cfd::Index cell_id = 0; cell_id < wrong_size_response.size(); ++cell_id)
    {
        wrong_size_response.u()[cell_id] = 17.0;
        wrong_size_response.v()[cell_id] = -19.0;
    }
    require_throws<std::invalid_argument>(
        [&]() { cfd::compute_momentum_pressure_response(mesh, u_system, v_system, wrong_size_response); },
        "Momentum pressure response accepted an output with incorrect cardinality.");
    require_response_unchanged(wrong_size_response, 17.0, -19.0,
                               "An output-cardinality error partially modified the response.");

    cfd::CellMomentumPressureResponse response{mesh.cell_count()};
    response.u()[0] = 23.0;
    response.v()[0] = -29.0;
    cfd::MeshBuildResult other_build_result{cfd::build_mesh(make_single_quadrilateral_raw_mesh())};
    cfd::ScalarLinearSystem other_system{other_build_result.mesh};
    other_system.diagonal()[0] = 3.0;
    require_throws<std::invalid_argument>(
        [&]() { cfd::compute_momentum_pressure_response(mesh, other_system, v_system, response); },
        "Momentum pressure response accepted a u system referencing another Mesh.");
    require_response_unchanged(response, 23.0, -29.0, "A u-system Mesh mismatch partially modified the response.");
    require_throws<std::invalid_argument>(
        [&]() { cfd::compute_momentum_pressure_response(mesh, u_system, other_system, response); },
        "Momentum pressure response accepted a v system referencing another Mesh.");
    require_response_unchanged(response, 23.0, -29.0, "A v-system Mesh mismatch partially modified the response.");
}

void test_rejects_invalid_diagonals_before_mutating_output()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_quadrilateral_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    cfd::ScalarLinearSystem u_system{mesh};
    cfd::ScalarLinearSystem v_system{mesh};
    cfd::CellMomentumPressureResponse response{mesh.cell_count()};
    response.u()[0] = 31.0;
    response.v()[0] = -37.0;

    const auto require_rejected = [&](const bool invalidate_u, const double invalid_diagonal) {
        u_system.diagonal()[0] = invalidate_u ? invalid_diagonal : 2.0;
        v_system.diagonal()[0] = invalidate_u ? 4.0 : invalid_diagonal;
        require_throws<std::runtime_error>(
            [&]() { cfd::compute_momentum_pressure_response(mesh, u_system, v_system, response); },
            "Momentum pressure response accepted an invalid momentum diagonal.");
        require_response_unchanged(response, 31.0, -37.0,
                                   "An invalid momentum diagonal partially modified the response.");
    };

    require_rejected(true, 0.0);
    require_rejected(false, 0.0);
    require_rejected(true, -1.0);
    require_rejected(false, -1.0);
    require_rejected(true, std::numeric_limits<double>::quiet_NaN());
    require_rejected(false, std::numeric_limits<double>::infinity());
    require_rejected(true, std::numeric_limits<double>::denorm_min());
}

} // namespace

int main()
{
    int failure_count{};

    failure_count +=
        cfd::test::run_test("momentum pressure-response field access", test_field_construction_and_component_access);
    failure_count += cfd::test::run_test("momentum pressure-response deep copy", test_copy_construction_is_deep);
    failure_count += cfd::test::run_test("one-cell momentum pressure response",
                                         test_computes_one_cell_response_from_independent_diagonals);
    failure_count += cfd::test::run_test("multi-cell repeated momentum pressure response",
                                         test_multiple_cells_use_current_diagonals_and_overwrite_output);
    failure_count += cfd::test::run_test("momentum pressure-response input immutability",
                                         test_valid_computation_does_not_modify_inputs);
    failure_count += cfd::test::run_test("momentum pressure-response API validation",
                                         test_rejects_api_mismatches_before_mutating_output);
    failure_count += cfd::test::run_test("momentum pressure-response numerical validation",
                                         test_rejects_invalid_diagonals_before_mutating_output);

    return cfd::test::finish_tests(failure_count, "cell momentum pressure response");
}
