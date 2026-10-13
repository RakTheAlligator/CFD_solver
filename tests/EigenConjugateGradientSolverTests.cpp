#include "cfd/linear_algebra/EigenConjugateGradientSolver.hpp"

#include "cfd/linear_algebra/ScalarLinearSystem.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/meshing/RawMeshData.hpp"

#include "support/MeshFixtures.hpp"
#include "support/TestUtils.hpp"

#include <Eigen/SparseCore>

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace
{

using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws;

using DefaultEigenSparseMatrix = Eigen::SparseMatrix<double>;
using EigenSparseStorageIndex = DefaultEigenSparseMatrix::StorageIndex;

static_assert(std::is_signed_v<EigenSparseStorageIndex>);
static_assert(std::numeric_limits<EigenSparseStorageIndex>::digits <= std::numeric_limits<Eigen::Index>::digits);

void test_distinguishes_iteration_and_sparse_storage_index_ranges()
{
    if constexpr (std::numeric_limits<EigenSparseStorageIndex>::digits < std::numeric_limits<Eigen::Index>::digits &&
                  std::numeric_limits<EigenSparseStorageIndex>::digits < std::numeric_limits<cfd::Index>::digits)
    {
        constexpr cfd::Index first_unsupported_sparse_index{
            static_cast<cfd::Index>(std::numeric_limits<EigenSparseStorageIndex>::max()) + 1};
        static_assert(!std::in_range<EigenSparseStorageIndex>(first_unsupported_sparse_index));
        static_assert(std::in_range<Eigen::Index>(first_unsupported_sparse_index));

        const cfd::EigenConjugateGradientSolver solver{{1.0e-10, first_unsupported_sparse_index}};
        static_cast<void>(solver);
    }
}

[[nodiscard]]
cfd::RawMeshData make_two_cell_mesh()
{
    constexpr cfd::BoundaryId wall_boundary_id{0};

    cfd::RawMeshData raw_mesh;
    raw_mesh.nodes = {
        {0.0, 0.0}, {1.0, 0.0}, {2.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}, {2.0, 1.0},
    };
    raw_mesh.cell_types = {cfd::CellType::Quadrilateral, cfd::CellType::Quadrilateral};
    raw_mesh.cell_nodes = {0, 1, 4, 3, 1, 2, 5, 4};
    raw_mesh.cell_node_offsets = {0, 4, 8};
    raw_mesh.boundary_groups = {{wall_boundary_id, "wall"}};
    raw_mesh.boundary_edges = {
        {{0, 1}, wall_boundary_id}, {{1, 2}, wall_boundary_id}, {{2, 5}, wall_boundary_id},
        {{5, 4}, wall_boundary_id}, {{4, 3}, wall_boundary_id}, {{3, 0}, wall_boundary_id},
    };
    return raw_mesh;
}

[[nodiscard]]
cfd::Index internal_face_id(const cfd::Mesh &mesh)
{
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (!mesh.face_adjacencies()[face_id].is_boundary())
        {
            return face_id;
        }
    }
    throw std::runtime_error("Two-cell Eigen solver fixture has no internal face.");
}

void set_spd_matrix(cfd::ScalarLinearSystem &system)
{
    const cfd::Index face_id{internal_face_id(system.mesh())};
    system.diagonal()[0] = 4.0;
    system.diagonal()[1] = 3.0;
    system.owner_neighbor_coefficients()[face_id] = -1.0;
    system.neighbor_owner_coefficients()[face_id] = -1.0;
}

double normalized_residual(const cfd::ScalarLinearSystem &system, const std::array<double, 2> &rhs,
                           const std::array<double, 2> &solution)
{
    std::array<double, 2> matrix_product{};
    system.apply_matrix(solution, matrix_product);
    const double residual_norm{std::hypot(rhs[0] - matrix_product[0], rhs[1] - matrix_product[1])};
    return residual_norm / std::hypot(rhs[0], rhs[1]);
}

void test_solves_known_spd_system_and_reuses_preparation()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_mesh())};
    cfd::ScalarLinearSystem system{build_result.mesh};
    set_spd_matrix(system);
    cfd::EigenConjugateGradientSolver solver{{1.0e-14, 20}};
    solver.compute_matrix(system);

    const std::array first_rhs{2.0, 5.0};
    std::array first_solution{0.0, 0.0};
    const cfd::LinearSolveResult first_result{solver.solve(first_rhs, first_solution)};

    require(first_result.converged, "Eigen conjugate gradient did not converge for the first SPD system.");
    require_near(first_solution[0], 1.0, 1.0e-13, "Eigen conjugate gradient returned an incorrect first value.");
    require_near(first_solution[1], 2.0, 1.0e-13, "Eigen conjugate gradient returned an incorrect second value.");
    require(normalized_residual(system, first_rhs, first_solution) < 1.0e-13,
            "The first Eigen solution has an excessive explicit algebraic residual.");

    const std::array second_rhs{-7.0, 10.0};
    std::array second_solution{first_solution};
    const cfd::LinearSolveResult second_result{solver.solve(second_rhs, second_solution)};

    require(second_result.converged, "Eigen conjugate gradient did not converge after matrix reuse.");
    require_near(second_solution[0], -1.0, 1.0e-13, "Repeated Eigen solve returned an incorrect first value.");
    require_near(second_solution[1], 3.0, 1.0e-13, "Repeated Eigen solve returned an incorrect second value.");
    require(normalized_residual(system, second_rhs, second_solution) < 1.0e-13,
            "The repeated Eigen solution has an excessive explicit algebraic residual.");
}

void test_reports_iteration_limit_without_throwing()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_mesh())};
    cfd::ScalarLinearSystem system{build_result.mesh};
    set_spd_matrix(system);
    cfd::EigenConjugateGradientSolver solver{{1.0e-15, 1}};
    solver.compute_matrix(system);
    const std::array rhs{1.0, 0.0};
    std::array solution{0.0, 0.0};

    const cfd::LinearSolveResult result{solver.solve(rhs, solution)};

    require(!result.converged, "Eigen conjugate gradient silently accepted iteration-limit non-convergence.");
    require(result.iteration_count == 1, "Eigen conjugate gradient reported an unexpected limited iteration count.");
    require(std::isfinite(result.estimated_relative_error),
            "Eigen conjugate gradient reported a non-finite non-convergence error estimate.");
}

void test_rejects_invalid_options_and_usage()
{
    const auto require_invalid_options = [](const cfd::ConjugateGradientOptions options) {
        require_throws<std::invalid_argument>([options]() { const cfd::EigenConjugateGradientSolver solver{options}; },
                                              "Eigen conjugate gradient accepted invalid options.");
    };
    require_invalid_options({0.0, 10});
    require_invalid_options({1.0, 10});
    require_invalid_options({std::numeric_limits<double>::quiet_NaN(), 10});
    require_invalid_options({1.0e-10, 0});

    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_mesh())};
    cfd::ScalarLinearSystem system{build_result.mesh};
    set_spd_matrix(system);
    cfd::EigenConjugateGradientSolver solver;
    const std::array valid_rhs{2.0, 5.0};
    std::array valid_solution{0.0, 0.0};

    require_throws<std::logic_error>(
        [&solver, &valid_rhs, &valid_solution]() { static_cast<void>(solver.solve(valid_rhs, valid_solution)); },
        "Eigen conjugate gradient accepted solve before matrix preparation.");

    system.diagonal()[0] = std::numeric_limits<double>::infinity();
    require_throws<std::invalid_argument>([&solver, &system]() { solver.compute_matrix(system); },
                                          "Eigen conjugate gradient accepted a non-finite matrix coefficient.");
    set_spd_matrix(system);
    solver.compute_matrix(system);

    require_throws<std::invalid_argument>(
        [&solver, &valid_solution]() {
            const std::array wrong_rhs{1.0, 2.0, 3.0};
            static_cast<void>(solver.solve(wrong_rhs, valid_solution));
        },
        "Eigen conjugate gradient accepted an RHS with incorrect cardinality.");
    require_throws<std::invalid_argument>(
        [&solver, &valid_rhs]() {
            std::array wrong_solution{0.0, 0.0, 0.0};
            static_cast<void>(solver.solve(valid_rhs, wrong_solution));
        },
        "Eigen conjugate gradient accepted a solution with incorrect cardinality.");
    require_throws<std::invalid_argument>(
        [&solver]() {
            std::array aliased{1.0, 2.0};
            static_cast<void>(solver.solve(aliased, aliased));
        },
        "Eigen conjugate gradient accepted overlapping RHS and solution spans.");
}

void test_internal_tolerance_margin_preserves_public_contract()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(cfd::test::make_single_quadrilateral_raw_mesh())};
    cfd::ScalarLinearSystem system{build_result.mesh};
    system.diagonal()[0] = 1.0;
    constexpr double requested_tolerance{0.125};
    cfd::EigenConjugateGradientSolver solver{{requested_tolerance, 20}};
    solver.compute_matrix(system);
    const std::array rhs{1.0};
    std::array solution{0.8828125};
    // The initial residual 0.1171875 already satisfies the public threshold,
    // but not Eigen's stricter internal target.
    const auto result{solver.solve(rhs, solution)};
    const double independent_residual{std::abs(rhs[0] - solution[0]) / std::abs(rhs[0])};
    require(result.converged && independent_residual <= requested_tolerance,
            "Internal tolerance margin failed the public explicit-residual contract.");
    require(solution[0] == 1.0, "Eigen did not target a stricter internal tolerance.");
    require_near(result.estimated_relative_error, independent_residual, 0.0,
                 "Internal margin changed the reported explicit residual.");

    for (const double tiny_tolerance :
         {std::numeric_limits<double>::denorm_min(), 2.0 * std::numeric_limits<double>::denorm_min(),
          std::numeric_limits<double>::min()})
    {
        cfd::EigenConjugateGradientSolver tiny_solver{{tiny_tolerance, 20}};
        tiny_solver.compute_matrix(system);
        solution[0] = 1.0;
        const auto exact_result{tiny_solver.solve(rhs, solution)};
        require(exact_result.converged && exact_result.estimated_relative_error == 0.0 &&
                    exact_result.iteration_count == 0 && solution[0] == 1.0,
                "Tiny valid tolerance did not preserve an exact initial guess.");
    }
}

void test_explicit_residual_at_extreme_rhs_scales()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(cfd::test::make_single_quadrilateral_raw_mesh())};
    cfd::ScalarLinearSystem system{build_result.mesh};
    system.diagonal()[0] = 1.0;
    constexpr double tolerance{1.0e-10};
    cfd::EigenConjugateGradientSolver solver{{tolerance, 20}};
    solver.compute_matrix(system);

    for (const double rhs_value : {1.0, 1.0e-200, 1.0e200, std::numeric_limits<double>::denorm_min()})
    {
        const std::array rhs{rhs_value};
        std::array solution{0.0};
        cfd::LinearSolveResult result;
        try
        {
            result = solver.solve(rhs, solution);
        }
        catch (const std::runtime_error &)
        {
            require(rhs_value != 1.0, "Unit-scale solve failed numerically.");
            continue; // Eigen may fail at extreme scales, but must never report false success.
        }

        require(std::isfinite(solution[0]), "Extreme-scale solve returned a non-finite solution.");
        // A scalar relative residual needs no squaring, even at extreme scales.
        const double explicit_error{
            static_cast<double>(std::abs((static_cast<long double>(rhs_value) - solution[0]) / rhs_value))};
        require(std::isfinite(result.estimated_relative_error), "Extreme-scale residual diagnostic is not finite.");
        require_near(result.estimated_relative_error, explicit_error, 1.0e-14,
                     "Returned error differs from the robust explicit residual.");
        require(result.converged == (explicit_error <= tolerance),
                "Extreme-scale solve reported convergence inconsistent with its explicit residual.");
    }
}

void test_explicit_residual_uses_prepared_matrix_and_initial_guess()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(cfd::test::make_single_quadrilateral_raw_mesh())};
    cfd::ScalarLinearSystem system{build_result.mesh};
    system.diagonal()[0] = 1.0;
    cfd::EigenConjugateGradientSolver solver{{1.0e-10, 20}};
    solver.compute_matrix(system);
    system.diagonal()[0] = 2.0; // Residual validation must not read the unprepared FV coefficients.

    for (const double rhs_value : {1.0, 1.0e-100, 1.0e100})
    {
        const std::array rhs{rhs_value};
        std::array solution{rhs_value};
        const auto result{solver.solve(rhs, solution)};
        require(result.converged && result.iteration_count == 0 && result.estimated_relative_error == 0.0,
                "Exact nonzero initial guess was not recognized against the prepared matrix.");
        require(solution[0] == rhs_value, "Solve changed an exact caller-provided initial guess.");
    }

    const std::array rhs{-0.0};
    std::array solution{7.0};
    const auto result{solver.solve(rhs, solution)};
    require(result.converged && result.iteration_count == 0 && result.estimated_relative_error == 0.0 &&
                solution[0] == 0.0,
            "Exactly zero RHS did not return an exact zero solution and diagnostic.");
}

void test_explicit_residual_respects_tolerance()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(cfd::test::make_single_quadrilateral_raw_mesh())};
    cfd::ScalarLinearSystem system{build_result.mesh};
    system.diagonal()[0] = 1.0;
    constexpr double tolerance{0.125};
    cfd::EigenConjugateGradientSolver solver{{tolerance, 20}};
    solver.compute_matrix(system);

    for (const double rhs_value : {1.0, 1.0e-160, 1.0e150})
    {
        for (const double initial_error :
             {0.0625, std::nextafter(tolerance, 0.0), tolerance, std::nextafter(tolerance, 1.0), 0.25})
        {
            const std::array rhs{rhs_value};
            std::array solution{(1.0 - initial_error) * rhs_value};
            cfd::LinearSolveResult result;
            try
            {
                result = solver.solve(rhs, solution);
            }
            catch (const std::runtime_error &)
            {
                require(rhs_value == 1.0e-160, "Tolerance-threshold solve unexpectedly failed numerically.");
                continue;
            }
            const long double explicit_error{std::abs((static_cast<long double>(rhs_value) - solution[0]) / rhs_value)};
            require_near(result.estimated_relative_error, static_cast<double>(explicit_error), 1.0e-14,
                         "Tolerance-threshold error is not the explicit residual.");
            require(!result.converged || (std::isfinite(explicit_error) && explicit_error <= tolerance),
                    "Solve accepted a residual above the configured tolerance.");
        }
    }
}

} // namespace

int main()
{
    int failure_count{};

    failure_count += cfd::test::run_test("Eigen ConjugateGradient internal margin and public tolerance",
                                         test_internal_tolerance_margin_preserves_public_contract);

    failure_count += cfd::test::run_test("Eigen ConjugateGradient extreme-scale explicit residual",
                                         test_explicit_residual_at_extreme_rhs_scales);
    failure_count += cfd::test::run_test("Eigen ConjugateGradient prepared matrix and exact initial guess",
                                         test_explicit_residual_uses_prepared_matrix_and_initial_guess);
    failure_count += cfd::test::run_test("Eigen ConjugateGradient explicit residual tolerance",
                                         test_explicit_residual_respects_tolerance);

    failure_count += cfd::test::run_test("Eigen conjugate gradient iteration and sparse storage index distinction",
                                         test_distinguishes_iteration_and_sparse_storage_index_ranges);
    failure_count += cfd::test::run_test("Eigen conjugate gradient known SPD systems",
                                         test_solves_known_spd_system_and_reuses_preparation);
    failure_count += cfd::test::run_test("Eigen conjugate gradient non-convergence reporting",
                                         test_reports_iteration_limit_without_throwing);
    failure_count += cfd::test::run_test("Eigen conjugate gradient validation", test_rejects_invalid_options_and_usage);

    return cfd::test::finish_tests(failure_count, "Eigen conjugate gradient solver");
}
