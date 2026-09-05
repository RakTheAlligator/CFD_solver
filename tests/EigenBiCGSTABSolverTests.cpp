#include "cfd/linear_algebra/EigenBiCGSTABSolver.hpp"

#include "cfd/linear_algebra/ScalarLinearSystem.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/meshing/RawMeshData.hpp"

#include "support/TestUtils.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <type_traits>

namespace
{

using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws;

static_assert(!std::is_copy_constructible_v<cfd::EigenBiCGSTABSolver>);
static_assert(!std::is_copy_assignable_v<cfd::EigenBiCGSTABSolver>);
static_assert(!std::is_move_constructible_v<cfd::EigenBiCGSTABSolver>);
static_assert(!std::is_move_assignable_v<cfd::EigenBiCGSTABSolver>);

[[nodiscard]]
cfd::RawMeshData make_three_cell_mesh()
{
    constexpr cfd::BoundaryId wall_boundary_id{0};

    cfd::RawMeshData raw_mesh;
    raw_mesh.nodes = {
        {0.0, 0.0}, {1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}, {2.0, 1.0}, {3.0, 1.0},
    };
    raw_mesh.cell_types = {
        cfd::CellType::Quadrilateral,
        cfd::CellType::Quadrilateral,
        cfd::CellType::Quadrilateral,
    };
    raw_mesh.cell_nodes = {
        0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6,
    };
    raw_mesh.cell_node_offsets = {0, 4, 8, 12};
    raw_mesh.boundary_groups = {{wall_boundary_id, "wall"}};
    raw_mesh.boundary_edges = {
        {{0, 1}, wall_boundary_id}, {{1, 2}, wall_boundary_id}, {{2, 3}, wall_boundary_id}, {{3, 7}, wall_boundary_id},
        {{7, 6}, wall_boundary_id}, {{6, 5}, wall_boundary_id}, {{5, 4}, wall_boundary_id}, {{4, 0}, wall_boundary_id},
    };
    return raw_mesh;
}

void set_off_diagonal(cfd::ScalarLinearSystem &system, const cfd::Index row, const cfd::Index column,
                      const double value)
{
    const auto face_adjacencies{system.mesh().face_adjacencies()};
    for (cfd::Index face_id = 0; face_id < system.face_count(); ++face_id)
    {
        const cfd::FaceAdjacency &adjacency{face_adjacencies[face_id]};
        if (adjacency.is_boundary())
        {
            continue;
        }
        if (adjacency.owner == row && adjacency.neighbor == column)
        {
            system.owner_neighbor_coefficients()[face_id] = value;
            return;
        }
        if (adjacency.neighbor == row && adjacency.owner == column)
        {
            system.neighbor_owner_coefficients()[face_id] = value;
            return;
        }
    }
    throw std::runtime_error("Three-cell Eigen BiCGSTAB fixture is missing an expected adjacency.");
}

[[nodiscard]]
double off_diagonal(const cfd::ScalarLinearSystem &system, const cfd::Index row, const cfd::Index column)
{
    const auto face_adjacencies{system.mesh().face_adjacencies()};
    for (cfd::Index face_id = 0; face_id < system.face_count(); ++face_id)
    {
        const cfd::FaceAdjacency &adjacency{face_adjacencies[face_id]};
        if (adjacency.is_boundary())
        {
            continue;
        }
        if (adjacency.owner == row && adjacency.neighbor == column)
        {
            return system.owner_neighbor_coefficients()[face_id];
        }
        if (adjacency.neighbor == row && adjacency.owner == column)
        {
            return system.neighbor_owner_coefficients()[face_id];
        }
    }
    throw std::runtime_error("Three-cell Eigen BiCGSTAB fixture is missing an expected adjacency.");
}

void set_nonsymmetric_matrix(cfd::ScalarLinearSystem &system)
{
    system.clear_matrix();
    system.diagonal()[0] = 4.0;
    system.diagonal()[1] = 5.0;
    system.diagonal()[2] = 6.0;
    set_off_diagonal(system, 0, 1, -1.0);
    set_off_diagonal(system, 1, 0, -2.0);
    set_off_diagonal(system, 1, 2, -1.0);
    set_off_diagonal(system, 2, 1, -3.0);
}

template <std::size_t Size>
double normalized_residual(const cfd::ScalarLinearSystem &system, const std::array<double, Size> &rhs,
                           const std::array<double, Size> &solution)
{
    std::array<double, Size> matrix_product{};
    system.apply_matrix(solution, matrix_product);

    double residual_squared{};
    double rhs_squared{};
    for (std::size_t index = 0; index < Size; ++index)
    {
        const double residual{matrix_product.at(index) - rhs.at(index)};
        residual_squared += residual * residual;
        rhs_squared += rhs.at(index) * rhs.at(index);
    }
    return std::sqrt(residual_squared / rhs_squared);
}

void test_solves_known_nonsymmetric_system_and_reuses_matrix()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_three_cell_mesh())};
    cfd::ScalarLinearSystem system{build_result.mesh};
    set_nonsymmetric_matrix(system);
    require(off_diagonal(system, 0, 1) != off_diagonal(system, 1, 0),
            "Eigen BiCGSTAB test matrix is unexpectedly symmetric.");
    require(off_diagonal(system, 1, 2) != off_diagonal(system, 2, 1),
            "Eigen BiCGSTAB test matrix does not contain directed coefficients.");

    cfd::EigenBiCGSTABSolver solver{{1.0e-14, 100}};
    solver.compute_matrix(system);

    const std::array first_rhs{2.0, 9.0, -12.0};
    std::array first_solution{0.0, 0.0, 0.0};
    const cfd::LinearSolveResult first_result{solver.solve(first_rhs, first_solution)};

    require(first_result.converged, "Eigen BiCGSTAB did not converge for the known nonsymmetric system.");
    require_near(first_solution[0], 1.0, 1.0e-12, "Eigen BiCGSTAB returned an incorrect first value.");
    require_near(first_solution[1], 2.0, 1.0e-12, "Eigen BiCGSTAB returned an incorrect second value.");
    require_near(first_solution[2], -1.0, 1.0e-12, "Eigen BiCGSTAB returned an incorrect third value.");
    require(normalized_residual(system, first_rhs, first_solution) < 1.0e-13,
            "Eigen BiCGSTAB solution has an excessive explicit normalized residual.");

    const std::array second_rhs{-8.5, 3.5, 16.5};
    std::array second_solution{first_solution};
    const cfd::LinearSolveResult second_result{solver.solve(second_rhs, second_solution)};

    require(second_result.converged, "Eigen BiCGSTAB did not converge after matrix reuse.");
    require_near(second_solution[0], -2.0, 1.0e-12, "Repeated BiCGSTAB solve returned an incorrect first value.");
    require_near(second_solution[1], 0.5, 1.0e-12, "Repeated BiCGSTAB solve returned an incorrect second value.");
    require_near(second_solution[2], 3.0, 1.0e-12, "Repeated BiCGSTAB solve returned an incorrect third value.");
    require(normalized_residual(system, second_rhs, second_solution) < 1.0e-13,
            "Repeated BiCGSTAB solution has an excessive explicit normalized residual.");
}

void test_uses_caller_initial_guess()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_three_cell_mesh())};
    cfd::ScalarLinearSystem system{build_result.mesh};
    set_nonsymmetric_matrix(system);
    cfd::EigenBiCGSTABSolver solver{{1.0e-14, 100}};
    solver.compute_matrix(system);
    const std::array rhs{3.0, -8.0, 15.0};
    std::array solution{0.5, -1.0, 2.0};

    const cfd::LinearSolveResult result{solver.solve(rhs, solution)};

    require(result.converged, "Eigen BiCGSTAB rejected an exact caller-provided initial guess.");
    require(result.iteration_count == 0, "Eigen BiCGSTAB did not use the caller-provided initial guess.");
    require_near(solution[0], 0.5, 0.0, "BiCGSTAB changed the exact initial guess's first value.");
    require_near(solution[1], -1.0, 0.0, "BiCGSTAB changed the exact initial guess's second value.");
    require_near(solution[2], 2.0, 0.0, "BiCGSTAB changed the exact initial guess's third value.");
}

void test_ignores_boundary_off_diagonal_storage()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_three_cell_mesh())};
    cfd::ScalarLinearSystem system{build_result.mesh};
    set_nonsymmetric_matrix(system);
    for (cfd::Index face_id = 0; face_id < system.face_count(); ++face_id)
    {
        if (system.mesh().face_adjacencies()[face_id].is_boundary())
        {
            system.owner_neighbor_coefficients()[face_id] = 1.0e6;
            system.neighbor_owner_coefficients()[face_id] = -2.0e6;
        }
    }

    cfd::EigenBiCGSTABSolver solver{{1.0e-14, 100}};
    solver.compute_matrix(system);
    const std::array rhs{-6.0, 11.5, -3.0};
    std::array solution{0.0, 0.0, 0.0};

    const cfd::LinearSolveResult result{solver.solve(rhs, solution)};

    require(result.converged, "Boundary off-diagonal storage prevented BiCGSTAB convergence.");
    require_near(solution[0], -1.0, 1.0e-12, "Boundary storage altered the first solution value.");
    require_near(solution[1], 2.0, 1.0e-12, "Boundary storage altered the second solution value.");
    require_near(solution[2], 0.5, 1.0e-12, "Boundary storage altered the third solution value.");
}

void test_reports_iteration_limit_without_throwing()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_three_cell_mesh())};
    cfd::ScalarLinearSystem system{build_result.mesh};
    set_nonsymmetric_matrix(system);
    cfd::EigenBiCGSTABSolver solver{{1.0e-15, 1}};
    solver.compute_matrix(system);
    const std::array rhs{1.0, 2.0, 3.0};
    std::array solution{0.0, 0.0, 0.0};

    const cfd::LinearSolveResult result{solver.solve(rhs, solution)};

    require(!result.converged, "Eigen BiCGSTAB silently accepted iteration-limit non-convergence.");
    require(result.iteration_count == 1, "Eigen BiCGSTAB reported an unexpected limited iteration count.");
    require(std::isfinite(result.estimated_relative_error),
            "Eigen BiCGSTAB reported a non-finite non-convergence error estimate.");
}

void test_rejects_invalid_options_and_usage()
{
    const auto require_invalid_options = [](const cfd::BiCGSTABOptions options) {
        require_throws<std::invalid_argument>([options]() { const cfd::EigenBiCGSTABSolver solver{options}; },
                                              "Eigen BiCGSTAB accepted invalid options.");
    };
    require_invalid_options({0.0, 10});
    require_invalid_options({-1.0e-10, 10});
    require_invalid_options({1.0, 10});
    require_invalid_options({1.5, 10});
    require_invalid_options({std::numeric_limits<double>::quiet_NaN(), 10});
    require_invalid_options({std::numeric_limits<double>::infinity(), 10});
    require_invalid_options({1.0e-10, 0});
    require_invalid_options(
        {1.0e-10, static_cast<cfd::Index>(std::numeric_limits<Eigen::Index>::max()) + cfd::Index{1}});

    cfd::MeshBuildResult build_result{cfd::build_mesh(make_three_cell_mesh())};
    cfd::ScalarLinearSystem system{build_result.mesh};
    set_nonsymmetric_matrix(system);
    cfd::EigenBiCGSTABSolver solver;
    const std::array valid_rhs{2.0, 9.0, -12.0};
    std::array valid_solution{0.0, 0.0, 0.0};

    require_throws<std::logic_error>(
        [&solver, &valid_rhs, &valid_solution]() { static_cast<void>(solver.solve(valid_rhs, valid_solution)); },
        "Eigen BiCGSTAB accepted solve before matrix preparation.");

    solver.compute_matrix(system);
    system.diagonal()[0] = std::numeric_limits<double>::infinity();
    require_throws<std::invalid_argument>([&solver, &system]() { solver.compute_matrix(system); },
                                          "Eigen BiCGSTAB accepted a non-finite diagonal coefficient.");
    require_throws<std::logic_error>(
        [&solver, &valid_rhs, &valid_solution]() { static_cast<void>(solver.solve(valid_rhs, valid_solution)); },
        "Eigen BiCGSTAB reused stale prepared state after matrix preparation failed.");
    set_nonsymmetric_matrix(system);
    set_off_diagonal(system, 0, 1, std::numeric_limits<double>::quiet_NaN());
    require_throws<std::invalid_argument>([&solver, &system]() { solver.compute_matrix(system); },
                                          "Eigen BiCGSTAB accepted a non-finite directed coefficient.");
    set_nonsymmetric_matrix(system);
    solver.compute_matrix(system);

    require_throws<std::invalid_argument>(
        [&solver, &valid_solution]() {
            const std::array wrong_rhs{1.0, 2.0, 3.0, 4.0};
            static_cast<void>(solver.solve(wrong_rhs, valid_solution));
        },
        "Eigen BiCGSTAB accepted an RHS with incorrect cardinality.");
    require_throws<std::invalid_argument>(
        [&solver, &valid_rhs]() {
            std::array wrong_solution{0.0, 0.0, 0.0, 0.0};
            static_cast<void>(solver.solve(valid_rhs, wrong_solution));
        },
        "Eigen BiCGSTAB accepted a solution with incorrect cardinality.");
    require_throws<std::invalid_argument>(
        [&solver]() {
            std::array storage{1.0, 2.0, 3.0, 4.0};
            const std::span<const double> rhs{storage.data(), 3};
            const std::span<double> solution{storage.data() + 1, 3};
            static_cast<void>(solver.solve(rhs, solution));
        },
        "Eigen BiCGSTAB accepted overlapping RHS and solution spans.");
    require_throws<std::invalid_argument>(
        [&solver, &valid_solution]() {
            const std::array rhs{2.0, std::numeric_limits<double>::quiet_NaN(), -12.0};
            static_cast<void>(solver.solve(rhs, valid_solution));
        },
        "Eigen BiCGSTAB accepted a non-finite RHS.");
    require_throws<std::invalid_argument>(
        [&solver, &valid_rhs]() {
            std::array solution{0.0, std::numeric_limits<double>::infinity(), 0.0};
            static_cast<void>(solver.solve(valid_rhs, solution));
        },
        "Eigen BiCGSTAB accepted a non-finite initial guess.");
}

} // namespace

int main()
{
    int failure_count{};

    failure_count += cfd::test::run_test("Eigen BiCGSTAB known nonsymmetric systems and matrix reuse",
                                         test_solves_known_nonsymmetric_system_and_reuses_matrix);
    failure_count += cfd::test::run_test("Eigen BiCGSTAB caller initial guess", test_uses_caller_initial_guess);
    failure_count += cfd::test::run_test("Eigen BiCGSTAB boundary storage", test_ignores_boundary_off_diagonal_storage);
    failure_count +=
        cfd::test::run_test("Eigen BiCGSTAB non-convergence reporting", test_reports_iteration_limit_without_throwing);
    failure_count += cfd::test::run_test("Eigen BiCGSTAB validation", test_rejects_invalid_options_and_usage);

    return cfd::test::finish_tests(failure_count, "Eigen BiCGSTAB solver");
}
