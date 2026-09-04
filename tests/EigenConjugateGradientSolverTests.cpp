#include "cfd/linear_algebra/EigenConjugateGradientSolver.hpp"

#include "cfd/linear_algebra/ScalarLinearSystem.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/meshing/RawMeshData.hpp"

#include "support/TestUtils.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace
{

using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws;

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

} // namespace

int main()
{
    int failure_count{};

    failure_count += cfd::test::run_test("Eigen conjugate gradient known SPD systems",
                                         test_solves_known_spd_system_and_reuses_preparation);
    failure_count += cfd::test::run_test("Eigen conjugate gradient non-convergence reporting",
                                         test_reports_iteration_limit_without_throwing);
    failure_count += cfd::test::run_test("Eigen conjugate gradient validation", test_rejects_invalid_options_and_usage);

    return cfd::test::finish_tests(failure_count, "Eigen conjugate gradient solver");
}
