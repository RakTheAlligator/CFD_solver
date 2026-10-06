#include "cfd/linear_algebra/EigenBiCGSTABSolver.hpp"
#include "cfd/linear_algebra/EigenConjugateGradientSolver.hpp"
#include "cfd/linear_algebra/EigenSparseMatrixPattern.hpp"
#include "cfd/linear_algebra/ScalarLinearSystem.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/meshing/RawMeshData.hpp"
#include "cfd/numerics/PressureCorrectionReference.hpp"

#include "support/TestUtils.hpp"

#include <algorithm>
#include <array>
#include <memory>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace
{

using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws;
using SparseMatrix = cfd::EigenSparseMatrixPattern::SparseMatrix;
using StorageIndex = SparseMatrix::StorageIndex;

cfd::RawMeshData make_two_cell_mesh(const bool triangles)
{
    cfd::RawMeshData raw;
    raw.boundary_groups = {{0, "wall"}};
    if (triangles)
    {
        raw.nodes = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
        raw.cell_types = {cfd::CellType::Triangle, cfd::CellType::Triangle};
        raw.cell_nodes = {0, 1, 2, 0, 2, 3};
        raw.cell_node_offsets = {0, 3, 6};
        raw.boundary_edges = {{{0, 1}, 0}, {{1, 2}, 0}, {{2, 3}, 0}, {{3, 0}, 0}};
    }
    else
    {
        raw.nodes = {{0, 0}, {1, 0}, {2, 0}, {0, 1}, {1, 1}, {2, 1}};
        raw.cell_types = {cfd::CellType::Quadrilateral, cfd::CellType::Quadrilateral};
        raw.cell_nodes = {0, 1, 4, 3, 1, 2, 5, 4};
        raw.cell_node_offsets = {0, 4, 8};
        raw.boundary_edges = {{{0, 1}, 0}, {{1, 2}, 0}, {{2, 5}, 0}, {{5, 4}, 0}, {{4, 3}, 0}, {{3, 0}, 0}};
    }
    return raw;
}

SparseMatrix reference_matrix(const cfd::ScalarLinearSystem &system)
{
    using Triplet = Eigen::Triplet<double, StorageIndex>;
    std::vector<Triplet> entries;
    for (cfd::Index cell = 0; cell < system.cell_count(); ++cell)
    {
        const auto index{static_cast<StorageIndex>(cell)};
        entries.emplace_back(index, index, system.diagonal()[cell]);
    }
    for (cfd::Index face = 0; face < system.face_count(); ++face)
    {
        const auto &adjacency{system.mesh().face_adjacencies()[face]};
        if (!adjacency.is_boundary())
        {
            entries.emplace_back(static_cast<StorageIndex>(adjacency.owner),
                                 static_cast<StorageIndex>(adjacency.neighbor),
                                 system.owner_neighbor_coefficients()[face]);
            entries.emplace_back(static_cast<StorageIndex>(adjacency.neighbor),
                                 static_cast<StorageIndex>(adjacency.owner),
                                 system.neighbor_owner_coefficients()[face]);
        }
    }
    SparseMatrix matrix{static_cast<Eigen::Index>(system.cell_count()), static_cast<Eigen::Index>(system.cell_count())};
    matrix.setFromTriplets(entries.begin(), entries.end());
    return matrix;
}

void require_same_pattern(const SparseMatrix &matrix, const SparseMatrix &reference)
{
    require(matrix.isCompressed() && matrix.rows() == reference.rows() && matrix.cols() == reference.cols() &&
                matrix.nonZeros() == reference.nonZeros(),
            "Sparse dimensions or compression differ from the triplet reference.");
    const std::span<const StorageIndex> offsets{matrix.outerIndexPtr(), static_cast<cfd::Index>(matrix.cols()) + 1};
    const std::span<const StorageIndex> reference_offsets{reference.outerIndexPtr(), offsets.size()};
    const std::span<const StorageIndex> rows{matrix.innerIndexPtr(), static_cast<cfd::Index>(matrix.nonZeros())};
    const std::span<const StorageIndex> reference_rows{reference.innerIndexPtr(), rows.size()};
    require(std::equal(offsets.begin(), offsets.end(), reference_offsets.begin()) &&
                std::equal(rows.begin(), rows.end(), reference_rows.begin()),
            "Compressed sparse pattern changed.");
}

void require_same_matrix(const SparseMatrix &matrix, const SparseMatrix &reference)
{
    require_same_pattern(matrix, reference);
    const std::span<const double> values{matrix.valuePtr(), static_cast<cfd::Index>(matrix.nonZeros())};
    const std::span<const double> reference_values{reference.valuePtr(), values.size()};
    require(std::equal(values.begin(), values.end(), reference_values.begin()),
            "Numerical matrix differs from the triplet reference.");
}

void set_matrix(cfd::ScalarLinearSystem &system, const double scale)
{
    system.diagonal()[0] = 4.0 * scale;
    system.diagonal()[1] = 3.0 * scale;
    for (cfd::Index face = 0; face < system.face_count(); ++face)
    {
        const double coefficient{system.mesh().face_adjacencies()[face].is_boundary() ? 123.0 : -scale};
        system.owner_neighbor_coefficients()[face] = coefficient;
        system.neighbor_owner_coefficients()[face] = coefficient;
    }
}

void test_triangles_quads_zeros_and_pressure_gauge()
{
    for (const bool triangles : {false, true})
    {
        auto build{cfd::build_mesh(make_two_cell_mesh(triangles))};
        cfd::ScalarLinearSystem system{build.mesh};
        const cfd::EigenSparseMatrixPattern pattern{build.mesh};
        const SparseMatrix zero_matrix{pattern.create_matrix(system)};
        require(zero_matrix.nonZeros() == 4, "Explicit zero coefficients were pruned.");
        require_same_matrix(zero_matrix, reference_matrix(system));

        set_matrix(system, 1.0);
        const SparseMatrix first{pattern.create_matrix(system)};
        require_same_matrix(first, reference_matrix(system));
        set_matrix(system, 2.0);
        const SparseMatrix second{pattern.create_matrix(system)};
        require_same_matrix(second, reference_matrix(system));
        require_same_pattern(first, second);
        require(first.coeff(0, 0) == 4.0, "Independent matrix values were overwritten.");

        cfd::apply_zero_pressure_correction_reference(0, system);
        const SparseMatrix gauged{pattern.create_matrix(system)};
        require_same_matrix(gauged, reference_matrix(system));
        require_same_pattern(second, gauged);
        require(gauged.coeff(0, 1) == 0.0 && gauged.coeff(1, 0) == 0.0, "Pressure gauge did not zero both couplings.");
    }
}

void test_shared_pattern_updates_and_independent_solvers()
{
    auto build{cfd::build_mesh(make_two_cell_mesh(false))};
    cfd::ScalarLinearSystem system{build.mesh};
    const auto pattern{std::make_shared<cfd::EigenSparseMatrixPattern>(build.mesh)};
    cfd::EigenBiCGSTABSolver u{pattern, {1.0e-14, 20}};
    cfd::EigenBiCGSTABSolver v{pattern, {1.0e-14, 20}};
    cfd::EigenConjugateGradientSolver pressure{pattern, {1.0e-14, 20}};
    set_matrix(system, 1.0);
    u.compute_matrix(system);
    set_matrix(system, 2.0);
    v.compute_matrix(system);
    set_matrix(system, 3.0);
    pressure.compute_matrix(system);
    require(pattern.use_count() == 4, "The three solvers did not retain one shared pattern.");

    const std::array rhs{2.0, 5.0};
    std::array<double, 2> solution{};
    require(u.solve(rhs, solution).converged, "First shared-pattern solve failed.");
    require_near(solution[0], 1.0, 1.0e-13, "u matrix lost its independent values.");
    require_near(solution[1], 2.0, 1.0e-13, "u matrix lost its independent values.");
    require(v.solve(rhs, solution).converged, "Second shared-pattern solve failed.");
    require_near(solution[0], 0.5, 1.0e-13, "v matrix lost its independent values.");
    require_near(solution[1], 1.0, 1.0e-13, "v matrix lost its independent values.");
    require(pressure.solve(rhs, solution).converged, "Pressure shared-pattern solve failed.");
    require_near(solution[0], 1.0 / 3.0, 1.0e-13, "Pressure matrix lost its independent values.");

    set_matrix(system, 4.0);
    u.compute_matrix(system);
    v.compute_matrix(system);
    pressure.compute_matrix(system);
    require(pattern.use_count() == 4, "Coefficient updates replaced the shared pattern.");
    require(u.solve(rhs, solution).converged && v.solve(rhs, solution).converged &&
                pressure.solve(rhs, solution).converged,
            "Repeated shared-pattern preparation/solve failed.");
    require_near(solution[0], 0.25, 1.0e-13, "Updated coefficient was not used.");
    require_near(solution[1], 0.5, 1.0e-13, "Updated coefficient was not used.");

    cfd::apply_zero_pressure_correction_reference(0, system);
    u.compute_matrix(system);
    v.compute_matrix(system);
    pressure.compute_matrix(system);
    const std::array gauge_rhs{0.0, 12.0};
    require(u.solve(gauge_rhs, solution).converged, "u explicit-zero update failed.");
    require_near(solution[0], 0.0, 1.0e-13, "u retained its old coupling.");
    require_near(solution[1], 1.0, 1.0e-13, "u retained its old coupling.");
    require(v.solve(gauge_rhs, solution).converged, "v explicit-zero update failed.");
    require_near(solution[0], 0.0, 1.0e-13, "v retained its old coupling.");
    require_near(solution[1], 1.0, 1.0e-13, "v retained its old coupling.");
    require(pressure.solve(gauge_rhs, solution).converged, "Gauge update solve failed.");
    require_near(solution[0], 0.0, 1.0e-13, "Gauge constraint was not applied.");
    require_near(solution[1], 1.0, 1.0e-13, "Zeroed coupling was not updated.");
}

void test_standalone_solvers_update_and_rebuild_for_another_mesh()
{
    auto quads{cfd::build_mesh(make_two_cell_mesh(false))};
    auto triangles{cfd::build_mesh(make_two_cell_mesh(true))};
    cfd::ScalarLinearSystem quad_system{quads.mesh};
    cfd::ScalarLinearSystem triangle_system{triangles.mesh};
    cfd::EigenBiCGSTABSolver momentum;
    cfd::EigenConjugateGradientSolver pressure;
    const std::array rhs{2.0, 5.0};
    std::array<double, 2> solution{};
    const auto check_solutions = [&](cfd::ScalarLinearSystem &system, const double scale) {
        set_matrix(system, scale);
        momentum.compute_matrix(system);
        pressure.compute_matrix(system);
        require(momentum.solve(rhs, solution).converged, "Standalone BiCGSTAB preparation failed.");
        require_near(solution[0], 1.0 / scale, 1.0e-12, "BiCGSTAB used stale matrix coefficients.");
        require_near(solution[1], 2.0 / scale, 1.0e-12, "BiCGSTAB used stale matrix coefficients.");
        require(pressure.solve(rhs, solution).converged, "Standalone CG preparation failed.");
        require_near(solution[0], 1.0 / scale, 1.0e-12, "CG used stale matrix coefficients.");
        require_near(solution[1], 2.0 / scale, 1.0e-12, "CG used stale matrix coefficients.");
    };
    check_solutions(quad_system, 1.0);
    check_solutions(quad_system, 2.0);
    check_solutions(triangle_system, 3.0);
}

void test_rejects_incompatible_or_replaced_mesh()
{
    auto first{cfd::build_mesh(make_two_cell_mesh(false))};
    auto other{cfd::build_mesh(make_two_cell_mesh(false))};
    cfd::ScalarLinearSystem first_system{first.mesh};
    cfd::ScalarLinearSystem other_system{other.mesh};
    const auto pattern{std::make_shared<cfd::EigenSparseMatrixPattern>(first.mesh)};
    cfd::EigenBiCGSTABSolver u{pattern, {}};
    cfd::EigenConjugateGradientSolver pressure{pattern, {}};
    set_matrix(first_system, 1.0);
    u.compute_matrix(first_system);
    pressure.compute_matrix(first_system);
    require_throws<std::invalid_argument>([&] { u.compute_matrix(other_system); }, "Another Mesh was accepted.");
    require_throws<std::invalid_argument>([&] { pressure.compute_matrix(other_system); }, "Another Mesh was accepted.");
    const std::array rhs{2.0, 5.0};
    std::array<double, 2> solution{};
    require_throws<std::logic_error>([&] { static_cast<void>(u.solve(rhs, solution)); },
                                     "Stale u preparation survived.");
    require_throws<std::logic_error>([&] { static_cast<void>(pressure.solve(rhs, solution)); },
                                     "Stale pressure preparation survived.");

    first.mesh = std::move(other.mesh);
    require_throws<std::invalid_argument>([&] { static_cast<void>(pattern->create_matrix(first_system)); },
                                          "Replaced topology storage was accepted.");
    require_throws<std::invalid_argument>([&] { u.compute_matrix(first_system); }, "Replaced Mesh was accepted.");
    require_throws<std::invalid_argument>([&] { pressure.compute_matrix(first_system); },
                                          "Replaced Mesh was accepted.");
    require_throws<std::invalid_argument>([] { const cfd::EigenBiCGSTABSolver solver{nullptr, {}}; },
                                          "Null shared BiCGSTAB pattern was accepted.");
    require_throws<std::invalid_argument>([] { const cfd::EigenConjugateGradientSolver solver{nullptr, {}}; },
                                          "Null shared CG pattern was accepted.");
}

} // namespace

int main()
{
    int failures{};
    failures += cfd::test::run_test("Sparse pattern triangles/quads, zeros and gauge",
                                    test_triangles_quads_zeros_and_pressure_gauge);
    failures += cfd::test::run_test("Shared pattern updates and independent numerical matrices",
                                    test_shared_pattern_updates_and_independent_solvers);
    failures += cfd::test::run_test("Standalone solvers update/rebuild patterns",
                                    test_standalone_solvers_update_and_rebuild_for_another_mesh);
    failures += cfd::test::run_test("Sparse pattern identity and storage validation",
                                    test_rejects_incompatible_or_replaced_mesh);
    return cfd::test::finish_tests(failures, "Eigen sparse matrix pattern");
}
