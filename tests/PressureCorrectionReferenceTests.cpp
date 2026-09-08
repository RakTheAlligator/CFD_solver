#include "cfd/numerics/PressureCorrectionReference.hpp"

#include "cfd/field/FaceFluxField.hpp"
#include "cfd/field/FacePressureResponseField.hpp"
#include "cfd/field/PressureCorrectionBoundaryConditions.hpp"
#include "cfd/linear_algebra/EigenConjugateGradientSolver.hpp"
#include "cfd/linear_algebra/ScalarLinearSystem.hpp"
#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/mesh/Types.hpp"
#include "cfd/meshing/RawMeshData.hpp"
#include "cfd/numerics/IncompressiblePressureCorrectionAssembler.hpp"

#include "support/TestUtils.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{

using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws;
using cfd::test::test_tolerance;

constexpr cfd::BoundaryId bottom_boundary_id{0};
constexpr cfd::BoundaryId right_boundary_id{1};
constexpr cfd::BoundaryId top_boundary_id{2};
constexpr cfd::BoundaryId left_boundary_id{3};

[[nodiscard]]
cfd::RawMeshData make_two_cell_rectangle_raw_mesh(const bool separate_boundary_groups)
{
    cfd::RawMeshData raw_mesh;
    raw_mesh.nodes = {
        {0.0, 0.0}, {1.0, 0.0}, {2.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}, {2.0, 1.0},
    };
    raw_mesh.cell_types = {cfd::CellType::Quadrilateral, cfd::CellType::Quadrilateral};
    raw_mesh.cell_nodes = {0, 1, 4, 3, 1, 2, 5, 4};
    raw_mesh.cell_node_offsets = {0, 4, 8};

    if (separate_boundary_groups)
    {
        raw_mesh.boundary_groups = {
            {bottom_boundary_id, "bottom"},
            {right_boundary_id, "right"},
            {top_boundary_id, "top"},
            {left_boundary_id, "left"},
        };
        raw_mesh.boundary_edges = {
            {{0, 1}, bottom_boundary_id}, {{1, 2}, bottom_boundary_id}, {{2, 5}, right_boundary_id},
            {{5, 4}, top_boundary_id},    {{4, 3}, top_boundary_id},    {{3, 0}, left_boundary_id},
        };
    }
    else
    {
        raw_mesh.boundary_groups = {{0, "wall"}};
        raw_mesh.boundary_edges = {
            {{0, 1}, 0}, {{1, 2}, 0}, {{2, 5}, 0}, {{5, 4}, 0}, {{4, 3}, 0}, {{3, 0}, 0},
        };
    }
    return raw_mesh;
}

[[nodiscard]]
cfd::RawMeshData make_three_cell_strip_raw_mesh()
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
    raw_mesh.cell_nodes = {0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6};
    raw_mesh.cell_node_offsets = {0, 4, 8, 12};
    raw_mesh.boundary_groups = {{wall_boundary_id, "wall"}};
    raw_mesh.boundary_edges = {
        {{0, 1}, wall_boundary_id}, {{1, 2}, wall_boundary_id}, {{2, 3}, wall_boundary_id}, {{3, 7}, wall_boundary_id},
        {{7, 6}, wall_boundary_id}, {{6, 5}, wall_boundary_id}, {{5, 4}, wall_boundary_id}, {{4, 0}, wall_boundary_id},
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
    throw std::runtime_error("Pressure-correction solve fixture has no internal face.");
}

[[nodiscard]]
cfd::Index boundary_face_id(const cfd::Mesh &mesh, const cfd::BoundaryId boundary_id)
{
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (mesh.face_boundary_ids()[face_id] == boundary_id)
        {
            return face_id;
        }
    }
    throw std::runtime_error("Pressure-correction solve fixture has no requested boundary face.");
}

[[nodiscard]]
double normalized_residual(const cfd::ScalarLinearSystem &system, const std::array<double, 2> &solution)
{
    std::array<double, 2> matrix_product{};
    system.apply_matrix(solution, matrix_product);
    const double residual_norm{std::hypot(system.rhs()[0] - matrix_product[0], system.rhs()[1] - matrix_product[1])};
    return residual_norm / std::hypot(system.rhs()[0], system.rhs()[1]);
}

void test_zero_reference_preserves_symmetry_and_unrelated_entries()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_three_cell_strip_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    cfd::ScalarLinearSystem system{mesh};
    std::array<cfd::Index, 2> internal_faces{cfd::invalid_index, cfd::invalid_index};
    cfd::Index internal_count{};
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};
        if (adjacency.is_boundary())
        {
            system.owner_neighbor_coefficients()[face_id] = 101.0 + static_cast<double>(face_id);
            system.neighbor_owner_coefficients()[face_id] = -201.0 - static_cast<double>(face_id);
            continue;
        }
        switch (internal_count)
        {
        case 0:
            internal_faces[0] = face_id;
            break;
        case 1:
            internal_faces[1] = face_id;
            break;
        default:
            cfd::test::fail("Reference fixture has too many internal faces.");
        }
        ++internal_count;
    }
    require(internal_count == internal_faces.size(), "Reference fixture does not have two internal faces.");
    system.diagonal()[0] = 2.0;
    system.diagonal()[1] = 5.0;
    system.diagonal()[2] = 3.0;
    for (const cfd::Index face_id : internal_faces)
    {
        const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};
        const double coefficient{adjacency.owner == 0 || adjacency.neighbor == 0 ? -2.0 : -3.0};
        system.owner_neighbor_coefficients()[face_id] = coefficient;
        system.neighbor_owner_coefficients()[face_id] = coefficient;
    }
    system.rhs()[0] = 7.0;
    system.rhs()[1] = 8.0;
    system.rhs()[2] = 9.0;
    const std::array constant_mode{1.0, 1.0, 1.0};
    std::array matrix_product{99.0, 99.0, 99.0};
    system.apply_matrix(constant_mode, matrix_product);
    for (const double value : matrix_product)
    {
        require_near(value, 0.0, 0.0, "Internal Laplacian matrix does not have the constant nullspace.");
    }

    cfd::apply_zero_pressure_correction_reference(0, system);

    require(system.diagonal()[0] == 2.0 && system.diagonal()[1] == 5.0 && system.diagonal()[2] == 3.0,
            "Zero reference changed a matrix diagonal.");
    require(system.rhs()[0] == 0.0 && system.rhs()[1] == 8.0 && system.rhs()[2] == 9.0,
            "Zero reference changed an unrelated RHS or failed to clear the reference RHS.");
    cfd::Index incident_face_count{};
    cfd::Index unrelated_internal_face_count{};
    for (const cfd::Index face_id : internal_faces)
    {
        const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};
        if (adjacency.owner == 0 || adjacency.neighbor == 0)
        {
            require(system.owner_neighbor_coefficients()[face_id] == 0.0 &&
                        system.neighbor_owner_coefficients()[face_id] == 0.0,
                    "Zero reference did not remove both directed incident-face couplings.");
            ++incident_face_count;
        }
        else
        {
            require(system.owner_neighbor_coefficients()[face_id] == -3.0 &&
                        system.neighbor_owner_coefficients()[face_id] == -3.0,
                    "Zero reference changed an unrelated internal-face coupling.");
            ++unrelated_internal_face_count;
        }
    }
    require(incident_face_count == 1 && unrelated_internal_face_count == 1,
            "Reference fixture does not exercise both incident and unrelated internal faces.");
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (!mesh.face_adjacencies()[face_id].is_boundary())
        {
            continue;
        }
        require(system.owner_neighbor_coefficients()[face_id] == 101.0 + static_cast<double>(face_id) &&
                    system.neighbor_owner_coefficients()[face_id] == -201.0 - static_cast<double>(face_id),
                "Zero reference changed an unrelated boundary coefficient.");
    }
    system.apply_matrix(constant_mode, matrix_product);
    require(std::hypot(matrix_product[0], std::hypot(matrix_product[1], matrix_product[2])) > 0.0,
            "Zero reference did not remove the constant nullspace.");
}

void test_zero_reference_rejects_invalid_inputs_before_mutation()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh(false))};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::Index face_id{internal_face_id(mesh)};

    {
        cfd::ScalarLinearSystem system{mesh};
        system.diagonal()[0] = 2.0;
        system.diagonal()[1] = 2.0;
        system.owner_neighbor_coefficients()[face_id] = -2.0;
        system.neighbor_owner_coefficients()[face_id] = -2.0;
        system.rhs()[0] = 3.0;
        system.rhs()[1] = -3.0;
        require_throws<std::invalid_argument>(
            [&]() { cfd::apply_zero_pressure_correction_reference(mesh.cell_count(), system); },
            "Zero pressure-correction reference accepted an out-of-range cell.");
        require(system.diagonal()[0] == 2.0 && system.diagonal()[1] == 2.0 &&
                    system.owner_neighbor_coefficients()[face_id] == -2.0 &&
                    system.neighbor_owner_coefficients()[face_id] == -2.0 && system.rhs()[0] == 3.0 &&
                    system.rhs()[1] == -3.0,
                "Invalid reference index partially modified the system.");
    }

    for (const double invalid_diagonal :
         {0.0, -1.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
          -std::numeric_limits<double>::infinity()})
    {
        cfd::ScalarLinearSystem system{mesh};
        system.diagonal()[0] = invalid_diagonal;
        system.diagonal()[1] = 2.0;
        system.owner_neighbor_coefficients()[face_id] = -2.0;
        system.neighbor_owner_coefficients()[face_id] = -2.0;
        system.rhs()[0] = 3.0;
        system.rhs()[1] = -3.0;
        require_throws<std::runtime_error>([&]() { cfd::apply_zero_pressure_correction_reference(0, system); },
                                           "Zero reference accepted an invalid reference diagonal.");
        const bool invalid_diagonal_unchanged{std::isnan(invalid_diagonal) ? std::isnan(system.diagonal()[0])
                                                                           : system.diagonal()[0] == invalid_diagonal};
        require(invalid_diagonal_unchanged && system.diagonal()[1] == 2.0 &&
                    system.owner_neighbor_coefficients()[face_id] == -2.0 &&
                    system.neighbor_owner_coefficients()[face_id] == -2.0 && system.rhs()[0] == 3.0 &&
                    system.rhs()[1] == -3.0,
                "Invalid reference diagonal partially modified the system.");
    }
}

void test_all_fixed_mass_flux_system_solves_after_zero_reference()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh(false))};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressiblePressureCorrectionAssembler assembler{mesh};
    cfd::FaceFluxField provisional_mass_flux{mesh.face_count()};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count()};
    const cfd::Index internal_id{internal_face_id(mesh)};
    face_pressure_response[internal_id] = 2.0;
    provisional_mass_flux[boundary_face_id(mesh, 0)] = 1.0;
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};
        if (adjacency.is_boundary() && adjacency.owner == 1)
        {
            provisional_mass_flux[face_id] = -1.0;
            break;
        }
    }
    const cfd::PressureCorrectionBoundaryConditions boundary_conditions{
        mesh.boundary_groups().size(),
        {cfd::PressureCorrectionBoundaryConditionType::FixedMassFlux},
    };
    cfd::ScalarLinearSystem system{mesh};
    assembler.add_internal_face_contributions(provisional_mass_flux, face_pressure_response, system);
    assembler.add_boundary_provisional_flux_rhs(provisional_mass_flux, system);
    assembler.add_boundary_pressure_response(boundary_conditions, face_pressure_response, system);
    const std::array constant_mode{1.0, 1.0};
    std::array matrix_product{99.0, 99.0};
    system.apply_matrix(constant_mode, matrix_product);
    require_near(matrix_product[0], 0.0, 0.0, "All-FixedMassFlux matrix lacks the constant nullspace.");
    require_near(matrix_product[1], 0.0, 0.0, "All-FixedMassFlux matrix lacks the constant nullspace.");

    cfd::apply_zero_pressure_correction_reference(0, system);
    require(system.owner_neighbor_coefficients()[internal_id] == 0.0 &&
                system.neighbor_owner_coefficients()[internal_id] == 0.0,
            "Gauge fixing did not preserve pressure-correction symmetry.");
    cfd::EigenConjugateGradientSolver solver{{1.0e-14, 20}};
    solver.compute_matrix(system);
    std::array solution{0.0, 0.0};
    const cfd::LinearSolveResult result{solver.solve(system.rhs(), solution)};

    require(result.converged, "CG did not converge for the gauge-fixed all-FixedMassFlux system.");
    require_near(solution[0], 0.0, 1.0e-13, "CG did not preserve the exact zero reference.");
    require_near(solution[1], 0.5, 1.0e-13, "CG returned an incorrect non-reference pressure correction.");
    require(normalized_residual(system, solution) < 1.0e-13,
            "Gauge-fixed pressure-correction solution has an excessive algebraic residual.");
}

void test_fixed_pressure_response_anchors_system_without_reference()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh(true))};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressiblePressureCorrectionAssembler assembler{mesh};
    cfd::FaceFluxField provisional_mass_flux{mesh.face_count()};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count()};
    const cfd::Index internal_id{internal_face_id(mesh)};
    const cfd::Index left_face_id{boundary_face_id(mesh, left_boundary_id)};
    const cfd::Index right_face_id{boundary_face_id(mesh, right_boundary_id)};
    face_pressure_response[internal_id] = 2.0;
    face_pressure_response[right_face_id] = 3.0;
    provisional_mass_flux[left_face_id] = -1.5;
    provisional_mass_flux[right_face_id] = 3.0;
    std::vector<cfd::PressureCorrectionBoundaryConditionType> conditions(
        mesh.boundary_groups().size(), cfd::PressureCorrectionBoundaryConditionType::FixedMassFlux);
    conditions[right_boundary_id] = cfd::PressureCorrectionBoundaryConditionType::FixedPressure;
    const cfd::PressureCorrectionBoundaryConditions boundary_conditions{mesh.boundary_groups().size(),
                                                                        std::move(conditions)};
    cfd::ScalarLinearSystem system{mesh};
    assembler.add_internal_face_contributions(provisional_mass_flux, face_pressure_response, system);
    assembler.add_boundary_provisional_flux_rhs(provisional_mass_flux, system);
    assembler.add_boundary_pressure_response(boundary_conditions, face_pressure_response, system);
    const std::array constant_mode{1.0, 1.0};
    std::array matrix_product{99.0, 99.0};
    system.apply_matrix(constant_mode, matrix_product);
    require(std::hypot(matrix_product[0], matrix_product[1]) > 0.0,
            "FixedPressure boundary response did not remove the constant nullspace.");
    require(system.owner_neighbor_coefficients()[internal_id] == system.neighbor_owner_coefficients()[internal_id],
            "FixedPressure pressure-correction matrix is not symmetric.");

    cfd::EigenConjugateGradientSolver solver{{1.0e-14, 20}};
    solver.compute_matrix(system);
    std::array solution{0.0, 0.0};
    const cfd::LinearSolveResult result{solver.solve(system.rhs(), solution)};

    require(result.converged, "CG did not converge for the FixedPressure-anchored system.");
    require_near(solution[0], 0.25, 1.0e-13, "CG returned an incorrect owner pressure correction.");
    require_near(solution[1], -0.5, 1.0e-13, "CG returned an incorrect neighbor pressure correction.");
    require(normalized_residual(system, solution) < 1.0e-13,
            "FixedPressure-anchored solution has an excessive algebraic residual.");
}

} // namespace

int main()
{
    int failure_count{};

    failure_count += cfd::test::run_test("zero pressure-correction reference structure",
                                         test_zero_reference_preserves_symmetry_and_unrelated_entries);
    failure_count += cfd::test::run_test("zero pressure-correction reference validation",
                                         test_zero_reference_rejects_invalid_inputs_before_mutation);
    failure_count += cfd::test::run_test("all-FixedMassFlux pressure-correction CG solve",
                                         test_all_fixed_mass_flux_system_solves_after_zero_reference);
    failure_count += cfd::test::run_test("FixedPressure-anchored pressure-correction CG solve",
                                         test_fixed_pressure_response_anchors_system_without_reference);

    return cfd::test::finish_tests(failure_count, "pressure-correction reference and solve");
}
