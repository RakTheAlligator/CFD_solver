#include "cfd/numerics/IncompressiblePressureCorrectionAssembler.hpp"

#include "cfd/field/FaceFluxField.hpp"
#include "cfd/field/FacePressureResponseField.hpp"
#include "cfd/field/PressureCorrectionBoundaryConditions.hpp"
#include "cfd/linear_algebra/ScalarLinearSystem.hpp"
#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/mesh/Types.hpp"
#include "cfd/meshing/RawMeshData.hpp"

#include "support/MeshFixtures.hpp"
#include "support/TestUtils.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace
{

using cfd::test::make_two_triangle_raw_mesh;
using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws;

static_assert(std::is_nothrow_constructible_v<cfd::IncompressiblePressureCorrectionAssembler, const cfd::Mesh &>);
static_assert(!std::is_copy_constructible_v<cfd::IncompressiblePressureCorrectionAssembler>);
static_assert(!std::is_copy_assignable_v<cfd::IncompressiblePressureCorrectionAssembler>);
static_assert(std::is_nothrow_move_constructible_v<cfd::IncompressiblePressureCorrectionAssembler>);
static_assert(!std::is_move_assignable_v<cfd::IncompressiblePressureCorrectionAssembler>);

static_assert(!std::is_default_constructible_v<cfd::PressureCorrectionBoundaryConditions>);
static_assert(std::is_copy_constructible_v<cfd::PressureCorrectionBoundaryConditions>);
static_assert(!std::is_copy_assignable_v<cfd::PressureCorrectionBoundaryConditions>);
static_assert(std::is_nothrow_move_constructible_v<cfd::PressureCorrectionBoundaryConditions>);
static_assert(!std::is_move_assignable_v<cfd::PressureCorrectionBoundaryConditions>);

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

[[nodiscard]]
cfd::RawMeshData make_four_boundary_group_two_triangle_raw_mesh()
{
    constexpr cfd::BoundaryId bottom_boundary_id{0};
    constexpr cfd::BoundaryId right_boundary_id{1};
    constexpr cfd::BoundaryId top_boundary_id{2};
    constexpr cfd::BoundaryId left_boundary_id{3};

    cfd::RawMeshData raw_mesh{make_two_triangle_raw_mesh()};
    raw_mesh.boundary_groups = {
        {bottom_boundary_id, "bottom"},
        {right_boundary_id, "right"},
        {top_boundary_id, "top"},
        {left_boundary_id, "left"},
    };
    raw_mesh.boundary_edges = {
        {{0, 1}, bottom_boundary_id},
        {{1, 2}, right_boundary_id},
        {{2, 3}, top_boundary_id},
        {{3, 0}, left_boundary_id},
    };
    return raw_mesh;
}

[[nodiscard]]
cfd::PressureCorrectionBoundaryConditions make_uniform_pressure_correction_boundary_conditions(
    const cfd::Mesh &mesh, const cfd::PressureCorrectionBoundaryConditionType condition)
{
    return {mesh.boundary_groups().size(),
            std::vector<cfd::PressureCorrectionBoundaryConditionType>(mesh.boundary_groups().size(), condition)};
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
    throw std::runtime_error("Pressure-correction test fixture has no internal face.");
}

[[nodiscard]]
std::array<cfd::Index, 2> two_internal_face_ids(const cfd::Mesh &mesh)
{
    std::array<cfd::Index, 2> face_ids{cfd::invalid_index, cfd::invalid_index};
    cfd::Index internal_count{};
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (mesh.face_adjacencies()[face_id].is_boundary())
        {
            continue;
        }
        switch (internal_count)
        {
        case 0:
            face_ids[0] = face_id;
            break;
        case 1:
            face_ids[1] = face_id;
            break;
        default:
            cfd::test::fail("Three-cell fixture has too many internal faces.");
        }
        ++internal_count;
    }
    require(internal_count == face_ids.size(), "Three-cell fixture does not have two internal faces.");
    return face_ids;
}

[[nodiscard]]
std::array<cfd::Index, 2> two_boundary_face_ids(const cfd::Mesh &mesh)
{
    std::array<cfd::Index, 2> face_ids{cfd::invalid_index, cfd::invalid_index};
    cfd::Index boundary_count{};
    for (cfd::Index face_id = 0; face_id < mesh.face_count() && boundary_count < face_ids.size(); ++face_id)
    {
        if (!mesh.face_adjacencies()[face_id].is_boundary())
        {
            continue;
        }
        switch (boundary_count)
        {
        case 0:
            face_ids[0] = face_id;
            break;
        case 1:
            face_ids[1] = face_id;
            break;
        default:
            cfd::test::fail("Boundary-face collection exceeded its requested size.");
        }
        ++boundary_count;
    }
    require(boundary_count == face_ids.size(), "Pressure-correction test fixture has fewer than two boundary faces.");
    return face_ids;
}

void seed_system(cfd::ScalarLinearSystem &system)
{
    for (cfd::Index cell_id = 0; cell_id < system.cell_count(); ++cell_id)
    {
        system.diagonal()[cell_id] = 11.0 + static_cast<double>(cell_id);
        system.rhs()[cell_id] = 21.0 + static_cast<double>(cell_id);
    }
    for (cfd::Index face_id = 0; face_id < system.face_count(); ++face_id)
    {
        system.owner_neighbor_coefficients()[face_id] = 31.0 + static_cast<double>(face_id);
        system.neighbor_owner_coefficients()[face_id] = 41.0 + static_cast<double>(face_id);
    }
}

void require_seeded_system_unchanged(const cfd::ScalarLinearSystem &system, const std::string &context)
{
    for (cfd::Index cell_id = 0; cell_id < system.cell_count(); ++cell_id)
    {
        require(system.diagonal()[cell_id] == 11.0 + static_cast<double>(cell_id), context + " diagonal changed.");
        require(system.rhs()[cell_id] == 21.0 + static_cast<double>(cell_id), context + " RHS changed.");
    }
    for (cfd::Index face_id = 0; face_id < system.face_count(); ++face_id)
    {
        require(system.owner_neighbor_coefficients()[face_id] == 31.0 + static_cast<double>(face_id),
                context + " owner-neighbor coefficient changed.");
        require(system.neighbor_owner_coefficients()[face_id] == 41.0 + static_cast<double>(face_id),
                context + " neighbor-owner coefficient changed.");
    }
}

void require_seeded_matrix_unchanged(const cfd::ScalarLinearSystem &system, const std::string &context)
{
    for (cfd::Index cell_id = 0; cell_id < system.cell_count(); ++cell_id)
    {
        require(system.diagonal()[cell_id] == 11.0 + static_cast<double>(cell_id), context + " diagonal changed.");
    }
    for (cfd::Index face_id = 0; face_id < system.face_count(); ++face_id)
    {
        require(system.owner_neighbor_coefficients()[face_id] == 31.0 + static_cast<double>(face_id),
                context + " owner-neighbor coefficient changed.");
        require(system.neighbor_owner_coefficients()[face_id] == 41.0 + static_cast<double>(face_id),
                context + " neighbor-owner coefficient changed.");
    }
}

void require_seeded_rhs_and_off_diagonal_coefficients_unchanged(const cfd::ScalarLinearSystem &system,
                                                                const std::string &context)
{
    for (cfd::Index cell_id = 0; cell_id < system.cell_count(); ++cell_id)
    {
        require(system.rhs()[cell_id] == 21.0 + static_cast<double>(cell_id), context + " RHS changed.");
    }
    for (cfd::Index face_id = 0; face_id < system.face_count(); ++face_id)
    {
        require(system.owner_neighbor_coefficients()[face_id] == 31.0 + static_cast<double>(face_id),
                context + " owner-neighbor coefficient changed.");
        require(system.neighbor_owner_coefficients()[face_id] == 41.0 + static_cast<double>(face_id),
                context + " neighbor-owner coefficient changed.");
    }
}

template <typename Exception, typename Function>
void require_rejected_without_mutation(Function &&function, const cfd::ScalarLinearSystem &system,
                                       const std::string &message)
{
    require_throws<Exception>(std::forward<Function>(function), message);
    require_seeded_system_unchanged(system, message + " The system was partially modified.");
}

void test_exact_internal_face_contributions_and_local_conservation()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_triangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressiblePressureCorrectionAssembler assembler{mesh};
    cfd::FaceFluxField provisional_mass_flux{mesh.face_count()};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count()};
    const cfd::Index face_id{internal_face_id(mesh)};
    const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};
    provisional_mass_flux[face_id] = 3.0;
    face_pressure_response[face_id] = 2.0;
    cfd::ScalarLinearSystem system{mesh};

    assembler.add_internal_face_contributions(provisional_mass_flux, face_pressure_response, system);

    require_near(system.diagonal()[adjacency.owner], 2.0, 0.0, "Owner diagonal contribution is incorrect.");
    require_near(system.owner_neighbor_coefficients()[face_id], -2.0, 0.0, "Owner-neighbor contribution is incorrect.");
    require_near(system.diagonal()[adjacency.neighbor], 2.0, 0.0, "Neighbor diagonal contribution is incorrect.");
    require_near(system.neighbor_owner_coefficients()[face_id], -2.0, 0.0, "Neighbor-owner contribution is incorrect.");
    require_near(system.rhs()[adjacency.owner], -3.0, 0.0, "Owner provisional-flux RHS is incorrect.");
    require_near(system.rhs()[adjacency.neighbor], 3.0, 0.0, "Neighbor provisional-flux RHS is incorrect.");

    require_near(system.diagonal()[adjacency.owner] + system.owner_neighbor_coefficients()[face_id], 0.0, 0.0,
                 "Owner matrix row is not locally conservative.");
    require_near(system.diagonal()[adjacency.neighbor] + system.neighbor_owner_coefficients()[face_id], 0.0, 0.0,
                 "Neighbor matrix row is not locally conservative.");
    require_near(system.rhs()[adjacency.owner] + system.rhs()[adjacency.neighbor], 0.0, 0.0,
                 "Internal-face RHS is not locally conservative.");
    require(system.owner_neighbor_coefficients()[face_id] == system.neighbor_owner_coefficients()[face_id],
            "Internal-face pressure-correction matrix is not symmetric.");
}

void test_negative_provisional_flux_reverses_rhs_contributions()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_triangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressiblePressureCorrectionAssembler assembler{mesh};
    cfd::FaceFluxField provisional_mass_flux{mesh.face_count()};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count(), 1.0};
    const cfd::Index face_id{internal_face_id(mesh)};
    const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};
    provisional_mass_flux[face_id] = -4.5;
    cfd::ScalarLinearSystem system{mesh};

    assembler.add_internal_face_contributions(provisional_mass_flux, face_pressure_response, system);

    require_near(system.rhs()[adjacency.owner], 4.5, 0.0, "Negative provisional flux produced an incorrect owner RHS.");
    require_near(system.rhs()[adjacency.neighbor], -4.5, 0.0,
                 "Negative provisional flux produced an incorrect neighbor RHS.");
}

void test_assembly_is_additive()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_triangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressiblePressureCorrectionAssembler assembler{mesh};
    cfd::FaceFluxField provisional_mass_flux{mesh.face_count()};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count(), 1.0};
    const cfd::Index face_id{internal_face_id(mesh)};
    const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};
    constexpr double provisional_flux{-1.25};
    constexpr double pressure_response{2.5};
    provisional_mass_flux[face_id] = provisional_flux;
    face_pressure_response[face_id] = pressure_response;
    cfd::ScalarLinearSystem system{mesh};
    seed_system(system);

    assembler.add_internal_face_contributions(provisional_mass_flux, face_pressure_response, system);

    require_near(system.diagonal()[adjacency.owner], 11.0 + static_cast<double>(adjacency.owner) + pressure_response,
                 0.0, "Assembly overwrote the owner diagonal.");
    require_near(system.diagonal()[adjacency.neighbor],
                 11.0 + static_cast<double>(adjacency.neighbor) + pressure_response, 0.0,
                 "Assembly overwrote the neighbor diagonal.");
    require_near(system.owner_neighbor_coefficients()[face_id], 31.0 + static_cast<double>(face_id) - pressure_response,
                 0.0, "Assembly overwrote the owner-neighbor coefficient.");
    require_near(system.neighbor_owner_coefficients()[face_id], 41.0 + static_cast<double>(face_id) - pressure_response,
                 0.0, "Assembly overwrote the neighbor-owner coefficient.");
    require_near(system.rhs()[adjacency.owner], 21.0 + static_cast<double>(adjacency.owner) - provisional_flux, 0.0,
                 "Assembly overwrote the owner RHS.");
    require_near(system.rhs()[adjacency.neighbor], 21.0 + static_cast<double>(adjacency.neighbor) + provisional_flux,
                 0.0, "Assembly overwrote the neighbor RHS.");
}

void test_boundary_values_are_ignored_and_boundary_coefficients_unchanged()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_triangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressiblePressureCorrectionAssembler assembler{mesh};
    cfd::FaceFluxField provisional_mass_flux{mesh.face_count()};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count()};
    const cfd::Index internal_id{internal_face_id(mesh)};
    const cfd::FaceAdjacency &internal_adjacency{mesh.face_adjacencies()[internal_id]};
    provisional_mass_flux[internal_id] = 1.0;
    face_pressure_response[internal_id] = 2.0;
    cfd::ScalarLinearSystem system{mesh};
    seed_system(system);

    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (mesh.face_adjacencies()[face_id].is_boundary())
        {
            provisional_mass_flux[face_id] = std::numeric_limits<double>::infinity();
            face_pressure_response[face_id] = -1.0;
        }
    }

    assembler.add_internal_face_contributions(provisional_mass_flux, face_pressure_response, system);

    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        const bool is_owner{cell_id == internal_adjacency.owner};
        const bool is_neighbor{cell_id == internal_adjacency.neighbor};
        const double expected_diagonal{11.0 + static_cast<double>(cell_id) + (is_owner || is_neighbor ? 2.0 : 0.0)};
        const double expected_rhs{21.0 + static_cast<double>(cell_id) + (is_owner ? -1.0 : is_neighbor ? 1.0 : 0.0)};
        require(system.diagonal()[cell_id] == expected_diagonal, "A boundary value contributed to a cell diagonal.");
        require(system.rhs()[cell_id] == expected_rhs, "A boundary provisional flux contributed to a cell RHS.");
    }

    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (!mesh.face_adjacencies()[face_id].is_boundary())
        {
            continue;
        }
        require(std::isinf(provisional_mass_flux[face_id]) && provisional_mass_flux[face_id] > 0.0,
                "Assembly modified a boundary provisional flux.");
        require(face_pressure_response[face_id] == -1.0, "Assembly modified a boundary pressure response.");
        require(system.owner_neighbor_coefficients()[face_id] == 31.0 + static_cast<double>(face_id),
                "Assembly modified a boundary owner-neighbor coefficient.");
        require(system.neighbor_owner_coefficients()[face_id] == 41.0 + static_cast<double>(face_id),
                "Assembly modified a boundary neighbor-owner coefficient.");
    }
}

void test_rejects_incompatible_inputs_before_mutation()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_triangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressiblePressureCorrectionAssembler assembler{mesh};
    const cfd::FaceFluxField provisional_mass_flux{mesh.face_count()};
    const cfd::FacePressureResponseField face_pressure_response{mesh.face_count(), 1.0};
    cfd::ScalarLinearSystem system{mesh};
    seed_system(system);

    const cfd::FaceFluxField wrong_flux{mesh.face_count() + 1};
    require_rejected_without_mutation<std::invalid_argument>(
        [&]() { assembler.add_internal_face_contributions(wrong_flux, face_pressure_response, system); }, system,
        "Pressure-correction assembly accepted a provisional flux with incorrect cardinality.");

    const cfd::FacePressureResponseField wrong_response{mesh.face_count() + 1};
    require_rejected_without_mutation<std::invalid_argument>(
        [&]() { assembler.add_internal_face_contributions(provisional_mass_flux, wrong_response, system); }, system,
        "Pressure-correction assembly accepted a face response with incorrect cardinality.");

    cfd::MeshBuildResult other_build_result{cfd::build_mesh(make_two_triangle_raw_mesh())};
    cfd::ScalarLinearSystem other_system{other_build_result.mesh};
    seed_system(other_system);
    require_rejected_without_mutation<std::invalid_argument>(
        [&]() {
            assembler.add_internal_face_contributions(provisional_mass_flux, face_pressure_response, other_system);
        },
        other_system, "Pressure-correction assembly accepted a system referencing another Mesh.");
}

void test_rejects_nonfinite_internal_flux_before_mutation()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_triangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressiblePressureCorrectionAssembler assembler{mesh};
    cfd::FaceFluxField provisional_mass_flux{mesh.face_count()};
    const cfd::FacePressureResponseField face_pressure_response{mesh.face_count(), 1.0};
    cfd::ScalarLinearSystem system{mesh};
    seed_system(system);
    const cfd::Index face_id{internal_face_id(mesh)};

    for (const double invalid_flux : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
                                      -std::numeric_limits<double>::infinity()})
    {
        provisional_mass_flux[face_id] = invalid_flux;
        require_rejected_without_mutation<std::runtime_error>(
            [&]() { assembler.add_internal_face_contributions(provisional_mass_flux, face_pressure_response, system); },
            system, "Pressure-correction assembly accepted a non-finite internal provisional flux.");
    }
}

void test_rejects_invalid_internal_pressure_response_before_mutation()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_triangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressiblePressureCorrectionAssembler assembler{mesh};
    const cfd::FaceFluxField provisional_mass_flux{mesh.face_count()};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count(), 1.0};
    cfd::ScalarLinearSystem system{mesh};
    seed_system(system);
    const cfd::Index face_id{internal_face_id(mesh)};

    for (const double invalid_response :
         {0.0, -1.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
          -std::numeric_limits<double>::infinity()})
    {
        face_pressure_response[face_id] = invalid_response;
        require_rejected_without_mutation<std::runtime_error>(
            [&]() { assembler.add_internal_face_contributions(provisional_mass_flux, face_pressure_response, system); },
            system, "Pressure-correction assembly accepted an invalid internal pressure response.");
    }
}

void test_later_invalid_internal_face_is_transactional()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_three_cell_strip_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressiblePressureCorrectionAssembler assembler{mesh};
    cfd::FaceFluxField provisional_mass_flux{mesh.face_count()};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count()};
    const std::array internal_ids{two_internal_face_ids(mesh)};
    provisional_mass_flux[internal_ids[0]] = 3.0;
    face_pressure_response[internal_ids[0]] = 2.0;
    provisional_mass_flux[internal_ids[1]] = -4.0;
    face_pressure_response[internal_ids[1]] = 0.0;
    cfd::ScalarLinearSystem system{mesh};
    seed_system(system);

    require_rejected_without_mutation<std::runtime_error>(
        [&]() { assembler.add_internal_face_contributions(provisional_mass_flux, face_pressure_response, system); },
        system, "Pressure-correction assembly mutated the system before rejecting a later invalid face.");
}

void test_boundary_flux_sign_convention()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_triangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressiblePressureCorrectionAssembler assembler{mesh};
    const cfd::Index face_id{two_boundary_face_ids(mesh)[0]};
    const cfd::Index owner{mesh.face_adjacencies()[face_id].owner};

    for (const double boundary_flux : {3.0, -3.0})
    {
        cfd::FaceFluxField provisional_mass_flux{mesh.face_count()};
        provisional_mass_flux[face_id] = boundary_flux;
        cfd::ScalarLinearSystem system{mesh};

        assembler.add_boundary_provisional_flux_rhs(provisional_mass_flux, system);

        require_near(system.rhs()[owner], -boundary_flux, 0.0,
                     "Boundary provisional flux produced an incorrect owner RHS contribution.");
    }
}

void test_boundary_fluxes_accumulate_additively_without_changing_matrix()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_triangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressiblePressureCorrectionAssembler assembler{mesh};
    cfd::FaceFluxField provisional_mass_flux{mesh.face_count()};
    cfd::ScalarLinearSystem system{mesh};
    seed_system(system);

    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (mesh.face_adjacencies()[face_id].is_boundary())
        {
            provisional_mass_flux[face_id] = 0.5 + static_cast<double>(face_id);
        }
    }

    assembler.add_boundary_provisional_flux_rhs(provisional_mass_flux, system);

    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        double expected_rhs{21.0 + static_cast<double>(cell_id)};
        for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
        {
            const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};
            if (adjacency.is_boundary() && adjacency.owner == cell_id)
            {
                expected_rhs -= provisional_mass_flux[face_id];
            }
        }
        require_near(system.rhs()[cell_id], expected_rhs, 0.0,
                     "Boundary provisional fluxes did not accumulate on the seeded RHS.");
    }
    require_seeded_matrix_unchanged(system, "Boundary provisional-flux assembly");
}

void test_internal_flux_values_are_ignored_by_boundary_assembly()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_triangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressiblePressureCorrectionAssembler assembler{mesh};
    cfd::FaceFluxField provisional_mass_flux{mesh.face_count()};
    cfd::ScalarLinearSystem system{mesh};
    seed_system(system);
    const cfd::Index face_id{internal_face_id(mesh)};

    for (const double invalid_flux : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
                                      -std::numeric_limits<double>::infinity()})
    {
        provisional_mass_flux[face_id] = invalid_flux;
        assembler.add_boundary_provisional_flux_rhs(provisional_mass_flux, system);
        require_seeded_system_unchanged(system, "Boundary assembly consumed an internal provisional flux.");
    }
}

void test_boundary_assembly_rejects_incompatible_inputs_before_mutation()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_triangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressiblePressureCorrectionAssembler assembler{mesh};
    const cfd::FaceFluxField provisional_mass_flux{mesh.face_count()};
    cfd::ScalarLinearSystem system{mesh};
    seed_system(system);

    const cfd::FaceFluxField wrong_flux{mesh.face_count() + 1};
    require_rejected_without_mutation<std::invalid_argument>(
        [&]() { assembler.add_boundary_provisional_flux_rhs(wrong_flux, system); }, system,
        "Boundary pressure-correction assembly accepted a provisional flux with incorrect cardinality.");

    cfd::MeshBuildResult other_build_result{cfd::build_mesh(make_two_triangle_raw_mesh())};
    cfd::ScalarLinearSystem other_system{other_build_result.mesh};
    seed_system(other_system);
    require_rejected_without_mutation<std::invalid_argument>(
        [&]() { assembler.add_boundary_provisional_flux_rhs(provisional_mass_flux, other_system); }, other_system,
        "Boundary pressure-correction assembly accepted a system referencing another Mesh.");
}

void test_rejects_nonfinite_boundary_flux_before_mutation()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_triangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressiblePressureCorrectionAssembler assembler{mesh};
    cfd::FaceFluxField provisional_mass_flux{mesh.face_count()};
    cfd::ScalarLinearSystem system{mesh};
    seed_system(system);
    const cfd::Index face_id{two_boundary_face_ids(mesh)[0]};

    for (const double invalid_flux : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
                                      -std::numeric_limits<double>::infinity()})
    {
        provisional_mass_flux[face_id] = invalid_flux;
        require_rejected_without_mutation<std::runtime_error>(
            [&]() { assembler.add_boundary_provisional_flux_rhs(provisional_mass_flux, system); }, system,
            "Boundary pressure-correction assembly accepted a non-finite boundary provisional flux.");
    }
}

void test_later_invalid_boundary_face_is_transactional()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_triangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressiblePressureCorrectionAssembler assembler{mesh};
    cfd::FaceFluxField provisional_mass_flux{mesh.face_count()};
    const std::array boundary_ids{two_boundary_face_ids(mesh)};
    provisional_mass_flux[boundary_ids[0]] = 3.0;
    provisional_mass_flux[boundary_ids[1]] = std::numeric_limits<double>::infinity();
    cfd::ScalarLinearSystem system{mesh};
    seed_system(system);

    require_rejected_without_mutation<std::runtime_error>(
        [&]() { assembler.add_boundary_provisional_flux_rhs(provisional_mass_flux, system); }, system,
        "Boundary pressure-correction assembly mutated the system before rejecting a later invalid face.");
}

void test_combined_provisional_flux_rhs_is_negative_cell_mass_imbalance()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_triangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressiblePressureCorrectionAssembler assembler{mesh};
    cfd::FaceFluxField provisional_mass_flux{mesh.face_count()};
    const cfd::FacePressureResponseField face_pressure_response{mesh.face_count(), 2.0};
    cfd::ScalarLinearSystem system{mesh};
    seed_system(system);

    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        provisional_mass_flux[face_id] = 0.25 + static_cast<double>(face_id);
    }

    assembler.add_internal_face_contributions(provisional_mass_flux, face_pressure_response, system);
    assembler.add_boundary_provisional_flux_rhs(provisional_mass_flux, system);

    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        double outward_mass_flux{};
        for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
        {
            const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};
            if (adjacency.owner == cell_id)
            {
                outward_mass_flux += provisional_mass_flux[face_id];
            }
            else if (!adjacency.is_boundary() && adjacency.neighbor == cell_id)
            {
                outward_mass_flux -= provisional_mass_flux[face_id];
            }
        }
        const double expected_rhs{21.0 + static_cast<double>(cell_id) - outward_mass_flux};
        require_near(system.rhs()[cell_id], expected_rhs, 0.0,
                     "Combined pressure-correction assembly did not produce the negative cell mass imbalance.");
    }
}

void test_pressure_correction_boundary_condition_collection_validation()
{
    using cfd::PressureCorrectionBoundaryConditionType;

    const cfd::PressureCorrectionBoundaryConditions boundary_conditions{
        2,
        {PressureCorrectionBoundaryConditionType::FixedMassFlux,
         PressureCorrectionBoundaryConditionType::FixedPressure},
    };
    require(boundary_conditions.size() == 2,
            "Pressure-correction boundary conditions contain an incorrect number of values.");
    require(boundary_conditions[0] == PressureCorrectionBoundaryConditionType::FixedMassFlux,
            "Pressure-correction boundary condition 0 has an incorrect type.");
    require(boundary_conditions[1] == PressureCorrectionBoundaryConditionType::FixedPressure,
            "Pressure-correction boundary condition 1 has an incorrect type.");

    require_throws<std::invalid_argument>(
        []() {
            static_cast<void>(cfd::PressureCorrectionBoundaryConditions{
                2,
                {cfd::PressureCorrectionBoundaryConditionType::FixedMassFlux},
            });
        },
        "Pressure-correction boundary conditions accepted an incorrect condition count.");

    require_throws<std::invalid_argument>(
        []() {
            // NOLINTNEXTLINE(clang-analyzer-optin.core.EnumCastOutOfRange)
            const auto invalid_type{static_cast<cfd::PressureCorrectionBoundaryConditionType>(255)};
            static_cast<void>(cfd::PressureCorrectionBoundaryConditions{1, {invalid_type}});
        },
        "Pressure-correction boundary conditions accepted an unsupported condition type.");
}

void test_fixed_mass_flux_boundaries_do_not_modify_system_or_read_responses()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_triangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressiblePressureCorrectionAssembler assembler{mesh};
    const cfd::PressureCorrectionBoundaryConditions boundary_conditions{
        make_uniform_pressure_correction_boundary_conditions(
            mesh, cfd::PressureCorrectionBoundaryConditionType::FixedMassFlux)};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count()};
    const std::array invalid_responses{
        0.0,
        -1.0,
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
    };
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        face_pressure_response[face_id] = invalid_responses.at(face_id % invalid_responses.size());
    }
    cfd::ScalarLinearSystem system{mesh};
    seed_system(system);

    assembler.add_boundary_pressure_response(boundary_conditions, face_pressure_response, system);

    require_seeded_system_unchanged(
        system, "FixedMassFlux boundary assembly consumed an unused boundary or internal response.");
}

void test_mixed_boundaries_add_only_fixed_pressure_response()
{
    constexpr cfd::BoundaryId fixed_pressure_boundary_id{1};

    cfd::MeshBuildResult build_result{cfd::build_mesh(make_four_boundary_group_two_triangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressiblePressureCorrectionAssembler assembler{mesh};
    std::vector<cfd::PressureCorrectionBoundaryConditionType> conditions(
        mesh.boundary_groups().size(), cfd::PressureCorrectionBoundaryConditionType::FixedMassFlux);
    conditions[fixed_pressure_boundary_id] = cfd::PressureCorrectionBoundaryConditionType::FixedPressure;
    const cfd::PressureCorrectionBoundaryConditions boundary_conditions{mesh.boundary_groups().size(),
                                                                        std::move(conditions)};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count(), std::numeric_limits<double>::quiet_NaN()};
    cfd::Index fixed_pressure_owner{cfd::invalid_index};
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (mesh.face_boundary_ids()[face_id] == fixed_pressure_boundary_id)
        {
            face_pressure_response[face_id] = 2.0;
            fixed_pressure_owner = mesh.face_adjacencies()[face_id].owner;
        }
    }
    require(fixed_pressure_owner != cfd::invalid_index, "Mixed-boundary fixture has no FixedPressure face.");
    cfd::ScalarLinearSystem system{mesh};
    seed_system(system);

    assembler.add_boundary_pressure_response(boundary_conditions, face_pressure_response, system);

    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        const double expected_diagonal{11.0 + static_cast<double>(cell_id) +
                                       (cell_id == fixed_pressure_owner ? 2.0 : 0.0)};
        require(system.diagonal()[cell_id] == expected_diagonal,
                "Mixed boundary assembly added an incorrect diagonal response.");
    }
    require_seeded_rhs_and_off_diagonal_coefficients_unchanged(system, "Mixed boundary pressure-response assembly");
}

void test_multiple_fixed_pressure_faces_accumulate_on_owner_diagonals()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_triangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressiblePressureCorrectionAssembler assembler{mesh};
    const cfd::PressureCorrectionBoundaryConditions boundary_conditions{
        make_uniform_pressure_correction_boundary_conditions(
            mesh, cfd::PressureCorrectionBoundaryConditionType::FixedPressure)};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count(), std::numeric_limits<double>::quiet_NaN()};
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (mesh.face_adjacencies()[face_id].is_boundary())
        {
            face_pressure_response[face_id] = 2.0;
        }
    }
    cfd::ScalarLinearSystem system{mesh};
    seed_system(system);

    assembler.add_boundary_pressure_response(boundary_conditions, face_pressure_response, system);

    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        double expected_diagonal{11.0 + static_cast<double>(cell_id)};
        cfd::Index contributing_face_count{};
        for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
        {
            const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};
            if (adjacency.is_boundary() && adjacency.owner == cell_id)
            {
                expected_diagonal += 2.0;
                ++contributing_face_count;
            }
        }
        require(contributing_face_count > 1,
                "Multiple-response fixture does not attach multiple boundary faces to the cell.");
        require(system.diagonal()[cell_id] == expected_diagonal,
                "Multiple FixedPressure responses did not accumulate on the owner diagonal.");
    }
    require_seeded_rhs_and_off_diagonal_coefficients_unchanged(system, "Multiple FixedPressure boundary assembly");
}

void test_rejects_invalid_fixed_pressure_response_before_mutation()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_triangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressiblePressureCorrectionAssembler assembler{mesh};
    const cfd::PressureCorrectionBoundaryConditions boundary_conditions{
        make_uniform_pressure_correction_boundary_conditions(
            mesh, cfd::PressureCorrectionBoundaryConditionType::FixedPressure)};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count(), 1.0};
    const cfd::Index face_id{two_boundary_face_ids(mesh)[0]};
    cfd::ScalarLinearSystem system{mesh};
    seed_system(system);

    for (const double invalid_response :
         {0.0, -1.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
          -std::numeric_limits<double>::infinity()})
    {
        face_pressure_response[face_id] = invalid_response;
        require_rejected_without_mutation<std::runtime_error>(
            [&]() { assembler.add_boundary_pressure_response(boundary_conditions, face_pressure_response, system); },
            system, "Boundary pressure-response assembly accepted an invalid FixedPressure response.");
    }
}

void test_boundary_pressure_response_rejects_incompatible_inputs_before_mutation()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_triangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressiblePressureCorrectionAssembler assembler{mesh};
    const cfd::PressureCorrectionBoundaryConditions boundary_conditions{
        make_uniform_pressure_correction_boundary_conditions(
            mesh, cfd::PressureCorrectionBoundaryConditionType::FixedMassFlux)};
    const cfd::FacePressureResponseField face_pressure_response{mesh.face_count()};
    cfd::ScalarLinearSystem system{mesh};
    seed_system(system);

    const cfd::FacePressureResponseField wrong_response{mesh.face_count() + 1};
    require_rejected_without_mutation<std::invalid_argument>(
        [&]() { assembler.add_boundary_pressure_response(boundary_conditions, wrong_response, system); }, system,
        "Boundary pressure-response assembly accepted a face response with incorrect cardinality.");

    const cfd::PressureCorrectionBoundaryConditions wrong_boundary_conditions{
        mesh.boundary_groups().size() + 1,
        std::vector<cfd::PressureCorrectionBoundaryConditionType>(
            mesh.boundary_groups().size() + 1, cfd::PressureCorrectionBoundaryConditionType::FixedMassFlux)};
    require_rejected_without_mutation<std::invalid_argument>(
        [&]() { assembler.add_boundary_pressure_response(wrong_boundary_conditions, face_pressure_response, system); },
        system, "Boundary pressure-response assembly accepted an incorrect boundary-condition cardinality.");

    cfd::MeshBuildResult other_build_result{cfd::build_mesh(make_two_triangle_raw_mesh())};
    cfd::ScalarLinearSystem other_system{other_build_result.mesh};
    seed_system(other_system);
    require_rejected_without_mutation<std::invalid_argument>(
        [&]() { assembler.add_boundary_pressure_response(boundary_conditions, face_pressure_response, other_system); },
        other_system, "Boundary pressure-response assembly accepted a system referencing another Mesh.");
}

void test_later_invalid_fixed_pressure_face_is_transactional()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_triangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressiblePressureCorrectionAssembler assembler{mesh};
    const cfd::PressureCorrectionBoundaryConditions boundary_conditions{
        make_uniform_pressure_correction_boundary_conditions(
            mesh, cfd::PressureCorrectionBoundaryConditionType::FixedPressure)};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count(), 1.0};
    const std::array boundary_ids{two_boundary_face_ids(mesh)};
    face_pressure_response[boundary_ids[0]] = 2.0;
    face_pressure_response[boundary_ids[1]] = std::numeric_limits<double>::infinity();
    cfd::ScalarLinearSystem system{mesh};
    seed_system(system);

    require_rejected_without_mutation<std::runtime_error>(
        [&]() { assembler.add_boundary_pressure_response(boundary_conditions, face_pressure_response, system); },
        system, "Boundary pressure-response assembly mutated the system before rejecting a later invalid face.");
}

} // namespace

int main()
{
    int failure_count{};

    failure_count += cfd::test::run_test("pressure-correction exact internal-face assembly",
                                         test_exact_internal_face_contributions_and_local_conservation);
    failure_count += cfd::test::run_test("pressure-correction negative provisional flux",
                                         test_negative_provisional_flux_reverses_rhs_contributions);
    failure_count += cfd::test::run_test("pressure-correction additive assembly", test_assembly_is_additive);
    failure_count += cfd::test::run_test("pressure-correction boundary isolation",
                                         test_boundary_values_are_ignored_and_boundary_coefficients_unchanged);
    failure_count +=
        cfd::test::run_test("pressure-correction API validation", test_rejects_incompatible_inputs_before_mutation);
    failure_count += cfd::test::run_test("pressure-correction provisional-flux validation",
                                         test_rejects_nonfinite_internal_flux_before_mutation);
    failure_count += cfd::test::run_test("pressure-correction face-response validation",
                                         test_rejects_invalid_internal_pressure_response_before_mutation);
    failure_count += cfd::test::run_test("pressure-correction transactional validation",
                                         test_later_invalid_internal_face_is_transactional);
    failure_count += cfd::test::run_test("pressure-correction boundary flux sign", test_boundary_flux_sign_convention);
    failure_count += cfd::test::run_test("pressure-correction boundary additive assembly",
                                         test_boundary_fluxes_accumulate_additively_without_changing_matrix);
    failure_count += cfd::test::run_test("pressure-correction boundary internal-face isolation",
                                         test_internal_flux_values_are_ignored_by_boundary_assembly);
    failure_count += cfd::test::run_test("pressure-correction boundary API validation",
                                         test_boundary_assembly_rejects_incompatible_inputs_before_mutation);
    failure_count += cfd::test::run_test("pressure-correction boundary flux validation",
                                         test_rejects_nonfinite_boundary_flux_before_mutation);
    failure_count += cfd::test::run_test("pressure-correction boundary transactional validation",
                                         test_later_invalid_boundary_face_is_transactional);
    failure_count += cfd::test::run_test("pressure-correction combined provisional-flux RHS",
                                         test_combined_provisional_flux_rhs_is_negative_cell_mass_imbalance);
    failure_count += cfd::test::run_test("pressure-correction boundary-condition collection",
                                         test_pressure_correction_boundary_condition_collection_validation);
    failure_count += cfd::test::run_test("pressure-correction FixedMassFlux response isolation",
                                         test_fixed_mass_flux_boundaries_do_not_modify_system_or_read_responses);
    failure_count += cfd::test::run_test("pressure-correction mixed boundary pressure response",
                                         test_mixed_boundaries_add_only_fixed_pressure_response);
    failure_count += cfd::test::run_test("pressure-correction multiple FixedPressure responses",
                                         test_multiple_fixed_pressure_faces_accumulate_on_owner_diagonals);
    failure_count += cfd::test::run_test("pressure-correction FixedPressure response validation",
                                         test_rejects_invalid_fixed_pressure_response_before_mutation);
    failure_count += cfd::test::run_test("pressure-correction boundary pressure-response API validation",
                                         test_boundary_pressure_response_rejects_incompatible_inputs_before_mutation);
    failure_count += cfd::test::run_test("pressure-correction boundary pressure-response transaction",
                                         test_later_invalid_fixed_pressure_face_is_transactional);

    return cfd::test::finish_tests(failure_count, "incompressible pressure correction assembler");
}
