#include "cfd/numerics/IncompressiblePressureVelocityCorrection.hpp"

#include "cfd/field/CellMomentumPressureResponse.hpp"
#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVectorField.hpp"
#include "cfd/field/CellVelocityField.hpp"
#include "cfd/field/FaceFluxField.hpp"
#include "cfd/field/FacePressureResponseField.hpp"
#include "cfd/field/PressureCorrectionBoundaryConditions.hpp"
#include "cfd/linear_algebra/EigenConjugateGradientSolver.hpp"
#include "cfd/linear_algebra/ScalarLinearSystem.hpp"
#include "cfd/math/Vector2.hpp"
#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/mesh/Types.hpp"
#include "cfd/meshing/RawMeshData.hpp"
#include "cfd/numerics/IncompressiblePressureCorrectionAssembler.hpp"
#include "cfd/numerics/MassFluxBalance.hpp"

#include "support/TestUtils.hpp"

#include <algorithm>
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

using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws;
using cfd::test::test_tolerance;

constexpr cfd::BoundaryId bottom_boundary_id{0};
constexpr cfd::BoundaryId right_boundary_id{1};
constexpr cfd::BoundaryId top_boundary_id{2};
constexpr cfd::BoundaryId left_boundary_id{3};

static_assert(std::is_nothrow_constructible_v<cfd::IncompressiblePressureVelocityCorrection, const cfd::Mesh &>);
static_assert(!std::is_copy_constructible_v<cfd::IncompressiblePressureVelocityCorrection>);
static_assert(!std::is_copy_assignable_v<cfd::IncompressiblePressureVelocityCorrection>);
static_assert(std::is_nothrow_move_constructible_v<cfd::IncompressiblePressureVelocityCorrection>);
static_assert(!std::is_move_assignable_v<cfd::IncompressiblePressureVelocityCorrection>);

[[nodiscard]]
cfd::RawMeshData make_two_cell_rectangle_raw_mesh()
{
    cfd::RawMeshData raw_mesh;
    raw_mesh.nodes = {
        {0.0, 0.0}, {1.0, 0.0}, {2.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}, {2.0, 1.0},
    };
    raw_mesh.cell_types = {cfd::CellType::Quadrilateral, cfd::CellType::Quadrilateral};
    raw_mesh.cell_nodes = {0, 1, 4, 3, 1, 2, 5, 4};
    raw_mesh.cell_node_offsets = {0, 4, 8};
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
cfd::RawMeshData make_single_oblique_triangle_raw_mesh()
{
    constexpr cfd::BoundaryId wall_boundary_id{0};

    cfd::RawMeshData raw_mesh;
    raw_mesh.nodes = {
        {0.0, 0.0},
        {2.0, 0.0},
        {0.0, 1.0},
    };
    raw_mesh.cell_types = {cfd::CellType::Triangle};
    raw_mesh.cell_nodes = {0, 1, 2};
    raw_mesh.cell_node_offsets = {0, 3};
    raw_mesh.boundary_groups = {{wall_boundary_id, "wall"}};
    raw_mesh.boundary_edges = {
        {{0, 1}, wall_boundary_id},
        {{1, 2}, wall_boundary_id},
        {{2, 0}, wall_boundary_id},
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
    throw std::runtime_error("Pressure-velocity correction fixture has no internal face.");
}

[[nodiscard]]
std::array<cfd::Index, 2> internal_face_ids(const cfd::Mesh &mesh)
{
    std::array<cfd::Index, 2> result{cfd::invalid_index, cfd::invalid_index};
    cfd::Index count{};
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (mesh.face_adjacencies()[face_id].is_boundary())
        {
            continue;
        }
        require(count < result.size(), "Three-cell correction fixture has too many internal faces.");
        result[count] = face_id;
        ++count;
    }
    require(count == result.size(), "Three-cell correction fixture must have two internal faces.");
    return result;
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
    throw std::runtime_error("Pressure-velocity correction fixture has no requested boundary face.");
}

[[nodiscard]]
cfd::Index oblique_boundary_face_id(const cfd::Mesh &mesh)
{
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        const cfd::Vector2 &area_vector{mesh.face_area_vectors()[face_id]};
        if (mesh.face_adjacencies()[face_id].is_boundary() && area_vector.x != 0.0 && area_vector.y != 0.0)
        {
            return face_id;
        }
    }
    throw std::runtime_error("Pressure-velocity correction fixture has no oblique boundary face.");
}

[[nodiscard]]
cfd::PressureCorrectionBoundaryConditions make_uniform_boundary_conditions(
    const cfd::Mesh &mesh, const cfd::PressureCorrectionBoundaryConditionType condition)
{
    return {mesh.boundary_groups().size(),
            std::vector<cfd::PressureCorrectionBoundaryConditionType>(mesh.boundary_groups().size(), condition)};
}

[[nodiscard]]
cfd::PressureCorrectionBoundaryConditions make_right_fixed_pressure_boundary_conditions(const cfd::Mesh &mesh)
{
    std::vector<cfd::PressureCorrectionBoundaryConditionType> conditions(
        mesh.boundary_groups().size(), cfd::PressureCorrectionBoundaryConditionType::FixedMassFlux);
    conditions[right_boundary_id] = cfd::PressureCorrectionBoundaryConditionType::FixedPressure;
    return {mesh.boundary_groups().size(), std::move(conditions)};
}

void seed_flux(cfd::FaceFluxField &mass_flux)
{
    for (cfd::Index face_id = 0; face_id < mass_flux.size(); ++face_id)
    {
        mass_flux[face_id] = 10.0 + static_cast<double>(face_id);
    }
}

void require_seeded_flux_unchanged(const cfd::FaceFluxField &mass_flux, const std::string &context)
{
    for (cfd::Index face_id = 0; face_id < mass_flux.size(); ++face_id)
    {
        require(mass_flux[face_id] == 10.0 + static_cast<double>(face_id), context);
    }
}

void test_internal_face_correction_sign_and_additive_semantics()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressiblePressureVelocityCorrection corrector{mesh};
    const cfd::PressureCorrectionBoundaryConditions boundary_conditions{
        make_uniform_boundary_conditions(mesh, cfd::PressureCorrectionBoundaryConditionType::FixedMassFlux)};
    const cfd::Index face_id{internal_face_id(mesh)};
    const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};
    cfd::CellScalarField pressure_correction{mesh.cell_count()};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count(), 1.0};
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    pressure_correction[adjacency.owner] = 3.0;
    pressure_correction[adjacency.neighbor] = 1.0;
    face_pressure_response[face_id] = 2.0;
    mass_flux[face_id] = 5.0;

    corrector.correct_face_mass_flux(pressure_correction, boundary_conditions, face_pressure_response, mass_flux);
    require_near(mass_flux[face_id], 9.0, 0.0, "Positive internal pressure difference has the wrong flux sign.");

    corrector.correct_face_mass_flux(pressure_correction, boundary_conditions, face_pressure_response, mass_flux);
    require_near(mass_flux[face_id], 13.0, 0.0, "Repeated face correction is not additive and in place.");

    pressure_correction[adjacency.owner] = 1.0;
    pressure_correction[adjacency.neighbor] = 3.0;
    mass_flux[face_id] = 5.0;
    corrector.correct_face_mass_flux(pressure_correction, boundary_conditions, face_pressure_response, mass_flux);
    require_near(mass_flux[face_id], 1.0, 0.0, "Negative internal pressure difference has the wrong flux sign.");
}

void test_boundary_face_correction_semantics_and_isolation()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressiblePressureVelocityCorrection corrector{mesh};
    const cfd::PressureCorrectionBoundaryConditions boundary_conditions{
        make_right_fixed_pressure_boundary_conditions(mesh)};
    cfd::CellScalarField pressure_correction{mesh.cell_count(), 1.0};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count(), 1.0};
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    seed_flux(mass_flux);
    const cfd::Index internal_id{internal_face_id(mesh)};
    const cfd::Index right_face_id{boundary_face_id(mesh, right_boundary_id)};
    const cfd::Index left_face_id{boundary_face_id(mesh, left_boundary_id)};
    const double right_initial_flux{mass_flux[right_face_id]};
    const double left_initial_flux{mass_flux[left_face_id]};
    face_pressure_response[right_face_id] = 2.0;
    face_pressure_response[left_face_id] = std::numeric_limits<double>::quiet_NaN();

    corrector.correct_face_mass_flux(pressure_correction, boundary_conditions, face_pressure_response, mass_flux);

    require_near(mass_flux[right_face_id], right_initial_flux + 2.0, 0.0,
                 "FixedPressure boundary correction is incorrect.");
    require(mass_flux[left_face_id] == left_initial_flux,
            "FixedMassFlux boundary changed or its NaN response was consumed.");
    require(mass_flux[internal_id] == 10.0 + static_cast<double>(internal_id),
            "Zero internal pressure difference changed the internal mass flux.");
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (!mesh.face_adjacencies()[face_id].is_boundary() || face_id == right_face_id)
        {
            continue;
        }
        require(mass_flux[face_id] == 10.0 + static_cast<double>(face_id),
                "A FixedMassFlux boundary entry was modified.");
    }
}

void test_face_correction_rejects_incompatible_cardinalities_before_mutation()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressiblePressureVelocityCorrection corrector{mesh};
    const cfd::PressureCorrectionBoundaryConditions boundary_conditions{
        make_right_fixed_pressure_boundary_conditions(mesh)};
    const cfd::PressureCorrectionBoundaryConditions wrong_boundary_conditions{
        1, {cfd::PressureCorrectionBoundaryConditionType::FixedMassFlux}};
    cfd::CellScalarField pressure_correction{mesh.cell_count()};
    cfd::CellScalarField wrong_pressure_correction{mesh.cell_count() - 1};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count(), 1.0};
    cfd::FacePressureResponseField wrong_face_pressure_response{mesh.face_count() - 1, 1.0};
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    cfd::FaceFluxField wrong_mass_flux{mesh.face_count() - 1};
    seed_flux(mass_flux);
    seed_flux(wrong_mass_flux);

    require_throws<std::invalid_argument>(
        [&]() {
            corrector.correct_face_mass_flux(wrong_pressure_correction, boundary_conditions, face_pressure_response,
                                             mass_flux);
        },
        "Face correction accepted the wrong pressure-correction cardinality.");
    require_seeded_flux_unchanged(mass_flux, "Wrong pressure-correction cardinality modified face fluxes.");
    require_throws<std::invalid_argument>(
        [&]() {
            corrector.correct_face_mass_flux(pressure_correction, wrong_boundary_conditions, face_pressure_response,
                                             mass_flux);
        },
        "Face correction accepted the wrong boundary-condition cardinality.");
    require_seeded_flux_unchanged(mass_flux, "Wrong boundary-condition cardinality modified face fluxes.");
    require_throws<std::invalid_argument>(
        [&]() {
            corrector.correct_face_mass_flux(pressure_correction, boundary_conditions, wrong_face_pressure_response,
                                             mass_flux);
        },
        "Face correction accepted the wrong response cardinality.");
    require_seeded_flux_unchanged(mass_flux, "Wrong response cardinality modified face fluxes.");
    require_throws<std::invalid_argument>(
        [&]() {
            corrector.correct_face_mass_flux(pressure_correction, boundary_conditions, face_pressure_response,
                                             wrong_mass_flux);
        },
        "Face correction accepted the wrong mass-flux cardinality.");
    require_seeded_flux_unchanged(wrong_mass_flux, "Wrong mass-flux cardinality modified face fluxes.");
}

void test_face_correction_rejects_invalid_used_values_before_mutation()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressiblePressureVelocityCorrection corrector{mesh};
    const cfd::PressureCorrectionBoundaryConditions boundary_conditions{
        make_right_fixed_pressure_boundary_conditions(mesh)};
    const cfd::Index internal_id{internal_face_id(mesh)};
    const cfd::Index right_face_id{boundary_face_id(mesh, right_boundary_id)};
    const std::array invalid_responses{
        0.0,
        -1.0,
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
    };

    for (const cfd::Index face_id : {internal_id, right_face_id})
    {
        for (const double invalid_response : invalid_responses)
        {
            cfd::CellScalarField pressure_correction{mesh.cell_count(), 1.0};
            cfd::FacePressureResponseField face_pressure_response{mesh.face_count(), 1.0};
            cfd::FaceFluxField mass_flux{mesh.face_count()};
            seed_flux(mass_flux);
            face_pressure_response[face_id] = invalid_response;
            require_throws<std::runtime_error>(
                [&]() {
                    corrector.correct_face_mass_flux(pressure_correction, boundary_conditions, face_pressure_response,
                                                     mass_flux);
                },
                "Face correction accepted an invalid used pressure response.");
            require_seeded_flux_unchanged(mass_flux, "Invalid pressure response partially modified face fluxes.");
        }
    }

    {
        cfd::CellScalarField pressure_correction{mesh.cell_count(), 1.0};
        cfd::FacePressureResponseField face_pressure_response{mesh.face_count(), 1.0};
        cfd::FaceFluxField mass_flux{mesh.face_count()};
        seed_flux(mass_flux);
        pressure_correction[mesh.face_adjacencies()[internal_id].neighbor] = std::numeric_limits<double>::quiet_NaN();
        require_throws<std::runtime_error>(
            [&]() {
                corrector.correct_face_mass_flux(pressure_correction, boundary_conditions, face_pressure_response,
                                                 mass_flux);
            },
            "Face correction accepted a non-finite used pressure correction.");
        require_seeded_flux_unchanged(mass_flux, "Non-finite pressure correction partially modified face fluxes.");
    }

    for (const cfd::Index face_id : {internal_id, right_face_id})
    {
        cfd::CellScalarField pressure_correction{mesh.cell_count(), 1.0};
        cfd::FacePressureResponseField face_pressure_response{mesh.face_count(), 1.0};
        cfd::FaceFluxField mass_flux{mesh.face_count()};
        seed_flux(mass_flux);
        mass_flux[face_id] = std::numeric_limits<double>::infinity();
        require_throws<std::runtime_error>(
            [&]() {
                corrector.correct_face_mass_flux(pressure_correction, boundary_conditions, face_pressure_response,
                                                 mass_flux);
            },
            "Face correction accepted a non-finite corrected-face input flux.");
        for (cfd::Index other_face_id = 0; other_face_id < mesh.face_count(); ++other_face_id)
        {
            if (other_face_id != face_id)
            {
                require(mass_flux[other_face_id] == 10.0 + static_cast<double>(other_face_id),
                        "Non-finite current flux partially modified another face.");
            }
        }
    }
}

void test_face_correction_is_transactional_for_later_invalid_face()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_three_cell_strip_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressiblePressureVelocityCorrection corrector{mesh};
    const cfd::PressureCorrectionBoundaryConditions boundary_conditions{
        make_uniform_boundary_conditions(mesh, cfd::PressureCorrectionBoundaryConditionType::FixedMassFlux)};
    const std::array face_ids{internal_face_ids(mesh)};
    cfd::CellScalarField pressure_correction{mesh.cell_count()};
    pressure_correction[0] = 3.0;
    pressure_correction[1] = 2.0;
    pressure_correction[2] = 1.0;
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count(), 1.0};
    face_pressure_response[face_ids[0]] = 2.0;
    face_pressure_response[face_ids[1]] = std::numeric_limits<double>::infinity();
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    seed_flux(mass_flux);

    require_throws<std::runtime_error>(
        [&]() {
            corrector.correct_face_mass_flux(pressure_correction, boundary_conditions, face_pressure_response,
                                             mass_flux);
        },
        "Face correction accepted a later invalid internal response.");
    require_seeded_flux_unchanged(mass_flux, "Later invalid internal face caused partial flux correction.");
}

void test_pressure_correction_with_full_and_partial_relaxation()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressiblePressureVelocityCorrection corrector{mesh};
    cfd::CellScalarField pressure_correction{mesh.cell_count()};
    pressure_correction[0] = 2.0;
    pressure_correction[1] = -4.0;

    cfd::CellScalarField full_pressure{mesh.cell_count()};
    full_pressure[0] = 10.0;
    full_pressure[1] = 20.0;
    corrector.correct_pressure(pressure_correction, 1.0, full_pressure);
    require_near(full_pressure[0], 12.0, 0.0, "Full pressure correction is incorrect for cell 0.");
    require_near(full_pressure[1], 16.0, 0.0, "Full pressure correction is incorrect for cell 1.");

    cfd::CellScalarField relaxed_pressure{mesh.cell_count()};
    relaxed_pressure[0] = 10.0;
    relaxed_pressure[1] = 20.0;
    corrector.correct_pressure(pressure_correction, 0.25, relaxed_pressure);
    require_near(relaxed_pressure[0], 10.5, 0.0, "Relaxed pressure correction is incorrect for cell 0.");
    require_near(relaxed_pressure[1], 19.0, 0.0, "Relaxed pressure correction is incorrect for cell 1.");
}

void test_pressure_correction_validation_is_transactional()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressiblePressureVelocityCorrection corrector{mesh};
    const std::array invalid_relaxation_factors{
        0.0, -0.5, 1.1, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
    };
    for (const double relaxation_factor : invalid_relaxation_factors)
    {
        cfd::CellScalarField pressure_correction{mesh.cell_count(), 1.0};
        cfd::CellScalarField pressure{mesh.cell_count()};
        pressure[0] = 10.0;
        pressure[1] = 20.0;
        require_throws<std::invalid_argument>(
            [&]() { corrector.correct_pressure(pressure_correction, relaxation_factor, pressure); },
            "Pressure correction accepted an invalid relaxation factor.");
        require(pressure[0] == 10.0 && pressure[1] == 20.0,
                "Invalid pressure relaxation factor partially modified pressure.");
    }

    {
        cfd::CellScalarField pressure_correction{mesh.cell_count(), 1.0};
        pressure_correction[1] = std::numeric_limits<double>::quiet_NaN();
        cfd::CellScalarField pressure{mesh.cell_count()};
        pressure[0] = 10.0;
        pressure[1] = 20.0;
        require_throws<std::runtime_error>([&]() { corrector.correct_pressure(pressure_correction, 0.5, pressure); },
                                           "Pressure correction accepted a non-finite correction value.");
        require(pressure[0] == 10.0 && pressure[1] == 20.0,
                "Non-finite pressure correction partially modified pressure.");
    }

    {
        cfd::CellScalarField pressure_correction{mesh.cell_count(), 1.0};
        cfd::CellScalarField pressure{mesh.cell_count()};
        pressure[0] = 10.0;
        pressure[1] = std::numeric_limits<double>::infinity();
        require_throws<std::runtime_error>([&]() { corrector.correct_pressure(pressure_correction, 0.5, pressure); },
                                           "Pressure correction accepted a non-finite pressure value.");
        require(pressure[0] == 10.0 && std::isinf(pressure[1]), "Non-finite pressure input caused partial mutation.");
    }

    cfd::CellScalarField aliased_pressure{mesh.cell_count(), 1.0};
    require_throws<std::invalid_argument>(
        [&]() { corrector.correct_pressure(aliased_pressure, 0.5, aliased_pressure); },
        "Pressure correction accepted aliased fields.");
    require(aliased_pressure[0] == 1.0 && aliased_pressure[1] == 1.0, "Aliasing rejection modified pressure.");

    cfd::CellScalarField wrong_pressure_correction{mesh.cell_count() - 1};
    cfd::CellScalarField pressure{mesh.cell_count(), 1.0};
    require_throws<std::invalid_argument>(
        [&]() { corrector.correct_pressure(wrong_pressure_correction, 0.5, pressure); },
        "Pressure correction accepted the wrong correction cardinality.");
    cfd::CellScalarField valid_pressure_correction{mesh.cell_count()};
    cfd::CellScalarField wrong_pressure{mesh.cell_count() - 1, 1.0};
    require_throws<std::invalid_argument>(
        [&]() { corrector.correct_pressure(valid_pressure_correction, 0.5, wrong_pressure); },
        "Pressure correction accepted the wrong pressure cardinality.");
}

void test_velocity_correction_is_componentwise_and_unrelaxed()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressiblePressureVelocityCorrection corrector{mesh};
    cfd::CellVectorField pressure_correction_gradient{mesh.cell_count()};
    pressure_correction_gradient[0] = {4.0, 5.0};
    pressure_correction_gradient[1] = {-2.0, 0.25};
    cfd::CellMomentumPressureResponse momentum_response{mesh.cell_count()};
    momentum_response.u()[0] = 2.0;
    momentum_response.v()[0] = 3.0;
    momentum_response.u()[1] = 0.5;
    momentum_response.v()[1] = 4.0;
    cfd::CellVelocityField velocity{mesh.cell_count()};
    velocity.u()[0] = 10.0;
    velocity.v()[0] = 20.0;
    velocity.u()[1] = -1.0;
    velocity.v()[1] = 2.0;

    corrector.correct_velocity(pressure_correction_gradient, momentum_response, velocity);

    require_near(velocity.u()[0], 2.0, 0.0, "Cell 0 u correction is incorrect.");
    require_near(velocity.v()[0], 5.0, 0.0, "Cell 0 v correction is incorrect.");
    require_near(velocity.u()[1], 0.0, 0.0, "Cell 1 u correction is incorrect.");
    require_near(velocity.v()[1], 1.0, 0.0, "Cell 1 v correction is incorrect.");
}

void test_fixed_mass_flux_is_independent_of_anisotropic_cell_velocity_correction()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_oblique_triangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::Index face_id{oblique_boundary_face_id(mesh)};
    const cfd::Vector2 &area_vector{mesh.face_area_vectors()[face_id]};
    const cfd::Vector2 pressure_correction_gradient{-area_vector.y, area_vector.x};

    // This is the gradient of a linear p' satisfying homogeneous scalar
    // Neumann data on the oblique face, but an anisotropic D_P does not map it
    // to a zero normal velocity response.
    require_near(pressure_correction_gradient.x * area_vector.x + pressure_correction_gradient.y * area_vector.y, 0.0,
                 test_tolerance, "The discriminating gradient is not tangent to the oblique boundary.");

    cfd::CellVectorField gradient{mesh.cell_count(), pressure_correction_gradient};
    cfd::CellMomentumPressureResponse momentum_response{mesh.cell_count()};
    momentum_response.u()[0] = 2.0;
    momentum_response.v()[0] = 3.0;
    const double response_flux{momentum_response.u()[0] * pressure_correction_gradient.x * area_vector.x +
                               momentum_response.v()[0] * pressure_correction_gradient.y * area_vector.y};
    require(std::abs(response_flux) > 1.0,
            "Anisotropic momentum response did not distinguish the two boundary constraints.");

    const cfd::IncompressiblePressureVelocityCorrection corrector{mesh};
    cfd::CellVelocityField velocity{mesh.cell_count()};
    corrector.correct_velocity(gradient, momentum_response, velocity);
    require_near(velocity.u()[0] * area_vector.x + velocity.v()[0] * area_vector.y, -response_flux, test_tolerance,
                 "Cell velocity correction does not match the anisotropic response flux.");

    const cfd::PressureCorrectionBoundaryConditions boundary_conditions{
        make_uniform_boundary_conditions(mesh, cfd::PressureCorrectionBoundaryConditionType::FixedMassFlux)};
    const cfd::CellScalarField pressure_correction{mesh.cell_count(), 7.0};
    const cfd::FacePressureResponseField face_pressure_response{mesh.face_count(),
                                                                std::numeric_limits<double>::quiet_NaN()};
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    seed_flux(mass_flux);
    corrector.correct_face_mass_flux(pressure_correction, boundary_conditions, face_pressure_response, mass_flux);
    require_seeded_flux_unchanged(mass_flux,
                                  "FixedMassFlux face correction depended on the anisotropic corrected cell velocity.");
}

void test_velocity_correction_validation_is_transactional()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::IncompressiblePressureVelocityCorrection corrector{mesh};
    const std::array invalid_responses{
        0.0,
        -1.0,
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::infinity(),
    };

    for (const double invalid_response : invalid_responses)
    {
        cfd::CellVectorField gradient{mesh.cell_count(), {1.0, 1.0}};
        cfd::CellMomentumPressureResponse response{mesh.cell_count()};
        response.u()[0] = 1.0;
        response.v()[0] = 1.0;
        response.u()[1] = invalid_response;
        response.v()[1] = 1.0;
        cfd::CellVelocityField velocity{mesh.cell_count(), {10.0, 20.0}};
        require_throws<std::runtime_error>([&]() { corrector.correct_velocity(gradient, response, velocity); },
                                           "Velocity correction accepted an invalid momentum response.");
        require(velocity.u()[0] == 10.0 && velocity.u()[1] == 10.0 && velocity.v()[0] == 20.0 &&
                    velocity.v()[1] == 20.0,
                "Invalid momentum response partially modified velocity.");
    }

    {
        cfd::CellVectorField gradient{mesh.cell_count(), {1.0, 1.0}};
        gradient[1].y = std::numeric_limits<double>::infinity();
        cfd::CellMomentumPressureResponse response{mesh.cell_count()};
        response.u()[0] = response.v()[0] = response.u()[1] = response.v()[1] = 1.0;
        cfd::CellVelocityField velocity{mesh.cell_count(), {10.0, 20.0}};
        require_throws<std::runtime_error>([&]() { corrector.correct_velocity(gradient, response, velocity); },
                                           "Velocity correction accepted a non-finite gradient.");
        require(velocity.u()[0] == 10.0 && velocity.u()[1] == 10.0 && velocity.v()[0] == 20.0 &&
                    velocity.v()[1] == 20.0,
                "Non-finite gradient partially modified velocity.");
    }

    {
        cfd::CellVectorField gradient{mesh.cell_count(), {1.0, 1.0}};
        cfd::CellMomentumPressureResponse response{mesh.cell_count()};
        response.u()[0] = response.v()[0] = response.u()[1] = response.v()[1] = 1.0;
        cfd::CellVelocityField velocity{mesh.cell_count(), {10.0, 20.0}};
        velocity.v()[1] = std::numeric_limits<double>::quiet_NaN();
        require_throws<std::runtime_error>([&]() { corrector.correct_velocity(gradient, response, velocity); },
                                           "Velocity correction accepted a non-finite velocity.");
        require(velocity.u()[0] == 10.0 && velocity.u()[1] == 10.0 && velocity.v()[0] == 20.0 &&
                    std::isnan(velocity.v()[1]),
                "Non-finite velocity input caused partial mutation.");
    }

    cfd::CellVectorField wrong_gradient{mesh.cell_count() - 1};
    cfd::CellMomentumPressureResponse response{mesh.cell_count()};
    cfd::CellVelocityField velocity{mesh.cell_count()};
    require_throws<std::invalid_argument>([&]() { corrector.correct_velocity(wrong_gradient, response, velocity); },
                                          "Velocity correction accepted the wrong gradient cardinality.");
    cfd::CellVectorField gradient{mesh.cell_count()};
    cfd::CellMomentumPressureResponse wrong_response{mesh.cell_count() - 1};
    require_throws<std::invalid_argument>([&]() { corrector.correct_velocity(gradient, wrong_response, velocity); },
                                          "Velocity correction accepted the wrong response cardinality.");
    cfd::CellVelocityField wrong_velocity{mesh.cell_count() - 1};
    require_throws<std::invalid_argument>([&]() { corrector.correct_velocity(gradient, response, wrong_velocity); },
                                          "Velocity correction accepted the wrong velocity cardinality.");
}

void test_solved_pressure_correction_enforces_cell_continuity()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::Index internal_id{internal_face_id(mesh)};
    const cfd::Index left_face_id{boundary_face_id(mesh, left_boundary_id)};
    const cfd::Index right_face_id{boundary_face_id(mesh, right_boundary_id)};
    const cfd::PressureCorrectionBoundaryConditions boundary_conditions{
        make_right_fixed_pressure_boundary_conditions(mesh)};
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    mass_flux[left_face_id] = -1.5;
    mass_flux[right_face_id] = 3.0;
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count()};
    face_pressure_response[internal_id] = 2.0;
    face_pressure_response[right_face_id] = 3.0;
    cfd::CellScalarField provisional_imbalance{mesh.cell_count()};
    cfd::compute_cell_mass_imbalance(mesh, mass_flux, provisional_imbalance);

    const cfd::IncompressiblePressureCorrectionAssembler assembler{mesh};
    cfd::ScalarLinearSystem system{mesh};
    assembler.add_internal_face_contributions(mass_flux, face_pressure_response, system);
    assembler.add_boundary_provisional_flux_rhs(mass_flux, system);
    assembler.add_boundary_pressure_response(boundary_conditions, face_pressure_response, system);

    cfd::EigenConjugateGradientSolver solver{{1.0e-14, 20}};
    solver.compute_matrix(system);
    cfd::CellScalarField pressure_correction{mesh.cell_count()};
    const cfd::LinearSolveResult solve_result{solver.solve(system.rhs(), pressure_correction.values())};
    require(solve_result.converged, "Pressure-correction integration solve did not converge.");
    require(std::max(std::abs(pressure_correction[0]), std::abs(pressure_correction[1])) > 0.1,
            "Pressure-correction integration solution is trivial.");

    cfd::CellScalarField matrix_pressure_correction{mesh.cell_count()};
    system.apply_matrix(pressure_correction.values(), matrix_pressure_correction.values());
    const cfd::IncompressiblePressureVelocityCorrection corrector{mesh};
    corrector.correct_face_mass_flux(pressure_correction, boundary_conditions, face_pressure_response, mass_flux);
    cfd::CellScalarField corrected_imbalance{mesh.cell_count()};
    cfd::compute_cell_mass_imbalance(mesh, mass_flux, corrected_imbalance);

    double maximum_provisional_imbalance{};
    double maximum_corrected_imbalance{};
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        require_near(corrected_imbalance[cell_id], provisional_imbalance[cell_id] + matrix_pressure_correction[cell_id],
                     test_tolerance,
                     "Corrected imbalance does not satisfy provisional imbalance plus A times p-prime.");
        maximum_provisional_imbalance =
            std::max(maximum_provisional_imbalance, std::abs(provisional_imbalance[cell_id]));
        maximum_corrected_imbalance = std::max(maximum_corrected_imbalance, std::abs(corrected_imbalance[cell_id]));
    }
    require(maximum_corrected_imbalance < 1.0e-12, "Solved pressure correction did not enforce cell continuity.");
    require(maximum_corrected_imbalance < 1.0e-10 * maximum_provisional_imbalance,
            "Corrected mass imbalance was not substantially reduced.");
}

} // namespace

int main()
{
    int failure_count{};

    failure_count += cfd::test::run_test("internal face pressure correction",
                                         test_internal_face_correction_sign_and_additive_semantics);
    failure_count +=
        cfd::test::run_test("boundary face pressure correction", test_boundary_face_correction_semantics_and_isolation);
    failure_count += cfd::test::run_test("face correction cardinality validation",
                                         test_face_correction_rejects_incompatible_cardinalities_before_mutation);
    failure_count += cfd::test::run_test("face correction value validation",
                                         test_face_correction_rejects_invalid_used_values_before_mutation);
    failure_count += cfd::test::run_test("face correction transactional validation",
                                         test_face_correction_is_transactional_for_later_invalid_face);
    failure_count +=
        cfd::test::run_test("pressure correction update", test_pressure_correction_with_full_and_partial_relaxation);
    failure_count +=
        cfd::test::run_test("pressure correction validation", test_pressure_correction_validation_is_transactional);
    failure_count +=
        cfd::test::run_test("velocity correction update", test_velocity_correction_is_componentwise_and_unrelaxed);
    failure_count += cfd::test::run_test("FixedMassFlux anisotropic cell-response isolation",
                                         test_fixed_mass_flux_is_independent_of_anisotropic_cell_velocity_correction);
    failure_count +=
        cfd::test::run_test("velocity correction validation", test_velocity_correction_validation_is_transactional);
    failure_count += cfd::test::run_test("solved pressure correction continuity",
                                         test_solved_pressure_correction_enforces_cell_continuity);

    return cfd::test::finish_tests(failure_count, "incompressible pressure-velocity correction");
}
