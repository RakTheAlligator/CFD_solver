#include "cfd/numerics/RhieChowBoundaryFaceInterpolation.hpp"

#include "cfd/field/CellMomentumPressureResponse.hpp"
#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVectorField.hpp"
#include "cfd/field/CellVelocityField.hpp"
#include "cfd/field/FaceFluxField.hpp"
#include "cfd/field/FacePressureResponseField.hpp"
#include "cfd/field/PressureCorrectionBoundaryConditions.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/math/Point2.hpp"
#include "cfd/math/Vector2.hpp"
#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/mesh/Types.hpp"
#include "cfd/meshing/RawMeshData.hpp"

#include "support/TestUtils.hpp"

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

static_assert(!std::is_copy_constructible_v<cfd::RhieChowBoundaryFaceInterpolation>);
static_assert(!std::is_copy_assignable_v<cfd::RhieChowBoundaryFaceInterpolation>);
static_assert(std::is_nothrow_move_constructible_v<cfd::RhieChowBoundaryFaceInterpolation>);
static_assert(!std::is_move_assignable_v<cfd::RhieChowBoundaryFaceInterpolation>);

[[nodiscard]]
cfd::RawMeshData make_single_cell_rectangle_raw_mesh()
{
    cfd::RawMeshData raw_mesh;
    raw_mesh.nodes = {{0.0, 0.0}, {2.0, 0.0}, {2.0, 1.0}, {0.0, 1.0}};
    raw_mesh.cell_types = {cfd::CellType::Quadrilateral};
    raw_mesh.cell_nodes = {0, 1, 2, 3};
    raw_mesh.cell_node_offsets = {0, 4};
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
cfd::RawMeshData make_single_cell_sheared_raw_mesh()
{
    cfd::RawMeshData raw_mesh{make_single_cell_rectangle_raw_mesh()};
    raw_mesh.nodes = {{0.0, 0.0}, {2.0, 0.0}, {2.5, 1.0}, {0.5, 1.0}};
    return raw_mesh;
}

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
cfd::RawMeshData make_unusable_boundary_projection_raw_mesh()
{
    constexpr cfd::BoundaryId fixed_pressure_boundary_id{0};
    constexpr cfd::BoundaryId fixed_mass_flux_boundary_id{1};

    cfd::RawMeshData raw_mesh;
    raw_mesh.nodes = {{0.0, 0.0}, {1.0, 0.0}, {1.0e15, 1.0}};
    raw_mesh.cell_types = {cfd::CellType::Triangle};
    raw_mesh.cell_nodes = {0, 1, 2};
    raw_mesh.cell_node_offsets = {0, 3};
    raw_mesh.boundary_groups = {
        {fixed_pressure_boundary_id, "fixedPressure"},
        {fixed_mass_flux_boundary_id, "fixedMassFlux"},
    };
    raw_mesh.boundary_edges = {
        {{0, 1}, fixed_pressure_boundary_id},
        {{1, 2}, fixed_mass_flux_boundary_id},
        {{2, 0}, fixed_mass_flux_boundary_id},
    };
    return raw_mesh;
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
    throw std::runtime_error("Boundary Rhie-Chow fixture has no face for the requested boundary group.");
}

[[nodiscard]]
double dot(const cfd::Vector2 &first, const cfd::Vector2 &second) noexcept
{
    return first.x * second.x + first.y * second.y;
}

[[nodiscard]]
cfd::ScalarBoundaryConditions make_pressure_boundary_conditions(const cfd::Mesh &mesh,
                                                                const cfd::BoundaryId dirichlet_boundary_id,
                                                                const double pressure_value)
{
    std::vector<cfd::ScalarBoundaryCondition> conditions(mesh.boundary_groups().size(),
                                                         {cfd::ScalarBoundaryConditionType::Neumann, 0.0});
    conditions[dirichlet_boundary_id] = {cfd::ScalarBoundaryConditionType::Dirichlet, pressure_value};
    return {mesh.boundary_groups().size(), std::move(conditions)};
}

[[nodiscard]]
cfd::PressureCorrectionBoundaryConditions make_pressure_correction_boundary_conditions(
    const cfd::Mesh &mesh, const cfd::BoundaryId fixed_pressure_boundary_id)
{
    std::vector<cfd::PressureCorrectionBoundaryConditionType> conditions(
        mesh.boundary_groups().size(), cfd::PressureCorrectionBoundaryConditionType::FixedMassFlux);
    conditions[fixed_pressure_boundary_id] = cfd::PressureCorrectionBoundaryConditionType::FixedPressure;
    return {mesh.boundary_groups().size(), std::move(conditions)};
}

void initialize_valid_inputs(cfd::CellVelocityField &velocity, cfd::CellScalarField &pressure,
                             cfd::CellVectorField &pressure_gradient,
                             cfd::CellMomentumPressureResponse &momentum_response)
{
    for (cfd::Index cell_id = 0; cell_id < velocity.size(); ++cell_id)
    {
        velocity.u()[cell_id] = 1.0 + 0.25 * static_cast<double>(cell_id);
        velocity.v()[cell_id] = -0.5 + 0.2 * static_cast<double>(cell_id);
        pressure[cell_id] = 3.0 + static_cast<double>(cell_id);
        pressure_gradient[cell_id] = {0.4, -0.3};
        momentum_response.u()[cell_id] = 0.5 + 0.1 * static_cast<double>(cell_id);
        momentum_response.v()[cell_id] = 1.2 + 0.2 * static_cast<double>(cell_id);
    }
}

void seed_outputs(cfd::FaceFluxField &mass_flux, cfd::FacePressureResponseField &face_pressure_response)
{
    for (cfd::Index face_id = 0; face_id < mass_flux.size(); ++face_id)
    {
        mass_flux[face_id] = 101.0 + static_cast<double>(face_id);
    }
    for (cfd::Index face_id = 0; face_id < face_pressure_response.size(); ++face_id)
    {
        face_pressure_response[face_id] = -201.0 - static_cast<double>(face_id);
    }
}

void require_seeded_outputs_unchanged(const cfd::FaceFluxField &mass_flux,
                                      const cfd::FacePressureResponseField &face_pressure_response,
                                      const std::string &context)
{
    for (cfd::Index face_id = 0; face_id < mass_flux.size(); ++face_id)
    {
        require(mass_flux[face_id] == 101.0 + static_cast<double>(face_id), context + " mass flux changed.");
    }
    for (cfd::Index face_id = 0; face_id < face_pressure_response.size(); ++face_id)
    {
        require(face_pressure_response[face_id] == -201.0 - static_cast<double>(face_id),
                context + " pressure response changed.");
    }
}

template <typename Exception, typename Function>
void require_rejected_without_output_mutation(Function &&function, const cfd::FaceFluxField &mass_flux,
                                              const cfd::FacePressureResponseField &face_pressure_response,
                                              const std::string &message)
{
    require_throws<Exception>(std::forward<Function>(function), message);
    require_seeded_outputs_unchanged(mass_flux, face_pressure_response, message + " Outputs were partially modified.");
}

void test_constructor_rejects_invalid_density()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};

    for (const double density : {0.0, -1.0, std::numeric_limits<double>::quiet_NaN(),
                                 std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()})
    {
        require_throws<std::invalid_argument>(
            [&mesh, density]() { const cfd::RhieChowBoundaryFaceInterpolation interpolation{mesh, density}; },
            "Boundary Rhie-Chow interpolation accepted an invalid density.");
    }
}

void test_constant_pressure_and_orthogonal_response()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    constexpr double density{1.7};
    const cfd::RhieChowBoundaryFaceInterpolation interpolation{mesh, density};
    cfd::CellVelocityField velocity{mesh.cell_count(), {1.25, -0.75}};
    const cfd::CellScalarField pressure{mesh.cell_count(), 12.0};
    const cfd::CellVectorField pressure_gradient{mesh.cell_count()};
    cfd::CellMomentumPressureResponse momentum_response{mesh.cell_count()};
    momentum_response.u()[0] = 0.8;
    momentum_response.v()[0] = 1.4;
    const cfd::ScalarBoundaryConditions pressure_boundary_conditions{
        make_pressure_boundary_conditions(mesh, right_boundary_id, 12.0)};
    const cfd::PressureCorrectionBoundaryConditions pressure_correction_boundary_conditions{
        make_pressure_correction_boundary_conditions(mesh, right_boundary_id)};
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count()};
    seed_outputs(mass_flux, face_pressure_response);
    const cfd::Index face_id{boundary_face_id(mesh, right_boundary_id)};

    interpolation.update_fixed_pressure_boundaries(
        velocity, pressure, pressure_gradient, momentum_response, pressure_boundary_conditions,
        pressure_correction_boundary_conditions, mass_flux, face_pressure_response);

    require_near(mass_flux[face_id], density * velocity.u()[0], test_tolerance,
                 "Constant pressure did not reduce to the owner normal velocity flux.");
    require_near(face_pressure_response[face_id], density * momentum_response.u()[0], test_tolerance,
                 "Orthogonal boundary interpolation produced an incorrect pressure response.");
}

void test_nonorthogonal_anisotropic_response()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_cell_sheared_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    constexpr double density{2.0};
    const cfd::RhieChowBoundaryFaceInterpolation interpolation{mesh, density};
    cfd::CellVelocityField velocity{mesh.cell_count(), {1.2, -0.4}};
    cfd::CellScalarField pressure{mesh.cell_count(), 5.0};
    cfd::CellVectorField pressure_gradient{mesh.cell_count(), {0.3, -0.2}};
    cfd::CellMomentumPressureResponse momentum_response{mesh.cell_count()};
    momentum_response.u()[0] = 0.4;
    momentum_response.v()[0] = 1.6;
    const cfd::ScalarBoundaryConditions pressure_boundary_conditions{
        make_pressure_boundary_conditions(mesh, right_boundary_id, 7.0)};
    const cfd::PressureCorrectionBoundaryConditions pressure_correction_boundary_conditions{
        make_pressure_correction_boundary_conditions(mesh, right_boundary_id)};
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count()};
    const cfd::Index face_id{boundary_face_id(mesh, right_boundary_id)};
    const cfd::Vector2 &area_vector{mesh.face_area_vectors()[face_id]};
    const cfd::Point2 &owner_center{mesh.cell_centers()[0]};
    const cfd::Point2 &face_center{mesh.face_centers()[face_id]};
    const cfd::Vector2 displacement{face_center.x - owner_center.x, face_center.y - owner_center.y};
    const cfd::Vector2 response_area_vector{momentum_response.u()[0] * area_vector.x,
                                            momentum_response.v()[0] * area_vector.y};
    const double pressure_coefficient{dot(area_vector, response_area_vector) / dot(area_vector, displacement)};
    const cfd::Vector2 tangential_response{response_area_vector.x - pressure_coefficient * displacement.x,
                                           response_area_vector.y - pressure_coefficient * displacement.y};
    const cfd::Vector2 pressure_free_velocity{
        velocity.u()[0] + momentum_response.u()[0] * pressure_gradient[0].x,
        velocity.v()[0] + momentum_response.v()[0] * pressure_gradient[0].y,
    };
    const double expected_flux{density *
                               (dot(pressure_free_velocity, area_vector) - pressure_coefficient * (7.0 - pressure[0]) -
                                dot(pressure_gradient[0], tangential_response))};

    interpolation.update_fixed_pressure_boundaries(
        velocity, pressure, pressure_gradient, momentum_response, pressure_boundary_conditions,
        pressure_correction_boundary_conditions, mass_flux, face_pressure_response);

    require(area_vector.x != 0.0 && area_vector.y != 0.0,
            "Non-orthogonal boundary fixture lacks both area-vector components.");
    require(momentum_response.u()[0] != momentum_response.v()[0],
            "Anisotropic boundary fixture has equal momentum responses.");
    require_near(dot(tangential_response, area_vector), 0.0, test_tolerance,
                 "Boundary tangential response is not orthogonal to the area vector.");
    require_near(face_pressure_response[face_id], density * pressure_coefficient, test_tolerance,
                 "Non-orthogonal anisotropic pressure response is incorrect.");
    require_near(mass_flux[face_id], expected_flux, test_tolerance,
                 "Non-orthogonal anisotropic boundary mass flux is incorrect.");
}

void test_linear_pressure_consistency_and_shift_invariance()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_cell_sheared_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    constexpr double density{1.4};
    constexpr double pressure_x_coefficient{2.25};
    constexpr double pressure_y_coefficient{-1.75};
    constexpr double pressure_constant{6.5};
    constexpr double pressure_shift{128.0};
    const cfd::RhieChowBoundaryFaceInterpolation interpolation{mesh, density};
    cfd::CellVelocityField velocity{mesh.cell_count(), {-0.8, 1.7}};
    cfd::CellScalarField pressure{mesh.cell_count()};
    const cfd::CellVectorField pressure_gradient{mesh.cell_count(), {pressure_x_coefficient, pressure_y_coefficient}};
    cfd::CellMomentumPressureResponse momentum_response{mesh.cell_count()};
    momentum_response.u()[0] = 0.3;
    momentum_response.v()[0] = 1.4;
    const cfd::Index face_id{boundary_face_id(mesh, right_boundary_id)};
    const cfd::Point2 &cell_center{mesh.cell_centers()[0]};
    const cfd::Point2 &face_center{mesh.face_centers()[face_id]};
    pressure[0] = pressure_x_coefficient * cell_center.x + pressure_y_coefficient * cell_center.y + pressure_constant;
    const double boundary_pressure{pressure_x_coefficient * face_center.x + pressure_y_coefficient * face_center.y +
                                   pressure_constant};
    const cfd::PressureCorrectionBoundaryConditions pressure_correction_boundary_conditions{
        make_pressure_correction_boundary_conditions(mesh, right_boundary_id)};
    const cfd::ScalarBoundaryConditions original_pressure_boundary_conditions{
        make_pressure_boundary_conditions(mesh, right_boundary_id, boundary_pressure)};
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count()};

    interpolation.update_fixed_pressure_boundaries(
        velocity, pressure, pressure_gradient, momentum_response, original_pressure_boundary_conditions,
        pressure_correction_boundary_conditions, mass_flux, face_pressure_response);
    const double original_flux{mass_flux[face_id]};
    const double original_response{face_pressure_response[face_id]};
    require_near(original_flux, density * dot({velocity.u()[0], velocity.v()[0]}, mesh.face_area_vectors()[face_id]),
                 test_tolerance, "A linear pressure field disturbed the owner normal velocity flux.");

    pressure[0] += pressure_shift;
    const cfd::ScalarBoundaryConditions shifted_pressure_boundary_conditions{
        make_pressure_boundary_conditions(mesh, right_boundary_id, boundary_pressure + pressure_shift)};
    interpolation.update_fixed_pressure_boundaries(
        velocity, pressure, pressure_gradient, momentum_response, shifted_pressure_boundary_conditions,
        pressure_correction_boundary_conditions, mass_flux, face_pressure_response);

    require_near(mass_flux[face_id], original_flux, test_tolerance,
                 "Adding a constant pressure changed the boundary Rhie-Chow mass flux.");
    require_near(face_pressure_response[face_id], original_response, 0.0,
                 "Adding a constant pressure changed the boundary pressure response.");
}

void test_only_fixed_pressure_outputs_are_modified()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::RhieChowBoundaryFaceInterpolation interpolation{mesh, 1.0};
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    cfd::CellVectorField pressure_gradient{mesh.cell_count()};
    cfd::CellMomentumPressureResponse momentum_response{mesh.cell_count()};
    initialize_valid_inputs(velocity, pressure, pressure_gradient, momentum_response);
    const cfd::ScalarBoundaryConditions pressure_boundary_conditions{
        make_pressure_boundary_conditions(mesh, right_boundary_id, 4.0)};
    const cfd::PressureCorrectionBoundaryConditions pressure_correction_boundary_conditions{
        make_pressure_correction_boundary_conditions(mesh, right_boundary_id)};
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count()};
    seed_outputs(mass_flux, face_pressure_response);
    const cfd::FaceFluxField original_mass_flux{mass_flux};
    const cfd::FacePressureResponseField original_face_pressure_response{face_pressure_response};
    const cfd::Index fixed_pressure_face_id{boundary_face_id(mesh, right_boundary_id)};

    interpolation.update_fixed_pressure_boundaries(
        velocity, pressure, pressure_gradient, momentum_response, pressure_boundary_conditions,
        pressure_correction_boundary_conditions, mass_flux, face_pressure_response);

    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (face_id == fixed_pressure_face_id)
        {
            require(mass_flux[face_id] != original_mass_flux[face_id] &&
                        face_pressure_response[face_id] != original_face_pressure_response[face_id],
                    "FixedPressure boundary outputs were not overwritten.");
            continue;
        }
        require(mass_flux[face_id] == original_mass_flux[face_id],
                "Boundary interpolation modified an internal or FixedMassFlux mass flux.");
        require(face_pressure_response[face_id] == original_face_pressure_response[face_id],
                "Boundary interpolation modified an internal or FixedMassFlux pressure response.");
    }
}

void test_fixed_pressure_requires_physical_dirichlet_pressure()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::RhieChowBoundaryFaceInterpolation interpolation{mesh, 1.0};
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    cfd::CellVectorField pressure_gradient{mesh.cell_count()};
    cfd::CellMomentumPressureResponse momentum_response{mesh.cell_count()};
    initialize_valid_inputs(velocity, pressure, pressure_gradient, momentum_response);
    const cfd::ScalarBoundaryConditions pressure_boundary_conditions{
        mesh.boundary_groups().size(),
        std::vector<cfd::ScalarBoundaryCondition>(mesh.boundary_groups().size(),
                                                  {cfd::ScalarBoundaryConditionType::Neumann, 0.0})};
    const cfd::PressureCorrectionBoundaryConditions pressure_correction_boundary_conditions{
        make_pressure_correction_boundary_conditions(mesh, right_boundary_id)};
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count()};
    seed_outputs(mass_flux, face_pressure_response);

    require_rejected_without_output_mutation<std::invalid_argument>(
        [&]() {
            interpolation.update_fixed_pressure_boundaries(
                velocity, pressure, pressure_gradient, momentum_response, pressure_boundary_conditions,
                pressure_correction_boundary_conditions, mass_flux, face_pressure_response);
        },
        mass_flux, face_pressure_response, "FixedPressure pressure correction accepted a physical Neumann pressure.");
}

void test_rejects_incompatible_cardinalities_before_mutation()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::RhieChowBoundaryFaceInterpolation interpolation{mesh, 1.0};
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    cfd::CellVectorField pressure_gradient{mesh.cell_count()};
    cfd::CellMomentumPressureResponse momentum_response{mesh.cell_count()};
    initialize_valid_inputs(velocity, pressure, pressure_gradient, momentum_response);
    const cfd::ScalarBoundaryConditions pressure_boundary_conditions{
        make_pressure_boundary_conditions(mesh, right_boundary_id, 4.0)};
    const cfd::PressureCorrectionBoundaryConditions pressure_correction_boundary_conditions{
        make_pressure_correction_boundary_conditions(mesh, right_boundary_id)};
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count()};
    seed_outputs(mass_flux, face_pressure_response);

    const auto require_rejected = [&](const auto &input_velocity, const auto &input_pressure,
                                      const auto &input_gradient, const auto &input_response,
                                      const auto &physical_conditions, const auto &correction_conditions,
                                      auto &output_flux, auto &output_response, const std::string &message) {
        require_rejected_without_output_mutation<std::invalid_argument>(
            [&]() {
                interpolation.update_fixed_pressure_boundaries(input_velocity, input_pressure, input_gradient,
                                                               input_response, physical_conditions,
                                                               correction_conditions, output_flux, output_response);
            },
            output_flux, output_response, message);
    };

    const cfd::CellVelocityField wrong_velocity{mesh.cell_count() + 1};
    require_rejected(wrong_velocity, pressure, pressure_gradient, momentum_response, pressure_boundary_conditions,
                     pressure_correction_boundary_conditions, mass_flux, face_pressure_response,
                     "Boundary interpolation accepted an incorrectly sized velocity.");
    const cfd::CellScalarField wrong_pressure{mesh.cell_count() + 1};
    require_rejected(velocity, wrong_pressure, pressure_gradient, momentum_response, pressure_boundary_conditions,
                     pressure_correction_boundary_conditions, mass_flux, face_pressure_response,
                     "Boundary interpolation accepted an incorrectly sized pressure.");
    const cfd::CellVectorField wrong_gradient{mesh.cell_count() + 1};
    require_rejected(velocity, pressure, wrong_gradient, momentum_response, pressure_boundary_conditions,
                     pressure_correction_boundary_conditions, mass_flux, face_pressure_response,
                     "Boundary interpolation accepted an incorrectly sized pressure gradient.");
    const cfd::CellMomentumPressureResponse wrong_momentum_response{mesh.cell_count() + 1};
    require_rejected(velocity, pressure, pressure_gradient, wrong_momentum_response, pressure_boundary_conditions,
                     pressure_correction_boundary_conditions, mass_flux, face_pressure_response,
                     "Boundary interpolation accepted an incorrectly sized momentum response.");

    const cfd::ScalarBoundaryConditions wrong_pressure_conditions{
        mesh.boundary_groups().size() + 1,
        std::vector<cfd::ScalarBoundaryCondition>(mesh.boundary_groups().size() + 1,
                                                  {cfd::ScalarBoundaryConditionType::Dirichlet, 4.0})};
    require_rejected(velocity, pressure, pressure_gradient, momentum_response, wrong_pressure_conditions,
                     pressure_correction_boundary_conditions, mass_flux, face_pressure_response,
                     "Boundary interpolation accepted an incorrect physical pressure-condition count.");
    const cfd::PressureCorrectionBoundaryConditions wrong_correction_conditions{
        mesh.boundary_groups().size() + 1,
        std::vector<cfd::PressureCorrectionBoundaryConditionType>(
            mesh.boundary_groups().size() + 1, cfd::PressureCorrectionBoundaryConditionType::FixedPressure)};
    require_rejected(velocity, pressure, pressure_gradient, momentum_response, pressure_boundary_conditions,
                     wrong_correction_conditions, mass_flux, face_pressure_response,
                     "Boundary interpolation accepted an incorrect pressure-correction condition count.");

    cfd::FaceFluxField wrong_mass_flux{mesh.face_count() + 1};
    seed_outputs(wrong_mass_flux, face_pressure_response);
    require_rejected(velocity, pressure, pressure_gradient, momentum_response, pressure_boundary_conditions,
                     pressure_correction_boundary_conditions, wrong_mass_flux, face_pressure_response,
                     "Boundary interpolation accepted an incorrectly sized mass flux.");
    cfd::FacePressureResponseField wrong_face_response{mesh.face_count() + 1};
    seed_outputs(mass_flux, wrong_face_response);
    require_rejected(velocity, pressure, pressure_gradient, momentum_response, pressure_boundary_conditions,
                     pressure_correction_boundary_conditions, mass_flux, wrong_face_response,
                     "Boundary interpolation accepted an incorrectly sized face pressure response.");
}

void test_rejects_nonfinite_used_inputs_before_mutation()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::RhieChowBoundaryFaceInterpolation interpolation{mesh, 1.0};
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    cfd::CellVectorField pressure_gradient{mesh.cell_count()};
    cfd::CellMomentumPressureResponse momentum_response{mesh.cell_count()};
    initialize_valid_inputs(velocity, pressure, pressure_gradient, momentum_response);
    const cfd::ScalarBoundaryConditions pressure_boundary_conditions{
        make_pressure_boundary_conditions(mesh, right_boundary_id, 4.0)};
    const cfd::PressureCorrectionBoundaryConditions pressure_correction_boundary_conditions{
        make_pressure_correction_boundary_conditions(mesh, right_boundary_id)};
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count()};
    seed_outputs(mass_flux, face_pressure_response);
    const auto require_rejected = [&](const std::string &message) {
        require_rejected_without_output_mutation<std::runtime_error>(
            [&]() {
                interpolation.update_fixed_pressure_boundaries(
                    velocity, pressure, pressure_gradient, momentum_response, pressure_boundary_conditions,
                    pressure_correction_boundary_conditions, mass_flux, face_pressure_response);
            },
            mass_flux, face_pressure_response, message);
    };

    velocity.u()[0] = std::numeric_limits<double>::quiet_NaN();
    require_rejected("Boundary interpolation accepted a NaN u velocity.");
    velocity.u()[0] = 1.0;
    velocity.v()[0] = std::numeric_limits<double>::infinity();
    require_rejected("Boundary interpolation accepted an infinite v velocity.");
    velocity.v()[0] = -0.5;
    pressure[0] = std::numeric_limits<double>::quiet_NaN();
    require_rejected("Boundary interpolation accepted a NaN pressure.");
    pressure[0] = 3.0;
    pressure_gradient[0].x = std::numeric_limits<double>::infinity();
    require_rejected("Boundary interpolation accepted an infinite pressure-gradient x component.");
    pressure_gradient[0].x = 0.4;
    pressure_gradient[0].y = -std::numeric_limits<double>::infinity();
    require_rejected("Boundary interpolation accepted an infinite pressure-gradient y component.");
}

void test_rejects_invalid_momentum_responses_and_overflow_before_mutation()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_single_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::RhieChowBoundaryFaceInterpolation interpolation{mesh, 1.0};
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    cfd::CellVectorField pressure_gradient{mesh.cell_count()};
    cfd::CellMomentumPressureResponse momentum_response{mesh.cell_count()};
    initialize_valid_inputs(velocity, pressure, pressure_gradient, momentum_response);
    const cfd::ScalarBoundaryConditions pressure_boundary_conditions{
        make_pressure_boundary_conditions(mesh, right_boundary_id, 4.0)};
    const cfd::PressureCorrectionBoundaryConditions pressure_correction_boundary_conditions{
        make_pressure_correction_boundary_conditions(mesh, right_boundary_id)};
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count()};
    seed_outputs(mass_flux, face_pressure_response);
    const auto require_rejected = [&](const std::string &message) {
        require_rejected_without_output_mutation<std::runtime_error>(
            [&]() {
                interpolation.update_fixed_pressure_boundaries(
                    velocity, pressure, pressure_gradient, momentum_response, pressure_boundary_conditions,
                    pressure_correction_boundary_conditions, mass_flux, face_pressure_response);
            },
            mass_flux, face_pressure_response, message);
    };

    for (const double invalid_response :
         {0.0, -1.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
          -std::numeric_limits<double>::infinity()})
    {
        momentum_response.u()[0] = invalid_response;
        require_rejected("Boundary interpolation accepted an invalid u-momentum response.");
    }
    momentum_response.u()[0] = 0.5;
    momentum_response.v()[0] = 0.0;
    require_rejected("Boundary interpolation accepted an invalid v-momentum response.");

    momentum_response.u()[0] = std::numeric_limits<double>::max();
    momentum_response.v()[0] = std::numeric_limits<double>::max();
    const cfd::RhieChowBoundaryFaceInterpolation high_density_interpolation{mesh, 2.0};
    require_rejected_without_output_mutation<std::runtime_error>(
        [&]() {
            high_density_interpolation.update_fixed_pressure_boundaries(
                velocity, pressure, pressure_gradient, momentum_response, pressure_boundary_conditions,
                pressure_correction_boundary_conditions, mass_flux, face_pressure_response);
        },
        mass_flux, face_pressure_response, "Boundary interpolation accepted an overflowing pressure response.");
}

void test_rejects_unusable_boundary_geometry()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_unusable_boundary_projection_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::RhieChowBoundaryFaceInterpolation interpolation{mesh, 1.0};
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    cfd::CellVectorField pressure_gradient{mesh.cell_count()};
    cfd::CellMomentumPressureResponse momentum_response{mesh.cell_count()};
    initialize_valid_inputs(velocity, pressure, pressure_gradient, momentum_response);
    const cfd::ScalarBoundaryConditions pressure_boundary_conditions{make_pressure_boundary_conditions(mesh, 0, 4.0)};
    const cfd::PressureCorrectionBoundaryConditions pressure_correction_boundary_conditions{
        make_pressure_correction_boundary_conditions(mesh, 0)};
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count()};
    seed_outputs(mass_flux, face_pressure_response);

    require_rejected_without_output_mutation<std::runtime_error>(
        [&]() {
            interpolation.update_fixed_pressure_boundaries(
                velocity, pressure, pressure_gradient, momentum_response, pressure_boundary_conditions,
                pressure_correction_boundary_conditions, mass_flux, face_pressure_response);
        },
        mass_flux, face_pressure_response, "Boundary interpolation accepted unusable projection geometry.");
}

void test_later_invalid_fixed_pressure_boundary_is_transactional()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::RhieChowBoundaryFaceInterpolation interpolation{mesh, 1.0};
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    cfd::CellVectorField pressure_gradient{mesh.cell_count()};
    cfd::CellMomentumPressureResponse momentum_response{mesh.cell_count()};
    initialize_valid_inputs(velocity, pressure, pressure_gradient, momentum_response);
    std::vector<cfd::ScalarBoundaryCondition> physical_conditions(mesh.boundary_groups().size(),
                                                                  {cfd::ScalarBoundaryConditionType::Neumann, 0.0});
    physical_conditions[left_boundary_id] = {cfd::ScalarBoundaryConditionType::Dirichlet, 3.0};
    physical_conditions[right_boundary_id] = {cfd::ScalarBoundaryConditionType::Dirichlet, 4.0};
    const cfd::ScalarBoundaryConditions pressure_boundary_conditions{mesh.boundary_groups().size(),
                                                                     std::move(physical_conditions)};
    std::vector<cfd::PressureCorrectionBoundaryConditionType> correction_conditions(
        mesh.boundary_groups().size(), cfd::PressureCorrectionBoundaryConditionType::FixedMassFlux);
    correction_conditions[left_boundary_id] = cfd::PressureCorrectionBoundaryConditionType::FixedPressure;
    correction_conditions[right_boundary_id] = cfd::PressureCorrectionBoundaryConditionType::FixedPressure;
    const cfd::PressureCorrectionBoundaryConditions pressure_correction_boundary_conditions{
        mesh.boundary_groups().size(), std::move(correction_conditions)};
    require(boundary_face_id(mesh, left_boundary_id) < boundary_face_id(mesh, right_boundary_id),
            "Transactional fixture does not visit the valid FixedPressure face first.");
    momentum_response.u()[1] = std::numeric_limits<double>::infinity();
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count()};
    seed_outputs(mass_flux, face_pressure_response);

    require_rejected_without_output_mutation<std::runtime_error>(
        [&]() {
            interpolation.update_fixed_pressure_boundaries(
                velocity, pressure, pressure_gradient, momentum_response, pressure_boundary_conditions,
                pressure_correction_boundary_conditions, mass_flux, face_pressure_response);
        },
        mass_flux, face_pressure_response,
        "Boundary interpolation modified outputs before rejecting a later FixedPressure face.");
}

} // namespace

int main()
{
    int failure_count{};

    failure_count +=
        cfd::test::run_test("boundary Rhie-Chow density validation", test_constructor_rejects_invalid_density);
    failure_count += cfd::test::run_test("boundary Rhie-Chow constant pressure and orthogonal response",
                                         test_constant_pressure_and_orthogonal_response);
    failure_count += cfd::test::run_test("boundary Rhie-Chow non-orthogonal anisotropic response",
                                         test_nonorthogonal_anisotropic_response);
    failure_count += cfd::test::run_test("boundary Rhie-Chow linear consistency and pressure-shift invariance",
                                         test_linear_pressure_consistency_and_shift_invariance);
    failure_count +=
        cfd::test::run_test("boundary Rhie-Chow output isolation", test_only_fixed_pressure_outputs_are_modified);
    failure_count += cfd::test::run_test("boundary Rhie-Chow physical pressure semantics",
                                         test_fixed_pressure_requires_physical_dirichlet_pressure);
    failure_count += cfd::test::run_test("boundary Rhie-Chow cardinality validation",
                                         test_rejects_incompatible_cardinalities_before_mutation);
    failure_count += cfd::test::run_test("boundary Rhie-Chow finite-input validation",
                                         test_rejects_nonfinite_used_inputs_before_mutation);
    failure_count += cfd::test::run_test("boundary Rhie-Chow momentum-response validation",
                                         test_rejects_invalid_momentum_responses_and_overflow_before_mutation);
    failure_count +=
        cfd::test::run_test("boundary Rhie-Chow geometry validation", test_rejects_unusable_boundary_geometry);
    failure_count += cfd::test::run_test("boundary Rhie-Chow transactional validation",
                                         test_later_invalid_fixed_pressure_boundary_is_transactional);

    return cfd::test::finish_tests(failure_count, "Rhie-Chow boundary-face interpolation");
}
