#include "cfd/numerics/RhieChowInternalFaceInterpolation.hpp"

#include "cfd/field/CellMomentumPressureResponse.hpp"
#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVectorField.hpp"
#include "cfd/field/CellVelocityField.hpp"
#include "cfd/field/FaceFluxField.hpp"
#include "cfd/field/FacePressureResponseField.hpp"
#include "cfd/math/Point2.hpp"
#include "cfd/math/Vector2.hpp"
#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/meshing/RawMeshData.hpp"

#include "support/TestUtils.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

namespace
{

using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws;
using cfd::test::test_tolerance;

static_assert(!std::is_copy_constructible_v<cfd::RhieChowInternalFaceInterpolation>);
static_assert(!std::is_copy_assignable_v<cfd::RhieChowInternalFaceInterpolation>);
static_assert(std::is_nothrow_move_constructible_v<cfd::RhieChowInternalFaceInterpolation>);
static_assert(!std::is_move_assignable_v<cfd::RhieChowInternalFaceInterpolation>);

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
    throw std::runtime_error("Rhie-Chow test fixture has no internal face.");
}

[[nodiscard]]
double dot(const cfd::Vector2 &first, const cfd::Vector2 &second) noexcept
{
    return first.x * second.x + first.y * second.y;
}

[[nodiscard]]
double interpolation_weight(const cfd::Mesh &mesh, const cfd::Index face_id) noexcept
{
    const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};
    const cfd::Point2 &owner_center{mesh.cell_centers()[adjacency.owner]};
    const cfd::Point2 &neighbor_center{mesh.cell_centers()[adjacency.neighbor]};
    const cfd::Point2 &face_center{mesh.face_centers()[face_id]};
    const cfd::Vector2 area_vector{mesh.face_area_vectors()[face_id]};
    const cfd::Vector2 center_displacement{
        neighbor_center.x - owner_center.x,
        neighbor_center.y - owner_center.y,
    };
    const cfd::Vector2 owner_to_face{
        face_center.x - owner_center.x,
        face_center.y - owner_center.y,
    };
    return dot(area_vector, owner_to_face) / dot(area_vector, center_displacement);
}

void initialize_valid_inputs(cfd::CellVelocityField &velocity, cfd::CellScalarField &pressure,
                             cfd::CellVectorField &pressure_gradient,
                             cfd::CellMomentumPressureResponse &momentum_response)
{
    for (cfd::Index cell_id = 0; cell_id < velocity.size(); ++cell_id)
    {
        velocity.u()[cell_id] = 1.0 + static_cast<double>(cell_id);
        velocity.v()[cell_id] = -2.0 + 0.5 * static_cast<double>(cell_id);
        pressure[cell_id] = 3.0 - static_cast<double>(cell_id);
        pressure_gradient[cell_id] = {0.25, -0.75};
        momentum_response.u()[cell_id] = 0.5 + 0.25 * static_cast<double>(cell_id);
        momentum_response.v()[cell_id] = 1.0 + 0.5 * static_cast<double>(cell_id);
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
        require(mass_flux[face_id] == 101.0 + static_cast<double>(face_id), context + " (mass flux changed). ");
    }
    for (cfd::Index face_id = 0; face_id < face_pressure_response.size(); ++face_id)
    {
        require(face_pressure_response[face_id] == -201.0 - static_cast<double>(face_id),
                context + " (face pressure response changed). ");
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
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};

    for (const double density : {0.0, -1.0, std::numeric_limits<double>::quiet_NaN(),
                                 std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()})
    {
        require_throws<std::invalid_argument>(
            [&mesh, density]() { const cfd::RhieChowInternalFaceInterpolation interpolation{mesh, density}; },
            "Rhie-Chow interpolation accepted an invalid density.");
    }
}

void test_constant_pressure_reduces_to_velocity_interpolation_and_preserves_boundaries()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    constexpr double density{1.7};
    const cfd::RhieChowInternalFaceInterpolation interpolation{mesh, density};
    cfd::CellVelocityField velocity{mesh.cell_count()};
    velocity.u()[0] = 1.25;
    velocity.v()[0] = -0.75;
    velocity.u()[1] = 3.5;
    velocity.v()[1] = 2.25;
    const cfd::CellScalarField pressure{mesh.cell_count(), 12.0};
    const cfd::CellVectorField pressure_gradient{mesh.cell_count()};
    cfd::CellMomentumPressureResponse momentum_response{mesh.cell_count()};
    momentum_response.u()[0] = 0.4;
    momentum_response.v()[0] = 1.1;
    momentum_response.u()[1] = 0.9;
    momentum_response.v()[1] = 1.8;
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count()};
    seed_outputs(mass_flux, face_pressure_response);
    const cfd::FaceFluxField original_mass_flux{mass_flux};
    const cfd::FacePressureResponseField original_face_pressure_response{face_pressure_response};

    interpolation.update_internal_faces(velocity, pressure, pressure_gradient, momentum_response, mass_flux,
                                        face_pressure_response);

    const cfd::Index face_id{internal_face_id(mesh)};
    const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};
    const double lambda{interpolation_weight(mesh, face_id)};
    const cfd::Vector2 interpolated_velocity{
        (1.0 - lambda) * velocity.u()[adjacency.owner] + lambda * velocity.u()[adjacency.neighbor],
        (1.0 - lambda) * velocity.v()[adjacency.owner] + lambda * velocity.v()[adjacency.neighbor],
    };
    require_near(mass_flux[face_id], density * dot(interpolated_velocity, mesh.face_area_vectors()[face_id]),
                 test_tolerance, "Constant pressure did not reduce to ordinary velocity interpolation.");
    require(std::isfinite(face_pressure_response[face_id]) && face_pressure_response[face_id] > 0.0,
            "Constant-pressure interpolation produced an invalid face pressure response.");

    for (cfd::Index boundary_face_id = 0; boundary_face_id < mesh.face_count(); ++boundary_face_id)
    {
        if (!mesh.face_adjacencies()[boundary_face_id].is_boundary())
        {
            continue;
        }
        require(mass_flux[boundary_face_id] == original_mass_flux[boundary_face_id],
                "Internal Rhie-Chow interpolation changed a boundary mass flux.");
        require(face_pressure_response[boundary_face_id] == original_face_pressure_response[boundary_face_id],
                "Internal Rhie-Chow interpolation changed a boundary pressure response.");
    }
}

void test_additive_pressure_shift_invariance()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_sheared_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::RhieChowInternalFaceInterpolation interpolation{mesh, 2.3};
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    cfd::CellVectorField pressure_gradient{mesh.cell_count()};
    cfd::CellMomentumPressureResponse momentum_response{mesh.cell_count()};
    initialize_valid_inputs(velocity, pressure, pressure_gradient, momentum_response);
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count()};

    interpolation.update_internal_faces(velocity, pressure, pressure_gradient, momentum_response, mass_flux,
                                        face_pressure_response);
    const cfd::Index face_id{internal_face_id(mesh)};
    const double original_flux{mass_flux[face_id]};
    const double original_response{face_pressure_response[face_id]};

    constexpr double pressure_shift{128.0};
    for (double &value : pressure.values())
    {
        value += pressure_shift;
    }
    interpolation.update_internal_faces(velocity, pressure, pressure_gradient, momentum_response, mass_flux,
                                        face_pressure_response);

    require_near(mass_flux[face_id], original_flux, test_tolerance,
                 "Adding a constant pressure changed the Rhie-Chow mass flux.");
    require_near(face_pressure_response[face_id], original_response, 0.0,
                 "Adding a constant pressure changed the face pressure response.");
}

void test_linear_pressure_is_exact_on_sheared_anisotropic_face()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_sheared_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    constexpr double density{1.4};
    constexpr double pressure_x_coefficient{2.25};
    constexpr double pressure_y_coefficient{-1.75};
    constexpr double pressure_constant{6.5};
    const cfd::RhieChowInternalFaceInterpolation interpolation{mesh, density};
    cfd::CellVelocityField velocity{mesh.cell_count()};
    velocity.u()[0] = -0.8;
    velocity.v()[0] = 1.7;
    velocity.u()[1] = 2.4;
    velocity.v()[1] = -1.1;
    cfd::CellScalarField pressure{mesh.cell_count()};
    cfd::CellVectorField pressure_gradient{mesh.cell_count(), {pressure_x_coefficient, pressure_y_coefficient}};
    cfd::CellMomentumPressureResponse momentum_response{mesh.cell_count()};
    momentum_response.u()[0] = 0.3;
    momentum_response.v()[0] = 1.4;
    momentum_response.u()[1] = 1.1;
    momentum_response.v()[1] = 0.6;
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        const cfd::Point2 &center{mesh.cell_centers()[cell_id]};
        pressure[cell_id] = pressure_x_coefficient * center.x + pressure_y_coefficient * center.y + pressure_constant;
    }
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count()};

    interpolation.update_internal_faces(velocity, pressure, pressure_gradient, momentum_response, mass_flux,
                                        face_pressure_response);

    const cfd::Index face_id{internal_face_id(mesh)};
    const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};
    const double lambda{interpolation_weight(mesh, face_id)};
    const cfd::Vector2 interpolated_velocity{
        (1.0 - lambda) * velocity.u()[adjacency.owner] + lambda * velocity.u()[adjacency.neighbor],
        (1.0 - lambda) * velocity.v()[adjacency.owner] + lambda * velocity.v()[adjacency.neighbor],
    };
    require_near(mass_flux[face_id], density * dot(interpolated_velocity, mesh.face_area_vectors()[face_id]),
                 test_tolerance,
                 "An exact linear pressure disturbed velocity interpolation on an anisotropic sheared face.");
}

void test_checkerboard_pressure_jump_sensitivity_and_owner_oriented_sign()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::RhieChowInternalFaceInterpolation interpolation{mesh, 1.0};
    const cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    const cfd::CellVectorField pressure_gradient{mesh.cell_count()};
    cfd::CellMomentumPressureResponse momentum_response{mesh.cell_count()};
    momentum_response.u()[0] = 0.8;
    momentum_response.v()[0] = 1.3;
    momentum_response.u()[1] = 1.2;
    momentum_response.v()[1] = 0.7;
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count()};
    const cfd::Index face_id{internal_face_id(mesh)};
    const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};

    // This alternating pair represents the local pressure jump hidden by a
    // zero supplied gradient; it tests sensitivity, not global mode removal.
    pressure[adjacency.owner] = -2.0;
    pressure[adjacency.neighbor] = 3.0;
    interpolation.update_internal_faces(velocity, pressure, pressure_gradient, momentum_response, mass_flux,
                                        face_pressure_response);
    require(mass_flux[face_id] < 0.0,
            "An owner-to-neighbor pressure increase did not produce a negative owner-oriented flux.");
    require_near(mass_flux[face_id],
                 -face_pressure_response[face_id] * (pressure[adjacency.neighbor] - pressure[adjacency.owner]),
                 test_tolerance, "The flux does not directly contain the expected adjacent pressure jump.");
    require(std::abs(mass_flux[face_id]) > test_tolerance,
            "The direct pressure-jump test produced a trivial Rhie-Chow flux.");

    pressure[adjacency.owner] = 4.0;
    pressure[adjacency.neighbor] = -1.0;
    interpolation.update_internal_faces(velocity, pressure, pressure_gradient, momentum_response, mass_flux,
                                        face_pressure_response);
    require(mass_flux[face_id] > 0.0,
            "A neighbor-to-owner pressure decrease did not reverse the owner-oriented flux sign.");
    require_near(mass_flux[face_id],
                 -face_pressure_response[face_id] * (pressure[adjacency.neighbor] - pressure[adjacency.owner]),
                 test_tolerance, "The reversed pressure jump produced an incorrect direct flux.");
}

void test_isotropic_orthogonal_reduction()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    constexpr double density{2.0};
    constexpr double scalar_response{0.75};
    const cfd::RhieChowInternalFaceInterpolation interpolation{mesh, density};
    const cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    const cfd::CellVectorField pressure_gradient{mesh.cell_count()};
    cfd::CellMomentumPressureResponse momentum_response{mesh.cell_count()};
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        momentum_response.u()[cell_id] = scalar_response;
        momentum_response.v()[cell_id] = scalar_response;
    }
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count()};
    const cfd::Index face_id{internal_face_id(mesh)};
    const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};
    const cfd::Vector2 &area_vector{mesh.face_area_vectors()[face_id]};
    const cfd::Point2 &owner_center{mesh.cell_centers()[adjacency.owner]};
    const cfd::Point2 &neighbor_center{mesh.cell_centers()[adjacency.neighbor]};
    const cfd::Vector2 displacement{
        neighbor_center.x - owner_center.x,
        neighbor_center.y - owner_center.y,
    };
    pressure[adjacency.owner] = 1.0;
    pressure[adjacency.neighbor] = 4.0;

    interpolation.update_internal_faces(velocity, pressure, pressure_gradient, momentum_response, mass_flux,
                                        face_pressure_response);

    const double expected_response{density * scalar_response * dot(area_vector, area_vector) /
                                   dot(area_vector, displacement)};
    require_near(face_pressure_response[face_id], expected_response, test_tolerance,
                 "Isotropic orthogonal interpolation produced an incorrect pressure coefficient.");
    require_near(mass_flux[face_id], -expected_response * (pressure[adjacency.neighbor] - pressure[adjacency.owner]),
                 test_tolerance,
                 "Isotropic orthogonal interpolation did not reduce to the conventional pressure-jump form.");
}

void test_anisotropic_face_pressure_response()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_sheared_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    constexpr double density{1.6};
    const cfd::RhieChowInternalFaceInterpolation interpolation{mesh, density};
    const cfd::CellVelocityField velocity{mesh.cell_count()};
    const cfd::CellScalarField pressure{mesh.cell_count()};
    const cfd::CellVectorField pressure_gradient{mesh.cell_count()};
    cfd::CellMomentumPressureResponse momentum_response{mesh.cell_count()};
    momentum_response.u()[0] = 0.2;
    momentum_response.v()[0] = 1.7;
    momentum_response.u()[1] = 1.0;
    momentum_response.v()[1] = 0.5;
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count()};

    interpolation.update_internal_faces(velocity, pressure, pressure_gradient, momentum_response, mass_flux,
                                        face_pressure_response);

    const cfd::Index face_id{internal_face_id(mesh)};
    const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};
    const cfd::Vector2 &area_vector{mesh.face_area_vectors()[face_id]};
    require(area_vector.x != 0.0 && area_vector.y != 0.0,
            "The anisotropic-response fixture does not have both area-vector components.");
    const cfd::Point2 &owner_center{mesh.cell_centers()[adjacency.owner]};
    const cfd::Point2 &neighbor_center{mesh.cell_centers()[adjacency.neighbor]};
    const cfd::Vector2 displacement{
        neighbor_center.x - owner_center.x,
        neighbor_center.y - owner_center.y,
    };
    const double lambda{interpolation_weight(mesh, face_id)};
    const double face_u_response{(1.0 - lambda) * momentum_response.u()[adjacency.owner] +
                                 lambda * momentum_response.u()[adjacency.neighbor]};
    const double face_v_response{(1.0 - lambda) * momentum_response.v()[adjacency.owner] +
                                 lambda * momentum_response.v()[adjacency.neighbor]};
    require(face_u_response != face_v_response,
            "The anisotropic-response fixture accidentally produced equal face responses.");
    const double expected_response{
        density * (face_u_response * area_vector.x * area_vector.x + face_v_response * area_vector.y * area_vector.y) /
        dot(area_vector, displacement)};
    require_near(face_pressure_response[face_id], expected_response, test_tolerance,
                 "Anisotropic interpolation collapsed the independent momentum responses.");
}

void test_valid_update_does_not_modify_inputs()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_sheared_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::RhieChowInternalFaceInterpolation interpolation{mesh, 1.2};
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    cfd::CellVectorField pressure_gradient{mesh.cell_count()};
    cfd::CellMomentumPressureResponse momentum_response{mesh.cell_count()};
    initialize_valid_inputs(velocity, pressure, pressure_gradient, momentum_response);
    const cfd::CellVelocityField velocity_before{velocity};
    const cfd::CellScalarField pressure_before{pressure};
    const cfd::CellVectorField pressure_gradient_before{pressure_gradient};
    const cfd::CellMomentumPressureResponse momentum_response_before{momentum_response};
    const cfd::Index node_count_before{mesh.node_count()};
    const cfd::Index cell_count_before{mesh.cell_count()};
    const cfd::Index face_count_before{mesh.face_count()};
    const cfd::Index face_id{internal_face_id(mesh)};
    const cfd::Point2 face_center_before{mesh.face_centers()[face_id]};
    const cfd::Vector2 area_vector_before{mesh.face_area_vectors()[face_id]};
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count()};

    interpolation.update_internal_faces(velocity, pressure, pressure_gradient, momentum_response, mass_flux,
                                        face_pressure_response);

    require(mesh.node_count() == node_count_before && mesh.cell_count() == cell_count_before &&
                mesh.face_count() == face_count_before,
            "Rhie-Chow interpolation changed Mesh cardinalities.");
    require(mesh.face_centers()[face_id].x == face_center_before.x &&
                mesh.face_centers()[face_id].y == face_center_before.y &&
                mesh.face_area_vectors()[face_id].x == area_vector_before.x &&
                mesh.face_area_vectors()[face_id].y == area_vector_before.y,
            "Rhie-Chow interpolation changed Mesh geometry.");
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        require(velocity.u()[cell_id] == velocity_before.u()[cell_id] &&
                    velocity.v()[cell_id] == velocity_before.v()[cell_id],
                "Rhie-Chow interpolation changed the velocity field.");
        require(pressure[cell_id] == pressure_before[cell_id], "Rhie-Chow interpolation changed pressure.");
        require(pressure_gradient[cell_id].x == pressure_gradient_before[cell_id].x &&
                    pressure_gradient[cell_id].y == pressure_gradient_before[cell_id].y,
                "Rhie-Chow interpolation changed the pressure gradient.");
        require(momentum_response.u()[cell_id] == momentum_response_before.u()[cell_id] &&
                    momentum_response.v()[cell_id] == momentum_response_before.v()[cell_id],
                "Rhie-Chow interpolation changed the momentum response.");
    }
}

void test_rejects_incompatible_cardinalities_before_mutation()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::RhieChowInternalFaceInterpolation interpolation{mesh, 1.0};
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    cfd::CellVectorField pressure_gradient{mesh.cell_count()};
    cfd::CellMomentumPressureResponse momentum_response{mesh.cell_count()};
    initialize_valid_inputs(velocity, pressure, pressure_gradient, momentum_response);
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count()};
    seed_outputs(mass_flux, face_pressure_response);

    const cfd::CellVelocityField wrong_velocity{mesh.cell_count() + 1};
    require_rejected_without_output_mutation<std::invalid_argument>(
        [&]() {
            interpolation.update_internal_faces(wrong_velocity, pressure, pressure_gradient, momentum_response,
                                                mass_flux, face_pressure_response);
        },
        mass_flux, face_pressure_response, "Rhie-Chow interpolation accepted an incorrectly sized velocity.");

    const cfd::CellScalarField wrong_pressure{mesh.cell_count() + 1};
    require_rejected_without_output_mutation<std::invalid_argument>(
        [&]() {
            interpolation.update_internal_faces(velocity, wrong_pressure, pressure_gradient, momentum_response,
                                                mass_flux, face_pressure_response);
        },
        mass_flux, face_pressure_response, "Rhie-Chow interpolation accepted an incorrectly sized pressure.");

    const cfd::CellVectorField wrong_gradient{mesh.cell_count() + 1};
    require_rejected_without_output_mutation<std::invalid_argument>(
        [&]() {
            interpolation.update_internal_faces(velocity, pressure, wrong_gradient, momentum_response, mass_flux,
                                                face_pressure_response);
        },
        mass_flux, face_pressure_response, "Rhie-Chow interpolation accepted an incorrectly sized pressure gradient.");

    const cfd::CellMomentumPressureResponse wrong_response{mesh.cell_count() + 1};
    require_rejected_without_output_mutation<std::invalid_argument>(
        [&]() {
            interpolation.update_internal_faces(velocity, pressure, pressure_gradient, wrong_response, mass_flux,
                                                face_pressure_response);
        },
        mass_flux, face_pressure_response, "Rhie-Chow interpolation accepted an incorrectly sized momentum response.");

    cfd::FaceFluxField wrong_mass_flux{mesh.face_count() + 1};
    seed_outputs(wrong_mass_flux, face_pressure_response);
    require_rejected_without_output_mutation<std::invalid_argument>(
        [&]() {
            interpolation.update_internal_faces(velocity, pressure, pressure_gradient, momentum_response,
                                                wrong_mass_flux, face_pressure_response);
        },
        wrong_mass_flux, face_pressure_response, "Rhie-Chow interpolation accepted an incorrectly sized mass flux.");

    cfd::FacePressureResponseField wrong_face_response{mesh.face_count() + 1};
    seed_outputs(mass_flux, wrong_face_response);
    require_rejected_without_output_mutation<std::invalid_argument>(
        [&]() {
            interpolation.update_internal_faces(velocity, pressure, pressure_gradient, momentum_response, mass_flux,
                                                wrong_face_response);
        },
        mass_flux, wrong_face_response,
        "Rhie-Chow interpolation accepted an incorrectly sized face pressure response.");
}

void test_rejects_nonfinite_inputs_before_mutation()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::RhieChowInternalFaceInterpolation interpolation{mesh, 1.0};
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    cfd::CellVectorField pressure_gradient{mesh.cell_count()};
    cfd::CellMomentumPressureResponse momentum_response{mesh.cell_count()};
    initialize_valid_inputs(velocity, pressure, pressure_gradient, momentum_response);
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count()};
    seed_outputs(mass_flux, face_pressure_response);

    const auto require_rejected = [&](const std::string &message) {
        require_rejected_without_output_mutation<std::runtime_error>(
            [&]() {
                interpolation.update_internal_faces(velocity, pressure, pressure_gradient, momentum_response, mass_flux,
                                                    face_pressure_response);
            },
            mass_flux, face_pressure_response, message);
    };

    velocity.u()[0] = std::numeric_limits<double>::quiet_NaN();
    require_rejected("Rhie-Chow interpolation accepted a NaN velocity.");
    velocity.u()[0] = 1.0;
    velocity.v()[1] = std::numeric_limits<double>::infinity();
    require_rejected("Rhie-Chow interpolation accepted an infinite velocity.");
    velocity.v()[1] = -1.5;

    pressure[0] = std::numeric_limits<double>::quiet_NaN();
    require_rejected("Rhie-Chow interpolation accepted a NaN pressure.");
    pressure[0] = 3.0;
    pressure[1] = std::numeric_limits<double>::infinity();
    require_rejected("Rhie-Chow interpolation accepted an infinite pressure.");
    pressure[1] = 2.0;

    pressure_gradient[0].x = std::numeric_limits<double>::quiet_NaN();
    require_rejected("Rhie-Chow interpolation accepted a NaN pressure gradient.");
    pressure_gradient[0].x = 0.25;
    pressure_gradient[1].y = -std::numeric_limits<double>::infinity();
    require_rejected("Rhie-Chow interpolation accepted an infinite pressure gradient.");
}

void test_rejects_invalid_or_overflowing_responses_before_mutation()
{
    cfd::MeshBuildResult build_result{cfd::build_mesh(make_two_cell_rectangle_raw_mesh())};
    const cfd::Mesh &mesh{build_result.mesh};
    const cfd::RhieChowInternalFaceInterpolation interpolation{mesh, 1.0};
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    cfd::CellVectorField pressure_gradient{mesh.cell_count()};
    cfd::CellMomentumPressureResponse momentum_response{mesh.cell_count()};
    initialize_valid_inputs(velocity, pressure, pressure_gradient, momentum_response);
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    cfd::FacePressureResponseField face_pressure_response{mesh.face_count()};
    seed_outputs(mass_flux, face_pressure_response);

    const auto require_rejected = [&](const std::string &message) {
        require_rejected_without_output_mutation<std::runtime_error>(
            [&]() {
                interpolation.update_internal_faces(velocity, pressure, pressure_gradient, momentum_response, mass_flux,
                                                    face_pressure_response);
            },
            mass_flux, face_pressure_response, message);
    };

    momentum_response.u()[0] = 0.0;
    require_rejected("Rhie-Chow interpolation accepted a zero momentum response.");
    momentum_response.u()[0] = 0.5;
    momentum_response.v()[0] = -1.0;
    require_rejected("Rhie-Chow interpolation accepted a negative momentum response.");
    momentum_response.v()[0] = 1.0;
    momentum_response.u()[1] = std::numeric_limits<double>::quiet_NaN();
    require_rejected("Rhie-Chow interpolation accepted a NaN momentum response.");
    momentum_response.u()[1] = 0.75;
    momentum_response.v()[1] = std::numeric_limits<double>::infinity();
    require_rejected("Rhie-Chow interpolation accepted an infinite momentum response.");

    momentum_response.v()[1] = std::numeric_limits<double>::max();
    momentum_response.v()[0] = std::numeric_limits<double>::max();
    momentum_response.u()[0] = std::numeric_limits<double>::max();
    momentum_response.u()[1] = std::numeric_limits<double>::max();
    const cfd::RhieChowInternalFaceInterpolation high_density_interpolation{mesh, 2.0};
    require_rejected_without_output_mutation<std::runtime_error>(
        [&]() {
            high_density_interpolation.update_internal_faces(velocity, pressure, pressure_gradient, momentum_response,
                                                             mass_flux, face_pressure_response);
        },
        mass_flux, face_pressure_response, "Rhie-Chow interpolation accepted an overflowing face pressure response.");
}

} // namespace

int main()
{
    int failure_count{};

    failure_count += cfd::test::run_test("Rhie-Chow density validation", test_constructor_rejects_invalid_density);
    failure_count +=
        cfd::test::run_test("Rhie-Chow constant-pressure and boundary semantics",
                            test_constant_pressure_reduces_to_velocity_interpolation_and_preserves_boundaries);
    failure_count +=
        cfd::test::run_test("Rhie-Chow additive pressure invariance", test_additive_pressure_shift_invariance);
    failure_count += cfd::test::run_test("Rhie-Chow linear pressure exactness",
                                         test_linear_pressure_is_exact_on_sheared_anisotropic_face);
    failure_count += cfd::test::run_test("Rhie-Chow direct pressure-jump sign",
                                         test_checkerboard_pressure_jump_sensitivity_and_owner_oriented_sign);
    failure_count +=
        cfd::test::run_test("Rhie-Chow isotropic orthogonal reduction", test_isotropic_orthogonal_reduction);
    failure_count +=
        cfd::test::run_test("Rhie-Chow anisotropic pressure response", test_anisotropic_face_pressure_response);
    failure_count += cfd::test::run_test("Rhie-Chow input immutability", test_valid_update_does_not_modify_inputs);
    failure_count += cfd::test::run_test("Rhie-Chow cardinality validation",
                                         test_rejects_incompatible_cardinalities_before_mutation);
    failure_count +=
        cfd::test::run_test("Rhie-Chow finite-input validation", test_rejects_nonfinite_inputs_before_mutation);
    failure_count += cfd::test::run_test("Rhie-Chow response validation",
                                         test_rejects_invalid_or_overflowing_responses_before_mutation);

    return cfd::test::finish_tests(failure_count, "Rhie-Chow internal-face interpolation");
}
