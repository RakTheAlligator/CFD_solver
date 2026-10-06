#include "app/BoundaryBasedInitialization.hpp"
#include "app/IncompressibleCaseInitialization.hpp"

#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVectorField.hpp"
#include "cfd/field/CellVelocityField.hpp"
#include "cfd/field/FaceFluxField.hpp"
#include "cfd/field/PressureCorrectionBoundaryConditions.hpp"
#include "cfd/input/OpenFOAMCaseReader.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/numerics/CellFieldTransfer.hpp"
#include "cfd/numerics/IncompressibleSimpleSolver.hpp"

#include "support/MeshFixtures.hpp"
#include "support/TestUtils.hpp"

#include <algorithm>
#include <limits>
#include <span>
#include <vector>

namespace
{
using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws_with_message;

constexpr double density{1.0};
constexpr double viscosity{0.1};

cfd::input::MeshInput mesh_input()
{
    return {.geometry = cfd::RectangleGeometry{2.0, 1.0},
            .generation_options = {.mesh_size = 0.2, .cell_type = cfd::CellType::Quadrilateral},
            .automatic_meshing = {},
            .backward_facing_step_meshing = {}};
}

cfd::Mesh build_case_mesh(const cfd::input::MeshInput &input)
{
    return cfd::build_mesh(cfd::generate_mesh(input.geometry, input.generation_options, input.automatic_meshing,
                                              input.backward_facing_step_meshing))
        .mesh;
}

cfd::input::ScalarFieldInput velocity_input()
{
    using Type = cfd::ScalarBoundaryConditionType;
    return {.object_name = "u",
            .dimensions = {},
            .internal_value = 0.0,
            .boundary_conditions = {
                {"inlet", {Type::Neumann, 0.0}}, {"wall", {Type::Dirichlet, 0.0}}, {"outlet", {Type::Neumann, 0.0}}}};
}

cfd::input::ScalarFieldInput pressure_input()
{
    using Type = cfd::ScalarBoundaryConditionType;
    return {.object_name = "p",
            .dimensions = {},
            .internal_value = 0.0,
            .boundary_conditions = {{"inlet", {Type::Dirichlet, 0.04}},
                                    {"wall", {Type::Neumann, 0.0}},
                                    {"outlet", {Type::Dirichlet, 0.0}}}};
}

cfd::IncompressibleSimpleOptions simple_options()
{
    return {.maximum_iterations = 2000,
            .momentum_relaxation_factor = 1.0,
            .pressure_relaxation_factor = 0.1,
            .rhie_chow_flux_relaxation_factor = 0.3,
            .velocity_relative_tolerance = 1.0e-10,
            .rhie_chow_flux_relative_tolerance = 1.0e-10,
            .continuity_relative_tolerance = 1.0e-10,
            .momentum_linear_solver = {.relative_tolerance = 1.0e-6, .maximum_iterations = 5000},
            .pressure_correction_linear_solver = {.relative_tolerance = 1.0e-3, .maximum_iterations = 5000}};
}

void require_values_near(std::span<const double> actual, std::span<const double> expected, double tolerance,
                         const std::string &message)
{
    require(actual.size() == expected.size(), message + " cardinality mismatch.");
    for (cfd::Index i = 0; i < actual.size(); ++i)
    {
        require_near(actual[i], expected[i], tolerance, message);
    }
}

void test_coarse_mesh_plan()
{
    const cfd::input::InitializationInput automatic;
    const auto plan{cfd::app::make_coarse_mesh_plan(100, 0.2, automatic)};
    require(plan.target_cell_count == 25, "Automatic target is not Nf/4.");
    require_near(plan.mesh_size, 0.4, 1.0e-12, "Automatic coarse size is not 2h.");
    require(cfd::app::make_coarse_mesh_plan(7, 0.2, automatic).target_cell_count == 1,
            "Automatic target did not use integer division.");
    require(cfd::app::make_coarse_mesh_plan(1, 0.2, automatic).target_cell_count == 1,
            "Automatic target did not retain the lower bound of one.");
    const cfd::input::InitializationInput explicit_target{.target_coarse_cell_count = 4};
    require_near(cfd::app::make_coarse_mesh_plan(100, 0.2, explicit_target).mesh_size, 1.0, 1.0e-12,
                 "Explicit target conversion is incorrect.");
    for (const cfd::Index target : {cfd::Index{0}, cfd::Index{100}, cfd::Index{101}})
    {
        require_throws_with_message<std::invalid_argument>(
            [target]() {
                static_cast<void>(cfd::app::make_coarse_mesh_plan(100, 0.2, {.target_coarse_cell_count = target}));
            },
            "smaller than the final", "Accepted invalid explicit target.");
    }
}

void test_zero_preserves_historical_initialization()
{
    const auto input{mesh_input()};
    const auto mesh{build_case_mesh(input)};
    auto u_input{velocity_input()};
    auto v_input{velocity_input()};
    auto p_input{pressure_input()};
    u_input.internal_value = 0.125;
    v_input.internal_value = -0.25;
    p_input.internal_value = 0.5;
    const auto u_boundary{cfd::input::resolve_boundary_conditions(mesh, u_input)};
    const auto v_boundary{cfd::input::resolve_boundary_conditions(mesh, v_input)};
    const auto p_boundary{cfd::input::resolve_boundary_conditions(mesh, p_input)};
    const auto correction{cfd::app::make_pressure_correction_boundary_conditions(p_boundary)};
    cfd::CellVelocityField velocity{mesh.cell_count(), {u_input.internal_value, v_input.internal_value}};
    cfd::CellScalarField pressure{mesh.cell_count(), p_input.internal_value};
    cfd::FaceFluxField flux{mesh.face_count()};
    cfd::app::initialize_fixed_mass_flux_boundaries(mesh, density, u_boundary, v_boundary, correction, flux);
    const cfd::FaceFluxField original_flux{flux};
    const auto result{cfd::app::initialize_from_coarse_mesh(
        mesh, input, {.type = cfd::input::InitializationType::Zero}, u_input, v_input, p_input, density, viscosity,
        cfd::ScalarConvectionScheme::Linear, simple_options(), velocity, pressure)};
    require(!result.used_coarse_mesh && !result.automatic_fallback && result.actual_cell_count == 0,
            "Zero mode generated a coarse mesh.");
    for (cfd::Index cell = 0; cell < mesh.cell_count(); ++cell)
    {
        require(velocity.u()[cell] == u_input.internal_value && velocity.v()[cell] == v_input.internal_value &&
                    pressure[cell] == p_input.internal_value,
                "Zero mode changed historical internalField values.");
    }
    require_values_near(flux.values(), original_flux.values(), 0.0, "Zero mode changed target flux.");
}

void test_automatic_fallback_and_explicit_target_errors()
{
    const auto input{mesh_input()};
    const auto single{cfd::build_mesh(cfd::test::make_single_quadrilateral_raw_mesh()).mesh};
    const auto u_input{velocity_input()};
    const auto p_input{pressure_input()};
    cfd::CellVelocityField velocity{single.cell_count(), {2.0, 3.0}};
    cfd::CellScalarField pressure{single.cell_count(), 4.0};
    const auto fallback{cfd::app::initialize_from_coarse_mesh(single, input, {}, u_input, u_input, p_input, density,
                                                              viscosity, cfd::ScalarConvectionScheme::Linear,
                                                              simple_options(), velocity, pressure)};
    require(fallback.automatic_fallback && !fallback.used_coarse_mesh, "One-cell mesh did not fall back.");
    require(velocity.u()[0] == 2.0 && velocity.v()[0] == 3.0 && pressure[0] == 4.0,
            "Fallback changed the initial fields.");
    require_throws_with_message<std::invalid_argument>(
        [&]() {
            static_cast<void>(cfd::app::initialize_from_coarse_mesh(
                single, input, {.target_coarse_cell_count = 1}, u_input, u_input, p_input, density, viscosity,
                cfd::ScalarConvectionScheme::Linear, simple_options(), velocity, pressure));
        },
        "smaller than the final", "Accepted target equal to final cell count.");

    // Deliberately inconsistent final resolution exercises the generated-size guard.
    const auto small{cfd::build_mesh(cfd::test::make_two_triangle_raw_mesh()).mesh};
    auto inconsistent_input{input};
    inconsistent_input.generation_options.cell_type = cfd::CellType::Triangle;
    cfd::CellVelocityField small_velocity{small.cell_count(), {2.0, 3.0}};
    cfd::CellScalarField small_pressure{small.cell_count(), 4.0};
    const auto generated_fallback{cfd::app::initialize_from_coarse_mesh(
        small, inconsistent_input, {}, u_input, u_input, p_input, density, viscosity,
        cfd::ScalarConvectionScheme::Linear, simple_options(), small_velocity, small_pressure)};
    require(generated_fallback.automatic_fallback && generated_fallback.actual_cell_count >= small.cell_count(),
            "Generated coarse mesh that is not smaller did not fall back.");
    require_throws_with_message<std::invalid_argument>(
        [&]() {
            static_cast<void>(cfd::app::initialize_from_coarse_mesh(
                small, inconsistent_input, {.target_coarse_cell_count = 1}, u_input, u_input, p_input, density,
                viscosity, cfd::ScalarConvectionScheme::Linear, simple_options(), small_velocity, small_pressure));
        },
        "did not produce a smaller", "Explicit target silently accepted a non-smaller mesh.");

    // Gmsh can reject very coarse all-QUAD requests. Automatic mode must not
    // make an otherwise usable final mesh depend on this auxiliary mesh.
    auto coarse_quad_request{input};
    coarse_quad_request.generation_options.mesh_size = 1.0;
    const auto coarse_failure{cfd::app::initialize_from_coarse_mesh(
        small, coarse_quad_request, {}, u_input, u_input, p_input, density, viscosity,
        cfd::ScalarConvectionScheme::Linear, simple_options(), small_velocity, small_pressure)};
    require(coarse_failure.automatic_fallback && !coarse_failure.used_coarse_mesh &&
                !coarse_failure.fallback_reason.empty(),
            "Unusable automatic coarse meshing did not fall back with an explanation.");
    require_throws_with_message<std::invalid_argument>(
        [&]() {
            static_cast<void>(cfd::app::initialize_from_coarse_mesh(
                small, coarse_quad_request, {.target_coarse_cell_count = 1}, u_input, u_input, p_input, density,
                viscosity, cfd::ScalarConvectionScheme::Linear, simple_options(), small_velocity, small_pressure));
        },
        "coarse mesh", "Unusable explicit coarse meshing was silently accepted.");
}

void test_single_coarse_solve_linear_transfer_and_final_solution()
{
    const auto input{mesh_input()};
    const auto target{build_case_mesh(input)};
    const auto u_input{velocity_input()};
    const auto v_input{velocity_input()};
    const auto p_input{pressure_input()};
    const auto options{simple_options()};
    const auto u_boundary{cfd::input::resolve_boundary_conditions(target, u_input)};
    const auto v_boundary{cfd::input::resolve_boundary_conditions(target, v_input)};
    const auto p_boundary{cfd::input::resolve_boundary_conditions(target, p_input)};
    const auto correction{cfd::app::make_pressure_correction_boundary_conditions(p_boundary)};
    cfd::CellVelocityField warm_velocity{target.cell_count()};
    cfd::CellScalarField warm_pressure{target.cell_count()};
    cfd::FaceFluxField warm_flux{target.face_count()};
    cfd::app::initialize_fixed_mass_flux_boundaries(target, density, u_boundary, v_boundary, correction, warm_flux);
    const cfd::FaceFluxField target_initial_flux{warm_flux};
    const auto initialized{cfd::app::initialize_from_coarse_mesh(target, input, {}, u_input, v_input, p_input, density,
                                                                 viscosity, cfd::ScalarConvectionScheme::Linear,
                                                                 options, warm_velocity, warm_pressure)};
    require(initialized.used_coarse_mesh && initialized.actual_cell_count < target.cell_count(),
            "Coarse initialization was not used.");

    // Independent single-Mesh solve gives the expected source state, with no sequencing.
    auto coarse_input{input};
    coarse_input.generation_options.mesh_size = initialized.plan.mesh_size;
    const auto source{build_case_mesh(coarse_input)};
    require(initialized.actual_cell_count == source.cell_count(), "Coarse meshing options changed.");
    const auto source_u_boundary{cfd::input::resolve_boundary_conditions(source, u_input)};
    const auto source_v_boundary{cfd::input::resolve_boundary_conditions(source, v_input)};
    const auto source_p_boundary{cfd::input::resolve_boundary_conditions(source, p_input)};
    const auto source_correction{cfd::app::make_pressure_correction_boundary_conditions(source_p_boundary)};
    cfd::CellVelocityField source_velocity{source.cell_count()};
    cfd::CellScalarField source_pressure{source.cell_count()};
    cfd::FaceFluxField source_flux{source.face_count()};
    cfd::app::initialize_from_boundary_conditions(source, source_u_boundary, source_v_boundary, source_p_boundary,
                                                  source_correction, {}, 0.0, source_velocity, source_pressure);
    cfd::app::initialize_fixed_mass_flux_boundaries(source, density, source_u_boundary, source_v_boundary,
                                                    source_correction, source_flux);
    {
        cfd::IncompressibleSimpleSolver source_solver{source, density, viscosity, cfd::ScalarConvectionScheme::Linear,
                                                      options};
        const auto result{source_solver.solve(source_u_boundary, source_v_boundary, source_p_boundary,
                                              source_correction, source_velocity, source_pressure, source_flux)};
        require(result.converged && result.iteration_count == initialized.iteration_count,
                "Initialization did not use the direct nonrecursive coarse solve.");
    }
    cfd::CellVelocityField expected_velocity{target.cell_count()};
    cfd::CellScalarField expected_pressure{target.cell_count()};
    cfd::CellVectorField workspace{source.cell_count()};
    const cfd::CellFieldTransfer mapping{source, target};
    mapping.apply_linear_reconstruction(source_velocity.u(), source_u_boundary, workspace, expected_velocity.u());
    mapping.apply_linear_reconstruction(source_velocity.v(), source_v_boundary, workspace, expected_velocity.v());
    mapping.apply_linear_reconstruction(source_pressure, source_p_boundary, workspace, expected_pressure);
    require_values_near(warm_velocity.u().values(), expected_velocity.u().values(), 1.0e-12, "Incorrect u transfer.");
    require_values_near(warm_velocity.v().values(), expected_velocity.v().values(), 1.0e-12, "Incorrect v transfer.");
    require_values_near(warm_pressure.values(), expected_pressure.values(), 1.0e-12, "Incorrect p transfer.");
    require_values_near(warm_flux.values(), target_initial_flux.values(), 0.0, "Source flux was transferred.");
    require(std::ranges::any_of(source_flux.values(), [](double value) { return value != 0.0; }),
            "Source flux test did not exercise nonzero converged flux.");

    cfd::CellVelocityField zero_velocity{target.cell_count()};
    cfd::CellScalarField zero_pressure{target.cell_count()};
    cfd::FaceFluxField zero_flux{target_initial_flux};
    {
        cfd::IncompressibleSimpleSolver zero_solver{target, density, viscosity, cfd::ScalarConvectionScheme::Linear,
                                                    options};
        const auto result{
            zero_solver.solve(u_boundary, v_boundary, p_boundary, correction, zero_velocity, zero_pressure, zero_flux)};
        require(result.converged, "Historical zero-start solution did not converge.");
    }
    {
        cfd::IncompressibleSimpleSolver warm_solver{target, density, viscosity, cfd::ScalarConvectionScheme::Linear,
                                                    options};
        const auto result{
            warm_solver.solve(u_boundary, v_boundary, p_boundary, correction, warm_velocity, warm_pressure, warm_flux)};
        require(result.converged, "Sequenced target solution did not converge.");
    }
    // Inner solves stop at 1e-6: compare physical fixed points, not bitwise trajectories.
    require_values_near(warm_velocity.u().values(), zero_velocity.u().values(), 2.0e-5, "Final u changed.");
    require_values_near(warm_velocity.v().values(), zero_velocity.v().values(), 2.0e-5, "Final v changed.");
    require_values_near(warm_pressure.values(), zero_pressure.values(), 2.0e-5, "Final p changed.");
    require_values_near(warm_flux.values(), zero_flux.values(), 2.0e-5, "Final flux changed.");
}

void test_nonconverged_coarse_solve_is_rejected()
{
    const auto input{mesh_input()};
    const auto mesh{build_case_mesh(input)};
    const auto u_input{velocity_input()};
    const auto p_input{pressure_input()};
    auto options{simple_options()};
    options.maximum_iterations = 1;
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    require_throws_with_message<std::runtime_error>(
        [&]() {
            static_cast<void>(cfd::app::initialize_from_coarse_mesh(mesh, input, {}, u_input, u_input, p_input, density,
                                                                    viscosity, cfd::ScalarConvectionScheme::Linear,
                                                                    options, velocity, pressure));
        },
        "Coarse SIMPLE did not converge", "Nonconverged coarse state was accepted.");
    for (cfd::Index cell = 0; cell < mesh.cell_count(); ++cell)
    {
        require(velocity.u()[cell] == 0.0 && velocity.v()[cell] == 0.0 && pressure[cell] == 0.0,
                "A failed coarse solve changed target fields.");
    }
}
} // namespace

int main()
{
    int failures{};
    failures += cfd::test::run_test("coarse mesh planning", test_coarse_mesh_plan);
    failures += cfd::test::run_test("historical zero initialization", test_zero_preserves_historical_initialization);
    failures +=
        cfd::test::run_test("coarse fallback and target errors", test_automatic_fallback_and_explicit_target_errors);
    failures += cfd::test::run_test("single coarse solve and linear sequencing",
                                    test_single_coarse_solve_linear_transfer_and_final_solution);
    failures += cfd::test::run_test("reject unconverged coarse solve", test_nonconverged_coarse_solve_is_rejected);
    return cfd::test::finish_tests(failures, "IncompressibleCaseInitialization");
}
