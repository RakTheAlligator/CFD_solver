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
#include <memory>
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
            .generation_options = {.mesh_size = 0.1, .cell_type = cfd::CellType::Quadrilateral},
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

void test_grid_level_report_completion_is_explicit()
{
    cfd::app::GridLevelReport report;
    require(!report.iteration_count.has_value(), "A new report claims a completed solve.");
    report.iteration_count.emplace(0);
    require(report.iteration_count.has_value() && *report.iteration_count == 0,
            "A completed zero-iteration solve is indistinguishable from an unexecuted level.");
}

void test_grid_sequencing_plan()
{
    const auto plan{cfd::app::make_grid_sequencing_plan(1600, 0.1, {})};
    require(plan.size() == 3, "Expected three levels including the final Mesh.");
    require(plan.at(0).target_cell_count == 100 && plan.at(1).target_cell_count == 400 &&
                plan.at(2).target_cell_count == 1600,
            "Wrong automatic cell targets.");
    require_near(plan.at(0).mesh_size, 0.4, 1e-12, "Wrong first automatic size.");
    require_near(plan.at(1).mesh_size, 0.2, 1e-12, "Wrong intermediate automatic size.");
    require_near(plan.at(2).mesh_size, 0.1, 1e-12, "Final size changed.");
    const auto small{cfd::app::make_grid_sequencing_plan(7, 0.1, {})};
    require(small.at(0).target_cell_count == 1 && small.at(1).target_cell_count == 1,
            "Cell targets did not retain the lower bound of one.");
    const auto explicit_plan{cfd::app::make_grid_sequencing_plan(1600, 0.1, {.target_coarse_cell_count = 25})};
    require(explicit_plan.at(0).target_cell_count == 25, "Explicit target is not the coarsest target.");
    require_near(explicit_plan.at(0).mesh_size, 0.8, 1e-12, "Wrong explicit first size.");
    require_near(explicit_plan.at(0).mesh_size / explicit_plan.at(1).mesh_size,
                 explicit_plan.at(1).mesh_size / explicit_plan.at(2).mesh_size, 1e-12,
                 "Explicit sizes are not geometrically spaced.");
    require(explicit_plan.at(1).target_cell_count == 200, "Wrong explicit intermediate target.");
    require(cfd::app::make_grid_sequencing_plan(1600, 0.1, {.type = cfd::input::InitializationType::Zero}).empty(),
            "Zero mode constructed a plan.");
    for (const cfd::Index target : {cfd::Index{0}, cfd::Index{1600}, cfd::Index{1601}})
    {
        require_throws_with_message<std::invalid_argument>(
            [target]() {
                static_cast<void>(cfd::app::make_grid_sequencing_plan(1600, 0.1, {.target_coarse_cell_count = target}));
            },
            "smaller than the final", "Accepted invalid explicit target.");
    }
    for (double size : {0.0, -1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
    {
        require_throws_with_message<std::invalid_argument>(
            [size]() { static_cast<void>(cfd::app::make_grid_sequencing_plan(1600, size, {})); }, "positive finite",
            "Accepted an invalid final mesh size.");
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
    const auto result{cfd::app::initialize_with_grid_sequencing(
        mesh, input, {.type = cfd::input::InitializationType::Zero}, u_input, v_input, p_input, density, viscosity,
        cfd::ScalarConvectionScheme::Linear, simple_options(), velocity, pressure)};
    require(!result.used_grid_sequencing && !result.automatic_fallback && result.levels.empty(),
            "Zero mode generated a coarse mesh.");
    for (cfd::Index cell = 0; cell < mesh.cell_count(); ++cell)
    {
        require(velocity.u()[cell] == u_input.internal_value && velocity.v()[cell] == v_input.internal_value &&
                    pressure[cell] == p_input.internal_value,
                "Zero mode changed historical internalField values.");
    }
    require_values_near(flux.values(), original_flux.values(), 0.0, "Zero mode changed target flux.");
}

cfd::Mesh make_small_mesh()
{
    cfd::RawMeshData raw;
    raw.nodes = {{0, 0}, {0.5, 0}, {1, 0}, {0, 0.5}, {0.5, 0.5}, {1, 0.5}, {0, 1}, {0.5, 1}, {1, 1}};
    raw.cell_types = std::vector<cfd::CellType>(4, cfd::CellType::Quadrilateral);
    raw.cell_nodes = {0, 1, 4, 3, 1, 2, 5, 4, 3, 4, 7, 6, 4, 5, 8, 7};
    raw.cell_node_offsets = {0, 4, 8, 12, 16};
    raw.boundary_groups = {{0, "wall"}};
    raw.boundary_edges = {{{0, 1}, 0}, {{1, 2}, 0}, {{2, 5}, 0}, {{5, 8}, 0},
                          {{8, 7}, 0}, {{7, 6}, 0}, {{6, 3}, 0}, {{3, 0}, 0}};
    return cfd::build_mesh(std::move(raw)).mesh;
}

void test_automatic_fallback_and_explicit_target_errors()
{
    const auto input{mesh_input()};
    const auto single{cfd::build_mesh(cfd::test::make_single_quadrilateral_raw_mesh()).mesh};
    const auto u_input{velocity_input()};
    const auto p_input{pressure_input()};
    cfd::CellVelocityField velocity{single.cell_count(), {2.0, 3.0}};
    cfd::CellScalarField pressure{single.cell_count(), 4.0};
    const auto fallback{cfd::app::initialize_with_grid_sequencing(single, input, {}, u_input, u_input, p_input, density,
                                                                  viscosity, cfd::ScalarConvectionScheme::Linear,
                                                                  simple_options(), velocity, pressure)};
    require(fallback.automatic_fallback && !fallback.used_grid_sequencing, "One-cell mesh did not fall back.");
    require(std::ranges::all_of(fallback.levels, [](const auto &level) { return !level.iteration_count.has_value(); }),
            "An early fallback reported a completed level.");
    require(velocity.u()[0] == 2.0 && velocity.v()[0] == 3.0 && pressure[0] == 4.0, "Fallback changed fields.");
    require_throws_with_message<std::invalid_argument>(
        [&]() {
            static_cast<void>(cfd::app::initialize_with_grid_sequencing(
                single, input, {.target_coarse_cell_count = 1}, u_input, u_input, p_input, density, viscosity,
                cfd::ScalarConvectionScheme::Linear, simple_options(), velocity, pressure));
        },
        "smaller than the final", "Accepted target equal to final count.");

    // An inconsistent supplied final resolution deterministically exercises the actual-count guard.
    const auto small{make_small_mesh()};
    auto inconsistent_input{input};
    inconsistent_input.generation_options.cell_type = cfd::CellType::Triangle;
    cfd::CellVelocityField small_velocity{small.cell_count(), {2.0, 3.0}};
    cfd::CellScalarField small_pressure{small.cell_count(), 4.0};
    const auto generated_fallback{cfd::app::initialize_with_grid_sequencing(
        small, inconsistent_input, {}, u_input, u_input, p_input, density, viscosity,
        cfd::ScalarConvectionScheme::Linear, simple_options(), small_velocity, small_pressure)};
    require(generated_fallback.automatic_fallback &&
                generated_fallback.levels.front().actual_cell_count >= small.cell_count(),
            "Actual-count failure did not fall back.");
    require(std::ranges::all_of(generated_fallback.levels,
                                [](const auto &level) { return !level.iteration_count.has_value(); }),
            "A constructed but rejected level was reported as completed.");
    require_throws_with_message<std::invalid_argument>(
        [&]() {
            static_cast<void>(cfd::app::initialize_with_grid_sequencing(
                small, inconsistent_input, {.target_coarse_cell_count = 1}, u_input, u_input, p_input, density,
                viscosity, cfd::ScalarConvectionScheme::Linear, simple_options(), small_velocity, small_pressure));
        },
        "increasing grid hierarchy", "An inconsistent explicit hierarchy was accepted.");

    auto coarse_quad_request{input};
    coarse_quad_request.generation_options.mesh_size = 1.0;
    const auto failed_meshing{cfd::app::initialize_with_grid_sequencing(
        small, coarse_quad_request, {}, u_input, u_input, p_input, density, viscosity,
        cfd::ScalarConvectionScheme::Linear, simple_options(), small_velocity, small_pressure)};
    require(failed_meshing.automatic_fallback && !failed_meshing.fallback_reason.empty(),
            "Unusable meshing did not fall back.");
    for (cfd::Index cell = 0; cell < small.cell_count(); ++cell)
        require(small_velocity.u()[cell] == 2.0 && small_velocity.v()[cell] == 3.0 && small_pressure[cell] == 4.0,
                "A hierarchy fallback changed the historical fields.");
}

struct ReferenceLevel
{
    cfd::Mesh mesh;
    cfd::ScalarBoundaryConditions u_boundary;
    cfd::ScalarBoundaryConditions v_boundary;
    cfd::ScalarBoundaryConditions p_boundary;
    cfd::PressureCorrectionBoundaryConditions correction;
    cfd::CellVelocityField velocity;
    cfd::CellScalarField pressure;

    ReferenceLevel(const cfd::input::MeshInput &input, const cfd::input::ScalarFieldInput &u,
                   const cfd::input::ScalarFieldInput &v, const cfd::input::ScalarFieldInput &p)
        : mesh(build_case_mesh(input)), u_boundary(cfd::input::resolve_boundary_conditions(mesh, u)),
          v_boundary(cfd::input::resolve_boundary_conditions(mesh, v)),
          p_boundary(cfd::input::resolve_boundary_conditions(mesh, p)),
          correction(cfd::app::make_pressure_correction_boundary_conditions(p_boundary)),
          velocity(mesh.cell_count(), {u.internal_value, v.internal_value}),
          pressure(mesh.cell_count(), p.internal_value)
    {
    }
};

void test_multilevel_linear_transfer_and_final_solution()
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
    const auto initialized{cfd::app::initialize_with_grid_sequencing(
        target, input, {}, u_input, v_input, p_input, density, viscosity, cfd::ScalarConvectionScheme::Linear, options,
        warm_velocity, warm_pressure)};
    require(initialized.used_grid_sequencing && initialized.levels.size() == 3,
            "Multilevel initialization was not used.");
    require(initialized.levels.at(0).actual_cell_count < initialized.levels.at(1).actual_cell_count &&
                initialized.levels.at(1).actual_cell_count < target.cell_count(),
            "Hierarchy is not strictly increasing.");
    require(initialized.levels.back().actual_cell_count == target.cell_count() &&
                !initialized.levels.back().iteration_count.has_value(),
            "The existing final Mesh was solved or replaced.");
    require_values_near(warm_flux.values(), target_initial_flux.values(), 0.0, "Final target flux changed.");

    // An independent level-by-level reference uses BC initialization only once
    // and creates fresh target fluxes at each solve. Counts and transferred
    // fields protect both that ordering and the absence of recursive sequencing.
    std::unique_ptr<ReferenceLevel> source;
    for (cfd::Index level = 0; level + 1 < initialized.levels.size(); ++level)
    {
        auto level_input{input};
        level_input.generation_options.mesh_size = initialized.levels.at(level).plan.mesh_size;
        auto state{std::make_unique<ReferenceLevel>(level_input, u_input, v_input, p_input)};
        require(state->mesh.cell_count() == initialized.levels.at(level).actual_cell_count,
                "Meshing settings changed.");
        if (!source)
        {
            cfd::app::initialize_from_boundary_conditions(state->mesh, state->u_boundary, state->v_boundary,
                                                          state->p_boundary, state->correction, {}, 0.0,
                                                          state->velocity, state->pressure);
        }
        else
        {
            const cfd::CellFieldTransfer mapping{source->mesh, state->mesh};
            cfd::CellVectorField workspace{source->mesh.cell_count()};
            mapping.apply_linear_reconstruction(source->velocity.u(), source->u_boundary, workspace,
                                                state->velocity.u());
            mapping.apply_linear_reconstruction(source->velocity.v(), source->v_boundary, workspace,
                                                state->velocity.v());
            mapping.apply_linear_reconstruction(source->pressure, source->p_boundary, workspace, state->pressure);
        }
        source.reset();
        cfd::FaceFluxField flux{state->mesh.face_count()};
        cfd::app::initialize_fixed_mass_flux_boundaries(state->mesh, density, state->u_boundary, state->v_boundary,
                                                        state->correction, flux);
        cfd::IncompressibleSimpleSolver solver{state->mesh, density, viscosity, cfd::ScalarConvectionScheme::Linear,
                                               options};
        const auto result{solver.solve(state->u_boundary, state->v_boundary, state->p_boundary, state->correction,
                                       state->velocity, state->pressure, flux)};
        const auto &reported_iterations{initialized.levels.at(level).iteration_count};
        require(result.converged && reported_iterations.has_value() && result.iteration_count == *reported_iterations,
                "Auxiliary solve did not use the expected single initialization and target flux.");
        source = std::move(state);
    }
    cfd::CellVelocityField expected_velocity{target.cell_count()};
    cfd::CellScalarField expected_pressure{target.cell_count()};
    {
        const cfd::CellFieldTransfer mapping{source->mesh, target};
        cfd::CellVectorField workspace{source->mesh.cell_count()};
        mapping.apply_linear_reconstruction(source->velocity.u(), source->u_boundary, workspace, expected_velocity.u());
        mapping.apply_linear_reconstruction(source->velocity.v(), source->v_boundary, workspace, expected_velocity.v());
        mapping.apply_linear_reconstruction(source->pressure, source->p_boundary, workspace, expected_pressure);
    }
    source.reset();
    require_values_near(warm_velocity.u().values(), expected_velocity.u().values(), 1e-12,
                        "Incorrect multilevel u transfer.");
    require_values_near(warm_velocity.v().values(), expected_velocity.v().values(), 1e-12,
                        "Incorrect multilevel v transfer.");
    require_values_near(warm_pressure.values(), expected_pressure.values(), 1e-12, "Incorrect multilevel p transfer.");

    cfd::CellVelocityField zero_velocity{target.cell_count()};
    cfd::CellScalarField zero_pressure{target.cell_count()};
    cfd::FaceFluxField zero_flux{target_initial_flux};
    {
        cfd::IncompressibleSimpleSolver solver{target, density, viscosity, cfd::ScalarConvectionScheme::Linear,
                                               options};
        const auto result{
            solver.solve(u_boundary, v_boundary, p_boundary, correction, zero_velocity, zero_pressure, zero_flux)};
        require(result.converged, "Historical zero-start did not converge.");
    }
    {
        cfd::IncompressibleSimpleSolver solver{target, density, viscosity, cfd::ScalarConvectionScheme::Linear,
                                               options};
        const auto result{
            solver.solve(u_boundary, v_boundary, p_boundary, correction, warm_velocity, warm_pressure, warm_flux)};
        require(result.converged, "Multilevel final solution did not converge.");
    }
    require_values_near(warm_velocity.u().values(), zero_velocity.u().values(), 2e-5, "Final u changed.");
    require_values_near(warm_velocity.v().values(), zero_velocity.v().values(), 2e-5, "Final v changed.");
    require_values_near(warm_pressure.values(), zero_pressure.values(), 2e-5, "Final p changed.");
    require_values_near(warm_flux.values(), zero_flux.values(), 2e-5, "Final flux changed.");
}

void test_late_coverage_failure_preserves_final_fields()
{
    const auto input{mesh_input()};
    auto raw{cfd::generate_mesh(input.geometry, input.generation_options, input.automatic_meshing,
                                input.backward_facing_step_meshing)};
    for (auto &point : raw.nodes)
        point.x += 10.0;
    const auto target{cfd::build_mesh(std::move(raw)).mesh};
    const auto u_input{velocity_input()};
    const auto p_input{pressure_input()};
    cfd::CellVelocityField velocity{target.cell_count(), {2.0, 3.0}};
    cfd::CellScalarField pressure{target.cell_count(), 4.0};
    const auto fallback{cfd::app::initialize_with_grid_sequencing(target, input, {}, u_input, u_input, p_input, density,
                                                                  viscosity, cfd::ScalarConvectionScheme::Linear,
                                                                  simple_options(), velocity, pressure)};
    require(fallback.automatic_fallback && !fallback.used_grid_sequencing, "Missing final coverage was accepted.");
    require(fallback.levels.at(0).iteration_count.has_value() && fallback.levels.at(1).iteration_count.has_value(),
            "The coverage test did not exercise a late failure.");
    require(!fallback.levels.back().iteration_count.has_value(),
            "A late fallback reported the unexecuted final solve as completed.");
    require(fallback.fallback_reason.find("Transfer to level 2") != std::string::npos, "Unclear coverage failure.");
    for (cfd::Index cell = 0; cell < target.cell_count(); ++cell)
        require(velocity.u()[cell] == 2.0 && velocity.v()[cell] == 3.0 && pressure[cell] == 4.0,
                "Coverage fallback partially changed final fields.");
    require_throws_with_message<std::invalid_argument>(
        [&]() {
            static_cast<void>(cfd::app::initialize_with_grid_sequencing(
                target, input, {.target_coarse_cell_count = target.cell_count() / 16}, u_input, u_input, p_input,
                density, viscosity, cfd::ScalarConvectionScheme::Linear, simple_options(), velocity, pressure));
        },
        "Transfer to level 2", "Explicit coverage failure silently fell back.");
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
            static_cast<void>(cfd::app::initialize_with_grid_sequencing(
                mesh, input, {}, u_input, u_input, p_input, density, viscosity, cfd::ScalarConvectionScheme::Linear,
                options, velocity, pressure));
        },
        "Grid-sequencing SIMPLE level 0 did not converge", "Nonconverged coarse state was accepted.");
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
    failures += cfd::test::run_test("explicit level completion", test_grid_level_report_completion_is_explicit);
    failures += cfd::test::run_test("grid sequencing plan", test_grid_sequencing_plan);
    failures += cfd::test::run_test("historical zero initialization", test_zero_preserves_historical_initialization);
    failures +=
        cfd::test::run_test("coarse fallback and target errors", test_automatic_fallback_and_explicit_target_errors);
    failures += cfd::test::run_test("multilevel linear sequencing", test_multilevel_linear_transfer_and_final_solution);
    failures += cfd::test::run_test("reject unconverged coarse solve", test_nonconverged_coarse_solve_is_rejected);
    failures +=
        cfd::test::run_test("late transfer coverage fallback", test_late_coverage_failure_preserves_final_fields);
    return cfd::test::finish_tests(failures, "IncompressibleCaseInitialization");
}
