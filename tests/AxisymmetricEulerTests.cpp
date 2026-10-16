#include "cfd/compressible/AxisymmetricEulerSolver.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "support/EulerTube.hpp"
#include "support/TestUtils.hpp"

#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <type_traits>

namespace
{
using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws_with_message;
using Type = cfd::EulerBoundaryConditionType;
const cfd::CaloricallyPerfectGas gas{1.4, 1.0};
constexpr double tolerance{2.0e-12};
static_assert(!std::is_copy_constructible_v<cfd::AxisymmetricEulerSolver>);
static_assert(std::is_nothrow_move_constructible_v<cfd::AxisymmetricEulerSolver>);

void require_balance(const cfd::AxisymmetricEulerResult &result)
{
    require_near(result.balance_error.mass, 0.0, tolerance, "Mass balance defect.");
    require_near(result.balance_error.momentum_x, 0.0, tolerance, "Axial momentum balance defect.");
    require_near(result.balance_error.momentum_r, 0.0, tolerance, "Radial balance must include pressure source.");
    require_near(result.balance_error.energy, 0.0, tolerance, "Energy balance defect.");
}

void require_uniform(const cfd::ConservativeField &field, const cfd::ConservativeState expected)
{
    for (cfd::Index cell = 0; cell < field.size(); ++cell)
    {
        const auto actual{field.state(cell)};
        require_near(actual.density, expected.density, tolerance, "Uniform density changed.");
        require_near(actual.momentum_x, expected.momentum_x, tolerance, "Uniform axial momentum changed.");
        require_near(actual.momentum_r, expected.momentum_r, tolerance, "Spurious radial momentum.");
        require_near(actual.total_energy_density, expected.total_energy_density, tolerance, "Uniform energy changed.");
    }
}

void test_uniform_rest_and_axial_flow()
{
    for (const auto type : {cfd::verification::EulerTubeCells::Quad, cfd::verification::EulerTubeCells::Triangle,
                            cfd::verification::EulerTubeCells::Mixed})
    {
        const auto built{cfd::build_mesh(cfd::verification::make_euler_tube(3, 2, type))};
        for (const double ux : {0.0, 0.7})
        {
            const auto initial{gas.to_conservative({1.0, ux, 0.0, 2.0})};
            const cfd::ConservativeField field{built.mesh.cell_count(), initial};
            cfd::AxisymmetricEulerSolver solver{built.mesh, field, gas, cfd::verification::euler_tube_conditions()};
            while (solver.diagnostics().time < 0.2)
            {
                require_balance(solver.step(0.2));
                require_uniform(solver.state(), initial);
            }
            require(solver.diagnostics().accepted_steps > 1, "Uniform test must exercise several steps.");
            require(solver.diagnostics().rejected_steps == 0, "Uniform flow should not reject steps.");
            require_near(solver.diagnostics().integrals.mass, std::numbers::pi, tolerance, "Wrong cylinder mass.");
        }
    }
}

void test_pressure_source_cfl_axis_and_final_time()
{
    const auto built{cfd::build_mesh(cfd::verification::make_euler_tube(1, 1))};
    const cfd::ConservativeField field{1, gas.to_conservative({1.0, 0.0, 0.0, 1.0})};
    cfd::AxisymmetricEulerSolver solver{built.mesh, field, gas, cfd::verification::euler_tube_conditions()};
    // Cylinder: V=pi, axial areas=pi each, radial wall=2*pi, axis=0.
    // Sum(A*lambda)=4*pi*sqrt(1.4); dt=CFL/(4*sqrt(1.4)).
    const double dt{0.2 / (4.0 * std::sqrt(1.4))};
    const auto first{solver.step(1.0)};
    require_near(first.time_step, dt, tolerance, "Incorrect analytical cylinder CFL step.");
    require_near(first.boundary_outflow_rate.momentum_r, 2.0 * std::numbers::pi, tolerance,
                 "Incorrect wall pressure flux.");
    require_near(first.radial_pressure_source_rate, 2.0 * std::numbers::pi, tolerance,
                 "Incorrect geometric pressure source.");
    require(first.boundary_outflow_rate.mass == 0.0 && first.boundary_outflow_rate.energy == 0.0,
            "Rest/axis must have exactly zero mass and energy outflow.");
    require_balance(first);
    const auto final{solver.advance_to(0.123)};
    require(final.time == 0.123 && final.time_step > 0.0, "Final time must be reached without overshooting.");
    require(solver.step(0.123).accepted_steps == final.accepted_steps, "Already reached target must be a no-op.");
    require_throws_with_message<std::invalid_argument>([&] { solver.step(0.1); }, "final time",
                                                       "Past target accepted.");
}

void test_internal_face_conservation_and_closed_walls()
{
    const auto built{cfd::build_mesh(cfd::verification::make_euler_tube(2, 1))};
    cfd::ConservativeField field{2, gas.to_conservative({1.0, 0.0, 0.0, 1.0})};
    field.set_state(1, gas.to_conservative({2.0, 0.0, 0.0, 1.0}));
    cfd::AxisymmetricEulerSolver solver{
        built.mesh, field, gas, {4, {Type::Axis, Type::SlipWall, Type::SlipWall, Type::SlipWall}}};
    const auto initial{solver.diagnostics().integrals};
    const auto first{solver.step(0.001)};
    // Independent 1D Rusanov mass flux: -sqrt(1.4)*(2-1)/2, area/volume=2.
    const double increment{0.001 * std::sqrt(1.4)};
    require_near(solver.state().density()[0], 1.0 + increment, tolerance, "Owner face sign or area incorrect.");
    require_near(solver.state().density()[1], 2.0 - increment, tolerance, "Neighbor flux must be opposite.");
    require_balance(first);
    while (solver.diagnostics().time < 0.1)
    {
        const auto &result{solver.step(0.1)};
        require_balance(result);
        require(result.boundary_outflow_rate.mass == 0.0 && result.boundary_outflow_rate.energy == 0.0,
                "Slip walls must have zero mass/energy flux.");
        require_near(result.integrals.mass, initial.mass, tolerance, "Closed-tube mass changed.");
        require_near(result.integrals.energy, initial.energy, tolerance, "Closed-tube energy changed.");
    }
}

void test_oblique_slip_walls_are_impermeable()
{
    cfd::RawMeshData raw;
    raw.nodes = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}};
    raw.cell_types = {cfd::CellType::Triangle};
    raw.cell_nodes = {0, 1, 2};
    raw.cell_node_offsets = {0, 3};
    raw.boundary_groups = {{0, "axis"}, {1, "wall"}};
    raw.boundary_edges = {{{0, 1}, 0}, {{1, 2}, 1}, {{2, 0}, 1}};
    const auto built{cfd::build_mesh(std::move(raw))};
    const cfd::ConservativeField initial{1, gas.to_conservative({1.0, 0.6, 0.4, 1.0})};
    cfd::AxisymmetricEulerSolver solver{built.mesh, initial, gas, {2, {Type::Axis, Type::SlipWall}}};
    const auto integrals{solver.diagnostics().integrals};
    while (solver.diagnostics().time < 0.1)
    {
        const auto result{solver.step(0.1)};
        require(result.boundary_outflow_rate.mass == 0.0 && result.boundary_outflow_rate.energy == 0.0,
                "Oblique reflected walls must be exactly impermeable.");
        require_near(result.integrals.mass, integrals.mass, tolerance, "Oblique wall leaked mass.");
        require_near(result.integrals.energy, integrals.energy, tolerance, "Oblique wall leaked energy.");
        require_balance(result);
    }
}

void test_input_and_boundary_validation()
{
    const auto built{cfd::build_mesh(cfd::verification::make_euler_tube(1, 1))};
    const cfd::ConservativeField field{1, gas.to_conservative({1.0, 0.0, 0.0, 1.0})};
    require_throws_with_message<std::invalid_argument>([] { const cfd::EulerBoundaryConditions bc{1, {}}; }, "count",
                                                       "Missing BC accepted.");
    require_throws_with_message<std::invalid_argument>(
        [] {
            const auto unsupported{std::bit_cast<Type>(std::uint8_t{255})};
            const cfd::EulerBoundaryConditions bc{1, {unsupported}};
        },
        "Unsupported", "Unknown BC accepted.");
    for (const auto conditions : {std::array{Type::SlipWall, Type::SlipWall, Type::SlipWall, Type::SlipWall},
                                  std::array{Type::Axis, Type::Axis, Type::SlipWall, Type::SlipWall}})
    {
        require_throws_with_message<std::invalid_argument>(
            [&] {
                const cfd::AxisymmetricEulerSolver invalid{
                    built.mesh, field, gas, {4, {conditions.begin(), conditions.end()}}};
            },
            "Axis BC", "Inconsistent axis boundary accepted.");
    }
    require_throws_with_message<std::invalid_argument>(
        [&] { const cfd::AxisymmetricEulerSolver invalid{built.mesh, field, gas, {1, {Type::SlipWall}}}; }, "boundary",
        "Wrong Mesh BC cardinality accepted.");
    require_throws_with_message<std::invalid_argument>(
        [&] {
            const cfd::AxisymmetricEulerSolver invalid{built.mesh, cfd::ConservativeField{2}, gas,
                                                       cfd::verification::euler_tube_conditions()};
        },
        "cardinality", "Wrong initial cardinality accepted.");
    for (const double cfl :
         {0.0, -1.0, 1.1, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
    {
        require_throws_with_message<std::invalid_argument>(
            [&] {
                const cfd::AxisymmetricEulerSolver invalid{
                    built.mesh, field, gas, cfd::verification::euler_tube_conditions(), {.cfl = cfl}};
            },
            "CFL", "Invalid CFL accepted.");
    }
    for (const auto options :
         {cfd::AxisymmetricEulerOptions{.maximum_step_attempts = 0}, cfd::AxisymmetricEulerOptions{.maximum_steps = 0}})
        require_throws_with_message<std::invalid_argument>(
            [&] {
                const cfd::AxisymmetricEulerSolver invalid{built.mesh, field, gas,
                                                           cfd::verification::euler_tube_conditions(), options};
            },
            "limits", "Invalid step/attempt count accepted.");
    for (const cfd::ConservativeState bad :
         {cfd::ConservativeState{}, cfd::ConservativeState{1.0, 2.0, 0.0, 1.0},
          cfd::ConservativeState{1.0, 0.0, 0.0, std::numeric_limits<double>::infinity()}})
        require_throws_with_message<std::runtime_error>(
            [&] {
                const cfd::AxisymmetricEulerSolver invalid{built.mesh, cfd::ConservativeField{1, bad}, gas,
                                                           cfd::verification::euler_tube_conditions()};
            },
            "Perfect gas", "Invalid initial gas state accepted.");
    cfd::AxisymmetricEulerSolver limited{
        built.mesh, field, gas, cfd::verification::euler_tube_conditions(), {.maximum_steps = 1}};
    limited.step(1.0);
    require_throws_with_message<std::runtime_error>([&] { limited.step(1.0); }, "maximum", "Step limit ignored.");
    require_throws_with_message<std::invalid_argument>([&] { limited.step(std::numeric_limits<double>::quiet_NaN()); },
                                                       "final time", "Non-finite target accepted.");
}

void test_candidate_backtracking_is_global_and_transactional()
{
    const auto built{cfd::build_mesh(cfd::verification::make_euler_tube(1, 1))};
    // Fast radial outflow next to the axis: the explicit geometric source and
    // energy update make the face-CFL candidate inadmissible, without clipping.
    const auto initial{gas.to_conservative({1.0, 0.0, 10.0, 0.001})};
    const cfd::ConservativeField field{1, initial};
    const cfd::EulerBoundaryConditions bc{4, {Type::Axis, Type::Transmissive, Type::Transmissive, Type::Transmissive}};
    cfd::AxisymmetricEulerSolver rejected{built.mesh, field, gas, bc, {.cfl = 1.0, .maximum_step_attempts = 1}};
    require_throws_with_message<std::runtime_error>(
        [&] { rejected.step(1.0); }, "cell=0, dt=", "Exhausted candidate attempts must identify cell and dt.");
    require(rejected.diagnostics().time == 0.0 && rejected.diagnostics().accepted_steps == 0 &&
                rejected.diagnostics().rejected_steps == 1,
            "Failed step changed time or counters incorrectly.");
    require_uniform(rejected.state(), initial);
    cfd::AxisymmetricEulerSolver retried{built.mesh, field, gas, bc, {.cfl = 1.0}};
    const auto first{retried.step(1.0)};
    const double sound{std::sqrt(1.4 * 0.001)};
    require(first.rejected_steps == 1 && first.accepted_steps == 1, "Expected one global step reduction.");
    require_near(first.time_step, 0.5 / (20.0 + 4.0 * sound), tolerance, "Retry must halve the global CFL step.");
    require(first.minimum_density > 0.0 && first.minimum_pressure > 0.0, "Accepted retry must be admissible.");
    require_balance(first);
    require_uniform(field, initial);
}

void test_state_buffers_alternate_without_copying()
{
    const auto built{cfd::build_mesh(cfd::verification::make_euler_tube(1, 1))};
    const cfd::ConservativeField field{1, gas.to_conservative({1.0, 0.0, 0.0, 1.0})};
    cfd::AxisymmetricEulerSolver solver{built.mesh, field, gas, cfd::verification::euler_tube_conditions()};
    const auto *first{solver.state().density().data()};
    solver.step(1.0);
    const auto *second{solver.state().density().data()};
    require(first != second, "An accepted step must switch buffers.");
    solver.step(1.0);
    require(solver.state().density().data() == first, "Second step must reuse the original buffer.");
    solver.step(1.0);
    require(solver.state().density().data() == second, "Third step must reuse candidate storage.");
}

void test_flux_errors_are_not_hidden_by_backtracking()
{
    const auto built{cfd::build_mesh(cfd::verification::make_euler_tube(1, 1))};
    // Finite admissible state, but (rhoE+p)*ux cannot be represented as double.
    const auto value{gas.to_conservative({1.0, 1.0e110, 0.0, 1.0e215})};
    const cfd::ConservativeField initial{1, value};
    cfd::AxisymmetricEulerSolver solver{built.mesh, initial, gas, cfd::verification::euler_tube_conditions()};
    require_throws_with_message<std::runtime_error>([&] { solver.step(1.0); }, "Euler flux",
                                                    "A flux failure must propagate without time-step retries.");
    require(solver.diagnostics().time == 0.0 && solver.diagnostics().accepted_steps == 0 &&
                solver.diagnostics().rejected_steps == 0,
            "Flux failure was treated as a candidate rejection.");
    require_uniform(solver.state(), value);
}

void test_receding_slip_wall_has_positive_pressure()
{
    const auto built{cfd::build_mesh(cfd::verification::make_euler_tube(1, 1))};
    const double ux{0.8 * std::sqrt(1.4)};
    const cfd::ConservativeField initial{1, gas.to_conservative({1.0, ux, 0.0, 1.0})};
    // Left wall normal=(-1,0): un=-0.8*a. Other non-axis faces extrapolate.
    cfd::AxisymmetricEulerSolver solver{
        built.mesh, initial, gas, {4, {Type::Axis, Type::Transmissive, Type::SlipWall, Type::Transmissive}}};
    const auto result{solver.step(1.0e-6)};
    // Right momentum outflow=(rho*ux^2+p)*pi; left outflow=-p_wall*pi.
    const double wall_pressure{1.0 + ux * ux - result.boundary_outflow_rate.momentum_x / std::numbers::pi};
    require(wall_pressure > 0.0, "Receding admissible flow produced negative wall pressure.");
    // Isentropic expansion to un_wall=0 gives p_wall/p=(1-0.2*0.8)^7.
    require_near(wall_pressure, 0.29509034655744, tolerance, "Wall pressure must satisfy the normal Riemann problem.");
    require_balance(result);
}

} // namespace

int main()
{
    int failures{};
    failures +=
        cfd::test::run_test("Euler uniform rest and axial flow TRI/QUAD/mixed", test_uniform_rest_and_axial_flow);
    failures +=
        cfd::test::run_test("Euler source, axis, CFL and final time", test_pressure_source_cfl_axis_and_final_time);
    failures += cfd::test::run_test("Euler internal faces and closed-wall conservation",
                                    test_internal_face_conservation_and_closed_walls);
    failures += cfd::test::run_test("Euler input and BC validation", test_input_and_boundary_validation);
    failures += cfd::test::run_test("Euler oblique slip-wall impermeability", test_oblique_slip_walls_are_impermeable);
    failures += cfd::test::run_test("Euler candidate rejection and global retry",
                                    test_candidate_backtracking_is_global_and_transactional);
    failures +=
        cfd::test::run_test("Euler persistent alternating buffers", test_state_buffers_alternate_without_copying);
    failures += cfd::test::run_test("Euler unrecoverable flux error", test_flux_errors_are_not_hidden_by_backtracking);
    failures += cfd::test::run_test("Euler receding slip-wall pressure", test_receding_slip_wall_has_positive_pressure);
    return cfd::test::finish_tests(failures, "axisymmetric Euler");
}
