#include "cfd/numerics/SteadyScalarTransportSolver.hpp"

#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/numerics/LeastSquaresGradient.hpp"

#include "support/ScalarTransportFixtures.hpp"
#include "support/TestUtils.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

namespace
{

using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws;

static_assert(!std::is_copy_constructible_v<cfd::SteadyScalarTransportSolver>);
static_assert(!std::is_copy_assignable_v<cfd::SteadyScalarTransportSolver>);
static_assert(!std::is_move_constructible_v<cfd::SteadyScalarTransportSolver>);
static_assert(!std::is_move_assignable_v<cfd::SteadyScalarTransportSolver>);

[[nodiscard]]
cfd::SteadyScalarTransportOptions test_options()
{
    cfd::SteadyScalarTransportOptions options;
    options.field_relative_tolerance = 1.0e-10;
    options.discrete_relative_tolerance = 1.0e-10;
    options.linear_solver.relative_tolerance = 1.0e-13;
    return options;
}

void require_converged(const cfd::SteadyScalarTransportResult &result, const cfd::SteadyScalarTransportOptions &options)
{
    require(result.converged && result.linear_solve.converged, "Transport solve did not converge.");
    require(result.field_relative_change <= options.field_relative_tolerance, "Field-change tolerance was not met.");
    require(result.discrete_relative_residual <= options.discrete_relative_tolerance,
            "Full discrete equation tolerance was not met.");
}

void test_affine_diffusion_and_converged_initial_guess()
{
    for (const bool sheared : {false, true})
    {
        auto build{cfd::build_mesh(cfd::test::make_transport_raw_mesh(sheared))};
        const cfd::Mesh &mesh{build.mesh};
        const auto conditions{cfd::test::transport_affine_conditions(mesh)};
        const cfd::FaceFluxField flux{mesh.face_count()};
        const auto options{test_options()};
        cfd::SteadyScalarTransportSolver solver{mesh, 0.25, cfd::ScalarConvectionScheme::FirstOrderUpwind,
                                                cfd::ScalarConvectionLimiter::None, options};
        cfd::CellScalarField field{mesh.cell_count()};
        const auto result{solver.solve(flux, conditions, field)};
        require_converged(result, options);
        require(result.iteration_count > 1,
                "A non-exact initial field unexpectedly met the change criterion immediately.");
        for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
        {
            require_near(field[cell_id], cfd::test::transport_affine_value(mesh.cell_centers()[cell_id]), 1.0e-9,
                         "Affine diffusion solution is incorrect.");
            field[cell_id] = cfd::test::transport_affine_value(mesh.cell_centers()[cell_id]);
        }
        const auto exact_result{solver.solve(flux, conditions, field)};
        require_converged(exact_result, options);
        require(exact_result.iteration_count == 1, "Exact initial field required extra outer iterations.");
    }
}

void test_zero_flow_zero_field_and_neumann_data()
{
    auto build{cfd::build_mesh(cfd::test::make_transport_raw_mesh())};
    const cfd::Mesh &mesh{build.mesh};
    std::vector<cfd::ScalarBoundaryCondition> boundary_data(mesh.boundary_groups().size(),
                                                            {cfd::ScalarBoundaryConditionType::Neumann, 0.0});
    boundary_data.front().type = cfd::ScalarBoundaryConditionType::Dirichlet;
    const cfd::ScalarBoundaryConditions conditions{mesh.boundary_groups().size(), std::move(boundary_data)};
    const cfd::FaceFluxField flux{mesh.face_count()};
    cfd::CellScalarField field{mesh.cell_count()};
    cfd::SteadyScalarTransportSolver solver{mesh, 0.25};
    // A zero Dirichlet boundary fixes the constant while retaining 0/0 diagnostics.
    const auto result{solver.solve(flux, conditions, field)};
    require(result.converged && result.iteration_count == 1, "Zero state did not remain a fixed point.");
    require_near(result.field_relative_change, 0.0, 0.0, "Zero field change was not normalized to zero.");
    require_near(result.discrete_relative_residual, 0.0, 0.0, "Zero discrete residual was not normalized to zero.");
}

void test_explicit_iteration_and_distinct_residuals()
{
    auto build{cfd::build_mesh(cfd::test::make_transport_raw_mesh())};
    const cfd::Mesh &mesh{build.mesh};
    const auto conditions{cfd::test::transport_channel_conditions()};
    const auto flux{cfd::test::transport_flux(mesh)};
    auto options{test_options()};
    options.maximum_iterations = 1;
    cfd::SteadyScalarTransportSolver one_step{mesh, 0.25, cfd::ScalarConvectionScheme::LinearUpwind,
                                              cfd::ScalarConvectionLimiter::None, options};
    cfd::CellScalarField field{mesh.cell_count()};
    const auto result{one_step.solve(flux, conditions, field)};
    require(!result.converged && result.iteration_count == 1, "Outer iteration limit was not reported normally.");
    require(result.linear_solve.converged && result.linear_solve.estimated_relative_error < 1.0e-12,
            "The frozen linear equation was not solved accurately.");
    require(result.discrete_relative_residual > options.discrete_relative_tolerance,
            "Explicit corrections did not distinguish full from inner residual.");

    cfd::ScalarTransportAssembler assembler{mesh, 0.25, cfd::ScalarConvectionScheme::LinearUpwind};
    cfd::ScalarLinearSystem system{mesh};
    cfd::CellVectorField gradient{mesh.cell_count()};
    cfd::compute_least_squares_gradient(mesh, field, conditions, gradient);
    assembler.assemble_matrix(conditions, flux, system);
    assembler.assemble_rhs(field, gradient, conditions, flux, system);
    std::vector<double> product(mesh.cell_count());
    system.apply_matrix(field.values(), product);
    double residual_norm{};
    double product_norm{};
    double rhs_norm{};
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        residual_norm = std::hypot(residual_norm, product[cell_id] - system.rhs()[cell_id]);
        product_norm = std::hypot(product_norm, product[cell_id]);
        rhs_norm = std::hypot(rhs_norm, system.rhs()[cell_id]);
    }
    require_near(result.discrete_relative_residual, residual_norm / std::max(product_norm, rhs_norm), 1.0e-12,
                 "Full residual was not rebuilt at the new field.");

    options.maximum_iterations = 1000;
    cfd::SteadyScalarTransportSolver complete{mesh, 0.25, cfd::ScalarConvectionScheme::LinearUpwind,
                                              cfd::ScalarConvectionLimiter::None, options};
    std::ranges::fill(field.values(), 0.0);
    const auto completed{complete.solve(flux, conditions, field)};
    require_converged(completed, options);
    require(completed.iteration_count > 2, "Deferred correction was not closed through multiple outer iterations.");
}

void test_field_change_is_also_required()
{
    auto build{cfd::build_mesh(cfd::test::make_transport_raw_mesh())};
    const cfd::Mesh &mesh{build.mesh};
    const auto conditions{cfd::test::transport_channel_conditions()};
    // Avoid cancellation of the right Dirichlet diffusion/convection terms
    // on this very coarse centered stencil, so the solution actually changes.
    const auto flux{cfd::test::transport_flux(mesh, 0.5)};
    auto options{test_options()};
    options.maximum_iterations = 1;
    cfd::CellScalarField field{mesh.cell_count()};
    cfd::SteadyScalarTransportSolver solver{mesh, 0.25, cfd::ScalarConvectionScheme::Linear,
                                            cfd::ScalarConvectionLimiter::None, options};
    const auto result{solver.solve(flux, conditions, field)};
    require(result.discrete_relative_residual <= options.discrete_relative_tolerance,
            "Implicit Cartesian problem has an unexpected full residual.");
    require(result.field_relative_change > options.field_relative_tolerance && !result.converged,
            "A large field change was accepted solely from a small equation residual.");
}

void test_linear_upwind_with_and_without_limiter()
{
    auto build{cfd::build_mesh(cfd::test::make_transport_raw_mesh())};
    const cfd::Mesh &mesh{build.mesh};
    const auto conditions{cfd::test::transport_channel_conditions()};
    const auto flux{cfd::test::transport_flux(mesh)};
    const std::vector<double> original_flux{flux.values().begin(), flux.values().end()};
    const auto options{test_options()};
    for (const auto limiter : {cfd::ScalarConvectionLimiter::None, cfd::ScalarConvectionLimiter::BarthJespersen})
    {
        cfd::CellScalarField field{mesh.cell_count()};
        cfd::SteadyScalarTransportSolver solver{mesh, 0.25, cfd::ScalarConvectionScheme::LinearUpwind, limiter,
                                                options};
        require_converged(solver.solve(flux, conditions, field), options);
        for (const double value : field.values())
        {
            require(std::isfinite(value), "LinearUpwind produced a non-finite field.");
            if (limiter == cfd::ScalarConvectionLimiter::BarthJespersen)
            {
                require(value >= 0.0 && value <= 1.0, "Limited channel fixture escaped its prescribed range.");
            }
        }
        for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
        {
            require_near(flux[face_id], original_flux[face_id], 0.0, "Solver changed caller face fluxes.");
        }
    }
}

void test_invalid_options_and_inputs()
{
    auto build{cfd::build_mesh(cfd::test::make_transport_raw_mesh())};
    const cfd::Mesh &mesh{build.mesh};
    const auto reject_options = [&](const cfd::SteadyScalarTransportOptions options) {
        require_throws<std::invalid_argument>(
            [&]() {
                const cfd::SteadyScalarTransportSolver solver{mesh, 0.25, cfd::ScalarConvectionScheme::FirstOrderUpwind,
                                                              cfd::ScalarConvectionLimiter::None, options};
            },
            "Invalid transport options were accepted.");
    };
    auto options{test_options()};
    options.maximum_iterations = 0;
    reject_options(options);
    for (const double invalid :
         {0.0, -1.0, 1.0, 2.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()})
    {
        options = test_options();
        options.field_relative_tolerance = invalid;
        reject_options(options);
        options = test_options();
        options.discrete_relative_tolerance = invalid;
        reject_options(options);
        options = test_options();
        options.linear_solver.relative_tolerance = invalid;
        reject_options(options);
    }
    options = test_options();
    options.linear_solver.maximum_iterations = 0;
    reject_options(options);

    cfd::SteadyScalarTransportSolver solver{mesh, 0.25};
    const auto conditions{cfd::test::transport_channel_conditions()};
    const cfd::ScalarBoundaryConditions wrong_conditions{0, {}};
    cfd::FaceFluxField flux{cfd::test::transport_flux(mesh)};
    const cfd::FaceFluxField wrong_flux{mesh.face_count() + 1};
    cfd::CellScalarField field{mesh.cell_count()};
    cfd::CellScalarField wrong_field{mesh.cell_count() + 1};
    require_throws<std::invalid_argument>([&]() { (void)solver.solve(flux, conditions, wrong_field); },
                                          "Solver accepted wrong field cardinality.");
    require_throws<std::invalid_argument>([&]() { (void)solver.solve(wrong_flux, conditions, field); },
                                          "Solver accepted wrong flux cardinality.");
    require_throws<std::invalid_argument>([&]() { (void)solver.solve(flux, wrong_conditions, field); },
                                          "Solver accepted wrong BC cardinality.");
    for (const double invalid : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
                                 -std::numeric_limits<double>::infinity()})
    {
        field[1] = invalid;
        require_throws<std::invalid_argument>([&]() { (void)solver.solve(flux, conditions, field); },
                                              "Solver accepted non-finite field.");
        field[1] = 0.0;
        flux[mesh.face_count() - 1] = invalid;
        require_throws<std::invalid_argument>([&]() { (void)solver.solve(flux, conditions, field); },
                                              "Solver accepted non-finite flux.");
        flux[mesh.face_count() - 1] = 0.0;
    }
}

} // namespace

int main()
{
    int failures{};
    failures += cfd::test::run_test("Transport affine diffusion and exact initial guess",
                                    test_affine_diffusion_and_converged_initial_guess);
    failures +=
        cfd::test::run_test("Transport zero state and Neumann data", test_zero_flow_zero_field_and_neumann_data);
    failures += cfd::test::run_test("Transport explicit iterations and residual distinction",
                                    test_explicit_iteration_and_distinct_residuals);
    failures +=
        cfd::test::run_test("Transport convergence also requires field change", test_field_change_is_also_required);
    failures +=
        cfd::test::run_test("Transport LinearUpwind limiter variants", test_linear_upwind_with_and_without_limiter);
    failures += cfd::test::run_test("Transport solver validation", test_invalid_options_and_inputs);
    return failures == 0 ? 0 : 1;
}
