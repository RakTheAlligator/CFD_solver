#include "cfd/numerics/ScalarTransportAssembler.hpp"

#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVectorField.hpp"
#include "cfd/linear_algebra/ScalarLinearSystem.hpp"
#include "cfd/mesh/MeshBuilder.hpp"

#include "support/ScalarTransportFixtures.hpp"
#include "support/TestUtils.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <span>
#include <type_traits>
#include <vector>

namespace
{

using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws;

static_assert(!std::is_copy_constructible_v<cfd::ScalarTransportAssembler>);
static_assert(!std::is_copy_assignable_v<cfd::ScalarTransportAssembler>);
static_assert(std::is_nothrow_move_constructible_v<cfd::ScalarTransportAssembler>);
static_assert(!std::is_move_assignable_v<cfd::ScalarTransportAssembler>);

struct SchemeCase
{
    cfd::ScalarConvectionScheme scheme;
    cfd::ScalarConvectionLimiter limiter;
};

constexpr std::array cases{
    SchemeCase{cfd::ScalarConvectionScheme::FirstOrderUpwind, cfd::ScalarConvectionLimiter::None},
    SchemeCase{cfd::ScalarConvectionScheme::Linear, cfd::ScalarConvectionLimiter::None},
    SchemeCase{cfd::ScalarConvectionScheme::Hybrid, cfd::ScalarConvectionLimiter::None},
    SchemeCase{cfd::ScalarConvectionScheme::LinearUpwind, cfd::ScalarConvectionLimiter::None},
    SchemeCase{cfd::ScalarConvectionScheme::LinearUpwind, cfd::ScalarConvectionLimiter::BarthJespersen}};

void require_values(const std::span<const double> actual, const std::span<const double> expected)
{
    require(actual.size() == expected.size(), "Transport array cardinalities differ.");
    for (cfd::Index index = 0; index < actual.size(); ++index)
    {
        require_near(actual[index], expected[index], 1.0e-12, "Transport composition differs from existing operators.");
    }
}

void require_matrix(const cfd::ScalarLinearSystem &actual, const cfd::ScalarLinearSystem &expected)
{
    require_values(actual.diagonal(), expected.diagonal());
    require_values(actual.owner_neighbor_coefficients(), expected.owner_neighbor_coefficients());
    require_values(actual.neighbor_owner_coefficients(), expected.neighbor_owner_coefficients());
}

void seed_system(cfd::ScalarLinearSystem &system)
{
    std::ranges::fill(system.diagonal(), 11.0);
    std::ranges::fill(system.owner_neighbor_coefficients(), 12.0);
    std::ranges::fill(system.neighbor_owner_coefficients(), 13.0);
    std::ranges::fill(system.rhs(), 14.0);
}

void require_seeded(const cfd::ScalarLinearSystem &system)
{
    for (const double value : system.diagonal())
    {
        require_near(value, 11.0, 0.0, "Invalid input changed diagonal.");
    }
    for (const double value : system.owner_neighbor_coefficients())
    {
        require_near(value, 12.0, 0.0, "Invalid input changed owner-neighbor coefficients.");
    }
    for (const double value : system.neighbor_owner_coefficients())
    {
        require_near(value, 13.0, 0.0, "Invalid input changed neighbor-owner coefficients.");
    }
    for (const double value : system.rhs())
    {
        require_near(value, 14.0, 0.0, "Invalid input changed RHS.");
    }
}

void test_composition_overwrite_and_input_immutability()
{
    for (const bool sheared : {false, true})
    {
        auto build{cfd::build_mesh(cfd::test::make_transport_raw_mesh(sheared))};
        const cfd::Mesh &mesh{build.mesh};
        const auto conditions{cfd::test::transport_affine_conditions(mesh)};
        cfd::CellScalarField field{mesh.cell_count()};
        field[0] = 0.2;
        field[1] = 0.8;
        const cfd::CellVectorField gradient{mesh.cell_count(), {3.0, -0.75}};
        for (const double speed : {0.0, 1.0, 10.0})
        {
            const auto flux{cfd::test::transport_flux(mesh, speed)};
            const std::vector<double> saved_flux{flux.values().begin(), flux.values().end()};
            for (const auto &study : cases)
            {
                cfd::ScalarTransportAssembler assembler{mesh, 0.25, study.scheme, study.limiter};
                const cfd::ScalarDiffusionOperator diffusion{mesh, 0.25};
                const cfd::ScalarConvectionOperator convection{mesh, study.scheme, study.limiter};
                cfd::ScalarLinearSystem expected{mesh};
                cfd::ScalarLinearSystem actual{mesh};
                std::vector<double> limiter(
                    study.limiter == cfd::ScalarConvectionLimiter::BarthJespersen ? mesh.cell_count() : 0);
                diffusion.add_matrix_contributions(conditions, expected);
                convection.add_matrix_contributions(conditions, flux, expected, diffusion.face_primary_coefficients());
                diffusion.add_boundary_rhs(conditions, expected.rhs());
                convection.add_boundary_rhs(conditions, flux, expected.rhs(), diffusion.face_primary_coefficients());
                convection.add_deferred_correction_rhs(field, conditions, gradient, flux, limiter, expected.rhs());
                diffusion.add_non_orthogonal_rhs(conditions, gradient, expected.rhs());

                seed_system(actual);
                assembler.assemble_matrix(conditions, flux, actual);
                require_matrix(actual, expected);
                for (const double value : actual.rhs())
                {
                    require_near(value, 14.0, 0.0, "Matrix assembly changed RHS.");
                }
                assembler.assemble_rhs(field, gradient, conditions, flux, actual);
                require_matrix(actual, expected);
                require_values(actual.rhs(), expected.rhs());
                assembler.assemble_matrix(conditions, flux, actual);
                require_values(actual.rhs(), expected.rhs());
                assembler.assemble_rhs(field, gradient, conditions, flux, actual);
                require_matrix(actual, expected);
                require_values(actual.rhs(), expected.rhs());
                require_values(flux.values(), saved_flux);
                require_near(field[0], 0.2, 0.0, "Assembly changed previous field.");
                require_near(field[1], 0.8, 0.0, "Assembly changed previous field.");
                for (const auto &value : gradient.values())
                {
                    require_near(value.x, 3.0, 0.0, "Assembly changed gradient x.");
                    require_near(value.y, -0.75, 0.0, "Assembly changed gradient y.");
                }
            }
        }
    }
}

void test_internal_conservation_and_direct_balance()
{
    auto build{cfd::build_mesh(cfd::test::make_transport_raw_mesh(true))};
    const cfd::Mesh &mesh{build.mesh};
    const cfd::ScalarBoundaryConditions conditions{
        mesh.boundary_groups().size(),
        std::vector<cfd::ScalarBoundaryCondition>(mesh.boundary_groups().size(),
                                                  {cfd::ScalarBoundaryConditionType::Neumann, 0.0})};
    cfd::FaceFluxField flux{mesh.face_count()};
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (!mesh.face_adjacencies()[face_id].is_boundary())
        {
            flux[face_id] = 3.0;
        }
    }
    cfd::CellScalarField field{mesh.cell_count()};
    field[0] = 0.2;
    field[1] = 0.8;
    const cfd::CellVectorField gradient{mesh.cell_count(), {0.8, -0.2}};
    for (const auto &study : cases)
    {
        cfd::ScalarTransportAssembler assembler{mesh, 0.25, study.scheme, study.limiter};
        cfd::ScalarLinearSystem system{mesh};
        assembler.assemble_matrix(conditions, flux, system);
        assembler.assemble_rhs(field, gradient, conditions, flux, system);
        std::vector<double> product(mesh.cell_count());
        system.apply_matrix(field.values(), product);
        const cfd::ScalarConvectionOperator convection{mesh, study.scheme, study.limiter};
        const cfd::ScalarDiffusionOperator diffusion{mesh, 0.25};
        cfd::CellScalarField convective{mesh.cell_count()};
        cfd::CellScalarField diffusive{mesh.cell_count()};
        std::vector<double> limiter(mesh.cell_count());
        if (study.scheme == cfd::ScalarConvectionScheme::LinearUpwind)
        {
            convection.compute_flux_balance(field, conditions, flux, gradient, limiter, convective);
        }
        else
        {
            convection.compute_flux_balance(field, conditions, flux, convective, diffusion.face_primary_coefficients());
        }
        diffusion.compute_flux_balance(field, conditions, gradient, diffusive);
        double total{};
        for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
        {
            const double residual{product[cell_id] - system.rhs()[cell_id]};
            require_near(residual, convective[cell_id] + diffusive[cell_id], 1.0e-12,
                         "Transport residual differs from direct outward flux balance.");
            total += residual;
        }
        require_near(total, 0.0, 1.0e-12, "Internal transport fluxes are not conservative.");
    }
}

void test_invalid_inputs_before_clear()
{
    auto build{cfd::build_mesh(cfd::test::make_transport_raw_mesh())};
    auto other_build{cfd::build_mesh(cfd::test::make_transport_raw_mesh())};
    const cfd::Mesh &mesh{build.mesh};
    const auto conditions{cfd::test::transport_channel_conditions()};
    const cfd::ScalarBoundaryConditions wrong_conditions{0, {}};
    cfd::FaceFluxField flux{cfd::test::transport_flux(mesh)};
    const cfd::FaceFluxField wrong_flux{mesh.face_count() + 1};
    cfd::CellScalarField field{mesh.cell_count()};
    const cfd::CellScalarField wrong_field{mesh.cell_count() + 1};
    cfd::CellVectorField gradient{mesh.cell_count()};
    const cfd::CellVectorField wrong_gradient{mesh.cell_count() + 1};
    cfd::ScalarLinearSystem system{mesh};
    cfd::ScalarLinearSystem other_system{other_build.mesh};
    seed_system(system);
    seed_system(other_system);
    cfd::ScalarTransportAssembler assembler{mesh, 0.25};
    require_throws<std::invalid_argument>([&]() { assembler.assemble_matrix(wrong_conditions, flux, system); },
                                          "Wrong BC cardinality was accepted.");
    require_throws<std::invalid_argument>([&]() { assembler.assemble_matrix(conditions, wrong_flux, system); },
                                          "Wrong flux cardinality was accepted.");
    require_throws<std::invalid_argument>([&]() { assembler.assemble_matrix(conditions, flux, other_system); },
                                          "Other Mesh identity was accepted.");
    require_throws<std::invalid_argument>(
        [&]() { assembler.assemble_rhs(field, gradient, conditions, flux, other_system); },
        "RHS assembly accepted another Mesh.");
    require_throws<std::invalid_argument>(
        [&]() { assembler.assemble_rhs(wrong_field, gradient, conditions, flux, system); },
        "Wrong field cardinality was accepted.");
    require_throws<std::invalid_argument>(
        [&]() { assembler.assemble_rhs(field, wrong_gradient, conditions, flux, system); },
        "Wrong gradient cardinality was accepted.");
    require_throws<std::invalid_argument>(
        [&]() { assembler.assemble_rhs(field, gradient, wrong_conditions, flux, system); },
        "RHS accepted wrong BC cardinality.");
    require_throws<std::invalid_argument>(
        [&]() { assembler.assemble_rhs(field, gradient, conditions, wrong_flux, system); },
        "RHS accepted wrong flux cardinality.");
    for (const double invalid : {0.0, -1.0, std::numeric_limits<double>::quiet_NaN(),
                                 std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()})
    {
        require_throws<std::invalid_argument>(
            [&]() { const cfd::ScalarTransportAssembler invalid_assembler{mesh, invalid}; },
            "Invalid diffusion coefficient was accepted.");
    }
    require_throws<std::invalid_argument>(
        [&]() {
            const cfd::ScalarTransportAssembler invalid_assembler{mesh, 0.25, cfd::ScalarConvectionScheme::Linear,
                                                                  cfd::ScalarConvectionLimiter::BarthJespersen};
        },
        "Unsupported scheme/limiter combination was accepted.");
    for (const double invalid : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
                                 -std::numeric_limits<double>::infinity()})
    {
        flux[mesh.face_count() - 1] = invalid;
        require_throws<std::invalid_argument>([&]() { assembler.assemble_matrix(conditions, flux, system); },
                                              "Non-finite flux was accepted.");
        require_throws<std::invalid_argument>(
            [&]() { assembler.assemble_rhs(field, gradient, conditions, flux, system); },
            "RHS accepted non-finite flux.");
        flux[mesh.face_count() - 1] = 0.0;
        field[1] = invalid;
        require_throws<std::invalid_argument>(
            [&]() { assembler.assemble_rhs(field, gradient, conditions, flux, system); },
            "Non-finite field was accepted.");
        field[1] = 0.0;
        gradient[1].y = invalid;
        require_throws<std::invalid_argument>(
            [&]() { assembler.assemble_rhs(field, gradient, conditions, flux, system); },
            "Non-finite gradient was accepted.");
        gradient[1].y = 0.0;
    }
    require_seeded(system);
    require_seeded(other_system);
}

} // namespace

int main()
{
    int failures{};
    failures += cfd::test::run_test("Transport composition, overwrite and immutability",
                                    test_composition_overwrite_and_input_immutability);
    failures +=
        cfd::test::run_test("Transport conservation and direct balance", test_internal_conservation_and_direct_balance);
    failures += cfd::test::run_test("Transport input validation", test_invalid_inputs_before_clear);
    return failures == 0 ? 0 : 1;
}
