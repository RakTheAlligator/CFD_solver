#include "app/BoundaryBasedInitialization.hpp"

#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVelocityField.hpp"
#include "cfd/field/FaceFluxField.hpp"
#include "cfd/field/PressureCorrectionBoundaryConditions.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"

#include "support/MeshFixtures.hpp"
#include "support/TestUtils.hpp"

#include <array>
#include <limits>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace
{
using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws_with_message;

constexpr cfd::BoundaryId bottom{0};
constexpr cfd::BoundaryId right{1};
constexpr cfd::BoundaryId top{2};
constexpr cfd::BoundaryId left{3};
constexpr cfd::Index boundary_count{4};

cfd::Mesh make_mesh()
{
    auto raw{cfd::test::make_two_triangle_raw_mesh()};
    for (auto &point : raw.nodes)
    {
        point.x *= 2.0;
    }
    raw.boundary_groups = {{bottom, "bottom"}, {right, "right"}, {top, "top"}, {left, "left"}};
    raw.boundary_edges = {{{0, 1}, bottom}, {{1, 2}, right}, {{2, 3}, top}, {{3, 0}, left}};
    return cfd::build_mesh(std::move(raw)).mesh;
}

struct BoundaryData
{
    cfd::ScalarBoundaryConditions u;
    cfd::ScalarBoundaryConditions v;
    cfd::ScalarBoundaryConditions p;
    cfd::PressureCorrectionBoundaryConditions correction;
};

BoundaryData conditions(const std::array<cfd::Vector2, boundary_count> &velocity = {},
                        const std::array<std::optional<double>, boundary_count> &pressure = {})
{
    std::vector<cfd::ScalarBoundaryCondition> u;
    std::vector<cfd::ScalarBoundaryCondition> v;
    std::vector<cfd::ScalarBoundaryCondition> p;
    std::vector<cfd::PressureCorrectionBoundaryConditionType> correction;
    for (cfd::BoundaryId boundary = 0; boundary < boundary_count; ++boundary)
    {
        u.emplace_back(cfd::ScalarBoundaryConditionType::Dirichlet, velocity.at(boundary).x);
        v.emplace_back(cfd::ScalarBoundaryConditionType::Dirichlet, velocity.at(boundary).y);
        const bool fixed{pressure.at(boundary).has_value()};
        p.emplace_back(fixed ? cfd::ScalarBoundaryConditionType::Dirichlet : cfd::ScalarBoundaryConditionType::Neumann,
                       pressure.at(boundary).value_or(0.0));
        correction.push_back(fixed ? cfd::PressureCorrectionBoundaryConditionType::FixedPressure
                                   : cfd::PressureCorrectionBoundaryConditionType::FixedMassFlux);
    }
    return {{boundary_count, std::move(u)},
            {boundary_count, std::move(v)},
            {boundary_count, std::move(p)},
            {boundary_count, std::move(correction)}};
}

void require_state(const cfd::CellVelocityField &velocity, const cfd::CellScalarField &pressure,
                   cfd::Vector2 expected_velocity, double expected_pressure)
{
    for (cfd::Index cell = 0; cell < velocity.size(); ++cell)
    {
        require_near(velocity.u()[cell], expected_velocity.x, 1.0e-14, "Wrong bulk u.");
        require_near(velocity.v()[cell], expected_velocity.y, 1.0e-14, "Wrong bulk v.");
        require_near(pressure[cell], expected_pressure, 1.0e-14, "Wrong bulk pressure.");
    }
}

void test_uniform_inlet_and_no_flux_changes()
{
    const auto mesh{make_mesh()};
    std::array<cfd::Vector2, boundary_count> values{};
    values.at(left) = {2.0, 0.5};
    const auto boundary{conditions(values)};
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    const cfd::FaceFluxField flux{mesh.face_count(), 7.0};
    cfd::app::initialize_from_boundary_conditions(mesh, boundary.u, boundary.v, boundary.p, boundary.correction,
                                                  {17.0, 19.0}, 23.0, velocity, pressure);
    require_state(velocity, pressure, {2.0, 0.5}, 23.0);
    for (double value : flux.values())
    {
        require(value == 7.0, "Boundary initialization changed a face flux.");
    }
}

void test_weighted_inflows_ignore_outflow()
{
    const auto mesh{make_mesh()};
    std::array<cfd::Vector2, boundary_count> values{};
    values.at(left) = {2.0, 0.0};
    values.at(bottom) = {0.0, 3.0};
    values.at(right) = {8.0, 0.0};
    values.at(top) = {0.0, 9.0};
    const auto boundary{conditions(values)};
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    cfd::app::initialize_from_boundary_conditions(mesh, boundary.u, boundary.v, boundary.p, boundary.correction,
                                                  {17.0, 19.0}, 23.0, velocity, pressure);
    // Left length is 1, bottom length is 2; positive outward fluxes contribute nothing.
    require_state(velocity, pressure, {2.0 / 3.0, 2.0}, 23.0);
}

void test_no_inflow_preserves_fallback()
{
    const auto mesh{make_mesh()};
    std::array<cfd::Vector2, boundary_count> values{};
    values.at(left) = {-2.0, 0.0};
    values.at(bottom) = {0.0, -3.0};
    values.at(right) = {8.0, 0.0};
    values.at(top) = {0.0, 9.0};
    const auto boundary{conditions(values)};
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    cfd::app::initialize_from_boundary_conditions(mesh, boundary.u, boundary.v, boundary.p, boundary.correction,
                                                  {17.0, 19.0}, 23.0, velocity, pressure);
    require_state(velocity, pressure, {17.0, 19.0}, 23.0);
}

void test_only_two_dirichlet_components_define_inflow()
{
    const auto mesh{make_mesh()};
    const auto boundary{conditions()};
    std::vector<cfd::ScalarBoundaryCondition> u(boundary_count, {cfd::ScalarBoundaryConditionType::Dirichlet, 0.0});
    u.at(left) = {cfd::ScalarBoundaryConditionType::Neumann, 100.0};
    const cfd::ScalarBoundaryConditions partial{boundary_count, std::move(u)};
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    cfd::app::initialize_from_boundary_conditions(mesh, partial, boundary.v, boundary.p, boundary.correction,
                                                  {17.0, 19.0}, 23.0, velocity, pressure);
    require_state(velocity, pressure, {17.0, 19.0}, 23.0);
}

void test_scale_aware_tangential_flux_guard()
{
    const auto mesh{make_mesh()};
    std::array<cfd::Vector2, boundary_count> values{};
    values.at(right) = {-std::numeric_limits<double>::epsilon(), 1.0};
    auto boundary{conditions(values)};
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    cfd::app::initialize_from_boundary_conditions(mesh, boundary.u, boundary.v, boundary.p, boundary.correction,
                                                  {17.0, 19.0}, 23.0, velocity, pressure);
    require_state(velocity, pressure, {17.0, 19.0}, 23.0);

    values = {};
    values.at(left) = {1.0e-20, 0.0};
    const auto tiny_inflow{conditions(values)};
    cfd::app::initialize_from_boundary_conditions(mesh, tiny_inflow.u, tiny_inflow.v, tiny_inflow.p,
                                                  tiny_inflow.correction, {17.0, 19.0}, 23.0, velocity, pressure);
    require_near(velocity.u()[0], 1.0e-20, 1.0e-30, "An absolute inflow floor discarded a genuine small inflow.");
}

void test_dirichlet_pressure_averages()
{
    const auto mesh{make_mesh()};
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    std::array<std::optional<double>, boundary_count> values{};
    values.at(left) = 10.0;
    const auto unique{conditions({}, values)};
    cfd::app::initialize_from_boundary_conditions(mesh, unique.u, unique.v, unique.p, unique.correction, {17.0, 19.0},
                                                  23.0, velocity, pressure);
    require_state(velocity, pressure, {17.0, 19.0}, 10.0);
    values.at(left) = 2.0;
    values.at(bottom) = 8.0;
    const auto multiple{conditions({}, values)};
    cfd::app::initialize_from_boundary_conditions(mesh, multiple.u, multiple.v, multiple.p, multiple.correction,
                                                  {17.0, 19.0}, 23.0, velocity, pressure);
    require_state(velocity, pressure, {17.0, 19.0}, 6.0);
}

void test_nonfinite_values_rejected_before_writes()
{
    const auto mesh{make_mesh()};
    const auto boundary{conditions()};
    cfd::CellVelocityField velocity{mesh.cell_count(), {17.0, 19.0}};
    cfd::CellScalarField pressure{mesh.cell_count(), 23.0};
    for (double invalid : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
                           -std::numeric_limits<double>::infinity()})
    {
        // A vector move preserves references to its elements. Retain a test-only
        // borrow to emulate invalid BC data after collection construction.
        std::vector<cfd::ScalarBoundaryCondition> values(boundary_count,
                                                         {cfd::ScalarBoundaryConditionType::Dirichlet, 0.0});
        const std::span<cfd::ScalarBoundaryCondition> borrow{values};
        const cfd::ScalarBoundaryConditions corrupted{boundary_count, std::move(values)};
        borrow[left].value = invalid;
        for (const bool corrupt_pressure : {false, true})
        {
            const cfd::PressureCorrectionBoundaryConditions fixed_pressure{
                boundary_count, std::vector<cfd::PressureCorrectionBoundaryConditionType>(
                                    boundary_count, cfd::PressureCorrectionBoundaryConditionType::FixedPressure)};
            require_throws_with_message<std::invalid_argument>(
                [&]() {
                    cfd::app::initialize_from_boundary_conditions(
                        mesh, corrupt_pressure ? boundary.u : corrupted, boundary.v,
                        corrupt_pressure ? corrupted : boundary.p,
                        corrupt_pressure ? fixed_pressure : boundary.correction, {}, 0.0, velocity, pressure);
                },
                "finite boundary", "Accepted a non-finite boundary value.");
        }
        require_throws_with_message<std::invalid_argument>(
            [&]() {
                cfd::app::initialize_from_boundary_conditions(mesh, boundary.u, boundary.v, boundary.p,
                                                              boundary.correction, {invalid, 0.0}, 0.0, velocity,
                                                              pressure);
            },
            "finite internalField", "Accepted a non-finite fallback.");
        require_state(velocity, pressure, {17.0, 19.0}, 23.0);
    }
}

void test_cardinalities_and_pressure_pairing()
{
    const auto mesh{make_mesh()};
    const auto boundary{conditions()};
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    const cfd::ScalarBoundaryConditions empty{0, {}};
    require_throws_with_message<std::invalid_argument>(
        [&]() {
            cfd::app::initialize_from_boundary_conditions(mesh, empty, boundary.v, boundary.p, boundary.correction, {},
                                                          0.0, velocity, pressure);
        },
        "cardinalities", "Accepted an incompatible boundary collection.");
    cfd::CellVelocityField short_velocity{1};
    require_throws_with_message<std::invalid_argument>(
        [&]() {
            cfd::app::initialize_from_boundary_conditions(mesh, boundary.u, boundary.v, boundary.p, boundary.correction,
                                                          {}, 0.0, short_velocity, pressure);
        },
        "cardinalities", "Accepted incompatible outputs.");
    const cfd::PressureCorrectionBoundaryConditions mismatch{
        boundary_count, std::vector<cfd::PressureCorrectionBoundaryConditionType>(
                            boundary_count, cfd::PressureCorrectionBoundaryConditionType::FixedPressure)};
    require_throws_with_message<std::invalid_argument>(
        [&]() {
            cfd::app::initialize_from_boundary_conditions(mesh, boundary.u, boundary.v, boundary.p, mismatch, {}, 0.0,
                                                          velocity, pressure);
        },
        "Dirichlet pressure", "Accepted inconsistent pressure conditions.");
}
} // namespace

int main()
{
    int failures{};
    failures += cfd::test::run_test("uniform inlet and unchanged flux", test_uniform_inlet_and_no_flux_changes);
    failures += cfd::test::run_test("weighted inflow excludes outflow", test_weighted_inflows_ignore_outflow);
    failures += cfd::test::run_test("fallback without inflow", test_no_inflow_preserves_fallback);
    failures += cfd::test::run_test("require both Dirichlet velocity components",
                                    test_only_two_dirichlet_components_define_inflow);
    failures += cfd::test::run_test("scale-aware tangential guard", test_scale_aware_tangential_flux_guard);
    failures += cfd::test::run_test("Dirichlet pressure averages", test_dirichlet_pressure_averages);
    failures += cfd::test::run_test("reject non-finite inputs", test_nonfinite_values_rejected_before_writes);
    failures += cfd::test::run_test("cardinalities and pressure pairing", test_cardinalities_and_pressure_pairing);
    return cfd::test::finish_tests(failures, "BoundaryBasedInitialization");
}
