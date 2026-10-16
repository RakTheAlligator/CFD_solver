#include "cfd/compressible/ConservativeField.hpp"
#include "support/TestUtils.hpp"

#include <array>
#include <span>
#include <type_traits>
#include <utility>

namespace
{
using cfd::test::require;

static_assert(std::is_trivially_copyable_v<cfd::ConservativeState>);
static_assert(std::is_copy_constructible_v<cfd::ConservativeField>);
static_assert(!std::is_copy_assignable_v<cfd::ConservativeField>);
static_assert(std::is_nothrow_move_constructible_v<cfd::ConservativeField>);
static_assert(!std::is_move_assignable_v<cfd::ConservativeField>);
// Every component exposes only a view, never a mutable reference to its owner.
static_assert(std::is_same_v<decltype(std::declval<cfd::ConservativeField &>().density()), std::span<double>>);
static_assert(
    std::is_same_v<decltype(std::declval<const cfd::ConservativeField &>().density()), std::span<const double>>);
static_assert(std::is_same_v<decltype(std::declval<cfd::ConservativeField &>().momentum_x()), std::span<double>>);
static_assert(
    std::is_same_v<decltype(std::declval<const cfd::ConservativeField &>().momentum_x()), std::span<const double>>);
static_assert(std::is_same_v<decltype(std::declval<cfd::ConservativeField &>().momentum_r()), std::span<double>>);
static_assert(
    std::is_same_v<decltype(std::declval<const cfd::ConservativeField &>().momentum_r()), std::span<const double>>);
static_assert(
    std::is_same_v<decltype(std::declval<cfd::ConservativeField &>().total_energy_density()), std::span<double>>);
static_assert(std::is_same_v<decltype(std::declval<const cfd::ConservativeField &>().total_energy_density()),
                             std::span<const double>>);
static_assert(
    !std::is_constructible_v<cfd::CellScalarField, decltype(std::declval<cfd::ConservativeField &>().density())>);

void test_cardinality_zero_and_uniform_initialization()
{
    const cfd::ConservativeField empty{0};
    require(empty.size() == 0 && empty.total_energy_density().size() == 0, "Empty field cardinality mismatch.");
    const cfd::ConservativeField zero{4};
    const cfd::ConservativeField uniform{4, {2.0, 3.0, -4.0, 20.0}};
    for (const auto *field : {&zero, &uniform})
    {
        require(field->size() == 4 && field->density().size() == 4 && field->momentum_x().size() == 4 &&
                    field->momentum_r().size() == 4 && field->total_energy_density().size() == 4,
                "Conservative component cardinalities differ.");
    }
    for (cfd::Index cell = 0; cell < zero.size(); ++cell)
    {
        const auto z{zero.state(cell)}, u{uniform.state(cell)};
        require(z.density == 0.0 && z.momentum_x == 0.0 && z.momentum_r == 0.0 && z.total_energy_density == 0.0,
                "Default values must be exactly zero.");
        require(u.density == 2.0 && u.momentum_x == 3.0 && u.momentum_r == -4.0 && u.total_energy_density == 20.0,
                "Uniform initialization changed a component.");
    }
}

void test_mutable_independent_contiguous_components()
{
    cfd::ConservativeField field{3};
    field.density()[0] = 1.0;
    field.momentum_x()[1] = 2.0;
    field.momentum_r()[2] = -3.0;
    field.total_energy_density()[0] = 4.0;
    const std::array<std::span<double>, 4> components{field.density(), field.momentum_x(), field.momentum_r(),
                                                      field.total_energy_density()};
    for (cfd::Index i = 0; i < components.size(); ++i)
    {
        const auto component{components.at(i)};
        require(&component[1] == component.subspan(1).data(), "Component is not contiguous.");
        for (cfd::Index j = i + 1; j < components.size(); ++j)
            require(component.data() != components.at(j).data(), "Component arrays alias.");
    }
    require(field.density()[1] == 0.0 && field.momentum_x()[0] == 0.0 && field.momentum_r()[0] == 0.0 &&
                field.total_energy_density()[1] == 0.0,
            "Mutation leaked to another component/cell.");
    field.set_state(1, {5.0, 6.0, 7.0, 8.0});
    const cfd::ConservativeField &read_only{field};
    const auto state{read_only.state(1)};
    require(state.density == 5.0 && state.momentum_x == 6.0 && state.momentum_r == 7.0 &&
                state.total_energy_density == 8.0,
            "Local state load/store mismatch.");
}

void test_rebinding_views_preserves_component_storage()
{
    cfd::ConservativeField field{3, {1.0, 2.0, 3.0, 4.0}};
    const std::array<std::span<double>, 4> components{field.density(), field.momentum_x(), field.momentum_r(),
                                                      field.total_energy_density()};
    for (const auto view : components)
    {
        auto extracted{view};
        extracted[0] = 5.0;
        extracted = extracted.subspan(1);
        require(extracted.size() == 2, "Rebinding a component view failed.");
    }
    const std::array<std::span<const double>, 4> current{field.density(), field.momentum_x(), field.momentum_r(),
                                                         field.total_energy_density()};
    require(field.size() == 3, "Rebinding a view changed field cardinality.");
    for (cfd::Index i = 0; i < current.size(); ++i)
    {
        require(current.at(i).size() == 3 && current.at(i).data() == components.at(i).data(),
                "Rebinding a view changed component ownership or cardinality.");
        require(current.at(i)[0] == 5.0, "A component view must still permit value mutation.");
    }
}

void test_deep_copy_and_storage_preserving_move()
{
    cfd::ConservativeField original{2, {1.0, 2.0, 3.0, 4.0}};
    cfd::ConservativeField copy{original};
    copy.set_state(0, {5.0, 6.0, 7.0, 8.0});
    require(original.density()[0] == 1.0 && original.momentum_x()[0] == 2.0 && original.momentum_r()[0] == 3.0 &&
                original.total_energy_density()[0] == 4.0,
            "Copy construction did not isolate all four arrays.");
    const auto *rho{copy.density().data()}, *mx{copy.momentum_x().data()}, *mr{copy.momentum_r().data()},
        *energy{copy.total_energy_density().data()};
    const cfd::ConservativeField moved{std::move(copy)};
    require(moved.size() == 2 && moved.density().data() == rho && moved.momentum_x().data() == mx &&
                moved.momentum_r().data() == mr && moved.total_energy_density().data() == energy,
            "Move construction must transfer, not copy, the four arrays.");
    require(moved.state(0).total_energy_density == 8.0, "Move construction lost values.");
}
} // namespace

int main()
{
    int failures{};
    failures += cfd::test::run_test("conservative cardinality and initialization",
                                    test_cardinality_zero_and_uniform_initialization);
    failures += cfd::test::run_test("conservative SoA access", test_mutable_independent_contiguous_components);
    failures += cfd::test::run_test("conservative view ownership", test_rebinding_views_preserves_component_storage);
    failures += cfd::test::run_test("conservative copy/move", test_deep_copy_and_storage_preserving_move);
    return cfd::test::finish_tests(failures, "conservative field");
}
