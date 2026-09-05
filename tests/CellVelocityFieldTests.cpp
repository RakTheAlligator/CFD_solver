#include "cfd/field/CellVelocityField.hpp"

#include "support/TestUtils.hpp"

#include <type_traits>
#include <utility>

namespace
{

using cfd::test::require;

static_assert(!std::is_default_constructible_v<cfd::CellVelocityField>);
static_assert(std::is_copy_constructible_v<cfd::CellVelocityField>);
static_assert(!std::is_copy_assignable_v<cfd::CellVelocityField>);
static_assert(std::is_nothrow_move_constructible_v<cfd::CellVelocityField>);
static_assert(!std::is_move_assignable_v<cfd::CellVelocityField>);

static_assert(std::is_same_v<decltype(std::declval<cfd::CellVelocityField &>().u()), cfd::CellScalarField &>);
static_assert(
    std::is_same_v<decltype(std::declval<const cfd::CellVelocityField &>().u()), const cfd::CellScalarField &>);
static_assert(std::is_same_v<decltype(std::declval<cfd::CellVelocityField &>().v()), cfd::CellScalarField &>);
static_assert(
    std::is_same_v<decltype(std::declval<const cfd::CellVelocityField &>().v()), const cfd::CellScalarField &>);

void test_zero_initialization()
{
    const cfd::CellVelocityField velocity{4};

    for (cfd::Index cell_id = 0; cell_id < velocity.size(); ++cell_id)
    {
        require(velocity.u()[cell_id] == 0.0, "Velocity u-component was not initialized to zero.");
        require(velocity.v()[cell_id] == 0.0, "Velocity v-component was not initialized to zero.");
    }
}

void test_uniform_initialization()
{
    const cfd::CellVelocityField velocity{3, {2.5, -1.5}};

    for (cfd::Index cell_id = 0; cell_id < velocity.size(); ++cell_id)
    {
        require(velocity.u()[cell_id] == 2.5, "Velocity u-component does not match the initial value.");
        require(velocity.v()[cell_id] == -1.5, "Velocity v-component does not match the initial value.");
    }
}

void test_mutable_component_access_is_independent()
{
    cfd::CellVelocityField velocity{3};
    velocity.u()[1] = 4.5;
    velocity.v()[2] = -3.5;

    require(velocity.u()[0] == 0.0 && velocity.u()[1] == 4.5 && velocity.u()[2] == 0.0,
            "Mutable u-component access changed an unrelated u value.");
    require(velocity.v()[0] == 0.0 && velocity.v()[1] == 0.0 && velocity.v()[2] == -3.5,
            "Mutable v-component access changed an unrelated v value.");
}

void test_const_component_access()
{
    cfd::CellVelocityField velocity{2};
    velocity.u()[1] = 7.0;
    velocity.v()[1] = -8.0;
    const cfd::CellVelocityField &const_velocity{velocity};

    require(const_velocity.u()[0] == 0.0 && const_velocity.u()[1] == 7.0,
            "Const u-component access returned incorrect values.");
    require(const_velocity.v()[0] == 0.0 && const_velocity.v()[1] == -8.0,
            "Const v-component access returned incorrect values.");
}

void test_component_cardinalities_agree()
{
    const cfd::CellVelocityField velocity{5};

    require(velocity.size() == 5, "Velocity field has an incorrect cardinality.");
    require(velocity.u().size() == velocity.size(), "Velocity u-component cardinality does not match the field.");
    require(velocity.v().size() == velocity.size(), "Velocity v-component cardinality does not match the field.");
}

void test_copy_construction_is_deep()
{
    cfd::CellVelocityField original{3};
    original.u()[0] = 1.0;
    original.u()[1] = 2.0;
    original.u()[2] = 3.0;
    original.v()[0] = -1.0;
    original.v()[1] = -2.0;
    original.v()[2] = -3.0;

    cfd::CellVelocityField copy{original};
    copy.u()[1] = 8.0;
    copy.v()[2] = -9.0;

    require(copy.size() == original.size(), "Copy-constructed velocity field has an incorrect cardinality.");
    require(copy.u()[0] == 1.0 && copy.u()[1] == 8.0 && copy.u()[2] == 3.0,
            "Copy-constructed velocity field has incorrect u-component values.");
    require(copy.v()[0] == -1.0 && copy.v()[1] == -2.0 && copy.v()[2] == -9.0,
            "Copy-constructed velocity field has incorrect v-component values.");
    require(original.u()[1] == 2.0, "Modifying a copied u-component changed the original velocity field.");
    require(original.v()[2] == -3.0, "Modifying a copied v-component changed the original velocity field.");
}

} // namespace

int main()
{
    int failure_count{};

    failure_count += cfd::test::run_test("cell velocity field zero initialization", test_zero_initialization);
    failure_count += cfd::test::run_test("cell velocity field uniform initialization", test_uniform_initialization);
    failure_count += cfd::test::run_test("cell velocity field mutable component access",
                                         test_mutable_component_access_is_independent);
    failure_count += cfd::test::run_test("cell velocity field const component access", test_const_component_access);
    failure_count +=
        cfd::test::run_test("cell velocity field component cardinalities", test_component_cardinalities_agree);
    failure_count += cfd::test::run_test("cell velocity field copy construction", test_copy_construction_is_deep);

    return cfd::test::finish_tests(failure_count, "cell velocity field");
}
