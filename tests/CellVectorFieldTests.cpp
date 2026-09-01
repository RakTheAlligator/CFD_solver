#include "cfd/field/CellVectorField.hpp"

#include "support/TestUtils.hpp"

#include <span>
#include <type_traits>
#include <utility>

namespace
{

using cfd::test::require;

static_assert(!std::is_default_constructible_v<cfd::CellVectorField>);

static_assert(std::is_copy_constructible_v<cfd::CellVectorField>);
static_assert(!std::is_copy_assignable_v<cfd::CellVectorField>);
static_assert(std::is_nothrow_move_constructible_v<cfd::CellVectorField>);
static_assert(!std::is_move_assignable_v<cfd::CellVectorField>);

static_assert(std::is_same_v<decltype(std::declval<cfd::CellVectorField &>()[cfd::Index{}]), cfd::Vector2 &>);
static_assert(
    std::is_same_v<decltype(std::declval<const cfd::CellVectorField &>()[cfd::Index{}]), const cfd::Vector2 &>);

static_assert(std::is_same_v<decltype(std::declval<cfd::CellVectorField &>().values()), std::span<cfd::Vector2>>);
static_assert(
    std::is_same_v<decltype(std::declval<const cfd::CellVectorField &>().values()), std::span<const cfd::Vector2>>);

void test_default_value_initialization()
{
    const cfd::CellVectorField field{4};

    require(field.size() == 4, "Cell vector field contains an incorrect number of values.");

    for (const cfd::Vector2 &value : field.values())
    {
        require(value.x == 0.0, "Cell vector field x-component was not initialized to zero.");
        require(value.y == 0.0, "Cell vector field y-component was not initialized to zero.");
    }
}

void test_uniform_value_initialization()
{
    const cfd::Vector2 initial_value{1.5, -2.5};
    const cfd::CellVectorField field{3, initial_value};

    require(field.size() == 3, "Uniformly initialized cell vector field has an incorrect size.");

    for (const cfd::Vector2 &value : field.values())
    {
        require(value.x == initial_value.x, "Cell vector field x-component does not match the initial value.");
        require(value.y == initial_value.y, "Cell vector field y-component does not match the initial value.");
    }
}

void test_mutable_indexed_access()
{
    cfd::CellVectorField field{3};

    field[1].x = 4.5;
    field[1].y = -3.5;

    require(field[0].x == 0.0 && field[0].y == 0.0, "Indexed write modified the preceding vector value.");
    require(field[1].x == 4.5, "Indexed write did not modify the requested x-component.");
    require(field[1].y == -3.5, "Indexed write did not modify the requested y-component.");
    require(field[2].x == 0.0 && field[2].y == 0.0, "Indexed write modified the following vector value.");
}

void test_const_indexed_access()
{
    cfd::CellVectorField field{2};
    field[1].x = 7.0;
    field[1].y = -8.0;

    const cfd::CellVectorField &const_field{field};

    require(const_field[0].x == 0.0 && const_field[0].y == 0.0,
            "Const indexed access returned an incorrect first vector.");
    require(const_field[1].x == 7.0, "Const indexed access returned an incorrect x-component.");
    require(const_field[1].y == -8.0, "Const indexed access returned an incorrect y-component.");
}

void test_mutable_values_view()
{
    cfd::CellVectorField field{3, {1.0, 2.0}};
    std::span<cfd::Vector2> values{field.values()};

    values[0].x = 3.0;
    values[0].y = 4.0;
    values[2].x = 5.0;
    values[2].y = 6.0;

    require(values.size() == field.size(), "Mutable vector values view has an incorrect size.");
    require(field[0].x == 3.0 && field[0].y == 4.0, "Mutable values view did not modify the first field vector.");
    require(field[1].x == 1.0 && field[1].y == 2.0, "Mutable values view modified an unrelated field vector.");
    require(field[2].x == 5.0 && field[2].y == 6.0, "Mutable values view did not modify the final field vector.");
}

void test_const_values_view()
{
    const cfd::CellVectorField field{3, {6.0, -7.0}};
    const std::span<const cfd::Vector2> values{field.values()};

    require(values.size() == field.size(), "Const vector values view has an incorrect size.");

    for (const cfd::Vector2 &value : values)
    {
        require(value.x == 6.0, "Const values view returned an incorrect x-component.");
        require(value.y == -7.0, "Const values view returned an incorrect y-component.");
    }
}

void test_copy_construction_is_deep()
{
    cfd::CellVectorField original{3};
    original[0] = {1.0, -1.0};
    original[1] = {2.0, -2.0};
    original[2] = {3.0, -3.0};

    cfd::CellVectorField copy{original};
    copy[1].x = 8.0;
    copy[1].y = -9.0;

    require(copy.size() == original.size(), "Copy-constructed cell vector field has an incorrect size.");
    require(copy[0].x == 1.0 && copy[0].y == -1.0, "Copy-constructed field lost the first vector value.");
    require(copy[1].x == 8.0 && copy[1].y == -9.0,
            "Copy-constructed field did not accept an independent modification.");
    require(copy[2].x == 3.0 && copy[2].y == -3.0, "Copy-constructed field lost the final vector value.");
    require(original[1].x == 2.0 && original[1].y == -2.0,
            "Modifying a copy-constructed field changed the original field.");
}

} // namespace

int main()
{
    int failure_count{};

    failure_count +=
        cfd::test::run_test("cell vector field default value initialization", test_default_value_initialization);
    failure_count +=
        cfd::test::run_test("cell vector field uniform value initialization", test_uniform_value_initialization);
    failure_count += cfd::test::run_test("cell vector field mutable indexed access", test_mutable_indexed_access);
    failure_count += cfd::test::run_test("cell vector field const indexed access", test_const_indexed_access);
    failure_count += cfd::test::run_test("cell vector field mutable values view", test_mutable_values_view);
    failure_count += cfd::test::run_test("cell vector field const values view", test_const_values_view);
    failure_count += cfd::test::run_test("cell vector field copy construction", test_copy_construction_is_deep);

    return cfd::test::finish_tests(failure_count, "cell vector field");
}
