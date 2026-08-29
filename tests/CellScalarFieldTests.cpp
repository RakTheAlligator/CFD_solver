#include "cfd/field/CellScalarField.hpp"

#include "support/TestUtils.hpp"

#include <span>
#include <type_traits>
#include <utility>

namespace
{

using cfd::test::require;

static_assert(!std::is_default_constructible_v<cfd::CellScalarField>);

static_assert(std::is_copy_constructible_v<cfd::CellScalarField>);
static_assert(std::is_move_constructible_v<cfd::CellScalarField>);
static_assert(!std::is_copy_assignable_v<cfd::CellScalarField>);
static_assert(!std::is_move_assignable_v<cfd::CellScalarField>);

static_assert(std::is_same_v<decltype(std::declval<cfd::CellScalarField &>()[cfd::Index{}]), double &>);
static_assert(std::is_same_v<decltype(std::declval<const cfd::CellScalarField &>()[cfd::Index{}]), const double &>);

static_assert(std::is_same_v<decltype(std::declval<cfd::CellScalarField &>().values()), std::span<double>>);
static_assert(std::is_same_v<decltype(std::declval<const cfd::CellScalarField &>().values()), std::span<const double>>);

void test_default_value_initialization()
{
    const cfd::CellScalarField field{4};

    require(field.size() == 4, "Cell scalar field contains an incorrect number of values.");

    for (const double value : field.values())
    {
        require(value == 0.0, "Cell scalar field value was not initialized to zero.");
    }
}

void test_uniform_value_initialization()
{
    const cfd::CellScalarField field{3, 2.5};

    require(field.size() == 3, "Uniformly initialized cell scalar field has an incorrect size.");

    for (const double value : field.values())
    {
        require(value == 2.5, "Cell scalar field value does not match the requested initial value.");
    }
}

void test_mutable_indexed_access()
{
    cfd::CellScalarField field{3};

    field[1] = 4.5;

    require(field[0] == 0.0, "Indexed write modified the preceding cell value.");
    require(field[1] == 4.5, "Indexed write did not modify the requested cell value.");
    require(field[2] == 0.0, "Indexed write modified the following cell value.");
}

void test_const_indexed_access()
{
    cfd::CellScalarField field{2};
    field[1] = 7.0;

    const cfd::CellScalarField &const_field{field};

    require(const_field[0] == 0.0, "Const indexed access returned an incorrect first value.");
    require(const_field[1] == 7.0, "Const indexed access returned an incorrect second value.");
}

void test_mutable_values_view()
{
    cfd::CellScalarField field{3, 1.0};
    std::span<double> values{field.values()};

    values[0] = 3.0;
    values[2] = 5.0;

    require(values.size() == field.size(), "Mutable values view has an incorrect size.");
    require(field[0] == 3.0, "Mutable values view did not modify the first field value.");
    require(field[1] == 1.0, "Mutable values view modified an unrelated field value.");
    require(field[2] == 5.0, "Mutable values view did not modify the final field value.");
}

void test_const_values_view()
{
    const cfd::CellScalarField field{3, 6.0};
    const std::span<const double> values{field.values()};

    require(values.size() == field.size(), "Const values view has an incorrect size.");

    for (const double value : values)
    {
        require(value == 6.0, "Const values view returned an incorrect field value.");
    }
}

void test_copy_construction_is_deep()
{
    cfd::CellScalarField original{3};
    original[0] = 1.0;
    original[1] = 2.0;
    original[2] = 3.0;

    cfd::CellScalarField copy{original};
    copy[1] = 8.0;

    require(copy.size() == original.size(), "Copy-constructed field has an incorrect size.");
    require(copy[0] == 1.0, "Copy-constructed field lost the first source value.");
    require(copy[1] == 8.0, "Copy-constructed field did not accept an independent modification.");
    require(copy[2] == 3.0, "Copy-constructed field lost the final source value.");
    require(original[1] == 2.0, "Modifying a copy-constructed field changed the original field.");
}

} // namespace

int main()
{
    int failure_count{};

    failure_count +=
        cfd::test::run_test("cell scalar field default value initialization", test_default_value_initialization);
    failure_count +=
        cfd::test::run_test("cell scalar field uniform value initialization", test_uniform_value_initialization);
    failure_count += cfd::test::run_test("cell scalar field mutable indexed access", test_mutable_indexed_access);
    failure_count += cfd::test::run_test("cell scalar field const indexed access", test_const_indexed_access);
    failure_count += cfd::test::run_test("cell scalar field mutable values view", test_mutable_values_view);
    failure_count += cfd::test::run_test("cell scalar field const values view", test_const_values_view);
    failure_count += cfd::test::run_test("cell scalar field copy construction", test_copy_construction_is_deep);

    return cfd::test::finish_tests(failure_count, "cell scalar field");
}
