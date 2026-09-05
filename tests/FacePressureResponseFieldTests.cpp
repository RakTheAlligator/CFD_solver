#include "cfd/field/FacePressureResponseField.hpp"

#include "support/TestUtils.hpp"

#include <span>
#include <type_traits>
#include <utility>

namespace
{

using cfd::test::require;

static_assert(!std::is_default_constructible_v<cfd::FacePressureResponseField>);

static_assert(std::is_copy_constructible_v<cfd::FacePressureResponseField>);
static_assert(std::is_nothrow_move_constructible_v<cfd::FacePressureResponseField>);
static_assert(!std::is_copy_assignable_v<cfd::FacePressureResponseField>);
static_assert(!std::is_move_assignable_v<cfd::FacePressureResponseField>);

static_assert(std::is_same_v<decltype(std::declval<cfd::FacePressureResponseField &>()[cfd::Index{}]), double &>);
static_assert(
    std::is_same_v<decltype(std::declval<const cfd::FacePressureResponseField &>()[cfd::Index{}]), const double &>);

static_assert(std::is_same_v<decltype(std::declval<cfd::FacePressureResponseField &>().values()), std::span<double>>);
static_assert(
    std::is_same_v<decltype(std::declval<const cfd::FacePressureResponseField &>().values()), std::span<const double>>);

void test_zero_initialization()
{
    const cfd::FacePressureResponseField field{4};

    require(field.size() == 4, "Face pressure-response field contains an incorrect number of values.");
    for (const double value : field.values())
    {
        require(value == 0.0, "Face pressure-response value was not initialized to zero.");
    }
}

void test_uniform_initialization()
{
    const cfd::FacePressureResponseField field{3, 2.5};

    require(field.size() == 3, "Uniformly initialized face pressure-response field has an incorrect size.");
    for (const double value : field.values())
    {
        require(value == 2.5, "Face pressure-response value does not match the requested initial value.");
    }
}

void test_mutable_and_const_access()
{
    cfd::FacePressureResponseField field{3};
    field[1] = 4.5;
    std::span<double> mutable_values{field.values()};
    mutable_values[2] = 6.5;

    const cfd::FacePressureResponseField &const_field{field};
    const std::span<const double> const_values{const_field.values()};

    require(field[0] == 0.0, "Indexed write modified the preceding face pressure response.");
    require(const_field[1] == 4.5, "Const indexed access returned an incorrect face pressure response.");
    require(const_values[2] == 6.5, "Const values view returned an incorrect face pressure response.");
    require(mutable_values.size() == field.size(), "Mutable face pressure-response view has an incorrect size.");
    require(const_values.size() == field.size(), "Const face pressure-response view has an incorrect size.");
}

void test_copy_construction_is_deep()
{
    cfd::FacePressureResponseField original{3};
    original[0] = 1.0;
    original[1] = 2.0;
    original[2] = 3.0;

    cfd::FacePressureResponseField copy{original};
    copy[1] = 8.0;

    require(copy.size() == original.size(), "Copy-constructed face pressure-response field has an incorrect size.");
    require(copy[0] == 1.0 && copy[1] == 8.0 && copy[2] == 3.0,
            "Copy-constructed face pressure-response field contains incorrect values.");
    require(original[1] == 2.0, "Modifying a copied face pressure-response field changed the original field.");
    require(copy.values().data() != original.values().data(),
            "Copied face pressure-response fields unexpectedly share storage.");
}

} // namespace

int main()
{
    int failure_count{};

    failure_count += cfd::test::run_test("face pressure-response zero initialization", test_zero_initialization);
    failure_count += cfd::test::run_test("face pressure-response uniform initialization", test_uniform_initialization);
    failure_count += cfd::test::run_test("face pressure-response access", test_mutable_and_const_access);
    failure_count += cfd::test::run_test("face pressure-response deep copy", test_copy_construction_is_deep);

    return cfd::test::finish_tests(failure_count, "face pressure response field");
}
