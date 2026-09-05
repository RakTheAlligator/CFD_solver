#include "cfd/field/FaceFluxField.hpp"

#include "support/TestUtils.hpp"

#include <span>
#include <type_traits>
#include <utility>

namespace
{

using cfd::test::require;

static_assert(!std::is_default_constructible_v<cfd::FaceFluxField>);

static_assert(std::is_copy_constructible_v<cfd::FaceFluxField>);
static_assert(std::is_nothrow_move_constructible_v<cfd::FaceFluxField>);
static_assert(!std::is_copy_assignable_v<cfd::FaceFluxField>);
static_assert(!std::is_move_assignable_v<cfd::FaceFluxField>);

static_assert(std::is_same_v<decltype(std::declval<cfd::FaceFluxField &>()[cfd::Index{}]), double &>);
static_assert(std::is_same_v<decltype(std::declval<const cfd::FaceFluxField &>()[cfd::Index{}]), const double &>);

static_assert(std::is_same_v<decltype(std::declval<cfd::FaceFluxField &>().values()), std::span<double>>);
static_assert(std::is_same_v<decltype(std::declval<const cfd::FaceFluxField &>().values()), std::span<const double>>);

void test_zero_initialization()
{
    const cfd::FaceFluxField field{4};

    require(field.size() == 4, "Face flux field contains an incorrect number of values.");

    for (const double value : field.values())
    {
        require(value == 0.0, "Face flux field value was not initialized to zero.");
    }
}

void test_uniform_initialization()
{
    const cfd::FaceFluxField field{3, -2.5};

    require(field.size() == 3, "Uniformly initialized face flux field has an incorrect size.");

    for (const double value : field.values())
    {
        require(value == -2.5, "Face flux field value does not match the requested initial value.");
    }
}

void test_mutable_and_const_indexed_access()
{
    cfd::FaceFluxField field{3};
    field[1] = 4.5;
    const cfd::FaceFluxField &const_field{field};

    require(field[0] == 0.0, "Indexed write modified the preceding face flux.");
    require(const_field[1] == 4.5, "Indexed access did not preserve the requested face flux.");
    require(field[2] == 0.0, "Indexed write modified the following face flux.");
}

void test_mutable_and_const_values_views()
{
    cfd::FaceFluxField field{3, 1.0};
    std::span<double> mutable_values{field.values()};
    mutable_values[0] = -3.0;
    mutable_values[2] = 5.0;

    const cfd::FaceFluxField &const_field{field};
    const std::span<const double> const_values{const_field.values()};

    require(mutable_values.size() == field.size(), "Mutable face flux view has an incorrect size.");
    require(const_values.size() == field.size(), "Const face flux view has an incorrect size.");
    require(const_values[0] == -3.0, "Mutable face flux view did not modify the first value.");
    require(const_values[1] == 1.0, "Mutable face flux view modified an unrelated value.");
    require(const_values[2] == 5.0, "Const face flux view returned an incorrect final value.");
}

void test_copy_construction_is_deep()
{
    cfd::FaceFluxField original{3};
    original[0] = -1.0;
    original[1] = 2.0;
    original[2] = -3.0;

    cfd::FaceFluxField copy{original};
    copy[1] = 8.0;

    require(copy.size() == original.size(), "Copy-constructed face flux field has an incorrect size.");
    require(copy[0] == -1.0, "Copy-constructed face flux field lost the first source value.");
    require(copy[1] == 8.0, "Copy-constructed face flux field did not accept an independent modification.");
    require(copy[2] == -3.0, "Copy-constructed face flux field lost the final source value.");
    require(original[1] == 2.0, "Modifying a copied face flux field changed the original field.");
}

} // namespace

int main()
{
    int failure_count{};

    failure_count += cfd::test::run_test("face flux field zero initialization", test_zero_initialization);
    failure_count += cfd::test::run_test("face flux field uniform initialization", test_uniform_initialization);
    failure_count += cfd::test::run_test("face flux field indexed access", test_mutable_and_const_indexed_access);
    failure_count += cfd::test::run_test("face flux field values views", test_mutable_and_const_values_views);
    failure_count += cfd::test::run_test("face flux field copy construction", test_copy_construction_is_deep);

    return cfd::test::finish_tests(failure_count, "face flux field");
}
