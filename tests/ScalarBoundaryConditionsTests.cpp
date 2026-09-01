#include "cfd/field/ScalarBoundaryConditions.hpp"

#include "support/TestUtils.hpp"

#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace
{

using cfd::test::require;
using cfd::test::require_throws;

static_assert(!std::is_default_constructible_v<cfd::ScalarBoundaryCondition>);
static_assert(!std::is_default_constructible_v<cfd::ScalarBoundaryConditions>);

static_assert(std::is_copy_constructible_v<cfd::ScalarBoundaryConditions>);
static_assert(std::is_move_constructible_v<cfd::ScalarBoundaryConditions>);
static_assert(!std::is_copy_assignable_v<cfd::ScalarBoundaryConditions>);
static_assert(!std::is_move_assignable_v<cfd::ScalarBoundaryConditions>);

static_assert(std::is_same_v<decltype(std::declval<const cfd::ScalarBoundaryConditions &>()[cfd::BoundaryId{}]),
                             const cfd::ScalarBoundaryCondition &>);

void test_valid_construction()
{
    using cfd::ScalarBoundaryCondition;
    using cfd::ScalarBoundaryConditionType;

    std::vector<ScalarBoundaryCondition> conditions{
        {ScalarBoundaryConditionType::Dirichlet, 1.0},
        {ScalarBoundaryConditionType::Dirichlet, 0.0},
        {ScalarBoundaryConditionType::Neumann, 0.0},
    };

    const cfd::ScalarBoundaryConditions boundary_conditions{3, std::move(conditions)};

    require(boundary_conditions.size() == 3, "Scalar boundary conditions contain an incorrect number of values.");
    require(boundary_conditions[0].type == ScalarBoundaryConditionType::Dirichlet,
            "Boundary 0 has an incorrect scalar condition type.");
    require(boundary_conditions[0].value == 1.0, "Boundary 0 has an incorrect scalar condition value.");
    require(boundary_conditions[1].type == ScalarBoundaryConditionType::Dirichlet,
            "Boundary 1 has an incorrect scalar condition type.");
    require(boundary_conditions[1].value == 0.0, "Boundary 1 has an incorrect scalar condition value.");
    require(boundary_conditions[2].type == ScalarBoundaryConditionType::Neumann,
            "Boundary 2 has an incorrect scalar condition type.");
    require(boundary_conditions[2].value == 0.0, "Boundary 2 has an incorrect scalar condition value.");
}

void test_accepts_negative_finite_value()
{
    using cfd::ScalarBoundaryCondition;
    using cfd::ScalarBoundaryConditionType;

    std::vector<ScalarBoundaryCondition> conditions{
        {ScalarBoundaryConditionType::Neumann, -2.5},
    };

    const cfd::ScalarBoundaryConditions boundary_conditions{1, std::move(conditions)};

    require(boundary_conditions[0].type == ScalarBoundaryConditionType::Neumann,
            "Negative-value boundary condition has an incorrect type.");
    require(boundary_conditions[0].value == -2.5, "A finite negative boundary condition value was not preserved.");
}

void test_rejects_incorrect_condition_count()
{
    require_throws<std::invalid_argument>(
        []() {
            std::vector<cfd::ScalarBoundaryCondition> conditions{
                {cfd::ScalarBoundaryConditionType::Dirichlet, 1.0},
            };

            static_cast<void>(cfd::ScalarBoundaryConditions{2, std::move(conditions)});
        },
        "Scalar boundary conditions accepted a condition count different from the boundary count.");
}

void test_rejects_nan_value()
{
    require_throws<std::invalid_argument>(
        []() {
            std::vector<cfd::ScalarBoundaryCondition> conditions{
                {cfd::ScalarBoundaryConditionType::Dirichlet, std::numeric_limits<double>::quiet_NaN()},
            };

            static_cast<void>(cfd::ScalarBoundaryConditions{1, std::move(conditions)});
        },
        "Scalar boundary conditions accepted a NaN value.");
}

void test_rejects_positive_infinity()
{
    require_throws<std::invalid_argument>(
        []() {
            std::vector<cfd::ScalarBoundaryCondition> conditions{
                {cfd::ScalarBoundaryConditionType::Neumann, std::numeric_limits<double>::infinity()},
            };

            static_cast<void>(cfd::ScalarBoundaryConditions{1, std::move(conditions)});
        },
        "Scalar boundary conditions accepted positive infinity.");
}

void test_rejects_invalid_condition_type()
{
    require_throws<std::invalid_argument>(
        []() {
            // NOLINTNEXTLINE(clang-analyzer-optin.core.EnumCastOutOfRange)
            const auto invalid_type{static_cast<cfd::ScalarBoundaryConditionType>(255)};

            std::vector<cfd::ScalarBoundaryCondition> conditions{
                {invalid_type, 0.0},
            };

            static_cast<void>(cfd::ScalarBoundaryConditions{1, std::move(conditions)});
        },
        "Scalar boundary conditions accepted an unsupported condition type.");
}

void test_copy_construction_preserves_conditions()
{
    using cfd::ScalarBoundaryCondition;
    using cfd::ScalarBoundaryConditionType;

    std::vector<ScalarBoundaryCondition> conditions{
        {ScalarBoundaryConditionType::Dirichlet, 4.0},
        {ScalarBoundaryConditionType::Neumann, -1.5},
    };

    const cfd::ScalarBoundaryConditions original{2, std::move(conditions)};

    // This copy is intentional: the test verifies copy construction.
    const cfd::ScalarBoundaryConditions copy{original}; // NOLINT(performance-unnecessary-copy-initialization)

    require(copy.size() == 2, "Copy-constructed boundary conditions have an incorrect size.");
    require(copy[0].type == ScalarBoundaryConditionType::Dirichlet,
            "Copy-constructed boundary conditions lost the first condition type.");
    require(copy[0].value == 4.0, "Copy-constructed boundary conditions lost the first condition value.");
    require(copy[1].type == ScalarBoundaryConditionType::Neumann,
            "Copy-constructed boundary conditions lost the second condition type.");
    require(copy[1].value == -1.5, "Copy-constructed boundary conditions lost the second condition value.");
}

} // namespace

int main()
{
    int failure_count{};

    failure_count += cfd::test::run_test("valid scalar boundary conditions", test_valid_construction);
    failure_count +=
        cfd::test::run_test("negative scalar boundary condition value", test_accepts_negative_finite_value);
    failure_count +=
        cfd::test::run_test("reject incorrect scalar boundary condition count", test_rejects_incorrect_condition_count);
    failure_count += cfd::test::run_test("reject NaN scalar boundary condition value", test_rejects_nan_value);
    failure_count +=
        cfd::test::run_test("reject infinite scalar boundary condition value", test_rejects_positive_infinity);
    failure_count +=
        cfd::test::run_test("reject unsupported scalar boundary condition type", test_rejects_invalid_condition_type);
    failure_count += cfd::test::run_test("scalar boundary conditions copy construction",
                                         test_copy_construction_preserves_conditions);

    return cfd::test::finish_tests(failure_count, "scalar boundary conditions");
}
