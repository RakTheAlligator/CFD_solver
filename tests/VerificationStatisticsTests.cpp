#include "support/VerificationStatistics.hpp"

#include "support/TestUtils.hpp"

#include <cmath>

namespace
{

using cfd::test::require;
using cfd::test::require_near;
using cfd::test::test_tolerance;

void test_accumulates_signed_errors()
{
    cfd::verification::ErrorAccumulator accumulator;
    accumulator.add(1.0, 2.0);
    accumulator.add(3.0, -5.0);

    const cfd::verification::ErrorStatistics statistics{accumulator.finish()};

    require(statistics.cell_count == 2, "Error statistics contain an incorrect cell count.");
    require_near(statistics.total_area, 4.0, test_tolerance, "Error statistics contain an incorrect total area.");
    require_near(statistics.area_weighted_squared_error, 79.0, test_tolerance,
                 "Error statistics contain an incorrect area-weighted squared error.");
    require_near(statistics.area_weighted_rms_error, std::sqrt(79.0 / 4.0), test_tolerance,
                 "Error statistics contain an incorrect area-weighted RMS error.");
    require_near(statistics.linf_error, 5.0, test_tolerance,
                 "Error statistics do not use the absolute error for the L-infinity norm.");
}

} // namespace

int main()
{
    int failure_count{};

    failure_count += cfd::test::run_test("verification statistics signed errors", test_accumulates_signed_errors);

    return cfd::test::finish_tests(failure_count, "verification statistics");
}
