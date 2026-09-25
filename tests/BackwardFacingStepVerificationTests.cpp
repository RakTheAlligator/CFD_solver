#include "support/BackwardFacingStepReattachment.hpp"
#include "support/GridConvergenceIndex.hpp"

#include "support/TestUtils.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>

namespace
{

using cfd::test::require;
using cfd::test::require_near;

void test_parabolic_inlet_face_averages_preserve_flow()
{
    constexpr double test_step_height{0.4851485148514851};
    constexpr double test_inlet_height{0.5148514851485149};
    constexpr double test_mean_velocity{1.0};
    constexpr std::array eta_boundaries{0.0, 0.07, 0.24, 0.51, 0.83, 1.0};

    double integrated_flow{};
    for (std::size_t face_id = 0; face_id + 1 < eta_boundaries.size(); ++face_id)
    {
        const double lower_y{test_step_height + test_inlet_height * eta_boundaries.at(face_id)};
        const double upper_y{test_step_height + test_inlet_height * eta_boundaries.at(face_id + 1)};
        const double average{cfd::verification::parabolic_inlet_face_average(lower_y, upper_y, test_step_height,
                                                                             test_inlet_height, test_mean_velocity)};
        integrated_flow += average * (upper_y - lower_y);
    }

    require_near(integrated_flow, test_mean_velocity * test_inlet_height, 1.0e-13,
                 "Analytical inlet face averages do not preserve the prescribed volumetric flow.");
}

void test_celik_first_table_example()
{
    constexpr double fine_cell_count{18'000.0};
    constexpr double medium_cell_count{8'000.0};
    constexpr double coarse_cell_count{4'500.0};
    const double r21{std::sqrt(fine_cell_count / medium_cell_count)};
    const double r32{std::sqrt(medium_cell_count / coarse_cell_count)};

    const std::optional<cfd::verification::GridConvergenceIndexResult> result{
        cfd::verification::compute_grid_convergence_index(6.063, 5.972, 5.863, r21, r32)};

    if (!result)
    {
        throw std::runtime_error("Celik's first table example did not produce a GCI result.");
    }
    const cfd::verification::GridConvergenceIndexResult &values{*result};
    require(values.behavior == cfd::verification::GridConvergenceBehavior::Monotonic,
            "Celik's first table example was not classified as monotonic.");
    require_near(values.apparent_order, 1.53, 0.02, "Celik apparent order differs from the published value.");
    require_near(values.extrapolated_value, 6.1685, 0.002,
                 "Celik extrapolated value differs from the published value.");
    require_near(values.approximate_relative_error, 0.015, 0.001,
                 "Celik approximate relative error differs from the published value.");
    require_near(values.extrapolated_relative_error, 0.017, 0.001,
                 "Celik extrapolated relative error differs from the published value.");
    require_near(values.fine_grid_convergence_index, 0.022, 0.001,
                 "Celik fine-grid GCI differs from the published value.");
}

void test_gci_reports_nearly_identical_values_as_unusable()
{
    const std::optional<cfd::verification::GridConvergenceIndexResult> result{
        cfd::verification::compute_grid_convergence_index(1.0, 1.0, 0.9, 2.0, 2.0)};

    require(!result.has_value(), "GCI accepted a zero fine-to-medium solution difference.");
}

void test_estimates_durable_primary_reattachment()
{
    constexpr std::array samples{
        cfd::verification::NearWallVelocitySample{1.0, -0.3},
        cfd::verification::NearWallVelocitySample{2.0, -0.1},
        cfd::verification::NearWallVelocitySample{3.0, 0.2},
        cfd::verification::NearWallVelocitySample{4.0, 0.3},
    };

    require_near(cfd::verification::estimate_primary_reattachment_x(samples), 7.0 / 3.0, 1.0e-14,
                 "Primary reattachment interpolation is incorrect.");
}

void test_rejects_missing_primary_recirculation()
{
    constexpr std::array samples{
        cfd::verification::NearWallVelocitySample{1.0, 0.1},
        cfd::verification::NearWallVelocitySample{2.0, 0.2},
    };

    bool threw{};
    try
    {
        static_cast<void>(cfd::verification::estimate_primary_reattachment_x(samples));
    }
    catch (const std::runtime_error &)
    {
        threw = true;
    }
    require(threw, "Reattachment estimation accepted samples without recirculation.");
}

} // namespace

int main()
{
    int failure_count{};

    failure_count += cfd::test::run_test("parabolic inlet face averages preserve flow",
                                         test_parabolic_inlet_face_averages_preserve_flow);
    failure_count += cfd::test::run_test("Celik first table GCI example", test_celik_first_table_example);
    failure_count += cfd::test::run_test("GCI zero difference", test_gci_reports_nearly_identical_values_as_unusable);
    failure_count += cfd::test::run_test("durable primary reattachment", test_estimates_durable_primary_reattachment);
    failure_count += cfd::test::run_test("missing primary recirculation", test_rejects_missing_primary_recirculation);

    return cfd::test::finish_tests(failure_count, "backward-facing-step verification support");
}
