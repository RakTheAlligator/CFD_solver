#include "cfd/io/SimpleConvergenceCsvWriter.hpp"

#include "cfd/numerics/IncompressibleSimpleSolver.hpp"

#include "support/TestUtils.hpp"

#include <filesystem>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

namespace
{

using cfd::test::read_text_file;
using cfd::test::require;

constexpr auto expected_header =
    "iteration,continuity,x_velocity,y_velocity,velocity_change,corrected_continuity,u_linear_residual,"
    "u_linear_iterations,v_linear_residual,v_linear_iterations,pressure_correction_linear_residual,"
    "pressure_correction_linear_iterations,maximum_pressure_correction";

class TemporaryDirectory
{
  public:
    TemporaryDirectory()
    {
        constexpr std::size_t maximum_attempt_count{64};

        const std::filesystem::path temporary_root{std::filesystem::temp_directory_path()};
        std::random_device random_source;
        std::uniform_int_distribution<unsigned long long> token_distribution;
        std::error_code last_error;

        for (std::size_t attempt = 0; attempt < maximum_attempt_count; ++attempt)
        {
            const std::filesystem::path candidate{temporary_root / ("cfd_simple_convergence_csv_writer_tests_" +
                                                                    std::to_string(token_distribution(random_source)))};
            std::error_code creation_error;
            if (std::filesystem::create_directory(candidate, creation_error))
            {
                path_ = candidate;
                return;
            }
            last_error = creation_error;
        }

        std::string message{"Unable to create a unique temporary directory for SIMPLE convergence CSV tests."};
        if (last_error)
        {
            message += " Last filesystem error: " + last_error.message();
        }
        throw std::runtime_error(message);
    }

    ~TemporaryDirectory() noexcept
    {
        std::error_code cleanup_error;
        std::filesystem::remove_all(path_, cleanup_error);
    }

    TemporaryDirectory(const TemporaryDirectory &) = delete;
    TemporaryDirectory &operator=(const TemporaryDirectory &) = delete;
    TemporaryDirectory(TemporaryDirectory &&) = delete;
    TemporaryDirectory &operator=(TemporaryDirectory &&) = delete;

    [[nodiscard]]
    const std::filesystem::path &path() const noexcept
    {
        return path_;
    }

  private:
    std::filesystem::path path_;
};

[[nodiscard]]
std::vector<std::string> lines(const std::string &text)
{
    std::istringstream input{text};
    std::vector<std::string> result;
    for (std::string line; std::getline(input, line);)
    {
        result.push_back(std::move(line));
    }
    return result;
}

void test_writes_and_flushes_stable_csv_rows()
{
    const TemporaryDirectory temporary_directory;
    const std::filesystem::path file_path{temporary_directory.path() / "convergence.csv"};
    const cfd::SimpleIterationInfo first_info{
        .iteration = 7,
        .u_solve = {true, 3, 0.015625},
        .v_solve = {true, 4, 0.0078125},
        .pressure_correction_solve = {true, 5, 0.00390625},
        .x_velocity_equation_residual = 0.25,
        .y_velocity_equation_residual = 0.5,
        .velocity_relative_change = 0.0625,
        .provisional_continuity_relative_residual = 0.125,
        .corrected_continuity_relative_residual = 0.03125,
        .maximum_pressure_correction = 2.5,
    };
    cfd::SimpleIterationInfo second_info{first_info};
    second_info.iteration = 8;
    second_info.provisional_continuity_relative_residual = 0.0625;

    {
        cfd::SimpleConvergenceCsvWriter writer{file_path};
        writer.write(first_info);

        const std::vector<std::string> first_lines{lines(read_text_file(file_path))};
        require(first_lines.size() == 2, "SIMPLE convergence CSV did not flush its first completed row.");
        require(first_lines[0] == expected_header, "SIMPLE convergence CSV header is not stable.");
        require(first_lines[1] == "7,0.125,0.25,0.5,0.0625,0.03125,0.015625,3,0.0078125,4,0.00390625,5,2.5",
                "SIMPLE convergence CSV first row contains incorrect values.");

        writer.write(second_info);
        const std::vector<std::string> second_lines{lines(read_text_file(file_path))};
        require(second_lines.size() == 3, "SIMPLE convergence CSV did not append its second row.");
        require(second_lines[0] == expected_header, "SIMPLE convergence CSV rewrote or duplicated its header.");
        require(second_lines[2] == "8,0.0625,0.25,0.5,0.0625,0.03125,0.015625,3,0.0078125,4,0.00390625,5,2.5",
                "SIMPLE convergence CSV second row contains incorrect values.");
    }

    require(std::filesystem::exists(file_path), "SIMPLE convergence CSV file was not created.");
}

} // namespace

int main()
{
    int failure_count{};

    failure_count += cfd::test::run_test("SIMPLE convergence CSV rows", test_writes_and_flushes_stable_csv_rows);

    return cfd::test::finish_tests(failure_count, "SIMPLE convergence CSV writer");
}
