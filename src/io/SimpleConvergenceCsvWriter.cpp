#include "cfd/io/SimpleConvergenceCsvWriter.hpp"

#include "cfd/numerics/IncompressibleSimpleSolver.hpp"

#include <iomanip>
#include <limits>
#include <stdexcept>

namespace cfd
{

namespace
{

constexpr auto csv_header =
    "iteration,continuity,x_velocity,y_velocity,velocity_change,corrected_continuity,u_linear_residual,"
    "u_linear_iterations,v_linear_residual,v_linear_iterations,pressure_correction_linear_residual,"
    "pressure_correction_linear_iterations,maximum_pressure_correction";

} // namespace

SimpleConvergenceCsvWriter::SimpleConvergenceCsvWriter(const std::filesystem::path &file_path)
    : output_(file_path, std::ios::out | std::ios::trunc)
{
    if (!output_.is_open())
    {
        throw std::runtime_error("Unable to open SIMPLE convergence CSV file: " + file_path.string());
    }

    output_ << csv_header << '\n';
    if (!output_)
    {
        throw std::runtime_error("Unable to write SIMPLE convergence CSV header.");
    }
    output_.flush();
    if (!output_)
    {
        throw std::runtime_error("Unable to flush SIMPLE convergence CSV header.");
    }
    output_ << std::setprecision(std::numeric_limits<double>::max_digits10);
}

void SimpleConvergenceCsvWriter::write(const SimpleIterationInfo &info)
{
    output_ << info.iteration << ',' << info.provisional_continuity_relative_residual << ','
            << info.x_velocity_equation_residual << ',' << info.y_velocity_equation_residual << ','
            << info.velocity_relative_change << ',' << info.corrected_continuity_relative_residual << ','
            << info.u_solve.estimated_relative_error << ',' << info.u_solve.iteration_count << ','
            << info.v_solve.estimated_relative_error << ',' << info.v_solve.iteration_count << ','
            << info.pressure_correction_solve.estimated_relative_error << ','
            << info.pressure_correction_solve.iteration_count << ',' << info.maximum_pressure_correction << '\n';
    if (!output_)
    {
        throw std::runtime_error("Unable to write SIMPLE convergence CSV row.");
    }

    output_.flush();
    if (!output_)
    {
        throw std::runtime_error("Unable to flush SIMPLE convergence CSV row.");
    }
}

} // namespace cfd
