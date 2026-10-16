#include "cfd/io/EulerSteadyCsvWriter.hpp"

#include "cfd/compressible/EulerSteadyMonitor.hpp"

#include <iomanip>
#include <limits>
#include <stdexcept>

namespace cfd
{
EulerSteadyCsvWriter::EulerSteadyCsvWriter(const std::filesystem::path &path) : output_(path)
{
    output_
        << "step,time,dt,CFL,residual_rho,residual_rhoux,residual_rhour,residual_rhoE,"
           "mass_flow_in,mass_flow_out,mass_flow_relative_imbalance,minimum_density,minimum_pressure,rejected_steps,"
           "absolute_residual_rho,absolute_residual_rhoux,absolute_residual_rhour,absolute_residual_rhoE,"
           "change_rho,change_rhoux,change_rhour,change_rhoE,mass_flow_change,status\n";
    output_ << std::setprecision(std::numeric_limits<double>::max_digits10);
    output_.flush();
    if (!output_)
        throw std::runtime_error("Cannot open/write Euler steady CSV header: " + path.string());
}

void EulerSteadyCsvWriter::write(const EulerSteadyObservation &observation)
{
    const auto &time{observation.temporal};
    output_ << time.accepted_steps << ',' << time.time << ',' << time.time_step << ',' << time.cfl;
    for (const double residual : observation.normalized_residual)
        output_ << ',' << residual;
    output_ << ',' << observation.current.mass_flow_in << ',' << observation.current.mass_flow_out << ','
            << observation.mass_flow_relative_imbalance << ',' << time.minimum_density << ',' << time.minimum_pressure
            << ',' << time.rejected_steps;
    for (const double residual : observation.current.absolute_residual)
        output_ << ',' << residual;
    for (const double change : observation.field_change)
        output_ << ',' << change;
    output_ << ',' << observation.mass_flow_change << ',' << static_cast<int>(observation.status) << '\n';
    output_.flush();
    if (!output_)
        throw std::runtime_error("Cannot write/flush Euler steady CSV row.");
}
} // namespace cfd
