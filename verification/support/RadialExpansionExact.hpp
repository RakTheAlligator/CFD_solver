#pragma once

#include "cfd/compressible/ConservativeField.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace cfd::verification
{

struct RadialExpansionParameters
{
    double density{1.0};
    double pressure{1.0};
    double gamma{1.4};
    double alpha0{1.0};
};

struct RadialExpansionValues
{
    double density{};
    double pressure{};
    double rate{}; ///< ur/r, including its regular limit at r=0.
};

struct RadialMoments
{
    double mean_radius{};
    double mean_squared_radius{};
    double radius_variance{};
};

/// Exact smooth radial expansion and rectangular volume projection.
/// Reference arithmetic is independent of the production gas/flux/geometry APIs.
/// See docs/AxisymmetricEulerRadialExpansion.md for the Euler proof.
class RadialExpansionExact
{
  public:
    explicit RadialExpansionExact(const RadialExpansionParameters parameters = {}) : parameters_(parameters)
    {
        if (!std::isfinite(parameters.density) || !(parameters.density > 0.0) || !std::isfinite(parameters.pressure) ||
            !(parameters.pressure > 0.0) || !std::isfinite(parameters.gamma) || !(parameters.gamma > 1.0) ||
            !std::isfinite(parameters.alpha0) || parameters.alpha0 < 0.0)
            throw std::invalid_argument("Radial expansion requires finite rho0,p0>0, gamma>1, alpha0>=0.");
    }

    [[nodiscard]] double gamma() const noexcept
    {
        return parameters_.gamma;
    }

    [[nodiscard]] RadialExpansionValues values(const double time) const
    {
        if (!std::isfinite(time) || time < 0.0)
            throw std::invalid_argument("Radial expansion requires finite non-negative time.");
        const long double s{1.0L + static_cast<long double>(parameters_.alpha0) * time};
        return {positive(parameters_.density / s / s),
                positive(std::exp(std::log(static_cast<long double>(parameters_.pressure)) -
                                  2.0L * parameters_.gamma * std::log(s))),
                finite(parameters_.alpha0 / s)};
    }

    [[nodiscard]] static RadialMoments moments(const double inner, const double outer)
    {
        if (!std::isfinite(inner) || !std::isfinite(outer) || inner < 0.0 || !(outer > inner))
            throw std::invalid_argument("Radial moments require finite 0<=ri<ro.");
        const long double h{static_cast<long double>(outer) - inner}, m{inner + h / 2.0L};
        const long double ratio{h / m};
        return {positive(m + ratio * (h / 12.0L)), positive(m * m + h * h / 4.0L),
                positive(h * h / 12.0L * (1.0L - ratio * ratio / 12.0L))};
    }

    [[nodiscard]] ConservativeState cell_average(const double inner, const double outer, const double time) const
    {
        const auto state{values(time)};
        const auto moment{moments(inner, outer)};
        const long double rho{state.density}, rate{state.rate};
        return {state.density, 0.0, finite(rho * rate * moment.mean_radius),
                positive(state.pressure / (static_cast<long double>(parameters_.gamma) - 1.0L) +
                         0.5L * rho * rate * rate * moment.mean_squared_radius)};
    }

    [[nodiscard]] ConservativeState point_state(const double radius, const double time) const
    {
        if (!std::isfinite(radius) || radius < 0.0)
            throw std::invalid_argument("Radial expansion requires finite non-negative radius.");
        const auto state{values(time)};
        const long double velocity{static_cast<long double>(state.rate) * radius};
        return {state.density, 0.0, finite(state.density * velocity),
                positive(state.pressure / (static_cast<long double>(parameters_.gamma) - 1.0L) +
                         0.5L * state.density * velocity * velocity)};
    }

    [[nodiscard]] double projected_pressure(const double inner, const double outer, const double time) const
    {
        const auto state{values(time)};
        return positive(state.pressure + 0.5L * (static_cast<long double>(parameters_.gamma) - 1.0L) * state.density *
                                             state.rate * state.rate * moments(inner, outer).radius_variance);
    }

  private:
    static double finite(const long double value)
    {
        if (!std::isfinite(value) || std::abs(value) > std::numeric_limits<double>::max())
            throw std::runtime_error("Radial reference is not representable as finite double.");
        const double result{static_cast<double>(value)};
        if (value != 0.0L && result == 0.0)
            throw std::runtime_error("Radial reference underflows to zero.");
        return result;
    }

    static double positive(const long double value)
    {
        const double result{finite(value)};
        if (!(result > 0.0))
            throw std::runtime_error("Radial reference requires a representable positive quantity.");
        return result;
    }

    RadialExpansionParameters parameters_;
};

} // namespace cfd::verification
