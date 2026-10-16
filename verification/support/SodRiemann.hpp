#pragma once

#include <cmath>
#include <stdexcept>

namespace cfd::verification
{

struct SodPrimitive
{
    double density{};
    double velocity{};
    double pressure{};
};

/// Independent exact classical Sod solution: gamma=1.4, left=(1,0,1),
/// right=(0.125,0,0.1), discontinuity x=0.5. No numerical flux or gas API.
/// Rarefaction uses isentropic Riemann invariants; shock uses Rankine-Hugoniot.
/// Wave-curve construction: https://www.clawpack.org/riemann_book/html/Euler.html
class SodRiemann
{
  public:
    SodRiemann()
    {
        // The unique star pressure lies between the two initial pressures.
        double lower{0.1}, upper{1.0};
        for (int iteration = 0; iteration < 100; ++iteration)
        {
            const double middle{0.5 * (lower + upper)};
            if (left_curve(middle) + right_curve(middle) > 0.0)
                upper = middle;
            else
                lower = middle;
        }
        pressure_ = 0.5 * (lower + upper);
        velocity_ = 0.5 * (right_curve(pressure_) - left_curve(pressure_));
        left_density_ = std::pow(pressure_, 1.0 / gamma);
        constexpr double beta{(gamma - 1.0) / (gamma + 1.0)};
        const double ratio{pressure_ / 0.1};
        right_density_ = 0.125 * (ratio + beta) / (beta * ratio + 1.0);
        shock_speed_ = std::sqrt(gamma * 0.1 / 0.125) *
                       std::sqrt((gamma + 1.0) / (2.0 * gamma) * ratio + (gamma - 1.0) / (2.0 * gamma));
    }

    [[nodiscard]] SodPrimitive sample(const double x, const double time) const
    {
        if (!std::isfinite(x) || !std::isfinite(time) || time < 0.0)
            throw std::invalid_argument("Sod reference requires finite x and non-negative time.");
        if (time == 0.0)
            return x < 0.5 ? SodPrimitive{1.0, 0.0, 1.0} : SodPrimitive{0.125, 0.0, 0.1};
        const double similarity{(x - 0.5) / time};
        const double sound_left{std::sqrt(gamma)};
        const double sound_star_left{sound_left * std::pow(pressure_, (gamma - 1.0) / (2.0 * gamma))};
        if (similarity < -sound_left)
            return {1.0, 0.0, 1.0};
        if (similarity < velocity_ - sound_star_left)
        {
            const double velocity{2.0 / (gamma + 1.0) * (sound_left + similarity)};
            const double sound{2.0 / (gamma + 1.0) * (sound_left - 0.5 * (gamma - 1.0) * similarity)};
            return {std::pow(sound / sound_left, 2.0 / (gamma - 1.0)), velocity,
                    std::pow(sound / sound_left, 2.0 * gamma / (gamma - 1.0))};
        }
        if (similarity < velocity_)
            return {left_density_, velocity_, pressure_};
        if (similarity < shock_speed_)
            return {right_density_, velocity_, pressure_};
        return {0.125, 0.0, 0.1};
    }

    [[nodiscard]] double shock_speed() const noexcept
    {
        return shock_speed_;
    }

  private:
    static constexpr double gamma{1.4};

    static double left_curve(const double pressure)
    {
        return 2.0 * std::sqrt(gamma) / (gamma - 1.0) * (std::pow(pressure, (gamma - 1.0) / (2.0 * gamma)) - 1.0);
    }

    static double right_curve(const double pressure)
    {
        constexpr double coefficient{2.0 / ((gamma + 1.0) * 0.125)};
        constexpr double offset{(gamma - 1.0) / (gamma + 1.0) * 0.1};
        return (pressure - 0.1) * std::sqrt(coefficient / (pressure + offset));
    }

    double pressure_{};
    double velocity_{};
    double left_density_{};
    double right_density_{};
    double shock_speed_{};
};

} // namespace cfd::verification
