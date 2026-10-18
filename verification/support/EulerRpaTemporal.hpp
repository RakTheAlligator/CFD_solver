#pragma once

#include "EulerRpaExploratory.hpp"
#include "cfd/compressible/EulerCharacteristicBoundaryFlux.hpp"
#include "cfd/compressible/EulerFlux.hpp"
#include "cfd/geometry/AxisymmetricGeometry.hpp"

#include <algorithm>
#include <limits>
#include <optional>

namespace cfd::verification
{

inline double rpa_nozzle_acoustic_time(const double length)
{
    const double time{length /
                      std::sqrt(EulerRpaModel::gamma * EulerRpaModel::gas_constant * EulerRpaModel::temperature)};
    if (!std::isfinite(length) || !(length > 0) || !std::isfinite(time) || !(time > 0))
        throw std::invalid_argument("RPA acoustic time requires a representable positive nozzle length.");
    return time;
}

/// Five physical output targets, optionally one extra comparison time, plus
/// an initial export handled by the driver. No allocation is required.
/// Clipping uses the existing solver.step(target); exported timestamps always
/// come from accepted diagnostics, never from this requested schedule.
class RpaPhysicalOutputs
{
  public:
    explicit RpaPhysicalOutputs(const double final_time, const std::optional<double> comparison_time = {})
    {
        constexpr std::array fractions{.1, .25, .5, .75, 1.0};
        double previous{};
        for (Index k = 0; k < fractions.size(); ++k)
        {
            times_.at(k) = fractions.at(k) * final_time;
            if (!std::isfinite(times_.at(k)) || !(times_.at(k) > previous))
                throw std::invalid_argument("RPA physical output times must be positive, finite and increasing.");
            previous = times_.at(k);
        }
        if (comparison_time)
        {
            if (!std::isfinite(*comparison_time) || !(*comparison_time > 0) || !(*comparison_time < final_time))
                throw std::invalid_argument("RPA comparison time must be finite and strictly inside the run interval.");
            const auto planned{times()};
            if (std::ranges::find(planned, *comparison_time) == planned.end())
            {
                times_.back() = *comparison_time;
                count_ = times_.size();
                std::ranges::sort(times_);
            }
        }
    }

    [[nodiscard]] std::span<const double> times() const noexcept
    {
        return {times_.data(), count_};
    }

    /// New output segment: past targets are neither re-exported nor relabelled.
    void resume_after(const double accepted_time)
    {
        if (!std::isfinite(accepted_time) || accepted_time < 0)
            throw std::invalid_argument("Invalid RPA restart time.");
        while (next_ < count_ && times_.at(next_) <= accepted_time)
            ++next_;
    }

    [[nodiscard]] Index next_index() const noexcept
    {
        return next_;
    }

    [[nodiscard]] double next_time() const
    {
        return times_.at(next_ < count_ ? next_ : count_ - 1);
    }

    [[nodiscard]] bool due(const double accepted_time) const
    {
        return next_ < count_ && accepted_time >= times_.at(next_);
    }

    void exported(const double accepted_time)
    {
        if (!due(accepted_time) || accepted_time != times_.at(next_))
            throw std::runtime_error("RPA physical export must match its accepted solver time exactly.");
        ++next_;
    }

  private:
    std::array<double, 6> times_{};
    Index count_{5};
    Index next_{};
};

struct RpaSectionSample
{
    double mach{}, axial_mach{}, minimum_mach{std::numeric_limits<double>::infinity()}, maximum_mach{};
    double pressure{}, temperature{}, mass_flow{}, area{};
};

/// Diagnostic for an actual face-aligned x-section only, not a general cutter.
/// Face IDs, Engine-side cells and revolution areas are selected once. Scalar
/// values are area-weighted Engine-cell estimates (two-sided at the throat,
/// one-sided at inlet/exit), not claimed pointwise physical traces. Mass flow
/// instead uses exactly the production Rusanov/characteristic flux, toward +x.
class RpaFaceSection
{
  public:
    RpaFaceSection(const Mesh &mesh, const AxisymmetricGeometry &geometry, const std::span<const double> regions,
                   const double x, const double radius)
        : mesh_(&mesh)
    {
        if (&geometry.mesh() != &mesh || regions.size() != mesh.cell_count() || !std::isfinite(x) ||
            !std::isfinite(radius) || !(radius > 0))
            throw std::invalid_argument("Invalid RPA face-section geometry/cardinality.");
        const double tolerance{128 * std::numeric_limits<double>::epsilon() * std::max(std::abs(x), radius)};
        for (Index face = 0; face < mesh.face_count(); ++face)
        {
            const auto nodes{mesh.faces()[face].node_ids};
            const auto a{mesh.nodes()[nodes[0]]}, b{mesh.nodes()[nodes[1]]};
            if (std::abs(a.x - x) > tolerance || std::abs(b.x - x) > tolerance ||
                std::max(a.y, b.y) > radius + tolerance || a.y == b.y)
                continue;
            const auto adjacency{mesh.face_adjacencies()[face]};
            const bool owner_engine{regions[adjacency.owner] == static_cast<double>(RpaRegion::Engine)};
            const bool neighbor_engine{!adjacency.is_boundary() &&
                                       regions[adjacency.neighbor] == static_cast<double>(RpaRegion::Engine)};
            if (!owner_engine && !neighbor_engine)
                continue;
            const auto area{geometry.face_area_vectors()[face]};
            faces_.push_back({face, owner_engine ? adjacency.owner : adjacency.neighbor,
                              owner_engine && neighbor_engine ? adjacency.neighbor : invalid_index,
                              std::hypot(area.x, area.y)});
        }
        // Exact mesh-face coverage: no gap, overlap or silently partial section.
        const auto lower = [&](const Entry &entry) {
            const auto ids{mesh.faces()[entry.face].node_ids};
            return std::min(mesh.nodes()[ids[0]].y, mesh.nodes()[ids[1]].y);
        };
        std::ranges::sort(faces_, {}, lower);
        double end{};
        for (const auto &entry : faces_)
        {
            const auto ids{mesh.faces()[entry.face].node_ids};
            if (std::abs(lower(entry) - end) > tolerance || !(entry.area > 0))
                throw std::runtime_error("RPA section has a radial gap, overlap or zero area.");
            end = std::max(mesh.nodes()[ids[0]].y, mesh.nodes()[ids[1]].y);
        }
        if (faces_.empty() || std::abs(end - radius) > tolerance)
            throw std::runtime_error("RPA x-section is not fully aligned with actual Engine mesh faces.");
    }

    [[nodiscard]] Index face_count() const noexcept
    {
        return faces_.size();
    }

    [[nodiscard]] RpaSectionSample sample(const Mesh &mesh, const ConservativeField &field,
                                          const CaloricallyPerfectGas &gas,
                                          const EulerBoundaryConditions &conditions) const
    {
        if (&mesh != mesh_ || field.size() != mesh.cell_count() || conditions.size() != mesh.boundary_groups().size())
            throw std::invalid_argument("RPA section requires its original Mesh and matching field/BC cardinalities.");
        RpaSectionSample result;
        for (const auto &entry : faces_)
        {
            const auto a{gas.evaluate(field.state(entry.first))};
            const auto b{entry.second == invalid_index ? a : gas.evaluate(field.state(entry.second))};
            const auto mach = [](const auto &state) {
                return std::hypot(state.primitive.velocity_x, state.primitive.velocity_r) /
                       state.thermodynamics.sound_speed;
            };
            const double m{.5 * (mach(a) + mach(b))};
            result.area += entry.area;
            result.mach += entry.area * m;
            result.axial_mach += entry.area * .5 *
                                 (a.primitive.velocity_x / a.thermodynamics.sound_speed +
                                  b.primitive.velocity_x / b.thermodynamics.sound_speed);
            result.minimum_mach = std::min(result.minimum_mach, m);
            result.maximum_mach = std::max(result.maximum_mach, m);
            result.pressure += entry.area * .5 * (a.primitive.pressure + b.primitive.pressure);
            result.temperature += entry.area * .5 * (a.thermodynamics.temperature + b.thermodynamics.temperature);
            const auto adjacency{mesh.face_adjacencies()[entry.face]};
            const auto sf{mesh.face_area_vectors()[entry.face]};
            const double length{mesh.face_lengths()[entry.face]};
            const Vector2 normal{sf.x / length, sf.y / length};
            const auto flux{
                adjacency.is_boundary()
                    ? compute_euler_characteristic_boundary_flux(
                          field.state(adjacency.owner), normal,
                          conditions.condition(mesh.face_boundary_ids()[entry.face]), gas)
                          .flux
                    : compute_rusanov_flux(field.state(adjacency.owner), field.state(adjacency.neighbor), normal, gas)
                          .flux};
            result.mass_flow += entry.area * (normal.x > 0 ? 1.0 : -1.0) * flux.mass;
        }
        result.mach /= result.area;
        result.axial_mach /= result.area;
        result.pressure /= result.area;
        result.temperature /= result.area;
        return result;
    }

  private:
    struct Entry
    {
        Index face, first, second;
        double area;
    };

    const Mesh *mesh_; ///< Non-owning; the original Mesh must outlive this diagnostic.
    std::vector<Entry> faces_;
};
} // namespace cfd::verification
