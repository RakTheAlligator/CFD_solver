#pragma once

#include "cfd/field/CellScalarField.hpp"

#include <span>

namespace cfd
{

/// Local 2D conservative state, with no dynamic storage.
struct ConservativeState
{
    double density{};              ///< rho [kg/m^3].
    double momentum_x{};           ///< rho*ux [kg/(m^2*s)].
    double momentum_r{};           ///< rho*ur [kg/(m^2*s)].
    double total_energy_density{}; ///< rho*E [J/m^3], not specific energy.
};

/// Fixed-cardinality, owning Structure-of-Arrays conservative field.
///
/// Owns four contiguous CellScalarField arrays (32 bytes/cell of values),
/// with four allocations for nonempty fields. Neither owns nor references Mesh.
/// Zero initialization is storage initialization, not an admissible gas state;
/// thermodynamic validation belongs to CaloricallyPerfectGas.
///
/// Copy construction intentionally performs a deep O(N) copy. Component
/// accessors return non-owning spans without allocation: values can be edited,
/// but owning arrays cannot be moved or replaced individually. Move construction
/// transfers the whole field's storage. Assignment is disabled to keep
/// cardinality fixed. Views must not outlive their underlying field storage.
class ConservativeField
{
  public:
    explicit ConservativeField(Index cell_count, const ConservativeState &initial_value = {});

    ConservativeField(const ConservativeField &) = default;
    ConservativeField &operator=(const ConservativeField &) = delete;
    ConservativeField(ConservativeField &&) noexcept = default;
    ConservativeField &operator=(ConservativeField &&) noexcept = delete;
    ~ConservativeField() = default;

    [[nodiscard]] Index size() const noexcept
    {
        return density_.size();
    }

    [[nodiscard]] std::span<double> density() noexcept
    {
        return density_.values();
    }

    [[nodiscard]] std::span<const double> density() const noexcept
    {
        return density_.values();
    }

    [[nodiscard]] std::span<double> momentum_x() noexcept
    {
        return momentum_x_.values();
    }

    [[nodiscard]] std::span<const double> momentum_x() const noexcept
    {
        return momentum_x_.values();
    }

    [[nodiscard]] std::span<double> momentum_r() noexcept
    {
        return momentum_r_.values();
    }

    [[nodiscard]] std::span<const double> momentum_r() const noexcept
    {
        return momentum_r_.values();
    }

    [[nodiscard]] std::span<double> total_energy_density() noexcept
    {
        return total_energy_density_.values();
    }

    [[nodiscard]] std::span<const double> total_energy_density() const noexcept
    {
        return total_energy_density_.values();
    }

    /// Loads a small local state without allocation. @pre cell < size().
    [[nodiscard]] ConservativeState state(Index cell) const noexcept;
    /// Stores a local state without validation/allocation. @pre cell < size().
    void set_state(Index cell, const ConservativeState &value) noexcept;

  private:
    CellScalarField density_;
    CellScalarField momentum_x_;
    CellScalarField momentum_r_;
    CellScalarField total_energy_density_;
};

} // namespace cfd
