#pragma once

#include "cfd/numerics/ScalarConvectionOperator.hpp"
#include "cfd/numerics/ScalarDiffusionOperator.hpp"

namespace cfd
{

class CellVectorField;
class CellVelocityField;
class FaceFluxField;
class Mesh;
class ScalarBoundaryConditions;
class ScalarLinearSystem;

/// Assembles steady incompressible momentum-predictor systems in two dimensions.
///
/// Per unit depth, the component equations are
/// `div(rho U u) - div(mu grad(u)) = -dp/dx` and
/// `div(rho U v) - div(mu grad(v)) = -dp/dy`. The supplied face flux already
/// contains density as the integrated mass flux `rho U_f . S_f`; this class
/// does not multiply convection by density again.
///
/// Pressure gradients and component gradients for explicit non-orthogonal
/// diffusion correction are supplied externally. Each call clears and fully
/// rebuilds both output systems, then algebraically under-relaxes their
/// diagonals. The final relaxed diagonals are intentionally retained for later
/// pressure-velocity coupling work.
///
/// @note The referenced Mesh is not owned and must outlive this assembler.
/// @note Repeated valid calls to `assemble()` perform no dynamic allocation.
class IncompressibleMomentumAssembler
{
  public:
    /// Constructs an assembler for a fixed Mesh and constant dynamic viscosity.
    ///
    /// @param mesh Validated Mesh referenced by all output systems.
    /// @param dynamic_viscosity Constant dynamic viscosity `mu`.
    /// @param convection_scheme Scalar convection interpolation scheme used for
    ///        both velocity components.
    /// @throws std::invalid_argument If `dynamic_viscosity` is non-finite or not
    ///         strictly positive, or if `convection_scheme` is unsupported.
    /// @throws std::runtime_error If face geometry is unusable by a composed
    ///         operator.
    IncompressibleMomentumAssembler(
        const Mesh &mesh, double dynamic_viscosity,
        ScalarConvectionScheme convection_scheme = ScalarConvectionScheme::FirstOrderUpwind);

    IncompressibleMomentumAssembler(const IncompressibleMomentumAssembler &) = delete;
    IncompressibleMomentumAssembler &operator=(const IncompressibleMomentumAssembler &) = delete;

    IncompressibleMomentumAssembler(IncompressibleMomentumAssembler &&) noexcept = default;
    IncompressibleMomentumAssembler &operator=(IncompressibleMomentumAssembler &&) noexcept = delete;

    ~IncompressibleMomentumAssembler() = default;

    /// Clears and assembles the complete relaxed `u` and `v` momentum systems.
    ///
    /// `pressure_gradient` is the physical cell-centered `grad(p)`. The
    /// pressure contributions are `-A_P dp/dx` and `-A_P dp/dy`. The supplied
    /// `u_gradient` and `v_gradient` are used only by explicit non-orthogonal
    /// diffusion correction. Equation relaxation leaves the final diagonal as
    /// `a_P / relaxation_factor`; no post-solve field blending is performed.
    ///
    /// @throws std::invalid_argument If a cardinality is incompatible,
    ///         `relaxation_factor` is non-finite or outside `(0, 1]`, an output
    ///         system references another Mesh, or both outputs are the same
    ///         object.
    void assemble(const CellVelocityField &previous_velocity, const CellVectorField &u_gradient,
                  const CellVectorField &v_gradient, const CellVectorField &pressure_gradient,
                  const ScalarBoundaryConditions &u_boundary_conditions,
                  const ScalarBoundaryConditions &v_boundary_conditions, const FaceFluxField &mass_flux,
                  double relaxation_factor, ScalarLinearSystem &u_system, ScalarLinearSystem &v_system) const;

  private:
    const Mesh *mesh_;
    ScalarDiffusionOperator diffusion_;
    ScalarConvectionOperator convection_;
};

} // namespace cfd
