#pragma once

#include "cfd/numerics/ScalarConvectionOperator.hpp"
#include "cfd/numerics/ScalarDiffusionOperator.hpp"

#include <vector>

namespace cfd
{

class CellScalarField;
class CellVectorField;
class FaceFluxField;
class Mesh;
class ScalarBoundaryConditions;
class ScalarLinearSystem;

/// Assembles integrated steady scalar transport: `div(F phi) - div(Gamma grad(phi)) = 0`.
///
/// The supplied face flux is already integrated and owner-oriented. Gamma is
/// a constant positive diffusion coefficient, not necessarily a kinematic
/// diffusivity. Boundary Neumann data prescribe the outward normal derivative
/// of phi, not its diffusive flux. No source or relaxation is added.
///
/// @note Mesh is not owned and must outlive this assembler.
/// @note Repeated valid calls allocate no storage. Instances reuse their own
///       limiter workspace and must not be used concurrently.
class ScalarTransportAssembler
{
  public:
    /// Constructs the composed operators and, only for Barth-Jespersen, a
    /// cell-sized limiter workspace.
    /// @throws std::invalid_argument If the diffusion coefficient is non-finite
    ///         or non-positive, or the scheme/limiter combination is unsupported.
    /// @throws std::runtime_error If face geometry is unusable by an operator.
    ScalarTransportAssembler(const Mesh &mesh, double diffusion_coefficient,
                             ScalarConvectionScheme scheme = ScalarConvectionScheme::FirstOrderUpwind,
                             ScalarConvectionLimiter limiter = ScalarConvectionLimiter::None);

    ScalarTransportAssembler(const ScalarTransportAssembler &) = delete;
    ScalarTransportAssembler &operator=(const ScalarTransportAssembler &) = delete;
    ScalarTransportAssembler(ScalarTransportAssembler &&) noexcept = default;
    ScalarTransportAssembler &operator=(ScalarTransportAssembler &&) noexcept = delete;
    ~ScalarTransportAssembler() = default;

    /// Replaces only the matrix with implicit convection plus primary diffusion.
    /// RHS is unchanged. For fixed flux and BC types, this matrix is independent
    /// of the scalar field and its explicit corrections.
    /// @throws std::invalid_argument If cardinalities or Mesh identity differ,
    ///         or a face flux is non-finite. Inputs are validated before clearing.
    void assemble_matrix(const ScalarBoundaryConditions &boundary_conditions, const FaceFluxField &face_flux,
                         ScalarLinearSystem &system) const;

    /// Replaces only RHS with boundary terms and explicit convection/diffusion
    /// corrections from previous_field and its raw WLS gradient. Matrix is unchanged.
    /// Hybrid borrows the composed diffusion operator's conductances; the limiter
    /// changes only the convective reconstruction, not the diffusion gradient.
    /// @throws std::invalid_argument If cardinalities or Mesh identity differ,
    ///         or field, gradient or face flux values are non-finite. Inputs are
    ///         validated before clearing.
    void assemble_rhs(const CellScalarField &previous_field, const CellVectorField &gradient,
                      const ScalarBoundaryConditions &boundary_conditions, const FaceFluxField &face_flux,
                      ScalarLinearSystem &system);

  private:
    const Mesh *mesh_;
    ScalarDiffusionOperator diffusion_;
    ScalarConvectionOperator convection_;
    std::vector<double> limiter_workspace_;
};

} // namespace cfd
