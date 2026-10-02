#include "cfd/numerics/ScalarTransportAssembler.hpp"

#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVectorField.hpp"
#include "cfd/field/FaceFluxField.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/linear_algebra/ScalarLinearSystem.hpp"
#include "cfd/mesh/Mesh.hpp"

#include <cmath>
#include <stdexcept>

namespace cfd
{
namespace
{

void validate_inputs(const Mesh &mesh, const ScalarBoundaryConditions &boundary_conditions,
                     const FaceFluxField &face_flux, const ScalarLinearSystem &system)
{
    if (boundary_conditions.size() != mesh.boundary_groups().size() || face_flux.size() != mesh.face_count())
    {
        throw std::invalid_argument("Scalar transport boundary/flux cardinalities must match the Mesh.");
    }
    if (&system.mesh() != &mesh || system.cell_count() != mesh.cell_count() || system.face_count() != mesh.face_count())
    {
        throw std::invalid_argument("Scalar transport system must reference the assembler Mesh.");
    }
    for (const double flux : face_flux.values())
    {
        if (!std::isfinite(flux))
        {
            throw std::invalid_argument("Scalar transport face fluxes must be finite.");
        }
    }
}

} // namespace

ScalarTransportAssembler::ScalarTransportAssembler(const Mesh &mesh, const double diffusion_coefficient,
                                                   const ScalarConvectionScheme scheme,
                                                   const ScalarConvectionLimiter limiter)
    : mesh_(&mesh), diffusion_(mesh, diffusion_coefficient), convection_(mesh, scheme, limiter),
      limiter_workspace_(limiter == ScalarConvectionLimiter::BarthJespersen ? mesh.cell_count() : 0)
{
}

void ScalarTransportAssembler::assemble_matrix(const ScalarBoundaryConditions &boundary_conditions,
                                               const FaceFluxField &face_flux, ScalarLinearSystem &system) const
{
    validate_inputs(*mesh_, boundary_conditions, face_flux, system);
    system.clear_matrix();
    diffusion_.add_matrix_contributions(boundary_conditions, system);
    convection_.add_matrix_contributions(boundary_conditions, face_flux, system,
                                         diffusion_.face_primary_coefficients());
}

void ScalarTransportAssembler::assemble_rhs(const CellScalarField &previous_field, const CellVectorField &gradient,
                                            const ScalarBoundaryConditions &boundary_conditions,
                                            const FaceFluxField &face_flux, ScalarLinearSystem &system)
{
    validate_inputs(*mesh_, boundary_conditions, face_flux, system);
    if (previous_field.size() != mesh_->cell_count() || gradient.size() != mesh_->cell_count())
    {
        throw std::invalid_argument("Scalar transport field/gradient cardinalities must match the Mesh.");
    }
    for (Index cell_id = 0; cell_id < mesh_->cell_count(); ++cell_id)
    {
        const Vector2 &value{gradient[cell_id]};
        if (!std::isfinite(previous_field[cell_id]) || !std::isfinite(value.x) || !std::isfinite(value.y))
        {
            throw std::invalid_argument("Scalar transport fields and gradients must be finite.");
        }
    }

    system.clear_rhs();
    diffusion_.add_boundary_rhs(boundary_conditions, system.rhs());
    convection_.add_boundary_rhs(boundary_conditions, face_flux, system.rhs(), diffusion_.face_primary_coefficients());
    convection_.add_deferred_correction_rhs(previous_field, boundary_conditions, gradient, face_flux,
                                            limiter_workspace_, system.rhs());
    diffusion_.add_non_orthogonal_rhs(boundary_conditions, gradient, system.rhs());
}

} // namespace cfd
