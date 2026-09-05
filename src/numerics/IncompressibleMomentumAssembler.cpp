#include "cfd/numerics/IncompressibleMomentumAssembler.hpp"

#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVectorField.hpp"
#include "cfd/field/CellVelocityField.hpp"
#include "cfd/field/FaceFluxField.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/linear_algebra/ScalarLinearSystem.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/Types.hpp"

#include <cmath>
#include <stdexcept>

namespace cfd
{
namespace
{

void validate_cell_field_cardinalities(const Mesh &mesh, const CellVelocityField &previous_velocity,
                                       const CellVectorField &u_gradient, const CellVectorField &v_gradient,
                                       const CellVectorField &pressure_gradient)
{
    const Index cell_count{mesh.cell_count()};
    if (previous_velocity.u().size() != cell_count || previous_velocity.v().size() != cell_count)
    {
        throw std::invalid_argument("Previous velocity cardinality must match the momentum Mesh cell count.");
    }
    if (u_gradient.size() != cell_count || v_gradient.size() != cell_count || pressure_gradient.size() != cell_count)
    {
        throw std::invalid_argument("Momentum gradient cardinalities must match the Mesh cell count.");
    }
}

void validate_boundary_condition_cardinalities(const Mesh &mesh, const ScalarBoundaryConditions &u_boundary_conditions,
                                               const ScalarBoundaryConditions &v_boundary_conditions)
{
    const Index boundary_count{mesh.boundary_groups().size()};
    if (u_boundary_conditions.size() != boundary_count || v_boundary_conditions.size() != boundary_count)
    {
        throw std::invalid_argument("Momentum boundary-condition counts must match the Mesh boundary count.");
    }
}

void validate_output_system(const Mesh &mesh, const ScalarLinearSystem &system)
{
    if (&system.mesh() != &mesh || system.cell_count() != mesh.cell_count() || system.face_count() != mesh.face_count())
    {
        throw std::invalid_argument("Momentum output systems must reference the assembler Mesh.");
    }
}

void apply_equation_relaxation(const CellScalarField &previous_field, const double relaxation_factor,
                               ScalarLinearSystem &system) noexcept
{
    if (relaxation_factor == 1.0)
    {
        return;
    }

    auto diagonal{system.diagonal()};
    auto rhs{system.rhs()};
    const auto previous_values{previous_field.values()};
    const double previous_value_coefficient{(1.0 - relaxation_factor) / relaxation_factor};

    for (Index cell_id = 0; cell_id < system.cell_count(); ++cell_id)
    {
        const double unrelaxed_diagonal{diagonal[cell_id]};
        diagonal[cell_id] = unrelaxed_diagonal / relaxation_factor;
        rhs[cell_id] += previous_value_coefficient * unrelaxed_diagonal * previous_values[cell_id];
    }
}

} // namespace

IncompressibleMomentumAssembler::IncompressibleMomentumAssembler(const Mesh &mesh, const double dynamic_viscosity,
                                                                 const ScalarConvectionScheme convection_scheme)
    : mesh_(&mesh), diffusion_(mesh, dynamic_viscosity), convection_(mesh, convection_scheme)
{
}

void IncompressibleMomentumAssembler::assemble(const CellVelocityField &previous_velocity,
                                               const CellVectorField &u_gradient, const CellVectorField &v_gradient,
                                               const CellVectorField &pressure_gradient,
                                               const ScalarBoundaryConditions &u_boundary_conditions,
                                               const ScalarBoundaryConditions &v_boundary_conditions,
                                               const FaceFluxField &mass_flux, const double relaxation_factor,
                                               ScalarLinearSystem &u_system, ScalarLinearSystem &v_system) const
{
    validate_cell_field_cardinalities(*mesh_, previous_velocity, u_gradient, v_gradient, pressure_gradient);
    if (mass_flux.size() != mesh_->face_count())
    {
        throw std::invalid_argument("Momentum mass-flux cardinality must match the Mesh face count.");
    }
    validate_boundary_condition_cardinalities(*mesh_, u_boundary_conditions, v_boundary_conditions);
    if (!std::isfinite(relaxation_factor) || !(relaxation_factor > 0.0) || relaxation_factor > 1.0)
    {
        throw std::invalid_argument("Momentum relaxation factor must be finite and in (0, 1].");
    }
    if (&u_system == &v_system)
    {
        throw std::invalid_argument("Momentum component output systems must be distinct objects.");
    }
    validate_output_system(*mesh_, u_system);
    validate_output_system(*mesh_, v_system);

    u_system.clear();
    v_system.clear();

    diffusion_.add_matrix_contributions(u_boundary_conditions, u_system);
    convection_.add_matrix_contributions(u_boundary_conditions, mass_flux, u_system);
    diffusion_.add_boundary_rhs(u_boundary_conditions, u_system.rhs());
    convection_.add_boundary_rhs(u_boundary_conditions, mass_flux, u_system.rhs());
    diffusion_.add_non_orthogonal_rhs(u_boundary_conditions, u_gradient, u_system.rhs());

    diffusion_.add_matrix_contributions(v_boundary_conditions, v_system);
    convection_.add_matrix_contributions(v_boundary_conditions, mass_flux, v_system);
    diffusion_.add_boundary_rhs(v_boundary_conditions, v_system.rhs());
    convection_.add_boundary_rhs(v_boundary_conditions, mass_flux, v_system.rhs());
    diffusion_.add_non_orthogonal_rhs(v_boundary_conditions, v_gradient, v_system.rhs());

    const auto cell_areas{mesh_->cell_areas()};
    const auto pressure_gradients{pressure_gradient.values()};
    auto u_rhs{u_system.rhs()};
    auto v_rhs{v_system.rhs()};
    for (Index cell_id = 0; cell_id < mesh_->cell_count(); ++cell_id)
    {
        u_rhs[cell_id] -= cell_areas[cell_id] * pressure_gradients[cell_id].x;
        v_rhs[cell_id] -= cell_areas[cell_id] * pressure_gradients[cell_id].y;
    }

    apply_equation_relaxation(previous_velocity.u(), relaxation_factor, u_system);
    apply_equation_relaxation(previous_velocity.v(), relaxation_factor, v_system);
}

} // namespace cfd
