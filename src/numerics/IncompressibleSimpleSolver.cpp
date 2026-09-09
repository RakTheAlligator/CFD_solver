#include "cfd/numerics/IncompressibleSimpleSolver.hpp"

#include "cfd/field/FaceFluxField.hpp"
#include "cfd/field/PressureCorrectionBoundaryConditions.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/linear_algebra/LinearSolveResult.hpp"
#include "cfd/math/Vector2.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/numerics/LeastSquaresGradient.hpp"
#include "cfd/numerics/MassFluxBalance.hpp"
#include "cfd/numerics/MomentumPressureResponse.hpp"
#include "cfd/numerics/PressureCorrectionReference.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cfd
{
namespace
{

void require_connected_cell_domain(const Mesh &mesh)
{
    const Index cell_count{mesh.cell_count()};
    if (cell_count == 0)
    {
        throw std::invalid_argument("SIMPLE v1 requires a single connected cell domain.");
    }

    std::vector<bool> reached(cell_count);
    std::vector<Index> pending_cells;
    pending_cells.reserve(cell_count);
    reached[0] = true;
    pending_cells.push_back(0);

    Index reached_cell_count{};
    const auto cell_offsets{mesh.cell_node_offsets()};
    const auto cell_faces{mesh.cell_faces()};
    const auto face_adjacencies{mesh.face_adjacencies()};
    while (!pending_cells.empty())
    {
        const Index cell_id{pending_cells.back()};
        pending_cells.pop_back();
        ++reached_cell_count;

        for (Index position = cell_offsets[cell_id]; position < cell_offsets[cell_id + 1]; ++position)
        {
            const FaceAdjacency &adjacency{face_adjacencies[cell_faces[position]]};
            if (adjacency.is_boundary())
            {
                continue;
            }

            const Index adjacent_cell{adjacency.owner == cell_id ? adjacency.neighbor : adjacency.owner};
            if (!reached[adjacent_cell])
            {
                reached[adjacent_cell] = true;
                pending_cells.push_back(adjacent_cell);
            }
        }
    }

    if (reached_cell_count != cell_count)
    {
        throw std::invalid_argument("SIMPLE v1 requires a single connected cell domain.");
    }
}

[[nodiscard]]
double validate_positive_physical_coefficient(const double value, const char *const name)
{
    if (!std::isfinite(value) || !(value > 0.0))
    {
        throw std::invalid_argument(std::string{name} + " must be finite and strictly positive.");
    }
    return value;
}

[[nodiscard]]
IncompressibleSimpleOptions validate_options(IncompressibleSimpleOptions options)
{
    if (options.maximum_iterations == 0)
    {
        throw std::invalid_argument("SIMPLE maximum iterations must be nonzero.");
    }
    if (options.momentum_relaxation_factor != 1.0)
    {
        throw std::invalid_argument(
            "SIMPLE v1 requires a momentum relaxation factor of exactly 1 until relaxation-consistent momentum "
            "interpolation is implemented.");
    }
    if (!std::isfinite(options.pressure_relaxation_factor) || !(options.pressure_relaxation_factor > 0.0) ||
        !(options.pressure_relaxation_factor <= 1.0))
    {
        throw std::invalid_argument("SIMPLE pressure relaxation factor must be finite and in (0, 1].");
    }
    if (!std::isfinite(options.rhie_chow_flux_relaxation_factor) || !(options.rhie_chow_flux_relaxation_factor > 0.0) ||
        !(options.rhie_chow_flux_relaxation_factor <= 1.0))
    {
        throw std::invalid_argument("SIMPLE Rhie-Chow flux relaxation factor must be finite and in (0, 1].");
    }
    if (!std::isfinite(options.velocity_relative_tolerance) || !(options.velocity_relative_tolerance > 0.0) ||
        !(options.velocity_relative_tolerance < 1.0))
    {
        throw std::invalid_argument("SIMPLE velocity relative tolerance must be finite and in (0, 1).");
    }
    if (!std::isfinite(options.continuity_relative_tolerance) || !(options.continuity_relative_tolerance > 0.0) ||
        !(options.continuity_relative_tolerance < 1.0))
    {
        throw std::invalid_argument("SIMPLE continuity relative tolerance must be finite and in (0, 1).");
    }
    return options;
}

void validate_field_cardinalities(const Mesh &mesh, const CellVelocityField &velocity, const CellScalarField &pressure,
                                  const FaceFluxField &mass_flux)
{
    const Index cell_count{mesh.cell_count()};
    if (velocity.size() != cell_count || velocity.u().size() != cell_count || velocity.v().size() != cell_count)
    {
        throw std::invalid_argument("SIMPLE velocity cardinality must match the Mesh cell count.");
    }
    if (pressure.size() != cell_count)
    {
        throw std::invalid_argument("SIMPLE pressure cardinality must match the Mesh cell count.");
    }
    if (mass_flux.size() != mesh.face_count())
    {
        throw std::invalid_argument("SIMPLE mass-flux cardinality must match the Mesh face count.");
    }
}

void validate_boundary_cardinalities(
    const Mesh &mesh, const ScalarBoundaryConditions &u_boundary_conditions,
    const ScalarBoundaryConditions &v_boundary_conditions, const ScalarBoundaryConditions &pressure_boundary_conditions,
    const PressureCorrectionBoundaryConditions &pressure_correction_boundary_conditions)
{
    const Index boundary_count{mesh.boundary_groups().size()};
    if (u_boundary_conditions.size() != boundary_count || v_boundary_conditions.size() != boundary_count ||
        pressure_boundary_conditions.size() != boundary_count ||
        pressure_correction_boundary_conditions.size() != boundary_count)
    {
        throw std::invalid_argument("SIMPLE boundary-condition counts must match the Mesh boundary count.");
    }
}

void validate_initial_values(const CellVelocityField &velocity, const CellScalarField &pressure,
                             const FaceFluxField &mass_flux)
{
    for (Index cell_id = 0; cell_id < velocity.size(); ++cell_id)
    {
        if (!std::isfinite(velocity.u()[cell_id]) || !std::isfinite(velocity.v()[cell_id]))
        {
            throw std::invalid_argument("SIMPLE initial velocity values must be finite.");
        }
        if (!std::isfinite(pressure[cell_id]))
        {
            throw std::invalid_argument("SIMPLE initial pressure values must be finite.");
        }
    }
    for (const double flux : mass_flux.values())
    {
        if (!std::isfinite(flux))
        {
            throw std::invalid_argument("SIMPLE initial mass-flux values must be finite.");
        }
    }
}

[[nodiscard]]
bool validate_pressure_boundary_semantics(
    const ScalarBoundaryConditions &pressure_boundary_conditions,
    const PressureCorrectionBoundaryConditions &pressure_correction_boundary_conditions)
{
    bool has_fixed_pressure{};
    for (BoundaryId boundary_id = 0; boundary_id < pressure_boundary_conditions.size(); ++boundary_id)
    {
        const ScalarBoundaryConditionType physical_type{pressure_boundary_conditions[boundary_id].type};
        switch (pressure_correction_boundary_conditions[boundary_id])
        {
        case PressureCorrectionBoundaryConditionType::FixedMassFlux:
            if (physical_type != ScalarBoundaryConditionType::Neumann)
            {
                throw std::invalid_argument(
                    "SIMPLE FixedMassFlux pressure correction requires a physical pressure Neumann condition.");
            }
            break;

        case PressureCorrectionBoundaryConditionType::FixedPressure:
            if (physical_type != ScalarBoundaryConditionType::Dirichlet)
            {
                throw std::invalid_argument(
                    "SIMPLE FixedPressure pressure correction requires a physical pressure Dirichlet condition.");
            }
            has_fixed_pressure = true;
            break;
        }
    }
    return has_fixed_pressure;
}

[[nodiscard]]
ScalarBoundaryConditions make_pressure_correction_gradient_boundary_conditions(
    const PressureCorrectionBoundaryConditions &pressure_correction_boundary_conditions)
{
    std::vector<ScalarBoundaryCondition> conditions;
    conditions.reserve(pressure_correction_boundary_conditions.size());
    for (BoundaryId boundary_id = 0; boundary_id < pressure_correction_boundary_conditions.size(); ++boundary_id)
    {
        switch (pressure_correction_boundary_conditions[boundary_id])
        {
        case PressureCorrectionBoundaryConditionType::FixedMassFlux:
            conditions.emplace_back(ScalarBoundaryConditionType::Neumann, 0.0);
            break;

        case PressureCorrectionBoundaryConditionType::FixedPressure:
            conditions.emplace_back(ScalarBoundaryConditionType::Dirichlet, 0.0);
            break;
        }
    }
    return {pressure_correction_boundary_conditions.size(), std::move(conditions)};
}

void validate_fixed_mass_flux_compatibility(const Mesh &mesh, const FaceFluxField &mass_flux)
{
    long double boundary_flux_sum{};
    long double boundary_flux_magnitude{};
    const auto face_adjacencies{mesh.face_adjacencies()};
    for (Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (!face_adjacencies[face_id].is_boundary())
        {
            continue;
        }
        const long double flux{mass_flux[face_id]};
        boundary_flux_sum += flux;
        boundary_flux_magnitude += std::abs(flux);
    }

    constexpr long double relative_tolerance{64.0L * std::numeric_limits<double>::epsilon()};
    const long double compatibility_tolerance{relative_tolerance * boundary_flux_magnitude};
    if (std::abs(boundary_flux_sum) > compatibility_tolerance)
    {
        throw std::runtime_error("SIMPLE all-FixedMassFlux boundary flow is globally incompatible.");
    }
}

[[nodiscard]]
double relative_velocity_change(const CellVelocityField &velocity, const CellVelocityField &previous_velocity)
{
    double maximum_change{};
    double velocity_scale{};
    for (Index cell_id = 0; cell_id < velocity.size(); ++cell_id)
    {
        const double new_magnitude{std::hypot(velocity.u()[cell_id], velocity.v()[cell_id])};
        const double old_magnitude{std::hypot(previous_velocity.u()[cell_id], previous_velocity.v()[cell_id])};
        const double change_magnitude{std::hypot(velocity.u()[cell_id] - previous_velocity.u()[cell_id],
                                                 velocity.v()[cell_id] - previous_velocity.v()[cell_id])};
        maximum_change = std::max(maximum_change, change_magnitude);
        velocity_scale = std::max({velocity_scale, new_magnitude, old_magnitude});
    }
    if (velocity_scale == 0.0)
    {
        return maximum_change == 0.0 ? 0.0 : std::numeric_limits<double>::infinity();
    }
    return maximum_change / velocity_scale;
}

struct ContinuityDiagnostics
{
    double relative_residual{};
    double maximum_imbalance{};
};

[[nodiscard]]
ContinuityDiagnostics compute_continuity_diagnostics(const std::span<const double> mass_imbalance,
                                                     const FaceFluxField &mass_flux)
{
    double imbalance_sum{};
    double maximum_imbalance{};
    for (const double imbalance : mass_imbalance)
    {
        const double magnitude{std::abs(imbalance)};
        imbalance_sum += magnitude;
        maximum_imbalance = std::max(maximum_imbalance, magnitude);
    }

    double flux_sum{};
    for (const double flux : mass_flux.values())
    {
        flux_sum += std::abs(flux);
    }
    if (flux_sum == 0.0)
    {
        return {imbalance_sum == 0.0 ? 0.0 : std::numeric_limits<double>::infinity(), maximum_imbalance};
    }
    return {imbalance_sum / flux_sum, maximum_imbalance};
}

[[nodiscard]]
double maximum_absolute_value(const CellScalarField &field)
{
    double maximum{};
    for (const double value : field.values())
    {
        maximum = std::max(maximum, std::abs(value));
    }
    return maximum;
}

void require_converged(const LinearSolveResult &result, const char *const system_name)
{
    if (!result.converged)
    {
        throw std::runtime_error(std::string{"SIMPLE "} + system_name + " linear solve did not converge.");
    }
}

} // namespace

IncompressibleSimpleSolver::IncompressibleSimpleSolver(const Mesh &mesh, const double density,
                                                       const double dynamic_viscosity,
                                                       const ScalarConvectionScheme convection_scheme,
                                                       IncompressibleSimpleOptions options)
    : mesh_(&mesh), density_(validate_positive_physical_coefficient(density, "SIMPLE density")),
      options_(validate_options(options)), u_momentum_solver_(options_.momentum_linear_solver),
      v_momentum_solver_(options_.momentum_linear_solver),
      pressure_correction_solver_(options_.pressure_correction_linear_solver),
      momentum_assembler_(mesh, validate_positive_physical_coefficient(dynamic_viscosity, "SIMPLE dynamic viscosity"),
                          convection_scheme),
      internal_face_interpolation_(mesh, density_), boundary_face_interpolation_(mesh, density_),
      pressure_correction_assembler_(mesh), pressure_velocity_corrector_(mesh), previous_velocity_(mesh.cell_count()),
      previous_mass_flux_(mesh.face_count()), u_gradient_(mesh.cell_count()), v_gradient_(mesh.cell_count()),
      pressure_gradient_(mesh.cell_count()), pressure_correction_gradient_(mesh.cell_count()),
      momentum_response_(mesh.cell_count()), face_pressure_response_(mesh.face_count()),
      pressure_correction_(mesh.cell_count()), mass_imbalance_(mesh.cell_count()), u_momentum_system_(mesh),
      v_momentum_system_(mesh), pressure_correction_system_(mesh)
{
    require_connected_cell_domain(mesh);
}

IncompressibleSimpleResult IncompressibleSimpleSolver::solve(
    const ScalarBoundaryConditions &u_boundary_conditions, const ScalarBoundaryConditions &v_boundary_conditions,
    const ScalarBoundaryConditions &pressure_boundary_conditions,
    const PressureCorrectionBoundaryConditions &pressure_correction_boundary_conditions, CellVelocityField &velocity,
    CellScalarField &pressure, FaceFluxField &mass_flux, SimpleIterationCallback iteration_callback)
{
    validate_field_cardinalities(*mesh_, velocity, pressure, mass_flux);
    validate_boundary_cardinalities(*mesh_, u_boundary_conditions, v_boundary_conditions, pressure_boundary_conditions,
                                    pressure_correction_boundary_conditions);
    validate_initial_values(velocity, pressure, mass_flux);
    const bool has_fixed_pressure{
        validate_pressure_boundary_semantics(pressure_boundary_conditions, pressure_correction_boundary_conditions)};
    const ScalarBoundaryConditions pressure_correction_gradient_boundary_conditions{
        make_pressure_correction_gradient_boundary_conditions(pressure_correction_boundary_conditions)};
    if (!has_fixed_pressure)
    {
        validate_fixed_mass_flux_compatibility(*mesh_, mass_flux);
    }

    IncompressibleSimpleResult result{};
    for (Index iteration = 0; iteration < options_.maximum_iterations; ++iteration)
    {
        const Index iteration_count{iteration + 1};
        std::copy(velocity.u().values().begin(), velocity.u().values().end(), previous_velocity_.u().values().begin());
        std::copy(velocity.v().values().begin(), velocity.v().values().end(), previous_velocity_.v().values().begin());

        compute_least_squares_gradient(*mesh_, previous_velocity_.u(), u_boundary_conditions, u_gradient_);
        compute_least_squares_gradient(*mesh_, previous_velocity_.v(), v_boundary_conditions, v_gradient_);
        compute_least_squares_gradient(*mesh_, pressure, pressure_boundary_conditions, pressure_gradient_);

        momentum_assembler_.assemble(previous_velocity_, u_gradient_, v_gradient_, pressure_gradient_,
                                     u_boundary_conditions, v_boundary_conditions, mass_flux,
                                     options_.momentum_relaxation_factor, u_momentum_system_, v_momentum_system_);
        u_momentum_solver_.compute_matrix(u_momentum_system_);
        const LinearSolveResult u_solve{u_momentum_solver_.solve(u_momentum_system_.rhs(), velocity.u().values())};
        require_converged(u_solve, "u-momentum");
        v_momentum_solver_.compute_matrix(v_momentum_system_);
        const LinearSolveResult v_solve{v_momentum_solver_.solve(v_momentum_system_.rhs(), velocity.v().values())};
        require_converged(v_solve, "v-momentum");

        compute_momentum_pressure_response(*mesh_, u_momentum_system_, v_momentum_system_, momentum_response_);
        const bool relax_rhie_chow_flux{options_.rhie_chow_flux_relaxation_factor != 1.0};

        if (relax_rhie_chow_flux)
        {
            std::copy(mass_flux.values().begin(), mass_flux.values().end(), previous_mass_flux_.values().begin());
        }
        internal_face_interpolation_.update_internal_faces(velocity, pressure, pressure_gradient_, momentum_response_,
                                                           mass_flux, face_pressure_response_);
        boundary_face_interpolation_.update_fixed_pressure_boundaries(
            velocity, pressure, pressure_gradient_, momentum_response_, pressure_boundary_conditions,
            pressure_correction_boundary_conditions, mass_flux, face_pressure_response_);
        if (relax_rhie_chow_flux)
        {
            const double relaxation_factor{options_.rhie_chow_flux_relaxation_factor};
            const auto face_adjacencies{mesh_->face_adjacencies()};
            const auto face_boundary_ids{mesh_->face_boundary_ids()};
            for (Index face_id = 0; face_id < mesh_->face_count(); ++face_id)
            {
                if (face_adjacencies[face_id].is_boundary())
                {
                    switch (pressure_correction_boundary_conditions[face_boundary_ids[face_id]])
                    {
                    case PressureCorrectionBoundaryConditionType::FixedMassFlux:
                        continue;

                    case PressureCorrectionBoundaryConditionType::FixedPressure:
                        break;
                    }
                }
                mass_flux[face_id] =
                    relaxation_factor * mass_flux[face_id] + (1.0 - relaxation_factor) * previous_mass_flux_[face_id];
            }
        }

        pressure_correction_system_.clear();
        pressure_correction_assembler_.add_internal_face_contributions(mass_flux, face_pressure_response_,
                                                                       pressure_correction_system_);
        pressure_correction_assembler_.add_boundary_provisional_flux_rhs(mass_flux, pressure_correction_system_);
        pressure_correction_assembler_.add_boundary_pressure_response(
            pressure_correction_boundary_conditions, face_pressure_response_, pressure_correction_system_);
        const ContinuityDiagnostics provisional_continuity{
            compute_continuity_diagnostics(pressure_correction_system_.rhs(), mass_flux)};
        if (!has_fixed_pressure)
        {
            apply_zero_pressure_correction_reference(0, pressure_correction_system_);
        }

        std::fill(pressure_correction_.values().begin(), pressure_correction_.values().end(), 0.0);
        pressure_correction_solver_.compute_matrix(pressure_correction_system_);
        const LinearSolveResult pressure_correction_solve{
            pressure_correction_solver_.solve(pressure_correction_system_.rhs(), pressure_correction_.values())};
        require_converged(pressure_correction_solve, "pressure-correction");
        compute_least_squares_gradient(*mesh_, pressure_correction_, pressure_correction_gradient_boundary_conditions,
                                       pressure_correction_gradient_);

        pressure_velocity_corrector_.correct_face_mass_flux(
            pressure_correction_, pressure_correction_boundary_conditions, face_pressure_response_, mass_flux);
        pressure_velocity_corrector_.correct_velocity(pressure_correction_gradient_, momentum_response_, velocity);
        pressure_velocity_corrector_.correct_pressure(pressure_correction_, options_.pressure_relaxation_factor,
                                                      pressure);
        compute_cell_mass_imbalance(*mesh_, mass_flux, mass_imbalance_);

        const ContinuityDiagnostics continuity{compute_continuity_diagnostics(mass_imbalance_.values(), mass_flux)};
        result = {
            false,
            iteration_count,
            relative_velocity_change(velocity, previous_velocity_),
            provisional_continuity.relative_residual,
            continuity.relative_residual,
            continuity.maximum_imbalance,
            maximum_absolute_value(pressure_correction_),
        };
        result.converged = result.velocity_relative_change <= options_.velocity_relative_tolerance &&
                           result.provisional_continuity_relative_residual <= options_.continuity_relative_tolerance;
        const SimpleIterationInfo iteration_info{
            .iteration = iteration_count,
            .u_solve = u_solve,
            .v_solve = v_solve,
            .pressure_correction_solve = pressure_correction_solve,
            .velocity_relative_change = result.velocity_relative_change,
            .provisional_continuity_relative_residual = result.provisional_continuity_relative_residual,
            .corrected_continuity_relative_residual = result.continuity_relative_residual,
            .maximum_pressure_correction = result.maximum_pressure_correction,
        };
        if (iteration_callback)
        {
            iteration_callback(iteration_info);
        }
        if (result.converged)
        {
            return result;
        }
    }
    return result;
}

} // namespace cfd
