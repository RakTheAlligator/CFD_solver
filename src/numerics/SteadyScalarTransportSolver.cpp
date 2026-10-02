#include "cfd/numerics/SteadyScalarTransportSolver.hpp"

#include "cfd/field/FaceFluxField.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/numerics/LeastSquaresGradient.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <span>
#include <sstream>
#include <stdexcept>

namespace cfd
{
namespace
{

[[nodiscard]]
SteadyScalarTransportOptions validated_options(const SteadyScalarTransportOptions options)
{
    if (options.maximum_iterations == 0 || !std::isfinite(options.field_relative_tolerance) ||
        !(options.field_relative_tolerance > 0.0) || !(options.field_relative_tolerance < 1.0) ||
        !std::isfinite(options.discrete_relative_tolerance) || !(options.discrete_relative_tolerance > 0.0) ||
        !(options.discrete_relative_tolerance < 1.0))
    {
        throw std::invalid_argument("Scalar transport requires positive iteration limit and tolerances in (0, 1).");
    }
    return options;
}

[[nodiscard]]
double relative_difference(const std::span<const double> first, const std::span<const double> second)
{
    double difference_norm{};
    double first_norm{};
    double second_norm{};
    for (Index index = 0; index < first.size(); ++index)
    {
        difference_norm = std::hypot(difference_norm, first[index] - second[index]);
        first_norm = std::hypot(first_norm, first[index]);
        second_norm = std::hypot(second_norm, second[index]);
    }
    if (!std::isfinite(difference_norm) || !std::isfinite(first_norm) || !std::isfinite(second_norm))
    {
        throw std::runtime_error("Scalar transport residual norms must be finite.");
    }
    const double scale{std::max(first_norm, second_norm)};
    if (scale == 0.0)
    {
        return difference_norm == 0.0 ? 0.0 : std::numeric_limits<double>::infinity();
    }
    return difference_norm / scale;
}

} // namespace

SteadyScalarTransportSolver::SteadyScalarTransportSolver(const Mesh &mesh, const double diffusion_coefficient,
                                                         const ScalarConvectionScheme scheme,
                                                         const ScalarConvectionLimiter limiter,
                                                         const SteadyScalarTransportOptions options)
    : mesh_(&mesh), options_(validated_options(options)), assembler_(mesh, diffusion_coefficient, scheme, limiter),
      system_(mesh), linear_solver_(options_.linear_solver), gradient_(mesh.cell_count()),
      previous_field_(mesh.cell_count()), matrix_product_(mesh.cell_count())
{
}

SteadyScalarTransportResult SteadyScalarTransportSolver::solve(const FaceFluxField &face_flux,
                                                               const ScalarBoundaryConditions &boundary_conditions,
                                                               CellScalarField &field)
{
    if (field.size() != mesh_->cell_count())
    {
        throw std::invalid_argument("Scalar transport field cardinality must match the Mesh.");
    }
    for (const double value : field.values())
    {
        if (!std::isfinite(value))
        {
            throw std::invalid_argument("Scalar transport initial field must be finite.");
        }
    }

    assembler_.assemble_matrix(boundary_conditions, face_flux, system_);
    linear_solver_.compute_matrix(system_);
    compute_least_squares_gradient(*mesh_, field, boundary_conditions, gradient_);
    assembler_.assemble_rhs(field, gradient_, boundary_conditions, face_flux, system_);

    SteadyScalarTransportResult result;
    for (Index iteration = 0; iteration < options_.maximum_iterations; ++iteration)
    {
        std::copy(field.values().begin(), field.values().end(), previous_field_.values().begin());
        result.linear_solve = linear_solver_.solve(system_.rhs(), field.values());
        if (!result.linear_solve.converged)
        {
            std::ostringstream message;
            message << "Scalar transport inner BiCGSTAB did not converge at outer iteration " << iteration + 1
                    << "; iterations=" << result.linear_solve.iteration_count
                    << ", estimated error=" << result.linear_solve.estimated_relative_error << '.';
            throw std::runtime_error(message.str());
        }

        result.iteration_count = iteration + 1;
        result.field_relative_change = relative_difference(field.values(), previous_field_.values());
        compute_least_squares_gradient(*mesh_, field, boundary_conditions, gradient_);
        assembler_.assemble_rhs(field, gradient_, boundary_conditions, face_flux, system_);
        system_.apply_matrix(field.values(), matrix_product_);
        result.discrete_relative_residual = relative_difference(matrix_product_, system_.rhs());
        result.converged = result.field_relative_change <= options_.field_relative_tolerance &&
                           result.discrete_relative_residual <= options_.discrete_relative_tolerance;
        if (result.converged)
        {
            return result;
        }
    }
    return result;
}

} // namespace cfd
