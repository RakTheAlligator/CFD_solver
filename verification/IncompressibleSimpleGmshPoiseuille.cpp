#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVelocityField.hpp"
#include "cfd/field/FaceFluxField.hpp"
#include "cfd/field/PressureCorrectionBoundaryConditions.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/mesh/MeshStatistics.hpp"
#include "cfd/mesh/Types.hpp"
#include "cfd/meshing/GmshMesher.hpp"
#include "cfd/meshing/RectangleGeometry.hpp"
#include "cfd/numerics/IncompressibleSimpleSolver.hpp"
#include "cfd/numerics/ScalarConvectionOperator.hpp"

#include "support/VerificationStatistics.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace
{

using cfd::verification::ErrorAccumulator;

constexpr double domain_length{4.0};
constexpr double domain_height{1.0};
constexpr double inlet_pressure{0.04};
constexpr double outlet_pressure{0.0};
constexpr double dynamic_viscosity{0.1};
constexpr double target_mesh_size{0.1};
constexpr double baseline_density{1.0};
constexpr double alternate_density{1.225};
constexpr double baseline_pressure_relaxation{0.1};
constexpr double baseline_flux_relaxation{0.3};
constexpr double alternate_pressure_relaxation{0.3};
constexpr double alternate_flux_relaxation{0.1};
constexpr double radians_to_degrees{180.0 / std::numbers::pi};
constexpr double corrected_continuity_tolerance{1.0e-10};
constexpr double maximum_mass_imbalance_tolerance{1.0e-12};
constexpr double relative_flow_mismatch_tolerance{1.0e-10};
constexpr double maximum_v_to_u_scale{5.0e-2};
constexpr double quad_relative_u_error_tolerance{2.0e-2};
constexpr double quad_relative_pressure_error_tolerance{2.0e-2};
constexpr double quad_relative_flow_error_tolerance{3.0e-2};
constexpr double tri_relative_u_error_tolerance{1.5e-2};
constexpr double tri_relative_pressure_error_tolerance{1.5e-2};
constexpr double tri_relative_flow_error_tolerance{2.0e-2};
constexpr double density_independence_tolerance{5.0e-5};
constexpr double relaxation_independence_tolerance{1.0e-8};

[[nodiscard]]
double pressure_gradient_magnitude() noexcept
{
    return (inlet_pressure - outlet_pressure) / domain_length;
}

[[nodiscard]]
double analytical_u(const double y) noexcept
{
    return pressure_gradient_magnitude() * y * (domain_height - y) / (2.0 * dynamic_viscosity);
}

[[nodiscard]]
double analytical_pressure(const double x) noexcept
{
    return inlet_pressure - pressure_gradient_magnitude() * x;
}

[[nodiscard]]
double analytical_mean_velocity() noexcept
{
    return pressure_gradient_magnitude() * domain_height * domain_height / (12.0 * dynamic_viscosity);
}

[[nodiscard]]
double analytical_maximum_velocity() noexcept
{
    return pressure_gradient_magnitude() * domain_height * domain_height / (8.0 * dynamic_viscosity);
}

[[nodiscard]]
double analytical_volumetric_flow() noexcept
{
    return analytical_mean_velocity() * domain_height;
}

struct RunOptions
{
    double density;
    double pressure_relaxation_factor;
    double flux_relaxation_factor;
};

struct MeshDiagnostics
{
    cfd::Index node_count{};
    cfd::Index cell_count{};
    cfd::Index face_count{};
    double minimum_cell_area{};
    double minimum_cell_quality{};
    double maximum_non_orthogonality_degrees{};
};

struct RunDiagnostics
{
    bool converged{};
    cfd::Index iteration_count{};
    double velocity_relative_change{};
    double provisional_continuity{};
    double corrected_continuity{};
    double relative_l2_u_error{};
    double relative_l2_pressure_error{};
    double maximum_absolute_v{};
    double maximum_mass_imbalance{};
    double inlet_mass_flow{};
    double outlet_mass_flow{};
    double relative_flow_mismatch{};
    double volumetric_flow{};
    double relative_volumetric_flow_error{};
};

struct CoupledRun
{
    RunOptions options;
    RunDiagnostics diagnostics;
    cfd::CellVelocityField velocity;
    cfd::CellScalarField pressure;
    cfd::FaceFluxField mass_flux;
};

struct FieldComparison
{
    double relative_u_difference{};
    double scaled_v_rms_difference{};
    double relative_pressure_difference{};
    double relative_face_flux_difference{};
    double relative_volumetric_flow_difference{};
};

struct DensityComparison
{
    FieldComparison fields;
    double mass_flow_ratio{};
    double expected_density_ratio{};
};

[[nodiscard]]
cfd::IncompressibleSimpleOptions simple_options(const RunOptions &options)
{
    return {
        .maximum_iterations = 4000,
        .momentum_relaxation_factor = 1.0,
        .pressure_relaxation_factor = options.pressure_relaxation_factor,
        .rhie_chow_flux_relaxation_factor = options.flux_relaxation_factor,
        .velocity_relative_tolerance = 1.0e-10,
        .continuity_relative_tolerance = 1.0e-10,
        .momentum_linear_solver = {.relative_tolerance = 1.0e-12, .maximum_iterations = 5000},
        .pressure_correction_linear_solver = {.relative_tolerance = 1.0e-12, .maximum_iterations = 5000},
    };
}

[[nodiscard]]
cfd::ScalarBoundaryConditions velocity_boundary_conditions(const cfd::Mesh &mesh)
{
    std::vector<cfd::ScalarBoundaryCondition> conditions;
    conditions.reserve(mesh.boundary_groups().size());
    for (const cfd::BoundaryGroup &group : mesh.boundary_groups())
    {
        if (group.name == "wall")
        {
            conditions.emplace_back(cfd::ScalarBoundaryConditionType::Dirichlet, 0.0);
        }
        else if (group.name == "inlet" || group.name == "outlet")
        {
            conditions.emplace_back(cfd::ScalarBoundaryConditionType::Neumann, 0.0);
        }
        else
        {
            throw std::runtime_error("Gmsh Poiseuille mesh has an unexpected boundary group.");
        }
    }
    return {mesh.boundary_groups().size(), std::move(conditions)};
}

[[nodiscard]]
cfd::ScalarBoundaryConditions pressure_boundary_conditions(const cfd::Mesh &mesh)
{
    std::vector<cfd::ScalarBoundaryCondition> conditions;
    conditions.reserve(mesh.boundary_groups().size());
    for (const cfd::BoundaryGroup &group : mesh.boundary_groups())
    {
        if (group.name == "inlet")
        {
            conditions.emplace_back(cfd::ScalarBoundaryConditionType::Dirichlet, inlet_pressure);
        }
        else if (group.name == "outlet")
        {
            conditions.emplace_back(cfd::ScalarBoundaryConditionType::Dirichlet, outlet_pressure);
        }
        else if (group.name == "wall")
        {
            conditions.emplace_back(cfd::ScalarBoundaryConditionType::Neumann, 0.0);
        }
        else
        {
            throw std::runtime_error("Gmsh Poiseuille mesh has an unexpected boundary group.");
        }
    }
    return {mesh.boundary_groups().size(), std::move(conditions)};
}

[[nodiscard]]
cfd::PressureCorrectionBoundaryConditions pressure_correction_boundary_conditions(const cfd::Mesh &mesh)
{
    std::vector<cfd::PressureCorrectionBoundaryConditionType> conditions;
    conditions.reserve(mesh.boundary_groups().size());
    for (const cfd::BoundaryGroup &group : mesh.boundary_groups())
    {
        if (group.name == "inlet" || group.name == "outlet")
        {
            conditions.push_back(cfd::PressureCorrectionBoundaryConditionType::FixedPressure);
        }
        else if (group.name == "wall")
        {
            conditions.push_back(cfd::PressureCorrectionBoundaryConditionType::FixedMassFlux);
        }
        else
        {
            throw std::runtime_error("Gmsh Poiseuille mesh has an unexpected boundary group.");
        }
    }
    return {mesh.boundary_groups().size(), std::move(conditions)};
}

[[nodiscard]]
// Area-weighted relative L2 difference for two cell fields on the same Mesh.
double relative_difference(const std::span<const double> reference, const std::span<const double> candidate,
                           const std::span<const double> weights)
{
    if (reference.size() != candidate.size() || reference.size() != weights.size())
    {
        throw std::runtime_error("Cannot compare fields with different cardinalities.");
    }

    double squared_difference{};
    double squared_reference{};
    for (cfd::Index value_id = 0; value_id < reference.size(); ++value_id)
    {
        const double difference{candidate[value_id] - reference[value_id]};
        squared_difference += weights[value_id] * difference * difference;
        squared_reference += weights[value_id] * reference[value_id] * reference[value_id];
    }
    if (squared_reference == 0.0)
    {
        return squared_difference == 0.0 ? 0.0 : std::numeric_limits<double>::infinity();
    }
    return std::sqrt(squared_difference / squared_reference);
}

[[nodiscard]]
// Area-weighted RMS difference normalized by a nonzero physical field scale.
double scaled_rms_difference(const std::span<const double> reference, const std::span<const double> candidate,
                             const std::span<const double> weights, const double scale)
{
    if (reference.size() != candidate.size() || reference.size() != weights.size())
    {
        throw std::runtime_error("Cannot compare fields with different cardinalities.");
    }

    double total_weight{};
    double squared_difference{};
    for (cfd::Index value_id = 0; value_id < reference.size(); ++value_id)
    {
        const double difference{candidate[value_id] - reference[value_id]};
        total_weight += weights[value_id];
        squared_difference += weights[value_id] * difference * difference;
    }
    if (!(total_weight > 0.0) || !(scale > 0.0))
    {
        throw std::runtime_error("Cannot normalize a field difference with a non-positive scale.");
    }
    return std::sqrt(squared_difference / total_weight) / scale;
}

[[nodiscard]]
// Relative L2 difference of corresponding integrated face fluxes on one Mesh.
double relative_unweighted_difference(const std::span<const double> reference, const std::span<const double> candidate)
{
    if (reference.size() != candidate.size())
    {
        throw std::runtime_error("Cannot compare face fields with different cardinalities.");
    }

    double squared_difference{};
    double squared_reference{};
    for (cfd::Index value_id = 0; value_id < reference.size(); ++value_id)
    {
        const double difference{candidate[value_id] - reference[value_id]};
        squared_difference += difference * difference;
        squared_reference += reference[value_id] * reference[value_id];
    }
    if (squared_reference == 0.0)
    {
        return squared_difference == 0.0 ? 0.0 : std::numeric_limits<double>::infinity();
    }
    return std::sqrt(squared_difference / squared_reference);
}

[[nodiscard]]
double relative_scalar_difference(const double reference, const double candidate) noexcept
{
    const double difference{std::abs(candidate - reference)};
    return reference == 0.0 ? (difference == 0.0 ? 0.0 : std::numeric_limits<double>::infinity())
                            : difference / std::abs(reference);
}

[[nodiscard]]
double relative_flow_mismatch(const double inlet_mass_flow, const double outlet_mass_flow) noexcept
{
    const double scale{std::max(std::abs(inlet_mass_flow), std::abs(outlet_mass_flow))};
    const double difference{std::abs(inlet_mass_flow - outlet_mass_flow)};
    return scale == 0.0 ? (difference == 0.0 ? 0.0 : std::numeric_limits<double>::infinity()) : difference / scale;
}

[[nodiscard]]
MeshDiagnostics compute_mesh_diagnostics(const cfd::Mesh &mesh, const cfd::CellType expected_cell_type)
{
    if (mesh.cell_count() == 0 ||
        std::find(mesh.cell_types().begin(), mesh.cell_types().end(), expected_cell_type) == mesh.cell_types().end() ||
        std::find_if(mesh.cell_types().begin(), mesh.cell_types().end(),
                     [expected_cell_type](const cfd::CellType type) { return type != expected_cell_type; }) !=
            mesh.cell_types().end())
    {
        throw std::runtime_error("Gmsh Poiseuille mesh does not contain the requested cell type exclusively.");
    }

    const cfd::MeshStatistics statistics{cfd::compute_mesh_statistics(mesh)};
    double maximum_non_orthogonality{};
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        const cfd::FaceAdjacency &adjacency{mesh.face_adjacencies()[face_id]};
        if (adjacency.is_boundary())
        {
            continue;
        }
        const cfd::Point2 &owner_center{mesh.cell_centers()[adjacency.owner]};
        const cfd::Point2 &neighbor_center{mesh.cell_centers()[adjacency.neighbor]};
        const cfd::Vector2 displacement{
            neighbor_center.x - owner_center.x,
            neighbor_center.y - owner_center.y,
        };
        const cfd::Vector2 &area_vector{mesh.face_area_vectors()[face_id]};
        const double normalized_projection{(area_vector.x * displacement.x + area_vector.y * displacement.y) /
                                           (mesh.face_lengths()[face_id] * std::hypot(displacement.x, displacement.y))};
        maximum_non_orthogonality = std::max(
            maximum_non_orthogonality, std::acos(std::clamp(normalized_projection, -1.0, 1.0)) * radians_to_degrees);
    }

    if (!(statistics.cell_areas.minimum > 0.0) || !std::isfinite(statistics.cell_quality.minimum) ||
        !std::isfinite(maximum_non_orthogonality))
    {
        throw std::runtime_error("Gmsh Poiseuille mesh has invalid geometry diagnostics.");
    }
    return {
        .node_count = mesh.node_count(),
        .cell_count = mesh.cell_count(),
        .face_count = mesh.face_count(),
        .minimum_cell_area = statistics.cell_areas.minimum,
        .minimum_cell_quality = statistics.cell_quality.minimum,
        .maximum_non_orthogonality_degrees = maximum_non_orthogonality,
    };
}

[[nodiscard]]
RunDiagnostics compute_run_diagnostics(const cfd::Mesh &mesh, const double density,
                                       const cfd::IncompressibleSimpleResult &result,
                                       const cfd::CellVelocityField &velocity, const cfd::CellScalarField &pressure,
                                       const cfd::FaceFluxField &mass_flux)
{
    ErrorAccumulator u_error;
    ErrorAccumulator u_reference;
    ErrorAccumulator pressure_error;
    ErrorAccumulator pressure_reference;
    double maximum_absolute_v{};
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        const double area{mesh.cell_areas()[cell_id]};
        const double exact_u{analytical_u(mesh.cell_centers()[cell_id].y)};
        const double exact_p{analytical_pressure(mesh.cell_centers()[cell_id].x)};
        u_error.add(area, velocity.u()[cell_id] - exact_u);
        u_reference.add(area, exact_u);
        pressure_error.add(area, pressure[cell_id] - exact_p);
        pressure_reference.add(area, exact_p);
        maximum_absolute_v = std::max(maximum_absolute_v, std::abs(velocity.v()[cell_id]));
    }

    double inlet_mass_flow{};
    double outlet_mass_flow{};
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (!mesh.face_adjacencies()[face_id].is_boundary())
        {
            continue;
        }
        const std::string &boundary_name{mesh.boundary_groups()[mesh.face_boundary_ids()[face_id]].name};
        if (boundary_name == "inlet")
        {
            inlet_mass_flow -= mass_flux[face_id];
        }
        else if (boundary_name == "outlet")
        {
            outlet_mass_flow += mass_flux[face_id];
        }
        else if (boundary_name == "wall" && mass_flux[face_id] != 0.0)
        {
            throw std::runtime_error("Gmsh Poiseuille wall FixedMassFlux value changed from zero.");
        }
    }

    // FaceFluxField contains mass flow; average the conservative inlet/outlet
    // magnitudes before dividing by density to obtain volumetric flow.
    const double volumetric_flow{0.5 * (inlet_mass_flow + outlet_mass_flow) / density};
    const double u_reference_rms{u_reference.finish().area_weighted_rms_error};
    const double pressure_reference_rms{pressure_reference.finish().area_weighted_rms_error};
    return {
        .converged = result.converged,
        .iteration_count = result.iteration_count,
        .velocity_relative_change = result.velocity_relative_change,
        .provisional_continuity = result.provisional_continuity_relative_residual,
        .corrected_continuity = result.continuity_relative_residual,
        .relative_l2_u_error = u_error.finish().area_weighted_rms_error / u_reference_rms,
        .relative_l2_pressure_error = pressure_error.finish().area_weighted_rms_error / pressure_reference_rms,
        .maximum_absolute_v = maximum_absolute_v,
        .maximum_mass_imbalance = result.maximum_mass_imbalance,
        .inlet_mass_flow = inlet_mass_flow,
        .outlet_mass_flow = outlet_mass_flow,
        .relative_flow_mismatch = relative_flow_mismatch(inlet_mass_flow, outlet_mass_flow),
        .volumetric_flow = volumetric_flow,
        .relative_volumetric_flow_error = relative_scalar_difference(analytical_volumetric_flow(), volumetric_flow),
    };
}

[[nodiscard]]
CoupledRun solve_case(const cfd::Mesh &mesh, const RunOptions &run_options, const std::string_view case_name)
{
    cfd::CellVelocityField velocity{mesh.cell_count()};
    cfd::CellScalarField pressure{mesh.cell_count()};
    cfd::FaceFluxField mass_flux{mesh.face_count()};
    const cfd::ScalarBoundaryConditions velocity_conditions{velocity_boundary_conditions(mesh)};
    const cfd::ScalarBoundaryConditions pressure_conditions{pressure_boundary_conditions(mesh)};
    const cfd::PressureCorrectionBoundaryConditions pressure_correction_conditions{
        pressure_correction_boundary_conditions(mesh)};
    const cfd::IncompressibleSimpleOptions options{simple_options(run_options)};
    cfd::IncompressibleSimpleSolver solver{mesh, run_options.density, dynamic_viscosity,
                                           cfd::ScalarConvectionScheme::Linear, options};
    cfd::SimpleIterationInfo last_iteration;
    try
    {
        const cfd::IncompressibleSimpleResult result{solver.solve(
            velocity_conditions, velocity_conditions, pressure_conditions, pressure_correction_conditions, velocity,
            pressure, mass_flux, [&last_iteration](const cfd::SimpleIterationInfo &info) { last_iteration = info; })};
        const RunDiagnostics diagnostics{
            compute_run_diagnostics(mesh, run_options.density, result, velocity, pressure, mass_flux)};
        return {
            .options = run_options,
            .diagnostics = diagnostics,
            .velocity = std::move(velocity),
            .pressure = std::move(pressure),
            .mass_flux = std::move(mass_flux),
        };
    }
    catch (const std::exception &error)
    {
        std::ostringstream message;
        message << case_name << " failed after " << last_iteration.iteration << " completed SIMPLE iterations";
        if (last_iteration.iteration != 0)
        {
            message << " (r_U=" << last_iteration.velocity_relative_change
                    << ", r_cont*=" << last_iteration.provisional_continuity_relative_residual
                    << ", max|p'|=" << last_iteration.maximum_pressure_correction << ')';
        }
        message << ": " << error.what();
        throw std::runtime_error(message.str());
    }
}

[[nodiscard]]
FieldComparison compare_fields(const cfd::Mesh &mesh, const CoupledRun &reference, const CoupledRun &candidate)
{
    return {
        .relative_u_difference =
            relative_difference(reference.velocity.u().values(), candidate.velocity.u().values(), mesh.cell_areas()),
        .scaled_v_rms_difference =
            scaled_rms_difference(reference.velocity.v().values(), candidate.velocity.v().values(), mesh.cell_areas(),
                                  analytical_maximum_velocity()),
        .relative_pressure_difference =
            relative_difference(reference.pressure.values(), candidate.pressure.values(), mesh.cell_areas()),
        .relative_face_flux_difference =
            relative_unweighted_difference(reference.mass_flux.values(), candidate.mass_flux.values()),
        .relative_volumetric_flow_difference =
            relative_scalar_difference(reference.diagnostics.volumetric_flow, candidate.diagnostics.volumetric_flow),
    };
}

[[nodiscard]]
DensityComparison compare_density_runs(const cfd::Mesh &mesh, const CoupledRun &reference, const CoupledRun &candidate)
{
    const double reference_mass_flow{0.5 *
                                     (reference.diagnostics.inlet_mass_flow + reference.diagnostics.outlet_mass_flow)};
    const double candidate_mass_flow{0.5 *
                                     (candidate.diagnostics.inlet_mass_flow + candidate.diagnostics.outlet_mass_flow)};
    return {
        .fields = compare_fields(mesh, reference, candidate),
        .mass_flow_ratio = reference_mass_flow == 0.0 ? std::numeric_limits<double>::infinity()
                                                      : candidate_mass_flow / reference_mass_flow,
        .expected_density_ratio = candidate.options.density / reference.options.density,
    };
}

void require_at_most(const std::string_view case_name, const std::string_view quantity, const double value,
                     const double threshold)
{
    if (!std::isfinite(value) || value > threshold)
    {
        std::ostringstream message;
        message << case_name << ' ' << quantity << " = " << value << " exceeds threshold " << threshold << '.';
        throw std::runtime_error(message.str());
    }
}

void validate_run(const std::string_view case_name, const CoupledRun &run, const double relative_u_error_tolerance,
                  const double relative_pressure_error_tolerance, const double relative_flow_error_tolerance)
{
    const RunDiagnostics &result{run.diagnostics};
    const cfd::IncompressibleSimpleOptions options{simple_options(run.options)};
    if (!result.converged)
    {
        throw std::runtime_error(std::string{case_name} + " did not converge within the SIMPLE iteration limit.");
    }
    if (!(result.inlet_mass_flow > 0.0) || !(result.outlet_mass_flow > 0.0))
    {
        throw std::runtime_error(std::string{case_name} + " did not produce positive inlet/outlet flow magnitudes.");
    }

    require_at_most(case_name, "velocity relative change", result.velocity_relative_change,
                    options.velocity_relative_tolerance);
    require_at_most(case_name, "provisional continuity residual", result.provisional_continuity,
                    options.continuity_relative_tolerance);
    require_at_most(case_name, "corrected continuity residual", result.corrected_continuity,
                    corrected_continuity_tolerance);
    require_at_most(case_name, "maximum mass imbalance", result.maximum_mass_imbalance,
                    maximum_mass_imbalance_tolerance);
    require_at_most(case_name, "relative L2 u error", result.relative_l2_u_error, relative_u_error_tolerance);
    require_at_most(case_name, "relative L2 pressure error", result.relative_l2_pressure_error,
                    relative_pressure_error_tolerance);
    require_at_most(case_name, "max |v| / analytical u_max", result.maximum_absolute_v / analytical_maximum_velocity(),
                    maximum_v_to_u_scale);
    require_at_most(case_name, "relative inlet/outlet flow mismatch", result.relative_flow_mismatch,
                    relative_flow_mismatch_tolerance);
    require_at_most(case_name, "relative volumetric-flow error", result.relative_volumetric_flow_error,
                    relative_flow_error_tolerance);
}

void validate_density_comparison(const DensityComparison &comparison)
{
    require_at_most("density comparison", "relative u difference", comparison.fields.relative_u_difference,
                    density_independence_tolerance);
    require_at_most("density comparison", "scaled v RMS difference", comparison.fields.scaled_v_rms_difference,
                    density_independence_tolerance);
    require_at_most("density comparison", "relative pressure difference",
                    comparison.fields.relative_pressure_difference, density_independence_tolerance);
    require_at_most("density comparison", "relative volumetric-flow difference",
                    comparison.fields.relative_volumetric_flow_difference, density_independence_tolerance);
    require_at_most("density comparison", "relative mass-flow ratio error",
                    relative_scalar_difference(comparison.expected_density_ratio, comparison.mass_flow_ratio),
                    density_independence_tolerance);
}

void validate_relaxation_comparison(const FieldComparison &comparison)
{
    require_at_most("relaxation comparison", "relative u difference", comparison.relative_u_difference,
                    relaxation_independence_tolerance);
    require_at_most("relaxation comparison", "scaled v RMS difference", comparison.scaled_v_rms_difference,
                    relaxation_independence_tolerance);
    require_at_most("relaxation comparison", "relative pressure difference", comparison.relative_pressure_difference,
                    relaxation_independence_tolerance);
    require_at_most("relaxation comparison", "relative face-flux difference", comparison.relative_face_flux_difference,
                    relaxation_independence_tolerance);
    require_at_most("relaxation comparison", "relative volumetric-flow difference",
                    comparison.relative_volumetric_flow_difference, relaxation_independence_tolerance);
}

void print_mesh(const std::string_view name, const MeshDiagnostics &mesh)
{
    std::cout << name << " mesh: nodes=" << mesh.node_count << ", cells=" << mesh.cell_count
              << ", faces=" << mesh.face_count << std::scientific << std::setprecision(6)
              << ", min_area=" << mesh.minimum_cell_area << ", min_quality=" << mesh.minimum_cell_quality
              << ", max_nonorth_deg=" << mesh.maximum_non_orthogonality_degrees << '\n';
}

void print_run_table(const CoupledRun &quad, const CoupledRun &tri, const CoupledRun &dense,
                     const CoupledRun &alternate)
{
    std::cout << "\nCoupled cold-start SIMPLE results\n"
              << std::left << std::setw(17) << "case" << std::right << std::setw(8) << "iters" << std::setw(13) << "r_U"
              << std::setw(13) << "r_cont*" << std::setw(13) << "r_cont" << std::setw(13) << "relL2(u)" << std::setw(13)
              << "relL2(p)" << std::setw(13) << "max|v|" << std::setw(13) << "flow_err" << '\n';
    const auto print_row = [](const std::string_view name, const CoupledRun &run) {
        const RunDiagnostics &result{run.diagnostics};
        std::cout << std::left << std::setw(17) << name << std::right << std::setw(8) << result.iteration_count
                  << std::scientific << std::setprecision(5) << std::setw(13) << result.velocity_relative_change
                  << std::setw(13) << result.provisional_continuity << std::setw(13) << result.corrected_continuity
                  << std::setw(13) << result.relative_l2_u_error << std::setw(13) << result.relative_l2_pressure_error
                  << std::setw(13) << result.maximum_absolute_v << std::setw(13)
                  << result.relative_volumetric_flow_error << '\n';
    };
    print_row("QUAD rho=1", quad);
    print_row("TRI rho=1", tri);
    print_row("QUAD rho=1.225", dense);
    print_row("QUAD alternate", alternate);
}

void print_flow(const std::string_view name, const CoupledRun &run)
{
    const RunDiagnostics &result{run.diagnostics};
    std::cout << name << " flow: inlet_mass=" << std::scientific << std::setprecision(8) << result.inlet_mass_flow
              << ", outlet_mass=" << result.outlet_mass_flow << ", Q=" << result.volumetric_flow
              << ", mismatch=" << result.relative_flow_mismatch << ", max|Rm|=" << result.maximum_mass_imbalance
              << '\n';
}

void print_comparisons(const CoupledRun &quad, const CoupledRun &dense, const DensityComparison &density_comparison,
                       const CoupledRun &alternate, const FieldComparison &relaxation_comparison)
{
    std::cout << "\nDensity comparison (rho=" << quad.options.density << " vs " << dense.options.density << ")\n"
              << std::scientific << std::setprecision(8) << "  iterations: " << quad.diagnostics.iteration_count
              << " vs " << dense.diagnostics.iteration_count
              << ", du_rel=" << density_comparison.fields.relative_u_difference
              << ", dv_rms/Umax=" << density_comparison.fields.scaled_v_rms_difference
              << ", dp_rel=" << density_comparison.fields.relative_pressure_difference
              << ", dQ_rel=" << density_comparison.fields.relative_volumetric_flow_difference
              << ", mass_ratio=" << density_comparison.mass_flow_ratio
              << ", expected=" << density_comparison.expected_density_ratio << ", ratio_rel_err="
              << relative_scalar_difference(density_comparison.expected_density_ratio,
                                            density_comparison.mass_flow_ratio)
              << '\n'
              << "\nRelaxation comparison\n"
              << "  baseline (alpha_p=" << quad.options.pressure_relaxation_factor
              << ", alpha_F=" << quad.options.flux_relaxation_factor << "): " << quad.diagnostics.iteration_count
              << " iterations\n"
              << "  alternate (alpha_p=" << alternate.options.pressure_relaxation_factor
              << ", alpha_F=" << alternate.options.flux_relaxation_factor
              << "): " << alternate.diagnostics.iteration_count << " iterations\n"
              << "  du_rel=" << relaxation_comparison.relative_u_difference
              << ", dv_rms/Umax=" << relaxation_comparison.scaled_v_rms_difference
              << ", dp_rel=" << relaxation_comparison.relative_pressure_difference
              << ", dF_rel=" << relaxation_comparison.relative_face_flux_difference
              << ", dQ_rel=" << relaxation_comparison.relative_volumetric_flow_difference << '\n';
}

void print_acceptance_thresholds()
{
    std::cout << "\nAcceptance thresholds\n"
              << std::scientific << std::setprecision(3)
              << "  common: r_U<=1e-10, r_cont*<=1e-10, r_cont<=" << corrected_continuity_tolerance
              << ", max|Rm|<=" << maximum_mass_imbalance_tolerance
              << ", flow_mismatch<=" << relative_flow_mismatch_tolerance << ", max|v|/u_max<=" << maximum_v_to_u_scale
              << '\n'
              << "  QUAD: relL2(u)<=" << quad_relative_u_error_tolerance
              << ", relL2(p)<=" << quad_relative_pressure_error_tolerance
              << ", flow_error<=" << quad_relative_flow_error_tolerance << '\n'
              << "  TRI:  relL2(u)<=" << tri_relative_u_error_tolerance
              << ", relL2(p)<=" << tri_relative_pressure_error_tolerance
              << ", flow_error<=" << tri_relative_flow_error_tolerance << '\n'
              << "  density field/flow differences and mass-ratio error<=" << density_independence_tolerance << '\n'
              << "  relaxation field/flux/flow differences<=" << relaxation_independence_tolerance << '\n';
}

} // namespace

int main()
{
    try
    {
        cfd::MeshBuildResult quad_build{cfd::build_mesh(
            cfd::generate_mesh({.length = domain_length, .height = domain_height},
                               {.mesh_size = target_mesh_size, .cell_type = cfd::CellType::Quadrilateral}))};
        cfd::MeshBuildResult tri_build{
            cfd::build_mesh(cfd::generate_mesh({.length = domain_length, .height = domain_height},
                                               {.mesh_size = target_mesh_size, .cell_type = cfd::CellType::Triangle}))};
        const cfd::Mesh &quad_mesh{quad_build.mesh};
        const cfd::Mesh &tri_mesh{tri_build.mesh};
        const MeshDiagnostics quad_mesh_diagnostics{compute_mesh_diagnostics(quad_mesh, cfd::CellType::Quadrilateral)};
        const MeshDiagnostics tri_mesh_diagnostics{compute_mesh_diagnostics(tri_mesh, cfd::CellType::Triangle)};

        print_mesh("QUAD", quad_mesh_diagnostics);
        print_mesh("TRI", tri_mesh_diagnostics);

        const RunOptions baseline_options{
            .density = baseline_density,
            .pressure_relaxation_factor = baseline_pressure_relaxation,
            .flux_relaxation_factor = baseline_flux_relaxation,
        };
        const RunOptions dense_options{
            .density = alternate_density,
            .pressure_relaxation_factor = baseline_pressure_relaxation,
            .flux_relaxation_factor = baseline_flux_relaxation,
        };
        const RunOptions tri_options{
            .density = baseline_density,
            .pressure_relaxation_factor = baseline_pressure_relaxation,
            // The TRI cold start needs more conservative provisional-flux relaxation.
            .flux_relaxation_factor = alternate_flux_relaxation,
        };
        const RunOptions alternate_relaxation_options{
            .density = baseline_density,
            .pressure_relaxation_factor = alternate_pressure_relaxation,
            .flux_relaxation_factor = alternate_flux_relaxation,
        };

        std::cout << "Running TRI rho=1 (alpha_p=0.1, alpha_F=0.1)...\n" << std::flush;
        const CoupledRun tri_run{solve_case(tri_mesh, tri_options, "TRI rho=1")};
        std::cout << "Running QUAD rho=1 baseline...\n" << std::flush;
        const CoupledRun quad_run{solve_case(quad_mesh, baseline_options, "QUAD rho=1 baseline")};
        std::cout << "Running QUAD rho=1.225...\n" << std::flush;
        const CoupledRun dense_run{solve_case(quad_mesh, dense_options, "QUAD rho=1.225")};
        std::cout << "Running QUAD alternate relaxation...\n" << std::flush;
        const CoupledRun alternate_run{
            solve_case(quad_mesh, alternate_relaxation_options, "QUAD alternate relaxation")};
        const DensityComparison density_comparison{compare_density_runs(quad_mesh, quad_run, dense_run)};
        const FieldComparison relaxation_comparison{compare_fields(quad_mesh, quad_run, alternate_run)};

        std::cout << "Gmsh coupled SIMPLE Poiseuille verification\n"
                  << "L=" << domain_length << ", H=" << domain_height << ", mu=" << dynamic_viscosity
                  << ", p_in=" << inlet_pressure << ", p_out=" << outlet_pressure
                  << ", Q_exact=" << analytical_volumetric_flow() << "\n\n";
        print_run_table(quad_run, tri_run, dense_run, alternate_run);
        print_flow("QUAD rho=1", quad_run);
        print_flow("TRI rho=1", tri_run);
        print_flow("QUAD rho=1.225", dense_run);
        print_comparisons(quad_run, dense_run, density_comparison, alternate_run, relaxation_comparison);
        print_acceptance_thresholds();

        validate_run("QUAD rho=1", quad_run, quad_relative_u_error_tolerance, quad_relative_pressure_error_tolerance,
                     quad_relative_flow_error_tolerance);
        validate_run("TRI rho=1", tri_run, tri_relative_u_error_tolerance, tri_relative_pressure_error_tolerance,
                     tri_relative_flow_error_tolerance);
        validate_run("QUAD rho=1.225", dense_run, quad_relative_u_error_tolerance,
                     quad_relative_pressure_error_tolerance, quad_relative_flow_error_tolerance);
        validate_run("QUAD alternate relaxation", alternate_run, quad_relative_u_error_tolerance,
                     quad_relative_pressure_error_tolerance, quad_relative_flow_error_tolerance);
        validate_density_comparison(density_comparison);
        validate_relaxation_comparison(relaxation_comparison);
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "Gmsh coupled SIMPLE Poiseuille verification failed: " << error.what() << '\n';
        return 1;
    }
}
