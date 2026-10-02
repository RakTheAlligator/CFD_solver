#include "cfd/numerics/SteadyScalarTransportSolver.hpp"

#include "cfd/field/FaceFluxField.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/meshing/RawMeshData.hpp"
#include "cfd/numerics/LeastSquaresGradient.hpp"

#include "support/VerificationStatistics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace
{

constexpr double diffusion_coefficient{0.25};
constexpr double carrier_speed{1.0};
constexpr double domain_length{1.0};
constexpr double domain_height{1.0};
constexpr std::array<cfd::Index, 3> grid_levels{16, 32, 64};

struct SchemeStudy
{
    const char *name;
    cfd::ScalarConvectionScheme scheme;
    cfd::ScalarConvectionLimiter limiter;
};

constexpr std::array studies{
    SchemeStudy{"Upwind", cfd::ScalarConvectionScheme::FirstOrderUpwind, cfd::ScalarConvectionLimiter::None},
    SchemeStudy{"Linear", cfd::ScalarConvectionScheme::Linear, cfd::ScalarConvectionLimiter::None},
    SchemeStudy{"LinearUpwind", cfd::ScalarConvectionScheme::LinearUpwind, cfd::ScalarConvectionLimiter::None},
    SchemeStudy{"LinearUpwind+BJ", cfd::ScalarConvectionScheme::LinearUpwind,
                cfd::ScalarConvectionLimiter::BarthJespersen}};

[[nodiscard]]
cfd::SteadyScalarTransportOptions verification_options()
{
    cfd::SteadyScalarTransportOptions options;
    options.maximum_iterations = 1000;
    options.field_relative_tolerance = 1.0e-11;
    options.discrete_relative_tolerance = 1.0e-11;
    options.linear_solver = {1.0e-13, 5000};
    return options;
}

void require(const bool valid, const char *message)
{
    if (!valid)
    {
        throw std::runtime_error(message);
    }
}

[[nodiscard]]
cfd::RawMeshData make_raw_mesh(const cfd::Index nx, const cfd::Index ny, const double shear = 0.0,
                               const bool triangles = false)
{
    cfd::RawMeshData raw;
    const auto node_id = [nx](const cfd::Index i, const cfd::Index j) { return j * (nx + 1) + i; };
    for (cfd::Index j = 0; j <= ny; ++j)
    {
        const double y{domain_height * static_cast<double>(j) / static_cast<double>(ny)};
        for (cfd::Index i = 0; i <= nx; ++i)
        {
            raw.nodes.push_back({domain_length * static_cast<double>(i) / static_cast<double>(nx) + shear * y, y});
        }
    }
    raw.cell_node_offsets.push_back(0);
    for (cfd::Index j = 0; j < ny; ++j)
    {
        for (cfd::Index i = 0; i < nx; ++i)
        {
            const cfd::Index bottom_left{node_id(i, j)};
            const cfd::Index bottom_right{node_id(i + 1, j)};
            const cfd::Index top_right{node_id(i + 1, j + 1)};
            const cfd::Index top_left{node_id(i, j + 1)};
            if (triangles)
            {
                raw.cell_types.push_back(cfd::CellType::Triangle);
                raw.cell_nodes.insert(raw.cell_nodes.end(), {bottom_left, bottom_right, top_right});
                raw.cell_node_offsets.push_back(raw.cell_nodes.size());
                raw.cell_types.push_back(cfd::CellType::Triangle);
                raw.cell_nodes.insert(raw.cell_nodes.end(), {bottom_left, top_right, top_left});
            }
            else
            {
                raw.cell_types.push_back(cfd::CellType::Quadrilateral);
                raw.cell_nodes.insert(raw.cell_nodes.end(), {bottom_left, bottom_right, top_right, top_left});
            }
            raw.cell_node_offsets.push_back(raw.cell_nodes.size());
        }
    }
    // One group per boundary face permits exact affine face-center values
    // without extending the uniform-per-group scalar BC abstraction.
    const auto add_boundary = [&raw](const cfd::Index first, const cfd::Index second) {
        const cfd::BoundaryId id{raw.boundary_groups.size()};
        raw.boundary_groups.push_back({id, "boundary_" + std::to_string(id)});
        raw.boundary_edges.push_back({{first, second}, id});
    };
    for (cfd::Index i = 0; i < nx; ++i)
    {
        add_boundary(node_id(i, 0), node_id(i + 1, 0));
        add_boundary(node_id(i + 1, ny), node_id(i, ny));
    }
    for (cfd::Index j = 0; j < ny; ++j)
    {
        add_boundary(node_id(0, j + 1), node_id(0, j));
        add_boundary(node_id(nx, j), node_id(nx, j + 1));
    }
    return raw;
}

[[nodiscard]]
double affine_value(const cfd::Point2 &point) noexcept
{
    return 2.0 * point.x - 3.0 * point.y + 1.5;
}

[[nodiscard]]
double exponential_value(const cfd::Point2 &point)
{
    const double exponent{carrier_speed / diffusion_coefficient};
    return std::expm1(exponent * point.x) / std::expm1(exponent * domain_length);
}

[[nodiscard]]
cfd::ScalarBoundaryConditions make_conditions(const cfd::Mesh &mesh, const bool affine)
{
    std::vector<cfd::ScalarBoundaryCondition> conditions(mesh.boundary_groups().size(),
                                                         {cfd::ScalarBoundaryConditionType::Neumann, 0.0});
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (!mesh.face_adjacencies()[face_id].is_boundary())
        {
            continue;
        }
        auto &condition{conditions[mesh.face_boundary_ids()[face_id]]};
        const auto &center{mesh.face_centers()[face_id]};
        if (affine)
        {
            condition = {cfd::ScalarBoundaryConditionType::Dirichlet, affine_value(center)};
        }
        else if (mesh.face_area_vectors()[face_id].x != 0.0)
        {
            condition = {cfd::ScalarBoundaryConditionType::Dirichlet, exponential_value(center)};
        }
    }
    return {conditions.size(), std::move(conditions)};
}

struct LevelResult
{
    cfd::verification::ErrorStatistics error;
    cfd::SteadyScalarTransportResult solve;
};

void check_conservation(const cfd::Mesh &mesh, const cfd::ScalarBoundaryConditions &conditions,
                        const cfd::FaceFluxField &flux, const cfd::CellScalarField &field, const SchemeStudy &study)
{
    cfd::CellVectorField gradient{mesh.cell_count()};
    cfd::compute_least_squares_gradient(mesh, field, conditions, gradient);
    const cfd::ScalarConvectionOperator convection{mesh, study.scheme, study.limiter};
    const cfd::ScalarDiffusionOperator diffusion{mesh, diffusion_coefficient};
    cfd::CellScalarField convective{mesh.cell_count()};
    cfd::CellScalarField diffusive{mesh.cell_count()};
    std::vector<double> limiter(mesh.cell_count());
    if (study.scheme == cfd::ScalarConvectionScheme::LinearUpwind)
    {
        convection.compute_flux_balance(field, conditions, flux, gradient, limiter, convective);
    }
    else
    {
        convection.compute_flux_balance(field, conditions, flux, convective, diffusion.face_primary_coefficients());
    }
    diffusion.compute_flux_balance(field, conditions, gradient, diffusive);
    cfd::ScalarTransportAssembler assembler{mesh, diffusion_coefficient, study.scheme, study.limiter};
    cfd::ScalarLinearSystem system{mesh};
    assembler.assemble_matrix(conditions, flux, system);
    assembler.assemble_rhs(field, gradient, conditions, flux, system);
    std::vector<double> product(mesh.cell_count());
    system.apply_matrix(field.values(), product);
    double total_balance{};
    double total_absolute_balance{};
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        const double balance{convective[cell_id] + diffusive[cell_id]};
        require(std::abs(balance - (product[cell_id] - system.rhs()[cell_id])) < 1.0e-11,
                "Direct transport balance disagrees with assembled residual.");
        require(std::abs(balance) < 1.0e-9, "Cell transport balance is not closed.");
        total_balance += balance;
        total_absolute_balance += std::abs(balance);
    }
    require(std::abs(total_balance) < 1.0e-9, "Global transport balance is not closed.");
    std::cout << "  conservation: sum=" << total_balance << ", sum_abs=" << total_absolute_balance;
    if (study.limiter == cfd::ScalarConvectionLimiter::BarthJespersen)
    {
        cfd::Index limited_count{};
        double minimum_limiter{1.0};
        for (const double coefficient : limiter)
        {
            limited_count += static_cast<cfd::Index>(coefficient < 1.0 - 1.0e-12);
            minimum_limiter = std::min(minimum_limiter, coefficient);
        }
        std::cout << ", limited_cells=" << limited_count << ", min_limiter=" << minimum_limiter;
    }
    std::cout << '\n';
}

[[nodiscard]]
LevelResult solve_level(const cfd::Index nx, const cfd::Index ny, const SchemeStudy &study, const bool affine,
                        const double shear = 0.0, const bool triangles = false)
{
    auto build{cfd::build_mesh(make_raw_mesh(nx, ny, shear, triangles))};
    const cfd::Mesh &mesh{build.mesh};
    const auto conditions{make_conditions(mesh, affine)};
    cfd::FaceFluxField flux{mesh.face_count()};
    if (!affine)
    {
        for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
        {
            flux[face_id] = carrier_speed * mesh.face_area_vectors()[face_id].x;
        }
    }
    const auto options{verification_options()};
    cfd::SteadyScalarTransportSolver solver{mesh, diffusion_coefficient, study.scheme, study.limiter, options};
    cfd::CellScalarField field{mesh.cell_count()};
    LevelResult result;
    result.solve = solver.solve(flux, conditions, field);
    require(result.solve.converged, "Scalar transport outer iteration did not converge.");
    require(result.solve.field_relative_change <= options.field_relative_tolerance &&
                result.solve.discrete_relative_residual <= options.discrete_relative_tolerance,
            "Scalar transport convergence diagnostics exceed tolerances.");
    cfd::verification::ErrorAccumulator error;
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        const auto &center{mesh.cell_centers()[cell_id]};
        error.add(mesh.cell_areas()[cell_id],
                  field[cell_id] - (affine ? affine_value(center) : exponential_value(center)));
    }
    result.error = error.finish();
    check_conservation(mesh, conditions, flux, field, study);
    return result;
}

void verify_limiter_reconstruction()
{
    auto build{cfd::build_mesh(make_raw_mesh(1, 1))};
    const cfd::Mesh &mesh{build.mesh};
    const auto conditions{make_conditions(mesh, false)};
    const cfd::ScalarConvectionOperator limited{mesh, cfd::ScalarConvectionScheme::LinearUpwind,
                                                cfd::ScalarConvectionLimiter::BarthJespersen};
    const cfd::ScalarConvectionOperator unlimited{mesh, cfd::ScalarConvectionScheme::LinearUpwind};
    const cfd::CellScalarField field{mesh.cell_count(), 0.5};
    cfd::CellVectorField gradient{mesh.cell_count(), {1.0, 0.0}};
    cfd::FaceFluxField flux{mesh.face_count()};
    // Only the right outflow is nonzero: its balance exposes the reconstructed
    // face value directly. This intentionally inconsistent flux is a limiter
    // unit experiment, not a steady transport problem.
    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (mesh.face_area_vectors()[face_id].x > 0.0)
        {
            flux[face_id] = 1.0;
        }
    }
    std::vector<double> workspace(mesh.cell_count());
    cfd::CellScalarField limited_balance{mesh.cell_count()};
    cfd::CellScalarField unlimited_balance{mesh.cell_count()};
    limited.compute_flux_balance(field, conditions, flux, gradient, workspace, limited_balance);
    require(std::abs(workspace[0] - 1.0) < 1.0e-12 && std::abs(limited_balance[0] - 1.0) < 1.0e-12,
            "Barth-Jespersen modified an admissible smooth affine reconstruction.");
    gradient[0] = {4.0, 0.0};
    unlimited.compute_flux_balance(field, conditions, flux, gradient, unlimited_balance);
    limited.compute_flux_balance(field, conditions, flux, gradient, workspace, limited_balance);
    require(unlimited_balance[0] > 1.0 && std::abs(limited_balance[0] - 1.0) < 1.0e-12,
            "Barth-Jespersen failed to bound the synthetic overshoot.");
    std::cout << "Limiter: affine coefficient=1; synthetic unlimited=" << unlimited_balance[0]
              << ", limited=" << limited_balance[0] << ", coefficient=" << workspace[0] << '\n';
}

void print_order(const std::optional<double> order)
{
    if (order)
    {
        std::cout << *order;
    }
    else
    {
        std::cout << "n/a";
    }
}

} // namespace

int main()
{
    try
    {
        std::cout << std::scientific << std::setprecision(6);
        std::cout << "Pure diffusion: affine solution, Gamma=" << diffusion_coefficient << '\n';
        for (const auto &[name, shear, triangles] :
             {std::tuple{"Cartesian QUAD", 0.0, false}, std::tuple{"Sheared QUAD", 0.25, false},
              std::tuple{"Sheared TRI", 0.25, true}})
        {
            const auto result{solve_level(4, 2, studies.front(), true, shear, triangles)};
            require(result.error.linf_error < 1.0e-9, "Affine diffusion is not exact within iterative accuracy.");
            std::cout << name << " cells=" << result.error.cell_count << " RMS=" << result.error.area_weighted_rms_error
                      << " Linf=" << result.error.linf_error << " outer=" << result.solve.iteration_count
                      << " r_phi=" << result.solve.field_relative_change
                      << " r_discrete=" << result.solve.discrete_relative_residual << '\n';
        }
        std::cout << "Convection-diffusion: exponential, Pe_global="
                  << carrier_speed * domain_length / diffusion_coefficient << '\n';
        for (const auto &study : studies)
        {
            std::optional<LevelResult> previous;
            for (const auto nx : grid_levels)
            {
                const auto result{solve_level(nx, nx / 4, study, false)};
                std::cout << study.name << " nx=" << nx << " cells=" << result.error.cell_count
                          << " RMS=" << result.error.area_weighted_rms_error << " Linf=" << result.error.linf_error
                          << " order_RMS=";
                if (previous)
                {
                    require(result.error.area_weighted_rms_error < previous->error.area_weighted_rms_error &&
                                result.error.linf_error < previous->error.linf_error,
                            "Exponential solution errors did not decrease under refinement.");
                    print_order(cfd::verification::observed_order(previous->error.area_weighted_rms_error,
                                                                  result.error.area_weighted_rms_error, 2.0, 1.0));
                    std::cout << " order_Linf=";
                    print_order(cfd::verification::observed_order(previous->error.linf_error, result.error.linf_error,
                                                                  2.0, 1.0));
                }
                else
                {
                    std::cout << "- order_Linf=-";
                }
                std::cout << " outer=" << result.solve.iteration_count
                          << " r_phi=" << result.solve.field_relative_change
                          << " r_discrete=" << result.solve.discrete_relative_residual << '\n';
                previous = result;
            }
        }
        verify_limiter_reconstruction();
        std::cout << "All scalar transport verification checks passed.\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "Scalar transport verification failed: " << error.what() << '\n';
        return 1;
    }
}
