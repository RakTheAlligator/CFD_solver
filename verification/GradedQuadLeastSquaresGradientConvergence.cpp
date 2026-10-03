#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVectorField.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/meshing/RawMeshData.hpp"
#include "cfd/numerics/LeastSquaresGradient.hpp"

#include "support/GradientVerification.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace
{

using cfd::verification::domain_height;
using cfd::verification::domain_length;
using cfd::verification::ErrorAccumulator;
using cfd::verification::ErrorStatistics;

constexpr std::array<cfd::Index, 4> levels{8, 16, 32, 64};
constexpr std::array width_pattern{2.0 / 3.0, 2.0 / 3.0, 4.0 / 3.0, 4.0 / 3.0};
constexpr double epsilon{std::numeric_limits<double>::epsilon()};
// Margin for coordinate differences, four normalized observations and the
// two-dimensional solve. Error tolerances also scale with phi and inverse size.
constexpr double roundoff_margin{64.0};
constexpr double length_tolerance{roundoff_margin * epsilon * domain_length};

void require(const bool condition, const std::string &message)
{
    if (!condition)
    {
        throw std::runtime_error(message);
    }
}

void require_close(const double actual, const double expected, const double tolerance, const std::string &message)
{
    if (!std::isfinite(actual) || !std::isfinite(expected) || std::abs(actual - expected) > tolerance)
    {
        throw std::runtime_error(message + ": actual=" + std::to_string(actual) +
                                 ", expected=" + std::to_string(expected));
    }
}

struct Range
{
    double minimum{std::numeric_limits<double>::infinity()};
    double maximum{};

    void add(const double value)
    {
        require(std::isfinite(value) && value >= 0.0, "Invalid geometry diagnostic");
        minimum = std::min(minimum, value);
        maximum = std::max(maximum, value);
    }
};

[[nodiscard]]
std::vector<double> make_coordinates(const cfd::Index count, const double h, const bool periodic)
{
    std::vector<double> coordinates;
    coordinates.reserve(count + 1);
    coordinates.push_back(0.0);
    for (cfd::Index index = 0; index < count; ++index)
    {
        const double width{h * (periodic ? width_pattern.at(index % width_pattern.size()) : 1.0)};
        coordinates.push_back(coordinates.back() + width);
    }
    return coordinates;
}

[[nodiscard]]
cfd::RawMeshData make_raw_mesh(const cfd::Index n, const bool periodic)
{
    require(n >= 8 && n % width_pattern.size() == 0, "Grid must contain complete width patterns");
    const cfd::Index nx{2 * n};
    const cfd::Index ny{n};
    const double h{1.0 / static_cast<double>(n)};
    const auto x{make_coordinates(nx, h, periodic)};
    const auto y{make_coordinates(ny, h, periodic)};
    // Cumulative addition can accumulate one rounding contribution per width;
    // do not snap the endpoint or rescale the grading after construction.
    require_close(x.back(), domain_length, length_tolerance * static_cast<double>(nx), "Incorrect x extent");
    require_close(y.back(), domain_height, length_tolerance * static_cast<double>(ny), "Incorrect y extent");

    cfd::RawMeshData raw;
    raw.nodes.reserve((nx + 1) * (ny + 1));
    raw.cell_types.reserve(nx * ny);
    raw.cell_nodes.reserve(4 * nx * ny);
    raw.cell_node_offsets.reserve(nx * ny + 1);
    raw.cell_node_offsets.push_back(0);
    const auto node_id = [nx](const cfd::Index i, const cfd::Index j) { return j * (nx + 1) + i; };
    for (const double yi : y)
    {
        for (const double xi : x)
        {
            raw.nodes.push_back({xi, yi});
        }
    }
    for (cfd::Index j = 0; j < ny; ++j)
    {
        for (cfd::Index i = 0; i < nx; ++i)
        {
            raw.cell_types.push_back(cfd::CellType::Quadrilateral);
            raw.cell_nodes.insert(raw.cell_nodes.end(),
                                  {node_id(i, j), node_id(i + 1, j), node_id(i + 1, j + 1), node_id(i, j + 1)});
            raw.cell_node_offsets.push_back(raw.cell_nodes.size());
        }
    }
    using namespace cfd::verification;
    raw.boundary_groups = {{left_boundary_id, "left"},
                           {right_boundary_id, "right"},
                           {bottom_boundary_id, "bottom"},
                           {top_boundary_id, "top"}};
    raw.boundary_edges.reserve(2 * (nx + ny));
    for (cfd::Index j = 0; j < ny; ++j)
    {
        raw.boundary_edges.push_back({{node_id(0, j), node_id(0, j + 1)}, left_boundary_id});
        raw.boundary_edges.push_back({{node_id(nx, j), node_id(nx, j + 1)}, right_boundary_id});
    }
    for (cfd::Index i = 0; i < nx; ++i)
    {
        raw.boundary_edges.push_back({{node_id(i, 0), node_id(i + 1, 0)}, bottom_boundary_id});
        raw.boundary_edges.push_back({{node_id(i, ny), node_id(i + 1, ny)}, top_boundary_id});
    }
    return raw;
}

[[nodiscard]]
Range measure_widths(const cfd::Mesh &mesh, const double h, const bool periodic)
{
    Range widths;
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        require(mesh.cell_types()[cell_id] == cfd::CellType::Quadrilateral, "Expected only QUAD cells");
        const cfd::Index begin{mesh.cell_node_offsets()[cell_id]};
        const cfd::Index end{mesh.cell_node_offsets()[cell_id + 1]};
        require(end - begin == 4, "Expected four nodes per QUAD");
        double xmin{std::numeric_limits<double>::infinity()};
        double xmax{-std::numeric_limits<double>::infinity()};
        double ymin{std::numeric_limits<double>::infinity()};
        double ymax{-std::numeric_limits<double>::infinity()};
        for (cfd::Index position = begin; position < end; ++position)
        {
            const auto &node{mesh.nodes()[mesh.cell_nodes()[position]]};
            xmin = std::min(xmin, node.x);
            xmax = std::max(xmax, node.x);
            ymin = std::min(ymin, node.y);
            ymax = std::max(ymax, node.y);
        }
        require(xmax > xmin && ymax > ymin, "Non-positive QUAD dimensions");
        widths.add(xmax - xmin);
        widths.add(ymax - ymin);
    }
    require_close(widths.minimum, periodic ? 2.0 * h / 3.0 : h, length_tolerance, "Incorrect minimum width");
    require_close(widths.maximum, periodic ? 4.0 * h / 3.0 : h, length_tolerance, "Incorrect maximum width");
    require_close(widths.maximum / widths.minimum, periodic ? 2.0 : 1.0, 4.0 * length_tolerance / widths.minimum,
                  "Incorrect width ratio");
    return widths;
}

struct CellStencil
{
    bool boundary{};
    double a_x{};
    double b_x{};
    double a_y{};
    double b_y{};
};

struct GeometryDiagnostics
{
    Range widths;
    Range relative_asymmetry;
    double maximum_angular_defect{};
    double maximum_distance_difference_defect{};
    cfd::Index interior_count{};
    cfd::Index boundary_count{};
};

void check_pair(const double a, const double b, const double h, const bool periodic, GeometryDiagnostics &diagnostics)
{
    require(a > 0.0 && b > 0.0, "Incomplete opposite-neighbor pair");
    const double difference{std::abs(a - b)};
    const double expected_difference{periodic ? h / 3.0 : 0.0};
    require_close(difference, expected_difference, length_tolerance, "Incorrect opposite distance difference");
    diagnostics.maximum_distance_difference_defect =
        std::max(diagnostics.maximum_distance_difference_defect, std::abs(difference - expected_difference));
    const double asymmetry{difference / (a + b)};
    const double tolerance{4.0 * length_tolerance / std::min(a, b)};
    if (periodic)
    {
        require(std::min(std::abs(asymmetry - 1.0 / 5.0), std::abs(asymmetry - 1.0 / 7.0)) <= tolerance,
                "Relative asymmetry is not a predicted pattern phase");
    }
    else
    {
        require_close(asymmetry, 0.0, tolerance, "Uniform neighbors are not equidistant");
    }
    diagnostics.relative_asymmetry.add(asymmetry);
}

[[nodiscard]]
CellStencil inspect_stencil(const cfd::Mesh &mesh, const cfd::Index cell_id, const double h, const bool periodic,
                            GeometryDiagnostics &diagnostics)
{
    CellStencil stencil;
    std::array<cfd::Index, 4> neighbors{};
    std::array<cfd::Vector2, 4> displacements{};
    cfd::Index count{};
    const auto &center{mesh.cell_centers()[cell_id]};
    const cfd::Index begin{mesh.cell_node_offsets()[cell_id]};
    const cfd::Index end{mesh.cell_node_offsets()[cell_id + 1]};
    require(end - begin == 4, "Expected exactly four cell faces");
    for (cfd::Index position = begin; position < end; ++position)
    {
        const auto &adjacency{mesh.face_adjacencies()[mesh.cell_faces()[position]]};
        if (adjacency.is_boundary())
        {
            stencil.boundary = true;
            continue;
        }
        const cfd::Index neighbor{adjacency.owner == cell_id ? adjacency.neighbor : adjacency.owner};
        for (cfd::Index previous = 0; previous < count; ++previous)
        {
            require(neighbors.at(previous) != neighbor, "Cell has duplicate face neighbors");
        }
        neighbors.at(count) = neighbor;
        const auto &other{mesh.cell_centers()[neighbor]};
        displacements.at(count) = {other.x - center.x, other.y - center.y};
        ++count;
    }
    if (stencil.boundary)
    {
        ++diagnostics.boundary_count;
        return stencil;
    }
    require(count == 4, "Controlled interior cell must have four distinct internal neighbors");
    ++diagnostics.interior_count;
    // Determine directions from actual Mesh centers, never from construction
    // indices or an assumed local face ordering. Axis checks ensure two
    // independent directions, in addition to the shared opposite-pair metric.
    for (const auto &d : displacements)
    {
        const double distance{cfd::verification::magnitude(d)};
        double *slot{};
        if (std::abs(d.y) <= length_tolerance && std::abs(d.x) > length_tolerance)
        {
            slot = d.x > 0.0 ? &stencil.a_x : &stencil.b_x;
        }
        else if (std::abs(d.x) <= length_tolerance && std::abs(d.y) > length_tolerance)
        {
            slot = d.y > 0.0 ? &stencil.a_y : &stencil.b_y;
        }
        else
        {
            throw std::runtime_error("Neighbor displacement is not aligned with a Cartesian axis");
        }
        require(*slot == 0.0, "Repeated neighbor direction");
        *slot = distance;
    }
    const auto pairs{cfd::verification::compute_opposite_pair_metrics(displacements)};
    require(pairs.valid && pairs.angular_defect <= roundoff_margin * epsilon, "Neighbors are not opposite pairs");
    diagnostics.maximum_angular_defect = std::max(diagnostics.maximum_angular_defect, pairs.angular_defect);
    check_pair(stencil.a_x, stencil.b_x, h, periodic, diagnostics);
    check_pair(stencil.a_y, stencil.b_y, h, periodic, diagnostics);
    return stencil;
}

struct AnalyticalField
{
    std::string_view name;
    double constant;
    cfd::Vector2 linear;
    double quadratic;
};

constexpr std::array analytical_fields{AnalyticalField{"constant", 3.25, {}, 0.0},
                                       AnalyticalField{"affine", 1.5, {2.0, -3.0}, 0.0},
                                       AnalyticalField{"quadratic", 0.0, {}, 1.0}};

[[nodiscard]]
cfd::ScalarBoundaryConditions make_conditions(const AnalyticalField &field)
{
    // Only gradient reconstruction is performed; phi is analytically supplied.
    // These exact Neumann observations therefore introduce no PDE nullspace.
    return {4,
            {{cfd::ScalarBoundaryConditionType::Neumann, -field.linear.x},
             {cfd::ScalarBoundaryConditionType::Neumann, field.linear.x + 2.0 * field.quadratic * domain_length},
             {cfd::ScalarBoundaryConditionType::Neumann, -field.linear.y},
             {cfd::ScalarBoundaryConditionType::Neumann, field.linear.y + 2.0 * field.quadratic * domain_height}}};
}

struct FieldErrors
{
    ErrorStatistics interior;
    ErrorStatistics boundary;
    ErrorStatistics global;
    double roundoff_tolerance{};
    double maximum_prediction_mismatch{};
};

[[nodiscard]]
FieldErrors verify_field(const cfd::Mesh &mesh, const std::vector<CellStencil> &stencils, const AnalyticalField &field,
                         const double minimum_width)
{
    cfd::CellScalarField phi{mesh.cell_count()};
    cfd::CellVectorField gradient{mesh.cell_count()};
    double scale_phi{1.0};
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        const auto &p{mesh.cell_centers()[cell_id]};
        phi[cell_id] =
            field.constant + field.linear.x * p.x + field.linear.y * p.y + field.quadratic * (p.x * p.x + p.y * p.y);
        scale_phi = std::max(scale_phi, std::abs(phi[cell_id]));
    }
    const auto conditions{make_conditions(field)};
    cfd::compute_least_squares_gradient(mesh, phi, conditions, gradient);
    FieldErrors result;
    // Subtraction of phi values followed by division by a neighbor distance
    // sets the absolute roundoff scale of the reconstructed gradient.
    result.roundoff_tolerance = roundoff_margin * epsilon * scale_phi / minimum_width;
    ErrorAccumulator interior;
    ErrorAccumulator boundary;
    ErrorAccumulator global;
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        const auto &p{mesh.cell_centers()[cell_id]};
        const cfd::Vector2 error{gradient[cell_id].x - (field.linear.x + 2.0 * field.quadratic * p.x),
                                 gradient[cell_id].y - (field.linear.y + 2.0 * field.quadratic * p.y)};
        const double magnitude{cfd::verification::magnitude(error)};
        const double area{mesh.cell_areas()[cell_id]};
        require(std::isfinite(magnitude) && std::isfinite(area) && area > 0.0, "Invalid gradient error data");
        global.add(area, magnitude);
        const auto &stencil{stencils[cell_id]};
        if (stencil.boundary)
        {
            boundary.add(area, magnitude);
            continue;
        }
        interior.add(area, magnitude);
        // For the current inverse-distance-normalized WLS on two orthogonal
        // pairs, Taylor is exact for this quadratic: e = (a-b)/2 per axis.
        // This is an error diagnostic, not an alternative reconstruction.
        const cfd::Vector2 predicted{0.5 * field.quadratic * (stencil.a_x - stencil.b_x),
                                     0.5 * field.quadratic * (stencil.a_y - stencil.b_y)};
        require_close(error.x, predicted.x, result.roundoff_tolerance, "Cellwise x error disagrees with prediction");
        require_close(error.y, predicted.y, result.roundoff_tolerance, "Cellwise y error disagrees with prediction");
        result.maximum_prediction_mismatch =
            std::max(result.maximum_prediction_mismatch, std::hypot(error.x - predicted.x, error.y - predicted.y));
    }
    result.interior = interior.finish();
    result.boundary = boundary.finish();
    result.global = global.finish();
    require(result.interior.cell_count > 0 && result.boundary.cell_count > 0, "Empty verification population");
    return result;
}

struct LevelResult
{
    cfd::Index n{};
    double h{};
    GeometryDiagnostics geometry;
    FieldErrors quadratic;
    double expected_error{};
};

[[nodiscard]]
LevelResult run_level(const cfd::Index n, const bool periodic)
{
    auto build{cfd::build_mesh(make_raw_mesh(n, periodic))};
    const cfd::Mesh &mesh{build.mesh};
    require(mesh.cell_count() == 2 * n * n && mesh.boundary_groups().size() == 4, "Incorrect Mesh cardinalities");
    LevelResult result;
    result.n = n;
    result.h = 1.0 / static_cast<double>(n);
    result.geometry.widths = measure_widths(mesh, result.h, periodic);
    std::vector<CellStencil> stencils;
    stencils.reserve(mesh.cell_count());
    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        stencils.push_back(inspect_stencil(mesh, cell_id, result.h, periodic, result.geometry));
    }
    require(result.geometry.interior_count == (2 * n - 2) * (n - 2) &&
                result.geometry.boundary_count + result.geometry.interior_count == mesh.cell_count(),
            "Incorrect interior/boundary classification");
    if (periodic)
    {
        const double tolerance{4.0 * length_tolerance / result.geometry.widths.minimum};
        require_close(result.geometry.relative_asymmetry.minimum, 1.0 / 7.0, tolerance, "Asymmetry minimum changed");
        require_close(result.geometry.relative_asymmetry.maximum, 1.0 / 5.0, tolerance, "Asymmetry maximum changed");
    }
    for (const auto &field : analytical_fields)
    {
        const auto errors{verify_field(mesh, stencils, field, result.geometry.widths.minimum)};
        if (field.quadratic == 0.0)
        {
            require(errors.interior.linf_error <= 2.0 * errors.roundoff_tolerance,
                    "Constant/affine interior gradient is not exact within roundoff");
            std::cout << "  " << field.name << " interior Linf=" << errors.interior.linf_error << " PASS\n";
        }
        else
        {
            result.quadratic = errors;
        }
    }
    result.expected_error = periodic ? std::sqrt(2.0) * result.h / 6.0 : 0.0;
    require_close(result.quadratic.interior.area_weighted_rms_error, result.expected_error,
                  2.0 * result.quadratic.roundoff_tolerance, "Interior RMS differs from expected error");
    require_close(result.quadratic.interior.linf_error, result.expected_error,
                  2.0 * result.quadratic.roundoff_tolerance, "Interior Linf differs from expected error");
    return result;
}

void print_level(const LevelResult &result, const bool periodic, const std::optional<LevelResult> &previous)
{
    const auto &geometry{result.geometry};
    const auto &error{result.quadratic};
    std::cout << (periodic ? "periodic" : "uniform") << " N=" << result.n << " cells=" << error.global.cell_count
              << " interior=" << geometry.interior_count << " boundary=" << geometry.boundary_count << " h=" << result.h
              << '\n'
              << "  geometry: min_width=" << geometry.widths.minimum << " max_width=" << geometry.widths.maximum
              << " ratio=" << geometry.widths.maximum / geometry.widths.minimum
              << " max_angular_defect=" << geometry.maximum_angular_defect
              << " asymmetry_min=" << geometry.relative_asymmetry.minimum
              << " asymmetry_max=" << geometry.relative_asymmetry.maximum
              << " max_distance_difference_defect=" << geometry.maximum_distance_difference_defect << '\n'
              << "  interior: RMS=" << error.interior.area_weighted_rms_error << " Linf=" << error.interior.linf_error
              << " predicted=" << result.expected_error << " RMS_order=";
    if (!periodic)
    {
        // The uniform quadratic is exact in arithmetic; a logarithmic ratio
        // of floating-point noise is not a spatial convergence order.
        std::cout << "roundoff / n/a Linf_order=roundoff / n/a";
    }
    else if (!previous)
    {
        std::cout << "- Linf_order=-";
    }
    else
    {
        const auto rms_order{cfd::verification::observed_order(previous->quadratic.interior.area_weighted_rms_error,
                                                               error.interior.area_weighted_rms_error, previous->h,
                                                               result.h)};
        const auto linf_order{cfd::verification::observed_order(previous->quadratic.interior.linf_error,
                                                                error.interior.linf_error, previous->h, result.h)};
        if (!rms_order || !linf_order)
        {
            throw std::runtime_error("Periodic convergence order is undefined");
        }
        // Propagate the roundoff bounds on the two errors through log(Ec/Ef).
        const double order_tolerance{8.0 *
                                     (error.roundoff_tolerance / result.expected_error +
                                      previous->quadratic.roundoff_tolerance / previous->expected_error) /
                                     std::log(previous->h / result.h)};
        require_close(*rms_order, 1.0, order_tolerance, "Periodic RMS order is not one");
        require_close(*linf_order, 1.0, order_tolerance, "Periodic Linf order is not one");
        std::cout << *rms_order << " Linf_order=" << *linf_order;
    }
    std::cout << '\n'
              << "  boundary: RMS=" << error.boundary.area_weighted_rms_error << " Linf=" << error.boundary.linf_error
              << '\n'
              << "  global: RMS=" << error.global.area_weighted_rms_error << " Linf=" << error.global.linf_error << '\n'
              << "  prediction: max_cellwise_mismatch=" << error.maximum_prediction_mismatch
              << " component_roundoff_tolerance=" << error.roundoff_tolerance << '\n';
}

} // namespace

int main()
{
    try
    {
        std::cout << std::scientific << std::setprecision(9)
                  << "Controlled QUAD WLS verification: [0,2] x [0,1], phi=x^2+y^2.\n"
                  << "Periodic widths: h * {2/3,2/3,4/3,4/3} on both axes. No PDE solve.\n";
        for (const bool periodic : {false, true})
        {
            std::optional<LevelResult> previous;
            for (const auto n : levels)
            {
                std::cout << "\nChecking " << (periodic ? "periodic" : "uniform") << " N=" << n << '\n';
                const auto result{run_level(n, periodic)};
                print_level(result, periodic, previous);
                previous = result;
            }
        }
        std::cout << "\nPASS: uniform interior quadratic error is at roundoff.\n"
                  << "PASS: periodically non-uniform interior error follows O(h).\n"
                  << "PASS: measured cellwise error matches (a-b)/2 from actual Mesh centers.\n"
                  << "Boundary/global results are descriptive only; conclusions apply to this controlled family.\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "Graded QUAD WLS verification failed: " << error.what() << '\n';
        return 1;
    }
}
