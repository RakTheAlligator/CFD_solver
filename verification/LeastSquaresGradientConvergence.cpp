#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVectorField.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/io/VtkWriter.hpp"
#include "cfd/math/Point2.hpp"
#include "cfd/math/Vector2.hpp"
#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/meshing/GmshMesher.hpp"
#include "cfd/meshing/RawMeshData.hpp"
#include "cfd/meshing/RectangleGeometry.hpp"
#include "cfd/numerics/LeastSquaresGradient.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace
{

constexpr double domain_length{2.0};
constexpr double domain_height{1.0};
constexpr double boundary_coordinate_tolerance{1.0e-12};

constexpr cfd::BoundaryId left_boundary_id{0};
constexpr cfd::BoundaryId right_boundary_id{1};
constexpr cfd::BoundaryId bottom_boundary_id{2};
constexpr cfd::BoundaryId top_boundary_id{3};

constexpr std::array target_mesh_sizes{0.20, 0.10, 0.05, 0.025};

struct MeshFamily
{
    cfd::CellType cell_type;
    std::string_view table_name;
    std::string_view file_prefix;
};

constexpr std::array mesh_families{
    MeshFamily{cfd::CellType::Triangle, "triangles", "triangle"},
    MeshFamily{cfd::CellType::Quadrilateral, "quadrilaterals", "quad"},
};

struct LevelResult
{
    double target_mesh_size{};
    cfd::Index cell_count{};
    double characteristic_mesh_length{};
    double l2_error{};
    double linf_error{};
    std::optional<double> l2_order;
    std::optional<double> linf_order;
};

[[nodiscard]]
bool matches_coordinate(const double first, const double second, const double target) noexcept
{
    return std::abs(first - target) <= boundary_coordinate_tolerance &&
           std::abs(second - target) <= boundary_coordinate_tolerance;
}

void split_rectangle_boundary_groups(cfd::RawMeshData &raw_mesh)
{
    // The production rectangle mesher combines top and bottom as `wall`.
    // This analytical case needs different exact Neumann data on those sides,
    // so the verification study relabels only its generated raw boundary edges.
    raw_mesh.boundary_groups = {
        {left_boundary_id, "left"},
        {right_boundary_id, "right"},
        {bottom_boundary_id, "bottom"},
        {top_boundary_id, "top"},
    };

    std::array<cfd::Index, 4> boundary_edge_counts{};

    for (cfd::BoundaryEdge &edge : raw_mesh.boundary_edges)
    {
        const cfd::Point2 &first{raw_mesh.nodes[edge.node_ids[0]]};
        const cfd::Point2 &second{raw_mesh.nodes[edge.node_ids[1]]};

        if (matches_coordinate(first.x, second.x, 0.0))
        {
            edge.boundary_id = left_boundary_id;
        }
        else if (matches_coordinate(first.x, second.x, domain_length))
        {
            edge.boundary_id = right_boundary_id;
        }
        else if (matches_coordinate(first.y, second.y, 0.0))
        {
            edge.boundary_id = bottom_boundary_id;
        }
        else if (matches_coordinate(first.y, second.y, domain_height))
        {
            edge.boundary_id = top_boundary_id;
        }
        else
        {
            throw std::runtime_error("Unable to classify a generated rectangle boundary edge.");
        }

        ++boundary_edge_counts.at(edge.boundary_id);
    }

    if (std::any_of(boundary_edge_counts.begin(), boundary_edge_counts.end(),
                    [](const cfd::Index count) { return count == 0; }))
    {
        throw std::runtime_error("Generated rectangle does not contain all four verification boundaries.");
    }
}

// This is code verification against a smooth analytical field, not validation
// against a physical experiment. Here phi = sin(x) + y^2 and
// grad(phi) = (cos(x), 2y).
[[nodiscard]]
double analytical_phi(const cfd::Point2 &point) noexcept
{
    return std::sin(point.x) + point.y * point.y;
}

[[nodiscard]]
cfd::Vector2 analytical_gradient(const cfd::Point2 &point) noexcept
{
    return {std::cos(point.x), 2.0 * point.y};
}

[[nodiscard]]
double outward_normal_derivative(const std::string_view boundary_name)
{
    if (boundary_name == "left")
    {
        return -1.0;
    }

    if (boundary_name == "right")
    {
        return std::cos(domain_length);
    }

    if (boundary_name == "bottom")
    {
        return 0.0;
    }

    if (boundary_name == "top")
    {
        return 2.0;
    }

    throw std::runtime_error("Unsupported boundary group in gradient verification: " + std::string{boundary_name});
}

[[nodiscard]]
cfd::ScalarBoundaryConditions make_boundary_conditions(const cfd::Mesh &mesh)
{
    std::vector<cfd::ScalarBoundaryCondition> conditions;
    conditions.reserve(mesh.boundary_groups().size());

    // Pure Neumann data are valid here because phi is supplied analytically at
    // cell centers; this study reconstructs its gradient and does not solve a PDE.
    for (const cfd::BoundaryGroup &group : mesh.boundary_groups())
    {
        if (group.id != conditions.size())
        {
            throw std::runtime_error("Mesh boundary groups are not indexed contiguously.");
        }

        conditions.emplace_back(cfd::ScalarBoundaryConditionType::Neumann, outward_normal_derivative(group.name));
    }

    return {mesh.boundary_groups().size(), std::move(conditions)};
}

void write_verification_vtu(const cfd::Mesh &mesh, const cfd::CellScalarField &phi,
                            const cfd::CellVectorField &numerical_gradient, const cfd::CellVectorField &exact_gradient,
                            const cfd::CellVectorField &gradient_error, const cfd::CellScalarField &error_magnitude,
                            const std::filesystem::path &output_path)
{
    const std::array scalar_fields{
        cfd::VtkCellScalarData{"phi", phi.values()},
        cfd::VtkCellScalarData{"grad_phi_error_magnitude", error_magnitude.values()},
    };
    const std::array vector_fields{
        cfd::VtkCellVectorData{"grad_phi", numerical_gradient.values()},
        cfd::VtkCellVectorData{"grad_phi_exact", exact_gradient.values()},
        cfd::VtkCellVectorData{"grad_phi_error", gradient_error.values()},
    };
    const cfd::VtkCellData cell_data{
        .scalars = scalar_fields,
        .vectors = vector_fields,
    };

    cfd::write_vtu(mesh, output_path, cell_data);
}

[[nodiscard]]
LevelResult run_level(const MeshFamily &family, const double target_mesh_size, const std::size_t level_index,
                      const bool write_vtu, const std::filesystem::path &output_directory,
                      std::vector<std::filesystem::path> &written_paths)
{
    const cfd::RectangleGeometry geometry{
        .length = domain_length,
        .height = domain_height,
    };
    const cfd::MeshGenerationOptions options{
        .mesh_size = target_mesh_size,
        .cell_type = family.cell_type,
    };

    cfd::RawMeshData raw_mesh{cfd::generate_mesh(geometry, options)};
    split_rectangle_boundary_groups(raw_mesh);

    cfd::MeshBuildResult build_result{cfd::build_mesh(std::move(raw_mesh))};
    const cfd::Mesh &mesh{build_result.mesh};

    if (mesh.cell_count() == 0)
    {
        throw std::runtime_error("Gradient verification requires a non-empty mesh.");
    }

    cfd::CellScalarField phi{mesh.cell_count()};
    cfd::CellVectorField exact_gradient{mesh.cell_count()};

    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        phi[cell_id] = analytical_phi(mesh.cell_centers()[cell_id]);
        exact_gradient[cell_id] = analytical_gradient(mesh.cell_centers()[cell_id]);
    }

    const cfd::ScalarBoundaryConditions boundary_conditions{make_boundary_conditions(mesh)};
    cfd::CellVectorField numerical_gradient{mesh.cell_count()};
    cfd::compute_least_squares_gradient(mesh, phi, boundary_conditions, numerical_gradient);

    cfd::CellVectorField gradient_error{mesh.cell_count()};
    cfd::CellScalarField error_magnitude{mesh.cell_count()};

    double total_cell_area{};
    double area_weighted_squared_error{};
    double linf_error{};

    for (cfd::Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        const double cell_area{mesh.cell_areas()[cell_id]};
        if (!std::isfinite(cell_area) || !(cell_area > 0.0))
        {
            throw std::runtime_error("Gradient verification encountered an invalid cell area.");
        }

        gradient_error[cell_id] = {
            numerical_gradient[cell_id].x - exact_gradient[cell_id].x,
            numerical_gradient[cell_id].y - exact_gradient[cell_id].y,
        };

        const double squared_error{gradient_error[cell_id].x * gradient_error[cell_id].x +
                                   gradient_error[cell_id].y * gradient_error[cell_id].y};
        if (!std::isfinite(squared_error) || squared_error < 0.0)
        {
            throw std::runtime_error("Gradient verification produced a non-finite error.");
        }

        error_magnitude[cell_id] = std::sqrt(squared_error);
        total_cell_area += cell_area;

        // Normalization by total_cell_area below makes this an area-weighted,
        // domain-normalized (RMS) L2 error measure.
        area_weighted_squared_error += cell_area * squared_error;
        linf_error = std::max(linf_error, error_magnitude[cell_id]);
    }

    if (!std::isfinite(total_cell_area) || !(total_cell_area > 0.0))
    {
        throw std::runtime_error("Gradient verification produced an invalid total cell area.");
    }

    if (!std::isfinite(area_weighted_squared_error) || area_weighted_squared_error < 0.0 || !std::isfinite(linf_error))
    {
        throw std::runtime_error("Gradient verification produced an invalid error norm.");
    }

    const double l2_error{std::sqrt(area_weighted_squared_error / total_cell_area)};

    // Observed order uses this measured geometric length because a requested
    // Gmsh target size is not the exact spacing of the generated cells.
    const double characteristic_mesh_length{std::sqrt(total_cell_area / static_cast<double>(mesh.cell_count()))};

    if (!std::isfinite(l2_error) || !std::isfinite(characteristic_mesh_length) || !(characteristic_mesh_length > 0.0))
    {
        throw std::runtime_error("Gradient verification produced an invalid L2 error or mesh length.");
    }

    if (write_vtu)
    {
        const std::filesystem::path output_path{
            output_directory / (std::string{family.file_prefix} + "_level_" + std::to_string(level_index) + ".vtu")};

        write_verification_vtu(mesh, phi, numerical_gradient, exact_gradient, gradient_error, error_magnitude,
                               output_path);
        written_paths.push_back(output_path);
    }

    return {
        .target_mesh_size = target_mesh_size,
        .cell_count = mesh.cell_count(),
        .characteristic_mesh_length = characteristic_mesh_length,
        .l2_error = l2_error,
        .linf_error = linf_error,
        .l2_order = std::nullopt,
        .linf_order = std::nullopt,
    };
}

[[nodiscard]]
std::optional<double> observed_order(const double coarse_error, const double fine_error, const double coarse_h,
                                     const double fine_h, const std::string_view norm_name)
{
    if (!std::isfinite(coarse_error) || !std::isfinite(fine_error) || coarse_error < 0.0 || fine_error < 0.0)
    {
        throw std::runtime_error("Cannot compute " + std::string{norm_name} + " order from invalid error values.");
    }

    if (!std::isfinite(coarse_h) || !std::isfinite(fine_h) || !(coarse_h > 0.0) || !(fine_h > 0.0))
    {
        throw std::runtime_error("Cannot compute " + std::string{norm_name} + " order from invalid mesh lengths.");
    }

    // An exact zero makes the logarithmic error ratio undefined and is
    // reported explicitly as `n/a` in the convergence table.
    if (coarse_error == 0.0 || fine_error == 0.0)
    {
        return std::nullopt;
    }

    const double logarithmic_mesh_ratio{std::log(coarse_h / fine_h)};
    if (!std::isfinite(logarithmic_mesh_ratio) || logarithmic_mesh_ratio == 0.0)
    {
        throw std::runtime_error("Cannot compute " + std::string{norm_name} +
                                 " order because the mesh-length ratio is invalid.");
    }

    const double order{std::log(coarse_error / fine_error) / logarithmic_mesh_ratio};
    if (!std::isfinite(order))
    {
        throw std::runtime_error("Computed " + std::string{norm_name} + " order is non-finite.");
    }

    return order;
}

void compute_observed_orders(std::vector<LevelResult> &results)
{
    for (std::size_t level_index = 1; level_index < results.size(); ++level_index)
    {
        const LevelResult &coarse{results[level_index - 1]};
        LevelResult &fine{results[level_index]};

        fine.l2_order = observed_order(coarse.l2_error, fine.l2_error, coarse.characteristic_mesh_length,
                                       fine.characteristic_mesh_length, "L2");
        fine.linf_order = observed_order(coarse.linf_error, fine.linf_error, coarse.characteristic_mesh_length,
                                         fine.characteristic_mesh_length, "Linf");
    }
}

void print_order(const std::optional<double> order, const bool first_level)
{
    if (order.has_value())
    {
        std::cout << std::fixed << std::setprecision(3) << std::setw(10) << *order;
        return;
    }

    std::cout << std::setw(10) << (first_level ? "-" : "n/a");
}

void print_convergence_table(const MeshFamily &family, const std::vector<LevelResult> &results)
{
    std::cout << "\nLeast-squares gradient convergence - " << family.table_name << "\n\n"
              << std::left << std::setw(12) << "target_h" << std::setw(12) << "cells" << std::setw(16) << "h_char"
              << std::setw(16) << "L2_error" << std::setw(16) << "Linf_error" << std::setw(10) << "p_L2"
              << std::setw(10) << "p_Linf" << '\n'
              << std::string(92, '-') << '\n';

    for (std::size_t level_index = 0; level_index < results.size(); ++level_index)
    {
        const LevelResult &result{results[level_index]};

        std::cout << std::fixed << std::setprecision(3) << std::setw(12) << result.target_mesh_size << std::setw(12)
                  << result.cell_count << std::scientific << std::setprecision(6) << std::setw(16)
                  << result.characteristic_mesh_length << std::setw(16) << result.l2_error << std::setw(16)
                  << result.linf_error;

        print_order(result.l2_order, level_index == 0);
        print_order(result.linf_order, level_index == 0);
        std::cout << '\n';
    }
}

[[nodiscard]]
bool parse_write_vtu_flag(const std::span<char *> arguments)
{
    if (arguments.size() == 1)
    {
        return false;
    }

    if (arguments.size() == 2 && std::string_view{arguments[1]} == "--write-vtu")
    {
        return true;
    }

    throw std::invalid_argument("Unsupported arguments.\nUsage: cfd_least_squares_gradient_convergence [--write-vtu]");
}

} // namespace

int main(const int argument_count, char **arguments)
{
    try
    {
        const std::span argument_view{arguments, static_cast<std::size_t>(argument_count)};
        const bool write_vtu{parse_write_vtu_flag(argument_view)};
        const std::filesystem::path output_directory{"output/verification/least_squares_gradient"};

        if (write_vtu)
        {
            std::filesystem::create_directories(output_directory);
        }

        std::vector<std::filesystem::path> written_paths;

        for (const MeshFamily &family : mesh_families)
        {
            std::vector<LevelResult> results;
            results.reserve(target_mesh_sizes.size());

            std::size_t level_index{};
            for (const double target_mesh_size : target_mesh_sizes)
            {
                results.push_back(
                    run_level(family, target_mesh_size, level_index, write_vtu, output_directory, written_paths));
                ++level_index;
            }

            compute_observed_orders(results);
            print_convergence_table(family, results);
        }

        if (write_vtu)
        {
            std::cout << "\nVTU output:\n";
            for (const std::filesystem::path &path : written_paths)
            {
                std::cout << "  " << path.string() << '\n';
            }
        }
    }
    catch (const std::exception &error)
    {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }

    return 0;
}
