#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVelocityField.hpp"
#include "cfd/field/FaceFluxField.hpp"
#include "cfd/field/PressureCorrectionBoundaryConditions.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/input/OpenFOAMCaseReader.hpp"
#include "cfd/io/MeshReport.hpp"
#include "cfd/io/VtkWriter.hpp"
#include "cfd/math/Vector2.hpp"
#include "cfd/mesh/Boundary.hpp"
#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/mesh/MeshStatistics.hpp"
#include "cfd/meshing/GmshMesher.hpp"
#include "cfd/meshing/RawMeshData.hpp"
#include "cfd/meshing/RectangleGeometry.hpp"
#include "cfd/numerics/IncompressibleSimpleSolver.hpp"
#include "cfd/numerics/MassFluxBalance.hpp"
#include "cfd/numerics/ScalarConvectionOperator.hpp"

#include <array>
#include <chrono>
#include <cmath>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace
{

constexpr double density{1.0};
constexpr double dynamic_viscosity{0.1};
constexpr std::array<double, 7> velocity_dimensions{0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0};
constexpr std::array<double, 7> physical_pressure_dimensions{1.0, -1.0, -2.0, 0.0, 0.0, 0.0, 0.0};

[[nodiscard]]
std::string_view cell_type_name(const cfd::CellType cell_type)
{
    switch (cell_type)
    {
    case cfd::CellType::Triangle:
        return "triangles";

    case cfd::CellType::Quadrilateral:
        return "quadrilaterals";

    default:
        return "unknown";
    }
}

void require_field_metadata(const cfd::input::ScalarFieldInput &field_input,
                            const std::string_view expected_object_name,
                            const std::array<double, 7> &expected_dimensions,
                            const std::string_view expected_dimensions_text)
{
    if (field_input.object_name != expected_object_name)
    {
        throw std::invalid_argument("Expected field object '" + std::string{expected_object_name} + "', got '" +
                                    field_input.object_name + "'.");
    }

    if (field_input.dimensions != expected_dimensions)
    {
        throw std::invalid_argument("Field '" + field_input.object_name + "' must have dimensions " +
                                    std::string{expected_dimensions_text} + '.');
    }
}

[[nodiscard]]
cfd::PressureCorrectionBoundaryConditions make_pressure_correction_boundary_conditions(
    const cfd::ScalarBoundaryConditions &pressure_boundary_conditions)
{
    std::vector<cfd::PressureCorrectionBoundaryConditionType> conditions;
    conditions.reserve(pressure_boundary_conditions.size());

    for (cfd::BoundaryId boundary_id = 0; boundary_id < pressure_boundary_conditions.size(); ++boundary_id)
    {
        switch (pressure_boundary_conditions[boundary_id].type)
        {
        case cfd::ScalarBoundaryConditionType::Dirichlet:
            conditions.push_back(cfd::PressureCorrectionBoundaryConditionType::FixedPressure);
            break;

        case cfd::ScalarBoundaryConditionType::Neumann:
            conditions.push_back(cfd::PressureCorrectionBoundaryConditionType::FixedMassFlux);
            break;
        }
    }

    return {pressure_boundary_conditions.size(), std::move(conditions)};
}

void initialize_fixed_mass_flux_boundaries(
    const cfd::Mesh &mesh, const double density, const cfd::ScalarBoundaryConditions &u_boundary_conditions,
    const cfd::ScalarBoundaryConditions &v_boundary_conditions,
    const cfd::PressureCorrectionBoundaryConditions &pressure_correction_boundary_conditions,
    cfd::FaceFluxField &mass_flux)
{
    const auto face_adjacencies{mesh.face_adjacencies()};
    const auto face_boundary_ids{mesh.face_boundary_ids()};
    const auto face_area_vectors{mesh.face_area_vectors()};
    const auto boundary_groups{mesh.boundary_groups()};

    for (cfd::Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (!face_adjacencies[face_id].is_boundary())
        {
            continue;
        }

        const cfd::BoundaryId boundary_id{face_boundary_ids[face_id]};
        switch (pressure_correction_boundary_conditions[boundary_id])
        {
        case cfd::PressureCorrectionBoundaryConditionType::FixedMassFlux: {
            const cfd::ScalarBoundaryCondition &u_condition{u_boundary_conditions[boundary_id]};
            const cfd::ScalarBoundaryCondition &v_condition{v_boundary_conditions[boundary_id]};
            if (u_condition.type != cfd::ScalarBoundaryConditionType::Dirichlet ||
                v_condition.type != cfd::ScalarBoundaryConditionType::Dirichlet)
            {
                throw std::invalid_argument("FixedMassFlux boundary '" + boundary_groups[boundary_id].name +
                                            "' requires fixedValue conditions for both u and v.");
            }

            const cfd::Vector2 &area_vector{face_area_vectors[face_id]};
            const double boundary_mass_flux{density *
                                            (u_condition.value * area_vector.x + v_condition.value * area_vector.y)};
            if (!std::isfinite(boundary_mass_flux))
            {
                throw std::runtime_error("Computed mass flux is non-finite on boundary '" +
                                         boundary_groups[boundary_id].name + "'.");
            }
            mass_flux[face_id] = boundary_mass_flux;
            break;
        }

        case cfd::PressureCorrectionBoundaryConditionType::FixedPressure:
            break;
        }
    }
}

} // namespace

int main(const int argc, char *argv[])
{
    const std::span<char *> arguments{argv, static_cast<std::size_t>(argc)};

    if (arguments.size() != 2)
    {
        std::cerr << "Usage: CFD_solver <case-directory>\n";
        return 1;
    }

    try
    {
        const std::filesystem::path case_directory{arguments[1]};
        const cfd::input::MeshInput mesh_input{cfd::input::read_mesh_dict(case_directory / "system" / "meshDict")};
        const cfd::input::ScalarFieldInput u_input{cfd::input::read_scalar_field(case_directory / "0" / "u")};
        const cfd::input::ScalarFieldInput v_input{cfd::input::read_scalar_field(case_directory / "0" / "v")};
        const cfd::input::ScalarFieldInput pressure_input{cfd::input::read_scalar_field(case_directory / "0" / "p")};

        require_field_metadata(u_input, "u", velocity_dimensions, "[0 1 -1 0 0 0 0]");

        require_field_metadata(v_input, "v", velocity_dimensions, "[0 1 -1 0 0 0 0]");

        require_field_metadata(pressure_input, "p", physical_pressure_dimensions, "[1 -1 -2 0 0 0 0]");

        const cfd::RectangleGeometry &geometry{mesh_input.geometry};
        const cfd::MeshGenerationOptions &options{mesh_input.generation_options};

        std::cout << "============================================================\n"
                  << " CFD Solver\n"
                  << "============================================================\n";

        const auto mesh_generation_start{std::chrono::steady_clock::now()};

        cfd::RawMeshData raw_mesh{cfd::generate_mesh(geometry, options)};

        const auto mesh_generation_end{std::chrono::steady_clock::now()};

        const auto mesh_generation_duration{
            std::chrono::duration<double, std::milli>(mesh_generation_end - mesh_generation_start)};

        std::cout << "\n[Mesh generation]\n";
        std::cout << std::fixed << std::setprecision(3) << "  Domain            : rectangle " << geometry.length
                  << " x " << geometry.height << " m\n"
                  << "  Target mesh size  : " << options.mesh_size << " m\n";

        std::cout << "  Cell type         : " << cell_type_name(options.cell_type) << '\n'
                  << "  Nodes             : " << raw_mesh.nodes.size() << '\n'
                  << "  Cells             : " << raw_mesh.cell_types.size() << '\n';

        std::cout << std::fixed << std::setprecision(2) << "  Time              : " << mesh_generation_duration.count()
                  << " ms\n";

        cfd::MeshBuildResult mesh_build_result{cfd::build_mesh(std::move(raw_mesh))};
        const cfd::Mesh &mesh{mesh_build_result.mesh};
        const cfd::ScalarBoundaryConditions u_boundary_conditions{
            cfd::input::resolve_boundary_conditions(mesh, u_input)};
        const cfd::ScalarBoundaryConditions v_boundary_conditions{
            cfd::input::resolve_boundary_conditions(mesh, v_input)};
        const cfd::ScalarBoundaryConditions pressure_boundary_conditions{
            cfd::input::resolve_boundary_conditions(mesh, pressure_input)};
        const cfd::PressureCorrectionBoundaryConditions pressure_correction_boundary_conditions{
            make_pressure_correction_boundary_conditions(pressure_boundary_conditions)};

        const cfd::MeshStatistics mesh_statistics{cfd::compute_mesh_statistics(mesh)};
        cfd::write_mesh_report(std::cout, mesh, mesh_statistics, mesh_build_result.timings);

        cfd::CellVelocityField velocity{mesh.cell_count(),
                                        cfd::Vector2{u_input.internal_value, v_input.internal_value}};
        cfd::CellScalarField pressure{mesh.cell_count(), pressure_input.internal_value};
        cfd::FaceFluxField mass_flux{mesh.face_count()};
        initialize_fixed_mass_flux_boundaries(mesh, density, u_boundary_conditions, v_boundary_conditions,
                                              pressure_correction_boundary_conditions, mass_flux);

        const cfd::IncompressibleSimpleOptions simple_options{
            .maximum_iterations = 2000,
            .momentum_relaxation_factor = 1.0,
            .pressure_relaxation_factor = 0.1,
            .rhie_chow_flux_relaxation_factor = 0.3,
            .velocity_relative_tolerance = 1.0e-10,
            .continuity_relative_tolerance = 1.0e-10,
            .momentum_linear_solver = {.relative_tolerance = 1.0e-12, .maximum_iterations = 5000},
            .pressure_correction_linear_solver = {.relative_tolerance = 1.0e-12, .maximum_iterations = 5000},
        };
        cfd::IncompressibleSimpleSolver solver{mesh, density, dynamic_viscosity, cfd::ScalarConvectionScheme::Linear,
                                               simple_options};
        const cfd::IncompressibleSimpleResult simple_result{
            solver.solve(u_boundary_conditions, v_boundary_conditions, pressure_boundary_conditions,
                         pressure_correction_boundary_conditions, velocity, pressure, mass_flux)};
        if (!simple_result.converged)
        {
            throw std::runtime_error("SIMPLE did not converge within " + std::to_string(simple_result.iteration_count) +
                                     " iterations.");
        }

        cfd::CellScalarField mass_imbalance{mesh.cell_count()};
        cfd::compute_cell_mass_imbalance(mesh, mass_flux, mass_imbalance);

        const std::filesystem::path output_directory{case_directory / "results"};
        const std::filesystem::path solution_output_file{output_directory / "solution.vtu"};

        std::filesystem::create_directories(output_directory);

        const std::array<cfd::VtkCellScalarData, 2> scalar_data{
            cfd::VtkCellScalarData{"p", pressure.values()},
            cfd::VtkCellScalarData{"mass_imbalance", mass_imbalance.values()},
        };
        const std::array<cfd::VtkCellVectorComponentData, 1> vector_data{
            cfd::VtkCellVectorComponentData{"U", velocity.u().values(), velocity.v().values()},
        };
        cfd::write_vtu(mesh, solution_output_file,
                       cfd::VtkCellData{.scalars = scalar_data, .component_vectors = vector_data});

        std::cout << "\n[Summary]\n"
                  << "  Mesh              : " << mesh.node_count() << " nodes | " << mesh.cell_count() << " cells | "
                  << mesh.face_count() << " faces\n";

        std::cout << std::scientific << std::setprecision(6) << "  SIMPLE converged  : yes\n"
                  << "  Iterations        : " << simple_result.iteration_count << '\n'
                  << "  Velocity change   : " << simple_result.velocity_relative_change << '\n'
                  << "  Continuity rel.   : " << simple_result.continuity_relative_residual << '\n'
                  << "  Max mass imbalance: " << simple_result.maximum_mass_imbalance << '\n'
                  << "  Max pressure corr.: " << simple_result.maximum_pressure_correction << '\n';

        std::cout << "\n[Output]\n"
                  << "  Solution          : " << solution_output_file.string() << '\n';

        std::cout << "\n============================================================\n"
                  << " Steady SIMPLE solution complete\n"
                  << "============================================================\n";
    }
    catch (const std::exception &error)
    {
        std::cerr << "\n[Error]\n"
                  << "  " << error.what() << '\n';

        return 1;
    }

    return 0;
}
