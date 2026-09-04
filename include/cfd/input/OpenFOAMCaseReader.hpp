#pragma once

#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/meshing/GmshMesher.hpp"
#include "cfd/meshing/RectangleGeometry.hpp"

#include <filesystem>
#include <string>
#include <vector>

namespace cfd
{

class Mesh;

} // namespace cfd

namespace cfd::input
{

/// Rectangle and mesh-generation settings read from `system/meshDict`.
struct MeshInput
{
    RectangleGeometry geometry;
    MeshGenerationOptions generation_options;
};

/// Scalar boundary condition associated with a user-facing boundary name.
struct NamedScalarBoundaryCondition
{
    std::string boundary_name;
    ScalarBoundaryCondition condition;
};

/// Cell-scalar field data read from an OpenFOAM-inspired field file.
struct ScalarFieldInput
{
    std::string object_name;
    double internal_value{};
    std::vector<NamedScalarBoundaryCondition> boundary_conditions;
};

/// Reads the supported CFD_solver `meshDict` subset.
///
/// @throws std::runtime_error If the file cannot be read or does not conform to
///         the supported syntax and value constraints.
[[nodiscard]]
MeshInput read_mesh_dict(const std::filesystem::path &file_path);

/// Reads the supported OpenFOAM-inspired scalar-field subset.
///
/// @throws std::runtime_error If the file cannot be read or does not conform to
///         the supported scalar-field syntax.
[[nodiscard]]
ScalarFieldInput read_scalar_field(const std::filesystem::path &file_path);

/// Resolves named conditions into the Mesh boundary-ID ordering.
///
/// The returned collection owns its conditions. Name lookup is temporary and
/// does not remain in the numerical data path.
///
/// @throws std::invalid_argument If a boundary name is unknown or duplicated,
///         or if a Mesh boundary has no condition.
[[nodiscard]]
ScalarBoundaryConditions resolve_boundary_conditions(const Mesh &mesh, const ScalarFieldInput &field_input);

} // namespace cfd::input
