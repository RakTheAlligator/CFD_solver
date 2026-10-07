#pragma once

#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/meshing/GeometryInput.hpp"
#include "cfd/meshing/GmshMesher.hpp"

#include <array>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace cfd
{

class Mesh;

} // namespace cfd

namespace cfd::input
{

/// User-case initialization policy, independent of the single-Mesh SIMPLE solver.
enum class InitializationType : std::uint8_t
{
    Zero,
    CoarseMesh
};

/// Optional `initialization` settings in `system/controlDict`.
///
/// An absent block selects coarseMesh. The optional positive cell-count target
/// applies to the coarsest grid-sequencing level, is not an exact meshing
/// constraint and is incompatible with `type zero`.
struct InitializationInput
{
    InitializationType type{InitializationType::CoarseMesh};
    std::optional<Index> target_coarse_cell_count{};
};

/// Application control settings read from `system/controlDict`.
struct ControlInput
{
    bool live_convergence{true};
    InitializationInput initialization{};
};

/// Geometry and mesh-generation settings read from `system/meshDict`.
struct MeshInput
{
    GeometryInput geometry;
    MeshGenerationOptions generation_options;
    AutomaticMeshingOptions automatic_meshing;
    BackwardFacingStepMeshingOptions backward_facing_step_meshing;
};

/// Scalar boundary condition associated with a user-facing boundary name.
struct NamedScalarBoundaryCondition
{
    std::string boundary_name;
    ScalarBoundaryCondition condition;
};

/// Cell-scalar field data read from an OpenFOAM-inspired field file.
///
/// Dimension exponents are retained exactly as parsed; the reader performs no
/// unit conversion.
struct ScalarFieldInput
{
    std::string object_name;
    std::array<double, 7> dimensions{};
    double internal_value{};
    std::vector<NamedScalarBoundaryCondition> boundary_conditions;
};

/// Reads the supported CFD_solver `controlDict` subset.
///
/// A missing file returns the default control settings. Live convergence is
/// enabled when its monitoring entry is absent. Initialization defaults to
/// `coarseMesh`; `initialization { type zero; }` preserves the historical
/// internalField/target-flux initialization without a coarse solve.
///
/// @throws std::runtime_error If an existing file cannot be read or does not
///         conform to the supported syntax.
[[nodiscard]]
ControlInput read_control_dict(const std::filesystem::path &file_path);

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
