#pragma once

#include "cfd/mesh/Types.hpp"

#include <cstdint>
#include <string>

namespace cfd
{
class Mesh;
class CellVelocityField;
class CellScalarField;
class FaceFluxField;
class ScalarBoundaryConditions;
class PressureCorrectionBoundaryConditions;
struct IncompressibleSimpleOptions;
enum class ScalarConvectionScheme : std::uint8_t;

namespace input
{
struct MeshInput;
struct InitializationInput;
struct ScalarFieldInput;
} // namespace input
} // namespace cfd

namespace cfd::app
{

struct CoarseMeshPlan
{
    Index target_cell_count{};
    double mesh_size{};
};

/// Computes the requested coarse resolution from the actual final cell count.
/// Explicit targets must be positive and strictly smaller than the final count.
[[nodiscard]]
CoarseMeshPlan make_coarse_mesh_plan(Index final_cell_count, double final_mesh_size,
                                     const input::InitializationInput &initialization);

struct CoarseInitializationResult
{
    bool used_coarse_mesh{};
    bool automatic_fallback{};
    std::string fallback_reason;
    CoarseMeshPlan plan;
    Index actual_cell_count{};
    Index iteration_count{};
    double solve_seconds{};
    /// Mapping construction plus WLS reconstruction/transfer of u, v and p.
    double transfer_seconds{};
};

/// Maps physical scalar pressure conditions to the current p' boundary semantics.
[[nodiscard]]
PressureCorrectionBoundaryConditions make_pressure_correction_boundary_conditions(
    const ScalarBoundaryConditions &pressure_boundary_conditions);

/// Applies the historical target initialization to prescribed boundary fluxes.
/// Internal and FixedPressure face entries are unchanged.
void initialize_fixed_mass_flux_boundaries(
    const Mesh &mesh, double density, const ScalarBoundaryConditions &u_boundary_conditions,
    const ScalarBoundaryConditions &v_boundary_conditions,
    const PressureCorrectionBoundaryConditions &pressure_correction_boundary_conditions, FaceFluxField &mass_flux);

/// Optionally overwrites final u/v/p from one converged, independently generated
/// coarse mesh, using source boundary data and unlimited linear WLS reconstruction.
///
/// Zero mode leaves caller-initialized fields unchanged. Automatic coarse sizing
/// falls back unchanged if coarse meshing fails or produces a mesh that is not
/// smaller; an explicit target with either outcome is rejected. Invalid input
/// options and a non-converged coarse solve remain errors.
/// No face flux, pressure correction or linear system is transferred.
/// All coarse resources are local and released before this function returns;
/// the caller then constructs a single-Mesh solver for the final solve.
[[nodiscard]]
CoarseInitializationResult initialize_from_coarse_mesh(const Mesh &final_mesh, const input::MeshInput &mesh_input,
                                                       const input::InitializationInput &initialization,
                                                       const input::ScalarFieldInput &u_input,
                                                       const input::ScalarFieldInput &v_input,
                                                       const input::ScalarFieldInput &pressure_input, double density,
                                                       double dynamic_viscosity, ScalarConvectionScheme scheme,
                                                       const IncompressibleSimpleOptions &simple_options,
                                                       CellVelocityField &velocity, CellScalarField &pressure);

} // namespace cfd::app
