#pragma once

#include "cfd/mesh/Types.hpp"

#include <cstdint>
#include <functional>
#include <optional>
#include <string>
#include <vector>

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

struct GridLevelPlan
{
    Index target_cell_count{};
    double mesh_size{};
};

/// Plans three geometrically spaced levels, including the existing final Mesh.
/// The default sizes are {4h, 2h, h}. An explicit positive target smaller than
/// the final count sets the coarsest size; cell counts are reporting targets.
/// Zero mode returns an empty plan without changing historical initialization.
[[nodiscard]]
std::vector<GridLevelPlan> make_grid_sequencing_plan(Index final_cell_count, double final_mesh_size,
                                                     const input::InitializationInput &initialization);

struct GridLevelReport
{
    GridLevelPlan plan;
    Index actual_cell_count{};
    /// Present only after a successful auxiliary SIMPLE solve, including zero iterations.
    std::optional<Index> iteration_count{};
    double solve_seconds{};
    /// Outgoing mapping construction plus WLS/transfer of u/v/p to the next level.
    double transfer_seconds{};
};

struct GridSequencingResult
{
    bool used_grid_sequencing{};
    bool automatic_fallback{};
    std::string fallback_reason;
    /// Includes the final level, whose SIMPLE solve remains the caller's responsibility.
    std::vector<GridLevelReport> levels;
};

/// Synchronous observer of a converged auxiliary level, after its report is filled.
/// All references are borrowed only for this call and must not be retained.
/// The flux is the level's corrected mass flux, not a quantity to transfer.
using GridLevelConvergedCallback =
    std::function<void(Index, const Mesh &, const CellVelocityField &, const CellScalarField &, const FaceFluxField &)>;

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

/// Initializes final u/v/p through a loop of converged auxiliary single-Mesh solves.
/// Only level zero uses boundary-based bulk initialization; subsequent levels
/// receive unlimited linear WLS reconstruction with the preceding level's BCs.
///
/// Zero mode leaves caller-initialized fields unchanged. Automatic mode falls
/// back unchanged if meshing, strictly increasing cell counts or geometric
/// transfer coverage fail; explicit targets turn these failures into errors.
/// Invalid input, non-converged solves and numerical transfer failures are errors.
/// No face flux, pressure correction or linear system is transferred.
/// Auxiliary solvers are destroyed before transfer. Only the converged source
/// and receiving target states coexist; source resources are released before
/// the next solve. The existing final Mesh is borrowed and never regenerated.
/// All auxiliary resources are released before returning to the final caller.
/// The optional observer runs once per converged auxiliary level, before any
/// outgoing transfer. It is never called for zero mode or the final Mesh.
/// Observer exceptions propagate without fallback or cleanup of prior exports.
[[nodiscard]]
GridSequencingResult initialize_with_grid_sequencing(
    const Mesh &final_mesh, const input::MeshInput &mesh_input, const input::InitializationInput &initialization,
    const input::ScalarFieldInput &u_input, const input::ScalarFieldInput &v_input,
    const input::ScalarFieldInput &pressure_input, double density, double dynamic_viscosity,
    ScalarConvectionScheme scheme, const IncompressibleSimpleOptions &simple_options, CellVelocityField &velocity,
    CellScalarField &pressure, const GridLevelConvergedCallback &level_callback = {});

} // namespace cfd::app
