#pragma once

namespace cfd
{

class CellScalarField;
class FaceFluxField;
class Mesh;

/// Overwrites a cell field with the integrated outward mass-flux imbalance.
///
/// For every owner-oriented face flux, the owner receives `+F_f`, an internal
/// neighbor receives `-F_f`, and a boundary face contributes only to its
/// owner. The result is not normalized by cell area or volume.
///
/// All face values and computed cell imbalances are validated before the
/// caller-owned output is modified. Repeated valid calls are O(Nfaces +
/// Ncells), allocate no memory, and reuse the fixed output storage.
///
/// @throws std::invalid_argument If a field cardinality is incompatible with
///         the Mesh.
/// @throws std::runtime_error If a face flux or computed cell imbalance is
///         non-finite.
void compute_cell_mass_imbalance(const Mesh &mesh, const FaceFluxField &mass_flux, CellScalarField &imbalance);

} // namespace cfd
