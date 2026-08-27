#pragma once

namespace cfd
{

struct RawMeshData;

/// Validates the structural consistency of raw mesh data.
///
/// Validation covers node coordinates, flattened cell connectivity, supported
/// cell cardinalities, referenced node IDs, boundary-group definitions, and
/// boundary-edge connectivity.
///
/// Geometric properties such as cell area, shape, and face orientation are
/// validated later during mesh construction.
///
/// @param raw_mesh Raw mesh representation to validate.
/// @throws std::runtime_error If any raw-mesh invariant is violated.
void validate_raw_mesh(const RawMeshData &raw_mesh);

} // namespace cfd