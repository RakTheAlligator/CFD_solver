#pragma once

namespace cfd
{

class CellScalarField;
class CellVectorField;
class Mesh;
class ScalarBoundaryConditions;

/// Reconstructs a cell-centered scalar gradient with inverse-distance weighted
/// least squares in two dimensions.
///
/// Each internal-cell and Dirichlet Taylor equation is normalized by its
/// center-to-center or center-to-face distance. Dirichlet values are applied at
/// boundary face centers. Neumann values prescribe the derivative along the
/// outward unit normal.
///
/// @param mesh Validated mesh providing topology and geometry.
/// @param field Cell-centered scalar values; its size must equal the mesh cell
///        count.
/// @param boundary_conditions One condition per mesh boundary group.
/// @param gradient Caller-owned output whose size must equal the mesh cell
///        count.
/// @throws std::invalid_argument If a field or boundary-condition cardinality
///         is incompatible with the mesh.
/// @throws std::runtime_error If a cell stencil cannot resolve a two-dimensional
///         gradient.
void compute_least_squares_gradient(const Mesh &mesh, const CellScalarField &field,
                                    const ScalarBoundaryConditions &boundary_conditions, CellVectorField &gradient);

} // namespace cfd
