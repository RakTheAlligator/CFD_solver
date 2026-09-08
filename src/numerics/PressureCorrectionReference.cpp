#include "cfd/numerics/PressureCorrectionReference.hpp"

#include "cfd/linear_algebra/ScalarLinearSystem.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"

#include <cmath>
#include <stdexcept>

namespace cfd
{

void apply_zero_pressure_correction_reference(const Index reference_cell, ScalarLinearSystem &system)
{
    if (reference_cell >= system.cell_count())
    {
        throw std::invalid_argument("Pressure-correction reference cell must be inside the linear system.");
    }
    if (!std::isfinite(system.diagonal()[reference_cell]) || !(system.diagonal()[reference_cell] > 0.0))
    {
        throw std::runtime_error("Pressure-correction reference diagonal must be finite and strictly positive.");
    }

    auto owner_neighbor_coefficients{system.owner_neighbor_coefficients()};
    auto neighbor_owner_coefficients{system.neighbor_owner_coefficients()};
    const auto face_adjacencies{system.mesh().face_adjacencies()};
    for (Index face_id = 0; face_id < system.face_count(); ++face_id)
    {
        const FaceAdjacency &adjacency{face_adjacencies[face_id]};
        if (adjacency.is_boundary() || (adjacency.owner != reference_cell && adjacency.neighbor != reference_cell))
        {
            continue;
        }

        owner_neighbor_coefficients[face_id] = 0.0;
        neighbor_owner_coefficients[face_id] = 0.0;
    }

    system.rhs()[reference_cell] = 0.0;
}

} // namespace cfd
