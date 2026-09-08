#include "cfd/numerics/MassFluxBalance.hpp"

#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/FaceFluxField.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/Types.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace cfd
{

void compute_cell_mass_imbalance(const Mesh &mesh, const FaceFluxField &mass_flux, CellScalarField &imbalance)
{
    if (mass_flux.size() != mesh.face_count())
    {
        throw std::invalid_argument("Mass-flux size must match the balance Mesh face count.");
    }
    if (imbalance.size() != mesh.cell_count())
    {
        throw std::invalid_argument("Mass-imbalance size must match the balance Mesh cell count.");
    }

    const auto mass_flux_values{mass_flux.values()};
    for (Index face_id = 0; face_id < mesh.face_count(); ++face_id)
    {
        if (!std::isfinite(mass_flux_values[face_id]))
        {
            throw std::runtime_error("Mass-flux balance rejected non-finite face " + std::to_string(face_id) + ".");
        }
    }

    const auto cell_offsets{mesh.cell_node_offsets()};
    const auto cell_faces{mesh.cell_faces()};
    const auto face_adjacencies{mesh.face_adjacencies()};
    const auto evaluate_cell = [&](const Index cell_id) {
        double cell_imbalance{};
        for (Index offset = cell_offsets[cell_id]; offset < cell_offsets[cell_id + 1]; ++offset)
        {
            const Index face_id{cell_faces[offset]};
            const FaceAdjacency &adjacency{face_adjacencies[face_id]};
            const double oriented_flux{adjacency.owner == cell_id ? mass_flux_values[face_id]
                                                                  : -mass_flux_values[face_id]};
            cell_imbalance += oriented_flux;
            if (!std::isfinite(cell_imbalance))
            {
                throw std::runtime_error("Mass-flux balance is non-finite for cell " + std::to_string(cell_id) + ".");
            }
        }
        return cell_imbalance;
    };

    // Validate every result before overwriting caller-owned output storage.
    for (Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        static_cast<void>(evaluate_cell(cell_id));
    }

    auto imbalance_values{imbalance.values()};
    for (Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        imbalance_values[cell_id] = evaluate_cell(cell_id);
    }
}

} // namespace cfd
