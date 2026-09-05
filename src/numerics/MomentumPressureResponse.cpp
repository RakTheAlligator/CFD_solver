#include "cfd/numerics/MomentumPressureResponse.hpp"

#include "cfd/field/CellMomentumPressureResponse.hpp"
#include "cfd/linear_algebra/ScalarLinearSystem.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/Types.hpp"

#include <cmath>
#include <stdexcept>

namespace cfd
{
namespace
{

void validate_system(const Mesh &mesh, const ScalarLinearSystem &system)
{
    if (&system.mesh() != &mesh)
    {
        throw std::invalid_argument("Momentum pressure-response systems must reference the supplied Mesh.");
    }
    if (system.cell_count() != mesh.cell_count() || system.face_count() != mesh.face_count())
    {
        throw std::invalid_argument("Momentum pressure-response system cardinalities must match the Mesh.");
    }
}

void validate_positive_finite(const double value, const char *const message)
{
    if (!std::isfinite(value) || !(value > 0.0))
    {
        throw std::runtime_error(message);
    }
}

} // namespace

void compute_momentum_pressure_response(const Mesh &mesh, const ScalarLinearSystem &u_momentum_system,
                                        const ScalarLinearSystem &v_momentum_system,
                                        CellMomentumPressureResponse &response)
{
    if (response.size() != mesh.cell_count())
    {
        throw std::invalid_argument("Momentum pressure-response cardinality must match the Mesh cell count.");
    }
    validate_system(mesh, u_momentum_system);
    validate_system(mesh, v_momentum_system);

    const auto cell_areas{mesh.cell_areas()};
    const auto u_diagonal{u_momentum_system.diagonal()};
    const auto v_diagonal{v_momentum_system.diagonal()};
    for (Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        const double area{cell_areas[cell_id]};
        const double u_diagonal_value{u_diagonal[cell_id]};
        const double v_diagonal_value{v_diagonal[cell_id]};
        validate_positive_finite(area, "Momentum pressure response requires finite positive cell areas.");
        validate_positive_finite(u_diagonal_value, "Momentum pressure response requires finite positive u diagonals.");
        validate_positive_finite(v_diagonal_value, "Momentum pressure response requires finite positive v diagonals.");
        validate_positive_finite(area / u_diagonal_value,
                                 "Computed u-momentum pressure response must be finite and positive.");
        validate_positive_finite(area / v_diagonal_value,
                                 "Computed v-momentum pressure response must be finite and positive.");
    }

    auto u_response{response.u().values()};
    auto v_response{response.v().values()};
    for (Index cell_id = 0; cell_id < mesh.cell_count(); ++cell_id)
    {
        u_response[cell_id] = cell_areas[cell_id] / u_diagonal[cell_id];
        v_response[cell_id] = cell_areas[cell_id] / v_diagonal[cell_id];
    }
}

} // namespace cfd
