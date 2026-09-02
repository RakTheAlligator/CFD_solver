#include "cfd/numerics/LeastSquaresGradient.hpp"

#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVectorField.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/math/Point2.hpp"
#include "cfd/math/Vector2.hpp"
#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/Types.hpp"

#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace cfd
{

namespace
{

[[noreturn]]
void throw_singular_stencil(const Index cell_id, const double m00, const double m01, const double m11,
                            const double determinant, const double trace)
{
    std::ostringstream message;
    message << "Least-squares gradient reconstruction failed for cell " << cell_id << ": singular matrix [m00=" << m00
            << ", m01=" << m01 << ", m11=" << m11 << ", determinant=" << determinant << ", trace=" << trace;

    const double trace_squared{trace * trace};
    if (trace > 0.0 && std::isfinite(trace_squared))
    {
        message << ", determinant_over_trace_squared=" << determinant / trace_squared;
    }

    message << "].";

    throw std::runtime_error(message.str());
}

} // namespace

void compute_least_squares_gradient(const Mesh &mesh, const CellScalarField &field,
                                    const ScalarBoundaryConditions &boundary_conditions, CellVectorField &gradient)
{
    const Index cell_count{mesh.cell_count()};

    if (field.size() != cell_count)
    {
        throw std::invalid_argument("Scalar field size must match the mesh cell count.");
    }

    if (gradient.size() != cell_count)
    {
        throw std::invalid_argument("Gradient field size must match the mesh cell count.");
    }

    if (boundary_conditions.size() != mesh.boundary_groups().size())
    {
        throw std::invalid_argument("Scalar boundary-condition count must match the mesh boundary-group count.");
    }

    // Cell-to-face connectivity shares these offsets by the current 2D Mesh invariant.
    const auto cell_offsets{mesh.cell_node_offsets()};
    const auto cell_faces{mesh.cell_faces()};
    const auto face_adjacencies{mesh.face_adjacencies()};
    const auto face_boundary_ids{mesh.face_boundary_ids()};
    const auto cell_centers{mesh.cell_centers()};
    const auto face_centers{mesh.face_centers()};
    const auto face_lengths{mesh.face_lengths()};
    const auto face_area_vectors{mesh.face_area_vectors()};
    const auto field_values{field.values()};
    auto gradient_values{gradient.values()};

    constexpr double relative_singularity_tolerance{64.0 * std::numeric_limits<double>::epsilon()};

    for (Index cell_id = 0; cell_id < cell_count; ++cell_id)
    {
        double m00{};
        double m01{};
        double m11{};
        double rhs_x{};
        double rhs_y{};

        const Point2 &cell_center{cell_centers[cell_id]};
        const double cell_value{field_values[cell_id]};

        const Index cell_face_begin_offset{cell_offsets[cell_id]};
        const Index cell_face_end_offset{cell_offsets[cell_id + 1]};

        for (Index cell_face_position = cell_face_begin_offset; cell_face_position < cell_face_end_offset;
             ++cell_face_position)
        {
            const Index face_id{cell_faces[cell_face_position]};
            const FaceAdjacency &adjacency{face_adjacencies[face_id]};

            if (!adjacency.is_boundary())
            {
                const Index other_cell_id{adjacency.owner == cell_id ? adjacency.neighbor : adjacency.owner};
                const Point2 &other_center{cell_centers[other_cell_id]};

                const double dx{other_center.x - cell_center.x};
                const double dy{other_center.y - cell_center.y};
                const double squared_distance{dx * dx + dy * dy};
                const double inverse_squared_distance{1.0 / squared_distance};
                const double value_difference{field_values[other_cell_id] - cell_value};

                m00 += dx * dx * inverse_squared_distance;
                m01 += dx * dy * inverse_squared_distance;
                m11 += dy * dy * inverse_squared_distance;

                rhs_x += dx * value_difference * inverse_squared_distance;
                rhs_y += dy * value_difference * inverse_squared_distance;

                continue;
            }

            const BoundaryId boundary_id{face_boundary_ids[face_id]};
            const ScalarBoundaryCondition &condition{boundary_conditions[boundary_id]};

            switch (condition.type)
            {
            case ScalarBoundaryConditionType::Dirichlet: {
                const Point2 &face_center{face_centers[face_id]};

                const double dx{face_center.x - cell_center.x};
                const double dy{face_center.y - cell_center.y};
                const double squared_distance{dx * dx + dy * dy};
                const double inverse_squared_distance{1.0 / squared_distance};
                const double value_difference{condition.value - cell_value};

                m00 += dx * dx * inverse_squared_distance;
                m01 += dx * dy * inverse_squared_distance;
                m11 += dy * dy * inverse_squared_distance;

                rhs_x += dx * value_difference * inverse_squared_distance;
                rhs_y += dy * value_difference * inverse_squared_distance;
                break;
            }

            case ScalarBoundaryConditionType::Neumann: {
                const Vector2 &area_vector{face_area_vectors[face_id]};
                const double inverse_face_length{1.0 / face_lengths[face_id]};
                const double normal_x{area_vector.x * inverse_face_length};
                const double normal_y{area_vector.y * inverse_face_length};

                m00 += normal_x * normal_x;
                m01 += normal_x * normal_y;
                m11 += normal_y * normal_y;

                rhs_x += normal_x * condition.value;
                rhs_y += normal_y * condition.value;
                break;
            }
            }
        }

        const double determinant{m00 * m11 - m01 * m01};
        const double trace{m00 + m11};
        const double determinant_tolerance{relative_singularity_tolerance * trace * trace};

        if (!(determinant > determinant_tolerance))
        {
            throw_singular_stencil(cell_id, m00, m01, m11, determinant, trace);
        }

        gradient_values[cell_id] = {
            (rhs_x * m11 - rhs_y * m01) / determinant,
            (m00 * rhs_y - m01 * rhs_x) / determinant,
        };
    }
}

} // namespace cfd
