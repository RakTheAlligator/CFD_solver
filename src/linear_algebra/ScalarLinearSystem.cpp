#include "cfd/linear_algebra/ScalarLinearSystem.hpp"

#include "cfd/mesh/Face.hpp"
#include "cfd/mesh/Mesh.hpp"

#include <algorithm>
#include <functional>
#include <stdexcept>

namespace cfd
{

namespace
{

[[nodiscard]]
bool spans_overlap(const std::span<const double> input, const std::span<double> output) noexcept
{
    if (input.empty() || output.empty())
    {
        return false;
    }

    const double *const input_begin{input.data()};
    const double *const input_end{input.data() + input.size()};
    const double *const output_begin{output.data()};
    const double *const output_end{output.data() + output.size()};
    const std::less<const double *> precedes;
    return precedes(input_begin, output_end) && precedes(output_begin, input_end);
}

} // namespace

ScalarLinearSystem::ScalarLinearSystem(const Mesh &mesh)
    : mesh_(&mesh), diagonal_(mesh.cell_count()), owner_neighbor_coefficients_(mesh.face_count()),
      neighbor_owner_coefficients_(mesh.face_count()), rhs_(mesh.cell_count())
{
}

Index ScalarLinearSystem::cell_count() const noexcept
{
    return diagonal_.size();
}

Index ScalarLinearSystem::face_count() const noexcept
{
    return owner_neighbor_coefficients_.size();
}

const Mesh &ScalarLinearSystem::mesh() const noexcept
{
    return *mesh_;
}

std::span<double> ScalarLinearSystem::diagonal() noexcept
{
    return diagonal_;
}

std::span<const double> ScalarLinearSystem::diagonal() const noexcept
{
    return diagonal_;
}

std::span<double> ScalarLinearSystem::owner_neighbor_coefficients() noexcept
{
    return owner_neighbor_coefficients_;
}

std::span<const double> ScalarLinearSystem::owner_neighbor_coefficients() const noexcept
{
    return owner_neighbor_coefficients_;
}

std::span<double> ScalarLinearSystem::neighbor_owner_coefficients() noexcept
{
    return neighbor_owner_coefficients_;
}

std::span<const double> ScalarLinearSystem::neighbor_owner_coefficients() const noexcept
{
    return neighbor_owner_coefficients_;
}

std::span<double> ScalarLinearSystem::rhs() noexcept
{
    return rhs_;
}

std::span<const double> ScalarLinearSystem::rhs() const noexcept
{
    return rhs_;
}

void ScalarLinearSystem::clear_matrix() noexcept
{
    std::fill(diagonal_.begin(), diagonal_.end(), 0.0);
    std::fill(owner_neighbor_coefficients_.begin(), owner_neighbor_coefficients_.end(), 0.0);
    std::fill(neighbor_owner_coefficients_.begin(), neighbor_owner_coefficients_.end(), 0.0);
}

void ScalarLinearSystem::clear_rhs() noexcept
{
    std::fill(rhs_.begin(), rhs_.end(), 0.0);
}

void ScalarLinearSystem::clear() noexcept
{
    clear_matrix();
    clear_rhs();
}

void ScalarLinearSystem::apply_matrix(const std::span<const double> input, const std::span<double> output) const
{
    if (input.size() != cell_count())
    {
        throw std::invalid_argument("Scalar linear system input size must match its cell count.");
    }
    if (output.size() != cell_count())
    {
        throw std::invalid_argument("Scalar linear system output size must match its cell count.");
    }
    if (spans_overlap(input, output))
    {
        throw std::invalid_argument("Scalar linear system input and output spans must not overlap.");
    }

    for (Index cell_id = 0; cell_id < cell_count(); ++cell_id)
    {
        output[cell_id] = diagonal_[cell_id] * input[cell_id];
    }

    const auto face_adjacencies{mesh_->face_adjacencies()};
    for (Index face_id = 0; face_id < face_count(); ++face_id)
    {
        const FaceAdjacency &adjacency{face_adjacencies[face_id]};
        if (adjacency.is_boundary())
        {
            continue;
        }

        output[adjacency.owner] += owner_neighbor_coefficients_[face_id] * input[adjacency.neighbor];
        output[adjacency.neighbor] += neighbor_owner_coefficients_[face_id] * input[adjacency.owner];
    }
}

} // namespace cfd
