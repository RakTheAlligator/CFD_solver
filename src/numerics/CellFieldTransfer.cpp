#include "cfd/numerics/CellFieldTransfer.hpp"

#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVectorField.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/numerics/LeastSquaresGradient.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <span>
#include <stdexcept>
#include <string>

namespace cfd
{
namespace
{

struct BoundingBox
{
    Point2 minimum{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
    Point2 maximum{-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};

    void include(const Point2 &point) noexcept
    {
        minimum.x = std::min(minimum.x, point.x);
        minimum.y = std::min(minimum.y, point.y);
        maximum.x = std::max(maximum.x, point.x);
        maximum.y = std::max(maximum.y, point.y);
    }

    [[nodiscard]]
    bool contains(const Point2 &point) const noexcept
    {
        return point.x >= minimum.x && point.x <= maximum.x && point.y >= minimum.y && point.y <= maximum.y;
    }
};

struct SourceCell
{
    Index id{};
    BoundingBox box;
    double tolerance{};
};

struct SearchNode
{
    BoundingBox box;
    Index begin{};
    Index end{};
    Index left{invalid_index};
    Index right{invalid_index};
};

// Construction-only spatial index; no search tree persists in the transfer object.
class SourceCellSearch
{
  public:
    explicit SourceCellSearch(const Mesh &mesh) : mesh_(&mesh)
    {
        cells_.reserve(mesh.cell_count());
        nodes_.reserve(mesh.cell_count());
        const auto offsets{mesh.cell_node_offsets()};
        const auto connectivity{mesh.cell_nodes()};
        const auto points{mesh.nodes()};
        for (Index cell = 0; cell < mesh.cell_count(); ++cell)
        {
            SourceCell entry;
            entry.id = cell;
            for (Index position = offsets[cell]; position < offsets[cell + 1]; ++position)
            {
                entry.box.include(points[connectivity[position]]);
            }
            constexpr double relative_tolerance{64.0 * std::numeric_limits<double>::epsilon()};
            const double scale{
                std::max({std::abs(entry.box.minimum.x), std::abs(entry.box.minimum.y), std::abs(entry.box.maximum.x),
                          std::abs(entry.box.maximum.y), entry.box.maximum.x - entry.box.minimum.x,
                          entry.box.maximum.y - entry.box.minimum.y})};
            entry.tolerance = relative_tolerance * scale;
            entry.box.minimum.x -= entry.tolerance;
            entry.box.minimum.y -= entry.tolerance;
            entry.box.maximum.x += entry.tolerance;
            entry.box.maximum.y += entry.tolerance;
            cells_.push_back(entry);
        }
        static_cast<void>(build(0, cells_.size()));
    }

    [[nodiscard]]
    Index locate(const Point2 &point) const noexcept
    {
        return locate_in_node(point, 0);
    }

  private:
    Index build(const Index begin, const Index end)
    {
        const Index node_id{nodes_.size()};
        nodes_.emplace_back();
        BoundingBox box;
        for (Index position = begin; position < end; ++position)
        {
            box.include(cells_[position].box.minimum);
            box.include(cells_[position].box.maximum);
        }
        nodes_[node_id].box = box;
        nodes_[node_id].begin = begin;
        nodes_[node_id].end = end;
        constexpr Index maximum_leaf_cells{8};
        if (end - begin > maximum_leaf_cells)
        {
            const bool split_x{box.maximum.x - box.minimum.x >= box.maximum.y - box.minimum.y};
            const Index middle{begin + (end - begin) / 2};
            const std::span<SourceCell> cells{cells_};
            auto range{cells.subspan(begin, end - begin)};
            std::nth_element(
                range.begin(), range.begin() + static_cast<std::ptrdiff_t>(middle - begin), range.end(),
                [split_x](const SourceCell &first, const SourceCell &second) {
                    const double first_position{split_x ? std::midpoint(first.box.minimum.x, first.box.maximum.x)
                                                        : std::midpoint(first.box.minimum.y, first.box.maximum.y)};
                    const double second_position{split_x ? std::midpoint(second.box.minimum.x, second.box.maximum.x)
                                                         : std::midpoint(second.box.minimum.y, second.box.maximum.y)};
                    return first_position == second_position ? first.id < second.id : first_position < second_position;
                });
            const Index left{build(begin, middle)};
            const Index right{build(middle, end)};
            nodes_[node_id].left = left;
            nodes_[node_id].right = right;
        }
        return node_id;
    }

    bool contains(const SourceCell &cell, const Point2 &point) const noexcept
    {
        const auto offsets{mesh_->cell_node_offsets()};
        const auto faces{mesh_->cell_faces()};
        const auto adjacencies{mesh_->face_adjacencies()};
        const auto centers{mesh_->face_centers()};
        const auto lengths{mesh_->face_lengths()};
        const auto vectors{mesh_->face_area_vectors()};
        // Source triangles/quads are convex. Every outward half-plane must contain the point.
        for (Index position = offsets[cell.id]; position < offsets[cell.id + 1]; ++position)
        {
            const Index face{faces[position]};
            const double orientation{adjacencies[face].owner == cell.id ? 1.0 : -1.0};
            const double normal_x{orientation * vectors[face].x / lengths[face]};
            const double normal_y{orientation * vectors[face].y / lengths[face]};
            const double distance{normal_x * (point.x - centers[face].x) + normal_y * (point.y - centers[face].y)};
            if (!(distance <= cell.tolerance))
            {
                return false;
            }
        }
        return true;
    }

    Index locate_in_node(const Point2 &point, const Index node_id) const noexcept
    {
        const SearchNode &node{nodes_[node_id]};
        if (!node.box.contains(point))
        {
            return invalid_index;
        }
        if (node.left != invalid_index)
        {
            return std::min(locate_in_node(point, node.left), locate_in_node(point, node.right));
        }
        Index result{invalid_index};
        for (Index position = node.begin; position < node.end; ++position)
        {
            const SourceCell &cell{cells_[position]};
            if (cell.box.contains(point) && contains(cell, point))
            {
                result = std::min(result, cell.id);
            }
        }
        return result;
    }

    const Mesh *mesh_;
    std::vector<SourceCell> cells_;
    std::vector<SearchNode> nodes_;
};

} // namespace

CellFieldTransfer::CellFieldTransfer(const Mesh &source_mesh, const Mesh &target_mesh)
    : source_mesh_(&source_mesh), target_mesh_(&target_mesh), source_centers_(source_mesh.cell_centers().data()),
      target_centers_(target_mesh.cell_centers().data()), source_cell_count_(source_mesh.cell_count())
{
    if (source_mesh.cell_count() == 0 || target_mesh.cell_count() == 0)
    {
        throw std::invalid_argument("Cell field transfer requires nonempty source and target meshes.");
    }
    const SourceCellSearch search{source_mesh};
    entries_.resize(target_mesh.cell_count());
    for (Index target = 0; target < entries_.size(); ++target)
    {
        const Point2 &point{target_mesh.cell_centers()[target]};
        const Index source{search.locate(point)};
        if (source == invalid_index)
        {
            throw std::runtime_error("Cell field transfer: target cell " + std::to_string(target) +
                                     " center is not covered by the source mesh.");
        }
        const Point2 &center{source_mesh.cell_centers()[source]};
        entries_[target] = {source, {point.x - center.x, point.y - center.y}};
    }
}

void CellFieldTransfer::validate_fields(const CellScalarField &source, const CellScalarField &target) const
{
    if (source_mesh_->cell_count() != source_cell_count_ || target_mesh_->cell_count() != entries_.size() ||
        source_mesh_->cell_centers().data() != source_centers_ ||
        target_mesh_->cell_centers().data() != target_centers_)
    {
        throw std::invalid_argument("Cell field transfer requires the original unchanged source and target meshes.");
    }
    if (source.size() != source_cell_count_ || target.size() != entries_.size())
    {
        throw std::invalid_argument("Cell field transfer fields must match source and target cell counts.");
    }
    if (&source == &target)
    {
        throw std::invalid_argument("Cell field transfer source and target fields must be distinct.");
    }
}

void CellFieldTransfer::apply_piecewise_constant(const CellScalarField &source, CellScalarField &target) const
{
    validate_fields(source, target);
    const auto values{source.values()};
    auto output{target.values()};
    for (const Entry &entry : entries_)
    {
        if (!std::isfinite(values[entry.source_cell]))
        {
            throw std::runtime_error("Cell field transfer used source values must be finite.");
        }
    }
    for (Index cell = 0; cell < entries_.size(); ++cell)
    {
        output[cell] = values[entries_[cell].source_cell];
    }
}

void CellFieldTransfer::apply_linear_reconstruction(const CellScalarField &source,
                                                    const ScalarBoundaryConditions &source_boundary_conditions,
                                                    CellVectorField &source_gradient_workspace,
                                                    CellScalarField &target) const
{
    validate_fields(source, target);
    if (source_gradient_workspace.size() != source_cell_count_ ||
        source_boundary_conditions.size() != source_mesh_->boundary_groups().size())
    {
        throw std::invalid_argument("Cell field transfer WLS workspace and boundary data must match the source mesh.");
    }
    const auto values{source.values()};
    for (const double value : values)
    {
        if (!std::isfinite(value))
        {
            throw std::runtime_error("Cell field transfer source values must be finite for WLS.");
        }
    }
    compute_least_squares_gradient(*source_mesh_, source, source_boundary_conditions, source_gradient_workspace);
    const auto gradients{source_gradient_workspace.values()};
    const auto reconstructed_value = [&](const Entry &entry) {
        const Vector2 &gradient{gradients[entry.source_cell]};
        return values[entry.source_cell] + gradient.x * entry.displacement.x + gradient.y * entry.displacement.y;
    };
    // Validate before writing without allocating a target-sized temporary buffer.
    for (const Entry &entry : entries_)
    {
        if (!std::isfinite(reconstructed_value(entry)))
        {
            throw std::runtime_error("Cell field transfer reconstructed values must be finite.");
        }
    }
    auto output{target.values()};
    for (Index cell = 0; cell < entries_.size(); ++cell)
    {
        output[cell] = reconstructed_value(entries_[cell]);
    }
}

} // namespace cfd
