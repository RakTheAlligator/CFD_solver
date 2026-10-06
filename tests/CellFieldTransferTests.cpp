#include "cfd/numerics/CellFieldTransfer.hpp"

#include "cfd/field/CellScalarField.hpp"
#include "cfd/field/CellVectorField.hpp"
#include "cfd/field/CellVelocityField.hpp"
#include "cfd/field/ScalarBoundaryConditions.hpp"
#include "cfd/mesh/Mesh.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/meshing/RawMeshData.hpp"

#include "support/TestUtils.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace
{

using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws;
using cfd::test::require_throws_with_message;

static_assert(!std::is_copy_constructible_v<cfd::CellFieldTransfer>);
static_assert(!std::is_copy_assignable_v<cfd::CellFieldTransfer>);
static_assert(std::is_nothrow_move_constructible_v<cfd::CellFieldTransfer>);
static_assert(!std::is_move_assignable_v<cfd::CellFieldTransfer>);

constexpr std::array source_coordinates{0.0, 0.5, 1.0};
constexpr std::array target_x{0.0, 0.2, 0.6, 1.0};
constexpr std::array target_y{0.0, 0.3, 0.7, 1.0};

cfd::RawMeshData make_grid(const std::span<const double> x, const std::span<const double> y,
                           const bool triangles = false, const double shear = 0.0)
{
    cfd::RawMeshData raw;
    const auto node = [columns = x.size()](const cfd::Index i, const cfd::Index j) { return j * columns + i; };
    for (const double py : y)
    {
        for (const double px : x)
        {
            raw.nodes.push_back({px + shear * py, py});
        }
    }
    raw.cell_node_offsets.push_back(0);
    for (cfd::Index j = 0; j + 1 < y.size(); ++j)
    {
        for (cfd::Index i = 0; i + 1 < x.size(); ++i)
        {
            const cfd::Index lower_left{node(i, j)};
            const cfd::Index lower_right{node(i + 1, j)};
            const cfd::Index upper_right{node(i + 1, j + 1)};
            const cfd::Index upper_left{node(i, j + 1)};
            if (triangles)
            {
                raw.cell_types.push_back(cfd::CellType::Triangle);
                raw.cell_nodes.insert(raw.cell_nodes.end(), {lower_left, lower_right, upper_right});
                raw.cell_node_offsets.push_back(raw.cell_nodes.size());
                raw.cell_types.push_back(cfd::CellType::Triangle);
                raw.cell_nodes.insert(raw.cell_nodes.end(), {lower_left, upper_right, upper_left});
            }
            else
            {
                raw.cell_types.push_back(cfd::CellType::Quadrilateral);
                raw.cell_nodes.insert(raw.cell_nodes.end(), {lower_left, lower_right, upper_right, upper_left});
            }
            raw.cell_node_offsets.push_back(raw.cell_nodes.size());
        }
    }
    const auto boundary = [&](const cfd::Index first, const cfd::Index second) {
        const cfd::BoundaryId id{raw.boundary_groups.size()};
        raw.boundary_groups.push_back({id, "edge_" + std::to_string(id)});
        raw.boundary_edges.push_back({{first, second}, id});
    };
    for (cfd::Index i = 0; i + 1 < x.size(); ++i)
    {
        boundary(node(i, 0), node(i + 1, 0));
        boundary(node(i, y.size() - 1), node(i + 1, y.size() - 1));
    }
    for (cfd::Index j = 0; j + 1 < y.size(); ++j)
    {
        boundary(node(0, j), node(0, j + 1));
        boundary(node(x.size() - 1, j), node(x.size() - 1, j + 1));
    }
    return raw;
}

double affine_value(const cfd::Point2 &point, const cfd::Vector2 gradient, const double offset)
{
    return gradient.x * point.x + gradient.y * point.y + offset;
}

void fill_affine(const cfd::Mesh &mesh, cfd::CellScalarField &field, const cfd::Vector2 gradient, const double offset)
{
    for (cfd::Index cell = 0; cell < mesh.cell_count(); ++cell)
    {
        field[cell] = affine_value(mesh.cell_centers()[cell], gradient, offset);
    }
}

cfd::ScalarBoundaryConditions affine_boundary_data(const cfd::Mesh &mesh, const cfd::Vector2 gradient,
                                                   const double offset)
{
    std::vector<cfd::ScalarBoundaryCondition> conditions(mesh.boundary_groups().size(),
                                                         {cfd::ScalarBoundaryConditionType::Dirichlet, 0.0});
    for (cfd::Index face = 0; face < mesh.face_count(); ++face)
    {
        if (mesh.face_adjacencies()[face].is_boundary())
        {
            conditions[mesh.face_boundary_ids()[face]].value =
                affine_value(mesh.face_centers()[face], gradient, offset);
        }
    }
    return {mesh.boundary_groups().size(), std::move(conditions)};
}

void require_affine(const cfd::Mesh &mesh, const cfd::CellScalarField &field, const cfd::Vector2 gradient,
                    const double offset)
{
    for (cfd::Index cell = 0; cell < mesh.cell_count(); ++cell)
    {
        require_near(field[cell], affine_value(mesh.cell_centers()[cell], gradient, offset), 1.0e-12,
                     "Linear reconstruction did not reproduce an affine field.");
    }
}

void test_constant_and_affine_on_independent_triangles_and_quads()
{
    for (const bool triangles : {false, true})
    {
        auto source{cfd::build_mesh(make_grid(source_coordinates, source_coordinates, triangles, 0.3))};
        auto target{cfd::build_mesh(make_grid(target_x, target_y, !triangles, 0.3))};
        require(source.mesh.nodes().data() != target.mesh.nodes().data(), "Meshes are not independent.");
        const cfd::CellFieldTransfer transfer{source.mesh, target.mesh};
        cfd::CellScalarField input{source.mesh.cell_count(), 7.0};
        cfd::CellScalarField output{target.mesh.cell_count(), -999.0};
        cfd::CellVectorField gradient{source.mesh.cell_count()};
        const auto constant_boundary{affine_boundary_data(source.mesh, {}, 7.0)};
        transfer.apply_piecewise_constant(input, output);
        for (const double value : output.values())
        {
            require(value == 7.0, "Piecewise constant transfer changed a constant field.");
        }
        transfer.apply_linear_reconstruction(input, constant_boundary, gradient, output);
        for (const double value : output.values())
        {
            require(value == 7.0, "Linear transfer changed a constant field.");
        }

        constexpr cfd::Vector2 exact_gradient{2.0, -3.0};
        fill_affine(source.mesh, input, exact_gradient, 4.0);
        const auto affine_boundary{affine_boundary_data(source.mesh, exact_gradient, 4.0)};
        transfer.apply_linear_reconstruction(input, affine_boundary, gradient, output);
        require_affine(target.mesh, output, exact_gradient, 4.0);
        for (const auto &value : gradient.values())
        {
            require_near(value.x, exact_gradient.x, 1.0e-12, "Source WLS x gradient changed.");
            require_near(value.y, exact_gradient.y, 1.0e-12, "Source WLS y gradient changed.");
        }
    }
}

void test_piecewise_constant_and_interface_tie()
{
    auto source{cfd::build_mesh(make_grid(source_coordinates, source_coordinates))};
    auto target{cfd::build_mesh(make_grid(target_x, target_y))};
    const cfd::CellFieldTransfer transfer{source.mesh, target.mesh};
    cfd::CellScalarField input{source.mesh.cell_count()};
    cfd::CellScalarField output{target.mesh.cell_count()};
    for (cfd::Index cell = 0; cell < input.size(); ++cell)
    {
        input[cell] = 10.0 + static_cast<double>(cell);
    }
    for (int repetition = 0; repetition < 3; ++repetition)
    {
        transfer.apply_piecewise_constant(input, output);
        for (cfd::Index cell = 0; cell < target.mesh.cell_count(); ++cell)
        {
            const auto &point{target.mesh.cell_centers()[cell]};
            // The middle target row lies exactly on the source interface y=0.5.
            const cfd::Index expected{static_cast<cfd::Index>(point.x > 0.5) +
                                      2 * static_cast<cfd::Index>(point.y > 0.5)};
            require(output[cell] == input[expected], "Containing cell or deterministic interface tie is incorrect.");
        }
    }
}

void test_reuses_mapping_and_workspace_for_velocity_and_pressure()
{
    auto source{cfd::build_mesh(make_grid(source_coordinates, source_coordinates))};
    auto target{cfd::build_mesh(make_grid(target_x, target_y))};
    const cfd::CellFieldTransfer transfer{source.mesh, target.mesh};
    cfd::CellVelocityField input_velocity{source.mesh.cell_count()};
    cfd::CellVelocityField output_velocity{target.mesh.cell_count()};
    cfd::CellScalarField input_pressure{source.mesh.cell_count()};
    cfd::CellScalarField output_pressure{target.mesh.cell_count()};
    cfd::CellVectorField workspace{source.mesh.cell_count()};
    const auto *workspace_storage{workspace.values().data()};
    const auto *output_storage{output_velocity.u().values().data()};
    fill_affine(source.mesh, input_velocity.u(), {2.0, -3.0}, 4.0);
    fill_affine(source.mesh, input_velocity.v(), {-0.5, 0.25}, 7.0);
    fill_affine(source.mesh, input_pressure, {1.0, 2.0}, -4.0);
    const auto u_boundary{affine_boundary_data(source.mesh, {2.0, -3.0}, 4.0)};
    const auto v_boundary{affine_boundary_data(source.mesh, {-0.5, 0.25}, 7.0)};
    const auto p_boundary{affine_boundary_data(source.mesh, {1.0, 2.0}, -4.0)};
    for (int repetition = 0; repetition < 3; ++repetition)
    {
        transfer.apply_linear_reconstruction(input_velocity.u(), u_boundary, workspace, output_velocity.u());
        transfer.apply_linear_reconstruction(input_velocity.v(), v_boundary, workspace, output_velocity.v());
        transfer.apply_linear_reconstruction(input_pressure, p_boundary, workspace, output_pressure);
        require_affine(target.mesh, output_velocity.u(), {2.0, -3.0}, 4.0);
        require_affine(target.mesh, output_velocity.v(), {-0.5, 0.25}, 7.0);
        require_affine(target.mesh, output_pressure, {1.0, 2.0}, -4.0);
        require(workspace.values().data() == workspace_storage && output_velocity.u().values().data() == output_storage,
                "Repeated application resized caller-owned buffers.");
    }
}

void test_tree_localization_on_multiple_leaves()
{
    std::vector<double> coarse(17);
    std::vector<double> fine(34);
    for (cfd::Index i = 0; i < coarse.size(); ++i)
    {
        coarse[i] = static_cast<double>(i) / static_cast<double>(coarse.size() - 1);
    }
    for (cfd::Index i = 0; i < fine.size(); ++i)
    {
        fine[i] = static_cast<double>(i) / static_cast<double>(fine.size() - 1);
    }
    auto source{cfd::build_mesh(make_grid(coarse, coarse))};
    auto target{cfd::build_mesh(make_grid(fine, fine))};
    const cfd::CellFieldTransfer transfer{source.mesh, target.mesh};
    cfd::CellScalarField input{source.mesh.cell_count()};
    cfd::CellScalarField output{target.mesh.cell_count()};
    for (cfd::Index i = 0; i < input.size(); ++i)
    {
        input[i] = static_cast<double>(i);
    }
    transfer.apply_piecewise_constant(input, output);
    for (cfd::Index i = 0; i < output.size(); ++i)
    {
        const auto &point{target.mesh.cell_centers()[i]};
        cfd::Index expected{cfd::invalid_index};
        // Small brute-force reference only in this test; production uses the tree.
        for (cfd::Index cell = 0; cell < source.mesh.cell_count(); ++cell)
        {
            const auto &center{source.mesh.cell_centers()[cell]};
            if (std::abs(point.x - center.x) <= 0.5 / 16.0 + 1.0e-14 &&
                std::abs(point.y - center.y) <= 0.5 / 16.0 + 1.0e-14)
            {
                expected = cell;
                break;
            }
        }
        require(expected != cfd::invalid_index && output[i] == input[expected],
                "Tree localization disagrees with the small exhaustive reference.");
    }
}

void test_boundary_roundoff_and_uncovered_centers()
{
    constexpr std::array unit{0.0, 1.0};
    auto source{cfd::build_mesh(make_grid(unit, unit))};
    // The target center is just outside x=1 by roundoff, not by a geometric gap.
    constexpr std::array edge_x{0.75, 1.25 + 8.0 * std::numeric_limits<double>::epsilon()};
    auto near_boundary{cfd::build_mesh(make_grid(edge_x, unit))};
    const cfd::CellFieldTransfer boundary_transfer{source.mesh, near_boundary.mesh};
    cfd::CellScalarField input{1, 8.0};
    cfd::CellScalarField output{1};
    boundary_transfer.apply_piecewise_constant(input, output);
    require(output[0] == 8.0, "Roundoff-distance boundary center was rejected.");

    constexpr std::array outside_x{1.0, 2.0};
    auto outside{cfd::build_mesh(make_grid(outside_x, unit))};
    require_throws_with_message<std::runtime_error>(
        [&] { const cfd::CellFieldTransfer transfer{source.mesh, outside.mesh}; }, "not covered",
        "An uncovered target center received a silent nearest-cell fallback.");

    cfd::RawMeshData triangle;
    triangle.nodes = {{0, 0}, {1, 0}, {0, 1}};
    triangle.cell_types = {cfd::CellType::Triangle};
    triangle.cell_nodes = {0, 1, 2};
    triangle.cell_node_offsets = {0, 3};
    triangle.boundary_groups = {{0, "wall"}};
    triangle.boundary_edges = {{{0, 1}, 0}, {{1, 2}, 0}, {{2, 0}, 0}};
    auto triangle_source{cfd::build_mesh(std::move(triangle))};
    constexpr std::array corner{0.6, 0.8};
    auto excluded_corner{cfd::build_mesh(make_grid(corner, corner))};
    require_throws_with_message<std::runtime_error>(
        [&] { const cfd::CellFieldTransfer transfer{triangle_source.mesh, excluded_corner.mesh}; }, "not covered",
        "Bounding-box membership was incorrectly used as polygon membership.");
}

void test_rejects_incompatible_fields_and_nonfinite_values()
{
    auto source{cfd::build_mesh(make_grid(source_coordinates, source_coordinates))};
    auto target{cfd::build_mesh(make_grid(target_x, target_y))};
    const cfd::CellFieldTransfer transfer{source.mesh, target.mesh};
    cfd::CellScalarField input{source.mesh.cell_count(), 2.0};
    cfd::CellScalarField output{target.mesh.cell_count(), -999.0};
    cfd::CellScalarField wrong{1};
    cfd::CellVectorField workspace{source.mesh.cell_count()};
    cfd::CellVectorField wrong_workspace{1};
    const auto boundary{affine_boundary_data(source.mesh, {}, 2.0)};
    const cfd::ScalarBoundaryConditions wrong_boundary{0, {}};
    require_throws<std::invalid_argument>([&] { transfer.apply_piecewise_constant(wrong, output); },
                                          "Wrong source cardinality was accepted.");
    require_throws<std::invalid_argument>([&] { transfer.apply_piecewise_constant(input, wrong); },
                                          "Wrong target cardinality was accepted.");
    require_throws<std::invalid_argument>(
        [&] { transfer.apply_linear_reconstruction(input, boundary, wrong_workspace, output); },
        "Wrong source gradient workspace was accepted.");
    require_throws<std::invalid_argument>(
        [&] { transfer.apply_linear_reconstruction(input, wrong_boundary, workspace, output); },
        "Wrong source boundary cardinality was accepted.");
    input[input.size() - 1] = std::numeric_limits<double>::quiet_NaN();
    require_throws<std::runtime_error>([&] { transfer.apply_piecewise_constant(input, output); },
                                       "A used non-finite source value was accepted.");
    require_throws<std::runtime_error>(
        [&] { transfer.apply_linear_reconstruction(input, boundary, workspace, output); },
        "Non-finite WLS input was accepted.");
    for (const double value : output.values())
    {
        require(value == -999.0, "Rejected transfer partially modified target values.");
    }
    const cfd::CellFieldTransfer self_transfer{source.mesh, source.mesh};
    require_throws<std::invalid_argument>([&] { self_transfer.apply_piecewise_constant(input, input); },
                                          "In-place constant transfer was accepted.");
    require_throws<std::invalid_argument>(
        [&] { self_transfer.apply_linear_reconstruction(input, boundary, workspace, input); },
        "In-place linear transfer was accepted.");
}

void test_rejects_replaced_source_and_target_meshes()
{
    auto source{cfd::build_mesh(make_grid(source_coordinates, source_coordinates))};
    auto target{cfd::build_mesh(make_grid(target_x, target_y))};
    cfd::CellScalarField input{source.mesh.cell_count(), 1.0};
    cfd::CellScalarField output{target.mesh.cell_count()};
    const cfd::CellFieldTransfer first{source.mesh, target.mesh};
    auto replacement_source{cfd::build_mesh(make_grid(source_coordinates, source_coordinates))};
    source.mesh = std::move(replacement_source.mesh);
    require_throws_with_message<std::invalid_argument>([&] { first.apply_piecewise_constant(input, output); },
                                                       "unchanged", "Replaced source mesh reused a stale mapping.");
    const cfd::CellFieldTransfer second{source.mesh, target.mesh};
    auto replacement_target{cfd::build_mesh(make_grid(target_x, target_y))};
    target.mesh = std::move(replacement_target.mesh);
    require_throws_with_message<std::invalid_argument>([&] { second.apply_piecewise_constant(input, output); },
                                                       "unchanged", "Replaced target mesh reused a stale mapping.");
}

} // namespace

int main()
{
    int failures{};
    failures += cfd::test::run_test("Independent triangles/quads: constant and affine transfer",
                                    test_constant_and_affine_on_independent_triangles_and_quads);
    failures += cfd::test::run_test("Piecewise constant localization and interface ties",
                                    test_piecewise_constant_and_interface_tie);
    failures += cfd::test::run_test("Reusable u/v/p transfer and WLS workspace",
                                    test_reuses_mapping_and_workspace_for_velocity_and_pressure);
    failures += cfd::test::run_test("Spatial tree multiple leaves", test_tree_localization_on_multiple_leaves);
    failures +=
        cfd::test::run_test("Boundary roundoff and uncovered centers", test_boundary_roundoff_and_uncovered_centers);
    failures += cfd::test::run_test("Field validation and non-finite transfer inputs",
                                    test_rejects_incompatible_fields_and_nonfinite_values);
    failures += cfd::test::run_test("Mesh replacement invalidation", test_rejects_replaced_source_and_target_meshes);
    return cfd::test::finish_tests(failures, "Cell field transfer");
}
