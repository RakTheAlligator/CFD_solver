#include "support/EulerBoundaryWave.hpp"
#include "support/EulerRpaExploratory.hpp"
#include "support/EulerRpaTemporal.hpp"
#include "support/TestUtils.hpp"

#include "cfd/io/VtkWriter.hpp"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <limits>
#include <numbers>
#include <sstream>

namespace
{
using cfd::test::require;
using cfd::test::require_near;
using cfd::test::require_throws;
using Model = cfd::verification::EulerRpaModel;

void test_nominal_derived_and_exploratory_states()
{
    const cfd::CaloricallyPerfectGas gas{1.177, 414.2};
    require_near(Model::total_temperature(), 3350.294882116409, 2e-14, "Derived T0.");
    require_near(Model::total_pressure(), 1185747.2501228245, 2e-14, "Derived P0.");
    const auto hot{gas.evaluate(gas.to_conservative(Model::hot()))};
    require_near(hot.primitive.pressure, 1170000, 2e-14, "Hot pressure.");
    require_near(hot.thermodynamics.temperature, 3343.5658, 2e-14, "Hot temperature.");
    const auto ambient{gas.evaluate(gas.to_conservative(Model::ambient(300, 0, 0)))};
    require(std::abs(ambient.primitive.pressure / 1e-4 - 1.0) < 2e-14, "Nominal ambient pressure was altered.");
    require_near(ambient.thermodynamics.temperature, 300.0, 2e-14, "Assumed ambient temperature.");
    for (const double invalid :
         {0.0, -1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
        require_throws<std::invalid_argument>([&] { (void)Model::ambient(invalid, 0, 0); }, "Invalid ambient T.");
    require_throws<std::invalid_argument>([&] { (void)Model::ambient(300, INFINITY, 0); }, "Invalid ambient velocity.");
}

void test_eight_boundary_groups_and_no_silent_mapping()
{
    std::vector<cfd::BoundaryGroup> groups{{0, "inlet"},           {1, "axis"},
                                           {2, "inner_wall"},      {3, "outer_wall"},
                                           {4, "nozzle_lip"},      {5, "farfield_upstream"},
                                           {6, "farfield_radial"}, {7, "farfield_downstream"}};
    const auto conditions{cfd::verification::rpa_euler_conditions(groups, Model::ambient(300, 0, 0))};
    using Type = cfd::EulerBoundaryConditionType;
    require(conditions[0] == Type::TotalConditionsInlet && conditions[1] == Type::Axis, "Inlet/axis mapping.");
    for (cfd::Index id = 2; id < 5; ++id)
        require(conditions[id] == Type::SlipWall, "Solid boundary mapping.");
    for (cfd::Index id = 5; id < 8; ++id)
        require(conditions[id] == Type::Farfield, "Open boundary mapped as wall/outlet.");
    groups.back().name = "farfield_radial";
    require_throws<std::invalid_argument>(
        [&] { (void)cfd::verification::rpa_euler_conditions(groups, Model::ambient(300, 0, 0)); }, "Duplicate group.");
    groups.back().name = "invented";
    require_throws<std::invalid_argument>(
        [&] { (void)cfd::verification::rpa_euler_conditions(groups, Model::ambient(300, 0, 0)); }, "Unknown group.");
    groups.pop_back();
    require_throws<std::invalid_argument>(
        [&] { (void)cfd::verification::rpa_euler_conditions(groups, Model::ambient(300, 0, 0)); }, "Missing group.");
}

void test_region_initialization_and_conservative_first_step()
{
    auto raw{cfd::verification::make_euler_tube(4, 1)};
    for (auto &node : raw.nodes)
    {
        node.x *= 2;
        node.y *= 0.02;
    }
    const auto built{cfd::build_mesh(std::move(raw))};
    cfd::verification::require_connected_rpa(built.mesh);
    const cfd::RpaNozzleDomain domain{.inner_wall = {{0, .1}, {.5, .05}, {1, .1}},
                                      .outer_wall = {{0, .12}, {1, .12}},
                                      .curve_junctions = {},
                                      .throat = {.5, .05},
                                      .downstream_x = 2,
                                      .farfield_radius = 1};
    const cfd::RpaRegionLayout layout{domain, {.band_thickness_factor = .5, .plume_outlet_radius_factor = 6}};
    cfd::ConservativeField field{built.mesh.cell_count()};
    std::vector<double> regions(field.size());
    const auto ambient{Model::ambient(300, 0, 0)};
    require(cfd::verification::initialize_rpa(built.mesh, layout, ambient, field, regions) == 2,
            "Engine/exterior partition.");
    const cfd::CaloricallyPerfectGas gas{Model::gamma, Model::gas_constant};
    for (cfd::Index cell = 0; cell < field.size(); ++cell)
    {
        const auto state{gas.evaluate(field.state(cell))};
        require(state.primitive.density > 0 && state.primitive.pressure > 0 && state.thermodynamics.temperature > 0,
                "Initial admissibility.");
        require(regions[cell] == (cell < 2 ? 1.0 : 3.0), "Actual layout regions.");
    }
    cfd::AxisymmetricEulerSolver solver{
        built.mesh,
        field,
        gas,
        {4,
         {cfd::EulerBoundaryConditionType::Axis, cfd::EulerBoundaryConditionType::SlipWall,
          cfd::EulerBoundaryCondition::total_conditions_inlet(Model::total_pressure(), Model::total_temperature(),
                                                              {1, 0}),
          cfd::EulerBoundaryCondition::farfield(ambient)}},
        {.cfl = .01}};
    (void)solver.evaluate_steady_residual(); // Includes initial characteristic traces and every physical area.
    cfd::verification::EulerBoundaryWaveBalance balance{.initial = solver.diagnostics().integrals};
    balance.add(solver.step(1e-8));
    require(solver.diagnostics().accepted_steps == 1 && solver.diagnostics().rejected_steps == 0,
            "One accepted update.");
    const auto final{solver.diagnostics().integrals}, first{balance.initial};
    const std::array initial_values{first.mass, first.momentum_x, first.momentum_r, first.energy};
    const std::array final_values{final.mass, final.momentum_x, final.momentum_r, final.energy};
    const auto defects{balance.defect(final)};
    for (cfd::Index k = 0; k < defects.size(); ++k)
    {
        const long double scale{std::abs(initial_values.at(k)) + std::abs(final_values.at(k)) +
                                std::abs(balance.integrated_boundary.at(k)) +
                                (k == 2 ? std::abs(balance.integrated_source) : 0.0L)};
        require(std::abs(defects.at(k)) <= 512 * std::numeric_limits<double>::epsilon() * scale,
                "Component-scaled four-equation balance including source.");
    }
}

void test_output_protection_including_symlink()
{
    const auto path{std::filesystem::temp_directory_path() /
                    std::filesystem::path{"cfd-rpa-protect-" +
                                          std::to_string(std::chrono::steady_clock::now().time_since_epoch().count())}};
    std::filesystem::create_directories(path / "private");
    const auto input{path / "private"};
    std::ofstream{input / "sentinel.txt"} << "preserve";
    require_throws<std::invalid_argument>(
        [&] { cfd::verification::require_external_empty_output(input, input / "results"); }, "Private child output.");
    std::filesystem::create_directory_symlink(input, path / "alias");
    require_throws<std::invalid_argument>(
        [&] { cfd::verification::require_external_empty_output(input, path / "alias" / "results"); },
        "Symlink private output.");
    std::filesystem::create_directories(path / "old");
    std::ofstream{path / "old" / "sentinel.txt"} << "previous run";
    require_throws<std::invalid_argument>(
        [&] { cfd::verification::require_external_empty_output(input, path / "old"); },
        "Old results must not be overwritten.");
    cfd::verification::require_external_empty_output(input, path / "new output with spaces");
    require(cfd::test::read_text_file(input / "sentinel.txt") == "preserve", "Private data changed.");
    require(cfd::test::read_text_file(path / "old" / "sentinel.txt") == "previous run", "Old results changed.");
    std::filesystem::remove_all(path); // Only this test's uniquely created temporary fixture.
}

void test_acoustic_scale_physical_outputs_and_limits()
{
    using cfd::verification::RpaPhysicalOutputs;
    require(std::abs(cfd::verification::rpa_nozzle_acoustic_time(.650175702) / .000509251737436 - 1) < 2e-14,
            "Independently recorded nozzle acoustic scale.");
    for (const double invalid : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                                 std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::denorm_min()})
    {
        require_throws<std::invalid_argument>([&] { (void)RpaPhysicalOutputs{invalid}; }, "Invalid output target.");
        require_throws<std::invalid_argument>([&] { (void)cfd::verification::rpa_nozzle_acoustic_time(invalid); },
                                              "Invalid/unrepresentable acoustic scale.");
    }
    RpaPhysicalOutputs schedule{1e-8};
    require(!schedule.due(0), "No duplicate initial physical export.");
    require_throws<std::runtime_error>([&] { schedule.exported(0); }, "Cannot label an early state with a target.");
    require_throws<std::runtime_error>([&] { schedule.exported(1e-8); }, "Cannot silently skip targets.");
    const auto built{cfd::verification::build_boundary_wave_mesh(4, 2, 1, .2)};
    const cfd::CaloricallyPerfectGas gas{Model::gamma, Model::gas_constant};
    const cfd::ConservativeField field{built.mesh.cell_count(), gas.to_conservative(Model::hot())};
    const auto conditions{cfd::verification::boundary_wave_conditions({.gamma = Model::gamma,
                                                                       .gas_constant = Model::gas_constant,
                                                                       .density = Model::hot().density,
                                                                       .pressure = Model::pressure,
                                                                       .mach = Model::mach})};
    cfd::AxisymmetricEulerSolver solver{built.mesh, field, gas, conditions, {.cfl = .01, .maximum_steps = 5}};
    double previous{};
    for (const double target : schedule.times())
    {
        const auto &accepted{solver.step(schedule.next_time())};
        require(accepted.time == target && accepted.time > previous, "Actual accepted, strictly ordered output time.");
        require(schedule.due(accepted.time), "Export occurs only at its target.");
        schedule.exported(accepted.time);
        previous = accepted.time;
    }
    require(!schedule.due(previous) && solver.diagnostics().accepted_steps == 5, "Final export is not duplicated.");
    RpaPhysicalOutputs comparison{1.0, .2};
    constexpr std::array expected{.1, .2, .25, .5, .75, 1.0};
    require(std::ranges::equal(comparison.times(), expected), "Extra comparison time is inserted in physical order.");
    for (const double time : expected)
    {
        require(comparison.next_time() == time, "Next clipped target includes the comparison time.");
        comparison.exported(time);
    }
    require(!comparison.due(1), "Additional export does not duplicate the last state.");
    require(RpaPhysicalOutputs{1.0, .25}.times().size() == 5, "Existing comparison time is not duplicated.");
    for (const double invalid :
         {0.0, -1.0, 1.0, 2.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
        require_throws<std::invalid_argument>([&] { (void)RpaPhysicalOutputs{1.0, invalid}; },
                                              "Invalid comparison time.");
    cfd::AxisymmetricEulerSolver capped{built.mesh, field, gas, conditions, {.cfl = .01, .maximum_steps = 1}};
    require(capped.step(1e-4).time < 1e-4, "A step cap must not be reported as reaching the physical target.");
    require_throws<std::runtime_error>([&] { capped.step(1e-4); }, "The production step cap remains enforced.");
}

void test_face_aligned_section_geometry_and_numerical_flow()
{
    const auto built{cfd::verification::build_boundary_wave_mesh(4, 3, 1, .2)};
    const cfd::AxisymmetricGeometry geometry{built.mesh};
    const cfd::CaloricallyPerfectGas gas{Model::gamma, Model::gas_constant};
    cfd::ConservativeField field{built.mesh.cell_count(), gas.to_conservative(Model::hot())};
    std::vector<double> regions(field.size(), 1);
    for (cfd::Index cell = 0; cell < field.size(); ++cell)
        if (built.mesh.cell_centers()[cell].x > .75)
            regions[cell] = 3;
    const auto conditions{cfd::verification::boundary_wave_conditions({.gamma = Model::gamma,
                                                                       .gas_constant = Model::gas_constant,
                                                                       .density = Model::hot().density,
                                                                       .pressure = Model::pressure,
                                                                       .mach = Model::mach})};
    const double exact_area{std::numbers::pi * .04};
    for (const double x : {0.0, .5, .75})
    {
        const cfd::verification::RpaFaceSection section{built.mesh, geometry, regions, x, .2};
        const auto value{section.sample(built.mesh, field, gas, conditions)};
        require(section.face_count() == 3, "Complete radial face coverage.");
        require_near(value.area, exact_area, 2e-14, "Analytical disk section, not meridional length.");
        require_near(value.mach, Model::mach, 2e-14, "Section Mach.");
        require_near(value.axial_mach, Model::mach, 2e-14, "Axial section Mach.");
        require_near(value.pressure, Model::pressure, 2e-14, "Section pressure.");
        require_near(value.temperature, Model::temperature, 2e-14, "Section temperature.");
        require_near(value.mass_flow, Model::hot().density * Model::hot().velocity_x * exact_area, 2e-14,
                     "Independent uniform mass flow, oriented +x including the inlet.");
    }
    const cfd::verification::RpaFaceSection exit{built.mesh, geometry, regions, .75, .2};
    const auto other{cfd::verification::build_boundary_wave_mesh(4, 3, 1, .2)};
    require_throws<std::invalid_argument>(
        [&] { (void)cfd::verification::RpaFaceSection{other.mesh, geometry, regions, .75, .2}; },
        "Section geometry must belong to the exact same Mesh.");
    require_throws<std::invalid_argument>([&] { (void)exit.sample(other.mesh, field, gas, conditions); },
                                          "Section cannot sample another Mesh.");
    const cfd::ConservativeField short_field{1};
    require_throws<std::invalid_argument>([&] { (void)exit.sample(built.mesh, short_field, gas, conditions); },
                                          "Section field cardinality.");
    for (cfd::Index cell = 0; cell < field.size(); ++cell)
        if (regions[cell] != 1)
            field.set_state(cell, gas.to_conservative(Model::ambient(300, 0, 0)));
    require_near(exit.sample(built.mesh, field, gas, conditions).pressure, Model::pressure, 2e-14,
                 "Exit scalar estimate uses the Engine side, not ambient averaging.");
    require_throws<std::runtime_error>(
        [&] { (void)cfd::verification::RpaFaceSection{built.mesh, geometry, regions, .125, .2}; },
        "A non-face-aligned section cannot silently use nearest cells.");
    require_throws<std::runtime_error>(
        [&] { (void)cfd::verification::RpaFaceSection{built.mesh, geometry, regions, .5, .21}; },
        "A partial section must be rejected.");
}

void test_hybrid_euler_sections_source_cfl_boundaries_and_export()
{
    auto raw{cfd::verification::make_euler_tube(8, 3, cfd::verification::EulerTubeCells::Mixed)};
    raw.boundary_groups = {{0, "inlet"},           {1, "axis"},
                           {2, "inner_wall"},      {3, "outer_wall"},
                           {4, "nozzle_lip"},      {5, "farfield_upstream"},
                           {6, "farfield_radial"}, {7, "farfield_downstream"}};
    cfd::Index top{};
    for (auto &edge : raw.boundary_edges)
    {
        switch (edge.boundary_id)
        {
        case 0:
            edge.boundary_id = 1;
            break;
        case 1:
            edge.boundary_id = 2 + top++ % 5;
            break;
        case 2:
            edge.boundary_id = 0;
            break;
        case 3:
            edge.boundary_id = 7;
            break;
        default:
            throw std::runtime_error("Unexpected tube boundary.");
        }
    }
    for (auto &node : raw.nodes)
        node.y *= .2;
    auto invalid{raw};
    invalid.cell_nodes[1] = invalid.cell_nodes[0];
    require_throws<std::runtime_error>([&] { (void)cfd::build_mesh(std::move(invalid)); },
                                       "Degenerate hybrid cell accepted.");
    const auto built{cfd::build_mesh(std::move(raw))};
    const auto &mesh{built.mesh};
    cfd::verification::require_connected_rpa(mesh);
    const cfd::AxisymmetricGeometry geometry{mesh};
    cfd::Index mixed_faces{};
    for (cfd::Index face = 0; face < mesh.face_count(); ++face)
    {
        const auto adj{mesh.face_adjacencies()[face]};
        if (!adj.is_boundary() && mesh.cell_types()[adj.owner] != mesh.cell_types()[adj.neighbor])
        {
            ++mixed_faces;
            const auto sf{mesh.face_area_vectors()[face]};
            require(sf.x * (mesh.cell_centers()[adj.neighbor].x - mesh.cell_centers()[adj.owner].x) +
                            sf.y * (mesh.cell_centers()[adj.neighbor].y - mesh.cell_centers()[adj.owner].y) >
                        0,
                    "Mixed interface orientation.");
        }
    }
    require(mixed_faces > 0, "Fixture must exercise TRI/QUAD interfaces.");
    const cfd::CaloricallyPerfectGas gas{Model::gamma, Model::gas_constant};
    for (const bool moving : {false, true})
    {
        const cfd::PrimitiveState primitive{moving ? Model::hot() : cfd::PrimitiveState{1, 0, 0, 1}};
        const auto state{gas.to_conservative(primitive)};
        const cfd::ConservativeField field{mesh.cell_count(), state};
        const auto mapped{cfd::verification::rpa_euler_conditions(mesh.boundary_groups(), primitive)};
        // Synthetic closed inlet at rest: no total-condition reservoir is required.
        std::vector<cfd::EulerBoundaryCondition> data;
        data.reserve(8);
        for (const auto &group : mesh.boundary_groups())
            data.push_back(!moving && group.id == 0
                               ? cfd::EulerBoundaryCondition{cfd::EulerBoundaryConditionType::SlipWall}
                               : mapped.condition(group.id));
        const cfd::EulerBoundaryConditions conditions{8, std::move(data)};
        cfd::AxisymmetricEulerSolver solver{mesh, field, gas, conditions};
        const auto result{solver.step(1)};
        require(result.accepted_steps == 1 && result.rejected_steps == 0 && result.minimum_density > 0 &&
                    result.minimum_pressure > 0,
                "Hybrid Euler accepted-state invariants.");
        require_near(result.radial_pressure_source_rate, 2 * std::numbers::pi * .2 * primitive.pressure, 2e-12,
                     "Analytical integrated radial pressure source on mixed cells.");
        for (cfd::Index cell = 0; cell < mesh.cell_count(); ++cell)
        {
            const auto actual{solver.state().state(cell)};
            require_near(actual.density, state.density, 2e-12, "Hybrid uniform density.");
            require_near(actual.momentum_x, state.momentum_x, 2e-12, "Hybrid uniform axial momentum.");
            require_near(actual.momentum_r, 0, 2e-12, "Hybrid pressure/source balance at rest/tangency.");
            require_near(actual.total_energy_density, state.total_energy_density, 2e-12, "Hybrid uniform energy.");
        }
        for (const double defect : {result.balance_error.mass, result.balance_error.momentum_x,
                                    result.balance_error.momentum_r, result.balance_error.energy})
            require(std::isfinite(defect), "Nonfinite hybrid global conservation diagnostic.");
        if (!moving)
        {
            double expected{std::numeric_limits<double>::infinity()};
            for (cfd::Index cell = 0; cell < mesh.cell_count(); ++cell)
            {
                double area{};
                for (auto j = mesh.cell_node_offsets()[cell]; j < mesh.cell_node_offsets()[cell + 1]; ++j)
                {
                    const auto sf{geometry.face_area_vectors()[mesh.cell_faces()[j]]};
                    area += std::hypot(sf.x, sf.y);
                }
                expected = std::min(expected, .2 * geometry.cell_volumes()[cell] / (std::sqrt(Model::gamma) * area));
            }
            require_near(result.time_step, expected, 2e-14, "Independent mixed-cell CFL at rest.");
        }
        std::vector<double> regions(mesh.cell_count(), 1);
        for (const double x : {0.0, .5, .75})
        {
            if (!moving && x == 0)
                continue; // The synthetic closed inlet is not an open RPA section.
            const cfd::verification::RpaFaceSection section{mesh, geometry, regions, x, .2};
            const auto sample{section.sample(mesh, field, gas, conditions)};
            require(section.face_count() == 3, "Mixed section coverage.");
            require_near(sample.area, std::numbers::pi * .04, 2e-14, "Mixed section analytical disk area.");
            require_near(sample.mass_flow, primitive.density * primitive.velocity_x * sample.area, 2e-12,
                         "Mixed section production numerical flux.");
        }
    }
    const auto directory{
        std::filesystem::temp_directory_path() /
        ("cfd-rpa-hybrid-vtu-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()))};
    std::filesystem::create_directory(directory);
    cfd::write_vtu(mesh, directory / "mixed.vtu");
    const auto xml{cfd::test::read_text_file(directory / "mixed.vtu")};
    const auto array = [&](const std::string &name) {
        const auto begin{xml.find('>', xml.find("Name=\"" + name + "\"")) + 1};
        std::istringstream stream{xml.substr(begin, xml.find("</DataArray>", begin) - begin)};
        std::vector<cfd::Index> values;
        cfd::Index value{};
        while (stream >> value)
            values.push_back(value);
        return values;
    };
    const auto offsets{array("offsets")}, types{array("types")};
    require(std::ranges::equal(offsets, mesh.cell_node_offsets().subspan(1)), "VTU mixed CSR offsets.");
    require(std::find(types.begin(), types.end(), 5) != types.end() &&
                std::find(types.begin(), types.end(), 9) != types.end(),
            "VTU lost a mixed cell type.");
    std::filesystem::remove_all(directory); // Only this test's unique output.
}
} // namespace

int main()
{
    int failures{};
    failures += cfd::test::run_test("nominal and exploratory gas", test_nominal_derived_and_exploratory_states);
    failures += cfd::test::run_test("mixed Euler geometry/CFL/BC/sections/VTU",
                                    test_hybrid_euler_sections_source_cfl_boundaries_and_export);
    failures += cfd::test::run_test("eight explicit BC groups", test_eight_boundary_groups_and_no_silent_mapping);
    failures += cfd::test::run_test("regional initial state and first-step conservation",
                                    test_region_initialization_and_conservative_first_step);
    failures += cfd::test::run_test("private/old output protection", test_output_protection_including_symlink);
    failures += cfd::test::run_test("acoustic scale, physical exports and limits",
                                    test_acoustic_scale_physical_outputs_and_limits);
    failures += cfd::test::run_test("face-aligned section geometry and conservative flow",
                                    test_face_aligned_section_geometry_and_numerical_flow);
    return cfd::test::finish_tests(failures, "RPA exploratory");
}
