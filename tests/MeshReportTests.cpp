#include "cfd/io/MeshReport.hpp"
#include "cfd/mesh/MeshBuilder.hpp"
#include "cfd/mesh/MeshStatistics.hpp"
#include "cfd/meshing/RawMeshData.hpp"

#include "support/MeshFixtures.hpp"
#include "support/TestUtils.hpp"

#include <chrono>
#include <iomanip>
#include <sstream>
#include <string>
#include <utility>

namespace
{

using cfd::test::make_single_triangle_raw_mesh;
using cfd::test::require;
using cfd::test::require_contains;

void test_single_triangle_mesh_report()
{
    cfd::RawMeshData raw_mesh{make_single_triangle_raw_mesh()};

    cfd::MeshBuildResult build_result{cfd::build_mesh(std::move(raw_mesh))};

    const cfd::MeshStatistics statistics{cfd::compute_mesh_statistics(build_result.mesh)};

    const cfd::MeshBuildTimings timings{
        .topology = std::chrono::duration<double, std::milli>{1.25},
        .geometry = std::chrono::duration<double, std::milli>{2.50},
    };

    std::ostringstream output;

    output << std::scientific << std::setprecision(7);

    const auto original_flags{output.flags()};
    const auto original_precision{output.precision()};

    cfd::write_mesh_report(output, build_result.mesh, statistics, timings);

    require(output.flags() == original_flags, "Mesh report modified the output stream formatting flags.");

    require(output.precision() == original_precision, "Mesh report modified the output stream precision.");

    const std::string report{output.str()};

    require_contains(report, "[Mesh topology]", "Mesh report does not contain the topology section.");

    require_contains(report, "Faces             : 3", "Mesh report contains an incorrect face count.");

    require_contains(report, "Internal         : 0", "Mesh report contains an incorrect internal-face count.");

    require_contains(report, "Boundary         : 3", "Mesh report contains an incorrect boundary-face count.");

    require_contains(report, "1.25 ms", "Mesh report contains an incorrect topology timing.");

    require_contains(report, "[Mesh geometry]", "Mesh report does not contain the geometry section.");

    require_contains(report, "Total area        : 0.5000 m^2", "Mesh report contains an incorrect total area.");

    require_contains(report, "2.50 ms", "Mesh report contains an incorrect geometry timing.");

    require_contains(report, "Cell area", "Mesh report does not contain cell-area statistics.");

    require_contains(report, "Cell size", "Mesh report does not contain cell-size statistics.");

    require_contains(report, "Face length", "Mesh report does not contain face-length statistics.");

    require_contains(report, "Cell quality", "Mesh report does not contain cell-quality statistics.");

    require_contains(report, "Worst cell    : cell 0", "Mesh report contains an incorrect worst-quality cell.");
}

} // namespace

int main()
{
    int failure_count{};

    failure_count += cfd::test::run_test("single triangle mesh report", test_single_triangle_mesh_report);

    return cfd::test::finish_tests(failure_count, "mesh report");
}