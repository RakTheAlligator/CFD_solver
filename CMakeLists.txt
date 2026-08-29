cmake_minimum_required(VERSION 3.20)

project(
    CFD_solver
    VERSION 0.1.0
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# -----------------------------------------------------------------------------
# Gmsh
# -----------------------------------------------------------------------------

find_path(
    GMSH_INCLUDE_DIR
    NAMES gmsh.h
    REQUIRED
)

find_library(
    GMSH_LIBRARY
    NAMES gmsh
    REQUIRED
)

add_library(Gmsh::Gmsh UNKNOWN IMPORTED)

set_target_properties(
    Gmsh::Gmsh
    PROPERTIES
        IMPORTED_LOCATION "${GMSH_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${GMSH_INCLUDE_DIR}"
)

# -----------------------------------------------------------------------------
# Compiler warnings
# -----------------------------------------------------------------------------

function(cfd_enable_warnings target_name)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(
            ${target_name}
            PRIVATE
                -Wall
                -Wextra
                -Wpedantic
        )
    endif()
endfunction()

# -----------------------------------------------------------------------------
# CFD core library
# -----------------------------------------------------------------------------

add_library(
    cfd_core
    STATIC

    src/io/MeshReport.cpp
    src/io/VtkWriter.cpp

    src/mesh/MeshBuilder.cpp
    src/mesh/MeshStatistics.cpp

    src/meshing/GmshMesher.cpp
    src/meshing/RawMeshValidation.cpp

    src/mesh_build/GeometryBuilder.cpp
    src/mesh_build/GeometryValidation.cpp
    src/mesh_build/TopologyBuilder.cpp
    src/mesh_build/TopologyValidation.cpp
)

target_include_directories(
    cfd_core

    PUBLIC
        ${PROJECT_SOURCE_DIR}/include

    PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(
    cfd_core
    PRIVATE
        Gmsh::Gmsh
)

cfd_enable_warnings(cfd_core)

# -----------------------------------------------------------------------------
# Main executable
# -----------------------------------------------------------------------------

add_executable(
    CFD_solver
    src/main.cpp
)

target_link_libraries(
    CFD_solver
    PRIVATE
        cfd_core
)

cfd_enable_warnings(CFD_solver)

# -----------------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------------

include(CTest)

if(BUILD_TESTING)
    function(cfd_add_test target_name test_name source_file)
        add_executable(
            ${target_name}
            ${source_file}
        )

        target_include_directories(
            ${target_name}
            PRIVATE
                ${PROJECT_SOURCE_DIR}/tests
        )

        target_link_libraries(
            ${target_name}
            PRIVATE
                cfd_core
        )

        cfd_enable_warnings(${target_name})

        add_test(
            NAME ${test_name}
            COMMAND ${target_name}
        )
    endfunction()

    cfd_add_test(
        cfd_cell_scalar_field_tests
        cell_scalar_field
        tests/CellScalarFieldTests.cpp
    )

    cfd_add_test(
        cfd_gmsh_mesher_tests
        gmsh_mesher
        tests/GmshMesherTests.cpp
    )

    cfd_add_test(
        cfd_mesh_topology_tests
        mesh_topology
        tests/MeshTopologyTests.cpp
    )

    cfd_add_test(
        cfd_mesh_geometry_tests
        mesh_geometry
        tests/MeshGeometryTests.cpp
    )

    cfd_add_test(
        cfd_mesh_validation_tests
        mesh_validation
        tests/MeshValidationTests.cpp
    )

    cfd_add_test(
        cfd_mesh_report_tests
        mesh_report
        tests/MeshReportTests.cpp
    )

    cfd_add_test(
        cfd_vtk_writer_tests
        vtk_writer
        tests/VtkWriterTests.cpp
    )
endif()
