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

    src/io/VtkWriter.cpp

    src/MeshBuilder.cpp

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
    add_executable(
        cfd_tests
        tests/MeshPreprocessingTests.cpp
    )

    target_link_libraries(
        cfd_tests
        PRIVATE
            cfd_core
    )

    cfd_enable_warnings(cfd_tests)

    add_test(
        NAME mesh_preprocessing
        COMMAND cfd_tests
    )
endif()