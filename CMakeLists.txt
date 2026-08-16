cmake_minimum_required(VERSION 3.20)

project(
    CFD_solver
    VERSION 0.1.0
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


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


add_executable(
    CFD_solver

    src/main.cpp
    src/Node.cpp
    src/Face.cpp
    src/Cell.cpp
    src/Mesh.cpp
    src/GmshMesher.cpp
)


target_include_directories(
    CFD_solver
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)


target_link_libraries(
    CFD_solver
    PRIVATE
        Gmsh::Gmsh
)


if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(
        CFD_solver
        PRIVATE
            -Wall
            -Wextra
            -Wpedantic
    )
endif()