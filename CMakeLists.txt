cmake_minimum_required(VERSION 3.20)

project(cfd_solver
    VERSION 0.1.0
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(CFD_solver
    src/main.cpp
    src/Node.cpp
    src/Face.cpp
    src/Cell.cpp
    src/Mesh.cpp
)

target_include_directories(CFD_solver
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(CFD_solver PRIVATE
        -Wall
        -Wextra
        -Wpedantic
    )
elseif(MSVC)
    target_compile_options(CFD_solver PRIVATE /W4)
endif()
