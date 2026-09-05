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
# Eigen
# -----------------------------------------------------------------------------

find_package(Eigen3 3.4 REQUIRED NO_MODULE)

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

    src/input/OpenFOAMCaseReader.cpp

    src/io/MeshReport.cpp
    src/io/VtkWriter.cpp

    src/linear_algebra/EigenBiCGSTABSolver.cpp
    src/linear_algebra/EigenConjugateGradientSolver.cpp
    src/linear_algebra/ScalarLinearSystem.cpp

    src/mesh/MeshBuilder.cpp
    src/mesh/MeshStatistics.cpp

    src/meshing/GmshMesher.cpp
    src/meshing/RawMeshValidation.cpp

    src/mesh_build/GeometryBuilder.cpp
    src/mesh_build/GeometryValidation.cpp
    src/mesh_build/TopologyBuilder.cpp
    src/mesh_build/TopologyValidation.cpp

    src/numerics/LeastSquaresGradient.cpp
    src/numerics/ScalarConvectionOperator.cpp
    src/numerics/ScalarDiffusionOperator.cpp
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
    PUBLIC
        Eigen3::Eigen

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
# Numerical verification
# -----------------------------------------------------------------------------

option(
    BUILD_VERIFICATION
    "Build numerical verification studies"
    OFF
)

if(BUILD_VERIFICATION)
    function(cfd_add_verification target_name source_file)
        add_executable(
            ${target_name}
            ${source_file}
        )

        target_link_libraries(
            ${target_name}
            PRIVATE
                cfd_core
        )

        cfd_enable_warnings(${target_name})
    endfunction()

    cfd_add_verification(
        cfd_least_squares_gradient_convergence
        verification/LeastSquaresGradientConvergence.cpp
    )

    cfd_add_verification(
        cfd_structured_least_squares_gradient_convergence
        verification/StructuredLeastSquaresGradientConvergence.cpp
    )

    cfd_add_verification(
        cfd_controlled_quad_transition_convergence
        verification/ControlledQuadTransitionConvergence.cpp
    )

    cfd_add_verification(
        cfd_scalar_diffusion_operator_convergence
        verification/ScalarDiffusionOperatorConvergence.cpp
    )

    cfd_add_verification(
        cfd_scalar_diffusion_solution_convergence
        verification/ScalarDiffusionSolutionConvergence.cpp
    )

    cfd_add_verification(
        cfd_unstructured_scalar_diffusion_convergence
        verification/UnstructuredScalarDiffusionConvergence.cpp
    )
endif()

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
        cfd_cell_vector_field_tests
        cell_vector_field
        tests/CellVectorFieldTests.cpp
    )

    cfd_add_test(
        cfd_face_flux_field_tests
        face_flux_field
        tests/FaceFluxFieldTests.cpp
    )

    cfd_add_test(
        cfd_scalar_boundary_conditions_tests
        scalar_boundary_conditions
        tests/ScalarBoundaryConditionsTests.cpp
    )

    cfd_add_test(
        cfd_openfoam_case_reader_tests
        openfoam_case_reader
        tests/OpenFOAMCaseReaderTests.cpp
    )

    cfd_add_test(
        cfd_least_squares_gradient_tests
        least_squares_gradient
        tests/LeastSquaresGradientTests.cpp
    )

    cfd_add_test(
        cfd_scalar_diffusion_operator_tests
        scalar_diffusion_operator
        tests/ScalarDiffusionOperatorTests.cpp
    )

    cfd_add_test(
        cfd_scalar_convection_operator_tests
        scalar_convection_operator
        tests/ScalarConvectionOperatorTests.cpp
    )

    cfd_add_test(
        cfd_scalar_linear_system_tests
        scalar_linear_system
        tests/ScalarLinearSystemTests.cpp
    )

    cfd_add_test(
        cfd_eigen_bicgstab_solver_tests
        eigen_bicgstab_solver
        tests/EigenBiCGSTABSolverTests.cpp
    )

    cfd_add_test(
        cfd_eigen_conjugate_gradient_solver_tests
        eigen_conjugate_gradient_solver
        tests/EigenConjugateGradientSolverTests.cpp
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
