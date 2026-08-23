# CFD Solver

A 2D finite-volume CFD solver written in modern C++.

This project is developed as a learning project for both CFD and scientific C++.
The current focus is on building a clean, validated and efficient mesh infrastructure
before implementing the flow solver itself.

## Status

The CFD equations are **not implemented yet**.

The current preprocessing pipeline is:

```text
Gmsh
  ↓
RawMeshData
  ↓
Raw mesh validation
  ↓
Topology construction
  ↓
Topology validation
  ↓
Geometry construction
  ↓
Geometry validation
  ↓
Mesh
```

The next stage is the numerical infrastructure required by the finite-volume solver:
boundary conditions, scalar/vector fields and gradient reconstruction.

## Current features

- 2D rectangular domain generation with the Gmsh C++ API
- triangular unstructured meshes
- compact zero-based internal indexing
- flattened cell connectivity
- face construction with owner/neighbor adjacency
- physical boundary groups
- cell areas and centroids
- face centers, lengths and oriented area vectors
- mesh topology and geometry validation
- basic mesh-quality diagnostics
- automated positive and negative preprocessing tests

The current triangle quality metric is:

```text
q = 4 * sqrt(3) * A / (l1² + l2² + l3²)
```

with `q = 1` for an equilateral triangle.

## Project structure

```text
CFD_solver/
├── include/cfd/
│   ├── mesh/
│   └── meshing/
├── src/
│   ├── meshing/
│   ├── mesh_build/
│   ├── MeshBuilder.cpp
│   └── main.cpp
├── tests/
│   └── MeshPreprocessingTests.cpp
├── .clang-format
├── .clang-tidy
├── CMakeLists.txt
└── README.md
```

The main CMake targets are:

- `cfd_core`: reusable mesh/preprocessing library
- `CFD_solver`: main executable
- `cfd_tests`: automated preprocessing tests

## Requirements

- C++20 compiler
- CMake >= 3.20
- Gmsh development library

On Ubuntu:

```bash
sudo apt install gmsh libgmsh-dev
```

Optional development tools:

```bash
sudo apt install clang-format clang-tidy
```

## Build

```bash
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build
```

Run the current executable:

```bash
./build/CFD_solver
```

## Tests

Run the automated tests with:

```bash
ctest --test-dir build --output-on-failure
```

The current tests cover both valid and invalid meshes, including:

- analytical single-triangle geometry
- reference rectangular mesh
- total cell area versus boundary Shoelace area
- face orientation and geometry consistency
- invalid connectivity
- duplicated nodes and boundary edges
- open boundaries
- zero-area cells
- non-manifold faces
- invalid boundary IDs

## Code quality

The project currently uses:

- `-Wall -Wextra -Wpedantic`
- `clang-format` with the Microsoft preset
- `clang-tidy` with selected `clang-analyzer`, `bugprone`,
  `cppcoreguidelines` and `performance` checks

The project follows the C++ Core Guidelines where they are useful for scientific
code, without adding unnecessary abstraction to numerical data structures.

## Reference case

Current rectangular test case:

```text
Length    = 5
Height    = 1
Mesh size = 0.2
```

Typical Gmsh output:

```text
183 nodes
304 cells
486 faces
```

The computed total cell area and the independently reconstructed boundary area
both equal the analytical area:

```text
A = 5
```

## Roadmap

1. boundary-condition representation
2. scalar and vector fields
3. gradient reconstruction
4. diffusion and convection operators
5. sparse matrix assembly
6. linear solver
7. momentum equations
8. pressure correction
9. SIMPLE coupling
10. Poiseuille-flow validation
11. backward-facing-step validation
12. profiling and targeted optimization

## Scope

The solver is currently limited to 2D.

Triangles are supported first. Quadrilateral support may be added later without
changing the fundamental mesh representation.
