# CFD Solver

A 2D finite-volume CFD solver written in modern C++.

This project is developed as a learning project in both CFD and scientific C++, with an emphasis on correctness, simple architecture, robust validation, contiguous data structures, automated testing, and performance-aware design.

> **Current status:** mesh preprocessing and visualization are implemented. The CFD equations are not implemented yet.

## Current pipeline

```text
Domain geometry
      ↓
Gmsh mesh generation
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
      ↓
VTU export
      ↓
ParaView
```

## Implemented features

- 2D rectangular domains
- Gmsh C++ API integration
- triangular unstructured meshes
- compact zero-based internal indexing
- flattened cell connectivity
- face construction with owner/neighbor adjacency
- physical boundary groups
- cell areas and centroids
- face centers, lengths and oriented area vectors
- triangle-quality computation
- topology and geometry validation
- VTU export for ParaView
- automated positive and negative tests
- `clang-format` and `clang-tidy` integration

The current VTU export includes:

```text
cell_id
cell_area
cell_quality
```

which can be visualized directly in ParaView.

## Reference case

Current rectangular test case:

```text
Length    = 5.0 m
Height    = 1.0 m
Mesh size = 0.2 m
```

Typical mesh:

```text
183 nodes
304 triangular cells
486 faces
├── 426 internal faces
└── 60 boundary faces
```

Typical preprocessing diagnostics:

```text
Total area        = 5.0000 m²
Triangle quality  ≈ 0.892 ... 1.000
```

The total cell area is also checked against an independent Shoelace reconstruction of the external boundary.

## Visualization

The solver writes:

```text
results/mesh.vtu
```

Open it with ParaView:

```bash
paraview results/mesh.vtu
```

Useful views:

- `Surface With Edges` for the mesh
- `cell_area` for cell-size distribution
- `cell_quality` for mesh-quality visualization
- `cell_id` for identifying individual cells

The visualization layer is intentionally separate from the numerical core. The C++ solver only writes a standard VTU file; ParaView handles interactive rendering.

## Project structure

```text
CFD_solver/
├── include/cfd/
│   ├── io/
│   │   └── VtkWriter.hpp
│   ├── mesh/
│   └── meshing/
├── src/
│   ├── io/
│   │   └── VtkWriter.cpp
│   ├── mesh_build/
│   ├── meshing/
│   ├── MeshBuilder.cpp
│   └── main.cpp
├── tests/
│   └── MeshPreprocessingTests.cpp
├── .clang-format
├── .clang-tidy
├── CMakeLists.txt
└── README.md
```

Main CMake targets:

- `cfd_core` — reusable mesh/preprocessing library
- `CFD_solver` — main executable
- `cfd_tests` — automated tests

## Requirements

- C++20 compiler
- CMake >= 3.20
- Gmsh development library

Ubuntu:

```bash
sudo apt install gmsh libgmsh-dev
```

Optional development/visualization tools:

```bash
sudo apt install clang-format clang-tidy paraview
```

ParaView is not linked to the solver.

## Build

```bash
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build
```

Run:

```bash
./build/CFD_solver
```

Run tests:

```bash
ctest --test-dir build --output-on-failure
```

## Testing

Current tests include:

- analytical single-triangle geometry
- rectangular Gmsh mesh
- cell-area sum versus boundary Shoelace area
- face orientation consistency
- invalid connectivity
- duplicated cell nodes
- duplicated boundary edges
- open boundaries
- zero-area cells
- non-manifold faces
- invalid boundary IDs
- VTU export of a reference triangle

## Code quality

The project currently uses:

```text
-Wall
-Wextra
-Wpedantic
clang-format
clang-tidy
CTest
```

Formatting is based on the Microsoft `clang-format` preset.

Static analysis uses selected checks from:

```text
clang-analyzer-*
bugprone-*
cppcoreguidelines-*
performance-*
```

The project follows the C++ Core Guidelines where they provide practical value for scientific code, while avoiding unnecessary abstraction in numerical data structures.

## Design principles

- correctness before optimization
- explicit ownership
- contiguous storage for frequently traversed data
- no unnecessary copies of large mesh arrays
- no dynamic allocation inside future numerical loops
- validate data early
- keep Gmsh and visualization outside the numerical core
- test both valid and invalid inputs
- measure performance before optimizing

## Next development stage

The next step is the numerical-field infrastructure required by the finite-volume solver:

1. cell-centered scalar fields
2. scalar boundary conditions
3. vector fields
4. gradient reconstruction
5. diffusion and convection operators
6. sparse matrix assembly
7. linear-system solution
8. momentum equations
9. pressure correction
10. SIMPLE pressure-velocity coupling

Planned CFD validation cases:

- 2D Poiseuille flow
- backward-facing step

## Scope

The solver is currently limited to 2D.

Triangles are supported first. Quadrilateral support may be added later without changing the fundamental mesh representation.

The project does not currently target general 3D polyhedral meshes.
