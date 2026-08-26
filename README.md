# CFD Solver

A 2D finite-volume CFD solver written in modern C++.

This project is developed as a learning project in both CFD and scientific C++, with an emphasis on mathematical correctness, robust validation, simple architecture, contiguous data structures, automated testing, and performance-aware design.

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
      ├── MeshStatistics / MeshReport
      └── VTU export
              ↓
           ParaView
```

## Implemented preprocessing features

- 2D rectangular domains
- Gmsh C++ API integration
- unstructured triangular meshes
- recombined quadrilateral meshes
- compact zero-based internal indexing
- flattened cell connectivity
- unique face construction
- owner/neighbor face adjacency
- physical boundary groups
- cell-to-face connectivity
- cell areas and centroids
- face centers and lengths
- oriented face-area vectors
- triangle and quadrilateral cell-quality metrics
- raw-mesh validation
- topology validation
- geometry validation
- per-cell face-area-vector closure validation
- preprocessing timings and mesh statistics
- VTU export for ParaView

The final `Mesh` owns contiguous arrays and exposes read-only views through `std::span`. It is move-only to avoid accidental copies of large mesh storage.

## Cell types

The preprocessing pipeline currently supports:

```text
Triangle
Quadrilateral
```

For quadrilateral generation, the Gmsh surface is recombined. The imported mesh is checked so that the generated cells match the requested cell type.

Concave or self-intersecting quadrilaterals are currently rejected.

## Geometry conventions

For a face shared by an owner cell `P` and a neighbor cell `N`, the stored face-area vector is oriented outward from the owner:

```text
Sf : owner → neighbor
```

For the neighbor cell, the corresponding outward vector is therefore `-Sf`.

For every closed cell, preprocessing verifies the finite-volume closure relation

```text
sum(Sf_cell) ≈ 0
```

using the appropriate sign for owner and neighbor incidences.

## Mesh quality

Triangle quality is based on

```text
q = 4 sqrt(3) A / (l1² + l2² + l3²)
```

with `q = 1` for an equilateral triangle.

Quadrilateral quality uses a local corner metric and keeps the worst corner. A square has quality `1`; skewness, poor angles, and large aspect ratios reduce the value.

The public mesh representation exposes the result generically as:

```text
cell_quality
```

## Physical boundaries

The reference rectangle uses the physical groups:

```text
inlet
wall
outlet
```

Automated tests verify that the groups survive the complete

```text
Gmsh → RawMeshData → Mesh
```

pipeline and that their accumulated boundary lengths match the analytical rectangle geometry.

## Reference case

The demonstration case uses a rectangular domain:

```text
Length    = 5.0 m
Height    = 1.0 m
Mesh size = 0.2 m
```

The expected total cell area is:

```text
5.0 m²
```

Exact node and cell counts are intentionally not treated as regression values because they may legitimately vary with the Gmsh version or meshing algorithm.

## Visualization

The solver writes:

```text
results/mesh.vtu
```

Open it with ParaView:

```bash
paraview results/mesh.vtu
```

Useful views include:

- `Surface With Edges`
- `cell_area`
- `cell_quality`
- `cell_id`

The VTU writer supports both VTK triangle type `5` and VTK quadrilateral type `9`.

ParaView is not linked to the solver. The numerical code only writes a standard VTU file.

## Project structure

```text
CFD_solver/
├── include/cfd/
│   ├── io/
│   │   ├── MeshReport.hpp
│   │   └── VtkWriter.hpp
│   ├── mesh/
│   └── meshing/
├── src/
│   ├── io/
│   ├── mesh/
│   ├── mesh_build/
│   ├── meshing/
│   ├── MeshBuilder.cpp
│   └── main.cpp
├── tests/
│   ├── GmshMesherTests.cpp
│   ├── MeshTopologyTests.cpp
│   ├── MeshGeometryTests.cpp
│   ├── MeshValidationTests.cpp
│   ├── MeshReportTests.cpp
│   ├── VtkWriterTests.cpp
│   └── support/
│       ├── MeshFixtures.hpp
│       └── TestUtils.hpp
├── .github/workflows/ci.yml
├── .clang-format
├── .clang-tidy
├── CMakeLists.txt
└── README.md
```

Main CMake targets:

- `cfd_core` — reusable preprocessing library
- `CFD_solver` — main executable
- dedicated CTest executables for meshing, topology, geometry, validation, reporting, and VTU output

## Requirements

- C++20 compiler
- CMake >= 3.20
- Gmsh development library

Ubuntu:

```bash
sudo apt install gmsh libgmsh-dev
```

Optional development and visualization tools:

```bash
sudo apt install clang-format clang-tidy paraview
```

## Build

Debug build:

```bash
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Debug \
    -DBUILD_TESTING=ON

cmake --build build
```

Run:

```bash
./build/CFD_solver
```

Run all tests:

```bash
ctest --test-dir build --output-on-failure
```

Release build:

```bash
cmake -S . -B build-release \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_TESTING=ON

cmake --build build-release --parallel
ctest --test-dir build-release --output-on-failure
```

## Testing

The test suite covers several levels of the preprocessing pipeline.

### Gmsh meshing

- invalid and non-finite domain dimensions
- invalid mesh sizes
- triangle generation
- quadrilateral generation
- physical boundary groups
- complete quadrilateral Gmsh-to-Mesh preprocessing

### Topology

- deterministic two-triangle connectivity
- unique shared face
- owner/neighbor assignment
- boundary-face assignment
- non-manifold rejection

### Geometry

- analytical triangle area and centroid
- translated small-cell numerical robustness
- analytical quadrilateral geometry
- clockwise quadrilateral connectivity
- face orientation
- face-area-vector norms
- cell-quality bounds
- per-cell face-area-vector closure
- total cell area versus independent boundary Shoelace area

### Raw-mesh and failure validation

- empty mandatory arrays
- non-finite node coordinates
- malformed offsets
- incorrect cell cardinality
- unsupported cell types
- invalid and duplicated node indices
- malformed boundary groups
- malformed boundary edges
- unused boundary groups
- open physical boundaries
- zero-area cells
- concave quadrilaterals

### Output and reporting

- mesh statistics
- report contents and stream-state preservation
- VTU triangle export
- VTU quadrilateral export

## Code quality

The project uses:

```text
-Wall
-Wextra
-Wpedantic
clang-format
clang-tidy
CTest
GitHub Actions
AddressSanitizer
UndefinedBehaviorSanitizer
```

Formatting is based on the Microsoft `clang-format` preset.

Static analysis uses selected checks from:

```text
clang-analyzer-*
bugprone-*
cppcoreguidelines-*
performance-*
```

CI builds and tests the project in both Debug and Release configurations, runs sanitizers, checks formatting, and runs `clang-tidy`.

Leak detection is disabled in the sanitizer CI job because Gmsh is an external library; AddressSanitizer still checks invalid memory accesses in the process.

## Design principles

- mathematical correctness before optimization
- software correctness before optimization
- explicit ownership
- contiguous storage for frequently traversed data
- no unnecessary copies of large mesh arrays
- no dynamic allocation inside future numerical loops
- validate imported data early
- separate construction, validation, statistics, reporting, and visualization
- keep Gmsh outside the future numerical hot loops
- test both valid and invalid inputs
- measure performance before optimizing

## Preprocessing v1 status

The preprocessing architecture is considered complete enough to serve as the stable base of the finite-volume solver.

Future geometric quantities should only be added when required by the numerical schemes. Likely examples include:

```text
owner-neighbor vectors
face interpolation weights
non-orthogonality measures
skewness corrections
```

These should be introduced when their mathematical use is defined, rather than precomputed speculatively.

## Next development stage

The next stage is the numerical-field infrastructure:

1. cell-centered scalar fields
2. boundary conditions
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

The solver currently targets 2D finite-volume CFD with triangular and quadrilateral cells.

General 3D polyhedral meshes are outside the current scope.
