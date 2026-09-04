# CFD Solver

A 2D finite-volume CFD solver written in modern C++.

This project is developed as a learning project in both computational fluid dynamics and scientific C++, with an emphasis on mathematical correctness, robust validation, maintainable architecture, contiguous data structures, automated testing, and performance-aware design.

> **Current status:** mesh preprocessing and VTU export are implemented and validated. The finite-volume equations are not implemented yet.

## Current pipeline

```text
RectangleGeometry
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

The current preprocessing pipeline provides:

* 2D rectangular domains
* Gmsh C++ API integration
* unstructured triangular meshes
* recombined quadrilateral meshes
* compact zero-based internal IDs
* flattened cell-to-node connectivity
* unique face construction
* owner/neighbor face adjacency
* physical boundary groups
* cell-to-face connectivity
* cell areas and centroids
* face centers and lengths
* oriented face-area vectors
* triangle and quadrilateral cell-quality metrics
* raw-mesh validation
* topology validation
* geometry validation
* per-cell face-area-vector closure validation
* preprocessing timings
* mesh statistics and reporting
* VTU export for ParaView

The final `Mesh` owns contiguous arrays and exposes read-only views through `std::span`.

`Mesh` is move-only to avoid accidental copies of potentially large mesh storage.

## Cell types

The preprocessing pipeline currently supports:

```text
Triangle
Quadrilateral
```

For quadrilateral generation, the Gmsh surface is recombined.

The imported mesh is validated to ensure that generated cells match the requested cell type.

Concave or self-intersecting quadrilateral cells are currently rejected.

## Geometry conventions

For a face shared by an owner cell `P` and a neighbor cell `N`, the stored face-area vector is normal to the face, has magnitude equal to the face length, and is oriented outward from the owner:

```text
|Sf| = Lf

Sf · (x_N - x_P) > 0
```

On a non-orthogonal mesh, `Sf` is not generally parallel to the owner-neighbor vector `x_N - x_P`.

For the neighbor cell, the corresponding outward face-area vector is therefore:

```text
-Sf
```

For every closed cell, preprocessing verifies the finite-volume closure relation:

```text
Σ Sf ≈ 0
```

with the appropriate sign depending on whether the cell is the owner or the neighbor of each face.

## Mesh quality

### Triangles

Triangle quality is defined as:

```text
q = 4 sqrt(3) A / (l1² + l2² + l3²)
```

where:

* `A` is the triangle area
* `l1`, `l2`, `l3` are the edge lengths

An equilateral triangle has:

```text
q = 1
```

### Quadrilaterals

Quadrilateral quality is based on a local corner metric.

The quality of the cell is taken as the worst quality among its corners.

A square has quality:

```text
q = 1
```

Poor angles, skewness, and large aspect ratios reduce the metric.

The public mesh representation exposes both cell types through the same quantity:

```text
cell_quality
```

## Physical boundaries

The reference rectangular domain uses the physical groups:

```text
inlet
wall
outlet
```

`BoundaryId` values are compact zero-based indices into boundary-group storage.

External Gmsh physical-group tags are converted to these solver-internal IDs during mesh import.

Automated tests verify that the groups survive the complete pipeline:

```text
Gmsh
  ↓
RawMeshData
  ↓
Mesh
```

The accumulated lengths of the physical boundaries are also checked against the analytical dimensions of the rectangle.

## Reference case

The demonstration case uses:

```text
Length    = 5.0 m
Height    = 1.0 m
Mesh size = 0.2 m
```

The analytical domain area is:

```text
5.0 m²
```

The preprocessing pipeline verifies that the sum of the cell areas is consistent with the domain area.

Exact node and cell counts are deliberately not used as regression values because they may legitimately vary between Gmsh versions or meshing algorithms.

## Visualization

The solver writes:

```text
results/mesh.vtu
```

The file can be opened directly in ParaView:

```bash
paraview results/mesh.vtu
```

Useful ParaView views and fields include:

* `Surface With Edges`
* `cell_id`
* `cell_area`
* `cell_quality`

The VTU writer supports:

```text
VTK_TRIANGLE      = 5
VTK_QUAD          = 9
```

ParaView is not linked to the numerical core. The project only writes a standard VTU file.

## Project structure

```text
CFD_solver/
├── include/
│   └── cfd/
│       ├── io/
│       │   ├── MeshReport.hpp
│       │   └── VtkWriter.hpp
│       ├── mesh/
│       │   └── ...
│       └── meshing/
│           ├── GmshMesher.hpp
│           └── RectangleGeometry.hpp
│
├── src/
│   ├── io/
│   │   ├── MeshReport.cpp
│   │   └── VtkWriter.cpp
│   ├── mesh/
│   │   ├── MeshBuilder.cpp
│   │   └── MeshStatistics.cpp
│   ├── mesh_build/
│   │   ├── MeshBuildData.hpp
│   │   ├── GeometryBuilder.cpp
│   │   ├── GeometryBuilder.hpp
│   │   ├── GeometryValidation.cpp
│   │   ├── GeometryValidation.hpp
│   │   ├── TopologyBuilder.cpp
│   │   ├── TopologyBuilder.hpp
│   │   ├── TopologyValidation.cpp
│   │   └── TopologyValidation.hpp
│   ├── meshing/
│   │   ├── GmshMesher.cpp
│   │   └── RawMeshValidation.cpp
│   └── main.cpp
│
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
│
├── .github/
│   └── workflows/
│       └── ci.yml
│
├── .clang-format
├── .clang-tidy
├── CMakeLists.txt
└── README.md
```

The main CMake targets are:

* `cfd_core` — reusable preprocessing library
* `CFD_solver` — main executable
* dedicated CTest executables for meshing, topology, geometry, validation, reporting, and VTU output

## Requirements

Required:

* C++20 compiler
* CMake >= 3.20
* Eigen 3.4 development headers
* Gmsh development library

On Ubuntu:

```bash
sudo apt install gmsh libeigen3-dev libgmsh-dev
```

Optional development and visualization tools:

```bash
sudo apt install clang-format clang-tidy paraview
```

## Build

### Debug

Configure:

```bash
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Debug \
    -DBUILD_TESTING=ON
```

Build:

```bash
cmake --build build --parallel
```

Run:

```bash
./build/CFD_solver
```

Run the complete test suite:

```bash
ctest --test-dir build --output-on-failure
```

### Release

Configure:

```bash
cmake -S . -B build-release \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_TESTING=ON
```

Build and test:

```bash
cmake --build build-release --parallel

ctest --test-dir build-release --output-on-failure
```

Run:

```bash
./build-release/CFD_solver
```

## Testing

The test suite covers several levels of the preprocessing pipeline.

### Gmsh meshing

Tests cover:

* invalid and non-finite rectangle dimensions
* invalid mesh sizes
* triangular mesh generation
* quadrilateral mesh generation
* physical boundary groups
* complete Gmsh-to-Mesh preprocessing

### Topology

Tests cover:

* deterministic two-triangle connectivity
* unique shared-face construction
* owner/neighbor assignment
* boundary-face assignment
* cell-to-face connectivity
* non-manifold rejection

### Geometry

Tests cover:

* analytical triangle area and centroid
* small-cell numerical robustness under large coordinate translations
* analytical quadrilateral geometry
* analytical non-square quadrilateral quality (`q = 0.8` for a `2 x 1` rectangle)
* clockwise quadrilateral connectivity
* face orientation
* face-area-vector norms
* triangle and quadrilateral cell quality
* per-cell face-area-vector closure
* total cell area versus independently computed boundary area

### Raw-mesh and failure validation

Tests cover:

* empty mandatory arrays
* non-finite node coordinates
* malformed offsets
* incorrect cell cardinality
* unsupported cell types
* invalid node IDs
* duplicated node IDs inside cells
* malformed boundary groups
* non-compact boundary-group IDs
* malformed boundary edges
* unused boundary groups
* open physical boundaries
* zero-area cells
* concave quadrilaterals

### Output and reporting

Tests cover:

* mesh statistics
* mesh report contents
* stream-state preservation
* VTU triangle export
* VTU quadrilateral export
* VTU connectivity and offset values
* replacement of existing VTU output files

## Code quality

The project currently uses:

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

The CI pipeline:

* builds and tests Debug
* builds and tests Release
* checks formatting
* runs `clang-tidy`
* runs AddressSanitizer
* runs UndefinedBehaviorSanitizer

Leak detection is disabled in the sanitizer CI job because the external Gmsh library reports allocations that remain live at process termination.

AddressSanitizer remains active for invalid memory accesses, and UndefinedBehaviorSanitizer remains active for undefined behavior.

## Naming conventions

The preprocessing code follows these naming conventions:

```text
Types / classes / structs     PascalCase

Functions                     snake_case
Variables                     snake_case
Private members               snake_case_

Internal mesh IDs             *_id
Gmsh tags                     *_tag
Local indices                 *_index
Flattened-array positions     *_position
Entity counts                 *_count
Offsets                       *_offset / *_offsets
Booleans                      is_* / has_*
Containers                    plural names
```

Domain-specific terms such as `owner` and `neighbor` are retained where they are standard finite-volume terminology.

## Design principles

The project follows several explicit design principles:

* mathematical correctness before optimization
* software correctness before optimization
* simple and explicit ownership
* contiguous storage for frequently traversed numerical data
* no unnecessary copies of large mesh arrays
* no dynamic allocation in future numerical hot loops unless justified
* validate imported data early
* separate construction from validation
* separate mesh data from statistics and reporting
* keep Gmsh outside future numerical hot loops
* test both valid and invalid inputs
* measure performance before optimizing

Performance-related decisions are intended to remain evidence-based rather than speculative.

## Preprocessing v1 status

The preprocessing architecture is considered sufficiently complete to serve as the stable foundation of the finite-volume solver.

Current preprocessing provides the topology and geometry required to begin implementing numerical fields and finite-volume operators.

Additional geometric quantities will be introduced only when required by a numerical scheme.

Likely future quantities include:

```text
owner-neighbor vectors
face interpolation weights
non-orthogonality measures
skewness corrections
```

These quantities should be added when their mathematical role is defined rather than precomputed speculatively.

## Next development stage

The next stage is the numerical-field infrastructure.

Planned progression:

1. cell-centered scalar fields
2. boundary conditions
3. vector fields
4. gradient reconstruction
5. diffusion operators
6. convection operators
7. sparse matrix assembly
8. linear-system solution
9. momentum equations
10. pressure correction
11. SIMPLE pressure-velocity coupling

Planned CFD validation cases include:

* 2D Poiseuille flow
* backward-facing step

## Scope

The solver currently targets 2D finite-volume CFD with triangular and quadrilateral cells.

General 3D polyhedral meshes are outside the current scope.
