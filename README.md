# CFD Solver

A 2D finite-volume CFD solver written in modern C++.

The project is developed primarily as a scientific-computing and C++ learning
project, with particular attention to:

- mathematical correctness;
- clear and maintainable software architecture;
- explicit data ownership;
- robust mesh validation;
- memory-efficient data structures;
- cache-friendly numerical loops;
- performance measurement before optimization.

The long-term objective is to solve steady incompressible 2D flows using the
finite-volume method, including pressure-velocity coupling with SIMPLE.

The numerical core is implemented in C++.
Gmsh is currently used only for mesh generation.

---

## Current status

The mesh preprocessing pipeline currently implements:

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
Mesh
      ↓
Discrete mesh geometry        [next step]
      ↓
Finite-volume discretization
      ↓
Flow solver
```

No CFD equations are solved yet.

The current development stage focuses on constructing a reliable and efficient
internal representation of an unstructured 2D mesh before implementing the
numerical solver.

---

## Implemented features

### Mesh generation

- 2D rectangular domain.
- Gmsh C++ API integration.
- Triangular mesh generation.
- Infrastructure prepared for quadrilateral cells.
- Physical boundary groups:
  - inlet;
  - wall;
  - outlet.

### Raw mesh import

Gmsh data is converted to an internal `RawMeshData` representation containing:

- node coordinates;
- cell types;
- compressed cell-to-node connectivity;
- physical boundary groups;
- boundary edges.

Gmsh node tags are converted to compact zero-based internal indices.

### Raw mesh validation

The imported mesh is checked before any derived connectivity is constructed.

Current checks include:

- non-empty required data;
- finite node coordinates;
- consistent compressed connectivity offsets;
- valid node indices;
- correct node count for each supported cell type;
- no repeated node inside one cell;
- valid and unique boundary groups;
- valid boundary-edge node indices;
- no degenerate boundary edge;
- valid boundary IDs;
- no duplicated boundary edge.

### Mesh topology

The internal topology builder constructs:

- unique faces;
- cell-to-face connectivity;
- face-to-cell connectivity;
- owner and neighbor cells;
- physical boundary assignment for external faces.

Faces are identified using canonical undirected edge keys during mesh
construction.

The temporary lookup structure is discarded once the final mesh topology has
been created.

### Topology validation

The constructed topology is validated before the final `Mesh` receives the
data.

Current checks include:

- consistent topology-array sizes;
- valid cell-to-face references;
- reciprocal cell-to-face / face-to-cell connectivity;
- correct nodes for every local cell face;
- valid owner and neighbor cells;
- at most two cells per face;
- external faces have exactly one physical boundary assignment;
- internal faces have no physical boundary assignment;
- every imported boundary edge corresponds to an external mesh face;
- no pair of cells shares more than one face;
- consistency between local, internal and boundary face counts.

A `Mesh` object is therefore created only after both raw data and constructed
topology have passed validation.

---

## Reference case

The current reference geometry is a rectangle:

```text
Length      = 5
Height      = 1
Mesh size   = 0.2
Cell type   = triangle
```

Typical generated mesh:

```text
183 nodes
304 cells
486 faces
├── 426 internal faces
└── 60 boundary faces
```

The topology satisfies:

```text
3 × 304 = 2 × 426 + 60
```

and, for this connected domain without holes, Euler's relation:

```text
183 - 486 + 304 = 1
```

Euler's relation is used as a reference-case consistency check rather than as
a universal mesh invariant.

---

## Project architecture

The code is separated according to responsibility.

```text
CFD_solver/
│
├── include/
│   └── cfd/
│       │
│       ├── mesh/
│       │   ├── Boundary.hpp
│       │   ├── Cell.hpp
│       │   ├── Face.hpp
│       │   ├── Mesh.hpp
│       │   ├── MeshBuilder.hpp
│       │   ├── Node.hpp
│       │   └── Types.hpp
│       │
│       └── meshing/
│           ├── DomainGeometry.hpp
│           ├── GmshMesher.hpp
│           ├── RawMeshData.hpp
│           └── RawMeshValidation.hpp
│
├── src/
│   ├── meshing/
│   │   ├── GmshMesher.cpp
│   │   └── RawMeshValidation.cpp
│   │
│   ├── mesh_build/
│   │   ├── BuildData.hpp
│   │   ├── TopologyBuilder.hpp
│   │   ├── TopologyBuilder.cpp
│   │   ├── TopologyValidation.hpp
│   │   └── TopologyValidation.cpp
│   │
│   ├── MeshBuilder.cpp
│   └── main.cpp
│
├── CMakeLists.txt
└── README.md
```

### `mesh`

Contains the persistent mesh representation used by the future CFD solver.

These data remain alive during the numerical simulation.

### `meshing`

Contains mesh-generation and raw-import infrastructure.

This layer produces `RawMeshData` but is not used inside CFD iteration loops.

Gmsh belongs exclusively to this part of the code.

### `mesh_build`

Contains internal implementation details used to transform validated raw mesh
data into the final `Mesh`.

These types and algorithms are not part of the public API.

---

## Data representation

The project avoids one dynamically allocated object per cell or face.

Cell connectivity uses contiguous flattened arrays:

```text
cell_nodes
cell_node_offsets
```

For example:

```text
cell_nodes:
[n0 n1 n2] [n3 n4 n5] [n6 n7 n8 n9]

cell_node_offsets:
0           3           6               10
```

The same offsets are reused for `cell_faces`, since a 2D polygon has the same
number of nodes and faces.

This representation avoids per-cell dynamic allocations and is intended to
remain suitable for large meshes.

The current internal index type is:

```cpp
using Index = std::size_t;
```

This choice may be revisited later if memory measurements justify using a
smaller integer representation.

---

## Design principles

The project follows several rules intended for scientific C++ development:

1. Correctness before optimization.
2. Validate imported and derived mesh data early.
3. Keep Gmsh outside the numerical core.
4. Prefer contiguous storage for frequently traversed data.
5. Avoid dynamic allocation inside numerical loops.
6. Avoid unnecessary copies of large arrays.
7. Use move semantics when transferring ownership of mesh data.
8. Keep temporary construction data out of the final `Mesh`.
9. Measure performance before changing algorithms for speed.
10. Keep abstractions simple unless they provide a concrete benefit.

The project intentionally avoids unnecessary inheritance, dynamic
polymorphism and shared ownership in the numerical mesh representation.

---

## Requirements

Current development environment:

- C++20 compiler;
- CMake >= 3.20;
- Gmsh development library.

On Ubuntu, Gmsh and its development files can be installed with:

```bash
sudo apt install gmsh libgmsh-dev
```

---

## Build

From the repository root:

```bash
cmake -S . -B build
cmake --build build
```

Run:

```bash
./build/CFD_solver
```

A typical run currently reports the Gmsh meshing process followed by a concise
summary of the constructed CFD topology:

```text
[CFD] Topology: 486 faces (426 internal, 60 boundary) [...]
Number of nodes: 183
Number of cells: 304
Number of faces: 486
```

Timing values are diagnostic only at this stage and should not be interpreted
as formal benchmarks unless the build configuration and benchmark conditions
are explicitly controlled.

---

## Next development step

The next stage is the construction and validation of the discrete mesh
geometry.

Planned persistent quantities include:

### Per cell

- area;
- centroid.

### Per face

- center;
- length;
- oriented area vector.

The geometry stage will also detect:

- degenerate cells;
- zero or near-zero face lengths;
- invalid orientations;
- inconsistent owner/neighbor geometry.

For the rectangular reference case, the computed cell areas will additionally
be checked against:

```text
sum(cell areas) ≈ 5
```

as an integration test.

---

## Planned solver development

After mesh geometry:

1. automated mesh and geometry tests;
2. boundary-condition representation;
3. scalar fields;
4. gradient reconstruction;
5. finite-volume diffusion and convection operators;
6. sparse matrix assembly;
7. linear-system solution;
8. momentum equations;
9. pressure correction;
10. SIMPLE pressure-velocity coupling;
11. Rhie-Chow interpolation if required;
12. convergence monitoring;
13. Poiseuille-flow validation;
14. backward-facing-step validation;
15. performance profiling and targeted optimization.

---

## Scope

The current solver is intentionally limited to 2D.

The architecture supports triangular and quadrilateral cells, but the project
does not currently target general 3D polyhedral meshes.

This constraint is deliberate: the goal is to build a correct and
well-structured 2D finite-volume solver before adding unnecessary
generality.
