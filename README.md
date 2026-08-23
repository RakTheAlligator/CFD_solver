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
Discrete geometry construction
      ↓
Geometry validation
      ↓
Mesh
      ↓
Finite-volume discretization       [next numerical stage]
      ↓
Flow solver
```

No CFD equations are solved yet.

The current code builds a validated internal representation of an unstructured
2D mesh, including the topology and the geometric quantities required by a
finite-volume solver.

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

The constructed topology is validated before being transferred to the final
`Mesh`.

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

### Mesh geometry

The geometry builder computes the quantities required by the future
finite-volume discretization.

Per cell:

- area;
- centroid.

Per face:

- center;
- length;
- oriented area vector.

For a face with tangent vector

```text
(dx, dy)
```

the 2D area vector is constructed from

```text
(dy, -dx)
```

and oriented so that it points outward from the owner cell.

Its norm therefore satisfies:

```text
|Sf| = face length
```

The cell area and centroid are computed using polygon formulas, allowing the
same implementation to support triangles and simple quadrilaterals.

### Geometry validation and diagnostics

The computed geometry is validated before being transferred to the final
`Mesh`.

Current checks include:

- consistent geometry-array sizes;
- finite and positive cell areas;
- finite cell centers;
- finite and positive face lengths;
- finite face centers;
- finite face area vectors;
- consistency between face length and area-vector norm;
- outward orientation of the face area vector relative to the owner cell;
- consistent owner-to-neighbor orientation for internal faces.

The preprocessing stage also computes global geometry statistics:

- minimum, mean and maximum cell area;
- minimum, mean and maximum characteristic cell size;
- minimum, mean and maximum face length;
- total cell area;
- minimum, mean and maximum triangle quality;
- ID of the worst-quality triangle.

The current characteristic cell size is defined as:

```text
h = sqrt(cell area)
```

This is a convenient length scale and should not be interpreted as an edge
length or a complete cell-shape metric.

Triangle quality is currently defined as:

```text
q = 4 * sqrt(3) * A / (l1² + l2² + l3²)
```

with:

```text
q = 1       equilateral triangle
q -> 0      degenerate triangle
```

Mesh quality is currently reported as a diagnostic rather than used with an
arbitrary rejection threshold.

---

## Reference case

The current reference geometry is a rectangle:

```text
Length      = 5
Height      = 1
Mesh size   = 0.2
Cell type   = triangle
```

A typical generated mesh contains:

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

The geometry currently gives:

```text
Total cell area = 5

Cell area:
    min  = 0.011484
    mean = 0.0164474
    max  = 0.0209765

Characteristic cell size:
    min  = 0.107164
    mean = 0.12803
    max  = 0.144833

Face length:
    min  = 0.153073
    mean = 0.195962
    max  = 0.239533

Triangle quality:
    min  = 0.892375
    mean = 0.984439
    max  = 1
```

The total cell area is consistent with the analytical rectangle area:

```text
5 × 1 = 5
```

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
│       │   ├── Types.hpp
│       │   └── Vector2.hpp
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
│   │   ├── TopologyValidation.cpp
│   │   ├── GeometryBuilder.hpp
│   │   ├── GeometryBuilder.cpp
│   │   ├── GeometryValidation.hpp
│   │   └── GeometryValidation.cpp
│   │
│   ├── MeshBuilder.cpp
│   └── main.cpp
│
├── CMakeLists.txt
└── README.md
```

### `mesh`

Contains the persistent mesh representation used by the future CFD solver.

This includes:

- nodes;
- cell connectivity;
- face connectivity;
- boundary information;
- cell geometry;
- face geometry.

These data remain alive during the numerical simulation.

### `meshing`

Contains mesh-generation and raw-import infrastructure.

This layer produces `RawMeshData` but is not used inside CFD iteration loops.

Gmsh belongs exclusively to this part of the code.

### `mesh_build`

Contains internal implementation details used to transform validated raw mesh
data into the final `Mesh`.

This includes:

- topology construction;
- topology validation;
- geometry construction;
- geometry validation;
- temporary build data;
- preprocessing statistics.

These types and algorithms are implementation details and are not part of the
public API.

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

The persistent geometry is also stored in contiguous arrays:

```text
cell_areas
cell_centers

face_centers
face_lengths
face_area_vectors
```

This representation avoids per-cell and per-face dynamic allocations and is
intended to remain suitable for large meshes.

The current internal index type is:

```cpp
using Index = std::size_t;
```

This choice may be revisited later if memory measurements justify using a
smaller integer representation.

---

## Data ownership

The preprocessing pipeline uses explicit ownership transfer.

Conceptually:

```text
RawMeshData
    ↓
validated raw data

TopologyBuildData
    ↓
validated topology

GeometryBuildData
    ↓
validated geometry

Mesh
```

Large `std::vector` buffers are moved into the final `Mesh` instead of being
deep-copied.

Temporary construction data such as the topology hash table and preprocessing
statistics are destroyed before the CFD iterations begin.

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

Current development requirements:

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

A typical run currently reports the Gmsh meshing process followed by concise
topology and geometry diagnostics:

```text
[CFD] Topology: 486 faces (426 internal, 60 boundary) [...]
[CFD] Geometry: total area=5 [...]
      Cell area:   min=0.011484, mean=0.0164474, max=0.0209765
      Cell size:   min=0.107164, mean=0.12803, max=0.144833
      Face length: min=0.153073, mean=0.195962, max=0.239533
      Triangle q:  min=0.892375, mean=0.984439, max=1, worst cell=33
Number of nodes: 183
Number of cells: 304
Number of faces: 486
```

Timing values are diagnostic only at this stage and should not be interpreted
as formal benchmarks unless the build configuration and benchmark conditions
are explicitly controlled.

---

## Next development step

The next step is to introduce automated tests for the mesh preprocessing
pipeline before adding the finite-volume numerical operators.

The reference rectangular case will be used to verify automatically:

- node, cell and face counts;
- internal and boundary face counts;
- topology invariants;
- total domain area;
- geometric consistency;
- face orientation;
- triangle-quality bounds.

Once this preprocessing foundation is covered by automated tests, development
will continue with the numerical infrastructure required by the
finite-volume solver.

---

## Planned solver development

The current roadmap is:

1. automated mesh, topology and geometry tests;
2. boundary-condition representation;
3. scalar and vector fields;
4. gradient reconstruction;
5. finite-volume diffusion operator;
6. finite-volume convection operator;
7. sparse matrix assembly;
8. linear-system solution;
9. momentum equations;
10. pressure correction;
11. SIMPLE pressure-velocity coupling;
12. Rhie-Chow interpolation if required;
13. convergence monitoring;
14. Poiseuille-flow validation;
15. backward-facing-step validation;
16. performance profiling and targeted optimization.

Performance work will distinguish at least:

- mesh generation;
- topology construction;
- geometry construction;
- matrix assembly;
- gradient computation;
- linear solver;
- SIMPLE iterations;
- result export.

---

## Scope

The current solver is intentionally limited to 2D.

The architecture supports triangular cells and is designed to accommodate
quadrilateral cells without changing the fundamental mesh representation.

The project does not currently target general 3D polyhedral meshes.

This constraint is deliberate: the objective is to build a correct,
well-structured and measurable 2D finite-volume solver before adding
unnecessary generality.
