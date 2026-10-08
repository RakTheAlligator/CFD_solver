# CFD Solver

CFD_solver is a modern C++20 finite-volume solver for steady, two-dimensional incompressible flow. The project emphasizes explicit ownership, contiguous numerical storage, input validation, numerical verification, and measured performance work.

The current solver supports triangular and quadrilateral cell-centered meshes and a case-based command-line workflow.

## Implemented capabilities

- Rectangle and backward-facing-step meshing through the Gmsh C++ API, including automatic grading/refinement controls
- Triangular and recombined quadrilateral cells
- Validated topology, owner/neighbor connectivity, boundary groups, cell geometry, and oriented face-area vectors
- Fixed-cardinality scalar, vector, velocity, face-flux, and pressure-response fields
- Inverse-distance-weighted least-squares gradient reconstruction
- First-order upwind, geometry-aware Linear, Hybrid, and LinearUpwind deferred convection correction
- Optional Barth-Jespersen limiting of LinearUpwind reconstructions
- Isotropic diffusion with non-orthogonal correction
- Momentum-weighted/Rhie-Chow face-flux interpolation
- Steady SIMPLE pressure-velocity coupling with momentum under-relaxation and Majumdar-consistent face interpolation
- Steady scalar transport with supplied face flux and constant diffusivity (numerical API and verification drivers)
- Scalar finite-volume sparse systems
- Eigen BiCGSTAB momentum solves and conjugate-gradient pressure-correction solves
- Reusable Eigen sparse patterns with numerical coefficient updates
- Independent-mesh cell-field transfer and multilevel coarse-grid initialization
- VTU cell-data export for ParaView
- Per-iteration convergence CSV output and optional live plotting
- Structured and Gmsh numerical-verification drivers

The detailed current SIMPLE contract is documented in [IncompressibleSimpleV1.md](docs/IncompressibleSimpleV1.md).
Development studies, measured decisions, and evidence limits are recorded in [DevelopmentHistoryAndEvidence.md](docs/DevelopmentHistoryAndEvidence.md).

## Requirements

- CMake 3.20 or newer
- C++20 compiler
- Eigen 3.4 development headers
- Gmsh development library

On Ubuntu:

```bash
sudo apt install gmsh libeigen3-dev libgmsh-dev
```

The optional live convergence window also requires Python 3 and Matplotlib:

```bash
python3 -m pip install matplotlib
```

## Build and run

Configure and build a Debug tree:

```bash
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Debug \
    -DBUILD_TESTING=ON

cmake --build build --parallel
```

Run the supplied pressure-driven Poiseuille case:

```bash
./build/CFD_solver cases/poiseuille
```

The executable accepts exactly one case-directory argument.

For optimized runs:

```bash
cmake -S . -B build-release \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_TESTING=ON

cmake --build build-release --parallel
./build-release/CFD_solver cases/poiseuille
```

## Case structure

The input syntax is a deliberately limited, OpenFOAM-inspired dictionary subset; it is not full OpenFOAM compatibility.

```text
cases/poiseuille/
├── 0/
│   ├── p
│   ├── u
│   └── v
└── system/
    ├── controlDict
    └── meshDict
```

### `system/meshDict`

`meshDict` defines rectangle or backward-facing-step geometry and Gmsh controls:

```text
geometry
{
    type rectangle;
    length 4.0;
    height 1.0;
}

mesh
{
    cellType quadrilateral;
    size 0.1;
}
```

`cellType` accepts the currently supported `triangle` and `quadrilateral` meshes. `size` is the requested Gmsh characteristic length.
The backward-facing-step example is in `cases/backwardFacingStep`; its mesh controls also support automatic meshing, growth limits, and wall/step refinement.

### Initial and boundary fields

The files `0/u`, `0/v`, and `0/p` provide:

- the field object name and dimensions;
- a uniform cell-centered initial value;
- one named condition for every mesh boundary group.

The supported scalar conditions are uniform `fixedValue`, `fixedGradient`, and `zeroGradient`. Velocity components are currently supplied as separate scalar files. Pressure `fixedValue` boundaries map to `FixedPressure`, while pressure `zeroGradient` and `fixedGradient` boundaries map to `FixedMassFlux` in SIMPLE. The current case runner requires both velocity components to be `fixedValue` on a fixed-mass-flux boundary so it can initialize the imposed boundary flux.

Density, dynamic viscosity, SIMPLE controls, and linear-solver tolerances are currently application defaults in `main.cpp`; they are not yet read from the case dictionaries.
The committed CLI uses Linear convection, momentum BiCGSTAB tolerance `1e-6`, and pressure CG tolerance `1e-3`. Other schemes are available through the numerical APIs and verification drivers.

### `system/controlDict`

Initialization and monitoring can be selected independently:

```text
initialization
{
    type coarseMesh;
}

monitoring
{
    liveConvergence true;
}
```

`coarseMesh` is the default initialization mode. It solves two independently generated auxiliary meshes with default sizes `4h` and `2h`, then initializes the existing final mesh of size `h` by linear WLS reconstruction. Only the first auxiliary level receives boundary-based bulk initialization. Optional `targetCoarseCellCount` sets the coarsest target, not an exact Gmsh cell count.

Use `type zero;` to disable sequencing and preserve the supplied `internalField` values; this does not force nonzero input fields to zero. Unusable automatic hierarchies or transfer coverage fall back to that historical initialization; an explicit target makes those failures errors. SIMPLE/numerical/I/O failures are not silently converted to fallback. Sequencing has no performance-based size cutoff and can cost more on small cases.

`liveConvergence` accepts `true` or `false` and defaults to enabled when the entry, or the complete `controlDict`, is absent.

When enabled, the executable launches `tools/live_convergence.py` as a detached process. A plotter-launch failure produces a warning and does not stop the CFD solve. The default plot shows provisional continuity and the two external momentum-equation residuals. The script can also be run manually:

```bash
python3 tools/live_convergence.py \
    cases/poiseuille/results/convergence.csv \
    --diagnostics
```

## Outputs

Before initialization, a run removes only previous `results/initialization/`, `results/solution.vtu`, and `results/convergence.csv`. Other files in `results/` are preserved; filesystem errors abort the run.

After successful coarseMesh initialization and a converged final solve:

```text
cases/poiseuille/results/
├── initialization/
│   ├── level_0.vtu
│   └── level_1.vtu
├── convergence.csv
└── solution.vtu
```

Each auxiliary VTU is written immediately after that level converges, before transfer to the next mesh. It contains scalar `u`, `v`, `p`, vector `U`, and `mass_imbalance`, plus intrinsic cell metadata. Completed-level files remain available after a later failure. Zero mode creates no initialization directory. These are diagnostic files, not restart files.

`convergence.csv` is created for the final SIMPLE solve after initialization, receives one flushed row per completed iteration, and contains:

- provisional and corrected continuity;
- external x/y momentum-equation residuals;
- velocity change;
- Eigen residual estimates and iteration counts for all three linear solves;
- maximum pressure correction.

`solution.vtu` contains cell-centered pressure `p`, velocity `U`, and cell mass imbalance, together with the writer's intrinsic cell metadata. Open it directly in ParaView:

```bash
paraview cases/poiseuille/results/solution.vtu
```

Use `Surface With Edges` to inspect the mesh and choose `U`, `p`, or `mass_imbalance` from the coloring menu.

## Numerical conventions

For an internal face shared by owner `P` and neighbor `N`, the stored face-area vector is owner-oriented:

```text
|S_f| = face length
S_f · (x_N - x_P) > 0
```

The same vector is outward for the owner and inward for the neighbor. Boundary face vectors are outward from their owner cell. Integrated face mass fluxes use this owner-oriented convention.

The `Linear` convection scheme is geometry-aware but does not make a boundedness claim. LinearUpwind uses an implicit upwind matrix and explicit reconstructed correction; Barth-Jespersen limiting is optional. Non-orthogonal diffusion corrections use reconstructed cell gradients. SIMPLE requires a single connected cell domain and supports momentum relaxation factors in `(0, 1]` with Majumdar correction.

## Tests and numerical verification

Run the fast correctness and regression suite with:

```bash
ctest --test-dir build --output-on-failure
```

Optional numerical-verification executables are enabled separately:

```bash
cmake -S . -B build-release \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_TESTING=ON \
    -DBUILD_VERIFICATION=ON

cmake --build build-release --parallel
```

`tests/` contains focused unit and regression tests. `verification/` contains manually run gradient, diffusion, convection-diffusion, momentum, coupled SIMPLE, and performance studies. Numerical verification compares against mathematical reference solutions; physical or experimental validation is not yet claimed.

## Repository layout

```text
include/cfd/        Public fields, mesh, I/O, linear algebra, and numerics
src/                Implementations and the main case runner
tests/              Fast CTest targets
verification/       Manual numerical-verification and performance drivers
tools/              Live convergence plotter
cases/              Supplied Poiseuille and backward-facing-step cases
```

## Design principles

- Preserve mathematical meaning in types and interfaces
- Validate imported data and public numerical inputs
- Keep large storage ownership explicit and contiguous
- Avoid allocations and field-sized copies in numerical hot loops
- Separate fast regression tests from numerical-verification studies
- Measure performance before optimizing
