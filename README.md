# CFD Solver

CFD_solver is a modern C++20 finite-volume solver for steady, two-dimensional incompressible flow. The project emphasizes explicit ownership, contiguous numerical storage, focused validation, numerical verification, and measured performance work.

The current solver supports triangular and quadrilateral cell-centered meshes and a case-based command-line workflow.

## Implemented capabilities

- Rectangle meshing through the Gmsh C++ API
- Triangular and recombined quadrilateral cells
- Validated topology, owner/neighbor connectivity, boundary groups, cell geometry, and oriented face-area vectors
- Fixed-cardinality scalar, vector, velocity, face-flux, and pressure-response fields
- Inverse-distance-weighted least-squares gradient reconstruction
- First-order upwind and geometry-aware linear convection
- Isotropic diffusion with non-orthogonal correction
- Momentum-weighted/Rhie-Chow face-flux interpolation
- Steady SIMPLE pressure-velocity coupling
- Scalar finite-volume sparse systems
- Eigen BiCGSTAB momentum solves and conjugate-gradient pressure-correction solves
- VTU cell-data export for ParaView
- Per-iteration convergence CSV output and optional live plotting
- Structured and Gmsh numerical-verification drivers

The detailed current SIMPLE contract is documented in [IncompressibleSimpleV1.md](docs/IncompressibleSimpleV1.md).

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

`meshDict` defines the rectangular geometry and Gmsh controls:

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

### Initial and boundary fields

The files `0/u`, `0/v`, and `0/p` provide:

- the field object name and dimensions;
- a uniform cell-centered initial value;
- one named condition for every mesh boundary group.

The supported scalar conditions are uniform `fixedValue`, `fixedGradient`, and `zeroGradient`. Velocity components are currently supplied as separate scalar files. Pressure `fixedValue` boundaries map to `FixedPressure`, while pressure `zeroGradient` and `fixedGradient` boundaries map to `FixedMassFlux` in SIMPLE. The current case runner requires both velocity components to be `fixedValue` on a fixed-mass-flux boundary so it can initialize the imposed boundary flux.

Density, dynamic viscosity, SIMPLE controls, and linear-solver tolerances are currently application defaults in `main.cpp`; they are not yet read from the case dictionaries.

### `system/controlDict`

The supported monitoring control is:

```text
monitoring
{
    liveConvergence true;
}
```

`liveConvergence` accepts `true` or `false` and defaults to enabled when the entry, or the complete `controlDict`, is absent.

When enabled, the executable launches `tools/live_convergence.py` as a detached process. A plotter-launch failure produces a warning and does not stop the CFD solve. The default plot shows provisional continuity and the two external momentum-equation residuals. The script can also be run manually:

```bash
python3 tools/live_convergence.py \
    cases/poiseuille/results/convergence.csv \
    --diagnostics
```

## Outputs

Each run creates the case-local `results/` directory and convergence CSV. After a successful converged solve, the directory contains:

```text
cases/poiseuille/results/
├── convergence.csv
└── solution.vtu
```

`convergence.csv` is truncated at startup, receives one flushed row per completed SIMPLE iteration, and contains:

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

The `Linear` convection scheme is geometry-aware but does not make a boundedness claim. Non-orthogonal diffusion corrections use reconstructed cell gradients. SIMPLE v1 currently requires a single connected cell domain and a momentum relaxation factor of exactly one.

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
cases/poiseuille/   Supplied steady incompressible case
```

## Design principles

- Preserve mathematical meaning in types and interfaces
- Validate imported data and public numerical inputs
- Keep large storage ownership explicit and contiguous
- Avoid allocations and field-sized copies in numerical hot loops
- Separate fast regression tests from numerical-verification studies
- Measure performance before optimizing
