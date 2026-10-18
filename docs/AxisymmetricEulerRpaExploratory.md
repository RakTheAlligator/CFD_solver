# RPA Euler exploration (F3b, approach A)

This is a transient numerical experiment, **not a validated prediction of a
space plume**. It adds a verification driver, not another Euler solver. The
F3a Farfield implementation is a prerequisite; its numerical responsibilities
remain separate from the F3b driver, support header, tests and output validator.
SIMPLE, the private case and production numerical formulas are not changed.

This document records the original F3b experiments. Their temporary output files
are no longer available. For current checkpoint/restart support and persistent
output storage, see [AxisymmetricEulerRestart.md](AxisymmetricEulerRestart.md).

## Data and explicit assumptions

The original `Rapport_propulsion_GOx-GCH4.pdf` (nominal environment, section
2.6.1), `04_Results_RPA.pdf` (table 2) and NASA-CEA input/output were inspected.
The design is Pc=12 bar, O/F=3.1, Ae/At=100, nominal mass flow 0.423 kg/s.
Mass flow is a comparison, not an additional imposed inlet condition.

The frozen equivalent gas uses gamma=1.177, R=414.2 J/(kg K) everywhere.
Static inlet p=1.1700 MPa, T=3343.5658 K, M=0.1508 give
P0=1185747.2501228245 Pa and T0=3350.294882116409 K using perfect-gas
isentropic relations. These derived totals are not exact thermochemical RPA
results. The equivalent sound speed is 1276.7275 m/s, not RPA's 1248.8984 m/s.

The nominal ambient pressure remains **1e-4 Pa**. Ambient T=300 K,
ux=ur=0 and the use of the same equivalent hot gas outside are explicitly
authorized **exploratory assumptions**, not measured conditions at 200–800 km.
The EOS gives rho=8.047642040882022e-10 kg/m3. No pressure/density floor or
exact-vacuum treatment is introduced. The 1000 Pa optimal-expansion example
is a different case and is not used here.

The model includes neither chemistry/species nor viscosity/turbulence.
Euler's continuum assumptions need not hold in this extremely rarefied
exterior. No molecular collision data or Knudsen-number validation is supplied.
Numerical admissibility and conservative balances do not establish continuum,
discretization or mission-level physical validity.

## Mesh, boundaries and initial state

Read the existing private `cases/rpaNozzleSimple/system/meshDict` and DXF
without modifying them. All RPA generation options are passed to the existing
generator. The measured mesh has **2240 QUAD, 4675 faces, eight groups**, one
connected fluid component and minimum stored shape quality 0.10959655.
That shape metric is not the angular-quality threshold. MeshBuilder and
AxisymmetricGeometry validate orientation, positive volumes and revolution
metrics; axis surfaces may vanish. Face fluxes use revolution areas, not
meridional lengths.

- `inlet`: existing TotalConditionsInlet with the derived P0/T0, direction +x.
- `axis`: Axis.
- `inner_wall`, `outer_wall`, `nozzle_lip`: SlipWall.
- `farfield_upstream`, `farfield_radial`, `farfield_downstream`: existing
  F3a full-state Farfield, each face independently sampled.

Names must match all eight groups exactly; missing, duplicated, unknown or
misindexed groups are rejected. Actual RpaRegionLayout cell-centre labels
initialize 984 Engine cells to the hot inlet state and 1256 other cells to
ambient. `rpa_region` is geometry, not a species fraction.
The uniform axial Engine state is not tangent to every wall or steady through
the varying area; wall transients and the hot/ambient discontinuity are part
of this deliberately simple initialization, not an established nozzle flow.

## Reproduction and safeguards

Configure CMake with an existing Python containing `tools/requirements-cad.txt`
(`ezdxf==1.4.3`). The run used Python 3.12.3 in `~/.venvs/cfd-rpa`, explicitly
selected with `Python3_EXECUTABLE`; no ephemeral Python path is committed.

```sh
cmake -U 'Python3_*' -S . -B build-release \
  -DCMAKE_BUILD_TYPE=Release -DBUILD_VERIFICATION=ON \
  -DPython3_EXECUTABLE="$HOME/.venvs/cfd-rpa/bin/python"
cmake --build build-release --target cfd_axisymmetric_euler_rpa_exploratory
./build-release/cfd_axisymmetric_euler_rpa_exploratory \
  cases/rpaNozzleSimple /tmp/cfd-rpa-smoke 300 0 0.01 200 25
./build-release/cfd_axisymmetric_euler_rpa_exploratory \
  cases/rpaNozzleSimple /tmp/cfd-rpa-run1 300 0 0.01 3000 100
python3 verification/validate_rpa_euler_output.py /tmp/cfd-rpa-run1
paraview /tmp/cfd-rpa-run1/rpa.pvd
```

Use fresh output paths if these already contain results. Both commands start
from t=0; the second is **not a restart**. Output inside the input case,
including through symlinks, and nonempty output directories are rejected.
No cleanup or overwrite of previous runs is performed.

Arguments after CASE/OUTPUT are T_inf, ux_inf, CFL, accepted-step cap,
save interval, optionally `ur_inf max_time diagnostic_interval max_rejected
max_wall_seconds` in that order. Defaults for the optional arguments are
0 m/s, 1 s, 25 accepted steps, 20 cumulative rejects, 600 s of loop wall time.
Limits are checked between updates; a single update/output is not interruptible.
Invalid/nonfinite numbers, nonpositive limits and excessive snapshot budgets
are rejected. The per-step reported CFL is the production **upper bound**,
not an independently measured effective CFL after subcycling/time clipping.

The existing EulerSteadyMonitor is used without relaxing its thresholds:
residual/change/flow-change tolerances 1e-7, flow-imbalance tolerance 1e-6,
three consecutive checks after at least 100 steps. Its hot-gas reference
scales and whole-domain normalization are not a local relative-error measure
for the very dilute ambient; small exterior disturbances can be masked by
those scales. Absolute residuals, positivity and boundary diagnostics remain
available. No new convergence criterion is invented.

Outputs include initial/scheduled/last-accepted VTU and chronological PVD,
`steps.csv` (every accepted step), `convergence.csv`, `cell_diagnostics.csv`
and `boundary_groups.csv` (sampled every 25 steps by default), `assumptions.txt`
and `status.txt`. Cell diagnostics supply min T/max Mach at observations;
min/max rho/p and conservative balances are recorded every accepted step.
Boundary traces use the same existing BC evaluators and physical areas, but
describe the current state, not the flux used for the preceding update.
Status distinguishes convergence, step/time/wall/reject limits and failures.
An exception preserves accepted data and reports its phase; a rejected state
is never exported as a successful final solution.

## Measured first runs

| Run | Accepted/rejected | Physical time (s) | Nozzle acoustic crossings | Max Mach | Max p (Pa) | Max rho (kg/m3) |
|---|---:|---:|---:|---:|---:|---:|
| Smoke | 200 / 0 | 3.78718455e-7 | 7.43676315e-4 | 1.364457 | 1215806.69 | 0.87258302 |
| Longer | 3000 / 0 | 5.80641245e-6 | 0.0114018510 | 2.229139 | 1291909.12 | 0.91879120 |

Both end at **STEP_LIMIT, not Converged**. At all accepted steps p>=1e-4 Pa,
rho>=8.04764204e-10 kg/m3. Observed min T remains 300 K to rounding;
final max T is 3394.7285 K. No invalid state or pressure clipping was needed.
The nozzle acoustic time L/a_hot=5.09251737e-4 s, inlet convective time
L/u_hot=3.37700091e-3 s. The longer run covers only 0.0017194 inlet convective
crossing. These intervals are insufficient to develop a stationary nozzle/jet.

Normalized residuals (rho, rhoux, rhour, rhoE) evolve from
(0.00241784, 0.00173294, 0.000263789, 0.0118755) initially to
(0.00256750, 0.00168634, 0.000375340, 0.0129207) at step 3000.
There is no simultaneous residual convergence. Final conservative changes
are (2.43429e-7, 1.59763e-7, 3.55173e-8, 1.22490e-6).

Final outward boundary mass flows [kg/s] are:

| Group | Smoke | Longer |
|---|---:|---:|
| inlet | -0.417436044721 | -0.417436044721 |
| farfield_upstream | 2.21865e-26 | 3.40158e-25 |
| farfield_radial | 2.06030e-25 | 3.15879e-24 |
| farfield_downstream | 2.28611e-26 | 3.49738e-25 |

The inlet differs from nominal 0.423 kg/s by -1.3154%; it is not independently
forced to that flow. Farfield trace pressures remain 1e-4 Pa and trace normal
Mach magnitudes below 5e-18. Their subsonic in/out counts (4/7 upstream,
1/19 radial, 5/11 downstream) reflect roundoff-level signs, **not meaningful
backflow or entrainment**. Inlet has eight subsonic inflow faces. No sonic or
supersonic farfield switch is observed. The whole-domain flow imbalance stays
one while mass accumulates; tiny outlet flow is not a conservation failure.

An independent polygon-integration/time-ledger check gives longer-run defects
of 3.84e-18 kg, 0 axial momentum, 1.11e-16 radial momentum, 8.28e-12 J.
It includes every boundary contribution and the integrated radial pressure
source. Component-scaled defects are <=6.5e-16. The driver's separate
long-double accumulation agrees to rounding. Passing conservation alone
does not demonstrate physical accuracy.

Area-weighted **cell-based estimates**, averaging Engine cells adjacent to
actual throat/exit section faces (not numerical boundary traces), give:

| Section | M at 200 | M at 3000 | p at 3000 (Pa) | T at 3000 (K) |
|---|---:|---:|---:|---:|
| throat | 0.150907 | 0.177964 | 1160810.96 | 3340.2993 |
| exit | 0.162658 | 0.321314 | 896029.98 | 3326.4049 |

The throat has not reached choking (maximum Mach over all Engine cells is
0.321575). The global maximum Mach occurs outside
Engine (Transition, around x=0.6853 m, r=0.2537 m), not at the throat or exit.
Views show an initial expansion localized near the lip/hot–ambient interface,
not a developed plume through the downstream farfield. RPA's steady throat
M=1 and exit M=4.563, T=1399.5206 K, p about 900 Pa are not meaningful
accuracy targets at this very early time. Initialization, insufficient time,
coarse spatial resolution, numerical diffusion and the frozen EOS are distinct
possible discrepancies; they cannot all be attributed to Rusanov diffusion.

Longer-run indicative timings [s]: mesh generation 0.57533, Mesh build 0.001704,
initialization 0.000136, solver construction 0.000662, loop 11.34397;
assembly 9.01222, BC subset 4.22278, update 1.42007, monitoring 0.57294,
sampled Farfield subset 0.16073, export 0.66457. BC time is nested in assembly;
sampled Farfield time is nested in monitoring. Observations can also prepare
the cached production assembly, so assembly and monitoring overlap.
Export includes the initial
snapshot outside the loop timer; connectivity checking and the extra CAD read
for layout are not separately timed. These are not additive disjoint phases
or controlled CPU performance benchmarks. Production BC timing combines all
BC types, so it does not isolate the true per-step Farfield cost.

## Verification and visualization

Dedicated short tests use synthetic public fixtures, never the private DXF:
hot/ambient admissibility, independently checked total conditions, eight-name
BC mapping and invalid groups, real-layout initialization, initial residual
evaluation, one-step four-equation balance and symlink/old-output protection.
Existing input/CAD tests cover meshDict reading; both real runs exercise it.
CLI rejection checks also cover NaN/Inf, invalid CFL/limits and private output.

The independent stdlib Python checker verifies every exported mesh, region
initialization, positive revolution volume, connectivity/orientation, EOS,
derived fields, physical times and four integrated balances. It passed all
9 smoke and 31 longer-run snapshots. ParaView PVDReader loaded every time;
offscreen images were rendered with fixed inter-time ranges for Mach,
log10(p), T, rho and speed. The sandbox-only render failed to access X;
the authorized local-display render succeeded. No interactive GUI opening
is claimed. PNGs and all CFD results remain outside the repository under
`/tmp/cfd-rpa-smoke` and `/tmp/cfd-rpa-run1`.

Next: assess a physically meaningful time horizon and its computational budget
before extending the transient. Do not claim stationarity, tune thresholds,
infer rarefied-plume validity or start viscous/reacting extensions from this
short run alone.
