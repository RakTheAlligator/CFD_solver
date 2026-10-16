# Axisymmetric Euler nozzle — lot F

## Theory established before implementation

Sources read:

- [NASA Glenn, Nozzle Design](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/nozzle-design/),
  continuity/momentum derivation and subsonic/supersonic acceleration.
- [NASA Glenn, Mass Flow Rate Equations](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/mass-flow-rate-equations/),
  total/static relations and critical mass flow. Some rendered intermediate
  equations are malformed; the exponents below follow independently from
  perfect-gas energy and entropy, and agree with the next source.
- [NASA Glenn, Thrust Equations Summary](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/thrust-equations-summary/),
  displayed critical-flow, area–Mach and static/total formulas.
- [NASA NPARC, Converging-Diverging Verification Nozzle](https://www.grc.nasa.gov/www/wind/valid/cdv/cdv.html),
  sections Flow Description, Geometry and Comparison Data; curator J. W. Slater.
  Its references identify J. D. Anderson, *Modern Compressible Flow* (1984),
  and M. S. Liou, AIAA 87-0355 / NASA TM 88926 (1987). These books/papers
  were not themselves read here. The archive illustrates conservation and
  quasi-1D comparisons, not exact equivalence of a 2D solution to quasi-1D.

Assume steady, inviscid, adiabatic, shock-free flow, uniform over each section,
constant gamma/R and negligible transverse kinetic energy. With
`theta=1+(gamma-1)*M^2/2`, energy and entropy give

`T=T0/theta`, `p=P0*theta^(-gamma/(gamma-1))`,
`rho=P0/(R*T0)*theta^(-1/(gamma-1))`, `u=M*sqrt(gamma*R*T)`.

Combining `mdot=rho*u*A` at arbitrary M and at M=1 gives

`A/Astar = (1/M)*[2*theta/(gamma+1)]^((gamma+1)/(2*(gamma-1)))`,

`mdot_star=Astar*P0/sqrt(T0)*sqrt(gamma/R)*[2/(gamma+1)]^((gamma+1)/(2*(gamma-1)))`.

Area ratios are dimensionless; the latter expression has units kg/s.
The area–Mach function has a unique minimum 1 at M=1 and separate monotone
subsonic/supersonic branches. A bracketed independent scalar solver selects
the requested branch; no CFD flux or iteration is used by this reference.
For an unchoked solution, Astar is smaller than the geometric throat area.

## F1: physical stationarity, not a small time increment

Use the **same** Euler assembly as the evolution solver on its current state:
`R_C=-sum(Af*Fhat)+G_C*(0,0,p_C,0)`. A diagnostic evaluation does not advance
time, alter the state, or overwrite accepted-step diagnostics. Its prepared
residual/CFL can serve the next step, because Mesh/gas/steady BCs and state
are unchanged; an accepted update invalidates it. There is no second assembler.

For each equation report `Rabs_k=sum_C abs(R_C,k)`: units kg/s, N, N, W.
Choose explicit fixed positive reference density rho_ref, velocity U_ref
(an acoustic speed is appropriate), and length L_ref. Let
`Qref=(rho_ref,rho_ref*U_ref,rho_ref*U_ref,rho_ref*U_ref^2)`,
`tref=L_ref/U_ref`, `Vtot=sum(V_C)` and

`Rnorm_k=Rabs_k/(Vtot*Qref_k/tref)`.

This is a volume-L1 norm of the local time derivative normalized by physical
scales, independent of cancellation between cells and of dt. The fixed
references transform with a change of units; there is no epsilon denominator
or observed dominant-equation floor. References must reflect the problem's
physical scales; these are not componentwise relative backward errors.

At scheduled checks, also report
`change_k=sum(V_C*abs(U_C,k-U_previous_check,C,k))/(Vtot*Qref_k)`.
This is a separate sampled field-change diagnostic, not a replacement for
stationarity. The previous-check buffer is allocated once and overwritten
in place only at checks. For inward/outward **actual numerical** boundary
mass fluxes, use positive `mdot_in/out`,
`imbalance=abs(out-in)/(out+in)` (zero if both vanish), and
`flow_change=max(abs(in-in_prev),abs(out-out_prev))/(rho_ref*U_ref*Vtot/L_ref)`.
All four residuals, all four changes, both mass-flow checks, minimum accepted
steps and consecutive successful checks must pass together. Limits on accepted
steps and physical time return a distinct non-converged status. NaN/Inf cannot
establish convergence. Histories contain observations of the current state,
not stale pre-update residuals. Per-step balance diagnostics retain their
original meaning and do not prove stationarity.

CSV uses a dedicated small writer because the existing SimpleConvergenceCsvWriter
has a SIMPLE-specific schema/types, not a reusable generic CSV engine. Its
checked stream/precision/flush approach is reused without changing SIMPLE.
The CSV residual columns are normalized; absolute residuals and field/flow
changes are separate columns. `dt` describes the last accepted step whereas
the residual and boundary flux describe the state at the observation time.

## F2: controlled analytic geometry and initial data

Use the NASA CDV cosine **area ratios**, scaled to L=1 and throat radius 0.1:
with z=x/L and `At=pi*Rt^2`,

`A/At=1.75-0.75*cos((2*z-1)*pi)` for z<=1/2,
`A/At=1.25-0.25*cos((2*z-1)*pi)` for z>1/2,
`R(x)=sqrt(A/pi)`.

This gives inlet/throat/exit ratios 2.5/1/1.5, zero slope at the throat
and endpoints, continuous slope at the throat (C1, not globally C2). The dedicated verification
RawMeshData maps a structured grid to radial fractions of R(x); cells and
walls are the resulting straight-sided QUAD discretization, not curved CAD.
The production generators and Mesh are not changed. Boundary groups are
explicitly inlet/axis/wall/outlet with owner-outward normals.

First use a fully subsonic reference with throat Mach 0.5, determining Astar
and compatible exit static pressure from that branch. Use total inlet,
SlipWall, Axis, PressureOutlet. Next use Astar=At and the subsonic/sonic/
supersonic branches on either side of the throat, with a SupersonicOutlet.
This outlet is used only if the actual interior remains normally supersonic;
no PressureOutlet rejection is suppressed to start the choked case.

Uniform initialization is an explicit primitive state converted by the EOS.
Quasi-1D initialization is a distinct selectable strategy, with ur=0; it
is not an exact slip-compatible 2D state on inclined walls. Conservative
cell averages use quadrature weighted by the **actual polygon's** revolution
section area, not point sampling at a meridional centre. The profile is never
injected into interior cells again during evolution.

At each axial mesh plane, integrate the actual numerical mass flux over
revolved rings. Other profile quantities are area-weighted primitive estimates
from adjacent cells (boundary traces at inlet/outlet). Export their quasi-1D
counterparts separately. The true numerical mass flux includes Rusanov
dissipation and must not be confused with an arithmetic mean of rho*u.
Refinement keeps the same domain and thermodynamic BCs. Report discrete-grid
sensitivity, stationarity and model differences separately; do not label
convergence to quasi-1D as an exact 2D spatial order.

## Reuse and scope

Direct reuse: MeshBuilder, revolution geometry, SoA fields, EOS, existing
Euler flux/BC/positivity/CFL logic, VTU writer and verification fixtures.
MassFluxBalance's owner/neighbor signs apply but its scalar face storage is
unnecessary when the Euler assembler already has all flux components. WLS is
not needed by this first-order method. CellFieldTransfer/grid sequencing is
deferred: it neither guarantees conservative integrals nor positivity of an
unlimited reconstruction. SIMPLE pressure corrections/Rhie–Chow and its
convergence criteria are mathematically incompatible and are not imported.

An external RPA domain requires actual ambient inflow/backflow data and a
validated farfield closure, not just a pressure number. F3 must report this
restriction rather than turn open boundaries into walls or invent a reservoir.
No reacting model, implicit solver, pseudo-time, HLLC or second-order scheme
is introduced in this lot.

## Reproduction and measured F1/F2 results

```sh
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release -DBUILD_VERIFICATION=ON
cmake --build build-release --target cfd_axisymmetric_euler_nozzle --parallel
./build-release/cfd_axisymmetric_euler_nozzle /tmp/euler-nozzle all quasi1d
# Optional independent uniform-start subsonic check, not a choked startup:
./build-release/cfd_axisymmetric_euler_nozzle /tmp/euler-nozzle-uniform subsonic uniform 24
paraview /tmp/euler-nozzle/choked/N96/nozzle.pvd
```

Optional arguments select `subsonic|choked|all`, `quasi1d|uniform`, and a
single supported axial count `24|48|96`. Invalid arguments fail before runs.
Uniform choked startup is explicitly rejected: one uniform state cannot
simultaneously satisfy this subsonic-inlet and supersonic-outlet startup.

These are **normalized ideal-gas verification units**, gamma=1.4, R=P0=T0=1,
L=1, Rt=0.1; the flow numbers below are not RPA SI predictions.
For dimensional interpretation choose P0_phys, T0_phys, R_phys and L_phys;
the scales are rho_ref=P0_phys/(R_phys*T0_phys),
velocity=sqrt(R_phys*T0_phys), time=L_phys/velocity and
mass_flow=rho_ref*velocity*L_phys^2. The monitor uses the reservoir acoustic
speed sqrt(gamma) in the normalized calculation, not the zero radial velocity.

Common settings: CFL=0.2, residual/field-change/flow-change tolerances=1e-7,
relative inlet/outlet imbalance tolerance=1e-6, check every 50 accepted steps,
at least 100 steps and three consecutive successful checks, at most 200000
steps and 80 acoustic crossing times. An initial 40-crossing-time bound
correctly returned TimeLimit on N48 (max residual=1.15935e-7); only the common
physical horizon was extended. No convergence tolerance was loosened.
CSV status codes: Running=0, Converged=1, StepLimit=2, TimeLimit=3.

Subsonic reference: Mthroat=0.5, Mexit=0.3042348830, pexit/P0=0.9378185384,
mdot=0.01605521027. Choked reference: Mthroat=1, Mexit=1.854123527,
pexit/P0=0.1601759817, mdot=0.02151147313 (pressure is **not** imposed at
the supersonic outlet).

All six computations converged with **zero rejected steps**:

| Case | Nx/Nr | Cells | Steps | Physical time | Inlet flow | Outlet flow | Mean throat Mach | Mean exit Mach | Flow error vs quasi-1D |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Subsonic | 24/3 | 72 | 23950 | 32.77450 | 0.01172118742 | 0.01172118619 | 0.334660 | 0.223061 | -26.9945% |
| Subsonic | 48/6 | 288 | 51050 | 34.30728 | 0.01338342707 | 0.01338342582 | 0.396554 | 0.254311 | -16.6412% |
| Subsonic | 96/12 | 1152 | 105450 | 35.03850 | 0.01452784314 | 0.01452784171 | 0.440318 | 0.275745 | -9.51323% |
| Choked | 24/3 | 72 | 5850 | 6.764117 | 0.02242901101 | 0.02242901478 | 0.827400 | 1.465001 | +4.26536% |
| Choked | 48/6 | 288 | 10850 | 6.263054 | 0.02196531871 | 0.02196532007 | 0.884366 | 1.640263 | +2.10979% |
| Choked | 96/12 | 1152 | 21200 | 6.126425 | 0.02172644684 | 0.02172644842 | 0.922713 | 1.741936 | +0.999352% |

| Case/Nx | Rrho | Rrhoux | Rrhour | RrhoE | Relative flow imbalance |
|---|---:|---:|---:|---:|---:|
| Sub/24 | 2.22140e-8 | 9.49309e-8 | 5.38621e-9 | 1.87553e-8 | 5.28264e-8 |
| Sub/48 | 2.23834e-8 | 9.72445e-8 | 5.32994e-9 | 2.39437e-8 | 4.66244e-8 |
| Sub/96 | 2.57267e-8 | 9.81512e-8 | 5.24827e-9 | 3.19381e-8 | 4.93689e-8 |
| Choked/24 | 6.75820e-8 | 3.26579e-8 | 1.58869e-9 | 3.15649e-8 | 8.39878e-8 |
| Choked/48 | 2.43190e-8 | 4.80036e-8 | 2.77326e-9 | 5.39008e-8 | 3.08648e-8 |
| Choked/96 | 2.83532e-8 | 6.68241e-8 | 3.96502e-9 | 6.21843e-8 | 3.63817e-8 |

Profile comparisons use area-weighted section primitives. The following RMS
differences use trapezoidal integration in x/L, normalized by 1 for Mach,
P0 for pressure, T0 for temperature, rho0 for density and a0 for velocity.
They are **differences from a reduced model**, not exact Euler 2D error norms.

| Case/Nx | RMS Mach | RMS p/P0 | RMS T/T0 | RMS rho/rho0 | RMS ux/a0 |
|---|---:|---:|---:|---:|---:|
| Sub/24 | 0.105446 | 0.0373641 | 0.0160948 | 0.0240840 | 0.101909 |
| Sub/48 | 0.0655836 | 0.0242954 | 0.0102755 | 0.0159757 | 0.0633372 |
| Sub/96 | 0.0376026 | 0.0143280 | 0.00601693 | 0.00952991 | 0.0364644 |
| Choked/24 | 0.206077 | 0.0570424 | 0.0672488 | 0.0353593 | 0.116299 |
| Choked/48 | 0.113546 | 0.0308965 | 0.0367528 | 0.0188423 | 0.0621646 |
| Choked/96 | 0.0607941 | 0.0170075 | 0.0195337 | 0.0105712 | 0.0333490 |

The section-mean Mach is monotone in all three choked profiles, with its
M=1 crossing at x=0.57485/0.54384/0.52731 rather than exactly x=0.5.
This diagnoses a resolved transsonic branch without a resolved normal shock
in the mean profile; it does not prove absence of every weak 2D wave.
The subsonic case still has a significant first-order grid bias. Do not claim
a mesh-independent accurate nozzle result or a spatial order toward quasi-1D.

In particular, the choked numerical flow exceeding the quasi-1D critical value
is **not** enhanced physical performance. Independent decomposition at the
throat of the Rusanov flux gives:

| Nx | Central advective mass flux | Rusanov mass diffusion | Their sum |
|---:|---:|---:|---:|
| 24 | 0.02016978803 | 0.00225922458 | 0.02242901261 |
| 48 | 0.02079064760 | 0.00117467172 | 0.02196531933 |
| 96 | 0.02112763778 | 0.00059880975 | 0.02172644753 |

The central contributions are below the theoretical critical flow; the positive
Rusanov density-difference term is measured and decreases with refinement.
The refined wall also changes its polygonal approximation and 2D transverse
motion remains absent from quasi-1D. These effects are not conflated with
incomplete stationarity: the residual/balance thresholds are much smaller.

Global ledgers integrate the actual accepted-step numerical boundary fluxes
and radial source, reusing the existing verification balance helper. Maximum
absolute cumulative defect among the four normalized-unit components over
the six runs is 1.54e-16; maximum per-step defect is 1.55e-17.
All section-flow ranges also obey the independently checked bound
`max(mdot_section)-min(mdot_section) <= sum(abs(Rrho_cell)) + summation_roundoff`.
On N96 their relative ranges are 8.93e-8 (subsonic) and 7.35e-8 (choked).

The independent uniform-start N24 subsonic run converges in 23700 steps
versus 23950 for quasi-1D. Final max differences (rho,ux,ur,p) are
3.54e-7, 1.29e-6, 1.18e-7, 4.20e-7: negligible compared with the coarse-grid
discretization discrepancy, not an assertion of bitwise identical steady states.

### Cost and output validation

One instrumented Release run per case, not a controlled repeated benchmark:
elapsed iteration-loop times are 2.27/16.52/124.25 s subsonic and
0.59/3.92/27.43 s choked. CPU contention can affect these wall times.
On N96, assembly including BC/source/CFL takes 98.43/21.67 s;
BC evaluation is a **subset**, 7.52/1.46 s; updates/admissibility/accepted-state
summary take 25.52/5.66 s. Observation time 2.06/0.45 s includes some assembly
and must not be added to it as a disjoint phase. VTU export takes 0.078/0.068 s.
Mesh generation/build/initialization and solver construction are reported
separately by the driver; construction groups revolution metrics and persistent
buffers, so **isolated revolution-metric construction was not timed**.
Clock instrumentation is optional, disabled by default, has boundary-face
overhead when enabled, and is not a basis for an optimization claim.

The six PVD series reference 55 VTU snapshots at actual accepted physical
times (no relabeling or artificial reference fields). Each case also writes
convergence.csv, profiles.csv and a separate quasi1d.csv. All ten computed
VTU fields, mesh orientations/connectivity, EOS, total quantities, chronology,
positivity, nonuniform evolution, convergence histories and section integration
were checked independently. Maximum derived-field discrepancy: 1.0e-15.
ParaView PVDReader loads all six series without GUI/rendering.

## F3: actual private RPA mesh, no artificial simulation

The four original propulsion/RPA/CEA documents were accessible and consulted
again; [the lot-E data specification](RpaCompressibleBoundaryConditions.md)
remains the source distinction for static/total, gauge/absolute and reactive
versus equivalent-gas inputs. The private meshDict was **read, not modified**:
QUAD, baseCellSize=0.02, minCellQuality=0.05, wallThicknessFactor=0.02,
upstreamFactor=auto, downstreamFactor=20, outerRadiusFactor=10,
bandThicknessFactor=0.5, plumeOutletRadiusFactor=6.
An external temporary inspection program generated exactly that mesh:
2240 QUAD cells, 4675 faces, one connected cell component, valid revolution
geometry. Only a clearly labeled mesh-only VTU was written outside the case.

| Actual group | Faces | Revolution area (approx.) | Expected normal regime / proposed BC | Feasibility |
|---|---:|---:|---|---|
| inlet | 8 | 0.00256640 | hot-gas normal subsonic inflow / TotalConditionsInlet | Local inlet data available |
| axis | 167 | 0 | Axis | Compatible |
| inner_wall | 123 | 0.281409 | SlipWall | Compatible for inviscid model only |
| outer_wall | 44 | 0.293475 | SlipWall | Same |
| nozzle_lip | 1 | 0.00259206 | SlipWall | Same |
| farfield_upstream | 11 | 6.41290 | Can be subsonic inflow (outward nx=-1) | PressureOutlet lacks incoming thermal/directional data |
| farfield_radial | 20 | 31.5021 | Tangential ambient limit, entrainment/backflow possible | Ambient reservoir/farfield closure missing |
| farfield_downstream | 16 | 6.41600 | Mixed plume supersonic outflow / ambient subsonic or backflow | Not a uniformly supersonic outlet |

Constant-property hot-gas preparation at gamma=1.1770, R=414.2 J/(kg K),
p=1.1700 MPa, T=3343.5658 K, M=0.1508 gives approximately
T0=3350.294882 K, P0=1.185747250 MPa. Nominal RPA mdot=0.4230 kg/s is
a comparison datum, not another simultaneously imposed inlet constraint.
The 1.2000 MPa chamber value, 1.44 MPa initial Fluent gauge value and
ambient choices 2 Pa / 0.0001 Pa / 1000 Pa have different meanings.

Feasibility decision: geometry **yes**; axis/solid/inlet groups **compatible**;
open farfield groups **not yet fully specified/supported**. Positive local
hot-gas primitives are admissible, but no global ambient-compatible initial
state/BC setup is validated. A physically defensible full Euler run on this
external RPA domain is therefore **not launched**. Missing items are an
explicit ambient thermal/reservoir state, intended ambient pressure scenario,
a separately verified mixed-regime/backflow farfield closure and an admissible
startup. Chemistry/species/viscosity are additionally needed for the original
reacting Fluent problem. No private configuration or generator was changed.

## Validation status

Full Debug/Release builds and full CTest: **52/52 each**. Combined ASan/UBSan
CTest: **52/52**, with `ASAN_OPTIONS=detect_leaks=0` for existing Gmsh constraints
(no LeakSanitizer claim). The short subsonic uniform-start and choked N24
simulations also pass ASan/UBSan; the full fine-grid study is Release only.
Sod, radial expansion and boundary-wave verification drivers pass in Release;
existing SIMPLE, gradient, transfer and mesh unit tests remain enabled/pass.
Targeted strict clang-tidy, clang-format and diff whitespace checks pass.

F1 and the first controlled F2 study are delivered. Quantitative nozzle
accuracy still needs finer grids/assessment of model discrepancy; F3 remains
explicitly blocked by physical boundary data/closure, not by mesh generation.
No full RPA simulation, reacting physics or next solver block is claimed.
