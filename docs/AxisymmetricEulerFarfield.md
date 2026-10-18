# Euler farfield: decision before implementation (lot F3)

## Sources actually inspected

- [SU2 v7 markers](https://su2code.github.io/docs_v7/Markers-and-BC/) and
  [physical definition](https://su2code.github.io/docs_v7/Physical-Definition/):
  MARKER_FAR uses a complete freestream, with only two independent thermal inputs.
  [SU2 v7.5.1 CEulerSolver.cpp](https://github.com/su2code/SU2/blob/v7.5.1/SU2_CFD/src/solvers/CEulerSolver.cpp),
  `BC_Far_Field`, lines 3938–4117: the code negates its geometric normal,
  selects acoustic invariants and entropy/tangent using freestream normal speed,
  then calls its configured numerical flux. The nearby incoming/outgoing prose
  is not adopted as our sign convention. This is not a copy of SU2.
- [NASA Cart3D, Full State with Riemann Solver](https://www.nas.nasa.gov/publications/software/docs/cart3d/pages/howto/samples_power/README.html):
  complete exterior primitives provide both incoming thermal and directional
  information. Its use of the interior numerical solver is compared below;
  the linked AIAA papers were not read and no equations are attributed to them.
- [Giles, CFDL-TR-88-1 (1988)](https://people.maths.ox.ac.uk/~gilesm/files/bcs.pdf),
  introduction and normal-characteristic discussion: ignoring tangential
  derivatives is not multidimensional nonreflection. The
  [1990 journal article](https://doi.org/10.2514/3.10521) was identified, not read.
- [Clawpack Riemann book, Euler chapter](https://www.clawpack.org/riemann_book/html/Euler.html),
  integral curves, jump conditions and exact-solution construction: intersect
  the two wave curves in pressure/normal velocity, then sample the entropy
  solution. The implementation below is independently derived, not copied.
- [NASA axisymmetric jet](https://www.grc.nasa.gov/www/wind/valid/axinoz/axinoz.html):
  a future independent jet comparison with specified ambient and plenum data.
  Its viscous/turbulent velocity profiles are not an exact Euler reference.

## Alternatives and choice

Outward normal n=(nx,nr), tangent=(-nr,nx). Eigenvalues are un-a,un,un,un+a;
negative means incoming. On a smooth isentropic acoustic wave, J+=un+2a/(g-1)
is constant through a left-family wave and J-=un-2a/(g-1) through a right-family
wave. These are not simultaneous invariants across arbitrary shocks.

1. Combine J+ interior/J- exterior, then pick entropy/tangent: cheap normal
   approximation, but switching the entropy at un=0 with fixed reconstructed
   a can jump pressure when the two entropies differ. An interior sonic switch
   can also jump when exterior data are incompatible. Not selected.
2. Full exterior state + Rusanov: uses the existing solver, but its scalar
   dissipation mixes incoming/outgoing families even at supersonic outflow;
   exterior pressure can affect a purely outgoing flux. Not selected.
3. **Selected:** full exterior state + local exact normal Euler Riemann trace,
   then physical trace flux, with no second Rusanov dissipation. This is a
   boundary-only perfect-gas calculation, not a new generic Riemann framework.
   Internal Rusanov, existing characteristic BCs, wall and axis are unchanged.

Let L=interior, R=exterior, g=gamma, c=(g-1)/(g+1). For side k, define

```
f_k(p) = (p-p_k)*sqrt(2/((g+1)*rho_k)/(p+c*p_k))      if p>p_k
       = 2*a_k/(g-1)*((p/p_k)^((g-1)/(2*g))-1)        otherwise.
f_L(p*) + f_R(p*) + un_R-un_L = 0
un* = (un_L+un_R+f_R(p*)-f_L(p*))/2
```

The pressure equation is monotone. Bracket/bisect in log pressure without a
pressure floor. Reject vacuum-generating data
`un_R-un_L >= 2*(a_L+a_R)/(g-1)` and nonrepresentable states explicitly.
Rarefactions preserve p/rho^g and tangential velocity; shocks use
Rankine–Hugoniot, rho*/rho_k=(p*/p_k+c)/(c*p*/p_k+1).
Sample x_normal/t=0, including sonic fan interiors. Entropy/tangent follow
the upstream side of the contact; at zero contact speed choose L. Density
can jump at a stationary contact, but pressure/velocity and **flux** agree.
No arbitrary clipping, sonic epsilon or group-wide regime switch.

## Regime and conservation contract

Regime diagnostics use the **sampled trace**, not total Mach or merely the
interior state: un_b/a_b<=-1 supersonic inflow, -1<...<0 subsonic inflow,
0<=...<1 subsonic outflow, >=1 supersonic outflow. Zero-speed characteristics
are stationary; endpoint names are deterministic reporting conventions.
Incoming counts off the endpoints are 4/3/1/0. Every face is independent.

If all waves leave, the sample is L and ambient changes cannot alter the flux.
If all enter, it is R. A very strong exterior compression can send a shock
into a supersonic *interior* state: this is no longer the purely outgoing
branch. Ignoring that shock would not be a full-state Riemann boundary.
This distinction is explicit, not an untested claim that interior Mach alone
determines nonlinear wave propagation. Sonic fan sampling gives continuous
fluxes when a wave crosses zero speed; contact crossing also has continuous
flux. Incoming entropy and tangential velocity come from R, not L.

Flux is `[rho*un,rho*ux*un+p*nx,rho*ur*un+p*nr,(rhoE+p)*un]` on the sampled
trace, multiplied once by the existing revolution area. The CFL bound covers
both initial states, both star states and shock speeds. One boundary flux
enters the FV residual and global ledger; radial pressure source is unchanged.
There is no promise of exact 2D/nonlinear nonreflection or preservation of
ambient pressure at an outflow. Transverse/source evolution remains in the
FV solver. One fixed gamma/R applies to both sides; no species, viscosity,
rarefied/vacuum model or physical RPA ambience is inferred.

## Implementation and verification

`EulerBoundaryCondition::farfield(PrimitiveState)` owns immutable exterior data.
The independent `compute_euler_farfield_boundary_flux()` returns trace, physical
flux, wave bound, normal Mach and `EulerFarfieldRegime`. It validates the EOS
using the same gas as the interior. No ghost fields, allocations on valid face
evaluations, per-step deep copies or dynamic dispatch are introduced. The Euler
assembler adds one switch branch; existing BC formulas and internal Rusanov
are untouched. The legacy characteristic evaluator explicitly rejects Farfield
so it cannot accidentally apply a different closure.

Pressure is solved with long-double wave curves and bounded log-pressure
bisection. Besides stagnation in log space, the search stops at adjacent
representable double pressures: a near-equilibrium root with log(p*) near zero
does not need an arbitrarily large number of log-space digits. No residual
multiplier, pressure floor or empirical sonic threshold is used. Nonzero values
underflowing to zero are rejected. Precision/range gains from long double are
platform dependent. At nearly tangential flow, floating-point projection can
give a tiny signed Mach; the reported sign is not silently clipped.

Six dedicated unit-test groups check independent primitive fluxes; uniform
rest, normal-sub/supersonic and oblique states; all-in/all-out causality;
incoming strong shocks; two-rarefaction acoustic invariants; upstream entropy
and tangential velocity; a stationary contact approached from both sides;
Sod star values; independent symmetric-compression pressure; both sonic fan
heads; normal reversal; near-vacuum expansion and explicit vacuum rejection;
invalid/EOS-unrepresentable states; uniform scaling at 1e-200/1e200 and the
smallest representable positive pressure endpoint; full Euler
uniform preservation and conservative updates on a group with four actual
mixed-regime mesh faces. References are not computed using the tested flux.

The existing boundary-wave helper has one optional farfield-outlet argument;
its default and previous PressureOutlet/Transmissive tests are unchanged.
The new driver reuses its mesh, pulse, extended-domain comparison and the
accepted-step four-equation conservation ledger, not another CFD assembler.
SIMPLE's field/BC and convergence representations are specific to momentum/
pressure coupling, not suitable for an Euler characteristic trace. No new
framework or duplicated general-purpose output writer was needed.

### Acoustic reflection measured against an extended domain

Same dx/dr, gamma=1.4, R=1, rho=p=1, base Mach=0.3, radius=0.1, right-running
isentropic pressure pulse of relative amplitude 1e-4 centred at x=0.65 with
width=0.04. Short length=1, reference length=2, Nr=3, CFL=0.2. Compare at
t=0.85/a0; remote-boundary continuum signals cannot return to the comparison
region by that time. Numerical diffusion tails are not claimed compact.
The reflected pressure estimate is max(abs(p_short-p_extended)); its ratio
uses the initial sampled peak, including its sampling and numerical diffusion.

| Nx | Boundary | Reflected pressure | Reflection / incident | Incoming acoustic amplitude |
|---:|---|---:|---:|---:|
| 60 | PressureOutlet | 1.3083044e-5 | 0.13663388 | 2.2114442e-5 |
| 60 | Transmissive | 1.2308488e-11 | 1.2854473e-7 | 2.0806859e-11 |
| 60 | Farfield | 6.8024475e-12 | 7.1041936e-8 | 1.1497797e-11 |
| 120 | PressureOutlet | 1.8190810e-5 | 0.18389268 | 3.0748199e-5 |
| 120 | Transmissive | 1.3329338e-11 | 1.3474758e-7 | 2.2536220e-11 |
| 120 | Farfield | 8.9811492e-12 | 9.0791315e-8 | 1.5178126e-11 |

All cases have zero rejected steps; cumulative conservation defects are below
8.0e-17 in these normalized units. The farfield guard requires a finite
reflection estimate smaller than PressureOutlet, not perfect nonreflection.
Both Farfield and Transmissive work extremely well on this aligned weak simple
wave. Their very small difference does not establish general superiority:
Transmissive lacks a complete incoming reservoir, while this test does not
measure oblique-wave or strong-shock reflection. Fixed static pressure clearly
reflects this outgoing acoustic perturbation. The reference and incoming
amplitude distinguish reflection from merely observing a diffused pulse.

### Nonuniform backflow demonstration and ParaView

A separate 120x3 tube uses the explicitly **normalized verification gas**,
ambient rho=p=1, ux=-1e-4*a0, ur=0, and a right-running pulse of amplitude
1e-3. Both ends are Farfield, outer wall SlipWall, axis Axis. The right
boundary starts with inflow, becomes outflow as the pulse arrives, then
returns to inflow. It is not an RPA ambience or a stationary solution.

Measured: 1282 accepted steps, zero rejects, minimum final rho=0.99999998359,
p=0.99999999666, right trace Mach range [-0.000100000000000059,
0.000243655248199465]. Cumulative maximum balance defect=1.028e-17;
maximum accepted-step defect=1.351e-17, including radial pressure source.
The solver enforces positivity at every accepted state, not merely at exports.

Eight computed VTU snapshots and a chronological `farfield.pvd` are written
at t*a0=0,0.1,0.2,0.3,0.4,0.5,0.65,0.85. Reusable output fields contain rho,
p,T,ux,ur,velocity_magnitude,Mach,rhoE. Per-face CSV diagnostics contain the
sample regime/Mach and integrated outward mass flux on the current accepted
state, using exactly the production boundary evaluator. They are not the
previous step's flux and are not exported as fictitious cell fields.

Independent XML checks cover all eight snapshots: 360 QUAD cells, 484 points,
valid orientation/connectivity, finite fields, positive rho/p/T, EOS and energy,
nonuniform evolution, chronology and both signs of right boundary mass flow.
Maximum derived-field discrepancy is 5.33e-16. ParaView PVDReader loads all
eight times/fields; no GUI or rendering was used.

One instrumented Release run (not a controlled performance campaign) gives:
elapsed=0.695 s, assembly=0.559 s, all-boundary evaluation=0.256 s **within
assembly**, update/admissibility/summary=0.0951 s, verification diagnostics=
0.0254 s, export=0.0156 s. The explicitly requested observation-only farfield
sample takes 0.02374 s for 816 evaluations (about 29.1 microseconds/call),
excluding CSV, **within diagnostics**. This is a sample, not a new timer in
the production face loop. Geometry/field setup is outside this loop timing.
Per-boundary clock overhead is included; instrumentation is disabled by
default. These times are observations under concurrent build/test load, not
an optimization claim or an A/B performance comparison.

## F3b feasibility: nominal near-space case, not the 1000 Pa alternative

The original `Rapport_propulsion_GOx-GCH4.pdf`, `04_Results_RPA.pdf`,
`05_Input_NASA-CEA.txt` and `06_Output_NASA-CEA.txt` were consulted again.
Report section 2.6.1/equation 2.31 sets the nominal environmental pressure
to **1e-4 Pa** for the 200–800 km mission. Section 3.3.2 at 1000 Pa is a
different optimum-expansion study; Fluent's 1 Pa gauge plus 1 Pa operating
pressure is another convention. None is silently substituted for another.
The user confirmed the nominal near-space scenario during F3.

Retain RPA chamber=12 bar, O/F=3.1, geometric area ratio=100 and nominal
mdot=0.4230 kg/s. With gamma=1.177 and R=414.2 J/(kg K) frozen at the nozzle
inlet, p=1.1700 MPa, T=3343.5658 K, M=0.1508 imply independently

```
rho = p/(R*T) = 0.84482331897 kg/m^3
T0 = T*(1+(gamma-1)*M^2/2) = 3350.294882116 K
P0 = p*(T0/T)^(gamma/(gamma-1)) = 1185747.250123 Pa
a = sqrt(gamma*R*T) = 1276.727508626 m/s
p_infinity/P0 = 8.4335004774e-11
```

The frozen sound speed differs from the tabulated RPA thermochemical sound
speed; do not substitute one into the other model. Inlet P0/T0 are derived
under this approximation, not extra original RPA outputs. mdot is a comparison
target, not prescribed alongside both P0 and T0. One gamma/R everywhere is an
explicit equivalent-gas approximation, not combustion or a composition model.

The exact private meshDict was read and regenerated without alteration:
QUAD, baseCellSize=0.02, quality=0.05, downstreamFactor=20, outerRadiusFactor=10,
bandThicknessFactor=0.5, plumeOutletRadiusFactor=6, wallThicknessFactor=0.02,
upstreamFactor=auto. The current mesh has 2240 QUAD, 4675 faces and one connected
cell component; revolution geometry is valid.

| Actual group | Faces | Outward nx | Candidate BC / required data | Data and startup status |
|---|---:|---|---|---|
| inlet | 8 | -1 | TotalConditionsInlet, derived P0/T0, axial direction | Local hot-gas data sufficient; compatible subsonic initial state still needed |
| axis | 167 | 0 | Axis, exactly zero surface | Compatible |
| inner_wall | 123 | -0.6274 to 0.5 | SlipWall | Compatible with inviscid model only |
| outer_wall | 44 | -0.5 to 0.5932 | SlipWall | Same |
| nozzle_lip | 1 | -1 | SlipWall | Same |
| farfield_upstream | 11 | -1 | Farfield, complete exterior state | All normal regimes supported by F3a; nominal exterior thermal/directional state unspecified |
| farfield_radial | 20 | 0 | Farfield, complete exterior state | Same; entrainment/backflow not replaced with a wall |
| farfield_downstream | 16 | +1 | Farfield, complete exterior state | Same; cannot assume the whole group is supersonic |

### Evaluation of a vacuum-compatible formulation (no invented reservoir)

The pressure 1e-4 Pa does not determine rho or incoming entropy/energy. Neither
300 K (cold propellant inlets) nor 2.7 K (Fluent backflow setting without a
complete physical ambience) is adopted. No atmospheric composition or ambient
velocity was invented. A finite-positive full-state Farfield therefore cannot
be populated from these sources alone.

A mathematically distinct alternative is **exact vacuum outside**, rho_R=p_R=0.
The [Clawpack authors' exact Euler solver](https://github.com/clawpack/riemann_book/blob/master/exact_solvers/euler.py)
explicitly distinguishes an initially vacuum side and a vacuum generated
between two rarefactions. This is not the positive-state branch implemented
above. Independently, for interior state L expanding normally into right vacuum,
the head and free-surface speeds are un_L-a_L and un_L+2*a_L/(gamma-1).
At the fixed boundary xi=0:

```
un_L >= a_L: interior trace (all outgoing)
un_L <= -2*a_L/(gamma-1): vacuum trace, zero flux
otherwise: sonic rarefaction trace,
  a_b = 2/(gamma+1)*(a_L+(gamma-1)*un_L/2), un_b=a_b
  rho_b/rho_L=(a_b/a_L)^(2/(gamma-1))
  p_b/p_L=(a_b/a_L)^(2*gamma/(gamma-1))
```

This requires no arbitrary vacuum temperature/composition and provides an
entropy-admissible Euler **mathematical limit**, not a specification of a finite
1e-4 Pa reservoir. Applying zero pressure to the existing PressureOutlet is
not equivalent: it violates its positive/subsonic contract and misses sonic
fan sampling. Applying Transmissive everywhere likewise is not a vacuum BC.
A tiny positive prescribed pressure cannot force the existing subsonic outlet
past its choking/backflow restrictions either; those validations stay enabled.

There are two unresolved, separate issues before an RPA experiment:

1. The current gas/state/time-step contracts require strictly positive cell
   rho/p. A true initially empty exterior is not admissible, and a zero-state
   face cannot be passed to the present EOS/Rusanov. Extending vacuum support
   requires explicit tests of flux/CFL/source/update and startup; silently
   installing small rho/p floors would be a different, unvalidated model.
   A small positive numerical background could only be an explicitly approved
   exploratory regularization with declared parameters and sensitivity study.
2. Euler is a continuum approximation. A space plume includes regions where
   local molecular mean free path is not small compared with geometric or
   gradient length. [NASA-CR-190515, section 1.1 and figure 1](https://ntrs.nasa.gov/api/citations/19920021737/downloads/19920021737.pdf)
   explicitly distinguishes continuum, transitional and free-molecular plume
   regions. It does not justify Euler in the whole rarefied backflow region.
   No quantitative Knudsen number is asserted without molecular collision
   data, density/temperature and a local length scale. Positive tiny p alone
   does not establish continuum validity. No rarefied solver is added here.

Startup candidates are therefore assessed, not silently selected: uniform
ambient lacks the required state and can violate the hot inlet's outgoing
invariant/root; hot nozzle plus ambient exterior and quasi-1D internal data
still need an admissible exterior/free interface; gradual startup needs a
separate justified time-dependent prescription. Existing inlet validation is
not disabled. The full-domain frozen-hot-gas fill is not presented as physical
near-space initialization.

**Decision: F3a is independently deliverable; nominal F3b is not launched.**
The next minimal decision is whether to authorize a separately verified
exact-vacuum Euler-limit experiment, including its startup/state support, or
to specify an explicitly illustrative positive background and its sensitivity
checks. Neither would be labeled a predictive validation of the rarefied
200–800 km external plume. No 1000 Pa substitute run was launched.

## Reproduction and software validation

```
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DBUILD_VERIFICATION=ON
cmake --build build-release --parallel
ctest --test-dir build-release --output-on-failure
./build-release/cfd_axisymmetric_euler_farfield /tmp/cfd-f3-farfield
paraview /tmp/cfd-f3-farfield/farfield.pvd
```

Use ParaView's time controls and color by Mach, then p/T/rho/speed. CSV regime
codes are SupersonicInflow=0, SubsonicInflow=1, SubsonicOutflow=2,
SupersonicOutflow=3. CSV describes faces; VTU describes cells.

Full builds and CTest pass **53/53 Debug, 53/53 Release, 53/53 ASan/UBSan**.
Leak detection is disabled for the existing Gmsh sanitizer configuration;
LeakSanitizer is not claimed. Existing tests remain enabled/unmodified.
Release Sod, radial expansion, legacy boundary-wave and both analytic N24
nozzle branches pass. Full prior fine nozzle campaigns were not repeated.
Targeted format, strict tidy and whitespace checks are reported with the task.

Private meshDict/DXF checksums are unchanged; no case outputs, Mesh, generator,
SIMPLE or main.cpp were changed. F3 output is outside the private directory.
