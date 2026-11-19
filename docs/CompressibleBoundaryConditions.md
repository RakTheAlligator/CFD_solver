# Compressible Euler boundary conditions (lot E)

## Sources inspected before implementation

- M. B. Giles, *Nonreflecting Boundary Conditions for Euler Equation
  Calculations*, AIAA Journal 28(12), 2050–2058 (1990),
  [DOI](https://doi.org/10.2514/3.10521). The journal PDF was not accessible in
  this session. The author's accessible earlier report, *Non-Reflecting Boundary
  Conditions for the Euler Equations*, CFDL-TR-88-1 (1988),
  [PDF](https://people.maths.ox.ac.uk/~gilesm/files/bcs.pdf), introduction,
  sections 2.4, 2.8.1 and 3.8 (normal characteristic variables),
  describes normal one-dimensional characteristic conditions and their distinction
  from multidimensional nonreflection. Its coordinate/sign convention is not
  copied: below normals are outward.
- T. J. Poinsot and S. K. Lele, *Boundary conditions for direct simulations of
  compressible viscous flows*, JCP 101(1), 104–129 (1992),
  [DOI](https://doi.org/10.1016/0021-9991(92)90046-2). Bibliographic identification
  was verified, but the full text was inaccessible; no equation is attributed
  to an unread section. This implementation is **not** NSCBC.
- [SU2 v7 boundary documentation](https://su2code.github.io/docs_v7/Markers-and-BC/),
  sections Total Conditions, Mass Flow Inlet, Supersonic Inlet and Pressure Outlet;
  [inviscid bump tutorial](https://su2code.github.io/tutorials/Inviscid_Bump/).
  [SU2 v7.5.1 CEulerSolver.cpp](https://github.com/su2code/SU2/blob/v7.5.1/SU2_CFD/src/solvers/CEulerSolver.cpp),
  `BC_Inlet`/`TOTAL_CONDITIONS` and `BC_Outlet`, directly confirms extrapolation
  of the outgoing acoustic invariant and the outlet entropy/tangential velocity.
  We do not copy its clipping or its total-Mach regime decision.
- [NASA Cart3D propulsive BC documentation](https://www.nas.nasa.gov/publications/software/docs/cart3d/pages/howto/samples_power/README.html),
  sections Stagnation Properties for Subsonic Inflow and Mass Flow Rate and Total
  Temperature for Subsonic Inflow: these are distinct sets of physical inputs.
  The uniform normal mass-flux distribution below is our explicit closure,
  not a claim to reproduce Cart3D's implementation.
- [NASA Glenn mass flow equations](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/mass-flow-rate-equations/),
  mass flow, total/static relations and choking. Some rendered intermediate
  equations on that page are malformed; the relations below are independently
  derived from the perfect-gas energy and entropy equations.

## Equations and sign convention (independent derivation)

Coordinates are `(x,r)`, without swirl. Let `U=(rho,rho*ux,rho*ur,rhoE)`.
The axisymmetric equation is

`dU/dt + dFx/dx + (1/r)*d(r*Fr)/dr = (0,0,p/r,0)`.

The cell-integrated source remains `G*(0,0,p,0)` with `G=2*pi*A_2D`.
For outward unit normal `n=(nx,nr)` and tangent `tau=(-nr,nx)`,
`un=u.n`, `ut=u.tau`. The normal physical flux is

`Fn=(rho*un, rho*ux*un+p*nx, rho*ur*un+p*nr, (rhoE+p)*un)`.

The normal Jacobian eigenvalues are `un-a, un, un, un+a`,
`a=sqrt(gamma*p/rho)`. Negative speeds enter the domain:

| Regime | Incoming | Outgoing | Physical inputs |
|---|---:|---:|---|
| `un < -a` | 4 | 0 | full static primitive state |
| `-a < un < 0` | 3 | 1 | `P0,T0,direction`, or `mdot,T0,direction` |
| `0 < un < a` | 1 | 3 | static absolute pressure |
| `un > a` | 0 | 4 | none |

At sonic/tangential limits some speeds are zero, not incoming or outgoing.
Total Mach cannot substitute for normal Mach. The source does not change
these normal eigenvalues. Normal simple-wave invariants are an approximation
when transverse/source effects or finite shocks matter, not global invariants
of arbitrary multidimensional flows.

## Subsonic inlet

Extrapolate `Jplus=un_i+2*a_i/(gamma-1)` from the interior. Normalize the
finite nonzero supplied direction `d` once; require `mu=d.n<0` on each face.
For speed `q>=0`, prescribed total temperature gives

`a(q)^2=gamma*R*T0-(gamma-1)*q^2/2`,

`f(q)=mu*q+2*a(q)/(gamma-1)-Jplus=0`.

`f` is strictly decreasing on the physical branch. The normal-sonic speed is
`qsonic=sqrt(gamma*R*T0/(mu^2+(gamma-1)/2))`.
Bracket the unique root between zero and `qsonic`; reject a missing root or
sonic inlet instead of clipping. Bounded scalar bisection has no allocation.
An initially stagnant interior is allowed as the continuous subsonic limit;
outward or normally supersonic interior flow is rejected for these inlet types.

For total conditions, `T=a^2/(gamma*R)`,
`p=P0*(T/T0)^(gamma/(gamma-1))`, `rho=p/(R*T)`, `u=q*d`.
No additional independent density/velocity/pressure is prescribed.

For mass flow, `mdot>0` means flow **into** the boundary group. Precompute
`Ag=sum_f |S_axi,f|` once. Prescribe uniform inward normal mass-flux density
`j=mdot/Ag` (not uniform mass per face). Set
`rho=j/(-mu*q)`, `p=rho*R*T`. Thus `sum_f Fn.mass*|S_axi,f|=-mdot`
to floating-point accuracy, even for unequal face areas/normals.
`P0` is an output, not another constraint. There is no universal maximum mdot
at fixed T0 alone when density/total pressure are free. Nevertheless the normal
subsonic root must exist: no promise is made outside this compatible branch.
There is no global iterative feedback/pressure steering; the only nonlocal
operation is the group-area sum at construction.

## Outlet and supersonic boundaries

For `0<=un_i<a_i`, prescribe `pb>0` and extrapolate entropy and tangential
velocity:

`rho_b=rho_i*(pb/pi)^(1/gamma)`,
`ab=sqrt(gamma*pb/rho_b)`, `un_b=Jplus-2*ab/(gamma-1)`,
`u_b=u_i+(un_b-un_i)*n`.

Reject reconstructed backflow or a sonic/supersonic trace: the one-incoming
characteristic closure no longer applies. A very low pb is not silently
clamped. If the **interior** has `un_i>=a_i`, PressureOutlet ignores pb and
extrapolates the interior. SupersonicOutlet requires this same regime.
At exact outward sonic speed extrapolation is the limiting no-incoming case.

SupersonicInlet prescribes a complete valid primitive state. Both its trace
and the interior must have `un<-a`; sonic and reverse-flow states fail with
an explicit diagnostic. Boundary regime is checked on every evaluation, so
supported transitions are immediate, without hysteresis. A subsonic outlet
turning into inflow needs thermal/directional data that this outlet does not
have: it fails rather than inventing incoming information.

## Actual numerical flux and CFL

For the **new** types the characteristic reconstruction is the boundary trace,
not an exterior ghost state. Assemble its physical flux `Fn(Ub)` directly.
Using Rusanov between `Ui` and that trace would generally change the prescribed
mass flux/static-pressure closure through a second dissipation term. The
direct flux is a first-order characteristic boundary discretization; it is
not an exact shock Riemann solver. Use
`max(abs(un_i)+a_i, abs(un_b)+a_b)` for the face wave bound, and multiply the
flux exactly once by the physical revolution area in the existing assembler.
Internal Rusanov, Axis, SlipWall and Transmissive are unchanged.

The face flux used by the update is also the one used by global conservation
diagnostics. No component is overwritten after computing a different flux.
Gas conversions/admissibility use the existing EOS; no pressure/temperature
clipping. Unrepresentable values and incompatible regimes throw. Only supported
characteristic regime/root incompatibilities use `EulerCharacteristicCompatibilityError`;
configuration, geometry, invalid gas states and arithmetic errors remain distinct.
Current-state BC failures propagate without retry. The local steady solver also
checks candidate boundary fluxes before commitment: typed characteristic
incompatibilities may then reject the candidate and halve its pseudo-step.
Other failures remain fatal. Physical time stepping is unchanged.

## Scope and reuse

Reuse directly: Mesh groups/connectivities, AxisymmetricGeometry, gas EOS,
ConservativeField, physical Euler flux, SlipWall, CFL/backtracking, VTU writer.
Adapt: add immutable BC values and a small local trace/flux routine, plus a
precomputed group area in the Euler solver. SIMPLE scalar/velocity BCs,
Rhie–Chow and pressure-correction semantics are incompatible with Euler
characteristic data and are neither reused nor modified. No generic framework,
time-dependent hidden state, per-face allocation or new species model.

FixedPressure is a physical pressure constraint, **not** a nonreflecting outlet.
In linear 1D acoustics its pressure perturbation is zero, so an outgoing
pressure wave reflects with opposite sign. The dynamic verification compares
against an extended domain and reports the reflection rather than claiming
perfect transmission. General nonreflecting multidimensional BCs, backflow
reservoirs, choked inlet switching and reacting thermodynamics remain outside
this implementation.

## Reproducible verification and measured results

```sh
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release -DBUILD_VERIFICATION=ON
cmake --build build-release --target cfd_euler_characteristic_boundary_tests cfd_axisymmetric_euler_boundary_demo --parallel
./build-release/cfd_euler_characteristic_boundary_tests
./build-release/cfd_axisymmetric_euler_boundary_demo /tmp/euler-boundary
paraview /tmp/euler-boundary/boundary_demo.pvd
```

The unit tests independently check all four flux components, outgoing invariant,
total/static conditions, tangential preservation, several oblique orientations,
total-supersonic/normal-subsonic flow, mass flow on 2/5/9 unequal radial face
distributions with varying normals, uniform tube preservation (M=0.05/0.3/0.8/2),
normal sonic limits, stagnation, invalid data and explicit backflow failures.
The mass-inlet test also measures the actual solver boundary flux and mass
increment with the other end closed. Its target is 0.7 kg/s, with a 3e-12
combined absolute/relative test tolerance, not meridional-length weighting.

### Acoustic interaction (measured, not a nonreflection theorem)

Base rho=p=R=1, gamma=1.4, M=0.3, radius=0.1, short length=1,
extended length=2, three radial layers. Initial right-running isentropic
Gaussian pressure amplitude=1e-4, centre x=0.65, width=0.04.
Both domains use the same dx/dr and initial data on their common part.
The reference advances to each short-domain step time, subcycling if its CFL
requires it. At `t*a0=0.85`, the outgoing pulse centre is near x=1.755 in
the extended tube; reflected waves from x=2 cannot reach the common domain
by continuum characteristics. Discrete diffusive tails are not exactly zero.
All common-domain cells contribute to the RMS norm.

| Nx (short) | Outlet | Accepted/rejected steps | Short-reference pressure RMS | Reflected pressure maximum | Reflected / initial incident pressure |
|---:|---|---:|---:|---:|---:|
| 60 | PressureOutlet | 919/0 | 5.947672e-6 | 1.308304e-5 | 0.136634 |
| 120 | PressureOutlet | 1582/0 | 7.012214e-6 | 1.819081e-5 | 0.183893 |
| 60 | Transmissive control | 919/0 | 6.402407e-12 | 1.230849e-11 | 1.285447e-7 |
| 120 | Transmissive control | 1582/0 | 6.878559e-12 | 1.332934e-11 | 1.347476e-7 |

The pressure difference is negative at the reflected peak, consistent with
the sign reversal of linear fixed-pressure acoustics. The measured incoming
characteristic amplitude `max|delta u-delta p/(rho0*a0)|` is 2.211444e-5 and
3.074820e-5 for PressureOutlet. Ratios include first-order Rusanov diffusion
along both travel paths; they are **not** intrinsic boundary reflection
coefficients. Higher reflected amplitude on the finer mesh is not an
instability: less wave attenuation is expected. The largest absolute global
balance defect among the four dimensionless-gas ledger components is 7.97e-17.
Transmissive's small reflection here is a result for this aligned, weak 1D
wave, not validation as a general compressible inlet/outlet.

### Hot equivalent-gas demonstration

The separate nonuniform tube demonstration freezes gamma=1.1770 and R=414.2
J/(kg K) at the RPA nozzle inlet, p=1.1700 MPa, T=3343.5658 K, M=0.1508.
It derives rho=0.844823319, a=1276.727509 m/s, ux=192.530508 m/s,
T0=3350.294882 K and P0=1185747.250 Pa. This is a constant-property surrogate,
not a fit to RPA's tabulated sound speed or a combustion calculation.
The pulse amplitude is 1e-4 relative pressure, L=1 m, R=0.1 m, Nx=120, Nr=3.
It runs to `t*a0=0.85` in 1433 accepted steps, no rejects, with final
rho_min=0.844805213 and p_min=1169970.487 Pa. The largest absolute ledger
defect over mass/momenta/energy is 1.39e-10 (these are dimensional quantities,
not a relative norm; the energy scale is much larger than the unit-gas test).

Six VTU snapshots at `t*a0=0,0.1,0.25,0.45,0.65,0.85` and a chronological
PVD contain rho,p,T,ux,ur,velocity_magnitude,Mach,rhoE. `boundary_faces.csv`
holds face-only normal Mach, BC type and integrated outward mass flux at
those instants; these data are not mislabeled as cell fields or step averages.
Snapshot clipping changes the sequence of time steps only for the output demo.

Independent XML parsing verifies all connectivity/field sizes, positive cell
orientations, positivity, EOS, speed, Mach, energy, chronology and nontrivial
evolution. Maximum relative derived-field discrepancy: 1.41e-16. ParaView's
PVDReader also loads all six time steps (no GUI/rendering claimed).
