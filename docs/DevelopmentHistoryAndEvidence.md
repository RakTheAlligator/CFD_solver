# Development history and evidence

## 1. Purpose, scope and evidence limitations

This record covers the audit of 10 September 2026 at `53bf506` through committed HEAD `75dc51a`: 24 commits. It preserves study results, rejected experiments and decision rationale that commit messages alone cannot establish. The current solver contract is in [IncompressibleSimpleV1.md](IncompressibleSimpleV1.md); the case workflow is in [README.md](../README.md).

Git and selected public contracts establish implementation status. Historical numbers below are supplied development-campaign records, not results rerun while writing this document. Raw logs, complete fields, compiler/CPU metadata and repeat counts are not available for every earlier study. Missing details are not inferred. Local uncommitted case/CLI edits are outside this committed-history scope.

| Evidence category | Meaning in this record |
|---|---|
| Implementation | Functionality present in committed production/application code |
| Unit/regression testing | Focused contracts and known failure modes exercised by tests |
| Numerical verification | Mathematical-reference, refinement or discrete-consistency checks |
| Experimental validation | Comparison with physical measurements; not established exhaustively here |
| Performance measurement | Timings and iteration counts for specified configurations, not complexity guarantees |
| Exploratory prototype | Temporary or separate-branch experiment; not an integrated production feature |

Outer SIMPLE convergence, spatial convergence and agreement with a physical experiment are distinct. A small corrected continuity residual alone does not establish a coupled fixed point. Gains from different campaigns must not be added, and solver-only times must not be confused with complete workflow costs.

## 2. Baseline audit (`53bf506`)

The supplied audit verdict was a sound architecture with no demonstrated critical defect. This is the historical verdict, not a new repository audit.

| Finding | Audit concern | Position at the end of this period |
|---|---|---|
| F1 | FixedMassFlux and anisotropic grad(p') reconstruction are not strictly equivalent | Approximation remains documented; exact discrete boundary flux correction is separately zero |
| F2 | Insufficient genuinely two-dimensional SIMPLE verification | Rotated Poiseuille and coupled Kovasznay strengthen coverage; neither proves robustness for every flow |
| F3 | Too many application responsibilities in main.cpp | Initialization orchestration has its own application module; main still owns case setup, monitoring and output |
| F4 | Repeated Eigen sparse-structure reconstruction | Fixed-pattern mapping integrated; numerical coefficients and diagonal preconditioners still updated |
| F5 | Repeated traversals for transactional guarantees | Guarantees and associated costs remain; individual hot-path checks improved, not eliminated globally |
| F6 | Partially duplicated Eigen wrapper interfaces/backends | Shared sparse pattern reduces duplication; separate CG/BiCGSTAB wrappers remain |
| F7 | Contracts, documentation and lifetime clarity | Borrowing and level lifetimes made explicit; continuing discipline, not a claim that all risks are removed |

## 3. Development timeline

Commits are ordered chronologically and grouped by feature. All 24 commits after the audit are represented.

| Commit(s) | Milestone | Evidence/status |
|---|---|---|
| `b7c0855`, `6622848` | Face-flux fixed-point convergence with alpha_rc relaxation; CSV regression adjustment | Implementation and regression testing |
| `12014b7` | Rotated Poiseuille | Numerical verification of coupled velocity components |
| `b48b27b` | Multi-mesh SIMPLE scaling benchmark | Performance measurement infrastructure |
| `0400e95` | Coupled Kovasznay Re=40 | Two-dimensional analytical/refinement verification |
| `33135a7`, `7ad97ed`, `9dc2e63` | Backward-facing-step geometry, automatic meshing and Armaly study | Implementation, meshing tests and numerical study; not exhaustive experimental validation |
| `9090d94`, `162e6e6` | Algebraic momentum relaxation and Majumdar interpolation; Armaly driver update | Implementation, regression tests and coupled study |
| `f605af0`, `103deb2` | Hybrid convection and explicitly requested face-branch diagnostics | Implementation/testing; diagnostics add no counters inside SIMPLE assembly |
| `a085824` | LinearUpwind deferred correction | Implementation, tests and numerical verification |
| `a3a5460` | Barth-Jespersen limiter | Implementation and focused bounded-reconstruction/aliasing tests |
| `fedec03` | Steady scalar transport | Implementation and numerical verification with supplied carrier flux |
| `589b77f`, `9507fce` | Relax CLI pressure and momentum linear tolerances | Case-specific performance/solution studies preceding the committed settings |
| `c33c2e5` | Exactly zero BiCGSTAB RHS reporting | Regression fix; zero iterations/error with an exact zero solution |
| `91a2d2a` | Fixed Eigen sparse patterns | Prototype evidence followed by production integration and equivalence checks |
| `fd4b3f9` | CellFieldTransfer | Implementation, geometric coverage/transfer tests and Armaly mapping smoke |
| `4f2d326` | Coarse-grid SIMPLE initialization | Application integration after grid-sequencing experiments |
| `a7df3a0` | BoundaryBasedInitialization | Bulk initial-guess heuristic and focused tests |
| `df56f2e` | Multilevel grid sequencing | Application orchestration, explicit completion reports and failure tests |
| `75dc51a` | Export converged initialization levels | Callback-based VTU output, artifact cleanup and persistence tests |

CompactQ2Local, ILUT, pressure warm-start and pressure IncompleteCholesky experiments must not be read as additional committed solver capabilities.

## 4. Mathematical and numerical verification

| Study | Retained result | Scope/interpretation |
|---|---|---|
| Rotated Poiseuille | Velocity convergence approximately order 2; transverse velocity near machine precision in studied configurations | Analytical verification, not a general oblique-boundary closure proof |
| Kovasznay Re=40 | Velocity RMS near order 2; pressure RMS tending toward order 2; pressure Linf near order 1 | No second-order pressure Linf claim |
| Armaly Re=100, historical QUAD series | 3,955 / 14,888 / 57,694 cells at h=0.10 / 0.05 / 0.025; reported xR/S 2.11402134 / 2.44205295 / 2.64078033 | Refinement trend; not exhaustive comparison with experimental Armaly data |
| Convection/transport | Hybrid, LinearUpwind deferred correction and Barth-Jespersen implemented/tested; steady scalar transport verified | No universal accuracy or performance gain inferred from scheme availability |

The initial Kovasznay N=64 BiCGSTAB breakdown was investigated rather than hidden by a looser linear tolerance. Temporary SparseLU momentum solves did not stabilize SIMPLE: the linear failure was a symptom of unstable nonlinear iteration. Lower provisional-flux relaxation allowed a common N=8/16/32/64 campaign with Linear, WLS, alpha_u=1, alpha_p=0.1 and alpha_rc=0.1, zero cell initial fields and unchanged analytical BCs. The committed driver retains strict outer tolerances `1e-10` and linear tolerances `1e-12`; temporary diagnostic backends were not retained.

Kovasznay acceptance checks require strictly decreasing RMS and Linf errors for u/v/p, final RMS orders above 1.5, and boundary flux error/imbalance within `1024 * epsilon`. They deliberately impose no second-order Linf(p) requirement. Gradient/error studies use area-weighted, domain-normalized RMS measures; Linf uses the absolute signed error.

The Armaly coarse/medium/fine series is systematic. A locally refined fine_D sensitivity mesh is not a fourth GCI level. Solver convergence, grid uncertainty and experimental uncertainty must remain separate.

## 5. Performance investigations

The Armaly tolerance checks used FirstOrderUpwind, alpha_u=0.7, alpha_p=0.3, alpha_rc=1, outer criteria `1e-8`, maximum outer 4,000 and maximum linear 5,000. The pressure sweep kept momentum at `1e-10`; the subsequent momentum sweep kept pressure at `1e-3`. Earlier benchmark metadata are not complete enough here to reproduce every timing.

| Investigation | Recorded observation | Decision/limit |
|---|---|---|
| SIMPLE scaling | Empirical total-time exponents approximately 1.34–1.48 versus cell count; Krylov work grows with refinement | Specific benchmark configurations, not an asymptotic complexity bound |
| Momentum ILUT, integrated experiment | Total SIMPLE time increased by 45.6%, 33.1%, 18.7% at h=0.20, 0.14, 0.10 despite much fewer BiCGSTAB iterations | Same outer counts/physical solutions; retain Diagonal. Default ILUT parameters, no tuning campaign |
| Pressure tolerance on Armaly fine | CG total 3,613,888 -> 2,943,466 (-18.55%); 1,182 outer in both runs; indicative total-time gain about 19% | Supports CLI pressure CG `1e-3`, after Kovasznay and Armaly coarse/medium/fine checks |
| Momentum tolerance on Armaly fine | Historical `1e-10` -> `1e-6`: 1,182 -> 679 outer; indicative total time -43.11%; xR/S relative change -0.001446% | Fields/flow practically unchanged in these checks; supports CLI momentum `1e-6`, not arbitrary tolerance relaxation |
| Pressure warm-start | Exploratory test, not adopted | No new production initial-guess policy; retained numerical/timing details insufficient here for a universal conclusion |
| Pressure IncompleteCholesky | Exploratory test, not adopted | Setup plus solve, not Krylov count alone, governs the choice; no general rejection of IC on other systems |
| Fixed sparse pattern | Prototype momentum/pressure preparation reductions about 77%/72%; integrated N64 total gain indicative about 16% | Fixed topology enables direct coefficient-slot updates; do not add these percentages to other campaign gains |
| Precomputed CompactQ2Local | About 0.733 times WLS reconstruction cost at approximately one million cells; about 64 MB per million cells | Separate exploratory work, absent from this HEAD; no demonstrated gain on a complete CFD solve |

Earlier Kovasznay N64 profiling attributed about 56.8% to pressure solving with tight linear settings. After relaxing pressure CG to `1e-3` while momentum was still tight, momentum u+v accounted for about 68.55%, pressure about 12.45%. These fractions belong to different historical settings, not one additive speedup model. They motivated measuring momentum tolerance next.

The CLI changes were pressure `1e-12` -> `1e-3` and momentum `1e-12` -> `1e-6`. Generic CG/BiCGSTAB option defaults remain `1e-10`; historical verification drivers may explicitly request tighter values. Linear tolerances are not interchangeable with outer SIMPLE tolerances. Less work is accepted only with external convergence, field and physical-observable checks.

The zero-RHS fix treats exact zero, including negative zero, without an epsilon. Eigen could report Success with an obsolete iteration counter; the wrapper now returns a zero solution, zero iterations and zero estimated error. Finite-RHS validation and exact-zero detection share one traversal. Nonzero-RHS numerics were not redesigned.

Sparse-pattern reuse retains default Eigen sparse storage indices, explicit zero slots and validation. SIMPLE shares the immutable mapping for u/v/p; numerical matrices and preconditioners remain independent. Coefficients and diagonal preconditioners are recomputed as needed: this is not frozen-matrix or preconditioner reuse.

## 6. Grid sequencing experiments

CellFieldTransfer locates target centers geometrically in independent TRI/QUAD source meshes. PiecewiseConstant and unlimited LinearReconstruction are available; the latter uses source WLS gradients and source BCs. Neither method guarantees integral conservation. No face flux, p', matrix or system is transferred. Every receiving mesh gets its own historical boundary-flux initialization.

| Historical experiment | Retained observation | Cost interpretation |
|---|---|---|
| Armaly coarse -> medium, Linear | Workflow -27.28% versus medium zero-start | Includes source solve, mapping and transfer, not just the warm medium solve |
| Armaly medium -> fine, Linear | Fine zero: 2,064.285 s / 679 outer; fine warm: 1,234.355 s / 467 outer | Complete medium-plus-transfer-plus-fine workflow -36.97%; do not compare these timings with the later 1e-5 campaign |
| Kovasznay N32 -> N64, Linear | Complete workflow +15.25% | Source solve overhead can outweigh target savings; no universal benefit |

The application decision was coarseMesh by default, with type zero as an explicit opt-out and no performance-based size threshold. The default sizes are h0=4hf, h1=2hf, h2=hf. Optional targetCoarseCellCount changes the coarsest target; Gmsh counts are not guaranteed exactly.

Gmsh constructs the meshes independently. This is grid sequencing for initial guesses, not AMR, a parent-child hierarchy, a multigrid pressure solver or a restart mechanism. Zero mode preserves the caller's internalField, which need not be zero. Automatic hierarchy/coverage failure falls back to that state; explicit targets make those failures errors. Actual SIMPLE/numerical failures are not masked.

## 7. Boundary-based initialization experiments

V1 writes a uniform bulk cell velocity from face-length-weighted admissible inflow velocities, and uniform pressure from Dirichlet pressure faces. Eligible inflow requires both velocity components Dirichlet and an inward flux exceeding the documented roundoff-scaled threshold. Missing eligible data uses the respective internalField fallback. This is not a solved quasi-1D flow model and does not change target fluxes.

| Isolated trial | Outer iterations | Interpretation |
|---|---|---|
| Armaly coarse | 97 -> 96 | Small first-solve benefit in that configuration |
| Kovasznay N32 | 1,582 -> 1,582 | No iteration reduction observed |

The complete tolerance/CPU settings of these isolated trials are not retained here; their counts must not be compared directly with the later complete campaign.

### Controlled complete Armaly campaign

Re=100, QUAD 3,955 -> 14,888 -> 57,694 cells; FirstOrderUpwind and WLS. rho=1, mu=0.0102970297029703, mean inlet velocity=1; historical integrated parabolic inlet, no-slip walls and p=0 outlet. alpha_u=0.7, alpha_p=0.3, alpha_rc=1; all three outer criteria `1e-5`, maximum outer 4,000. Momentum `1e-6`, pressure `1e-3`, both maximum linear 5,000.

Geometry: upstream length `20/1.01`, downstream length 25, downstream height 1, inlet height `0.52/1.01`. Nominal h=0.10/0.05/0.025; automatic meshing, maximum growth 1.2, wall/step refinement factors 1, wall layers 4 and step layers 6.

All meshes/mappings were prepared once and reused. Order D -> A -> B -> C; CPU2, performance governor, turbo disabled, one measurement per solve, no throttle-counter increase observed. No old timing was used as the principal baseline. Totals below are conceptual autonomous costs including required auxiliary preprocessing, solves and transfers; common final-mesh generation is excluded.

| Workflow | Coarse outer | Medium outer | Final outer | Complete cost (s) | Gain versus D |
|---|---:|---:|---:|---:|---:|
| A: internalField -> coarse -> medium -> fine | 84 | 144 | 359 | 882.200 | 30.11% |
| B: BoundaryBased -> coarse -> medium -> fine | 83 | 144 | 359 | 906.263 | 28.20% |
| C: BoundaryBased -> medium -> fine | — | 202 | 359 | 919.779 | 27.13% |
| D: internalField -> fine | — | — | 549 | 1,262.214 | reference |

| Comparison | Measurement | Interpretation |
|---|---|---|
| A -> B | Coarse 84 -> 83; medium/final outer unchanged; first coarse provisional continuity reduced by about 91.7 times | Better first start, practically erased after the coarse solve converges |
| A -> B timings | B +24.063 s (+2.73%), with virtually identical later Krylov work | One sample cannot attribute the slowdown causally to BoundaryBased |
| C -> B | Medium 202 -> 144; complete workflow saves 13.516 s (1.47%) | Extra coarse modestly amortized in the observed run; small time margin |
| A/B/C -> D | Final 549 -> 359 outer; complete workflow saves 27–30% | Sequencing beneficial for this Armaly fine configuration |

Final xR/S was about 2.64074035 for A/B, 2.64074030 for C and 2.64058324 for D. Relative flow mismatch was at most about `1.01e-10` on the final levels. Against D, A/B/C maximum differences were approximately `1.58e-3` in u (0.105% of the reference maximum), `6.40e-4` in v (0.333%), `1.49e-4` in p (0.00142%) and `3.86e-5` in face flux (0.105%). These are not machine-identical fields.

D differed from the strict historical fine xR/S=2.64078033 by about -0.00019709 (-0.00746%), only 0.099% of the historical medium-to-fine difference 0.19872738. The `1e-5` criterion was sufficient for this performance study, not for validating fine field accuracy; no archived strict-field comparison was available. All four starts had initial flux continuity residual 1 because internal/outlet target fluxes were initially zero. The transferred fields do not transfer continuity or flux history.

Conclusions are restricted to this case/configuration: grid sequencing helps Armaly fine; BoundaryBased mainly helps the first start; the additional coarse level has a modest benefit. They do not establish an optimal universal initialization policy.

## 8. Architectural decisions and limitations

| Decision | Rationale and remaining limit |
|---|---|
| Single-Mesh numerical solvers; sequencing in src/app | Reuses established assembly/solve components without teaching SIMPLE about hierarchy or file output |
| Stable, explicitly owned auxiliary states | Borrowed Mesh addresses stay valid; only source and receiving target auxiliary states coexist during transfer, in addition to the already existing final mesh |
| Release resources in lifecycle order | Auxiliary solver/matrices die before transfer; mapping/workspace die before source; source dies before the next solve. Reports retain metadata, not all meshes/fields |
| Explicit optional completion count | nullopt means no successful auxiliary solve; a present zero would still be a completed solve, not a sentinel |
| Synchronous converged-level observer | solve -> convergence check -> report -> export -> transfer; caller writes VTU while borrowed data is alive, without copying u/v/p or retaining references |
| Persist completed-level VTUs | Later level/final failure does not remove levels already exported in the current run; unconverged levels are not exported |
| Remove stale run artifacts once before sequencing | Remove results/initialization/, solution.vtu and convergence.csv, not the whole results directory; propagate filesystem errors rather than fallback |
| Diagnostic VTU, not restart | Mesh, u/v/p, U and mass imbalance are visible; no solver internals, p', matrices or complete restart state |
| Preserve FixedMassFlux closure distinction | Discrete F'_b=0 is exact; auxiliary scalar grad(p') . n=0 is not generally equivalent to (D_P grad(p')) . S_b=0 on oblique faces with anisotropic response |
| Minimal numerical abstractions | No retained generic preconditioner strategy or dynamic solver hierarchy solely for rejected experiments; shared sparse mapping has a measured purpose |

The current model remains steady, 2D and constant-property incompressible flow on one connected cell domain. Case BC data are uniform per group, not arbitrary spatial functions. Linear is unbounded; limiting LinearUpwind reconstructions is not a proof of global coupled-solution boundedness. Scalar transport requires a caller-supplied well-posed problem and is not automatically coupled to SIMPLE. BoundaryBased is a bulk heuristic; transfer is nonconservative; numerical failures need explicit handling rather than a successful-looking fallback.

## 9. Outstanding questions

- How robust is the FixedMassFlux reconstruction approximation on strongly oblique, anisotropic-response stencils?
- Which parts of genuinely 2D/nonlinear and physical experimental validation remain insufficient beyond the current analytical and Armaly studies?
- When do auxiliary solves pay for themselves across mesh sizes, flows and schemes? The Kovasznay counterexample and single-measurement margins prevent a universal size/policy claim.
- How much initialization trajectory dependence is acceptable for a chosen accuracy study? Outer tolerance is not a field-error bound.
- Do the remaining transactional traversals, backend duplication and application responsibilities warrant changes after measurement, rather than on stylistic grounds?
- Does exploratory CompactQ2Local help a complete CFD solver under controlled repeated timings? A reconstruction microbenchmark alone does not answer this.

## 10. Next development phase

The intended direction is axisymmetric finite-volume metrics, then compressible Euler with conservative variables and a cold-nozzle case. Energy, viscous extensions, species, turbulence, thermochemistry and combustion follow as later phases. These steps are not implemented or validated at this HEAD. Existing incompressible verification and diagnostic contracts are a foundation, not evidence that those new equations, metrics or boundary treatments will work without their own verification.
