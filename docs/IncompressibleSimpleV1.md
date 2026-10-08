# Incompressible SIMPLE: current production contract

This document records the current production contract of `IncompressibleSimpleSolver`. It describes implemented behavior only. The historical filename is retained; the contract now includes momentum under-relaxation and Majumdar interpolation.

## Outer iteration

One successfully completed SIMPLE iteration performs the following sequence:

1. Save the iteration-start cell velocity.
2. Reconstruct the cell gradients of `u`, `v`, and physical pressure with inverse-distance-weighted least squares.
3. Assemble and algebraically under-relax the two convection-diffusion momentum systems, then evaluate their external equation residuals at the iteration-start velocity.
4. Prepare and solve the nonsymmetric `u` and `v` systems with Eigen BiCGSTAB, producing the provisional cell velocity `U_star`.
5. Derive the cell momentum pressure response from the assembled momentum diagonals.
6. Interpolate internal and `FixedPressure` boundary mass fluxes with momentum-weighted/Rhie-Chow interpolation and Majumdar correction. Measure the face-flux fixed-point residual before applying the independent provisional face-flux relaxation.
7. Assemble the pressure-correction system and evaluate provisional continuity before any pressure-reference modification.
8. Apply a zero pressure-correction reference only when no `FixedPressure` boundary anchors the system.
9. Solve the pressure-correction system with Eigen conjugate gradient and reconstruct `grad(p')`.
10. Correct face mass flux, cell velocity, and physical pressure.
11. Evaluate corrected continuity, velocity change, and the three fixed-point convergence criteria.
12. On a convergence candidate, reconstruct the current corrected gradients and reassemble both unrelaxed physical momentum equations. Accept convergence only if each component passes its relative or absolute momentum check; report every successfully completed iteration, including rejected candidates and the converged one.

An inner linear-solver failure interrupts the iteration before it is reported as completed.

## Momentum response and corrections

For cell area `A_P` and the assembled, algebraically relaxed component momentum diagonals `a_P,u` and `a_P,v`, the response coefficients are

```text
d_P,u = A_P / a_P,u
d_P,v = A_P / a_P,v
```

They form the diagonal response tensor

```text
D_P = diag(d_P,u, d_P,v)
```

used by the cell-velocity correction:

```text
U_new = U_star - D_P grad(p')
```

Physical pressure is updated with pressure under-relaxation:

```text
p_new = p_old + alpha_p p'
```

The momentum relaxation factor satisfies `0 < alpha_u <= 1`. For each component, algebraic relaxation divides the unrelaxed diagonal by `alpha_u` and adds `(1 - alpha_u) / alpha_u * a_unrelaxed * U_previous` to the RHS. The response uses this relaxed diagonal.

For an internal face, momentum-weighted interpolation adds the Majumdar term

```text
F_Majumdar = F_RhieChow,naive
           + (1 - alpha_u) * (F_previous - rho * I(U_previous) · S_f)
```

`I(U_previous)` uses the face interpolation geometry. FixedPressure boundary faces use the corresponding owner-velocity term. `F_previous` is the preceding corrected mass flux; at `alpha_u == 1`, this correction vanishes and the unrelaxed path is preserved.

The independent factor `alpha_rc` relaxes only computed provisional Rhie-Chow face fluxes:

```text
F_provisional = alpha_rc F_Majumdar + (1 - alpha_rc) F_previous
```

This blend applies to internal and `FixedPressure` faces. It does not change imposed `FixedMassFlux` values or face pressure-response coefficients, and it is not momentum under-relaxation or Majumdar correction.

## Momentum residuals

The existing `x_velocity_equation_residual` and `y_velocity_equation_residual` monitors use the algebraically relaxed systems before the corresponding linear solve:

```text
sum_P |b_P - (A U_start)_P|
---------------------------------
sum_P |b_P| + sum_P |(A U_start)_P|
```

These require an additional matrix-vector product per momentum component. They measure the relaxed-equation imbalance at the iteration-start velocity, not physical momentum convergence; their denominator includes artificial relaxation terms.

The separate `x_momentum_relative_residual` and `y_momentum_relative_residual` convergence checks use the current corrected velocity, pressure, face flux and corresponding freshly reconstructed gradients. Assembly uses `relaxation_factor = 1`, including the existing non-orthogonal diffusion and deferred convection terms, so the residual is `b_physical - A_physical * U_current`.

```text
S_u = sum_P (|b_u,P| + |(A_u u)_P|)
S_v = sum_P (|b_v,P| + |(A_v v)_P|)
R_abs,u = sum_P |b_u,P - (A_u u)_P|
R_abs,v = sum_P |b_v,P - (A_v v)_P|
R_rel,u = R_abs,u / S_u
R_rel,v = R_abs,v / S_v
```

The two scales are independent of algebraic relaxation and of the other momentum component. A common scale could hide a significant relative error in a weak component. With zero component scale, a zero imbalance gives zero relative residual; a nonzero imbalance gives infinity. Every nonzero scale is used as-is: there is no epsilon or dominant-component floor. Non-finite equation values or accumulated imbalances are rejected; a non-finite residual cannot establish convergence.

Each component passes if `R_rel,k <= velocity_relative_tolerance` **or** `R_abs,k <= momentum_absolute_tolerance`. Both components must pass independently. The relative branch reuses the external velocity accuracy setting, not the inner BiCGSTAB tolerance. Absolute and relative diagnostics are reported separately, even when only the absolute branch passes.

The explicit absolute tolerance is finite and non-negative; zero removes the allowance for a nonzero absolute imbalance. Its default is `1e-12`, in the units of the integrated component equations: force per unit out-of-plane depth (`N/m` with SI inputs). `R_abs` is a domain sum of already integrated cell imbalances, not an area-normalized RMS or a per-cell bound. Refinement changes the local integrals and accumulated numerical errors; the sum is not guaranteed mesh-independent or simply proportional to cell count. The tolerance must be adapted to the domain, physical units, force scales and desired accuracy; `1e-12` is not a universal CFD accuracy guarantee.

In the tested Cartesian Poiseuille grids (64/256/1024 cells), the analytically zero transverse component had absolute imbalances about `6.13e-17` to `3.31e-13` at the convergence candidate, despite relative residuals of 0.124 to 0.803. The absolute default admits these negligible imbalances without changing the existing analytical acceptance thresholds. It still rejects the two-scale regression `4v=4e-10`, `v=1.1e-10`, whose absolute imbalance is `4e-11`. For data of amplitude `1e-20`, a similar 10% relative discrepancy gives about `4e-21 N/m` and is deliberately acceptable under the default absolute budget. Set a smaller explicit budget, including zero, if that physical scale requires relative accuracy.

This full check is performed only when the three existing criteria pass. Existing field/system workspace is reused without field-sized allocations. The optional result/callback diagnostics are absent on iterations where no check occurred, rather than retaining a stale residual. Their gradient, assembly and matrix-product costs enter the existing corresponding timing counters. The check does not solve a new system or change relaxation, Rhie-Chow or correction formulas. An inner linear-solver failure remains an error regardless of the absolute allowance.

Eigen separately reports its own estimated relative residual and iteration count for each BiCGSTAB or conjugate-gradient solve. The Eigen values describe the internal iterative linear solve and are not interchangeable with the external momentum residuals.

## Continuity and convergence

After provisional flux assembly, the pressure-correction RHS is the negative provisional cell mass imbalance before gauge fixing:

```text
rhs_P = -R_star,P
```

The provisional relative continuity residual is

```text
sum_P |R_star,P| / sum_f |F_star,f|
```

After pressure, velocity, and face-flux correction, corrected continuity is evaluated as

```text
sum_P |R_P| / sum_f |F_f|
```

If the flux denominator is zero, the relative residual is zero when the corresponding imbalance sum is also zero and infinity otherwise.

When momentum or provisional face-flux relaxation is active, the solver monitors
the fixed-point residual of the momentum-weighted face flux before the alpha_rc blend:

```text
max_f |F_Majumdar,f - F_previous,f|
-----------------------------------
max_f max(|F_Majumdar,f|, |F_previous,f|)
```

The maximum is evaluated over internal and `FixedPressure` faces, i.e. the
faces on which provisional Rhie-Chow flux relaxation is applied.
`FixedMassFlux` boundaries are excluded because their imposed flux is not
relaxed.

With `alpha_u < 1` and `alpha_rc == 1`, this residual remains active. Only the
fully unrelaxed path (`alpha_u == 1` and `alpha_rc == 1`) reports zero for this
diagnostic in the current implementation.

Outer SIMPLE convergence therefore requires:

```text
velocity_relative_change <= velocity_relative_tolerance
rhie_chow_flux_relative_residual <= rhie_chow_flux_relative_tolerance
provisional_continuity_relative_residual <= continuity_relative_tolerance
(x_momentum_relative_residual <= velocity_relative_tolerance
 OR x_momentum_absolute_residual <= momentum_absolute_tolerance)
(y_momentum_relative_residual <= velocity_relative_tolerance
 OR y_momentum_absolute_residual <= momentum_absolute_tolerance)
```

The last two conditions require an actual current-state check. If they fail, iterations continue; reaching the limit returns `converged == false`. Stagnant variables or an inner solver accepting its initial guess are not sufficient evidence of physical momentum convergence.

The face-flux criterion prevents convergence from being reported while the
momentum-weighted flux is still moving toward its fixed point, independently
of whether an additional alpha_rc blend is applied.

Corrected continuity and maximum cell mass imbalance remain diagnostics.
Corrected continuity is not the outer stopping criterion because pressure
correction can make it nearly zero while the provisional state is still far
from the coupled fixed point.

## Pressure-correction boundaries

`FixedPressure` corresponds to prescribed physical pressure and imposes

```text
p'_b = 0
```

Its boundary pressure response contributes to the owner diagonal and its corrected boundary flux is computed from `p'_P`.

`FixedMassFlux` imposes an unmodified boundary mass flux. Its discrete face-flux correction is exactly

```text
F'_b = 0
```

For the separate cell reconstruction of `grad(p')`, SIMPLE maps this condition to the auxiliary scalar condition

```text
grad(p') · n = 0
```

This is only an approximation to the anisotropic response condition

```text
(D_P grad(p')) · S_b = 0
```

on general oblique faces with unequal Cartesian responses. The exact zero `FixedMassFlux` face correction is enforced independently by the face-flux equations and is not recomputed from the corrected cell velocity.

## Current limitations

- Steady, two-dimensional, constant-property incompressible flow only
- One connected fluid-cell component per solve
- Case input accepts uniform scalar boundary values per group, not arbitrary spatial boundary functions
- Auxiliary scalar approximation for `FixedMassFlux` pressure-correction gradient reconstruction
- FirstOrderUpwind, Linear, Hybrid, and LinearUpwind are available; LinearUpwind uses deferred correction with optional Barth-Jespersen limiting, not a general coupled-solution boundedness guarantee
- No transient, turbulence, energy, multiphase, or three-dimensional equations
