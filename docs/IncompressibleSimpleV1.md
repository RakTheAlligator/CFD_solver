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
11. Evaluate corrected continuity, velocity change, and final iteration diagnostics; test outer convergence and report every successfully completed iteration, including the converged one.

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

The reported x- and y-momentum residuals are external assembled-equation diagnostics evaluated before the corresponding linear solve:

```text
sum_P |b_P - (A U_start)_P|
---------------------------------
sum_P |b_P| + sum_P |(A U_start)_P|
```

These require an additional matrix-vector product per momentum component. They measure the outer momentum-equation imbalance at the iteration-start velocity.

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
```

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
