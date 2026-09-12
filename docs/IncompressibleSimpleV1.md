# Incompressible SIMPLE v1

This document records the current production contract of `IncompressibleSimpleSolver`. It describes implemented behavior only.

## Outer iteration

One successfully completed SIMPLE iteration performs the following sequence:

1. Save the iteration-start cell velocity.
2. Reconstruct the cell gradients of `u`, `v`, and physical pressure with inverse-distance-weighted least squares.
3. Assemble the two convection-diffusion momentum systems and evaluate their external equation residuals at the iteration-start velocity.
4. Prepare and solve the nonsymmetric `u` and `v` systems with Eigen BiCGSTAB, producing the provisional cell velocity `U_star`.
5. Derive the cell momentum pressure response from the assembled momentum diagonals.
6. Interpolate provisional internal and `FixedPressure` boundary mass fluxes with momentum-weighted/Rhie-Chow interpolation, then apply provisional face-flux relaxation where enabled.
7. Assemble the pressure-correction system and evaluate provisional continuity before any pressure-reference modification.
8. Apply a zero pressure-correction reference only when no `FixedPressure` boundary anchors the system.
9. Solve the pressure-correction system with Eigen conjugate gradient and reconstruct `grad(p')`.
10. Correct face mass flux, cell velocity, and physical pressure.
11. Evaluate corrected continuity, velocity change, the Rhie-Chow flux fixed-point residual when applicable, and final iteration diagnostics; report the completed iteration and test outer convergence.

An inner linear-solver failure interrupts the iteration before it is reported as completed.

## Momentum response and corrections

For cell area `A_P` and the assembled component momentum diagonals `a_P,u` and `a_P,v`, the response coefficients are

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

SIMPLE v1 requires the momentum relaxation factor `alpha_u` to equal exactly `1`. The momentum-weighted interpolation uses the final assembled momentum diagonal and does not implement a relaxation-consistent Majumdar correction.

The independent factor `alpha_rc` relaxes only computed provisional Rhie-Chow face fluxes:

```text
F_provisional = alpha_rc F_RhieChow + (1 - alpha_rc) F_previous
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

When provisional Rhie-Chow face-flux relaxation is enabled, the solver also
monitors the fixed-point residual of the computed face flux:

```text
max_f |F_RhieChow,f - F_previous,f|
-----------------------------------
max_f max(|F_RhieChow,f|, |F_previous,f|)
```

The maximum is evaluated over internal and `FixedPressure` faces, i.e. the
faces on which provisional Rhie-Chow flux relaxation is applied.
`FixedMassFlux` boundaries are excluded because their imposed flux is not
relaxed.

When `alpha_rc == 1`, this diagnostic is defined as zero and does not constrain
convergence.

Outer SIMPLE convergence therefore requires:

```text
velocity_relative_change <= velocity_relative_tolerance
rhie_chow_flux_relative_residual <= rhie_chow_flux_relative_tolerance
provisional_continuity_relative_residual <= continuity_relative_tolerance
```

The Rhie-Chow flux criterion prevents convergence from being reported while a
relaxed divergence-free face flux is still moving toward its fixed point.

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

For the separate cell reconstruction of `grad(p')`, SIMPLE v1 maps this condition to the auxiliary scalar condition

```text
grad(p') · n = 0
```

This is only an approximation to the anisotropic response condition

```text
(D_P grad(p')) · S_b = 0
```

on general oblique faces or when `d_P,u != d_P,v`. The exact zero `FixedMassFlux` face correction is enforced independently by the face-flux equations and is not recomputed from the corrected cell velocity.

## Current limitations

- Steady, two-dimensional, constant-property incompressible flow only
- One connected fluid-cell component per solve
- `alpha_u` fixed to exactly `1`
- No relaxation-consistent Majumdar momentum interpolation
- Uniform scalar boundary data only; spatially varying boundary values are not yet supported.
- Auxiliary scalar approximation for `FixedMassFlux` pressure-correction gradient reconstruction
- First-order upwind or unbounded geometry-aware linear convection only
- No transient, turbulence, energy, multiphase, or three-dimensional equations
