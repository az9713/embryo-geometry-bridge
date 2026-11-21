# Numerical integration of the actual activated P2 model

**Completed:** all six cases were re-equilibrated at higher quadrature and pass a still higher-rule force check. Eight synthetic checks and 68 independent replay checks pass. The largest shortening correction is 0.019584 µm; the largest rigorously bounded apical-curve change is approximately 0.025714 µm.

This is a target-free integration-accuracy study on the **same 12×6 finite-dimensional P2 space**, not a new constitutive model. The frozen P1, prestress, mixed, full-pointwise P2 and robust-controller sources are preserved. The separate implementation is `scripts/quadratic_quadrature_study.py`; source hashes, original-state hashes and actual integration policies are recorded in `evidence/activated_quadrature_*.json`.

## 1. Why a solved field must be tested under higher quadrature

The default six-point degree-four triangle rule integrates some polynomial terms exactly, but the hyperelastic potential contains \(\log J\), \((\log J)^2\), and the surface terms contain curved-edge lengths. These are not generally degree-four polynomials of the reference coordinates. Stationarity of one quadrature approximation does not imply stationarity of a more accurately integrated potential on the same displacement space.

For a fixed finite-element coefficient vector \(u_h\), define the actual rule-dependent potential and residual

\[
 E_Q(u_h)=\sum_{T,q}w_{Tq}W(F_h(x_{Tq})) + E_{s,Q}(u_h)+E_{a,Q}(u_h),
 \qquad R_Q=\nabla_{u_h}\frac{E_Q}{\mu V_0}.
\]

The original check was \(\|R_{Q_0}\|_\infty<10^{-7}\). Here we first hold \(u_h\) fixed and evaluate \(R_{Q_8}\) and \(R_{Q_{12}}\). A larger residual is integration sensitivity, even if the total energies differ by only a small relative amount. A constant prestress reference energy dominates the absolute total in these runs; agreement of total energy alone is especially weak evidence of force accuracy.

## 2. Declared higher-order rules

On the unit triangle let

\[
 (\lambda_0,\lambda_1,\lambda_2)=(1-u,u(1-v),uv),\quad 0<u,v<1.
\]

The Jacobian of \((u,v)\mapsto(\xi,\eta)=(u(1-v),uv)\) is \(u\). We use tensor Gauss–Legendre integration with 8×8 and 12×12 interior points and the physical axisymmetric weight \(2\pi r\det D_m\,u\,du\,dv\). Surface rules use 12 and 16 Gauss points on each exact quadratic edge, respectively.

Every cache is rebuilt explicitly: P2 basis values and gradients, reference radii, volume weights, reference volume, pressure mass matrix and diagnostic projection basis, and quadratic-edge basis/derivative/surface weights. No global quadrature function is monkeypatched. The inherited parameter fields still describe the frozen nominal construction; every new saved result separately records its **actual** integration policy.

The constitutive potential remains

\[
 W(F)=\frac\mu2(I_1-3-2\log J)+\frac\lambda2(\log J)^2,
 \quad E_a=\int\frac f2\chi_A|F e_z|^2dV_0,
 \quad E_s=\int\left[\frac{K_s}{2}(J_s-1)^2+T_sJ_s\right]dA_0.
\]

Continuous P1 projected pressure is auxiliary diagnostic output only. It does not enter the physical potential, so no hidden projected-volume cancellation changes the derivatives. Element chunks of size 12 bound vectorized-Hessian memory without changing the summed potential; chunk sizes one and twelve are separately compared in verification.

## 3. Existing equilibria fail the higher-rule force check

The actual geometry and material settings match the robust small-region runs: height scale 33.5 µm, measured diameter 18.71024074 µm, activation 5–25 µm below the apical surface, \(K_b=2\), \(K_a=0.2\), \(T_b=2\), \(T_a=0\) in the declared normalized units. No gain or material fit is performed.

| \(\lambda/\mu\) | \(f/\mu\) | Higher-rule force residual, Duffy8/surface12 | Duffy12/surface16 |
|---:|---:|---:|---:|
| 500 | 1 | 3.30997e−5 | 3.30997e−5 |
| 500 | 5 | 1.30289e−4 | 1.30289e−4 |
| 500 | 10 | 1.98022e−4 | 1.98022e−4 |
| 5000 | 1 | 1.73241e−4 | 1.73241e−4 |
| 5000 | 5 | 6.73553e−4 | 6.73553e−4 |
| 5000 | 10 | 9.94272e−4 | 9.94272e−4 |

All six default equilibria fail the original \(10^{-7}\) force threshold under both higher rules. This does not invalidate the earlier statement that they solve the original quadrature equations. It shows that the default quadrature equations are insufficiently stationary approximations to the higher-rule equations at this strict force tolerance.

The two higher rules agree closely on the fixed fields, while their bulk energies differ from the default rule. Surface terms agree essentially to rounding in this particular study. Full component and residual receipts are in `evidence/activated_quadrature_evaluation.json`.

## 4. Re-equilibration and verification

Higher-rule fixed-load relaxations are retained separately under `runs/activated_quadrature`. An original-rule field is only an **initial guess** under the new integration policy. It is never passed off as a verified previous equilibrium for that new policy. The procedure is same-strength numerical relaxation, not a physical activation trajectory.

Before re-equilibration, eight synthetic checks pass: zero-reference force, exact affine isochoric energy, positive integration weights, reference volume, positive pressure mass matrix, automatic derivative versus directional finite difference, and chunk-size invariance of both potential and gradient. The directional derivative relative discrepancy is approximately 2.33e−10. These are component checks, not continuum error estimates.

All six new Duffy8/surface12 equilibria pass the Duffy12/surface16 force check. The same new field is evaluated under that second rule; an independently optimized second field is unnecessary for this declared residual gate and has not been claimed.

| \(\lambda/\mu\) | \(f/\mu\) | New-rule force residual | Duffy12 cross residual | Change in maximum shortening (µm) | Apical curve change bound (µm) |
|---:|---:|---:|---:|---:|---:|
| 500 | 1 | 3.01e−13 | 7.68e−12 | −0.000435 | 0.001739 |
| 500 | 5 | 4.75e−13 | 3.61e−11 | +0.005188 | 0.007102 |
| 500 | 10 | 2.79e−13 | 3.21e−11 | +0.019584 | 0.025714 |
| 5000 | 1 | 3.25e−12 | 2.50e−11 | −0.004122 | 0.005230 |
| 5000 | 5 | 3.14e−12 | 6.65e−11 | +0.000837 | 0.011921 |
| 5000 | 10 | 2.81e−12 | 5.16e−11 | +0.006522 | 0.014481 |

Curve-bound displays are rounded upward here. The exact bounds are stored as rational squared lengths; they compare the old and new **discrete curves**, not either curve to a continuum solution. The maximum basal-curve bound over all cases is less than 0.003130 µm.

![Quadrature residuals, shortening corrections and exact curve-difference bounds](../runs/activated_quadrature/quadrature-overview.png)

## 5. Independent isolation and exact geometry receipts

The post-run reviewer independently evaluates four policies on each original field: volume/surface orders (8,8), (8,12), (12,12), (12,16). Isolating surface order 8→12 changes the force vector by at most 9.97e−14; isolating volume order 8→12 changes it by at most 6.85e−11. These are tiny compared with the original default-to-higher residuals. Thus the observed failure is attributable to the default volume integration in this tested setting, rather than the surface rule.

All 864 new quadratic triangles pass exact rational positive-J certificates. All six bodies pass exact global embedding checks. All source and input-state hashes match. The independent post-run receipt is `evidence/activated_quadrature_verification.json`; full exact geometry coefficients are in `evidence/activated_quadrature_reequilibration.json`.

For boundary differences, let \(d_0,d_m,d_1\) be exact rational differences of old and new quadratic nodal positions on the same material edge. The Bézier difference controls are

\[
 c_0=d_0,\quad c_1=2d_m-(d_0+d_1)/2,\quad c_2=d_1.
\]

Nonnegative Bernstein weights sum to one, so
\(\sup_t|y_{new}(t)-y_{old}(t)|\le\max_i|c_i|\).
The maximum squared norm is rational and retained exactly. Pairing corresponding material parameters also bounds the Hausdorff distance between the two curve images. This is a rigorous **same-space discrete shape-difference bound**. It is not a bound on the vertical Eulerian graph difference, equilibrium error, quadrature remainder, or continuum uncertainty.

![Higher-rule and original actual profiles, with their small difference resolved separately](../runs/activated_quadrature/quadrature-profiles.png)

The plots derive exclusively from retained fields and receipts. `runs/activated_quadrature/make_visual.py` reproduces them; `evidence/activated_quadrature_visuals.json` binds the figure and generator hashes.

## 6. Interpretation and remaining accuracy requirements

The default rule's fields were excellent solutions of their own algebraic equations, yet failed a higher-rule force test. After re-equilibration the same P2 space passes both higher rules, with a small but resolved geometric correction. This illustrates why force, geometry and integration checks must be reported together. A small total-energy change alone can conceal a meaningful force discrepancy.

Agreement of Duffy8 and Duffy12 is empirical integration evidence, not a rigorous quadrature-remainder certificate. Neither it nor the exact discrete geometry certificates establish continuum accuracy, uniform penalty robustness, a global energy minimum or biological validity. The separate mesh study intentionally holds the default rule fixed, so its mesh effect must be kept distinct from this same-space integration effect. A final combined higher-rule/finer-mesh check remains appropriate before a physical prediction freeze.

Reproduce with:

```powershell
python -X utf8 controlled_furrowing/scripts/quadratic_quadrature_study.py --self-test
python -X utf8 controlled_furrowing/scripts/quadratic_quadrature_study.py --evaluate
python -X utf8 controlled_furrowing/scripts/quadratic_quadrature_study.py --solve
python -X utf8 controlled_furrowing/runs/activated_quadrature/make_visual.py
```

The `--solve` command re-runs synthetic verification and original-field evaluation before the six fixed-load solves. Original robust fields are preserved. No biological labels, response outcomes or fitted gains are read.
