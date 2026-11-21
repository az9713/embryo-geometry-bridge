# Restoring pointwise volume control and testing a forced continuum solution

The projected-pressure experiment exposed an important failure: small resolved pressure residuals and small total volume change can coexist with large local volume changes. The separate P2 penalty model restores the original pointwise constitutive energy while retaining quadratic displacement, full finite-height coupling, true surface geometry and optional baseline tensions. Neither earlier solver is modified.

`scripts/axisymmetric_quadratic.py` is this separate model. It makes no claim that P2 penalty elements are locking-free or that one benchmark certifies biological predictions. Its forced manufactured experiment goes beyond interpolation: the solver is given a known body force and approximate Dirichlet data, and its computed equilibrium is compared with an explicitly known, coefficient-independent nonlinear continuum solution.

## The missing volume modes

For the positive quadrature inner product used in [the mixed formulation](MIXED-MECHANICS.md), set \(\ell=\log J\) and let \(\Pi_h\ell\) be its continuous P1 projection. Orthogonality gives

\[
 \frac{\lambda}{2}\|\ell\|_Q^2
 =\frac{\lambda}{2}\|\Pi_h\ell\|_Q^2
 +\frac{\lambda}{2}\|\ell-\Pi_h\ell\|_Q^2.
\]

The projected model retained only the first term. The new module adds the nonnegative second term, thereby restoring exactly the full pointwise penalty energy at the same quadrature points:

\[
 E_v=\int_{\Omega_0}\frac{\lambda}{2}(\log J)^2\,dV_0.
\]

This is a P2 displacement penalty method, not a claim that a continuous P1 pressure constraint controls every local volume mode. The governing volumetric stress coefficient is now \(p_{\rm vol}=\lambda\log J\) at quadrature points. The projected pressure coefficients remain separately labeled auxiliary diagnostics. As before, the inherited parameter name `bulk` denotes \(\lambda\), not the infinitesimal bulk modulus \(\lambda+2\mu/3\).

The numerical identity between the added kernel energy and the explicitly integrated projection residual agrees to \(1.08\times10^{-19}\) in the synthetic check. The full objective's AD/finite-difference derivative relative error is \(3.87\times10^{-8}\). Zero-dose and activated reference checks pass. Protocol, source/parent hashes and values are retained in `evidence/quadratic_protocol.json` and `quadratic_verification.json`.

## A nonlinear manufactured equilibrium with exactly J=1

On the reference cylinder \(0\le r\le R_o\), \(0\le z\le H\), prescribe

\[
 R(r,z)=r e^{\epsilon z},\qquad
 Z(r,z)=\frac{1-e^{-2\epsilon z}}{2\epsilon},\qquad\epsilon\ne0.
\]

Write \(A=e^{\epsilon z}\), \(B=e^{-2\epsilon z}\). Its full axisymmetric deformation gradient is

\[
 F=\begin{pmatrix}A&0&\epsilon rA\\0&A&0\\0&0&B\end{pmatrix},\qquad J=A^2B=1.
\]

The volumetric stress vanishes for every \(\lambda>0\). With no activation and no surface energies, the exact first Piola stress is \(P=\mu(F-F^{-T})\), with nonzero cylindrical components

\[
 P_{rr}=P_{\theta\theta}=\mu(A-A^{-1}),\quad
 P_{rz}=\mu\epsilon rA,\quad
 P_{zr}=\mu\epsilon r/B,\quad
 P_{zz}=\mu(B-B^{-1}).
\]

The reference cylindrical divergence has components

\[
 (\operatorname{Div}P)_r=\partial_rP_{rr}+\partial_zP_{rz}
                          +(P_{rr}-P_{\theta\theta})/r,
 \qquad
 (\operatorname{Div}P)_z=\partial_rP_{zr}+\partial_zP_{zz}+P_{zr}/r.
\]

Substitution gives the independently derived body force

\[
 f_r=-\mu\epsilon^2r e^{\epsilon z},\qquad
 f_z=2\mu\epsilon e^{-2\epsilon z},\qquad
 -\operatorname{Div}P=f.
\]

In the vertical equation the two \(\mu\epsilon/B\) terms cancel the \(-2\mu\epsilon/B\) term. Omitting the cylindrical \(P_{zr}/r\) contribution would give an incorrect force.

The solver includes the conservative force potential \(-\int_{\Omega_0}f\cdot y\,dV_0\). Its first variation is \(-\int f\cdot v\,dV_0\), so the weak equilibrium is \(\int P:\delta F[v]\,dV_0=\int f\cdot v\,dV_0\), with Dirichlet boundary variations zero. Every boundary displacement node is prescribed by the exact map. Along an edge the numerical Dirichlet curve is its P2 interpolant, so boundary approximation error remains part of the measured discretization error.

Crucially this is the same exact continuum solution and the same force for all tested volume coefficients. Changes with \(\lambda\) cannot be explained by a changing manufactured continuum deformation. They expose numerical coefficient sensitivity in this benchmark. Stationarity of the exact continuum solution is proved; global uniqueness or global minimality is not asserted.

## Predeclared numerical test

The retained protocol `evidence/quadratic_manufactured_protocol.json` fixes \(\epsilon=0.2\), \(\mu=H=1\), \(R_o=3\), meshes 6×3, 8×4 and 12×6, and \(\lambda\in\{5,50,500\}\). Activation, area elasticity and baseline surface tensions are all zero. The starting displacement is the analytic P2 nodal interpolant. Thus the test has no biological data, inferred constitutive parameter or post hoc fit.

For numerical displacement \(y_h\), errors use the reference volume weight \(2\pi r\):

\[
 e_{L^2}=\left(V_0^{-1}\int|y_h-y_*|^2dV_0\right)^{1/2},
\]

\[
 e_{H^1}=\left(V_0^{-1}\int\left[|F_{m,h}-F_{m,*}|^2
              +|R_h/r-e^{\epsilon z}|^2\right]dV_0\right)^{1/2}.
\]

The second expression includes the hoop gradient, giving the complete axisymmetric three-dimensional displacement-gradient error. The records also report RMS \(J-1\) and maximum sampled \(|J-1|\). Errors are evaluated by the stated degree-4 quadrature, not rigorously enclosed continuous norms.

All nine solves passed the normalized nodal force threshold \(10^{-7}\) and the built-in sufficient positive-orientation test. The H1 RMS errors are:

| Mesh | \(\lambda=5\) | \(\lambda=50\) | \(\lambda=500\) |
|---|---:|---:|---:|
| 6×3 | 0.002219 | 0.002907 | 0.008260 |
| 8×4 | 0.0006513 | 0.0007484 | 0.001760 |
| 12×6 | 0.0002821 | 0.0003134 | 0.0008006 |

The displacement L2 RMS errors are:

| Mesh | \(\lambda=5\) | \(\lambda=50\) | \(\lambda=500\) |
|---|---:|---:|---:|
| 6×3 | 0.0001768 | 0.0002157 | 0.0005859 |
| 8×4 | 0.00002427 | 0.00002882 | 0.00008493 |
| 12×6 | 0.000006954 | 0.000007990 | 0.00002489 |

Refinement improves every tested coefficient, but higher volume penalty still raises displacement/gradient error. This is evidence of coefficient sensitivity consistent with penalty locking; it is not a uniform-in-\(\lambda\) accuracy proof. Fine-mesh \(\lambda=500\) sampled J lies approximately between 0.998272 and 1.000305; its conservative Bernstein lower bound is about 0.956732. The latter is not an attained compression and must not be plotted as such.

Complete initial-interpolation and solved errors, forces, bounds, runtimes, fields, protocol and source hash are in `evidence/quadratic_manufactured.json` and `runs/quadratic_checks/manufactured/`. Numerical outputs remain reproducible with:

```powershell
python -B controlled_furrowing/scripts/axisymmetric_quadratic.py
python -B controlled_furrowing/scripts/axisymmetric_quadratic.py --manufactured
```

## What this establishes and what remains

The test verifies a solved nonlinear equilibrium against a known continuum solution, independently of the larger-region experiment. It distinguishes implementation error, finite-element approximation and volume-penalty sensitivity more effectively than checking conservation or replaying one optimized figure. These methods and the manufactured solution construction are standard numerical verification, not a new biological mechanism or a novelty claim.

It does not establish an a posteriori continuum error bound, asymptotic convergence order, locking-free behavior, realistic material parameters or validated transient bending. Quadrature convergence and independent interval-enclosed residual/error estimates would be stronger numerical receipts. Biological prediction still requires the separate small-only calibration, geometry/constitutive uncertainty, independently defined observation/time map and withheld larger-region scoring.

![Forced manufactured accuracy and remaining penalty sensitivity](../runs/quadratic_checks/quadratic_accuracy.png)

The visual compares solved equilibria with the same exact continuum map, rather than comparing only interpolated functions. Figure/input hashes are retained in `evidence/quadratic_visual.json`.
