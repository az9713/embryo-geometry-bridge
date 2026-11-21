# Quadratic displacement and a projected logarithmic volume constraint

The earlier P1 displacement model had substantial mesh sensitivity and appreciable local volume change in strongly activated states. This separate component tests a higher-order displacement and mixed-pressure alternative. It does not modify either frozen P1 source. It does not presume that a mixed method fixes locking or guarantees local incompressibility.

`scripts/axisymmetric_mixed.py` implements a finite-height axisymmetric body with shared quadratic displacement nodes and continuous linear pressure on the same material triangulation. This is a genuine coupled solid. Its finite-penalty volumetric energy differs from the earlier pointwise penalty formulation. The constitutive and geometric assumptions remain exploratory, and no experimental outcomes are read by the module or its synthetic checks.

## The mixed problem and its exact condensation

Use the flat reference meridian, three-dimensional axisymmetric deformation gradient, passive shear energy, prescribed active fibers, area elasticity and true baseline tensions described in [the original model](AXISYMMETRIC-MODEL.md) and [surface-tension revision](BASELINE-TENSION.md). Write the logarithmic volume coefficient as \(\lambda\); the inherited API calls this parameter `bulk`. It is a coefficient, not the infinitesimal bulk modulus, which is \(\lambda+2\mu/3\).

Let \(V_h\) be continuous P2 vector displacement with axis/rim constraints and \(Q_h\) continuous P1 scalar pressure. Define the discrete mixed potential

\[
 \mathcal L_h(y_h,p_h)=E_{\rm shear}(y_h)+E_{\rm surface}(y_h)+E_a(y_h)
 +\int_{\Omega_0}\left[p_h\log J_h-\frac{p_h^2}{2\lambda}\right]dV_0,
\]

where \(E_{\rm shear}=\int\mu(I_1-3-2\log J_h)/2\,dV_0\). The pressure stationarity equation is

\[
 \int_{\Omega_0}q_h\left(\log J_h-\frac{p_h}{\lambda}\right)dV_0=0
 \quad(q_h\in Q_h).
\]

For pressure basis \(\psi_i\), let

\[
 M_{ij}=\int_{\Omega_0}\psi_i\psi_j\,dV_0,\qquad
 b_i(y_h)=\int_{\Omega_0}\psi_i\log J_h\,dV_0.
\]

The mass matrix is positive definite: \(c^TMc=\int(\sum_i c_i\psi_i)^2dV_0>0\) for nonzero \(c\), since the axis has zero volume and \(r>0\) almost everywhere in the body. Thus the unique stationary pressure at any admissible displacement is \(p=\lambda M^{-1}b\). Substitution gives the condensed energy

\[
 E_h(y_h)=E_{\rm shear}+E_{\rm surface}+E_a
          +\frac{\lambda}{2}b(y_h)^TM^{-1}b(y_h).
\]

This is elimination of the concave quadratic pressure variable, not minimization in both variables. Automatic differentiation of the condensed energy gives the correct mixed displacement equation. The pressure contribution is \(\int p_h\,\delta\log J_h\,dV_0\), with Piola stress \(p_hF_h^{-T}\). Positive pressure coefficients here mean positive volumetric stress; compression generally produces negative coefficients. This generalized material pressure must not be confused with a compression-positive fluid-pressure convention.

With exact integration, \(p_h=\lambda\Pi_{Q_h}\log J_h\) and the energy is \(\lambda\|\Pi_{Q_h}\log J_h\|_{L^2}^2/2\). Projection controls resolved moments of \(\log J\); it does not control its unresolved spatial oscillations. A small total volume change or a tiny pressure constraint residual does not imply \(J\approx1\) pointwise. No weighted inf-sup estimate, uniform-in-\(\lambda\) error estimate or locking-free theorem has been proved here.

This distinction is structural. If \(\log J_h\) has a component orthogonal to \(Q_h\), increasing \(\lambda\) does not directly penalize that component. The passive shear term still opposes compression: for positive principal stretches \(s_i\), its density is \(\mu\sum_i(s_i^2-1-2\log s_i)/2\ge0\). But this resistance is set by \(\mu\), rather than the large projected-volume coefficient. Consequently neither exact pressure condensation nor positive mass-matrix eigenvalues proves local near-incompressibility. This is an established limitation of weak volume constraints, not a novelty claim.

## Actual spaces and quadrature

Every triangle has the six displacement nodes in order \([v_0,v_1,v_2,e_{01},e_{12},e_{20}]\). Edge midpoints are shared across adjacent triangles, so the field is continuous. Original vertices occupy the first indices; pressure lives on those vertices. The reference element mapping remains affine. Quadratic Lagrange shape functions are

\[
 N_i=\ell_i(2\ell_i-1),\qquad N_{ij}=4\ell_i\ell_j,
\]

where \(\ell_i\) are barycentric coordinates. The hoop stretch \(R/r\), all meridional gradient terms, true surface geometry and through-height deformation are retained. Flat-reference curvature restriction and the clamped outer rim are inherited.

Positive degree-4 Dunavant quadrature is used in the volume. It integrates the reference mass matrix exactly in real arithmetic because the weighted integrand \(r\psi_i\psi_j\) has degree 3. Logarithms and the nonlinear elastic energy are not polynomial; their integrals have quadrature error. The implemented pressure equation is exact for this quadrature model, within floating arithmetic. Surface integrals use 8-point Gauss by default, with a minimum of 5 points.

For this P2 motion, \(\det F_m\) has degree 2 and \(R\) degree 2. Hence

\[
 V=2\pi\sum_T\int_T R\det F_m\,dr\,dz
\]

has degree-4 integrands and is integrated exactly in real arithmetic by the chosen volume quadrature. Its floating implementation is not an interval volume certificate.

## Positivity is checked beyond quadrature samples

The code separately reports sampled \(J\) and sufficient Bernstein coefficient bounds. A conservative lower bound is not an attained minimum; a low bound is not evidence that actual tissue compression reaches that value.

Express \(\det F_m\) in the nonnegative degree-2 Bernstein basis. All-positive coefficients imply positive meridional determinant throughout the triangle. Express the quadratic radius \(R\) and the degree-elevated reference radius \(r\) in the same basis. The coefficients of \(r\) are nonnegative. If every zero denominator coefficient also has exactly zero numerator coefficient, then

\[
 \frac{R}{r}=\sum_{i:r_i^B>0}\frac{r_i^BB_i}{\sum_jr_j^BB_j}\frac{R_i^B}{r_i^B}.
\]

It is therefore bounded by the nonzero coefficient ratios. Axis boundary nodes are constrained so that zero-radius coefficients vanish. On an entire axis edge the quadratic \(R\) is divisible by the linear coordinate \(r\), making the axis limit well defined. Positive determinant and hoop bounds give a sufficient positive \(J\) bound. Failure of the coefficient test is reported as uncertified; it is not called inversion, and sampled positivity does not override it.

The solver's built-in calculation uses floating coefficients. The independent [rational geometry certificate](MIXED-GEOMETRY-PROOF.md), `scripts/mixed_geometry_certificate.py`, derives exact rational bounds from stored binary floating coordinates, independently verifies shared edges and a sufficient embedded boundary, and reports a discrete embedding certificate. It does not turn the computed field into a continuum PDE error certificate. The first zero, activated and scaled synthetic fields all passed that independent certificate.

## Quadratic surfaces and observables

Current top/bottom boundaries are quadratic curves, not node polylines. The module returns every segment's three controls. Its reported profiles include endpoints, midpoints and every interior root of \(dZ/dt=0\), evaluated without a curvature cutoff. Radial monotonicity is checked using the endpoint values of the linear derivative \(dR/dt\).

When the two boundaries are radially monotone, their common Eulerian graph grid includes every boundary Z-extremum. Values on that grid are evaluated by inverting each exact quadratic radial map, using floating bisection. The global height extrema are therefore represented to floating accuracy. Straight interpolation between those reported values is still an approximation to the true P2 graph. The independent rational certificate separately verifies extrema, including arbitrarily small nonzero curvature. The module does not label a sampled graph-thickness minimum as a rigorous minimum; it also supplies the sufficient lower bound \(\min a-\max b\).

Material-column shortening is quadratic on each paired top/bottom reference edge. Its maximum inside the activated patch is computed from endpoints and interior derivative roots, rather than only displacement nodes. It is still a static mechanical observable, not the experimental maximum over time. The source transient and the missing time-to-equilibrium map remain unresolved.

## API and evidence

```python
from axisymmetric_mixed import Parameters, build_model
model = build_model(Parameters(nr=8, nz=4, bulk=50,
                              basal_surface_tension=2,
                              apical_surface_tension=0.2))
first = model.relax(0.2)
second = model.relax(1.0, first['free_displacement'])
```

The warm-start and continuation methods are inherited from the frozen optimizer. `converged` requires verified force stationarity and the built-in sufficient orientation certificate. `stationary_discrete_equations` is reported separately, as are sampled J, conservative bounds, pressure constraint residual and boundary status. Unknown or uncertified geometry is not silently passed.

`evidence/mixed_protocol.json` and `mixed_verification.json` retain eleven synthetic checks: reference equilibrium, mass positivity, isochoric affine J/volume/energy, derivative checks, Bernstein orientation, pressure constraint and uniform material scaling. The activated synthetic force residual is about \(3.6\times10^{-10}\); affine energy relative error \(1.8\times10^{-16}\), AD/finite-difference relative error \(6.8\times10^{-8}\), and pressure constraint residual \(2.7\times10^{-20}\). These quantify implementation consistency, not biological accuracy.

The first successful implementation is retained under `runs/mixed_checks/initial_verified_version/`. Adding complete quadratic output diagnostics exposed an indexing error; `evidence/mixed_initial_failure.json` retains it. The corrected intermediate source is retained under `runs/mixed_checks/interim_profile_version/`. The final output-only correction removes a curvature tolerance from the extrema search; eleven checks were replayed without changing the energy, solver thresholds or either parent source.

## A target-independent nonaffine check

The smooth axisymmetric deformation

\[
 R=r e^{\epsilon z},\qquad
 Z=\frac{1-e^{-2\epsilon z}}{2\epsilon}
\]

has \(\det F_m=e^{-\epsilon z}\), hoop stretch \(e^{\epsilon z}\), and exactly \(J=1\). It is a kinematic check; it need not satisfy the clamped-rim boundary condition and is not an equilibrium test. The exact passive shear energy is

\[
 E_{\rm shear}=\pi\mu\left[
 \frac{R_o^2}{2}(2I_2+I_{-4}-3H)
 +\frac{\epsilon^2R_o^4}{4}I_2\right],
\]

where \(I_2=(e^{2\epsilon H}-1)/(2\epsilon)\) and \(I_{-4}=(1-e^{-4\epsilon H})/(4\epsilon)\). Interpolating this map at P2 nodes for \(\epsilon=0.2\) yields:

| Mesh | Maximum sampled J error | Relative shear-energy error |
|---|---:|---:|
| 6×3 | 0.004631 | 0.00002263 |
| 8×4 | 0.001142 | 0.000001614 |
| 12×6 | 0.0005051 | 0.0000003232 |
| 16×8 | 0.0002835 | 0.00000009817 |

The interpolation errors decrease, supporting the shape-gradient and energy implementation. They do not measure locking in a solved forced problem. A separate nonisochoric affine deformation checks \(p_h=\lambda\log J\) and \(E_v=\lambda V_0(\log J)^2/2\), obtaining pressure error \(1.84\times10^{-13}\) and energy error \(5.73\times10^{-14}\). Protocol, complete values and source hash are in `evidence/mixed_manufactured_protocol.json` and `mixed_manufactured.json`.

## Fixed activation and material sensitivity

`runs/mixed_checks/run_fixed10.py` predeclares gain 10, logarithmic volume coefficients 5/50/500, meshes 8×4/12×6/16×8, ten continuation increments and the original no-tension geometry/materials. No calibration or target observations enter this test. `evidence/mixed_fixed10_protocol.json` is saved before the runs; `mixed_fixed10.json` records each completed case, solver/boundary status, sampled J and conservative bounds. Completion of that matrix is a numerical sensitivity study, not a scientific PASS or a continuum error bound.

All nine final fixed-gain fields were computed and their final force residuals replayed below \(10^{-7}\). Only the first two cases passed every intermediate continuation stationarity gate; the other seven retain their failed earlier stages and overall `converged=false`. A stationary final field does not retrospectively certify its entire relaxation path. `evidence/mixed_kernel_audit.json` retains every failed gain and residual.

The local-volume outcome is negative. At gain 10 on 8×4 cells, increasing \(\lambda\) from 5 to 50 to 500 gives sampled J ranges approximately 0.480–1.903, 0.518–2.039 and 0.521–2.048. On 16×8 cells at \(\lambda=500\), sampled J ranges approximately 0.424–1.975 while total volume change is only +0.0002103. These are actual quadrature evaluations, distinct from conservative Bernstein coefficient bounds. Near-zero resolved pressure residual and small total volume change therefore did not achieve local near-incompressibility in this calculation.

For example the 8×4, \(\lambda=500\) field has RMS logJ approximately 0.03044, but RMS projected logJ only 0.0002838. Its unresolved component has RMS approximately 0.03044 and would contribute about 6.55 energy units if penalized at the same \(\lambda\). The projection-kernel audit computes these quantities directly from the retained field. This is concrete evidence of the limitation identified by the projection formula, not a new theorem about biology. The separate [P2 full pointwise penalty model and forced accuracy test](QUADRATIC-ACCURACY.md) addresses that constitutive omission while retaining all failures.

The next substantive accuracy gate remains a matched forced manufactured equilibrium with exactly the same continuum solution across volume coefficients, checked under independent refinement and quadrature. Only that comparison can separate discretization error, locking, physical finite compressibility and a changed volumetric formulation. A higher-order space alone establishes none of them.

![Negative projected-pressure numerical result](../runs/mixed_checks/mixed_summary.png)

This visual uses actual sampled J values and separately displays global-volume error and unresolved log-volume modes. It does not plot conservative Bernstein bounds as attained local compression. Figure/input hashes are retained in `evidence/mixed_visual.json`.
