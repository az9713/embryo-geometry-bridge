# Numerical accuracy before adding material complexity

Status: independent numerical-methods assessment of the frozen v1 model. Three additional target-free computations have been executed; the proposed mixed solver and the full locking benchmark have **not** been implemented here. This work does not close the original geometric validation gates.

## What the current model actually enforces

The stored energy in `scripts/axisymmetric_mechanics.py`, `energy_components`, is

\[
 W(F)=\frac\mu2(I_1-3-2\log J)+\frac\kappa2(\log J)^2,
 \qquad J=\det F>0.
\]

The source parameter named `bulk` is \(\kappa\), not the infinitesimal bulk modulus. Expanding about \(F=I\) gives

\[
 W=\mu\,\varepsilon:\varepsilon
   +\frac\kappa2(\operatorname{tr}\varepsilon)^2+O(|\nabla u|^3),
 \quad K=\kappa+2\mu/3,
 \quad \nu=\frac\kappa{2(\kappa+\mu)}.
\]

Thus \(\kappa/\mu=50\) implies small-strain \(\nu=0.490196\), not an exact constraint \(J=1\). The active term is \(f|Fn_0|^2/2\); relevant nondimensional load ratios include both \(f/\mu\) and \(f/\kappa\). Large activation can cause substantial compression even when the infinitesimal Poisson ratio is close to one half. At \(J=0.8\), the penalty pressure conjugate to \(\log J\) is \(\kappa\log J=-11.16\mu\) for \(\kappa=50\mu\), comparable to prescribed activation strengths around ten. Positive \(J\) checks orientation; it does not check biological near-incompressibility.

There are two different unresolved effects: finite compressibility in the stated continuum constitutive law, and discretization error that can include volumetric locking. Their contributions cannot be separated by terminology, force residuals, or global volume change.

## Executed diagnostic and retained receipts

The driver `runs/numerical_accuracy/run_diagnostic.py` imports the frozen solver only. Its protocol was written before these computations. It reads no biological outcome files and performs no fitting. The small-region geometry is fixed at measured diameter 18.71024074 µm, height 33.5 µm, and activation 5–25 µm below the apical surface; \(H=1\), \(R_o/H=3\), \(\mu=1\), \(K_b=2\), \(K_a=0.2\). The same \(12\times6\) mesh and fixed dimensionless activation \(f=10\) are used for all cases, with ten continuation increments.

| \(\kappa/\mu\) | min \(J\) | max \(J\) | global volume change | max column shortening / \(H\) | scaled force residual |
|---|---:|---:|---:|---:|---:|
| 5 | 0.630287 | 1.289477 | −0.00310827 | 0.356179 | 9.42e−13 |
| 50 | 0.830797 | 1.118653 | −0.000644038 | 0.297379 | 6.30e−11 |
| 500 | 0.942053 | 1.050965 | −0.0000531218 | 0.147463 | 2.11e−13 |

All accepted local orientations and every continuation stage passed the solver's stationarity criterion. Nevertheless, none is certified as a continuum-accurate solution or as the global energy minimum. The bulk sweep shows strong dependence of the **discrete** geometry on compressibility. It does not prove locking: the continuum problem itself changes when \(\kappa\) changes. A global volume loss of only 0.0644% coexists with 16.9% local volume loss at \(\kappa=50\). Global volume accounting therefore cannot replace a local \(J\) audit.

Separately, the existing synthetic suite reports maximum shortening 0.0109952, 0.0135872 and about 0.01654 \(H\) on 8×4, 12×6 and 16×8 meshes at \(f=0.2\), \(\kappa=50\). Those runs use the module's default activation geometry, not the precise source geometry above. They demonstrate mesh sensitivity but cannot be combined with this sweep to form a matched 3×3 study. Neither sequence supplies a rigorous error bound.

Full fields, strain diagnostics, continuation histories, source/protocol hashes and timings are retained in `runs/numerical_accuracy/bulk_{5,50,500}.json`. `summary.json` and `evidence/numerical_accuracy_plan.json` index them.

![Fixed-load compressibility diagnostic](../runs/numerical_accuracy/compressibility-diagnostic.png)

## Recommended upgrade: mixed displacement–pressure, quadratic displacement

Use one finite-height axisymmetric implementation of the **same energy and activation**, with continuous quadratic displacement and continuous linear pressure on a shape-regular triangulation. Do not obtain a more favorable shape by silently switching constitutive law. An exact algebraic mixed potential for the present volumetric penalty is

\[
 \mathcal L(u,q)=\int_{\Omega_0}
 \left[\frac\mu2(I_1-3-2\log J)
       +q\log J-\frac{q^2}{2\kappa}
       +\frac f2\chi_A|Fn_0|^2\right]dV_0+E_s.
\]

Stationarity in \(q\) gives \(q=\kappa\log J\), recovering the frozen continuum energy pointwise before discretization. The displacement stress contribution is

\[
 P=\mu(F-F^{-T})+qF^{-T}+f\chi_A(Fn_0)\otimes n_0.
\]

The pressure residual is \(\int s(\log J-q/\kappa)dV_0=0\); displacement equilibrium is \(\int P:\nabla v\,dV_0+\delta E_s[v]=0\). This is a saddle-point Newton problem, not an unconstrained minimization in both fields. Axisymmetric integration must use \(2\pi r\,dr\,dz\), including hoop strain \(u_r/r\) and regularity \(u_r(0,z)=0\). A pressure mean constraint is appropriate only for an actual constant-pressure nullspace; it must not be imposed blindly on a traction-boundary problem.

Taylor–Hood is a well established starting point, not a blanket proof of stability for every finite-strain state or weighted axisymmetric mesh. Verify its tangent inf-sup behavior on the actual mesh family and monitor pressure oscillation. A genuine stabilized or macroelement pressure method is an alternative if its stability and constitutive consistency are established. Merely adding P1 pressure to P1 displacement is insufficient. A three-field \(u,p,J\) formulation is a defensible alternative, especially with an explicit isochoric/volumetric split, but changing that split changes the constitutive law and requires a separate retained model version.

The method is standard: the [Ferrite nearly incompressible tutorial](https://ferrite-fem.github.io/Ferrite.jl/stable/gallery/quasi_incompressible_hyperelasticity/) demonstrates mixed Taylor–Hood finite-strain mechanics; its volumetric energy differs from ours and must not be copied without adaptation. The [deal.II step-44 tutorial](https://dealii.org/9.7.0/doxygen/deal.II/step_44.html) documents a three-field approach to avoiding overly stiff displacement discretizations and traces it to Simo, Taylor and Pister (1985), DOI [10.1016/0045-7825(85)90033-7](https://doi.org/10.1016/0045-7825(85)90033-7). These sources support the method choice, not a claim that locking has already been diagnosed in this project.

## A diagnostic that can distinguish locking from physical incompressibility

First use a manufactured **small-strain, zero-surface-tension** axisymmetric benchmark with an exact nonaffine divergence-free displacement. On \(0\le r\le R\), \(0\le z\le H\), choose

\[
 \psi=r^2 f(r)g(z),\quad
 f(r)=(1-r^2/R^2)^2,\quad g(z)=z^2(H-z)^2,
 \qquad
 u_r=-\frac1r\partial_z\psi=-rf g',\quad
 u_z=\frac1r\partial_r\psi=(2f+rf')g.
\]

Then \(\operatorname{div}_{axi}u=r^{-1}\partial_r(ru_r)+\partial_z u_z=0\) identically; the smooth axis limits exist, and this field vanishes at the rim and at \(z=0,H\). Set pressure zero and manufacture body force \(b=-\operatorname{div}_{axi}(2\mu\varepsilon(u))\), with the exact displacement boundary conditions. The exact linear-elastic displacement is now **the same for every** \(\kappa\), because its divergence vanishes. Compare energy- and \(L^2\)-norm displacement errors against that known solution for both P1 and mixed discretizations. Deterioration with increasing \(\kappa\) at fixed mesh, relieved by the stable mixed method, directly identifies lack of locking robustness. An affine patch test alone is inadequate: affine fields are exactly representable by P1 and may hide this problem. This benchmark validates the tangent regime; it does not prove finite-strain convergence.

Then, without a new biological fit, execute a matched nine-case matrix on the actual small geometry: \(\kappa/\mu\in\{5,50,500\}\), \((n_r,n_z)\in\{(8,4),(12,6),(16,8)\}\), initially \(f=0.2\), then separately \(f=10\). Keep geometry, material surface terms, load, boundary and quadrature fixed. Run both formulations. Exact uniform refinement or a declared mesh generator must avoid changing the activation support. Retain failed solves. Compare full surfaces on a common physical radial domain, local \(J\), pressure, shortening, total energy and free-force residual. Quadrature must be independently raised, including the surface integrals. Extend refinement when the initial matrix is insufficient.

## Acceptance before returning to biological transfer

1. Pass zero-load, affine/isochoric and manufactured nonaffine benchmarks; verify consistent mixed tangents, pressure conventions, axis regularity and absence of a spurious pressure mode.
2. At fixed finite \(\kappa\), establish repeatable surface convergence and cross-formulation agreement. A small Newton residual verifies algebraic solution, not discretization accuracy.
3. Separate a **model choice** (finite compressibility or exact incompressibility) from mesh effects. Report local volume distortion and its localization. A biologically admissible volume tolerance requires independent evidence; do not infer it from a desired answer. If exact incompressibility is selected, verify its mixed constraint directly.
4. Recheck branch stability, outer-boundary dependence, contact/global injectivity and activation-interface refinement. P1 pointwise extrema cannot be reused unchanged for curved P2 elements; use certified polynomial/interval bounds if rigorous element-wide \(J>0\) is claimed.
5. Report empirical mesh uncertainty honestly. To call a continuum feature rigorously certified requires an applicable a posteriori bound or validated discretization and nonlinear-solution bound. Coarse/fine differences alone do not provide one, particularly for a plug-slice observer whose gap can be exactly zero.
6. Only after numerical adequacy is passed, repeat small-only calibration under the retained revised protocol, freeze all material scenarios, and evaluate larger geometry. Keep static peak-proxy and unknown slice/time qualifications.

Adding baseline tension is physically motivated but does not solve these numerical gates. Its flat-reference first variation is balanced by the clamped rim, and \(T A\) contributes \(\tfrac T2\int|\nabla w|^2\) at small normal slope, whereas a stress-free area spring begins at quartic normal slope. This explains the different linear normal stiffness; it does not establish either biological adequacy or numerical accuracy. No novel mechanics principle or confirmed novel biological result is claimed here.
