# The seventeen selected states are strict local energy minima

**All seventeen selected root families are now independently verified strict local minima of the fixed finite-dimensional energy on the same shared parameter cell.** A separate 192-bit replay passed every record after the producer calculation, with its own exact matrix construction and norm arithmetic.

This strengthens the existing stationary-root theorem without changing a material parameter, root, source dose, observation, or discrepancy bound. It concerns the contraction-to-geometry model motivated by [Countryman et al. (2025)](https://doi.org/10.1038/s41467-025-62483-6), not a new molecular mechanism for embryonic development.

![Exact matrix bounds and uniform normalized energy curvature](../runs/cohort_energy_stability_v1/cohort-energy-stability.png)

Each left-panel point is a rigorous upper bound for a matrix-certificate residual. Each right-panel bar is a rigorous lower bound for energy curvature in the fixed normalized free coordinates. They are not attained norms or measured eigenvalues. F/K are highlighted because their already verified root families use separately corrected numerical proof centers; this calculation makes no further root correction.

## What is being certified

Let \(E_i(u,\theta)\) be the exact reduced energy for source record \(i\), with its prescribed activation \(gI_i\). It is the scalar finite-quadrature potential defined by the frozen binary cache, restricted through the existing linear map to 611 allowed displacement coordinates. Write

\[
F_i(u,\theta)=\nabla_u E_i(u,\theta),
\qquad H_i(u,\theta)=\nabla_u^2 E_i(u,\theta).
\]

The previously reviewed root theorem supplies \(F_i(u_i(\theta),\theta)=0\) and full predictor/corrector tubes for every \(\theta\in C\). Its analytic-domain checks keep the regularized Jacobian penalty inactive and the surface terms nonsingular. Thus the exact scalar potential is \(C^2\), and the reduced Hessian \(H_i\) is symmetric throughout those tubes.

The new producer certificates establish

\[
H_i(u,\theta)\succeq m_i I,
\qquad m_i>3.94\times10^{-4}>0,
\]

for all seventeen records, every \(\theta\in C\), and every \(u\) in the corresponding declared full root tube. These numbers are normalized free-coordinate energy curvatures. They are not physical moduli or relaxation rates and depend on the coordinate basis and normalization.

At fixed activation and parameters, sufficiently small allowed displacements from the selected root therefore raise the model energy. This does not identify how the embryo reaches that state.

The cell is unchanged: relative half-width \(10^{-8}\) for the four positive material/gain coordinates, and \(T_a\in[0,2\times10^{-10}]\). The geometry, 12×6 quadratic axisymmetric mesh, finite volumetric coefficient, binary integration policy and source observations also remain unchanged. This tiny computational cell is not a physiological prior or the full compatible material domain.

## 1. The numerical preconditioner is a proof witness

Each existing root certificate uses a frozen binary matrix \(R_i\), numerically constructed as an approximate inverse Hessian. It is a **numerical proof witness**, not a new constitutive coefficient, physical mobility, or biological dynamical law. It need not be exactly symmetric.

The independent root audit already bounds, uniformly on the full tube,

\[
\|I-R_iH_i(u,\theta)\|_\infty\le Z_i<1.
\]

This is the full-tube bound \(Z_i\), not the smaller center-only defect. A positive-stable preconditioned matrix alone would not prove a symmetric Hessian positive. The additional premise needed here is

\[
P_i:=\frac{R_i+R_i^T}{2}\succ0.
\]

Every entry of this target is formed from the **exact binary rational entries** of \(R_i\), followed by addition and division by two in rigorous Arb arithmetic. A separately rounded floating symmetric matrix does not define \(P_i\).

Floating Cholesky and a triangular solve propose a square binary matrix \(K_i\). They are untrusted proposals: failure produces an unresolved result, not a theorem about eigenvalue signs. The actual certificate is the rigorous 160-bit bound

\[
q_i:=\|I-K_iP_iK_i^T\|_\infty
\le \overline q_i<8\times10^{-10}<1.
\]

The exact matrix \(K_iP_iK_i^T\) is symmetric. Its eigenvalues lie in \([1-\overline q_i,1+\overline q_i]\), so it is positive definite. Consequently square \(K_i\) is invertible, and congruence implies \(P_i\succ0\). No independent floating rank test or assumption that \(R_i\) is symmetric is used.

## 2. Why these facts imply the actual energy Hessian is positive

Fix a record, parameter and point in its tube; suppress indices. Set \(A=RH\). From

\[
\|I-A\|_\infty\le Z<1
\]

every eigenvalue of \(A\) has positive real part. In particular \(A\) is invertible. Since \(\operatorname{sym}(R)\succ0\), \(R\) is invertible, hence so is \(H\).

Because \(H=H^T\),

\[
A^TH+HA
=H(R^T+R)H
=2H\,\operatorname{sym}(R)\,H
=:Q\succ0.
\]

The positive-stability Lyapunov identity then gives

\[
H=\int_0^\infty e^{-A^Tt}Qe^{-At}\,dt\succ0.
\]

Equivalently, differentiating \(e^{-A^Tt}He^{-At}\) and integrating proves the identity directly. The exponential is an auxiliary algebraic proof device; it is not a postulated trajectory or evolution law for tissue.

This argument is applied pointwise throughout each already reviewed full tube and throughout the same \(C\). It therefore proves a uniform positive reduced Hessian, rather than positivity only at sampled states.

## 3. A quantitative lower bound without a new Hessian assembly

Neumann's bound yields

\[
H^{-1}=(RH)^{-1}R,\qquad
\|H^{-1}\|_\infty\le
\frac{\|R\|_\infty}{1-Z}.
\]

After positivity of \(H\) has been established,

\[
\lambda_{\min}(H)
=\frac{1}{\rho(H^{-1})}
\ge\frac{1}{\|H^{-1}\|_\infty}
\ge
m:=\frac{1-\overline Z}{\|R\|_\infty}>0.
\]

The matrix infinity norm of \(R\) is computed as exact rational absolute row sums of its binary entries. The numerator uses the independently audited outward upper endpoint \(\overline Z\). No dimension factor, sampled eigenvalue, or fresh Hessian assembly is hidden in this estimate.

Across the actual producer receipts, the smallest lower bound is approximately \(0.0003948742051692295\); the conservative display \(m_i>3.94\times10^{-4}\) is checked by exact rational comparisons. It is a lower bound, not the actual smallest Hessian eigenvalue.

For a selected root \(u_*\), Taylor's theorem gives the local energy-gap inequality

\[
E_i(u,\theta)-E_i(u_*,\theta)
\ge\frac{m_i}{2}\|u-u_*\|_2^2
\]

when the segment from \(u_*\) to \(u\) stays in the certified tube. The predictor-centered coordinate cube is convex, so this holds for points in that cube at the fixed parameter. In particular the selected root is a strict local minimum.

## Verification and receipts

The [theory and source-applicability review](../evidence/cohort_energy_stability_theory_review.json) verifies the nonsymmetric-preconditioner lemma, the quantitative bound, the scalar \(C^2\) target and applicability to the actual seventeen independently reviewed root tubes. It deliberately held the energy-minimum conclusion pending the separate positivity calculations.

The [new frozen protocol](../runs/cohort_energy_stability_v1/protocol.json) binds all seventeen selected caches, root certificates, specific independent root-audit capsules, source modules and the common cell. F/K use their new corrected-center caches. M uses the immutable nested independent bounds and root-audit binding. Full corrector radii are matched; center-only defects and unrelated weak-load certificates are not substituted.

The [guard tests](../runs/cohort_energy_stability_v1/guard_tests.json) include nonsymmetric matrices with positive symmetric parts, indefinite symmetric parts despite positive eigenvalues, nearly singular exact-binary matrices, rank-deficient proposals, the exact \(q=1\) boundary, and tampered/contradictory cache hashes. Authoritative expected hashes cannot be overwritten by a current-file snapshot. The initial unsolved implementation and the caught binding gap are retained under the versioned prelaunch archive.

The [producer result](../runs/cohort_energy_stability_v1/results.json), SHA-256 5f590ed1bea078e2e88d17bb39f8030d144a95fb73bf96ab10b5823951883164, contains seventeen successful matrix certificates. Per-record files retain \(K_i\), rowwise residual bounds, exact \(R_i\) norms, independent \(\overline Z_i\), and exact curvature lower bounds. The [terminal receipt](../runs/cohort_energy_stability_v1/terminal.json) records actual completion of the single authorized job. No equilibrium solve, parameter fitting or fresh Hessian assembly occurred.

The final [independent replay](../evidence/cohort_energy_stability_independent_review.json), SHA-256 d4f4c4d30d16bbbf2141cabf1cee18b47fcc35ab958a2a8a203000f0a7acaf95, passed all seventeen records at 192-bit precision: eight checks per record and five summary checks. The reviewer formed the exact symmetric part independently, used a different matrix-product association, and computed the binary row norm through integer dyadic summation. It did not call the producer certifier or assemble a new Hessian. The frozen floating \(K_i\) remained only a witness proposal. The exact \(R_i\) norms and curvature lower bounds were independently reproduced.

The review also has [twenty separate adversarial guard checks](../runs/cohort_energy_stability_review_v1/guard_tests.json). Existing full-tube root/Hessian bounds and the independently reviewed scalar-energy applicability are explicitly reused; the new contribution of this replay is the exact-preconditioner positivity and energy-minimum composition.

## What changes scientifically—and what remains open

The [68-observation result](SHARED-CELL-COHORT.md) is unchanged:

\[
5.18638714\le D(\theta)\le5.18638737\,\mu\mathrm m
\qquad(\theta\in C).
\]

This computation changes the mathematical qualification of the already selected roots, not their observational fit. The independently accepted result shows that the earlier stationary-tuple witness on \(C\) also consists of strict local energy minima of this fixed model. A declared allowance of \(5.2\,\mu\mathrm m\) still is not a measurement-noise estimate or validation of the source observation contract.

The claim is restricted to allowed displacements in the 611-coordinate constrained axisymmetric discretization. It neither tests every nonaxisymmetric or continuum mode nor proves an exact integral/PDE error bound. It does not certify a global energy minimum, a shared stable loading history, actual biological dynamics, chemical realization, source temporal peaks/reference/ROI choices, the full material/geometry identified set, or successful prediction at larger activated diameters.

Local energy positivity closes a specific mathematical gap in the selected model families. It is not proof of a new biological mechanism, a novel general concept or completion of the [original five-gate study](../PLAN.md). The methods—congruence certificates, norm bounds, Lyapunov identities and the second-derivative test—are standard; the retained model-specific evidence is what makes the conclusion assessable.

## References

1. Countryman et al. (2025), [Endogenous OptoRhoGEFs reveal biophysical principles of epithelial tissue furrowing](https://doi.org/10.1038/s41467-025-62483-6). Experimental context for the prescribed contraction-to-geometry problem.
2. Hank Yang, [Optimal Control and Estimation, Appendix C, §C.1.3, Theorem C.3](https://hankyang.seas.harvard.edu/OptimalControlEstimation/app-lti-system-theory.html). The appendix credits Shucheng Kang. The standard Hurwitz-sign Lyapunov theorem and integral apply here after replacing its system matrix by \(-A\); the self-contained argument above uses positive-stable \(A=RH\). This reference supports the mathematical identity, not a biological dynamics claim.
