# Four exact-spatial finite-element minima

For each of two frozen material candidates, an independent calculation certifies one stationary root and a strict local minimum of the exactly integrated finite-element energy at small record L and at the matched large input. Each result applies to one closed cube in 611 free displacement coordinates at one fixed input. All four roots lie strictly inside their cubes.

## Numerical certificates

The common radius is \(r=2^{-34}\approx5.821\times10^{-11}\) in normalized displacement coordinates. Candidate + maximizes the frozen affine contrast direction and candidate − maximizes its negative. The final material points for this test were both frozen at half of the chosen affine steps. The table rounds upper bounds upward and the curvature lower bound downward.

![Four stacked horizontal bars give bounds on the center residual divided by the radius plus the contraction bound. Totals from top to bottom are at most 0.732, 0.381, 0.677 and 0.341. All totals fall below the dashed threshold at one.](selfmap_budget.png)

| Case | \(z\) | Self-map ratio | Error ×10⁻¹¹ | Curvature ×10⁻⁴ |
|:--|--:|--:|--:|--:|
| + · small L | 0.3976 | 0.7316 | 3.2275 | 3.6360 |
| + · large | 0.3418 | 0.3801 | 0.3389 | 5.9185 |
| − · small L | 0.3766 | 0.6767 | 2.8024 | 3.9158 |
| − · large | 0.3194 | 0.3410 | 0.1855 | 6.5236 |

The self-map ratio bounds \((\eta+zr)/r\). Error bounds the free-coordinate infinity distance from the frozen center to the exact-spatial root. Curvature bounds the smallest eigenvalue of the **normalized** Hessian on the entire cube in frozen free coordinates. It has no direct interpretation as a measured modulus. Exact values and each observed native proof are linked below.

## Why the certificates hold

Let \(E_I\) be the normalized scalar energy evaluated by exact spatial integration on the fixed quadratic axisymmetric mesh. Its gradient and Hessian are \(F_I\) and \(H_I\). The matching cached operator has derivatives \(F_Q,H_Q\) and serves as an auxiliary calculation. The same case-specific geometry, linear free map, singleton material parameters, exact gain times dose and energy normalization bind both calculations. Independently reviewed body and trace domains provide admissibility throughout each cube. Strict positive domain margins and smooth integrands on the compact reference domains give a \(C^2\) energy near the cube.

The proofs separately enclose the signed point defect \(F_I(u_0)-F_Q(u_0)\) and the symmetric whole-cube derivative defect \(H_I-H_Q\). The independent reduction checks the 144 bulk and 24 surface component maps and uses separate exact endpoint scatter and directed interval products. It adds the already-preconditioned cached center force once, **before** taking norms. With the frozen binary matrix \(R\), the resulting bounds are

\[
\begin{gathered}
\eta\ge\|RF_I(u_0)\|_\infty,\\
z\ge\sup_{u\text{ in cube}}\|I-RH_I(u)\|_\infty.
\end{gathered}
\]

The derivative bound combines the cached defect bound with the preconditioned Hessian defect bound. Since \(z<1\), \(RH_I\) and the square matrix \(R\) are invertible. The strict inequality \(\eta+zr<r\) makes \(u\mapsto u-RF_I(u)\) a contraction into the cube interior. Banach's fixed-point theorem gives a unique stationary point on the cube and

\[
\|u_*-u_0\|_\infty\le e=\frac{\eta}{1-z}<r .
\]

For the minimum, the separately certified symmetric part \(P=(R+R^{\mathsf T})/2\) is positive definite. All eigenvalues of \(A=RH_I\) have positive real part. With \(S=H_I^{-1}\), symmetry gives \(AS+SA^{\mathsf T}=2P\). Integrating the derivative of \(e^{-tA}Se^{-tA^{\mathsf T}}\) from zero to infinity gives a positive definite integral representation of \(S\). Thus \(H_I\) is positive definite. The Neumann bound on the inverse then gives

\[
\lambda_{\min}(H_I)\ge m=\frac{1-z}{\|R\|_\infty}>0 .
\]

Taylor integration along each segment from the interior stationary point gives an energy gap of at least \(m\|u-u_*\|_2^2/2\). The point is a strict local minimum and the unique minimum on its declared cube.

## Coverage and scientific meaning

These four proofs cover **4 of 36** final outcome fields. They cover the two large contrast inputs and the L field for each candidate. The two L fields span **8 of 136** primary small comparisons. The other 32 fields and their **128** small comparisons still require exact-spatial root certificates. The root errors pass the prior state-error design budgets for these four fields. A reported observation bound also requires transfer through the chosen observer, including stored-coordinate rounding and moving evaluation points.

The separate stored-field report finds a large contrast gap of about 0.154755 µm under the predeclared scenario. The frozen 5.38983567 µm allowance is a conditional model residual budget. Both saved candidates fit it under the primary disk/spatial-maximum convention. The same candidates exceed it under the other three audited static conventions. Experimental compatibility, dynamics along a stable loading path, the physical identified set, continuum and nonaxisymmetric accuracy, and literature novelty remain unresolved. The original five-gate objective remains active.

## Evidence

- [Consolidated independent four-case certificate](../bridge_reduction_v1/independent_actual_preparation/all_four_readback.json) · [Exact spatial reduction method](../bridge_reduction_v1/independent_actual_preparation/METHOD.md)
- [+ · small L proof](../bridge_reduction_v1/independent_actual_preparation/actual_replay/max_contrast_change_small_17/receipt.json) · [+ · large proof](../bridge_reduction_v1/independent_actual_preparation/actual_replay/max_contrast_change_large_03/receipt.json) · [− · small L proof](../bridge_reduction_v1/independent_actual_preparation/actual_replay/max_negative_contrast_change_small_17/receipt.json) · [− · large proof](../bridge_reduction_v1/independent_actual_preparation/actual_replay/max_negative_contrast_change_large_03/receipt.json)
- [+ · small L OS](../bridge_reduction_v1/independent_actual_preparation/actual_replay/max_contrast_change_small_17/os_exit_observation.json) · [+ · large OS](../bridge_reduction_v1/independent_actual_preparation/actual_replay/max_contrast_change_large_03/os_exit_observation.json) · [− · small L OS](../bridge_reduction_v1/independent_actual_preparation/actual_replay/max_negative_contrast_change_small_17/os_exit_observation.json) · [− · large OS](../bridge_reduction_v1/independent_actual_preparation/actual_replay/max_negative_contrast_change_large_03/os_exit_observation.json)
- [Whole-cube geometry replay](../independent_actual/replay/receipt.json) · [Positive symmetric parts of the four preconditioners](../preconditioner_spd_v1/independent_actual_preparation/actual_replay/receipt.json)
- [Fixed-field coverage ledger](../../full_outcome_coverage_v1/COVERAGE.md) · [Conditional nonlinear report](../../../report/REPORT.html) · [Observer sensitivity replay](../../../observer_robustness/independent_actual/receipt.json)
- [Earlier cached-energy minima](../cached_minimum_application_v1/MINIMA.html) · [Originality audit](../../../originality_followup/receipt.json) · [Five-gate research plan](../../../../../PLAN.md)
