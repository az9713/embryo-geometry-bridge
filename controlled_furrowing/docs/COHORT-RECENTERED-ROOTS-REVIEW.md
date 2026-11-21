# Independent F/K proof-center refinement review

**Independent replay PASS for both new F/K local root certificates.** The [new frozen review receipt](../evidence/cohort_recentered_roots_independent_review.json) retains the exact outcome, enclosures and hashes. The [earlier 15-of-17 result](COHORT-UNIFORM-ROOTS-REVIEW.md) remains immutable. This completes the local root-family component on one fixed computational cell; the original five-gate research objective remains incomplete.

## A numerical proof center changed; the physical equation did not

The old F and K tests could not prove their self-map inequalities. That is a limitation of those sufficient interval bounds, not a nonexistence theorem. Simply enlarging the tested correction cubes made the derivative majorants worse. The new calculation proposes one closer proof center while preserving the physical model and the same five-coordinate parameter cell.

For an original binary seed \(u_0\), exact binary preconditioner \(R\), and a stored dyadic midpoint \(m_i\) of the producer's interval \((RF(u_0,\theta_0))_i\), the numerical proposal is

\[
u_{1,i}=\operatorname{RN}_{64}\!\left(\operatorname{exact}(u_{0,i})-m_i\right).
\]

There is one nearest-even binary64 rounding of the exact subtraction. Rounding the correction to binary64 first can produce a different answer; the independent guard includes an explicit counterexample. All 611 proposed coordinates and their exact actual translations \(u_1-u_0\) are checked.

This is a change of proof center, not a material fit, a new biological observation, a changed activation dose or another nonlinear equilibrium optimization. The original seed state stays identified separately from the explicitly labelled auxiliary proof-center state. The latter is not falsely marked as a converged experimental or solver result.

## Exact target and cache audit

The old and new cache arrays must have identical names, shapes, dtypes and bytes **except for the center vector `u`**. This includes all reference and quadrature data, active masks, force normalization, mesh, free indices and the approximate inverse \(R\). The affine predictor file \(S\) is bit-identical. Model coefficients, exact common-gain–dose law, baseline parameter vector and shared cell also stay fixed.

The independent verifier recomputes the array and full-potential fingerprints. It separately checks three identities:

- the original tangent archive's `free_vector` belongs to the old seed \(u_0\);
- the new cache and auxiliary state's `free_displacement` both equal \(u_1\);
- the exact recorded center translation is \(u_1-u_0\).

The old float Hessian/tangent archive is not relabelled as a newly evaluated archive at \(u_1\). Its preserved preconditioner and predictor are proposals. Their adequacy must be demonstrated by fresh interval force, Hessian and parameter-defect calculations at the new center.

The review reproduces the exact single-rounding recipe from the recorded midpoints and checks their inclusion in the recorded force-correction balls. It does not independently repeat the old-center force assembly merely to establish that the proposal was a good Newton step. Root existence is instead proved from scratch at \(u_1\); that theorem is valid for any proposed center passing its fresh tests.

## Fresh uniform root theorem

The new predictor is

\[
\widehat u(\theta)=u_1+S(\theta-\theta_0),\qquad\theta\in C.
\]

The independent implementation uses the reviewed local strain algebra and separately assembles the global force, Hessian, all five force-parameter columns and Taylor remainder at \(u_1\). It does not call the producer's point or cell routines and performs no equilibrium optimization. Each record requires all 9,504 terms and 611 free coordinates at 192-bit Arb precision.

With \(H_1=D_uF(u_1,\theta_0)\) and \(B_1=D_\theta F(u_1,\theta_0)\), the exact residual decomposition is

\[
F(u_1+Sd,\theta_0+d)=F(u_1,\theta_0)+(H_1S+B_1)d
+\int_0^1[H(u_1+tSd,\theta_0+d)-H_1]Sd\,dt.
\]

It remains valid even though \(S\) was computed at the old seed. Its accuracy is tested through the new explicit defect and remainder, rather than assumed.

The proof then bounds, uniformly over the unchanged cell and the declared corrector cube,

\[
\eta\ge\|RF(\widehat u(\theta),\theta)\|_\infty,
\qquad
Z\ge\|I-RD_uF(\widehat u(\theta)+v,\theta)\|_\infty.
\]

Strict \(Z<1\) and \(\eta+Zr<r\) imply one unique local stationary root for every parameter through the Banach argument. Positive Jacobian and nondegenerate surface domains are checked throughout every admitted strain box. A center-defect bound below one establishes nonsingularity of \(R\), so the fixed point is a root of \(F\).

The normal corrector radius \(r\) is measured about \(\widehat u(\theta)\). It is **not** the center translation \(u_1-u_0\). The fresh proof uses \(u_1\) as its center; an observer using the old saved geometry must add the exact center translation separately to its geometric error envelope, or explicitly rebase to the auxiliary geometry with its own serialization gap.

## Executed independent result

The sole independent process exited zero after 32.63 seconds. It freshly assembled both 9,504-term, 611-coordinate systems at 192-bit precision. Both records passed their 19 root/metadata checks and 16 explicit lineage checks; all six final checks passed. These counts include logical summaries and overlapping provenance checks, not independent scientific experiments.

| Record | η upper, approximate | Z upper, approximate | Full radius r | (η+Zr) upper, approximate | Exact center-translation norm, approximate |
|---|---:|---:|---:|---:|---:|
| F | 7.4237296711e-13 | 0.0754081646 | 1e−11 | 1.4964546131e-12 | 8.5844803488e-11 |
| K | 2.4987031192e-12 | 0.1289284033 | 1e−11 | 3.7879871524e-12 | 5.0061621515e-11 |

The displayed values are approximate; acceptance used exact Arb endpoint fractions. Displacement coordinates and radii here are normalized solver coordinates, while Z is dimensionless. The translation column is the exact binary-center difference displayed approximately. It is larger than the normal corrector radius in both examples and has not been confused with that radius.

Fresh positive-J analytic-domain bounds also passed: the full strain tubes have quadrature-J lower bounds above approximately 0.988968597 for F and 0.964298642 for K. Whole-body injectivity and output transfer require their own separate composition proof; quadrature positivity alone is not substituted for it.

These two certificates combine with the immutable 15 earlier certificates to provide a selected local stationary-root family for every one of the 17 source inputs on the same parameter cell. They do not rewrite the earlier failed tests or automatically produce eight new observer bounds. Those require a separately reviewed center-offset or rebased geometry/observation composition.

The lesson is numerical: the old sufficient bounds failed around less accurate proof centers, while one arithmetic center update permitted a rigorous proof for the unchanged equations. This is standard validated-numerics machinery, not an established new biological principle or a novelty claim.

Even all 17 root families and all 68 static output enclosures would supply only a conditional stationary-root-tuple calculation for this fixed finite model. A physical loading trajectory, experimental temporal maxima/ROI alignment, uncertainty in data/geometry/material domain, other branches and parameter tails, continuum-PDE error and novelty remain separate original requirements.
