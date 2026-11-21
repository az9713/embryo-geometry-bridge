# Independent review of the shared-cell cohort root bounds

**Completed independent replay; mathematical outcome PARTIAL_UNRESOLVED.** Fifteen of the 17 selected stationary-root families are certified. F and K remain unresolved by the tested sufficient bounds. The [frozen receipt](../evidence/cohort_uniform_roots_independent_review.json) lists every record, exact bounds, checks and hashes. The earlier [actual-M review](PARAMETRIC-ACTUAL-M-REVIEW.md) remains immutable. The original five-gate research objective remains incomplete.

## Target and separation from the producer

The target is the exact-real stationary equation of the declared finite-quadrature, frozen-binary-coefficient P2 energy. Every record uses the same five-coordinate parameter cell, fixed reference geometry and compression coefficient, and its own exact source dose. The input strength is the exact product of the common gain and that dose; it is not silently replaced by the rounded strength recorded by the floating-point solver.

The [reviewer](../scripts/review_cohort_uniform_roots.py) restores each saved cache, approximate inverse and affine predictor. It does not optimize an equilibrium or call the producer's `build_point` or `evaluate_cell`. It uses the already reviewed local strain-coordinate derivative formulas, then separately assembles the global force, Hessian, all five force-parameter columns, Taylor remainder and contraction bounds. This is an independent assembly and proof replay; it is not an independent rederivation of every shared local algebra helper.

Record M reuses its already completed, hash-bound independent mechanical proof. Its dose is compared as an exact rational number, rather than by decimal versus fraction spelling. A cohort wrapper does not supply a new biological observation proof.

The remaining 16 records each require fresh assembly of all 9,504 quadrature terms and all 611 free coordinates at 192-bit Arb precision. No arbitrary equality of numerical upper majorants across different precisions is required. Two sound enclosures can differ while certifying the same theorem.

## Uniform root theorem

Write the free-displacement equation as

\[
F(u,\theta)=0,\qquad \theta\in C,
\]

and let \(R\) be the exact binary matrix used as a preconditioner. For each fixed parameter, use the affine predictor

\[
\widehat u(\theta)=u_0+S(\theta-\theta_0).
\]

The review bounds, uniformly over the closed parameter cell and the full corrector cube,

\[
\eta\ge\sup_{\theta\in C}\|RF(\widehat u(\theta),\theta)\|_\infty,
\quad
Z\ge\sup_{\substack{\theta\in C\\\|v\|_\infty\le r}}
\|I-RD_uF(\widehat u(\theta)+v,\theta)\|_\infty.
\]

If \(Z<1\) and \(\eta+Zr<r\), then

\[
v\longmapsto v-RF(\widehat u(\theta)+v,\theta)
\]

maps the closed corrector cube into itself and is a contraction. Banach's theorem gives exactly one fixed point in that cube for each parameter. The separate center-defect bound \(\|I-RH_0\|_\infty<1\) implies that \(R\) is nonsingular, so its fixed point solves \(F=0\). This establishes a selected local stationary root for every parameter. It does not exclude roots outside the declared branch domain.

Every independent acceptance uses the producer's **full declared corrector radius**. Downstream geometry and observers can therefore use that common cube without depending on tiny crossprecision differences in the sharper posterior estimate \(\eta/(1-Z)\).

## Force-parameter columns and residual remainder

For fixed displacement, the force is affine in the varied coefficients. At the center let \(H_0=D_uF(u_0,\theta_0)\), \(B_0=D_\theta F(u_0,\theta_0)\), \(d=\theta-\theta_0\). The exact identity is

\[
F(u_0+Sd,\theta_0+d)
=F_0+(H_0S+B_0)d
+\int_0^1[H(u_0+tSd,\theta_0+d)-H_0]Sd\,dt.
\]

Thus predictor accuracy is tested by the explicit linear defect and interval remainder; a floating-point tangent is not assumed exact. A hull containing \(td\), \(0\le t\le1\), encloses the path in the last term. The projection \(B_{\rm strain}S\) is formed before absolute values to retain useful cancellation.

For a surface contribution \(K(J_s-1)^2/2+TJ_s\), the parameter columns are \((J_s-1)\nabla J_s\) and \(\nabla J_s\). For the active bulk term \(gI\,m(b^2+d^2)/2\), the gain column is \(Im(0,b,0,d,0)^T\). The independent lightweight checks differentiate these identities symbolically.

A local Hessian difference is projected with a rigorous absolute row bound. With a strain map \(B\), a local change \(D\), and positive quadrature weight \(w\),

\[
\sum_j|(B^TDB)_{ij}|
\le\sum_a|B_{ai}|\sum_b|D_{ab}|\sum_j|B_{bj}|.
\]

Summing positive weighted contributions and multiplying their row majorants by \(|R|\) yields a valid global infinity-norm bound. Every logarithm, square root and division must remain in the analytic positive-J, positive-radius and nondegenerate-tangent region throughout the admitted box.

## Evidence and failure policy

Before expensive replay, the independent source, provider manifest, numerical versions, exact per-record result/cache/predictor/state hashes, producer protocols and any separately reviewed radius-only supplement are frozen. A parent prelaunch receipt must bind the reviewer source and its final protocol/preflight. No output is called a completed all-17 review until the batch actually terminates.

The original failed producer run and all unresolved radius tests remain retained. Failure of a sufficient contraction/self-map test means **UNRESOLVED**, not nonexistence. A radius-only retry may change the proof neighborhood while retaining exactly the same parameter cell, target equation and saved center. It must be recorded as a separate version with its own prelaunch review.

Each final per-record audit capsule binds the exact reviewed root result, shared cell and observation-contract hash. Uniform root acceptance is distinct from whole-body injectivity, uniform output bounds and all-68-row calibration compatibility. Those require their own composition proofs.

## Executed result and retained failures

The sole replay process exited successfully after 262.7 seconds. It freshly assembled 16 records at 192-bit precision: 152,064 center quadrature-term assemblies, plus the predictor-remainder and full-tube evaluations, and 16 complete 611-coordinate force/Hessian/parameter-column systems. Record M reused its immutable prior independent review through seven explicit mechanical-cache and context checks. All frozen sources and inputs still matched at the end.

There are 263 per-record checks: 257 true and six false. The false entries are exactly the three expected entries at each of F and K: producer status is not certified; the strict self-map inequality fails; and overall Banach acceptance therefore fails. In both cases the contraction bound itself is below one, and the analytic positive-J domain check passes. These are unresolved sufficient tests, not proofs that physical or discrete equilibria do not exist.

| Source column | Independent status | η upper, approximate | Z upper, approximate | Full corrector radius, approximate | (η+Zr)/r, approximate |
|---|---|---:|---:|---:|---:|
| B | PASS | 5.28746e-13 | 0.0486824 | 1.00000e-11 | 0.1015570 |
| C | PASS | 4.74721e-13 | 0.0649802 | 1.00000e-11 | 0.1124524 |
| D | PASS | 1.19872e-12 | 0.0910061 | 1.00000e-11 | 0.2108783 |
| E | PASS | 2.28339e-12 | 0.1092921 | 1.00000e-11 | 0.3376309 |
| F | UNRESOLVED | 8.65398e-11 | 0.6903874 | 1.73080e-10 | 1.1903874 |
| G | PASS | 5.53799e-12 | 0.1185511 | 1.10760e-11 | 0.6185511 |
| H | PASS | 1.45340e-12 | 0.0994095 | 1.00000e-11 | 0.2447497 |
| I | PASS | 7.57192e-13 | 0.0759205 | 1.00000e-11 | 0.1516397 |
| J | PASS | 1.79982e-12 | 0.1100145 | 1.00000e-11 | 0.2899967 |
| K | UNRESOLVED | 5.22615e-11 | 0.5653228 | 1.04523e-10 | 1.0653228 |
| L | PASS | 3.00586e-11 | 0.3831477 | 6.01172e-11 | 0.8831477 |
| M | PASS, immutable review reused | — | — | 1e−11 | — |
| N | PASS | 1.31838e-12 | 0.0925103 | 1.00000e-11 | 0.2243483 |
| O | PASS | 3.78161e-13 | 0.0606550 | 1.00000e-11 | 0.0984711 |
| P | PASS | 1.46762e-12 | 0.0983929 | 1.00000e-11 | 0.2451548 |
| Q | PASS | 1.31378e-12 | 0.0951385 | 1.00000e-11 | 0.2265164 |
| R | PASS | 1.63098e-12 | 0.1012818 | 1.00000e-11 | 0.2643797 |

The table is a diagnostic display, not the proof arithmetic. Acceptance used exact rational endpoints of Arb enclosures, retained separately for every record. Every certified record satisfies strict contraction and self-map inequalities at the exact declared full radius.

The independent batch reviewed the final v2 wrappers at their selected original radii. The failed v1 run, the decimal-versus-fraction dose guard correction and the subsequent failed radius-only supplement are retained by exact hash in the [provider manifest](../runs/cohort_uniform_root_review_v1/provider_manifest.json). The six larger-radius supplement attempts were not independently reassembled by this batch; their immutable history is distinguished from the fresh 16-record proof calculation. No parameter cell was narrowed to turn F or K into a pass.

Before replay, a reviewer-only serialization defect was caught: one guard returned a NumPy boolean rather than a native boolean. The entire unsolved source/protocol/preflight revision was preserved under `runs/cohort_uniform_root_review_v1/preflight_revision_0`, then the native-bool cast and a regression check were added. The final frozen reviewer passed 24 lightweight guards; the separate parent prelaunch review passed 18 checks. No numerical replay had started before that correction.

The 15 successful root families can be composed with separately reviewed geometry and observation enclosures. They do not supply roots or output bounds for F and K. In particular, they do not justify an all-68-row upper discrepancy bound, an IN classification or a common stable loading path. A later proof-center refinement, if attempted, is a separate versioned mathematical calculation and cannot silently replace these retained partial results.

## Limits of this result

Even a complete all-17 root review would establish only local stationary-root families of this fixed finite model. It would not prove a common stable or quasistatic trajectory, an experimental temporal maximum, ideal-quadrature or continuum-PDE accuracy, biological material identification, a full-domain identified set, or novelty. The original five-gate research objective remains incomplete.
