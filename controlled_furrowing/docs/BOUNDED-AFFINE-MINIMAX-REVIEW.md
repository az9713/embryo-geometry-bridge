# Independent review of bounded affine minimax calibration

**PASS for the defined affine optimization problems and their source-to-coefficient dataflow.** All four actual cohort certificates were reconstructed independently, including the historical and hardened-helper extractions. This is not a certificate of mechanical derivative accuracy, nonlinear improvement, a physical optimum or the full compatible parameter set.

The [independent receipt](../evidence/bounded_affine_minimax_independent_review.json) retains executed replay programs, exact certificates, source hashes, adversarial checks and findings. It binds the current [helper](../scripts/bounded_affine_minimax.py), the frozen [cohort results](../runs/material_cohort_descent_v1/results.json), and the authoritative [version-2 exact candidates](../runs/material_cohort_descent_v1/exact_extraction_v2/exact_candidates.json). Mechanical tangent verification is a [separate independent review](MATERIAL-COHORT-DESCENT-REVIEW.md).

## What the certificate proves

For explicitly supplied residual intercepts \(r_i\), derivative coefficients \(J_{ij}\) and a bounded step box, the problem is

\[
\min_{\ell\leq s\leq u}\max_i|r_i+J_i s|.
\]

Introduce \(t\geq0\) and write all residual and coordinate constraints as \(Ax\leq b\), where \(x=(s,t)\) and \(c=(0,\ldots,0,1)\). A certificate supplies a feasible \(x_*\) and nonnegative multipliers \(\lambda\) satisfying

\[
A^T\lambda=-c,\qquad -\lambda^Tb=c^Tx_*.
\]

For every feasible \(x\), multiplication of \(Ax\leq b\) by \(\lambda\geq0\) gives

\[
-c^Tx=\lambda^TAx\leq\lambda^Tb,
\qquad c^Tx\geq-\lambda^Tb=c^Tx_*.
\]

Thus feasibility, nonnegative dual multipliers, stationarity and objective equality prove a global optimum of this exact rational LP. No floating solver status is needed for the proof once that certificate is available.

The implementation uses floating HiGHS to propose an active basis, then solves and checks exact rational systems. Bounded basis recovery can fail and return `UNRESOLVED`; it is not a complete symbolic LP algorithm. Failure to recover a certificate is not evidence of infeasibility. An unrepresentable floating proposal also returns `UNRESOLVED`. No numerical proposal is accepted solely because its floating residual is small.

## Independent helper checks

The independent test program uses determinant/Cramer's-rule vertex enumeration rather than the helper's elimination/basis routine. It covers 29 analytic or exhaustively enumerated cases, 194 feasible enumerated vertices, and 145 solver calls including residual permutations, row sign flips, coefficient scaling and coordinate permutations. Cases include fixed boxes, zero sensitivities, duplicate constraints, interior rational optima, one-sided bounds and a five-coordinate analytic example.

Four adversarial certificates each fail exactly one necessary obligation: primal feasibility, dual nonnegativity, stationarity or objective equality. Independent arithmetic and the root verifier agree on all four failures. The test program records 206 passing checks. A further 15 guard checks cover malformed sequences/certificates, dimensions, NaN/infinity, valid NumPy arrays and finite exact coefficients outside floating range.

Two implementation findings were corrected before the authoritative extraction: string containers could be interpreted as digit sequences, and an enormous finite rational could crash floating conversion. Neither produced a false mathematical optimum for the valid decoded numbers. The original source bytes and its extraction remain preserved; the hardened helper is linked to a new versioned extraction. For receipt chronology, `initial_source_snapshot` captures the already hardened `a0156ac…` revision; `reported_findings_source_sha256` and `historical_source_archive` identify the original `ac2c3bf…` bytes on which the two findings were reproduced. The independent review verifies both sets of actual certificates without rewriting either history.

## Actual cohort association and rounding checks

The cohort contains 17 biological records and 68 observed summaries. Four observer conventions reuse the same 17 baseline states at the primary disk/spatial-max gain; there are 272 candidate scalar residuals, not 272 independent observations. This common baseline differs from evaluating each convention at its earlier separately optimized gain.

The review independently checked every source-column identity/dose, baseline-field hash, actual gain-dose association, observer-record hash and mechanical archive/verification hash. All 17 baseline graph observers were re-executed. Every coordinate order and finite-difference step was checked, and all 2,720 central-difference interval calculations were reconstructed from their retained plus/minus endpoint intervals.

The actual affine tables were rebuilt without calling the root minimax or certificate checker: all 272 intercepts, all 1,360 selected derivative coefficients, all four 68-row label orders and every trust bound agree exactly. Intercepts use exact baseline-interval midpoints minus the canonical source values. Derivatives use the smaller-step interval midpoint converted to binary64, then retain that binary value as an exact fraction. The trust bounds similarly retain their binary representations. These are definitions of the affine surrogate; they are not rounding-error enclosures for an exact mechanical derivative.

Both extractions use these saved fractions verbatim. Reinterpreting them through decimal floating strings would define a slightly different LP. The hardened helper's generic floating-input convention does not change already rational inputs.

The independent arithmetic reconstructed the producer's weak-dual lower bounds and feasible upper bounds, and then checked exact primal/dual equality for both historical and authoritative extractions. The exponential/additive parameter proposal mapping was also reproduced.

| Convention | Common-baseline worst residual (µm) | Exact affine optimum, rounded (µm) | Predicted affine reduction (µm) |
|---|---:|---:|---:|
| Disk / spatial-max components | 5.379836 | 5.183863 | 0.195973 |
| Disk / regional-mean components | 6.467716 | 6.076690 | 0.391026 |
| Line / spatial-max components | 6.183270 | 5.553978 | 0.629291 |
| Line / regional-mean components | 6.183270 | 6.026441 | 0.156829 |

For the predeclared primary disk/spatial-max convention, the proposed display parameters are \(K_a=0.2210341836\), \(K_b=1.8096748361\), \(T_a=0\), \(T_b=2.2103418362\), and gain approximately 7.8126545043 in the declared normalized solver units. The compression coefficient remains fixed at 5000 by assumption. The gain display is a convenience rounding; a subsequent direct protocol must bind its exact gain-transformation and target-rounding convention explicitly.

## Scientific scope of this pass

The five coordinates are logarithmic changes in apical/basal area coefficients, an additive apical tension, a logarithmic basal tension and a logarithmic common gain. Their trust box is a numerical local-search choice. A logarithmic step of 0.1 multiplies a coefficient by \(e^{0.1}\); it is not exactly a 10% physical uncertainty. Fixing compression and geometry does not mean those quantities were identified by the data.

The saved perturbed graph coefficients are diagnostic fields, not new equilibria. Finite-difference agreement at two steps does not bound truncation error or establish differentiability of a spatial maximum. Exact smooth-mean differentiation does not make the floating mechanical tangents or nonsmooth maxima exact. The separate tangent review supplies additional numerical checks, with its own stated limits.

No equilibrium under the proposed new material parameters was solved by this cohort calculation. Consequently all four reductions in the table remain **affine predictions** until a separately frozen nonlinear replay tests them. An exact LP optimum can still propose a step that worsens the nonlinear response. A successful later replay would establish one local improvement under a conditional static observer; it would not prove a nonlinear global optimum, biological validity, source temporal-operator equivalence or completion of the original identified-set requirement.
