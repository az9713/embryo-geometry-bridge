# From one height fit to a joint mechanical calibration

**Direct nonlinear replay of all 68 candidate strengths is now complete and independently verified.** The same 68 small-region endpoint entries are analysed under four retained observation conventions in one fixed mechanical scenario. The expanded strength20 table gives interior, independently certified surrogate candidates. Actual candidate discrepancy requirements are 5.3798–6.2451 µm, with maximum sampled surrogate/direct differences of 0.4175–0.4346 µm. This completes candidate verification within Gate3; it does not complete its physical compatible-set requirement or any new gate.

## What is fixed, and what is being inferred?

The conditional mechanical scenario uses the full pointwise P2 energy, 12×6 mesh, Duffy8/surface12 integration, compression coefficient/shear ratio 5000, basal/apical area coefficients 2/0.2 and basal/apical tensions 2/0. The flat reference height is 33.5 µm and the reported activation diameter is 18.71024074 µm. These are one declared numerical/material scenario, not independently measured material constants or an exhaustive parameter domain.

The fitted variable is one common nonnegative gain \(g\), mapping measured activation \(I_i\) to model strength \(f_i=gI_i\). For each candidate, the intended loss is a common micrometre minimax discrepancy over all four summaries:

\[
D(g)=\max_{i=1,\ldots,17}\max_{k=1,\ldots,4}
|\mathcal O_k(q(gI_i))-y_{ik}|.
\]

This objective is a transparent deterministic choice. It does not assume that the four summaries are independent, estimate measurement noise, produce a confidence interval or calibrate biological force units. Per-channel residuals remain available so one large mismatch is not hidden by an average score.

The present static-endpoint closure is conditional: \(q(gI_i)\) is one selected stationary discrete field, not a demonstrated experimental temporal maximum, unique equilibrium or certified energy minimizer. The [source contract](JOINT-OBSERVATION-CONTRACT.md) and [observer](QUADRATIC-ENDPOINT-OBSERVER.md) retain growth-reference, ROI transport, time/average ordering and component-reference limits. This calculation does not silently resolve them.

## Four observation alternatives are retained

The source-plausible alternatives combine line versus projected-disk spatial weighting with spatial-maximum versus regional-mean components. Each combination has its own gain candidate. The canonical common source-column association is explicit; the [joint constraint module](JOINT-CALIBRATION-CONSTRAINTS.md) also implements a separate unknown-pairing sensitivity, but that is not silently substituted into the current gain fit.

These alternatives are not selected using larger-region outcomes. No held-out response file enters the preparation, fitting or direct-calibration driver.

## What the exact surrogate optimizer proves

The first table contains verified mechanical fields at strengths \(0,1,5,10\). The endpoint observer encloses the four static quantities at each node. Taking each interval midpoint defines an **exactly specified numerical surrogate**; linearly interpolating those values does not establish an interpolation-error bound for the actual mechanics.

For this defined surrogate, \(\widehat{\mathcal O}_k(gI_i)\) changes affine formula only at gain values

\[
g=\frac{f_j}{I_i}.
\]

Partitioning the supported gain interval by all such breakpoints makes every signed residual affine on each cell. Its worst absolute residual is a convex upper envelope

\[
\widehat D(g)=\max_m(a_mg+b_m).
\]

A cell minimum has an exact certificate: at an interior point zero belongs to the interval between the smallest and largest active slopes; at the left boundary the right derivative is nonnegative; at the right boundary the left derivative is nonpositive. A floating-point linear program only proposes candidates. Rational replay of every affine inequality and the active-slope condition decides whether a certificate is accepted. Comparing the certified cell minima gives the global optimum **for this surrogate on the supported interval**.

The [independent optimizer review](GAIN-SURROGATE-REVIEW.md) verifies all 108 actual cells across four variants and 14,796 affine lines against the source-linked table/data, plus independent small examples. It does not transform midpoint interpolation into a physical forward bound.

## First table: all four candidates hit the boundary

The largest observed dose is 1.51057344. Avoiding extrapolation beyond strength ten therefore restricts the initial search to

\[
0\le g\le\frac{10}{1.51057344}
=\frac{15625000}{2360271}\approx6.620002534.
\]

This upper limit comes from available computation, not physiology or data-based identification. The exact initial surrogate results are:

| Weighting | Component convention | Surrogate minimax discrepancy | Gain location |
|---|---|---:|---|
|Projected disk|Spatial maxima|5.93810437 µm|Upper support boundary|
|Radial line|Spatial maxima|5.93810437 µm|Upper support boundary|
|Projected disk|Regional means|7.27136527 µm|Upper support boundary|
|Radial line|Regional means|6.78228355 µm|Upper support boundary|

These are neither noise estimates nor nonlinear physical minima. They show that the initial support does not close the gain search. Reporting the boundary value as an inferred best biological gain would be unjustified.

The first source/protocol/table/candidate files remain unchanged under `runs/joint_gain_calibration/`. Its direct replay was deliberately not started before addressing this support issue. The [new versioned calculation](JOINT-GAIN-CALIBRATION.md), under `runs/joint_gain_extension_v2/`, first extends the unchanged 12×6 higher-rule table to strength20. It may extend to40 and80 only if a declared support boundary remains active. These are computational support budgets, not measured prior limits. A failed or unresolved extension remains a reported limitation.

Once candidate support is interior and the independent surrogate gate passes, direct nonlinear solves at the distinct \(g_{\rm variant}I_i\) values test all 68 endpoint predictions per variant. The driver retains every failed attempt, last-verified warm-start chain, observer interval, state/source hash and surrogate-to-direct discrepancy. A small direct-surrogate difference would support candidate replay at those doses; it would not prove a global nonlinear optimum or an entire identified set.

### Extended table: candidate gains are now interior

The same 12×6 mechanics has verified the added strength20 state, including its higher-rule force and exact-geometry checks. All four defined surrogate optimizers are now interior to the extended support, so no strength40 or80 run is needed for this stage:

| Weighting | Component convention | Surrogate candidate gain | Surrogate minimax discrepancy |
|---|---|---:|---:|
|Projected disk|Spatial maxima|8.062972|5.246851 µm|
|Projected disk|Regional means|9.509957|6.024626 µm|
|Radial line|Spatial maxima|7.327159|5.599342 µm|
|Radial line|Regional means|8.160246|6.076291 µm|

The exact rational candidates, affine certificates, observer intervals and state/source hashes are in the versioned run directory. These entries remain surrogate results, and their selected inputs have now been directly replayed. The [independent hash-bound gate](GAIN-SURROGATE-V2-REVIEW.md) has passed: 168 cell certificates and 23,016 source-linked affine lines were verified. The sole direct job finished all 68 distinct candidate strengths from the verified zero state, with 68 accepted attempts and no rejected attempts. The [independent direct review](JOINT-DIRECT-REVIEW.md) passes all 16 final checks, including every field’s Q8/Q12 force, exact geometry, observer and residual reconstruction. These verification audits are not newly completed research gates.

Even in direct replay, the exact spatial intervals refer to the stored finite-element coefficient vector. The force-tolerance check does not itself provide a validated enclosure of an exact nonlinear algebraic root, much less the continuum solution. The algebraic-solution, interpolation, mesh, source-observation and model-discrepancy errors remain distinct.

## Completed direct candidate results

| Convention | Surrogate-selected gain | Direct candidate discrepancy | Largest sampled surrogate/direct difference |
|---|---:|---:|---:|
| Disk / spatial-max components | 8.06297202 | 5.37983567 µm | 0.43354056 µm |
| Disk / regional-mean components | 9.50995731 | 6.18013536 µm | 0.43295265 µm |
| Line / spatial-max components | 7.32715933 | 5.69666515 µm | 0.41749944 µm |
| Line / regional-mean components | 8.16024598 | 6.24513544 µm | 0.43462117 µm |

The source contains 68 observations analysed four ways, not 272 independent observations. All 272 residual entries and sampled surrogate/direct differences are retained. All four worst residuals are regional-mean shortening; numerical verification does not turn these substantial discrepancies into a successful biological fit. The [complete study](JOINT-GAIN-CALIBRATION.md) supplies exact interval formulas, per-channel results, three inspected figures and receipt links.

There are 9,792 exact stored-map triangle certificates. Maximum Q8/Q12 scaled force residuals are about \(1.8641\times10^{-10}\), below \(10^{-7}\); the combined exact Jacobian enclosure is about \([0.955379,1.026142]\). Positive orientation is certified, but one-percent incompressibility is not. These are stationary stored fields, not continuum-error certificates or proven stable/global energy minima.

The [final direct result](../runs/joint_gain_extension_v2/direct_results.json) SHA256 is `78c779112966adaa6819778637d3abf573afbc5af61ee74ea1ec28b8b5967dea`. Exact surrogate optimality does not certify optimality for nonlinear mechanics. Sampled surrogate/direct differences are not uniform interpolation bounds and may include branch/solver-path effects.

## Why the finer numerical study still matters

The separate 24×12 sweep has reached strengths1,5,10. At strength10 its material-column maximum changes only approximately0.03394 µm from the matched 16×8 field, whereas the exact discrete Eulerian apical-surface discrepancy is enclosed in approximately[0.670904,0.680214] µm. Thus a stable-looking scalar statistic can coexist with substantial unresolved shape sensitivity. These are inter-mesh comparisons, not continuum error bounds.

The joint-calibration calculation remains a conditional discrete experiment while that error issue and the source-observation assumptions are unresolved. It cannot close the original shape-prediction gate merely by returning a gain and passing force tolerance.

## What remains under the original objective

The full material/geometry compatible set is still uncomputed; no discrepancy law is independently measured; source-matched geometry/time evaluation remains incomplete; and revised uncertainty-propagated larger-region predictions have not been frozen and scored. Numerical readiness, a one-scenario gain fit, a matching certificate and a source audit address different requirements. No confirmed novelty follows from these standard methods or from the current partial results.

## Completed local parameter-inference step

The planned **multivariate local descent using all 17 records and all four endpoints** has now executed. It starts from the lowest-discrepancy existing conditional candidate, disk weighting with spatial-max components, selected using small-region calibration alone. At its same 17 base states, mechanical displacement tangents cover five material/input coordinates: \(\log K_a,\log K_b,T_a,\log T_b,\log g\). The logarithmic compression coefficient stays at 5000 as the declared numerical penalty and geometry stays fixed.

All four observation alternatives use those same base fields, rather than each variant's separate gain-fit branch. The [local study](MATERIAL-COHORT-DESCENT.md) proposes one minimax step per convention within a declared numerical trust region; exact primal/dual certificates cover the defined affine problems. Two graph-difference steps and independent smooth-mean derivatives are checked, while nonsmooth maxima, tangent truncation and branch invertibility remain explicit limitations.

The predeclared primary joint material/gain proposal then received a [complete 17-dose nonlinear replay](MATERIAL-COHORT-DIRECT.md) and [435-check independent review](MATERIAL-COHORT-DIRECT-REVIEW.md). Its primary discrete discrepancy decreases from 5.37984 to 5.18639 µm, but the disk/regional-mean comparison worsens from 6.46772 to 6.52849 µm. Gain changes jointly with material ratios; no optimized gain-only control establishes a material-specific benefit. The old 68-input and new 17-input continuation histories differ, so branch identity is unproved. This is a better tested joint candidate for one static operator, not a globally optimal or observation-robust biological model.

This advances inference beyond one fitted gain without pretending that a local trust region is a physical prior or an exhaustive parameter set. The subsequent rigorous target remains a tolerance-indexed compatibility calculation with uniform forward enclosures. For a declared parameter cell \(C\), every endpoint needs a validated enclosure \(\mathcal O_k(q(\theta,I_i))\in\mathcal B_{ik}(C)\) throughout the cell. Cells can be excluded when some response enclosure is disjoint from its allowed observation band, accepted wholly only when every enclosure lies within its band, and otherwise remain unresolved. Point solves, tangent predictions and midpoint interpolation alone do not supply those uniform certificates.

Separately, an [independently validated root certificate](ALGEBRAIC-ROOT-INDEPENDENT-REVIEW.md) and [actual observer-error transfer](ALGEBRAIC-ROOT-OBSERVER-ERROR.md) now bound the algebraic-error contribution for one old archived weak-load state. This is not certification of all seventeen new fields or the whole admitted parameter domain, and it gives no continuum error bound.

The **full physical identified set**, **uniform validated physical/continuum forward bound**, and **source-matched frozen larger-region predictions** remain **INCOMPLETE** under the original five-gate contract. No new gate becomes complete through these discrete component results. The [current original-scope audit](CURRENT-FIVE-GATE-AUDIT.md) lists the remaining requirements.
