# Independent review of direct joint-gain candidate replay

**Independent computational review: PASS for all 68 completed direct targets.** The sole forward job exited successfully, its process ended, and its terminal record reports all 68 targets verified. The [independent receipt](../evidence/joint_direct_independent_review.json) retains executed replay programs, exact fractions and per-field hashes. This verifies the declared discrete candidate calculations; the original research objective remains incomplete.

The reviewed candidate procedure is [the frozen v2 protocol](../runs/joint_gain_extension_v2/protocol.json), [direct request list](../runs/joint_gain_extension_v2/direct_protocol.json), [surrogate candidates](../runs/joint_gain_extension_v2/surrogate_candidates.json), and [runner](../scripts/joint_gain_extension.py). The [surrogate review](GAIN-SURROGATE-V2-REVIEW.md) certifies the preceding finite table optimization, not the nonlinear responses evaluated here.

## Completed results and coverage

The final [direct results](../runs/joint_gain_extension_v2/direct_results.json) have SHA256 `78c779112966adaa6819778637d3abf573afbc5af61ee74ea1ec28b8b5967dea`. The independently reconstructed worst residual requirements and maximum sampled surrogate differences are:

| Observation convention | Selected surrogate gain | Worst direct residual (µm) | Maximum sampled surrogate difference (µm) |
|---|---:|---:|---:|
| Disk mean / spatial-max components | 8.062972 | 5.379836 | 0.433541 |
| Disk mean / regional-mean components | 9.509957 | 6.180135 | 0.432953 |
| Line mean / spatial-max components | 7.327159 | 5.696665 | 0.417499 |
| Line mean / regional-mean components | 8.160246 | 6.245135 | 0.434621 |

The exact lower and upper residual bounds are retained as rational numbers; they coincide at the displayed precision. Each convention's largest residual is the regional-average endpoint for source column L. All four direct worst residuals exceed their corresponding table-surrogate optima, by approximately 0.0973–0.1688 µm. The selected gains therefore remain **surrogate-selected candidates**, rather than certified nonlinear optima. The roughly 0.42–0.44 µm sampled differences establish that surrogate interpolation cannot simply be treated as exact nonlinear mechanics. They are not a uniform response-table error bound.

Coverage was exhaustive over the declared replay: 68 completed fields, 136 Q8/Q12 force evaluations, 68 exact-geometry replays, and all 272 used observer intervals. All 272 source residual rows and 272 surrogate-versus-direct error rows were independently reconstructed with rational arithmetic. The fields used 68 continuation attempts, with no rejected attempts or intermediate bisection states. The largest replayed force norms were 1.8640755e-10 (Q8) and 1.8641613e-10 (Q12), in the model's declared force normalization.

The maximum exact requested-to-binary input difference was approximately 8.539131e-16. This numerical rounding is retained explicitly; no response-sensitivity bound at the exact rational request is inferred.

The narrowest numerical certificates concern only the functionals of the stored polynomial curves. Across all 272 outputs the largest interval width was approximately 1.55×10⁻⁶ µm. This is neither a continuum discretization-error estimate nor an experimental error bar. The sharp worst-residual bounds can be much narrower because their maximizing channels need not be those with the widest individual interval.

## Exact data association and gate binding

There are 17 small-region source-column records and four observation conventions. Each convention uses one fixed common surrogate-selected gain and all four endpoints. This gives 272 predicted scalar outputs compared with the **same 68 observations**, not 272 independent biological observations.

For variant \(v\), record \(i\), the requested strength is the exact rational product

\[
A_{vi}^{\rm req}=g_v I_i.
\]

The preflight independently reconstructed all 68 distinct products and their variant/record uses. The immutable request list is ascending, begins from the verified zero state, and is bound to the candidate, response table, source data, protocol, optimizer and independent-review hashes before computation. The runner rejects duplicate launch and any changed binding. No v1 direct replay was started.

The actual numerical input is the binary floating-point value of each rational request. Its exact binary rational and difference from the requested value are retained. Subsequent stored-field intervals apply at that actual numerical strength. No response-sensitivity bound for this rounding difference is assumed.

## What is checked for each completed target

The independent prefix verifier processes only targets with completed verification records. It never takes a Newton step or changes a model field. Its retained code and per-target hashes allow the checks to be replayed.

1. Verify common physical parameters, requested versus actual strength, reached-target and discrete-stationarity flags, and immutable field/verification hashes. Reconstruct nodal coordinates from the displacement coefficients and check exact binary equality with the exported mesh and both observed boundary-control arrays, so force and geometry checks concern the same field.
2. Follow every continuation attempt from the previous target's verified vector. Independently hash the binary displacement coefficients, check each attempt's input against the last accepted vector, and advance that reference only after acceptance. Verify returned attempt artifacts and accepted-increment records.
3. Re-evaluate the frozen potential gradient at each completed target under both Duffy8/surface12 and Duffy12/surface16 integration. Separately replay the exact positive-J and embedding verifier.
4. Re-evaluate accepted internal continuation states under their declared Q8/local-orientation gate when they differ from the completed target. Those internal numerical warm starts are not silently claimed to have target-level Q12/global-embedding certificates.
5. Recompute the exact stored-curve observation intervals and bind the four used outputs to the correct variant and source-column record. The [observation contract](JOINT-OBSERVATION-CONTRACT.md) remains conditional on its ROI, aggregation, reference and static-peak conventions.

Exception attempts may have error records without returned optimizer fields. Such records must remain distinguishable from a retained failed field. A false failure or a missing result is never replaced by zero, and no failed vector is promoted into the warm-start chain.

Re-executing the frozen geometry/observer verifiers checks the new inputs and retained outputs. It is not an independent proof of those algorithms; their separate component reviews provide that evidence.

## Exact residual and surrogate-error comparisons

For a stored-field observation enclosed by \(q\in[\ell,u]\) and source value \(y\), the sharp scalar box bounds are

\[
d_{\min}=\max\{\ell-y,y-u,0\},\qquad
d_{\max}=\max\{|\ell-y|,|u-y|\}.
\]

Taking the maximum over the 68 channels/records of one convention encloses that candidate's worst residual:

\[
\max_{i,k}d^{\min}_{ik}\le\max_{i,k}|q_{ik}-y_{ik}|
\le\max_{i,k}d^{\max}_{ik}.
\]

This is a residual requirement for the supplied stored fields, not a noise estimate or a global optimum over gain/materials. The lower endpoint need not be jointly attained because the outputs arise from common mechanics.

The same interval-distance calculation compares each direct field output with its exact midpoint-surrogate prediction at the requested gain-dose product. That is an observed discrepancy at the candidate points, including the distinction between requested and actual numerical strength. It is not a uniform interpolation-error bound between all response-table nodes.

## Source observation caveat retained at final review

The primary article's Methods, section `Sec16` (Furrow quantification), describes averaging local shortening peaks over the activation region. A literal reading therefore favors the **mean of local temporal maxima**, rather than the temporal maximum of the regional mean. For a fixed tracked domain and nonnegative normalized weights,

\[
\sup_t\langle s(\cdot,t)\rangle
\leq\langle\sup_t s(\cdot,t)\rangle.
\]

The operations need not agree. The frozen direct calculation uses a single stored mechanical state and declared line/disk and component-aggregation alternatives. It therefore does not reproduce the source temporal observer merely by passing its static interval checks. In particular, a separate zero-residual construction for the peak of a regional mean would not establish a source-aligned fit. An exact implementation crosswalk and time-dependent observations remain necessary.

The same primary record states that no custom code was generated, while primary microscopy is available by author request. The unresolved request is for the exact analysis workflow/settings and derived tracked time courses, rather than for allegedly declared but missing custom code. No author contact was made in this review.

## Limits that remain after any successful replay

The observer intervals bound functionals of stored finite-element curves with binary coefficients. They do not enclose the exact algebraic equilibrium, continuum solution or experimental response merely because a force tolerance is passed. A force residual needs a stability/error argument before it becomes a displacement-error radius.

An embedded positive-J stationary field is not by itself a proof of energetic stability or uniqueness. This review supplies no new stability/inertia or branch-uniqueness theorem.

This calculation concerns one fixed material/geometry scenario and four declared static observation conventions. It does not recover the full compatible parameter set, settle the source component/reference ambiguity, establish experimental temporal peak coverage, fit a biological error law or provide a novel concept. A completed direct replay can test the surrogate's candidate predictions; it cannot by itself certify the surrogate gain as the nonlinear global optimum.

The four conventions are different measurement hypotheses, not independent datasets or a pair of mechanisms compatible under one fixed observation contract. Different selected gains across those conventions do not establish the proposed transfer-ambiguity theorem. Surrogate-to-direct differences may also include numerical branch/path effects; they are not automatically attributable solely to interpolation.
