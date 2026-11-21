# Direct nonlinear verification of one material-calibration step

**The primary multivariate proposal was tested with 17 newly computed nonlinear fields.** In the fixed 12×6 numerical model, its disk/spatial-max discrepancy decreases from **5.37983567 to 5.18638725 µm**. The affine prediction was 5.18386299 µm. This is a verified numerical step under a declared continuation path, not a satisfactory biological fit, nonlinear global optimum, full parameter identified set or continuum result.

The same material update was evaluated under all four retained observation conventions. The disk/regional-mean discrepancy **increases**, from 6.46771627 to 6.52848857 µm. Thus this is not an improvement under every source-plausible operator. The observation ambiguity remains scientifically consequential.

The [direct runner](../scripts/material_cohort_direct.py), [frozen protocol](../runs/material_cohort_direct_v1/protocol.json), [complete results](../runs/material_cohort_direct_v1/results.json) and [producer verification/comparison](../runs/material_cohort_direct_v1/verification_and_comparison.json) retain every field, input, attempt and residual.

## 1. What was fixed before the new computation

The preceding [cohort tangent calculation](MATERIAL-COHORT-DESCENT.md) derived four locally optimal affine proposals from the same 17 baseline fields and all 68 endpoints. Only the **primary disk/spatial-max proposal** was selected for this direct test. That continuation of the primary baseline was declared before any new nonlinear outputs; neither held-out larger-region outcomes nor the new results selected the proposal.

The normalized material/input values change as follows:

| Quantity | Old baseline | Proposed value used for all 17 new fields |
|---|---:|---:|
| Apical area modulus \(K_a\) | 0.2 | 0.2210341836 |
| Basal area modulus \(K_b\) | 2 | 1.8096748361 |
| Apical tension \(T_a\) | 0 | 0 |
| Basal tension \(T_b\) | 2 | 2.2103418362 |
| Common gain \(g\) | 8.0629720200 | 7.8126545043 |
| Compression coefficient \(\kappa\) | 5000 | 5000 |

The shear scale, fixed physical length 33.5 µm, geometry, activation support, boundary conditions, constitutive law, 12×6 finite-element space and Duffy8/surface12 integration stay fixed. Compression remains the declared numerical penalty, not a fitted physiological compressibility. Absolute physical force units are not inferred; [the force-scale quotient](PHYSICAL-PARAMETER-DOMAIN.md) still applies.

The log steps are applied using recorded binary64 exponential/material multiplication values. Gain is defined by the exact old gain Fraction times the binary64 exponential of its step. Consequently its display differs by one ulp from a convenience expression that first rounds the old gain. The exact-product convention and every requested-versus-binary64 target rounding are retained; neither is a new fit.

## 2. A new-material reference and an independently checked continuation

The new material starts from its own verified zero-activation state. An old-material state is never relabelled as a verified new-material state. The 17 measured-dose inputs are then processed in increasing order, using

\[
 f_i=g_{\rm new}I_i.
\]

The last accepted state supplies the next warm start. A failed trial may be retained and bisected but cannot become the next verified starting vector. This particular run needed **17 accepted nonzero increments and no rejected attempts**.

Every target must pass the unchanged Q8 force test, a Q12/surface16 cross-force test, exact stored-map positive-orientation and global-embedding tests, and the current-radius endpoint observer. All four conventions are evaluated at the same new fields. The other three convention-specific affine proposals have not been directly solved or verified.

Gain changed jointly with the material ratios, and the old gain was only surrogate-selected, not proved nonlinear-optimal. There is no optimized gain-only control. The decrease therefore does **not** prove that varying material ratios adds benefit beyond gain reoptimization. It establishes a better tested joint candidate for the primary convention, not a discovered material effect.

The earlier baseline fields were selected from a union of 68 ascending inputs; the new path visits 17 targets. Without uniqueness or stability proofs, loading-history differences can affect which stationary branch is reached. The test checks proposal predictions under the declared new path. It does not isolate a causal material-change contribution; no additional matched baseline-control job is silently claimed.

## 3. Actual results and what the affine predictor got right

| Observation convention | Shared old baseline, µm | Affine prediction at the primary step, µm | New direct discrepancy, µm | Old minus new, µm |
|---|---:|---:|---:|---:|
| Disk / spatial-max components | 5.37983567 | 5.18386299 | 5.18638725 | +0.19344842 |
| Disk / regional-mean components | 6.46771627 | 6.52653959 | 6.52848857 | −0.06077230 |
| Line / spatial-max components | 6.18326952 | 5.97629596 | 5.97824258 | +0.20502694 |
| Line / regional-mean components | 6.18326952 | 6.06214476 | 6.06426585 | +0.11900367 |

Every affine column in this table uses the **same primary step**. It is not the table of four separately optimized affine proposals. Every direct column uses the **same new material and gain**.

The largest sampled endpoint difference between the primary-step affine prediction and the new direct fields is approximately 0.00526180 µm. The primary uniform discrepancy differs from its affine prediction by approximately 0.00252427 µm. These are successful local numerical prediction checks at the tested inputs. They are not uniform bounds over the trust box, derivative remainder theorems, experimental uncertainty estimates or continuum error bars.

![Shared baseline, primary affine prediction and actual new fields under all four operators](../runs/material_cohort_direct_v1/material-primary-step-direct-comparison.png)

The figure has no invented physical error bars. The exact stored-field observer intervals and all 272 affine/direct comparison intervals are retained separately. The [visual receipt](../runs/material_cohort_direct_v1/visual_verification.json) binds the generator, complete result, comparison JSON and exported PNG; the PNG was opened and visually inspected after the run finished.

## 4. Exact meaning of the computed discrepancy

For convention \(v\), let \(y_{ik}\) be the 68 recorded endpoints, \(u_i^0\) the saved baseline coefficients and \(u_i^1\) the new saved coefficients. The statistic is

\[
 D_v(u)=\max_{i,k}|\mathcal O^v_k(u_i)-y_{ik}|.
\]

If the observer encloses one prediction in \([a_{ik},b_{ik}]\), then

\[
 \max(a_{ik}-y_{ik},y_{ik}-b_{ik},0)
 \le|\mathcal O^v_k(u_i)-y_{ik}|
 \le\max(|a_{ik}-y_{ik}|,|b_{ik}-y_{ik}|).
\]

Taking maxima gives an enclosure of \(D_v\) for the supplied fields. Subtracting the new upper bound from the old lower bound and the new lower bound from the old upper bound encloses their discrepancy reduction. Exact rational values are in the comparison JSON; displayed endpoints can coincide without making the physical model exact.

These bounds concern functionals of **serialized finite-element coefficient vectors**. Passing a force tolerance is not a validated enclosure of an exact nonlinear algebraic root, and it does not certify the continuum solution. No bound here controls the change between discretizations. Earlier inter-mesh shape differences also cannot be repurposed as common error bars on this new comparison.

The source contract remains conditional: a static endpoint is used as the experimental temporal-peak proxy, with flat component references, a fixed registered current-radius ROI and canonical source-column pairing. Line/disk spatial weighting, component aggregation, growth correction and temporal ordering are not resolved by obtaining a smaller residual. The [source observation contract](JOINT-OBSERVATION-CONTRACT.md) states those gaps precisely.

## 5. Verification receipts and chronology

Producer verification passes all **15 aggregate checks**, together with per-record hash, parameter, dose, force, geometry and continuation checks:

- 17 new nonzero fields and one new-material zero field;
- 17 accepted nonzero increments, zero rejected attempts;
- 2,592 exact triangle certificates including the zero reference;
- maximum Q8 scaled force residual approximately \(1.16194\times10^{-10}\);
- maximum Q12 scaled force residual approximately \(1.16075\times10^{-10}\), below the declared \(10^{-7}\) threshold;
- nonzero-field exact Jacobian enclosures jointly contained in approximately \([0.958291001,1.023404301]\).

The Jacobian interval is a conservative enclosure, not a claim that its endpoints are attained. It proves positive orientation of the stored maps, not exact or one-percent incompressibility. The global embedding certificate rules out self-overlap for those stored maps under its verified hypotheses. Neither geometry nor force checks prove stability, uniqueness or global energy minimization.

Before launch, [the concrete prelaunch review](../evidence/material_cohort_direct_prelaunch_review.json) passed ten aggregate checks and rejected 13 mutated binding/false/truthy guard cases. Recursive source checking verifies all 17 runtime local modules inside 21 declared source hashes. Two unsolved source/protocol revisions are preserved after transitive dependency omissions were caught; no numerical job preceded the corrected review.

The review file gained archive-hash metadata during authorized launch. Its exact initial bytes were recovered and retained, with the original launch hash; only that documentary field differs, while PASS checks and source/protocol/candidate bindings are identical. Both versions and their chronology are recorded in [the provenance addendum](../evidence/material_cohort_direct_prelaunch_provenance_addendum.json) and [the launch provenance](../runs/material_cohort_direct_v1/launch_receipt_provenance.json). No numerical source, parameter, policy or job was changed or restarted.

Frozen driver SHA256:

    448675cc6f59973c07a372f119aea5dd9be0f915398f57a3f9bd7d5cfe91074e

Frozen direct protocol SHA256:

    1ed458c8f0f7b9c9bacfdf0c04b907d474d06e5607ced5956ecc2d2f1adac953

Complete direct result SHA256:

    638a95bc97055e33a10e4f6c34f4750629b23f53688d19b3c8f9e7ca1ecce2b7

The [independent all-field review](MATERIAL-COHORT-DIRECT-REVIEW.md) passes **435 checks**. It replays all 17 actual new fields and the zero reference, source/parameter bindings, Q8/Q12 forces, exact geometry, nodal boundary traces, all operator outputs, retained attempts and last-verified warm-start chains. The [independent receipt](../evidence/material_cohort_direct_independent_review.json) has SHA256 `162ddbf46cf9347683dbeaf56bc29bb51b6d5c22bd255ae375516a7052ceb61f`. This review remains a conditional discrete numerical check, not experimental or continuum validation.

## 6. What this advances and what remains open

This step moves beyond a single fitted gain: a five-coordinate local material/input proposal has now been followed by new nonlinear mechanical computations for every calibration record. Its primary discrepancy decreases by approximately 0.19345 µm in this fixed numerical experiment, and the local affine prediction is close to the sampled direct response.

That does not establish robust biological improvement. One retained observation convention worsens, the primary discrepancy remains above five micrometres, the continuation paths differ, source-operator assumptions remain unresolved, and numerical/continuum uncertainty is not enclosed. The test establishes no new scientific concept or confirmed novelty.

The next requirements are validated nonlinear/parameter-cell forward bounds, externally justified discrepancy allowances or tolerance-indexed compatibility sets, material and geometry uncertainty propagation, and source-matched frozen larger-region predictions. A local optimization step and a successful replay cannot replace the full physical identified set. No incomplete original research gate is promoted to complete here.
