# Joint calibration candidate using all four small-region endpoints

This is a **conditional candidate study in one fixed mechanical scenario**, not a completed physical identified-set calculation. It uses all 17 source-column records and all four recorded endpoints—68 numbers—with four explicitly retained observation conventions. The same observations are analysed four ways; this is not four independent data sets.

The first response-table surrogate selected its computational gain cap for every convention. That boundary value is not reported as an inferred biological gain. A separately versioned extension to strength 20 has now produced four interior **surrogate candidates**. Their independent 22-check certificate audit passes. The sole 68-target direct nonlinear replay has now finished successfully. Its four candidate discrepancy requirements are approximately 5.38–6.25 µm. These are substantial residuals in this fixed scenario, not a successful reconstruction or a globally optimal physical fit.

The implementation is `scripts/joint_gain_extension.py`, with protocol, states and progress retained under `runs/joint_gain_extension_v2`. Parent preparation under `runs/joint_gain_calibration` remains unchanged. No held-out larger-region response file enters this study, and no 24×12 job is launched or duplicated here.

## 1. Fixed model and explicit observation choices

The mechanical scenario is the existing 12×6 higher-quadrature P2 model: height scale 33.5 µm, small-patch diameter 18.71024074 µm, activation 5–25 µm below the apical surface, shear scale one, logarithmic-volume coefficient 5000, basal/apical area moduli 2/0.2 and basal/apical tensions 2/0. Its actual integration is Duffy8/surface12; each new accepted state is also checked with Duffy12/surface16. Material, geometry and numerical policies are unchanged.

The observer supplies current-radius vertical height, independently maximized or regional-mean component motions, and line/disk regional means. Its source ambiguities and mathematical certificates are documented in `QUADRATIC-ENDPOINT-OBSERVER.md` and `JOINT-OBSERVATION-CONTRACT.md`. All four conventions remain:

- current-disk mean with spatial-maximum component endpoints;
- current-line mean with spatial-maximum component endpoints;
- current-disk mean with regional-mean component endpoints;
- current-line mean with regional-mean component endpoints.

Every convention additionally assumes a static endpoint represents the measured peak, a flat initial reference with component heights H/0, a registered current-radius ROI equal to the initial patch radius, and canonical source-column pairing. These are declared candidate assumptions, not newly established facts about the experiment. Temporal coverage, source counterfactual growth and component reference splitting remain unresolved.

## 2. The defined exact surrogate problem

For dose \(I_i\), four observed endpoints \(y_{ik}\), and observation convention \(v\), let \(\widetilde r_{vk}(f)\) be piecewise-linear interpolation of **interval-midpoint** response-table values. The candidate gain minimizes

\[
 L_v(g)=\max_{i=1,\ldots,17;\,k=1,\ldots,4}
 \left|\widetilde r_{vk}(gI_i)-y_{ik}\right|.
\]

All channels receive the same micrometre minimax loss. This is a declared discrepancy objective, not a likelihood or estimated error distribution.

On gain cells separated by \(f_j/I_i\), every residual is affine in gain, so the loss is the upper envelope of its positive and negative affine forms and zero. `piecewise_gain_minimax.py` uses a floating LP only to suggest candidate locations, then verifies an exact rational optimality condition for every cell: an interior zero subgradient, or the appropriate one-sided boundary derivative. All cells are covered and their certified minima compared. The resulting certificate is global for this **defined finite response-table surrogate and gain interval only**.

It does not enclose physical nonlinear interpolation errors, determine the nonlinear globally optimal gain, or cover the full material/geometry domain. An observer interval midpoint is a defined surrogate number, not an exact unknown physical endpoint.

## 3. Why the original cap must be extended

The original verified strength nodes were 0, 1, 5 and 10. Forbidding table extrapolation imposed

\[
 0\le g\le\frac{10}{\max_i I_i}
 =\frac{15625000}{2360271}\approx6.62000253.
\]

All four exact surrogate optima hit that cap. The spatial-component variants had loss approximately 5.9381 µm; the disk-regional-mean and line-regional-mean component variants had losses approximately 7.2714 and 6.7823 µm. These are truncated-surrogate results, not evidence that the cap is physiological or that the gain is identified.

Version 2 freezes a support schedule **before new computation**: add strength 20, then 40 and 80 only if at least one retained convention still selects a declared gain boundary. These limits are computational support limits, not physical priors. Every expansion refits all four conventions without extrapolation. Reaching an interior optimum only removes the visible boundary-selection problem on that expanded surrogate; it does not prove that still larger support or the true nonlinear response cannot change the optimum.

The first extension starts from the verified strength-10 higher-rule field. Adaptive continuation must reuse only the last verified state, retain every failed trial and pass force, higher-rule cross-force, exact geometry and observer-resolution gates before a node enters the table. A target that cannot be reached or observed is retained as unresolved, not silently dropped or converted into a response value.

## 4. Independent candidate gate and direct replay

When all four expanded-surrogate candidates are interior, an independent audit must bind the exact candidate, response-table, protocol and optimizer-source hashes. It reconstructs all gain cuts, all 68 residual affine forms per cell, every signed subgradient/boundary certificate and the global cell minimum. Direct replay cannot start on an unmatched or incomplete review receipt.

The direct stage then forms the union of all distinct exact strengths \(g_vI_i\), orders them increasingly, and starts from the separately verified zero state. At most 68 distinct target strengths are needed. Coincident targets are reused while retaining each convention/record association. The solver receives binary64 strengths; the requested exact value, actual binary64 value and their rational difference are recorded separately.

Each completed target retains its full field, failed and accepted continuation attempts, state hash, higher-rule force evaluation, exact discrete geometry and all four observation alternatives. Every convention must receive all 17 records and all four endpoints. No high-residual record is skipped. A failed target prevents a complete-candidate claim.

For each convention, the final report must include:

- all 68 direct residual intervals;
- the uniform discrepancy requirement at this candidate, with observer interval uncertainty;
- all 68 surrogate-to-direct error intervals at the actual requested strengths;
- the largest sampled surrogate error and complete source/field provenance.

These discrepancy requirements are not estimated biological noise. Surrogate-to-direct checks at selected strengths are not a uniform nonlinear interpolation bound. A direct replay validates candidate predictions at those inputs; it does not make the selected surrogate gain a nonlinear global optimum or compute \(\Theta_\delta\).

Ascending continuation is a declared numerical branch-selection convention, not evidence that the separate experiments followed that loading path. The seed table and direct replay can also differ in nonlinear initialization history. Reported surrogate-to-direct discrepancies may contain interpolation and branch/solver-path effects; they do not isolate either contribution without an additional comparison.

## 5. Retained status and reproducibility

The frozen version-2 source is

```text
b54b46251fc4df0c057765c55dcc64305e359645b97fbda34fd54ecaaff99e37
```

and its pre-computation protocol hash is

```text
3b72bcf4583c0b136d01c420028a20a9324dc9c0b3b4ef212b14ec76a62e2031
```

Six pre-launch guard tests pass, including exact linear surrogate interpolation, no extrapolation, output-path resolution and preservation of the original unrun direct stage. An earlier preparation revision is retained under `pre_launch_revision_0`; no mechanical extension preceded that source-hardening revision.

The strength-20 support node was reached in one accepted increment from the verified strength-10 field, with no rejected attempts. Its Q8 force residual is approximately 5.515e−12 and Q12 cross residual 5.452e−11. Exact discrete geometry passes; its correlated whole-body J interval is contained in [0.950008, 1.031191], which is a positive-orientation certificate rather than an exact incompressibility claim. The original fields remain intact.

All four candidates became interior after this extension, so strengths 40 and 80 were **not** computed:

| Observation convention | Exact-surrogate gain, displayed | Exact-surrogate minimax loss, displayed µm |
|---|---:|---:|
| Disk / spatial-max components | 8.06297202 | 5.24685117 |
| Disk / regional-mean components | 9.50995731 | 6.02462625 |
| Line / spatial-max components | 7.32715933 | 5.59934156 |
| Line / regional-mean components | 8.16024598 | 6.07629139 |

These numbers are exact optima only of the defined midpoint-table surrogate on its expanded computational interval. Eight numerical-support/provenance checks pass in `evidence/joint_gain_extension_support.json`. They do not replace the independent LP audit, direct forward checks or physical inference.

Current extension status is recorded in `runs/joint_gain_extension_v2/extension_job.json`, with node progress in `extension_progress.json`. The independent LP audit `evidence/piecewise_gain_minimax_v2_independent_review.json` passes all 22 checks, reconstructing 43 gain breakpoints, 168 cell certificates and 23,016 affine coefficients. It binds all candidate/table/protocol/optimizer hashes and preserves the original truncated v1 record.

The direct stage has a separate immutable `direct_protocol.json` and live `direct_job.json`. It contains 68 distinct targets, ranging from approximately 1.084496 to 14.365489 in normalized strength, all inside the verified table support. The original v1 direct stage was never run. The sole direct session exited successfully, and `direct_job.json` records `FINISHED` with all 68 targets verified. Final fields, residuals and hashes are retained in `direct_results.json`; the completed results follow below.

```powershell
python -X utf8 controlled_furrowing/scripts/joint_gain_extension.py prepare
python -X utf8 controlled_furrowing/scripts/joint_gain_extension.py extend
python -X utf8 controlled_furrowing/scripts/joint_gain_extension.py direct --review controlled_furrowing/evidence/piecewise_gain_minimax_v2_independent_review.json
```

Existing preparation/job records prevent accidental duplication. Inspect the current live handle and terminal record before invoking a stage. This conditional single-scenario exercise does not complete the original numerical, source-observation or full-parameter inference gates and establishes no confirmed novelty.


## 6. Completed direct nonlinear replay

All **68 distinct target strengths** were reached in the one ascending continuation job. This required 68 accepted attempts and no rejected attempts. Every field retains its full nodal deformation, source/parameter information, continuation record, Q8 evaluation, Q12 cross-evaluation and exact discrete geometry certificate. Every convention receives all 17 source records and four endpoints; the final result contains 272 residual entries, reflecting four analyses of the same 68 observations.

| Retained observation convention | Surrogate-selected gain | Defined surrogate loss, µm | Actual candidate discrepancy requirement, µm | Largest sampled surrogate/direct difference, µm |
|---|---:|---:|---:|---:|
| Disk / spatial-max components | 8.06297202 | 5.24685117 | 5.37983567 | 0.43354056 |
| Disk / regional-mean components | 9.50995731 | 6.02462625 | 6.18013536 | 0.43295265 |
| Line / spatial-max components | 7.32715933 | 5.59934156 | 5.69666515 | 0.41749944 |
| Line / regional-mean components | 8.16024598 | 6.07629139 | 6.24513544 | 0.43462117 |

The direct evaluation therefore changes the surrogate predictions measurably. The displayed direct discrepancy endpoints agree at this precision; exact rational lower and upper bounds are retained in the JSON. Agreement of displayed endpoints must not be interpreted as exact physical knowledge. The largest endpoint-observer interval width over all four candidate replays is approximately \(1.55\times10^{-6}\) µm, conditional on the stored finite-dimensional fields. This small observer uncertainty does not bound mechanical discretization error, experimental uncertainty or source-operator ambiguity.

### What the discrepancy interval proves

For one observed endpoint \(y\) and a certified stored-field prediction interval \([a,b]\), the sharp interval-wise residual bounds are

\[
 d^-=\max(a-y,y-b,0),\qquad
 d^+=\max(|a-y|,|b-y|).
\]

For all 68 endpoint constraints in one convention,

\[
 \delta^-=\max_{i,k}d^-_{ik}
 \ \le\ \max_{i,k}|p_{ik}-y_{ik}|
 \ \le\ \delta^+=\max_{i,k}d^+_{ik}.
\]

Thus any common tolerance smaller than \(\delta^-\) excludes this **specific candidate** under its declared observation contract; a tolerance at least \(\delta^+\) covers every point in the supplied intervals. These inequalities are direct consequences of the minimum and maximum of absolute distance on a closed interval. They neither assign a probability distribution to error nor prove that another gain or mechanical parameter choice cannot perform better.

The worst discrepancy in all four replays is regional-mean shortening. The per-channel maximum absolute requirements are:

| Convention | Maximum height, µm | Regional mean, µm | Apical endpoint, µm | Basal endpoint, µm |
|---|---:|---:|---:|---:|
| Disk / spatial-max | 4.45327 | 5.37984 | 5.11462 | 3.13650 |
| Disk / regional-mean | 5.14141 | 6.18014 | 5.89431 | 3.32252 |
| Line / spatial-max | 4.36937 | 5.69667 | 5.50916 | 3.19809 |
| Line / regional-mean | 4.50176 | 6.24514 | 5.95470 | 3.31484 |

The maximum-height, mean-height and component observations therefore constrain different features. Fitting one scalar height response would conceal the mean and component discrepancies exposed here. At the tested candidates, the model tends to predict regional means above many recorded means and basal motions below many recorded basal endpoints. This is evidence against interpreting these particular fixed-scenario candidates as quantitatively adequate descriptions of all four observables. It does not falsify the full material family or establish which missing mechanism or observation assumption causes the mismatch.

## 7. Visual evidence

The plots are generated exclusively from the completed direct receipt and the small-region calibration CSV. Lines joining directly solved strengths are display guides; they do not add intermediate solved states or continuous-response certificates. The source-column ordering and four retained operator conventions remain visible.

![Four directly replayed endpoint responses](../runs/joint_gain_extension_v2/joint-direct-endpoints.png)

![All 68 residuals for each observation convention](../runs/joint_gain_extension_v2/joint-direct-residuals.png)

![Surrogate and direct candidate discrepancy requirements](../runs/joint_gain_extension_v2/joint-direct-surrogate-comparison.png)

The figure hashes and generator/input hashes are retained in `evidence/joint_gain_extension_visuals.json`. All three exported figures were opened and visually inspected after generation.

## 8. Numerical and provenance receipts

The completed direct result has SHA256

    78c779112966adaa6819778637d3abf573afbc5af61ee74ea1ec28b8b5967dea

`evidence/joint_gain_extension_verification.json` records **34 passing aggregate checks**, together with per-state receipt checks. It verifies frozen sources and protocols, unchanged parent/seed artifacts, every common gain–dose association, actual binary64 input rounding, every field/verification hash and the last-verified-only continuation vector chain. Its all-state summary is:

- 68 verified stationary discrete fields, 68 accepted increments and zero rejected attempts;
- 144 triangles per field, hence **9,792 exact stored-map triangle certificates**;
- maximum scaled Q8 force residual \(1.86408\times10^{-10}\);
- maximum scaled Q12 cross-force residual \(1.86416\times10^{-10}\), against the declared \(10^{-7}\) threshold;
- exact rational whole-field Jacobian enclosures jointly contained in approximately \([0.955378588,1.026141882]\);
- positive local orientation and certified global embedding of every stored field.

The Jacobian interval is a conservative enclosure of the stored P2 maps, not a claim that either endpoint is attained. It certifies positive orientation, **not** exact incompressibility or a one-percent volume-error guarantee. Small stationarity residuals also do not establish uniqueness, energy-global optimality, stability or a PDE error bound. Existing higher-mesh comparisons already show appreciable shape sensitivity.

The independent all-cell surrogate audit is `evidence/piecewise_gain_minimax_v2_independent_review.json`; the separate independent full direct replay is documented in `JOINT-DIRECT-REVIEW.md` and `evidence/joint_direct_independent_review.json`. The independent direct review passes all **16 final checks**: every actual field is replayed under Q8/Q12, exact geometry and the endpoint observer, and all 272 residuals plus all 272 sampled surrogate/direct errors are reconstructed. It binds the final direct-result hash above. Its receipt SHA256 is `c37960e959b55e001b87ec2b5683d087b8419bc957436931399a93152f7efc41`. This is a distinct review, not replaced by the producer’s bookkeeping.

Reproducible final postprocessing, which refuses partial runs:

    python -X utf8 controlled_furrowing/runs/joint_gain_extension_v2/final_verify.py
    python -X utf8 controlled_furrowing/runs/joint_gain_extension_v2/make_visual.py

## 9. What has and has not advanced

This step completes an actual all-observable candidate replay rather than another single-height calibration. It replaces a visibly truncated gain selection with an independently certified interior **surrogate** selection and checks every selected input by a fresh stationary nonlinear solve. It also measures, rather than assumes away, the surrogate/direct disagreement.

The result does not establish a satisfactory biological fit. A single deterministic mechanical scenario with one common gain and the retained static observation proxies still requires roughly 5.4–6.3 µm of uniform discrepancy. The discrepancy may combine biological heterogeneity, missing dynamics, incorrect material/geometry assumptions, observation-reference or aggregation mismatch, and remaining numerical error. This experiment does not identify their separate contributions.

The next defensible inference step is to resolve the source observation contract and obtain numerical response enclosures across a declared parameter family, then characterize which parameters satisfy externally justified endpoint tolerances. Uniform nonlinear response/error control would be needed to upgrade the surrogate gain certificate into a physical gain-optimality result. Simply rerunning the same optimization on a denser table does not provide that proof.

No held-out larger-region outcomes were used to pick a convention or fit a parameter. No out-of-sample geometry-transfer claim, full \(\Theta_\delta\), continuum-accuracy certificate or novel general scientific principle is established here. The original physical and observation-alignment gates remain open.
