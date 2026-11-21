# Independent review of frozen spatial experiment v1

**Arithmetic, source hashes and the declared computational data flow pass this review. Successful biological shape prediction is not established.** This reviews the original stress-free, finite-height spatial experiment only; the prestress extension is separate ongoing work. It updates the earlier interim implementation status: an actual coupled discrete solver and93frozen geometry predictions now exist. Their existence does not close the original five-gate study.

This internal AI review uses separate standard-library/NumPy calculations and does not import or rerun the project solver or slice-observer implementation. The reviewer authored the data extraction previously; the present arithmetic is independent of the model/scorer implementation, not external peer review.

## What was independently replayed

* All **93 prediction-result SHA-256 hashes**:31larger-region records in each of3predeclared basal-material scenarios.
* The6source/protocol/input hashes in every frozen manifest, matching record identities, source dosage/diameter inputs, and the scored-response hashes.
* The33response-table entries per branch and all51direct small-region calibration cases. A separate `bisect` piecewise-linear interpolator and scalar summation reproduce each501-point gain-grid minimizer. Direct predictions and the declared0.05µm selected-calibration interpolation tolerance are independently checked.
* All9scalar prediction/measurement pairs per branch, their residuals, RMSE, MAE and signed bias.
* All31possible-label sets per branch at assumed surface errors0,0.05,0.2µm, using independently written gap algebra.
* Local volume-Jacobian extrema and global revolved volume from retained reference/current nodes and triangles, without using the stored Jacobian values as inputs.
* The training-only apical/basal motion-partition diagnostic and the simpler calibration-only linear comparator.

The exact verification program is retained inside [the independent receipt](../evidence/spatial_v1_independent_review.json), together with snapshot and intermediate-file hashes. To replay from the embryo-project root:

```powershell
python -c "import json;exec(json.load(open('controlled_furrowing/evidence/spatial_v1_independent_review.json'))['verification_code'])"
```

The program refuses to silently accept changed audited inputs. A replay-code bookkeeping correction is documented in the receipt; its previous code is retained. No model, protocol, source measurement or frozen prediction was changed during this review.

## Data-flow finding

The preparation stage reads the public source table to remove outcomes, producing a31-record file with exactly `record_id`, `reported_diameter_um`, and `activation_magnitude`. Calibration reads the17small-region response records and model response table. Prediction reads the frozen calibration and larger-region input descriptors; neither callable path reads larger-region shortening or bending labels. The scorer checks every branch and prediction hash before loading responses.

This supports the declared **computational withholding**. It is a source/call-path and artifact-dependency audit, not an independent operating-system I/O trace or immutable proof of execution chronology. Public outcomes were already known to the analysts. Hashes identify bytes; they do not establish biological blinding or scientific truth.

All three predeclared material branches remain in the results. Their fitted gains—8.66,8.88,9.50—are reproduced from calibration only and are not at the declared gain-grid endpoints. These are minima of the specified interpolated finite grid, not proofs of an exact nonlinear or continuous-parameter global optimum. The material scenarios are assumptions, not a complete calibration identified set.

## Scalar results

| Basal area-modulus scenario | Fitted gain | Large-region RMSE (µm) | MAE (µm) | Prediction-minus-observation bias (µm) |
|---:|---:|---:|---:|---:|
| 0.2 | 8.66 | 3.088137 | 2.552934 | +2.446180 |
| 2 | 8.88 | 3.047690 | 2.527290 | +2.376537 |
| 20 | 9.50 | 3.150300 | 2.601510 | +2.482836 |

The calibration-only linear comparator gives **RMSE2.227042µm** on the same9large-region observations. The zero-shortening baseline gives6.892817µm. Thus all three coupled v1 scenarios outperform the zero response empirically but have larger errors than the simple linear comparator. This is descriptive performance, not a significance or equivalence claim.

The model endpoint is maximum shortening of paired initial material columns. The source measurement subtracts a cellularization reference and compares apical/basal surfaces in image coordinates. Those are not automatically the same observation operator when surfaces acquire different radial displacements. The numerical comparison is therefore explicitly a proxy diagnostic.

The independent small-region partition check also finds **zero of17predicted cases apical-dominant in every branch**, whereas the source's separate apical/basal maxima are apical-dominant in16of17. Those source maxima can occur at different times/locations, and the model is static, so this is a warning about the chosen peak-equilibrium proxy rather than a universal biological falsification. Nevertheless a scalar-height fit does not validate the predicted partition of motion.

## Why zero label exclusions proves no binary accuracy here

For every one of the31records in every branch, the unknown-slice possible-label set is **{0,1}**, even at zero assumed surface error. Both published labels are admitted. All31sets remain ambiguous at0.05and0.2µm assumptions.

Therefore the reported **zero excluded labels is vacuous compatibility**, not31correct binary predictions or100%predictive accuracy. There is no discriminative binary prediction in these sets. The lack of actual acquisition slice/time correspondence remains even if the discrete mechanical geometry were exact.

The0/0.05/0.2µm numbers are user-declared error assumptions for a conditional operator, not measured image errors, numerical error estimates or certified continuum bounds. Enlarging such assumptions cannot provide missing optical observations or certify statistical coverage.

The independently checked gap is

\[
G=\min\{z_a(0),z_a(R)\}
-\max\{z_b(0),z_b(R),\min_r z_a(r)\}.
\]

For the stated graph/body assumptions, the nominal unknown-slice rule returns{0}whenG≤0and{0,1}otherwise. A uniform surface-error assumptionε enlarges the gap by at most2ε. This finite geometry calculation is distinct from proving thatεactually bounds the finite-element error or that the source used the corresponding ideal slice.

## Local volume loss is substantial despite small global change

| Basal scenario | Minimum local J in larger predictions | Maximum local volume loss | Largest absolute global volume change |
|---:|---:|---:|---:|
| 0.2 | 0.812454 | 18.75% | 0.389% |
| 2 | 0.809285 | 19.07% | 0.388% |
| 20 | 0.797252 | 20.27% | 0.415% |

These values were independently reconstructed. Small global volume change can coexist with local compression and compensating expansion. Consequently these results do **not** validate biologically near-incompressible behavior. They remain solutions of the declared finite-compressibility penalty model.

For each linear meridional triangle, the meridional deformation gradient is constant. The hoop stretch isr/R, a weighted average of non-axis nodal radius ratios. This gives the element's Jacobian extrema from nodal ratios times the meridional determinant. The revolved volume is

\[
V=\sum_T2\pi |T|\,\det F_{\rm mer}\,
\frac{r_1+r_2+r_3}{3}.
\]

The receipt verifies these identities against all retained geometries in floating point. They justify the reported discrete geometric diagnostics; they are not interval-arithmetic certificates, continuum-error bounds, proofs of a global energy minimizer, or proofs of global injectivity/self-contact avoidance. Positive localJand small force residuals do not establish those stronger statements.

## Original five-gate assessment after v1

| Gate | What v1 adds | What remains |
|---|---|---|
| 1. Originality | Bounded primary-source comparison exists | No novelty established |
| 2. Data feasibility | Actual numerical cohorts and public example images audited | Specimen geometry, slice/time rules and calibrated uncertainties remain incomplete |
| 3. Model and identifiability | Genuine coupled finite-height **discrete** implementation, calibration and checks now exist | Small-motion partition warning, no complete material identified set, unvalidated specimen/constitutive assumptions and continuum accuracy |
| 4. Frozen prediction | All93conditional spatial outputs frozen and independently hash/arithmetic checked | Original uncertainty/observation requirements not fully satisfied; binary label sets uninformative |
| 5. Contribution | Reproducible stress test exposes deficiencies of this v1 model | No established successful biological shape transfer or novel concept/theorem |

The original study is therefore **not complete merely because v1 executes correctly**. The new prestress model must retain this failed/limited branch, use an explicit versioned protocol and avoid choosing material or observation parameters to match the already inspected larger responses. Its performance and scientific status require their own audit.
