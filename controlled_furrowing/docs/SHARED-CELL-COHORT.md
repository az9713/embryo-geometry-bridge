# One shared parameter cell, seventeen roots, sixty-eight observations

**All seventeen source records now have individually audited uniform stationary-root and static-observation certificates on the same five-parameter cell.** The original fifteen provide 60 outputs; separately verified proof-center corrections close the F/K gaps and supply the remaining eight. No physical parameter, source dose or shared cell was changed to obtain those two proofs.

The source is Countryman et al., [*Endogenous OptoRhoGEFs reveal biophysical principles of epithelial tissue furrowing*](https://doi.org/10.1038/s41467-025-62483-6) (2025). The experiments use light-controlled recruitment of molecular regulators to perturb contraction in fly-embryo epithelial tissue. Our 17 small-region source entries supply activation doses and four deformation summaries. This model studies the contraction-to-shape part of that causal chain; it does not derive the molecular signaling from chemistry. Larger-region prediction remains a separate scientific gate.

For the selected stationary-root family, the maximum discrepancy from the 68 published response values is uniformly enclosed. A rigorously outward-rounded display is

\[
\boxed{\quad
5.18638714\le D(\theta)\le5.18638737\;\mu\mathrm m,
\qquad\theta\in C.
\quad}
\]

The [receipt's exact rational endpoints](../evidence/shared_cell_cohort_inference.json) are authoritative. Their numerical approximations are \(L\approx5.186387141964336\) and \(U\approx5.186387360305146\,\mu\mathrm m\); the **underlying exact enclosure** has width less than \(2.184\times10^{-7}\,\mu\mathrm m\). This is a pointwise uniform bound throughout \(C\), not a proof that either endpoint is attained, a sharp range, or a globally optimal calibration.

That narrow width concerns the frozen finite-dimensional mathematical target on a minute parameter cell. It is neither experimental precision nor continuum accuracy, and it does not bound the separate modeling, observation-law or PDE errors.

![All 68 conditional output enclosures](../runs/shared_cell_cohort_v1/shared-cell-68-observations.png)

White crosses are published source values. Cyan marks are the original fifteen records; orange marks are F/K, whose fresh root certificates use corrected numerical proof centers. The intervals are too narrow to resolve at this plot scale. They are mathematical bounds, not measurement error bars. The visible discrepancies did not disappear when the proofs became rigorous.

## The exact claim and its two tolerance queries

Let \(F_i(u_i,\theta)=0\) denote stationary force balance for source record \(i\), with its fixed dose \(I_i\) and active coefficient \(gI_i\). Let \(O_{ik}\) be the four stated static observations and \(y_{ik}\) the published values. Define

\[
D(\theta)=\max_{1\le i\le17,\;1\le k\le4}
|O_{ik}(u_i(\theta),\theta)-y_{ik}|.
\]

The root certificates select one locally unique stationary state in a declared nearby domain for each record. The joined theorem is

\[
\forall\theta\in C\quad
\exists!\,(u_1,\ldots,u_{17})\in\prod_i\mathcal D_i(\theta):
\quad F_i(u_i,\theta)=0,\quad
L\le D(\theta)\le U.
\]

Uniqueness is restricted to that product of nearby domains. It is not uniqueness among every possible mechanical state. The seventeen states are independent stationary responses to seventeen inputs; a stable common loading trajectory or experimental time path has not been established.

For clarity, define the restricted stationary compatibility set

\[
\Theta^{\rm stat}_\delta=
\left\{\theta:
\exists(u_1,\ldots,u_{17})\text{ admissible stationary states},\quad
\max_{i,k}|O_{ik}-y_{ik}|\le\delta\right\}.
\]

This definition uses the fixed geometry, normalization, numerical energy and conditional static observation contract of this calculation. It includes no energetic-stability or source-time requirement.

| Declared common allowance | What follows |
|---|---|
| \(\delta=5.1\,\mu\mathrm m<L\) | The **selected root tuple** fails throughout \(C\). Other branches remain unresolved; this does not prove \(C\cap\Theta^{\rm stat}_{5.1}=\varnothing\). |
| \(\delta=5.2\,\mu\mathrm m>U\) | Every parameter in \(C\) has the certified stationary tuple satisfying all 68 constraints. Thus \(C\subseteq\Theta^{\rm stat}_{5.2}\): a nonempty inner **stationary-tuple witness** under this restricted contract. |

These are illustrative conditional queries made after the M threshold was known, as the [frozen combination protocol](../runs/shared_cell_cohort_v1/protocol.json) states. They are not estimated measurement noise, a prospective hypothesis test, or evidence that a \(5.2\,\mu\mathrm m\) discrepancy is biologically acceptable.

![Declared allowance and certified residual bounds](../runs/shared_cell_cohort_v1/shared-cell-residual-threshold.png)

For \(\delta<L\), this selected family is excluded. For \(\delta\ge U\), it supplies a stationary compatibility witness throughout \(C\). Between the exact endpoints, these bounds do not decide the question. No claim is made that every value between the endpoints occurs.

## Why the shared cell matters—and how small it is

The same parameter vector

\[
\theta=(K_a,K_b,T_a,T_b,g)
\]

is used for **every** record: normalized apical/basal area moduli, apical/basal surface tensions and common activation gain. The four positive coordinates have relative half-width \(10^{-8}\); \(T_a\in[0,2\times10^{-10}]\). The cell has nonzero width in every coordinate. It is a tiny computational box around one candidate, not an inferred physiological prior, broad robustness region or complete material/geometry identified set.

The finite-height axisymmetric model, 12×6 quadratic mesh, 611 free coordinates, fixed \(\mu=1\), finite volumetric coefficient \(\kappa=5000\), geometry and binary Duffy-8/surface-12 integration cache remain fixed. The exact source-dose law is \(f_i=gI_i\); archived rounded floating solver drives are not substituted for that law.

The common contract and cell hashes are respectively `2498493…942f82` and `64f93fe…394dbf`. The [combined receipt](../evidence/shared_cell_cohort_inference.json), SHA-256 `1f8c4dc1ad8924f8fda13e6748d5e73aac426b5771b02404803d27c1955aa677`, binds all 17 proof pairs. Its protocol locks 71 input files before combination. There is no separate fitted material vector or chosen cell for F or K.

## How the theorem is obtained

### 1. Retain the original fifteen proofs

The [first independent cohort review](../evidence/cohort_uniform_roots_independent_review.json) verified 15 records on \(C\) and retained F/K as unresolved. For each successful record, interval arithmetic established a predictor-centered root tube

\[
u_i(\theta)\in u_{i,0}+S_i(\theta-\theta_0)+[-r_i,r_i]^{611}
\]

using the Banach conditions \(Z_i<1\) and \(\eta_i+Z_i r_i<r_i\). The [first observation application review](../evidence/cohort_static_observers_independent_review.json) then verified 60 output intervals. Those files are unchanged. The [60/8 checkpoint](COHORT-INFERENCE-CHECKPOINT.md) remains the historical record of that partial phase.

### 2. Correct a proof center, not a physical parameter

Larger error cubes had failed to certify F/K on the same \(C\). The [bound diagnostics](COHORT-BOUND-DIAGNOSTICS.md) showed that the original-center residual dominated, while the parameter-only derivative bounds were modest. A separately frozen calculation therefore proposed one numerical correction

\[
u_{i,1}=\operatorname{fl}_{64}\!\left(
u_{i,0}-\operatorname{mid}_{\rm Arb}[R_iF_i(u_{i,0},\theta_0)]
\right).
\]

It retained one nearest-even binary64 rounding of the exact subtraction. Original saved fields remained immutable. All cached arrays except `u` were verified unchanged, **including \(R_i\)**; the binary predictor \(S_i\), parameters, exact gain-times-dose law, geometry, quadrature and cell were also unchanged. Auxiliary center files explicitly say they are numerical proof centers, not experimental data or converged nonlinear-solver states.

Fresh interval proofs succeeded for both records with the normal corrector radius \(r_i=10^{-11}\). The [independent corrected-root review](../evidence/cohort_recentered_roots_independent_review.json), SHA-256 `9eb54489…111ee34f`, verified those new certificates. A numerical correction alone was never treated as a proof. The earlier failures remain valid records of an insufficient enclosure, rather than physical nonexistence results.

### 3. Map the new root tubes back to the original graphs

For F/K, observations still use the **original saved seed graphs**. A new adapter handles the center offset explicitly. Its coordinate bound is

\[
\epsilon_{ik}=
|X_{ik}+u_{i,1,k}-x^{\rm seed}_{ik}|
+\max_{\theta\in C}|(S_i(\theta-\theta_0))_k|+r_i
\]

for a free coordinate, with only the exact gap on constrained coordinates. All quantities are exact rationals of frozen binary constants and exact cell endpoints. The center translation is not added to the Banach radius, its residual, or its derivative test. The original strict same-center observer was preserved and still rejects a new center passed off as the old numerical state.

Each complete nodal box then passes the reviewed whole-body embedding certificate. Exact Bernstein bounds, axis conditions and strong monotonicity apply over entire elements, not only sampled quadrature points. The fixed radial endpoint coordinates and coverage of the same current-radius ROI are verified separately.

### 4. Transfer geometry errors to the four source-named outputs

On the same fixed current radius, define \(a=A_0-A\), \(b=B-B_0\) and \(s=a+b\). The primary observations are

\[
\left(\sup_R s,\quad
\int_0^{R_*}s(R)\frac{2R}{R_*^2}\,dR,\quad
\sup_R a,\quad\sup_R b\right).
\]

For a quadratic boundary trace with coordinate error \(\epsilon\), interpolation gives a value bound \(5\epsilon/4\) and derivative perturbation bound \(8\epsilon\). A strict radial margin and the base graph slope bound give separate Eulerian errors \(e_A,e_B\). Spatial maxima and normalized nonnegative means are 1-Lipschitz in the uniform norm, so the four output errors are

\[
(e_A+e_B,\;e_A+e_B,\;e_A,\;e_B).
\]

Exact base graph intervals are expanded by those errors. Height is evaluated between the two surfaces at the **same current physical radius**; it is not a material-fibre length or a sum of separately located component maxima.

The [F/K application review](../evidence/recentered_observers_independent_review.json), SHA-256 `f0e99c43…60fcd97`, rebuilt the offset arithmetic, hulls, full normal radii, raw boundary connectivity, transport constants and eight constraints. Whole-body geometry and base-observer checks are disclosed replays of previously reviewed components, not independently rewritten implementations. A [metadata-only scope correction](../runs/cohort_recentered_observers_v1/scope_correction.json) explicitly retained F/K-only coverage as eight channels; the global join counts actual 17×4 rows and does not trust the old helper's misleading all-68 flag.

### 5. Combine all 68 exact residual inequalities

For each uniform prediction interval \([l_{ik},h_{ik}]\),

\[
d^-_{ik}=\max(l_{ik}-y_{ik},y_{ik}-h_{ik},0),\qquad
d^+_{ik}=\max(|l_{ik}-y_{ik}|,|h_{ik}-y_{ik}|).
\]

Every selected root satisfies \(d^-_{ik}\le|O_{ik}-y_{ik}|\le d^+_{ik}\). Therefore

\[
L:=\max_{i,k}d^-_{ik}\le D(\theta)
\le\max_{i,k}d^+_{ik}=:U,
\qquad\forall\theta\in C.
\]

This step needs all 68 channels for the sufficient upper bound. The two newly audited F/K records now supply the eight channels missing from the historical checkpoint. The M/apical constraint still determines the reported enclosure; completing the proof did not optimize away its discrepancy.

The final [independent combination review](../evidence/shared_cell_cohort_independent_review.json), SHA-256 `2f7da9f2bee4c1b7e4fdd7f3a88ed90dd914fb218269b3d65bf12bb8a6ab9114`, passed 180 checks. It verified all 71 locked inputs, 17 proof pairs, 68 constraint rows, residual bounds, conditional queries, source controls and retained archives. The [review explanation](SHARED-CELL-COHORT-REVIEW.md) records its scope. This is separate from the independent root and application proofs above; no additional equilibrium calculation or parameter optimization is performed by the combination.

## What has advanced, and what has not

The advance is a **complete, parameter-uniform, finite-dimensional stationary compatibility witness or exclusion test on one shared cell**, with individually validated geometry and observations for every record. It closes the two previously missing root/observation proofs and supplies an all-68 sufficient bound under the declared static contract.

This remains one local part of Gate 3 in the [original five-gate plan](../PLAN.md). It does not infer the full material/geometry compatible set, exclude other branches when this one fails, validate energetic stability or a common physical trajectory, resolve source temporal peaks/reference/ROI/record-pairing ambiguities, bound quadrature or continuum PDE errors, or establish a larger-region prediction. A declared \(5.2\,\mu\mathrm m\) allowance is not a measurement-error estimate. The old fixed-material weak-load energy-norm gain family remains a distinct result and is not used as a substitute for these cohort proofs.

No conceptual novelty or publication priority follows from employing Banach bounds, exact graph observers, proof-center refinement or interval residual arithmetic. The project has become substantially more assessable: claims of compatibility, exclusion and missing evidence are now tied to explicit quantifiers and retained receipts. The wider scientific question—whether a source-matched mechanically adequate model predicts embryonic geometry beyond its calibration conditions—remains open.
