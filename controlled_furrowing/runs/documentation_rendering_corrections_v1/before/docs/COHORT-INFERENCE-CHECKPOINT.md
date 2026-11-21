# Cohort inference checkpoint: 60 verified outputs, eight unresolved

**Fifteen of the seventeen small-region records now have independently verified uniform stationary-root and static-observation bounds on the same five-parameter cell.** This yields 60 verified conditional output intervals. Records F and K still lack the required uniform root certificates, so their eight output channels remain unresolved. There is no complete 17-record root tuple, no all-68 sufficient fit bound, and no completion of the original five-gate project.

This is progress toward an assessable chemistry–mechanics–geometry bridge: a prescribed optogenetic contraction input is connected to a mechanically consistent tissue deformation and to explicit experimental observables. It does not yet explain how an embryo generates its entire body plan, infer a molecular regulatory mechanism, or validate a prediction of the larger activated tissue regions.

![Verified partial cohort: 60 enclosures and eight unresolved channels](../runs/cohort_static_observers_v1/composition_final/cohort-static-observations.png)

White crosses are published values; cyan marks are uniform model intervals across the declared computational cell. These intervals are too narrow to resolve at the plot scale. They are not experimental error bars. F and K are displayed as unresolved rather than supplied with unverified numerical predictions.

## What was actually verified

The fixed model is a finite-height, axisymmetric, nonlinear elastic tissue with surface mechanics and an imposed contraction drive. Its current implementation has 611 free displacement coordinates on a 12×6 quadratic finite-element mesh. The exact mathematical target is the finite-dimensional energy using its frozen binary integration cache, with 9,504 local volume/surface contributions. This is a precise discrete target; it is not a continuum PDE error certificate.

Every record uses the **same** closed cell

\[
C\subset(K_a,K_b,T_a,T_b,g),
\]

where (K_a,K_b) are apical/basal area moduli, (T_a,T_b) are surface tensions and (g) is the shared activation gain. Its relative half-width is (10^{-8}) in the four positive coordinates; (T_a\in[0,2\times10^{-10}]). Geometry, \(\mu=1\), finite volumetric coefficient \(\kappa=5000\), and the integration policy remain fixed. This tiny cell is computational support around one candidate, not a physiological prior or the full compatible material/geometry set.

For canonical source dose (I_i), the drive is exactly (f_i=gI_i). The earlier floating solver's rounded drive is recorded separately. A proof must handle the exact gain-times-dose target rather than substituting that rounded input.

For each of the 15 successful records, the validated statement is

\[
\forall\theta\in C\quad
\exists!\,u_i(\theta)\in
u_i^0+S_i(\theta-\theta_0)+[-r_i,r_i]^{611},
\qquad F_i(u_i(\theta),\theta)=0.
\]

Uniqueness applies only inside the declared nearby domain. It does not rule out other stationary states, establish energetic stability, or connect different records along a common loading or experimental time trajectory.

The [independent root review](../evidence/cohort_uniform_roots_independent_review.json), SHA-256 `67c9d5195c36fc11814d08517070fe53ea0f10f67d45911df517bd900e098297`, reviewed all 17 records at 192-bit interval precision: 15 passed and F/K remained unresolved. Its overall status is deliberately `PARTIAL_UNRESOLVED`, rather than a global PASS that would conceal the two missing proofs. M's earlier immutable independent proof is explicitly reused.

For every successful record, the observation composition then checks the entire root-tube body geometry, fixed radial endpoints, coverage of the common current-radius ROI and strict radial graph monotonicity. It transfers the full corrector cube, exact affine predictor hull and saved-node serialization gap into four output intervals:

\[
\sup_R s(R),\qquad
\int_0^{R_*}s(R)\frac{2R}{R_*^2}\,dR,\qquad
\sup_R a(R),\qquad\sup_R b(R),
\]

with (a=A_0-A), (b=B-B_0), (s=a+b), all evaluated at the **same current physical radius**. Height shortening is not a material-fibre length or the sum of separately located component maxima.

The [independent application review](../evidence/cohort_static_observers_independent_review.json), SHA-256 `5e71814a8a266ff309f990dfcd07b7ba8078abee39c5db8ac10ea04aaa2af8ad`, has status `PASS_PARTIAL_APPLICATION`. Its 18 finalization checks bind the 15 passing independent root reviews to the corresponding observer outputs and retain the two unresolved records. The underlying [arithmetic replay](../runs/cohort_observer_review_v1/arithmetic.json) passed 278 checks. Predictor hulls, serialization gaps, full radii, boundary extraction, transfer constants and residual arithmetic were independently reconstructed. Whole-body geometry and the base static observer were **disclosed replays of previously reviewed components**, not independently rewritten implementations.

The [observation proof and receipts](COHORT-STATIC-OBSERVERS.md) give the detailed transfer argument. No additional geometry or observer failure occurred among the 15 root-certified records. The final coverage is therefore exactly 60 outputs plus eight missing channels.

## What the residual bound means

If (y_{ik}) is a published response and the uniform model output is enclosed in ([l_{ik},h_{ik}]), then

\[
\max(l_{ik}-y_{ik},y_{ik}-h_{ik},0)
\le |O_{ik}(\theta)-y_{ik}|.
\]

Taking the maximum over the 60 covered channels gives the exact necessary lower bound retained in the [terminal result](../runs/cohort_static_observers_v1/composition_final/results.json). A conservative displayed consequence is

\[
D(\theta):=\max_{i,k}|O_{ik}(\theta)-y_{ik}|
\ge 5.18638714\,\mu\mathrm m
\]

for any hypothetical complete static root tuple extending the covered selected branches. The controlling constraint is M's apical response.

This rules out such an extension for a **declared** common discrepancy allowance below that value. It does not establish a measured noise level or imply that (5.19\,\mu\mathrm m) is biologically acceptable. It does not reject every branch or every material/geometry choice. Source-reference and aggregation assumptions remain conditional.

Most importantly, the missing F/K channels prevent any all-68 sufficient upper bound. Even a large allowance cannot be certified to fit the full cohort from this packet, because root existence and observations for those records have not been proved. A numerical point prediction or a passing matrix calculation cannot fill this gap.

## Why F and K failed the present proof

The Banach test requires both

\[
Z<1,\qquad \eta+Zr<r.
\]

The initial radii were approximately (1.73080\times10^{-10}) for F and (1.04523\times10^{-10}) for K. Their derivative bounds were below one, but the self-map inequality failed. The separately retained [radius supplement](../runs/cohort_radius_supplement_v1/results.json) tried 2×, 4× and 8× those radii on the **same** parameter cell, with the same physical constants, predictors, caches and preconditioners. All six tests remained unresolved; the larger cubes made the conservative derivative bound exceed one. No record-specific cell shrink was used.

The [fresh diagnostic decomposition](COHORT-BOUND-DIAGNOSTICS.md) explains the immediate obstruction:

| Validated bound | F | K |
|---|---:|---:|
| Original-center preconditioned residual | \(8.58449\times10^{-11}\) | \(5.00617\times10^{-11}\) |
| Predictor-only derivative bound, zero error radius | 0.0376960 | 0.0827620 |
| Affine tangent-defect contribution over the cell | \(2.890\times10^{-21}\) | \(5.336\times10^{-21}\) |

The original-center residual dominates the predictor residual. The parameter hull alone gives much smaller derivative bounds than the full error cubes. These are diagnoses of this sufficient proof method. They do **not** prove root absence, instability, physical collapse or biological incompatibility. Nor is a preconditioned residual itself a certified distance to an unknown root.

## Next bounded computation: refine a proof center, then prove again

The proposed next step is one separately declared correction of the **auxiliary numerical proof center**, not a change to tissue parameters, dose, data, cell or energy. With the same stored preconditioner (R), the proposal is

\[
u_1=\operatorname{fl}_{64}\!\left(u_0-
\operatorname{mid}_{\rm Arb}[R F(u_0,\theta_0)]\right).
\]

The intended benefit is a smaller residual around the same exact physical target. Success is not assumed. The new center is a frozen binary proposal and must receive a **fresh** interval residual, derivative, analytic-domain and self-map proof on the unchanged cell (C). Keeping the old (R) and (S) does not make the old certificate valid at (u_1).

The source/preflight for this correction are being prepared; **no corrected-center result is claimed here**. The original saved numerical field and all unsuccessful proof attempts remain immutable. Metadata must separate `original_seed_state` from the auxiliary `proof_center` and retain the exact translation (u_1-u_0). The correction is not an experimental observation or a new physiological trajectory.

If a corrected-center root certificate succeeds, observation transfer needs a separately reviewed center-offset adapter. The present strict same-center observer should continue rejecting that new cache. Relative to the original saved nodes, the new coordinate error must explicitly include

\[
|X+u_1-x_{\rm saved}|+
\max_{\theta\in C}|S(\theta-\theta_0)|+r_1.
\]

It must then recheck whole-body geometry, radial graph/ROI conditions and all four observations. A changed proof center cannot silently become a changed calibration datum or erase the original failures. Even a successful F/K refinement would finish only this local discrete-root/observer component of Gate 3.

## A separate successful family must not be conflated with this one

The earlier [energy-norm gain family](ENERGY-GAIN-CELL.md) has an independently reviewed uniform stationary-root/local-minimum certificate for an **original-material weak-load case**, with only relative gain varying over \([-0.001,0.001]\). Its [independent review](../evidence/energy_gain_cell_independent_review.json) passed 71 checks; a separate [whole-body geometry composition review](../evidence/energy_gain_geometry_independent_review.json) passed 1,858 checks.

That is an existing result in this project, not a newly published theorem and not a substitute for F/K. It fixes material coefficients, uses a different center and a different energy-norm domain, and does not cover all 17 new-material source doses or their 68 observations. Its local-minimum conclusion must not be transferred to the present stationary-root boxes. It supplies a verified methodological precedent, not a completed cohort result or a novelty claim.

## The original five gates remain the project

The [original execution contract](../PLAN.md) still asks whether one mechanically adequate model, calibrated only on the small-region data, predicts larger-region deformation without condition-specific refitting.

| Original gate | Current position and remaining work |
|---|---|
| 1. Exact originality audit | Executed. Scalar-transfer novelty was rejected; originality of the specific stronger spatial-transfer claim remains unconfirmed. Standard validation machinery is not itself a new concept. |
| 2. Actual data feasibility | Inventory executed: 17 calibration records, nine larger-region scalar records and 31 larger-region bending labels. Individual acquisition/geometry/classification mapping and aspects of the source observation law remain unresolved. |
| 3. Adequate model and identifiability | In progress. This checkpoint verifies 15 local root/observer families. It does not supply all 17 roots, the full compatible material/geometry/branch set, a source-justified discrepancy law, or uniform integration/continuum error bounds. |
| 4. Frozen larger-region prediction | A conditional earlier version was executed, but the original shape-transfer target remains incomplete. A revised source-matched prediction must carry the unresolved parameter/observation and numerical uncertainties before scoring. No new successful larger-region prediction follows from this checkpoint. |
| 5. Contribution assessment | Earlier assessments found no demonstrated predictive increment or confirmed novel concept. A revised conclusion must compare against simpler models, source limits, numerical uncertainty and the closest research after the stronger prediction is actually evaluated. |

In parallel with repairing the two numerical proofs, the scientific workflow still needs source clarification for temporal peaks, component references, ROI transport/weighting and record linkage; a defensible discrepancy contract; and calibration-set inference that retains other parameters, geometry choices, branches and tails. Contacting authors or collecting new data is a separate action, not something this checkpoint claims has happened.

The useful lesson is precise: an approximate nonlinear field can now be connected to a verified parameter-uniform body and to experimentally named observables, with explicit failure propagation. That makes the project more auditable. It does not turn an unresolved physical inference problem into a completed fit, a continuum theorem, or a novel explanation of embryonic shape.
