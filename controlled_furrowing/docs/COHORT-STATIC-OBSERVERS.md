# Uniform static observation of the 17-record stationary-root family

The independently reviewed partial composition encloses **60 outputs for 15 records**, with **eight channels for records F and K unresolved**, over one shared five-parameter cell. All 68 requested source constraints remain in the covered/missing ledger. There is a rigorous necessary discrepancy lower bound conditional on the underlying root certificates; there is no sufficient upper bound or all-68 IN statement. This is a partial checkpoint, not completion of the original inference goal.

The composition performs no equilibrium optimization. The mechanical root proofs and observation composition retain separate independent reviews. Passing the partial application review does not supply the missing eight channels.

![Uniform static outputs versus published source values](../runs/cohort_static_observers_v1/composition_final/cohort-static-observations.png)

White crosses are published values; cyan marks are the conditional uniform model intervals. Their widths are too small to resolve on this scale. They are mathematical root/parameter-cell enclosures, not measurement error bars. Source-column labels F and K remain explicitly unresolved.

## Actual partial result

The [terminal result](../runs/cohort_static_observers_v1/composition_final/results.json), SHA-256 `ffa23af7660364a3d4d40460024ba3565c086cef1662cb2cb55365bb8ed27e0c`, retains every record. All 15 root-certified records passed the additional whole-body geometry, fixed-endpoint/current-ROI, radial graph and endpoint-composition checks. No additional records were lost at these stages.

The largest certified necessary residual is the M/apical constraint. Its exact rational lower endpoint is stored in the result; an outward conservative decimal statement is

\[
D(\theta)\ge 5.18638714\,\mu\mathrm m
\]

for any hypothetical complete static root tuple extending the covered selected branches. This bound already excludes such an extension for a declared common allowance below that value. It does not assert that the unresolved F and K roots exist, exclude other branches, establish a biological noise level, or reject every model/material parameter choice.

There is **no** all-68 sufficient upper bound: eight missing constraints remain missing. The v1 dose-encoding failure, immutable v2 completion and unsuccessful supplemental F/K radius attempts are preserved through the [provider selection](../runs/cohort_static_observers_v1/provider_selection.json) and its lineage. Their failed sufficient proofs are not evidence of nonexistence or physical collapse. The next mathematical gap is to obtain sound uniform root bounds for these two records before any all-17 tuple statement.

All 15 base maximum-height intervals happened to have exact width zero under the retained exact rational adaptive signed-supremum certificates; all requested reporting widths were therefore achieved. The general code still preserves wider valid enclosures if a future reporting budget is reached.

The common contract is [model_observation_contract.json](../runs/cohort_static_observers_v1/model_observation_contract.json), SHA-256 `2498493bf50049295eb2ad0d2360b1defaf387ee8db8828ceb576f3d7e942f82`. It contains the fixed geometry, material-coordinate conventions, integration policy, exact shared cell, physical length, primary observer, source pairing, and limits. Per-record dose, saved state, predictor, root domain, and provenance are deliberately separate.

## Four observations at a common current radius

Let \(A(R)\) and \(B(R)\) be the current apical and basal graph heights. The flat reference is \(A_0=33.5\,\mu\mathrm m\), \(B_0=0\). On the same fixed registered current-radius disk \(0\le R\le R_*\), define

\[
a(R)=A_0-A(R),\qquad b(R)=B(R)-B_0,\qquad s(R)=a(R)+b(R).
\]

The primary four-component output is

\[
\mathcal O(A,B)=\left(
\sup_R s(R),\quad
\int_0^{R_*}s(R)\frac{2R}{R_*^2}\,dR,\quad
\sup_R a(R),\quad\sup_R b(R)
\right).
\]

Height is measured between the two surfaces at the **same current physical radius**. A material fibre connecting two corresponding mesh nodes generally has a different observable. Separate apical and basal maxima are not added to obtain maximum shortening. All outputs are signed; there is no local positive-part clipping.

These are explicit conditional choices, not a resolution of the source's registration, component-reference, disk-versus-section weighting, or temporal-aggregation ambiguities. The source's reported peak quantities are compared with one static endpoint only under the stated endpoint-as-peak proxy. Published source-column order supplies the component pairing; it is not an independently verified specimen-ID crosswalk. See [the source observation audit](JOINT-OBSERVATION-CONTRACT.md).

## From a uniform root tube to an exact nodal box

For record \(i\), the independent mechanical proof must establish, for every parameter \(\theta\in C\), a unique stationary root in

\[
u_i(\theta)\in u_i^0+S_i(\theta-\theta_0)+[-r_i,r_i]^{n_i}.
\]

The proof uses the **full accepted corrector cube** \(r_i\), not the smaller a posteriori distance estimate. A strict Banach self-map/contraction bound is checked again with exact rational endpoints:

\[
0\le Z_i<1,\qquad \eta_i+Z_i r_i<r_i.
\]

Every binary64 predictor coefficient is treated as its exact binary rational. For each free coordinate \(k\), interval linear arithmetic gives the exact range

\[
[p^-_{ik},p^+_{ik}]
=\sum_{j=1}^{5}S_{ikj}[C_j^--\theta_{0j},C_j^+-\theta_{0j}].
\]

The serialized saved node is not silently identified with the exact reference-plus-displacement centre. If \(x_{ik}^{\rm save}\) denotes that saved coordinate and \(X_{ik}+u^0_{ik}\) its exact arithmetic centre, the nodal coordinate error is bounded by

\[
\epsilon_{ik}^{\rm norm}=
|X_{ik}+u^0_{ik}-x_{ik}^{\rm save}|
+\max(|p^-_{ik}|,|p^+_{ik}|)+r_i
\]

for a free coordinate. For a fixed coordinate only the serialization gap remains. This records every term instead of estimating error from force residuals, mesh differences, or a posteriori sampling.

The complete nodal box is passed through the exact whole-body geometry verifier: correlated Bernstein determinant/hoop bounds, axis conditions, and strong monotonicity imply an admissible embedded body for the entire box. A failed sufficient condition leaves that record unresolved; it is not evidence that the physical body fails.

## Boundary-specific transfer to the Eulerian graphs

Let \(\epsilon_A\) and \(\epsilon_B\) be \(33.5\,\mu\mathrm m\) times the largest coordinate radius on each corresponding surface's boundary nodes. The quadratic interpolation value bound is \(5\epsilon/4\), and the derivative perturbation bound is \(8\epsilon\) on the unit edge parameter. Strict radial monotonicity is therefore verified separately for each surface.

If \(L_A,L_B\) bound the saved graph slopes, then

\[
\|A_\theta-A_{\rm save}\|_\infty\le
e_A=\frac54\epsilon_A(1+L_A),\qquad
\|B_\theta-B_{\rm save}\|_\infty\le
e_B=\frac54\epsilon_B(1+L_B).
\]

This comparison is on a common current-radius domain. The code verifies that the axis and outer radial endpoint coordinates are fixed DOFs with **exactly zero** predictor/corrector/serialization radii. It then checks that both graph domains cover the prescribed positive-width ROI. It does not infer fixed endpoints from a flag passed to the transfer helper.

The map \(f\mapsto\sup f\) is 1-Lipschitz in the uniform norm. A normalized nonnegative weighted mean also has operator norm one. Consequently the four endpoint error bounds are

\[
(e_A+e_B,\quad e_A+e_B,\quad e_A,\quad e_B).
\]

The base-state endpoint intervals come from exact rational quadratic graph inversion and integration. These intervals are expanded by the appropriate channel-specific transfer error. Both the finite base-observer enclosure and the root-tube transfer error are retained.

The base maximum computation requests a reporting width of \(10^{-7}\,\mu\mathrm m\), with an 8192-cell budget. A certified interval remains a valid bound if that requested width is not reached. The composition retains the full maximum certificate, actual width, success flag and cell counts, and uses the actual interval in every residual calculation. It does not label an unmet reporting target as achieved or discard an otherwise rigorous conservative bound.

## All 68 constraints and the tolerance family

For observed value \(y_{ik}\) and a uniform prediction enclosure \([l_{ik},h_{ik}]\), exact arithmetic gives

\[
d^-_{ik}=\max(l_{ik}-y_{ik},y_{ik}-h_{ik},0),\qquad
d^+_{ik}=\max(|l_{ik}-y_{ik}|,|h_{ik}-y_{ik}|).
\]

Thus, when all 68 channels are covered,

\[
D^-:=\max_{i,k}d^-_{ik}
\le D(\theta):=\max_{i,k}|\mathcal O_{ik}(\theta)-y_{ik}|
\le D^+:=\max_{i,k}d^+_{ik},\qquad \theta\in C.
\]

For a **declared** common discrepancy allowance \(\delta\), the selected static root tuple is excluded throughout \(C\) if \(\delta<D^-\), and every selected tuple in \(C\) satisfies all constraints if \(\delta\ge D^+\). Between these cutoffs the enclosure is undecided. No measurement-noise level or biologically justified tolerance is fabricated.

If any channel is missing, its absence is retained. Covered channels can still establish a necessary lower bound for any hypothetical complete tuple extending those selected branches. There is then **no sufficient upper bound and no all-68 IN statement**.

## Quantifier and branch limits

The conclusion is: for every shared parameter value in \(C\), independently for each record there exists its selected stationary root with the stated static outputs. Separate dose-specific root certificates do not establish a coherent loading trajectory, dynamic or energetic stability, a common experimental time, or absence of other roots. Excluding this selected root tuple does not exclude every mechanical branch or the full physical parameter domain. The cell is a small computational support region, not a material prior or the full identified set.

The target is the exact finite-dimensional energy with the frozen binary integration cache. This composition certifies neither quadrature error relative to an ideal integration rule nor finite-element error relative to a continuum PDE, and it does not establish biological validity or novelty.

## Verification and retained artifacts

The [implementation](../scripts/cohort_static_observers.py) has separate synthetic-test, protocol-freeze, and composition stages. The [preflight](../runs/cohort_static_observers_v1/preflight.json) verifies all 17 saved graph inputs without root composition. The [synthetic receipt](../runs/cohort_static_observers_v1/synthetic_checks.json) checks exact affine hulls against all 32 box corners, residual formulas, incomplete-coverage handling, the full-radius rule, and rejection of nonuniform/mismatched inputs.

The [frozen protocol](../runs/cohort_static_observers_v1/protocol.json) binds recursive local runtime dependencies, the common contract, source CSV, all saved fields, the mechanical cohort protocol, and its independent prelaunch review. A separate immutable provider selection manifest can retain v1 root results together with an explicitly declared supplemental corrector-radius proof. It must cover every source record and lock each selected root result and protocol before observation composition. Supplemental certificates cannot silently extend a parameter cell, change a dose, replace a physical state, or erase an unresolved v1 attempt.

The [prelaunch independent review](../runs/cohort_static_observers_v1/prelaunch_independent_review.json) passed 24 checks before the single composition process. The [job receipt](../runs/cohort_static_observers_v1/composition_final/job.json) records true terminal completion and binds the exact source, protocol, review and provider selection. Per-record outputs retain the full root radius, predictor hulls, serialization gaps, whole-body result, transfer constants and endpoint intervals. The [visual receipt](../runs/cohort_static_observers_v1/composition_final/visual_receipt.json) and [inspection receipt](../runs/cohort_static_observers_v1/composition_final/visual_inspection.json) bind the displayed figure to the numerical result.

The subsequent [independent root review](../evidence/cohort_uniform_roots_independent_review.json) verified 15 record families and retained F/K as unresolved. The [application arithmetic replay](../runs/cohort_observer_review_v1/arithmetic.json) passed 278 checks; the final [independent partial application receipt](../evidence/cohort_static_observers_independent_review.json), SHA-256 `5e71814a8a266ff309f990dfcd07b7ba8078abee39c5db8ac10ea04aaa2af8ad`, passed 18 finalization checks and binds those 15 independently reviewed roots to 60 output intervals. Whole-body and base-observer checks are disclosed replays of previously reviewed components, not independently rewritten implementations. Status `PASS_PARTIAL_APPLICATION` certifies this limited composition; it does not declare the complete inference study passed. The [cohort checkpoint](COHORT-INFERENCE-CHECKPOINT.md) explains the remaining proof and scientific gaps.
