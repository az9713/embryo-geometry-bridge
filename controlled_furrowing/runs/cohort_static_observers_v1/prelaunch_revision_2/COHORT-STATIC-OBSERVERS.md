# Uniform static observation of the 17-record stationary-root family

This composition implements the four primary observations for every small-region source record over **one shared five-parameter cell**. It performs no equilibrium optimization. The mechanical root proofs and the observation composition are separate steps, and a producer composition is not an independently audited inference result.

The common contract is [model_observation_contract.json](../runs/cohort_static_observers_v1/model_observation_contract.json), SHA-256 `2498493bf50049295eb2ad0d2360b1defaf387ee8db8828ceb576f3d7e942f82`. It contains the fixed geometry, material-coordinate conventions, integration policy, exact shared cell, physical length, primary observer, source pairing, and limits. Per-record dose, saved state, predictor, root domain, and provenance are deliberately separate.

## Four observations at a common current radius

Let (A(R)) and (B(R)) be the current apical and basal graph heights. The flat reference is (A_0=33.5\,\mu\mathrm m), (B_0=0). On the same fixed registered current-radius disk (0\le R\le R_*), define

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

For record (i), the independent mechanical proof must establish, for every parameter \(\theta\in C\), a unique stationary root in

\[
u_i(\theta)\in u_i^0+S_i(\theta-\theta_0)+[-r_i,r_i]^{n_i}.
\]

The proof uses the **full accepted corrector cube** (r_i), not the smaller a posteriori distance estimate. A strict Banach self-map/contraction bound is checked again with exact rational endpoints:

\[
0\le Z_i<1,\qquad \eta_i+Z_i r_i<r_i.
\]

Every binary64 predictor coefficient is treated as its exact binary rational. For each free coordinate (k), interval linear arithmetic gives the exact range

\[
[p^-_{ik},p^+_{ik}]
=\sum_{j=1}^{5}S_{ikj}[C_j^--\theta_{0j},C_j^+-\theta_{0j}].
\]

The serialized saved node is not silently identified with the exact reference-plus-displacement centre. If (x_{ik}^{\rm save}) denotes that saved coordinate and (X_{ik}+u^0_{ik}) its exact arithmetic centre, the nodal coordinate error is bounded by

\[
\epsilon_{ik}^{\rm norm}=
|X_{ik}+u^0_{ik}-x_{ik}^{\rm save}|
+\max(|p^-_{ik}|,|p^+_{ik}|)+r_i
\]

for a free coordinate. For a fixed coordinate only the serialization gap remains. This records every term instead of estimating error from force residuals, mesh differences, or a posteriori sampling.

The complete nodal box is passed through the exact whole-body geometry verifier: correlated Bernstein determinant/hoop bounds, axis conditions, and strong monotonicity imply an admissible embedded body for the entire box. A failed sufficient condition leaves that record unresolved; it is not evidence that the physical body fails.

## Boundary-specific transfer to the Eulerian graphs

Let \(\epsilon_A\) and \(\epsilon_B\) be (33.5\,\mu\mathrm m) times the largest coordinate radius on each corresponding surface's boundary nodes. The quadratic interpolation value bound is (5\epsilon/4), and the derivative perturbation bound is (8\epsilon) on the unit edge parameter. Strict radial monotonicity is therefore verified separately for each surface.

If (L_A,L_B) bound the saved graph slopes, then

\[
\|A_\theta-A_{\rm save}\|_\infty\le
e_A=\frac54\epsilon_A(1+L_A),\qquad
\|B_\theta-B_{\rm save}\|_\infty\le
e_B=\frac54\epsilon_B(1+L_B).
\]

This comparison is on a common current-radius domain. The code verifies that the axis and outer radial endpoint coordinates are fixed DOFs with **exactly zero** predictor/corrector/serialization radii. It then checks that both graph domains cover the prescribed positive-width ROI. It does not infer fixed endpoints from a flag passed to the transfer helper.

The map (f\mapsto\sup f) is 1-Lipschitz in the uniform norm. A normalized nonnegative weighted mean also has operator norm one. Consequently the four endpoint error bounds are

\[
(e_A+e_B,\quad e_A+e_B,\quad e_A,\quad e_B).
\]

The base-state endpoint intervals come from exact rational quadratic graph inversion and integration. These intervals are expanded by the appropriate channel-specific transfer error. Both the finite base-observer enclosure and the root-tube transfer error are retained.

## All 68 constraints and the tolerance family

For observed value (y_{ik}) and a uniform prediction enclosure ([l_{ik},h_{ik}]), exact arithmetic gives

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

For a **declared** common discrepancy allowance \(\delta\), the selected static root tuple is excluded throughout (C) if \(\delta<D^-\), and every selected tuple in (C) satisfies all constraints if \(\delta\ge D^+\). Between these cutoffs the enclosure is undecided. No measurement-noise level or biologically justified tolerance is fabricated.

If any channel is missing, its absence is retained. Covered channels can still establish a necessary lower bound for any hypothetical complete tuple extending those selected branches. There is then **no sufficient upper bound and no all-68 IN statement**.

## Quantifier and branch limits

The conclusion is: for every shared parameter value in (C), independently for each record there exists its selected stationary root with the stated static outputs. Separate dose-specific root certificates do not establish a coherent loading trajectory, dynamic or energetic stability, a common experimental time, or absence of other roots. Excluding this selected root tuple does not exclude every mechanical branch or the full physical parameter domain. The cell is a small computational support region, not a material prior or the full identified set.

The target is the exact finite-dimensional energy with the frozen binary integration cache. This composition certifies neither quadrature error relative to an ideal integration rule nor finite-element error relative to a continuum PDE, and it does not establish biological validity or novelty.

## Verification and retained artifacts

The [implementation](../scripts/cohort_static_observers.py) has separate synthetic-test, protocol-freeze, and composition stages. The [preflight](../runs/cohort_static_observers_v1/preflight.json) verifies all 17 saved graph inputs without root composition. The [synthetic receipt](../runs/cohort_static_observers_v1/synthetic_checks.json) checks exact affine hulls against all 32 box corners, residual formulas, incomplete-coverage handling, the full-radius rule, and rejection of nonuniform/mismatched inputs.

The [frozen protocol](../runs/cohort_static_observers_v1/protocol.json) binds recursive local runtime dependencies, the common contract, source CSV, all saved fields, the mechanical cohort protocol, and its independent prelaunch review. A separate immutable provider selection manifest can retain v1 root results together with an explicitly declared supplemental corrector-radius proof. It must cover every source record and lock each selected root result and protocol before observation composition. Supplemental certificates cannot silently extend a parameter cell, change a dose, replace a physical state, or erase an unresolved v1 attempt.

At this document's initial checkpoint, source/protocol and preflight are frozen; observation composition awaits concrete prelaunch review and a terminal, hash-bound provider selection. Final outputs remain **producer** results until independent review of both the root family and observation composition.
