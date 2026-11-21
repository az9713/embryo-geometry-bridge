# Independent review of the F/K center-offset observers

**PASS for the eight conditional static outputs of F and K.** The review independently reconstructs the translated node boxes, predictor hulls, trace derivatives, graph-error constants and all eight residuals. It also replays the previously reviewed whole-body and baseline-observer components. Finalization binds the **separate fresh F/K root audit**, not the producer root calculations alone.

The [final receipt](../evidence/recentered_observers_independent_review.json) records 62 application/scope checks and 11 independent-root binding checks. The [arithmetic receipt](../runs/recentered_observer_review_v1/arithmetic.json), [frozen review protocol](../runs/recentered_observer_review_v1/protocol.json) and [review source](../scripts/review_recentered_observers.py) retain the calculation and its hashes. This review covers **8 outputs, not all 68**, and makes no global-fit claim.

## Why a new adapter was necessary

The original saved solver state supplies nominal P2 tissue graphs. The root proof now uses a nearby numerical center \(u_1\), obtained by one declared proof-center correction from \(u_0\). It is an auxiliary computational state, not new experimental data and not a claimed converged nonlinear-solver result.

The review checks every old/new cached array except `u` is byte-identical, **including the preconditioner \(R\)**; \(S\) is unchanged too. Independent fingerprints bind the reference geometry, quadrature, parameters, base coefficients, exact dose law and shared parameter cell. The old seed and new auxiliary-center identities stay distinct. Neither root proof nor observer silently replaces the seed data.

The fresh independent root audit is [this separate receipt](../evidence/cohort_recentered_roots_independent_review.json). Its two audit capsules bind the actual new root-certificate bytes, selected domains, shared cell and common model contract. Both reassemble all 611 free coordinates and 9,504 terms at the new centers with 192-bit interval arithmetic. The application review checks those bindings and their per-record files before issuing PASS; it does not recompute those root proofs itself.

## Exact node-box reconstruction

For free coordinate \(j\), a root in the reviewed family has

\[
 x_j=X_j+(u_1)_j+(S\delta\theta)_j+e_j,
 \qquad |e_j|\le r.
\]

Let \(\widetilde x_j\) be the **original serialized seed coordinate**. The correct node radius is

\[
 \nu_j=|X_j+(u_1)_j-\widetilde x_j|
       +\sup_{\theta\in C}|(S\delta\theta)_j|+r.
\]

The first term charges both center translation and serialization directly, with exact rational arithmetic. It is not the old seed serialization gap alone. At constrained coordinates there is only the exact fixed-coordinate gap; no corrector radius is assigned. The review reconstructs all 611 predictor hulls, all 650 coordinate radii and the exact translation vector.

The normal Banach radius remains \(10^{-11}\), centered at the **new** affine predictor. Translation is not added to the radius used by the root theorem. It enters the observer's original-seed-centered node box instead. Keeping those two roles separate is essential.

## Same-current-radius trace error

The review independently extracts top and bottom P2 edges from raw reference coordinates and six-node connectivity. For each trace, it checks the radial derivative, fixed radial endpoints, the same positive-width physical ROI and the actual surface-node error \(\varepsilon\), with physical scale \(L_0=33.5\) µm applied once.

P2 cardinal interpolation gives coordinate error at most \(5\varepsilon/4\) and radial-derivative error at most \(8\varepsilon\). Every edge satisfies the required strict radial margin. If the nominal piecewise graph has global Lipschitz upper bound \(L\), its same-current-radius vertical error is

\[
 E=\frac54\varepsilon(1+L).
\]

The global slope bound remains valid across edge junctions. The radial term matters because a point's material parameter generally moves to a different current radius.

With surface errors \(E_a,E_b\), the errors assigned to maximum height shortening and disk-mean height shortening are \(E_a+E_b\). Apical and basal component maxima receive \(E_a,E_b\), respectively. This does not identify a maximum of a sum with the sum of separate maxima. All four assignments and every interval endpoint are checked exactly.

## Disclosed component replays and residuals

The whole-body nodal-box certificate and baseline static observer are replays of already reviewed implementations. This review does **not** claim independent reimplementations of those algorithms. The body replay confirms the translated box is admissible, including the axis. The baseline replay evaluates the original saved seed graphs, with the fixed current-radius ROI, flat reference levels, signed spatial-maximum components and projected-disk mean. Actual maximum-enclosure widths are retained rather than assuming a requested reporting tolerance was achieved.

For observed value \(y\) and enclosed prediction \([a,b]\), the independently recomputed sharp residual range is

\[
 [\max(a-y,y-b,0),\ \max(|a-y|,|b-y|)].
\]

For orientation, the following are approximate interval midpoints and published source values, in µm. Exact prediction/residual intervals are in the arithmetic receipt.

| Record | Static quantity | Prediction midpoint | Source value |
|---|---|---:|---:|
| F | Maximum same-radius height shortening | 5.617206 | 7.648084 |
| F | Disk-mean height shortening | 4.487594 | 2.646384 |
| F | Apical downward maximum | 4.833357 | 6 |
| F | Basal upward maximum | 0.783848 | 3.648084 |
| K | Maximum same-radius height shortening | 10.490407 | 12.620972 |
| K | Disk-mean height shortening | 8.387087 | 8.383192 |
| K | Apical downward maximum | 9.137541 | 13 |
| K | Basal upward maximum | 1.352866 | 2.862915 |

Tiny computational intervals are not experimental precision or continuum error bars. Their narrowness does not validate the static proxy, source-time interpretation, canonical cross-sheet pairing or fixed geometry assumptions.

## Scope finding and preserved correction

The first producer aggregate correctly contained eight rows but inherited a generic field named `all68_uniform_static_tuple_enclosed=True`. That was an incorrect global-coverage label. The mathematical intervals were unaffected.

The original result, source, protocol and numerical record packets remain unchanged. The authoritative [scoped version-two wrapper](../runs/cohort_recentered_observers_v1/results_scoped_v2.json) and [scope-correction receipt](../runs/cohort_recentered_observers_v1/scope_correction.json) explicitly set the global all-68 flag false, declare F/K 8-of-8 coverage and forbid global IN. Reversing **only** those declared metadata changes recovers the entire original packet exactly; the review verifies this and every numerical row/hash link.

The eight F/K bounds can be combined with the separately audited 60 outputs by a distinct, hash-bound join. This receipt alone does not certify an all-68 fit, a full identified material set, other stationary branches, a stable/common experimental path, biological observation semantics, continuum accuracy or novelty. This is an independent project-agent check, not external peer review.
