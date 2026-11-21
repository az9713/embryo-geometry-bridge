# Independent review of the all-68 shared-cell stationary calculation

**The mathematical and provenance audit passes.** All 17 source inputs and four channels per input are accounted for, with the same parameter cell and observation contract and an independently reviewed root/observer proof pair for every record. The [receipt](../evidence/shared_cell_cohort_independent_review.json) binds the actual [combined result](../evidence/shared_cell_cohort_inference.json), its source and the publication-control snapshot reviewed alongside it. No force, Hessian, equilibrium, geometry or baseline-observer computation is repeated by this audit.

This completes a bounded conditional component: a uniform stationary-root tuple and 68 static output enclosures on one very small parameter cell. The original five-gate research objective is still incomplete.

## The object actually certified

The five coordinates are apical and basal area moduli, apical and basal prestresses, and a common activation gain. Every source dose uses exactly the same closed cell \(C\). Its four positive baselines vary by relative half-width \(10^{-8}\); apical prestress varies from zero to \(2\times10^{-10}\) in the declared normalized units. These are computational support bounds, not measured physiological uncertainties or a confidence region.

At each of the 17 exact source doses, a selected local stationary root \(u_i(\theta)\) exists for every \(\theta\in C\). The first 15 root families retain their original reviewed certificates. F and K use separately reviewed proof-center refinements for the unchanged physical target. They are not silently substituted into the old failed calculation.

The current output law has four quantities: maximum same-current-radius height shortening, disk-averaged height shortening, maximum apical downward displacement and maximum basal upward displacement. The fixed ROI and flat reference heights are part of the declared static contract. They are not assumed to reconstruct the experimental temporal observation pipeline.

Every output enclosure \([\ell_{ij},v_{ij}]\) is uniform over the cell and is derived from a full root/predictor/geometry error envelope. The audit compares these bounds with each observer's `prediction_intervals_exact`; it does not substitute the baseline point predictions. All 68 actual intervals have positive width.

## The exact discrepancy theorem

Define the maximum absolute static discrepancy

\[
D(\theta)=\max_{\substack{1\le i\le17\\1\le j\le4}}
\left|O_j(u_i(\theta))-y_{ij}\right|.
\]

For each observed value \(y_{ij}\), interval arithmetic gives

\[
d^-_{ij}=\max(\ell_{ij}-y_{ij},\;y_{ij}-v_{ij},\;0),
\qquad
d^+_{ij}=\max(|\ell_{ij}-y_{ij}|,\;|v_{ij}-y_{ij}|).
\]

Since the actual output lies inside its uniform interval, its absolute residual lies between these two bounds. Taking maxima preserves both inequalities. Therefore

\[
L:=\max_{ij}d^-_{ij}\le D(\theta)\le
U:=\max_{ij}d^+_{ij}
\qquad\text{for every }\theta\in C.
\]

The exact rational endpoints are retained in the combined result. Their decimal displays are

\[
L\approx5.186387141964336\;\mu\mathrm m,
\qquad
U\approx5.186387360305146\;\mu\mathrm m.
\]

Both bounding maxima are supplied by M's apical channel. The audit independently recomputes every one of the 68 residual pairs from the source CSV and the selected observer intervals, then recomputes the maxima and range width with exact fractions.

This is a pointwise conservative enclosure of \(D\) throughout the specified cell. It is not a proof that the interval endpoints are attained, a sharp minimax theorem, or a global material optimum. The small enclosure width does not include mesh/integration/continuum error, source processing uncertainty or physiological parameter uncertainty.

## What the two tolerance queries mean

For a declared common allowance \(\delta\), the theorem implies:

| Query | Certified implication on this cell |
|---|---|
| \(\delta<L\) | Every parameter in the cell fails for the supplied selected stationary tuple. Other branches remain unresolved. |
| \(\delta\ge U\) | Every parameter in the cell has the supplied stationary tuple satisfying all 68 declared static comparisons. |
| \(L\le\delta<U\) | The interval calculation alone does not decide compatibility. |

The recorded 5.1 micrometre query is below \(L\), so the selected-tuple decision is OUT. That does not prove that all stationary branches or the whole physical parameter domain are excluded.

The recorded 5.2 micrometre query is above \(U\), so every point of this nonempty cell has a stationary-tuple witness for that allowance. It supplies an inner component of the **declared conditional** compatibility problem. It does not assert that every other branch passes or that the full identified set has been computed.

These illustrative queries were selected after the earlier M threshold was known. They are neither prospective statistical tests nor estimates of experimental noise. Increasing an allowance changes the question being asked; it does not demonstrate that the model matches the source's actual measurement uncertainty.

## The 60-plus-8 evidence chain

The combination selects 15 earlier audited records, giving 60 outputs, and exactly F and K's two newly audited records, giving eight outputs. It rejects duplicate selected records and requires four channels per selected record. The independent audit checks:

- exact source-dose association and the same nonempty cell and contract across all records;
- each root certificate's bytes, passing uniform-root capsule and local branch identifier;
- each observer certificate's bytes, passing application checks and its matching root proof;
- the new F/K application's explicit uniform-quantifier and admissible-family flags;
- all frozen combination inputs and both exact tolerance-query decisions.

The raw two-record F/K output originally inherited a misleading nested all-68 flag. Its separately retained metadata correction makes the nested scope explicitly eight outputs and leaves all record files, numeric intervals and aggregate numeric bounds unchanged. The combined 68-output conclusion is made only in the new outer calculation. The audit verifies the corrected wrapper, correction receipt and independent F/K application review together; it does not treat the original nested flag as a valid global result.

The underlying body and baseline-observer implementations have their disclosed reviewed-component replays. This audit does not relabel those as another independent implementation. Neither does a stationary tuple constitute one stable physical loading path through all activation doses.

## What remains outside this result

The original contract still requires adequate source-compatible physical inference and uncertainty-propagated larger-region prediction. One tiny fixed-geometry cell is not the full material/geometry domain. Other branches, zero-coefficient faces, unbounded parameter tails and independent physiological constraints are not resolved by this calculation.

The experimental temporal peak/reference/ROI/aggregation procedure, specimen-level source pairing and measurement/discrepancy law remain conditional. The finite model uses frozen binary quadrature and geometry constants; its exact discrete-root theorem is not an ideal-integral or continuum-PDE error theorem. No revised larger-region biological validation or confirmed scientific novelty follows from the compatibility witness alone.

The meaningful next operations concern those remaining inference, forward-error, source and transfer requirements. Another local PASS or another editorial audit cannot substitute for them.

## Publication-control snapshot

The original PLAN remains byte-identical. The current gate ledger still marks Gate 3 in progress, Gate 4's original target incomplete and the full objective unachieved. Its latest shared-cell entry binds the exact 17-record/68-channel result; the earlier partial phase is explicitly historical and superseded. All eight prior control snapshots in the parameter-family archive match their retained hashes.

The README still states that the original 15-step project is incomplete and claims no new GitHub Pages deployment. Both landing-page generators lead with the completed conditional calculation and retain its scientific limits. Their shortened public enclosure, [5.18638714, 5.18638737] micrometres, is rounded safely outward from the exact endpoints.

This audit identified one older parent-landing sentence that still called the all-68 calculation incomplete. It was corrected to distinguish that historical checkpoint from the current conditional result. The child single-root heading and limitation were likewise clarified as historical and specific to that earlier one-field proof.

This receipt binds publication **sources**. Final generated HTML, equation rendering, link layout and any live external service are outside this particular audit and need their separately retained build/browser receipts. The final HTML can include this review without a circular source-versus-generated-page hash claim.
