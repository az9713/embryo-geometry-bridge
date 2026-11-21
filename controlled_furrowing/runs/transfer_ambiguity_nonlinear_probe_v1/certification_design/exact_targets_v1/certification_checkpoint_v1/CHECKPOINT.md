# Four cached roots and the exact-spatial test

The two frozen material candidates retain a strictly positive large-patch contrast difference in the independently verified **cached-quadrature model**. The signed gap is **0.15475515644140 to 0.15475515688952 µm**. The exact-spatial finite-element root test remains pending in this checkpoint.

## What the four proofs cover

Each candidate contributes one small state, record L, and one large state at the common matched input. The proofs certify four local roots of the **611-coordinate binary Q8/S12 quadrature potential**, with the same frozen discrete model and each state's exact requested optical drive. Q8/S12 refers to the quadrature orders for bulk and surface terms. The drive in the exact arithmetic is the candidate gain times its exact dose. The rounded drive that the numerical solver used is retained separately.

The four independent root proofs use closed cubes of radius \(2^{-36}\) in dimensionless free coordinates. A separate checked geometry certificate covers the entire larger cube of radius \(r=2^{-34}\) for each center. A second independent interval replay encloses the cached derivative over that whole larger cube with the same binary preconditioner \(R\). Define \(\overline{Z}_Q\) as the resulting upper bound, with \(H_Q\) the cached-quadrature Hessian.

\[
\begin{aligned}
\|I-RH_Q(u)\|_\infty &\le \overline Z_Q,\\
\|u-u_0\|_\infty &\le 2^{-34}.
\end{aligned}
\]

| Frozen candidate and outcome | \(\overline Z_Q\) | Room below one |
|:--|--:|--:|
| Candidate + · small L | 0.320883 | 0.679117 |
| Candidate + · large | 0.340514 | 0.659486 |
| Candidate − · small L | 0.306458 | 0.693542 |
| Candidate − · large | 0.318140 | 0.681860 |

The table and figure display rounded dimensionless bounds. The receipts retain exact rationals. The room below one leaves headroom for the spatial Hessian defect. The signed point-force condition is a separate test.

![Four reviewed cached derivative bounds at radius 2 to the minus 34. The bars show the room left below the contraction threshold of one.](cached_derivative_headroom.png)

[Four cached-root proofs](../../cached_pilots_v1/independent_actual_preparation/actual_replay/receipt.json) · [Whole-cube geometry replay](../independent_actual/replay/receipt.json) · [Larger cached derivative replay](../cached_derivative_extension_v1/independent_actual_preparation/actual_replay_v2/receipt.json) · [Figure provenance](figure_receipt.json)

## What the certified cached roots say about the outputs

The large-patch observable \(Q\) is the current apical height on the axis minus its height at current radius **24.94698765 µm**. The reported interval encloses \(Q_{+}-Q_{-}\), with the candidate order frozen before the nonlinear comparison. The independent transfer uses each root displacement bound plus the exact binary addition gap between the analytic center and serialized coordinates. A monotone quadratic trace on the common current-radius interval then bounds the apical shoulder value. This gap measures the response of two fixed model candidates. It gives no estimate of measurement error or experimental resolution.

For the **two corresponding L states**, the cached-root observer application independently checked **32** inequalities across four static conventions. **28 were IN and 4 were OUT** at the original allowance \(\delta \simeq 5.3898356748\,\mu\mathrm{m}\). Every OUT comparison was the regional mean of total shortening under radial-line weighting. IN and OUT refer to whether a rigorous absolute residual interval lies within or above the frozen allowance. The allowance is a model-discrepancy threshold chosen for this diagnostic.

The full saved-field comparison is sensitive to observation convention. These worst-case common required allowances come from the separate, independently reviewed sweep of the **34 saved small fields**, using the same two material candidates and recorded values.

| Static convention | Common required allowance (µm) | At the frozen allowance |
|:--|--:|:--|
| Disk · separate component maxima | 5.385023 | Within |
| Disk · component means | 6.576090 | Exceeds |
| Radial line · separate component maxima | 6.191999 | Exceeds |
| Radial line · component means | 6.191999 | Exceeds |

Only the declared primary convention, disk weighting with separate apical and basal spatial maxima, fits the same allowance for both saved-field candidates. The stored-field sweep leaves the source timing, initial profiles, specimen-level pairing and ROI registration unresolved.

[Exact cached observer readback](../../cached_pilots_v1/observer_application_v1/independent_readback/DERIVATION.md) · [Exact observer intervals](../../cached_pilots_v1/observer_application_v1/independent_readback/receipt.json) · [Independent saved-field observer sweep](../../../observer_robustness/independent_actual/receipt.json) · [Full nonlinear report](../../../report/REPORT.html)

## What exact spatial integration must still prove

The next test keeps the same closed material singleton and binary \(R\). Let \(F_I\) and \(F_Q\) denote the exact-spatial and cached-quadrature reduced forces. At the exact center \(u_0\), the spatial force defect \(D_0\) and cached anchor must be combined with their signs intact.

\[
\begin{aligned}
D_0 &= F_I(u_0)-F_Q(u_0),\\
R F_I(u_0) &= R F_Q(u_0)+R D_0 .
\end{aligned}
\]

Over the **entire** matching cube, a bound \(b\) on \(\|R(H_I-H_Q)\|_\infty\) gives a spatial derivative bound \(z=\overline Z_Q+b\). Let \(\eta\) bound the infinity norm of the signed full center force \(R F_I(u_0)\). A sufficient strict local root test is

\[
z<1,\qquad \eta+zr<r .
\]

The plan calls for **144 bulk and 24 surface components** per outcome, with independent replay and signed global reduction. The retained source and synthetic checks verify the arithmetic interface. This checkpoint has no admitted exact-spatial root result for these witnesses. A spatial-root certificate would still leave the energy-minimum proof and continuum error estimate open.

[Independent component-replay method](../bridge_v1/independent_actual_preparation/ROOT-SOURCE-READBACK.md) · [Signed reduction derivation](../bridge_reduction_v1/PROOF.md) · [Pure reduction source review](../bridge_reduction_v1/independent_review/final_binding.json)

## Coverage of the original five gates

The frozen comparison has **36 outcomes**, comprising **34 small final fields and two large final fields**. The other eight saved checkpoints document numerical loading. At the four-pilot snapshot, **32 outcomes** lay outside the cached proof set. Even if all four exact-spatial pilots succeed, they cover just **8 of 136** primary small-region inequalities and the two large operands. **128 primary inequalities** still need exact-root observation transfer.

**Cached coverage update.** A separate follow-up now certifies the two M states in the cached-quadrature model, giving **6 of 36** cached roots and leaving **30** at that stage. Candidate + at M passed on radius \(2^{-36}\). Candidate − at M was unresolved on radius \(2^{-32}\). Its independent replay passed on the narrower \(2^{-34}\) cube at the same center, drive and binary preconditioner. [Positive M proof](../../full_outcome_coverage_v1/remaining_cached_wave1/independent_actual_preparation/actual_replay/max_contrast_change_small_09.json) · [Negative M narrow-cube proof](../../full_outcome_coverage_v1/remaining_cached_wave1/negative_M_narrow_v1/actual/receipt.json)

1. **Originality.** Existing work already establishes related small-to-large transfer and mechanical discrimination strategies. The bounded audit records the prior art. A new general concept has not been established.
2. **Data feasibility.** Public summaries and observation conventions are documented. A matched geometry, reference-profile, timing, transported ROI and measurement-discrepancy contract still needs evidence.
3. **Model and identifiability.** The two fixed models supply a conditional comparison. Full physical parameter and geometry coverage, branch admissibility, and continuum error control remain open.
4. **Frozen prediction.** The pair and its static diagnostic were frozen. Revised uncertainty-propagated predictions and matched larger-region outcome scoring remain open.
5. **Contribution.** The positive cached contrast gap is quantitative evidence under a declared model and observer. It does not establish biological transfer, experimental detectability or novelty. The study's full contribution assessment remains open.

[36-outcome coverage ledger](../../full_outcome_coverage_v1/COVERAGE.md) · [Full checkpoint readback](receipt.json) · [Source observation audit](../../../../../docs/JOINT-OBSERVATION-CONTRACT.html) · [Originality follow-up](../../../originality_followup/README.md) · [Five-gate contract](../../../../../PLAN.md)
