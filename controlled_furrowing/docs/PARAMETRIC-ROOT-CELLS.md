# Uniform equilibrium enclosures over parameter cells

**The new calculation covers whole five-parameter boxes, rather than isolated coefficient choices.** A local branch of discrete roots is certified for every parameter in selected small cells. One cell is attached to actual calibration record M and produces a uniform static-apical prediction. The certified five-dimensional widths remain extremely small; these cells are neither measured physiological priors nor the full compatible parameter set.

The [implementation](../scripts/parametric_root_cells.py), [pilot protocol/results](../runs/parametric_root_cells_v1/protocol.json), [actual-M protocol](../runs/parametric_root_cells_v1/actual_M/protocol.json), and [actual-M results](../runs/parametric_root_cells_v1/actual_M/results.json) retain the target, exact endpoints, predictor, every failed sufficient test and an explicit cover ledger. The [actual-M independent review](PARAMETRIC-ACTUAL-M-REVIEW.md) distinguishes the uniform-root proof from its composition with geometry and observations.

## The family being certified

The five coordinates are

\[
\theta=(K_a,K_b,T_a,T_b,g).
\]

Here apical means the outer cellular surface and basal the inner boundary. The two \(K\)'s are normalized surface area-elastic coefficients, the two \(T\)'s are surface tensions, and \(g\) converts the reported normalized fluorescence input into a model loading coefficient. That proportionality is a modeling assumption, not an independent force measurement.

The shear normalization is one, compression coefficient is fixed at 5000, and geometry, finite-element space and quadrature are fixed. A trial cell is a closed rectangle \(C\) in these coordinates. For a fixed input dose \(I\), activation enters as the **exact product** \(gI\). A stored numerical center need not have been solved at exactly that product: its residual is freshly enclosed for this declared target.

The pilot uses the previously reviewed weak field at representative dose one. It covers no actual source-record input. The second calculation uses the saved new-primary-material field for `optoGEF2_small_col_M`, whose source dose is exactly 0.78208913 in the canonical extract. Its solver strength was approximately 6.1101921642. The requested product and its binary solver rounding remain distinct; the uniform proof does not inherit a root claim at the wrong strength.

All five-dimensional cells have positive width. The apical tension begins at zero and varies on its nonnegative side. The other coordinates use relative half-widths around the stated base values. This is an explicitly chosen numerical family, with no inference that zero faces, other coefficient ratios, geometries or unbounded tails have become inadmissible.

## A parameter-dependent predictor

Let \(u_0\) be the archived 611-variable coefficient vector, and choose any fixed binary matrix \(S\in\mathbb R^{611\times5}\). Define

\[
\widehat u(\theta)=u_0+S(\theta-\theta_0).
\]

The predictor is constructed from stored implicit tangent columns, converted to absolute coefficient coordinates. It is an approximation, not an exact derivative assumption. The proof explicitly retains its first-order defect. The constant-center comparison uses \(S=0\).

At fixed geometry and normalization, the residual is affine in the five material/input coefficients. Write

\[
f_0=f(u_0,\theta_0),\quad H_0=D_uf(u_0,\theta_0),\quad
B_0=D_\theta f(u_0,\theta_0),\quad d=\theta-\theta_0.
\]

Exact integration along the displacement segment gives

\[
f(u_0+Sd,\theta_0+d)=f_0+(H_0S+B_0)d
+\int_0^1[H(u_0+tSd,\theta_0+d)-H_0]Sd\,dt.
\]

This identity includes mixed parameter/displacement effects. Dropping the integral or replacing \(H_0S+B_0\) by zero because a floating tangent equation was approximately solved would invalidate the certificate.

For each quadrature point, local strain is affine in displacement. The predictor variation is formed as \((B_qS)d\) **before** taking interval bounds. This preserves cancellations that would be lost by first replacing every nodal variation with one global absolute radius.

The existing independently reviewed bulk/surface derivative formulas bound local Hessian changes throughout both the predictor segment and the error tube. Elementwise row bounds are assembled and preconditioned with the actual stored matrix \(R\). Every parameter endpoint, normalized weight, denominator and positive-domain check is carried through Arb enclosures.

## The uniform theorem

Suppose the local analytic domains hold for all \(\theta\in C\) and \(\|e\|_\infty\le r\). Assume validated bounds satisfy

\[
\eta(C)\ge\sup_{\theta\in C}\|Rf(\widehat u(\theta),\theta)\|_\infty,
\]

\[
Z(C,r)\ge\sup_{\theta\in C,\|e\|_\infty\le r}
\|I-RH(\widehat u(\theta)+e,\theta)\|_\infty,
\]

\[
Z(C,r)<1,\qquad \eta(C)+Z(C,r)r<r.
\]

For each parameter, the map

\[
T_\theta(e)=e-Rf(\widehat u(\theta)+e,\theta)
\]

is a contraction from the same closed error cube into itself. The independently bounded center defect proves \(R\) nonsingular. Banach's theorem therefore gives, for **every** \(\theta\in C\), a unique root in \(\widehat u(\theta)+B_r\). Its corrector is bounded by \(\eta/(1-Z)\). Continuous parameter dependence follows from uniform contraction and continuity of the maps.

This is a nearby stationary branch. It neither excludes roots outside the declared domain nor proves energetic stability or a common quasistatic loading history across experimental doses.

## Declared attempts and their limits

Before evaluating cells, the protocols froze six relative half-widths from \(10^{-12}\) to \(10^{-2}\), constant and affine predictors, a radius-selection rule, and one bounded bisection demonstration. Each target produced fourteen attempts: twelve initial method/cell combinations and two retained children of the unresolved cover root.

| Target | Largest tested certified constant-center width | Largest tested certified affine-predictor width | Certified / unresolved attempts |
|---|---:|---:|---:|
| Old weak representative field | \(10^{-10}\) | \(10^{-8}\) | 5 / 9 |
| Actual source-record M, new material | \(10^{-10}\) | \(10^{-8}\) | 5 / 9 |

![Declared local cell certification statuses](../runs/parametric_root_cells_v1/cell_status.svg)

These are largest **tested successful** widths, not optimal cutoffs. The cells are nested and methods overlap; the counts are not disjoint coverage volumes or probabilities. At larger widths the sufficient tests remained unresolved. Some proposed radii also exceed the declared test budget or analytic-domain enclosure; neither condition means the parameter cell has no equilibrium.

For actual M's selected affine cell, the author enclosures give approximately

\[
\eta=1.74424215\times10^{-12},\qquad
Z=0.096862311,\qquad r=10^{-11}.
\]

The producer's posterior corrector bound is about \(1.931314\times10^{-12}\). The forward observer deliberately uses the **full common certified cube** \(r=10^{-11}\), together with exact rational predictor extrema. This avoids requiring upper majorants obtained at different Arb precisions to be numerically identical.

The selected parameter ranges, rounded only for display, are:

| Coordinate | Lower | Upper |
|---|---:|---:|
| \(K_a\) | 0.2210341814 | 0.2210341859 |
| \(K_b\) | 1.8096748179 | 1.8096748542 |
| \(T_a\) | 0 | \(2\times10^{-10}\) |
| \(T_b\) | 2.2103418140 | 2.2103418583 |
| \(g\) | 7.8126544261 | 7.8126545824 |

All exact endpoints are retained. The small width is a limitation of this first infinity-norm enclosure, not evidence of biological parameter precision. A separate energy-metric investigation addresses wider gain-only intervals; it must not be presented as an already proved five-dimensional or actual-M result.

## An actual-M uniform forward statement

The [model/observation contract](../runs/parametric_root_cells_v1/actual_M/model_observation_contract.json) fixes the primary static operator: projected-disk weighting for means, spatial-max signed components, flat reference profiles, the same current-radius ROI on both surfaces, and canonical source-column pairing. A line-weighted or temporal local-peak observer is a different contract.

The [root certificate](../runs/parametric_root_cells_v1/actual_M/root_certificate.json) is composed with a whole-body nodal-box geometry proof and the reviewed nodal-to-current-radius observer theorem. Exact binary addition gaps are included. The geometry proof covers every field in the outer box, preserving the axis and strong monotonicity; it is not just a check at volume quadrature points.

The [uniform observer certificate](../runs/parametric_root_cells_v1/actual_M/observer_certificate.json) supplies only the apical endpoint. Its interval is approximately

\[
A_M(\theta)\in[6.813612639694854,\ 6.813612858035664]\ \mu\mathrm m
\quad\text{for all }\theta\in C.
\]

The canonical source value is 12 µm. Consequently this one constraint has a residual bracket approximately

\[
|A_M(\theta)-12|\in
[5.186387141964336,\ 5.186387360305146]\ \mu\mathrm m.
\]

The exact [tolerance-family diagnostic](../runs/parametric_root_cells_v1/actual_M/contradiction_diagnostic.json) therefore excludes this selected-branch cell for every apical allowance below the lower bound. An allowance at least the upper bound satisfies this **one** constraint throughout the cell; it does not establish compatibility with the other 67 constraints. The illustrative allowance 5.1 µm is declared, not estimated measurement noise. The interval endpoints are conservative bounds, not proved attained thresholds.

![Conditional apical tolerance-family cut for actual M](../runs/parametric_root_cells_v1/delta_profile.svg)

This is conditional on the recorded static operator and nearby branch. Other branches, source-operator alternatives and the rest of the parameter domain are not excluded. Independent review must bind the root, observation, data, cell and model-contract hashes before a cover engine consumes the provider.

That review now passes: the [frozen independent receipt](../evidence/parametric_actual_M_independent_review.json) retains 19 root-replay checks, 22 composition checks and six packaging checks. It independently replays the selected actual-M cell, checks the exact gain–dose law, reconstructs rational predictor and nodal-gap bounds, replays the body box, and reconstructs the apical extrema against the source workbook. It binds the current full-corrector observer, not the preserved earlier composition. This is independent verification of the selected actual-M provider; it is not a claim that every one of the 28 retained assemblies was independently rerun.

The [author's retained-record verifier](../runs/parametric_root_cells_v1/verification.json) separately checks all 28 attempts, 64 radius tests, exact self-map/error inequalities, artifact bindings and the bisection ledger. Its 26 checks are arithmetic/provenance checks of the retained records, with the narrower coverage stated explicitly.

## Adaptive coverage without silently dropping uncertainty

The included primitive bisects the widest normalized coordinate interval exactly. Its two closed children have union equal to the parent, with only a shared boundary. After the declared budget, every terminal leaf is retained as certified or unresolved. A failed contraction never becomes `OUT`.

That ledger covers only its stated root box. It does not cover the whole ratio domain. Even the chart's positive-coordinate assumption is local: globally allowed zero faces and unbounded material ratios remain outside this cover. Geometry, other source doses and physical observation uncertainty require their own enclosures and quantifiers.

## Reproducibility and what advances

The historical pilot source/protocol/results remain unchanged and hash-bound in `runs/parametric_root_cells_v1/pilot_source_snapshot/`. Generalization to actual M used a separately frozen subrun and freshly evaluated Hessian/preconditioner at the saved state, with **no new equilibrium optimization**. The first observer composition and source remain in `actual_M/observer_revision_0/`; version 2 uses exact predictor hulls and the conservative full corrector radius. No failed cell or earlier artifact was erased from the record.

The advance is the quantified forward implication over a continuum of coefficients, and its first connection to an actual calibration constraint. The bounds remain local and narrow. This is standard validated-numerics machinery, not a new concept, a physical identified set or a validated chemistry-to-shape explanation. The full five-gate requirement still includes the remaining doses, domain tails, geometry, quadrature/continuum error, source temporal semantics and revised frozen larger-region predictions.
