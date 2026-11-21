# Joint calibration requires all four observations and a matching operator

**A valid height-only fit need not fit the other three reported quantities.** The new exact constraint evaluator processes all 68 entries in the 17-record small-region table, requires four explicit discrepancy allowances, and distinguishes certified acceptance, certified rejection and unresolved forward intervals. It is an inference prerequisite; it does not by itself compute the physical compatible set.

The [implementation](../scripts/joint_calibration_contract.py), [30-check verification](../evidence/joint_calibration_contract_verification.json), and [actual all-four zero-activation baseline](../evidence/joint_zero_activation_baseline.json) are retained. No larger-region outcome enters this work. The evaluator does not validate the physical model, observation semantics, interval provenance or continuum accuracy supplied by its caller.

## Exact interval membership predicate

For small-region record \(i\), observable \(k\), let \(y_{ik}\) be the recorded value, \(\delta_k\ge0\) an explicitly declared common discrepancy allowance, and suppose the chosen forward model's true output is enclosed by

\[
q_{ik}(\theta)\in[\ell_{ik},u_{ik}].
\]

The forward bound must refer to the same quantity as the data. An exact interval for a stored discrete static curve is not automatically an interval for a continuum or time-dependent experimental quantity.

The least and greatest possible residual on this interval are

\[
d_{ik}^{\min}=\max\{\ell_{ik}-y_{ik},\ y_{ik}-u_{ik},\ 0\},
\qquad
d_{ik}^{\max}=\max\{|\ell_{ik}-y_{ik}|,|u_{ik}-y_{ik}|\}.
\]

The first formula is the distance to a closed interval. The second follows because the convex absolute-value function attains its maximum at an endpoint. Therefore:

* If some \(d_{ik}^{\min}>\delta_k\), the parameter candidate is **OUT**: that prediction interval is disjoint from its permitted observation band.
* If every \(d_{ik}^{\max}\le\delta_k\), the candidate is **IN**, conditional on valid model-output enclosures and matching observation definitions.
* Otherwise the supplied intervals leave membership **UNRESOLVED**.

Overlapping an allowed band is not proof that a physically realizable response lies in the overlap. Different output coordinates can be coupled by the same mechanics. In particular, choosing the best point independently in each interval does not construct one common compatible parameter.

The per-observable requirements

\[
L_k=\max_i d_{ik}^{\min},\qquad U_k=\max_i d_{ik}^{\max}
\]

give necessary and sufficient thresholds for rejection/robust interval acceptance in the preceding sense. A discrepancy below \(L_k\) excludes the candidate; discrepancy at least \(U_k\) for every observable accepts every response enclosed by the supplied intervals. Values between them require tighter forward bounds or another argument. All calculations in this module use exact rationals for the supplied decimals.

These are standard set-membership inequalities, not a new inference theorem. No missing observation is dropped or replaced by zero: the exact record IDs and all four keys are mandatory.

## An optional unknown specimen-pairing contract

The [source audit](JOINT-OBSERVATION-CONTRACT.md) verifies the dose/maximum/mean cross-panel association, and that the apical and basal values within each Fig2e column form a pair. Linking those component pairs to the other panels uses the author-provided common B:R column order, without a separate specimen identifier. Canonical ordering is a reasonable explicit assumption. It is not upgraded to an independently verified specimen crosswalk.

The evaluator also implements a conservative alternative: retain fixed dose/maximum/mean associations but allow **one unknown permutation of the 17 intact apical/basal pairs**. The two components are never permuted independently. This is a different, weaker inference contract; it is not an assertion that the original columns were unordered or that biological samples were dependent.

For a fixed mechanical candidate, make a bipartite graph from predicted dose-records to observed component-pairs. An outer edge exists when both component intervals intersect their allowed observation bands; an inner edge exists when both intervals are wholly inside those bands. Every inner edge is an outer edge. Then:

* If a fixed maximum/mean constraint is OUT, or the outer graph has no perfect matching, the candidate is OUT even allowing an unknown pairing.
* If every fixed maximum/mean constraint is IN and the inner graph has a perfect matching, that one pairing certifies IN under this explicitly weakened contract.
* Otherwise the available bounds leave the decision UNRESOLVED.

In the second case, the same permutation works for every component value enclosed by the supplied intervals. The first case does not rely on choosing independent values in each interval: absence of even a possible pairing excludes every realizable one. The intermediate case makes no existence claim.

Failure of matching is retained as a directly checkable Hall obstruction: a set \(S\) of predicted records with

\[
|N(S)|<|S|.
\]

No injective assignment can send all of \(S\) into its smaller neighbor set. The implementation obtains the obstruction from alternating paths starting at unmatched left vertices. It also retains complete inner/outer edge lists and any matching, so the certificate can be replayed without trusting the search algorithm.

All 512 possible 3×3 bipartite graphs were checked against independent brute-force assignments; every deficient case's Hall witness was checked directly. A swapped-pair example is correctly rejected by canonical pairing and accepted by the optional permutation contract. Integer vertex/graph-size and closed-interval input guards were strengthened after independent review found permissive type handling. Matching and Hall's theorem are established mathematics, not a project novelty.

## A real but deliberately uninformative mechanical baseline

With common activation gain \(g=0\), the flat reference is an equilibrium of the stated fixed-rim model for every admitted positive shear/compression coefficient and nonnegative membrane coefficients. It is also a global energy minimizer under the usual admissible embedded-disk surfaces and fixed circular rim. To see this, each bulk singular-value contribution is nonnegative:

\[
s^2-1-2\log s\ge0\quad(s>0),
\]

and the logarithmic-volume and area-spring terms are nonnegative. Each deformed membrane has area at least the projected disk area \(\pi R_o^2\), by the fixed boundary and the area/projection inequality. Nonnegative surface tension therefore cannot reduce its energy below the flat reference value. The reference attains all these lower bounds. This is an elementary variational baseline, not a physiological claim.

Under the declared zero component reference, all four response observables vanish, independent of whether the spatial averaging convention is line or disk area. Direct evaluation against the actual cohort gives the following exact minimum discrepancy vector for accepting this zero-response subfamily:

| Observable | Required discrepancy for the zero-response baseline |
|---|---:|
| Maximum shortening |12.73324306 µm|
| Region-averaged shortening |8.38319196 µm|
| Separate apical entry |13 µm|
| Separate basal entry |4.36972161 µm|

At these four allowances every one of the 68 constraints accepts the zero response. Tightening any one allowance below its displayed maximum excludes this zero-response family. **These large values are not inferred error bars and are not proposed as acceptable biological discrepancies.** They simply give an exact executable baseline and expose why silently generous tolerances would make compatibility uninformative.

This establishes one conditional subfamily, not an exhaustive identified set or a fitted nonzero mechanical response. It cannot support a larger-region folding prediction or novelty claim.

## Source inequalities and their limits

In all 17 records, the published numbers satisfy

\[
\text{regional mean}\le\text{maximum shortening}
\le\text{apical entry}+\text{basal entry}.
\]

The first inequality follows for several plausible mean/maximum conventions. The second would follow from pointwise additive component shortening with separate suprema over a shared domain. However, other aggregation/reference conventions can also satisfy the observed inequality. Their satisfaction does **not** identify the source algorithm. The source component aggregation, growth-reference allocation and optical/ROI transport need their own audit rather than inference from these numerical coincidences.

The exact table of 17 differences is retained in the baseline receipt. Four records have equality between maximum shortening and the component sum; the others have a strict gap. Separate maxima must not be added as an identity for the measured maximum.

## What the joint fit must still supply

A real candidate must supply a common material/gain law, allowed geometry and observation nuisance variables, and bounded predictions for each included record/observable. The source meanings of the four quantities must be carried explicitly; a central material-node displacement is not automatically a maximum or a regional mean. A quasistatic path model must specify its path coverage and max/average ordering rather than treating one final field as a measured temporal maximum.

The exact static observer is developed separately, and the [source audit](JOINT-OBSERVATION-CONTRACT.md) records the remaining aggregation, reference and registration limits. The evaluator can test explicitly declared alternatives, but cannot decide which alternative generated the published data. It also cannot turn an exploratory scenario collection into the full physical compatible set. The original identified-set and uncertainty-propagation requirement remains open.
