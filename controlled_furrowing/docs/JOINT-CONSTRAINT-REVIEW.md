# Independent review of the joint interval-constraint evaluator

**PASS for exact interval logic and the declared paired-permutation sensitivity.** The review independently checked all 65,536 bipartite graphs with four vertices on each side, 131 smaller rectangular/empty graphs, 560 interval cases, 21 invalid-input cases and the actual 68-entry zero-response baseline. This does not establish a fitted mechanical model, valid forward enclosures, source-operator equivalence or a physical identified set.

Reviewed files are [the evaluator](../scripts/joint_calibration_contract.py), [its 30-check receipt](../evidence/joint_calibration_contract_verification.json), [the baseline](../evidence/joint_zero_activation_baseline.json), and [the mathematical contract](JOINT-CALIBRATION-CONSTRAINTS.md). The [independent receipt](../evidence/joint_calibration_contract_independent_review.json) retains results, exact source hashes, the original guard findings and a read-only replay program. The [source-observation audit](JOINT-OBSERVATION-CONTRACT.md) remains necessary: exact constraints cannot repair an undefined measurement operator.

## Interval acceptance is conditional and conservative

For one recorded value \(y\), a supplied closed prediction interval \([\ell,u]\), and an explicit allowance \(\delta\ge0\), the evaluator returns:

- **IN** when the entire interval lies in \([y-\delta,y+\delta]\).
- **OUT** when the two intervals are disjoint.
- **UNRESOLVED** otherwise.

The review checked 560 combinations of signed observations, ordered integer endpoint pairs and rational allowances, including equality at boundaries. This independently reproduces both the minimum and maximum possible absolute residuals. Joint canonical membership is OUT if any included constraint is OUT, IN if every constraint is IN, and unresolved otherwise.

IN is sound **if** a fixed candidate's actual output is enclosed by each supplied interval and the measured quantities match. It is not a test that those premises hold. For example, the coupled output vector \((t,1-t)\), \(0\le t\le1\), has both scalar coordinate intervals overlapping zero, but no jointly zero response. The evaluator correctly leaves zero-tolerance membership unresolved rather than declaring a realizable fit.

The generic function requires prediction IDs and observable keys to match the observation mapping supplied by its caller. It does not itself prove that an arbitrary caller supplied the authoritative 17-record cohort. The source loader and actual baseline test establish that the recorded application includes all 17 records and all 68 entries.

## Unknown pairing preserves intact component pairs

The optional alternative keeps the dose/maximum/mean association fixed and permits one bijection of **paired** apical/basal source columns. A possible edge means both component intervals intersect their respective allowed bands. A guaranteed edge means both intervals are wholly inside those bands.

If fixed maximum/mean constraints are IN and a guaranteed-edge perfect matching exists, one fixed permutation works for every output enclosed by the intervals: conditional IN is justified. If a fixed constraint is OUT, or even the possible-edge graph has no perfect matching, OUT is justified. A possible matching without a guaranteed one remains unresolved. Unknown measurement identities do not justify independently permuting apical and basal values.

The review explicitly verified:

- A complete swap of two intact component pairs is canonical OUT but optional-permutation IN.
- Observed pairs \((0,1)\) and \((1,0)\) cannot be used to manufacture the pair \((0,0)\) by separate component shuffles: OUT.
- Two predicted records competing for one possible observed pair yield a genuine Hall deficit: OUT.
- Wide component intervals with only possible matches remain unresolved.
- Fixed maximum/mean OUT overrides a component matching; fixed maximum/mean unresolved prevents IN.

This is a conservative identity-nuisance sensitivity, not verification of specimen identity or evidence that the author's column order was wrong. The canonical assumption remains distinct.

## Independent graph verification and Hall witnesses

The root implementation uses augmenting paths. The review instead computed maximum cardinality from the largest Hall deficiency by enumerating all subsets of left vertices:

\[
\nu(G)=|L|-\max_{S\subseteq L}(|S|-|N(S)|).
\]

For each of the 65,536 four-by-four graphs, the resulting cardinality agrees with the root implementation. Every returned matching uses distinct right vertices on valid edges. For all 27,713 deficient graphs, the reported Hall set has exactly the stated neighborhood and positive deficit. Another 131 rectangular/empty graphs also agree. This supplies a separate finite reference calculation, rather than rerunning the implementation's own search as its oracle.

The underlying matching and Hall-obstruction principles are standard. Exhaustive small-graph tests are implementation evidence; they are not experimental validation or a novelty claim.

## Input guard findings were corrected

The first review found two input-domain problems: a string such as `'12'` could be interpreted as a two-endpoint interval, and booleans or in-range noninteger vertices could pass the matching range check. A Boolean could silently act as an integer vertex; a floating vertex could fail later with an indexing error.

The corrected code requires an interval list/tuple of two ordered finite rational values, a non-Boolean integral right-vertex count, and non-Boolean integral vertices checked before sorting/deduplication and range validation. The review tested these corrections, negative/out-of-range vertices, malformed rows, nonfinite/Boolean values, missing/extra record IDs, missing observable keys, missing allowances and negative allowances. All 21 cases reject with ValueError. Negative **predicted responses** remain allowed; signed motion is meaningful and must not be rejected merely to make a source inequality hold.

## Actual source fractions and the zero-response subfamily

The independent calculation reconstructs the 17 records from the frozen small CSV using exact fractions. Its zero-response thresholds agree exactly with the retained baseline:

| Observable | Exact decimal threshold (µm) |
|---|---:|
|Maximum shortening|12.73324306|
|Region-averaged shortening|8.38319196|
|Apical entry|13|
|Basal entry|4.36972161|

At these allowances, all 68 zero-response constraints are IN under either pairing convention. Tightening any one allowance by 10⁻⁸ µm makes the result OUT, independently checked for both conventions. These are the maximum observed values and a computational baseline—not estimated errors or recommended tolerances.

The baseline's mechanical interpretation is conditional on the declared flat fixed-rim, zero-gain, zero-reference model. Its variational argument is consistent: the bulk singular-value contributions and area-spring penalties are nonnegative; fixed-rim membrane area cannot fall below the projected disk; nonnegative tensions therefore cannot lower the total below the flat reference. This does not make zero response an adequate biological model or produce a nonzero transfer prediction.

The module performs no optimization, chooses no discrepancy law and validates no continuum response. It supplies a checked membership predicate that a subsequent physical inference calculation may use after its observation definitions and forward bounds are justified.
