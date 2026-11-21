# Independent review: hypothetical calibration exceptions

**PASS for the declared finite scalar problem.** All 131,072 subset costs and all 17 exception-budget optima were independently verified. The minimum uniform discrepancy falls only from 3.131943715 µm to 3.084159745 µm when one arbitrary record may hypothetically be excepted. No record is classified as a biological outlier, and the actual mechanical calibrations retain all 17 records.

This reviews [the implementation](../scripts/calibration_contamination.py), [root results](../evidence/calibration_contamination.json), and [scientific interpretation](OUTLIER-SENSITIVITY.md). The independent [review receipt](../evidence/calibration_contamination_review.json) records source hashes, every budget result, adversarial guard tests, and a self-contained replay program. No larger-region outcomes were read. The figures and primary observation definitions are documented separately in [the calibration uncertainty analysis](CALIBRATION-UNCERTAINTY.md).

## Exact problem and scope

The verified inputs are the 17 frozen CSV records, sorted by strictly distinct **positive** measured doses, with nonnegative maximum-shortening responses. Numerical exactness concerns those recorded decimals, not biological measurement precision. Dose order is fixed. For an exception budget \(k\), retain at least \(17-k\) records and minimize the largest residual to one nonnegative nondecreasing curve. A zero-response origin can be added without cost because all retained doses are positive: interpolate monotonically from \((0,0)\) to the first fitted value, between retained values, and then use a constant tail.

No concavity, dose uncertainty, specimen-specific material variation or physical constitutive law is imposed. The common monotone curve is a conditional null class, not an established exact biological law. The nested-set peak theorem can motivate this class only under its separately stated common-branch and path-coverage assumptions.

## Independent methods

The root implementation builds each subset's maximum violating drop from a smaller prefix subset, then independently checks interval-subsequence feasibility with a dynamic program. This review did not use either implementation to establish the reported optima.

First, the review generated every subset mask independently. For each ordered source pair, it located **all** masks containing both observations and updated their costs with that pair's positive response drop. The resulting 131,072-entry integer array matches the retained array entry for entry. At the verified scale of \(10^8\) units per micrometer, all integer costs fit safely in the recorded 64-bit representation.

Second, for each candidate doubled discrepancy \(t=2\delta\), it constructed a conflict graph: records \(i<j\) conflict exactly when \(y_i-y_j>t\). A retained set admits a monotone fit within \(\delta\) precisely when it is an independent set. The necessary pair condition follows from monotonicity; sufficiency follows by taking the nonnegative prefix maximum of the lower observation-interval endpoints. The review computes the largest independent set by a separate exact include/exclude recursion:

\[
\alpha(G)=\max\{\alpha(G-v),\;1+\alpha(G-v-N(v))\}.
\]

This algorithm branches on a graph vertex, rather than propagating minimal interval tails. It agrees with every reported feasible retained count and every predecessor-threshold rejection.

Finally, every reported matching construction was checked against the source record, exact dose, exact response, nonnegativity, retained order, residual bound, and retained/omitted partition. These are finite scalar constructions, not solutions of the tissue-mechanical model.

## Verified results

| Hypothetical exceptions | Retained records | Sharp uniform discrepancy (µm) |
|---:|---:|---:|
|0|17|3.131943715|
|1|16|3.084159745|
|2|15|1.695976145|
|3|14|1.502609595|
|4|13|1.387728410|
|5|12|1.241829440|
|6|11|0.807538065|
|7|10|0.634125125|
|8|9|0.346913030|
|9|8|0.064503065|
|10|7|0.056135720|
|11–16|6–1|0|

For every positive optimum, the preceding distinct critical threshold admits fewer than the required number of records. At zero discrepancy, the largest compatible subset contains six records. That combinatorial fact does **not** imply that the other 11 embryos are erroneous; ordinary biological variability or measurement noise can create order reversals.

The root's one-exception matching example omits source-column I. Record M then has response 12 µm at dose 0.78208913, preceding record P with response 5.83168051 µm at dose 0.80580697. Their reversal alone forces half their 6.16831949 µm response drop, or 3.084159745 µm. Exhaustive verification establishes that no other single omission produces a smaller optimum. This is a mathematical sensitivity witness, not grounds for deleting I.

## Why the predecessor check proves sharpness

Every subset's optimal doubled discrepancy is zero or one of its positive ordered response drops. Hence a global optimum lies in the finite set of all such critical values. Between consecutive critical values, the conflict graph does not change. If the exact maximum independent-set size is sufficient at a claimed threshold and insufficient at the previous critical threshold, no intermediate value can achieve the required count. When the optimum is zero, nonnegative discrepancy itself supplies the lower bound. This covers all budgets, including plateaus and the zero-error tail.

## Guard checks and limitations

The implementation now explicitly rejects nonpositive source doses, duplicate measured doses, negative response inputs, and negative discrepancy arguments to its feasibility helper. Seven adversarial tests passed: three helper tests and four synthetic-input tests of the main entry point. The invalid synthetic data were supplied in memory; original data were unchanged. This review does not certify a duplicate-dose algorithm or an unrestricted negative-response formulation.

The conclusion is narrowly useful: the all-data incompatibility with a common monotone curve is not eliminated by one exceptional record. It neither estimates a true contamination fraction nor establishes measurement noise, confidence coverage, biological mechanisms, novelty, or a physical parameter identified set. All actual mechanical fits must continue to use every record unless a separate, scientifically justified exclusion protocol is adopted.

Reviewed implementation SHA256: `493bfff54c97c17e8886fbe000cd773d15176d0f2436ebe7459c83b03a4ddce0`. The receipt identifies the exact source, root-output and subset-array hashes, so later implementation changes require a fresh review rather than inheriting this PASS automatically.
