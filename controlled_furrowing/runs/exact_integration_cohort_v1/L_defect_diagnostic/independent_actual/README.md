# Independent postrun review of L element 21

The producer finished with exit code 0 and the exact result digest in its terminal receipt. The review checked 47 frozen source and input bindings. It reconstructed the interval arithmetic for all 12 local force coordinates with exact rational endpoints. [Review](../actual_independent_review.json) · [Endpoint replay](endpoint_replay.json)

| Quantity | Largest across 12 coordinates |
|---|---:|
| Old whole-box interval width | 5.23196888e-8 |
| New direct whole-box interval width | 5.23196888e-8 |
| Centered whole-box interval width | 9.32225561e-12 |
| Fixed-center point interval width | 3.13136626e-19 |
| Fixed-center absolute defect upper bound | 2.17060600e-14 |
| Absolute state-increment upper bound | 4.66112764e-12 |
| Absolute gain-increment upper bound | 1.98976229e-25 |

Values are rounded displays of normalized force bounds. The largest whole-box width decreases by about 5,612×, with the maximum at free coordinate 287 in both enclosures. All 12 centered intervals lie inside both direct whole-box enclosures and equal the reported intersection of the new direct and centered intervals. The original direct interval contains the new direct interval in only one coordinate. All 12 old/new pairs overlap.

At the exact proof center and base gain, the actual point intervals certify five positive and seven negative local defect components. These signs come from the signed endpoint enclosures. Every centered whole-box interval still contains zero. The direct whole-box width is dominated by the volumetric truncation contribution's interval width. At the point, the finite cache correction has the largest absolute bound among the three signed terms. These two comparisons measure different things.

The centered check uses the point defect, the previously bound whole-box Hessian defect at the base gain along the segment to the current node vector, and the affine-in-gain increment at the current node vector. The review recomputed endpoint products and sums and replayed local derivative norms for the two-mass scalar tails. It relies on the old local Hessian packet and the prelaunch-reviewed moment, cache, and gain derivative primitives. It did not evaluate a new force integral, Hessian moment, weighted cache derivative sum, or global reduction. The result is local to element 21 and establishes no whole-record force or root conclusion.
