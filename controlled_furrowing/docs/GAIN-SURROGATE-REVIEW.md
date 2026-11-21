# Independent review of the first joint-gain surrogate

**PASS for the exact defined surrogate; its gain support is exhausted.** All four observation variants attain their optimum at the table-imposed cap \(g=15625000/2360271\approx6.6200025336\). This is not a physiological bound or a physical gain estimate. Extending the response table is necessary to investigate the boundary candidate, while the exact result on the declared interval remains valid.

The reviewed [protocol](../runs/joint_gain_calibration/protocol.json), [response table](../runs/joint_gain_calibration/response_table.json), [candidate certificates](../runs/joint_gain_calibration/surrogate_candidates.json), [optimizer](../scripts/piecewise_gain_minimax.py), and [independent receipt](../evidence/piecewise_gain_minimax_independent_review.json) are retained separately. No actual fitting run or root writer was rerun for this review. Pure synthetic examples were checked against a different exhaustive algorithm; stored-state observer extraction was replayed only to verify table provenance. No held-out response was inspected.

## What mathematical object was optimized?

For each of four static observation conventions, the table gives closed endpoint intervals at activation strengths 0, 1, 5, 10. Their exact rational midpoints define numbers \(v_{jk}\). Linear interpolation through those numbers defines a **new, explicitly specified surrogate** \(\widetilde q_k(f)\). It is not a proved interpolant or enclosure of the nonlinear mechanical response between strength nodes.

The objective uses all 17 small-region records and all four reported channels:

\[
L(g)=\max_{i=1,\ldots,17;\ k=1,\ldots,4}
|\widetilde q_k(gI_i)-y_{ik}|,
\qquad 0\le g\le\frac{10}{\max_i I_i}.
\]

Thus all 68 entries participate with a common micrometer worst-residual objective. This is a diagnostic optimization criterion, not a measurement-noise model or a chosen biological error tolerance. The gain cap prevents table extrapolation; it is determined by the largest recorded dose 1.51057344, not by independent physiology.

The alternatives are current-radius disk versus radial-line averaging, each with either spatial-maximum or regional-mean component aggregation. They retain a flat no-growth reference, canonical source-column pairing and a static endpoint used as a peak proxy. The [source observation contract](JOINT-OBSERVATION-CONTRACT.md) explains why these remain assumptions, rather than a uniquely recovered source algorithm.

## Exact global certification over the declared interval

Every possible change of interpolation segment occurs at a gain \(f_j/I_i\). Together with the two bounds, these yield 28 distinct breakpoints and 27 positive-width closed cells. The review reconstructs all of them from the actual source doses; no cell is omitted.

On one open cell, each record/channel prediction has fixed affine dependence on gain. Therefore

\[
L(g)=\max_m(a_mg+b_m)
\]

on the closed cell, with 136 signed residual lines and the zero line. Its restriction is convex even though the full piecewise objective need not be globally convex across cells. The review independently reconstructed **14,796 exact line coefficients across 108 cells** and matched them to the retained certificates. Every coefficient uses the correct source record, channel, strength segment and interval midpoint.

At a candidate \(g_*\), let \(A_*\) be the set of active lines. For a maximum of affine functions, the left and right derivatives are respectively the minimum and maximum active slopes. A sufficient and necessary closed-cell minimum condition is:

\[
\begin{cases}
\max_{m\in A_*}a_m\ge0,&g_*\text{ at the left boundary},\\
\min_{m\in A_*}a_m\le0,&g_*\text{ at the right boundary},\\
\min_{m\in A_*}a_m\le0\le\max_{m\in A_*}a_m,&g_*\text{ in the interior}.
\end{cases}
\]

Every retained candidate, loss, active set and slope condition passes exact rational replay. Taking the smallest certified cell loss therefore proves a global minimum on the union of all cells. Floating-point linear programming proposes a candidate; the exact certificate, not the floating optimum, establishes the result. This is standard piecewise-affine convex certification, not a new optimization theorem.

## Results and active constraints

| Mean weighting | Component aggregation | Exact-surrogate minimum residual (µm) | Active limiting source entry at the cap |
|---|---|---:|---|
|Current disk|Spatial maximum|5.938104366|Column M, apical entry|
|Current disk|Regional mean|7.271365271|Column G, apical entry|
|Current radial line|Spatial maximum|5.938104366|Column M, apical entry|
|Current radial line|Regional mean|6.782283545|Column G, apical entry|

All four gains equal the same upper cap. The active residual is observed apical shortening minus the surrogate prediction; its slope is strictly negative at the cap from the left. This identifies which small-data summary limits this fixed-scenario surrogate fit. It does not establish that the original biological component measurement uses that convention.

A boundary optimum does **not** prove that the unrestricted physical optimum lies beyond the cap: the uncomputed response could change behavior. It does establish that the current finite table cannot support an interior gain estimate. The table extension and later direct evaluations must remain separately versioned; they must not overwrite this result or reinterpret its restricted domain.

## Provenance and independent checks

The protocol, table and candidate hashes match. All 11 frozen source hashes and the small CSV hash match current files. The strength 1, 5, 10 state hashes match the already retained higher-quadrature equilibria. Replaying the frozen observer on those three states and the zero state reproduces all 64 strength/variant/channel endpoint intervals. This binds table values to their source fields; it is not an independent re-proof of the observer's geometry. A separate 21-check observer review supplies that component's independent evidence.

Six synthetic cases were independently solved by enumerating every cell endpoint and every pairwise line intersection. They include a nonmonotone response, nonlinear sampled multichannel values, a flat minimizer interval, multiple equal minima and rational breakpoints. All reported optimum values agree exactly. Seven invalid-input cases are rejected. Synthetic tests invoke no root receipt writer and no biological fit.

The exact midpoint definition is important. Small endpoint intervals describe stored discrete-curve quantities at the sampled strengths. Their midpoints do not bound nonlinear interpolation error, continuum error, experimental temporal peaks or unresolved observation-reference choices. Even interpolating interval endpoints would not automatically enclose the actual nonlinear response between nodes.

## What remains before a physical inference claim

Extend table support under a new frozen version, independently certify the resulting surrogate candidates, and evaluate the actual nonlinear fields at every requested candidate-dose product. Preserve requested exact strengths, actual numerical strengths, failed attempts and all four observation intervals. Compare direct outputs against the surrogate rather than promoting the table optimum to a nonlinear optimum.

Even successful direct verification establishes only responses at those candidate gains in one fixed mechanical scenario. It does not prove global optimality over nonlinear gain, material/geometry combinations or observation conventions. It supplies neither a discrepancy law nor the full compatible set required by the original plan. At this reviewed version, the surrogate calculation is complete; support for an interior or physical gain estimate is incomplete.
