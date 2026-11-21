# Observer bounds at four exact-spatial roots

At the matched large input, the two frozen finite-element models have a positive signed difference in the apical contrast. The independently verified gap is **0.15475515613497 to 0.15475515719595 µm**. Both individual contrasts are negative. The two verified small-L roots account for eight primary comparisons and 24 alternative comparisons. They provide **28 IN and 4 OUT** at the frozen allowance of approximately **5.389835675 µm**.

## Same sign, different magnitude

The large observable is \(Q=Z_a(0)-Z_a(R_s)\), with the axis at zero and the shoulder at current radius \(R_s=24.94698765\,\mu\mathrm m\) in both fields. A negative value means the axis lies below the shoulder in this axial coordinate. The labels + and − identify the two original affine search directions. The contrast result gives a difference in magnitude and no evidence of opposite folding directions or biological significance.

![On a magnified range of the negative Q axis, the minus candidate is at minus 4.290758 micrometers and the plus candidate is at minus 4.136003 micrometers. The difference is about plus 0.154755 micrometers.](exact_spatial_contrast.png)

| Large-state value | Lower bound (µm) | Upper bound (µm) |
|:--|--:|--:|
| Candidate + | −4.136003261986 | −4.136003261301 |
| Candidate − | −4.290758418498 | −4.290758418120 |
| Signed gap | 0.154755156134 | 0.154755157196 |

The interval widths are smaller than the plot resolution. The numeric bounds round outward. The signed difference adds the uncertainty from both root observables through interval subtraction:

\[
[Q_{+,l}-Q_{-,u},\ Q_{+,u}-Q_{-,l}] .
\]

## Comparison at the two small L states

Each line includes both candidates and four small-L observables. IN requires the entire residual interval to stay at or below the allowance. OUT requires its lower bound to exceed the allowance.

| Static convention | IN | OUT |
|:--|--:|--:|
| Disk · component maxima | 8 | 0 |
| Disk · component means | 8 | 0 |
| Radial line · component maxima | 6 | 2 |
| Radial line · component means | 6 | 2 |

All four OUT rows are the regional mean of total shortening under radial-line weighting. This mean is shared by the two radial-line conventions. The exact residual lower bounds at L are **6.07609207 µm** for + and **6.19199947 µm** for −. The allowance is a conditional model discrepancy threshold. The experiment’s measurement uncertainty, time aggregation and reference registration remain open.

## How the root certificate reaches an observation

Let \(e\) bound the free-coordinate distance of the spatial root from the exact binary center and let \(\xi\) bound the coordinate addition error in the saved field. The independent readback reconstructed both quadratic traces from the stored nodes, confirmed fixed axis and rim radii, and checked the radial derivative against \(8\epsilon\) on every edge for the new \(\epsilon\). A face with stored slope bound \(M_s\) then has Eulerian height error \(E_s\) given by

\[
\epsilon=33.5(e+\xi)\,\mu\mathrm m,
\qquad E_s=\tfrac54(1+M_s)\epsilon .
\]

For total shortening, the two face errors add. A component maximum or mean uses its own face bound. The axis-node error \(\epsilon\) adds to the apical shoulder error \(E_a\) for \(Q\). This retains the fixed observer, including evaluation at the same current radius as the curve changes.

## What the gap would imply under full compatibility

If both models belong to the genuinely calibration-compatible set, every single-point prediction \(q\) for their large observable incurs worst-case absolute error at least half their gap. The triangle inequality gives the classical two-value bound

\[
\max\{|q-Q_+|,|q-Q_-|\}\ge\tfrac12|Q_+-Q_-| .
\]

With these two Q bounds that conditional lower bound would be **at least 0.07737757806748 µm**. The four spatial certificates alone cover only eight of 136 primary comparisons. The other 128 and their associated 32 exact-spatial roots remain unverified at this scope. The separate 36-state result uses cached quadrature. Source compatibility, a stable loading path, continuum adequacy and scientific novelty remain unresolved.

## Evidence

- [Independent four-root observer receipt](independent_actual_readback/receipt.json) · [Observed independent replay exit](independent_actual_readback/observed_os_exit.json)
- [Four exact-spatial root and minimum certificates](../four_spatial_checkpoint_v1/SPATIAL.html) · [Spatial proof readback](../bridge_reduction_v1/independent_actual_preparation/all_four_readback.json)
- [Full 36-field cached-quadrature observer checkpoint](../../full_outcome_coverage_v1/observer_checkpoint_v1/OBSERVERS.html) · [Nonlinear stored-field report](../../../report/REPORT.html) · [Independent stored-field convention sweep](../../../observer_robustness/independent_actual/receipt.json)
