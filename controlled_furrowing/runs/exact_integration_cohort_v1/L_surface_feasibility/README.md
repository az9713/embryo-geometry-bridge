# L surface feasibility on the original proof box

All 24 closed surface edges have certified positive hoop ratio and squared surface Jacobian G on the coefficientwise box built from L's original binary proof center, its affine predictor over the unchanged parameter cell, and its original normal corrector radius of about 6.01172e-11. The smallest edge hoop lower bound is 0.9211355, and the G enclosures span approximately [0.8742944, 1.6128926] across all edges.

Every edge attains the B selector target of rho ≤ 1/8 with at most five Bernstein subdivisions in depth. The largest certified global rho is 0.12444073 at zero-based apical edge 11. The order-28 selector passes both two-mass tail limits on all 24 edges. The most conservative order-28 force and Hessian-entry bounds across all edges are approximately 1.001e-27 and 1.106e-24, respectively.

| Face | Edge | Subdivision depth | Global rho bound | Minimum force order | Minimum Hessian order | Minimum joint order |
|---|---:|---:|---:|---:|---:|---:|
| Basal | 0 | 0 | 0.00809666 | 9 | 8 | 9 |
| Basal | 1 | 4 | 0.11926254 | 22 | 20 | 22 |
| Basal | 2 | 3 | 0.12154835 | 22 | 20 | 22 |
| Basal | 3 | 4 | 0.09903134 | 20 | 18 | 20 |
| Basal | 4 | 5 | 0.08270245 | 19 | 17 | 19 |
| Basal | 5 | 4 | 0.06729188 | 18 | 15 | 18 |
| Basal | 6 | 3 | 0.08653202 | 20 | 17 | 20 |
| Basal | 7 | 3 | 0.06424961 | 17 | 15 | 17 |
| Basal | 8 | 2 | 0.10040086 | 21 | 18 | 21 |
| Basal | 9 | 2 | 0.08368767 | 19 | 16 | 19 |
| Basal | 10 | 2 | 0.07184603 | 18 | 16 | 18 |
| Basal | 11 | 2 | 0.06525574 | 18 | 15 | 18 |
| Apical | 0 | 0 | 0.01186889 | 10 | 8 | 10 |
| Apical | 1 | 4 | 0.12192896 | 22 | 20 | 22 |
| Apical | 2 | 3 | 0.12105751 | 22 | 19 | 22 |
| Apical | 3 | 3 | 0.09894889 | 20 | 17 | 20 |
| Apical | 4 | 5 | 0.09233862 | 20 | 17 | 20 |
| Apical | 5 | 4 | 0.07119174 | 18 | 15 | 18 |
| Apical | 6 | 3 | 0.09517694 | 20 | 17 | 20 |
| Apical | 7 | 3 | 0.08234684 | 19 | 16 | 19 |
| Apical | 8 | 2 | 0.10571638 | 21 | 18 | 21 |
| Apical | 9 | 2 | 0.08583248 | 19 | 16 | 19 |
| Apical | 10 | 2 | 0.07222139 | 18 | 15 | 18 |
| Apical | 11 | 1 | 0.12444073 | 23 | 20 | 23 |

The selector searched orders 0–128 against a force limit of 1e-22 and a Hessian-entry limit of 1e-16. It used the exact normalized auxiliary zeroth mass and an enclosure of the analytic normalized material zeroth mass. For each edge it used the frozen physical-radius derivative polynomials and bounds for the same face's full modulus and tension intervals. Minimum orders refer to those sufficient tail bounds. At the limiting apical edge, order 22 yields a force bound of 3.02249e-22, and order 23 lowers it to 3.68284e-23.

The run took 5.32 seconds and terminated with exit code 0 and zero failed edges. That timing covers geometry and scalar tails. A future integration adapter could keep order 28 on this domain or use the per-edge passing orders to reduce series work. Any translated predictor or common cube needs its own domain checks. Spatial series, moment defects, finite-cache corrections, global matrices and exact-integral root certification are outside this receipt.

The [result](results.json) has the exact per-edge domains, mass and derivative inputs, pass/fail comparisons, and timings. The [protocol](protocol.json) binds 18 source and input files, including the frozen target-preparation record, selected L cache and predictor, and prior L proof context. The [start](started.json) and [terminal](terminal.json) receipts bind the execution and result.
