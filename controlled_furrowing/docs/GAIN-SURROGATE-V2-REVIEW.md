# Independent gate for the extended gain surrogate

**PASS for the frozen v2 defined surrogate.** Extending the verified strength table from 10 to 20 produces interior optima for all four observation alternatives. The review certifies the complete finite surrogate calculation and its source linkage; it does not certify a nonlinear physical gain optimum or complete identified set. Direct nonlinear candidate evaluation may now proceed against the exact reviewed hashes.

The [frozen protocol](../runs/joint_gain_extension_v2/protocol.json), [extended response table](../runs/joint_gain_extension_v2/response_table.json), [candidate certificates](../runs/joint_gain_extension_v2/surrogate_candidates.json), and [independent gate receipt](../evidence/piecewise_gain_minimax_v2_independent_review.json) are separate from the preserved [v1 review](GAIN-SURROGATE-REVIEW.md). No biological fit or solver optimization was rerun by this reviewer.

## Frozen inputs and new-node checks

The table now has strengths 0, 1, 5, 10, 20. Its first four rows are identical to v1, and all parent-artifact, seed-field, source-code and small-data hashes match. The new strength 20 state has the same physical parameters, reaches its requested target, and is bound to the new table by SHA256.

The reviewer re-evaluated the frozen physical potential's gradient at that saved field under both integration policies, without taking an optimization step. The normalized force norms were 5.5149149×10⁻¹² under Duffy8/surface12 and 5.4514844×10⁻¹¹ under Duffy12/surface16; both pass the declared 10⁻⁷ gate. Replaying the frozen exact geometry verifier reproduces the retained positive-J and global-embedding certificates. The observer reproduces all 16 new strength/variant/channel intervals. These checks bind one stored discrete state to the table; they are not continuum error bounds or independent proofs of the verifier algorithms.

The cap is now \(20/\max I_i=31250000/2360271\approx13.2400050672\). The protocol's conditional extension rule stopped at 20 because every variant's restricted surrogate optimum is interior. No strength 40 or 80 field was needed or promoted into this reviewed table.

## Complete exact certificate replay

All 17 records and all four channels remain in the objective. The reviewer reconstructed 43 breakpoints, 42 gain cells per variant, and 23,016 signed affine line coefficients across 168 cell certificates. Every coefficient agrees exactly with the actual source response, dose, channel, interpolation segment and rational midpoint of the node interval.

Each candidate loss and active set passes the appropriate exact convex subgradient or cell-boundary condition. The cells cover the complete declared gain interval, and the reported global result is the minimum over every cell. The optimizer is byte-identical to the one independently checked using exhaustive vertex enumeration on six synthetic cases in v1. No optimizer was invoked to reproduce the actual v2 results: the retained certificates were independently replayed.

| Spatial mean weight | Component aggregation | Interior surrogate gain | Minimum worst residual (µm) |
|---|---|---:|---:|
|Current disk|Spatial maximum|8.062972020|5.246851170|
|Current disk|Regional mean|9.509957310|6.024626246|
|Current radial line|Spatial maximum|7.327159326|5.599341560|
|Current radial line|Regional mean|8.160245983|6.076291390|

All four optima balance two active residuals: overprediction of column L's regional-mean shortening and underprediction of column M's apical entry. This identifies a quantitative tradeoff within the current fixed-scenario surrogate. It is not a demonstrated biological inconsistency: source observation semantics, heterogeneous specimens, nonlinear interpolation, time dependence and numerical error remain separate issues.

## What this gate authorizes—and what it does not establish

The receipt's PASS is bound to the exact protocol, table, candidate and optimizer hashes. It permits direct evaluation of the reviewed candidate-dose strengths under the frozen v2 procedure, retaining failures and all four endpoint intervals. Changing any reviewed artifact requires another audit or a separately versioned result.

An interior surrogate optimum resolves v1's immediate table-edge problem. It does not prove that no better gain exists outside the declared range, that the actual nonlinear response follows the midpoint interpolant, or that a static endpoint equals the experimental temporal maximum. The exact long rational gain strings describe this finite mathematical surrogate, not biological measurement precision.

Direct replay remains necessary. Its completed outputs can assess interpolation discrepancies and candidate residual requirements, but point checks alone will not establish a uniform interpolation bound or global nonlinear optimum. The four conventions retain the [source-observation limitations](JOINT-OBSERVATION-CONTRACT.md). No discrepancy law, noise estimate, full material/geometry identified set, larger-region predictive success or novel concept follows from this gate.
