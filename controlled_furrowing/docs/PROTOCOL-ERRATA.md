# Protocol wording corrections, without changing frozen calculations

Frozen inputs and code are preserved so historical hashes remain verifiable. These corrections clarify interpretation; they do not change fitted values, numerical thresholds or predictions.

## Surface-tension revision: reference versus prestress

`SPATIAL-TENSION-PROTOCOL.json` was derived from v1 and retains the phrase “flat stress-free idealization” in its geometry description. For nonzero basal tension, the correct description is **flat equilibrated reference with a stress-free bulk and a prestressed basal membrane, balanced by the clamped rim's reactions**. The constant surface tension has zero interior first variation on a flat membrane, but it is not zero stress. The actual executed energy contains \(T_b A_b\), as proved and checked in [BASELINE-TENSION](BASELINE-TENSION.html).

The zero-tension member is stress-free in the original sense. Curved references with nonzero tension are rejected by this revision until balanced loads are specified. Source calculations, generation receipts and frozen protocol hashes remain unchanged; subsequent protocols must use the corrected wording.

## Volume-change coefficient is not an incompressibility certificate

The parameter called `bulk` in solver APIs multiplies \((\log J)^2/2\). Near identity it is Lamé \(\lambda\); the infinitesimal bulk modulus is \(\lambda+2\mu/3\). A large ratio to shear modulus alone does not establish small volume change under large active loads. Report actual local Jacobians as well as global volume, and distinguish conservative Bernstein bounds from evaluated Jacobian values.

## Prescribed potential work is not calibrated optical energy

The active-potential parameter is externally prescribed. Its work ledger is mechanically consistent for the declared potential gauge. It does not quantify photons, ATP use, biochemical conversion efficiency or a finite energy reserve. Adding a parameter-dependent constant changes that ledger without changing the mechanical forces. None of these embryo models inherits the finite-battery claims of the earlier, separate survival project.

## Observational versus numerical uncertainty

The 0.05-µm tolerance in the scalar-to-spatial calibration workflow checks a surrogate against direct solves at selected training points. It is not a biological measurement uncertainty or a uniform continuum error bound. The 0/0.05/0.2-µm surface-error scenarios in the v1 label diagnostic are assumptions, not certified errors. Unknown slice/time mapping remains a separate source of uncertainty.
