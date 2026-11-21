# Independent review: complete B-center bulk assembly

**PASS157checks.** After the producer's true144-element terminal result, the independent reduction terminated with exit0. The frozen main receipt is `evidence/exact_bulk_assembly_B_independent_review.json`, SHA2565e32293697e1bc398c2d5037dcb0438d5d11ae95905e152cde2290ff3dfc06ae. A separate eight-check followup binds the prior four independent replays directly to the new adaptive-order element packets.

## What was checked completely

All144local packet hashes, element IDs, truncation orders, rho values, local free maps and force/Hessian tail targets were checked. Exact Fraction arithmetic sums their signed interval endpoints into all611global force intervals and6676unique upper-triangular Hessian slots. The producer's global intervals contain these exact endpoint sums. Unsorted local-to-global indices and shared nodes are retained; local upper Hessian entries are sorted into global pairs and included once in each symmetric row.

The611×611preconditioner is decoded from its hash-bound binary cache. Every entry becomes its exact binary rational, including its sign. All611signed R*DeltaFintervals are independently computed by rational interval matrix multiplication; every producer interval contains the corresponding independent enclosure. The force norm, preconditioned force norm, symmetric Hessian row norm and conservative product ||R||infinity*||DeltaH||infinity are bounded by the reported producer quantities. There is no floating tolerance or midpoint-only comparison.

| Exact-rational upper-bound calculation | Approximate display |
|---|---:|
| ||R||infinity |1686.93121545|
| Bulk force-defect infinity norm |2.76009685×10^-11|
| Signed preconditioned bulk force-defect infinity norm |1.67368904×10^-11|
| Bulk Hessian-defect row-sum upper bound |6.62364706×10^-8|
| Conservative preconditioned bulk Hessian-defect upper bound |1.11736370×10^-4|

Exact fractions and every assembled interval endpoint are retained in `exact_reduction.json`. “Exact reduction” means exact rational interval arithmetic applied to already validated local enclosures. It does not mean exact evaluation of the integrals or a sharp attainable norm over their correlations.

## Independent local integration coverage

Two additional elements,48and143, were fixed before reading their producer results and independently integrated at order18/640bits under a frozen protocol. Their22force entries,133upper Hessian entries and7040exact coherent-basis comparisons passed. They use independently constructed cubic-Jacobian jets, formal coefficient convolution, direct barycentric coherent maps and manual constitutive lifts, with explicit reuse of the reviewed moment/tail/dual/constitutive kernels.

The earlier four independent replays, elements0,6,72,66, were then explicitly compared with the new adaptive-order outputs. All38force and210Hessian interval overlaps passed and the original source/cache contexts remained unchanged. Hence six distinct local elements have independent numerical replay coverage tied to this assembly. The remaining138were not independently reintegrated: they use the mathematically reviewed generic local engine, complete packet/tail/hash/free-map checks and the independent full endpoint reduction. This distinction is part of the receipt.

## What the result does not establish

This is the bulk contribution at one saved center of one record. It includes no surface discrepancy, parameter/root-tube uniformity, full17record bridge or exact-integral root/minimum theorem. In particular, the bulk preconditioned-force upper bound exceeding an old radius does **not** prove that the old cube fails: an upper bound is not a displacement lower bound, surface terms can cancel bulk terms, and the complete residual/self-map condition has not been evaluated. Point Hessian discrepancy is also not a uniform tube bound.

The series-tail targets constrain omitted-series errors only; total cache/integration defects and arithmetic width are retained separately. No biological error model, continuum/PDE estimate, source temporal validation or scientific novelty is inferred. Original project completion remains outstanding.
