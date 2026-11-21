# Independent all144 uniform-bulk wrapper prelaunch review

The revised preparation passes1238 checks. This is permission for one frozen producer invocation, not a statement that the137 new local integrations have completed. No nonlinear local integral or64-site cache force/Hessian sum was run by this reviewer.

The original unsolved preparation consumed137 additional corrected point diagnostic files without binding them. The parent preserved that preparation, added all144 consumed point-file hashes with conflict-strict merging, and strengthened the review schema check before this review. Those point packets supply overlap/map diagnostics, not the numerical uniform bounds themselves. The original finding is retained in initial_findings.json.

## Local-to-global argument

For each material element, the reviewed SAFE V2 component encloses its force discrepancy and the upper triangle of its symmetric Hessian discrepancy on the same fresh common nodal box and the same parameter cell. Let its Boolean projection into the611 free coordinates be \(E_e\). The bulk sum is

\[
\Delta F_{\rm bulk}=\sum_eE_e^T\Delta F_e,\qquad
\Delta H_{\rm bulk}=\sum_eE_e^T\Delta H_eE_e.
\]

The wrapper sums signed interval endpoints in these common coordinates. Sorting each global index pair handles reversed local/global ordering; inserting each assembled off-diagonal in both symmetric positions and each diagonal once reconstructs the symmetric matrix. Sharing uncertain nodes between elements can make the interval sum wider, but cannot invalidate its containment.

For the fixed exact-binary preconditioner \(R\), interval matrix multiplication encloses every entry of \(R\Delta H\). Thus

\[
b_1=\max_i\sum_j\sup\lvert(R\Delta H)_{ij}\rvert
\]

is an upper bound for its infinity norm. A second bound is

\[
b_2=\max_i\sum_k\lvert R_{ik}\rvert
        \sum_j\sup\lvert\Delta H_{kj}\rvert.
\]

Both calculations use nonnegative absolute upper bounds; therefore \(\min(b_1,b_2)\) remains a valid upper bound. The same signed interval multiplication encloses the preconditioned force discrepancy. These are bulk-only bounds. Surface contributions, the exact-integral root-transfer conditions, a PDE error estimate, and biological/source-time validation remain separate obligations.

## Verification receipts

- All144 selectors are unique and exhaustive, and every fresh-geometry \(J\) lower bound exceeds the inactive-penalty threshold.
- Each scalar analytic-domain radius is independently recomputed by outward dyadic rounding. Every first adequate derivative-series order is independently reconstructed from exact harmonic-number remainder formulas, both positive measure masses, and polynomial coefficient majorants.
- Every corrected point packet has the expected complete free map and upper-triangle indices. The seven reused uniform packets have exactly matching selectors, domain lineage, and bytes.
- The actual local assembly and norm-calculation AST blocks were executed on a three-coordinate toy (only the constant611 was replaced by3). All64 exact rational Hessian endpoint choices and16 force endpoint choices are contained; reversed mappings, shared diagonal addition, signed cancellation, symmetry and both norm bounds are covered.
- The actual prelaunch guard is tested with matching hashes and rejected false/truthy/mismatched cases. Every frozen source/input hash is checked.

Authoritative receipt: [producer prelaunch review](../uniform_bulk_assembly_B_safe_v2/prelaunch_independent_review.json), SHA25699a8fa459fd331fe78392703d4d71b9cb0c9a9c08a07fc8be6c124c397648690. The retained review implementation and duplicate review output are in this folder.
