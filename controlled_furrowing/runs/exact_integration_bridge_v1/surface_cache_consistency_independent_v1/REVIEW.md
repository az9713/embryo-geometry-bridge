# Independent local surface cache-consistency review

**PASS for the inspected unweighted local formula and two common-box sites.** This does not certify all 24 surface edges, an integrated error, an assembled Hessian, a root, or completion of the five-gate research objective.

Frozen producer: `scripts/exact_surface_cache_consistency.py`, SHA-256 `02c39fa9dfa459a1cc174211a4731dc54966e6aeb12a0133fd9129b100644d2e`. Independent receipt: `evidence/exact_surface_cache_consistency_independent_review.json`, SHA-256 `a40a0b4833044cc6c8d04957edfa831506f4ad3416b3d5c5e88aabcd5dfd8b3e`.

The single completed independent review run exited zero, retained 1,960 checks and 393 source/input bindings, and took approximately 8.23 seconds. Two fresh actual sites were checked: basal face 0 / edge 0 / auxiliary point 0, and apical face 1 / edge 11 / auxiliary point 11. Each was tested against 256 combinations of local nodal-box vertices and shared material endpoints. These endpoint tests supplement the analytic interval argument; they are not, alone, a proof covering the continuum of inputs.

## Mathematical checks

For canonical horizontal reference length \(L>0\), the transformation
\[
D=\operatorname{diag}(\alpha L,L^{-1},L^{-1}),\qquad
\alpha=(r_{0c}L_{0c})^{-1}
\]
preserves the old scalar surface Jacobian exactly:
\[
(Dq)_0\sqrt{(Dq)_1^2+(Dq)_2^2}
=\alpha q_0\sqrt{q_1^2+q_2^2}.
\]
The same fixed matrix must transform the raw nodal map, so the pulled-back force and Hessian are unchanged. The two actual reference edges were independently verified horizontal. This review does not extend this horizontal-reference recipe to arbitrary curved material boundaries.

Write \(q=(a,v)\), \(s=|v|>0\), and \(d=T-K\). Independent explicit derivatives of
\[
W(q)=\tfrac12Ka^2|v|^2+da|v|+\tfrac12K
\]
were checked against exact symbolic rational evaluations and the producer's formal jet. In particular,
\[
H_{aa}=Ks^2,\quad H_{av}=(2Ka+d/s)v,\quad
H_{vv}=Ka^2I+da(I/s-vv^\mathsf{T}/s^3).
\]
The independent third derivatives include
\[
W_{av_iv_j}=2Ka\delta_{ij}+d(\delta_{ij}/s-v_iv_j/s^3)
\]
and
\[
W_{v_iv_jv_k}=da\left[-\frac{\delta_{ij}v_k+\delta_{ik}v_j+\delta_{jk}v_i}{s^3}
+\frac{3v_iv_jv_k}{s^5}\right].
\]
Numeric Torch automatic differentiation at eight additional points agreed to relative error at most \(2.67\times10^{-15}\); this is a numerical diagnostic, separately labeled from the exact/interval checks.

With full six-coordinate current nodal vector \(Y\), first form the rational map difference \(\Delta B=B_*-B_c\), then enclose \(\delta q=\Delta B Y\). The fundamental theorem of calculus gives \(\delta g\in[H]_{\rm line}\delta q\) and \(\delta H_{ij}\in\sum_k[W_{ijk}]_{\rm line}\delta q_k\) on the entire convex line box, provided its hoop and tangent norm are positive. The exact algebra is
\[
\Delta F=\Delta B^\mathsf{T}g_c+B_*^\mathsf{T}\delta g,
\]
\[
\Delta H=\Delta B^\mathsf{T}H_cB_c+B_*^\mathsf{T}H_c\Delta B
+B_*^\mathsf{T}\delta H B_*.
\]
The second Hessian term contains the quadratic \(\Delta B^\mathsf{T}H_c\Delta B\) contribution because \(B_*=B_c+\Delta B\); it is not missing. A fixed-coordinate adversary changes only an excluded nodal column and produces nonzero force/Hessian effects. Those effects were enclosed, verifying that fixed-coordinate contributions are retained before selecting final free columns.

## Safe arithmetic and retained reviewer failure

The direct jet cleanup is `not(v == 0)`, not the invalid `v != 0`. Arb equality with exact zero is false for a nonexact zero-crossing ball, so this predicate retains it. Tests cover zero-crossing and near-cancelled coefficients through products at radii near \(10^{-1}\), \(10^{-80}\), and \(10^{-300}\), plus nonfinite rejection.

The module imports the legacy target, but the inspected local path uses its rational conversion, immutable binary decoder/container, and binding utilities. Both actual site replays also succeeded with legacy sparse `clean` replaced in the review process by a raising sentinel. No old integration arithmetic is approved by that limited dependency result.

The first review-harness run reached actual vertex checks but failed a byte-for-byte producer replay assertion after changing Arb precision from 512 to 256 without clearing the precision-dependent `core.A` conversion cache. A fresh 256-bit process reproduced the original fixture exactly. The failed harness/protocol/log are retained in `precision_cache_failure_revision_1`; the corrected harness clears that cache at precision switches. Both fixtures then reproduced exactly at the producer's 256-bit precision, independently of the 512-bit analytic checks. No producer source or result was modified.

## Remaining work

A future aggregate must multiply every local contribution by its actual old-normalized auxiliary weight exactly once and assemble the correct global indices. The exact-integral-minus-auxiliary part remains separate. These local results do not establish any whole-surface integral/Hessian error or root-persistence inequality and do not restore the withdrawn unsafe point-integration certificate.
