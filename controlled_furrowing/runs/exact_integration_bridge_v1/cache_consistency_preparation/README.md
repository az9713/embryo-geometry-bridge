# Bulk cache consistency: local uniform primitive

`../../../scripts/exact_cache_consistency.py` supplies a bulk-only local Hessian consistency bound over a supplied nodal/parameter box. Preparation includes synthetic proofs/tests and **two actual B auxiliary-site fixtures**. No full144-element/611-coordinate cache job, quadrature sum, integration, surface bound, solve or exact-integral root transfer has been run.

## Preserving the small difference

Bstar and Bc map the SAME exact current nodal vector Y to the canonical and cached five-component strains. Set

\[
\Delta B=B_\star-B_c,\qquad \delta q=\Delta B_{\rm full}Y.
\]

Full maps include all12 nodal coordinates, including fixed coordinates and retained binary offsets. Only the final Hessian projection selects local free columns. Using free displacements alone for delta-q would omit reference/fixed-coordinate contributions.

Exact coefficient subtraction occurs before interval multiplication. Two synthetic fixtures have nodal halfwidth0.001 and drive in[1,2]. Their coefficientwise delta-q upper bounds are about2e-12, whereas direct interval subtraction has an upper bound about0.002. These are interval-conditioning diagnostics, not physical model differences.

## Signed decomposition and proof

For shared physical parameters on both sides,

\[
B_\star^TH_\star B_\star-B_c^TH_cB_c
=\Delta B^TH_cB_c+B_\star^TH_c\Delta B
+B_\star^T\delta H B_\star.
\]

The production bound never subtracts wide constitutive-Hessian interval matrices. Instead,

\[
\delta H_{\alpha\beta}=\int_0^1\sum_\gamma
W_{\alpha\beta\gamma}(q_c+t\delta q)\delta q_\gamma\,dt.
\]

The convex-line enclosure is the cache strain box plus the interval hull of zero and coefficientwise delta-q. Their common-variable dependence is conservatively overenclosed; independence is not assumed. The third tensor is bounded on the entire resulting box.

Write J=(ad-bc)h and

\[
\phi(J)=-\mu\log J+\frac\kappa2(\log J)^2.
\]

Quadratic neo-Hookean and active-drive terms have zero third strain derivative. Hence

\[
W_{\alpha\beta\gamma}
=\phi^{(3)}J_\alpha J_\beta J_\gamma
+\phi^{(2)}(J_{\alpha\beta}J_\gamma+J_{\alpha\gamma}J_\beta+J_{\beta\gamma}J_\alpha)
+\phi^{(1)}J_{\alpha\beta\gamma},
\]

\[
\phi^{(1)}=\frac{\kappa\log J-\mu}{J},\quad
\phi^{(2)}=\frac{\kappa+\mu-\kappa\log J}{J^2},\quad
\phi^{(3)}=\frac{2\kappa\log J-3\kappa-2\mu}{J^3}.
\]

The only nonzero third coefficients of cubic J are permutations of(a,d,h) with sign+1 and(b,c,h) with sign-1. The entire convex strain box must certify detFm>0,h>0,J>1e-6, keeping the numerical trial penalty inactive. Failure of this sufficient test is unresolved, rather than evidence of physical invalidity.

## API

- `jacobian_tensors(q)` returns J and its first three strain derivatives with the analytic guard.
- `bulk_third_tensor(q,mu,kappa)` returns the symmetric constitutive third tensor and scalar/domain intervals.
- `exact_bulk_maps(target,element,quadrature_index)` constructs full5-by-12 rational maps at one derived Qstar point and supplies local/free associations and the old normalized weight.
- `consistency_bound(Bstar_full,Bcache_full,SAME_Y,selected_local,mu,kappa,drive)` returns the three signed terms, their sum, coefficientwise delta-q, local delta-H and the certified convex-line domain.

B coefficients are exact rational/binary values; Y and shared parameters may be interval enclosures. Outputs are **unweighted local free-coordinate contributions**. A caller must apply each exact old normalized weight once and assemble every required sparse global entry. A whole-field bound also needs all auxiliary sites and a surface counterpart.

## Verification

`synthetic_tests.json` retains329 true checks, including35 exact Sympy identities for all distinct third-derivative index multisets of fullW,155 J-tensor checks, three scalar-phi derivative identities, permutation-envelope consistency, analytic guards, zero-map difference,128 endpoint witnesses in two interval nodal/drive fixtures, and narrow-delta-q checks. Endpoint tests supplement the uniform interval proof mechanism.

The protocol was frozen before two actual local fixtures. `local_element0_q0.json` and `local_element6_q63.json` cover the selected B old-root nodal cube and common parameter cell. Both sites are inactive; parameter dependence still enters their nodal predictor tube. Synthetic fixtures separately exercise an interval active coefficient.

Local constitutive delta-H bounds are approximately8.95e-13 and2.26e-12. Individual projected signed-term bounds are approximately2.1e-9 through4.6e-9 **before weights**. These are conservative local bounds, not an assembled norm or root-transfer certificate; exact endpoints are in the packets.

`engineering_revision_0` retains a test that incorrectly used Arb provable-point equality for interval envelopes, corrected to exact lower/upper endpoint comparison. `engineering_revision_1` retains a synthetic[0,2] Arb hull crossing zero outward and conservatively rejected by the strict drive guard. The fixture now uses[1,2]; no physical guard was weakened. Actual inactive sites use exact zero.

Independent source/protocol review is next. Sources and local receipts are frozen before any full-cache job. Old physics, roots, public documents and author-contact state are unchanged. Surface consistency, a uniform assembled bound and exact-integral root transfer remain outstanding.
