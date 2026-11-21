# Uniform local surface cache consistency

`../../../scripts/exact_surface_cache_consistency.py` supplies an unweighted local force/Hessian difference enclosure over the SAME full six-coordinate nodal box and shared K/T parameter box. Preparation includes synthetic identities and two actual local sites over the NEW common B geometry boxes. It does not sum all24 surface edges or establish an integral/root bound.

## Exact normalization of the old cache

The old cache uses raw coordinates qraw=(r,r-prime,z-prime), old positive radius r0c and length L0c, and alpha-c=1/(r0c*L0c). Let LI be the EXACT canonical horizontal endpoint distance. Set

\[
D_c=\operatorname{diag}(\alpha_c L_I,1/L_I,1/L_I),\qquad
q_{cn}=D_cq_{\rm raw}.
\]

Because LI>0,

\[
(q_{cn})_0\sqrt{(q_{cn})_1^2+(q_{cn})_2^2}
=\alpha_c r\sqrt{(r')^2+(z')^2}=J_{s,\rm old}.
\]

Dc is a fixed exact geometry matrix, independent of nodal Y and K/T. Consequently replacing Braw by Bcn=Dc*Braw gives exactly the old nodal gradient/Hessian pullbacks. Normalization does not replace old cached radius/length values by canonical ones; their differences remain explicitly in Dc.

The canonical map Bstar uses exact P2 shape functions at the derived interior Qstar point, radius Re=Rlo+LI*s-star and derivatives divided by LI. Both full3-by-6 maps act on the SAME Y, including fixed coordinates. The proof uses

\[
\Delta B=B_\star-B_{cn},\qquad \delta q=\Delta B_{\rm full}Y
\]

with exact rational coefficient subtraction before interval evaluation. Full fixed-coordinate contributions are retained; only the final projection selects free columns. The axis-first-edge point is interior, so its positive Re causes no0/0 operation.

## Local energy and formal third-order jet

For normalized q=(a,b,c),

\[
W(q)=\frac K2a^2(b^2+c^2)+(T-K)a\sqrt{b^2+c^2}+\frac K2.
\]

The entire convex q-line box must certify a>0 and b squared plus c squared>0. The implementation uses a three-variable Arb polynomial modulo total degree4: exactly20 monomials of degree at most3. The square-root jet is

\[
\sqrt{C+N}=\sqrt C\left(1+\frac N{2C}-\frac{N^2}{8C^2}+\frac{N^3}{16C^3}\right)
\pmod{\text{degree }4}.
\]

C is a positive constant coefficient and N has its constant coefficient REMOVED exactly. N is formally nilpotent. No numerical smallness or spatial Taylor convergence assumption on N/C is needed to extract derivatives through order3. Interval coefficients conservatively cover all base q and shared K/T values in their boxes.

For each actual fixture, physical K>0,T>=0 is checked on the exact declared cell before conversion to outward Arb balls. The analytic polynomial/square-root potential remains valid on slightly enlarged finite K/T enclosures, so an outward ball touching a small negative value at the T=0 boundary is not silently clipped.

## Signed small-difference proof

With cache local gradient gc and Hessian Hc,

\[
\delta g=\int_0^1 H(q_c+t\delta q)\delta q\,dt,
\qquad
\delta H_{ij}=\int_0^1\sum_kW_{ijk}(q_c+t\delta q)\delta q_k\,dt.
\]

Uniform local derivatives on the entire convex line box give the integral enclosures. Then

\[
\Delta F=\Delta B^Tg_c+B_\star^T\delta g,
\]

\[
\Delta H_{\rm nodal}
=\Delta B^TH_cB_{cn}+B_\star^TH_c\Delta B
+B_\star^T\delta HB_\star.
\]

Both models use the SAME physical K/T. No K/T difference term is introduced. Signed interval terms are retained, avoiding subtraction of wide gradient/Hessian boxes.

Outputs are UNWEIGHTED local point contributions. A future whole-cache runner must apply each actual old-normalized auxiliary weight exactly once and assemble the correct global free indices. This component does not bound the exact-material-integral-minus-Qstar part.

## API

- `surface_derivatives(q,K,T)` returns value, gradient, Hessian, full third tensor, formal jet and analytic-domain bounds.
- `canonical_maps` constructs the exact normalization and paired maps from declared reference/cache data.
- `actual_maps(target,face,edge,iq)` derives the point from the frozen old shape coefficients and returns canonical/normalized-cache full maps, free indices and old normalized weight. All indices reject booleans.
- `consistency_bound(Bstar_full,Bcache_normalized_full,SAME_Y,selected_local,K,T)` returns force/H differences, each signed term, coefficientwise delta-q and convex-line domain.

## Retained verification and scope

`synthetic_tests.json` records409 true checks:19 general Sympy formal-jet derivative identities;120 exact point value/gradient/Hessian/third checks; square-root/nilpotent/zero-overlap guards; analytic/index guards; normalized/raw pullback equivalence; matched-map zero difference;256 mismatched-map interval nodal/K/T endpoint witnesses; and a small-difference diagnostic. The interval fixture bounds delta-q near1.71e-13, whereas direct interval subtraction gives about2e-5. Endpoint witnesses supplement the uniform interval mechanism, rather than replacing it.

The protocol was frozen BEFORE the two actual local fixtures: basal face edge0/q0 and apical face edge11/q11. They use the NEW old-centered common geometry boxes, unchanged C and shared surface coefficients. The largest individual unweighted force/H signed-term bounds are approximately3.83e-14/4.05e-12 at the first site and2.88e-16/6.67e-15 at the outer site. Exact rational endpoints are retained in the JSON packets. These are local bounds, not global assembled norms.

`engineering_revision_0` retains a test-conversion issue: reconstructing an exact rational as an Arb ball can produce a wider test input than the derivative enclosure. Exact Sympy targets are now checked directly against exported rational endpoints.

`engineering_revision_1` retains a substantive interval cleanup bug caught by the mismatched-box fixture: Arb `v!=0` means provably unequal and can be false for a nonzero-width ball overlapping0. The corrected jet cleanup drops ONLY exact zero (`not(v==0)` retains all uncertain coefficients). A dedicated guard prevents recurrence. The analogous pattern in the frozen target polynomial helper was reported to the parent for a separate versioned audit; this module does not modify or assess the old integration receipts and does not use that helper's sparse operations for interval jets.

Independent review is pending. No whole-surface job, full geometry/H assembly, integration, Newton solve, physical parameter change, source-observation validation, new root, PDE bound, publication edit or author contact is supplied by this preparation.
