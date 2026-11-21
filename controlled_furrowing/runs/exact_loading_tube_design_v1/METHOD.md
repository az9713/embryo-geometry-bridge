# Certifying a loading branch between exact finite-element roots

The archived nonlinear run has 44 stored checkpoints. Each material has one 18-state small-geometry schedule from zero through 17 doses and a separate four-state large-geometry schedule from zero through one-third, two-thirds and the full large dose. The saved continuation histories record numerical warm starts and rounded solver inputs. The independently certified endpoint minima apply at their exact requested loads. A branch between them needs a uniform certificate and a proof that its endpoints are the certified roots. [Frozen probe](../transfer_ambiguity_nonlinear_probe_v1/preparation/protocol.json) · [All checkpoint identities](../transfer_ambiguity_nonlinear_probe_v1/field_manifest.json) · [Nonlinear report](../transfer_ambiguity_nonlinear_probe_v1/report/REPORT.md)

This design concerns the exactly integrated, normalized, finite-element energy at fixed candidate material, fixed reference mesh, fixed support and fixed activation mask. Each of the four material–geometry combinations needs its own chain. The previous gain-family theorem covers a different near-unit weak-load cached-quadrature benchmark. Its proof architecture is reusable. The old weak-load bounds and the new singleton endpoint bounds have no uniform authority for the intervening dose slabs. [Existing gain-family theorem](../../docs/ENERGY-GAIN-CELL.md) · [Its geometric composition](../../docs/ENERGY-GAIN-GEOMETRY.md) · [Identified-set branch contract](../../docs/VALIDATED-IDENTIFIED-SET.md)

## The dose terms are simple on a fixed reference

Use \(I\) for the exact optical dose and \(u\in\mathbb R^{611}\) for the free displacement. Write \(F(u,I)=\partial_u E(u,I)\) and \(H=\partial_uF\). The active energy is \(gI/2\) times the integral of \(\chi\lVert F_{\mathrm{def}}e_Z\rVert^2\) on the fixed active reference elements, with the same normalization as the passive terms. On an element, \(F_{\mathrm{def}}e_Z=a+B_zu\) is affine in \(u\), including the exact stored reference offset \(a\). Let \(A\) and \(b\) be the assembled integrals of \(B_z^TB_z\) and \(B_z^Ta\) with that reference measure and normalization. For a fixed gain \(g\),

\[
C=gA,\quad d=gb,\qquad
\partial_I F(u,I)=Cu+d,\quad
\partial_I H=C,\quad \partial_I^2F=0.
\]

Thus \(F(u,I_c+t)=F(u,I_c)+t(Cu+d)\) and \(H(u,I_c+t)=H(u,I_c)+tC\) exactly. The active term has no force contribution from the surface terms. The active meridional integrands have degree at most three including the reference radial weight on the affine triangle. Exact rational simplex moments suffice before multiplying by the common rigorous \(2\pi/S_{\rm binary}\) factor. The cached operator uses its own frozen binary maps and normalized weights. Its active coefficients may differ. One can assemble the exact spatial coefficients directly, or include the signed active-coefficient difference once with the cache coefficients. The identity assumes a fixed mask in reference coordinates and a fixed dose-to-strength law \(gI\). [Exact target](../../scripts/exact_fe_target_safe.py) · [Active implementation](../../scripts/axisymmetric_mixed.py#L100) · [Reviewed force derivative wrt gain](../../scripts/exact_centered_bulk_L.py#L133)

## A closed-slab root and minimum test

Let \(J=[I_-,I_+]\) contain \(I_c\), set \(t=I-I_c\), and take \(h=\max_{I\in J}|t|\). Choose a continuous affine predictor \(p(I)=u_c+t v\), a fixed binary square preconditioner \(R\), and a radius \(r>0\). The fiber is \(X_I=p(I)+[-r,r]^{611}\). Certify a convex nodal envelope \(Y\) containing every fiber and every segment from \(u_c\) to \(p(I)\). A simple centered envelope gives free coordinate \(i\) radius \(r+h|v_i|\) and every constrained coordinate radius zero.

On that whole envelope, prove admissible finite-element geometry, positive body Jacobian and trace domains for the exact energy and its derivatives. Then obtain outward bounds

\[
Z_c\ge \sup_{u\in Y}\lVert \mathrm{Id}-RH(u,I_c)\rVert_\infty,
\qquad M\ge\lVert RC\rVert_\infty,
\qquad Z=Z_c+hM.
\]

The old cached-plus-defect decomposition can supply \(Z_c\) at the single central dose, using new whole-envelope bounds. The exact dose identity supplies the entire interval contribution \(hM\). Require a signed predictor residual enclosure with \(\eta\ge\sup_{I\in J}\lVert RF(p(I),I)\rVert_\infty\). A directly executable first test uses \(v=0\) and takes \(\eta=\eta_c+hB\), with \(\eta_c\ge\lVert RF(u_c,I_c)\rVert_\infty\), \(B\ge\lVert R(Cu_c+d)\rVert_\infty\). The sufficient step inequality is

\[
\eta_c+Z_c r+h(B+rM)<r.
\]

If \(B+rM>0\) and \((1-Z_c)r-\eta_c>0\), the ratio \(((1-Z_c)r-\eta_c)/(B+rM)\) is a strict upper limit for this sufficient half-width. If the denominator is zero, this bound imposes no additional dose step limit. The other domain guards still apply.

For an affine predictor, put \(q=RF(u_c,I_c)\) and enclose the signed slopes \(s=R[H(u_c,I_c)v+Cu_c+d]\) and \(k=RCv\). Along the whole predictor segment the exact identity is

\[
RF(p(I),I)=q+ts+t^2 k+
t\int_0^1 R\bigl[H(u_c+\lambda t v,I_c)-H(u_c,I_c)\bigr]v\,d\lambda.
\]

Use an interval difference of the central-dose Hessian bound and the center Hessian to bound the directional remainder, preserving signs before taking its row maximum. This requires no new mixed-parameter AD and no derivative of order three. A Lipschitz bound \(L_H\) for \(u\mapsto RH(u,I_c)\) gives the sharper remainder \(L_H t^2\lVert v\rVert_\infty^2/2\) if proved. The active \(t^2k\) term is explicit. Moving the state first at \(I_c\) yields this identity. An interval sum with the same \(t\) is already safe. Directed scalar polynomial range bounds can keep more cancellation.

If \(Z<1\) and \(\eta+Zr<r\), the map \(u\mapsto u-RF(u,I)\) contracts every fiber into its interior and has a unique fixed point there. Since \(RH\) is invertible, the fixed point is a force root with \(\lVert u_*(I)-p(I)\rVert_\infty\le e=\eta/(1-Z)<r\). The translated map on the fixed corrector cube is continuous in \(I\) and uniformly contractive. The fixed point is therefore continuous across the slab. The implicit-function theorem upgrades it to a smooth local branch as long as the energy remains smooth on the admitted domain. Its first derivative obeys \(u_*'(I)=-H^{-1}(Cu_*+d)\), with a bound from \(\sup_Y\lVert R(Cu+d)\rVert_\infty/(1-Z)\) when needed.

Strict minima also hold throughout the slab if the very same binary \(R\) has certified positive symmetric part. Indeed \(RH\) has spectrum in the right half-plane and \((RH)H^{-1}+H^{-1}(RH)^T=R+R^T\). The Lyapunov integral is positive definite, so the symmetric \(H\) is positive definite. The uniform normalized curvature floor is \((1-Z)/\lVert R\rVert_\infty\). A previously verified symmetric-part proof applies only if the exact matrix bytes and their normalization are unchanged. [Existing stability proof](../../docs/COHORT-ENERGY-STABILITY.md) · [Exact signed residual and bridge](../../docs/L-CENTERED-ROOT.md)

## Glue actual roots rather than geometric boxes

Adjacent closed slabs must cover the exact dose interval in order with no gap. At a common dose \(I_*\), they use the same exact force operator, and one certified root enclosure must lie in the other's uniqueness cube. A sufficient directional check is

\[
\lVert p_j(I_*)-p_{j+1}(I_*)\rVert_\infty+e_j(I_*)<r_{j+1}.
\]

If the other direction works, use it. A separately certified small junction root enclosure contained in both cubes also works. A further option is uniform injectivity of the force map on a convex bridge containing both root enclosures at \(I_*\). Overlap of the two cube domains alone does not identify their roots. The same enclosure-to-uniqueness check attaches original endpoint certificates. When an archived endpoint proof uses a recentered binary vector, use that proof center and its certified error for gluing, and retain its original shift and storage gap for observer transfer.

Intermediate radii are chosen by geometric and contraction headroom and junction inclusion. The old small endpoint observer goals need apply only at their scored endpoints. A later continuous-dose peak observer would need an enclosure on each complete slab plus any extremum switches and declared observer/reference/window semantics. These are separate from proving a continuous path of stable stationary states. [Nodal observer transfer](../../docs/NODAL-OBSERVER-ERROR-TRANSFER.md) · [Quasistatic peak scope](../../docs/JOINT-OBSERVATION-CONTRACT.md)

## Bounded first pilot and reusable premises

Use the negative-direction small-geometry adjacent archived endpoints K and L, respectively `small_16` and `small_17`. Their exact optical doses are \(128898317/10^8\) and \(2360271/1562500\), a gap of \(0.22159027\). Both already have independently verified exact-spatial local minima, under one material, the same small reference and active element support, and the same normalization. Start with K’s exact recentered proof center and L’s actual proof center. Propose a binary secant predictor and a rational slab deck. An actual first trial can test one slab joining the two anchors, then bisect the dose interval if the envelope test or gluing fails. If a junction needs a new numerical seed, freeze it as a proposal and independently verify it. A failed sufficient bound remains unresolved. [K exact-spatial anchor](../transfer_ambiguity_nonlinear_probe_v1/certification_design/exact_targets_v1/remaining32_bridge_v1/first_wave_reduction_v1/independent_actual_preparation/actual_replay/receipt.json) · [L exact-spatial anchor](../transfer_ambiguity_nonlinear_probe_v1/certification_design/exact_targets_v1/bridge_reduction_v1/independent_actual_preparation/all_four_readback.json)

The anchor crosswalk compares the exact material singleton, reference topology, free-coordinate map, both surface routes, reference mask and cached integration arrays. The activation axis \(e_Z\) and the pointwise logarithmic-volume energy match the same exact target implementation. Both loads are the same exact gain times their respective nominal doses. The common positive normalization matches exactly here. In general, a different positive normalization of the same physical energy preserves its zeros and strict minima. Any reused preconditioned proof bounds would need the corresponding normalization of the preconditioner. [Anchor crosswalk and source hashes](receipt.json)

| Needed bound | Existing reviewed building block and fresh obligation |
|:--|:--|
| Whole-envelope body embedding, trace domains, and cache analytic guard | Use `root_body_geometry.certify_box`, fixed axis/free map and existing surface/trace domain methods on the larger envelope. Prove every relevant slope/domain on the supplied box. Endpoint checks alone cannot be reused across it. |
| Signed center force, uniform central-dose derivative, and optional directional residual | Use cached point force and fresh variation bounds plus the reviewed same-geometry exact spatial defect reduction at \(I_c\). If a cached variation routine takes an isotropic radius, use the containing radius \(r+h\lVert v\rVert_\infty\) and certify its whole domain, or independently verify per-coordinate bounds. The singleton driver and its receipts remain fixed. The B/L signed predictor method illustrates the optional signed slope before taking norms. |
| \(R(Cu_c+d)\), \(RC\), and if used \(R(H(u_c,I_c)v+Cu_c+d)\) | Assemble the fixed mask’s axial moments and the prescribed fixed-coordinate offset at the same normalization. The existing centered active gain derivative verifies the moment/map pattern for the defect; it is not a ready-made numerical coefficient for this K–L slab. Recheck sparse projections and exact/cache identities. For the constant predictor no central point-Hessian action is required. |
| Minimum and overlap | Preserve exact \(R\) where possible, with its reviewed positive symmetric part and exact infinity norm. Certify new \(R\) if changed. Compute overlap and endpoint inclusion using exact nominal doses, centers and directed radii. |

The pilot would establish one finite-dimensional branch segment if every slab and junction passes. Extending it to zero and the full schedules requires fresh zero anchors and coverage of all intervals. It would still leave experimental time, reference growth, dose calibration, branch selection from dynamics, continuum accuracy and full physical-parameter inference open. The source’s own published maxima also challenge one common nested-peak model at literal zero discrepancy under fixed geometry, a separate already-audited constraint. No path was computed or certified here and the continuation tools provide no novelty claim. [Source limits and conditional peak incompatibility](../../docs/PHYSICAL-INFERENCE-READINESS.md)
