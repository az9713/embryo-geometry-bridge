# Uniform surface bridge: preparation checkpoint

**Preparation only. No actual edge integration or cache sum has run.** The frozen protocol selects all 24 B surface edges, the same physical parameter cell C, and the explicit existing common nodal box. Root membership in that box is not assumed. The box remains an explicit set; no withdrawn unsafe point-integration conclusion is used as a premise.

Source: `../../../scripts/exact_uniform_surface_bridge_B_safe.py`.

Current preparation: 48 passing checks, 161 source/input bindings, 24 whole-edge interval domains and 288 auxiliary-site convex-line domains. Every whole-edge proof covers interval coefficients from all supplied nodal values; it is not a midpoint or point-only positivity check. The maximum expansion ratio is approximately 0.123731, below the declared goal 1/8; spatial subdivision reaches depth five. Decimal values are descriptive; the protocol retains exact rational bounds and every Bernstein coefficient.

## Whole-box domain proof

For each closed chart leaf, expand a power polynomial after the exact affine change `s=lo+(hi-lo)t`, then convert degree-n power coefficients `d_j` to Bernstein coefficients

\[
b_i=\sum_{j=0}^i\frac{\binom{i}{j}}{\binom{n}{j}}d_j.
\]

All interval coefficient endpoints are converted outward to exact rationals before those operations. The convex hull of the interval Bernstein coefficients encloses the polynomial for every spatial point and every allowed coefficient value. Positive denominator division bounds the hoop stretch. Exact interval squaring bounds

\[
G=a^2(b^2+c^2),\qquad
q_0=(G_-+G_+)/2,\qquad
\rho=(G_+-G_-)/(G_++G_-).
\]

Every leaf must collectively yield positive hoop, positive G and rho below one. Axis-edge radial division is removed in the safe field constructor before the closed endpoint is included. Synthetic tests reject positive endpoint samples with an interior collapse and a positive midpoint with a hoop interval crossing zero.

## Frozen integration/cache algorithm awaiting approval

With the safe polynomial jet of G, and the same shared K/T intervals, the chain coefficients are

\[
A=K/2+(T-K)/(2\sqrt G),\qquad B=-(T-K)/(4G^{3/2}).
\]

The run will use order-28 inverse-half and inverse-three-half series about each exact q0 at 640-bit precision. Signed material-integral-minus-auxiliary monomial defects are contracted before taking upper bounds; binary rule moment defects remain nonzero. Uniform tails use both positive material and auxiliary masses, the supremum of `|T-K|`, and rational upper bounds on the interval-coefficient derivative polynomials. This controls the actual nodal/material family even though coefficient interval arithmetic may include additional uncorrelated polynomials.

At every one of the twelve auxiliary sites per edge, the separately reviewed local consistency primitive forms

\[
\delta q=(B_*-B_c)_{\rm full}Y
\]

before interval evaluation, retaining fixed nodal coordinates. Its first/second/third derivatives on the full convex line give the two force terms and three Hessian terms. Each unweighted local bound is multiplied by its actual old-normalized positive auxiliary weight exactly once and assembled through the frozen free-column map. No naive subtraction of broad gradient/Hessian boxes replaces that small-difference identity.

The final decomposition is

\[
I[f]-Q_{\rm cache}[f_{\rm cache}]
=(I-Q_*)f+(Q_*f-Q_{\rm cache}[f_{\rm cache}]).
\]

The first term will receive the signed polynomial moment contraction and two-mass tail. The second receives the full-Y local cache-consistency sum. This remains a surface-only component, not a complete bulk+surface error or a root-persistence test.

## Execution gate and retained versions

`run` requires `prelaunch_independent_review.json` with status exactly PASS, nonempty literal-True checks, and the exact current source, protocol and preflight hashes. Source/input hashes are rechecked before and after execution. Outputs are write-once. The original preparation revision is retained under `prelaunch_revision_0`; no actual integration/cache sums were run in that revision either.

No physical parameters, geometry priors, source observations or frozen kernels were edited. No equilibrium optimization, complete root theorem, continuum validation, novelty claim or five-gate completion is asserted.
