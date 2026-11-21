# Point surface integration pilot: preparation only

The source `../../../scripts/exact_surface_bridge_pilot.py` contains a gated integration engine. **It has not been run on any actual surface field.** The current preparation certifies point-field geometry domains for B edges 0,1,11 on both faces and freezes a 512-bit, order-24 numerical policy. It supplies no full surface sum, uniform parameter/root-tube error, exact-integral root transfer, biological fit or continuum/PDE result.

## Target and measure

All fields use the same exact coordinates Y=Xbinary+free-u as the selected B point cache, retaining fixed coordinates and their small reference-midpoint offsets. The material edge is the canonical horizontal interval Rlo to Rhi from the exact endpoints. The normalized integration functional is

\[
I(f)=\frac{2\pi}{S_{\rm binary}}\int_{R_{\rm lo}}^{R_{\rm hi}} Rf(R)\,dR.
\]

The reviewed `surface_moment` already includes the factor R. It is included once. The auxiliary points are the exact canonical positions Rlo+L*s-star derived from the old cached shape functions. Their weights are the actual old positive binary weights divided by the exact binary energy scale; no Gaussian exactness is assumed.

The finite cache correction evaluates the exact canonical strain/B matrices and the old raw cached N,dN,r0,l0 strain/B matrices at their paired auxiliary points, using the SAME exact Y. The old reference radius and tangent length are explicitly retained. It computes auxiliary-minus-cache in signed Arb arithmetic.

## Positive branch and derivative identity

Let

\[
a=r/R,\qquad b=dr/dR,\qquad c=dz/dR,
\qquad G=a^2(b^2+c^2)=J_s^2.
\]

Preparation encloses the rational chart functions on every closed material edge, subdividing s and centering their exact polynomial coefficients on each leaf. The first axis-edge hoop factor is cancelled symbolically before this evaluation. All leaf hoop lower bounds and the global G lower bound are positive. This matters: G>0 alone would allow a<0, for which sqrt(G) would not equal Js.

On the certified positive branch,

\[
E_s=\frac K2G+(T-K)\sqrt G+\frac K2,
\]

\[
A(G)=\frac K2+\frac{T-K}{2\sqrt G},\qquad
B(G)=-\frac{T-K}{4G^{3/2}},
\]

\[
F_i=A(G)G_i,\qquad H_{ij}=A(G)G_{ij}+B(G)G_iG_j.
\]

The helper's generic moment-dual contraction therefore applies to G and its exact polynomial jets. Constant K/2 is part of A and is not added a second time.

## Binomial tail proof and the two masses

Choose the exact rational center q0=(Glo+Ghi)/2 and set x=G/q0-1. The certified bound is |x|<=rho=(Ghi-Glo)/(Ghi+Glo)<1. This is a numerical expansion coordinate, not a physical parameter choice.

For alpha=-1/2 or -3/2, exact coefficients satisfy

\[
c_0=1,\qquad c_{n+1}=c_n\frac{\alpha-n}{n+1}.
\]

Writing m=N+1, the inverse-half tail is bounded by

\[
\tau_{1/2}=|c_m|\frac{\rho^m}{1-\rho}.
\]

For inverse-threehalf, |c_n|=(2n+1)|binom(-1/2,n)|<=2n+1 gives the universally valid majorant

\[
\tau_{3/2}\le\rho^m\left(
\frac{2m+1}{1-\rho}+\frac{2\rho}{(1-\rho)^2}\right).
\]

When the first omitted-tail ratio bound rho*(m+3/2)/(m+1) is below one, the implementation also uses the sharper first-omitted-term geometric bound and takes the smaller of the two valid bounds.

Multiplying by q0 to the appropriate negative half power and by the absolute constitutive prefactors gives uniform bounds tauA,tauB. For the signed moment defect I-Q-star, force tails are bounded by

\[
(M_I+M_Q)\,\tau_A\,\|G_i\|_\infty,
\]

and Hessian tails by

\[
(M_I+M_Q)\left(
\tau_A\|G_{ij}\|_\infty+
\tau_B\|G_i\|_\infty\|G_j\|_\infty\right).
\]

Here MI is the exact positive material mass and MQ is the sum of the ACTUAL positive auxiliary weights. Both are required. A retained synthetic signed-point example explicitly fails if the auxiliary remainder mass is omitted.

## Preparation evidence

`guard_tests.json` and `preflight.json` record 116 true checks. They cover exact binomial identities through order32, rigorous tail endpoint checks, closed-form synthetic integration, deliberately nonzero quadrature defects, the R weight counted once, full-edge positive-domain fixtures, negative-hoop rejection, local surface force/Hessian identities, matched and mismatched finite-cache toys, and strict review/helper bindings. The toy finite-cache calculations are not actual B surface integrals.

`protocol.json` contains all six full-edge domain certificates and the exact q0/rho values. The nearest-decimal rho values are approximately 0.003556,0.121648,0.123738 on the basal face and 0.002556,0.120936,0.123534 on the apical face. These are conservative point-field expansion bounds; they are not parameter/root-tube bounds or physical uncertainty estimates.

## Frozen interfaces and next action

- `positive_G_bounds` prepares the whole-edge positive branch and expansion domain using exact rational coefficient enclosures.
- `binomial_coefficients` and `binomial_tail` implement the scalar theorem above.
- `SurfaceFunctional` provides signed 1D Laurent material-minus-auxiliary moments.
- `normalized_series`, `cache_correction`, and `run_edge` implement the unexecuted polynomial/tail/finite-cache composition at 512 bits.
- `run` requires a NEW exact matching `prelaunch_independent_review.json` with status PASS, nonempty literal-True checks and matching source, protocol and preflight hashes.

The next action is independent source/protocol/preflight review. Only after approval may the actual six-edge integration pilot run. At that point interval widths, signed defect values and all tails must be retained. Fixed 512-bit precision does not guarantee a useful final width; a wide valid enclosure must not be reported as sharp. The integration engine's moment-cache count is cumulative over live functional instances.

The initial G-only positive-branch omission and the later empty-review/helper-binding guard hardening are retained in engineering_revision_0 and engineering_revision_1. No actual integration occurred under either archived revision. Old physics, caches, proofs and public documents were not edited. No author message was sent.
