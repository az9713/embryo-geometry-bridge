# Independent verification of one algebraic equilibrium certificate

**A separate 256-bit interval calculation confirms a unique local stationary root for one archived 611-degree-of-freedom finite-quadrature problem.** It rebuilds all 9,504 terms and their second derivatives independently. At normalized coefficient radius \(r=10^{-11}\), the independent calculation gives

\[
 \eta\approx2.5277098028\times10^{-12},\qquad
 Z\approx0.037226870681<1,
\]

and verifies the self-mapping inequality with exact rational bounds. The resulting coefficient-distance bound is

\[
 \|u^\dagger-u_0\|_\infty<2.626\times10^{-12}.
\]

This closes an algebraic-root gap for **this one fixed target**. It does not bound quadrature-integration error, finite-element/PDE error, biological error, a parameter family, or the new 17-field material replay.

The [independent implementation](../scripts/review_algebraic_root.py), [final PASS receipt](../evidence/algebraic_root_independent_review.json), [frozen review protocol](../runs/algebraic_root_independent_review/protocol.json) and [full higher-precision output](../runs/algebraic_root_independent_review/higher_precision_replay.json) retain the calculation. The [producer certificate and proof path](ALGEBRAIC-ROOT-CERTIFICATION-PATH.md) are separate artifacts and were never overwritten by this review.

## 1. Which root is actually certified?

The target is the archived small-region strength-one field, with a 12×6 P2 mesh, 611 free displacement coefficients, compression coefficient 5000, fixed material/geometry and Duffy8/surface12 quadrature. It is not one of the newly changed-material cohort states.

Every serialized binary64 cache constant is interpreted as an exact real number: reference coordinates, shape values and gradients, quadrature radii and weights, activation mask, surface constants and the proposed preconditioner. Coordinates are formed as exact \(X+Pu\), where \(P\) inserts the free displacement coefficients. The certificate therefore concerns a precise smooth finite sum with fixed binary constants.

That mathematical potential differs from performing every intermediate operation in rounded Torch arithmetic, and from a finite-element integral with ideal exact Gauss nodes. Those distinctions remain explicit. The review verifies the embedded cache checksum, all current-model cache arrays, state/parameter/strength association, archived free vector and free-index ordering. The 611×611 binary preconditioner is finite and square; its usefulness and invertibility are verified through the bounds below rather than assumed from numerical inversion.

## 2. Independent derivative implementation

The reviewer does not call the producer point assembler or Hessian-variation routine. It implements a separate second-order forward automatic-differentiation object with Arb values, gradients and Hessians.

For products \(w=uv\), it propagates

\[
 Dw=v\,Du+u\,Dv,\qquad
 D^2w=vD^2u+uD^2v+Du\otimes Dv+Dv\otimes Du.
\]

For a scalar composition \(w=\phi(u)\), it uses

\[
 Dw=\phi'(u)Du,\qquad
 D^2w=\phi'(u)D^2u+\phi''(u)Du\otimes Du,
\]

with interval-enclosed logarithm, inverse and square-root operations.

For meridional/hoop strain coordinates \(q=(a,b,c,d,h)\), the energy differentiated is

\[
 J=(ad-bc)h,\qquad
 W(q)=\frac{\mu}{2}
 (a^2+b^2+c^2+d^2+h^2-3-2\log J)
 +\frac{\kappa}{2}(\log J)^2
 +\frac{f\chi}{2}(b^2+d^2).
\]

This directly differentiates \(\log[(ad-bc)h]\), independently of the producer's manually written cofactor/Hessian expressions. The active terms use the vertical reference-fibre components \(b,d\).

For each surface term it differentiates

\[
 J_s=cR\sqrt{x^2+z^2},\qquad
 W_s=\frac K2(J_s-1)^2+TJ_s.
\]

Sixty-three local checks exercise positive/negative off-diagonal strain and tangent components, exact point comparisons, small interval boxes, independent interior-point containment and rejected analytic-domain violations. The checks explicitly include negative off-diagonal Hessian couplings; they do not test only diagonal or identity states.

## 3. Rebuilding the actual finite-dimensional function

The reviewer independently reconstructs every linear strain map from the frozen arrays, including fixed-node contributions and the free-index selection. In each term,

\[
 q_\ell(u)=q_{\ell,0}+B_\ell(u-u_0),\qquad
 F(u)=\nabla E(u)
 =\sum_\ell w_\ell B_\ell^\mathsf T
        \nabla W_\ell(q_\ell(u)),
\]

and

\[
 DF(u)=\sum_\ell w_\ell B_\ell^\mathsf T
       D^2W_\ell(q_\ell(u))B_\ell.
\]

All normalized weights are strictly positive. The full independent assembly contains **9,216 bulk quadrature terms and 288 surface terms**.

As a separate numerical mapping diagnostic, the independent exact-real-cache point values were compared with the original frozen Torch automatic differentiation:

| Diagnostic | Measured difference |
|---|---:|
| Point gradient, absolute maximum | \(2.12\times10^{-12}\) |
| Hessian entries, absolute maximum | \(5.46\times10^{-11}\) |
| Hessian entries, relative maximum | \(2.43\times10^{-15}\) |

The comparison measures mapping/rounding consistency. It is not a claim that the two arithmetic evaluations are identical or that a small floating gradient alone proves a root.

## 4. Full-ball Hessian bound

Let \(B_\infty(u_0,r)\) be the closed coefficient ball. The independent strain enclosure is

\[
 |q_{\ell,a}(u)-q_{\ell,a}(u_0)|
 \le r\,s_{\ell,a},\qquad
 s_{\ell,a}=\sum_j |B_{\ell,aj}|.
\]

Independent interval AD encloses the local Hessian throughout each strain box. If

\[
 \Delta_{\ell,ab}\ge
 \sup_{u\in B_\infty}
 |D^2W_{\ell,ab}(q_\ell(u))
   -D^2W_{\ell,ab}(q_\ell(u_0))|,
\]

then the row sum of the assembled Hessian variation is bounded by

\[
 v_i=\sum_\ell w_\ell\sum_a|B_{\ell,ai}|
                       \sum_b\Delta_{\ell,ab}s_{\ell,b}.
\]

This follows by expanding the assembled Hessian, applying the triangle inequality and summing over the second global index. It retains negative cross-couplings through absolute bounds; no diagonal-only approximation is made. For the fixed binary matrix \(R\),

\[
 \|R(DF(u)-DF(u_0))\|_\infty
 \le\||R|v\|_\infty.
\]

The independent ball calculation checks positive meridional determinant, positive hoop stretch and \(J>10^{-6}\) at every cached bulk quadrature point. It also checks positive surface radius, tangent norm, surface Jacobian and normalized weights. Therefore the logarithms and square roots remain analytic and the model's trial inversion penalty is inactive throughout the certified coefficient ball. This is a finite-quadrature analytic-domain statement, not a whole-body geometric certificate for every nearby deformation.

## 5. Banach proof and why the preconditioner suffices

Define

\[
 G(u)=u-RF(u),\quad
 \eta\ge\|RF(u_0)\|_\infty,\quad
 Z_0\ge\|I-RDF(u_0)\|_\infty.
\]

Combining the point defect and full-ball variation gives

\[
 \sup_{u\in B_\infty}\|DG(u)\|_\infty
 \le Z_0+\||R|v\|_\infty=:Z.
\]

The independent point defect is approximately \(6.59123\times10^{-7}\). The total independent bound is approximately 0.03722687. All acceptance inequalities use outward interval bounds or their exact rational endpoints, not these decimal displays.

For any \(u\) in the ball,

\[
 \|G(u)-u_0\|_\infty
 \le\eta+Zr<r,\qquad Z<1.
\]

Thus \(G\) is a contraction mapping the complete closed ball into itself. Banach's theorem gives one fixed point \(u^\dagger\) in that ball.

The fixed point satisfies \(RF(u^\dagger)=0\). Moreover \(Z_0<1\) implies \(RDF(u_0)\) is invertible by the Neumann-series criterion; in a square finite-dimensional problem this forces \(R\) to be invertible. Consequently \(F(u^\dagger)=0\), and any other root in the ball would be another fixed point, contradicting uniqueness. No energetic-minimum or stability assertion follows.

Finally,

\[
 \|u^\dagger-u_0\|_\infty
 \le\frac{\eta}{1-Z}
 \approx2.62544697794\times10^{-12}.
\]

The independent interval-AD bound is slightly looser than the producer's manual-derivative bound, but independently proves the same radius-\(10^{-11}\) ball contains one local root. Multiplication by the fixed 33.5 µm length unit gives a coefficient-coordinate distance below \(8.80\times10^{-11}\) µm. That is only a unit conversion, not a continuum or experimental accuracy estimate.

## 6. Receipts and limits

The final independent receipt passes **20 actual-target checks** in addition to the 63 local derivative checks. It verifies the exact rational producer inequalities and distance bound, repeats the actual full finite-sum proof at higher precision using separate derivative/assembly/variation code, checks source/cache/state identities, and confirms the producer certificate stayed unchanged.

Producer certificate SHA256:

    132bd7695619b1ff67b276f7b65c547d980a74ab76df9c736e221e888d186721

Producer source SHA256:

    7321d38bf17cf354cc03be2871f3fbaa45bcebbcbfca7be90d0f867b80078acf

Independent review source SHA256:

    6af2311d0667a28af263b6fb6237063481aa8e65209ae9722c2b0246a8cff27a

Independent PASS receipt SHA256:

    f327f2fab76f9fdcab3a83882adbb1a4ae3f0cf1e5c269e7335fb909691085d5

A separately justified observer propagation can turn this coefficient-root enclosure into endpoint enclosures. It must account for the distinction between exact \(X+Pu_0\) and rounded saved current nodes, graph/ROI admissibility, and the observation functional. The coefficient bound by itself does not provide those results.

The remaining scientific gaps are separate: parameter-uniform bounds, the other calibration states, source-matched temporal/ROI observations, quadrature/PDE accuracy, full material/geometry compatibility and frozen larger-region predictions. This certificate establishes neither a new biological mechanism nor confirmed novelty.
