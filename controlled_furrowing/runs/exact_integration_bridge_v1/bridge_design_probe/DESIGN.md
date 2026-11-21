# Local algebra probe and executable integration design

This probe uses the **actual B proof-center field** from the frozen all-17 target-preparation protocol. It does not integrate a complete field, expand order-16 coefficients, solve equilibrium equations, or establish a bridge from the old finite rule to the exact-integral target.

## Measured sparsity

The following counts are exact nonzero sparse-polynomial counts for the sampled fields. The last column is a **Minkowski-support upper bound**, obtained without computing the order-16 coefficients.

| Bulk chart | Local free coordinates | J support | Sum of Ju supports | Sum of upper Juv supports | Order-16 series-times-factor support upper bound |
|---|---:|---:|---:|---:|---:|
| Axis vertex | 11 | 14 | 146 | 500 | 2530 |
| Off axis near the patch | 12 | 15 | 175 | 751 | 2701 |
| Axis edge | 9 | 16 | 108 | 240 | 3025 |
| Off axis at the outer rim | 6 | 11 | 58 | 131 | 1864 |

In physical off-axis coordinates, J has radial exponents -1 through 3 and axial degree at most 4. In axis-vertex coordinates it has s-degree 3, t-degree 4 and total degree 7; the axis-edge example has degree 3 in both variables and total degree 6. These chart-degree increases come from the bilinear Duffy coordinate map, rather than an increase in the degree of the physical P2 displacement.

The union support of Ju_i Ju_j has 45 terms in the interior off-axis sample, 42 at the axis vertex, and 49 at the axis edge. In contrast, independently expanding every constitutive-series-times-Hessian-entry polynomial would repeatedly manipulate thousands of coefficients.

For surfaces use G=Js squared. The first axis edge has five G coefficients with radial powers 0 through 4. The sampled positive-radius edges have seven coefficients with powers -2 through 4. Their order-16 series-times-factor support bounds are only 73 or 109 terms. Both apical and basal samples were inspected.

The complete unchanged finite-dimensional connectivity has 611 free coordinates but only **6676 global upper-triangular structural entries**, or 12741 in the full symmetric pattern. There are 10952 local element/surface upper-triangle slots before assembly. These are connectivity counts; no Hessian values were assembled by the sparsity probe.

## Why the radial coordinate matters

Write the off-axis determinant defect exactly as

\[
J-1=\frac{P(R,\zeta)}{R},\qquad
R=R_c+h_R x.
\]

The helper retains this denominator. It never replaces 1/R with a finite polynomial. It shifts and scales only the polynomial numerator, keeping zeta unchanged. Exact recovery of the original physical-R Laurent polynomial is tested.

For the actual B near-patch sample, a direct physical-R coefficient absolute-sum bound on the enclosing coordinate box is approximately 0.051835. Shifting/scaling the numerator before dividing its absolute-sum bound by Rmin reduces it to approximately 0.006658, a factor of 7.8. In the outer-rim sample the same bounds are approximately 0.018224 and 0.000037653, a factor of 484.

These are valid coefficient bounds for these center polynomials, not attained J ranges, full-tube bounds, or tail acceptance results. A stronger independently established geometric J bound may avoid needing this coefficient enclosure for convergence. Translation is nevertheless valuable for avoiding a badly inflated coefficient-based convergence test.

Back-expanding a degree-4 numerator at order 16 has a conservative radial binomial-cancellation scale of about 128 bits near the patch and 266 bits at the rim. This is a conditioning budget, not a measured error or a lower bound on required precision. A fixed 256-bit evaluation should therefore not be accepted solely because its inputs are exact: the final Arb interval widths must decide whether the chosen precision is sufficient.

## Bulk gradient and Hessian structure

Let

\[
\phi(J)=-\mu\log J+\frac{\kappa}{2}(\log J)^2.
\]

The exact volumetric contributions are

\[
F_i=\phi'(J)J_i,\qquad
H_{ij}=\phi''(J)J_iJ_j+\phi'(J)J_{ij}.
\]

The remaining neo-Hookean and active terms are polynomial quadratic forms in q. Their gradients and Hessians can be handled directly by the same signed moment functional without a nonlinear constitutive series.

For a truncated expansion let A approximate phi-prime and B approximate phi-second, each formed **once per element**. Let Delta-m be the signed exact-material moment minus the moment of the declared auxiliary rule. Define the dual vectors

\[
d_A(k)=\sum_l A_l\,\Delta m_{k+l},\qquad
d_B(k)=\sum_l B_l\,\Delta m_{k+l}.
\]

Then all local entries follow from small contractions:

\[
\Delta F_i^{\rm trunc}=\sum_k (J_i)_k d_A(k),
\]

\[
\Delta H_{ij}^{\rm trunc}
=\sum_k(J_{ij})_k d_A(k)
 +\sum_k(J_iJ_j)_k d_B(k).
\]

Only the union factor supports need dual entries: at most 15 A shifts and 45 B shifts in the sampled interior off-axis case. The local polynomial products Ju_i Ju_j remain small. This avoids separately expanding each of 78 large Hessian-entry polynomials. Geometry/rule moment tables can be memoized and reused across all 17 fields, because the frozen geometry and auxiliary rules are identical.

`polynomial_jet.py` implements this contraction. It returns truncated contributions only. Inputs must be exact rational or outward Arb coefficients/moments. The caller supplies the constitutive expansion, signed moment function, admissibility proof, remainder bounds and finite-cache correction.

## Axis and surface moment conventions

For an axis-edge chart recover s from the off-axis barycentric coordinate, and recover t from the second axis barycentric coordinate divided by 1-s. For an axis-vertex chart use s=1-lambda_axis and t=lambda_second_offaxis/s. Auxiliary points are strictly interior, so both inversions are nonsingular. The axis moment already includes R times the absolute affine determinant and the Duffy Jacobian. Do not multiply those factors a second time.

For surface G=Js squared,

\[
E_s=\frac K2G+(T-K)\sqrt G+\frac K2.
\]

Thus its gradient/Hessian require G_i, G_ij and series or validated analytic bounds for G to the powers -1/2 and -3/2. `surface_squared_jacobian_jet` supplies the polynomial derivatives. Positivity and a suitable expansion range for G must be proved before using such a series; the B-center probe alone does not prove them for every root tube.

## Required remainder accounting

For the decomposition

\[
(I-Q_*)f+(Q_*f-Q_{\rm cache}f_{\rm cache}),
\]

and a polynomial approximation f=Pn+remainder, the first bracket is

\[
\Delta(P_n)+I(\text{remainder})-Q_*(\text{remainder}).
\]

A uniform remainder bound therefore carries **both the exact positive material mass and the sum of the actual positive auxiliary weights**. Counting only the material integral would omit the auxiliary remainder. An equivalent alternative is Delta(Pn)+(Q-star Pn minus Q-cache f-cache) plus the material-integral remainder; in that form the finite correction must use Pn consistently. These alternatives must not be mixed.

Keep signed moments and signed assembled defects until applying the chosen proof norm/preconditioner. Early absolute values discard the cancellation on which a sharp integration-defect bound may depend. A full-tube proof additionally needs uniform parameter/field coefficient bounds, analytic-domain positivity, and the applicable constitutive tail theorem.

## Proposed bounded execution order

1. Freeze the local pilot element, expansion order, moment/callback source hashes, precision ladder and exact auxiliary-coordinate/weight metadata.
2. Test low-degree signed material-minus-auxiliary moments, including deliberately nonzero defects and exact axis chart normalization.
3. Form one shared constitutive polynomial per pilot element and dual-contract its gradient/Hessian defects. Increase Arb precision only through the declared ladder, retaining failures and final widths.
4. Add both remainder contributions and the actual finite auxiliary-minus-cache correction. Report the pilot result without extrapolating to all elements or 17 records.
5. After independent review, assemble the sparse field defect with geometric moment reuse. Only a completed bound may feed a new exact-integral root certificate.

The present task stops before step 3. The helper tests use series polynomials of degree at most two and exact signed point-functionals to verify algebra; those are not physical integration results.

## Verification files

- `results.json`: actual supports/degrees, support-only high-order estimates, center conditioning bounds.
- `derivative_verification.json`: generic second-order polynomial automatic differentiation reproduces all representative derivative coefficient dictionaries.
- `polynomial_jet_tests.json`: 25 exact checks of the reusable helper, radial round trips and equality of dual/direct contractions for every tested local entry.
- `sparsity.json`: exact connectivity-only counts.

The frozen exact target, all old roots, physical parameters and scientific documents remain unchanged. No full field integral, exact-target root, continuum/PDE validation, new biological fit, author contact or novelty claim is supplied by this probe.
