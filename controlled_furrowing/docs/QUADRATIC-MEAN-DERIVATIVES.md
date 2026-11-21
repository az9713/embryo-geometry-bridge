# Exact derivatives of fixed-ROI quadratic regional means

**Completed scope: interval-enclosed node derivatives of the smooth conditional regional means.** The implementation uses saved surface geometry only. It performs no mechanical optimization, material fit, maximum differentiation, implicit equilibrium solve or local-identifiability claim.

The source `scripts/quadratic_mean_derivatives.py` is frozen at `681b69335a80ae0468dd98856c5562b24524827809ece465747b95454d152bff`. It uses the unchanged [quadratic endpoint observer](QUADRATIC-ENDPOINT-OBSERVER.md) only for validated saved-state curve extraction; the derivative formula and polynomial assembly are implemented separately. **35 synthetic checks and 49 checks on one saved 12x6 high-quadrature gain-10 state pass.**

## What is held fixed

Let a continuous piecewise-quadratic boundary have physical coordinates

\[
 R(t)=\sum_{j=0}^2N_j(t)R_j,\qquad
 Z(t)=\sum_{j=0}^2N_j(t)Z_j,\qquad 0\le t\le1,
\]

with the three P2 nodal functions

\[
 N_0=1-3t+2t^2,\qquad N_1=4t-4t^2,\qquad N_2=-t+2t^2.
\]

Assume the current radial derivative is strictly positive on every visible edge section, so current radius has a smooth local inverse. Fix a physical ROI \(I=[c,d]\), component references, and either \(w(R)=1\) or \(w(R)=R\). Its positive normalizer is

\[
 m_w=\int_c^d w(R)dR
 =\begin{cases}d-c,&w=1,\\\frac{d^2-c^2}{2},&w=R.\end{cases}
\]

The boundary mean is \(\overline Z_w=m_w^{-1}\int_c^d w(R)Z(R)dR\). The derivative holds **c, d, the reference heights, the weight rule and the physical length unit fixed**. It is not a derivative at fixed material radius or fixed edge parameter.

## Why radial transport and moving cut endpoints matter

On a visible edge interval \([\tau_-,\tau_+]\), the numerator is

\[
 I_e=\int_{\tau_-}^{\tau_+}w(R(t))Z(t)R'(t)dt.
\]

For a fixed physical cut \(R(\tau)=c\) or d, differentiation gives

\[
 \delta\tau=-\frac{\delta R(\tau)}{R'(\tau)}.
\]

The integrand variation, followed by integration by parts, is

\[
\begin{aligned}
 \delta\{w(R)ZR'\}
 &=w'(R)\delta R\,ZR'+w(R)\delta Z\,R'+w(R)Z\,\delta R',\\
 \int\delta\{w(R)ZR'\}dt
 &=\left[w(R)Z\,\delta R\right]_{\tau_-}^{\tau_+}
   +\int w(R)\left(R'\delta Z-Z'\delta R\right)dt.
\end{aligned}
\]

Leibniz endpoint motion contributes the negative of the displayed endpoint term at each fixed physical cut. At boundaries between full adjacent edges, the remaining endpoint terms cancel because the current trace is continuous and the shared nodal displacement has a single value. Thus the assembled variation is exactly

\[
 \boxed{\quad
 \delta\int_c^d w(R)Z(R)dR
 =\sum_e\int_{\tau_{e,-}}^{\tau_{e,+}}
 w(R)\left(R'\delta Z-Z'\delta R\right)dt.
 \quad}
\]

The \(-Z'\delta R\) term is horizontal transport through a sloped surface. Dropping it incorrectly treats the measured height as following a material point. Holding the parameter cuts fixed also gives the wrong answer when they represent fixed physical ROI boundaries.

## Explicit per-node derivatives

Substituting \(\delta R=N_j\delta R_j\), \(\delta Z=N_j\delta Z_j\) yields

\[
 \frac{\partial\overline Z_w}{\partial R_j}
 =-\frac1{m_w}\int w(R)Z'N_jdt,
 \qquad
 \frac{\partial\overline Z_w}{\partial Z_j}
 =\frac1{m_w}\int w(R)R'N_jdt.
\]

Contributions from edges sharing a node are added. The polynomial integrands have degree at most three for line means and five for disk means. Their primitives have degrees at most four and six. Exact rational coefficients and rational brackets for inverse-cut endpoints therefore enclose each integral and each derivative. Signed coefficients and negative P2 basis values are retained; individual nodal derivatives need not be positive even though the spatial averaging weight is nonnegative.

The code differentiates the physical functional analytically and then encloses those derivative integrals. It does not differentiate a finite number of bisection iterations or treat a branch of interval code as the physical derivative. The exact enclosures refer to the supplied stored piecewise-quadratic curves, not to a continuum PDE or experimental surface.

For the actual conditional component observables,

\[
 \overline d_a=a_0-\overline a_w,\qquad
 \overline d_b=\overline b_w-b_0,\qquad
 \overline s=\overline d_a+\overline d_b,
\]

the apical node derivative has a minus sign, the basal one a plus sign, and the height-mean gradient is their sum.

## Boundary-domain admissibility

The formula defines a directional derivative on continuous trace variations for which the same physical ROI remains covered. If the ROI begins at the physical axis, the axis radial displacement must remain zero in the actual admissible mechanics space. Moving that domain endpoint into the ROI can destroy coverage and need not admit an independent two-sided physical partial derivative.

The full array is a representation of the derivative on **admissible directions**, including formal entries for constrained coordinates. Those constrained entries are not promoted to independent physical partials. The mechanics caller must select the true free displacement coordinates. The current implementation checks strict positive radial derivative on each visible segment; a folded graph or unresolved zero inverse slope is rejected.

## API and units

```python
from quadratic_mean_derivatives import mean_state_gradients, midpoint_node_gradient

result = mean_state_gradients(
    saved_state,
    roi_um=None,                  # declared fixed default [0, patch_radius*33.5]
    length_unit_um="33.5",
    weighting="current_disk",   # or "current_radial_line"
    root_steps=80,
)
gradient_intervals = result["node_gradients"]["height_regional_mean_shortening"]
gradient_centres = midpoint_node_gradient(result, "height_regional_mean_shortening")
```

`node_gradients[key][node][coordinate]` is an exact interval for coordinate R or Z of the **physical** current node. Nonboundary nodes have exactly zero entries. The three keys are:

- `apical_regional_mean_downward`
- `basal_regional_mean_upward`
- `height_regional_mean_shortening`

The mean has units of length, its physical node coordinate has units of length, so these derivative entries are dimensionless. If the mechanical free variable is a dimensionless nodal displacement and \(x_{\rm phys}=L_0(X+u)\), then

\[
 \frac{\partial\overline s}{\partial u_i}
 =L_0\,\frac{\partial\overline s}{\partial x_{{\rm phys},i}}.
\]

Multiply by the physical length unit, flatten in the mechanics node/coordinate order and select the true free coordinates. Float centres are a convenience for a later numerical implicit calculation; they do not replace the retained exact intervals.

If a future parameter changes the reference, length unit, ROI, weight rule or normalization, its explicit chain-rule terms must be added. For example, this API supplies no derivative of the ROI radius with respect to activation radius or H. It also supplies no material derivative through a mechanically solved equilibrium. Those are separate operations.

## Analytical example that exposes both mistakes

Take \(R(t)=t+t^2\), \(Z(t)=t^2\) and the fixed physical ROI

\[
 [c,d]=[5/16,21/16],\qquad[\tau_-,\tau_+]=[1/4,3/4].
\]

The line mass is one and disk mass is 13/16. Exact derivative coefficients are:

| Measure | Node | Derivative in R | Derivative in Z |
|---|---:|---:|---:|
| Line | 0 | 0 | 1/48 |
| Line | 1 | -11/24 | 11/12 |
| Line | 2 | -1/24 | 1/16 |
| Disk | 0 | 17/780 | -17/624 |
| Disk | 1 | -149/312 | 237/260 |
| Disk | 2 | -109/1560 | 361/3120 |

For the line mean, varying the middle radial node has derivative -11/24. Omitting radial transport predicts zero. Differentiating with fixed parameter cuts instead predicts -1/12: it misses a boundary cancellation of 3/8. Independent numerical inverse/quadrature finite differences recover the correct fixed-current-ROI value.

A constant-height graph supplies another check: every radial derivative is exactly zero on a fixed current ROI, despite changes in its radial parameterization. Summing all vertical-node derivatives gives one, as required by rigid vertical translation of the whole boundary. No assumption that every individual P2 coefficient is positive is made.

## Verification receipts and retained scope

`evidence/quadratic_mean_derivatives_synthetic.json` retains 35 checks: the twelve rational coefficients above, twelve separately evaluated finite differences, translation identities, zero radial derivatives for a constant graph, the two missing-term counterexamples and malformed-input guards. Synthetic finite differences use physical step 1e-5 and absolute comparison tolerance 2e-8. These finite differences corroborate the implementation; they are not the proof.

`evidence/quadratic_mean_derivatives_saved.json` retains the complete full-node gradient arrays for line and disk means on **one** saved 12x6 high-quadrature gain-10 field, its hash, and all 49 checks. Six independently generated admissible physical directions per measure test all three means, giving 36 numerical directional finite differences. Physical steps are 2e-5 micrometres and the absolute derivative tolerance is 2e-8. Independent float polynomial inversion and scalar quadrature are used for the perturbed means, rather than the exact derivative primitive. Mean values are also checked against the frozen observer, unused-node gradients are exactly zero, and the source state hash is unchanged.

The three derivative identities, signed assembly and exact rational integral enclosures are the mathematical claims. No derivative of a spatial maximum at ties, unsampled temporal peak, full parameter family, biological mechanism or local-rank/inverse-uniqueness theorem is claimed. The calculus is standard; its role is to provide a checked bridge from the physical observation contract to a later, separately verified implicit material-inference calculation.
