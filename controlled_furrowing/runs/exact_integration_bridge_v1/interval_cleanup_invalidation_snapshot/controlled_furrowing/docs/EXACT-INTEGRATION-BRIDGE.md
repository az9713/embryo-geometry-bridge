# From cached quadrature to an exactly integrated mechanical energy

**The original five-gate study remains incomplete.** The complete bulk-plus-surface integration discrepancy is now bounded at one saved B-center field. It is large enough to prevent direct transfer of the old self-map certificate with its unchanged cube and preconditioner. This is not an equilibrium-nonexistence theorem. A uniform equilibrium-transfer theorem has not yet been verified. The model remains an idealized finite-height axisymmetric finite-element model, with no confirmed new biological mechanism or novelty.

## Why this step is necessary

The previous seventeen stationary-state and local-minimum certificates concern the precise finite-dimensional energy defined by its frozen binary quadrature cache. They do not automatically certify the energy obtained by integrating the same constitutive law over the actual material elements. This distinction matters especially when the self-map margins in a root proof are extremely small.

There are three separate targets: the cached-quadrature finite-element model, the exactly integrated finite-element model, and the continuum boundary-value problem. This chapter addresses the first transition. It does not treat finite-element interpolation error or prove continuum stability.

B labels the selected column-B small-region calibration record. Its saved displacement vector is a computational proof center. Lengths use a fixed reference-height normalization: this numerical tissue cylinder has height 1, outer radius 3 and an activated footprint of radius approximately 0.279257, with activation also restricted in depth. These are declared idealized geometry assumptions; matched per-specimen initial geometry and uncertainty remain unavailable.

The previously certified parameter cell is an extremely small computational box around one candidate. Its cached-model static outputs have a worst-case residual near 5.19 micrometres. It does not represent the full physically compatible parameter set or establish measurement uncertainty. The present chapter checks numerical target consistency for that candidate; no improved fit is claimed.

## 1. Define the exact target without silently changing the stored field

Each material triangle is the exact affine triangle determined by its three binary corner coordinates. Its P2 shape functions are the mathematical Lagrange polynomials. The current nodal coefficients remain the exact values of stored reference coordinates plus injected free displacement, including fixed coordinates. Stored midpoint roundings are retained; the stored reference interpolant therefore need not be precisely the affine identity. There are 611 free displacement coordinates.

With meridional coordinates \((R,Z)\), current coordinates \((r,z)\), and

\[
F_m=\begin{pmatrix}a&b\\c&d\end{pmatrix}
=\begin{pmatrix}r_R&r_Z\\z_R&z_Z\end{pmatrix},
\qquad h=r/R,\qquad J=(ad-bc)h,
\]

the bulk density is

\[
W=\frac\mu2(a^2+b^2+c^2+d^2+h^2-3)
-\mu\log J+\frac\kappa2(\log J)^2
+\frac{\chi gI}{2}(b^2+d^2).
\]

Here \(\chi\) is the frozen elementwise activation mask, \(I\) is the source dose and \(g\) the common gain. The current calculations use \(\mu=1\), \(\kappa=5000\), and the exact common-gain-times-dose value rather than the archived solver's rounded activation strength.

For a horizontal reference boundary, define

\[
J_s=\frac rR\sqrt{r_R^2+z_R^2},\qquad
\psi(J_s)=\frac K2(J_s-1)^2+T J_s.
\]

The exact energy is

\[
E_I=\frac{2\pi}{S_{\rm binary}}
\left[\sum_{\mathcal T}\int_{\mathcal T}R W\,dR\,dZ
+\sum_{\Gamma}\int_{\Gamma}R\psi\,dR\right].
\]

The positive normalization is the same frozen binary energy scale as before; the numerator uses mathematical pi. The logarithm requires positive \(J\); the surface formulas additionally require positive hoop ratio and nonzero tangent. These are proof obligations, not numerical clipping conventions.

**Receipt:** [target preparation](../runs/exact_integration_bridge_v1/target_preparation/preflight.json), 56 checks and 147 input bindings; [independent geometry alignment](../runs/exact_integration_bridge_v1/independent_design/geometry_alignment.json), all seventeen canonical fields enclosed by the existing body boxes.

## 2. Remove apparent singularities at the axis

The mesh has six triangles with an entire axis edge, six with only an axis vertex, and 132 positive-radius triangles. On an axis-vertex triangle use the chart

\[
\lambda_0=1-s,\qquad \lambda_1=s(1-t),\qquad\lambda_2=st.
\]

The equal non-axis radii in this structured mesh give \(R=sR_1\), while the current radial field also contains a factor \(s\). Cancel it symbolically before evaluation. The area Jacobian contributes another \(s\). An axis-edge triangle instead has chart Jacobian \(1-s\). The resulting weighted monomial integrals are respectively

\[
\int R\,s^a t^b\,|\det D_M|\,s\,ds\,dt
=\frac{R_1|\det D_M|}{(a+3)(b+1)},
\]

\[
\int R\,s^a t^b\,|\det D_M|(1-s)\,ds\,dt
=\frac{R_1|\det D_M|}{(a+2)(a+3)(b+1)}.
\]

A vertex-only axis field may have a direction-dependent corner limit in physical coordinates. The chart handles that dependence; it does not assume a unique limit. The first surface edge is regularized by the same exact cancellation principle.

**Receipt:** [moment component independent review](../evidence/exact_fe_moments_independent_review.json), including 48 axis-measure checks and independent simplex and Laurent comparisons.

## 3. Separate integration error from cache inconsistency

Define an auxiliary rule \(Q_*\) using exact rational canonical points reconstructed from the stored shape rows, with the actual old positive binary weights. Its points are not claimed to be exact Gaussian abscissae, and its moments are not assumed exact.

For any gradient or Hessian component,

\[
I[f]-Q_{\rm cache}[f_{\rm cache}]
=(I-Q_*)[f]+\bigl(Q_*[f]-Q_{\rm cache}[f_{\rm cache}]\bigr).
\]

The second term retains the actual differences in basis derivatives, reference radius and surface metric at the same exact nodal field. Direct pointwise signed subtraction at high precision is valid here. Extending that subtraction naively to independent wide intervals over a parameter tube would lose essential dependence.

The first term is computed from signed monomial moment defects, before taking absolute values. On positive-radius triangles these are Laurent moments with rational coefficients and rigorously enclosed logarithms. On the axis charts they are polynomial moments. Binary moment errors are retained: the constant weighted-measure defect is certifiably nonzero on all 144 triangles. This is a numerical fact, not a new biological principle.

**Receipts:** [144-check component producer](../runs/exact_integration_bridge_v1/moment_component/verification_v2.json); [1,028-check independent component review](../evidence/exact_fe_moments_independent_review.json). The latter includes 15 rigorous complex-integration cross-checks and independent exact barycentric formulas.

## 4. Bound derivatives, not merely the energy remainder

Set \(\varphi(J)=-\mu\log J+\kappa(\log J)^2/2\), and let \(J_i,J_{ij}\) be exact nodal derivatives. The volumetric contributions are

\[
F_i=\varphi'(J)J_i,\qquad
H_{ij}=\varphi''(J)J_iJ_j+\varphi'(J)J_{ij}.
\]

For \(e=J-1\), the convergent series are

\[
\varphi'(1+e)=\sum_{n\ge0}(-1)^{n+1}(\mu+\kappa H_n)e^n,
\]

\[
\varphi''(1+e)=\sum_{n\ge0}(-1)^n(n+1)(\mu+\kappa H_{n+1})e^n,
\qquad H_n=\sum_{k=1}^{n}\frac1k.
\]

The independent body-domain review licenses \(|e|<21/500\) on the existing tubes. Exact harmonic-number and geometric-series majorants bound the omitted derivative tails. If those scalar bounds are \(\tau_1,\tau_2\), the integral-minus-auxiliary force remainder is at most

\[
(M_I+M_Q)\tau_1\|J_i\|_\infty,
\]

and the Hessian remainder at most

\[
(M_I+M_Q)\left[\tau_2\|J_i\|_\infty\|J_j\|_\infty
+\tau_1\|J_{ij}\|_\infty\right].
\]

Both positive masses are essential: the omitted series contributes to the material integral and the finite auxiliary rule. A small scalar-energy remainder alone would not prove a force or Hessian bound.

Sparse moment-dual contraction reuses each constitutive polynomial across the local entries. It is an efficient evaluation of the same finite sums, not an additional approximation. Exact symbolic second-order differentiation checks its Jacobian jets independently.

**Receipts:** [theory and scalar-tail verification](../evidence/exact_integration_target_theory.json); [sparse-jet helper checks](../runs/exact_integration_bridge_v1/bridge_design_probe/polynomial_jet_tests.json).

## 5. Actual independently replayed element results

![Local force and Hessian discrepancies, with included series-tail bounds](../runs/exact_integration_bridge_v1/visuals/bulk-element-integration-defects.png)

At the saved B-center field, conservative rounded upper bounds are:

| Element | Chart | Maximum absolute local force-entry defect | Maximum absolute local Hessian-entry defect |
|---|---|---:|---:|
| 0 | Axis vertex | 1.92e-15 | 1.42e-13 |
| 6 | Positive radius near axis | 2.14e-13 | 3.39e-9 |
| 72 | Axis edge | 1.27e-16 | 1.03e-13 |
| 66 | Outer positive-radius field | 4.95e-13 | 4.22e-12 |

These quantities use normalized energy coordinates. They are not micrometre errors, physical force measurements, or equilibrium-position errors. The exact outward intervals and remainder accounting are retained in the [four-element producer](../runs/exact_integration_bridge_v1/bulk_pilot/results.json). The plot contrasts total bounds with the included truncation tails; increasing series order cannot erase the underlying integration/cache discrepancy.

An independent implementation used order 18 and 640-bit arithmetic, independently reconstructed cubic Jacobian jets and canonical nodal strain maps, and checked all 38 force and 210 upper-Hessian entries. All intervals overlap. The reviewed moment, constitutive-primitive and moment-dual components were reused explicitly; this is not a wholly separate software stack. A further 12,160 canonical basis-entry comparisons passed exactly.

**Receipt:** [independent actual-field replay](../evidence/exact_bulk_pilot_independent_review.json), with [method and accuracy details](../runs/exact_integration_bridge_v1/bulk_pilot_independent_replay/REVIEW.md).

## 6. Complete point result: the old cube is not invariant for the integral map

All 144 bulk elements and all 24 surface edges have now been evaluated and assembled at the same saved B-center field, including both integration and finite-cache discrepancies. The full signed free-coordinate force has 611 entries and the symmetric Hessian has 6,676 upper-triangle positions. Scalar absolute-error bounds were not substituted for the signed sum before applying the preconditioner.

The complete point calculation gives the conservative display enclosure

\[
1.67369209\times10^{-11}
\le \|R(F_I-F_Q)(u_0,\theta_0)\|_\infty
\le 1.67369210\times10^{-11}.
\]

Here \(\theta_0\) is checked to lie in the original five-parameter cell, and the old affine predictor at \(\theta_0\) is precisely the saved \(u_0\). The old uniform cached-force residual bound applies at this point. Using a conservative upper bound that covers both producer and independent residual evaluations,

\[
\|RF_Q(u_0,\theta_0)\|_\infty
\le 5.28745505\times10^{-13}.
\]

Consequently the reverse triangle inequality yields

\[
\|RF_I(u_0,\theta_0)\|_\infty
>1.62081753\times10^{-11}
>r,\qquad r=10^{-11}.
\]

Define the new integral map using the old preconditioner,

\[
T_I(u)=u-RF_I(u,\theta_0).
\]

An explicit witness is free coordinate 588, the axial coordinate of finite-element node 313: its map displacement lies between approximately -1.727e-11 and -1.620e-11, wholly outside the old coordinate interval [-1e-11, 1e-11]. The exact [coordinate witness](../runs/exact_integration_bridge_v1/complete_point_B_independent_review/coordinate_witness.json) uses the independently evaluated old residual bound.

Then \(T_I(u_0)\) is outside the closed old cube centered at \(u_0\). Thus **the new integral map cannot be a self-map of that unchanged cube**. This is stronger than merely observing that a sufficient majorant failed: it exhibits a point whose image leaves the cube. The old quadrature map \(T_Q\) and its existing certificates remain valid for their original target.

![Integral-map displacement versus the old proof-cube radius](../runs/exact_integration_bridge_v1/visuals/complete-point-integral-map-failure.png)

The upper bound for the preconditioned Hessian discrepancy at this point is approximately 1.118e-4. It is a pointwise quantity, not the full-tube bound required for contraction. The failure above does not prove that an exact-integral equilibrium is absent, that it cannot be certified after recentering, or that the physical tissue model is inaccurate by a measurable amount. These are normalized map-displacement scales, not experimental shape-error estimates. They matter because the old numerical proof cube is extraordinarily small. No new experimental precision is inferred.

**Receipts:** [complete point producer](../runs/exact_integration_bridge_v1/complete_point_B/results.json); [complete point independent review](../evidence/exact_complete_point_B_independent_review.json); [144-element bulk assembly review](../evidence/exact_bulk_assembly_B_independent_review.json). The bulk review independently reduces all 611 force and 6,676 Hessian entries and all 611 signed preconditioned force entries. Six bulk elements were separately reintegrated; the remaining 138 use the generic reviewed element engine plus packet checks. Six surface edges were independently reintegrated, with component reuse explicitly disclosed. The full-domain and assembly checks are distinct from reintegrating every element in a second implementation.

## 7. What would actually transfer an equilibrium theorem?

Let the old uniform root proof provide an affine predictor \(p_\theta\), cube radius \(r\), preconditioner \(R\), and bounds

\[
\|RF_Q(p_\theta)\|_\infty\le\eta,\qquad
\sup\|I-RH_Q\|_\infty\le Z<1.
\]

The missing quantities are validated, uniform bounds

\[
a\ge\sup_\theta\|R(F_I-F_Q)(p_\theta)\|_\infty,
\qquad b\ge\sup_{\theta,u}\|R(H_I-H_Q)(u)\|_\infty.
\]

If, on a genuinely verified proof domain,

\[
Z+b<1,\qquad \eta+a+(Z+b)r<r,
\]

then \(u\mapsto u-RF_I(u)\) maps the cube strictly into itself and contracts. Banach's theorem supplies a unique exact-integral FE root there. Together with the previously verified positive symmetric part of \(R\), the same Lyapunov argument supplies a strict local minimum. This statement is a conditional theorem; the actual inequalities have not yet been verified.

The complete point calculation does not supply \(a\) or \(b\). Parameter dependence and full-tube cache correction remain necessary. In this B case the unchanged self-map is already ruled out at its center, so direct transfer requires a newly justified map/preconditioner, numerical center or proof domain. Fresh geometry, root, energy and observation bounds must accompany that change; the physical material/gain cell cannot be silently narrowed to make the test pass.

A local bulk cache-consistency primitive has also been prepared and independently checked. It uses the small coefficient difference first,

\[
\delta q=(B_*-B_c)Y,
\]

and the mean-value integral of the constitutive third tensor to bound \(H_*-H_c\), preserving the same nodal field and shared parameters on both sides. Its independent review has 150 checks. Two actual unweighted auxiliary sites are covered, both inactive; their 16 endpoint evaluation checks represent eight distinct nodal witnesses, with the duplicated zero-drive endpoints disclosed. Active-drive dependence is separately exercised by synthetic interval fixtures. No full-cache or full-tube assembly follows merely from these local tests.

**Receipt:** [local uniform cache-consistency review](../evidence/exact_cache_consistency_independent_review.json).

## Remaining scientific obligations

The exact-integration bridge is numerical-mechanics groundwork. It does not close the continuum discretization error, common stable loading path, full physical parameter/geometry identified set, source-matched reference/ROI/time contract, larger-region prediction test, or novelty assessment. The source-data request remains unsent pending the user's authorization and sender details. The original five-gate contract remains unchanged.

The scientific target remains a defensible test of whether local calibration determines larger-scale tissue geometry, with a certified ambiguity and a separating measurement if it does not. Neither ambiguous mechanical fitting in general nor validated quadrature is claimed as a new concept.


Independent verification here means internal, separately implemented AI checks, with component reuse disclosed. No outside expert review or new biological experiment is claimed. No author contact or GitHub publication occurred in this step.
