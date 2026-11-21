# Extending the integral proof: the difficult L case

**The L calculation is complete and independently checked, but its original-cube root test is unresolved.** All 144 volume elements and 24 surface edges were integrated with validated local enclosures on the original L nodal box and the unchanged small parameter cell. Independent exact endpoint assembly confirms the producer's bounds. This extends the computational coverage; it does not extend the successful [B root theorem](EXACT-INTEGRAL-ROOT.md) to L.

![Contraction passes but the L self-map enclosure fails, including a hypothetical zero-force comparison](../runs/exact_integration_cohort_v1/checkpoint_L_v1/root_test_limits.png)

## 1. Same model, explicitly different record

The source record is `optoGEF2_small_col_L`, with measured dimensionless activation input \(I=1.51057344\). Its 611-dimensional predictor, binary preconditioner, old cached-quadrature proof and exact nodal domain are loaded from that record's selected artifacts. F and K require their own corrected centers; M has a separate reviewed root provider. A generic loader was checked against all 17 selected records before this L job. Reusing B's dose, center or force packet would be invalid.

The L box has positive whole-body Jacobian bounds, approximately \([0.95829,1.02340]\). All 144 volume domains, 24 surface domains and 288 surface-cache sites were checked before force integration. Geometry admissibility is a premise for evaluating this model, not proof that the model describes the tissue.

The scalar tail budgets were derived before execution from the exact preconditioner norm, the old root radius, the old contraction bound and the actual finite-element incidence counts. Bulk expansion orders ranged from 5 to 15; surface orders ranged from 7 to 19. These orders certify the declared truncation tails. They do not by themselves control interval dependency or finite-element error relative to the continuum.

## 2. The sufficient theorem and its actual outcome

Let \(F_Q\) denote the force for the frozen binary quadrature functional and \(F_I\) the exactly spatially integrated finite-element force. On the declared domain set \(D=F_I-F_Q\). The old parameter-matched predictor \(p(\theta)\), its radius \(r\), and preconditioner \(R\) satisfy the reviewed cached bounds

\[
 \|R F_Q(p(\theta);\theta)\|_\infty\le\eta_0,
 \qquad \|I-RH_Q(u;\theta)\|_\infty\le Z_0.
\]

The new component enclosures give bounds \(d\ge\sup\|RD\|_\infty\) and \(b\ge\sup\|R(H_I-H_Q)\|_\infty\). Thus

\[
 \eta=\eta_0+d,\qquad z=Z_0+b.
\]

For fixed \(\theta\), the map \(T(u)=u-RF_I(u;\theta)\) is a contraction and maps the predictor cube into itself if

\[
 z<1,\qquad \eta+zr<r.
\]

Indeed the derivative bound gives Lipschitz constant \(z\), and the center displacement plus that constant times the radius bounds every image displacement. The closed finite-dimensional cube is complete, so Banach's theorem gives a unique fixed point there. The bound \(\|I-RH_I\|_\infty<1\) also makes \(RH_I\) invertible by the Neumann series; since these matrices are square, \(R\) is invertible. Thus a fixed point is exactly a force root. This is a sufficient existence test: its failure does not imply absence of a root. The verified contraction still implies at most one root per parameter-matched old cube, without establishing existence.

| Independently replayed quantity | Rounded display |
|---|---:|
| Derivative-defect bound \(b\) | 0.1851092833 |
| Combined contraction bound \(z\) | 0.5682570024 |
| Original radius \(r\) | \(6.011722185\times10^{-11}\) |
| Self-map left-hand side \(\eta+zr\) | \(6.563471613\times10^{-7}\) |
| Left-hand side divided by radius | 10917.7893 |

Here \(r,\eta_0,d\) use normalized nodal-displacement coordinates, not micrometres. Therefore contraction passes and self-mapping is unresolved. Even setting the predictor-force contribution to the residual bound to zero while retaining the same derivative bound gives \((\eta_0+zr)/r\approx1.068257\). Tightening force alone cannot finish this particular proof. For that hypothetical zero contribution, the strict test would require \(b<1-Z_0-\eta_0/r\approx0.11685228\).

The producer used outward Arb arithmetic. The reviewer independently scattered exact rational endpoints, contracted with the exact binary entries of \(R\), and checked signed enclosure containment. Absolute Hessian row bounds are checked by upper domination, because widening an Arb ball can also increase the lower endpoint of a subsequently computed upper bound. The resulting reduction covers all 168 packets, 611 force coordinates and 6676 symmetric upper Hessian entries. It shares the already reviewed local enclosure and geometry engines; it is not a separate physical implementation.

## 3. A wide enclosure is not a measured error

Every one of the 611 signed preconditioned defect intervals contains zero. Their certified norm lower bound is therefore zero. The widest interval, at free coordinate 190 (node 105, axial component), has width approximately \(1.312548304\times10^{-6}\), almost entirely from volume terms. Its midpoint is not a force evaluation at the proof center.

Element 21 is the largest local contributor to that width. A separately frozen, prelaunch-reviewed force-only diagnostic evaluated that element at the exact proof center and over the whole nodal box. Their scalar tail bounds, both about \(1.57\times10^{-19}\), are much smaller than the original whole-box interval width. Reducing the scalar tail alone would leave the dominant interval width in the volumetric truncation contribution.

## 4. Centering recovers useful local information

For this element let \(G(Y,g)\) be its local integral-minus-cache force, where \(Y\) is the local nodal vector and \(g\) the common activation gain. A path first from \(Y_0\) to \(Y\) at fixed \(g_0\), then from \(g_0\) to \(g\) at fixed \(Y\), gives the exact identity

\[
 G(Y,g)=G(Y_0,g_0)
 +\int_0^1 H_D(Y_0+t(Y-Y_0),g_0)(Y-Y_0)\,dt
 +(g-g_0)\,\partial_gG(Y,g).
\]

The last equality uses the affine gain dependence of this bulk energy. The entire nodal segment is inside the same convex box. Its positive-J guarantee and the frozen whole-box Hessian packet justify the derivative bound. The exact gain derivative uses the dose-scaled difference between integral and cache axial stiffness operators on the same geometry. The point, state-increment and gain-increment terms are enclosed separately; intersecting this bound with a direct enclosure preserves validity.

| Local element-21 quantity | Largest value across 12 force coordinates |
|---|---:|
| Direct whole-box interval width | \(5.23196888\times10^{-8}\) |
| Centered whole-box interval width | \(9.32225561\times10^{-12}\) |
| Fixed-center point interval width | \(3.13136626\times10^{-19}\) |
| Fixed-center absolute defect upper bound | \(2.17060600\times10^{-14}\) |

The table uses normalized local generalized-force coordinates. The maximum width improves by about 5612 times. All 12 centered intervals lie within both direct enclosures and still contain zero. At the fixed center, the much narrower point intervals certify five positive and seven negative local defects. Those local signs do not determine the assembled global force or its root. No new root test is inferred from this one-element calculation.

## 5. An equivalent scaled chart tightens the derivative enclosure

A second, separately reviewed pilot changed only the exact algebraic representation of element 21. Write physical material coordinates as \(R=R_c+h x\), \(Z=Z_c+k y\), with exact rational constants chosen from the element's coordinate ranges. The dimensionless radial denominator \(\widehat R=R/R_c\) stays positive. Exact quadratic basis functions and their physical derivatives are constructed in \((x,y)\) before nodal intervals are inserted. The hoop strain retains its radial denominator. Strain, Jacobian and first/second derivative identities were checked exactly; they represent the same force and Hessian, not a modified constitutive law.

Signed moments of \(x^a y^b/\widehat R^d\) are pulled back to the reviewed physical rational/log moments. The physical radial measure is counted once. The original auxiliary quadrature sites, binary weights, cache correction, parameter cell, nodal box, order 12 and tail budgets are retained. Thus the two representations enclose the same quantities, and their interval intersection is valid.

| Largest local bound | Original representation | Scaled representation |
|---|---:|---:|
| Force interval width | \(5.23197\times10^{-8}\) | \(1.72151\times10^{-11}\) |
| Hessian entry interval width | \(3.00324\times10^{-3}\) | \(2.95026\times10^{-8}\) |
| Absolute Hessian row bound | \(8.12987\times10^{-3}\) | \(5.67727\times10^{-8}\) |

Independent endpoint replay confirms that all 12 force intervals and 78 upper-Hessian intervals lie inside their original counterparts. The largest absolute Hessian row bound decreases by about 143200 times. This demonstrates representation-dependent enclosure looseness, not increased physical precision. The [actual independent review](../runs/exact_integration_cohort_v1/L_defect_diagnostic/scaled_chart_v1/actual_independent_review.json) retains every intersection and tail check. No whole-L reduction is implied by this local result.

## 6. Scientific significance and the next obligation

This result distinguishes a failed enclosure from a demonstrated physical discrepancy and identifies a concrete numerical remedy. It is an application of standard interval and mean-value methods, not an established new biological concept. The next mathematical obligation is to tighten both force and derivative enclosures without changing the model, parameter cell or data contract. The next scientific obligation remains source-compatible calibration and a frozen, uncertainty-aware larger-region prediction. More numerical digits do not supply missing timing, references, ROI registration or measurement uncertainty.

## Retained receipts and replay sources

| Step | Evidence |
|---|---|
| Explicit record loader | [Independent context review](../evidence/exact_cohort_context_independent_review.json) |
| Local volume/surface integration | [Independent integrator review](../evidence/exact_cohort_integrator_independent_review.json) |
| Fresh L domain and expansion selection | [Domain/integration prelaunch review](../runs/exact_integration_cohort_v1/records/optoGEF2_small_col_L/original_cube_v1/integration_prelaunch_independent_review.json) |
| True terminal integrations | [144-volume terminal](../runs/exact_integration_cohort_v1/records/optoGEF2_small_col_L/original_cube_v1/bulk_terminal.json), [24-edge terminal](../runs/exact_integration_cohort_v1/records/optoGEF2_small_col_L/original_cube_v1/surface_terminal.json) |
| Independent exact endpoint reduction | [Actual result](../runs/exact_integration_cohort_v1/reduction_preparation/optoGEF2_small_col_L_reduction_review.json), [reducer prelaunch](../evidence/exact_cohort_reducer_prelaunch_independent_review.json) |
| Signed defect diagnosis | [Exact endpoints](../runs/exact_integration_cohort_v1/reduction_preparation/optoGEF2_small_col_L_signed_endpoint_diagnostics.json) |
| Centered element-21 diagnostic | [Frozen protocol](../runs/exact_integration_cohort_v1/L_defect_diagnostic/protocol.json), [independent prelaunch](../runs/exact_integration_cohort_v1/L_defect_diagnostic/prelaunch_independent_review.json), [actual independent review](../runs/exact_integration_cohort_v1/L_defect_diagnostic/actual_independent_review.json) |
| Exact scaled material chart | [Frozen protocol](../runs/exact_integration_cohort_v1/L_defect_diagnostic/scaled_chart_v1/protocol.json), [independent prelaunch](../runs/exact_integration_cohort_v1/L_defect_diagnostic/scaled_chart_v1/prelaunch_independent_review.json), [actual independent review](../runs/exact_integration_cohort_v1/L_defect_diagnostic/scaled_chart_v1/actual_independent_review.json) |
| Chapter numerical composition and figure | [223-binding receipt](../evidence/exact_cohort_L_checkpoint.json), [composition source](../scripts/cohort_L_checkpoint.py) |

The B proof remains the one completed exact-integral root result. The five-gate project remains incomplete; no continuum, source-compatible witness pair, larger-region success or confirmed novelty follows from this L checkpoint.
