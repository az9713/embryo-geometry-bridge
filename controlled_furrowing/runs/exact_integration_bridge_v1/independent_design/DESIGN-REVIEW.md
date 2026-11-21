# Exact-integration bridge: independent target and feasibility review

This is a design and read-only scale audit. No exact-integration error has yet been bounded. No new equilibrium, Hessian assembly or numerical proof is produced here. The existing17root/stability/body certificates concern the frozen binary-cache finite-quadrature potential. They do not already certify the exactly integrated finite-element energy, let alone a continuum equilibrium.

## 1. Freeze an actual new mathematical target

For each reference triangle use its three stored binary vertex coordinates as exact rationals, the exact affine reference map, and exact barycentric P2 Lagrange basis in the existing six-node order. The current coordinate field has the stored binary nodal X plus exact free displacement, including all stored midpoint roundings. Do not silently replace those midpoint coordinates with ideal affine midpoints. Existing whole-body proofs already use an affine reference and exact P2 current map, so their field interpretation is compatible with this choice. The stored X interpolant need not be exactly the affine identity map: that tiny fact belongs in the target definition.

Fix the existing free/fixed coordinates, mesh connectivity, stored element activation mask, common material/gain cell, specimen dose and physical coefficients. Use a declared positive normalization; retaining the frozen positive binary energy scale is legitimate, but the exactly integrated numerator must use actual2π, reference Jacobians and exact basis maps. Different normalization changes all gradient/Hessian budgets and cannot be silently substituted. The prior clipped log and inversion penalty vanish only after the **new target's** whole-tube positivity is established; they are not part of an analytic continuation through inversion.

The exact potential is the sum of the axisymmetric volume integrals and apical/basal surface integrals, with the same logarithmic-volume, neo-Hookean, active-fiber and membrane constitutive terms. Its variations must be taken with respect to the same611free coordinates. Fixing a new exactly integrated FE model still leaves interpolation/discretization/PDE, source observation and biological errors separate.

## 2. Computable persistence and minimum theorem

Write FQ=∇EQ for the already certified normalized cache potential, FI=∇EI for the proposed exactly integrated FE potential, and R for the same exact binary preconditioner. For each parameter θ let pθ be its frozen affine predictor. On the closed cube Bθ={u:||u−pθ||∞≤r}, existing independent proofs provide

    ||R FQ(pθ)||∞ ≤ η,    sup_Bθ ||I−R HQ(u)||∞ ≤ Z <1.

Define new, genuinely validated uniform defect bounds

    a ≥ sup_θ ||R(FI−FQ)(pθ)||∞,
    b ≥ sup_θ,u∈Bθ ||R(HI−HQ)(u)||∞.

If

    Z+b <1,     η+a+(Z+b)r <r,

then Tθ(u)=u−R FI(u) is a strict self-map and contraction of Bθ. Banach gives a unique fixed point for every θ. Since ||I−R HI||<1 implies R HI invertible, R is nonsingular and the fixed point is an actual FI=0 root. Its distance from the predictor is at most (η+a)/(1−Z−b), and the root lies in the same certified geometric cube. A weaker non-strict self-map can suffice for existence, but this review proposes retaining strict margins for robust replay.

The separately independently proved sym(R)>0 transfers directly. Since HI is a real symmetric energy Hessian and A=R HI has spectrum in Re(z)>0, HI^-1 solves A X+X Aᵀ=R+Rᵀ>0. The convergent Lyapunov integral makes X positive definite. Thus HI is positive definite throughout the tube and each root is a strict local energy minimum. Quantitatively

    HI ≥ [(1−Z−b)/||R||∞] I.

This avoids an unnecessarily severe independent unpreconditioned Hessian comparison. An alternative sufficient criterion is ||HI−HQ||2<m, where m is the existing uniform eigenvalue lower bound; a symmetric row-sum norm bounds that spectral norm. Neither route proves non-axisymmetric stability, a continuum minimum, dynamics or a common parameter-dependent history.

A convenient **sufficient, non-optimal** split is a≤s/4 and b≤s/(4r), with s=r−η−Zr>0. Then the self-map retains at least s/2margin. Using raw defects eF=||FI−FQ||∞ and eH=||HI−HQ||∞ instead gives a≤||R||∞eF and b≤||R||∞eH. These absolute-value reductions may be much too pessimistic; the useful target is a directly preconditioned, correlation-preserving error bound.

## 3. Actual scales, not aspirational tolerances

`audit_scales.py` reads all17 hash-linked root, observer, independent-root and independent-stability records. It does not regenerate a model, assemble a Hessian or solve. Exact rational arithmetic reproduces all17m=(1−Z)/||R||∞ values and all positive self-map margins. `scales.json` retains every exact fraction and hash binding; displayed figures below are rounded diagnostics, not certificate endpoints.

| Quantity | Record | Approximate bound |
|---|---|---:|
| Smallest remaining self-map margin s | G |4.2249217×10^-12 |
| Smallest sufficient raw-force budget s/(4||R||∞) | G |6.7877707×10^-16 |
| Smallest sufficient preconditioned Hessian budget s/(4r) | L |0.02921307 |
| Largest full-tube Z | L |0.38314772 |
| Smallest normalized Hessian lower bound m | L |0.00039487421 |
| Largest ||R||∞ | B |1686.9312 |
| Smallest whole-body J lower bound | L |0.95829092 |
| Largest whole-body J upper bound | L |1.02340445 |
| Smallest geometric strong-monotonicity lower bound | L |0.27276309 |

G is most restrictive for force-defect preservation of the **existing** tube; L is the geometric/conditioning and preconditioned-Hessian-budget case. The smallest raw-force budget is already comparable to binary64 rounding of order-one quantities and far below many uncancelled κ=5000intermediates. More floating precision shrinks arithmetic enclosures but cannot shrink true quadrature or cached-map defects. It is entirely possible that no unchanged old tube meets these inequalities. Failure would mean only that this sufficient bridge failed: explicitly enlarge/recenter a proof tube or improve the norm, then revalidate geometry and observation transport. It would not prove nonexistence of an exactly integrated equilibrium.

## 4. Exact geometry/cache consistency cannot be waved away

Actual policy is8-point tensor-Duffy Gauss in volume (64points/triangle),12-point surface Gauss,144P2triangles and325nodes. Legacy parameter labels quadrature4/surface8 do **not** describe the active rule. There are6axis-edge triangles,6triangles touching the axis only at one vertex,132off-axis triangles and24active triangles. Exact binary-coordinate checks show that the specified radial/depth interfaces cut no triangle interior and the frozen mask agrees with element containment in that rectangle. Retaining the cached mask is still the most explicit target; future geometric changes require a fresh interface audit.

Three exact cache discrepancies were checked, rather than presumed negligible:

* Stored midpoint coordinates differ from their exact binary-vertex averages by as much as2^-53;192element-coordinate occurrences are nonzero.
* The cached shape partition has max|ΣN−1|=45/2^57≈3.1225×10^-16.
* The cached derivative partition has max|Σ∇N|≈9.1184×10^-15.

These are exact rational facts about frozen binary arrays, not yet propagated force/Hessian errors. A bridge must include shape, derivative, radius, reference inverse/Jacobian, quadrature-node/weight, surface-tangent and normalization consistency. Replacing the cached arrays by exact P2 formulas changes the mathematical target; existing Arb proof of the old arrays does not bound that change automatically.

## 5. Axis and interface treatment

No quadrature error theorem using uniform log/inverse derivatives may directly divide by r on an axis-containing closed triangle. On an axis-edge triangle, fixed radial axis nodes make the P2 radial polynomial divisible by the reference radial linear factor. Cancel that factor symbolically in R/r **and in every free-coordinate variation** before interval differentiation.

For an axis-vertex-only triangle, R/r generally is neither a polynomial nor a function with one spatial corner limit. Example: an affine map R=r+z on a triangle with vertices(0,0),(1,0),(1,1) gives R/r=1+t after the Duffy substitution r=s,z=st. The corner s=0 has direction-dependent t but the transformed function is smooth on the square. In general choose the zero-radius vertex as λ0=1−s, λ1=s(1−t), λ2=st. Then r=s rbar(t), R=s Rtilde(s,t), with rbar strictly positive. Cancel s; retain the Duffy Jacobian s and prove denominator/tangent/J positivity on the transformed closed square. A single-axis-limit assumption would exclude real elements in this mesh.

Likewise cancel the surface axis factor before AD: in the weighted membrane energy, r0l0(K/2)(Js−1)^2 can be written K/2[R²|y′|²/(r0l0)−2R|y′|+r0l0], and the tension part is T R|y′|, apart from2π. The axis-edge radial factor is removable, but tangent nondegeneracy must still be proved. Existing body positivity and surface graph checks supply useful domain information; a new interval integrator must consume their exact bounds, not merely a sampled J.

The activation jump must be split along true element interfaces. Here element alignment was checked exactly. It is invalid to apply a smooth high-order remainder across an uncut material/activation jump in a future mesh.

## 6. A validated, computationally credible integration route

On each axis-regularized element/edge, build spatial polynomial or Taylor models for the **gradient and Hessian integrands**, including the actual uniform displacement/parameter tube. Exact reference-polynomial pieces and moments can be integrated rationally. Bound log, inverse and square-root remainders on positive compact domains, using outward Arb arithmetic and adaptive degree/subdivision only with retained accepted/rejected cells and complete coverage. A scalar energy error alone does not control its derivatives; derivative models require their own justified remainders.

Apply the quadrature-error functional L(f)=∫f−Q(f) to each polynomial **before** coefficient interval absolute values. Frozen binary Gauss nodes and weights are not an exactly degree15/23 Gaussian rule as rationals. Every monomial must retain its exact residual moment, rather than being set to zero from nominal Gaussian exactness. Alternatively certify ideal algebraic Gauss nodes/weights, then explicitly bound ideal-rule-to-binary-rule rounding. Refined spatial partitions do not inherit exact cancellation from the original global rule: either retain a global error functional or account for all refined-rule differences.

Separate the target-consistent rule Qc from frozen EQ and prove

    EI−EQ = (EI−EQc)+(EQc−EQ),

at the gradient and Hessian levels. The second term is a local, same-coordinate consistency comparison of exact P2/radius/weight maps with the cached maps, not subtraction of unrelated interval totals.

Naively integrating interval FI over the full parameter/predictor box and subtracting an independent interval FQ loses shared dependence: widths of order HΔu can dominate a10^-12self-map margin. Preserve parameter/predictor dependence in the polynomial error functional, or use a point defect plus rigorously bounded derivatives of the **defect** over the full tube. Sparse local12-coordinate bulk and6-coordinate surface maps can be contracted with the appropriate columns of dense R before final absolute values. Bound611gradient outputs and the induced Hessian matrix norm with a documented assembly identity. Independent interval differences atQ8/Q12, or high-order sampled agreement, do not establish this bound.

The model author's independently proposed route agrees on exact affine-reference P2 with stored nodal X, fixed χ, axis-adapted Duffy cancellation and spatial series/exact moments. In discussion the author explicitly confirmed that binary-rule low-degree moment defects must remain nonzero and that unchanged tiny tubes are not promised to pass. No implementation or numerical integral certificate has been independently reviewed yet.

## 7. Acceptance packet and stopping criterion

Each record needs a versioned target/protocol manifest binding old cache, predictor, full parameter cell, dose, R, exact geometry/interface choice, normalization, source versions, analytic domains and chosen error norm. Retain point/predictor a and full-tube b separately, all quadrature subdivisions/remainder budgets, binary consistency error and directed rounding evidence. Pass only exact rational/outward-ball inequalities above, with untouched source packets and an independent implementation/replay of the decisive error calculation. Precision is adaptive: there is no justified universal claim that160/192bits suffice until the actual remainder budgets close.

If the existing tube fails, its failure remains; a new proof domain must separately carry root, body geometry, surface/ROI and all68observation-error implications. One successful field would establish one integration bridge, not the17record cohort. Passing every record would advance finite-element integration consistency, while full material/geometry identified sets, tails, continuum discretization error, source temporal interpretation, causal prediction and novelty would still require their original work.
