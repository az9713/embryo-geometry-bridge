# Which force ratios affect the smooth mean responses?

**The saved model has strongly unequal local response sensitivities.** A conditional tangent calculation at strengths1,5,10 produces a six-by-six smooth-mean response matrix whose largest and smallest singular values differ by roughly four orders of magnitude in the declared coordinates. This is a numerical local conditioning diagnostic—not a fitted material confidence region, a global identifiability theorem or a measured biological property.

The [mechanical tangent implementation](../scripts/material_tangent_inference.py), [force-derivative/linear-system verification](../evidence/material_tangent_mechanics_verification.json), and [mean-response diagnostic](../evidence/material_tangent_mean_diagnostic.json) are retained. No mechanical optimization, new equilibrium or source-outcome fit is performed by this module. It evaluates three already saved 12×6 higher-rule fields.

An [independent 111-check review](MATERIAL-TANGENT-REVIEW.md) verifies the current force components, all 18 implicit tangent columns, mean chain rule and source/input/archive hashes. It identified a future-reuse precondition gap in the low-level `means()` function. Use the separate [guarded entry point](../scripts/material_tangent_guarded.py) for reuse: it verifies all retained source/input/archive hashes before a replay, checks them afterward, and archives both old and new outputs. Its default read-only check and simulated stale-dependency probes are retained in [the guard receipt](../evidence/material_tangent_guard_verification.json). The original audited mechanical source and protocols remain unchanged.

## Conditional implicit differentiation

Let \(u\) be the free displacement coefficients, \(\theta\) the force-ratio coordinates, and let

\[
R(u,\theta)=\partial_u\frac{E(u,\theta)}{\mu V_0}=0.
\]

If an **exact** equilibrium \(u_*\) is regular, so that \(H=\partial_uR(u_*,\theta)\) is invertible, the implicit-function theorem supplies a local branch and

\[
H\frac{\partial u_*}{\partial\theta_j}
=-\partial_{\theta_j}R.
\]

For an observable \(O\), the chain rule gives

\[
\frac{dO}{d\theta_j}
=O_u\frac{\partial u_*}{\partial\theta_j}+O_{\theta_j}.
\]

The current material changes hold geometry, references, current-radius ROI and length unit fixed, so the smooth mean observer has no explicit parameter term. Future height, ROI, reference, activation-support or geometry changes generally require one; this module does not omit such terms while claiming to differentiate those different problems.

The actual computation uses floating Hessians and residuals at tolerance-verified stored fields. Small force residuals and a successful linear solve are not a validated enclosure of an exact root or its inverse Hessian. Thus the theorem above is conditional, and the output is a diagnostic approximation to its derivative. No rigorous algebraic/PDE stability claim is made.

## The six declared coordinates

The shear coefficient is fixed as the force unit. The coordinates are log compression coefficient, log apical and basal area coefficients, additive apical tension, log basal tension, and log common activation gain. All lengths use the fixed \(L_0=33.5\) µm convention. Absolute common force scale is invisible to static shape and was removed by the [force-scale quotient](PHYSICAL-PARAMETER-DOMAIN.md).

Apical baseline tension is zero at these fields. Its physical variation is one-sided and additive; taking its logarithm would be undefined. An unconstrained singular-value decomposition does not enforce that nonnegative tangent cone. The other log coordinates describe relative coefficient changes. Their mixed normalization is explicit and affects singular values/vectors; no coordinate-invariant statistical meaning is attached to them.

The finite compression coefficient is a parameter of this declared law. Its sensitivity is not a measured physiological compressibility, nor permission to tune a numerical incompressibility penalty against held-out outcomes.

## Computing parameter-force derivatives without new solves

Each physical coefficient enters the energy linearly on the accepted positive-J states. For a positive coefficient \(c\), doubling it at a fixed field gives

\[
\partial_{\log c}R=R(u;2c)-R(u;c).
\]

For an additive tension, a unit increment gives its linear coefficient derivative. At fixed input dose, a relative gain change gives \(\partial_{\log g}R=f\,\partial_fR\), computed by doubling the active strength. These are fixed-field force evaluations, not refitted equilibria.

The source computes six right-hand sides \(B\), factors the unchanged sparse Hessian, and solves \(HS=-B\). All three states pass replayed force checks, six independent fixed-state finite-difference checks and the linear-system residual checks. The observed central log-difference error is approximately \(1.67\times10^{-7}\), consistent with \(\sinh(10^{-3})/10^{-3}-1\) for a coefficient-linear force. Linear solve relative residuals are below \(10^{-11}\).

The common-force scaling direction also cancels after including the normalization's shear derivative. This verifies the expected gauge numerically; it is not an additional physical observation.

## Differentiating the mean at fixed current radius

The [separate geometric derivative module](QUADRATIC-MEAN-DERIVATIVES.md) gives exact interval gradients with respect to physical surface-node positions. For one monotone parametric edge \((R(t),Z(t))\), a fixed-current-ROI weighted integral has variation

\[
\delta\int w(R)Z\,dR
=\int w(R)\bigl(R'\delta Z-Z'\delta R\bigr)\,dt.
\]

Moving endpoint contributions cancel because the physical ROI, rather than its material parameter, is held fixed. Neglecting the \(-Z'\delta R\) term would differentiate the wrong observation. With linear or disk weights and quadratic edges, the integrands are polynomial; rational inverse brackets handle the ROI cuts.

The root module takes gradient interval midpoints, multiplies by \(L_0\) to obtain derivatives with respect to dimensionless solver displacements, and selects the actual free degrees of freedom before multiplying by \(S\). The exact geometric-gradient bounds do not make the floating mechanical inverse or complete composed sensitivity rigorously enclosed.

## Numerical diagnostic and interpretation

There are two smooth mean channels per strength: apical and basal. Their height mean is their sum and is not counted as a third independent channel. Stacking these at three diagnostic strengths gives six rows. The singular values are approximately:

| Mean weighting | Singular values, descending |
|---|---|
|Projected disk|11.06445,1.85661,0.13353,0.04522,0.01216,0.001439|
|Radial line|12.12863,1.97683,0.14870,0.05767,0.01229,0.001310|

The weakest computed right-singular directions are predominantly basal-area-coefficient changes mixed with other coefficients. This says the specified local smooth-mean map responds weakly to some combinations in these coordinates. It does not prove exact nonidentifiability or that the basal area coefficient is biologically unmeasurable. Full spatial profiles, different inputs, maxima, time dependence and geometric measurements can carry different information.

These six rows are not the experimental 68-summary Jacobian or the completed17-dose cohort. Spatial maxima may be nonsmooth at competing extrema, and this module does not invent their derivatives. The source's component and time-aggregation interpretations also remain unresolved. Even numerical full rank cannot justify a physical confidence interval without validated forward error, correct observables, a discrepancy model and global analysis.

This component supplies a route toward material-ratio inference. It does not satisfy the original full compatible-set/propagation requirement on its own and establishes no novelty.
