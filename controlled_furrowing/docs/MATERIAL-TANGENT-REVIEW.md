# Independent review of the local material-tangent diagnostic

**Verdict: PASS for the current floating, fixed-geometry diagnostic.** The source, archived matrices, current saved inputs, composed smooth-mean derivatives and current explanatory document were independently checked. This is an internal AI review, not external scientific validation or a completed experimental inverse problem.

The review performs no optimization, new equilibrium calculation or fit. Its retained program and results are inside `evidence/material_tangent_independent_review.json`, allowing the numerical checks to be replayed without changing frozen sources. All **107 numerical/provenance checks pass** before the final documentation checks described below.

## What mathematical derivative is being approximated?

For the fixed free-displacement space, references, geometry, quadrature and activation support, define

\[
 R(u,\theta)=\nabla_u\frac{E(u,\theta)}{\mu V_0},\qquad
 H=\partial_uR,\qquad B=\partial_\theta R.
\]

At an **exact regular root**, differentiability and invertibility of H imply a local equilibrium branch with

\[
 H S=-B,\qquad S=\frac{\partial u_*}{\partial\theta}.
\]

The accepted positive-J states lie on the smooth branch of the implemented potential; the reference activation indicator and free-coordinate constraints stay fixed under these material variations. The additive apical-tension coordinate at zero is physically one-sided. One may use the smooth mathematical extension of the potential to formulate the local theorem, then restrict the branch to nonnegative tension; an unconstrained parameter SVD does not impose that restriction.

Actual calculations use small-residual saved fields and floating matrices. Neither a successful sparse factorization nor a small force residual encloses an exact regular root. The displayed S is therefore a **floating linearized response at the stored field**, not a rigorously enclosed equilibrium derivative. This qualification is explicit in the originating source and document.

## Independently assembled parameter-force derivatives

The six coordinates, in order, are log finite compression coefficient, log apical area coefficient, log basal area coefficient, additive apical tension, log basal tension and log common activation gain. Shear coefficient mu is fixed as the force unit. The compression coefficient is the finite coefficient of the logarithmic-volume term, not an independently measured physiological bulk modulus.

At the saved 12x6 high-quadrature gain-5 field, the review independently assembles six scalar component potentials and differentiates them by automatic differentiation. It shares the verified geometric P2 quadrature cache but never calls the production component-force assembly to form these columns:

\[
\begin{aligned}
 E_\kappa&=\int\frac\kappa2(\log J)^2\,dV_0,\\
 E_{K_a}&=\int_{\Gamma_a}\frac{K_a}2(J_s-1)^2\,dA_0,
 &E_{K_b}&=\int_{\Gamma_b}\frac{K_b}2(J_s-1)^2\,dA_0,\\
 E_{T_a,\mathrm{add}}&=\int_{\Gamma_a}dA,
 &E_{T_b,\log}&=T_b\int_{\Gamma_b}dA,\\
 E_{g,\log}&=\int\frac{f\chi}{2}|F n_0|^2\,dV_0.
\end{aligned}
\]

For each component, the column is \(\nabla_uE_j/(\mu V_0)\). The additive apical-tension column is the derivative of true area, not the derivative of a zero energy term at zero tension. Separate apical/basal Gauss loops use current revolved area \(2\pi R|\partial_t y|dt\) and reference area to form the surface terms.

All six independent columns match the archived B: the largest relative column difference is **1.21e-13**, and the largest absolute entry difference is **4.58e-16**. This corroborates the production doubling identity for coefficient-linear energies. Central log finite differences have the predictable factor

\[
 \frac{R(ce^h)-R(ce^{-h})}{2h}
 =\frac{\sinh h}{h}\,\partial_{\log c}R,
\]

so their reported relative error near 1.67e-7 at h=0.001 is expected truncation, rather than unexplained derivative error. The additive coefficient check is one-sided because its baseline value is zero.

The independently assembled shear contribution also reproduces the common-force gauge cancellation. If all physical coefficients, active strength and mu scale together, the normalized force is unchanged. With inactive trial-inversion penalty, the mu column contributes the shear gradient minus the normalized full force. Adding that term to the six coefficient columns, with apical tension multiplying its additive column, leaves a maximum residual **4.00e-16**. This is a force-scale symmetry check, not new biological information.

## Fresh Hessian-vector checks of all three saved tangents

For each saved gain 1, 5 and 10, a fresh automatic-differentiation Hessian-vector product differentiates the actual frozen potential along every archived column of S. It is independent of sparse coloring, Hessian recovery and the archived sparse entries. All eighteen columns satisfy the retained relative tolerance 1e-9 for \(HS+B\); the largest observed relative residual is **9.85e-12**. Their fresh HVPs also agree with multiplication by the archived sparse H.

The review separately checks that every archived free vector equals the saved field, every free-index array equals the current mechanical model, the actual integration rule is Duffy-8/surface-12, fresh forces pass the unchanged 1e-7 gate, and all current input and NPZ hashes equal the locked protocol/receipts. No stored field or production archive was modified.

## Units and the composed observation derivative

The [exact smooth geometric gradients](QUADRATIC-MEAN-DERIVATIVES.md) are derivatives of a mean measured in micrometres with respect to physical node coordinates measured in micrometres. If the mechanical displacement is dimensionless and \(x_{\rm phys}=L_0(X+u)\), then

\[
 D_uO=L_0D_{x_{\rm phys}}O,
 \qquad D_\theta O=(L_0D_{x_{\rm phys}}O)_{\rm free}\,S,
 \qquad L_0=33.5\,\mu\mathrm m.
\]

The root code uses the correct scale, flattening order and free-coordinate restriction. The current parameter variations hold ROI/reference/length unit fixed, so no explicit observation-parameter term is missing for this declared calculation. Varying height, ROI, references or activation support would require their separate explicit derivatives.

For each of the eighteen material directions and each of line/disk weighting, the review evaluates the mean directly on the displaced geometries \(u\pm10^{-4}S_j\). It calls only low-level quadratic graph integrals: it does **not** fabricate `converged` flags or describe these geometric perturbations as newly equilibrated material states. The largest absolute difference between these mean finite differences and the two-channel chain-rule derivative is **1.05e-10 micrometres per declared coordinate**, below the predeclared 2e-6 crosscheck threshold.

Both signs of a geometric displacement direction may be used to check the observable derivative even when the associated material coordinate has only one physical side. This does not evaluate a negative physical apical tension or assert existence of that branch. The apical and basal means are two algebraically distinct channels; no statistical noise independence is assumed. Height mean is their exact sum and is correctly excluded as a third independent row.

## What the SVD does—and does not—say

Recomputing the six-by-six matrices and singular values reproduces both reported diagnostics. Their condition ratios are approximately 7690 for disk weighting and 9261 for line weighting. The weakest reported right-singular vectors have the largest magnitude in the basal-area coordinate, consistent with the originating document's bounded interpretation.

These singular values mix five logarithmic coordinates and one additive coordinate. Their magnitudes and weakest directions depend on the chosen scaling, strengths, mean measure and row convention. They are not a parameterization-invariant uncertainty, a noise-weighted information matrix, the full experimental 68-summary Jacobian, the completed dose cohort or a global identifiability theorem. Full spatial profiles, maxima and dynamics can carry different information. No derivatives of competing maxima have been supplied here.

## One provenance weakness for future reuse

`means()` checks the mechanical source protocol but does not independently enforce the current saved-input and NPZ hashes before composing its outputs. If those files changed between stages, source equality alone would not prevent mixing a new observation gradient with an old mechanical tangent.

**Current snapshot:** all current input, archive, vector and index checks pass in this independent review, so this is not evidence that the retained result is incorrect. Before future reuse, the runner should check each input hash against `mechanics_protocol.json`, each NPZ hash against its mechanics row and the vector/index correspondence before evaluating means. The review records this limited missing runtime guard; it does not modify the frozen source or run a new production stage.

The first independent-report write encountered a NumPy-boolean JSON serialization error. Its original program hash and failure are retained in the final receipt. Casting report flags to ordinary booleans fixed reporting without changing a model, arithmetic method or threshold. The harmless existing Tensor-to-float warning in the production source is not hidden by editing that source.

The completed result is a checked local sensitivity component. Validated roots and forward errors, an empirically justified observation contract, experimental discrepancy accounting and the full parameter-domain analysis remain necessary before material inference or biological validation can be claimed.
