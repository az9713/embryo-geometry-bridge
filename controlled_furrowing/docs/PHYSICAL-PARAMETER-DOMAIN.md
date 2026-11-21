# Which mechanical parameters can shape determine?

**Absolute force scale is invisible to this quasistatic shape experiment.** The full nonlinear model therefore needs inference on force ratios, with one fixed physical length unit. The implementation now supplies that mapping and an exact rational force-scale quotient. This is standard dimensional analysis and an identifiability symmetry, not a new concept or a completed biological parameter fit.

The [implementation](../scripts/physical_parameter_domain.py) and [29-check receipt](../evidence/physical_parameter_domain_verification.json) are separate from every frozen physical source and prediction protocol. No biological records or larger-region outcomes enter these checks. The full identified-set propagation required by [the original plan](../PLAN.md) remains incomplete.

## Physical units before fitting

Write the physical constitutive law, on positive-J deformations, as

\[
E[y;I]=\int_\Omega\left[
\frac\mu2\bigl(\lVert F\rVert^2-3-2\log J\bigr)
+\frac\kappa2(\log J)^2
+\frac{gI}{2}\chi_{\rm act}\lVert Fn_0\rVert^2
\right]dV
+\sum_{s=a,b}\int_{\Gamma_s}
\left[\frac{K_s}{2}(J_s-1)^2+T_sJ_s\right]dA.
\]

Here \(\mu,\kappa,g\) have pressure units, with \(g\) understood as pressure per unit of the declared dimensionless fluorescence input. The surface coefficients \(K_s,T_s\) have force-per-length units. The source parameter `bulk` denotes \(\kappa\); the infinitesimal bulk modulus is \(\kappa+2\mu/3\), not \(\kappa\). Geometry, activation support, boundary conditions, observables and the branch-selection rule must also be specified.

Fix one reference length \(L_0=33.5\,\mu\mathrm m\), independent of the specimen height \(H\). Under \(x=L_0\hat x\), the energy divided by \(L_0^3\) has pressure coefficients

\[
c=\left(\mu,\kappa,\frac{K_a}{L_0},\frac{K_b}{L_0},
\frac{T_a}{L_0},\frac{T_b}{L_0},g\right).
\]

The dimensionless height is \(H/L_0\), and the radii and activated depths are divided by the same \(L_0\). Dividing this energy by \(\mu\) gives the solver's shear coefficient one and surface coefficients \(K_s/(\mu L_0),T_s/(\mu L_0)\). Thus changing \(H\), at fixed physical material, changes the geometry but leaves these fixed-unit material coefficients unchanged.

If a different implementation uses \(H\) itself as its length unit, it must instead update the surface ratios to \(K_s/(\mu H),T_s/(\mu H)\). Keeping these ratios fixed while changing \(H\) changes the physical surface material. The two conventions are equivalent only when all coefficients and geometry are converted consistently. This mapping does not accuse earlier fixed-height computations of a unit error; it prevents an error in the proposed height-uncertainty propagation.

The helper accepts pressures in Pa and surface coefficients in Pa·µm, so its volume-energy unit is Pa·µm³. A conversion to joules requires the factor \(10^{-18}\); no absolute biological energy is inferred here.

## Exact nonlinear force-scale symmetry

**Proposition.** Hold geometry, input \(I\), admissible deformations and the fixed-displacement boundary conditions constant. Assume there is no separately fixed applied force, pressure or body load. For every \(\alpha>0\), replacing all seven coefficients by \(\alpha c\) leaves the stationary deformations, their local-minimum character and the set of global minimizers unchanged.

**Proof.** Each term of the physical energy is linear in its coefficient, so

\[
E_{\alpha c}[y;I]=\alpha E_c[y;I],\quad
DE_{\alpha c}=\alpha DE_c,\quad
D^2E_{\alpha c}=\alpha D^2E_c.
\]

Stationarity is unchanged, and multiplication by a positive scalar preserves the ordering of energies and the sign of every second variation. The same holds for the stationary branch as the input varies, provided the branch-selection convention itself is invariant under positive energy scaling. Consequently every shape observable on that branch is unchanged. This proves a structural nonidentifiability of absolute force scale even with noiseless full-shape measurements.

The statement is not uniqueness of equilibrium. It also does not imply identical physical relaxation times with fixed mobility or drag. Independently known traction, pressure, force or dynamical calibration can break the symmetry. A measured activation fluorescence is not by itself a calibrated mechanical force.

The optimizer has a fixed-coefficient inversion penalty outside its accepted physical regime. That numerical extension is not part of this theorem. The retained kernel checks use positive-J fields with zero trial penalty; they verify both raw energy scaling and invariance of the objective and gradient normalized by \(\mu V_0\). No symmetry of invalid trial fields is claimed.

## An unbounded ratio domain without an arbitrary finite box

For \(\mu,\kappa>0\) and the other five coefficients nonnegative, define

\[
w_j=\frac{c_j}{\sum_{k=1}^7c_k}.
\]

Then \(w\) lies in the simplex \(\sum_j w_j=1\), with \(w_\mu,w_\kappa>0\). Conversely, \(c_j/c_\mu=w_j/w_\mu\); thus these normalized vectors represent exactly the positive common-scale equivalence classes. The exact implementation uses rational arithmetic to verify normalization, invariance and the inverse ratios.

This avoids pretending that a few material scenarios, or an arbitrary bounded ratio grid, exhaust the physical uncertainty. It does **not** make the admitted parameter domain compact: its closure is compact, but contains excluded faces such as \(\mu=0\) or \(\kappa=0\). Ratios can diverge near those faces. Existence, stability, continuity, attained prediction extrema and uniform numerical accuracy do not follow merely from this normalization. Those properties need separate hypotheses or proofs.

The rational representation covers rational inputs without floating-point underflow. The convenience floating-point quotient and solver conversion cover only representable ratios: they explicitly reject underflow that would replace a positive coefficient by zero. This guard was added after independent review found that an extreme finite input could otherwise be mapped onto an excluded simplex face. No finite-precision helper is claimed to cover the entire unbounded ratio domain.

If \(\kappa\) is selected as a numerical approximation to incompressibility rather than inferred as a physical compressibility coefficient, it is a method parameter. For a declared ratio \(r=\kappa/\mu\), impose \(w_\kappa=r w_\mu\) and study numerical/penalty convergence separately. Do not tune it to held-out folding outcomes or count numerical regularization as measured material information.

## The actual inference object still needed

Let \(\gamma\) denote independently justified geometry and observation assumptions, \(\mathcal O_k(w,\gamma,I_i)\) the model's kth observable on the declared branch or quasistatic peak path, and \(y_{ik}\) the small-region observations. With fixed measured doses and a declared vector of nonnegative discrepancy allowances \(\delta_k\), the conditional identified set is

\[
\Theta_\delta=\{(w,\gamma): w\in\mathcal W,\ \gamma\in\mathcal G,
\ |\mathcal O_k(w,\gamma,I_i)-y_{ik}|\le\delta_k
\ \text{for all available }i,k\}.
\]

The separate component maxima must remain separate observables; their maxima may occur at different times or positions and must not be summed. If input errors are allowed, their independently specified intervals enter this set explicitly. If specimen-specific geometry varies, its common constraints and specimen-specific nuisance variables must be declared before optimization; adding arbitrary per-embryo mechanics is not a free explanation of every residual.

There is no measured discrepancy law here that selects a unique \(\delta\). The [sharp monotone residual floor](CALIBRATION-UNCERTAINTY.md) is a necessary compatibility condition under its common-path assumptions. It is not an estimate of the error allowances and cannot supply them by fiat. At present the scientifically honest object is a **conditional family of identified sets**, not an invented confidence region.

The desired larger-region prediction envelope is

\[
\mathcal P_\delta(e)=\{\mathcal Q_e(w,\gamma):(w,\gamma)\in\Theta_\delta\}.
\]

This new module neither computes \(\Theta_\delta\) nor \(\mathcal P_\delta\). It makes their units, scale symmetry and domain explicit so that an actual inference implementation can be judged against the original requirement. A validated pair of compatible models can establish ambiguity; certifying agreement requires covering the whole declared set. An empty identified set requires rejection or revision of the assumptions, not reporting a best optimizer as compatible.

## Retained verification

The 29 checks include exact rational scale invariance, recovered coefficient ratios, dimensional conversion, a changed-height test holding physical material fixed, actual high-quadrature energy/gradient scaling on a nonzero admissible P2 deformation, and eleven invalid-input cases including extreme-ratio underflow. A receipt links every checked physical-kernel source by SHA-256. This is implementation verification and a standard symmetry proof, not independent experimental validation or novelty evidence.
