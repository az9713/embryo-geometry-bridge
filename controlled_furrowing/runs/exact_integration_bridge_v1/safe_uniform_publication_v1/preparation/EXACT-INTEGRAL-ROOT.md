# From a cached equilibrium to an exactly integrated finite-element equilibrium

**The scalar composition passed its final independent 1,500-check review. A second direct uniform-residual route also passed.** The result below concerns one source-dose record B, every parameter in the same small closed cell C, and a fixed 611-coordinate axisymmetric model. It does not yet concern all 17 records or the continuum elasticity PDE.

The broader project asks whether a mechanically consistent model can predict tissue geometry across perturbations. Localized light-induced contraction in embryonic tissue motivates the model; apical and basal mean its upper and lower surfaces. This chapter establishes a numerical mathematics prerequisite: an equilibrium of the exactly integrated finite-element energy, rather than only an equilibrium of a cached numerical approximation. [The biological source study](https://doi.org/10.1038/s41467-025-62483-6) and the [source-contract audit](PUBLIC-SUPPLEMENT-AUDIT.md) remain separate from this mathematical result.

## Three distinct targets

The model is a finite-height, axisymmetric solid with quadratic meridional displacement coefficients, fixed remote rim, and exactly constrained radial motion on the axis. For deformation gradient \(F\) and \(J=\det F\), its normalized bulk density is

\[
W=\tfrac12\bigl(\operatorname{tr}(F^{\mathsf T}F)-3-2\log J\bigr)
+\tfrac{5000}{2}(\log J)^2
+\tfrac{g I_B\chi}{2}|F e_z|^2.
\]

Here \(\chi\) is the unchanged elementwise illuminated mask and the exact source dose is \(I_B=0.14801043\). Upper/lower surface densities are \(K_s(J_s-1)^2/2+T_sJ_s\). The five varying coefficients are \(\theta=(K_a,K_b,T_a,T_b,g)\). Their original tiny closed cell C is retained, including its one-sided face at zero apical tension; it is not a full physical identified set. [Exact parameter bounds](../runs/exact_integration_bridge_v1/integral_root_B_safe_v2/protocol.json) are retained as fractions.

* \(E_Q\) is the scalar energy defined by the actual frozen binary cache: its stored strain maps, nodes, weights and normalization. Its gradient is \(F_Q\).
* \(E_I\) integrates the coherent P2 field exactly over the same affine-reference material triangles and surfaces. Its gradient is \(F_I\). The full current nodal coefficients and stored reference-midpoint offsets are retained; nothing is silently repaired.
* The continuum energy allows functions outside this fixed finite-element space. Neither existence/error transfer to its PDE nor nonaxisymmetric stability is established here.

Both discrete energies use the same positive normalization and physical constraints. Proving a root of \(F_Q\) therefore does not automatically prove a root of \(F_I\).

## The claim and its domains

Let \(p_0(\theta)=u_0+S(\theta-\theta_0)\) be the fixed affine predictor, with the stored binary coefficients interpreted exactly. Put

\[
r=10^{-11},\qquad p_1(\theta)=p_0(\theta)-d,
\qquad r_c=r+\|d\|_\infty.
\]

The decrement d is the exact difference of two retained binary candidate vectors, not an assumed equilibrium. Its error relative to the preconditioned integration defect is bounded by m, defined below. No nonlinear equilibrium solve or material refit supplies the result.

For each theta, the common cube
\[
X_\theta=p_0(\theta)+B_\infty(r_c)
\]
contains the old predictor cube, translated candidate cube and their connecting segments. A larger, explicit independent nodal box X contains every such cube and the segments from u0 to each predictor. This distinction matters: the cached derivative bound is used on the coupled cubes \(X_\theta\), while integration-defect and parameter-variation bounds cover the independent domain \(X\times C\).

The successful sufficient tests establish, uniformly for \(\theta\in C\), a unique integral root in the translated radius-r cube, a strict local energy minimum, and no integral root in the old radius-r cube. Uniqueness here extends only to the common cube at fixed theta, not to every configuration or parameter.

## Residual transfer, including parameter variation

Use the unchanged square preconditioner R and define

\[
D(u,\theta)=R(F_I-F_Q)(u,\theta),\qquad
\|R(H_I-H_Q)\|_\infty\le b\quad\text{on }X\times C.
\]

The existing proofs give \(\|RF_Q(p_0,\theta)\|_\infty\le\eta_Q\), cached derivative bound \(\|I-RH_Q\|_\infty\le z_Q\) on the common cubes, and

\[
\|D(u_0,\theta_0)-d\|_\infty\le m.
\]

Because the energy coefficients enter affinely at fixed u, separately bounded base-field parameter derivatives supply a coefficient-only variation bound A. With
\(s\ge\sup_{\theta\in C}\|S(\theta-\theta_0)\|_\infty\), the fundamental theorem of calculus gives

\[
\|D(p_0(\theta),\theta)-D(u_0,\theta_0)\|_\infty
\le V:=bs+A.
\]

This term cannot be dropped merely because C is small. Applying the same theorem along \(p_0-td\) gives the exact identity

\[
\begin{aligned}
RF_I(p_1,\theta)
={}&RF_Q(p_0,\theta)+[D(u_0,\theta_0)-d]\\
&+[D(p_0,\theta)-D(u_0,\theta_0)]
+\int_0^1[I-RH_I(p_0-td,\theta)]d\,dt.
\end{aligned}
\]

Consequently,

\[
z=z_Q+b,\qquad
\eta_I=\eta_Q+m+V+z\|d\|_\infty.
\]

No favorable cancellation between separately enclosed errors is assumed.

## Why the first representation failed, and what changed

The first safe uniform assembly had z below one but \((\eta_I+zr)/r\approx1.230797\). Its sufficient self-map test was unresolved. This did not prove that the integral energy lacked an equilibrium.

Six elements have a full edge on the rotation axis. Write \(R_{\rm ref},Z_{\rm ref}\) for material coordinates and \(q=(r_{R_{\rm ref}},r_{Z_{\rm ref}},z_{R_{\rm ref}},z_{Z_{\rm ref}},r/R_{\rm ref})\) for the five strain variables. Their three constrained radial nodal values vanish exactly, so the quadratic radial polynomial is divisible by \(R_{\rm ref}\). Dividing the exact rational basis before intervalization makes all five strain variables affine physical polynomials and J cubic; the Duffy expression has degree six. Exact coefficient identities verify the same target and free derivatives. The revised assembly replaces those six enclosures, retaining the other 138 bulk contributions and all 24 surfaces. [Representation proof](../runs/exact_integration_bridge_v1/axis_edge_physical_preparation/PROOF.md), [six-element execution](../runs/exact_integration_bridge_v1/uniform_axis_edge_B_safe_v2/RESULTS.md).

The arithmetic also required a separate repair: the earlier unsafe sparse cleanup dropped nonexact Arb coefficients overlapping zero. Those affected point claims were withdrawn. The present chain uses the separately tested exact-zero-only kernel and safe replays; agreement between two implementations sharing the former bug was insufficient verification. [Safe kernel audit](../evidence/exact_safe_polynomial_kernel_independent_review.json).

![Strict contraction and self-map margins](../runs/exact_integration_bridge_v1/safe_uniform_publication_v1/preparation/root_margins.png)

## Existence and strict local minimality

For \(T_I(u)=u-RF_I(u,\theta)\), the successful exact-rational tests are

\[
z\approx0.1050417212<1,\qquad
\eta_I+zr\approx3.337978164\times10^{-12}<10^{-11}.
\]

Banach's theorem yields a unique fixed point in the translated cube, with

\[
\|u_I-p_1\|_\infty\le e_I:=\frac{\eta_I}{1-z}
\approx2.556053177\times10^{-12}.
\]

The cached base bound \(\|I-RH_0\|_\infty<1\) makes R nonsingular, so a fixed point is an actual zero of FI.

An independently proved symmetric positive-definite cached base Hessian H0 supplies the minimum argument. For any integral Hessian H on a common cube, the symmetric matrix homotopy \(H_t=(1-t)H_0+tH\) satisfies

\[
\|I-RH_t\|_\infty\le(1-t)z_0+tz<1.
\]

No eigenvalue can cross zero, hence H is positive definite. This does not require symmetric R. Neumann inversion and the spectral-radius bound give

\[
\lambda_{\min}(H)\ge\frac{1-z}{\|R\|_\infty}
\approx5.305244640\times10^{-4}>0.
\]

Thus energy rises at least quadratically near the root while the connecting segment stays in the common cube. This is curvature of the normalized discrete energy in its declared coordinates, not a physical frequency, dynamical-stability theorem or global minimum.

A second [direct uniform-force residual route](../evidence/exact_uniform_direct_root_B_crosscheck.json), avoiding the base-point-defect and parameter-derivative premises, gives

\[
\eta_{\rm dir}\approx2.299666995\times10^{-12},\qquad
\eta_{\rm dir}+zr\approx3.350084207\times10^{-12}<r.
\]

It yields a slightly wider error bound, approximately \(2.569580113\times10^{-12}\), with the same contraction constant. This is a distinct successful residual composition sharing the validated uniform force/Hessian kernels; it is not a wholly independent physics implementation. The sharper primary bounds above remain the ones used below.

## Root location and the cached/integral branch gap

The retained decrement in free coordinate 588 is positive. The existence/error theorem implies

\[
|u_{I,588}-p_{0,588}|\ge d_{588}-e_I
\approx1.418087117\times10^{-11}>r.
\]

The integral root is outside the old cube. Moreover, two integral zeros in the same common convex cube would obey \(\|x-y\|_\infty\le z\|x-y\|_\infty\), so must coincide. Since the newly established root is inside that common cube, there is **no integral root in the old cube**. This conclusion needs existence plus common-cube uniqueness; the older point-map escape alone could not establish it.

The already proved cached root has \(\|u_Q-p_0\|_\infty\le e_Q\approx5.908037852\times10^{-13}\). Therefore,

\[
|u_{I,588}-u_{Q,588}|\ge d_{588}-e_I-e_Q
\approx1.359006738\times10^{-11}>0.
\]

![Coordinate-588 proof enclosures](../runs/exact_integration_bridge_v1/safe_uniform_publication_v1/preparation/root_location.png)

The diamonds mark predictors, not equilibria. The displayed intervals are proof enclosures, not observations. This positive coordinate gap demonstrates that the two mathematical targets have distinguishable nearby stationary branches at the certificate scale; it does not show experimentally detectable tissue displacement.

## Receipts and unresolved requirements

The [producer scalar composition](../runs/exact_integration_bridge_v1/integral_root_B_safe_v2/results.json) retains exact fractions, while [its protocol](../runs/exact_integration_bridge_v1/integral_root_B_safe_v2/protocol.json) binds the premises. Supporting receipts include the [safe point reduction](../evidence/exact_complete_point_B_independent_review_safe_v2.json), [fresh common-body geometry](../evidence/exact_retube_B_geometry_independent_review.json), [cached common-cube reduction](../evidence/exact_cached_Q_common_B_reduction_review.json), [complete physical-axis uniform reduction](../runs/exact_integration_bridge_v1/uniform_full_B_physical_axis_independent_reduction_v2/results.json), [parameter variation](../evidence/exact_parameter_defect_B_independent_reduction.json), and [base energy-stability review](../evidence/cohort_energy_stability_independent_review.json). The [final independent scalar-composition review](../evidence/exact_integral_root_B_safe_independent_review.json) passed 1,500 checks. Independent project code/reduction checks are not external peer review.

What remains: integral certificates for the other 16 records; a full physically justified material/geometry identified set; stable common loading paths and dynamics; finite-element/PDE error and nonaxisymmetric perturbations; source-matched reference, ROI and temporal observation operators with defensible data uncertainty; revised larger-region predictions; and a claim-level originality audit. This strengthens the mechanics-to-geometry mathematics prerequisite. It does not establish a chemistry-to-biology mechanism, a general explanation of embryonic shape, confirmed novelty, or completion of the original five gates.
