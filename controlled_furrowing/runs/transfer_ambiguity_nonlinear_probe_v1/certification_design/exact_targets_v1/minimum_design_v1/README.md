# Strict local minimum after a same-target exact-spatial root

Design only. The four transfer targets currently have separate finite-cache root and geometry evidence. This memo specifies the additional inputs for a strict local minimum of the exact spatially integrated finite-element energy. No new matrix, Hessian, integration, equilibrium, or curvature calculation ran.

## What the root proof supplies

For one preset material tuple and exact rational gain–dose strength, write \(E_I\) for the exact spatially integrated finite-element energy and \(E_Q\) for the exact-real energy of the frozen binary quadrature cache. Both use the same free coordinate order, the same sign convention \(F_I=\nabla E_I\), the same positive binary energy normalization, and the same exact deformation \(Y=X+\operatorname{inj}(u)\). Denote \(H_I=\nabla^2E_I\) and \(H_Q=\nabla^2E_Q\). These are real symmetric Hessians on the admissible 611-dimensional free coordinate domain.

The bridge must provide a strict self-map test for the same admitted cube \(\mathcal B_r\) with \(r=2^{-34}\), \(z<1\) and \(\eta+zr<r\). Its signed force is evaluated at the exact cube center. The separate whole-cube evidence bounds
\[
\|I-RH_Q(u)\|_\infty\le Z_Q,\qquad
\|R[H_I(u)-H_Q(u)]\|_\infty\le b,\qquad
z=Z_Q+b<1 .
\]
Here \(R\) is the actual exact binary-rational preconditioner. It may be nonsymmetric. The \(b\) bound comes from complete signed upper-Hessian defect packets on \(\mathcal B_r\), assembled symmetrically and contracted with that same \(R\) before absolute row sums. The normalized cache/auxiliary identity must cover the full operator \(H_I-H_Q\), including the finite same-deformation cache correction. A cache-only derivative, center-only derivative, argument with rounded drive, or operator evaluated in a different set of free coordinates cannot supply these premises.

The defect bound implies invertibility, but it contains no inertia information by itself. For example, \(R=H=\operatorname{diag}(-1,1)\) gives \(I-RH=0\); the critical point of the quadratic energy with Hessian \(H\) is a saddle. Likewise symmetry of \(R\) alone is insufficient.

## Preferred route with one additional matrix certificate

Certify the positive definiteness of the exact symmetric part
\[
P=\tfrac12(R+R^\mathsf T).
\]
Build \(P\) by exact binary-rational addition or outward interval addition of exact binary entries. Rounded floating symmetrization can propose a certificate but cannot define the certified \(P\). For this route only \(R\) needs a new matrix certificate. No new Hessian integration or quadrature assembly is required once the existing admitted \(z\) covers the right cube.

A reviewed local routine already supports this operation. In `scripts/cohort_energy_stability.py:76`, `certify_sympart(R,K=None,precision_bits=...)` accepts an explicit finite binary64 matrix and can propose an inverse Cholesky factor \(K\) from floating arithmetic. Freeze the actual binary \(K\) and accept only an outward bound
\[
q\ge\|I-KPK^\mathsf T\|_\infty,\quad q<1.
\]
The exact congruence matrix is symmetric with eigenvalues strictly between 0 and 2. It is therefore SPD; since \(K\) is square, the congruence is nonsingular only if \(K\) is invertible, which proves \(P\succ0\). Use `scripts/review_cohort_energy_stability.py:62,67`, `exact_P` and `congruence_bounds(R,K,precision)`, for an independent exact-binary reconstruction and a different multiplication association. Independently recompute the exact dyadic row-sum norm \(\|R\|_\infty\) with `norm_R_integer` at line 51 or use an independently certified upper norm for the identical \(R\). The historical drivers around these pure routines select original-cohort records, so a new explicit singleton caller and prelaunch review are required. The routines set Arb precision; isolate them in their own process and pin its precision and sources. Floating Cholesky failure or a bound \(q\ge1\) is unresolved and supplies no negative inertia claim. For example, take \(R=\bigl(\begin{smallmatrix}1&4\\0&1\end{smallmatrix}\bigr)\) and \(H=\operatorname{diag}(1,1/10)\). The symmetric part of \(R\) is indefinite although \(H\) is SPD and \(\|I-RH\|_\infty=9/10\).

**Why this suffices for the exact spatial Hessian.** Fix any \(u\in\mathcal B_r\) and write \(H=H_I(u)\), \(A=RH\). The bound \(\|I-A\|_\infty\le z<1\) puts every eigenvalue of \(A\) in the open right half-plane and makes \(R,H\) nonsingular. Because \(H=H^\mathsf T\),
\[
A^\mathsf T H+HA=2HPH\succ0.
\]
The equality follows by expanding \(A=RH\). Positive stability makes \(e^{-tA}\to0\). Differentiate \(e^{-tA^\mathsf T}He^{-tA}\), integrate on \(t\ge0\), and obtain
\[
H=\int_0^\infty e^{-tA^\mathsf T}(2HPH)e^{-tA}\,dt\succ0 .
\]
Thus the entire admitted exact-spatial cube is SPD, even when \(R\) is nonsymmetric. This is the same symmetric-part and Lyapunov route independently reviewed for the earlier cohort.

## Alternative with an independently certified SPD reference

If direct symmetric-part positivity does not pass, use an exact symmetric SPD reference \(H_{\rm ref}\) in the same coordinates with a verified \(q_{\rm ref}\ge\|I-RH_{\rm ref}\|_\infty<1\). For every \(u\in\mathcal B_r\), form \(H_t=(1-t)H_{\rm ref}+tH_I(u)\). It is symmetric and
\[
\|I-RH_t\|_\infty\le (1-t)q_{\rm ref}+tz<1.
\]
All \(H_t\) are nonsingular. Their ordered real eigenvalues vary continuously, so no eigenvalue can cross zero. The inertia of the SPD endpoint transfers to \(H_I(u)\). A nonsymmetric comparison such as an unchecked \(R^{-1}\) would invalidate this inertia step.

Three possible ways to obtain such an endpoint are:

- **Same-target exact cached center.** Certify \(H_Q(u_0)\succ0\) by verified factorization or exact congruence. Bind the already reviewed \(Z_0\ge\|I-RH_Q(u_0)\|_\infty<1\) to this exact center and cache. It is sufficient as a comparison endpoint even though this center-only bound supplies no whole-cube derivative control. Existing cached pilots certify a root and a derivative bound. Those facts alone do not certify \(H_Q(u_0)\succ0\).
- **Same-target cached whole cube.** If a fresh independent proof establishes \(H_Q(u)\succ0\) for all \(u\) on this larger cube, the pointwise path \(H_Q(u)+t(H_I(u)-H_Q(u))\) has defect bounded by \(Z_Q+tb\le z<1\). A cached SPD certificate on a smaller cube would not support this pointwise path at all \(u\); its correctly bound center can still serve as the first option.
- **Explicit SPD comparison matrix from an existing proposal.** Reuse the archived sparse binary Hessian proposal as an exact comparison after forming its exact symmetric part in the verified free order. Certify that comparison's positivity and its own residual \(q_{\rm ref}\) against the same binary \(R\). Its original evaluation used rounded drive. It is only a comparison matrix, whose separate exact certificates are what make the homotopy valid. Neither a floating inertia calculation nor the exact-drive \(Z_0\) for a different matrix can certify this comparison.

The older B implementation compares a previously certified cached base \(H_0\) with \(H_I\), using \(Z_0<1\) and the spatial bound. Its new caller, not its target data, is the template. The older L implementation uses the same-record cached SPD premise and symmetric pointwise path. See `scripts/exact_integral_root_B_safe.py:55,91,122,142`, `scripts/review_cohort_exact_record.py:435`, and their independent actual readbacks. **The old cohort L certificate does not transfer to either new small L material candidate.** Their exact material tuples, drives, caches, centers, preconditioners and radius authorities differ even when the physical source record ID agrees.

## Curvature and strict minimum

Once positivity is established, use \(H_I^{-1}=(RH_I)^{-1}R\) and the Neumann bound to obtain, pointwise and uniformly,
\[
\|H_I^{-1}\|_\infty\le \frac{N_R}{1-z},
\qquad
\lambda_{\min}(H_I)\ge m:=\frac{1-z}{N_R}>0,\quad
N_R\ge\|R\|_\infty .
\]
The second inequality follows by bounding the spectral radius of the symmetric SPD inverse by its induced infinity norm. It does not require \(R\) itself to be symmetric or to represent a true inverse. \(m\) is a lower bound in the normalized energy and the existing unscaled free-coordinate basis, not a physical modulus. If a new certificate reports a differently scaled energy, explicitly transform the bound instead of equating the two normalizations.

The strict map inequality gives \(e=\eta/(1-z)<r\), putting the proved stationary zero strictly inside the cube. Analytic domain evidence must show the underlying energies are \(C^2\) on a neighborhood of the cube in the admissible free space. This follows for the stated target from the reviewed positive bulk \(J\), regularized axis chart, strictly positive trace Jacobian/tangent domains, fixed material and mask, compatible exact derivative maps and uniform margins supporting differentiation under the integral. Body embedding and fixed constraints ensure these neighboring FE deformations have the stated physical domain. These properties do not substitute for the separate SPD result. For any \(v\) in the same cube the segment from \(u_*\) to \(v\) stays in the cube and the second-order integral formula gives
\[
E_I(v)-E_I(u_*)\ge \tfrac12 m\|v-u_*\|_2^2 .
\]
This proves a strict local minimum among the admissible free FE variations and establishes strict separation of the energy from that of any other point in the certified cube. A point-SPD Hessian at the interior root plus verified \(C^2\) would suffice for a strict local minimum by continuity, but would not justify the displayed uniform curvature claim without a surrounding bound.

## Required same-target input inventory

1. **Target and domain.** Unique pilot identity and singleton parameter cell, exact prescribed rational drive and frozen binary material coefficients, cache and embedded payload hashes, exact center \(X+\operatorname{inj}(u)\), free-map and matrix shape/order, precise binary preconditioner hash, energy units/normalization, current independently reviewed body/trace domains and the same closed \(2^{-34}\) cube. Retain both exact small and larger geometry review bindings, while using the larger one only for this derivative.
2. **Cache derivative and energy identity.** Native same-exact-drive \(Z_Q(2^{-34})\) and its checked site and domain witnesses, reviewed full terminal provenance, and the previously reviewed signed \(RF_Q(u_0)\) if needed to bind the root result. Preserve the distinction between \(Z_0\), \(Z_Q(2^{-36})\) and \(Z_Q(2^{-34})\). The old cached pilots explicitly left energetic stability open.
3. **Spatial derivative and root.** Admit complete 144 bulk and 24 trace packets, terminal worker evidence, independent local F/H containment, the exact singleton and point/whole-cube scopes, complete projected upper-Hessian occupation without duplicate entries, and independent global reductions reproducing \(b,z,\eta,e\) and the strict existence result. Record all source/phase and endpoint array hashes, including the signed \(RF_I\) and row sums. No minimum can be promoted from a raw packet, a plan, or an unreviewed contraction candidate. The existing reduction caller and pure reduction kernel intentionally make no minimum claim.
4. **Inertia certificate and exact scalar review.** Choose either the new \(P=\mathrm{Sym}(R)\) congruence or a verified SPD endpoint with its valid residual. Record actual \(R\) and, if applicable, binary \(K\) or \(H_{\rm ref}\) byte hashes, exact construction and all interval row bounds, proof precision, current producer terminal and independent replay receipts. Recompute \(m>0\) using an admitted exact norm or outward upper norm. Confirm the interior root, actual \(C^2\) energy/geometry hypotheses, correct Eulerian spatial bridge target, and the precise scope in the final minimum receipt. Run no actual proof until a standalone source/protocol/preflight and independent launch gate are approved.

The resulting minimum, if certified, would be for one selected fixed input and exact-spatial FE neighborhood. Stable loading would require a verified continuous root branch through the loading parameter with positive, matched derivative and inertia bounds over a connected domain. Time evolution, dynamical stability, global minima, continuum/PDE guarantees, source compatibility, and other calibration states require separate evidence.
