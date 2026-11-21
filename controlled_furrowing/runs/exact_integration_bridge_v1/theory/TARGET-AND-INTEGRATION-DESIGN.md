# Exact-integration target and a feasible derivative-error bridge

**Status: mathematical target/design, not an executed integral certificate.** No new 611-coordinate integration, Hessian assembly, equilibrium or optimization has been run for this document. The objective is a bridge from the already verified binary-cache energy to a coherent, exactly integrated P2 finite-element energy, for all seventeen selected families on their existing shared cell. This is not a continuum/PDE, full physical identified-set, source-observer, loading-path or biological theorem.

## 1. The new mathematical target is explicit

Let the material triangulation consist of the existing triangles with their three stored binary corner coordinates interpreted as exact rationals. Each material triangle uses its exact affine map and the exact P2 Lagrange functions in the established order: three vertices followed by the 01, 12 and 20 edge midpoints. Material interpolation nodes are the mathematical midpoints of those corners. Define the **current** field by

\[
y_h(x;u)=\sum_{j=1}^{6}\phi_j(x)
\bigl(X^{\rm binary}_j+(P u)_j\bigr),
\]

where the stored binary nodal coordinates, global connectivity, free-coordinate selection and fixed current coefficients are retained. The 611-vector \(u\) is the same coordinate vector as in the existing proofs. We do not replace the current nodal coefficients by repaired midpoint positions. Consequently \(u=0\) need not be exactly the identity material map. Equivalently, canonical displacement contains a fixed shift \(X^{\rm binary}-X^{\rm canonical}\), including prescribed shifts at fixed nodes. Calling that shift zero would change the target.

The existing cache is not one exactly coherent basis evaluation: exact binary partition-of-unity and derivative-sum defects were found. Therefore this target is a **new coherent integration target of the same stated constitutive/finite-element model**, not an assertion that the old cached scalar function was already exactly this Galerkin energy. The bridge must bound the difference.

Freeze the existing element mask \(\chi_T\in\{0,1\}\), each exact source dose \(I_i\), and the common five-coordinate parameter cell. The actual mask is constant on each of 144 triangles, with 24 active. Its inspected support aligns with material mesh faces for this model. Defining it as the union of those exact triangles is sufficient; no new biological activation geometry is inferred. All seventeen fields share this geometry, mask and integration-cache structure; their current coordinates and preconditioners differ.

Write \(y_h=(R_h,Z_h)\), \(F_m=\nabla_{r,z}y_h\), \(h=R_h/r\), and \(J=\det(F_m)h\). The bulk density remains

\[
W=\frac\mu2\left(\|F_m\|_F^2+h^2-3-2\log J\right)
+\frac\kappa2(\log J)^2
+\frac{gI_i\chi_T}{2}|F_m e_z|^2.
\]

For a material top/bottom edge parameterized by \(t\in[0,1]\), use its exact affine corner radius \(r_e(t)\) and exact endpoint length \(L_e\). Here the edges are horizontal, so \(L_e\) is an exact rational difference of endpoint radii. With the exact P2 current edge \(y_e(t)=(R_e(t),Z_e(t))\), set

\[
J_s=\frac{R_e(t)|y_e'(t)|}{r_e(t)L_e},
\qquad W_s=\frac{K_e}{2}(J_s-1)^2+T_eJ_s.
\]

The material surface radius is **not** the old rounded P2 reference-radius samples; their difference belongs to the cache-consistency error. Retain the positive exact-binary normalization \(\sigma_c\) from the old cache and define

\[
E_I(u,\theta)=\frac{2\pi}{\sigma_c}
\left[\sum_T\int_T rW\,dr\,dz
+\sum_e\int_0^1 r_e(t)L_eW_s\,dt\right].
\]

Actual mathematical \(2\pi\), not a silently substituted floating approximation, appears here. Cached floating \(2\pi\), areas, weights and radius multiplications are compared to this target. Retaining \(\sigma_c\) makes its normalization mismatch identically zero **by definition**, not by an unverified equality to an exact volume. If exact physical normalization \(\sigma_*\) is also reported, compute the positive factor \(c=\sigma_c/\sigma_*\) rigorously: \(E_*=cE_I\), roots/minimum signs are unchanged, Hessians scale by \(c\), and \(R_*=R/c\) preserves every preconditioned equation. Do not transfer normalized curvature numbers without that factor.

The target has no analytic extension through inverted elements. Its gradient and Hessian are defined only after proving positive meridional determinant, hoop factor, \(J\), and nondegenerate surface tangents on the entire proposed domain. The old cache's clamp and inversion penalty are inactive only on their separately proved old-cache analytic domain. No numerical regularizer is imported into the new physical integral.

## 2. Axis cancellation must precede differentiation

There are six axis-edge triangles, six axis-vertex-only triangles and 132 off-axis triangles in this mesh.

On an **axis-edge** triangle all radial current coefficients on that edge are fixed to zero. The quadratic radial field is divisible by its affine radial coordinate. Divide symbolically in both the field and every allowed variation; the hoop factor is then regular.

An **axis-vertex-only** triangle is different. Vanishing at one vertex does not imply divisibility by \(r\). For example, on vertices \((0,0),(1,0),(1,1)\), the field \(R=r+z\) has \(R/r=1+z/r\) and no single corner limit. Put

\[
\lambda_0=1-s,\quad\lambda_1=s(1-t),\quad\lambda_2=st.
\]

Then \(r=s\bar r(t)\), \(R=s\widetilde R(s,t)\), and the Duffy Jacobian contains \(s\). Since both non-axis vertices have positive radii, \(\bar r(t)>0\). Cancel \(s\) before interval evaluation or AD. On this structured first column the two non-axis radii are equal, so \(\bar r\) is constant and the remaining hoop expression is polynomial. Direction dependence in physical coordinates is not a singularity of this transformed integrand. An axis-edge chart similarly retains its appropriate Jacobian and cancels its common radial factor.

For the first surface edge, \(r_e=L_et\), and the fixed axis radial node makes \(R=t\widetilde R\). Cancel \(t\) in \(R/r_e\), including derivatives with respect to free nodes. Never bound singular \(1/t\) terms separately and hope interval arithmetic will recover their cancellation.

These transformations make each closed integration chart regular. They justify moving first and second free-coordinate derivatives under the integral when the stated positive denominator bounds hold on the compact parameter/tube domain.

## 3. Use the near-unit Jacobian, not separate logarithms

The reviewed old geometry boxes have whole-body \(J\) bounds approximately \([0.95829,1.02341]\). These are useful for the canonical map above, whose nodal interpretation matches the exact body proofs. Exact per-record bounds, not this decimal display, must be consumed. An enlarged proof tube requires fresh domain bounds.

Set

\[
\varphi(J)=-\mu\log J+\frac\kappa2(\log J)^2.
\]

Then

\[
\varphi'(J)=\frac{\kappa\log J-\mu}{J},
\qquad
\varphi''(J)=\frac{\kappa+\mu-\kappa\log J}{J^2}.
\]

For local free-coordinate indices \(a,b\), the volumetric/logarithmic gradient and Hessian contributions are

\[
\varphi'(J)J_a,
\qquad
\varphi''(J)J_aJ_b+\varphi'(J)J_{ab}.
\]

The remaining shear and active contributions are rational-polynomial expressions with known linear strain maps. The derivatives \(J_a,J_{ab}\) are computed from the exact polynomial/rational Jacobian, not from finite differences.

For \(e=J-1\), \(|e|\le\rho<1\), let \(H_n=\sum_{k=1}^n1/k\), \(H_0=0\). Directly combining log and inverse series gives

\[
\varphi'(1+e)=\sum_{n\ge0}(-1)^n(-\mu-\kappa H_n)e^n,
\]

\[
\varphi''(1+e)=\sum_{n\ge0}(-1)^n(n+1)(\mu+\kappa H_{n+1})e^n.
\]

All coefficients are exact rationals when the material coefficients are rational. This avoids large artificial ranges from independently expanding \(\log\det F_m\) and \(\log h\), or from forming an uncorrelated difference between their intervals.

After retaining terms through degree \(N\), valid absolute tails are

\[
B_{1,N}=\rho^{N+1}\left[
\frac{\mu+\kappa H_{N+1}}{1-\rho}
+\frac{\kappa\rho}{(N+2)(1-\rho)^2}\right],
\]

\[
\begin{aligned}
B_{2,N}=\rho^{N+1}\left[
\frac{(N+2)C}{1-\rho}
+\frac{((N+2)B+C)\rho}{(1-\rho)^2}
+\frac{B\rho(1+\rho)}{(1-\rho)^3}\right],\\
C=\mu+\kappa H_{N+2},\qquad B=\frac\kappa{N+3}.
\end{aligned}
\]

The proof uses \(H_{N+1+k}\le H_{N+1}+k/(N+2)\), and the geometric sums for \(1,k,k^2\). For the actual positive \(\mu,\kappa\), these tails bound the **derivative primitives themselves**. A small scalar energy tail alone would not control force or Hessian errors. The local force tail is bounded using \(B_{1,N}|J_a|\); the Hessian tail uses \(B_{2,N}|J_aJ_b|+B_{1,N}|J_{ab}|\), with the integration measure and free-coordinate lift included.

### Surface primitive

Put \(q=J_s^2>0\). Then

\[
W_s=\frac K2 q+(T-K)\sqrt q+\frac K2,
\quad W_s'(q)=\frac K2+\frac{T-K}{2\sqrt q},
\quad W_s''(q)=-\frac{T-K}{4q^{3/2}}.
\]

Apply those derivatives to \(q_a,q_{ab}\) as above. Choose a positive rational center \(q_0\) and prove \(|q/q_0-1|\le\rho_s<1\). The expansions of \((1+e)^{-p}\), for \(p=1/2,3/2\), have coefficients \((-1)^n(p)_n/n!\). Their tails can be bounded by a geometric majorant. For \(p=1/2\), use the ratio bound \(\rho_s\). For \(p=3/2\), after the first omitted term \(n=N+1\), use \(\rho_s(N+1+p)/(N+2)\), provided it is strictly below one. Multiply by the appropriate power of \(q_0\) and by \(|T-K|/2\) or \(|T-K|/4\). A surface with a poor relative range needs degree adaptation or a differently justified integration chart; no successful degree is assumed in advance.

## 4. Radial Laurent moments avoid an expensive inverse-radius expansion

On an off-axis element, the meridional derivatives are polynomials in material \(r,z\), while the hoop factor and Jacobian contain powers of \(1/r\). Thus finite truncations above are finite **Laurent-polynomial** sums \(\sum c_{ab}r^az^b\), with integer radial exponent \(a\) and nonnegative \(b\). Their moments can be evaluated exactly up to rigorously enclosed endpoint logarithms.

For these structured triangles, write the lower and upper edges as affine functions \(z_l(r),z_u(r)\), with \(r_l>0\). Then

\[
\int_T r^az^b\,dr\,dz
=\frac1{b+1}\int_{r_l}^{r_u}
r^a\left[z_u(r)^{b+1}-z_l(r)^{b+1}\right]dr.
\]

Expand the bracket with exact rational coefficients. Each final integral is a difference of integer powers divided by its exponent, or \(\log(r_u/r_l)\) for exponent \(-1\). The latter is enclosed by Arb rather than called rational. Surface moments with denominator \((r_l+L_et)^k\) reduce similarly after an affine substitution.

On axis-containing charts, cancel radial factors **before** expansion. Integrate the resulting ordinary polynomial moments on the unit square, including the Duffy/axisymmetric weight. Do not integrate individually divergent Laurent terms at \(r=0\), even if their formal sum is finite.

This shared moment engine is preferable to expanding an unrelated \(1/r\) series whose convergence is governed by radial element aspect rather than by the small physical variation in \(J\). It can reuse the same geometry moments for every record and derivative entry. Use shifted/scaled local polynomial coordinates to reduce coefficient growth. Large Laurent moments and small coefficients may still cancel; exact rational preprocessing and sufficiently high-precision outward arithmetic must control that loss. Neither the near-unit pointwise range nor increased truncation degree alone guarantees that interval coefficient widths remain useful. The bounded element pilot is required to assess this issue.

## 5. Subtract the quadrature functionals before losing dependence

Let \(Q_c^I\) evaluate the **coherent target** at exact-rational auxiliary nodes with the same positive exact-binary cache weights and normalization. Freeze the following intrinsic choice, which requires no regeneration of floating Gauss abscissae. From each cached bulk shape row in the six-node order, define

\[
\xi=N_1+\tfrac12(N_{01}+N_{12}),\qquad
\eta=N_2+\tfrac12(N_{12}+N_{20}),\qquad
\lambda_0=1-\xi-\eta.
\]

All additions use the exact rationals represented by those binary entries. Prove \(\xi,\eta,\lambda_0>0\) for every row. On a surface row ordered left/middle/right, take \(t=N_{\rm right}+N_{\rm middle}/2\), and prove \(0<t<1\). These identities would reconstruct coordinates exactly from a coherent P2 basis; here they merely define an auxiliary rational rule. They do not assert that these are the actual algebraic Gauss nodes or the original floating generation history. Any such interior-node rule is legitimate in the identity below; proximity to the old nodes affects tightness only. Node witnesses and their selected cache hashes belong in the protocol.

At the force and Hessian levels, split

\[
E_I-E_Q=(E_I-Q_c^I)+(Q_c^I-E_Q).
\]

For each finite spatial polynomial/Laurent model, apply the first difference symbolically:

\[
\mathcal L(P)=\sum_\alpha c_\alpha
\left(M_\alpha-\sum_q w_q^c m_\alpha(x_q)\right).
\]

Only after forming the moment defects take interval absolute values. The coefficients may be interval enclosures over the full parameter/predictor/error tube: their dependence is multiplied by small **moment defects**, rather than by two unrelated large totals. Every binary-rule residual moment is retained, including low-degree and constant moments. Nominal Gauss-8 degree-15 or surface-Gauss-12 degree-23 exactness does **not** make these exact rational residuals zero.

The analytic tails contribute at most the positive integral measure plus the positive discrete mass times the appropriate derivative-tail bound. This is a bound on the error functional, not merely an energy approximation.

The second difference compares the exact target strain maps, radius, reference metric and tangent data at a node with the frozen cache's linear maps. Keep the same \(u,\theta\). Form exact coefficient differences first, for example

\[
q_I-q_c=(B_I-B_c)u+(q_{I,0}-q_{c,0}).
\]

Use local mean-value bounds or a shared derivative Taylor representation to propagate these tiny differences into local force/Hessian differences. Force comparison uses the constitutive Hessian; Hessian comparison also needs controlled third derivatives or an equivalent exact difference identity. Independently interval-evaluating two full Hessians and subtracting them is not the recommended consistency bound. Surface differences include canonical affine \(r_e,L_e\) versus cached P2-radius/length samples. Cached weight and \(2\pi\) rounding already enter the retained moment defects.

For the surface consistency term, an exact algebraic rescaling reuses the existing reviewed local primitive. Let \(\alpha_c=1/(r_{0,c}L_c)\) be its fixed cache scale, and \(\alpha_I=1/(r_eL_e)\) the target scale at that auxiliary point. Evaluate the cache primitive with \(R^*=R_I\alpha_I/\alpha_c\) and the coherent target tangents. Then \(\alpha_cR^*|y_I'|=J_{s,I}\) exactly. The rescaled local map remains linear in the same free coordinates. Reference-metric differences are now included in \(\Delta q,\Delta B\), without treating an independently rounded reference scale as equal.

At each node with unchanged weight \(w\), write \(B^*=B_c+\Delta B\), \(g^*=g(q^*)\), \(H^*=H(q^*)\). Exact identities are

\[
\delta F=w\left[\Delta B^Tg_c+(B^*)^T(g^*-g_c)\right],
\]

\[
\delta H=w\left[\Delta B^TH_cB_c+(B^*)^TH_c\Delta B
+(B^*)^T(H^*-H_c)B^*\right].
\]

Bound the two constitutive differences along the segment joining the two local strain vectors, with all common parameters retained. Third derivatives needed for this comparison are analytic; useful scalar primitives are

\[
\varphi'''(J)=\frac{2\kappa\log J-3\kappa-2\mu}{J^3},
\qquad W_s'''(q)=\frac{3(T-K)}{8q^{5/2}}.
\]

Chain them with the exact polynomial derivatives of \(J\) or \(q\), including their nonzero higher coordinate derivatives. This is a proposed computable consistency bound; no actual coefficient-wise force/Hessian bound has yet been evaluated.

All spatial representations must preserve the parameter and predictor dependence needed by the error functional. Naive interval integration of \(F_I(\Theta,U)\), followed by subtraction of an independently interval-evaluated \(F_Q(\Theta,U)\), produces widths of order \(H\Delta u\) even when the true quadrature error is tiny. Likewise, splitting an element and treating its original global Gauss rule as locally moment-exact is invalid. Prefer one global element Laurent polynomial; if subdivision is used, prove the complete refined-rule-to-original-rule difference or retain a shared parameter Taylor model. Coverage of every spatial cell and all tails remains mandatory.

## 6. Assembly and a decisive persistence theorem

Local force and Hessian defects have at most twelve local displacement coordinates per bulk triangle and six per surface edge. Lift through the frozen global free-coordinate selectors. Directly precondition signed polynomial defect coefficients with the corresponding columns of \(R\) before the final absolute values when that improves the bound. A conservative fallback lifts absolute local component bounds with \(|R|\); its possible pessimism must be reported.

For each record and \(\theta\in C\), let \(p_\theta=u_0+S(\theta-\theta_0)\), and let the old proof supply

\[
\|RF_Q(p_\theta)\|_\infty\le\eta,
\quad\sup_{u\in B_\infty(p_\theta,r)}\|I-RH_Q(u,\theta)\|_\infty\le Z<1.
\]

New integrations must produce

\[
a_p\ge\sup_\theta\|R(F_I-F_Q)(p_\theta,\theta)\|_\infty,
\quad b\ge\sup_{\theta,u\in B}\|R(H_I-H_Q)(u,\theta)\|_\infty.
\]

If

\[
Z+b<1,\qquad \eta+a_p+(Z+b)r<r,
\]

then \(u\mapsto u-RF_I(u,\theta)\) is a strict self-map and contraction of the same cube for every \(\theta\). Banach gives a unique exact-integral root there. The earlier independently verified \(\operatorname{sym}R\succ0\), the symmetric target Hessian and the Lyapunov argument imply

\[
H_I\succeq\frac{1-Z-b}{\|R\|_\infty}I>0.
\]

It is therefore a strict local minimum of the target energy, not merely an algebraic root. If \(e_Q\) bounds the old root's distance to its predictor, the two roots additionally satisfy

\[
\|u_I-u_Q\|_\infty
\le \frac{a_p+b e_Q}{1-Z-b}.
\]

The numerator follows by applying the defect Hessian bound along the segment from the predictor to the old root. A force defect known only at the predictor must not silently be used as a uniform force defect on the full tube.

**The old tiny radii are not promised to pass.** The independent scale audit gives the smallest old self-map margin as approximately \(4.2249\times10^{-12}\) (record G). A sufficient, deliberately nonoptimal split \(a_p\le s/4\), \(b\le s/(4r)\), with \(s=r-\eta-Zr\), preserves half that margin. Using a raw force bound and a separate \(\|R\|_\infty\) factor can require errors near \(10^{-15}\); this is a computational warning, not a physical measurement tolerance.

### If the old cube fails

Retain that failure. Do not shrink the shared physical parameter cell or relabel integration error as biological noise. A separate protocol may enlarge/recenter the **numerical proof domain**, using exact-integral force information as an untrusted center proposal. It must freshly bound both cache and integral Hessians on the new domain, prove geometry/analyticity, and propagate the resulting actual root displacement to all four observations. Translation and new corrector radius remain separate.

If the infinity-norm bounds become too pessimistic, a concrete alternative is a frozen binary square witness \(W\), with \(u=p_\theta+Wz\), \(\|z\|_2\le\rho\). Prove

\[
\eta_E\ge\sup_\theta\|W^TF_Q(p_\theta)\|_2,
\quad a_E\ge\sup_\theta\|W^T(F_I-F_Q)(p_\theta)\|_2,
\]

\[
Z_E\ge\sup_{\theta,z}\|I-W^TH_QW\|_2,
\quad b_E\ge\sup_{\theta,z}\|W^T(H_I-H_Q)W\|_2.
\]

Then \(\eta_E+a_E+(Z_E+b_E)\rho<\rho\), \(Z_E+b_E<1\), gives the same contraction/minimum argument in \(z\). Enclose the ellipsoid by nodal radii \(\rho\|W_{j,:}\|_2\) plus predictor spans for all geometry checks. These new metric bounds must be actually validated; one cannot transfer an infinity-norm \(Z\) by relabeling it. This route has demonstrated relevance in earlier separate model calculations but is not yet proved for these seventeen integration targets.

## 7. Method choice and bounded implementation sequence

**Primary method:** exact spatial Laurent/polynomial moments plus closed, derivative-specific series tails. It exploits near-unit \(J\), keeps low-degree binary-rule errors, shares geometry across records, and exposes where width is lost. Allocate effort from the actual preconditioned root slack, not a guessed decimal tolerance.

**Secondary tool:** rigorous complex quadrature for selected scalar 1D primitives, surface crosschecks or moment verification. The inspected [python-flint 0.9.0 `acb.integral` documentation](https://python-flint.readthedocs.io/en/latest/acb.html) requires forwarding the analytic flag through log/sqrt/branch-sensitive callbacks. It explicitly demonstrates an incorrect result if the flag is ignored. A nested 2D callback must also prove joint analytic validity and enclose inner-integration error for every complex outer ball. Separate real point evaluations are insufficient. Naively invoking a nested scalar integrator for every one of the 611 force and many Hessian entries is not the primary architecture.

**Not recommended as the main proof:** interval Riemann sums, high-order agreement without a remainder, or Gauss–Peano derivative constants applied across activation jumps/axis singular expressions. Riemann bounds converge too slowly for the observed budgets; Peano or complex-analytic bounds require explicit high derivative/analytic-neighborhood control and are reasonable only when those constants are less expensive than the series method.

The bounded execution order is:

1. Freeze the target/normalization, per-record selected root/cache/predictor lineage, coherent node witness and every exact geometry/mask/source binding. No physical source observation is needed for this numerical target.
2. Test radial moments, axis charts, derivative identities, low-degree binary-rule defects and tail formulas on independent analytic fixtures. A scalar energy-only test cannot pass derivative integration.
3. On representative axis-edge, axis-vertex and off-axis elements, and one edge, validate local derivative integrals at a saved center before any 611-job launch. Retain error components separately: analytic tail, arithmetic width, moment/rule rounding and cache-map consistency.
4. On one predeclared record, bound the actual predictor force defect and whole-tube Hessian defect. Compare with its exact old slack; preserve failure. Choose no adaptive radius or changed proof center without a separately frozen attempt policy.
5. Only after a verified nonvacuous one-record bridge, run all seventeen under a common target and retained policies. A successful record is not all-record coverage. Independent replay must implement/check the decisive integral remainder and consistency bounds, not merely read producer status.
6. If roots move, use existing reviewed observer-transfer lemmas with the new numerical radii and fresh physical/ROI domains. The old 68 output bounds and discrepancy enclosure remain historical until this transfer is proved; do not claim unchanged values from an expectation that quadrature error is small.

No actual integral bound, transferred root/minimum, updated observation, computational timing, successful truncation degree, or new original-goal gate completion is claimed by this design. Every failure of a sufficient condition is **UNRESOLVED**, not nonexistence. Even all seventeen successful exact-integration bridges would still leave mesh/PDE error, nonaxisymmetric modes, full physical parameter coverage, source observation alignment, physical paths and larger-region prediction to their original gates.
