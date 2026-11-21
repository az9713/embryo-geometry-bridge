# Operator assembly and local truncation policy

This is a design addendum to the frozen target/theory receipt. The surface prelaunch and independent six-edge replay are recorded separately. Nothing here changes a frozen integration policy, reruns an actual field, proves a new root, or transfers a point calculation to a uniform parameter/tube claim.

## 1. The vector/Hessian path is explicit

`exact_fe_target.element` supplies five scalar strain polynomials \(q=(a,b,c,d,h)\), their state-independent local free-coordinate map \(B\), and global indices. Let

\[
D=ad-bc,\quad J=Dh,\quad
D_i=B_{ai}d+aB_{di}-B_{bi}c-bB_{ci},
\]

\[
D_{ij}=B_{ai}B_{dj}+B_{aj}B_{di}-B_{bi}B_{cj}-B_{bj}B_{ci},
\quad J_i=D_i h+DB_{hi},
\quad J_{ij}=D_{ij}h+D_iB_{hj}+D_jB_{hi}.
\]

For the scalar derivative series \(A_N\approx\varphi'(J)\), \(B_N\approx\varphi''(J)\), form the local force and Hessian polynomials

\[
G_i=\mu\sum_kq_kB_{ki}+f(bB_{bi}+dB_{di})+A_NJ_i,
\]

\[
H_{ij}=\mu\sum_kB_{ki}B_{kj}+f(B_{bi}B_{bj}+B_{di}B_{dj})
+A_NJ_{ij}+B_NJ_iJ_j.
\]

Cache powers of \(J-1\) once per element. Do not expand the full high-degree product separately for every coordinate pair. The already reviewed `polynomial_jet.contract_volumetric_defect` implements a more efficient moment dual. For a high-degree polynomial \(A=\sum_\alpha A_\alpha x^\alpha\), precompute

\[
d_A(\beta)=\sum_\alpha A_\alpha\,\Delta M(\alpha+\beta),
\quad \mathcal L(Af)=\sum_\beta f_\beta d_A(\beta).
\]

Only the supports of the low-degree factors \(J_i,J_{ij},J_iJ_j\) are needed. Reuse their duals across all local force/Hessian entries. The equivalent surface path uses \(G=J_s^2\), \(G_i,G_{ij}\) and the two surface derivative series. None of these truncated polynomials alone is the complete error bound: analytic tails and finite-cache consistency are mandatory.

The moment measure is an easy place to introduce a factor error. `normalized_bulk_moment` already includes normalized axisymmetric \(r\,dr\,dz\). `axis_chart_moment` already includes \(r\) and the affine/Duffy Jacobian. `surface_moment` includes \(r\,dr\), with \(2\pi/\sigma_c\) supplied separately. Do not multiply those moments by `fields.measure_polynomial` again. If ordinary unweighted square moments are used instead, explicitly multiply the measure there, exactly once.

Geometry/rule moments can be shared across the seventeen fields only with a complete key: exact corners/chart permutation, axial shift, normalization, auxiliary points and binary weights. An element number alone is not a safe key for a future different geometry. All nonzero constant/low-degree moment defects remain in every contraction.

## 2. A local Jacobian range can reduce order without weakening the proof

For an off-axis element, multiply the Laurent polynomial \(e=J-1\) by \(r^d\), where \(d\) removes every negative radial exponent. In local coordinates

\[
r=r_c+h_rx,\quad \zeta=\zeta_c+h_zy,\quad |x|,|y|\le1,
\]

expand its numerator as \(\sum c_{ab}x^ay^b\). Then the computable bound

\[
\sup_T|J-1|\le
\frac{\sum_{a,b}|c_{ab}|}{r_{\min}^d}
\]

is valid. It uses a containing rectangle and can be conservative. Exact Bernstein bounds on the actual triangle may improve it. On regularized axis charts use ordinary polynomial bounds; never divide by an axis minimum of zero. A separately validated whole-body range provides another bound. Taking the minimum of valid upper bounds is safe.

**The coefficient domain determines the scope.** Point coefficients give only a point-field range. Uniform predictor coefficients must include every \(\theta\in C\). Full-tube coefficients must also include the entire numerical corrector domain. A small point-field outer-element range cannot select orders for a uniform root-tube proof unless it is re-established on that full domain. Record the parameter cell, proof-domain identifier and coefficient provenance beside each \(\rho_T\).

Require an outward rational \(0\le\rho_T<1\). Floating values can propose a degree, but only exact rational/outward-ball inequalities accept it. A locally poor bound can fall back to the validated global range. Subdivision or a different expansion center must retain complete coverage and the original-rule difference. Increasing spatial subdivisions does not make the old binary rule locally Gaussian-exact.

## 3. Choose degree from operator budgets

Let \(M_T\) be the positive normalized material mass plus the positive auxiliary-rule mass. Let \(R_T\) denote the columns of the fixed preconditioner corresponding to the element's local coordinates. Valid preconditioned tail prefactors are

\[
\begin{aligned}
A_T&=M_T\sup\|R_TJ_u\|_\infty,\\
B_T&=M_T\sup\|R_TJ_u\|_\infty\,\sup\|J_u\|_1,\\
C_T&=M_T\sup\|R_TJ_{uu}\|_\infty.
\end{aligned}
\]

The suprema cover the same spatial/parameter/nodal domain as the local range. For the scalar primitive tails \(b_{1,T}(N),b_{2,T}(N)\),

\[
a_T^{\rm tail}\le A_Tb_{1,T}(N),\qquad
b_T^{\rm tail}\le B_Tb_{2,T}(N)+C_Tb_{1,T}(N).
\]

These are force and induced Hessian-row-norm bounds after preconditioning, not scalar-density errors. Cheap conservative prefactors can use \(\|R_T\|_\infty\) times local derivative norms; direct signed polynomial contraction before absolute values may improve them substantially. The surface formula is identical with \(G=J_s^2\) replacing \(J\) and its corresponding derivative tails.

For an old cube of radius \(r\), its exact self-map slack is \(s=r-\eta-Zr\). One concrete policy reserves, for example, one quarter of this **mathematical** slack for the complete truncation contribution:

\[
\sum_T\bigl(a_T^{\rm tail}+r b_T^{\rm tail}\bigr)
+\sum_e\bigl(a_e^{\rm tail}+r b_e^{\rm tail}\bigr)
\le s/4.
\]

This split is a sufficient computational budget, not a measurement tolerance. It does not guarantee that the signed polynomial defects or cache consistency fit in the remaining slack. The final proof must still check the actual complete inequalities \(Z+b<1\) and \(\eta+a_p+(Z+b)r<r\).

An executable bounded selection policy is: freeze a finite degree list and an element/edge budget allocation; choose the first degree whose exact operator-tail bound fits its allocation; retain every attempted degree/bound. Alternatively, raise the degree of the largest current tail contributor until the global budget fits, with a declared maximum degree/work budget. The complete final sum is the acceptance test; no unproved monotonic behavior of an intermediate floating estimate is needed. If the finite list fails, report an unresolved truncation budget. For fixed \(\rho<1\) and finite prefactors the analytic tails tend to zero, but arithmetic/interval wrapping and the actual target difference may prevent a practical persistence proof.

## 4. Uniform finite-cache force correction needs compensated differences

At an auxiliary node, form exact coefficient differences before multiplying the common nodal coordinates. Point subtraction at sufficiently high Arb precision is valid for the existing point pilots. For a full parameter/tube box, a useful signed force correction is

\[
g(q_c+\delta q)-g(q_c)=H(q_c)\delta q+r_g,
\qquad
|(r_g)_a|\le\frac12\sum_{b,c}
\sup|W_{abc}|\,|\delta q_b\delta q_c|.
\]

Keep the signed linear term through global assembly/preconditioning; take absolute values only on the controlled remainder and final interval. A first-order absolute third-derivative bound often suffices for the Hessian difference; a compensated Hessian difference would additionally require validated fourth derivatives. Directly subtracting two broad interval gradients would discard the small coefficient difference.

For `exact_fe_target.surface`, the canonical strain vector is already normalized: \(q_I=(R/r_e,R'/L_I,Z'/L_I)\). Normalize the old raw cache in the same coordinates using

\[
D_c=\operatorname{diag}(\alpha_cL_I,1/L_I,1/L_I),
\quad \alpha_c=(r_{0,c}L_{0,c})^{-1},
\quad q_{c,n}=D_cq_{c,\rm raw},\quad B_{c,n}=D_cB_{c,\rm raw}.
\]

Then the scale-one surface primitive has **exactly** the old scalar \(J_s\), free force and Hessian after the chain rule, while \(q_I-q_{c,n}\) is a genuinely small map-consistency difference. Rescaling only one factor while leaving the tangent normalization inconsistent can create large artificial \(\delta q\) even when \(J_s\) is identical. This normalization identity changes no old energy or physical parameter.

The point-pilot and independent surface-replay results remain point statements. No order adaptation, compensated uniform correction, enlarged tube, new exact-integral root, or complete seventeen-record bridge is executed by this addendum.
