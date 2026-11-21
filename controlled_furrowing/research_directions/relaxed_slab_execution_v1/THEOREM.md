# A constructive finite-height slab recovery theorem

**Shape convention:** use the [canonical deformation-map and surface guide](../../../docs/SHAPE-DESCRIPTION.html). This theorem concerns an explicit three-dimensional material embedding, not a force-balanced embryo or a certified minimizer. All lengths and energies below are nondimensionalized by one fixed length and bulk shear scale. Physical surface bending coefficients are converted by \(b_s=B_s/(\mu L_0^3)\).

## Statement and scope

Let \(\Omega=(0,L_1)\times(0,L_2)\times(-H,0)\), with affine-periodic lateral boundary conditions, and let

\[
y_0(X)=F_0X,\quad F_0=\operatorname{diag}(\lambda_1,\lambda_2,\lambda_3),\quad \lambda_i>0.
\]

On each free horizontal face take \(\psi_s(J)=K_s(J-1)^2/2+T_sJ\), where \(K_s>0\), \(T_s\ge0\). Set \(j_0=\lambda_1\lambda_2\) and \(j_s^*=\max(0,1-T_s/K_s)\). Define the candidate relaxed surface density by \(\psi_{s,\mathrm{rel}}(J)=\psi_s(\max(J,j_s^*))\). Only free faces with \(j_0<j_s^*\) are corrugated. A prescribed/clamped horizontal face is unchanged and its energy is not relaxed in the target.

Assume the bulk density \(W(X,F)\) is a Caratheodory integrand: measurable in \(X\) for every fixed \(F\), and continuous in \(F\) for almost every \(X\). Assume it is uniformly bounded on the compact positive-determinant gradient set constructed below. These assumptions make the composed bulk density measurable and integrable. They hold for the project's log-volume, quadratic shear and prescribed active-fiber density with bounded coefficients and bounded activation. No convexity of the bulk is needed for this upper-bound construction.

Let

\[
E_{\mathbf b}[y]=\int_\Omega W(X,\nabla y)\,dX+
 \sum_s\int_{\Gamma_{s,0}}\psi_s(J_s[y])\,dA_0+
 \sum_s\frac{b_s}{2}\int_{\Gamma_s[y]}(k_1^2+k_2^2)\,da.
\]

The bending norm is the invariant curvature norm, computed with the current surface metric. Let \(E_{\mathrm{rel}}[y_0]\) replace the free-face densities by \(\psi_{s,\mathrm{rel}}\) at the affine field, leaving any prescribed face unchanged.

**Theorem.** For fixed admissible material and geometry parameters, there are explicit constants and, for all sufficiently small \(b=\max_s b_s>0\), globally injective \(C^{1,1}\) embeddings \(y_{\mathbf b}\) satisfying the affine periods such that

\[
\det\nabla y_{\mathbf b}\ge\frac14\det F_0>0,
\qquad
\|y_{\mathbf b}-y_0\|_\infty\le C_y b^{1/3},
\]

\[
|E_{\mathbf b}[y_{\mathbf b}]-E_{\mathrm{rel}}[y_0]|
 \le C_E b^{1/3}.
\]

For every finite \(p\ge1\), the bulk gradients also converge strongly in \(L^p\), at the explicit construction rate \(O(b^{1/(3p)})\). Surface gradients need not converge strongly. If there are no compressed free faces, simply take \(y_{\mathbf b}=y_0\).

This is a constructive recovery bound at a diagonal affine field in a periodic slab of **fixed positive height**. It is not a sharp-rate theorem, a full Gamma-convergence theorem, a minimizer-convergence theorem, a nonlinear clamped-cylinder theorem, or a biological prediction.

## 1. A periodic angle with exactly the right mean

For each compressed free face write

\[
v_s=j_s^*/\lambda_2,\quad q_s=\lambda_1/v_s\in(0,1),\quad \delta_s=q_s/4.
\]

Temporarily suppress the face index. On one period, define a continuous piecewise-linear angle \(\theta\): it dwells at \(+\alpha\) for length \((1-\delta)/2\), falls linearly to \(-\alpha\) over length \(\delta/2\), dwells at \(-\alpha\) for length \((1-\delta)/2\), then rises linearly to \(+\alpha\) over length \(\delta/2\). Extend periodically. It is Lipschitz and satisfies \(\theta(t+1/2)=-\theta(t)\).

Its means and squared derivative integral are exactly

\[
\int_0^1\sin\theta\,dt=0,\quad
M(\alpha):=\int_0^1\cos\theta\,dt
 =(1-\delta)\cos\alpha+\delta\frac{\sin\alpha}{\alpha},
\]

\[
\int_0^1|\theta'|^2dt=\frac{16\alpha^2}{\delta}.
\]

Interpret the quotient by continuity at zero. There is a unique

\[
0<\alpha<\alpha_{\max}=\arccos(q/2)
\quad\text{such that}\quad M(\alpha)=q.
\]

Indeed, \(M(0)=1>q\), while

\[
M(\alpha_{\max})\le q/2+\delta=3q/4<q.
\]

On \((0,\pi/2)\), both \(\cos\alpha\) and \(\sin\alpha/\alpha\) decrease strictly. For the latter, \(\sin\alpha-\alpha\cos\alpha>0\), since its derivative is \(\alpha\sin\alpha>0\) and its value at zero is zero. Continuity and strict monotonicity prove existence and uniqueness. In particular

\[
v\cos\theta(t)>v q/2=\lambda_1/2.
\]

## 2. A constant-area corrugated face

Define periodic primitives by

\[
p'(t)=v\cos\theta(t)-\lambda_1,\qquad
h'(t)=v\sin\theta(t),\qquad p(0)=h(0)=0.
\]

Both derivatives have zero mean, so \(p,h\) are periodic. Crude but explicit bounds are

\[
|p|\le2v,\quad |h|\le v,\quad |p'|\le2v,\quad |h'|\le v.
\]

For a period \(\varepsilon=L_1/N\), place the face at

\[
x(X_1,X_2)=(\lambda_1X_1+\varepsilon p(X_1/\varepsilon),\
 \lambda_2X_2,\ z_{\mathrm{face}}+\varepsilon h(X_1/\varepsilon)).
\]

Its two tangent vectors are \((v\cos\theta,0,v\sin\theta)\) and \((0,\lambda_2,0)\). Hence

\[
J_s=v\lambda_2=j_s^*
\]

pointwise, including the transition regions. There is no transition-area penalty. Its intrinsic stretch does not oscillate, although its tangent direction does. Its first tangent component is strictly positive, so the face is a graph in the physical horizontal coordinate and has no self-overlap.

This improves on an abrupt two-slope laminate: smoothing a jump in its tangent would add an area-energy error. Here the tangent turns along a constant-speed arc. Lipschitz \(\theta\) makes the face \(C^{1,1}\), with finite squared curvature even though \(\theta'\) has jumps.

## 3. Extend into the finite-height body

Use a collar cutoff \(\chi\) with values in \([0,1]\), equal to zero below the collar and one at the face. The cubic smoothstep gives

\[
\max|\chi'|=M_\chi=3/2,
\]

where the derivative is with respect to the dimensionless collar coordinate. Put \(d=C\varepsilon\), with

\[
C\ge\max\left(1,\frac{12M_\chi v_{\max}^2}{\lambda_1\lambda_3}\right),
\qquad v_{\max}=\max_s v_s.
\]

In a top-face collar set

\[
y_1=\lambda_1X_1+\varepsilon p(X_1/\varepsilon)\chi(Z/d),\quad
y_2=\lambda_2X_2,\quad
y_3=\lambda_3Z+\varepsilon h(X_1/\varepsilon)\chi(Z/d).
\]

Here \(\chi\) runs from zero at \(Z=-d\) to one at \(Z=0\). Use the reflected cutoff for a bottom-face collar. Choose \(2d<H\) if both faces are corrugated, so the collars are disjoint. Away from them use \(y_0\). The cutoff has zero derivative at its joining endpoints, so the complete embedding is \(C^{1,1}\).

## 4. Positive determinant and global injectivity

Within either collar define \(A=\partial_1y_1\). Since it is a convex combination of \(\lambda_1\) and \(v\cos\theta\),

\[
\lambda_1/2\le A\le v.
\]

The determinant of the \((X_1,Z)\)-submap is

\[
D=\lambda_3 A+\frac\varepsilon d\chi'
 [A h-\chi h'p].
\]

The bracket has absolute value at most \(3v^2\). The chosen collar depth therefore gives

\[
D\ge\lambda_1\lambda_3/2-
 \frac{3M_\chi v_{\max}^2}{C}
 \ge\lambda_1\lambda_3/4.
\]

Multiplying by \(\lambda_2\) proves the determinant bound. It also gives the useful upper bound

\[
\det\nabla y\le\lambda_2\lambda_3(v_{\max}+\lambda_1/4).
\]

Local positivity is not used as a substitute for global injectivity. Extend the periodic formula to all \(X_1\in\mathbb R\). At every fixed \(Z\), the map \(X_1\mapsto y_1\) is strictly increasing and onto, because its derivative is at least \(\lambda_1/2\) and its correction to an affine map is bounded. For each physical horizontal coordinate \(u\), solve \(y_1(X_1,Z)=u\) uniquely. Along that curve,

\[
\frac{d y_3}{dZ}\Big|_{y_1=u}=\frac{D}{A}>0.
\]

Thus equal \(y_1,y_3\) force equal \(Z,X_1\). Equal \(y_2\) forces equal \(X_2\). The map on the infinite periodic slab is injective, and so is its restriction to a period cell. This argument also excludes contact between the two free faces. No new lumen or topology change is introduced.

## 5. Explicit bulk and bending costs

Let \(n\in\{1,2\}\) be the number of corrugated faces and \(A_0=L_1L_2\). The gradient differs from \(F_0\) only in volume \(nA_0d\). One uniform squared-gradient bound is

\[
F_{\max}^2=2v_{\max}^2+\lambda_2^2+
 (2M_\chi v_{\max}/C)^2+
 (\lambda_3+M_\chi v_{\max}/C)^2.
\]

Together with the determinant bounds this defines a compact subset of positive-J gradients. If \(|W|\le M_W\) on that set, including \(F_0\), then

\[
\left|\int_\Omega[W(X,\nabla y)-W(X,F_0)]dX\right|
 \le 2nA_0 C M_W\varepsilon.
\]

For the project's bulk law, an explicit possible bound is

\[
M_W=\frac\mu2(F_{\max}^2+3+2L_J)
 +\frac\kappa2L_J^2+\frac{a_{\max}}2F_{\max}^2,
\]

where \(L_J\) bounds \(|\log J|\) on the determinant interval and \(0\le a(X)\le a_{\max}\) bounds the active coefficient for a prescribed unit fiber \(n_0\). This bound is deliberately conservative and is not an experimentally inferred constant.

The corrugated face is an extruded constant-speed planar curve. Its principal curvatures are \(\theta'/(\varepsilon v)\) and zero, while its area element is \(v\lambda_2dA_0\). Thus its bending cost is exactly

\[
\frac{b_s}{2}\int(k_1^2+k_2^2)da
 =\frac{8A_0b_s\lambda_2\alpha_s^2}{v_s\delta_s\varepsilon^2}.
\]

Its area energy equals the relaxed target exactly. Unchanged noncompressed faces remain flat and already have their relaxed density. Put

\[
C_{\rm bulk}=2n C M_W,\qquad
C_{\rm bend}=\sum_s\frac{8\lambda_2\alpha_s^2}{v_s\delta_s}.
\]

The triangle inequality gives the complete construction estimate

\[
|E_{\mathbf b}[y_\varepsilon]-E_{\mathrm{rel}}[y_0]|
 \le A_0\left(C_{\rm bulk}\varepsilon+
            C_{\rm bend}\frac b{\varepsilon^2}\right).
\]

This controls the absolute difference for the constructed fields. It is not a lower bound over all competing deformations.

## 6. Select integer periods and obtain the rate

Let \(\ell=b^{1/3}\), \(N=\lceil L_1/\ell\rceil\) and \(\varepsilon=L_1/N\). If \(\ell\le L_1\), then \(\ell/2\le\varepsilon\le\ell\). If also \(\ell<H/(2C)\), the collars fit without overlap. Consequently

\[
C_E=A_0(C_{\rm bulk}+4C_{\rm bend})
\]

is a valid energy constant. The primitive bounds give

\[
C_y=\sqrt5\,v_{\max},\qquad
\|y_\varepsilon-y_0\|_\infty\le C_y\varepsilon\le C_y b^{1/3}.
\]

Finally, the bounded gradient correction is supported in volume \(nA_0C\varepsilon\). This proves the stated finite-p bulk-gradient convergence. It does not imply strong convergence of the surface tangents, which keep rotating through angles of nonvanishing size.

For a fixed material-ROI mean-height contrast, the direct construction also gives \(|Q(y_\varepsilon)-Q(y_0)|\le2v_{\max}\varepsilon\). This only concerns the constructed sequence. It does not establish that an energy minimizer or a loading trajectory retains that coarse observable.

## 7. What has and has not been established

The construction answers the periodic finite-slab question affirmatively: a preferred-area free surface can be realized through fine folds while preserving a quantitative positive determinant, global injectivity, finite curvature energy and fixed positive body height. Both free faces can be treated with disjoint collars. A bottom prescribed to remain affine can instead be left untouched.

There is no matching lower bound proving the exponent optimal. No full liminf theorem, compactness theorem for arbitrary minimizing sequences, variable-base injective construction, wall-clamped cylinder recovery, observer-growth constant, global minimizer or developmental-path selection has been proved here. Existing curved finite-element fields are not made curvature-admissible by this theorem. A biological application still requires a valid cell/continuum scale window and justified physical regularization.

The mathematical ingredients have substantial prior art. Smooth periodic folds and boundary matching appear in [Conti, DeSimone and Müller (2005), Section 3.1](https://people.sissa.it/~desimone/Publications/CDM2005.pdf). Negative surface stress and bending/bulk regularization are treated in [Binysh, Wilks and Souslov (2022)](https://arxiv.org/html/2101.04006v2). This is an explicitly scoped constructive adaptation with quantitative guards; no claim of a new general concept or publication priority is made.
