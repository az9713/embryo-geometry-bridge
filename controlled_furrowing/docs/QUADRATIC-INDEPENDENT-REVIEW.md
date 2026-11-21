# Independent review of the forced quadratic mechanical benchmark

**Result:** the frozen P2 pointwise-penalty implementation and all nine saved **forced solutions** pass 94 independent checks. This is stronger evidence than checking interpolation of a known map. It remains numerical verification on a synthetic problem, not a certified continuum-error estimate or biological validation. The newer exact correlated Bernstein check additionally proves that all six saved solutions on meshes 8×4 and 12×6 have less than 1% local volume distortion everywhere in their discrete bodies.

Primary receipt: `evidence/quadratic_independent_review.json`. It contains source hashes, every state hash, the replay program, independent norm evaluations, exact coefficient certificates and strong-monotonicity certificates. Solver source SHA256:

```text
5bb64d6b441f6db8ada9da0f1bcb44716c86ed5d816e44f75bd4e9e021fdf410
```

No biological response files, labels or fits were read for this review.

## 1. Independently derived exact nonlinear solution

Let material cylindrical coordinates be \((r,\theta,z)\), and choose \(\epsilon=0.2\). Define

\[
 A=e^{\epsilon z},\qquad B=e^{-2\epsilon z},\qquad
 R=rA,\quad\Theta=\theta,\quad Z=\frac{1-B}{2\epsilon}.
\]

In orthonormal cylindrical bases ordered \((r,\theta,z)\),

\[
 F=\begin{pmatrix}A&0&\epsilon rA\\0&A&0\\0&0&B\end{pmatrix},
 \qquad J=A^2B=1.
\]

The original pointwise material law is

\[
 W=\frac\mu2(I_1-3-2\log J)+\frac\lambda2(\log J)^2,
 \qquad P=\mu(F-F^{-T})+\lambda\log J\,F^{-T}.
\]

Since this map has \(\log J=0\), its stress and required body force are independent of \(\lambda\). Its nonzero first Piola components are

\[
 P_{rr}=P_{\theta\theta}=\mu(A-A^{-1}),\quad
 P_{rz}=\mu\epsilon rA,\quad
 P_{zr}=\mu\epsilon r/B,\quad
 P_{zz}=\mu(B-B^{-1}).
\]

The reference cylindrical divergence includes geometric terms:

\[
 (\operatorname{Div}P)_r=\partial_rP_{rr}+\partial_zP_{rz}
                         +(P_{rr}-P_{\theta\theta})/r,
 \quad
 (\operatorname{Div}P)_z=\partial_rP_{zr}+\partial_zP_{zz}+P_{zr}/r.
\]

Substitution gives the independent force

\[
 f=-\operatorname{Div}P
 =\left(-\mu\epsilon^2r e^{\epsilon z},\quad
          2\mu\epsilon e^{-2\epsilon z}\right),
\]

exactly matching `ManufacturedModel`. Symbolic computation verifies \(J=1\) and both force components. This is an exact stationary solution of the continuum bulk-energy-plus-body-force problem under its own Dirichlet trace, for every finite positive \(\lambda\). No claim of unique global nonlinear energy minimization is needed or established.

## 2. What is actually solved

`scripts/axisymmetric_quadratic.py` contains a real forced finite-element problem. Both components of every boundary displacement node—including the axis—are set to the analytic map. Interior displacement degrees of freedom are optimized against the bulk material energy and body-force potential, with no activation or surface energy. This stronger-than-biological boundary constraint is appropriate to this manufactured benchmark.

The potential is \(-\int f\cdot y\,dV_0\). Replacing \(y\) by displacement \(u=y-x\) changes it by the reference constant \(-\int f\cdot x\,dV_0\), so the implemented force sign is correct. Each saved field was reconstructed independently from its free vector and prescribed boundary and compared to saved current nodes. All boundary nodes match the prescribed map; all interior solutions differ from their initial analytic interpolants. Replay of the force residual confirms actual forced stationarity for all nine states.

The exact Dirichlet trace is nonpolynomial. P2 boundary edges interpolate its nodal values and therefore approximate its continuous trace. The numerical PDE has this boundary interpolation and degree-four volume quadrature. These approximation errors belong in the numerical accuracy discussion; exact nodal boundary agreement does not remove them.

## 3. Why restoring the pointwise penalty matters

Let \(g=\log J\) evaluated at positive-weight quadrature points, let \(W_q\) be the quadrature weight matrix, and let \(N\) interpolate the continuous P1 pressure basis. Set

\[
 M=N^TW_qN,\quad b=N^TW_qg,\quad
 \Pi g=NM^{-1}b.
\]

The earlier condensed mixed volumetric term was \(\lambda b^TM^{-1}b/2\). Orthogonality in the **same quadrature inner product** gives

\[
 \frac\lambda2g^TW_qg
 =\frac\lambda2b^TM^{-1}b
 +\frac\lambda2(g-\Pi g)^TW_q(g-\Pi g).
\]

The final nonnegative term is precisely what a pressure projection does not penalize. A very small projected pressure residual does not bound local \(J\). The new subclass adds this unresolved term back, restoring the original pointwise logarithmic-volume penalty at quadrature points. Our independent assembly of \(M,b\), full energy, projected energy and remainder verifies this identity on all nine saved fields. This algebra is a finite-quadrature identity; it is not a proof of exact incompressibility.

## 4. Independent reconstruction of three-dimensional errors

The reviewer reconstructs the six-node P2 field from its polynomial coefficients, rather than calling the implementation's error routine. With \(dV_0=2\pi r\,dr\,dz\), its normalized displacement and gradient errors are

\[
 E_0^2=\frac1{V_0}\int|y_h-y_*|^2dV_0,
 \qquad
 E_1^2=\frac1{V_0}\int
 \left(|F_{m,h}-F_{m,*}|_F^2+
        |R_h/r-A|^2\right)dV_0.
\]

The hoop term is required for the full three-dimensional gradient norm. The label H1 RMS in the receipts denotes this gradient **seminorm**, without adding the displacement term. We also recompute \(E_J^2=V_0^{-1}\int(J_h-1)^2dV_0\).

Independent degree-four quadrature reproduces the recorded norms to rounding error. Higher 8×8 and 12×12 tensor-Gauss/Duffy quadrature evaluates the same saved fields more densely. It changes H1 RMS very little, but the original degree-four displacement L2 RMS can underestimate its higher-quadrature value by approximately 1–16%. Sampling-based maximum \(|J-1|\) also increases when more points are evaluated, as expected; it was never an upper bound.

Using the independently computed 12×12 quadrature:

| Mesh | \(\lambda/\mu\) | Displacement RMS | Gradient H1 RMS | \(J-1\) RMS |
|---|---:|---:|---:|---:|
| 6×3 | 5 | 2.01278e−4 | 2.21985e−3 | 1.83490e−3 |
| 8×4 | 5 | 2.81636e−5 | 6.51380e−4 | 4.90860e−4 |
| 12×6 | 5 | 8.10214e−6 | 2.82097e−4 | 2.20581e−4 |
| 6×3 | 50 | 2.34459e−4 | 2.90729e−3 | 1.79511e−3 |
| 8×4 | 50 | 3.19749e−5 | 7.48476e−4 | 4.85378e−4 |
| 12×6 | 50 | 8.96501e−6 | 3.13402e−4 | 2.18925e−4 |
| 6×3 | 500 | 5.91889e−4 | 8.26070e−3 | 1.70937e−3 |
| 8×4 | 500 | 8.58371e−5 | 1.76046e−3 | 4.72260e−4 |
| 12×6 | 500 | 2.51634e−5 | 8.00613e−4 | 2.13342e−4 |

These values improve with refinement for every tested \(\lambda\). Higher \(\lambda\) causes larger displacement and gradient errors even though the exact solution is unchanged. That is evidence of remaining parameter sensitivity, not a locking-free result. Three nonuniform mesh configurations and sampled norms do not establish a formal convergence order, uniform-in-\(\lambda\) stability theorem or a posteriori PDE-error certificate.

## 5. Exact correlated Bernstein bounds for local volume

The original geometry verifier bounded the product \(J=\det F_m\,(R/r)\) by separate factor bounds. This loses correlations between the factors. The new independent helper `scripts/mixed_geometry_certificate_v2.py` leaves that verifier unchanged and combines the numerator first.

Let \(d=\sum_{|\alpha|=2}d_\alpha B_\alpha^2\) be the meridional determinant, and \(R=\sum_{|\beta|=2}a_\beta B_\beta^2\). Their product has degree-four coefficients

\[
 n_\gamma=\sum_{\alpha+\beta=\gamma}
 \frac{\binom2\alpha\binom2\beta}{\binom4\gamma}
 d_\alpha a_\beta.
\]

The affine reference radius elevated to degree four has coefficient
\(c_\gamma=(\gamma_0r_0+\gamma_1r_1+\gamma_2r_2)/4\ge0\).
Require \(n_\gamma=0\) wherever \(c_\gamma=0\). Nonnegative Bernstein weights then prove, for \(r>0\) and its axis limits,

\[
 \min_{c_\gamma>0}\frac{n_\gamma}{c_\gamma}
 \le J=\frac{\sum n_\gamma B_\gamma^4}{\sum c_\gamma B_\gamma^4}
 \le\max_{c_\gamma>0}\frac{n_\gamma}{c_\gamma}.
\]

All arithmetic is exact rational arithmetic on saved binary64 node values. The retained coefficients establish the following whole-body **discrete** intervals; displayed decimals are approximations of exact rational receipt bounds:

| Mesh | \(\lambda/\mu\) | Certified lower \(J\) | Certified upper \(J\) |
|---|---:|---:|---:|
| 6×3 | 5 | 0.97781355 | 1.00781698 |
| 6×3 | 50 | 0.97205401 | 1.00782806 |
| 6×3 | 500 | 0.97099587 | 1.00822198 |
| 8×4 | 5 | 0.99529826 | 1.00204807 |
| 8×4 | 50 | 0.99404499 | 1.00195240 |
| 8×4 | 500 | 0.99311832 | 1.00187583 |
| 12×6 | 5 | 0.99803254 | 1.00088647 |
| 12×6 | 50 | 0.99751008 | 1.00083211 |
| 12×6 | 500 | 0.99693388 | 1.00077150 |

For example, the 12×6, \(\lambda=500\) separate-product interval was approximately [0.95673,1.04311]; the correlated interval is [0.996934,1.000772]. This sharpening is an improvement in certification, not a change in the solved field. The stronger bounds do not extend to the unknown continuum solution or to a different activated biological scenario. Rational Bernstein bounding is an established technique, not a new conceptual discovery; see the source discussion in `MIXED-GEOMETRY-PROOF.md`.

## 6. Exact global injectivity for this manufactured family

The manufactured rim is curved, so the prior vertical-rim Jordan-boundary routine does not apply. Instead, we verify a sufficient strong-monotonicity condition. On each affine reference triangle, \(F_{m,h}\) is affine. Its symmetric part \(S_h\) is therefore an affine combination of its vertex values. At each vertex we verify exactly

\[
 S_{11}>0,\qquad\det S>0.
\]

These are Sylvester's conditions for a positive definite symmetric 2×2 matrix. Moreover, \(\lambda_{\min}(S)\ge\det(S)/\operatorname{tr}(S)>0\). Take the minimum of those exact rational lower bounds over all element vertices. Convexity of the positive-definite cone yields \(S_h\ge mI\) everywhere for one \(m>0\). Shared P2 edge nodes verify continuity of the complete quadratic traces; disk topology and the rectangular reference boundary are checked separately.

The reference meridian is convex. Integrating along the segment from \(x\) to \(x'\), across finitely many element pieces, gives

\[
 (y_h(x')-y_h(x))\cdot(x'-x)\ge m|x'-x|^2.
\]

Thus distinct reference points have distinct images. All nine discrete maps pass this exact sufficient global meridional injectivity check, preserve the axis, and have positive certified axisymmetric \(J\); their revolutions are injective bodies. The condition is sufficient and restrictive. Failure would be `UNCERTIFIED`, not proof of self-intersection.

## 7. What has and has not been established

We have verified correct forcing, genuine discrete forced equilibrium, correct three-dimensional error definitions, empirical refinement improvement, exact whole-element orientation/volume bounds, and global injectivity for nine saved synthetic bodies. The helper's separate exact synthetic suite passes 45 checks, including independent analytic product identities and a constant-J affine case; see `evidence/mixed_geometry_v2_synthetic.json`.

We have not proved that the solver is locking-free, that the active large-deformation biological experiment is converged, or that any fitted geometry is accurate. The manufactured map is smooth, externally forced and boundary prescribed. Returning to active furrowing still requires matched-load/mesh/quadrature checks and local-volume control on those actual activated states. No new biological finding or novel general mechanics theorem follows from this benchmark.
