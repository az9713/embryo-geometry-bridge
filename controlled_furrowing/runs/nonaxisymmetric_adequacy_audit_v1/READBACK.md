# Transverse stability of the flat finite-height model

**Finding.** The certified minima control a finite axisymmetric displacement space with zero swirl. They supply no positivity bound for any angular mode \(m\ge1\), the axisymmetric swirl sector, or the meridional continuum outside the P2 space. There is a direct falsification test of the three-dimensional energy restricted to one admissible transverse direction. Strict total compression on a resolved open surface subedge is enough to construct a negative high-\(m\) direction for the declared surface law with area elasticity and tension. No such compression or negative direction at an archived target has been evaluated in this audit. [Scope of the cached minima](../../docs/COHORT-ENERGY-STABILITY.md#what-changes-scientificallyand-what-remains-open) · [Scope of the exact integral minimum](../../docs/EXACT-INTEGRAL-ROOT.md#three-distinct-targets)

## Model and admissible directions

Use reference material coordinates \((R,\theta,Z)\) on the flat cylinder, and write the certified field as \(y_0=r(R,Z)e_r(\theta)+z(R,Z)e_z\). All fields here use the normalized P2 model with reference radius \(R\), height \(H\), frozen material mask \(\chi(R,Z)\), and flat prescribed \(n_0=e_Z\). The exact 3D energy uses \(dV_0=R\,dR\,d\theta\,dZ\), \(dA_0=R\,dR\,d\theta\) on the disks, and the common positive divisor \(S_{\rm binary}\). The \(2\pi\) factor in the meridional code is the already-integrated angular measure. The material cofactor and trace formulas below refer to the exact integral, as opposed to quadrature with rounded sites and weights. The full cylinder wall is clamped, including its circumferential component when the displacement space is enlarged. The two disks are free to move and carry surface energies. There is no imposed axial clamp on the basal disk. [Exact P2 target and normalization](../../scripts/exact_fe_target_safe.py#L397) · [Current finite-height energy with surface tension](../../scripts/quadratic_quadrature_study.py#L62) · [Reference mask](../../scripts/axisymmetric_mechanics.py#L93) · [P2 constraints](../../scripts/axisymmetric_mixed.py#L66) · [Physical energy and units](../../docs/PHYSICAL-PARAMETER-DOMAIN.md#physical-units-before-fitting)

For an angular harmonic \(e^{im\theta}\) with cylindrical amplitudes \((u,v,w)\), the gradient amplitude is
\[
 G_m =
 \begin{pmatrix}
 u_R &(im u-v)/R &u_Z\\
 v_R &(im v+u)/R &v_Z\\
 w_R &imw/R     &w_Z
 \end{pmatrix}.
\]
Take real fields in the variations below. Equivalently, use the paired cosine and sine phases or the Hermitian form of the real Hessian. The \(1/R\) terms include derivatives of the cylindrical basis. The 3D base coefficients in this basis are independent of \(\theta\), so the integrated real bilinear form is orthogonal across different Fourier modes. For any \(m\ge1\) the first variation integrates to zero on the full circle, even for a revolved P2 state that is only stationary in its meridional finite-dimensional space. The constrained P2 Hessian covers the \(m=0,\ v=0\) block. The \(m=0\) swirl block is absent from that space, and reflection \(\theta\mapsto-\theta\) makes its cross term with the zero-swirl \(m=0\) block vanish.

At the axis, impose regularity on the Cartesian vector. Its transverse combinations satisfy \(v_x\pm i v_y=e^{\pm i\theta}(v_r\pm i v_\theta)\), so their angular orders are \(m\pm1\) and the axial order is \(m\). Smooth amplitudes vanish to the corresponding orders for each nonzero Cartesian harmonic. In particular imposing the old axisymmetric radial rule on every amplitude would suppress allowed \(m=1\) displacements. The annular test below vanishes in a neighborhood of the axis, so it sidesteps singular weights and axis matching. The full wall clamp removes global displacement translations and rotations. The common positive force-scale symmetry is in material-parameter space, and creates no null direction in the clamped displacement Hessian. Angular phase equivalence of a nonzero mode comes from rotational symmetry, and also gives no zero-energy direction by itself. [Axis-aware base geometry and embedding](../../docs/ROOT-BODY-GEOMETRY.md) · [Material force-scale symmetry](../../docs/PHYSICAL-PARAMETER-DOMAIN.md#exact-nonlinear-force-scale-symmetry)

The explicit 3D constitutive density in the physical-domain document supplies the natural extension tested here. It uses the same Frobenius norm of \(F\), determinant \(J\), axial stretch norm \(\lVert Fe_Z\rVert\), physical disk area ratio \(J_s\), fixed material support and full wall clamp in three dimensions. The axisymmetric solver alone does not uniquely specify that extension. For example, an objective term proportional to \((F e_\theta^{\rm ref}\cdot F e_Z)^2\), with its material coefficient supported strictly off-axis, vanishes in the no-swirl axisymmetric space and can change transverse curvature. The test concerns the stated extension with no such additional term. Its volumetric term is the **local** squared log penalty. The production Q8/Q12 override only computes projected pressure as a diagnostic, so the older mixed parent is not the volumetric authority. [Three-dimensional constitutive density](../../docs/PHYSICAL-PARAMETER-DOMAIN.md#physical-units-before-fitting) · [Local volume term and pressure diagnostic](../../scripts/quadratic_quadrature_study.py#L73)

## Second variation at the same field

On the physical domain \(J>0\), write \(F=\nabla y_0\), \(G=\nabla v\), \(M=F^{-1}G\), and \(a=\chi gI\ge0\). Since the first and second derivatives of \(\log\det(F+tG)\) at zero are \(\operatorname{tr}M\) and \(-\operatorname{tr}(M^2)\), the bulk integrand of the second variation is
\[
 {\cal Q}_{\rm bulk}
 =\mu\lVert G\rVert_F^2+a\lVert Ge_Z\rVert^2
   +\kappa(\operatorname{tr}M)^2-(\kappa\log J-\mu)\operatorname{tr}(M^2).
\]
In particular \(\operatorname{tr}(M^2)\) is a signed contraction, distinct from the Frobenius norm squared. The active term is positive semidefinite and independent of the variation of \(J\). The log clamp and inversion penalty in the trial solver are inactive on the certified analytic domain. [No-penalty domain and diagonal active Hessian](../../scripts/algebraic_root_feasibility.py#L63)

For \(G=p\otimes q\), \(\operatorname{tr}(M^2)=(q\cdot F^{-1}p)^2\). The coefficient of that square in \({\cal Q}_{\rm bulk}\) is \(\kappa(1-\log J)+\mu\). Thus uniform \(J\le2\) is a simple sufficient condition for strong rank-one ellipticity with lower bound \(\mu |p|^2|q|^2\), as \(\log 2<1\) and \(\kappa>0\). This can be checked from the same-target whole-body \(J\) upper bound. The constitutive law does allow loss of rank-one ellipticity on wider configurations. With \(F=2\,{\rm Id}\), \(\mu=1\), \(\kappa=5000\), and \(p=q=e_R\), the active term vanishes and the rank-one value is \(1+[1+5000(1-3\log2)]/4<0\) because \(\log2>2/3\). This example describes the constitutive density only. Rank-one ellipticity on the actual field still gives no sign theorem for the whole constrained second variation, which includes signed bulk stress and free surface terms.

For one disk let \(a_s=\partial_R y_0,\ b_s=R^{-1}\partial_\theta y_0,\ p_s=\partial_R v,\ q_s=R^{-1}\partial_\theta v\), \(c_s=a_s\times b_s,\ J_s=|c_s|\), and \(D_s=p_s\times b_s+a_s\times q_s\). In the approved nondegenerate surface domain,
\[
 J_s'=\widehat c_s\cdot D_s,\qquad
 J_s''=\frac{|D_s|^2-(\widehat c_s\cdot D_s)^2}{J_s}
             +2\widehat c_s\cdot(p_s\times q_s).
\]
The surface contribution is \(K_s(J_s')^2+\tau_s J_s''\), with total surface stress coefficient \(\tau_s=K_s(J_s-1)+T_s\). A nonnegative baseline \(T_s\) does not ensure \(\tau_s\ge0\) under surface compression, since \(J_s\) measures current area relative to reference area. A full transverse analysis must retain both geometric terms in \(J_s''\), including the signed cross product, and the bulk signed stress term. [Exact edge trace law](../../scripts/exact_fe_target_safe.py#L148) · [Surface scalar density and derivative implementation](../../scripts/algebraic_root_feasibility.py#L88)

## Smallest falsification test

Take one real axial mode
\[
 v=\psi(R,Z)\cos(m\theta)e_z,\quad m\in\mathbb N,\ m\ge1,
\]
with material support in an annular subregion disjoint from the axis and wall. It touches one disk and vanishes near the opposite disk. The angular cross terms integrate to zero and both squared phases integrate to \(\pi\). Let \(h=r/R>0,\ \ell=\sqrt{r_R^2+z_R^2}\), so \(J_s=h\ell\) on each disk. For \(\gamma=\mu+\kappa(1-\log J)\) and \(s=F^{-1}e_z\), the form becomes
\[
 D^2E[y_0](v,v)=\frac{\pi}{S_{\rm binary}}(C_0+m^2 C_2),
\]
\[
\begin{aligned}
 C_0={}&\int_{\Omega_2}
  \left[
   \mu(\psi_R^2+\psi_Z^2)+a\psi_Z^2
   +\gamma(s_R\psi_R+s_Z\psi_Z)^2
  \right]R\,dR\,dZ\\
 &+\sum_s\int\left[K_s(hz_R/\ell)^2+\tau_s h r_R^2/\ell^3\right]\psi_R^2 R\,dR,\\
 C_2={}&\int_{\Omega_2}\mu\frac{\psi^2}{R}\,dR\,dZ
 +\sum_s\int \tau_s \frac{r_R^2}{J_s}\frac{\psi^2}{R}\,dR .
\end{aligned}
\]
The determinant of the axial-perturbed gradient is independent of its azimuthal slot because the unperturbed circumferential row contains only the hoop stretch. The active stretch also has no azimuthal term. The bulk rank-one determinant term contributes only to \(C_0\) because \((F^{-1}e_z)_\theta=0\). For the disk formula \(p_s\) and \(q_s\) are parallel, so \(p_s\times q_s=0\), and the circumferential perturbation changes the norm of the area vector at second order. With \(z_R=0\), the surface formula reduces to \(J_s''=(h/r_R)(\psi_R\cos m\theta)^2+(r_R/h)(m\psi\sin m\theta/R)^2\) for \(r_R>0\).

**One scalar sufficient screen.** On a predeclared strictly off-axis, interior apical subedge, use the same certified posterior nodal intervals and exact P2 tangents to bound
\[
 \tau_a=K_a\left[\frac rR\sqrt{r_R^2+z_R^2}-1\right]+T_a .
\]
Require a strict negative upper bound on a closed nonzero-width interval within one edge, plus positive \(r_R\) and nondegenerate \(J_s\) bounds there. A valid enclosure gives a uniform \(-\tau_a r_R^2/J_s\ge c>0\) on a smaller interior radial interval. The existing exact whole-box trace or a local exact Bernstein calculation supplies the geometric quantities. For \(K_a>0\) the inequality can also be tested directly by bounding \(J_s^2<(1-T_a/K_a)^2\) with positive \(1-T_a/K_a\). This avoids taking a square root for the sign decision.

Choose a nonzero smooth radial cutoff \(f(R)\) supported in that interval and write \(\psi=f(R)w((H-Z)/d)\), with \(w(0)=1\), \(|w|\le1\), and smooth cutoff to zero at depth \(d\). The active set and geometry remain frozen. Then
\[
 C_2\le (\mu d-c)\int f(R)^2/R\,dR<0
\]
by choosing \(d<c/\mu\) and \(d<H\). Place the compact radial support strictly inside one surface element and further decrease the depth so this support lies inside its adjacent material triangle. Hold the resulting smooth function fixed as the angular frequency increases. This thin radial–axial trial belongs to the enlarged continuum test space, which is richer than the fixed P2 meridian with added Fourier factors. For each fixed such \(d\), \(C_0\) is finite on the certified physical box. Large enough \(m\) then gives a strictly negative second variation even if its bulk contribution stays rank-one elliptic. This surface localization argument uses the absence of bending or other higher-gradient terms in the actual surface law. Under a uniform lower bound on the symmetric meridional gradient and hoop stretch, the 3D symmetric gradient of the base is positive definite almost everywhere in the cylinder. The chosen annular mode has bounded gradient at every fixed \(m\). Sufficiently small amplitude preserves that strict bound, hence injectivity by strong monotonicity in the convex reference cylinder. It preserves the wall clamp and the local determinant and tangent guards as well. The zero first variation on that mode and negative second variation then give admissible states of lower energy. A sign decision at one uncertified sampled root or isolated uncertain surface point has no such force. [Whole-body strong monotonicity](../../docs/ROOT-BODY-GEOMETRY.md#3-prove-global-injectivity-for-every-field-in-the-box)

The next bounded computation would check this scalar surface sign before assembling any transverse matrix. If the sufficient screen is unresolved, an inexpensive finite-mode alternative is one annular axial P2 trial with \(m=1\), and a rigorous signed upper bound on \(C_0+C_2\) over the same root enclosure. A strictly negative upper bound establishes a lower-energy perturbation in that direction. A nonnegative value in a single direction, or nonnegative surface stress on all tested edges, leaves the general question open. A complete finite-mode test requires the radial, swirl and axial components together, with exact boundary constraints, axis regularity, certified root-dependent bulk and surface coefficients, and a treatment of higher modes. A continuum stability claim would further require coercivity in the enriched function space and control of discretization and model validity.

The proposed high-mode implication concerns the declared finite-height energy evaluated at a certified revolved field and enlarged admissible directions. It makes no claim about a continuum equilibrium, a physically resolved high-frequency tissue instability, or which model additional bending or contact would support. Neither screen has been run on a scientific target in this audit.



## Relationship to prior work

Instability from negative surface stress and the role of bending in regularizing short wavelengths are established phenomena. Binysh, Wilks and Souslov study a bulk elastic solid with imposed negative surface tension and examine bending and bulk dispersion as large-wavenumber regularizations. The scalar in this audit is the derivative of the declared area-elastic energy at a specific nonlinear finite-height field. Andrenšek, Ziherl and Krajnc study a different epithelial wrinkling mechanism in an unsupported cellular sheet with differential apico-basal tension. These papers give context for the mechanism, while the audited root enclosures and conditional test concern this specific model. The analysis makes no general novelty claim. [Binysh, Wilks and Souslov, 2022](https://arxiv.org/pdf/2101.04006) · [Andrenšek, Ziherl and Krajnc, 2023](https://arxiv.org/abs/2210.11883)
