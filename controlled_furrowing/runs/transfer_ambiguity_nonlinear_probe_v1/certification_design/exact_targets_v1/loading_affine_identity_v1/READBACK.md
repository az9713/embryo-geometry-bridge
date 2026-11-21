# Exact dose affinity and a loading-path certificate

**Scope.** This is a source readback for one fixed flat axisymmetric P2 target, with fixed material, reference mesh, free-coordinate injection, reference activation mask and positive energy normalization. The same formula holds separately for each such target. It applies to the exact FE integral and to the real-valued potential of the frozen binary quadrature, with different coefficients. No loading interval or path is certified here.

The weak-load cached theorem already proves a smooth, stable finite-dimensional branch on one short gain interval. Its predictor residual retains active Hessian variation in the Hessian integral. The independent exact-spatial L proof also uses the pointwise gain derivative of the bulk integral/cache defect. The additional implementation consequence here is that the activation stiffness is an exactly integrable constant Gram matrix even at the axis. It can be separated from all nonlinear spatial variation and retained as one exact quadratic-in-dose predictor term. [Prior gain-cell theorem](../../../../../docs/ENERGY-GAIN-CELL.md) · [Whole-body composition](../../../../../docs/ENERGY-GAIN-GEOMETRY.md) · [L defect proof](../../../../../docs/COHORT-INTEGRATION-L.md)

## Force identity

Let \(I\) be the requested optical dose, \(g\ge0\) the fixed common gain, \(Y=X+Pu\) the full normalized vector of P2 coordinates and \(P\) the injection of the free components. Let \(D_Z\) take the two material axial derivatives of the P2 interpolant on each element. With \(\chi_e\in\{0,1\}\) the fixed reference activation flag, write \(B_e=D_ZP\), \(c_e=D_ZX\), and \(M=2\pi/S_{\rm binary}>0\). The frozen active direction in these flat models is material \(n_0=e_Z\). The active energy per unit drive is
\[
 V(u)=\frac{M}{2}\sum_e\chi_e\int_{T_e}R_{\rm ref}\,\|c_e+B_eu\|_2^2\,dR_{\rm ref}\,dZ.
\]
Here \(F=\nabla_u E\) is the code's free residual sign convention. Differentiating gives
\[
 F(u,I)=F_{\rm pass}(u)+gI(Au+b),\qquad
 A=M\sum_e\chi_e\int_{T_e}R_{\rm ref}B_e^TB_e,\quad
 b=M\sum_e\chi_e\int_{T_e}R_{\rm ref}B_e^Tc_e.
\]
The integrals in the definitions have measure \(dR_{\rm ref}\,dZ\). Symmetry follows termwise, and \(v^TAv=M\sum_e\chi_e\int R_{\rm ref}\|B_ev\|^2\ge0\). The kernel may be nontrivial. The active term has no cross-coupling between radial and axial displacement components in this flat law. Fixed coordinates still contribute to \(c_e\) and therefore to \(b\). The literal binary midside coordinates in \(X\) are authoritative. Replacing them by ideal exact midsides would change the target.

- [Exact target basis and gradient maps](../../../../../scripts/exact_fe_target_safe.py#L93), [full nodal vector and stored mask](../../../../../scripts/exact_fe_target_safe.py#L184), [energy and reference normalization definition](../../../../../scripts/exact_fe_target_safe.py#L397).
- [Production Q8/Q12 volume and surface energy](../../../../../scripts/quadratic_quadrature_study.py#L35), [full pointwise energy and active axial column](../../../../../scripts/quadratic_quadrature_study.py#L62), [normalized AD residual](../../../../../scripts/axisymmetric_mechanics.py#L144), [sparse production Hessian uses that same energy](../../../../../scripts/sparse_quadratic.py#L76).
- [Reference-fitted radial/depth grids and activation mask](../../../../../scripts/axisymmetric_mechanics.py#L68), [P2 axis/rim constraints and lower radial trial bounds](../../../../../scripts/axisymmetric_mixed.py#L53), [frozen exact-dose and rounded solver strength distinction](../../../../../scripts/transfer_ambiguity_exact_targets.py#L232).
- [Exact integrated bulk uses activation only on two axial entries](../../../../../scripts/exact_cohort_integrator.py#L227), [stored-cache local derivatives](../../../../../scripts/algebraic_root_feasibility.py#L63), [exact arithmetic cache weight and scatter rules](../../../../../scripts/algebraic_root_feasibility.py#L144).

The hoop stretch \(r/R_{\rm ref}\) belongs to the passive shear, Jacobian/log-volume and surface terms. The factor in the active integral is the **reference** radius, which stays fixed. The older curved P1 kernel uses a normalized reference orientation field and is still quadratic in the nodal vector because the orientation is prescribed. Its orientation field and coefficients differ from the flat target. A current-shape normal, current-volume weighting, or a dose-/shape-dependent mask would require a new derivation. Changing patch radius or depth changes the mesh/support here, so different targets need separate coefficients. [Older reference orientation and active term](../../../../../scripts/axisymmetric_mechanics.py#L93)

The solver's log clamp and inversion penalty depend on \(u\) and are dose-independent, so they preserve additive dose affinity wherever the relevant derivatives exist. They belong to its trial passive force. The physically certified analytic operator employs the ordinary log at positive hoop and meridional determinant and a strict \(J>10^{-6}\) margin, with the trial penalty zero. Surface graph/tangent and global embedding guards are separate. A full loading proof must stay inside those domains and away from numerical inequality bounds if inferring unconstrained stationarity. At an active bound, constrained equilibrium and its stability require a separate constrained argument. The fixed axis and rim constraints are eliminated by \(P\). [Trial guards in the production law](../../../../../scripts/quadratic_quadrature_study.py#L66) · [Analytic no-penalty guard](../../../../../scripts/algebraic_root_feasibility.py#L63) · [Solver bounds and verification](../../../../../docs/QUADRATIC-SOLVE-POLICY.md)

## Exact integration and cached arithmetic

Each \(B_e\) and \(c_e\) entry is a degree-at-most-one polynomial on a P1 material triangle. Thus every active force coefficient involves only degree-at-most-three polynomial moments after multiplying by \(R_{\rm ref}\). This contribution has no inverse radius, log, square-root or axis singularity. It is valid to integrate it on the unsplit material triangles, even for axis-touching elements. Expand in barycentric coordinates and use
\[
 \int_{T_e}\lambda_0^i\lambda_1^j\lambda_2^k
   = 2 |T_e|\,\frac{i!\,j!\,k!}{(i+j+k+2)!}.
\]
The measure is \(dR_{\rm ref}\,dZ\). The exact binary vertices, full binary reference vector and positive binary scale make \(A_{\rm int}/\pi\) and \(b_{\rm int}/\pi\) rational. A later implementation can apply directed intervals only at the common mathematical-\(\pi\) factor. The existing exact moment and axial stiffness machinery already supplies the material polynomial building blocks. [Moment primitive](../../../../../scripts/exact_fe_moments.py#L102) · [Current exact centered axial derivative](../../../../../scripts/exact_centered_bulk_L.py#L120)

For the real \(Q_{\rm cache}\) potential, the coefficients are exact binary site maps and positive binary weights divided by the same \(S_{\rm binary}\). Replace the exact integral by the finite sum \(A_Q=\sum_i(w_i/S_{\rm binary})\,\chi_{e(i)} B_i^TB_i\) and use the corresponding \(b_Q\). These coefficients are rational and enter the same identity with mathematical drive \(gI\). Binary quadrature weights already contain their stored, rounded cylindrical factor, including the numerical \(\pi\). Applying another \(2\pi R_{\rm ref}\) would double-count it. Exact \(Q_*\) site maps with positive old weights yield their own Gram matrix. The defect \(A_{\rm int}-A_Q\) is symmetric and need not be PSD. The finite \(Q_*-Q_{\rm cache}\) contribution can also be assembled exactly in this part. The real potential with fixed stored tensors is the exact object of the cached proof. Torch sums, coordinates and strength multiplication round in the solver. Sending the rounded strength \(\operatorname{RN}(gI)\) on individual solver calls therefore does not implement an algebraically affine map of exact optical dose at the bit level. [Exact auxiliary rule](../../../../../scripts/exact_fe_target_safe.py#L219) · [Candidate drive consumer](../../../../../scripts/transfer_ambiguity_cached_pilots.py#L314) · [Cached Gram assembly in prior gain certificate](../../../../../scripts/energy_gain_cell.py#L79)

## Consequences for validated continuation

The exact identities are \(\partial_I F=g(Au+b)\), \(\partial_I H=gA\), and \(\partial_I^2 F=0\) at fixed \(u\). All third displacement derivatives of the activation energy vanish. Thus all spatially nonlinear derivative variation and associated approximation tails can be treated in the passive terms. At one fixed \(u\), \(H(u,I)\succeq H(u,I_-)\) for \(I\ge I_-\) and \(g\ge0\). A verified uniform positive lower bound at the left dose over a common nodal enclosure would then hold throughout the product region. Changing \(u\) along an unknown path changes the passive Hessian, so pointwise dose monotonicity gives no automatic lower bound at later states.

On a candidate dose slab \(I=I_c+\delta\), let \(u_{\rm p}=u_c+p\delta\). An exact residual expansion separates the dose dependence:
\[
 F(u_{\rm p},I_c+\delta)=F(u_c,I_c)
   +[H(u_c,I_c)p+g(Au_c+b)]\,\delta
   +gAp\,\delta^2
   +\mathcal R_{\rm pass}(\delta),
\]
where \(\mathcal R_{\rm pass}(\delta)=\int_0^1[H_{\rm pass}(u_c+t p\delta)-H_{\rm pass}(u_c)]p\delta\,dt\). An approximate tangent \(p\) is sufficient when its linear defect is enclosed. For \(|\delta|\le h\), a preconditioner \(P_0\) gives the whole-tube derivative bound by enclosing \(I_n-P_0[H_{\rm pass}(U)+g I_cA]\) and adding \(gh\|P_0A\|\). Signed evaluation of the affine derivative term at dose endpoints may be sharper. Use the same cache or exact-integral operator consistently when forming \(A,b,H,F\), including the corresponding normalization.

To certify one stable quasistatic finite-element path:

1. Anchor a certified root at the starting dose of the specified target. Treat zero displacement as a proposal until exact stationarity or an enclosed correction is proved.
2. Cover the complete dose interval by adjacent slabs. In each slab certify a uniform self-map about the predictor, with force bound \(\eta\), whole-tube derivative contraction \(z<1\), and corrector radius \(r\) satisfying \(\eta+zr<r\). Certify physical geometry and a positive energy Hessian over the appropriate tube. The older energy norm construction is one way. A nonsymmetric approximate inverse requires the separate same-operator positive symmetric-part/Hessian argument if using a contraction proof to infer minima.
3. At a shared dose, show one slab's posterior root enclosure is contained in the other slab's uniqueness domain, or place both root enclosures in a third common uniqueness domain. Mere overlap of cubes or a shared dose is insufficient. Continuity follows from uniform contraction for the continuously moving predictor, and smoothness follows locally from nonsingular \(H\).
4. Match the final certified root to the archived endpoint certificate with the same common-uniqueness argument if it will carry that observer claim. Different geometry/activation targets must each have their own path and normalization.

On such a smooth stable branch \(H u'(I)=-g(Au+b)\). The unit-drive energy \(V\) decreases weakly on the branch because \(dV(u(I))/dI=-g(Au+b)^TH^{-1}(Au+b)\le0\). This concerns a material axial energy, whereas pointwise surface shortening and a maximum over the observed region are different observables. The current accepted discrete checkpoints and singleton root/minimum certificates do not supply the required slab coverage or branch matching. The completed path would concern quasi-static local energy minima of a specified FE model. Time-dependent tissue behavior and continuum approximation require other evidence. [Discrete accepted-state policy](../../../../../docs/QUADRATIC-SOLVE-POLICY.md) · [Previously proved short branch and limits](../../../../../docs/ENERGY-GAIN-GEOMETRY.md#3-a-local-parameter-branch-not-a-biological-activation-history) · [Peak-observer scope](../../../../../docs/CALIBRATION-UNCERTAINTY.md#why-the-floor-can-apply-to-a-quasistatic-peak-model)

