# Gate 3: mechanically defensible model candidates

**Decision:** a finite-height independent-cell null can be specified and tested quantitatively now. It cannot predict patch-scale bending. A coupled finite-thickness model is the appropriate shape candidate, but its large-region point prediction is not identified by a scalar small-region shortening curve alone. No mechanical simulation, fit or held-out score is performed in this design document. None of the standard energy, scaling or invariance arguments below is claimed novel.

## Experimental and input contract

Countryman et al., *Nature Communications* 16:7665 (2025), [primary article](https://doi.org/10.1038/s41467-025-62483-6), Figs. 2–4, reports lateral illumination of an 18-micrometre-diameter region at cell heights 32–35 micrometres. Illumination occupies 5–25 micrometres below the apical surface. Its dose is an optical recruitment proxy, not a force measurement. The source distinguishes cell shortening from apical-surface furrowing and from internalization. Its supplementary information pp. 9–10 and 13 discusses optical controls, larger-region morphology and tool limitations.

Allowed training outcomes are the small-region observations only. Larger-region diameter and measured optical dose are prescribed intervention inputs; larger-region shortening, bending labels, surface profiles and representative post-activation images are evaluation outcomes. They cannot select stiffness, active gain, a threshold, an initial shape or an error tolerance. Published outcomes are already known, so the test is retrospective, not blind.

The source-data audit reports 17 small-region maximum-shortening/dose pairs and nine nominal-50-micrometre pairs in Fig. 4h. The repeated small pairs in Fig. 2d/4h are not independent training cohorts. The continuous target must be **maximum shortening in both groups**. Fig. 2c average shortening is a different observable. Separate apical/basal maxima need not occur at the same place or time; their sum is not constrained to equal the maximum cell-height shortening. Report measured diameters separately from nominal labels. No specimen identity should be invented from row order.

The published activation normalizer is calculated across experiments. Using supplied normalized doses is a declared source optical calibration, not a strictly training-only preprocessing claim. Re-estimating a training-only normalizer requires raw activation intensities. No normalization may use held-out shape outcomes.

## Thickness and curvature invalidate an automatic plate reduction

Let activation radius be \(a=D/2\) and tissue height \(h\). At nominal diameter 18 micrometres,

\[
\frac h a=\frac{32\text{--}35}{9}=3.56\text{--}3.89.
\]

At nominal diameter 50 micrometres, \(h/a=1.28\text{--}1.40\). Neither is a thin-tissue limit \(h/a\ll1\). A thin basal cortical layer could have thickness much smaller than \(a\); that does not make the entire epithelium a thin plate. A plate law needs a separately justified reduction and comparison against finite-thickness mechanics.

Initial AP and DV curvatures must be measured before activation. They are independent geometric inputs, not parameters selected from the large-region bending response. Axisymmetry is a conditional circular-patch benchmark; anisotropic initial curvature or elliptical stimulation requires a genuinely non-axisymmetric model. Cell size relative to activation radius must also be checked before treating a few-cell patch as a continuum.

## M0: exact finite-height column null

Use an ideal circular column solely as a homogeneous local null. Its volume \(V\) is fixed, cross-sectional area is \(A>0\), height is \(H=V/A\), and lateral area is \(2V\sqrt\pi A^{-1/2}\). For nonnegative basal area modulus \(K_b\), positive lateral tension \(\sigma_l\) and positive combined apical/basal baseline tension \(\Sigma\), define

\[
E(A;I)=\Sigma A+\frac{K_b}{2A_0}(A-A_0)^2
+2\sigma_l(I)V\sqrt\pi A^{-1/2},
\qquad \sigma_l(I)=\sigma_{l0}+gI.
\]

The positive affine gain \(g\) is an explicit small-dose constitutive hypothesis, not a measured conversion from fluorescence to force. Positive tension imposes a restricted dose range. This energy structure follows existing cell-surface mechanics; it does not propose a new mechanism.

The equilibrium law is

\[
\Sigma+K_b(A/A_0-1)-\sigma_l(I)V\sqrt\pi A^{-3/2}=0,
\qquad S(I)=H_0-\frac V{A(I)}.
\]

It is well defined: energy diverges at zero and infinity under the stated assumptions, and

\[
E''(A)=\frac{K_b}{A_0}+\frac{3\sigma_l(I)V\sqrt\pi}{2A^{5/2}}>0.
\]

Thus exactly one positive minimizer exists. A positive-mobility overdamped relaxation is permissible if kinetics are explicitly needed, but endpoint maxima do not determine a unique activation time course or viscosity. The fixed-volume null is conditional, not a proof that cellularizing tissue has constant volume. Cellularization reference growth belongs in the observation operator rather than being misclassified as activation-induced shortening. The source's post-response reference correction may define the evaluated target; it cannot silently become a future-time forecast input.

### Independent-column invariance theorem

Suppose a patch contains \(N\) identical, uncoupled columns, all with the same dose and no diameter-dependent boundary load. Its energy is \(\sum_{i=1}^N E(A_i;I)\). Strict convexity gives \(A_i=A(I)\) for every cell, for every positive \(N\). Consequently both the mean and the maximum shortening equal \(S(I)\), independently of patch diameter. Multiplying the number of columns cannot generate a bending transition or reverse the displacement direction.

This is a rigorous statement about M0, not a theorem that all local biological models are size invariant. In particular, for independent **heterogeneous** cells with shortening CDF \(G_I\), the maximum satisfies

\[
\Pr(M_N\le s\mid I)=G_I(s)^N.
\]

Maximum shortening may therefore change with patch size even without mechanical coupling. A claim that a size effect proves mechanical propagation must exclude this order-statistics alternative. The workbook does not supply the cell-level distribution and activated cell counts needed to calibrate it. Scalar maximum shortening is a bounded response test, not a substitute for predicting surface geometry.

The M0 energy has no tissue embedding or center-of-mass constraint. It determines height, not the individual apical/basal positions in the laboratory frame. It cannot predict an outward bulge or an inward furrow. A classifier appended to its shortening curve would not be a mechanical derivation.

## M1: finite-height asymmetric-cell energy

The closest theoretical baseline is Guerra Santillán, Dahmann and Fischer-Friedrich, *PRX Life* 2:013004 (2024), [main paper](https://doi.org/10.1103/PRXLife.2.013004), Eq. (1), Fig. 6 and [Supplement §I, S1–S7](https://journals.aps.org/prxlife/supplemental/10.1103/PRXLife.2.013004/Supplementary_GuerraSantillan.pdf). It already couples lateral tension, basal elasticity and cell-shape asymmetry. Its wing-disc basement-membrane parameters cannot be assigned to embryonic basal actomyosin. The material structure is a candidate analogy requiring validation.

An exact geometric extension, not a validated new biological law, uses a regular hexagonal frustum. Let top/bottom edge lengths be \(w+\delta\) and \(w-\delta\), with \(w>|\delta|\), and \(C=3\sqrt3/2\). Then

\[
A_a=C(w+\delta)^2,\quad A_b=C(w-\delta)^2,
\quad H=\frac V{C(w^2+\delta^2/3)},
\quad A_l=6w\sqrt{H^2+3\delta^2}.
\]

The volume and lateral area are exact frustum formulas, including finite slant. There is no \(H/w\ll1\) assumption. Define

\[
E(w,\delta;I)=\sigma_a A_a+\sigma_{b0} A_b
+\frac{K_b}{2A_{b0}}(A_b-A_{b0})^2+\sigma_l(I)A_l.
\]

This is still a local cell model. It does not make independently optimized frusta into a compatible tissue: shared faces, cell tilts and neighboring-cell forces must be solved jointly. Basal/apical stiffness differences are capable of inducing shape asymmetry, but that principle and its weak-curvature analysis are already in the cited theory.

### An explicit local identifiability calculation

For an illustrative symmetric reference, normalize \(w_0=1\), set \(A_{b0}=C\), \(\sigma_a=\sigma_{b0}=T>0\), and let \(H_0=q\). The resting lateral tension is \(\sigma_{l0}=\sqrt3T/q\). The exact Hessian at the resting configuration is

\[
\mathcal H_{ww}=4C(3T+K_b),\qquad
\mathcal H_{w\delta}=-4CK_b,\qquad
\mathcal H_{\delta\delta}=2\sqrt3(T+3K_b+9T/q^2).
\]

Small lateral-tension increments produce

\[
\frac{\Delta\delta}{\Delta w}
=\frac{3K_b}{T+3K_b+9T/q^2}\in[0,1).
\]

The height-response slope fixes a combination of active gain and effective stiffness; it does not, by itself, fix the asymmetry ratio. For any admissible \(K_b/T\), active gain can be adjusted to reproduce the same infinitesimal shortening slope while changing predicted asymmetry. This is a local small-signal non-identifiability result for this parameterization, not a proof that every nonlinear dose curve or every additional apical/basal observation fails to identify it.

No global tissue curvature formula may be inferred merely by treating these cells as independent wedges. A parallel-surface relation such as \(\kappa\simeq2\delta/(wH)\) additionally assumes a compatible weak-curvature embedding. It is a diagnostic limit, not a thin-plate model or a large-region shape predictor.

## M2: minimal genuine shape candidate

The defensible mechanical extension is a coupled finite-thickness patch with explicit apical and basal interfaces. Two implementations are admissible:

* A 3D cell/vertex model with fixed cell volumes, shared lateral interfaces, apical/basal tensions, basal area elasticity and surrounding nonactivated cells. Activation changes specified lateral face tensions within the illuminated 5–25-micrometre depth range. Cell topology and initial surface geometry come from pre-activation observations. Shared faces must be counted once with an explicitly defined interface tension.
* A 3D nearly incompressible hyperelastic/viscoelastic body with resolved depth-dependent active stress and apical/basal surface constitutive laws. This is a homogenized approximation; it must be checked against activation radius/cell size and cannot silently become a Kirchhoff–Love plate. Mixed incompressibility, volume conservation, baseline prestress, contact and outer-boundary sensitivity need independent checks.

Both must embed the pre-activation tissue curvature, mechanically couple activated cells to an unactivated surrounding region, and account for vitelline contact if relevant. The outer domain and boundary conditions are fixed from the acquisition geometry or a declared remote-boundary sensitivity analysis. They cannot be tuned against the large-region morphology. Basal elasticity is a material assumption, not a generic clamp: displacement support, area stiffness and active prestress are different quantities.

No arbitrary diameter threshold is allowed. Continuous changes in curvature, apical depth and basal position must result from equilibrium/dynamics. Any bending label is applied afterward using an independent, predeclared observation operator; a label threshold chosen from held-out Fig. 4f is forbidden.

## Dimensionless combinations and what can actually be fitted

For the local surface-energy model, using \(w_0\) and \(\sigma_a w_0^2\) as scales leaves the ratios

\[
V/w_0^3,\quad A_{b0}/w_0^2,\quad
\sigma_{b0}/\sigma_a,\quad K_b/\sigma_a,
\quad\sigma_{l0}/\sigma_a,\quad g/\sigma_a.
\]

For a bulk 3D model with shear modulus \(\mu\), height \(h\) and time scale \(t_*\), relevant groups include

\[
a/h,\quad h/R_{AP},\quad h/R_{DV},\quad\ell_{cell}/h,
\quad K_b/(\mu h),\quad T_{a,b}/(\mu h),\quad
\sigma_{active}/\mu,\quad\eta/(\mu t_*).
\]

Illumination depths divided by \(h\), outer-domain size divided by \(h\), contact properties and activation turnover are additional inputs if used. A common positive scaling of all force/energy parameters, with matching viscous/friction parameters for dynamics, leaves shape trajectories invariant. Shape data therefore constrain ratios, not an absolute force scale, unless an independent mechanical measurement breaks that symmetry.

## Identified-set prediction contract

For a fully specified candidate, let \(Y_{18}(\theta)\) denote training observables and \(\mathcal U_{18}\) their declared uncertainty region. Use

\[
\Theta_{18}=\{\theta:\ Y_{18}(\theta)\in\mathcal U_{18}\}
\cap\Theta_{\rm physical},
\qquad
\mathcal P_D=\{Y_D(\theta):\theta\in\Theta_{18}\}.
\]

Only after freezing \(\Theta_{18}\), the observation operator and numerical uncertainty may a larger-diameter experiment be scored. A point estimate with a narrow optimizer covariance is not proof that \(\mathcal P_D\) is narrow. If several parameter sets fit the small-region data but give materially different large-region surface shapes, the correct prediction is that set or interval. If even its envelope cannot be computed reliably, the mechanical point-prediction gate remains unresolved.

The scalar workbook can support a no-refit homogeneous-column response null. A credible 3D prediction additionally requires initial cell geometry and curvature, a justified basal constitutive law, compatible neighboring-tissue constraints, and spatial apical/basal validation data or a sufficiently informative identified set. Do not invent these quantities to force a desired large-region outcome.

## Required verification before a coupled prediction claim

1. Check exact cell volume, surface areas, rigid-motion invariance and zero-activation equilibrium.
2. Check stable passive energy and nonnegative dissipation for the selected rheology; active work is explicitly supplied, not free equilibrium energy generation.
3. Verify finite-thickness solutions and coupled-face forces; compare any reduced limit with the full model.
4. Recover known parameters/combinations synthetically and exhibit alternative fits along unresolved directions.
5. Freeze the small-region calibration, activation observation model and numerical error protocol before large-region scoring.
6. Report scalar shortening, surface displacement and bending-classification outcomes separately. A scalar response success cannot pass the 3D shape-prediction gate.

The accompanying `model_design_*.json` receipts retain source hashes, exact-symbolic checks and the design boundary. Modest originality, if any, must lie in a genuinely new and independently checked no-refit prediction or nontrivial information limit, not in relabeling this standard mechanics.
