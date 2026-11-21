# The exact-integral root for record L

**The original predictor cube for record L contains a unique root of the exactly integrated finite-element force for every parameter vector in the same closed cell. The root is a strict local minimum of the normalized finite-element energy.** The signed residual calculation establishes this result within the original cube of radius \(6.011722184963558\times10^{-11}\) in normalized nodal coordinates. The displayed numbers are rounded. The strict comparisons use exact fractions. [Independent 256-bit certificate](../runs/exact_integration_cohort_v1/L_defect_diagnostic/signed_predictor_L_v1/actual_independent_review.json) · [Final mathematical readback](../runs/exact_integration_cohort_v1/L_defect_diagnostic/signed_predictor_L_v1/final_mathematical_readback/receipt.json)

## The model and the cube

The broader experiment uses light to induce local contraction in embryonic epithelial tissue. Record L, the workbook column identified as **optoGEF2_small_col_L**, has a dimensionless source dose of **1.51057344**. The calculation concerns a fixed axisymmetric quadratic, or P2, finite-element model of the tissue. Here \(u\in\mathbb R^{611}\) is the vector of free normalized nodal displacements, after imposing the axis and outer-rim constraints. Its parameter vector is \(\theta=(K_a,K_b,T_a,T_b,g)\). These are apical and basal area moduli, apical and basal surface tensions, and common activation gain. Apical and basal name the tissue's two opposing faces; their distinct elastic and contractile properties can affect bending. The active coefficient is the exact product of gain and source dose. [Original L protocol](../runs/exact_integration_cohort_v1/records/optoGEF2_small_col_L/original_cube_v1/protocol.json) · [Shared root cohort](COHORT-UNIFORM-ROOTS.md)

Write \(F_Q=\nabla E_Q\) for the force of the energy evaluated by the frozen binary quadrature cache, and \(F_I=\nabla E_I\) for the force of the energy integrated over the same fixed finite-element space. The force is the energy gradient under the sign convention used here. Both share a positive energy normalization. Write \(H_Q\) and \(H_I\) for their nodal derivatives, the respective energy Hessians.

The cached proof supplied the affine predictor \(p(\theta)=u_0+S(\theta-\theta_0)\) and a square preconditioner \(R\) that approximately inverts the force derivative. A predictor estimates the change in displacement as the parameters move away from the base point \(\theta_0\). The exact same closed parameter cell \(C\), stored \(S\) and \(R\), and parameter-matched cube \(X_\theta=p(\theta)+[-r,r]^{611}\) enter this proof. The cell includes the one-sided apical tension interval starting at zero. A larger independent nodal box covers all these cubes and the segments from \(u_0\) to each predictor. The positive-Jacobian and surface-domain bounds permit the integral and derivative enclosures along those segments. [Original nodal box](../runs/exact_integration_cohort_v1/records/optoGEF2_small_col_L/original_cube_v1/coordinate_box.json) · [Earlier L integration and domain checks](COHORT-INTEGRATION-L.md)

## Center the bulk integration defect

Split the difference \(F_I-F_Q=D_B+D_S\) into bulk and surface force terms in the same free coordinate order. Set \(\delta\theta=\theta-\theta_0\) and let \(\delta g=g-g_0\) be its gain component. The bulk defect \(D_B(u,g)\) is affine in \(g\) at fixed \(u\). Its gain derivative \(G_B(u)=\partial_gD_B(u,g)\) uses the source-dose-scaled difference of the integral and cache axial stiffness actions.

First move from \(u_0\) to \(p(\theta)\) at \(g_0\), then change gain at the final displacement. The fundamental theorem of calculus and the exact affine gain step give

\[
\begin{aligned}
R D_B(p(\theta),g)
  &= b_0 + R\overline H_B(\theta) S\,\delta\theta+R G_B(p(\theta))\,\delta g,\\
b_0&=R D_B(u_0,g_0),\\
\overline H_B(\theta)&=\int_0^1 \partial_u D_B(u_0+tS\,\delta\theta,g_0)\,dt .
\end{aligned}
\]

The computation centers all **144 bulk elements** using a point-force enclosure, a Hessian enclosure applied to the displacement increment, and the exact gain dependence with rigorously bounded derivative. For the **24 surface edges**, it retains the signed force intervals over the whole box. The combined reduction assembles **611 force rows** and **6,676 symmetric upper-triangular Hessian entries**, in the same normalization and coordinate order. The surface term may vary over the box, so its addition to the point bulk term supplies an interval enclosure rather than a surface force evaluated at a unique parameter point. [Centered bulk protocol](../runs/exact_integration_cohort_v1/L_defect_diagnostic/centered_bulk_L_v1/preparation/protocol.json) · [Centered bulk packet review](../runs/exact_integration_cohort_v1/L_defect_diagnostic/centered_bulk_L_v1/actual_independent_review.json) · [Centered reduction review](../runs/exact_integration_cohort_v1/L_defect_diagnostic/centered_bulk_L_v1/reduction_preparation/actual_independent_review.json)

## Preserve signs in the full predictor residual

At fixed \(u_0\), \(F_Q\) is affine in the five coefficients. To build its predictor identity, move the parameters first at \(u_0\), then move the displacement from \(u_0\) to \(p(\theta)\) at the **final** parameter \(\theta\). With \(q_0=R F_Q(u_0,\theta_0)\),

\[
\begin{aligned}
R F_Q(p(\theta),\theta) &= q_0+A_Q\,\delta\theta+\varepsilon_Q(\theta),\\
A_Q &= R\bigl[\partial_\theta F_Q(u_0,\theta_0)+H_Q(u_0,\theta_0)S\bigr],\\
\varepsilon_Q(\theta)
&=R\int_0^1\bigl[H_Q(u_0+tS\,\delta\theta,\theta)-H_Q(u_0,\theta_0)\bigr]S\,\delta\theta\,dt .
\end{aligned}
\]

This order makes the parameter step exact and accounts for both parameter and displacement dependence in the nodal remainder. Each closed parameter-offset interval contains zero, and the verified bounds cover the straight nodal segment at final \(\theta\). The cached-force audit reconstructs the signed anchor and five slopes from **9,504 cache sites** and gives a symmetric error interval for the remainder in each row. [Frozen signed residual protocol](../runs/exact_integration_cohort_v1/L_defect_diagnostic/signed_predictor_L_v1/preparation/protocol.json) · [Prelaunch review](../runs/exact_integration_cohort_v1/L_defect_diagnostic/signed_predictor_L_v1/preparation/prelaunch_independent_review.json) · [Independent cached replay source](../runs/exact_integration_cohort_v1/L_defect_diagnostic/signed_predictor_L_v1/independent_actual/replay.py)

To enclose \(RF_I(p(\theta),\theta)\), add \(q_0\), \(b_0\), and the bounded surface term \(RD_S\) in each row with their signs intact. Form each of the five interval columns by adding the cached slope to \(R\overline H_B S\), and include \(RG_B\) in the gain column. Each resulting column multiplies one shared parameter offset. Add the cached remainder enclosure, then take the largest absolute endpoint over all rows. Keeping the signs until that last operation preserves cancellations between the cached and integration-defect contributions. Evaluation over the closed offset intervals may relax correlations, so the result remains a uniform enclosure over \(C\). [Independent signed rows and exact bounds](../runs/exact_integration_cohort_v1/L_defect_diagnostic/signed_predictor_L_v1/independent_actual/results.json)

![Uniform residual bounds obtained with separate force norms and by adding signed rows, relative to the exact force budget.](../runs/exact_integration_cohort_v1/L_defect_diagnostic/signed_predictor_L_v1/publication/L_signed_force_budget.png)

## Apply the contraction test

Use the infinity norm, taking the largest component in magnitude for a vector and the maximum absolute row sum for a matrix. The proof uses a uniform bound \(\eta\) on \(\|RF_I(p(\theta),\theta)\|_\infty\) and a uniform bound \(z\) on \(\|I-RH_I(u,\theta)\|_\infty\) for \(u\in X_\theta\). The whole-cube derivative bound is \(z=Z_Q+b\), where \(Z_Q\) bounds \(\|I-RH_Q\|_\infty\) and \(b\) bounds \(\|R(H_I-H_Q)\|_\infty\). The signed calculation sharpens \(\eta\) while using the previously checked whole-cube derivative bound. [Scaled reduction review](../runs/exact_integration_cohort_v1/L_defect_diagnostic/scaled_offaxis_L_v1/reduction_preparation/actual_independent_review.json) · [Independent signed proof](../runs/exact_integration_cohort_v1/L_defect_diagnostic/signed_predictor_L_v1/actual_independent_review.json)

For fixed \(\theta\), take \(T_\theta(u)=u-RF_I(u,\theta)\). The derivative bound gives Lipschitz constant at most \(z\) on the convex cube. At its center the map moves the predictor by at most \(\eta\). Hence every image lies within distance \(\eta+zr\) of the predictor. The decision requires \(z<1\) and \(\eta+zr<r\). The original direct integral enclosure leaves this test unresolved. Rewriting the off-axis integrals in scaled material coordinates and then centering all bulk force contributions further tighten the bounds. The signed full residual is the first of the four stages to satisfy the strict inequality. [Direct reduction](../runs/exact_integration_cohort_v1/reduction_preparation/optoGEF2_small_col_L_reduction_review.json) · [Centered reduction](../runs/exact_integration_cohort_v1/L_defect_diagnostic/centered_bulk_L_v1/reduction_preparation/actual_independent_review.json) · [Exact arithmetic for both figures](../runs/exact_integration_cohort_v1/L_defect_diagnostic/signed_predictor_L_v1/publication/receipt.json)

| Enclosure stage | \((\eta+zr)/r\) |
|---|---:|
| Direct integration defect | 10,917.7892967 |
| Scaled coordinates on 132 off-axis elements | 25.4985717175 |
| Centered bulk with separate cached and defect norms | 1.19214056427 |
| Signed full predictor residual | **0.765587093559** |

![The ratio of the self-map bound to the original cube radius at each of the four reviewed stages, with the pass threshold at one.](../runs/exact_integration_cohort_v1/L_defect_diagnostic/signed_predictor_L_v1/publication/L_self_map_stages.png)

The final numbers are from the **independent 256-bit replay**. The predictor coordinate bounds, force-map residual and map radius use normalized nodal coordinates.

| Quantity | Rounded value |
|---|---:|
| Signed predictor residual bound \(\eta\) | \(2.298807859837408\times10^{-11}\) |
| Whole-cube Lipschitz bound \(z\) | \(0.38319951989736134\) |
| Original radius \(r\) | \(6.011722184963558\times10^{-11}\) |
| Self-map bound \(\eta+zr\) | \(4.602496914871760\times10^{-11}\) |
| Strict self-map margin \(r-(\eta+zr)\) | \(1.409225270091799\times10^{-11}\) |
| Predictor-to-root error bound \(\eta/(1-z)\) | \(3.726987792640620\times10^{-11}\) |
| Normalized Hessian curvature lower bound \(m\) | \(3.948410451880172\times10^{-4}\) |

The closed cube is complete, so Banach’s theorem supplies a unique fixed point there. Since \(\|I-RH_I\|_\infty<1\), the matrix \(RH_I\) is invertible by the Neumann series. Both factors are square, so \(R\) is invertible and a fixed point is exactly a force root. If \(u_I\) is the root, the contraction estimate gives \(\|u_I-p(\theta)\|_\infty\le\eta+z\|u_I-p(\theta)\|_\infty\), yielding the tabulated error bound. Every bound applies uniformly on the same closed \(C\), with one selected cube for each parameter.

For local minimality, the cached Hessian \(H_Q\) is positive definite on the same tube by the reviewed energy certificate. At a fixed \(u,\theta\) in the tube, the symmetric homotopy \(H_t=(1-t)H_Q+tH_I\) satisfies \(\|I-RH_t\|_\infty\le (1-t)Z_Q+tz<1\) for all \(0\le t\le1\). Thus \(H_t\) is nonsingular and its inertia, the counts of positive and negative eigenvalues, stays constant. This carries positive definiteness to \(H_I\). The Neumann bound on \((RH_I)^{-1}\) gives \(\|H_I^{-1}\|_\infty\le\|R\|_\infty/(1-z)\). Since the spectral radius is bounded by the induced norm, \(\lambda_{\min}(H_I)\ge m=(1-z)/\|R\|_\infty>0\). The interior root is therefore a strict local energy minimum. [Cached energy certificate](../evidence/cohort_energy_stability_independent_review.json) · [Final proof readback](../runs/exact_integration_cohort_v1/L_defect_diagnostic/signed_predictor_L_v1/final_mathematical_readback/receipt.json)

## Scope and reproducibility

The earlier B theorem places its integral root outside the original B cube, whereas this theorem puts the L root inside the original L cube. These comparisons depend on each record’s numerical predictor and chosen radius. The physical dose effect cannot be inferred from cube membership. B and L now have exact-integral root results for their respective selected branches. For the other fifteen records, integral root certification remains open, as do continuum PDE approximation, source-compatible physical fit, physical loading paths, and novelty. [B exact-integral theorem](EXACT-INTEGRAL-ROOT.md) · [Five-gate study contract](../PLAN.md)

The independent replay reconstructs the cached 9,504-site calculation and the 611 signed force rows. It shares the reviewed local derivative, integral-defect and box-contraction premises, which the final mathematical readback binds explicitly. The retained visualization source rehashes those bindings and derives both plots from the exact certificates. [Frozen source and preflight](../runs/exact_integration_cohort_v1/L_defect_diagnostic/signed_predictor_L_v1/preparation/preflight.json) · [Independent terminal receipt](../runs/exact_integration_cohort_v1/L_defect_diagnostic/signed_predictor_L_v1/independent_actual/terminal.json) · [Figure source](../runs/exact_integration_cohort_v1/L_defect_diagnostic/signed_predictor_L_v1/publication/figures.py) · [Figure protocol and source hashes](../runs/exact_integration_cohort_v1/L_defect_diagnostic/signed_predictor_L_v1/publication/protocol.json)
