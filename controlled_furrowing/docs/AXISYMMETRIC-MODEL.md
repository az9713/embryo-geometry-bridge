# A coupled finite-height axisymmetric solid

This component implements a three-dimensional elastic body restricted to axisymmetric motion. It is a conditional mathematical model of a circular activation patch coupled to surrounding tissue. It is not a plate, a specimen reconstruction, a calibrated biological constitutive law, or a prediction of the observed transient bending label. No experimental outcome is read by the mechanics module.

The construction, constitutive choices, automatic differentiation and continuation are standard methods. This implementation makes no originality claim.

## Reference geometry and the actual three-dimensional deformation

Let the material meridian have cylindrical coordinates \((r,z)\), with \(0\le r\le R_o\). A flat reference body has \(0\le z\le H\). Optional reference curvature replaces its boundaries by \(z=k r^2/2\) and \(z=H+k r^2/2\), represented piecewise linearly. Both geometries are explicitly stress-free idealizations. They do not supply measured prestress or reconstruct the embryo's initial geometry.

Axisymmetric motion without swirl maps \((r,\theta,z)\) to \((R(r,z),\theta,Z(r,z))\). In the cylindrical orthonormal bases,

\[
 F=\begin{pmatrix}R_{,r}&0&R_{,z}\\0&R/r&0\\Z_{,r}&0&Z_{,z}\end{pmatrix},\qquad
 F_m=\begin{pmatrix}R_{,r}&R_{,z}\\Z_{,r}&Z_{,z}\end{pmatrix},\qquad
 J=\det F_m\,R/r.
\]

In particular the hoop stretch \(R/r\) is retained. The integration measure is the three-dimensional reference volume \(dV_0=2\pi r\,dr\,dz\). The displacement is solved through the whole finite height and throughout the nonactivated annulus; these regions share finite-element nodes. There is no independent-cell or prescribed-displacement shortcut.

At the symmetry axis, \(R=0\). At the outer cylindrical boundary, both components of displacement vanish. This is a remote full clamp, not an experimentally established basal vertical clamp. Both apical and basal surfaces move freely subject to their elastic surface energies. Axisymmetry cannot predict an elliptical patch or distinguish the embryo's AP and DV curvatures.

## Passive and active energies

The passive bulk energy density is

\[
 W(F)=\frac{\mu}{2}\bigl(\operatorname{tr}(F^TF)-3-2\log J\bigr)
       +\frac{K}{2}(\log J)^2,\qquad \mu,K>0,
\]

with \(J>0\). Here \(K\) is the coefficient of the logarithmic volume penalty, not an independently measured embryo bulk modulus. Large \(K/\mu\) penalizes volume change but does not impose exact incompressibility. Around identity the infinitesimal bulk modulus is \(K+2\mu/3\).

For reference top/bottom surface area \(dA_0\), let \(J_s=dA/dA_0\). Each surface contributes

\[
 E_s=\int_{\Gamma_{s,0}}\frac{K_s}{2}(J_s-1)^2\,dA_0.
\]

For a straight meridional edge, \(J_s=(R/r)(\ell/\ell_0)\). Thus the basal modulus changes resistance to surface-area change; it does not pin the basal surface vertically. The model has no baseline surface tension, spontaneous curvature or prestress.

Let \(I(r,z)\) be the material activation indicator, supported inside radius \(a\) and a specified fraction of the reference height. The prescribed reference fiber is vertical in the flat case and \(n_0=(-kr,0,1)/\sqrt{1+k^2r^2}\) in the curved idealization. It is a declared constitutive direction, not an inferred measured cell orientation. For an externally controlled nonnegative gain \(g\),

\[
 E_a[y;g]=\int_{\Omega_0}\frac{gI}{2}|Fn_0|^2\,dV_0,
 \qquad
 P_a=gI(Fn_0)\otimes n_0,
 \qquad
 \sigma_a=\frac{gI}{J}(Fn_0)\otimes(Fn_0).
\]

The positive fiber tension favors shorter activated fibers and is invariant under a superposed rigid rotation. It is one possible idealized active-stress closure. It does not model measured Rho activation, ATP consumption or photon-to-mechanical conversion. Its tension changes with stretch; it is not an imposed strain-independent Cauchy stress.

The total mechanical potential is \(\mathcal E=\int W\,dV_0+E_a+E_{\rm basal}+E_{\rm apical}\). At each fixed gain the solver seeks a stationary admissible state by minimizing this discrete potential. No global minimum, uniqueness or global stability theorem is asserted.

For completeness the passive first Piola stress is

\[
 P_p=\mu(F-F^{-T})+K\log J\,F^{-T}.
\]

The weak equilibrium condition for every admissible axisymmetric virtual displacement \(v\), vanishing on the clamp and respecting the axis, is

\[
 \int_{\Omega_0}(P_p+P_a):\delta F[v]\,dV_0
 +\sum_{s\in\{a,b\}}\int_{\Gamma_{s,0}}K_s(J_s-1)\,\delta J_s[v]\,dA_0=0.
\]

The hoop component of \(\delta F\) is \(v_R/r\); omitting it would change the equilibrium. Automatic differentiation takes the variation of the discretized three-dimensional potential, including this term. The free nodal force residual is exactly this discrete first variation in the free coordinate basis, normalized by \(\mu V_0\).

## Discretization and a useful exact orientation check

The meridian is triangulated with continuous piecewise affine displacement. The radial grid contains \(a\), and the height grid contains both activation-depth boundaries, so an activated triangle does not straddle the prescribed discontinuity. The meridional gradient is constant on a triangle. The hoop stretch varies as a ratio of affine functions. Triangle quadrature can use degrees 1, 2 or 4; independent Gauss surface quadrature uses 1, 2, 3, 5 or 8 points. The default uses degree 2 in the volume and 3 surface points.

For barycentric coordinates \(\lambda_i\) in a triangle with \(r>0\), axis nodes satisfy \(R_i=r_i=0\), and

\[
 \frac{R}{r}=\frac{\sum_i\lambda_i R_i}{\sum_i\lambda_i r_i}
 =\sum_{i:r_i>0}\frac{\lambda_i r_i}{\sum_j\lambda_j r_j}\frac{R_i}{r_i}.
\]

The coefficients are nonnegative and sum to one. Consequently the elementwise extrema of \(R/r\) lie between its non-axis vertex ratios. Multiplying both endpoint ratios by the constant \(\det F_m\), and taking their minimum/maximum, gives exact elementwise \(J\) bounds in real arithmetic, including axis limits. The code tests these bounds rather than asserting orientation from quadrature samples alone. Floating point calculations are not interval certificates. Positive local orientation does not prove global injectivity or absence of self-contact.

The current volume is \(\int J\,dV_0=2\pi\sum_T\det(F_{m,T})\int_T R\,dr\,dz\). Its affine integrand makes the discrete volume exact in real arithmetic for all supported triangle quadratures. The nonlinear elastic energy and surface energy have quadrature error.

The optimizer uses bounded L-BFGS followed by automatic-differentiation Hessian/Newton polishing with backtracking. A temporary extension penalizes inverted trial states; its penalty must be exactly zero on every reported physical solution. The accepted convergence gate requires the normalized free nodal force infinity norm below \(10^{-7}\), positive elementwise orientation and stationary continuation stages. Raw L-BFGS success flags remain in every record: a short L-BFGS stage can hit its iteration limit and then pass the independent force gate after Newton polishing.

This is a displacement-only P1 penalty method. Near-incompressible locking and insufficient local resolution are real limitations. Small residuals show that the chosen discrete equations were solved; they do not bound continuum discretization error.

## Outputs, observation limits and units

The result retains every reference/current node and triangle, free displacement vector, material activation mask, full apical/basal \((R,Z)\) boundary profiles, signed central displacements, volume, elementwise \(J\) bounds, sampled stretches, energies and continuation diagnostics. If each boundary is radially monotone, the code constructs a common Eulerian radial grid containing every top/bottom breakpoint. The two piecewise-linear graphs can then be passed without arbitrary crest thresholds to the separate slice observer.

The reported material-column vertical height \(Z_{\rm top}(r)-Z_{\rm bottom}(r)\), the distance between those material endpoints, and the Eulerian graph thickness \(a(R)-b(R)\) are different observables. The code labels them separately. `maximum_activated_column_shortening` is a static equilibrium proxy over material columns. The experimental statistic is a maximum over time; equating these requires an independent observation/kinetic argument and is not established here.

Default units have \(H=1\), \(\mu=1\). Length outputs multiply by a declared physical height, stress by a declared shear modulus, surface moduli have units \(\mu H\), and energy has units \(\mu H^3\). No time unit is supplied. For the parent's exploratory transfer family, \(H=33.5\,\mu\mathrm m\), \(a/H=9.35512037/33.5\), the illuminated material-depth fractions are \((H-25)/H\) to \((H-5)/H\), \(R_o/H=3\), \(K/\mu=50\), \(K_a/(\mu H)=0.2\) and \(K_b/(\mu H)\in\{0.2,2,20\}\). These are declared scenarios, not a measured confidence region or identified parameter set. In particular published discrete basal/apical spring ratios cannot be copied directly into these continuum ratios.

An equilibrium graph may permit an apical slice with a disconnected annulus, or exclude that topology within this model. Without the actual slice plane, specimen geometry and a time-to-state map this is a conditional existence test. The source movie's early plug and later recovery cannot be predicted by a single static equilibrium. Adding a post hoc timing choice would not repair this gap.

## External activation work is a declared bookkeeping convention

At an activation increment with the old configuration held fixed,

\[
 W_j=\mathcal E(y_{j-1},g_j)-\mathcal E(y_{j-1},g_{j-1}),\qquad
 D_j=\mathcal E(y_{j-1},g_j)-\mathcal E(y_j,g_j).
\]

Telescoping gives \(\mathcal E(y_n,g_n)-\mathcal E(y_0,g_0)=\sum_jW_j-\sum_jD_j\). The receipt reports each term and the balance residual. \(D_j\ge0\) is checked rather than presumed for the numerical path. This is an identity for the chosen potential and idealized energy-decreasing relaxation, not a measured dissipative dynamics. Adding an arbitrary gain-dependent constant to the potential changes \(W_j\) without changing forces. Therefore this ledger cannot calibrate chemical or optical energy, establish a finite energy reserve, or supply thermodynamic efficiency.

## API and replay

```python
from axisymmetric_mechanics import Parameters, build_model
model = build_model(Parameters(nr=8, nz=4))
first = model.relax(0.2)  # one fixed-strength solve, no zero reset
second = model.relax(0.5, first['free_displacement'])
# For ramp-work bookkeeping with a declared previous gain:
ramped = model.solve(0.5, initial=first['free_displacement'],
                     continuation_start_strength=0.2,
                     continuation_steps=1, maxiter=100)
```

`relax` starts its work ledger at the same fixed strength; it does not pretend to include the previous activation increment. A normal `solve` starts at zero gain by default and performs continuation. All returned vectors are directly reusable only for the same mesh/material reference. A caller should set Torch's thread count appropriately; the CLI self-test uses one thread.

```powershell
python -B controlled_furrowing/scripts/axisymmetric_mechanics.py --self-test
```

## Verification receipts and limitations

The premeasurement thresholds are retained in `evidence/axisymmetric_protocol.json`. Checks include zero-dose identity, an exact isochoric affine deformation, its analytic bulk/surface energy, directional derivatives checked by finite differences, positive-J activated equilibria, an external-potential ledger, uniform material/activation scaling and independent mesh/quadrature/domain/basal-modulus/reference-curvature sensitivity runs.

The first run failed its force tolerance because L-BFGS stopped on relative objective reduction. It remains intact in `runs/axisymmetric/initial_failure/`. Newton polishing was added; no force threshold was relaxed. The intermediate successful version and receipt remain in `runs/axisymmetric/newton_verified_v1/`. `evidence/axisymmetric_verification.json` records the final replay and source hash. `evidence/axisymmetric_warm_start.json` records the fast warm-start test; its exact older source hash is retained, and subsequent changes exposed surface quadrature without changing defaults.

At gain 0.2 in the synthetic default geometry, shortening changes from approximately 0.010995H on 8×4 cells to 0.013587H on 12×6 and 0.016541H on 16×8. All are stationary to about \(10^{-11}\) in normalized nodal force. This is substantial unresolved mesh sensitivity, not numerical convergence evidence. It blocks a rigorous continuum prediction or a narrow error bar from this implementation alone. Quadrature, outer clamp, basal modulus and curvature sensitivity must also accompany any exploratory calibration.

The fast warm sequence on the parent's small-radius geometry reached gain 4 with positive J and force residuals at most \(4.6\times10^{-12}\), taking approximately 0.3–0.7 seconds per solve on the test host. A long cold continuation and short warm solve at gain 0.5 agreed to \(1.28\times10^{-11}H\). This is a local-branch implementation check, not a theorem ruling out other equilibria.

The component is a genuine coupled finite-height mechanical calculation. Completion of its synthetic checks does not complete biological calibration, a withheld perturbation prediction, the five research gates, or the whole embryo project.

![Synthetic coupled mechanics, dose response, unresolved mesh sensitivity and declared potential work](../runs/axisymmetric/axisymmetric_summary.png)

The upper-left panel displays the actual computed deformation without displacement magnification. The other panels distinguish a stationary discrete solution from continuum convergence and physical energy calibration. Input hashes for this figure are retained in `evidence/axisymmetric_visual.json`. The final-source warm replay is retained separately in `evidence/axisymmetric_warm_start_final.json`.
