# A separate flat-reference surface-tension extension

The initial area-elastic model failed the small-region directional adequacy check: the three explored basal-modulus branches predicted predominantly basal motion, whereas the source's small-region observations predominantly showed apical motion. This failure is retained. The present revision is motivated by that small-region discrepancy alone. It does not use larger-region outcomes to choose a constitutive change.

`scripts/axisymmetric_prestress.py` subclasses the frozen original mechanics module. It adds actual basal and apical surface tensions while retaining the full finite-height bulk, area elasticity, activation potential, mechanical coupling, mesh, outputs and warm-start API. The original source is unchanged; its required SHA256 is `ca431e30a2574a02f4a21c98935c50217aa79bc5512cf7279143f22ce1a580e0`.

## Area elasticity and baseline tension do different things

For each material surface \(\Gamma_{s,0}\), retain its original area-elastic term and add

\[
 E_s=\int_{\Gamma_{s,0}}\left[\frac{K_s}{2}(J_s-1)^2+T_sJ_s\right]\,dA_0,
 \qquad T_s\ge0.
\]

Thus the new term is exactly \(T_s A_s\), where \(A_s\) is the current surface's true area. It is not a vertical spring or an imposed basal clamp. \(K_s\) controls the change of tension with area, while \(T_s\) supplies baseline surface tension even when \(J_s=1\). Its first variation contributes \(T_s\delta A_s\). Uniform surface tension resists out-of-plane slopes even when the planar reference area is initially unchanged.

The stored energy is the nonnegative actual-area convention \(T_s A_s\), including \(T_s A_{s,0}\). Subtracting this activation-independent reference constant would leave all forces and activation-work increments unchanged. This is distinct from adding an activation-dependent energy gauge, which changes the externally controlled potential-work ledger without changing mechanical forces. Neither convention calibrates ATP or optical energy.

## Why the flat initial state is still an equilibrium

Consider either initially flat surface, parameterized by material radius \(0\le r\le R_o\). An axisymmetric variation has radial component \(v(r)\) and normal component \(w(r)\). Its current area under variation amplitude \(\varepsilon\) is

\[
 A(\varepsilon)=2\pi\int_0^{R_o}
 (r+\varepsilon v)\sqrt{(1+\varepsilon v')^2+\varepsilon^2(w')^2}\,dr.
\]

The symmetry axis gives \(v(0)=0\); the fixed outer rim gives \(v(R_o)=w(R_o)=0\). Therefore

\[
 A'(0)=2\pi\int_0^{R_o}(v+rv')\,dr
       =2\pi[rv]_0^{R_o}=0.
\]

The initial plane has zero curvature. Its uniform surface tension has zero distributed normal force and its boundary line force is balanced by the clamped rim. The passive bulk and original area-elastic energies also have zero first variation at the identity, so the total zero-activation reference is stationary. This is an equilibrium statement under the declared clamp, not a claim that the embryo is unstressed.

Expanding one order further gives

\[
 A''(0)=2\pi\int_0^{R_o}\left(2vv'+r(w')^2\right)\,dr
       =2\pi\int_0^{R_o}r(w')^2\,dr\ge0.
\]

The endpoint term \(2\pi[v^2]_0^{R_o}\) vanishes. Hence true baseline tension adds precisely a nonnegative surface-normal gradient stiffness around the flat state. It does not add a spurious tangential area stiffness for planar rearrangements with the same rim. This also holds for the P1 discrete boundary: edge contributions telescope at shared vertices.

A curved reference would have nonzero curvature force from uniform tension. Without balancing pressure, bulk prestress or external tractions, it would generally fail zero-dose equilibrium. The new parameter class therefore rejects all nonzero reference curvature. It does not silently relax a curved initial shape and call that state measured or stress-free.

## Exact area implementation and scaling

For a straight meridional edge with endpoint radii \(R_1,R_2\) and current length \(\ell\), the revolved frustum area is

\[
 A_{\rm edge}=\pi(R_1+R_2)\ell.
\]

The implementation integrates \(J_s\,dA_0=2\pi R\ell\,ds\) with the inherited edge quadrature. Its integrand is affine in the edge parameter, so any supported Gauss rule computes this tension contribution exactly in real arithmetic. The self-test independently checks it with the frustum formula. The original nonlinear area-elastic contribution still has quadrature error.

If every stress or surface-modulus parameter and the activation strength are multiplied by \(c>0\), including \(T_a,T_b\), then the entire energy is multiplied by \(c\). Its stationary geometry is unchanged. The normalized objective also stays unchanged because its denominator is \(\mu V_0\). The synthetic scaling test multiplies all these quantities together; scaling tensions alone changes the physical model and need not preserve the shape.

Dimensionless tension is \(T_s/(\mu H)\). The explored values \(T_b/(\mu H)\in\{0.2,2,20\}\) are sensitivity scenarios. They are not experimentally established bounds, an identified set, or a measured confidence interval. An independently measured tension/modulus relation would be needed for that interpretation.

## API and retained evidence

```python
from axisymmetric_prestress import Parameters, build_model
p = Parameters(nr=8, nz=4,
               basal_surface_tension=2.0,
               apical_surface_tension=0.2)
model = build_model(p)
first = model.relax(0.2)
second = model.relax(1.0, first['free_displacement'])
```

The defaults of both new tensions are zero, reproducing the original model's energy and forces. Current full boundary graphs, elementwise orientation checks, strain diagnostics and signed displacement outputs are inherited. Additional energy fields identify `basal_tension`, `apical_tension` and the reference tension energy separately.

Run the synthetic verification with:

```powershell
python -B controlled_furrowing/scripts/axisymmetric_prestress.py
```

`evidence/prestress_protocol.json` retains the checks and thresholds before numerical evaluation. `evidence/prestress_verification.json` retains the source hashes, zero-dose first variation, exact areas, derivative, scaling, warm/cold, work and stationary-scenario checks. Complete outputs remain in `runs/prestress_checks/`. The first attempt completed numerical calculations but failed receipt serialization because a NumPy boolean was passed to the JSON encoder; `evidence/prestress_initial_failure.json` records this failure. Explicit boolean conversion fixed serialization without changing the numerical solver or thresholds.

The separate replay `runs/prestress_checks/check_second_variation.py` compares the automatic-differentiation normal Hessian with the independently integrated edge formula above. Both give 8.126858310618834 for the specified synthetic direction. It also verifies bit-for-bit equality of the original and zero-tension extension's objective and gradient on a nontrivial displaced state. Its receipt is `evidence/prestress_second_variation.json`.

At activation gain 1 in the synthetic small patch with basal area modulus 2 and apical tension 0.2, the explored basal tensions yield:

| Basal tension / \(\mu H\) | Central apical displacement / \(H\) | Central basal displacement / \(H\) |
|---|---:|---:|
| 0.2 | −0.01258 | +0.02681 |
| 2 | −0.01843 | +0.007965 |
| 20 | −0.02095 | +0.001002 |

The sign convention is positive apically/outward. Stronger basal tension suppresses basal vertical motion, allowing apical motion to dominate in this idealized calculation. This demonstrates the mechanical relevance of a missing constitutive ingredient; it does not establish the embryo's tension, calibrate a directional ratio, validate a larger-patch prediction or prove a novel biological mechanism. All predeclared scenarios and their failures must remain visible in the subsequent calibration.

The original limitations remain: P1 penalty locking and mesh sensitivity, idealized reference/rim geometry, local equilibrium rather than a global minimum theorem, local positive determinant rather than global contact control, absent kinetic/peak observation map and no finite ATP/optical budget. The new term does not remove any of these gaps.

![Synthetic effect of basal surface tension on finite-height deformation](../runs/prestress_checks/prestress_summary.png)

The figure uses only the synthetic equilibria above. It displays both surface responses, the change in apical/basal motion sharing, and nonlinear dose response; no displacement exaggeration is used. Figure/input hashes are retained in `evidence/prestress_visual.json`.
