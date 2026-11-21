# Matched refinement of the activated finite-height body

This target-free study compares the existing 12×6 mechanical states with a 16×8 meridional mesh, keeping the constitutive law, activation strength, physical geometry, boundary conditions, integration rule and solver policy fixed. It does not read response data or refit activation gains. The purpose is numerical readiness: measure how much the computed activated tissue changes under refinement before treating its geometry as a biological prediction.

The source was frozen before the single job was launched. `scripts/run_activated_mesh_study.py` has SHA256 `78295655c236d465a4c4ac9d324c9d810d23359855e7d54761496535ffe0d8a2`. It checks six frozen dependencies before solving. The [robust continuation policy](QUADRATIC-SOLVE-POLICY.md) accepts only force-and-geometry-verified states, retains rejected trials, and starts every retry from the last verified vector. Failure remains a legitimate outcome.

## What is matched

The reference thickness is \(H=1\), with the declared dimensional conversion \(H=33.5\,\mu\mathrm m\). The physical activation radius is 9.35512037 µm and its depth is 5–25 µm below the top, corresponding to reference height fractions \(8.5/33.5\) through \(28.5/33.5\). The outer radius remains \(3H\); the reference is flat and its outer rim is clamped.

The shear modulus is 1, basal/apical area moduli are 2/0.2, and basal/apical baseline tensions are 2/0. Volume coefficients are 500 and 5000; each branch attempts activation strengths 1, 5 and 10 in sequence. No strength is changed to make the refined model match an earlier shortening value. All parameter fields besides the radial/depth cell counts are checked for exact equality before comparison.

The governing constitutive law is the frozen full pointwise logarithmic-volume penalty, implemented with P2 displacement. Both meshes use the same six-point degree-4 Dunavant volume rule and eight-point surface Gauss rule. A separate high-quadrature study examines integration error. Keeping this rule fixed prevents silently conflating a mesh change with an integration-policy change; it does not prove quadrature error is negligible.

## Eulerian geometry, not corresponding material node indices

The deformed top and bottom have quadratic meridional segments

\[
 R(t)=R_0+R_1t+R_2t^2,\qquad
 Z(t)=Z_0+Z_1t+Z_2t^2,\qquad 0\le t\le1.
\]

Different meshes have different material nodes and radial displacement. Comparing two nodes with the same index or the same reference radius does not generally compare the two surfaces at the same current physical position. The whole-curve comparison instead fixes a common **Eulerian physical radius** \(\rho\) and evaluates both current graphs there.

Stored controls are interpreted as exact rationals and multiplied by the same declared 33.5-µm scale. The radial derivative is linear, so its endpoint values verify monotonicity of each quadratic radial map. Exact rational bisection encloses the root of \(R(t)=\rho\). The corresponding height range is bounded using the bracket endpoints and any quadratic height critical point. No floating-point root approximation is promoted into an exact curve location.

Both models span the same physical radial domain, 0–100.5 µm. The code still explicitly computes and records the common covered domain; incomplete coverage or a nonmonotone surface raises an error rather than silently reinterpreting the surface as a graph.

## A rigorous discrepancy between two discrete curves

For fine/coarse graphs \(a_f,a_c\), partition the common physical radius interval into cells \(I_j\). Rational inversion supplies enclosures

\[
 a_f(I_j)\subseteq[F_j^-,F_j^+],\qquad
 a_c(I_j)\subseteq[C_j^-,C_j^+].
\]

For every \(\rho\in I_j\), their difference belongs to

\[
 [F_j^- - C_j^+,\;F_j^+ - C_j^-].
\]

The maximum absolute endpoint over all cells is therefore an upper bound \(U\) on the supremum discrepancy between the two **stored discrete curves**. At a sampled radius, an interval enclosing the point difference supplies a lower bound on its absolute value: zero if the interval contains zero, otherwise its distance from zero. The maximum of these pointwise lower bounds is \(L\). Thus

\[
 L\le\sup_\rho|a_f(\rho)-a_c(\rho)|\le U.
\]

The same construction is applied to the basal surface. Refining a radial cell preserves complete domain coverage, while every new point remains a valid lower-bound witness. The code refines the cell with the largest remaining upper bound until \(U-L\le0.01\,\mu\mathrm m\) or 4096 active cells are reached. Any unresolved enclosure width is reported. The 0.01-µm value is a requested reporting resolution, not a biological acceptance tolerance.

This enclosure separates two issues: uncertainty in the geometric comparison of saved discrete fields, and their unknown error relative to the continuum mechanical solution. It rigorously encloses the former. It is **not** a finite-element PDE error bound, a Richardson estimate, physical parameter uncertainty or a confidence interval. A small mesh discrepancy alone does not establish convergence to the correct continuum solution; shared underintegration or locking can affect both meshes.

The plot also uses a common 257-point physical radial grid for readable full-curve overlays. Those display samples do not replace the interval-cell supremum enclosure.

## Other measurements and retained evidence

The report separately compares central signed apical/basal displacement and maximum activated material-column shortening. Central displacements refer to the same symmetry axis on both meshes. The maximum column statistic is computed from each model's material-column quadratic extrema; it is not mislabeled an Eulerian thickness maximum or a maximum over biological time.

Each full field retains solver residuals, actual quadrature J values, conservative Bernstein J bounds, volume change, quadratic surface controls, pressure diagnostics and every accepted/rejected attempt. Conservative J bounds are not attained compression measurements. The run protocol records physical scaling, integration policy, solver policy, source hashes and the absence of target responses. Its live process record prevents accidental duplicate jobs.

`evidence/activated_mesh_graph_tests.json` verifies the comparison with synthetic exact tests. In particular, two different material parameterizations of the same Eulerian straight graph produce zero enclosed physical discrepancy, despite differing midpoint material values. A constant 0.1-unit height offset has an exact supremum of 0.1. The full common physical domain is checked explicitly.

The separate retained replay `runs/activated_mesh/verify_graph_bounds.py` compares twelve general quadratic graph differences with their independently computed analytic endpoint/critical-point maxima. Every true maximum lies in the rational interval-cell enclosure, and all twelve requested reporting widths are achieved. Its receipt is `evidence/activated_mesh_graph_adversarial.json`. These are implementation tests for the geometric comparison, not biological evidence.

The principal receipts are `activated_mesh_protocol.json`, `activated_mesh_results.json` and `activated_mesh_comparison.json` under `evidence/`. Complete new states and attempt ledgers are under `runs/activated_mesh/`; existing 12×6 states remain unchanged under `runs/quadratic_robust/`.

## Completed result: scalar agreement does not establish shape convergence

All six 16×8 targets were reached. The retained path contains 22 attempted increments: fourteen accepted and eight rejected. Every accepted field passes the unchanged force/positive-geometry gate. A separate replay verifies all 22 force calculations, all six cross-target input-hash chains and that all six frozen dependencies remain unchanged. All twelve apical/basal supremum enclosures achieve the requested 0.01-µm reporting width. The receipt status is deliberately `PASS_NUMERICAL_EXECUTION_NOT_MESH_ADEQUACY`.

The following intervals enclose the supremum discrepancy between the two stored discrete surfaces over the full common physical radial domain. They do not enclose either mesh's continuum error.

| Volume coefficient | Gain | Apical curve discrepancy (µm) | Basal curve discrepancy (µm) | Fine minus coarse maximum column shortening (µm) |
|---|---:|---|---|---:|
| 500 | 1 | [0.12359, 0.13307] | [0.00692, 0.01582] | +0.08830 |
| 500 | 5 | [0.45461, 0.46459] | [0.02802, 0.03786] | +0.04381 |
| 500 | 10 | [0.70222, 0.71216] | [0.04219, 0.05201] | −0.15177 |
| 5000 | 1 | [0.18652, 0.19622] | [0.00789, 0.01629] | +0.18514 |
| 5000 | 5 | [0.51174, 0.52068] | [0.04860, 0.05815] | +0.23601 |
| 5000 | 10 | [0.74786, 0.75551] | [0.09616, 0.10594] | +0.15060 |

Table intervals are rounded outward from the exact rational receipts. For example at coefficient 500 and gain 5, the scalar shortening changes by only about 0.044 µm, while the apical profile differs somewhere by at least 0.4546 µm. Good agreement in a central or maximum shortening statistic would therefore miss a substantially larger shape difference.

Even a qualitative parameter ordering is not stable in these discrete calculations. At gains 5 and 10, the 12×6 model gives greater shortening for coefficient 500 than 5000, whereas the 16×8 model reverses this ordering. This is a warning about numerical readiness, not a robust physical effect or a new constitutive discovery.

There is a further limitation beyond ordinary discretization error: the adaptive controller accepted different intermediate gains on different meshes. The energy is nonlinear and no uniqueness/global-minimum theorem has been proved. Consequently a local-equilibrium branch-selection effect has not been excluded. The differences are accurately described as matched numerical mesh-study discrepancies; they must not all be attributed to a rigorously isolated approximation-error mechanism.

![Whole Eulerian curve differences, unchanged-gain shortening and discrete supremum enclosures](../runs/activated_mesh/activated_mesh_comparison.png)

The top panels compare surfaces at the same current physical radius over the entire 0–100.5-µm domain. The lower-right bars are exact-arithmetic enclosures between discrete curves, not confidence intervals or PDE error bars. Figure/input hashes are retained in `evidence/activated_mesh_visual.json`.

## Consequence for the next accuracy gate

This calculation does not pass mesh adequacy. Further refinement or a verified adaptive mesh, branch checks, and integration-error control are needed before assigning a justified continuum uncertainty to the activated shape. No 24×12 job was launched as part of this first 16×8 delivery.

The separate [activated quadrature study](ACTIVATED-QUADRATURE.md) changes integration in the same 12×6 space and reequilibrates those fields. Its same-material/Hausdorff curve bound is not the same metric as the Eulerian vertical supremum used here; a numerical ratio between the two would require an appropriate conversion. Neither study alone supplies a rigorous continuum error bound, and mesh/integration interactions remain to be checked.

This is established numerical verification practice applied transparently to this model, not a new biological concept. Original study gates still require numerical adequacy, defensible material/geometry uncertainty, a justified observation/time operator and genuine larger-region testing.
