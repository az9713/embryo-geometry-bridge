# One actual multivariate calibration step from all 17 records

**Five material/input coordinates now enter a joint local calibration using all 68 small-region endpoint observations.** The calculation has run on the full cohort, with all four observation conventions retained at the same 17 baseline fields. Four local affine minimax proposals are exactly certified. Their apparent improvements remain predictions of the affine approximation until fresh nonlinear fields verify them. This is neither a full physical identified set nor a completed larger-region prediction.

The [runner](../scripts/material_cohort_descent.py), [frozen protocol](../runs/material_cohort_descent_v1/protocol.json), [completed cohort result](../runs/material_cohort_descent_v1/results.json) and [authoritative exact candidates](../runs/material_cohort_descent_v1/exact_extraction_v2/exact_candidates.json) retain the inputs, matrices, graph observations and arithmetic receipts. No new mechanical equilibrium was solved in this tangent step.

## 1. The shared baseline and the inferential question

The preceding [joint gain study](JOINT-GAIN-CALIBRATION.md) directly checked four gain candidates. The disk-weighted, spatial-maximum-component convention had the smallest existing calibration discrepancy, approximately 5.37984 µm. Its gain \(g_0=8.0629720200\ldots\) and 17 corresponding fields define this step's common baseline. Selection uses small-region calibration alone, without larger-region outcomes.

All four conventions are now evaluated at **these same 17 fields and this same gain**. Consequently the other three baseline discrepancies differ from those in the preceding study, which gave each convention its own gain. There is one cohort of 68 observations viewed four ways, not four independent cohorts.

Geometry and the normalized shear coefficient stay fixed. The compression coefficient is fixed at 5000 as the declared numerical penalty; it is not inferred here as a physiological modulus. The five local coordinates are

\[
 z=(\Delta\log K_a,\Delta\log K_b,\Delta T_a,
       \Delta\log T_b,\Delta\log g).
\]

The initial \(T_a=0\) makes an additive nonnegative coordinate appropriate. The computational trust region is

\[
 |z_1|,|z_2|,|z_4|,|z_5|\le0.1,\qquad 0\le z_3\le0.02.
\]

These are bounds on one numerical step, not measured material priors, confidence limits or an exhaustive uncertainty domain. Surface coefficients and geometry use the fixed physical length \(L_0=33.5\) µm described in [the parameter-domain analysis](PHYSICAL-PARAMETER-DOMAIN.md).

## 2. Mechanical sensitivity: what is proved and what is only checked

For equilibrium equations \(R(u,\eta)=D_uE(u,\eta)=0\), differentiation gives

\[
 H S_j+B_j=0,\qquad
 H=D_uR,\quad B_j=D_{\eta_j}R,\quad
 S_j=\frac{\partial u}{\partial\eta_j}.
\]

This derivative identity follows from the implicit-function theorem **if** an exact regular root and a smooth chosen branch exist. The present matrices are floating evaluations at saved force-tolerance states; no validated algebraic-root enclosure or certified invertibility constant is supplied. Their linear solutions are therefore numerical tangent diagnostics, not a completed application of a validated implicit-function theorem.

The [mechanical tangent implementation](MATERIAL-TANGENT-INFERENCE.md) computes six columns, including \(\log\kappa\). All six are retained and checked for every record; only the other five enter the step. Thus there are 102 mechanical linear systems, but \(\kappa\) does not become a fitted coordinate. Each record retains its sparse Hessian, force derivatives, displacement tangents, finite-difference force checks, source hashes and linear-system residuals.

## 3. Observing perturbed geometry without inventing equilibria

For a tangent column \(S_j\), the code evaluates boundary graphs from coefficient vectors

\[
 u^\pm=u\pm h S_j,\qquad
 D^{(h)}_{ikj}=
 \frac{\mathcal O_k(u^+)-\mathcal O_k(u^-)}{2h}.
\]

The two predeclared steps are \(h=10^{-3}\) and \(5\times10^{-4}\). Stale stored boundary controls are removed; fresh quadratic traces are extracted from the perturbed coordinates and passed directly to the graph observer. These perturbed coefficient vectors are **not** marked or promoted as newly solved equilibria. A local orientation screen and graph-admissibility checks precede endpoint evaluation. The protocol permits retaining failed graph trials and shrinking both steps; no shrink was needed in this cohort.

The observer retains the current-radius vertical height, line/disk means and spatial-max/regional-mean component alternatives. It preserves signs and uses the same conditional flat component references and registered current ROI. The [source observation ambiguities](JOINT-OBSERVATION-CONTRACT.md), temporal peak closure and canonical source-column pairing remain unresolved.

For a smooth regional mean, the separate exact trace-gradient formula also supplies

\[
 D\mathcal O_k(u)[S_j]
 =\nabla_{\text{physical nodes}}\mathcal O_k
       \cdot(L_0 S_j).
\]

Exact rational trace-gradient intervals are contracted with the serialized floating tangent coefficients. This checks the graph transport terms independently of the central-difference formula. It does not validate the mechanical tangent itself.

There are 340 graph perturbations: 17 records × five coordinates × two steps × two signs. The largest two-step midpoint derivative difference is approximately \(1.27\times10^{-8}\) µm per coordinate; the largest fine-step mean-difference distance from its exact trace-gradient contraction is approximately \(4.23\times10^{-9}\). The largest central-difference interval width is approximately \(1.32\times10^{-19}\) µm per coordinate. These are retained numerical diagnostics, not bounds on truncation error or proof that a maximum is differentiable.

A maximum can be nonsmooth when its maximizing location changes. For example, \(\max(z,-z)=|z|\) has zero symmetric difference at zero but a nonzero right derivative. The synthetic tests explicitly retain this warning. Agreement between two chosen steps cannot eliminate it for arbitrary parameters.

## 4. The affine minimax problem and its exact certificate

For each convention separately, let \(r_{ik}\) be the baseline interval-midpoint prediction minus the recorded endpoint, and let \(J_{ikj}\) be the smaller-step central-difference midpoint. The declared affine problem is

\[
 \min_{z\in B,t\ge0}t,\qquad
 -t\le r_{ik}+\sum_jJ_{ikj}z_j\le t
 \quad(i=1,\ldots,17;\ k=1,\ldots,4).
\]

All channels have the same micrometre minimax objective. This is a discrepancy criterion, not a fitted statistical noise model. Every record enters; none is discarded for having a large residual.

The serialized binary64 derivative midpoints define exact rational coefficients for this finite affine problem. The root [bounded-affine helper](BOUNDED-AFFINE-MINIMAX.md) accepts those Fractions without re-rounding. A floating LP proposes an active basis; exact rational linear algebra then reconstructs a feasible primal vector and nonnegative dual multipliers.

Writing the full endpoint, trust-box and \(t\ge0\) inequalities as \(Ax\le b\), with \(x=(z,t)\), the certificate verifies

\[
 Ax\le b,\qquad\lambda\ge0,\qquad
 A^\mathsf T\lambda=-e_t,\qquad
 t=-b^\mathsf T\lambda.
\]

Indeed, for any feasible \(x'\),

\[
 t'=e_t^\mathsf Tx'
   =-\lambda^\mathsf TAx'
   \ge-\lambda^\mathsf Tb=t.
\]

This proves global optimality of the **defined bounded affine LP**. It does not prove the derivative approximation, its trust-region error, nonlinear improvement or physical parameter compatibility. The producer's independent weak-dual interval also contains each exact optimum.

The initial extraction is preserved unchanged. After helper input-validation hardening, a separate `exact_extraction_v2` protocol reran only the exact arithmetic on the **same** rational coefficient tables; no tangent or equilibrium calculation was repeated. The old helper source and provenance are archived under `runs/bounded_affine_minimax_review/initial_revision/`.

## 5. Completed local proposals

| Observation convention | Common baseline discrepancy, µm | Exact affine optimum, displayed µm | Predicted reduction, µm |
|---|---:|---:|---:|
| Disk / spatial-max components | 5.37983567 | 5.18386299 | 0.19597269 |
| Disk / regional-mean components | 6.46771627 | 6.07669043 | 0.39102585 |
| Line / spatial-max components | 6.18326952 | 5.55397817 | 0.62929135 |
| Line / regional-mean components | 6.18326952 | 6.02644077 | 0.15682875 |

The proposed steps are:

| Convention | Δlog Ka | Δlog Kb | ΔTa | Δlog Tb | Δlog gain |
|---|---:|---:|---:|---:|---:|
| Disk / spatial max | 0.1 | −0.1 | 0 | 0.1 | −0.03153743 |
| Disk / regional mean | −0.1 | −0.1 | 0 | 0.1 | 0.1 |
| Line / spatial max | 0.1 | 0.1 | 0.00791817 | 0.1 | −0.1 |
| Line / regional mean | 0.1 | −0.1 | 0 | 0.1 | −0.02180194 |

Several coordinates reach computational trust limits. This signals the limitation of this local step; it does not identify biological moduli or justify expanding the physical domain by fiat. Even the affine losses remain several micrometres.

![Common baseline versus four affine proposals, with trust-step diagnostics](../runs/material_cohort_descent_v1/material-cohort-affine-descent.png)

The figure has no invented error bars. The bar comparison is between directly verified baseline values and an exactly solved local approximation; the second bars are not new nonlinear results. Its hashes and generator are retained in [the visual receipt](../runs/material_cohort_descent_v1/visual_verification.json).

## 6. Retained receipts and the next physical check

The tangent driver SHA256 is

    8343f4908a76a47f8c83159e7befbbaa23afdc7d4a5f4689bda9ec901e3957c1

The frozen cohort protocol SHA256 is

    dbec2723481dc9a3a66c6dbde0648721b7283b734b33d532477e2d3357d47c47

The complete cohort result SHA256 is

    70c3b02bd4a26e57181875715af5ad3865688a76374ca26de67bc7ac10d89259

The authoritative v2 exact candidate file SHA256 is

    a87d92a1b9216e5b94ae2ee68944e10aac46711e5c854c78d8da7e19842d7758

Ten producer synthetic checks pass. The [independent mechanical cohort review](MATERIAL-COHORT-DESCENT-REVIEW.md) passes 598 numerical checks. The separate [exact-affine independent receipt](../evidence/bounded_affine_minimax_independent_review.json) passes all 16 aggregate checks, including full cohort coefficient/observation reconstruction and all four exact certificates, plus 206 independent analytic/vertex/permutation checks and 15 input guards. These reviews remain conditional numerical/affine assessments, not physical derivative or nonlinear-fit certificates.

The next physical test is a separately versioned [direct runner](../scripts/material_cohort_direct.py) and [predeclared protocol](../runs/material_cohort_direct_v1/protocol.json) for the **primary disk/spatial-max proposal only**. It transforms the material to approximately \(K_a=0.22103418,K_b=1.80967484,T_a=0,T_b=2.21034184,g=7.81265450\), leaving geometry and compression coefficient fixed. It must create a new-material verified zero state and solve all 17 new inputs, with Q8/Q12 force, exact geometry and all four observation conventions. The other three variant-specific affine proposals are not silently treated as directly verified.

There is a branch-history caveat: the old 17 baseline fields were selected from a union of 68 ascending target strengths, whereas the new replay visits 17 targets. Without uniqueness and stability proofs, a changed result can include continuation/branch-history effects. This test checks the proposal under its declared new path; it does not isolate a causal material-change contribution. No matched additional baseline-control job is silently added.

No full material/geometry \(\Theta_\delta\), uniform nonlinear/continuum response bound, externally measured discrepancy law or source-matched larger-region prediction follows from this local calculation. None of the original incomplete research gates is marked complete.
