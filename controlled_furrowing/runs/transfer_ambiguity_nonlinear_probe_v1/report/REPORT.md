# Nonlinear transfer ambiguity in one epithelial model

Two frozen parameter choices satisfy all 68 small-region entries per choice under the declared **disk-weighted, spatial-maximum** observer. At the same large optical input, their saved apical curves differ in signed center-to-shoulder contrast by **0.154755156665 µm**. The result supports a narrow claim about two numerical fields under one explicit observation and residual contract.

**Observer choice changes the conclusion.** The common allowance is **5.38983567480 µm**, and the smaller remaining margin is **0.00481253906 µm**. Each of the other three predeclared static observer conventions needs a larger allowance for the same pair. The allowance is a conditional model discrepancy threshold derived from an affine residual. Experimental noise, timing, reference and registration remain separately unresolved. [Primary result](../actual/results.json) · [Independent primary replay](../independent_actual_preparation/actual_replay/receipt.json) · [Independent observer sweep](../observer_robustness/independent_actual/receipt.json)

![Two zoomed panels compare the signed nonlinear and affine gaps, and put both endpoint residual requirements below the frozen allowance.](contrast_and_allowance.png)

## Intervention and parameter selection

Both candidates use the same axisymmetric P2 finite-element law and geometric conventions. They have distinct material coefficients and one common gain *within each candidate* across its 17 small doses and its large intervention. The matched large optical dose is **0.89767374**, the diameter is **49.8939753 µm**, and the current-radius shoulder is **24.94698765 µm**. Optical dose is held common between candidates, while their mechanical activation strength is the exact candidate gain times that dose, rounded once for the numerical solver.

The frozen protocol selected two ends of the local affine contrast frontier at \(\delta=D_{\mathrm{aff},0}+0.01\,\mu\mathrm m\), where \(D_{\mathrm{aff},0}\approx 5.37983567480\,\mu\mathrm m\) is the archived baseline maximum residual. Each exact selected perturbation was halved before any new nonlinear calculation. Labels \(+\) and \(-\) identify the affine contrast directions. [Affine derivation](../../transfer_ambiguity_linear_v1/REPORT.html) · [Frozen nonlinear protocol](../preparation/protocol.json)

| Endpoint | \(\Delta\log K_a\) | \(\Delta\log K_b\) | \(\Delta T_a\) | \(\Delta\log T_b\) | \(\Delta\log g\) |
|---|---:|---:|---:|---:|---:|
| + | -0.05000000 | +0.05000000 | +0.00000000 | -0.05000000 | -0.02523642 |
| − | +0.05000000 | -0.05000000 | +0.00000000 | +0.05000000 | +0.00597390 |

The small cohort was computed in ascending dose order from a fresh zero equilibrium for each candidate. The large cases each began from a separate new zero and used one-third, two-thirds and full dose checkpoints. This procedure produced 34 small final fields, two large final fields and eight supplementary loading checkpoints. No large experimental response entered selection or scoring. [Whole field manifest](../field_manifest.json)

## Exact comparisons of the saved curves

For record \(i\) and static observable \(k\), let \(y_{ik}\) be the recorded small-region value and \(I_{ik,j}^{s}\) the exact interval for the saved curve of candidate \(j\in\{+,-\}\). Compute

\[
D^{s}_{j}=\max_{i,k}\ \sup_{x\in I_{ik,j}^{s}} |x-y_{ik}|,\qquad
D^{s}_{\mathrm{pair}}=\max(D^{s}_{+},D^{s}_{-}).
\]

The supremum of absolute residual on each closed interval is its larger endpoint residual. The decision uses all 68 rows per candidate and exact rational endpoints. The table displays rounded midpoints of the large contrast intervals and rounded upper bounds on the required allowances. The margin column is \(\delta-D^{s}_{j}\).

| Endpoint | Stored large \(Q\), µm | Required allowance upper bound, µm | Remaining margin, µm |
|---|---:|---:|---:|
| + | -4.13600326164 | 5.28396210625 | 0.10587356855 |
| − | -4.29075841831 | 5.38502313574 | 0.00481253906 |

For each large saved curve, write \(z_a(R)\) for apical height at current physical radius \(R\). The signed statistic is

\[
Q=z_a(0)-z_a(a),\qquad a=24.94698765\,\mu\mathrm m .
\]

A negative \(Q\) indicates a center below the shoulder. The actual positive signed gap persists at **0.154755156665 µm**, compared with the frozen affine half-step figure **0.154900794334 µm**. The actual minus affine gap is **-0.000145637669 µm**. This is agreement at the two selected endpoints. The experiment supplies no bound on the affine approximation over the rest of the parameter domain. [Exact gap and difference intervals](../independent_actual_preparation/actual_replay/receipt.json)

![Stored apical profiles and their pointwise height difference, both on current radius with the common shoulder marked.](large_current_radius_profiles.png)

The graph uses midpoints of the retained exact curve evaluations. If \(\Delta z(R)=z_{a,+}(R)-z_{a,-}(R)\), then \(Q_+-Q_-=\Delta z(0)-\Delta z(a)\). The lower panel illustrates that the gap is a difference between two registered heights, each of which also changes across the pair.

## Dependence on the static observer

All four conventions use the same saved small fields, recorded values, flat reference, current-radius region and total maximum shortening. For the regional mean, with \(0\le R\le R_0\), the normalized disk and radial-line weights are \(2R/R_0^2\) and \(1/R_0\), respectively. The two separate apical and basal component values are either spatial maxima or regional means. Separate maxima can occur at different radii. The maximum of the total shortening therefore has its own calculation.

| Static observer | Required common allowance, µm | Excess over frozen \(\delta\), µm | Pair at frozen \(\delta\) |
|---|---:|---:|---|
| Disk · component maxima | 5.3850231357 | -0.0048125391 | Both satisfy |
| Disk · component means | 6.5760903956 | +1.1862547208 | At least one exceeds |
| Radial line · component maxima | 6.1919994826 | +0.8021638078 | At least one exceeds |
| Radial line · component means | 6.1919994826 | +0.8021638078 | At least one exceeds |

![Signed excess above the frozen allowance for four observers. Only disk weighting with component maxima has negative excess.](observer_dependence.png)

An independent sweep reproduced all 544 scalar comparisons from the same 34 fields. The all-four common requirement is **6.5760903956 µm** in this sweep. Passing the primary observer provides limited evidence when the source-to-model observation rule remains uncertain. The alternatives here affect static weighting and aggregation. The unresolved source contract also includes timing of reported extrema, the reference for growth and the mapping of the region of interest. [Exact secondary results](../observer_robustness/actual/results.json) · [Independent secondary audit](../observer_robustness/independent_actual/receipt.json) · [Observation contract](../../../docs/JOINT-OBSERVATION-CONTRACT.html)

## Numerical evidence and the missing error accounts

The independent replay verified all 44 checkpoints, recomputed 88 numerical force evaluations and checked 6,336 discrete element geometries. The largest replayed scaled force infinity norms were **6.469e-09** for Q8 and **6.469e-09** for Q12, under the \(10^{-7}\) gate. The smallest exact stored lower Jacobian bound was **0.957775**. The stored fields passed exact positive orientation and embedding checks. Q8 and Q12 refer to the main and cross-check integration policies, respectively. They assess residual at the saved coordinates in two quadrature rules.

The stored signed gap interval has width about **5.48e-24 µm**. This resolves evaluation of the serialized quadratic curves. To bound a stationary point of the exactly integrated finite-element problem, a separate proof must enclose both the signed residual and a valid inverse or contraction on a neighborhood. An observer transfer theorem can then turn a certified nodal displacement bound, including binary addition gaps, into an error on the plotted contrast and all 136 small comparisons. The linked design allocates targets from the stored slacks and includes strict radial-derivative margins. Those targets constrain future certificates; actual root-position, continuum and source errors remain open for this pair. [Pointwise certification design](../certification_design/README.md) · [Nodal observer theorem](../../../docs/NODAL-OBSERVER-ERROR-TRANSFER.html)

Practical discrimination requires a registered measurement of the same \(Q\) plus a budget for the experimental, geometric, temporal and model errors. For two certified response intervals separated by \(\gamma>0\), any single scalar estimate has conditional worst-case absolute error at least \(\gamma/2\). This follows from the triangle inequality. Additional deterministic error radii \(\eta_+,\eta_-\) preserve separation if \(\eta_++\eta_-<\gamma\). Here the stored gap establishes neither a bound on those radii nor a necessary instrument specification.

## Originality and scope

The general same-law discrimination construction is established. Liu, Maini and Baker develop interventions that distinguish candidate parameter sets and give a logistic-model example. Kumar et al. find multiple epithelial parameter clusters that reproduce similar control morphology. Countryman et al. already establish activation-size-dependent furrowing and examine basal mechanics. Deterministic, bounded-error estimation of a quantity of interest has an existing optimal-recovery theory. [Liu et al. (2026)](https://doi.org/10.1016/j.mbs.2026.109710) · [Kumar et al. (2024)](https://www.nature.com/articles/s41540-024-00375-z.pdf) · [Countryman et al. (2025)](https://www.nature.com/articles/s41467-025-62483-6) · [Foucart and Liao](https://arxiv.org/abs/2111.02601) · [DeVore et al.](https://foucart.github.io/publi/DFPW_recovery_final.pdf)

The narrow contribution is a reproducible, predeclared nonlinear pair for this discrete epithelial model with documented dependence on the observation convention. A bounded literature search cannot establish priority. A source-compatible mechanical witness, exact equilibrium error, a continuous loading branch and practically attainable measurement discrimination remain open. [Originality follow-up and search locators](../originality_followup/README.md) · [Originality audit receipt](../originality_followup/receipt.json) · [Build provenance](report_build_receipt.json)
