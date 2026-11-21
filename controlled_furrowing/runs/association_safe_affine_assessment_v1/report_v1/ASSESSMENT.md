# How source linkage changes the affine transfer range

The modeled large-region contrast has a **40.19% wider affine range** under the primary disk-weighted observer when calibration keeps only the 34 dose-linked channels. Under the radial-line alternative with component maxima, the range widens **49.61%**. The comparison uses independently checked rational linear-programming bounds on archived tangents.

![Four paired width markers. The two disk observers share the 34-row width 0.08471 micrometres, and the two radial-line observers share 0.09363 micrometres. Every 68-row width is smaller.](affine_widths.png)

| Static observer | 68 ordered rows (µm) | 34 dose-linked rows (µm) | Increase |
|:--|--:|--:|--:|
| Disk, component maxima | 0.06042515 | 0.08470768 | 40.19% |
| Disk, component means | 0.06340250 | 0.08470768 | 33.60% |
| Radial line, component maxima | 0.06258378 | 0.09363464 | 49.61% |
| Radial line, component means | 0.06533135 | 0.09363464 | 43.32% |

For each spatial weighting, the two 34-row cases are exactly identical. Changing the component aggregation affects only the excluded component channels. Markers show the midpoint of each certified range-width interval. The exact rational brackets, all shorter than \(1.9\times10^{-15}\) µm, are narrower than a plotted pixel. The display rounding applies only to the archived affine optimization. [Exact width and ratio CSV](widths.csv) · [Independent mathematical replay](../independent_readback/actual/receipt.json) · [Acceptance and process evidence](../independent_readback/root_acceptance.json)

## What the workbook joins

![Duplicate scalar entries link 17 doses to 17 maximum and 17 mean summaries. Fig. 2e pairs apical and basal within columns and the dose association assumes the same B through R order. No verified specimen identifier crosses these panels.](source_associations.png)

The 34-row group retains the maximum and total regional mean for each dose. Adding the 17 intact apical–basal pairs yields the 68-row group under common B–R column order. The published scalar joins give numerical associations only. The specimen crosswalk, registered boundary profiles and transported region, weights and aggregation order, per-specimen timing and growth references, uncertainty, and a matching large-region measurement remain open. [Observation contract](../../../docs/JOINT-OBSERVATION-CONTRACT.md) · [Public data readback](../../source_contract_public_refresh_v1/PUBLIC-REFRESH.md)

## What the affine range measures

For archived small-channel derivative \(J\), large-contrast derivative \(c\), and parameter step \(d\) in the same trust box \(B\), the two endpoint optimizations are

\[
\max_{d\in B,\;\lVert Jd\rVert_\infty\le\epsilon}\;\pm c^T d,
\qquad \epsilon=0.001\,\mathrm{\mu m}.
\]

Here \(Q\) is the modeled apical center height minus its value at the fixed current radius of 24.94698765 µm. The allowance constrains change from the model baseline in each small channel. Its interpretation is numerical and conditional. The exact stored matrices all have rank five in the archived coordinates. That rank is specific to the five selected parameters, fixed geometry, selected observation conventions, and stored coefficients. A finite perturbation can move the spatial maximum or leave the local linear regime. No experimental error bound or physical identified set follows.

## Exact force scale and conditional dose scale

Multiply all seven force coefficients by \(\alpha>0\). On the positive physical branch, the model energy \(E\), first variation \(DE\), and second variation \(D^2E\) all scale by \(\alpha\). Stationary shapes and their local minimum character are preserved at fixed geometry, activation support, input and displacement constraints, with no independently fixed load. Branch predictions also require scale-invariant selection. The numerical trial penalty is excluded where inactive. With \(L_0=33.5\,\mathrm{\mu m}\), the force groups available to shape inference are

\[
\frac{\kappa}{\mu},\quad \frac{K_s}{\mu L_0},\quad
\frac{T_s}{\mu L_0},\quad \frac{g}{\mu}
\qquad(s=a,b).
\]

The archived tangent fixes \(\kappa/\mu=5000\). The modeled activation enters as \(gI\). If a common unknown conversion \(\beta\) links recorded dose to modeled dose, there is an additional exact invariance

\[
\begin{aligned}
I &= \beta I_{\mathrm{rec}},\\
(g,\beta) &\mapsto(g/a,a\beta),\qquad a>0.
\end{aligned}
\]

Under that assumption shape can constrain the effective ratio \(g\beta/\mu\). The archived calculation takes the dose scale as given. [Constitutive law and scope](../../../docs/PHYSICAL-PARAMETER-DOMAIN.md)

## Strongest next test

After the fixed-field exact-spatial certificates, bound the 34-channel transfer over one predeclared physical parameter cell. Carry a verified nonlinear remainder and allow the maximizing radius to switch. That would test whether the affine ambiguity persists for the modeled branch. Matching the actual experimental transfer still requires a registered larger-region observation with the same operator and an error budget. Compare simpler models under that same operator before assessing detectability or a scientific contribution. Prior work establishes the general same-law ambiguity construction. [Five-gate plan](../../../PLAN.md) · [Originality follow-up](../../transfer_ambiguity_nonlinear_probe_v1/originality_followup/README.md) · [Full analytical assessment](../ASSESSMENT.md)
