# Why the F and K root bounds remain unresolved

**The retained center-residual bound dominates the predictor residual; widening the error cube then makes the componentwise derivative bound too large.** The parameter-only derivative bounds are below 0.083. This diagnoses the failed sufficient tests, not root nonexistence or biological incompatibility.

The [protocol](../runs/cohort_bound_diagnostics_v1/protocol.json) froze both cached fields, predictors, shared cell, source dependencies and historical radius receipts before calculation. The [read-only diagnostic script](../scripts/cohort_bound_diagnostics.py) freshly assembled the cached point and the predictor-only \(r=0\) remainder/variation. Its [terminal results](../runs/cohort_bound_diagnostics_v1/results.json) reproduce both original 160-bit residual upper bounds exactly. No new root, displacement center, nonlinear optimization, parameter fit or cell shrink was used.

## Residual decomposition

For \(d=\theta-\theta_0\), the reviewed affine predictor gives

\[
 RF(\widehat u(\theta),\theta)=
 RF_0+R(H_0S+B_0)d+R\mathcal R(d).
\]

The kernel bounds the three terms with exact Arb interval arithmetic. The full \(\eta\) is the maximum of the **rowwise sums** of the combined linear enclosure and remainder bound. It is not necessarily the sum of their separate maximum norms.

| Validated sufficient bound | F | K |
|---|---:|---:|
| \(\|RF_0\|_\infty\) | \(8.58449\times10^{-11}\) | \(5.00617\times10^{-11}\) |
| Affine linear-defect bound over the cell | \(2.890\times10^{-21}\) | \(5.336\times10^{-21}\) |
| Predictor Taylor-remainder bound | \(7.424\times10^{-13}\) | \(2.499\times10^{-12}\) |
| Full predictor residual \(\eta\) | \(8.65398\times10^{-11}\) | \(5.22616\times10^{-11}\) |
| Center preconditioned-Hessian defect \(Z_0\) | \(7.269\times10^{-7}\) | \(1.535\times10^{-6}\) |
| Full predictor-only \(Z\), error radius zero | 0.0376960 | 0.0827620 |

Displayed entries are safe rounded upper summaries; exact rational endpoints and the worst-row decomposition are retained. The explicit tangent defect is negligible at this cell width. Most of the residual bound remains even at the original parameter center. This is an algebraic-numerical statement: the preconditioned residual itself is not a certified distance to a root.

## Why a larger error radius failed

The initial tested error cubes had radii approximately \(1.73080\times10^{-10}\) and \(1.04523\times10^{-10}\). Their full \(Z\) bounds were 0.69039 and 0.56533: contraction alone held, but the required self-map inequality

\[
 \eta+Zr<r
\]

did not. The additional 2×, 4× and 8× radii produced \(Z>1\) for both records. The detailed [cohort account](COHORT-UNIFORM-ROOTS.md) and [radius supplement](../runs/cohort_radius_supplement_v1/results.json) retain every attempted radius.

The fresh \(r=0\) bounds show that the shared parameter hull alone does not produce those large \(Z\) values. Adding the error cube drives the conservative bound upward. Subtracting two sufficient upper bounds would not, by itself, certify a separate operator contribution; no such difference is used as a theorem here.

## Floating matrix diagnostics

The archived coefficient-Euclidean Hessians admit floating Cholesky factors. Their approximate eigenvalue extrema are \((0.001175,43828.5)\) for F and \((0.001167,43746.4)\) for K, with spectral condition diagnostics around \(3.7\times10^7\). The products of archived-H and cached-R infinity norms are around \(9.6\times10^7\) and \(9.3\times10^7\).

These are floating diagnostics at the archived rounded solver strengths, not interval results for the exact gain-times-dose target or unknown roots. Eigenvalue magnitudes and conditioning depend on the coefficient basis and norm; they are not mass-normalized physical modes. They do not prove positive definiteness, energetic stability or a continuum stability property.

## Consequence for the next method

The evidence supports examining the center-residual treatment and the norm used for the error enclosure, while keeping the same physical parameters and shared cell. It gives no reason to silently shrink that cell or interpret these failures as rejection of the mechanical model. A changed center or a validated energy-norm tube would be a separately declared computation with its own residual, domain and self-map proof. Neither has been implemented for F or K at this checkpoint.

At this diagnostic checkpoint the cohort had **15 certified records and two unresolved records**, with no all-seventeen stationary-root tuple certificate. The original same-center observation packets consequently cover only 60 channels. The subsequent [separately declared single proof-center refinement](COHORT-UNIFORM-ROOTS.md) produced successful F/K root receipts without changing the model or cell; its independent audit and an explicit center-offset observer adapter remain separate gates. The original failures and this diagnostic calculation are retained unchanged.
