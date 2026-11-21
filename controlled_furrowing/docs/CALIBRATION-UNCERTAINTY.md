# What uncertainty does the small-region calibration require?

**A single nondecreasing curve cannot fit all 17 small-region maximum-shortening records with a uniform residual below 3.131943715 µm.** Adding nonnegative concavity and a zero-response origin raises the sharp floor to approximately 3.648762214 µm. A homogeneous linear response needs approximately 4.587984160 µm. These are exact deterministic compatibility statements for declared curve classes and the published numerical records. They are **not** estimates of segmentation noise, biological variation, or 95% confidence intervals.

This work uses only the 17 small-region OptoGEF2 observations. No larger-region outcomes or spatial-solver files enter the calculation. The [script](../scripts/calibration_uncertainty.py), [exact bounds and witnesses](../evidence/calibration_uncertainty_bounds.json), and [verification receipt](../evidence/calibration_uncertainty_verification.json) retain the calculations. The [original visualization](../runs/calibration_uncertainty/compatibility_envelopes.png) displays all four observables and their admissible envelopes; the shaded regions are not statistical confidence bands.

## Observation unit and measurement definition

The primary source is Countryman et al., Nature Communications 16, 7665 (2025), [DOI](https://doi.org/10.1038/s41467-025-62483-6), [article](https://pmc.ncbi.nlm.nih.gov/articles/PMC12361397/). Methods Sec22 explicitly says, “All measurements were taken from distinct samples.” It also identifies 17 OptoGEF2 embryos, including one without visible shortening. The reporting summary, page 2, describes random selection of individual embryos, genotype-based group allocation, automated apical segmentation blinded to condition, and manual basal segmentation blinded to activation magnitude but not genotype. The Fig2 caption does not itself contain an explicit sentence equating every marker with an individual embryo; the experiment/sample descriptions nevertheless support treating the 17 records as distinct biological samples.

Missing specimen IDs limit cross-file identity verification. They are **not evidence of dependent samples**. Conversely, distinct embryos do not establish independent, identically distributed measurement errors. Sec17 normalizes activation fluorescence against a reference curve averaged across experiments; a shared calibration can couple measurement errors even when specimens are distinct.

The response is not simply a cell's final height. Methods Sec16 constructs an embryo-specific counterfactual height curve using nonactivated tissue and its ongoing cellularization. Maximum shortening is the largest reference-relative reduction across **both spatial location and observation time**. Region-averaged shortening is a separate temporally maximized, spatially averaged summary; maximum and averaging must not be assumed to commute without a source-compatible implementation. Apical and basal displacement maxima are separate observables; they need not occur at the same position or time and cannot be added to reconstruct maximum shortening.

Sec14 gives sampling and protocol information, including nominal 2 µm z-plane spacing and 22 s acquisition intervals. Sec15 describes manual correction and manual basal annotation. Neither spatial sampling nor the integer-valued apical entries supplies a validated annotation-error bound. No numerical per-record segmentation, coordinate, dose-error interval, or complete measurement-noise law was found in the inspected methods, source worksheets, supplementary captions, and reporting summary. This is a bounded source audit, not a claim that such errors are intrinsically unknowable.

## A standard sharp monotone compatibility theorem

Let the fixed measured doses be \(I_i\), responses be \(y_i\), and consider the contract

\[
f\text{ nondecreasing},\qquad |f(I_i)-y_i|\le\delta\quad(i=1,\ldots,n).
\]

There is no stochastic error assumption. On the observed dose interval, the smallest compatible uniform residual is

\[
\boxed{\delta_* = \frac12\max_{I_i\le I_j}(y_i-y_j)_+.}
\]

**Necessity.** For every ordered pair, monotonicity and the two residual inequalities imply

\[
y_i-\delta\le f(I_i)\le f(I_j)\le y_j+\delta,
\]

so \(2\delta\ge y_i-y_j\).

**Sufficiency and exact envelopes.** For \(\delta\ge\delta_*\), define

\[
L_\delta(I)=\max_{I_i\le I}(y_i-\delta),\qquad
U_\delta(I)=\min_{I_j\ge I}(y_j+\delta).
\]

Both are nondecreasing and \(L_\delta\le U_\delta\). At each observed dose they lie within every applicable observation interval. Thus each endpoint envelope, and their midpoint, gives an admissible fitted vector. Between observations, monotone interpolation gives an admissible curve. Any admissible curve must satisfy \(L_\delta\le f\le U_\delta\); each endpoint value can be attained by an admissible extension, so these are sharp pointwise bounds. At coincident doses, all records at that dose must share one fitted value; the implementation explicitly handles this case. Outside the observed dose interval, one side of the unconstrained envelope is unbounded unless additional assumptions are supplied.

This is standard unweighted \(L^\infty\) isotonic regression, not a new theorem. See Stout, [*Weighted L∞ isotonic regression*, JCSS 91 (2018), 69–81](https://doi.org/10.1016/j.jcss.2017.09.001), especially the violating-pair lower bound and minimum/maximum regressions. A [primary author-hosted manuscript](https://web.eecs.umich.edu/~qstout/pap/LinfinityIsoReg.pdf) was checked in the literature audit.

## Exact small-data obstruction and its meaning

The strongest maximum-shortening witness is:

| Source-column record | Dose | Maximum shortening |
|---|---:|---:|
| F |0.46402059|7.64808362 µm|
| I |0.47133397|1.38419619 µm|

The slightly larger measured dose has a 6.26388743 µm smaller response. Therefore \(\delta_* = 3.131943715\) µm. This says that **any common monotone deterministic curve has at least that much worst-case residual**. It does not identify how much comes from measurement error, dose error, specimen heterogeneity, missing mechanics, or a wrong common-curve assumption.

Exactness refers to the numerical representation. The frozen CSV uses canonical shortest decimal representations of the source workbook's floating-point values. The original OOXML spells these two responses as 7.6480836200000004 and 1.3841961899999999, giving the exact token-based floor 3.13194371500000025 µm. The 2.5×10⁻¹⁶ µm difference is serialization, not meaningful biological precision. Both calculations and all 85 source-cell comparisons are retained.

The witness also depends on measured dose order. The dose separation is 0.00731338. If each dose has an unknown absolute error at most 0.00365669, the witness order cannot reverse, so its lower bound survives. A larger dose-error allowance can reverse this particular pair; it does not establish that every other obstruction disappears. No such dose-error allowance has been independently measured here.

Applying the same calculation separately gives:

| Declared monotone observable | Sharp CSV residual floor | One strongest witness |
|---|---:|---|
| Maximum shortening |3.131943715 µm|F → I|
| Region-averaged shortening |2.228529505 µm|K → L|
| Separate apical maximum |3 µm|F → I; M → P also attains the floor|
| Separate basal maximum |1.24138052 µm|F → Q|

Monotonicity of a biological component response is a conditional modeling restriction, not established by the table itself. The separate floors cannot be summed to infer a height error, because the component maxima are not additive.

## Why the floor can apply to a quasistatic peak model

A common instantaneous finite-element shortening curve need not be monotone. There is, however, a different elementary monotonicity statement for a correctly defined peak observable. Suppose one common quasistatic branch \(q(A;\theta)\) is available over activation amplitudes \(0\le A\le gI\), with \(g\ge0\), and define

\[
M(I;\theta)=\sup_{0\le A\le gI}\mathcal O(q(A;\theta)).
\]

If \(I_1\le I_2\), the first amplitude interval is contained in the second, so \(M(I_1;\theta)\le M(I_2;\theta)\), regardless of whether the instantaneous observable is monotone. Including a spatial maximum inside \(\mathcal O\), or using a separate apical/basal observable, leaves this argument unchanged. This is the standard maximum-over-nested-sets observation, not new mathematics.

Consequently, the 3.131943715 µm floor is necessary for a **common-branch, scaled-cap, full-path quasistatic peak model** with a uniform residual contract and fixed doses. The source's spatial-temporal maximum motivates considering such an operator. It does not prove its assumptions: finite observation windows, dose-dependent pulse histories or branches, hysteresis, reference-growth differences, specimen geometry and material heterogeneity can destroy the common nested-set representation. The current audit does not prove instantaneous FE monotonicity, path coverage, or equivalence between an equilibrium peak and the experimental temporal maximum.

## A sharper conditional concavity obstruction

For the relaxed class of nonnegative nondecreasing concave curves with \(f(0)=0\), the sharp floor is

\[
\delta_{\rm conc}
=\frac{15244809486711951}{4178077000000000}
\approx3.648762214461809\;\mu\mathrm m.
\]

Let record I have dose \(x_j=0.47133397\), response \(y_j=1.38419619\), and record M have \(x_k=0.78208913\), \(y_k=12\). With \(\beta=x_j/x_k\), origin-anchored concavity forces \(f(x_j)\ge\beta f(x_k)\). Combining the two observation intervals gives

\[
y_j+\delta\ge\beta(y_k-\delta),\qquad
\delta\ge\frac{\beta y_k-y_j}{1+\beta}=\delta_{\rm conc}.
\]

This lower bound is attained for all 17 records by the exact rational curve

\[
f(x)=\min(gx,C),\qquad
g=\frac{y_k-\delta_{\rm conc}}{x_k}\approx10.6781151472,
\qquad C=\max_i(y_i-\delta_{\rm conc})\approx9.08448084554\;\mu\mathrm m.
\]

It is nonnegative, nondecreasing and concave. All 68 rational linear-program constraints—34 observation inequalities, 17 monotonicity inequalities, 16 slope-concavity inequalities and nonnegative residual—are verified exactly. An exact dual witness combines the two observation inequalities and the origin-concavity consequence to produce \(-\delta\le-\delta_{\rm conc}\). A separate floating-point linear program agrees; the rational primal/dual witnesses, rather than the floating-point optimizer, establish sharpness. Repeating the exact check using original OOXML numeric tokens gives the same conclusion up to serialization differences, retained in the receipt.

This relaxed class is relevant to the separately derived column null, conditionally. Write

\[
I(S)=\frac{(1-2\nu)x^{3/2}+\nu x^{5/2}-(1-\nu)}{a},
\qquad x=(1-S/h)^{-1},
\]

with \(a,h>0\), \(0\le\nu\le1\), and \(0\le S<h\). Direct differentiation gives

\[
I'(S)=\frac{x^{5/2}}{2ah}[3(1-2\nu)+5\nu x]>0,
\]

\[
I''(S)=\frac{5x^{7/2}}{4ah^2}[3(1-2\nu)+7\nu x]>0.
\]

For \(x\ge1\), the two bracketed factors are bounded below by \(3-\nu>0\) and \(3+\nu>0\), respectively. Hence the inverse response has \(S'>0\), \(S''=-I''/(I')^3<0\), and \(S(0)=0\). The relaxed concavity floor is therefore a necessary compatibility bound for that column family. The piecewise-linear attaining curve is not asserted to belong to the more restrictive column family. Concavity has not been established for the full spatial FE maximum or for a general peak-over-branches operator.

For comparison, the stricter homogeneous linear null \(f(I)=gI\), \(g\ge0\), has exact minimum uniform error approximately 4.587984159785218 µm. Its admissible gain interval at a declared \(\delta\) is

\[
\left[\max\!\left(0,\max_i\frac{y_i-\delta}{I_i}\right),
\min_i\frac{y_i+\delta}{I_i}\right].
\]

The interval first becomes nonempty at that floor, with gain approximately 9.47720094283 µm per normalized dose. This response gain is a null-comparator coefficient, not a direct material modulus or the active-stress parameter of the coupled solver.

## Concrete uncertainty contracts and the remaining physical inference

A defensible report can present the family of deterministic contracts \(\mathcal C_\delta\), rather than invent one physical error bar. For a common monotone curve, \(\mathcal C_\delta\) is empty below the sharp floor; at and above it, the receipt gives exact pointwise envelopes at \(\delta_*\), \(\delta_*+0.5\) µm, \(\delta_*+1\) µm and \(2\delta_*\). Nonnegativity can additionally clip the lower envelope at zero without changing these data's minimum floor. This is an extra physical restriction, not an observed zero-dose datum.

For an actual physical family, the desired object is instead

\[
\Theta_\delta=\{\theta\in\Theta_{\rm prior}:
|\mathcal O_i(q(\theta,I_i))-y_i|\le\delta\;\forall i,
\ \text{declared geometry and solver-validity conditions hold}\}.
\]

The curve bounds can exclude a noise contract for that family only when the family's required shape property has been proved. They do not determine \(\Theta_\delta\), identify basal stiffness or prestress, give valid probabilities, or justify an arbitrary finite parameter grid. Completing physical inference still requires a justified observation operator, specimen/reference variation or a common-specimen assumption, dose uncertainty, a declared discrepancy contract, independently constrained geometry/material ranges, accurate branch solutions, and a validated search or enclosure of every remaining parameter combination. One least-squares optimum or a few carried branches is not that set.

The practical advance is a quantitative rejection threshold: a common-path peak model claiming uniformly smaller residual is already incompatible with the published small-only records under its own fixed-dose assumptions. Its originality lies neither in isotonic regression nor in this elementary application. It improves the project's uncertainty accounting and clarifies which new measurements or physical assumptions would actually sharpen inference.
