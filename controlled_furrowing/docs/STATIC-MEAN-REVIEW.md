# Independent review: the conditional static-mean compatibility floor

**Scope verdict: PASS for the source calculations, stated kinematic theorem and actual originating exposition.** The optional exact-isochoric strengthening has also been independently checked, as recorded below. This is an internal AI review, not external expert validation, a new biological result or a proof that the published measurements use the disputed static convention.

The review reads the 17 canonical source-column records directly from `data/derived/calibration_small_18um.csv` as exact decimal fractions. It imports none of the production data-loading, correction or matching functions. `evidence/static_mean_compatibility_independent_review.json` retains **149 passing checks**, all 17 independently formulated linear programs, exact constructive checks, independent matching scans and source hashes.

## Lower bound and sharp correction

With common spatial ROI, references and normalized measure at one static state, the signed regional means satisfy the identity

\[
 V=A+B.
\]

If all four observed channels have a common absolute discrepancy allowance \(\delta\), any compatible triple must obey

\[
 |A_{\rm obs}+B_{\rm obs}-V_{\rm obs}|
 \le |A_{\rm obs}-A|+|B_{\rm obs}-B|+|V_{\rm obs}-V|
 \le3\delta.
\]

Thus the per-record lower bound is \(|A_{\rm obs}+B_{\rm obs}-V_{\rm obs}|/3\). More generally the necessary bound is the sum of the three individual channel allowances; no equal-noise assumption is inferred.

Set \(t=(A_{\rm obs}+B_{\rm obs}-V_{\rm obs})/3\) and define

\[
 A'=A_{\rm obs}-t,\qquad B'=B_{\rm obs}-t,\qquad
 V'=V_{\rm obs}+t,\qquad M'=\max(M_{\rm obs},V').
\]

Then \(V'=A'+B'\). If the observations satisfy \(V_{\rm obs}\le M_{\rm obs}\), the M correction is zero for nonpositive t, and otherwise is at most t. All four residuals are therefore at most \(|t|\), with equality in the three mean channels. This exactly solves the four-channel L-infinity feasibility problem, provided the corrected maximum also satisfies the height requirement of the kinematic witness below.

Independent HiGHS linear programs use variables \((M,V,A,B,\delta)\), constraints \(|x_k-y_k|\le\delta\), \(V=A+B\), \(V\le M\), signed A/B/V, nonnegative delta and \(M\le33.5-10^{-6}\). All 17 optima match the exact per-record fractions to less than 1e-9 micrometres. The exact constructive fractions, rather than the floating LP solutions, supply the sharpness claim.

The cohort floor is

\[
 \delta_{\rm canonical}=\frac{43438757}{12500000}
 =3.47510056\;\mu\mathrm m.
\]

The unique worst record is canonical column G: \(A_{\rm obs}=13\), \(B_{\rm obs}=3.46648612\), \(V_{\rm obs}=6.04118444\), \(M_{\rm obs}=12.73324306\), all in micrometres. Its corrected triple is \((A',B',V')=(9.52489944,-0.00861444,9.516285)\). Signed basal motion is essential to this particular sharp witness; imposing an additional nonnegativity restriction would define a different problem.

## A genuine kinematic witness, with an explicit smooth extension

Fix a declared physical ROI radius R and height H. For either weighting separately, choose constant apical downward displacement A' and

\[
 S_0(r)=M'-c(r/R)^2,
 \qquad
 c=\begin{cases}3(M'-V'),&\text{line mean},\\2(M'-V'),&\text{disk mean}.\end{cases}
\]

The normalized second moments of r/R are 1/3 and 1/2 respectively. Hence the height mean is V', its maximum is M', and the basal upward field \(S_0-A'\) has mean B'. Since \(M'\ge V'\), c is nonnegative and \(S_0\le M'<H\).

To extend beyond the ROI into a larger clamped disk, choose a smooth radial cutoff \(0\le\chi\le1\), equal to one through the ROI and zero in a neighborhood of the outer rim. Continue the same quadratic \(S_0\) outside the ROI, and set

\[
 S=\chi S_0,\qquad D_a=\chi A',\qquad D_b=S-D_a.
\]

Because \(S_0\le M'\), multiplication by chi gives \(S\le\max(M',0)<H\) everywhere, even if \(S_0\) becomes negative. The fields are smooth at the axis because chi is constant there and the core dependence is on \(r^2\). The finite-height axisymmetric map

\[
 R_{\rm current}=r,\qquad
 Z_{\rm current}=D_b(r)+\left(1-\frac{S(r)}H\right)z
\]

has determinant \(J=1-S/H>0\), is injective because radius is fixed and each vertical fibre map is strictly increasing, and equals the identity near the outer rim. This supplies the missing global positive-thickness/clamped-rim argument without a mechanical solve.

For the actual canonical corrections at H=33.5 micrometres, the largest corrected M is 12.73324306 micrometres, and the smallest core J lower bound is approximately 0.619903. These original radial-identity witnesses are not nearly incompressible equilibria, force-balanced solutions, common-material fits or admissible optogenetic activation histories. The separately verified isochoric strengthening below removes the compression caveat for the wider kinematic class. Their purpose remains sharpness in the **larger class of independent kinematic fields**. The lower bound still constrains any smaller mechanical class sharing the static identity, but that class may require a larger discrepancy or be infeasible.

The displayed line and disk witnesses are different fields. The construction proves attainment for each fixed weighting separately; it does not display one identical field simultaneously attaining both observation conventions.

## Unknown intact component-pair permutation

The alternative identity contract permits one cohort-wide permutation of the paired apical/basal columns relative to the maximum/mean records. It does not permit separate shuffling of apical and basal channels. Let \(x_i=A_i+B_i\) and \(y_j=V_j\). Its bottleneck objective is

\[
 \min_\pi\max_i\frac{|x_{\pi(i)}-y_i|}{3}.
\]

Sorting both one-dimensional lists is optimal. To check the usual uncrossing argument, suppose \(x\le x'\), \(y\le y'\), and crossing pairs have maximum distance D. Then \(x-y\le x'-y\le D\), \(y-x\le y'-x\le D\); similarly both signs of \(x'-y'\) are bounded by D. Removing a crossing cannot increase the bottleneck. Repeated uncrossing yields sorted matching.

The independently computed optimum is

\[
 \delta_{\rm permutation}=\frac{895512391}{300000000}
 \simeq2.98504130333\;\mu\mathrm m.
\]

An independent compiled Hopcroft-Karp algorithm scanned **206 critical exact thresholds** through the first perfect matching. Every lower critical threshold was infeasible. The retained Hall witness immediately below the optimum has 16 left records but only 15 possible neighboring component pairs; its neighbor set was recomputed directly from exact costs. The preceding threshold is 74425479/25000000 micrometres. Twenty additional exact fixtures of sizes two through six independently compare sorted matching against enumeration of every permutation.

Each selected pairing also passes the exact correction and positive-height primal construction. Again, this weakens the source identity contract; it neither diagnoses an error in the authors' ordering nor provides a common mechanical parameter vector.

## Essential interpretation limits

The identity fails for separately maximized temporal means. At two times let the component means be (1,0) and (0,1): their separate maxima sum to two, while the maximum total mean is one. Pointwise component maxima likewise need not add to the total maximum. The source observation audit explicitly retains these unresolved alternatives.

Therefore this result excludes a **conditional static regional-mean interpretation below its discrepancy floor**. It does not refute the measurements, establish an experimental noise level or refute every dynamic/mechanical explanation. All 17 canonical records satisfy the necessary pointwise-summary inequalities mean <= maximum <= apical maximum + basal maximum; this alone does not determine the correct source convention.

The mathematics is established additive compatibility, convex feasibility and one-dimensional bottleneck matching. This is a useful exact audit of the current interpretation, with no confirmed novelty claim.


## Independently checked exact-isochoric strengthening

The actual `STATIC-MEAN-COMPATIBILITY.md`, `scripts/static_isochoric_construction.py` and its 204-check receipt were read. The wording explicitly distinguishes static interpretation from temporal peaks, canonical pairing from crosswalk sensitivity, independent configurations from common mechanical equilibria, and different fields for line/disk weighting. The requested precision refinement has been incorporated: the originating exposition now states determinant one **almost everywhere**, with no classical gradient asserted at the core/annulus interface.

The original 149-check review receipt is retained unchanged as a historical record. The new `evidence/static_isochoric_independent_review.json` contains **751 passing independent checks** across 68 configurations. It reconstructs every corrected record from the CSV, expands and integrates the polynomials directly, checks endpoint/critical-point extrema and derives the squared radial-transport polynomials independently. None of the production correction, matching or isochoric-construction functions is imported.

For ROI radius R, outer radius \(R_o\), \(d=R_o-R\), core edge value e and annular coordinate \(u=(r-R)/d\), write

\[
 S_{\rm ann}(u)=e(1-u)-Ku(1-u).
\]

Direct expansion of \(r=R+du\) gives

\[
 \int_R^{R_o}rS_{\rm ann}(r)dr
 =e\,d\left(\frac R2+\frac d6\right)
 -K\,d\left(\frac R6+\frac d{12}\right).
\]

The core integral is \(I_c=R^2(M/2-c/4)\). The specified K therefore gives exact zero total \(\int_0^{R_o}rS\,dr\). Independent quadratic extrema check positive height throughout each core and annulus; the smallest height across all 68 is exactly 20.76675694 micrometres.

Put \(h=H-S\) and define

\[
 \rho(r)^2=\frac2H\int_0^r xh(x)dx.
\]

For positive r, h>0 gives \(\rho'(r)=rh/(H\rho)>0\), so it has a strictly increasing inverse \(r=\mathcal R(\rho)\). Total zero shortening integral implies \(\rho(R_o)^2=R_o^2\). Near the axis, the core has

\[
 \rho(r)^2=\frac{H-M}{H}r^2+\frac{c}{2HR^2}r^4.
\]

Its leading coefficient is strictly positive. Thus the inverse radial scale tends to \(\sqrt{H/(H-M)}\), the radial Cartesian map is regular at the axis, and no singular volume concentration is introduced there.

For current basal height B(r), use the reference-cylinder map

\[
 (\rho,\theta,z)\mapsto
 \left(\mathcal R(\rho),\theta,
 B(\mathcal R(\rho))+\frac{h(\mathcal R(\rho))}{H}z\right).
\]

On every smooth piece its cylindrical Jacobian is

\[
 J=\frac{\mathcal R}{\rho}\mathcal R'\frac hH
 =\frac r\rho\frac{H\rho}{rh}\frac hH=1.
\]

The core/annulus fields are continuous, with bounded piecewise derivatives and height bounded away from zero. Their inverse radial map and full three-dimensional map are Lipschitz; determinant one holds almost everywhere. Radial injectivity and positive vertical scaling give global injectivity. At the outer rim, \(\rho=R_o\), B=0 and h=H, so every rim point is fixed. Gradient identity or traction balance at that rim is not claimed.

Independent Brent inversion plus central finite differences of the actual composed map check four noninterface radii in each configuration. The largest numerical determinant error is about 1.33e-10 against a declared 2e-7 tolerance. This numerical corroboration is separate from the analytic Jacobian cancellation and exact polynomial certificates.

This strengthening provides exactly volume-preserving **continuum kinematic** witnesses. They have finite neo-Hookean/surface energy on the bounded cylinder: the deformation gradient is bounded, J is one almost everywhere and surface slopes are bounded. They are not force-balanced configurations, exact P2 fields or realizations of the project's activation mask and material parameters. In particular, radial transport generally maps a different reference radius onto the fixed current ROI. Separate fields are still used for each record and weighting; a shared equilibrium/material/activation construction has not been supplied.

The originating exposition's main scope is sound. The result excludes one conditional static-mean interpretation below its exact discrepancy floor, while its volume-preserving sharpness belongs to a broader kinematic class. It leaves mechanical compatibility, source temporal semantics, joint material inference and biological validation unresolved.
