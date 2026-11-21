# Independent review: robustness to the two height-mean orderings

**Verdict: PASS for sixteen exact kinematic witnesses and the explicitly conditional nonattainment classification of record O.** The construction is robust to two temporal/spatial orders for the **height mean**. Its apical/basal summaries remain separately maximized regional means. It is not a common material fit, mechanical dynamics or a reconstruction of the source's observation code.

`evidence/temporal_order_robust_independent_review.json` retains the independent replay program and **493 passing algebraic/provenance checks** before the final document/figure checks. The review directly parses all 17 CSV records as exact fractions, imports no production loading/construction function, and independently checks all 32 positive-record/weight profiles and eight approximation terms for O. The original failed all-positive-mean attempt is retained; no zero observation is replaced or deleted.

## The two operators and the correct quantifiers

For one fixed current-radius ROI and normalized nonnegative line or disk measure, define

\[
 V_1=\sup_t\int S(t,r)d\nu(r),\qquad
 V_2=\int\sup_t S(t,r)d\nu(r).
\]

Since each instantaneous S is bounded by its pointwise temporal supremum, \(V_1\le V_2\). This inequality does not identify which operator produced a published number.

For each of the sixteen positive-mean records and for **each weighting separately**, the new construction provides one path with \(V_1=V_2=V\). It does not display one identical path realizing both line and disk weighting at once. Components remain

\[
 A=\sup_t\int a(t,r)d\nu(r),\qquad
 B=\sup_t\int b(t,r)d\nu(r),\qquad S=a+b.
\]

The ordering robustness is not transferred to the component channels. In fact, exact baseline checks identify 15 of the 32 constructed profiles with negative final basal motion near the ROI edge. Starting from zero clips those local temporal basal peaks, so their mean exceeds the declared B. This confirms that the originating document's component-order exclusion is necessary.

## Exact nonnegative profiles for positive V

Write \(u=r/R\). The basis bump \(f_p(u)=(1-u^2)^p\) lies in [0,1], equals one at the axis and zero at the edge. For line weighting its moment satisfies

\[
 F_p=\int_0^1(1-u^2)^pdu
 =\frac{2p}{2p+1}F_{p-1}
 =\frac{4^p(p!)^2}{(2p+1)!},\qquad F_0=1.
\]

The recurrence follows from integration by parts. For projected disk weighting, substituting \(v=u^2\) gives \(F_p=1/(p+1)\). Both tend to zero; dominated convergence also proves the limit because the measures are atomless and the axis is a set of measure zero.

For \(0<V\le M<H\), choose \(p\ge2\) with \(F_p\le V/M\) and put

\[
 C=\frac{V-MF_p}{1-F_p},\qquad
 S_*(r)=C+(M-C)f_p(r/R).
\]

The denominator is positive, selection gives C>=0, and V<=M gives C<=M. Hence

\[
 0\le C\le S_*\le M,\qquad S_*(0)=M,
 \qquad\int S_*d\nu=C+(M-C)F_p=V.
\]

Independent factorial moments, binomial coefficient expansion and direct coefficient integration reproduce every retained polynomial and mean. The chosen p is minimal subject to the declared p>=2 rule; its largest value for the actual positive-mean records is 20.

The original two-segment component path is retained. During the first segment S scales from zero to this nonnegative profile; during the second it stays fixed while the apical/basal partition changes. Consequently the pointwise height peak is exactly S*, and both height operators equal V. The all-time component mean peaks follow from the exact endpoint inequalities V-B<=A and V-A<=B; affine segments have no hidden interior extrema.

## Volume, geometry and the boundary outside the ROI

The core is nonnegative only **within the measured ROI**. A fixed-rim isochoric body needs compensating lengthening outside it. The physical radial-volume integral uses r dr, even when the declared observation is a line mean. Independent expansion checks

\[
 I_c=\frac{R^2}{2}\left(C+\frac{M-C}{p+1}\right).
\]

With d=Ro-R and annular coordinate x, set

\[
 S_{\rm ann}=C(1-x)-Kx(1-x),\quad
 I_L=d(R/2+d/6),\quad I_B=d(R/6+d/12),\quad
 K=\frac{I_c+CI_L}{I_B}>0.
\]

Then \(\int_0^{R_o}rSdr=0\). Also \(S_{\rm ann}=(1-x)(C-Kx)\le C\le M\), so the **global** height obeys H-S>=H-M>0. The core/annulus traces agree at the ROI and vanish at the outer rim; their derivative jump is allowed in the stated piecewise-smooth Lipschitz class.

The compensated radial coordinate and body map from the [independently checked isochoric proof](STATIC-MEAN-REVIEW.md) therefore have a regular axis, fixed rim, global injectivity and determinant one almost everywhere. On the initial time segment, compensation scales linearly and height is a convex combination of H and the positive final height. On the second segment the height/radial map is fixed and only the vertical offset changes. These are all-time statements, not a temporal sampling certificate.

Each finite-degree configuration has finite gradients and energy. No common constitutive equilibrium, force history, activation mask, physical timescale, biological spatial cutoff or uniform regularity bound is supplied. Compensation outside the ROI must not be silently included in the mean whose nonnegativity/order equality was proved only on the declared ROI.

## Record O: what is impossible, and what remains possible

Record O stores

\[
 V=0,\qquad M=2.41505646\,\mu\mathrm m,\qquad
 A=2\,\mu\mathrm m,\qquad B=0.41505646\,\mu\mathrm m.
\]

For nonattainment of V2, assume a common fixed ROI/reference, pointwise-zero initial shortening, a continuous spatial field and a nonnegative normalized measure with positive mass on every nonempty relative open interval. Assume the temporal supremum is measurable; the constructed jointly continuous finite-time paths certainly satisfy that condition.

Then \(P(r)=\sup_t S(t,r)\ge0\). If M>0, some time and point have shortening greater than M/2. Spatial continuity yields a nonempty relative open neighborhood U with P>M/2. Full support gives

\[
 V_2=\int P\,d\nu\ge\frac M2\nu(U)>0.
\]

Thus O cannot be attained under V2 in this class, and cannot have one exact path satisfying both V operators simultaneously. This proof also applies to disk weighting at the axis: its density vanishes at one point, but every relative open neighborhood still has positive measure.

It does **not** exclude V1=0<M for signed instantaneous profiles. The earlier parabolic path supplies exactly that counterexample. If nonnegative shortening is imposed at every time as an additional restriction, O is impossible for V1 as well; that is a stronger class, not a consequence of initial zero alone. Likewise, without initial zero, full support, continuity or the same ROI, the stated V2 nonattainment need not follow.

The source audit has not verified these pipeline assumptions or whether a recorded zero is rounded, censored or exact. Accordingly this is a mathematical classification of retained decimal records, not a contradiction of a biological measurement.

## Zero infimum is not exact attainment

For O use \(S_p=M(1-(r/R)^2)^p\). Its maximum stays M, but the two corresponding means \(\varepsilon_p=MF_p\) tend to zero. Each finite p has nonnegative ROI shortening and the exact compensated volume/radial construction. Replace V by epsilon in the two component states: as long as epsilon<=A+B, the A, B and M peaks remain exact, so the entire four-channel maximum residual is epsilon.

The independent review verifies all eight retained terms, at p=8,32,128,512 for line and disk weighting, including their exact moments, component peaks, compensating annulus and positive global height. Dominated convergence establishes the infinite-sequence claim; a finite table alone does not. This limit uses the atomless line/disk measures, not arbitrary full-support measures containing an atom at the peak.

Thus the attainable summary image has a nonclosed boundary in the chosen summary topology: the O vector is approached arbitrarily closely but is not attained in this observation class. Continuity alone supplies no positive discrepancy floor. The proof imposes no uniform spatial-gradient or concentration bound; no claim that a bounded mechanical parameter family can realize this concentrating sequence is made. Nor is any epsilon reinterpreted as an observed biological value or measured noise allowance.

## Review outcome

The current originating exposition correctly distinguishes height-order robustness from component ordering, V2 nonattainment from signed-V1 feasibility, a zero infimum from an attained zero, and kinematic witnesses from common material dynamics. No substantive mathematical error was found in the current code or exact receipt. The initial failure is retained and its zero-mean boundary is now explicitly classified.

The useful conclusion is a limit on what summary compatibility can establish before mechanics and the observation pipeline are specified. It proves neither a common model-compatible parameter set nor a new scientific mechanism or confirmed novelty.
