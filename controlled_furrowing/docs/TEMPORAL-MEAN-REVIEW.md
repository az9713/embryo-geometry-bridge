# Independent review: separate temporal regional-mean peaks

**Verdict: PASS for the stated continuous kinematic existence construction.** Under the explicitly defined peak-of-regional-mean operator, all 17 canonical source vectors have exact, volume-preserving geometric paths with the recorded four summaries. This does not construct mechanical dynamics or identify the experimental observation algorithm.

The independent review reconstructs all **68 endpoint configurations** from the CSV with exact fractions, without importing the production loading or construction functions. `evidence/temporal_mean_independent_review.json` retains 942 passing numerical/algebraic/provenance checks, the independent replay program, exact path coefficients and source hashes. The continuous-time statement follows from affine extrema and uniform geometric bounds, rather than sampling a few times.

## Exact component peaks along a continuous path

For each source vector, the recorded values satisfy

\[
 M,V,A,B\ge0,\qquad V\le M<H,\qquad V\le A+B.
\]

At normalized path times 0, 1 and 2, use instantaneous component means

\[
 (0,0),\qquad(A,V-A),\qquad(V-B,B).
\]

The corresponding total means are 0, V, V. The instantaneous identity is respected at every state. Interpolate these means affinely on each time segment. Then

\[
 \sup_t\overline a=A,\qquad
 \sup_t\overline b=B,\qquad
 \sup_t\overline s=V,
\]

because \(V-B\le A\), \(V-A\le B\), and A, B, V are nonnegative. Each mean is affine on each closed segment, so its continuous-time extrema occur at an endpoint. The two component peaks need not coincide in time. Negative values at the other component's peak are allowed; imposing nonnegative signed motion at every time would change this admissibility class.

For either fixed weighting separately, put

\[
 S_*(r)=M-c(r/R)^2,
 \qquad c=3(M-V)\;\text{(line)},\qquad
 c=2(M-V)\;\text{(disk)}.
\]

Inside the fixed current-radius ROI, let the apical downward motion be spatially constant with the prescribed instantaneous mean, and set basal upward motion to \(b=S-a\). On the first segment, \(S(t,r)=tS_*(r)\); on the second, \(S=S_*\). Since M is nonnegative, the peak spatial height shortening is exactly M. The weighting-specific polynomial moments give mean V at both nonzero endpoints. Thus all four recorded summaries are matched with zero discrepancy under this particular operator.

## Exact volume and all-time geometric admissibility

The compensated annulus is computed independently for both endpoint component partitions. Its shortening polynomial and radial transport are identical at the two nonzero states; only their vertical offset differs. Expanded polynomial integration checks zero total \(\int rS\,dr\), and exact endpoint/critical-point extrema check positive height for all 68 configurations.

For every actual source construction, the **global** maximum shortening, including the compensation annulus, equals the ROI value M. The smallest global height is exactly 20.76675694 micrometres. On the first segment, core and annular shortening scale together, so

\[
 \int rS_t\,dr=t\int rS_*\,dr=0,
 \qquad
 h_t=H-tS_*=(1-t)H+t h_*>0
 \quad(0\le t\le1).
\]

The uniform lower bound is \(\min(H,\inf h_*)\), so positivity is proved at **every time**, not inferred from a finite temporal grid. The radial volume coordinate

\[
 \rho_t(r)^2=\frac2H\int_0^r xh_t(x)dx
\]

is strictly increasing, continuous in t, regular at the axis and equal to \(R_o\) at the outer rim. Its inverse gives an embedded body map

\[
 (\rho,\theta,z)\mapsto
 \left(R_t(\rho),\theta,
 b_t(R_t(\rho))+\frac{h_t(R_t(\rho))}{H}z\right)
\]

with determinant one almost everywhere, as in the [independently checked isochoric proof](STATIC-MEAN-REVIEW.md). At t=0 the map is the undeformed reference. At t=1 the two time pieces agree. For \(1\le t\le2\), height and radial transport are fixed while the vertical offset varies continuously and affinely; that offset does not change the triangular cylindrical Jacobian. Rim displacement remains zero at every time.

The spatial fields are continuous with bounded piecewise derivatives; their interface gradient jumps are permitted in this Lipschitz kinematic class. The axis limit remains regular throughout the path because \(H-tM\) stays uniformly positive. No exact P2 representation, energy minimization, viscous force balance, acceleration law or physical timescale is supplied.

## The order of spatial averaging and temporal maximization matters

The zero-error result uses

\[
 V=\sup_t\left\langle S(t,\cdot)\right\rangle.
\]

It does not generally reproduce \(\langle\sup_t S(t,\cdot)\rangle\). Because the path starts from zero and then scales to \(S_*\), its pointwise temporal peak is \(\max(0,S_*)\). For 29 of the 34 source-record/weighting height paths, the final profile is negative near the ROI edge. Their mean of local temporal peaks is strictly greater than V.

For example, when \(c>M\), direct integration of the positive part gives

\[
 \left\langle\max(0,M-cu^2)\right\rangle
 =\begin{cases}
 \dfrac{2M}{3}\sqrt{M/c},&\text{line},\\
 \dfrac{M^2}{2c},&\text{disk}.
 \end{cases}
\]

These exceed the corresponding \(V=M-c/3\) or \(M-c/2\). The exact rational condition c>M determines the count; the receipt retains the analytic expressions and numerical illustrations. This reinforces the need to resolve the source's time/mean ordering rather than selecting whichever operator fits.

## Interpretation and review outcome

The actual originating exposition distinguishes peak regional means from mean local peaks, declares arbitrary normalized time, and disclaims a common constitutive law, gain, activation history or experiment-specific dynamics. Separate paths are used for each record and weighting; a single identical path is not claimed to realize both weights.

The static 3.47510056 micrometre floor and this zero-error temporal witness are compatible because they constrain different observation contracts. The latter proves that the static additive exclusion cannot simply be transferred to separately timed regional-mean peaks. It does not establish that the source used this algorithm or that any chosen mechanochemical model can produce the paths. The canonical cross-panel ordering remains an explicit identity assumption.

This is an exact kinematic counterexample and an observation-contract audit. It supplies neither a common fitted material system nor a new biological mechanism, physical uncertainty set or confirmed novel concept.
