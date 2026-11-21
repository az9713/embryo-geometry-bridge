# Independent review: nodal coordinate errors to static endpoint errors

**Verdict: PASS for the current conditional transfer theorem and its revised API. No actual mechanical error radius is certified by this review.** A verified or explicitly assumed nodal bound is required as input. Force residuals, agreement between meshes and higher-quadrature checks are not substituted for that bound.

The proof constants, common-domain conditions, cross-edge comparison and endpoint composition were checked independently. `evidence/nodal_observer_error_independent_review.json` retains **80 passing exact/control checks** before final documentation checks, including 64 independent radial/vertical sign corners, a deliberate edge-crossing example and the revised ROI guards. The mathematical proof covers all points; the finite control test grid is corroboration rather than that proof.

## Exact trace constants

On each corresponding P2 edge, use local parameter t in [0,1] and nodal functions

\[
 N_0=1-3t+2t^2,\qquad N_1=4t-4t^2,\qquad N_2=-t+2t^2.
\]

Assume every corresponding physical nodal coordinate differs by at most epsilon, separately in R and Z. For t<=1/2, the basis signs give

\[
 \sum_i|N_i|=1+2t-4t^2\le\frac54.
\]

Reflection handles the other half. At t=1/4 the signs (+,+,-) attain 5/4. Also \(\sum_i|N_i'|\) is a convex piecewise-linear function of t; its endpoint values are eight, hence its maximum is eight. The radial error signs (-,+,-) attain the derivative constant at t=0. Therefore

\[
 \|\widetilde R-R\|_\infty,
 \|\widetilde Z-Z\|_\infty\le e:=\frac54\varepsilon,
 \qquad\|\widetilde R'-R'\|_\infty\le8\varepsilon.
\]

These two constants are separately sharp for unconstrained edge nodal errors. The final composed observer bound is conservative and is not claimed sharp. Fixed-node constraints can improve constants but are not exploited here.

If every base edge has minimum normalized-parameter radial derivative m>8epsilon, every perturbed edge remains strictly radially increasing. Shared global nodes must retain conformity, so the assembled perturbed trace is continuous and increasing. These are assumptions about corresponding P2 fields; a force residual alone does not prove them.

## Why crossing an edge does not invalidate the bound

Let G be the continuous base Eulerian graph on its full radial domain. The exact code uses the global bound

\[
 L=\max_{\rm edges}
 \frac{\max_{t\in[0,1]}|Z'(t)|}{\min_{t\in[0,1]}R'(t)}.
\]

Each numerator and denominator extremum occurs at an endpoint of its linear derivative. The resulting piecewise bound gives one global Lipschitz constant across the continuous joined graph.

For any r in a common covered physical ROI, select the perturbed edge and local parameter with \(\widetilde R(t)=r\). The corresponding base R(t) is in the full base domain, although it may lie outside the ROI or on the other side of a graph knot. Then

\[
 |\widetilde Z(t)-G(r)|
 \le |\widetilde Z(t)-Z(t)|+|G(R(t))-G(r)|
 \le e+Le.
\]

The proof does not assume that the two inverse maps choose the same edge. Restricting L to the selected base edge or an unexpanded ROI would not be justified in general.

An independent exact example uses two base edges: flat height on R in [0,1], then slope four on [1,2]. Move their shared radial knot to 1+epsilon and each adjacent radial midpoint by epsilon/2, preserving outer radii and all heights. At r=1+epsilon/2, the perturbed inverse lies on the first edge while the base inverse lies on the second. The actual height error is 2epsilon. A same-edge slope-zero estimate of 5epsilon/4 would fail, but the implemented global-Lipschitz bound succeeds.

## Common domain and positive mean mass

If the two global endpoint radii are independently known to remain fixed, the whole base radial domain is inherited. The boolean API flag asserts this external fact; the function does not verify it from unknown perturbed nodes.

Without that assertion, the nodal error bound supplies a guaranteed common interior. Using the conservative e, an ROI contained in \([R_{\min}+e,R_{\max}-e]\) is sufficient. A full-domain request without fixed endpoints is correctly left unresolved. A boundary inverse that merely exists for the base trace is insufficient.

Apical and basal endpoint errors may be combined only on one identical fixed ROI. The current API rejects composition of differing default full-surface ROIs, while allowing an explicitly supplied common covered interval. The mean interface also rejects zero-width common ROIs, because line/disk normalization would have zero mass. A single-curve pointwise bound on a singleton remains meaningful and is kept available separately.

The original same-ROI composition gap was preserved in `runs/nodal_error_transfer/initial_revision/`. The reviewer then exposed the zero-width mean-mass case; its source/tests/finding are preserved in `runs/nodal_error_transfer/pre_mean_mass_guard/`, and the current version has thirteen passing source checks. Neither historical gap invalidated the current project's positive-width common ROI, but both generic API cases required explicit handling.

## Endpoint composition and reference/weight assumptions

For two static surfaces with Eulerian errors Ea and Eb on the same ROI, unchanged reference profiles imply

\[
 |\delta d_a|\le E_a,\qquad |\delta d_b|\le E_b,
 \qquad|\delta s|\le E_a+E_b.
\]

For any bounded functions, \(|\sup f-\sup g|\le\|f-g\|_\infty\). A fixed normalized nonnegative measure has operator norm one on bounded functions. Thus component maxima/means obey Ea or Eb, while height maximum/mean obey their sum. Signed shortening is permitted; no positive-part clipping is assumed.

The physical ROI, references and normalized measure must be identical in the comparison. Changing a counterfactual reference, registration, spatial mask or geometry-dependent surface-area weight adds other errors not covered here. In particular, unchanged surfaces with zero nodal error cannot bound an independently shifted reference height. The line/disk weights are fixed functions of current physical radius, which meets the stated condition.

## A future algebraic root cube requires a rounding bridge

Suppose an independently validated root satisfies \(\|u_*-u_0\|_\infty\le r\) in the fixed free-displacement coordinates. The exact reference-plus-centre sum and stored rounded `current_nodes` need not coincide. A sufficient physical-coordinate input is

\[
 \varepsilon=L_0 r+
 \max_i\left|L_0(X_i+u_{0,i})-x_{i,\rm stored}\right|,
\]

with fixed coordinates treated consistently and the last term evaluated with exact arithmetic or rigorous intervals. Omitting it is unjustified even if the numerical gap is very small. The root cube must also preserve the actual radial boundary constraints and admissible geometry.

A finite-dimensional certificate for one quadrature-defined root would still leave quadrature-model differences, continuum discretization, branch/stability issues, temporal observation semantics and experimental uncertainties separate. This component supplies the conditional last link from a nodal bound to static observations; it does not create the missing nodal bound or a physical confidence interval.
