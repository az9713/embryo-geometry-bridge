# From a verified nodal error to the measured static endpoints

**A bound on nodal coordinates can be transferred rigorously to the current-radius observer. No actual mechanical error radius is supplied by this result.** It is one conditional link needed between a future algebraic-root certificate and observable uncertainty. Quadrature, continuum, temporal and experimental errors remain separate.

The [implementation](../scripts/nodal_observer_error_transfer.py) and [thirteen exact synthetic checks](../evidence/nodal_observer_error_transfer_tests.json) use rational arithmetic. An independent review is still required before this component is used to support a physical error claim. An initial generic API gap allowed composition of different default surface ROIs; the corrected implementation requires a shared ROI. Independent review also required positive ROI width before claiming bounds for normalized means. The preceding sources, tests and findings are retained under `runs/nodal_error_transfer/initial_revision/` and `pre_mean_mass_guard/`, before any application to an actual mechanical error bound.

## Trace and derivative bounds

For one P2 surface edge, parameterized by \(0\le t\le1\), the nodal basis is

\[
 N_0=(1-t)(1-2t),\quad N_1=4t(1-t),\quad N_2=t(2t-1).
\]

Suppose corresponding nodes differ by at most \(\varepsilon\) in each physical coordinate. Then

\[
 \|\widetilde R-R\|_\infty,\ \|\widetilde Z-Z\|_\infty
 \le\frac54\varepsilon=:e,\qquad
 \|\widetilde R'-R'\|_\infty\le8\varepsilon.
\]

For \(t\le1/2\), \(\sum|N_i|=1+2t-4t^2\le5/4\); symmetry handles the other half. Equality holds at \(t=1/4\) with suitable nodal signs. The sum of absolute derivatives is convex in t and has maximum eight at an endpoint. Both constants are sharp for unconstrained edge nodal errors; additional fixed-node constraints can improve them.

Let the exact base edge have \(\min R'>8\varepsilon\). The perturbed edge remains strictly radially increasing. Shared nodes preserve conformity. Fixed endpoint radii, or a sufficiently interior ROI, guarantee a common radial graph domain. These are verified prerequisites or explicit assumptions—not conclusions obtained from a force residual.

## Comparing at the same current radius

Let \(L\) bound the Lipschitz constant of the base graph \(Z=G(R)\). One conservative exact choice is the maximum over edges of \(\max|Z'|/\min R'\). For a perturbed point at radius r, write \(\widetilde R(t)=r\). Then

\[
 |\widetilde Z(t)-G(r)|
 \le |\widetilde Z(t)-Z(t)|
 +|G(R(t))-G(r)|
 \le e+Le=e(1+L).
\]

This accounts for radial displacement as well as vertical displacement; using only vertical nodal error would compare different physical points. The estimate uses the base graph's global Lipschitz bound, so it remains valid when the inverse radius crosses an edge boundary.

For apical and basal errors \(E_a,E_b\), the static component maxima and regional means have errors at most \(E_a,E_b\), respectively. Same-radius height, its maximum and its regional mean have error at most \(E_a+E_b\). The mean bound holds for any fixed normalized nonnegative weight; it therefore covers the declared line and disk weights. Maxima obey \(|\sup f-\sup g|\le\|f-g\|_\infty\).

## Connecting a future root ball

If a validated finite-dimensional root lies in a free-displacement cube of radius r around the archived vector, fixed physical scale \(L_0\) gives nodal error at most \(L_0r\) relative to the exact sum of the reference nodes and that vector. The stored `current_nodes` are rounded binary sums, so comparison with the existing observer additionally requires the exact reference-plus-displacement rounding gap. It is small but cannot be silently dropped in a rigorous certificate.

The root ball must separately preserve the fixed radial boundary conditions and admissible positive geometry. An algebraic root certificate for a fixed quadrature potential would still not establish continuum discretization accuracy, physical stability, experimental observation semantics or a uniform parameter-domain bound. No hypothetical chosen \(\varepsilon\) is presented as an actual error estimate in this module.
