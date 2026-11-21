# A corrected point certificate: what survives the interval-arithmetic failure

The corrected computation again proves a restricted result: for one saved
finite-element field and its original parameters, the **exact-spatial-integral
Newton-like map sends the old center outside the old certified cube**. The
earlier point-integration receipts remain withdrawn. This statement now rests
on newly computed safe intervals and new independent checks; agreement with
the earlier displayed numbers is not a proof premise.

It does **not** prove that an equilibrium fails to exist, that an equilibrium
lies outside that cube, or that the biological model has been validated.

![Coordinate projection and reverse-triangle bound, with the corrected proof chain.](../runs/exact_integration_bridge_v1/safe_v2_publication_preparation/safe_point_restoration.png)

*The teal interval is one coordinate projection of the old radius cube. The
gold interval encloses that coordinate of the new map image of the old center;
it is not an equilibrium enclosure. The right panel shows the certified lower
displacement relative to the old radius. All plotted coordinates are normalized
displacements, scaled by the old radius, with decimal display approximations.
The underlying comparisons use retained exact fractions. The diagram is
original and generated from the corrected receipts.*

## The two mathematical targets

This concerns the project's model of controlled tissue folding. Bulk
elasticity, energies on the two opposed boundary surfaces, and a prescribed
local contractile drive enter a mechanical energy. The present update checks
the numerical mathematics of that model; it adds no microscopy observation or
molecular mechanism.

The controlled-furrowing model represents an axisymmetric epithelial body by
quadratic finite elements. At the selected resolution there are 611 free
displacement coordinates. Fixed material geometry, piecewise element
activation, material constants and boundary conditions define a finite
dimensional energy. Two implementations of its spatial integration must be
distinguished:

- \(E_Q\) is the previously certified finite sum, using the actual stored binary
  quadrature nodes, weights and strain-map coefficients.
- \(E_I\) integrates the coherent P2 displacement field over the reference
  elements and surfaces. Stored binary input coefficients are interpreted as
  exact rational constants. The integral itself is bounded with validated
  moments and constitutive-series remainder bounds.

An exact spatial integral of a fixed finite-element ansatz is still a finite
dimensional target. It is not the un-discretized elasticity PDE.

Write \(F_Q=\nabla E_Q\), \(F_I=\nabla E_I\), and let \(R\) be the same stored
binary approximate inverse, also interpreted exactly. The original field is
\(u_0\), at the original base parameter \(\theta_0\) for source record B. The
old cube is

\[
\mathcal B_r(u_0)=\{u:\|u-u_0\|_\infty\le r\},
\qquad r=10^{-11}.
\]

The original, independently checked dense-quadrature proof bounds
\(\|R F_Q(u_0)\|_\infty\). A different spatial target needs a new error bridge;
a small floating difference between quadrature orders is insufficient.

## The arithmetic flaw and its repair

Sparse polynomials store only coefficients thought to be nonzero. The earlier
helper retained a coefficient with `if v != 0`. For an Arb ball containing zero
and nonzero values, both `v == 0` and `v != 0` can be false. These tests certify
relations between uncertain quantities; they are not logical complements.
Consequently, the old cleanup could discard a coefficient whose true value
was nonzero. The simplest example is a ball centered at zero containing 1/2:
turning that coefficient into an omitted term incorrectly replaces its allowed
value 1/2 by zero.

The failure was observed in the actual new uniform-field construction:
ten uncertain coefficients were dropped across nine cleanup calls. The audit
did not establish observed corruption of every historical point result.
Nevertheless, the old full point certificates could not be relied upon: their
local enclosure engine was unsound, and their independent integrations reused
that same arithmetic. Hash agreement and two implementations sharing a bug do
not establish enclosure validity.

The new version removes only a coefficient **certified exactly zero**, using
`is_zero()` for finite Arb balls, and rejects nonfinite coefficients. All old
sources and receipts are preserved. An independently checked safe target and
jet helper pass 992 checks, including adversarial interval-coefficient cases.
The additional signed-cancellation fixture sets

\[
e_0=\frac1{5001}+2^{-1000}.
\]

For the degree-one expansion of the relevant volumetric derivative, the
constant coefficient is exactly \(5001\,2^{-1000}\). At 128-bit working precision
its computed ball can cross zero. The corrected polynomial retains it and
encloses its exact value. Silently deleting it is forbidden regardless of its
small size.

## Fresh calculations and their independent checks

The new producer recomputes all 144 bulk elements and 24 surface edges at the
same original field and parameters. It does not reuse historical point-error
outputs as validated numerical contributions. The normalization, actual rules,
stored physical inputs and constitutive laws remain fixed.

The integration method contracts **signed** integral-minus-rule moments of
the truncated constitutive polynomials. It retains nonzero low-degree defects
of the actual binary rule rather than assuming ideal Gaussian exactness. Force
and Hessian remainders are bounded separately with both positive masses: the
exact reference-integral mass and the sum of actual positive auxiliary weights.
The finite canonical-map-versus-cached-map correction is evaluated at the same
exact stored nodal field.

The independent verification has two distinct layers:

- **Independent local integration.** Six predeclared bulk elements are
  reintegrated at degree 18 and 640 bits. The reviewer builds the Jacobian jet
  from two signed cubic monomials, generates constitutive coefficients by a
  separate convolution, rebuilds canonical barycentric maps, and performs
  manual point lifts. All 60 force and 343 upper-Hessian interval comparisons
  pass, together with 19,200 canonical basis-entry checks. The runtime receipt
  contains 32 checks. Separately, six surface edges at degree 28/640 bits give
 30 force and92 upper-Hessian comparisons and a 349-check review receipt.
- **Complete exact reduction.** Every one of the 168 new local packets is
  checked for identity, hashes, free-coordinate maps, domains and order/tail
  policy. Serialized interval endpoints are added with exact rational
  arithmetic to reconstruct 611 force coordinates and 6676 symmetric
  upper-Hessian slots. All 611 signed products by the exact binary \(R\) are
  reduced independently. This layer passes 208 checks and retains 355 input
  bindings.

These are not 168 independent reintegrations. The remaining local packets rely
on the reviewed generic safe engines, their domain/tail checks, and the
complete exact reduction. The local replays explicitly reuse the reviewed
moment, tail and low dimensional constitutive primitives. Evidence counts
describe different checks; they are not independent statistical trials.

## The new bound and proof

Let \(D=R(F_I-F_Q)(u_0)\). The corrected intervals prove the following
conservatively rounded bounds:

\[
1.673692093\times10^{-11}
\le \|D\|_\infty
\le 1.673692097\times10^{-11}.
\]

The old dense-quadrature residual premise is still valid. The review uses the
larger of the author and independent upper bounds, not the more favorable
number:

\[
\|R F_Q(u_0)\|_\infty\le\eta,
\qquad \eta\le5.287455044\times10^{-13}.
\]

The base parameter is verified to belong to the same original closed parameter
cell. Its affine predictor equals \(u_0\) there, so this residual premise
applies at exactly the point used by the new calculation. The reverse triangle
inequality gives

\[
\begin{aligned}
\|R F_I(u_0)\|_\infty
&\ge \|D\|_\infty-\|R F_Q(u_0)\|_\infty\\
&\ge L-\eta\\
&\ge1.6208175\times10^{-11}>r,
\end{aligned}
\]

where \(L\) and \(\eta\) are the retained exact-fraction bounds. Therefore

\[
T_I(u)=u-R F_I(u)
\quad\Longrightarrow\quad
T_I(u_0)\notin\mathcal B_r(u_0).
\]

A fresh coordinate witness is free index 588, node 313, axial component. Its map
displacement is enclosed approximately by
\([-1.72656665,-1.62081754]\times10^{-11}\); the source receipt stores the exact
endpoints. That coordinate alone establishes escape from the old cube.

The conclusion is about the **self-map property**, a premise in a Banach
fixed-point proof. It does not establish root exclusion. For example,
\(T(x)=6/5-x/2\) sends the center 0 of \([-1,1]\) outside the interval, but has
the fixed point \(4/5\) inside it. A different center, radius, preconditioner,
or proof norm may succeed for the integral target; that remains further work.

## What remains valid, and what remains unproved

The 17 dense finite-quadrature root certificates and their 611-coordinate
axisymmetric local-energy-minimum proofs concern \(E_Q\). They do not use the
faulty sparse polynomial path and remain valid within their declared scope.
The separately checked exact rational moment formulas and pure interval
whole-body geometry results also remain valid. Geometry of a box does not
imply that an integral equilibrium belongs to that box.

The new point certificate supplies no whole-tube bounds on
\(F_I-F_Q\) or \(\nabla F_I-\nabla F_Q\), no integral-root existence or stability,
and no inference over the full physical parameter/geometry domain. It does not
close the source temporal-observation/registration gaps, establish a shared
dynamic path between doses, prove a continuum PDE result, or validate the
model biologically. No novel biological concept or new general mathematical
principle is claimed. The methods are standard validated numerical analysis
applied to this explicitly defined model; their value here is a reproducible
statement with precise limits.

The next mathematical operation is to use the corrected arithmetic for uniform
integration derivative bounds on a freshly verified comparison domain, then
attempt an appropriately centered integral-root enclosure. The original
five-gate objective remains incomplete. New point certification does not
replace full parameter-set inference, geometry/branch/tail coverage,
continuum analysis or source-compatible prediction and validation.

## Retained receipts

| Evidence | New receipt or artifact | Hash prefix |
| --- | --- | --- |
| Arithmetic failure and scoped withdrawal | `evidence/interval_polynomial_cleanup_audit.json` | `fe0f5414` |
| Independently reviewed safe polynomial kernel | `evidence/exact_safe_polynomial_kernel_independent_review.json` | `a2f6b6c4` |
| New six-element bulk replay | `evidence/exact_bulk_independent_replay_safe_v2.json` | `37515a8d` |
| New six-edge surface replay | `evidence/exact_surface_safe_independent_review.json` | `248ba500` |
| Corrected complete point reduction | `evidence/exact_complete_point_B_independent_review_safe_v2.json` | `433125c7` |
| Fresh coordinate witness | `runs/exact_integration_bridge_v1/bulk_independent_replay_safe_v2/coordinate_witness.json` | `98a70c84` |
| New component delivery bindings | `runs/exact_integration_bridge_v1/bulk_independent_replay_safe_v2/delivery.json` | `68383583` |

The preparation receipt binds full hashes, exact plotted values, current source
files and figure-generation code. Historical unsafe point receipts remain
withdrawn and are retained for the development record.

This corrected point result is now independently reviewed. The subsequent [uniform integral-root theorem](EXACT-INTEGRAL-ROOT.md) adds full-domain derivative and residual bounds; it proves more than the point certificate alone. Historical unsafe receipts remain withdrawn.
