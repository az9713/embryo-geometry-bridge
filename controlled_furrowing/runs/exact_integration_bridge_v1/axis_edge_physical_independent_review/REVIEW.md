# Independent physical-axis representation review

All188 review checks pass. The audit binds the frozen source ebe0cbe0,
protocol f880b96b and preflight ef20fac6, preserves all existing inputs, and
performs no moment integration, cache sum or mechanical solve.

The independent rational polynomial layer derives the barycentric coordinates
from the inverse affine corner matrix and constructs the six quadratic shape
functions directly. It accepts only rational coefficients, so its sparse
cleanup cannot repeat the uncertain-Arb zero-removal bug. On each actual
full-axis-edge triangle72–77, three radial nodal values are exactly zero and
their radial degrees of freedom are constrained. Each remaining radial basis
polynomial is exactly divisible by R. Division is performed before multiplying
by uncertain nodal values.

The independently obtained physical q components and all selected free
derivative columns match the producer coefficient-by-coefficient. They also
match the old safe Duffy representation under an independently constructed
Duffy substitution. These are identities in exact rational polynomial
arithmetic, not floating comparisons. Since q is affine in the free nodal
coefficients, equality at the saved point together with equality of its full
constant free-derivative map proves equivalence throughout the same affine
family, with fixed coordinates retained.

For every component and every monomial in the union of all rational basis
supports, the audit independently computes the exact signed extrema of that
linear coefficient over the complete local nodal box. Each interval lies in
the producer's safe Arb coefficient enclosure. This closes a support omission
loophole that sampling or iterating only over the observed output support
would not address. The audit does not rely on the producer's33 sampled fields
per case for uniformity. All q/B polynomials are ordinary affine physical
polynomials with nonnegative powers. Their determinant J is consequently
cubic at most, and its free-coordinate derivatives remain ordinary
nonnegative-power polynomials through the reviewed generic polynomial jet.

The representation proof justifies a separately reviewed physical-triangle
integration adapter. Such an adapter must use the same original affine
reference triangle and the normalized measure(2π/σ)R dR dζ. Nonnegative
physical moments are integrable including R=0; that radial factor is included
once by normalized_bulk_moment. No additional Duffy Jacobian or inverse-R
positive-radius norm may be introduced. A power-weighted supremum norm must
use the actual physical-coordinate ranges.

This is not a proof that the future adapter gives smaller error bounds. The
first full uniform bound remains retained independently. Any replacement of
its six packets needs a separately frozen same-target adapter and reduction
revision. These identities do not prove an integral root or a PDE result, nor
do they impose additional Cartesian smoothness of the existing axisymmetric
ansatz at the axis. Cells touching the axis at only one vertex are explicitly
excluded; their hoop expression need not share this polynomial divisibility.

Source and receipt are frozen after the successful preparation-only tool
session35474 returned exit0. No producer file or existing proof was changed.
