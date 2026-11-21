# Signed singleton bridge: reduction theorem and limitations

This is a mathematical kernel preparation, not an admission of actual spatial
integration packets. It proves no new physical root until its input hypotheses
are supplied and independently checked. The radius proposed by the current
bridge is `r = 2^-34`; the theorem itself applies to any positive radius.

Let `u0` be the exact-real stored binary displacement, with exact addition to the
reference coordinates. Let `B` be its closed infinity-norm radius-r cube in the
611 free coordinates. The physical singleton parameters, requested rational
drive, constraints and energy normalization must be the same in every operand.
Let `FI` be the gradient of the spatially integrated finite-element potential,
and `FQ` the gradient of the exact-real potential defined by the frozen binary
quadrature arrays. Assume both are continuously differentiable on an open
neighborhood of B. The independently certified analytic domains establish the
premises needed by each admitted local formula; this document does not establish
those domains itself.

Write `DF = FI(u0)-FQ(u0)` and `DH(v) = DFI(v)-DFQ(v)`. Suppose authenticated local
signed intervals, including polynomial truncation tails and canonical-to-binary
cache corrections, enclose these quantities at their stated point/whole-box
scopes. Assemble local contributions by exact addition in their actual free-DOF
maps. Every off-diagonal upper Hessian entry contributes once to each symmetric
global position, including when the local free ordering is permuted.

For the same fixed square binary preconditioner R, assume independently supplied
signed rows enclose `R FQ(u0)` and the whole-box bound satisfies

    sup_{v in B} ||I - R DFQ(v)||_infinity <= ZQ.

The reduction first assembles signed DF and DH, then contracts with R, and only
then takes absolute row sums. Outward Arb matrix arithmetic encloses the exact
products, including all rational-to-ball rounding. Define

    eta >= ||R FQ(u0) + R DF||_infinity,
    b   >= sup_{v in B} ||R DH(v)||_infinity,
    z   = ZQ + b.

The triangle inequality implies `sup_B ||I-R DFI||_infinity <= z`. Bounds obtained
by independent interval entries can be conservative; they remain valid without
assuming that entry extrema can occur simultaneously.

## Existence, uniqueness and posterior error

If `z < 1` and `eta + z*r < r`, let `T(v)=v-R FI(v)`. The segment integral formula
gives `||T(v)-T(w)|| <= z||v-w||` on the convex cube. Moreover

    ||T(v)-u0|| <= ||R FI(u0)|| + z||v-u0|| <= eta+zr < r.

Thus T maps the complete closed cube into itself and is a strict contraction.
Banach's theorem gives one fixed point u*. At any point in B the matrix
`R DFI = I-(I-R DFI)` is invertible by its Neumann series, since its defect norm
is below one. Both square factors, including R, are therefore invertible. The
fixed-point equation consequently implies `FI(u*)=0`. Every zero in B is a fixed
point, so the zero is unique in B. Finally,

    ||u*-u0|| <= eta + z||u*-u0||,
    ||u*-u0|| <= eta/(1-z) < r.

This proves a stationary point, not a local minimum. A same-target positive
definiteness argument is an additional requirement for a minimum claim.

## Excluding this specific cube

Independently of `z<1`, the same derivative defect bound gives
`||R DFI|| <= 1+z`. If a zero v* lies in B, the segment formula implies

    ||R FI(u0)|| = ||R FI(u0)-R FI(v*)|| <= (1+z)r.

A signed interval row that proves a lower bound strictly above `(1+z)r` therefore
excludes all zeros in this cube. This implication does not require R to be
invertible. It says nothing about roots outside this cube. Failure of either
sufficient test is unresolved, not evidence for the opposite conclusion.

## Input authority and verification

`scripts/transfer_ambiguity_exact_bridge_reduction.py` independently implements
signed scatter, rational boundary decisions and outward matrix contraction. Its
eight cold fixtures compare a permuted local scatter and signed matrix product
against direct rational arithmetic, test nondegenerate interval containment,
check strict posterior/exclusion inequalities and reject mixed scopes, missing
Hessian entries, duplicate packets and implicit float/bool rational inputs.

These fixtures do not verify any actual 168-component integral computation.
Before an actual application, a separate caller must bind complete component
coverage, true process exits, source/protocol/plan hashes, independent packet
reviews, exact singleton/radius/free-map identities, the same R, signed cached
anchor and whole-box ZQ. An independent review must then check both reduction
source and actual numerical operands/results.

Even a passing application remains a local finite-element equilibrium result.
It supplies no continuum discretization error, global branch selection, full
compatible parameter set, source-matched observation/noise contract, predictive
success or novelty by itself.
