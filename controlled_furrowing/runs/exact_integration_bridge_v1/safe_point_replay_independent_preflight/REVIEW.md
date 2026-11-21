# SAFE V2 full-point replay: independent launch review

Both prelaunch receipts pass. The bulk receipt contains178 literal-True checks
and the surface receipt56. They each include24 shared semantic/safety checks,
so these totals are not disjoint evidence counts.

The review is difference-focused. The new point libraries import the separately
reviewed safe target and polynomial-jet helper. Their constitutive series,
two-positive-mass tail bounds, exact reference moments, canonical-versus-cache
point correction and normalization formulas are unchanged. Source-level AST
comparison verifies the unchanged functions. The changed Functional and output
functions only add per-instance moment-key diagnostics. They do not substitute
old local point results. The standalone old execution paths are removed.

The dedicated cancellation fixture uses128bits and exact
e0=1/5001+2^(-1000), with N=1. The constant coefficient of the truncated first
volumetric derivative is exactly5001·2^(-1000). Its rounded Arb enclosure may
cross zero, but the coefficient remains in the safe polynomial and encloses
the exact value. Additional addition/product guards retain uncertain
zero-crossing coefficients. Bulk and surface libraries share the same safe
target/helper objects. The prior independent992-check safe-kernel receipt is
bound to the current exact module hashes. This is substantially stronger than
merely checking file integrity.

For bulk, the reviewer independently reconstructs the cubic determinant jet
and shifted/axis-centered coefficient majorant on all144 actual elements.
The independently reused jet implementation is explicitly rebound to the safe
target before execution. Its input and the complete coefficient-selector
calculation use exact Fractions. All144 smallest permitted degrees, local free
maps, rho values and declared exact force/Hessian tail fractions reproduce.
The selector ranges from degree4 through24 and the actual chosen degrees
remain5–16. These are bounds on the omitted constitutive tails, not the whole
integration discrepancy. Sparse scatter tests cover unsorted local indices,
shared nodes, signs, symmetric row sums and signed R·F reduction.

For surfaces, an independent rational Bernstein calculation proves the supplied
positive-hoop and G bounds over all24 complete edge covers:222 stored leaves,
666 inequality subinterval proofs. With q0=(Glo+Ghi)/2 and
rho=(Ghi−Glo)/(Ghi+Glo)<1, the chosen analytic square-root branch and the
relative-binomial series/tails apply to each entire edge. The retained actual
policy is degree28 at512bits. Unsorted scatter and strict literal-True receipt
guards also pass.

The pointwise decomposition remains
FI−Fcache=(FI−Fstar)+(Fstar−Fcache), with the analogous Hessian identity.
The first term contracts signed exact-integral-minus-actual-rule moments of
the truncated constitutive polynomial, then adds both positive-mass remainder
bounds. The second evaluates the canonical and cached strain maps at the same
exact stored field and parameters. No nominal Gaussian exactness is assumed;
no low-degree moment defect is silently set to zero. None of this is a
parameter/tube-uniform integration-error bound.

Both protocols bind the immutable original B cache/proof center and same
physical parameters/rules/normalization. The source closure includes the safe
target, helper, moments, dense constitutive core, point library and wrapper.
No old element/edge result is bound as a reused numerical contribution.

This is launch permission for fresh144-bulk and24-surface POINT computations.
No actual field integration was run by this reviewer. The old negative
center-escape result and old point certificates remain withdrawn. Fresh
independent local replays, complete signed reductions and explicit point/root
scope review are still required before any scientific claim is restored.
