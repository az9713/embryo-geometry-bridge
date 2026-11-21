# Safe independent point replay and complete signed reduction

This version preserves the historical implementations and replaces their
shared unsafe sparse-interval dependency. The six declared bulk elements are
0,6,72,66,48,143. Replays use degree18,640bits and the reviewed pointwise
|J−1|≤21/500 domain. They are not uniform in the parameter or nodal tube.

For the volume term, q=(a,b,c,d,h) is linear in the current nodal coefficients,
and J=adh−bch. The independent jet differentiates these two signed cubic
monomials directly. The producer first differentiates the meridional
determinant and multiplies h. The independent constitutive coefficients follow
from convolution of the log and reciprocal series, rather than calling the
producer coefficient generator. For μ=1,κ=5000, the coefficient of e^n in
φ′(1+e) is (−1)^(n+1)(1+5000 H_n), and the φ″ coefficient is(n+1) times the
next φ′ coefficient. H_0=0. The same independently reviewed tail formulas are
explicitly reused.

Let L=I−Q*, with I the normalized coherent physical element integral and Q*
the actual binary-node/binary-weight auxiliary rule evaluated on coherent P2
fields. For each monomial, L is evaluated with a signed exact reference moment
minus its actual rule moment. No nominal Gaussian exactness is assumed.
The reference moment includes the radial measure once. The signed polynomial
contraction happens before absolute values. Retained Arb coefficients must
not be discarded merely because they contain zero.

If t1,t2 bound the omitted φ′ and φ″ tails and M=I[1]+Q*[1]>0, the force-tail
bound is M t1 ||J_j||∞. The Hessian-tail bound is
M(t2 ||J_j||∞ ||J_k||∞ + t1 ||J_jk||∞).
Both positive masses are necessary: L of the remainder is bounded by the
sum of its integral and discrete-rule absolute bounds. A scalar-energy
remainder alone would not establish these derivative bounds.

At all64 auxiliary sites, the independent replay rebuilds canonical P2
barycentric values and physical derivatives, lifts both star and frozen-cache
strains from the same exact binary X+u nodal coefficients, and contracts the
constitutive force/Hessian entries with explicit sums. The reviewed local
constitutive primitives are reused; the producer's local matrix-lift routine
is not called. Every canonical strain/free-basis coefficient is checked
against the polynomial field. The complete point discrepancy adds the signed
finite Q*−Qcache correction to L and its tails. The independent direct
gradient0 and Hessian00 polynomial contractions also cross-check the shared
reviewed moment-dual accelerator.

The full reduction is a different verification layer. It reads all144 newly
computed bulk and24 newly computed surface local packets and adds every
serialized rational interval endpoint exactly. Upper-triangle local Hessian
slots are mapped to sorted global pairs; off-diagonal entries enter both row
sums once. The611-by611 preconditioner entries are interpreted as their exact
binary rational values. For each signed row, interval multiplication reverses
the endpoints when the preconditioner coefficient is negative. This produces
the complete611 signed R·δF intervals without replacing δF by its norm first.

Let L0 be the maximum of the lower absolute bounds of those611 intervals.
The independently valid old finite-quadrature certificate supplies
||R FQ(u0)||∞≤η at the exact same base parameters and u0. Use the larger of
the author and independent η bounds. Then
||R FI(u0)||∞≥max(0,L0−η).
Only if this exceeds the old exact radius r does the *new integral map*
u↦u−R FI(u) send the old center outside the old cube. That conclusion is
computed from the new safe results; it is not assumed from the withdrawn
version. It is not a root-nonexistence or root-location theorem. For example,
T(x)=6/5−x/2 maps0 outside[-1,1] but has the fixed point4/5 inside it.

The old quadrature-root and energy proofs address a different, unchanged
finite-sum target. The new point discrepancy does not yet provide a uniform
integration-error bound, a new integral root, PDE convergence, a time-dependent
observation operator, or a biological validation. Six of144 bulk elements and
six of24 surface edges are independently reintegrated; all remaining local
packets receive engine/domain/tail/hash/map and exact-reduction checks, not
individual independent reintegration. Historical unsafe point receipts remain
withdrawn regardless of whether a corrected conclusion ultimately agrees.
