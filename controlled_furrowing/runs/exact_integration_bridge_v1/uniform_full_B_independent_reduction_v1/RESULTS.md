# First complete uniform integration-error bound

The first combined reduction terminated successfully (session33739, exit0)
in84.91seconds. All348 review checks pass. This is an independently reduced
matrix-error component; no local moment integration or mechanical solve was
rerun by this reviewer.

The frozen producer records cover all144 bulk elements and24 surfaces on the
same explicit common B nodal box, the same five-coordinate C, dose, geometry,
normalization and physical parameters. Seven completed safe-pilot bulk packets
are reused explicitly by hash/selector, while137 bulk elements are newly
computed. No physical-axis replacement was silently inserted into this first
result. Every component's free-coordinate map is independently rebuilt from
the original connectivity and frozen free-index map.

The exact reduction retains all611 signed force endpoints and6676 symmetric
upper-Hessian endpoint sums. All local metadata, runtime positive-domain
checks and derivative-tail policies are checked against the frozen selectors
and prelaunch domain proofs. The original safe point u0/theta0 belongs to the
common box/C in all650 coordinates. The new safe point enclosures overlap all
local uniform intervals (and in this dataset are all contained); that is a
consistency check, not the proof of uniformity.

Let Hbox enclose HI−HQ entrywise for every field/parameter in this supplied
family. The exact symmetric absolute row bounds are
h_j=Σk sup|Hbox_jk|. Treating each original binary R entry as its exact rational
value gives the independently computed bound

    b_row=max_i Σj |Rij| h_j ≈0.01507322434861234.

A separate256bit Arb matrix product uses the independently assembled Hbox,
and its611 row sums are then reduced from exact outward endpoints. The full
611-by611 signed matrix endpoint output is retained as compressed NDJSON.
That bound is

    b_signed≈0.01507198211105738.

Both are valid upper bounds. Their minimum remains a valid upper bound, so a
readable outward statement is

    sup ||R(HI−HQ)||∞ ≤0.015071982112.

No float comparison is used to accept the norm inequality. The exact endpoints
and sums are in exact_reduction.json. The synthetic tests distinguish exact
submultiplicativity from possible small extra outward Arb rounding at equality;
neither rounded bound is assumed always to dominate the other numerically.

On the same common affine/free-coordinate tube, the separately computed
cached-Q derivative bound is ZQ≈0.10504042519531853. The triangle inequality
therefore gives a derivative bound for the integral map of approximately

    ZQ+b ≤0.1201124073063759 <1.

This is not a self-map or root-existence theorem. A residual bound at an
appropriate center must still fit its radius, uniformly over C. In particular,
the coordinator's coarse translated-residual estimate is about1.11e−11 and
can exceed the old candidate radius1e−11. The first valid derivative bound must
be retained even when that separate sufficient test is inconclusive. No
equilibrium nonexistence follows. A new same-target coordinate representation
may tighten local enclosures only through a separately reviewed adapter and
frozen packet-selection/reduction revision.

The direct uniform force box is also retained; it is much broader than a
derivative-controlled difference from the known point and is not a useful
replacement for the needed translated/predictor residual analysis. One must
not silently substitute the old point discrepancy for that uniform quantity.

Reused premises are the mathematically reviewed generic safe integration and
cache-consistency engines, the explicit whole-body and surface analytic-domain
proofs, and their frozen source/prelaunch reviews. This reduction does not
independently reintegrate all168 local moments. It supplies no integral root,
PDE convergence, biological validation, source-observation alignment, global
physical-domain inference, novelty or original-five-gate completion.
