# Exact endpoint refinement of the point-force tail

For each point-force coordinate, the raw independent replay computes a polynomial moment interval and a signed same-point finite-cache interval. Call the sum of these intervals (X), where the finite-cache interval has the same exact nodal vector in both routes. The omitted polynomial remainder lies in ([-t,t]) for the independently bounded scalar (t). Bulk uses an exact rational scalar tail. Surface uses the exact upper endpoint of an outward Arb enclosure for its nonnegative scalar tail. The same material and auxiliary measures, scalar series, coefficient norms, frozen occupancy, and exact activation are used throughout.

The raw code manufactures an Arb ball (S) around the exact zero using the magnitude conversion of the previously computed scalar tail. The radius field of that ball represents a nonnegative exact dyadic number (s) once converted to a rational. We require that field to be exact, that the zero midpoint be exact, that the two returned endpoints be (-s,+s), and that (tle s). The companion reconstructs that very scalar input at 640 bits from the admitted exact point and frozen plan. For the bulk tail this uses the exact-rational coefficient norm and the shared geometry mass bound. For the surface tail this uses the same symbolic rational zeroth measure moment, derivative polynomial norm, binomial tail, and rounded scalar operations as the raw replay. No higher moments or finite-site force sums are reevaluated.

Write the raw serialized enclosure as (J=[ell,h]). The bulk grouping is a previously calculated interval sum followed by addition of (S). The surface grouping adds (S) to the polynomial moment interval before adding the cache interval. Outward interval addition covers the complete Minkowski sum in either grouping. Thus (X+[-s,s]subseteq J) even when the point components are interval-valued. This is stronger than having (J) enclose an unspecified physical value. It implies (Xsubseteq[ell+s,h-s]), since for every (xin X), both (x-s) and (x+s) lie in (J). Require (ell+sle h-s). Adding the mathematically sufficient tail exactly gives the valid enclosure

[
[ell+s-t, h-s+t].
]

Compare its exact rational endpoints with the original unchanged producer interval. Its signed inclusion must hold in both directions. No overlap, arbitrary tolerance, intersection of competing intervals, or alteration of the producer packet meets this check. Persist the raw endpoints, exact (s,t), refinements, and any remaining misses as separate evidence.

The companion requires a truly completed 168-packet independent replay and its observed OS exit. The raw replay must already contain all whole-box Hessian and direct point cache terms. Neither is narrowed here. A failed containment remains unresolved. All other input and proof scope stays unchanged. The cold companion runs only a synthetic scalar fixture.
