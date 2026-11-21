# Direct uniform-force crosscheck of the integral-root family

This is a separately executed exact-rational scalar composition. It shares the reviewed integration/component kernels and their enclosures with the principal proof. It is not an independent physical implementation. No integration, Hessian assembly, matrix product, optimization or equilibrium solve was repeated.

Write p0(theta)=u0+S(theta-theta0), p1(theta)=p0(theta)-d, and let R be the fixed nonsingular preconditioner. The old finite-quadrature result supplies ||R FQ(p0(theta),theta)||∞ ≤ etaQ for every theta in the same closed C. The new uniform component enclosure supplies Dbox containing R(FI-FQ)(u,theta) on the entire reviewed common nodal box and C. The independent endpoint reduction retains signs of all 611 preconditioned force coordinates. Use b from the fully exact rational |R|-times-Hessian-row majorant, rather than the tighter signed-Arb matrix bound.

Define m=max_i sup |Dbox_i-d_i|, z=zQ+b, and dn=||d||∞. All bounds are exact rational numbers. The common radius is rc=r+dn. Thus both old/new radius-r cubes and the complete line between their centers are in the common domain.

The fundamental theorem of calculus gives

    R FI(p1)=R FQ(p0)+[R(FI-FQ)(p0)-d]
                + integral_0^1 [I-R HI(p0-t*d)] d dt.

Consequently eta_direct=etaQ+m+z*dn bounds ||R FI(p1)||∞ uniformly in theta. On the whole new radius-r cube, T(u)=u-R FI(u) has derivative norm at most z. Strict eta_direct+z*r<r and z<1 prove existence and uniqueness in that cube by Banach, with error ≤eta_direct/(1-z). R is nonsingular by the already verified center defect ||I-R HQ||<1. Thus a fixed point is an FI zero. The geometry and analytic-domain certificates ensure differentiability throughout these cubes; they do not assert a continuum PDE solution.

Actual scalar bounds (decimal displays only): m=1.284615006e-14, z=0.1050417212, eta_direct=2.299666996e-12. The self-map left side is below 3.350084208e-12, strictly less than r=1e-11. The resulting error is below 2.569580114e-12. Acceptance used the complete rational fractions in results.json, not these displays.

There is also at most one FI zero anywhere in the convex common cube: two zeros are fixed points and their distance would be at most z times itself. This uniqueness argument does not require a self-map of the common cube. Coordinate588 satisfies |d588|-error>r, so the newly certified root lies outside the old closed cube. Since the old cube is inside the common cube, there is no exact-integral FI zero in the old cube. This is a new consequence of existence plus common-domain uniqueness plus location, not a consequence of the earlier map-center escape alone. The known old Q root has error at most etaQ/(1-zQ) about p0; the distance between the Q and integral roots is at least max(0,|d588|-error-etaQ/(1-zQ)). Exact values are retained in location_consequences.json.

The result concerns one record B, the declared tiny five-parameter C, fixed geometry/mesh and axisymmetric finite-P2 energy with exact spatial integration. It establishes no PDE convergence, other-record integral roots, global parameter-domain cover, dynamic path, biological observation alignment, energetic stability, or novelty. The historical unsafe point receipts remain withdrawn. This route uses corrected safe v2 uniform packets only. The arbitrary candidate d is read from exact differences of the two stored binary vectors; its historical proposal recipe is not a validated residual premise here.

Nine synthetic checks preceded the frozen protocol and twelve actual binding/arithmetic checks passed. One initially mistaken synthetic subtraction example was corrected before the protocol: .3 and .2 lie within the exact-subtraction range of Sterbenz's lemma. The original source and failure are preserved in engineering_revision_0. There was no failed or repeated physical calculation.
