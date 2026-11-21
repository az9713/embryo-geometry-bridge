# Complete B-point bridge: rigorous negative self-map certificate

**PASS195checks.** The independent exact reduction terminated with exit0. It covers all144bulk elements and24surface edges, all611force entries,6676upper-Hessian slots and611signed preconditioned force entries. Main receipt: `evidence/exact_complete_point_B_independent_review.json`, SHA256bc22a23e7af96d0d1a99836401112378600c14be61d5cba957e6f15a12d699c5.

## The two equations and the common point

Let EQbe the old binary-cache quadrature potential and EIthe coherent exactly integrated P2-on-affine-reference potential. Both use the same611free coordinates u, exact stored current nodal coefficients X_binary+u, fixed coordinates, element mask, material coefficients, exact source-dose drive and positive frozen binary normalization. Let FQ=gradient(EQ), FI=gradient(EI), and Rbe the old exact binary preconditioner.

The independently reviewed old predictor is p_theta=u0+S(theta−theta0). Exact fraction checks show that the complete calculation's base parameters equal the original proof's theta0and lie in its **closed** common cell C. In particular the zero apical prestress is an allowed boundary value. Thus p_theta0=u0and the old uniform residual bound applies at the actual center used here. The exact common gain times dose is used; the archived rounded solver-strength value is not substituted. Cache, saved-vector, reference, free-map, preconditioner and normalization bindings agree.

The old producer's validated eta is approximately5.287455036108341×10^-13, while the independent192-bit replay gives5.287455043740370×10^-13. These are different valid outward majorants. This independent certificate uses their **maximum**, avoiding an unjustified equality of cross-precision bounds. The exact selected radius r=10^-11is contained in the independent radius enclosure; it is not falsely equated to that enclosure's lower endpoint.

## Complete error enclosure and proof

All168local packet hashes, policies, free maps and local Hessian coverage were checked. Exact rational endpoint sums reconstruct the global force and symmetric Hessian. Signed multiplication by the exact binary entries of R reconstructs every preconditioned error component. The producer's combined intervals contain these independent enclosures. Every surface domain is linked to the separate exact Bernstein proof, and all six independent640-bit surface replays overlap the current24-edge outputs with their immutable contexts verified. The earlier bulk assembly audit supplies full packet checks plus six distinct actual local numerical replays; this review does not pretend that all168local integrals were independently reintegrated.

Writing DeltaF=FI(u0,theta0)−FQ(u0,theta0), the complete error intervals give approximately

    1.6736920935997396e−11 ≤ ||R DeltaF||infinity
                            ≤1.6736920960160624e−11.

All authoritative endpoints are retained as exact fractions. Reverse triangle inequality and the old residual theorem give

    ||R FI(u0,theta0)||infinity
      ≥ ||R DeltaF||infinity−||R FQ(u0,theta0)||infinity
      ≥1.620817543162336e−11 >r=1e−11.

Therefore the **new integral map** TI(u)=u−R FI(u,theta0), using the old R, sends the old proof center outside the old radius-r cube. It cannot be a self-map of that cube. This rigorously rules out transferring the unchanged Banach self-map certificate to the new energy. The old quadrature map TQand its proofs remain valid for their original target.

An explicit witness strengthens assessability: free coordinate588is the axial coordinate of node313. The interval for TI(u0)−u0at that coordinate is approximately

    [−1.72656664645e−11, −1.62081754316e−11],

wholly below the allowed coordinate displacement interval[−1e−11,+1e−11]. `coordinate_witness.json` retains exact endpoints, the node/free-coordinate mapping and source hashes. Its scale is a numerical proof radius, not an observed biological deformation.

## What is and is not disproved

This does **not** prove that an exactly integrated equilibrium is absent, or even outside the old cube. A simple counterexample is T(x)=1.2−0.5x on[−1,1]: T(0)is outside, but its fixed point0.8is inside. This is also a gradient iteration for a strictly convex quadratic energy. Non-selfmapping of one chosen operator/domain is a different statement from equilibrium nonexistence.

No uniform error bounds a or b over the parameter/root tube have been established by this point calculation. The small point Hessian discrepancy is not automatically a uniform bound. A future uniform Lipschitz constant ZIwould yield the necessary condition L_center≤(1+ZI)rfor a root inside that cube; applying it would require actually proving that new uniform constant.

The meaningful next computation is a new numerical center/domain or preconditioner with validated uniform derivative errors, fresh geometry and observer transport. Source measurement alignment, full physical parameter coverage, continuum/PDE error, loading paths and biological validation remain separate original requirements. This negative certificate is a precise numerical-analysis result for the declared model, not a new biological mechanism or a completed project.
