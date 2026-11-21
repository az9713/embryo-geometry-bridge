# Exact finite-element target preparation

This is a frozen preparation checkpoint, not a validated integration or a proof that the old finite-quadrature roots solve the new exact-integral equations. No equilibrium solve or full force/Hessian integration was performed.

The implementation is `../../../scripts/exact_fe_target.py`. `protocol.json` freezes two sources and 147 input bindings for all 17 independently reviewed selected root families, including the recentered F/K proof centers. The original physical fields, caches, protocols and scientific documents are unchanged.

## API

- Sparse polynomials use `{(i,j): coefficient}`. `add`, `sub`, `scale`, `mul`, `power`, `diff`, `shift`, `compose`, and `evaluate` support exact rational coefficients and Arb coefficients through explicit exact coercion. Integer Laurent powers are allowed; negative powers of nonmonomials are rejected.
- `Target(point_cache_packet)` decodes and validates the frozen binary cache. `nodes()` returns **exact binary reference coordinates plus injected free displacements**, including fixed coordinates without alteration. `node_tube(S, cell, base, radius)` encloses the exact affine predictor and full normal corrector cube. This describes the old certified root tube; it does not certify a root of the new integral target.
- `Target.element(index, nodes=None)` exposes `q`, `B`, `global_indices`, `active`, `z_center`, `detDM`, and `chart`. Off-axis variables are `(R,zeta)` with `zeta=Z-z_center`. Axis variables are `(s,t)` after exact cancellation of radial factors. `chart.kind` is `off_axis`, `axis_edge`, or `axis_vertex`; axis metadata include `R1` and the nonnegative Duffy Jacobian. `measure_polynomial` includes `R*detDM*DuffyJac` on axis, and is `R` in physical coordinates off axis.
- `Target.surface(k,e,nodes=None)` exposes canonical horizontal `Rref(s)=Rlo+L*s`, `L`, `q_numerators`, `B_numerators`, and `hoop_denominator`. The first axis-edge hoop factor is cancelled symbolically. `q_physical_R` and `B_physical_R` provide the same functions as Laurent polynomials in physical R for one-dimensional moments. `evaluate_surface` applies the explicit denominator.
- `Target.normalization()` encloses the **mathematical** scalar `2*pi/S_binary` with Arb; it does not replace pi with a floating constant.
- `Target.auxiliary_rule()` derives exact Q-star barycentric coordinates from the frozen shape coefficients and surface coordinates from the frozen edge coefficients. All points are strictly interior. Its weights are the old exact normalized binary weights; they are not asserted Gaussian-exact or equal to the canonical material measure.
- `core.bulk_local` and `core.surface_local` may evaluate the resulting local strain values after independently checking their analytic domains. `q`/`B` remain available before evaluation for later polynomial moment and remainder calculations.

## Receipts and scope

`guard_tests.json` and `preflight.json` contain 56 true checks: exact partition and coordinate reproduction, canonical identity and isochoric deformation, exact bulk/surface nodal derivatives, independent symbolic energy directional checks, axis endpoint regularity, interval-coefficient operation, strict source/capsule bindings, and actual selected-cache preparation. The selected current fields retain two fixed reference-midpoint Z offsets: DOF587 is -2^-56 and DOF599 is +2^-55.

`cache_interpolation_defects.json` reports exact discrepancies across all 384 cached volume shape coefficients and all 110,592 physical-gradient coefficients at 9,216 volume rule sites. These are interpolation discrepancies, not integrated force bounds. Surface numerator/denominator representations and auxiliary coordinates are supplied, but no full surface correction or energy integration is claimed at this checkpoint. Three actual local center samples passed the local analytic-domain check; that sample test is not a whole-body proof.

A future bridge must bound both the exact-integral-minus-Q-star moment defect and the finite Q-star-minus-old-cache correction. It must also validate nonlinear remainders over the full tubes. None of those results follows merely from this preparation.

`engineering_revision_0` preserves a mistaken synthetic example (whose corner limit was actually unique). `engineering_revision_1` preserves the initial Fraction/Arb operator incompatibility caught by the interval-coefficient synthetic test. Both were corrected before the frozen protocol and did not affect old scientific artifacts.

No author message was sent. The original full research goal, continuum/PDE accuracy, source-time validation, complete physical parameter inference and confirmed novelty remain unresolved.
