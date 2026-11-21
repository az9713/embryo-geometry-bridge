# L surface force enclosure on the original box

## First reduction

Retain the 24 reviewed uniform surface force and Hessian packets for the first centered bulk reduction. Exact signed assembly of those 24 saved force vectors followed by the original binary preconditioner gives a surface-only upper bound of **2.9193098488439825e−15**. The strict original-cube force allowance at the latest reviewed derivative bound is **7.021720358707882e−12**. The entire retained surface contribution uses **0.0415754%** of that allowance. The surface-only absolute Hessian row projection `|| |R| h_surface ||∞` is at most **6.31367211138398e−10**. All 611 surface-only preconditioned force intervals contain zero.

This is an exact arithmetic check on saved enclosures. Preserve their signed endpoints when assembling the 144 centered bulk terms and 24 retained surfaces. Compare the resulting full preconditioned bound and the strict self-map test on the same L cube. This surface budget calculation does not establish that the centered bulk force fits or that an exact-integral root exists.

## Smallest centered surface implementation if needed

For one edge, let `D(Y, θ)` be its signed integral-minus-cache force in normalized free nodal coordinates, and let `B` be the reviewed original nodal box. With `Y0` the exact cached proof-center field, use

```text
D(Y, θ) = D(Y0, θ) + integral_0^1 H_D(Y0+t(Y−Y0), θ) (Y−Y0) dt.
```

1. Construct `Y0 = context.target.nodes()` and verify it equals the frozen unshifted coordinate centers exactly. Use a separate in-memory singleton box with the same material reference, triangles, and center and zero radii. The state increment uses only the selected local free columns, with exact signed differences or the original symmetric radii. Constrained coordinates have zero increment. Keep all six current coordinates when evaluating the edge fields and cache maps.
2. For each of the 24 edge cases, form `target.surface(face, edge, nodes=Y0)`. Build a fresh singleton `SUR.uniform_G_bounds` and 12 `SUR.local_site_domain` guards. Require a positive hoop, positive squared surface Jacobian, and convergent relative binomial radius on the full closed regularized edge. The cache guards need a positive hoop and nondegenerate tangent over the same-node cache-to-canonical strain segments. Certify a force-tail budget at the new expansion center and order. If reusing the original expansion center/radius, independently prove the singleton polynomial and its whole-edge values fit those bounds.
3. Call `exact_cohort_integrator.surface_component(ctx, Y0, rule, point_case)`. Pass the original record row and full closed parameter cell. This evaluates `D(Y0, θ)` uniformly over the same material intervals and keeps all context/rule checks. A point field can cancel polynomial terms exactly, so prepare fresh coefficient support. `verify_frozen_domain` requires matching supports when handed an old case. The returned Hessian at `Y0` is incidental.
4. Apply the already reviewed signed edge Hessian interval over `B × C` to the interval of free displacements `Y−Y0`. The state paths are contained in the convex original box, and the same parameter is held fixed along each path. Intersect the resulting signed force interval with the retained direct whole-box interval. The original signed Hessian packet remains the derivative premise. Assemble by the original global free map before preconditioning.

This ordering avoids a separate material derivative or mixed parameter-state remainder. Using full-cell component intervals at the singleton may lose cancellation through the repeated coefficients `K` and `T−K`. The affine method below preserves those parameter dependencies and can be tried if the singleton enclosure is still wide.

## Explicit affine material derivatives

For face 0 use the basal pair `K_b,T_b` and for face 1 use the apical pair `K_a,T_a`. The surface energy has zero derivative with respect to the common active gain. Its normalized strains in physical reference radius are `a=r/R, b=dr/dR, c=dz/dR`; at the axis the exact zero radial trace removes the hoop denominator symbolically. Let `J_s=a sqrt(b²+c²)>0`, `L=J_s²` and let `P` and `Q_s` be the full force defects for energies `L/2` and `sqrt(L)`. Then

```text
W = (K/2)L + (T−K)sqrt(L) + K/2
D(Y, K, T) = K P(Y) + (T−K) Q_s(Y)
∂K D = P − Q_s
∂T D = Q_s
D(Y, K, T) − D(Y, K0, T0) = δK P + (δT−δK) Q_s.
```

The parameter contribution for a state-first path needs `P,Q_s` on the whole nodal box, but they are independent of material parameters and the difference identity is exact at fixed Y. Holding θ variable in the singleton instead restricts coefficient evaluations to `P(Y0),Q_s(Y0)` and uses the retained full-cell Hessian for the state step. A point evaluation at only `K0,T0` without one of these parameter contributions omits parts of C.

For local force column i, (P_i) is the signed moment defect of `D_i L / 2` plus its finite cache correction. (Q_{s,i}) uses the signed moment defect of `D_i L/(2 sqrt(L))` plus its finite correction. `D_i L` comes from `JET.surface_squared_jacobian_jet`. `SurfaceFunctional.defect` subtracts Qstar at exact inferred sites from the canonical radial integral with original (2π/σ) normalization. The finite correction is Qstar minus the frozen binary rule. `exact_surface_cache_consistency.actual_maps` returns the exact 3×6 canonical and normalized cached strain maps, with the shared free-column map. Both are applied to the same six coordinates. The cached normalization is `diag(L_ref/(r0c L0c), 1/L_ref, 1/L_ref)`. Apply each positive old normalized binary weight once. For the two affine basis forces, explicit pairs `(K,T)=(1,1)` and `(0,1)` in the cache kernel yield (P) and (Q_s) respectively. The public physical-cell guard requires positive K, so call the explicit kernel for the zero-K derivative basis. The strain positivity checks still apply.

For a first-derivative truncation of degree `N`, (P) is exact polynomial/Laurent. With (q_0>0), (|(L-q_0)/q_0|≤ρ<1), and (μ_+=μ_I+μ_*), use per-column

```text
tail_Qs[i] ≤ (μ_+ / (2 sqrt(q0))) binomial_tail(−1/2,N,ρ) sup|D_i L|.
μ_I = (2π/σ) (Rhi²−Rlo²)/2, μ_* = sum of old normalized weights.
```

For base-pair force, multiply by `|T0−K0|`. For a direct parameter-box force, use an outward upper bound on `|T−K|`. The exponent `−3/2` is relevant to a fresh Hessian series and is unnecessary for a force-only enclosure using retained Hessian intervals. Do not truncate `K/2` or add another normalization to finite corrections.

The previous B helper `surface_poly_defect` is explicitly point-only. Its `split_surface` recovers `Q_s` by dividing by `T−K`. Explicit basis evaluation avoids that division when the gap is small or vanishes. Use explicit target, nodes, face, edge, order, domain, sites, and material coefficients in new code. Leave the frozen B and element-21 launch globals untouched.

## Source anchors

- `scripts/exact_fe_target_safe.py:148` builds the physical-radius P2 strains and free derivatives, with exact axis regularization.
- `scripts/exact_cohort_integrator.py:294` checks the frozen surface rule, and line 305 forwards point or interval nodes to the reviewed uniform kernel.
- `scripts/exact_uniform_surface_bridge_B_safe.py:136,186,205,219,238,271,304,330` provides Bernstein whole-edge positivity, domain checking, norms, binomial coefficients, face mapping, site guards, finite correction and edge assembly.
- `scripts/exact_surface_bridge_safe.py:77,116,130,150` provides binomial tails, signed surface moments and the point-specific reference implementation.
- `scripts/exact_surface_cache_consistency.py:75,86,98,110` supplies derivatives and exact same-node cache maps.
- `scripts/exact_parameter_defect_B_safe.py:43,78,82,94` records the point-only affine decomposition and polynomial cache recipe.

The accompanying receipt binds the source bytes, retained surface packet manifest and exact arithmetic. It contains no new integration, cache sum, or nonlinear solve.
