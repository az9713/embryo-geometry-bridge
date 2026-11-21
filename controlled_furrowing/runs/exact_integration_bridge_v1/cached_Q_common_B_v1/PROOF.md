# Fresh cached-Q derivative bound on the common B tube

This is a prerequisite calculation for the unchanged finite-quadrature energy.
It is not an exact-spatial-integral error bound, an integral-root certificate,
or a continuum PDE claim. Preparation checks exact inputs and containment but
does not assemble the611-coordinate Hessian or run its contraction bound.

## Domain and exact containment

The original affine predictor is p0(θ)=u0+S(θ−θ0), with the stored binary S
interpreted exactly and θ in the same closed five-coordinate cell C. The
original normal radius is r=1e−11. Let u1 be the already stored arbitrary
candidate vector, and define d=u0−u1 by exact rational subtraction of its
binary coordinates. Its infinity norm is

    ||d||∞ =2412045/144115188075855872.

The requested common normal radius is exactly

    rc=r+||d||∞=188144378943289/7036874417766400000000000.

For every θ in C and every ||e||∞≤r, the old cube p0(θ)+e belongs to
p0(θ)+B∞(rc). So does the candidate cube p0(θ)−d+e. More generally,

    ||−td+e||∞≤t||d||∞+r≤rc,   0≤t≤1.

Thus the full swept family p0(θ)−td+e is contained, including the translation
paths e=0. No accuracy or Newton-step interpretation of d is needed. The
historical origin of that candidate supplies no trusted integral residual
premise here.

Preparation independently computes exact rational affine-predictor coordinate
hulls and all650 physical-coordinate intervals. It verifies3250 containments:
old cube, candidate cube, translation path, swept cube, and inclusion in the
same newly certified canonical geometry box. Fixed coordinates remain fixed.
The normalization, original u0, binary R, S, geometry, cache, source dose,
common cell and exact gain-times-dose drive are unchanged.

## Analytic cached-Q family

For fixed geometry, each cached local strain has an affine dependence on u.
The actual body and surface quadrature maps/weights remain their stored binary
values. The calculation never regenerates quadrature or invokes a mechanical
optimizer. The physical activation strength is exact common_gain × source_dose,
rather than the separately rounded solver-display strength.

The reused pure functions are parametric_root_cells.build_point, cell_context
and variation, together with dense local gradient/Hessian formulas in
algebraic_root_feasibility. The sparse interval-polynomial helper with the
withdrawn cleanup path is not imported or used. The build_point input resolves
only its predictor path to an absolute path in memory; the original module's
OUT and all frozen files remain unchanged.

For each bulk site the entire enclosing cached-strain box must satisfy
det(Fm)>0, hoop>0 and J>1e−6. For each surface site its cached radial coordinate,
tangent squared norm and scale must be strictly positive. These checks cover
the convex boxes, not merely the center values. They are distinct from the
canonical whole-body geometry proof: slightly different frozen strain maps
require their own analytic-domain checks. Parameter and predictor segments to
the base are enclosed because θ0 lies in the convex cell C and the predictor
is affine. Every normalized local weight must be strictly positive.

## Row-majorant proof

Let H0=HQ(u0,θ0), and let hT be a local constitutive Hessian with local-to-free
strain map BT and positive normalized quadrature weight wT. The fresh point
assembly gives a rigorous Z0≥||I−R H0||∞. On the common tube, local interval
evaluation bounds DT,ab≥|hT,ab−hT,ab(base)|. Put sT,b=Σk|BT,bk|. Then

    Σk |(BTᵀ (hT−hT0) BT)jk|
    ≤ Σa |BT,aj| Σb DT,ab sT,b.

Summing these positive bounds over sites and applying entrywise |R| gives
611 preconditioned row majorants vi. Consequently

    ||I−R HQ(u,θ)||∞ ≤ Z0+max_i vi = ZQ

for every θ in C and u=p0(θ)+e with ||e||∞≤rc. The retained per-row bound is
Z0+vi; it uses the global point bound Z0, not an uncomputed sharper individual
point-row defect. Final maxima and additions are taken on exact rational
outward endpoints from the Arb enclosures. Nothing is accepted using a float
tolerance.

If ZQ<1, the cached-Q Newton-like map has a uniform derivative contraction
bound on the common tube for each fixed θ. This does not supply a new self-map
bound or an integral-root proof. In particular the whole-tube derivative
defect between HI and HQ is a separate missing quantity; a pointwise H defect
cannot replace it.

## Execution gate and retained evidence

The source, inputs, common radius and containment ledger are frozen in
protocol.json. preflight.json records preparation-only checks. The run refuses
to assemble the611-dimensional point matrix without a separate nonempty
literal-True parent review binding the exact source, protocol and preflight
hashes. Only one versioned job may start. It retains the fresh point summary,
all9504 cached-domain bounds,611 variation and uniform row bounds, final ZQ,
PID/job/terminal records and any failure. A failed analytic or sufficient
condition does not imply that an equilibrium does not exist.

The first preparation attempt was stopped by the strict Boolean gate because
one NumPy finite-array check returned numpy.bool_. The versioned source/failure
record is preserved; explicit conversion to a Python Boolean changed no
arithmetic or physical input. No611-coordinate calculation occurred then.
