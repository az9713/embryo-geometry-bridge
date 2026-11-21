# Same-target physical-axis reduction, version2

The single version2 reduction completed with exit0 in58.07seconds. All347
checks pass. Its777 frozen input bindings specify exactly138 unchanged
original bulk packets, six separately reviewed physical-axis packets72–77,
and all24 unchanged surface packets. The first complete bound remains
unchanged and hash-bound; no packet was silently substituted into that result.

This is the same finite-P2 potential, same C, same common nodal box, source
dose, R, physical constants, cached quadrature and boundary conditions. The six
new integrals express the same full-axis-edge strain/basis polynomials in
physical(R,zeta) coordinates, with ordinary nonnegative powers. The exact
representation proof applies throughout the whole nodal box. It does not
claim polynomial divisibility on the six different axis-vertex-only cells.

The newly reduced bounds are

    b_row   =1.2959939358739153e−6,
    b_signed=5.6125894396497984e−8,

where b_row comes from exact rational |R| times absolute Hessian-discrepancy
rows, and b_signed comes from a separate256bit interval R·H product whose
row sums are computed from retained exact outward endpoints. Conservative
readable statements are b_row≤1.295993936e−6 and
b_signed≤5.612589440e−8. Either is a valid upper bound; the chosen minimum
does not discard a failed condition or assert a lower bound.

The first chosen bound was approximately0.01507198211105738. Its ratio to the
new chosen upper bound is about268539. This sharpens a conservative numerical
enclosure, not the underlying physics or the actual approximation error by
that factor. A lower-degree representation retains useful structure before
interval wrapping; the potential and admissible field family are unchanged.

The complete signed force-error box is also much narrower: its preconditioned
infinity-norm upper bound is approximately1.674936347176935e−11, versus the
first bound1.857898967426724e−6. Those are uniform error-box bounds, not
translated residuals or root-displacement bounds. Subtracting an actual
corrector must preserve signs and coordinate intervals.

The new chosen bound and the same fresh cached-Q bound give

    ZQ+b ≤0.10504048132121292 <1.

The root composer may instead use the fully exact-rational b_row bound, which
is already small enough to make b*s about7.37e−16 for s≈5.684e−10. This is a
uniform derivative prerequisite only. A fresh actual residual/parameter-path
bound, self-map inequality and corrector/domain inclusion must be checked
before any integral-root existence, uniqueness or location conclusion.

The reducer independently reconstructs all611 force sums,6676 upper-Hessian
slots,611 signed R·F coordinates and611 matrix row bounds. Its compressed
NDJSON retains every signed611-by611 R·H interval endpoint. New physical
packet maps, domain/tail policies, old-normalized cache data and safe point
consistency are checked. Every other packet is byte-identical to the first
selection. Generic reviewed local integration and representation algorithms
are explicit reuse; this stage performs no independent local reintegration
or mechanical solve. It does not establish PDE convergence, biological
validation, global parameter-domain inference, novelty or original-five-gate
completion.
