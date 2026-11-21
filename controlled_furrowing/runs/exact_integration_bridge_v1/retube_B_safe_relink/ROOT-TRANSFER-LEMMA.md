# Conditional transfer to the exactly integrated finite-element energy

This lemma specifies the remaining proof obligations. It is not a statement
that their numerical inequalities have passed. The target is finite
dimensional, not the continuum elasticity PDE.

Let C be the unchanged closed material-parameter cell, and let
p(theta)=u0+S(theta-theta0) be the exact affine predictor determined by the
stored binary coefficients. Write FQ and FI for the gradients of the cached
and exactly integrated energies. R is the same square exact-binary
preconditioner. Set D(u,theta)=R(FI-FQ)(u,theta).

Use an explicitly convex nodal box X containing every p(theta), the base u0,
the candidate cubes p(theta)-d+B_infinity(r), and the swept cubes joining them
to the old predictors. Parameters and nodal fields are independent inputs in
X times C. This is stronger than bounding only the coupled family
p(theta)+B_infinity(rc). It is needed for the fixed-theta segment from u0 to
p(theta) in the variation argument below.

Assume all of the following have been proved on their stated domains:

1. At the old predictors, ||R FQ(p(theta),theta)|| <= etaQ uniformly on C.
2. On the swept cubes, ||I-R HQ|| <= zQ, with zQ<1.
3. On X times C, ||R(HI-HQ)|| <= b.
4. At the base point, ||D(u0,theta0)-d|| <= m.
5. ||D(p(theta),theta)-D(u0,theta0)|| <= V on C.
6. Both energies are C2 on these domains; canonical geometry and cached
   analytic-domain guards have been established separately.

Here all norms are induced infinity norms or the corresponding vector norm.
The field u1=p(theta)-d need not have been obtained by a nonlinear solve.

## Residual at the translated predictor

At a fixed theta the segment p(theta)-t*d, 0<=t<=1, lies in the swept domain.
The fundamental theorem of calculus gives

    R FI(p-d,theta)
    = R FQ(p,theta) + D(p,theta) - integral_0^1 R HI(p-td,theta)d dt
    = R FQ(p,theta) + [D(u0,theta0)-d]
      + [D(p,theta)-D(u0,theta0)]
      + integral_0^1 [I-R HI(p-td,theta)]d dt.

Consequently a valid translated-predictor residual upper bound is

    etaI = etaQ + m + V + (zQ+b)*||d||.

No cancellation between independently enclosed uncertain quantities is
assumed. Replacing V by zero without a separate theorem is invalid.

## A sufficient bound for V

The energy is affine in the five material/activation coefficients at fixed u.
Let aj bound ||R partial_theta_j(FI-FQ)(u0)||. Let Delta_j be the largest exact
distance of coordinate theta_j in C from theta0_j, and let
s bound sup_theta ||S(theta-theta0)||. The fixed-theta segment from u0 to
p(theta) stays in convex X, so

    V <= b*s + sum_j aj*Delta_j.

The partial derivatives are evaluated at the same base field, geometry,
normalization and physical force law. They are not derivatives through an
equilibrium solve. The activation coefficient is common gain times the exact
source dose. Force-scale symmetry has no role in this numerical bound.

## Fixed point and mechanical equilibrium

For each theta define TI(u)=u-R FI(u,theta). On the candidate cube,

    ||D_u TI|| <= z=zQ+b.

If z<1 and etaI+z*r<r, the map sends the closed candidate cube into itself and
is a contraction. Banach's theorem gives a unique fixed point there. Moreover
zQ<1 at a cached base point makes R HQ invertible by the Neumann-series
criterion. Since the matrices are square, R is nonsingular. Thus a fixed
point has FI=0, not merely R FI=0. The error from the translated predictor is
at most etaI/(1-z). Uniqueness is local to this cube.

## Conditional transfer of strict local stability

A small residual alone supplies no positive-definiteness theorem. However,
an independently established positive-definite cached base Hessian H0, at
exactly the same u0/theta0, provides an economical additional premise here.
Assume H0 is symmetric positive definite, and let z0 bound ||I-R H0||<1.
For any symmetric target Hessian H with ||I-R H||<=z<1, set
Ht=(1-t)H0+tH. Then

    ||I-R Ht|| <= (1-t)z0+t*z < 1,  0<=t<=1.

Every Ht is nonsingular, because R Ht is invertible. Eigenvalues of symmetric
matrices vary continuously; none can change sign along this path without
passing through zero. Thus H is positive definite. This argument needs no
assumption that R itself is symmetric. It must bind an actual independently
proved SPD H0, not merely a numerically positive sampled spectrum.

Moreover H^(-1)=(R H)^(-1)R and the induced-norm Neumann bound imply

    ||H^(-1)||_infinity <= ||R||_infinity/(1-z).

Since the positive-definite inverse is symmetric, its spectral radius is at
most its infinity norm. Consequently

    lambda_min(H) >= (1-z)/||R||_infinity > 0.

If the uniform integral Hessian bound above passes and the old base-SPD
premise is correctly linked, this certifies uniform energy curvature and
strict local minimality of the new integral equilibrium. No such numerical
composition is asserted by this preparation note before its premises pass.
It proves no common loading path, global
minimum, other-branch exclusion, finite-element/PDE convergence, experimental
fit, full physical identified set or originality.

## Current numerical premises

## A location consequence requiring the NEW existence proof

Suppose the new existence/error test passes with error eI from p(theta)-d.
For a fixed free coordinate j,

    |uI_j-p_j(theta)| >= |d_j|-eI.

If this lower bound exceeds the old radius r, the new integral root is outside
the closed old predictor cube. This follows from the new root/error theorem,
not from the earlier point map-center escape.

There is also an exclusion statement with additional domain bookkeeping. The
derivative bound zQ+b<1 holds on the entire common convex cube at fixed theta.
Two integral zeros there would be fixed points of TI and obey
||x-y||<=z*||x-y||, so must coincide. This at-most-one-root argument does not
require TI to map the common cube into itself. The candidate cube belongs to
the common cube because ||d||+r<=rc. Thus the newly proved root is its unique
integral root; if that root is outside the old cube, there is no integral
root in the old cube. All constants are uniform on C. No roots outside the
common cube or parameters outside C are excluded.

For the already known cached root uQ in the old cube, the fresh common-cube
Q derivative bound gives ||uQ-p(theta)||<=etaQ/(1-zQ)=eQ. Therefore

    |uI_j-uQ_j| >= max(0, |d_j|-eI-eQ).

These are conditional consequences to evaluate after the new exact tests,
not location conclusions from the historical unsafe or corrected point
certificate alone.

## Current numerical premises (still conditional until composition)

The new independently reviewed safe point calculation can supply premise 4
after exact subtraction of the actual unchanged candidate decrement; see
results.json in this directory. The fresh cached-Q common-tube calculation
supplies premise 2 within its precisely declared coupled swept family; see
the separate cached_Q_common_B_v1 proof and endpoint/domain audit. Premises
3 and 5 still require completed, reviewed uniform component calculations and
the parameter-derivative bridge. No integral-root conclusion is asserted by
this preparation note.
