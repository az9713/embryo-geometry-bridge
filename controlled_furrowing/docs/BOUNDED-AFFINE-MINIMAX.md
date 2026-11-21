# Exact certificates for a local multivariate calibration proposal

**The bounded affine minimax subproblem can be certified exactly even when its coefficients originated in approximate mechanical sensitivities.** These are different claims: an exact LP optimum does not certify a nonlinear mechanical optimum or the accuracy of a derivative. The [implementation](../scripts/bounded_affine_minimax.py), [32-check synthetic receipt](../evidence/bounded_affine_minimax_tests.json) and [independent review](BOUNDED-AFFINE-MINIMAX-REVIEW.md) retain that distinction.

## Defined local problem

For residual vector \(r\), response matrix \(D\), and a declared numerical step box \(\ell\le h\le u\), solve

\[
 \min_{h,t}t\quad\text{such that}\quad
 -t\le r_i+D_i h\le t\quad\text{for all }i,\qquad
 \ell\le h\le u,\quad t\ge0.
\]

All 68 endpoint constraints can enter this problem. A step bound limits one local proposal; it does not restrict the scientific parameter domain or certify that the linear approximation holds throughout the box. The actual material inference must still evaluate the proposed material and gain through nonlinear mechanics, preserve its physical constraints, and propagate the remaining compatible set.

Coefficient conversion is explicit. The helper preserves supplied exact `Fraction` values or rational strings. Floating inputs alone become the exact decimal values of `str(float(x))`. To certify an existing binary-rounded surrogate without changing its definition, callers pass `Fraction.from_float(x)` rather than silently converting its coefficients to a different convention.

## Proof by primal and dual witnesses

Collect the signed-residual, box and nonnegative-loss inequalities into

\[
 Ax\le b,\qquad x=(h,t),\qquad c=(0,\ldots,0,1).
\]

For any \(\lambda\ge0\) with \(A^T\lambda=-c\), multiplication of the inequalities gives

\[
 -c^Tx=\lambda^TAx\le\lambda^Tb,
 \qquad\text{hence}\quad c^Tx\ge-\lambda^Tb.
\]

If a feasible \(x_*\) also satisfies \(c^Tx_*=-\lambda^Tb\), it attains that lower bound and is globally optimal for this defined LP. This argument is elementary linear-programming duality, not a new optimization method.

HiGHS supplies floating candidate active constraints and multipliers. A separate rational elimination reconstructs a proposed primal vertex and dual support. The code then checks **every** primal inequality, nonnegativity of **every** multiplier, exact stationarity and equality of the two objective values. Unsuccessful basis recovery returns UNRESOLVED rather than accepting the floating optimizer. The finite basis search is a computational recovery budget, not a mathematical coverage assertion; a complete valid certificate needs only one verified matching primal/dual pair.

The independent certificate verifier takes the original coefficients and the stored exact witnesses, without trusting a floating optimal status. No tolerance appears in its final optimality predicates. Floating thresholds are used only to propose which exact systems to try.

## Checks and limits

Analytically solvable cases cover interior balancing, positive and negative boundary solutions, a one-sided tension step, an exact-zero fit, an unused coordinate, independent channels and a fixed coordinate box. Adversarial tests lower each certified objective and verify rejection. Invalid dimensions, reversed bounds, malformed string vectors and nonfinite coefficients are rejected; exact coefficients outside floating proposal range return UNRESOLVED. Independent review found the malformed-string and overflow handling gaps before this hardening. The exact pre-hardening source is archived at `runs/bounded_affine_minimax_review/initial_revision/bounded_affine_minimax.py`, with a hash matching the original cohort extraction protocol. The valid earlier rational certificates remain retained; a versioned re-extraction uses the hardened helper. These checks establish the affine component's behavior; the independent review supplies a separate assessment.

For an actual cohort, the mechanical equation \(R(u,\theta)=0\) suggests

\[
 R_u\,u_{\theta_j}=-R_{\theta_j},\qquad
 O_{\theta_j}=O_u u_{\theta_j}+O_{\theta_j}^{\rm explicit}.
\]

At a regular exact root these identities follow from the implicit-function theorem. The current tangent calculations use floating saved fields and inverses. Spatial maxima can be nonsmooth, and two finite-difference step sizes do not prove differentiability or a uniform approximation bound. Certifying the LP therefore does not make its coefficients exact physical derivatives, establish root existence, bound continuum error, validate a temporal observation contract or compute the full identified set. See [material sensitivities](MATERIAL-TANGENT-INFERENCE.md), [physical inference prerequisites](PHYSICAL-INFERENCE-READINESS.md) and the [original five-gate plan](../PLAN.md).
