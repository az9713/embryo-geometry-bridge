# Independent mathematical review of the auxiliary scalar null

**Scope:** derivation and numerical-root verification of the frozen scalar response. This is not completion of the coupled geometric prediction gates. A calibrated maximum-shortening curve is a homogeneous maximum proxy; it does not validate equilibrium kinetics, predict individual surface positions or establish a biological mechanism.

## Derivation from the column energy

For the conditional volume-conserving column energy in MODEL-CANDIDATES,

\[
E(A;I)=\Sigma A+\frac{K_b}{2A_0}(A-A_0)^2
+2(\sigma_{l0}+gI)V\sqrt\pi A^{-1/2},
\]

equilibrium requires

\[
\Sigma+K_b(A/A_0-1)
-(\sigma_{l0}+gI)V\sqrt\pi A^{-3/2}=0.
\]

Let \(A_*\) be its unactivated equilibrium area, \(h=V/A_*\), \(x=A/A_*\), and introduce positive baseline effective tension and nonnegative stiffness

\[
L_0=\Sigma+K_b(A_*/A_0-1)
=\sigma_{l0}V\sqrt\pi A_*^{-3/2},
\quad B=K_bA_*/A_0,\quad Q=V\sqrt\pi A_*^{-3/2}.
\]

Subtract the unactivated equilibrium after multiplication by \(x^{3/2}\):

\[
L_0(x^{3/2}-1)+B(x-1)x^{3/2}=gQI.
\]

Therefore, with

\[
\nu=\frac B{L_0+B},\qquad a=\frac{gQ}{L_0+B},
\]

the implemented equation is exactly

\[
F_\nu(x)=(1-\nu)(x^{3/2}-1)+\nu(x-1)x^{3/2}=aI,
\qquad S=h(1-x^{-1}).
\]

If \(A_*=A_0\), then \(L_0=\Sigma\), \(B=K_b\), giving the simpler \(\nu=K_b/(\Sigma+K_b)\). In general, fitted dimensionless parameters absorb baseline prestretch; they do not uniquely recover a physical basal modulus or active force.

Positive baseline lateral tension implies \(L_0>0\) and hence \(0\le\nu<1\). The declared grid endpoint \(\nu=1\) is a limiting zero-baseline-tension case, not the interior of that strictly positive class. Its mathematical response is well defined and may be used as a clearly labeled limiting null/sensitivity case. No physical basal stiffness can be inferred merely because that endpoint gives the best scalar fit.

## Root existence, sign and bracket

For \(I,a\ge0\), seek \(x\ge1\). We have \(F_\nu(1)=0\), unbounded growth, and

\[
F'_\nu(x)=\sqrt x\left[\frac32(1-\nu)+\frac\nu2(5x-3)\right]
\ge\frac32-\frac\nu2\ge1.
\]

There is consequently exactly one root on this physical branch, with area expansion and shortening satisfying

\[
x\ge1,\qquad 0\le S<h.
\]

The implementation's upper bracket \(1+2aI\) is valid. A stronger elementary bound is

\[
F_\nu(x)\ge(1-\nu)(x-1)+\nu(x-1)=x-1,
\qquad 1\le x_*\le1+aI.
\]

Thus the looser implemented upper bound is conservative, including \(aI=0\), where both endpoints are exactly one. Bisection uses monotonicity to select the bracket half. In exact arithmetic, after \(n\) steps the interval has width \(2aI/2^n\). Since \(dS/dx=h/x^2\le h\), a root-interval width also bounds the range of shortening by \(h\) times that width.

These last statements are exact-arithmetic facts. Double-precision midpoints eventually stagnate and computed signs have rounding error; 60 iterations are **not** an interval-arithmetic certificate. The independent receipt compares production values against high-precision roots of a different polynomial representation,

\[
\nu y^5+(1-2\nu)y^3-(1-\nu)-aI=0,\qquad y=\sqrt x\ge1,
\]

and checks the closed form \(x=(1+aI)^{2/3}\) at \(\nu=0\). Agreement is numerical verification over the recorded test range, not proof of a uniform floating-point enclosure.

## Interpretation, fitting and observational equivalence

The predicted quantity is conditional equilibrium height loss for identical cells under an effective affine optical-dose-to-lateral-tension law. The source outcome is a maximum over space and time after cellularization-reference correction. Equating these quantities is an explicitly testable auxiliary null; it is not derived from the dynamic experiment. Fixed heights 32, 33.5 and 35 micrometres form a declared sensitivity analysis, not embryo-specific height measurements or a confidence set.

The finite-grid fit selects a minimum only over the declared grid. It does not prove continuum optimality or uniqueness. At small dose,

\[
S(I)=\frac{ha}{3/2-\nu/2}I+O(I^2).
\]

A linear training response therefore constrains this combination, not \(a\) and \(\nu\) independently. Nonlinear data may provide additional information; no general impossibility theorem for that fit is asserted. More fundamentally, all physical parameter combinations yielding the same \((h,a,\nu)\) produce the same scalar curve. A common positive scaling of all energy/active-force parameters is one such equivalence.

The homogeneous model contains no diameter variable, tissue embedding, compatible neighboring-cell forces or curvature. Conditional on dose, its response is identical at both activation sizes. This is the independent-column null being tested, not a successful derivation of size-dependent shape. Independent heterogeneous cells can also give diameter-dependent maxima through ordinary order statistics; rejecting the homogeneous null would not uniquely demonstrate mechanical coupling.

The optical normalizer is supplied by the publication across experiments. Published outcomes were previously inspected. Prediction-file separation verifies computational withholding only; it is not human blinding or external preregistration. Both training and evaluation must use maximum shortening, and the repeated 18-micrometre panels cannot be counted as independent replicates.

## Code-review boundary

Initial review requested rejection of nonfinite heights, intensities and gains rather than silently generating NaNs, explicit labeling of the limiting material endpoint, and propagation of the protocol iteration count to evaluation calls. The corrected version 2 implements these requirements; its NaN guard tests pass. At version 1's value of 60, the original default and declared iteration count already coincided, so that correction changes no valid version-1 response. Earlier files/locks must be retained, and version-2 predictions must be freshly frozen before evaluation scoring, with final source/protocol hashes recorded.

The independent review receipt records these findings, 135 high-precision test cases, numerical discrepancies and the exact reviewed source/protocol hashes. Its maximum absolute shortening difference is less than 8.56e-15 micrometres across the declared cases. No held-out response was used in this mathematical review and no experimental score is calculated here. These checks do not validate a tissue shape, a biochemical mechanism or equilibrium-dynamics equivalence.
