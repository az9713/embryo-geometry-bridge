# B proof-center translation and fresh-tube obligations

This package prepares one explicit numerical proof-center candidate. It performs **no Newton solve, residual evaluation, Hessian assembly, integration, root proof or full611-coordinate proof job**. All physical coefficients, the exact energy target, C, R, S, source dose and binary normalization are preserved.

The input is the complete point discrepancy D0=R(FI-FQ)(u0,theta0), independently reviewed in `evidence/exact_complete_point_B_independent_review.json`. Here FQ is the old finite binary-rule force and FI is the canonical exact-integral finite-element force. A proof about either discrete finite-dimensional model is still distinct from a PDE or biological validation.

## Exact two-stage rounding

For every one of611 coordinates the recipe is

\[
d_{{\rm requested},j}=\operatorname{fl}_{64}
\frac{D_{0,j}^{\rm lo}+D_{0,j}^{\rm hi}}2,
\qquad u_{1,j}=\operatorname{fl}_{64}(u_{0,j}-d_{{\rm requested},j}).
\]

The displacement used by every subsequent bound is the ACTUAL exact-binary difference

\[
d_j=Q(u_{0,j})-Q(u_{1,j}),\qquad
p_1(\theta)=p_0(\theta)-d,
\quad p_0(\theta)=u_0+S(\theta-\theta_0).
\]

The requested decrement and applied decrement are not assumed identical. The stored exact mismatch bound covers the full D0 interval and the final subtraction rounding.

For orientation, nearest-decimal approximations are:

- applied decrement infinity norm:1.673692434644991e-11;
- norm of D0 minus applied decrement, upper bound:3.4527682162304635e-18;
- proposed old-centered common radius r-old plus decrement norm:2.6736924346449912e-11.

Exact rational values in `rounding_ledger.json` and `preparation_receipt.json` are authoritative. `candidate_vectors.npz` retains original u0, candidate u1 and the requested binary decrement. `candidate_exact_nodes.json` gives exact Xbinary+P*u1, including fixed coordinates. `candidate_state.json` is explicitly marked NOT SOLVED; its current-nodes field is rounded serialization only.

## Conditional reuse lemma

Assume all of the following have been certified on the required common domains:

\[
\sup_{\theta\in C}\|R F_Q(p_0(\theta),\theta)\|_\infty
\le\eta_{\rm old},
\]

\[
\|I-RH_Q\|_\infty\le Z_Q,
\qquad \|R(H_I-H_Q)\|_\infty\le b,
\]

on every translation path from p0(theta) to p1(theta), and define

\[
V\ge\sup_{\theta\in C}
\|R(F_I-F_Q)(p_0(\theta),\theta)-D_0\|_\infty.
\]

Then

\[
\eta_{\rm new}\le
\eta_{\rm old}+\|D_0-d\|_\infty
+(Z_Q+b)\|d\|_\infty+V.
\]

To prove it, insert the old force and the point discrepancy at p0. The exact translation identity is

\[
RF_I(p_1,\theta)-RF_I(p_0,\theta)+d
=\int_0^1[I-RH_I(p_0-t d,\theta)]d\,dt.
\]

The Hessian assumptions give norm at most(ZQ+b) times norm(d). The remaining two terms are RFQ(p0,theta) and the discrepancy variation relative to D0. Taking norms proves the inequality.

The variation V is essential. For example FQ(u,theta)=u-theta and FI=FQ+epsilon+gamma*theta, with theta0=0 and d=epsilon, have ZQ=b=eta-old=point mismatch=0, yet the new predictor residual is gamma*theta. Omitting V would assert a false zero bound.

A route to V uses a joint path in the convex parameter cell, theta(t)=theta0+t*delta-theta and u(t)=u0+t*S*delta-theta:

\[
R\Delta F(p_0(\theta),\theta)-D_0
=\int_0^1 R[\Delta H\,S+\partial_\theta\Delta F]
\,\delta\theta\,dt.
\]

Thus parameter-force derivative differences, or another direct uniform discrepancy-variation proof, are needed. A base-point Hessian error alone does not establish V. Differentiability and full path/domain coverage are part of this conditional statement.

The package leaves ZQ-fresh,b-fresh,V and eta-new **unset** in `conditional_reuse_lemma.json`. It does not substitute the available old-domain ZQ or the pointwise Hessian discrepancy for these missing uniform bounds.

## The old domains do not cover the candidate tube

The proposed candidate root cube keeps the old normal radius1e-11 around p1(theta). It is not contained in the old scalar root-proof cube: the decrement norm itself exceeds that radius, and adding the proposed new radius exceeds it further.

The exact predictor hull plus candidate cube was separately compared against the WHOLE-coordinate boxes used by the old B body certificate. **396 of650 coordinate intervals escape those boxes.** `domain_obligations.json` records every failure and exact endpoint excess. The old boxes were first verified to enclose their original old root tube, so this is a containment failure of the proposed reuse, not a reconstruction mismatch.

A simple old-predictor-centered common radius is

\[
r_{\rm common}=r_{\rm old}+\|d\|_\infty.
\]

It encloses both the translation paths and the new same-radius cubes. This is a proposed new proof domain, not a certified enlargement of any old theorem.

## Minimum fresh checks

1. Uniform HQ contraction/preconditioned-Hessian bounds on a domain covering all old predictors, translation paths and proposed new cubes.
2. Uniform HI-HQ error b on that domain, with all bulk and surface integral remainders and finite-cache differences included.
3. The parameter/predictor discrepancy variation V, including parameter-force differences.
4. Whole canonical P2 positive-J/no-penalty, axis/fixed constraints and global embedding over the fresh coordinate boxes.
5. Positive surface hoop and nondegenerate tangents/G over all fresh boundary tubes.
6. Fresh graph monotonicity, top/bottom separation, fixed current-radius ROI and observer error transport before reusing measurement endpoints. The center translation must be included explicitly.
7. A separate uniform HI positivity/contraction argument for any exact-integral energy-minimum claim. The unchanged positive symmetric-R witness can remain a witness, but proves nothing about the new Hessian without the missing Hessian bounds.

The scalar synthetic fixtures test the conditional inequality, demonstrate the necessity of V, guard against unverified-domain substitution, and expose subtraction-rounding/intercept differences. They do not establish any missing physical-domain bound.

## Files and provenance

`protocol.json` freezes the arithmetic recipe and all input/source bindings before the candidate is formed. `unchanged_target.json` fingerprints R, S and every unchanged geometry/auxiliary array and records exact parameters, dose and normalization. `rounding_ledger.json` retains all611 exact decrements and predictor-intercept changes, plus nodal serialization gaps. `preparation_receipt.json` binds the candidate files and states that no root was certified.

Old roots, proofs, source data, scientific documents, production README/gates and publication state were not edited. The original full project goal remains active; this is preparation for a fresh proof, not completion of it.
