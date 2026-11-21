# Independent numeric review of the all-record material descent

**Verdict: PASS for the current local tangent and observation plumbing.** All seventeen source-column records are covered under one declared common baseline gain. The saved forces, material right-hand sides, tangent systems, smooth-mean chain rule and source/input/archive hashes independently replay. The output remains a proposed step for a finite affine surrogate; no new updated-material equilibrium has been solved or validated.

The retained independent program and evidence are in `evidence/material_cohort_descent_independent_review.json`. It contains **594 passing numeric/provenance checks** before the final document checks. This review does not duplicate the separate exact affine primal/dual certificate audit.

## Baseline and parameter coordinates

The protocol starts from the previously verified calibration-only candidate with lowest discrepancy among four already tested observation/gain candidates. That comparison does not establish the globally optimal nonlinear gain. Its common gain is approximately 8.062972020048665. For every record, the actual activation strength is exactly the floating representation of the retained rational gain times the retained decimal dose.

All seventeen unique canonical source-column IDs match the CSV and saved fields, including record I with zero apical entry and record O with zero regional mean. No record is omitted or assigned a replacement value. This verifies source-record coverage; it does not independently establish biological specimen identities or remove the documented cross-panel pairing assumption.

The same seventeen baseline fields are reobserved under all four declared observation conventions. Alternate conventions are therefore compared locally at this common baseline, rather than at their own globally optimized nonlinear material/gain solutions.

Mechanical archives retain six columns:

\[
 (\log\kappa,\log K_a,\log K_b,T_a,\log T_b,\log g).
\]

The affine step uses only columns 1 through 5 in zero-based indexing:

\[
 (\log K_a,\log K_b,T_a,\log T_b,\log g).
\]

The finite compression coefficient kappa=5000, geometry, reference, quadrature and length unit remain fixed. Its computed diagnostic tangent is **not** an active fitting direction. The additive apical-tension step is constrained to [0,0.02] because the baseline tension is zero; the other coordinates have a numerical +/-0.1 trust interval. These intervals are numerical step limits, not measured material priors or an exhaustive physical domain.

## Independent mechanical replay on every field

At each saved baseline, let

\[
 R(u,\theta)=\partial_u\frac{E(u,\theta)}{\mu V_0},\qquad
 H=\partial_uR,\qquad B=\partial_\theta R.
\]

All seventeen fresh force evaluations pass the unchanged 1e-7 threshold; the largest normalized force is **1.86e-10**. Every archive free vector and index array matches the saved field and current model. All actual integration policies are Duffy-8/surface-12; nominal parent fields 4/8 are not substituted for the actual integration policy.

For every field, an independent scalar-energy assembly differentiates all six parameter components: pointwise logarithmic volume, separate apical/basal area elasticity, additive apical true area, basal tension times true area and active fibre energy. It shares the checked geometric P2 cache but does not call the production coefficient-doubling force routine. Separate surface Gauss loops form the current revolved area directly. The largest relative difference from archived B across all 102 columns is **2.06e-11**.

All 102 archived linear systems satisfy

\[
 H S_j+B_j\simeq0,
\]

with maximum relative residual **9.19e-12**. Thirty fresh automatic-differentiation Hessian-vector products at five selected records independently check the actual energy along the archived tangent directions. The selection covers weak, intermediate and strong activation and the I/O source edge cases. Their largest relative residual is **9.46e-12**. This bypasses the archived sparse-Hessian multiplication rather than checking only that an archive is internally consistent.

The implicit-function interpretation remains conditional on an exact regular equilibrium. Floating small residuals, positive numerical geometry screens and successful factorizations do not enclose an exact root or its inverse Hessian. The source does not claim that stronger result.

## Exact smooth-mean contraction and graph finite differences

The exact physical-node mean gradients from the [geometric derivative module](QUADRATIC-MEAN-DERIVATIVES.md) are contracted independently on every record. The correct chain rule is

\[
 D_\theta O
 =\left(L_0D_{x_{\rm phys}}O\right)_{\rm free}S,
 \qquad L_0=33.5\,\mu\mathrm m.
\]

All **510** interval contractions match the saved analytic tangents: seventeen records, two spatial measures, three smooth means and five active directions. They preserve signs, physical units and the free-coordinate restriction. These are exact geometric derivative bounds along serialized floating tangent vectors; the floating equilibrium inverse does not thereby become rigorously enclosed.

On five selected records, fresh low-level quadratic graph integrals evaluate both measures along \(u\pm hS_j\). They use no fabricated convergence metadata and do not re-equilibrate a changed material. The largest smooth-mean finite-difference discrepancy is **4.24e-9 micrometres per declared coordinate**, below the retained 2e-6 threshold. References, current-radius ROI and length unit stay fixed throughout.

The source also retains two central-difference steps for all five directions and all four observation conventions. This review independently checks their interval-subtraction algebra, step order and positive local geometry-screen labels. The affine matrix uses exactly the declared smaller-step midpoint for all 68 rows; no record-wise choice of a convenient derivative is introduced.

## What the maximum differences do not prove

Spatial maxima can be nonsmooth when extrema compete. Agreement of two central-difference steps is only a numerical diagnostic. For example,

\[
 f(x)=|x|\quad\Longrightarrow\quad
 \frac{f(h)-f(-h)}{2h}=0
\]

for every positive h, even though the derivative at zero does not exist. The retained synthetic check explicitly exposes this issue. Neither small two-step discrepancies nor the exact interval arithmetic for each displaced graph proves differentiability of the maximum or bounds the error of a nonlinear material update.

The displaced graph vectors are labeled as nonstationary geometric probes. Their local orientation screens and graph validity checks do not turn them into newly solved equilibria. The one-sided physical apical-tension constraint is imposed on the proposed material step; using both displacement signs to check a geometric derivative does not assert a negative physical-tension branch.

## Affine descent is a proposal, not a new mechanical fit

Each of the four local affine problems has one common five-coordinate step and 68 constraints in the fixed seventeen-record/four-observable order. The independently checked plumbing selects the declared smaller-step derivative midpoints. A separate reviewer checks the rational primal/weak-dual certificates for those finite rounded affine problems.

Such a certificate can tightly bound the optimum of the specified affine model. It does not bound derivative approximation error, maximum nonsmoothness, exact-root error, nonlinear remainder, finite-element continuum error or a changed-material equilibrium. In particular, an affine decrease is not evidence that the actual nonlinear common-material discrepancy decreased. A directly solved and independently verified material update is still needed for that claim.

The current result performs no new nonlinear equilibrium solves, reads no larger-region response for this step, fixes kappa rather than fitting numerical compressibility, and leaves the original full compatible-set/propagation goal open. Trust-box limits and four conditional observation variants are not a physical posterior. No biological mechanism, novel concept or global parameter-domain result follows from this local proposal.

## Retained verification scope

The receipt hashes all seventeen baseline fields, mechanical archives, mechanical verification records and full observer-tangent files, as well as the protocol, terminal job, results and frozen runner. Current bytes still match the predeclared source/input references. The original runner and physical sources were not edited, and no job was relaunched by the reviewer.

The completed component is an independently checked cohort-wide local numerical direction. Its strength is coverage and audited derivative plumbing; its limits are exact-root/maximum/nonlinear/continuum uncertainty and unresolved experimental observation semantics.
