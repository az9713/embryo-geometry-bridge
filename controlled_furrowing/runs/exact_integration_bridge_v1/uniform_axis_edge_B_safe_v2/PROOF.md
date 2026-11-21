# Preparation of six same-target physical-axis uniform contributions

**No actual six-element force/Hessian integration or cache sum has run.** The source `../../../scripts/exact_uniform_axis_edge_bridge_B_safe.py` prepares an integration adapter for full-axis-edge elements 72–77. It leaves all frozen kernels, domains, coordinates, material parameters and observations unchanged.

## Target and measure

The physical representation from `../axis_edge_physical_preparation/PROOF.md` gives the same P2 current map and the same free-coordinate derivatives, with affine ordinary physical polynomials q/B and cubic J. Only the spatial polynomial coordinates change. The common nodal box, C, fixed mu=1/kappa=5000, and shared gain-times-dose rule are identical to the existing uniform bulk calculation.

`PhysicalFunctional` inherits the reviewed finite-rule constructor, point-power evaluation, signed moment defect and polynomial application. For the explicit `axis_edge_physical` chart, the inherited constructor takes its physical-coordinate branch. The overridden `integral(a,b)` requires nonnegative integer exponents and evaluates

\[
\frac{2\pi}{\sigma}\int_{T_e}R^{a+1}(Z-Z_c)^b\,dR\,dZ.
\]

This includes the radial measure R exactly once, through `normalized_bulk_moment`. No Duffy Jacobian, second radial factor, inverse-R chart assumption or negative exponent is used at the axis. The reference triangle itself includes R=0; ordinary polynomial moments are integrable there.

Synthetic checks compare the physical mass with the old Duffy mass and verify exact rational agreement of all sixteen monomials with a,b from zero through three on an identity triangle. Actual field-defect moments are not needed for preparation: order selection uses only the positive geometry formula

\[
\int_{T_e}R\,dR\,dZ
=\frac{\det D_m}{6}\sum_{j=0}^2R_j
\]

and the sum of the original positive normalized auxiliary weights.

## Uniform norms and remainder policy

For every ordinary polynomial with interval coefficients,

\[
\sup_{T_e}|p|\le\sum_{i,j}|c_{ij}|_{\rm upper}
R_{\max}^i\zeta_{\max}^j.
\]

Coefficient upper bounds are exact rational values obtained from outward Arb absolute upper bounds. The implementation rejects negative exponents and nonfinite coefficients. It does not route the physical axis triangle through an off-axis inverse-R norm.

The scalar log-J derivative series and harmonic-coefficient tail bounds are reused unchanged. Each selected element uses the same freshly certified body-J interval and the same outward radius `ceil(rawrho*2^48)/2^48`. The declared order policy remains degrees 2–18, with force-tail target 1e-22 and individual Hessian-tail target 1e-16. These are numerical component targets; they are not measurement-error assumptions or a guarantee that the assembled contraction inequality will pass.

The first and second volumetric chain derivatives are contracted with signed material-integral-minus-auxiliary moments. Both positive measure masses multiply the scalar remainder bounds. The shear/activation polynomial terms use the identical U-series algebra and gain rule. All sparse interval coefficients, including zero-crossing balls, remain under the reviewed safe cleanup.

## Cache component and diagnostics

The cache part calls the reviewed `U.uniform_cache` without monkeypatching it. It retains the exact full-Y coefficientwise strain-map difference, all 64 original auxiliary sites, the three signed Hessian terms, and the same old-normalized weights and free-coordinate maps. The new representation changes no cache target.

All six consumed SAFE point packets are separately bound in the protocol. The saved point field and its drive are checked to lie in the same common nodal/material box. After a future run, force/Hessian packet overlaps are recorded as consistency diagnostics. Those comparisons do not prove a uniform bound; the analytic interval construction and tails supply that argument. A failed diagnostic remains explicit and cannot be promoted to nonexistence or root failure.

## Gate and retained development history

The actual `run` command requires an independent review whose PASS status, nonempty literal-True checks and source/protocol/preflight digests match exactly. All bound inputs are checked before and after execution; outputs are write-once. Geometry-review digests are also checked against the exact body/result/protocol files, and representation-test digests against the exact physical representation source/protocol.

A development synthetic test initially used Arb `==` to compare two identical nonexact auxiliary-mass enclosures. That is not endpoint identity; the test now compares packed exact endpoints. The source and failure were retained under `development_failure_0`. The first preparation-only revision is retained under `prelaunch_revision_0` before semantic receipt-to-file guards were tightened. Neither version ran an actual six-element integration or cache job.

The adapter prepares six possible replacements in an otherwise unchanged 144-element assembly. It establishes neither an improved numerical error bound before execution/review, nor any complete bulk-plus-surface/root, continuum, biological or novelty claim.
