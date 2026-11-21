# Independent review: energy-norm gain tube to whole-body geometry

**The audited one-parameter energy-norm root family now has a verified whole-body geometric enclosure.** For the old weak-load material and relative linear gain increment \(s\in[-10^{-3},10^{-3}]\), every root supplied by that theorem lies in a nodal box with

\[
 0.992464257<J<1.006576818.
\]

The body is axis-preserving and globally injective throughout that box. This is not the five-parameter actual-column-M family, a statement about all 17 calibration inputs, or a continuum-accuracy result.

The [composition implementation](../scripts/compose_energy_geometry.py), [producer geometry receipt](../evidence/energy_gain_geometry.json) and [independent review](../evidence/energy_gain_geometry_independent_review.json) retain all coefficients, bounds and hashes. The independent composition review passes **1,858 checks**, including exact coordinate-radius tests for all 611 free coefficients and a fresh 144-element geometry replay. No root was optimized and no producer receipt was overwritten.

## 1. The prerequisite theorem and its coordinates

The [audited energy-gain theorem](ENERGY-GAIN-CELL-REVIEW.md) concerns a fixed old weak-load model, centered at normalized active coefficient \(f_0=1\). Its parameter is a **relative linear increment**,

\[
 f=f_0(1+s),\qquad |s|\le h=10^{-3}.
\]

It is not a finite log-gain interval. All other material parameters and geometry remain fixed. The archive's log-gain tangent equals the relative-gain derivative at the center, but the proof does not assume that its floating coefficients are an exact derivative: the frozen binary predictor's residual defect is explicitly bounded.

For each s, the theorem supplies a root of the form

\[
 u(s)=u_0+sS+Wz,\qquad
 \|z\|_2\le\rho=10^{-6}.
\]

Here \(u_0\), \(S\) and the invertible 611×611 matrix W are the **frozen binary** arrays bound by that theorem. The radius is an auxiliary metric-coordinate radius, not a physical displacement or an energy in joules.

The composition uses the full declared \(\rho\), rather than replacing it by a smaller approximate root-distance display. Its correctness therefore depends on no unreviewed tightening of the metric tube.

The 71-check energy review verifies uniform self-mapping/contraction and strict transformed-Hessian positivity for its finite binary-cache problem. Its local-energy-minimum conclusion is inherited from that theorem; the geometry composition adds no global minimum, dynamical stability or biological path-selection claim.

## 2. Exact coordinate radii

For free coordinate i,

\[
 |u_i(s)-u_{0,i}|
 \le h|S_i|+|(Wz)_i|
 \le h|S_i|+\rho
       \sqrt{\sum_jW_{ij}^2}
 =:\varepsilon_i.
\]

The second inequality is Cauchy–Schwarz. It uses the Euclidean norm of the **row** of W, not its column, its diagonal alone, an infinity norm or an assumed orthogonality.

Every squared row norm is evaluated as an exact rational sum of squares of serialized binary64 values. The square root is rounded outward to a rational value. The review independently constructs 100-decimal integer-square-root brackets and checks

\[
 c_i\ge0,\qquad c_i^2\ge\sum_jW_{ij}^2
\]

for each producer upper root \(c_i\). Thus each \(\varepsilon_i\) is a proved upper bound, not a nearest-rounded display.

Fixed coordinates have exactly zero radius. Free radii are mapped through the archived free-index array; axis radial coordinates therefore remain exactly fixed. Centers are exact reference coordinates plus the archived binary displacement vector, rather than separately rounded stored current nodes.

The largest coordinate deviation bound is less than

\[
 5.918405\times10^{-5}
\]

in normalized units, or 0.001983 µm after multiplication by the fixed length 33.5 µm. This combines **parameter-induced predictor motion and algebraic correction** over the family. It is not a solver-error estimate or experimental resolution.

## 3. Composition with the whole-body proof

The coordinate boxes contain the entire metric root tube:

\[
 \{u_0+sS+Wz:|s|\le h,\ \|z\|_2\le\rho\}
 \subseteq
 \prod_i[u_{0,i}-\varepsilon_i,u_{0,i}+\varepsilon_i].
\]

The independently reviewed [P2 geometry primitive](ROOT-BODY-GEOMETRY-REVIEW.md) is then applied to this new, larger box. It bounds the correlated degree-four Jacobian coefficients on every whole element, preserves the axis exactly, and certifies positive symmetric meridional gradient.

The resulting strong-monotonicity argument proves meridional injectivity. Axis preservation and positive interior radius then imply injectivity after revolution. The implementation verifies these properties for **every nodal choice in the box**; it need not pretend that every such choice is itself a root.

The actual J enclosure is approximately

\[
 [0.992464257604187,\ 1.0065768174032732].
\]

It lies inside \([0.99,1.01]\) for this family. This is neither exact incompressibility nor a result at larger loads or different materials. Interval endpoints need not be attained. The proof concerns off-axis Jacobians and their bounded limiting values; classical Cartesian \(C^1\) regularity on the piecewise-polynomial axis is not inferred.

## 4. The certificates must not be mixed across families

| Certificate | Parameters/center | Corrector coordinates | What it covers |
|---|---|---|---|
| Old weak point root | Old fixed material, \(f=1\) | Nodal infinity ball | One local algebraic root |
| Energy-gain family composed here | Same old material, \(f=1+s\), \(|s|\le0.001\) | Fixed-W Euclidean metric ball | One nearby root family and its body geometry |
| Actual column-M cell | New primary material, actual measured \(I_M\), five material/gain intervals | Its separately audited predictor/corrector tube | One source-matched local parameter cell |

The third row is established through [a different uniform-cell proof](PARAMETRIC-ROOT-CELLS.md). Its source dose, material center, gain parameterization and branch domain cannot be replaced by the old weak field's representative input one. Nor can these separate certificates be combined into a full five-parameter or 17-record identified set.

## 5. Independent verification and provenance

The review checks:

- the exact energy producer and 71-check independent-review bindings;
- every frozen energy input and source hash, W-cache and model-cache checksums;
- the archived predictor column, actual free-index map and old weak input state;
- all 611 rational Cauchy–Schwarz radius inequalities and independent square-root brackets;
- zero radii on all 39 fixed coordinates and exact radial-axis preservation;
- equality of a fresh 144-element geometry replay with the retained composition;
- the audited geometry packet, its dependencies and all bound inputs again at the end.

The first composition-review version had the same numerical result; three additional explicit dependency/end-of-run guards were added before finalization. That review is preserved under `runs/energy_gain_geometry_independent_review/revision_0/`. No underlying geometry, energy input or producer receipt changed.

Producer composition SHA256:

    f1b747fc05d6fe55e7fc1365c09e0b44c8d05e56b36e75174e48a610723224f7

Independent final PASS receipt SHA256:

    0f1fe93804504a0e7964102c520d89a72b183763343619890a61a7ab94b23fe7

The receipt retains the full replay program and exact per-coordinate bounds. The validated finite-sum theorem and geometric composition establish no quadrature-integration/PDE error, physical force calibration, experimental temporal law, all-branch exhaustiveness or confirmed novelty.
