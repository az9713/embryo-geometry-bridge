# Independent review of the actual-M parameter cell

**The uniform local-root proof and conservative observation composition pass independent review.** For every parameter in the declared small box, the selected source-M root family has a certified apical static prediction in approximately \([6.813612639694854,6.813612858035664]\) µm. The corresponding source value is 12 µm, giving a residual greater than 5.1863871 µm. This excludes that selected family at the explicitly declared 5.1 µm tolerance. It is not a complete identified material set, exclusion of other branches, or a result for all seventeen source columns.

The [independent script](../scripts/review_parametric_actual_M.py) reconstructs the calculation without calling the producer's `build_point`, `variation`, or `evaluate_cell`. The [receipt](../evidence/parametric_actual_M_independent_review.json) retains both a fresh 192-bit proof and a same-160-bit reproduction. Neither performs an equilibrium optimization or changes producer outputs.

## Exact target, parameters and source association

The five independent coordinates are absolute normalized apical area modulus, basal area modulus, apical surface tension, basal surface tension, and common gain. The first, second, fourth and fifth have relative half-width \(10^{-8}\); apical tension ranges from zero to \(2\times10^{-10}\). Every coordinate has positive width. The shear normalization \(\mu=1\), numerical volumetric coefficient \(\kappa=5000\), flat reference geometry, activation support, physical length scale and boundary constraints are fixed.

The center uses the revised-material saved state for `optoGEF2_small_col_M`. Its source dose is exactly 0.78208913. The parametric family uses the **exact rational product** of common gain and source dose. That product is not identified with the binary-rounded activation strength passed to the earlier floating solver. The independent replay verifies both quantities and their nonzero difference. The archived state is only a proposed center and preconditioner; the uniform proof validates the actual exact-product target.

The independent review checks the direct-run protocol, saved geometry/material state, gain, dose, cache, predictor archive and all bound source hashes. Source-column M is a workbook locator, not an independently established biological embryo identifier.

## Uniform predictor theorem

Let \(\theta\) lie anywhere in the declared closed five-dimensional cell, \(d=\theta-\theta_0\), and

\[
 \widehat u(\theta)=u_0+Sd.
\]

The stored binary predictor matrix is reproduced from the archived tangent columns and the documented binary coefficient divisions. Its numerical accuracy is not assumed. With \(F=\nabla E\), \(H_0=F_u(u_0,\theta_0)\), and \(B_0=F_\theta(u_0,\theta_0)\), the force is affine in the five coefficients at fixed geometry. Hence

\[
 F(\widehat u(\theta),\theta)=F_0+(H_0S+B_0)d+
 \int_0^1[H(u_0+tSd,\theta)-H_0]Sd\,dt.
\]

The parameter in this integral is the full \(\theta\), rather than \(\theta_0+td\). The producer's hull of \(d\) with zero encloses the entire geometric segment, including the one-sided apical-tension direction. The independent reconstruction retains that same quantifier.

All 9,504 local terms are freshly assembled. The five force-derivative columns are independently derived: the active drive contributes a dose-weighted quadratic derivative in the two fiber coordinates; surface coefficient derivatives are \((J_s-1)\nabla J_s\) and \(\nabla J_s\). They are projected through the actual local strain maps, with positive normalized quadrature weights.

For the exact binary preconditioner \(R\), interval bounds enclose both the explicit center/linear defect and the integrated Hessian remainder. A separate whole-tube calculation gives

\[
 \|RF(\widehat u(\theta),\theta)\|_\infty\le\eta,
 \quad
 \sup_{\|e\|_\infty\le r}
 \|I-RH(\widehat u(\theta)+e,\theta)\|_\infty\le Z.
\]

At \(r=10^{-11}\), the fresh 192-bit calculation gives approximately

\[
 \eta\le1.744242152\times10^{-12},\quad
 Z\le0.096862311,\quad
 \eta+Zr<2.713\times10^{-12}<r.
\]

These are safe rounded summaries. The exact rational interval endpoints in the receipt prove the strict inequalities. They hold **for every common parameter vector in the entire cell**; the use of interval overapproximations does not require different physical parameters at different quadrature points. The map \(u\mapsto u-RF(u,\theta)\) therefore has a unique fixed point in each predictor-centered cube. The derivative bound makes \(RH\), hence \(R\), invertible, so its fixed point is a stationary point. No other stationary branches outside that domain are excluded.

## Two precisions and a retained review failure

The first independent run used 192-bit Arb arithmetic and passed its contraction, self-map, domain and input-provenance gates. Four additional checks incorrectly demanded that interval **majorants** agree with a 160-bit run to an arbitrary absolute tolerance of \(10^{-25}\). Those checks failed and remain retained, with the original review source.

Different valid majorants need not be equal or overlap. In particular, Arb radius bounds use finite mantissa radius arithmetic, so increasing midpoint precision need not change each outward majorant monotonically. A separate 160-bit replay reproduced the producer's values. This diagnoses reproducibility; it is not the logical basis of the theorem. Acceptance depends on the same exact target/cell/predictor, a valid independent derivation, positive analytic domains, and strict self-map/contraction inequalities.

The initial fresh proof's corrector upper bound is slightly larger than the producer's tighter corrector upper bound. For downstream observer/body composition the safe approach is to use **exact rational affine-predictor coordinate extrema plus the full certified \(10^{-11}\) corrector cube**, and add the exact serialized-coordinate gap. This dominates both independently established root enclosures without relying on cross-precision equality of their tighter estimates.

## Actual conservative observation composition

The revised adapter uses exact rational extrema of every affine predictor coordinate over the exact parameter endpoints. It then adds the **full certified corrector cube** \(10^{-11}\), and each coordinate's exact difference between binary reference-plus-displacement and its rounded archived value. All 611 predictor bounds and all 650 physical-coordinate entries are independently reconstructed. This box contains both the 160- and 192-bit root enclosures. The earlier tighter composition and its source are retained under `observer_revision_0/`.

The actual whole-body nodal-box certificate is replayed exactly. Its generic geometry proof has a separate independent review. This establishes admissibility of every field in this box, including the axis; it is stronger than positivity only at the quadrature points. It does not establish continuum geometry accuracy.

The review extracts both P2 boundary traces directly from raw reference coordinates and triangle connectivity. It independently checks strict radial monotonicity, fixed radial endpoints, the same positive-width current physical ROI, and fixed reference levels. For the archived apical maximum, it reconstructs every quadratic polynomial, brackets inverse-radius endpoints with 120 exact rational bisections, and includes every admissible internal stationary point. Its independent graph-extremum enclosure lies inside the published base interval. The exact same-radius transport bound then encloses all parameter-dependent roots.

Source workbook cells `Fig 2c!M4` and `Fig 2e!M4` independently confirm the dose 0.78208913 and apical shortening 12 µm. The header says apical shortening; interpreting that quantity as the model's static spatial maximum remains a conditional observation contract, not an experimental fact proved by this calculation.

The final root and observer audits bind the canonical parameter-cell hash, model/observation contract, specific root/observer certificate bytes, selected branch domain and the universal parameter quantifier. The observer audit records admissibility of the full forward family. They may be consumed by a guarded provider; execution of a larger covering engine is a separate operation.

The certified interval gives

\[
 \inf_{\theta\text{ in this cell}}
 |A_{\rm model}(\theta)-12|
 \ge12-6.813612858035664>5.1863871\ {\rm µm}.
\]

The displayed decimal is illustrative; the actual inequality uses exact rational bounds. The declared tolerance 5.1 µm is not an estimated experimental noise level. One apical-channel contradiction is sufficient to reject this selected family at that tolerance without supplying the other 67 channels. It cannot declare a complete 68-channel fit or exclude unreviewed stationary branches.

## Scientific limits

The source temporal observation operator, canonical cross-sheet pairing, numerical compressibility, chosen flat geometry and fixed experimental dose remain explicit conditional assumptions. A local finite-binary-quadrature root certificate is not a continuum error certificate, experimental validation, a physically selected time-dependent path or a novelty claim.
