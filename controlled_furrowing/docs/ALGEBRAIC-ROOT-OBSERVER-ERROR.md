# An actual root-error contribution to the static endpoints

**For one archived weak-load field, a nearby exact stationary point of the specified finite-quadrature potential now exists with a validated error bound.** Independent interval automatic differentiation reconstructs all 9,504 terms and verifies the contraction. Transferring the conservative root bound to the current-radius observer bounds its algebraic-error contribution by less than \(2.438\times10^{-10}\) µm for height maximum or height mean.

This tiny number measures numerical agreement with one mathematical target. It is not experimental precision or a continuum accuracy claim. Integration, spatial discretization, material/model error and source observation assumptions remain separate.

The evidence chain is [the producer proof](ALGEBRAIC-ROOT-CERTIFICATION-PATH.md), [independent interval-AD review](ALGEBRAIC-ROOT-INDEPENDENT-REVIEW.md), [exact transfer theorem](NODAL-OBSERVER-ERROR-TRANSFER.md), [transfer review](NODAL-OBSERVER-ERROR-REVIEW.md), and [actual application receipt](../evidence/algebraic_root_observer_error.json). The [application code](../scripts/apply_root_observer_error.py) requires both independent reviews and their source/certificate hashes before producing an output. The [application review](ROOT-OBSERVER-APPLICATION-REVIEW.md) checks the final composition separately.

## 1. The exact target and contraction

The target is the real-valued smooth energy using the frozen binary64 shape functions, quadrature caches, reference coordinates, parameters and activation strength of the archived 12×6, strength-one case. It has 611 free displacement coordinates. The center is the exact sum of the frozen reference and the archived free vector, rather than a second rounded serialization of that sum. It is neither ideal exact Gauss integration nor the continuum variational problem.

For gradient \(F(u)\), center \(u_0\), fixed binary preconditioner \(R\), and cube \(B_r(u_0)\), validated bounds give

\[
 \|RF(u_0)\|_\infty\le\eta,\qquad
 \sup_{u\in B_r}\|I-RF'(u)\|_\infty\le Z<1,
 \qquad\eta+Zr<r.
\]

Thus \(u\mapsto u-RF(u)\) maps the cube into itself and is a contraction. Banach's theorem supplies a unique fixed point there. The norm bound also makes \(RF'(u_0)\) invertible, hence R is nonsingular, so the fixed point satisfies \(F(u_*)=0\). The root distance is at most \(\eta/(1-Z)\). No local-minimum or global-uniqueness conclusion is inferred.

Two calculations verify the same target and cube of radius \(10^{-11}\):

| Calculation | Precision | Approximate \(\eta\) | Approximate Z | Root-distance upper bound, dimensionless |
|---|---:|---:|---:|---:|
| Analytic local derivatives with Arb balls | 160 bits | \(2.527710\times10^{-12}\) | 0.0351883 | \(2.619900\times10^{-12}\) |
| Independent second-order interval AD and cache reconstruction | 256 bits | \(2.527710\times10^{-12}\) | 0.0372269 | \(2.625447\times10^{-12}\) |

The values displayed here are rounded summaries. Exact outward bounds and strict acceptance inequalities are in the receipts. The downstream application takes the **larger** of the two exact error bounds. The two local roots coincide because both certificates concern the same center, target and uniqueness cube.

## 2. Add the serialization gap before converting units

Let \(r_*\) be that conservative dimensionless root bound. The largest difference between an exact reference-plus-displacement coordinate and its stored rounded `current_nodes` value is exactly

\[
 g_{\rm add}=\frac{8131}{36893488147419103232}
 \approx2.204\times10^{-16}.
\]

With fixed physical scale \(L_0=33.5\) µm, the nodal-coordinate bound relative to the stored observer field is

\[
 \varepsilon=L_0(r_*+g_{\rm add}).
\]

The application decodes and hashes the frozen cache, checks its reference and free vector against the archived state, and computes the addition gap with exact fractions. It verifies that the actual free-index set fixes the radial endpoints of both surfaces and that those endpoints agree exactly with the stored observer controls.

## 3. Transfer to the same current radius

For each surface, P2 interpolation bounds the coordinate error by \(5\varepsilon/4\) and the normalized-edge radial derivative error by \(8\varepsilon\). The actual base curves retain a strict positive radial margin after this perturbation. Fixed radial endpoints and a common positive-width ROI give common graph domains.

If the base graph has certified Lipschitz bound L, its Eulerian vertical error is at most

\[
 E=\frac54\varepsilon(1+L).
\]

The radial contribution is essential: simply comparing vertical nodal displacements would compare different current radii. Component maxima or means inherit the surface bound; height maximum and mean inherit the sum of the two surface bounds. Actual results are:

| Static endpoint | Safe rounded-up algebraic-error contribution, µm |
|---|---:|
| Apical component, spatial maximum or regional mean | \(<1.312\times10^{-10}\) |
| Basal component, spatial maximum or regional mean | \(<1.126\times10^{-10}\) |
| Maximum same-radius height shortening | \(<2.438\times10^{-10}\) |
| Line- or disk-mean same-radius height shortening | \(<2.438\times10^{-10}\) |

The JSON retains exact rational bounds and outward-rounded binary display values. A pre-application review found that nearest-rounded floats could lie below an exact upper bound; the corrected display explicitly rounds outward. The earlier source and finding are retained under `runs/root_observer_application/pre_outward_display/`.

## 4. What this closes, and what it leaves open

Previously a small force residual supplied no theorem relating the stored field to an exact stationary point. This calculation closes that link for **one** frozen finite-quadrature state and transfers the result to its static observer. It is a standard validated-numerics result applied to the project; no new general concept or originality claim follows.

It does not certify the new joint material/gain proposal at all seventeen doses, every compatible material combination, stability, a physical continuation branch, positive body Jacobian between all quadrature points for the unknown root, numerical integration accuracy, continuum spatial error, biological observation semantics or cross-size prediction. The original [five-gate objective](CURRENT-FIVE-GATE-AUDIT.md) remains incomplete.
