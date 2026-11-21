# Mesh sensitivity of the actual four endpoint functionals

**Scope: a target-free numerical audit of three saved high-quadrature equilibria.** No equilibrium, activation gain, source response or material parameter was fitted or recomputed. This component distinguishes the mesh sensitivity of conditional calibration observables from the full-domain shape discrepancy. It does not establish a continuum error bound or that biological inference is now adequate.

The [quadratic endpoint observer](QUADRATIC-ENDPOINT-OBSERVER.md) is frozen at `075dac2148ed26a3bfe83fadf1b97602e4b08536f2be38cb95c2054776c0b926`. The [exact graph-comparison helper](ACTIVATED-MESH-STUDY.md) is also unchanged. `scripts/endpoint_mesh_functionals.py` applies them to the saved 12x6, 16x8 and 24x12 states at activation 10 and volume coefficient 5000, with identical physical parameters and actual Duffy-8/surface-12 integration policies.

## The fixed observation contract

All three states use the same **current physical radius** interval

\[
 I=[0,A],\qquad A=33.5\,\mu\mathrm m\times\operatorname{Q}(\text{stored patch radius})
 \simeq9.35512037\,\mu\mathrm m.
\]

Here Q converts the stored binary number to its exact rational value; it does not replace it with a rounded decimal. The physical length unit is exactly the declared 33.5 micrometres. Flat references are \(a_0=33.5\,\mu\mathrm m\), \(b_0=0\).

This is a **conditional fixed registered current ROI**, equal in numerical radius to the initial activation patch. It is not a material-coordinate interval, a recovered PIV transport rule or a measured source ROI. The restriction wrapper changes only the graph object's declared domain; inherited exact quadratic inversion and range bounds remain unchanged. All comparison cells and display samples lie in this same interval. An analytical test verifies that deformation outside the ROI does not contaminate its reported supremum.

For top and bottom graphs \(a,b\), define signed fields

\[
 d_a=a_0-a,\qquad d_b=b-b_0,\qquad s=d_a+d_b.
\]

All four declared observer variants are retained: line versus disk weighting, crossed with spatial-maximum versus regional-mean components. In each variant the height endpoints remain \(M=\sup_I s\) and \(\overline s_w\). Component endpoints are either \((\sup_I d_a,\sup_I d_b)\) or their corresponding regional means. No positive-part clipping or temporal maximum is invented. These saved static states do not supply the source's measured response window or counterfactual growth correction.

## Functional error bounds on a fixed ROI

Let two stored graphs be \((a_1,b_1)\), \((a_2,b_2)\), with the same references, current ROI and weights. Suppose

\[
 \|a_2-a_1\|_{L^\infty(I)}\le\varepsilon_a,
 \qquad
 \|b_2-b_1\|_{L^\infty(I)}\le\varepsilon_b.
\]

Then \(\|s_2-s_1\|_\infty\le\varepsilon_a+\varepsilon_b\) by the triangle inequality. For any bounded functions \(f,g\),

\[
 f\le g+\|f-g\|_\infty
 \quad\Longrightarrow\quad
 \left|\sup_I f-\sup_I g\right|\le\|f-g\|_\infty,
\]

where the reverse inequality follows by exchanging f and g. No attained maximum or matching location of extrema is required.

For a nonnegative weight of positive finite mass, write \(d\mu_w=w(R)dR/\int_Iw\). Positivity and normalization imply

\[
 \left|\int_I f\,d\mu_w-\int_I g\,d\mu_w\right|
 \le\int_I|f-g|\,d\mu_w
 \le\|f-g\|_\infty.
\]

Therefore every declared variant satisfies

\[
 \begin{aligned}
 |A_2-A_1|&\le\varepsilon_a, & |B_2-B_1|&\le\varepsilon_b,\\
 |M_2-M_1|&\le\varepsilon_a+\varepsilon_b,
 & |\overline s_{w,2}-\overline s_{w,1}|&\le\varepsilon_a+\varepsilon_b.
 \end{aligned}
\]

The constants are sharp: the constant shifts \(a_2=a_1-\varepsilon_a\), \(b_2=b_1+\varepsilon_b\) shift the corresponding components by \(\varepsilon_a,\varepsilon_b\) and both height functionals by their sum. This is a standard supremum/positive-mean stability fact, not a new concept.

These bounds apply here to **differences between stored numerical fields**. Turning them into error bars for a true continuum solution would require independently established continuum surface-error bounds. Three meshes and small functional differences do not provide those bounds, justify Richardson extrapolation or prove identical local-branch selection.

## Measurement and convention errors are separate

If flat references also change by \(\delta a_0,\delta b_0\), the component budgets become \(\varepsilon_a+|\delta a_0|\), \(\varepsilon_b+|\delta b_0|\), and the height budget gains \(|\delta a_0-\delta b_0|\). Those errors are not estimated by mesh refinement. A source growth-correction convention can change the reference split while leaving a height measurement fixed.

Changing ROI, registration, temporal ordering or the spatial probability measure also falls outside the fixed-contract inequality. For two measures on a common domain, the standard total-variation bound is

\[
 \left|\int f\,d\mu-\int f\,d\nu\right|
 \le\operatorname{osc}(f)\,\|\mu-\nu\|_{\rm TV},
 \qquad\operatorname{osc}(f)=\sup f-\inf f,
\]

with \(\|\mu-\nu\|_{\rm TV}=\sup_E|\mu(E)-\nu(E)|\). To see this, subtract \(\inf f\), divide by its oscillation and use the layer-cake representation of a function between zero and one. Each level-set probability difference is bounded by total variation.

For centered line and disk averages on \([0,A]\), the densities become 1 and \(2t\) after \(t=R/A\). Thus

\[
 \|\mu_{\rm line}-\mu_{\rm disk}\|_{\rm TV}
 =\frac12\int_0^1|1-2t|dt=\frac14.
\]

The bound \(|\overline f_{\rm line}-\overline f_{\rm disk}|\le\operatorname{osc}(f)/4\) is sharp, for example for the indicator of \(t>1/2\); continuous functions approximate that extremizer. This is uncertainty about the observation convention, not a mesh error, and the source does not currently justify identifying one convention as canonical.

## Exact arithmetic, verification and limits

Observer intervals enclose the actual piecewise-quadratic discrete graphs using rational inversion, analytic polynomial antiderivatives and adaptive supremum bounds. Signed fine-minus-coarse intervals are obtained by interval subtraction, and absolute differences retain possible sign ambiguity. Curve suprema use the same exact current ROI, with requested report width 0.001 micrometres. That width is reporting resolution, not a biological acceptance tolerance.

Twenty-seven synthetic checks passed before processing the saved fields. They cover restricted constant-offset suprema, exclusion of an exterior deformation bump, radial reparameterization, malformed ROI rejection, sharp Lipschitz shifts for all four variants and the total-variation examples. Input, observer and comparison-helper hashes are checked before and after processing. All saved fields and the full-domain sensitivity receipts remain unchanged.

## What the saved-field comparisons show

The full-domain apical discrepancy is not used as the calibration-functional error. Restricting both surfaces to the identical current ROI reduces it, but does not remove component sensitivity. All intervals below are in micrometres and are rounded **outward** from the exact receipts.

| Mesh pair | Apical ROI sup difference | Basal ROI sup difference | Height-functional Lipschitz upper bound |
|---|---:|---:|---:|
| 12 to 16 | [0.225907, 0.226898] | [0.050587, 0.050588] | 0.277485 |
| 16 to 24 | [0.326502, 0.327122] | [0.123578, 0.124518] | 0.451640 |
| 12 to 24 | [0.515546, 0.516076] | [0.174166, 0.174167] | 0.690242 |

In particular, the 16-to-24 apical ROI sup difference is about 0.327 micrometres, while the full-domain value remains about 0.68 micrometres. Both are relevant to different questions. The former bounds fixed-ROI endpoint changes; the latter remains evidence of unresolved full-shape sensitivity. No field or initial numerical branch history was altered to obtain the narrower ROI.

### Spatial-maximum components and maximum height

These component values are identical under line and disk weighting, since no mean is taken for them. Approximate component values are shown to six decimals; height maxima retain outward intervals.

| Mesh | Apical spatial maximum | Basal spatial maximum | Maximum Eulerian height shortening |
|---|---:|---:|---:|
| 12x6 | 9.015839 | 1.474422 | [10.490260, 10.490261] |
| 16x8 | 9.223679 | 1.424029 | [10.647708, 10.647709] |
| 24x12 | 9.381109 | 1.300773 | [10.681602, 10.681883] |

### Regional-mean components and mean height

Together with the preceding table, these rows show all four declared variants. The height maximum in each row remains the preceding mesh-specific height maximum. Means below are approximate displays; exact rational endpoint intervals for every variant remain in the three observation files.

| Mesh | Measure | Apical regional mean | Basal regional mean | Regional mean height shortening |
|---|---|---:|---:|---:|
| 12x6 | Line | 7.853661 | 1.266960 | 9.120621 |
| 12x6 | Disk | 7.228724 | 1.168643 | 8.397368 |
| 16x8 | Line | 7.996534 | 1.227514 | 9.224048 |
| 16x8 | Disk | 7.335429 | 1.132550 | 8.467979 |
| 24x12 | Line | 8.175277 | 1.125618 | 9.300895 |
| 24x12 | Disk | 7.486893 | 1.041584 | 8.528477 |

![All four static endpoint conventions across the three saved meshes](../runs/endpoint_mesh_functionals/endpoint_mesh_values.png)

### Opposing component changes can hide sensitivity

The following are **signed** 24-minus-16 changes. A positive value means greater predicted shortening under the finer mesh. The component maximum values can occur at different radii, so their sum is not asserted to equal the maximum-height change.

| Component aggregation | Weight | Apical change | Basal change | Maximum-height change | Mean-height change |
|---|---|---:|---:|---:|---:|
| Spatial maximum | Line | [0.157430, 0.157431] | [-0.123257, -0.123256] | [0.033893, 0.034174] | [0.076847, 0.076848] |
| Regional mean | Line | [0.178742, 0.178743] | [-0.101896, -0.101895] | [0.033893, 0.034174] | [0.076847, 0.076848] |
| Spatial maximum | Disk | [0.157430, 0.157431] | [-0.123257, -0.123256] | [0.033893, 0.034174] | [0.060498, 0.060499] |
| Regional mean | Disk | [0.151464, 0.151465] | [-0.090966, -0.090965] | [0.033893, 0.034174] | [0.060498, 0.060499] |

The maximum-height change is only about 0.034 micrometres, while the apical and basal spatial-maximum components change by about +0.157 and -0.123 micrometres. The regional means likewise move in opposite directions. Thus stability of a combined endpoint can coexist with sensitivity in the joint component response. This matters when calibrating four observables together.

At 24x12 the line-versus-disk mean-height difference is approximately 0.772418 micrometres, compared with a 16-to-24 mesh change of approximately 0.076848 micrometres for the line mean and 0.060498 micrometres for the disk mean. The convention effect is not a statistical confidence interval, a measured source uncertainty distribution or a reason to choose the convention that fits better. It shows why the missing observation definition must be resolved before claiming a quantitative joint calibration.

![Restricted versus full-domain sensitivity and signed endpoint changes](../runs/endpoint_mesh_functionals/endpoint_mesh_sensitivity.png)

## Receipts and independently checked arithmetic

- `evidence/endpoint_mesh_functionals_protocol.json` locks the exact common ROI, references, four variants, reporting resolution and input/source hashes.
- `evidence/endpoint_mesh_functionals_synthetic.json`: all 27 analytical/guard checks pass.
- `evidence/endpoint_mesh_functionals_results.json`: all three mesh-pair comparisons, six surface sup enclosures and 48 endpoint differences, including their Lipschitz bounds.
- `runs/endpoint_mesh_functionals/mesh_{12,16,24}_observations.json`: twelve complete static observer results with exact rational bounds.
- `evidence/endpoint_mesh_functionals_verification.json`: all seven saved-field result gates pass; every requested supremum width was achieved and inputs/frozen sources are unchanged.
- `evidence/endpoint_mesh_functionals_independent.json`: 54 independent floating-arithmetic checks pass. Separate nodal polynomials, Brent inversion, SciPy quadrature and bounded numerical optimization check all 48 endpoints and six ROI sup intervals without calling the production observer or exact comparison helper. A 2e-8 micrometre crosscheck tolerance is declared; this is a sanity check of the exact certificates, not a second rigorous continuum bound.
- `evidence/endpoint_mesh_functionals_initial_report_failure.json`: an initial independent-report JSON serialization failure from a NumPy boolean is retained with the initial script. Casting the comparison to an ordinary boolean fixed reporting; no numerical method, tolerance, field or model changed.
- `evidence/endpoint_mesh_functionals_visuals.json`: plotting-source, result and image hashes.

The completed deliverable is a **verified conditional discrete error-budget component**. It supplies neither continuum calibration error bars nor a completed biological validation gate. Material uncertainty, source ROI/reference/weight conventions, transient response coverage, local branch identity and full spatial convergence remain distinct unresolved requirements. No novelty claim follows from these standard functional bounds.

