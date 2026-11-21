# Four observable endpoints from the actual quadratic tissue geometry

The new observer computes vertical tissue height at the **same current physical radius** on the apical and basal surfaces. It does not replace that height by a material-fibre length, a difference between differently located material nodes, or the sum of two independently maximized surface motions.

`scripts/quadratic_endpoint_observer.py` is independent of the mechanical solver. **23 analytical/guard tests and 115 saved-state consistency checks pass.** Eight existing higher-quadrature fields have been processed without any new optimization, data fit or stronger-load job. The source is frozen at

```text
075dac2148ed26a3bfe83fadf1b97602e4b08536f2be38cb95c2054776c0b926
```

The observation operator remains **conditional**: the published source does not uniquely settle component aggregation, component reference split, ROI transport, spatial weights or temporal ordering of the mean. The implementation retains alternatives instead of silently selecting a convenient one.

## 1. What the source establishes and leaves unresolved

The authoritative audit is `JOINT-OBSERVATION-CONTRACT.md` and `evidence/joint_observation_contract_audit.json`. In [Countryman et al. (2025)](https://www.nature.com/articles/s41467-025-62483-6), Methods Sec15 describes apical/basal boundary extraction and PIV position correction. Sec16 defines vertical height using apical and basal positions at the same tissue/image location and constructs reference growth from nonactivated tissue. This supports a vertical-height operator, not a normal thickness or material-column length. PIV registration exists; its exact ROI transport and averaging measure are not fully specified by the prose.

Sec16 establishes a maximum over the available response window and ROI. It does not fully disambiguate maximum-over-time of a spatial mean from spatial mean of local temporal maxima. Actual illustrative trace windows differ, so the schematic five-minute protocol cannot supply one exact shared observation horizon.

The component interpretation is also unresolved. Fig2b explicitly presents **regional-average** apical/basal positions, while Fig2e labels maximum component shortening. Its endpoint columns could represent maxima of regional-mean motion rather than maxima over both space and time. The Fig2e workbook headers identify apical and basal shortening without a full aggregation formula. Therefore the implementation exposes both component conventions.

The source uses preactivation anchors and growth correction; height information alone does not uniquely determine an apical/basal counterfactual split. Our default flat initial component references are a declared model hypothesis, not a reconstruction of that source correction. Joining all four endpoints by specimen also retains the source-column pairing assumption identified in the audit.

## 2. Instantaneous quantities on a common physical ROI

Let \(R\) be current physical radius, with apical height \(a(R)\) and basal height \(b(R)\), positive apically. On a declared registered current-radius ROI \(I=[R_0,R_1]\), choose flat component references \(a_0>b_0\). Define the signed quantities

\[
 d_a(R)=a_0-a(R),\qquad d_b(R)=b(R)-b_0,
 \qquad s(R)=(a_0-b_0)-[a(R)-b(R)]=d_a(R)+d_b(R).
\]

The static maximum-height endpoint is

\[
 M=\sup_{R\in I}s(R).
\]

For component aggregation `spatial_max`, the other component endpoints are \(A=\sup_I d_a\), \(B=\sup_I d_b\). For `regional_mean`, they are the corresponding spatial means. No positive-part operation or clipping to zero is applied. Cells can lengthen, and a signed component or mean may be negative in a supplied state.

For a common positive weight \(w\), the fourth endpoint is

\[
 \overline s_w=\frac{\int_I s(R)w(R)dR}{\int_Iw(R)dR}.
\]

Two explicitly retained measures are:

- `current_radial_line`: \(w=1\), equivalent to a symmetric diameter-line average under the axisymmetric idealization.
- `current_disk`: \(w=R\), equivalent to current annular/disk area weighting; the constant factor \(2\pi\) cancels.

Thus the four-endpoint family has **four declared variants**: component spatial maximum or regional mean, crossed with current line or disk weighting. The source does not justify choosing a variant using larger-region outcomes. Disk weighting itself invokes the axisymmetric area idealization; it is not automatically the averaging measure of a coverslip-adjacent image section.

At one instant, \(\overline d_a+\overline d_b=\overline s\) under the same ROI and weights. In contrast, \(M\) generally does not equal \(A+B\). For example, \(d_a(R)=R,d_b(R)=1-R\) on [0,1] gives separate maxima one and one, while \(M=1\). Negative basal motion also prevents a universal assertion \(M\ge A\).

## 3. Exact quadratic inversion and supremum enclosures

Each boundary edge is represented by its three actual P2 nodal positions. Its current coordinates have the form

\[
 R(t)=r_0+r_1t+r_2t^2,\qquad Z(t)=z_0+z_1t+z_2t^2,
 \quad 0\le t\le1.
\]

The code checks exact endpoint continuity and monotone radius, then inverts \(R(t)\) with rational bisection. A physical-radius interval is enlarged to an enclosing parameter interval. Its height range is bounded by endpoint values and the quadratic critical point when present. Arithmetic uses exact rationals of supplied serialized values, so these are certificates for the discrete curves.

Component extrema are enclosed using all relevant segment endpoints and height critical points, including root-bracket uncertainty at ROI cuts. For \(s=a_0-b_0-a+b\), an interval cell uses the two graphs' independently enclosed height ranges. Attained point enclosures supply a lower bound on the supremum; the maximum cell upper bound supplies an upper bound. Adaptive bisection tightens the enclosure. If the two surfaces share exactly the same radial parameterization, their difference is itself quadratic and exact critical-point analysis is used directly.

The default requested maximum-endpoint width is 0.001 µm. The return object states whether it was achieved. A valid but wider enclosure is not silently presented as a sharp numerical value. This reporting tolerance is not a statistical discrepancy allowance or biological accuracy threshold.

The observer also certifies positive Eulerian thickness on the ROI. A radial fold, out-of-domain ROI, inconsistent reference height, mismatched saved controls/current nodes, nonconforming trace or unverified state is rejected. Failure to certify a very small positive thickness is distinguished from a proven nonpositive thickness; neither is replaced by the material-column proxy.

## 4. Eulerian means from polynomial antiderivatives

No approximate sampled Eulerian quadrature is needed. Because \(dR=R'(t)dt\),

\[
 \int Z(R)dR=\int Z(t)R'(t)dt,\qquad
 \int RZ(R)dR=\int R(t)Z(t)R'(t)dt.
\]

The integrands have degrees at most three and five, respectively. Their primitives therefore have degrees at most four and six, with exact rational coefficients. Primitive values at the inverse-root endpoint brackets are enclosed by interval Horner evaluation and summed over segments. This remains valid when heights or displacement components are negative.

The normalizers are \(R_1-R_0\) for the line and \((R_1^2-R_0^2)/2\) for the disk. Component means and mean height use the same operator, preserving the instantaneous additive identity. When ROI endpoints coincide with exact graph knots, many mean intervals reduce to exact rational values. Otherwise the only integration uncertainty comes from the rigorously bracketed algebraic inverse endpoints; this is not uncertainty about the continuum mechanics or the experiment.

## 5. Why radial motion makes a material-column difference wrong in general

An analytical test uses reference material radius \(r\in[0,1]\), flat reference height ten, and

\[
 R_a(r)=1.5r-0.5r^2,\quad Z_a(r)=10-r,
 \qquad R_b(r)=r,\quad Z_b(r)=r.
\]

Both surfaces span the same full current radial domain, but their radial motions differ. In the registered current ROI \(0\le R\le5/8\), the outer apical point comes from material \(r=1/2\), whereas the basal point comes from \(r=5/8\). The correct maximum vertical shortening is \(9/8\). Pairing material columns with reference \(r\le5/8\) instead gives \(5/4\). The analytical tests certify the discrepancy.

In several existing numerical fields, the height maximum happens at the axis, where both surfaces share current radius zero. Their legacy material-column maxima can consequently coincide with the Eulerian maximum. That numerical coincidence does not justify using the proxy for another field or ROI.

## 6. Finite supplied paths are not experimental time coverage

`observe_graph_path` and `observe_state_path` retain maxima over **only the explicitly supplied finite states**. The saved-state version requires verified states with identical model/mesh parameters and a shared registered ROI/reference. It reports the supplied strengths, but does not infer physical times, unsampled strengths, complete quasistatic coverage or the source's experimental observation windows.

The output separates

\[
 \max_j\overline{s_j}\quad\hbox{from}\quad
 \overline{\max_j s_j}.
\]

The latter is rigorously bounded using radial-cell dominance and exact integrals when one supplied field dominates, with interval integral bounds elsewhere. A test with \(s_1(R)=R\), \(s_2(R)=1/5\) on [0,1] gives line-weighted peak-of-mean \(1/2\), but mean-of-peaks \(13/25\). No commutation is assumed. The static API computes neither temporal maximum by default.

## 7. Existing-state diagnostics and retained evidence

The fixed-current-radius diagnostic ROI equals the initial patch radius, with flat initial apical/basal references. This is a documented candidate registration/reference convention, not a recovered PIV transport rule or arbitrary spatially varying counterfactual growth field. The following figures are diagnostic model outputs, not fitted biological measurements:

| Mesh | \(\kappa/\mu\) | \(f/\mu\) | Maximum Eulerian shortening, approximate µm | Disk mean µm | Line mean µm |
|---|---:|---:|---:|---:|---:|
| 12×6 | 500 | 1 | [2.207782, 2.208159] | 1.717544 | 1.883652 |
| 12×6 | 500 | 5 | 7.181539 | 5.700515 | 6.216024 |
| 12×6 | 500 | 10 | 10.650784 | 8.485103 | 9.250522 |
| 12×6 | 5000 | 1 | 1.991137 | 1.590343 | 1.724562 |
| 12×6 | 5000 | 5 | 7.005392 | 5.599190 | 6.077914 |
| 12×6 | 5000 | 10 | 10.490261 | 8.397368 | 9.120621 |
| 16×8 | 5000 | 10 | [10.647708, 10.647709] | 8.467979 | 9.224048 |
| 24×12 | 5000 | 1 | [2.330113, 2.330147] | 1.797690 | 1.976288 |

Displayed intervals are rounded outward. Singleton displays and mean values are rounded summaries; complete exact intervals are in the receipts. These diagnostics do not certify the pending stronger-load 24×12 targets.

At the 12×6, \(\kappa=5000,f=10\) state, spatial-maximum components are approximately (9.015839, 1.474422) µm, while disk-regional-mean components are approximately (7.228724, 1.168643) µm. The observation convention can therefore matter far more than the 0.001 µm observer resolution. It must be declared before inference.

Receipts are `quadratic_endpoint_observer_synthetic.json`, `quadratic_endpoint_observer_diagnostics.json` and `quadratic_endpoint_observer_verification.json` under `evidence/`. The latter replays all four variants on all eight fields, binds their input hashes, checks mean/component identities and requested resolution, and records the source-contract hash. No simulation is launched.

```python
from quadratic_endpoint_observer import observe_state
out = observe_state(saved_result,
    weighting="current_disk",             # alternative: current_radial_line
    component_aggregation="spatial_max",   # alternative: regional_mean
    tolerance_um="0.001")
endpoint_intervals = out["four_endpoints"]  # exact Fraction .lo and .hi
```

```powershell
python -X utf8 controlled_furrowing/scripts/quadratic_endpoint_observer.py --test --diagnostics
```

The observer makes the geometric alignment more faithful and the remaining measurement assumptions assessable. It does not establish a canonical source operator, experimental time coverage, a computed compatible parameter set, statistical confidence, continuum accuracy or novelty. Any joint fit must retain these conditional choices and numerical endpoint intervals.
