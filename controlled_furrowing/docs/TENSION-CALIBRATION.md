# Small-region adequacy audit of the prestress revision

**The numerical checks pass. Stronger basal surface tension improves the small-region motion-partition proxy, while substantial local compression remains.** This audit uses the17small-region observations only. No old or new larger-region response, prediction or score file was opened. All four declared scenarios are retained; no larger-outcome branch selection is performed.

The revision adds ordinary membrane surface tension to the existing finite-height model after the small-region v1 partition failure. It changes neither the original frozen v1 files nor the experimental data. Tension values are exploratory assumptions, not independently measured embryo properties.

## Inputs and reproducibility

The height is33.5µm; normalized apical/basal area moduli are0.2and2. The four basal tension scenarios are

\[
\widehat T_b=T_b/(\mu h)\in\{0,0.2,2,20\},
\]

with zero apical tension. Each scenario's single activation gain was fitted to the same17maximum-shortening measurements. The table-based gain fit and its direct-solve checks were independently replayed from small-region data.

The audit verifies all68direct profile files, their source/protocol provenance hashes, recorded identities and input strengths. It retains their SHA-256 hashes as the audit snapshot. The original calibration records did not contain an independent historical checksum manifest for every direct profile; this review captures and replays those profile bytes rather than claiming such a prior manifest existed.

The only larger-condition descriptor file accessed was hash-read because the calculation's source manifest references it. It contains input descriptors, not outcomes; the audit never parses that file. An explicit path guard rejects evaluation, prediction and score files. All accessed paths and the complete verification program are retained in [the audit receipt](../evidence/tension_calibration_adequacy.json).

Replay from the embryo-project root:

```powershell
python -c "import json;exec(json.load(open('controlled_furrowing/evidence/tension_calibration_adequacy.json'))['verification_code'])"
```

The replay refuses changed source/profile snapshots. It imports neither mechanical solver nor project observer.

## Numerical stationarity independently checked

A separate NumPy implementation reconstructs the bulk, area-elastic, active-fibre and surface-tension energies. The tension contribution uses the exact revolved straight-edge area, independently of the implementation's surface quadrature. On accepted positive-Jstates these energies have an analytic local extension. Complex-step differentiation checks every free coordinate:

\[
\frac{\partial\mathcal E}{\partial q_j}
\simeq\frac{\operatorname{Im}E(q+i\delta e_j)}{\delta\mu V_0},
\qquad\delta=10^{-30},\qquad\mathcal E=E/(\mu V_0).
\]

All **10,948 free-coordinate derivatives** across68profiles were checked. The largest independently recomputed scaled force norm is **3.63×10⁻¹⁰**, below the declared10⁻⁷criterion. The maximum difference between independently recomputed and reported force norms is **1.26×10⁻¹⁵**.

This is strong floating-point verification of the discrete energy/gradient and stationarity. It does not certify a global minimum, unique equilibrium, global injectivity or continuum accuracy. Local-Jand exact P1revolved-volume calculations were independently reconstructed as well.

## Same observation proxy as the v1 small-region audit

Within the initial activated material radius, the independently recomputed component observables are

\[
A=\max(-u_{z,a})h,\qquad B=\max(u_{z,b})h,
\qquad S=\max\{1-[z_a-z_b]\}h.
\]

These are respectively apical-down maximum, basal-up maximum and vertical paired-material-column shortening, using dimensionless model coordinates. The experimental comparison uses the corresponding separate source maxima. **A+B is not substituted for S**: the individual maxima may occur at different locations or times. The source data are temporal/image-derived; these model outputs are static material-column proxies.

The source is apical-dominant in16of17small records. “Direction agreement” below compares the orderingA>Bwith the source ordering of its two maxima. It is a calibration adequacy diagnostic, not a validated biological classification accuracy.

| Basal tension scenario | Gain | Height RMSE (µm) | Apical RMSE (µm) | Basal RMSE (µm) | Predicted apical-dominant | Direction agreement |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 8.90 | 2.255813 | 5.567822 | 3.727482 | 0/17 | 1/17 |
| 0.2 | 10.40 | 2.248892 | 4.575639 | 2.517067 | 0/17 | 1/17 |
| 2 | 14.30 | 2.240081 | 2.673956 | 1.316931 | 17/17 | 16/17 |
| 20 | 15.65 | 2.235579 | 2.622288 | 2.479043 | 17/17 | 16/17 |

The independently reproduced height residual sums of squares are86.507731,85.977765,85.305407and84.962817µm². Gains are interior to the declared finite grid. The selected direct/table discrepancies are at most0.01099µm, below the declared0.05µm tolerance; that is a check at selected calibration points, not a global interpolation or continuous-parameter optimum certificate.

The four height fits are numerically similar, yet their predicted motion partitions differ substantially. This demonstrates limited discrimination by the scalar target **within this tested family**. It is not a theorem describing the full calibration identified set. The component diagnostics may inform a later declared small-data calibration, but this audit does not introduce new weights, fit a tension, or select a winning branch.

## Why a surface-tension revision is mechanically sensible

For a flat surface displaced normally by a small graphw, letq=|∇w|². Then the area ratio is

\[
J_s=\sqrt{1+q}=1+\frac q2+O(q^2).
\]

An initially unstressed area-elastic penalty contributes

\[
\frac{K_s}{2}(J_s-1)^2=\frac{K_s}{8}|\nabla w|^4
+O(|\nabla w|^6),
\]

whereas a baseline surface tension contributes

\[
T_s(J_s-1)=\frac{T_s}{2}|\nabla w|^2+O(|\nabla w|^4).
\]

Thus these two standard terms supply different leading-order normal restoring responses. Increasing an unstressed area modulus is not the same as adding tension. This local expansion motivates the revision; it does not by itself prove the full tissue response or identify the embryo's prestress. The mechanism is established membrane mechanics, not a new-concept claim.

One editorial qualification is necessary: the frozen protocol inherits the phrase “flat stress-free idealization.” Nonzero surface-tension cases are **flat equilibrated prestressed reference states**, with rim reactions supporting the tension. They are not wholly stress-free. This note clarifies the wording without altering the frozen protocol or its hashes.

## Compression remains a major physical limitation

| Basal tension scenario | Minimum local J | Maximum local volume reduction | Largest absolute global volume change |
|---:|---:|---:|---:|
| 0 | 0.805932 | 19.41% | 0.0769% |
| 0.2 | 0.782001 | 21.80% | 0.0898% |
| 2 | 0.725455 | 27.45% | 0.1163% |
| 20 | 0.708673 | 29.13% | 0.1238% |

Global volume changes remain small while local compression increases. Surface tension has therefore **not resolved the incompressibility problem**. Positive localJand stationary energy do not justify calling these predictions biologically near-incompressible.

## What this establishes

Within the four tested discrete assumptions, stronger basal tension changes the small-region apical/basal partition in the observed majority direction. Similar scalar shortening fits alone would obscure that distinction. This is useful model-development evidence from calibration data, with no larger-outcome fitting in this audit.

It does not identify a physiological tension value, establish a full admissible parameter set, resolve the static-versus-temporal observation mismatch, demonstrate a successful larger-region prediction, or complete the original five gates. All scenarios and the failed/limited v1 remain retained. Further comparisons require the declared geometry/material/observation assumptions and their numerical limits to remain explicit.
