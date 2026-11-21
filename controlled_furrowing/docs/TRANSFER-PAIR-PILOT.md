# Small-to-large transfer in two frozen P2 models

At the same optical dose, the two stored large-patch fields differ by **0.040641 µm** in the prespecified apical center-to-shoulder contrast. The exact intervals for the stored curves are disjoint. The result resolves this numerical comparison for two fields in the finite model. The corresponding gap between equilibria, source compatibility, and scientific novelty remain open. [Exact analysis receipt](../runs/transfer_pair_pilot_v1/analysis/receipt.json)

## What was held fixed

The transfer question is whether two material and activation-gain choices informed by the small-patch calibration produce different shapes when the stimulated diameter grows. One candidate comes from the earlier gain-only calibration. The second changes the two area moduli, basal surface tension, and gain jointly. The dimensionless gains are approximately 8.062972 and 7.812655, respectively. The calculation compares the two complete parameter vectors. Each vector is fixed across diameters, so the observed differences combine material and gain effects. [Gain-only calibration](../runs/joint_gain_extension_v2/direct_results.json) · [Joint calibration protocol](../runs/material_cohort_direct_v1/protocol.json) · [Joint calibration results](../runs/material_cohort_direct_v1/results.json)

Among the 31 larger input records, the selection rule first takes the largest reported diameter. It then selects the dose nearest the median of the 17 small-calibration doses, 0.78208913, with record ID breaking a tie. This selects `fig4f_col_I`, diameter **49.8939753 µm**, at optical dose **0.89767374**. The same dose at **18.71024074 µm** is a matched counterfactual. The larger response values were withheld from the calculation and no source score was computed. [Input-only table](../data/derived/spatial_inputs_only.json) · [Selection receipt](../runs/exact_integration_cohort_v1/scientific_next_action/receipt.json) · [Frozen pilot protocol](../runs/transfer_pair_pilot_v1/preparation/protocol.json)

All four cases start independently from the zero-displacement reference and use the same optical-dose fractions (0,1/3,2/3,1) of the final dose. At each checkpoint the solver uses the product of the archived gain and requested dose, rounded to binary64. The common continuation rule can bisect a failed step in activation strength, which would add intermediate implied optical doses to the record. All twelve requested nonzero steps were accepted without bisection. These checkpoints describe a discrete numerical loading sequence. Continuous branch stability and the experimental time history remain unverified. [Driver](../scripts/transfer_pair_pilot.py) · [Preflight](../runs/transfer_pair_pilot_v1/preparation/preflight.json) · [Independent launch gate](../runs/transfer_pair_pilot_v1/prelaunch_independent_review.json) · [Checkpoint inventory](../runs/transfer_pair_pilot_v1/analysis/receipt.json)

P2 means the meridional displacement uses quadratic shape functions on each triangle. Both diameters use the same 12-by-6 mesh construction, six-node triangle connectivity, activation-depth layers, boundary constraints, and solver policy. The patch radius changes the radial coordinates so the activation edge stays aligned to the mesh. At the symmetry axis radial motion is fixed, and both components are fixed at the outer rim. [Geometry and numerical protocol](../runs/transfer_pair_pilot_v1/preparation/protocol.json)

## The signed observation

For current physical radius \(R\), let \(z_a(R)\) be the apical height, with positive height toward the apical side. The contrast in micrometres is

\[
Q=z_a(0)-z_a(a),\qquad a=24.94698765\,\mu\mathrm{m}.
\]

The shoulder radius \(a\) is half the reported larger diameter and is common to all four fields. It is a location in the deformed field. A negative \(Q\) means the center lies below that shoulder. [Frozen observation definition](../runs/transfer_pair_pilot_v1/preparation/protocol.json)

On a P2 boundary edge, the current radius \(R(t)\) and height \(z(t)\) are quadratic interpolants of the start, midpoint, and end nodes, with \(0\le t\le1\). The calculation treats the archived binary64 node values as exact rationals. A positive radial derivative along each edge makes the trace a single-valued function of current radius. The observer locates the edge containing \(a\), brackets the inverse of \(R(t)=a\) by 80 exact rational bisections, and bounds \(z(t)\) on that bracket from its endpoints and any interior quadratic extremum. The same graph supplies the center value at the axis. For interval values \(X=[x_-,x_+]\) and \(Y=[y_-,y_+]\), the enclosure for a difference is

\[
x_- - y_+ \;\leq\; x-y \;\leq\; x_+ - y_-\qquad (x\in X,\;y\in Y).
\]

Applying this rule to the center and shoulder yields the four \(Q\) intervals. Applying it again yields the diameter effects and the difference between those effects. Individual \(Q\) interval widths range from approximately \(4.7\times10^{-25}\) to \(2.8\times10^{-24}\,\mu\mathrm m\). This width comes from the inverse-radius query on the stored P2 curves. Bounds on the distance to an equilibrium root, integration and continuum error, and measurement uncertainty for \(Q\) remain unavailable. [Observer source](../scripts/quadratic_endpoint_observer.py) · [Exact endpoints and arithmetic](../runs/transfer_pair_pilot_v1/analysis/receipt.json)

## Results at the matched dose

The diameter effect is \(Q_{\mathrm{large}}-Q_{\mathrm{small}}\). All entries below are in micrometres and rounded for display.

| Candidate | \(Q\), 18.71024074 µm patch | \(Q\), 49.8939753 µm patch | Diameter effect |
|---|---:|---:|---:|
| Gain-only | −7.414204699092435 | −4.231347036069552 | +3.182857663022884 |
| Joint material and gain | −7.322311450872728 | −4.271987657501987 | +3.050323793370740 |
| Joint minus gain-only | +0.091893248219708 | −0.040640621432436 | −0.132533869652144 |

The larger patch brings \(Q\) closer to zero in both stored models. The joint candidate’s diameter effect is about **0.132534 µm smaller** than that of the gain-only candidate. The difference of diameter effects includes the small-diameter gap as well. [Exact contrast intervals](../runs/transfer_pair_pilot_v1/analysis/receipt.json)

![Apical and basal height against current radius for all four stored fields, showing the common shoulder radius.](../runs/transfer_pair_pilot_v1/analysis/current_radius_profiles.png)

![Signed center-to-shoulder contrast at both diameters and the diameter effect for each candidate.](../runs/transfer_pair_pilot_v1/analysis/shape_contrast_and_diameter_effect.png)

The profile figure joins retained samples of both surfaces over the full current-radius domain. The exact piecewise-quadratic controls are retained with each field. A bound on the maximum discrepancy over the complete profiles remains open. [Gain-only small field](../runs/transfer_pair_pilot_v1/jobs/baseline_small/checkpoint_03/field.json) · [Gain-only large field](../runs/transfer_pair_pilot_v1/jobs/baseline_large/checkpoint_03/field.json) · [Joint small field](../runs/transfer_pair_pilot_v1/jobs/joint_small/checkpoint_03/field.json) · [Joint large field](../runs/transfer_pair_pilot_v1/jobs/joint_large/checkpoint_03/field.json)

## What the checks establish

All four solver processes returned OS exit code zero. Across the final fields, the maximum normalized force residual is approximately \(3.88\times10^{-12}\) under the solve quadrature and \(6.77\times10^{-11}\) under the higher-order cross-evaluation, both below \(10^{-7}\). The stored-field checks also establish positive Jacobian, global embedding, monotone radial surface graphs, and positive thickness. Small force residuals test numerical stationarity at the returned states. Validating a nearby root and bounding its displacement would require additional estimates. The independent replay checked all 16 stored checkpoints, including 2,304 stored-element geometry checks and 32 force evaluations using the reviewed kernels. It independently reconstructed the P2 boundary graphs and contrast. [Supervisor terminal receipt](../runs/transfer_pair_pilot_v1/job.json) · [Terminal results](../runs/transfer_pair_pilot_v1/results.json) · [Independent saved-field replay](../runs/transfer_pair_pilot_v1/independent_actual/actual_review.json) · [Independent analysis readback](../runs/transfer_pair_pilot_v1/independent_actual/analysis_readback.json)

The earlier small-region checks compare 68 static shortening summaries per candidate. Their largest conditional stored-curve discrepancy requirements are about 5.37984 µm for gain-only and 5.18639 µm for the joint candidate. These are residual requirements for different observables. They leave both measurement uncertainty and a detectable threshold for \(Q\) unspecified. Establishing a source-compatible pair would require a common observation and uncertainty contract across all 68 small-region constraints, plus a comparable registered profile for the larger-response test. [Small-region residual receipt](../runs/exact_integration_cohort_v1/scientific_next_action/receipt.json) · [Source observation contract](JOINT-OBSERVATION-CONTRACT.md)

The tested pair therefore supplies a local, stored-field separation for one observable and one matched-dose comparison. A strong transfer-ambiguity result and a novelty finding still lack the necessary evidence. [Transfer claim criteria](TRANSFER-CLAIM-AUDIT.md)
