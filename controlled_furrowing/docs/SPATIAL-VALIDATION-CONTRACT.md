# Conditional spatial validation contract

**Version 2, synchronized to SPATIAL-PROTOCOL.json before spatial scoring.** This defines a test of an idealized finite-height axisymmetric mechanical model. It does not replace the five-gate PLAN, make the scalar null a shape solver, or assert specimen-level validation. Publication outcomes have already been seen. The additional computational separation is calibration-only fitting, saved predictions, then label scoring. This revision supersedes the unexecuted basal-modulus proposal {0,1,10}; the authoritative scenarios are {0.2,2,20}, not an outcome-selected subset.

## What the published bending label actually observes

Countryman et al. (2025), [Results, large-region activation, Fig. 4d–f](https://doi.org/10.1038/s41467-025-62483-6), classifies an apical-slice pattern: central tissue is separated from surrounding tissue by an acellular perivitelline annulus. Positive displacement, cell shortening and positive incremental curvature are not automatically that pattern. The image data do not provide a calibrated slice-height/time/point-spread-function operator for every record.

The retained supplementary Movie 3 is one illustrative 50-micrometre activation. The reviewed 04:40 frame has an apical-plane central island and surrounding dark annulus. The data audit also records recovery by 09:00. This is evidence that the displayed phenotype can be transient, not validation of an equilibrium observer for all 31 larger-region labels. The 2-micrometre sub-apical label is a depth descriptor, not a freely tunable bending threshold. Playback frame rate is not biological acquisition time. Sources and decoded-image hashes are in `evidence/data_movie3_inspection.json`.

## Coordinates, material positions and optical geometry

Use physical micrometres, with z increasing apically/outward. Initially flat idealized tissue occupies \(0\le z\le h\). Let reference material radius be \(R\); the deformed apical and basal parametric boundaries are \((r_a(R),z_a(R))\) and \((r_b(R),z_b(R))\). They must be converted to single-valued graphs on the **same deformed physical radial coordinate** \(r\):

\[
b(r)<z<a(r),\qquad 0\le r\le R_o.
\]

Using top and bottom values at the same original material R as if they were the same spatial column is generally incorrect. If both deformed radial maps are strictly increasing, their piecewise-linear graph representations can be merged on the union of radial breakpoints. Their minima and thickness bounds are then exact for those graph interpolants. A nonmonotone surface, crossed top/bottom graph or absent common radial domain makes this observer unevaluable. It is not assigned a negative bending label.

Positive local Jacobians alone do not prove global injectivity or absence of self-contact. Flat reference geometry, a remote clamped boundary and the absence of vitelline contact are explicit idealizations, not reconstructed specimen properties. If a vitelline contact claim is made, the independently supplied envelope and gap must be included; they cannot be inferred from a held-out bending image.

## Exact conditional apical-void observer

For continuous, positive-thickness graph surfaces define

\[
U=\min\{a(0),a(R_o)\},\qquad
L=\max\{b(0),b(R_o),\inf_{0<r<R_o}a(r)\},\qquad g=U-L.
\]

There exists a nondegenerate horizontal slice containing central and far-field tissue, with an intervening **apical-void** annulus, precisely when \(g>0\). Every \(z_*\in(L,U)\) is a witness. This is a conditional observation result, not an original scientific principle.

**Proof.** Central and far-field interior occupancy requires

\[
\max\{b(0),b(R_o)\}<z_*<\min\{a(0),a(R_o)\}.
\]

An apical void at an interior radius requires \(a(r)<z_*\), equivalent to \(z_*>\inf_{0<r<R_o}a(r)\). Continuity turns a strict interior inequality into an open radial interval, hence an annulus after revolution. Combining the inequalities gives \(L<z_*<U\). Conversely a slice with those properties satisfies every bound defining that interval. Strict inequalities exclude tangent-only, zero-margin cases. The minimum over the closed radial domain is interchangeable with the interior infimum here: a strict witness cannot have its decisive apical minimum only at an occupied endpoint.

The operator concerns voids above the apical surface. It does not cover basal cavities, nonaxisymmetric tissue plugs, optical membrane-intensity thresholds or overhanging surfaces. With exact piecewise-linear graph profiles, compute endpoint values and the apical minimum at vertices; do not substitute a center-minus-rim contrast.

## Missing slice metadata gives sets, not a binary point forecast

If the actual slice height is unknown, positive existence permits both labels \(\{0,1\}\): the acquisition could miss the interval. Exact absence permits only \(\{0\}\), conditional on this idealized operator. With an independently supplied admissible slice interval \(Z\), the possible positive heights are \(Z\cap(L,U)\); a positive point prediction requires an appropriate known acquisition plane and time, not merely a nonempty intersection.

A static negative equilibrium result does not rule out a transient positive phenotype. An observed positive label can conflict with the declared equilibrium/observer combination, but cannot by itself falsify the biological mechanism without temporal and optical correspondence. Report such comparisons as **conditional equilibrium-observer compatibility**, not validated bending-classification accuracy.

Required reporting includes, for every material family and record: whether geometry is evaluable; \(L,U,g\) in micrometres; a possible slice interval; the predicted label set; the source label only in the scorer; any conflict; and every missing temporal/optical assumption. Report singleton fraction, ambiguous fraction and invalid-geometry fraction. A high set-compatibility fraction with mostly \(\{0,1\}\) predictions is uninformative and must not be advertised as high accuracy.

## Numerical error and zero-margin cases

If true and computed surfaces on the same physical radial domain have uniform error at most \(\epsilon\), minima/maxima are 1-Lipschitz, giving

\[
|\widetilde L-L|\le\epsilon,\qquad
|\widetilde U-U|\le\epsilon,\qquad
|\widetilde g-g|\le2\epsilon.
\]

Consequently \(\widetilde g>2\epsilon\) certifies existence, and \(\widetilde g<-2\epsilon\) certifies absence, conditional on that actual uniform error bound. Inside the unresolved band neither is certified. At \(g=0\), absence of an open witness is exact for the numerical interpolant but is not robust to perturbations.

The difference between two meshes, or a small force residual, is not itself an error enclosure. Record mesh/domain/quadrature sensitivity separately. If only empirical sensitivity is available, label the sign as a discrete-model result; do not invent a tolerance from observed labels. An apical pit commonly has exactly zero computed observer gap, so a nominal negative discrete label cannot robustly exclude an arbitrarily small unseen central wiggle under a nonzero uniform error budget. Robust exclusion needs a strictly negative margin under that bound, or additional justified derivative/shape control. Source optical resolution is not converted into a threshold unless independently characterized.

For an explicit instability example, take \(a(r)=h\), \(b(r)=0\), so \(g=0\). For any \(0<\delta<h\), replace the apical surface by \(a_\delta(r)=h-\delta\sin^2(\pi r/R_o)\). Its uniform change is only \(\delta\), its endpoints stay at h, and its observer gap is \(\delta>0\). Arbitrarily small uniform perturbations therefore change exact absence into possible existence. A zero-gap discrete negative cannot be promoted to a robust continuum negative merely by testing several nearby meshes.

## Calibration outcome is a separate proxy

For paired material boundary markers the model can report

\[
S_{\rm col}^{eq}=\max_{R\le a_{act}}
\{[z_{a0}(R)-z_{b0}(R)]-[z_a(R)-z_b(R)]\}.
\]

This is a static vertical material-column shortening proxy. The source maximum is over spatial locations and time, relative to a cellularization reference. Different current radii of corresponding top and bottom markers, spatial tracking, cell-axis tilts and post-activation reference growth prevent automatic observational equality. Also return graph thickness, material-fibre length, central apical/basal displacements and their locations. Do not require separate measured maxima to add to the maximum height loss.

Fitting this proxy to the 17 small-region maximum-shortening records is a declared approximation. It neither validates quasistatic peak-load dynamics nor makes the larger-region surface prediction specimen-specific. A full source-matched shortening observer would additionally need the actual spatial tracking and temporal reference construction.

## Predeclared finite material family

Use the supplied optical dose I and one fitted nonnegative active-stress gain per material scenario. The primary proposed conditional family is

\[
h=33.5\ \mu{\rm m},\quad \mu=1,\quad
K_{\log J}/\mu=50,\quad K_a/(\mu h)=0.2,\quad
K_b/(\mu h)\in\{0.2,2,20\}.
\]

Here \(K_{\log J}\) is the coefficient of \((\log J)^2/2\), not an independently measured biological bulk modulus. The three basal moduli are equal to, ten times and one hundred times the declared apical area modulus. A zero-basal run can be a separate synthetic sensitivity/negative control, not a quantitative simulation of scraps RNAi or an extra calibrated family member. Additional heights, compressibility choices and initial curvatures are separately declared sensitivity scenarios, not confidence limits. All choices, the physical outer radius, mesh and loading geometry must be copied into a hashed execution protocol before fitting/prediction. A later material-family revision is separately versioned and cannot use larger labels.

For a thickness-one implementation, convert each recorded diameter D and source activation depth explicitly:

\[
a_{act}=\frac D{2h},\qquad
z_{min}=1-25/h,\qquad z_{max}=1-5/h.
\]

Use reported diameters where supplied, preserving nominal labels as metadata; no field is inferred from the held-out morphology. Reference curvature is zero in this conditional family. Actual embryo AP/DV curvature is unknown here and cannot be fitted from an activation outcome. A dimensionless outer radius of three corresponds to \(3h\) physically; domain-size sensitivity must be disclosed, including that convention for height sensitivity runs.

For each material member, fit only gain to the small-region scalar proxy and disclose its calibration residuals and numerical failures. The authoritative gain grid is 0 to 10 in steps of 0.02, with nonlinear response-strength tabulation initially 0 to 16 in steps of 0.5. The selected surrogate candidate is checked against direct nonlinear solves using the predeclared 0.05-micrometre interpolation tolerance. This numerical interpolation tolerance is not a biological measurement error. Boundaries, failed nodes and missing strength coverage remain explicit; extrapolation beyond a response table is not silently accepted. No continuum optimum is claimed.

A calibration-only revision is allowed only as an explicitly retained new protocol before the new predictions are locked. No member is chosen because it agrees with larger labels. The entire family, including poor calibration fits or unresolved solves, is reported. It is a finite set of material scenarios, **not a complete identified set, physical prior, confidence interval or uncertainty coverage guarantee**.

Mesh refinement can change the inferred gain. Report both fixed-gain numerical sensitivity and, if performed, the complete recalibration-on-small-data prediction sensitivity. Recalibration never uses larger outcomes and is frozen before scoring; it is not a large-condition-specific refit. Agreement after recalibration alone does not erase under-resolved displacement or stress fields.

The solver is finite-strain nonlinear. Its predictions are not assumed linear in gain. In a tangent linear limit, a scalar maximum proxy can have exact gain/material confounding after slope calibration; that observation does not prove identical fits to the full nonlinear 17-record curve.

## Freeze, evaluate and preserve the original gates

The prediction process may read only calibration observations and larger-record inputs: IDs, diameter and optical dose. It saves all fitted gains, family predictions, complete profiles, solver diagnostics and protocol/code/input hashes before a separate scorer reads larger labels. Published outcomes are already known, so this is data-flow separation only.

At minimum retain: material/reference radial grids; actual top/bottom and initial boundary coordinates; graph profiles or invalidity reason; activation geometry; all units/scales; minimum J and graph thickness; volume/work/free-force residuals; continuation/optimization status; mesh, quadrature and domain sensitivities; scalar calibration definition; and every family's observer interval and label set. A failed solve is not a negative phenotype. A single local equilibrium is not a global energy-minimum proof.

This test can establish a mathematically implemented finite-height conditional mechanics prediction and its compatibility limits. The absence of specimen surfaces, timing and optical correspondence remains a limitation of the original geometry gate. No status in PLAN.md is changed by this document, and no novelty claim follows from a successful numerical solve or an ambiguous compatible label set.
