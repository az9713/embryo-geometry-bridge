# The four small-region observations: source definitions and model contract

**The public source supports four distinct mechanical summaries, but does not completely specify one reproducible four-output operator.** Maximum shortening has a clear spatial-and-temporal definition. The ordering of averaging and temporal maximization, component-reference construction, spatial weights and ROI transport require explicit assumptions or alternatives. A model can implement these alternatives precisely; numerical precision does not by itself establish equivalence to the original analysis.

This is a source/observation audit, not a fit or a choice of discrepancy allowances. It reads no held-out response values and changes no frozen protocol. The [audit receipt](../evidence/joint_observation_contract_audit.json) retains exact locators, source hashes, small-block workbook checks, illustrative trace comparisons and the distinctions below.

## What is being imaged?

Countryman et al., [Nature Communications 16, 7665 (2025)](https://doi.org/10.1038/s41467-025-62483-6), activate lateral membranes in a small region of the dorsal Drosophila embryonic epithelium. “Apical” is the outer surface near the embryo's envelope; “basal” is the inner cellular boundary. During cellularization, cells continue extending basally, so height may increase even without activation. A response relative to the initial image therefore differs from shortening relative to an unactivated counterfactual.

Primary locators are Methods Sec15, **Furrow segmentation and tracking**; Sec16, **Furrow quantification**; Sec17, **Activation quantification**; Results Sec5; and Fig2a–e. Sec15 specifies automated apical segmentation, manual confirmation/correction, manual basal segmentation and PIVlab motion correction. Sec16 restricts height measurements to coverslip-adjacent tissue to reduce curvature artifacts. These facts do not provide the actual segmentation masks, registration fields or source-analysis code.

## Definitions established by the source

Use a registered image coordinate \(x\) and time \(t\). With the upward/apical-positive z convention shown in Fig2b, write the apical and basal surface heights as \(z_a(x,t)\) and \(z_b(x,t)\). The source-compatible vertical height is

\[
H(x,t)=z_a(x,t)-z_b(x,t).
\]

This is a difference in z-position. It is not a surface-normal thickness, a curved cell-fiber length, or the distance between unrelated material points. Let \(H_{\rm ref}(x,t)\) denote the embryo-specific no-activation reference height, constructed using cellularization in nonactivated tissue. Positive shortening is

\[
s(x,t)=H_{\rm ref}(x,t)-H(x,t).
\]

The article describes an embryo-specific reference **curve** and cellularization rate, rather than publishing a full spatial reference field; the \(x\)-dependence above permits a declared spatial implementation. A spatially uniform curve is one possible convention, not an additional measured height map.

| Quantity | Source-supported content | Remaining implementation detail |
|---|---|---|
|Maximum shortening \(M\)|Sec16 maximizes shortening over both time and position within the tracked activation region. Fig2d contrasts this with region-averaged shortening.|Exact masks, discrete/interpolated time treatment, ROI transport, treatment of missing boundary pixels and the reference-growth fit.|
|Region-averaged shortening \(V\)|Sec16 describes a temporally maximal shortening summary averaged over points in the activation region; Fig2c/d and the workbook explicitly label it region-averaged.|The prose does not unambiguously order maximum and averaging. It does not give area-versus-section weights or a full ROI sampling rule.|
|Apical shortening \(A\)|Fig2e caption calls it maximum shortening of the apical surface. Fig2b explicitly plots **regional mean** surface positions over time.|Whether Fig2e takes a temporal maximum of the mean motion or a pointwise spatial/temporal maximum; exact apical reference and any drift/growth treatment.|
|Basal shortening \(B\)|Fig2e gives the basal counterpart; Sec16 recognizes ongoing basal extension and uses a reference-height correction.|Its component-specific counterfactual, aggregation and timing; the height reference alone does not specify these.|

For a specified tracked ROI \(D(t)\), the clearest total-maximum operator is

\[
M=\sup_{t\in\mathcal T}\ \sup_{x\in D(t)}s(x,t).
\]

The experimental window is finite: Fig2a indicates preactivation acquisition, a 40 s activation step and approximately 5 min response imaging; Sec5 describes stopping before later developmental compression. These schematic labels do not define one exact common end time. The published weak, moderate and strong example traces span respectively −44.3 to 402.6 s, −71.7 to 517.9 s and −44.5 to 321 s. Per-record cohort windows and reference fits are not tabulated. An equilibrium at one selected load is not automatically this temporal maximum.

## The averaging ambiguity is mathematically consequential

For a fixed registered domain and normalized nonnegative spatial weight \(w\), two possible readings are

\[
V_1=\sup_t\int_Ds(x,t)w(x)\,dx,
\qquad
V_2=\int_D\sup_ts(x,t)w(x)\,dx.
\]

In general \(V_1\le V_2\le M\), not \(V_1=V_2\). For two equally weighted locations, suppose shortening is \((2,0)\) at one time and \((0,2)\) at another. Then \(V_1=1\), \(V_2=2\), and \(M=2\). This exact counterexample shows why a wording ambiguity cannot be resolved by swapping the operations. Moving domains require an explicit registration/material chart before \(V_2\) can even be defined consistently across time.

The public table alone does not distinguish these cases. A static-state observer has no temporal ordering, so both reduce to the same spatial mean at that state. Taking the maximum of static-state means over a supplied path implements \(V_1\), not necessarily \(V_2\).

## Component references and aggregation

If component reference surfaces \(z_{a,\rm ref}\), \(z_{b,\rm ref}\) are declared with difference \(H_{\rm ref}\), define signed motions

\[
a=z_{a,\rm ref}-z_a,\qquad b=z_b-z_{b,\rm ref},\qquad s=a+b.
\]

Apical inward/downward motion and basal outward/upward motion then contribute positively to shortening. Local values may be negative; the source does not justify silently clipping every local displacement or shortening to its positive part.

The height reference alone does not determine the component references: adding a common function \(c(x,t)\) to both references leaves \(H_{\rm ref}\) and \(s\) unchanged, while changing \(a\) and \(b\) oppositely. This is a statement about the information in a height difference, **not** evidence that the authors failed to register their images. The published Fig2b coordinate frame and preactivation traces provide meaningful anchors. A fixed preactivation apical reference with the growth correction allocated basally is a plausible model convention, but its exact use in Fig2e is not stated in the inspected methods/caption/workbook. The audit does not assume that \(c\) is constant in time.

Two component-aggregation alternatives must be distinguished:

\[
A_{\rm point}=\sup_{t,x}a(x,t),\quad B_{\rm point}=\sup_{t,x}b(x,t),
\]

\[
A_{\rm mean}=\sup_t\int_Da(x,t)w(x)\,dx,\quad
B_{\rm mean}=\sup_t\int_Db(x,t)w(x)\,dx.
\]

They are generally different. Under the pointwise convention, \(M\le A_{\rm point}+B_{\rm point}\), but equality need not hold because the extrema can occur at different places/times. Nor may one assume \(M\ge A_{\rm point}\): the basal contribution at the apical maximum can be negative. Under the regional-mean convention the appropriate sum bounds \(V_1\), not necessarily the spatial maximum \(M\).

There is a useful but limited source clue. The published weak/moderate/strong Fig2b **regional-mean apical traces** have max-minus-min excursions 3, 6, 13 µm. These match Fig2e entries in source columns B,F,K; the trace captions' rounded doses 0.15, 0.46, 1.29 also match those columns' doses. This supports investigating a regional-mean component interpretation. It does not prove specimen identity or the exact reference/aggregation rule. In particular, max-minus-min of a trace is not generally the same as change from its chosen preactivation reference.

## What the workbook verifies—and what its column order assumes

Only the small-region blocks were read in this audit:

| Sheet/block | Source rows and columns | Verified content |
|---|---|---|
|Fig2c|A4 dose label; A5 regional-mean label; B:R observations|17 dose/mean pairs.|
|Fig2d|A4 maximum label; A5 regional-mean label; B:R|The same 17 regional means, exactly, with paired maximum shortening.|
|Fig2e|A4 apical label; A5 basal label; B:R|17 paired apical/basal entries. No dose or specimen-ID crosskey is provided in this block.|
|Fig4h, small block only|A3 identifies 18 µm region; A4 dose; A5 shortening; B:R|Dose/max entries exactly duplicate the small Fig2c/d data. Rows containing larger-region responses were not read.|

Fig2c and 2d link through the same 17 unique regional-mean values; Fig4h's small block independently repeats the dose/maximum arrays. Fig2e's component pairs have no such dose/mean/maximum key. Their linkage to the other panels relies on the author's common B:R ordering. The canonical source-column table retains that assumption; it is not upgraded to a verified biological-ID crosswalk.

The apical/basal pairing **within** each Fig2e column is author-supported: it is one paired coordinate of that panel's scatter plot. Preserving the common B:R ordering across panels is a reasonable documented default. A conservative identity-sensitivity alternative is one unknown permutation of the 17 **paired** component columns relative to the verified dose/mean/maximum records. This does not permit apical and basal values to be shuffled independently or claim the authors' ordering was wrong. Under that alternative, exact bipartite matching can test whether one cohort-wide permutation satisfies both component constraints. A perfect matching using guaranteed-compatible interval edges certifies existence; a Hall-deficient set even in the possible-compatible graph certifies impossibility. Intermediate cases remain unresolved. This is standard matching used to represent missing identity information, not a new inference principle or a selected fitting convention.

All 17 records satisfy mean≤maximum and maximum≤apical+basal in the canonical table. Four component sums equal the maximum; the others are strictly larger. These checks establish consistency with the pointwise convention, but do not prove it over other aggregation/reference choices. The four summaries are also not four statistically independent measurements: they share images, registration and reference calculations.

## A precise conditional static observer for the current model

Let an axisymmetric deformation produce radially single-valued apical/basal graphs \(Z_a(r),Z_b(r)\) over a common **current physical radius**. Let \(Z_a^0(r),Z_b^0(r)\) be declared reference graphs. For a declared current-radius ROI \(0\le r\le R\), define

\[
a(r)=Z_a^0(r)-Z_a(r),\quad
b(r)=Z_b(r)-Z_b^0(r),\quad
s(r)=a(r)+b(r).
\]

For the present flat static reference, choosing \(Z_a^0=H\), \(Z_b^0=0\) means no cellularization growth and no added reference drift. That is an explicit idealization. Compute \(M=\max_rs(r)\), and choose a declared spatial mean:

\[
\langle s\rangle_{\rm line}=\frac1R\int_0^Rs(r)\,dr,
\qquad
\langle s\rangle_{\rm disk}=\frac2{R^2}\int_0^Rr\,s(r)\,dr.
\]

The first is the mean along a symmetric diameter section; the second is the projected-area mean over a disk. Neither is automatically the source's sampling weight. Averaging over actual curved surface area is a third, different convention and is not supplied by the phrase region-averaged. Component outputs can similarly be either maxima or regional means. One observer variant therefore consists of an explicit ROI convention, spatial weighting, component aggregation and reference convention.

The old `maximum_activated_column_shortening` pairs top and bottom material nodes sharing an **initial** radial label. Their current radii can differ. It is not generally the difference between the two graphs at one current radius. Likewise, `central_apical_displacement` and `central_basal_displacement` measure axis motion, not extrema or regional means. Even an old maximum over initial-radius boundary nodes remains a conditional material-node proxy. The current observer must invert the separate radial surface maps and evaluate both boundaries at the same current radius. Sampled piecewise-linear interpolation of a quadratic boundary is not its exact geometry.

## Quasistatic peak extension and limits

For a declared common branch indexed by activation amplitude \(\alpha\in[0,gI]\), a conditional peak operator replaces experimental time maximization with maximization over that branch. For example,

\[
M_{\rm qs}(I)=\sup_{0\le\alpha\le gI}\max_rs(r;\alpha),\qquad
V_{1,\rm qs}(I)=\sup_{0\le\alpha\le gI}\langle s(\cdot;\alpha)\rangle.
\]

Regional-mean components require separate maxima over \(\alpha\) of their means; pointwise components require separate maxima over radius and \(\alpha\). Summing those separately maximized components does not reconstruct \(M_{\rm qs}\). A local-peak-then-average variant needs a shared radial chart across states.

If only a finite set of saved states is supplied, the output is a maximum over that **finite set**. It does not certify the continuous-branch supremum or experimental time-window coverage. Dynamics, pulse history, relaxation, growing reference heights, branches and specimen geometry can change the actual temporal response. None is established by renaming a static state a peak.

For joint inference, freeze one cohort-wide operator or explicitly carry a family of alternatives. “Compatible under at least one common operator variant” and “compatible under every allowed variant” have different quantifiers. Selecting a different convenient operator for each observation after seeing its residual would not be the same contract. This audit chooses no fitting tolerance, noise law, favored variant or parameter set; it supplies the definitions and unresolved source limits needed to make that subsequent calculation assessable.

The unmodified Fig2 source image already retained at `literature/material_constraints_countryman_fig2.jpg` was inspected locally and checked against the official checksum. No source pixels are embedded here or authorized for adapted public publication.
