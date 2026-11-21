# What all sixteen linked public supplements supply for a geometry test

**The sixteen linked supplementary assets have now been acquired and checked. The next data prerequisite is a source-aligned observation contract.** These conclusions concern the inspected public version of Countryman et al., [*Endogenous OptoRhoGEFs reveal biophysical principles of epithelial tissue furrowing*](https://doi.org/10.1038/s41467-025-62483-6), *Nature Communications* 16, 7665 (2025). They do not establish that additional processed data are unavailable from the authors or elsewhere.

This chapter explains the experimental source before separating file integrity, image interpretation, and model validation. It adds no model fit, noise allowance, anatomical segmentation, or novelty claim. The existing shared-cell mechanical results and their discrepancy bounds are unchanged.

![Original schematic: public files, measurement bridge, and conditional mechanics](../runs/source_contract_supplements_v1/public-source-contract.png)

The schematic is drawn from this audit's logical structure. It contains no source microscopy pixels or adapted publication figures. Extracted source frames and contact sheets remain local inspection assets.

## The biological experiment in mathematical terms

An epithelium is a layer of connected cells. Here the cells form part of a fly embryo. Their **apical** surfaces face outward and their **basal** surfaces face inward. Light-controlled molecular recruitment changes cellular contraction. A mechanical model can ask how that contraction changes tissue shape, but it needs to know exactly which shapes the experimental measurements represent.

In the small-region experiments, the authors activated a nominal 18-micrometre-wide cylinder, spanning 5–25 micrometres below the apical surface. The published protocol selected embryos with cells 32–35 micrometres long, used a 40-second activation scan, and followed the response. The tissue is still developing: unactivated cells continue lengthening. Therefore a reduction relative to expected unactivated growth is different from the raw difference between two photographed heights. These details are supplied in the source; they are not missing physical ideas invented by our model. [Primary text, Results “RhoGEF-driven cell shortening is dose-dependent,” paragraph Par18; Methods “Live imaging and photomanipulation,” Par56; “Furrow quantification,” Par59.](https://pmc.ncbi.nlm.nih.gov/articles/PMC12361397/)

If reconstructed apical and basal positions are \(z_a(x,t)\) and \(z_b(x,t)\), a declared sign convention gives thickness

\[
h(x,t)=z_a(x,t)-z_b(x,t).
\]

Let \(h_{\mathrm{ref}}(x,t)\) denote the counterfactual height without activation, and choose positive shortening as

\[
S(x,t)=h_{\mathrm{ref}}(x,t)-h(x,t).
\]

This notation specifies the comparison we would replay; it does not assert that the unpublished per-acquisition reference arrays or their complete allocation to the two surface-component summaries are known. The paper describes an embryo-specific reference derived from nonactivated tissue. The reference curve, the tracked region, and their connection to each reported scalar are part of the measurement, not optional numerical details.

## What was actually verified

The [all-assets manifest](../evidence/public_supplement_asset_manifest.json) contains **16 files, 111,690,608 bytes: four PDFs, two XLSX workbooks, and ten AVI movies**. Each file matches the official object size and MD5, and has a retained SHA256. Matching a checksum identifies a file; it does not validate an anatomical interpretation of its pixels.

The [workbook package audit](../evidence/source_contract_workbook_package_audit.json) checked both archives and XML structures. They contain **40 visible worksheets** and **1,000 formula-cell elements**. Those formulas belong to the Fig. 2f scaling-axis construction; they are not an unreported segmentation pipeline or complete timecourse reconstruction. No hidden worksheets, hidden rows/columns, external relationships, comments, or embedded scientific files were identified by that structural inspection. A **922-byte Google workbook-metadata part remains semantically undecoded**. We retain that limit: uninterpreted bytes are not evidence that they contain no scientific information.

The [small-movie audit](../evidence/source_contract_small_movies_audit.json), [complete Movie 3 decode receipt](../evidence/movie_geometry_audit.json), and [Movies 4–10 audit](../evidence/source_contract_other_movies_audit.json) account for **903 fully decoded video frames**: 159 in Movie 1, 20 in Movie 2, 26 in Movie 3, and 698 in Movies 4–10. First/middle/last frames and local contact sheets were inspected under their recorded policies. All frames were decoded; this is not a claim that every frame was manually annotated. Movie 1's small cross-decoder RGB differences were quantified and retained rather than relabeled identical. Decoder agreement, where present, remains a software consistency check, not a biological measurement-error bound.

The supplied movies are prepared views, not the complete primary microscopy volumes. They show cell membranes, selected slices or projections, display clocks, scale bars and activation overlays. Some relevant boundaries are visible in Movie 2. However, complete, registered, calibrated initial and peak apical/basal profiles for the relevant source entries have **not been verified** from those views. We do not claim all contours are clipped or that every displayed boundary is unusable.

### What each movie contributes

| Movie | Published example and view | What it does not establish for the current 17-entry calibration |
|---|---|---|
| 1 | Unactivated gastrulation, apical average projection | A matched activated specimen or full apical/basal surfaces |
| 2 | OptoGEF2 and OptoCysts small activation; transverse and two in-plane views | Full cohort identity map, registered profiles and reference arrays |
| 3 | Large-region activation; transverse and selected in-plane slices | Certified anatomical annulus boundaries or matched fullcohort geometry |
| 4 | Large-region activation followed by **vitelline-membrane micro-perforation** | Basal/front-ablation experiment or unmodified small-dose response |
| 5 | Natural cephalic furrow in control and scra-RNAi embryos; two depths | Current small-activation cohort timecourses |
| 6 | Early cellularization, circular activation, one subapical slice | The 32–35-micrometre gastrulation-stage geometry |
| 7 | Early dorsal closure, amnioserosa activation, two projections; activation ellipse redefined | Stored time-dependent response masks or complete surfaces |
| 8 | Late dorsal closure, bilateral epidermal activation, subapical slice | Paired apical/basal trajectories for the calibration entries |
| 9 | Single amnioserosa-cell activation, apical projection | Multicellular patch calibration or a full height field |
| 10 | Five sequential dorsal-closure interventions, two projections | Five independent specimens or an unperturbed single-load history |

The source locators are the [description-of-supplementary-files PDF](../data/raw/41467_2025_62483_MOESM3_ESM.pdf), pages 1–2, and the individual audit receipts. The captions specify nominal biological conditions. File-container metadata specifies codec, dimensions and playback timing. These are different evidence types: the latter does not independently calibrate the former. A drawn activation outline is also not automatically the tracked region used to average the response.

## The source already tells us several important protocol facts

**Staging is constrained.** The 32–35-micrometre cell-length criterion is an experimentally used selection criterion, not permission to assign every source entry our idealized 33.5-micrometre slab. Per-specimen measured initial profiles are still needed to compare that idealization with observations.

**A timing protocol is supplied.** Fig. 2 gives the preactivation → activation scan → response sequence and a cyan 40-second activation interval. The explicit sentence identifying time zero with activation initiation belongs to the Fig. 6 caption for its experiments; it must not be transplanted into a purported Fig. 2 quotation. A model may declare activation start as its nominal zero, but the matching of each Fig. 2 acquisition's frames, clocks and stages to that convention still needs to be verified. The protocol is not wholly undefined. Contraction can already occur during illumination. A state at the end of the 40-second scan cannot be assumed undeformed, including when constructing a reference for responding specimens.

**Sampling is not accuracy.** Methods Par56 specifies 44-micrometre stacks with 2-micrometre slice spacing, taking 22 seconds per frame for the cited membrane/myosin recordings. Neither the axial step nor the frame interval is a certified error bar on reconstructed anatomy. A 2-micrometre axial sampling step does not imply a 2-micrometre contour-error bound; pixel or voxel pitch alone cannot supply such a statement.

**A substantial analysis workflow is described.** Methods Par58 specifies trained ilastik apical segmentation, object classification, manual confirmation/correction, manual basal tracing in ImageJ, and PIV registration. The corresponding project/training state, corrected masks and per-specimen registrations are needed for replay. “No custom code was generated” does not mean no image-processing workflow was used. Methods Par59 describes cellularization correction and shortening summaries; Par60 describes normalization of first-activation fluorescence into the reported dose.

**Temporal aggregation needs an implementation-level crosswalk.** The written definition of average shortening is consistent with averaging pointwise temporal peaks. It should not be dismissed as though the publication supplied no definition. Nevertheless, implementing that definition requires point tracking, the corrected field, the actual region and time window; surface-component summaries need their corresponding aggregation/reference conventions. In general,

\[
\int_\Omega \max_{t\in T}S(x,t)\,d\nu(x)
\ne
\max_{t\in T}\int_\Omega S(x,t)\,d\nu(x).
\]

Here \(\nu\) is the declared normalized spatial measure. Two spatial locations with equal weights can each peak at one at different times while their summed shortening never exceeds one. This gives left side one and right side one-half. This elementary difference explains why a static identity or static model residual cannot automatically reject independently time-aggregated source measurements.

The [peer-review audit](../evidence/source_contract_peer_review_audit.json) read all 36 pages with two text extractors and inspected five source-figure pages. Three author responses sharpen the constraints: measured 32–35-micrometre selection restricted developmental timing (PDF p. 17, first rebuttal, Reviewer 1 minor point 11); shortening could begin before the approximately 40-second scan ended (p. 31, second rebuttal, Reviewer 1 point 2); roughly 2-micrometre basal motion was difficult to resolve with 2-micrometre axial sampling (pp. 35–36, second rebuttal, Reviewer 2 point 3d). These clarify protocol and limitations; they do not supply the matched trajectory/mask/reference data, establish a segmentation error bound, or justify our illustrative 5.2-micrometre allowance. No peer-review text is quoted.

## The remaining observation bridge

| Requirement | Already supplied or verified | Missing or not yet verified for the matched source cohort |
|---|---|---|
| File provenance | All 16 linked assets, official checksums, local hashes | No claim that these exhaust primary microscopy or all possible repositories |
| Specimen identity | Canonical published table entries; captioned representative examples | Stable IDs linking each dose, scalar column, raw acquisition, processed series and representative view |
| Initial and response geometry | Stage/height selection, representative views, calibrated units in methods | Per-specimen registered apical/basal profiles or surfaces before activation and throughout the response |
| Acquisition time | Nominal activation onset, 40-second scan, stated sampling | Per-acquisition frame times, activation start/end and mapping of derived maxima to source frames |
| Region and registration | Activation geometry and PIV/segmentation descriptions | Stored activation masks, response masks, transforms and reference coordinates, with changing ROIs tracked |
| Developing-tissue reference | Embryo-specific nonactivated cellularization correction described | The chosen nonactivated region, fitted/reference trace and its use in each component summary |
| Observation operator | Written mean/maximum definitions and representative traces | Executable settings or transparent intermediate tables reproducing all four scalar summaries from each processed trajectory |
| Anatomical uncertainty | Nominal imaging spacing and manual-correction description | Repeated annotations or other justified contour/registration/reference error estimates; no invented bound |
| Mechanical interpretation | Qualitative perturbations and separate model calculations | Independent material inputs, common loading history, continuum accuracy and larger-region validation remain separate gaps |

The matrix says what the inspected package supports. It is not an accusation of unavailable author data, a proof that the opaque metadata contains nothing useful, or a requirement for a new experiment.

## A compact processed-data request that would close specific gaps

No message has been sent. A useful request can ask for existing processed outputs, not an indiscriminate transfer of huge image volumes:

1. A specimen/acquisition crosswalk for the small and larger activation cohorts: genotype, stage, source-table location, dose, activation geometry and matching image-series IDs.
2. Physical frame times, activation start/end, voxel calibration and the registration/coordinate transforms used for each acquisition.
3. Corrected apical and basal boundary coordinates or masks at each recorded time, including pre-activation frames and any field-of-view exclusions. Profiles are sufficient when those were the original analysis objects; full volumes are not intrinsically required.
4. Activation and measurement-region masks, tracked through time; clarify which spatial measure produced each regional mean.
5. The nonactivated reference region, measured cellularization trajectory, fitted/reference height curve and how its correction was allocated to the two component summaries.
6. A small intermediate table showing the trajectory-to-summary calculation, including the spatial/temporal aggregation order, selected time window, maximum locations, and any clipping, rounding, thresholding or exclusions.
7. Existing ilastik project/training settings, manual-correction record and uncertainty/replicate-annotation information where available. If no such uncertainty estimate exists, that absence should be stated rather than replaced by a guessed pixel envelope.

A minimal deliverable could be four compact tables—`specimens`, `frames`, `boundaries`, and `reference_and_summary`—plus masks/transforms where coordinates alone are insufficient. The purpose is to reproduce the published measurements under the authors' actual definitions before judging a mechanistic prediction.

## What remains true about the mathematics

The [shared-cell result](SHARED-CELL-COHORT.md) retains 17 independently audited nearby stationary roots and 68 conditional static output enclosures on one small shared parameter cell. The [separate completed energy-stability component](../evidence/cohort_energy_stability_final_delivery.json), supported by its [independent review](../evidence/cohort_energy_stability_independent_review.json), additionally certifies that the selected roots are strict local minima of the same frozen discrete energy throughout that cell. This is local energetic stability in the stated free-coordinate space, not a physical time-path theorem. The outward-rounded discrepancy enclosure remains

\[
5.18638714\le D(\theta)\le5.18638737\;\mu\mathrm m,
\qquad\theta\in C.
\]

The exact rational bounds in its receipt are authoritative. These are bounds for the declared finite-dimensional energy and static observation contract. Neither more complete source-file coverage nor the stronger local-minimum proof identifies the source observation law or makes this discrepancy an experimental uncertainty estimate. Local energy minima do not supply an experimentally validated time path.

Full parameter-domain coverage, branch/loading-path selection, PDE and quadrature/discretization errors, source-aligned temporal observations, and the larger-region prediction remain distinct tasks. Completing this source inventory does not complete the original five scientific gates and does not establish a new biological or mathematical concept.

## Receipts and limits of this chapter

The [binding receipt](../runs/source_contract_supplements_v1/public_source_audit_binding.json) records the exact final audit, draft and original-visual hashes. The previous availability refresh remains a historical bounded search receipt; the new download/inspection inventory supersedes its previously uninspected supplement scope without altering its archived bytes. No original source file, physical model, proof packet, measurement-noise allowance or biological result was changed by preparing this chapter.
