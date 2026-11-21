# What the published movie can measure

**Supplementary Movie 3 verifies a qualitative bending example and supports its correspondence with Fig 4c/d. It does not provide a calibrated, uncensored initial-to-peak apical/basal geometry, or certified plug/void area bounds.** This is an audit of one published example, not a substitute for the 31 larger-region records, a material identified set, or raw-volume validation.

The source is Countryman et al., *Endogenous OptoRhoGEFs reveal biophysical principles of epithelial tissue furrowing*, Nature Communications 16, 7665 (2025), [article](https://pmc.ncbi.nlm.nih.gov/articles/PMC12361397/), [DOI](https://doi.org/10.1038/s41467-025-62483-6). Relevant locators are Results Sec 7, Fig 4c/d/f, Methods Sec 14–15 and 17, and the [Supplementary Movie 3 description](https://pmc-oa-opendata.s3.amazonaws.com/PMC12361397.1/41467_2025_62483_MOESM3_ESM.pdf). The [original movie](https://pmc-oa-opendata.s3.amazonaws.com/PMC12361397.1/41467_2025_62483_MOESM7_ESM.avi) remains unchanged locally. No author contact or model fit was performed.

## Verification record

The replayable [audit script](../runs/movie_geometry/audit_movie.py) independently reads all 26 frames, checks source checksums against official NCBI metadata, records decoded-frame hashes, checks the displayed scale bar at four thresholds, compares published image fields, and verifies that all source hashes are unchanged afterward. The [receipt](../evidence/movie_geometry_audit.json) retains results, source locators, software versions, exact crops and the distinction between automated checks and visual judgments.

| Item | Verified observation | Limit |
|---|---|---|
| Original movie |4,547,582 bytes; SHA256 `260555e647c56a74b157ef624fe263568cc47fbda553d3e8aa7fe0ad05e14be7`; official MD5 agrees | Published composite movie, not microscopy stack |
| Decoding |26 frames, 1908 × 826 pixels, 4 playback fps; ffmpeg decodes without error |4 fps is playback, not acquisition time |
| Biological labels | All 26 burned-in times transcribed: 0, 21, 43, …, 540 seconds; exactly floor(21.6 i) for frame i | Methods gives nominal 22 s/frame; no raw timestamp or exact activation-origin record |
| Horizontal scale | Caption 20 µm; white bar 160 × 15 display pixels at thresholds 128, 200, 230, 250 | Does not establish every panel's physical transformation |
| Transverse crop | Half-open source-pixel rectangle x [20, 576), y [170, 538) | Upper rectangular edge clips lateral membrane structures; it is not an apical contour |
| Source correspondence | Fig 4d 2 min best matches Movie 3 frame 5; grayscale NCC 0.915119, high-pass NCC 0.682942 | Evidence for the same published field; not an acquisition/specimen-ID crosswalk |

The exact timestamp list is in the receipt. Frame 5 displays 01:48, frame 13 displays 04:40, and frame 25 displays 09:00. Nominal figure labels such as 2 min are not precise acquisition timestamps. Fig 4d's nominal 0 min crop best matches movie frame 1 (00:21); this is another reason not to equate a displayed zero with a precisely registered preactivation state.

## Which optical plane is shown?

The Movie 3 description identifies a DV transverse view, a horizontal section 2 µm below the apical surface, and a horizontal section 30 µm below it. It does **not** specify whether the horizontal reference follows the initial apical location, a current global reference, or another registered stack convention. Sec 15 describes correction of tissue motion using PIVlab, but does not supply the actual transformations or a general z-plane tracking rule for this movie. Sec 17 describes following an apical myosin-intensity peak for a separate myosin-spreading analysis; transferring that rule to Movie 3 would be an unsupported inference.

A horizontal optical plane through a bent tissue is not the same object as a surface that remains 2 µm below every local apical point. The latter would generally be curved. A model comparison therefore needs an explicitly declared horizontal-plane observation operator and an unresolved plane reference, rather than silently treating the movie as a local surface-following slice.

## What is calibrated, and what is only consistent?

The displayed bar gives a nominal horizontal scale

\[
s_x=\frac{20}{160}=0.125\;\mu\mathrm m/\text{display pixel}.
\]

Allowing each raster endpoint an assumed one-pixel uncertainty yields a width interval [158, 162] and the conditional scale interval

\[
s_x\in\left[\frac{20}{162},\frac{20}{158}\right]
=[0.12345679,0.12658228]\;\mu\mathrm m/\text{display pixel}.
\]

This is a **digitization sensitivity allowance**, not a confidence interval or a complete physical calibration error bound. The movie container reports square display samples (SAR 1:1), but that does not establish biological x/z isotropy after transverse reslicing. No per-panel resampling, crop or registration matrix is supplied. Rotation and square-pixel reslicing are plausible for the horizontal panels, but remain explicit assumptions.

For consistency only, the transverse crop is 368 pixels tall, giving 46 µm if the same scale applies in z. The Methods 44 µm stack at 2 µm plane spacing would have 23 plane centers; 23 pixel-bin edge widths would span 46 µm. The agreement supports a plausible reslicing convention; it is not an independent calibration of that convention. The 0.12 µm spacing in Sec 14 describes the activation scan and must not silently replace the published display scale.

## Apical censoring prevents a height reconstruction

At movie frames 0 and 13, bright lateral membrane structures intersect the upper border of the transverse image rectangle. A continuous apical tissue boundary is not visible there. The lower bright band is compatible with the basal membrane region, but heterogeneous signal and adjacent structures prevent a validated anatomical contour/error envelope from being assigned by a brightest-pixel rule. The primary Methods themselves used manual basal segmentation because of poor signal.

Consequently, neither tracing the top crop edge nor treating the 32–35 µm protocol height as an observed per-embryo boundary is justified. No peak profile is extracted: an image at an illustrative late time is not independently established as the time of maximum shortening, and the missing apical contour persists.

An elementary censoring calculation makes the limit precise. Let image y increase downwards, let the crop start at y_c, and suppose an independently annotated basal ordinate satisfies b∈[b_0−δ,b_0+δ], while the apical ordinate is a≤y_c. If a calibrated positive z-scale satisfies s_z≥s_min, then

\[
H=s_z(b-a)\ge s_{\min}(b_0-\delta-y_c)
\]

when the right side is positive. No finite upper bound follows from this crop alone. Here neither the basal anatomical error δ nor a physical z-scale interval is independently validated, so the formula does **not** produce a rigorous biological bound in micrometers from the present movie.

More formally, the ideal operator that retains only boundary segments inside the crop has the same output for two distinct smooth apical curves lying entirely above it and one fixed visible basal curve. Their heights differ. This is nonidentifiability for that boundary-cropping operator, not a theorem about every possible microscopy measurement or a new scientific concept. A prior mechanical model could supply an extrapolation, but would no longer be an independent observed geometry.

The [independently generated schematic](../runs/movie_geometry/censoring_nonidentifiability_schematic.png) illustrates this argument. It contains no traced or adapted source frames.

## Can the central plug and dark moat supply areas?

The source defines its bending phenotype using a central tissue plug surrounded by acellular perivitelline space (Sec 7). The horizontal movie slice visibly supports that qualitative example. The published raster does not, by itself, establish solid occupancy from fluorescence intensity: dim cell membranes cross parts of the apparent gap, the central plug boundary is not one sharply resolved continuous contour, and no PSF, raw-volume segmentation, independent annotation, or anatomical error bound accompanies these frames. Cell interiors are dark even when cells are present.

We therefore do **not** report certified inner/outer plug-void areas or fabricate a ±pixel envelope as if it were ground truth. A future explicitly exploratory annotation/segmentation study could report how image proxies vary across declared thresholds, operators and annotators. Such variation is a sensitivity analysis, not a validated anatomical uncertainty interval. Black polygonal display corners must also be masked rather than counted as perivitelline space.

If validated masks become available, let A_i be the plug area and A_o the area enclosed by the outer void boundary; the gap area is A_o−A_i. The area-equivalent radii

\[
r_i=\sqrt{A_i/\pi},\qquad r_o=\sqrt{A_o/\pi}
\]

are invariant under translation and avoid choosing a future-output-derived activation center. They do not imply that a noncircular observed contour is a circle. Converting pixel areas to physical area requires the horizontal transformation determinant; under an assumed shared isotropic scale, areas multiply by s_x² and radii by s_x. The ratio A_i/A_o cancels a shared constant area scale, but not the missing semantic mask validity or plane/time ambiguity.

A horizontal-plane annulus test on an axisymmetric model can thus be a useful **conditional necessary check** once observational masks and plane uncertainty are supplied. It cannot currently be claimed as a scored, certified biological prediction from this movie.

## What is known about specimen identity?

Sec 7 discusses Fig 4c together with Movie 3; the Fig 4 caption explicitly makes Fig 4d the in-plane view of the embryo in Fig 4c and gives rounded activation magnitude 1.02. Image correspondence provides an additional, independent check on this association. Direct resizing of the movie middle panel to the Fig 4d 2 min rectangle produces NCC 0.915119 at frame 5; a high-pass comparison, which suppresses shared low-frequency shading, gives 0.682942. The same high-pass template compared with the movie's distinct 30 µm panel remains below 0.1 in all 26 frames. The exact algorithm and all values are retained; no feature-based warp was optimized against the figure.

The canonical 50 µm numerical table has exactly one dose rounding to 1.02: 1.02397148, row locator `optoGEF2_large_col_F`, source Fig 4h!F8/F9. The same dose appears in Fig 4f with diameter 49.8939753 µm. This is a **unique rounded-dose association**, not a verified biological specimen ID. The movie and worksheets do not provide an acquisition-ID crosswalk. The observed geometry is not used to choose material parameters or refit the model.

## Completion consequence

This audit exhausts the specified public example's support without converting publication-display assumptions into measurements. It improves the source correspondence and makes the observation gaps specific: unclipped apical boundaries, physical z/per-panel transformations, time/plane reference, independently validated masks/errors and a specimen crosswalk. It supplies neither the missing 31-record shape cohort nor a full coupled identified-set prediction. Those requirements remain distinct from the completed provenance/display audit.

Local source-image contact sheets and the downloaded unmodified article figure are retained only as inspection records in `runs/movie_geometry/`. They must not be embedded as adapted figures in the public Space/site. The source license is CC BY-NC-ND 4.0; publication of this project's original schematic does not imply permission to redistribute adapted source images.
