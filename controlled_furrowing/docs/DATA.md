# Data gate: controlled epithelial furrowing

**Verdict: GO for a retrospective small-to-large transfer test of maximum cell shortening.** The public source workbook also permits a larger-diameter bending-label check. It does not supply a complete individual-resolved 3D shape/time-series benchmark or independent material calibration. This is a data audit, not a model fit, predictive score, or claim of scientific novelty.

The source is Countryman et al., *Endogenous OptoRhoGEFs reveal biophysical principles of epithelial tissue furrowing*, Nature Communications 16, 7665 (18 August 2025), [DOI](https://doi.org/10.1038/s41467-025-62483-6), [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC12361397/). The inspected experiments activate lateral membranes in the dorsal primary epithelium of Drosophila at gastrulation onset. Methods specify cell heights of 32–35 µm, a cylindrical activation volume spanning 5–25 µm below the apical surface, and nominal 18-µm small regions. These are protocol-level descriptors, not per-embryo measured initial geometries in the workbook.

## Actual retrieval and provenance

Direct Springer CDN addresses failed DNS; PMC attachment requests returned HTML challenges despite HTTP 200; EuropePMC supplementary-archive requests timed out or returned HTTP 503. Those failures are retained, but they did **not** establish data unavailability. The official public NCBI S3 listing succeeded:

`https://pmc-oa-opendata.s3.amazonaws.com/?list-type=2&prefix=PMC12361397&max-keys=50`

It contains 26 objects: 6 figure images, 10 supplementary movies, 4 supplementary PDFs, 2 workbooks, and article XML/PDF/text/metadata. Object names, sizes, ETags and listing bytes are retained. Only selected relevant objects were downloaded; listing an object is distinguished from decoding it. The original full-volume microscopy is author-request according to the article; no author contact was made. Published numerical source data and an illustrative movie are nevertheless accessible and inspected.

Downloaded numerical sources:

| Original | Meaning | Bytes | SHA-256 |
|---|---|---:|---|
| [MOESM16](https://pmc-oa-opendata.s3.amazonaws.com/PMC12361397.1/41467_2025_62483_MOESM16_ESM.xlsx) | Source Data, 36 worksheets | 256659 | `11659c641ace017f2848e2418bd30a5b762b65ffb509b96cda08243ebc0f210d` |
| [MOESM4](https://pmc-oa-opendata.s3.amazonaws.com/PMC12361397.1/41467_2025_62483_MOESM4_ESM.xlsx) | Supplementary Data1, 4 primer-design worksheets | 11438 | `6c97e2046808a9e02d0958eff947f983d51437caaae361a31b83bfa8fdf3ddd7` |

The supplement PDF (MOESM1), file-description PDF (MOESM3), reporting summary (MOESM15), official object metadata, and Supplementary Movie3 (MOESM7) were also downloaded. The first is available to the originality/model audit. Sources remain byte-for-byte originals in `data/raw/`; derived extracts are separate. [Download receipts](../evidence/data_downloads.json), [additional downloads](../evidence/data_additional_downloads.json), [object inventory](../evidence/data_oa_object_inventory.json), and [integrity checks](../evidence/data_integrity.json) give sizes, hashes and original-source checksum comparisons.

The article and official metadata state **CC BY-NC-ND 4.0**. No permissive code/data license is inferred from “open access”; no custom code is provided. Originals retain attribution and licensing. This audit's numerical extractions are for reproducible analysis; the article's license is not represented as permission to redistribute adapted article figures. The movie inspection image is a local diagnostic, not an original project visualization or an independently acquired observation.

## All worksheets inspected and independently parsed

Every nonempty cell in all **40 worksheets** was compared between openpyxl and a separate standard-library ZIP/XML/shared-string parser. Cell coordinates, values and full worksheet inventories are retained in `data/derived/all_workbook_cells.json` and [the audit receipt](../evidence/data_workbook_audit.json). Formatting can inflate maximum worksheet dimensions; data counts use actual numeric cells.

The 36 Source Data sheets cover:

| Sheets | What they contain | Use/limitation |
|---|---|---|
| Fig2b | Three illustrated time/apical/basal position traces | 63 timepoints in total; not whole-cohort trajectories |
| Fig2c,2d,2e | Dose/average shortening, maximum shortening, separate apical/basal maxima | 17 OptoGEF2 and 11 OptoCysts small-region observations |
| Fig2f,2g | Published scaling sweep and expression/activity ratios | Derived summaries, not new embryos |
| Fig3b | One-copy/two-copy dose-response summaries | Different genotype comparison; not silently pooled |
| Fig3c,3d,3e,3f | Furrowing/recovery times versus dose/depth | Summaries with response-based eligibility; not complete traces |
| Fig4f | Diameter, activation magnitude, bending Yes/No | 48 source records across six diameters, including the reused 17 small observations |
| Fig4h | Dose and maximum shortening in small/large regions | Same17 small, plus9 nominal50-µm records |
| Fig6b,6d,6f | Area/opening-width trajectories | Other developmental tissues/activation protocols; not this calibration/test cohort |
| SupplementaryFig1g,1h,1i | Expression levels/variation and cellularization rates | No per-record initial-state map into Fig4h |
| SupplementaryFig2a-c,2e,2g,2h,2i,2j | Viability and myosin intensity/spreading/timing | Related tool characterization, not large-region shape coordinates |
| SupplementaryFig3a | Protein-expression comparison | Molecular abundance summary |
| SupplementaryFig7a,7b,7c,7d,7e,7g,7i | Alternative dose summaries/control responses | Some duplicate main-response observations; not independent validation cohorts |
| SupplementaryFig8f | Groove/crest cell-height comparison | Five values each; no full surface/curvature map |
| SupplementaryFig9c,9d | Ablation/retraction curves and initial distances | Basal-mechanics perturbation information, not an absolute modulus calibration |

Supplementary Data1's four sheets are named `RhoGEF2-tgRFPt-SspB`, `RhoGEF2-SspB`, `SspB-tgRFPt-Cysts`, `SspB-Cysts` and contain cloning primers/templates. They do not supply missing geometric trajectories.

## Endpoint and duplicate identity checks

The small-region calibration endpoint must be **maximum shortening**, because that is the large-region endpoint in Fig4h. Region-averaged shortening from Fig2c is a different observable. The script verifies exact dose and maximum-shortening equality between Fig2/Fig4h small records; these are **one set of17**, not34 independent samples. It also verifies equality of the repeated mean-shortening array and matching dose sets in Fig4f.

Each canonical CSV record has an explicit source-cell locator such as `Fig 2d!B4`. Key worksheets lack biological embryo IDs. IDs such as `optoGEF2_small_col_B` are **analysis locators**, not invented specimen identifiers. Same-column endpoint associations follow the source-table organization; exact duplicate dose/response checks support reuse across relevant panels. They do not reconstruct missing raw-file IDs or permit arbitrary joins to other sheets.

The primary Methods statistical section explicitly reports distinct samples, and its OptoGEF2 response discussion identifies 17 embryos. Missing table IDs therefore limit record crosswalks; they are **not evidence that the embryos are dependent samples**. Distinct biological units do not by themselves establish an identical additive-error law, independent dose-normalization errors, or a confidence-coverage guarantee. The new [calibration uncertainty analysis](CALIBRATION-UNCERTAINTY.html) retains the exact source locators and this distinction.

Separate maximum apical and basal displacements may occur at different sites/times. Their sum need not equal maximum height shortening. No such identity is imposed or used to repair data.

## Numerical diameter audit and feasible split

Fig4f stores the following diameters:

| Reported diameter (µm) | Records | Proposed role |
|---:|---:|---|
| 18.71024074 | 17 | Small calibration set, reused in main dose-response panels |
| 24.94698765 | 6 | Larger-region bending evaluation |
| 31.18373456 | 5 | Larger-region bending evaluation |
| 37.42048148 | 5 | Larger-region bending evaluation |
| 43.65722839 | 6 | Larger-region bending evaluation |
| 49.89397530 | 9 | Large-region maximum-shortening and bending evaluation |

The paper labels small/large groups nominally18/50 µm. Both nominal and workbook diameters are preserved; no silent rescaling or rounding is performed. The article discusses experiments up to75 µm, but this numerical source table ends at49.8939753 µm. A quantitative75-µm cohort must not be fabricated from that textual statement.

Canonical exports, with no fitting or scoring:

* `data/derived/calibration_small_18um.csv`:17 OptoGEF2 dose/maximum-shortening records, plus average/apical/basal calibration outcomes and source locators.
* `evaluation_large_50um.csv`:9 dose/maximum-shortening records; held-out responses are separated from calibration files.
* `evaluation_larger_bending.csv`:31 larger-region diameter/dose/bending records.
* `all_bending_records.csv`:48 records including the17 reused small observations, explicitly tagged.
* `example_apical_basal_trajectories.csv`:63 points from three illustrative traces; no unverified join to a cohort embryo ID.

These supply a defensible **retrospective withheld-diameter computation**. Public target values were inspected during data auditing; this is not a prospective or analyst-blinded experiment. The predictive model and selection rules must be frozen before test-score inspection. A data extractor's successful parse is not evidence of predictive accuracy.

## Inputs, outcomes, and what cannot be inferred

**Allowed test inputs:** intervention diameter, OptoGEF2 genotype/protocol and measured activation magnitude. The latter is normalized fluorescence during the first activation scan, not prescribed laser power, an absolute GEF concentration, or a known mechanical stress. The source standard curve is averaged across experiments; raw fluorescence is unavailable to reconstruct a training-only normalization. Consequently the honest target is prediction **conditional on the published realized-dose descriptor**, not prediction solely from pre-exposure laser commands.

**Held-out outcomes:** maximum shortening and bending classification. The source shortening measure subtracts a cellularization-reference trajectory inferred from non-activated tissue. This is an outcome-processing definition, not an admissible future trajectory to feed into a predictor. Calibration may use small-region outcomes; the large-region outcome/reference curve must not initialize a simulation.

Initial cell height is protocol-level32–35 µm; individual heights, full initial surfaces, embryo curvature, force/stress calibration and individual material moduli are not supplied in these tables. The 3 small-region example traces do not yield full-cohort dynamics. Bending is an image-based categorical label from a continuous phenotype, not a numerical curvature, sign-calibrated furrow-depth field, or 3D surface.

Supplementary Movie3 is an accessible50-µm example with transverse and two in-plane panels. It decodes as26 frames,1908×826 RGB,4 playback frames/s. The inspected frame carries biological time and a20-µm scale bar; playback fps is **not acquisition timing**. The publication includes crop/registration/display operations and no verified workbook-record crosswalk. It is useful qualitative geometry evidence, not a substitute for a complete paired raw-volume cohort. Its source hash and first/middle/last decoded hashes are in [the movie receipt](../evidence/data_movie3_inspection.json).

## Gate conclusion and replay

Run `python controlled_furrowing/data/audit_workbooks.py` from the project root after retrieving the unchanged workbooks. The script compares every cell independently, proves duplicate checks, and exports the separated tables. It performs no optimization, learning, parameter fitting, or model evaluation.

**Gate2 passes for actual-data inventory and the specified scalar/label test.** A successful maximum-shortening transfer would be meaningful for that endpoint. It would not alone establish that one mechanical model predicts full larger-region deformation, bending geometry or dynamics. Those stronger claims need additional identifiable mechanical parameters and shape/time-resolved observations; limited public raw access does not invalidate the available numerical comparison.
