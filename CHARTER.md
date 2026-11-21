# P4 Research Charter — version 1.1
**Date:** 18 November 2026. **Step:** 1. **Status:** PASS — independently reviewed as a research contract; data, model, novelty and biological claims remain unverified. This is not a novel result or biological validation.

## Scope and question
Primary system: quail early gastrulation in Caldarelli et al., Nature 2024, doi:10.1038/s41586-024-07934-8. Primary baseline: intact marginal/tissue flow; first perturbation: a documented cut/reattachment condition. Chick models are separate comparisons, not pooled training data.

We will test whether available observations discriminate specified mechanochemical closures and whether an identifiable closure predicts withheld perturbation responses beyond the published baseline. Unpaired observations restrict conclusions to ensemble prediction; equivalent mechanisms lead to non-identifiability conclusions.

## Prediction boundary and admissible inputs
Step 4 must lock numerical times t0 and T from source acquisition/staging records. Inputs are measured initial geometry, molecular/motor fields legitimately observed at or before t0, initial asymmetry, documented boundary conditions, and specified intervention schedules. Future velocity, final streak shape, and later marker patterns cannot be used to reconstruct initial inputs. Any baseline using supplied future fields is labeled conditional/descriptive and receives an explicit input-information comparison.

Published outcomes are known; all tests on them are retrospective. Split by biological embryo and intervention family, never adjacent frames. Prospective evidence requires genuinely untouched/new data.

## Observation operator and dependence
Primary measurements are image intensities, pixel/time calibration, documented interventions, and independently acquired molecular/force observations where present. Velocity, strain rate, trajectories, and deformation reconstructed from the same movie are dependent summaries, not independent corroboration. Scoring preserves their dependence and does not count them as four independent datasets.

Fluorescence is not motor activity, active stress, or absolute force. The observation model records reporter kinetics, background, calibration, and inferred quantities. Force/viscosity ratios may be identifiable when absolute force is not.

Paired marker and motion data may support within-embryo joint prediction. Unpaired endpoint staining supports only separately labeled ensemble-level comparisons with justified staging, registration, and between-embryo variation. It is not a time-series measurement.

## Candidate model classes
Published contractility/strain-rate feedback baseline; local tension feedback; accumulated-deformation/history feedback; and appropriately constrained molecular-to-mechanical coupling. These are candidate classes to audit, not asserted novel or biological facts.

Step 10 must determine distinguishability under the chosen rheology, observation operator, parameter symmetries, and interventions before a fit is interpreted as a mechanism identification. Statistical preference alone cannot distinguish exactly equivalent closures.

## Falsification and inconclusiveness table
| Target | Required observations | Prediction/comparator | Failure or inconclusiveness | Tolerance owner |
| --- | --- | --- | --- | --- |
| Published baseline reproduction | Reference arrays/profiles or documented figure extraction | Intact and one perturbation fields, profiles and amplitude | Unresolved discrepancies outside justified error; unavailable quantitative reference limits claim | Step 5, locked before pass |
| Measurement recovery | Calibrated images and synthetic known transformations | Known motion, deformation and marker fields | Error outside declared operating range; unknown units or missing planes | Step 6 |
| Geometry prediction | Independent withheld embryos and documented inputs | Spatial profiles, amplitudes, motion and streak geometry vs fixed baseline | No improvement or predictive mismatch; low-power data may be inconclusive | Steps 4/6/11 before Step 13 scoring |
| Molecular–mechanical joint response | Paired measurements, or justified matched ensembles | Marker and motion responses under perturbation | No pairing or unsupported temporal mapping; predictions outside uncertainty | Steps 4/11 |
| Closure discrimination | Constrained rheology and informative interventions | Distinct prediction sets after nuisance uncertainty | Equivalent/overlapping predictions or unsupported force scale → non-identifiable | Step 10 |
| 3D biology | Actual surface/thickness/folding observations | Curvature and shape from supported material/active inputs | Missing observations or imposed desired curvature → no 3D validation | Step 15 |

No numerical threshold is selected to favor an observed result. Later gate versions retain original criteria, reasons for revision, and affected-result invalidation.

## Approximations and exclusions
Continuum scale separation, low inertia, rheology, density balance, and planar validity are hypotheses to check in Step 7. The initial body-axis asymmetry is supplied, not explained. Planar tissue maps do not derive layer topology changes. No complete embryo, genome-to-shape derivation, quantum-chemistry derivation, or new universal principle is claimed.

## Dependencies and completion
Complete/verify this charter before Step 2 workspace creation. Steps 3 originality and 4 data audit may then proceed independently under this scope. Step 5 reproduction follows baseline/data feasibility classification. Formal step verification and scientific outcomes are separate. A negative audit can be verified; it does not establish a successful biological bridge. Missing molecular evidence limits claims; absent 3D validation leaves Step 15 and the complete 15-step project incomplete.

## Step 1 evidence and review
The initial internal AI review returned FAIL pending seven clarifications: prediction boundary, correlated observables, closure equivalence, paired vs unpaired observations, falsification criteria, observation-to-force mapping, and full-project completion limits. Version 1.1 explicitly incorporates each. Internal AI review is not external expert endorsement.

