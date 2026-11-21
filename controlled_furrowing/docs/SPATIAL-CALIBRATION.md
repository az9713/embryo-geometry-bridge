# Independent small-region spatial adequacy check

**This check uses the 17 small-region records only.** It reads no larger-region response or bending label, changes no model parameter and introduces no numerical threshold chosen from outcomes. It is a secondary calibration-data adequacy check, not held-out validation.

## Separate observation operators

The nonlinear model was fitted to maximum vertical material-column shortening. Its source counterpart is a maximum over space and time relative to a cellularization reference. The source additionally provides separate apical and basal maxima. Those maxima need not occur at the same location or time, so their sum is not constrained to equal the maximum cell-height loss.

For each stored equilibrium, coordinates increase apically/outward. Within the initially activated material radius, the independent calculations are

\[
A_{\rm pred}=\max_R[-u_{z,a}(R)],\qquad
B_{\rm pred}=\max_R[u_{z,b}(R)].
\]

Values are converted from thickness-one coordinates to micrometres using the locked height 33.5 micrometres. The numerical profiles are static; comparison to separate experimental time maxima is consequently descriptive and conditional. No sum-of-maxima assertion is used.

## What the fitted models actually predict

All 51 stored equilibria—17 records in each of three basal-area-modulus scenarios—are converged and locally orientation-preserving. Yet every scenario predicts a larger basal-up maximum than apical-down maximum in every record. The source maxima show the opposite ordering in 16 of 17 records.

| Basal area modulus / (shear modulus × height) | Height-proxy training RMSE (µm) | Apical-component descriptive RMSE (µm) | Basal-component descriptive RMSE (µm) | Records with predicted apical dominance | Records with source apical dominance |
|---|---:|---:|---:|---:|---:|
| 0.2 | 2.2472 | 5.6414 | 3.7962 | 0 / 17 | 16 / 17 |
| 2 | 2.2558 | 5.5707 | 3.7209 | 0 / 17 | 16 / 17 |
| 20 | 2.2694 | 5.2535 | 3.3550 | 0 / 17 | 16 / 17 |

The strict ordering comparison uses no amplitude cutoff. It agrees with the source order for only one record in each scenario. These records are source-column identifiers, not independently verified biological specimen IDs; no sampling-independence or significance claim is made.

This is a concrete warning that a maximum-height calibration can hide an incorrect partition of apical and basal motion. It does not show that every finite-height model, every basal modulus, or every dynamic observation operator fails. The three scenarios are a finite declared family, not a complete identified set. The results should not be upgraded into a universal negative theorem or a novel biological principle.

## Why area stiffness is not automatically basal normal support

For a flat surface with purely transverse graph displacement \(w(r)\) and no tangential displacement, its area stretch is

\[
J_s=\sqrt{1+(\partial_r w)^2}
=1+\tfrac12(\partial_r w)^2+O((\partial_r w)^4).
\]

Hence a stress-free area-elastic term has leading transverse energy

\[
\frac{K_b}{2}(J_s-1)^2
=\frac{K_b}{8}(\partial_r w)^4+O((\partial_r w)^6).
\]

It has no quadratic normal stiffness in this particular motion. In contrast, a pre-existing tension term gives

\[
T_b(J_s-1)=\frac{T_b}{2}(\partial_r w)^2+O((\partial_r w)^4).
\]

This is standard geometric stiffness, not a new result. It explains why raising a stress-free area modulus is not equivalent to adding prestress, a transverse clamp or a real basal load-bearing network. It does not prove that any one missing ingredient will fix the observed discrepancy. Prestress, initial curvature, tissue surroundings, contact and dynamic mapping remain distinct hypotheses.

## Decision boundary

The declared family has not demonstrated the source small-region apical/basal morphology under its static maximum proxy. Larger-region simulations can still be retained as conditional diagnostics or negative controls, but should not be advertised as an already validated geometric forecast. Any constitutive or observation-model revision must preserve the current results, use only the small-region information for fitting, declare new assumptions and freeze new predictions before larger-region scoring.

The machine-readable receipt `evidence/spatial_calibration_adequacy.json` retains all 51 source-file hashes, record IDs, direct profile-derived values and summary metrics. No larger outcomes were used in this check.
