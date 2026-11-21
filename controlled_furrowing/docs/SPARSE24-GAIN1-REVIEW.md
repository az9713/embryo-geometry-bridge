# Independent review of the completed 24×12 strength-1 state

**31 independent checks pass for one completed activation strength only.** The 24×12 higher-quadrature state at \(f/\mu=1\) is a verified discrete equilibrium with exact positive-J and global embedding certificates. This review does not certify completion of the strength-5 or strength-10 targets, and performs no new optimization.

Receipt: `evidence/sparse24_gain1_independent_review.json`, including the replay program, all exact element certificates, continuation provenance, source hashes and two-mesh comparison. Reviewed field:

```text
runs/sparse_quadratic/24x12_lambda5000/gain_1.0/final.json
SHA256 f1d03340c76d760730e62fe2b961e4751d0791411ae5df15604abe02124ffd3d
```

The folder's word “gain” denotes this runner's activation-strength target; it is not a newly fitted biological gain.

## Model and scope

The model is an idealized flat, finite-height, axisymmetric tissue solid with a localized vertical active-fibre contraction potential. The height scale is 33.5 µm, patch diameter 18.71024074 µm, and active depth 5–25 µm below the apical surface. The shear scale is one, logarithmic-volume penalty \(\kappa/\mu=5000\), basal/apical area moduli are 2 and 0.2, and basal/apical tensions are 2 and zero in the declared normalized units. The outer meridional rim is clamped.

The displacement field is continuous P2 on 24×12 meridional cells: **576 quadratic triangles and 2,375 free displacement degrees of freedom**. The actual integration is Duffy8 volume quadrature and surface Gauss12. The nominal inherited parameter fields `quadrature=4` and `surface_quadrature=8` do not describe the actual overridden integration policy; the separate saved policy does.

No biological response observations, target labels, parameter fits or new material laws enter this review. Maximum material-column shortening is an explicit model observable; it is not automatically apical furrow depth or a transient experimental maximum.

## Saved-state force, energy and geometry replay

The sparse source remains

```text
d4833318151c7b94b4c56d40a96158d55a06258619ef080548d489e2be0394f9
```

and every physical-kernel parent hash matches the original 24×12 protocol. Reconstructed energy and gradient are bit-identical to the non-sparse higher-quadrature parent kernel. The saved total energy is reproduced with zero absolute difference in this replay.

| Check | Independent result |
|---|---:|
| Duffy8/surface12 scaled free-force residual | 2.983e−12 |
| Duffy12/surface16 cross-evaluated residual | 3.216e−12 |
| Required force threshold | 1e−7 |
| Original sampled quadrature \(J\) range | [0.997521, 1.001901] |
| Exact whole-discrete-body \(J\) enclosure, outward rounded | [0.99737057, 1.00242578] |
| Exact global embedding | Pass |

The correlated degree-four rational Bernstein construction combines \(\det F_m R\) before division by reference radius \(r\), retaining correlations that separate factor bounds discard. It proves the stated whole-element bounds for every triangle, including valid axis limits. The stored rational endpoints also prove \(|J-1|<1\%\) everywhere in this particular discrete body. One percent is a numerical diagnostic, not an independently measured biological compressibility limit.

The earlier separate-factor bound [approximately 0.96256, 1.03824] was conservative; it was not evidence that the body actually attained 3.8% local distortion. The sharper interval certifies the same saved field without changing its solution. The construction and global-injectivity arguments are explained in `QUADRATIC-INDEPENDENT-REVIEW.md` and `MIXED-GEOMETRY-PROOF.md`. These are discrete geometry certificates, not a continuum approximation-error bound.

## Three attempts, one rejection, and no rejected-state reuse

The complete target-1 ledger is replayed, including the rejected attempt:

| Attempt | Last verified strength → attempted strength | Decision | Replayed force residual | State used next |
|---:|---|---|---:|---|
| 0 | 0 → 1 | Reject | approximately 1.550e−3 | Verified strength-0 vector retained |
| 1 | 0 → 0.5 | Accept | 3.143e−12 | Verified strength-0.5 vector |
| 2 | 0.5 → 1 | Accept | 2.983e−12 | Final verified strength-1 vector |

Every attempt's initial-vector hash matches the most recent accepted vector. The final vector matches the last accepted attempt. The first rejected vector is retained in its own file and never admitted as the next initial state. Its geometry was positive, but its force residual failed; “rejected” here describes numerical nonconvergence, not physical inversion.

Accepted-path work and relaxation drops are independently recomputed using the chosen mechanical potential. Their energy-balance residual is zero to the recorded arithmetic in this replay. The ledger is numerical continuation bookkeeping; it is not a biological time trajectory, calorimetric measurement or ATP account.

## Matched comparison with the verified 12×6 higher-rule state

The comparator is `runs/activated_quadrature/lambda_5000.0_strength_1.0/reequilibrated.json`, bound to its prior independent-verification hash. Its physical parameters, strength and actual integration policy match; only the displacement mesh differs. The reviewer does **not** substitute a 16×8 default-rule state. No verified matched 16×8 higher-rule strength-1 field is included, so no three-level convergence trend is claimed.

| Quantity | Verified higher-rule 12×6 → 24×12 comparison |
|---|---:|
| Maximum material-column shortening | 1.99113690 → 2.33011795 µm |
| Shortening change | +0.33898105 µm |
| Eulerian apical vertical-discrepancy sup | [0.35326005, 0.35848985] µm |
| Eulerian basal vertical-discrepancy sup | [0.01403053, 0.01852789] µm |
| Same-material apical vector-displacement sup enclosure | approximately [0.35297418, 0.35542345] µm |

The Eulerian intervals are exact rational enclosures displayed with outward rounding. They compare heights at the same **deformed physical radius**, using quadratic radial inversion and adaptive interval cells. Their requested reporting width 0.01 µm is achieved for both surfaces. It is not a scientific accuracy threshold.

The same-material metric instead pairs equal reference radii and bounds a vector displacement; its upper bound also bounds Hausdorff distance. Its lower bound is not a Hausdorff lower bound. These different metrics must not be substituted for one another. All are comparisons of two stored discrete curves, not uncertainty intervals for a continuum solution.

## What this establishes

One previously completed numerical target has independently verified source provenance, discrete stationarity, integration consistency, exact geometry and last-verified continuation behavior. A substantial two-mesh observable difference remains. The check does not establish mesh convergence, uniform penalty robustness, a unique equilibrium branch, a global energy minimum, biological accuracy, a compatible parameter set, novelty, or completion of the remaining 24×12 targets.

The live continuation to other strengths was neither restarted nor duplicated. Its outcomes require their own receipts and review.
