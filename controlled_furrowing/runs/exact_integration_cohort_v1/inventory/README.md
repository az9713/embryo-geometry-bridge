# Exact-integration cohort inventory

The 17 cases share the exact five-parameter cell with SHA-256 `64f93fefa034f4c004b93451d98b2fd04e0072855be8c5e3bba9babb8f394dbf`. All 17 have independently reviewed cached-quadrature roots, whole-body geometry on their original proof boxes, and strictly positive cached energy curvature. The verified exact-spatial-integral FE root covers case **B**. The other 16 cases need exact-integral bounds and root proofs.

The [machine-readable inventory](inventory.json) records exact doses, corrector radii, base parameters, proof context and artifact hashes for each full record ID. It includes every selected predictor, embedded point cache and proof-center path, plus old root, observer and energy receipts with their independent reviews. The [read-only checker](inspect.py) verified 209 input files, including all 147 bindings in the frozen preparation protocol. It uses NumPy to compare geometric arrays and read the 611×5 predictors. The 611×611 preconditioners are bound through the embedded cache hashes, with no matrix arithmetic. There were no new integrations, solves or proof replays.

## Old proof domains

The table gives rounded displays. Consult the exact rational values in the inventory for arithmetic. The contraction bound `Z_Q`, quadrature Jacobian bound `J_Q`, and energy Hessian bound `m_Q` apply to the selected cached-root tube. The field `J` and strong monotonicity bound `m` apply throughout the old observer nodal box, whose receipt supplies per-element bounds and boundary observer transfers. A new center and common box require fresh enclosures or verified containment.

| Case | Dose | Corrector radius | `Z_Q` upper | `J_Q` lower | Old box `J` | Old box `m` lower | `m_Q` lower |
|---|---:|---:|---:|---:|---|---:|---:|
| B | 0.14801043 | 1e-11 | 0.048682 | 0.997472 | [0.997288, 1.003034] | 0.480755 | 0.0005639 |
| C | 0.32718904 | 1e-11 | 0.064980 | 0.993328 | [0.992678, 1.005797] | 0.456807 | 0.0005689 |
| D | 0.69240310 | 1e-11 | 0.091006 | 0.980751 | [0.979447, 1.011614] | 0.410684 | 0.0005749 |
| E | 0.97607953 | 1e-11 | 0.109292 | 0.971690 | [0.969818, 1.016509] | 0.360826 | 0.0005713 |
| F | 0.46402059 | 1e-11 | 0.075408 | 0.988969 | [0.988202, 1.007429] | 0.439831 | 0.0005727 |
| G | 1.04731044 | 1.1076e-11 | 0.118551 | 0.969745 | [0.967754, 1.017618] | 0.348422 | 0.0005665 |
| H | 0.82160217 | 1e-11 | 0.099409 | 0.976394 | [0.974814, 1.013928] | 0.388031 | 0.0005739 |
| I | 0.47133397 | 1e-11 | 0.075921 | 0.988704 | [0.987919, 1.007504] | 0.438905 | 0.0005729 |
| J | 0.98739488 | 1e-11 | 0.110014 | 0.971690 | [0.969479, 1.016689] | 0.358847 | 0.0005710 |
| K | 1.28898317 | 1e-11 | 0.128928 | 0.964299 | [0.961983, 1.020955] | 0.307673 | 0.0005577 |
| L | 1.51057344 | 6.0117e-11 | 0.383148 | 0.960805 | [0.958291, 1.023404] | 0.272763 | 0.0003949 |
| M | 0.78208913 | 1e-11 | 0.096862 | 0.977691 | [0.976193, 1.013233] | 0.394997 | 0.0005743 |
| N | 0.71523649 | 1e-11 | 0.092510 | 0.979958 | [0.978604, 1.012031] | 0.406714 | 0.0005748 |
| O | 0.27583779 | 1e-11 | 0.060655 | 0.994727 | [0.994213, 1.005069] | 0.463513 | 0.0005673 |
| P | 0.80580697 | 1e-11 | 0.098393 | 0.976908 | [0.975360, 1.013652] | 0.390818 | 0.0005740 |
| Q | 0.75547653 | 1e-11 | 0.095138 | 0.978583 | [0.977141, 1.012759] | 0.399675 | 0.0005745 |
| R | 0.85073705 | 1e-11 | 0.101282 | 0.975460 | [0.973821, 1.014432] | 0.382888 | 0.0005735 |

All cases use the same reference nodes, connectivity, mask, free map, normalization, and cached quadrature arrays. Each has 325 P2 nodes, 144 triangles, 611 free coordinates among 650 interleaved radial and axial coordinates, and two sets of 12 surface edges. Zero-based axis-edge elements are 72–77, with 74–76 active. Axis-vertex elements are 0–5, with 2–4 active. The 132 off-axis elements contain 18 active elements. The binary cache lists all these indices and the exact common energy normalization. Its displacement array varies with the case.

All 17 selected proof centers and predictors matched their authoritative receipts. F and K use auxiliary recentered proof states in their corresponding `cohort_recentered_roots_v1/records` folders. Their cache matrix-archive pointers lead to original seed tangents. For the exact-target center, use the embedded cache's `u`, the linked auxiliary state, and the hash-bound predictor. Their `R` and `S` survived the earlier cached-root recentering unchanged. M uses the independently reviewed `actual_M` point cache and predictor through the cohort wrapper. G and L use supplemented radii. The other 15 use `1e-11`.

## Pilot selection

**L** is the strongest available challenge. It has the highest dose, widest old radius, largest old whole-box Jacobian deviation from one, largest cached `Z_Q`, and lowest old box monotonicity and cached energy curvature bounds. The old whole-box rho is about 0.04171 for L and 0.00303 for B. This ranking uses the existing bounds. Difficulty of the exact-integral defect at L remains an inference until evaluated.

For L, use the hash-bound `cohort_uniform_roots_v2/records/10_optoGEF2_small_col_L/result.json`, `cohort_uniform_roots_v1/records/10_optoGEF2_small_col_L/preparation_attempts/000` point cache and predictor, and `material_cohort_direct_v1/direct/16/final.json` center. This stresses both the deformation bounds and the cached contraction margin on the same mesh and parameter cell.
