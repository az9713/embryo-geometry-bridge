# Exact active force bounds for the K–L loading trial

This producer supplies the dose-dependent part of the negative K–L continuation test. Both stored models use the same normalized finite-element energy, binary reference, free-coordinate injection, 24 reference activation flags, and material axis \(e_Z\). The only load dependence is the exact common gain \(g\) times the requested optical dose. The program reads both frozen proof centers from their binary packets, verifies exact \(X+Pu\) at both anchors, and uses K’s actual recentered \(u_K\) and binary preconditioner \(R\). The direction \(v\) is the exact rational secant from those centers and the exact dose gap. [Prior physical crosswalk](../exact_loading_tube_design_v1/receipt.json) · [Frozen secant plan](../exact_loading_tube_pilot_preparation_v1/plan.json) · [Independent affine-law review](../transfer_ambiguity_nonlinear_probe_v1/certification_design/exact_targets_v1/loading_affine_identity_v1/independent_readback/review_receipt.json)

## Rational element coefficients

On an active reference triangle take \(\lambda_0,\lambda_1,\lambda_2\) as its barycentric coordinates and \(D_Z\) as its axial material derivative. Positive orientation is verified using the exact binary vertices and their exact determinant. In the P2 vertex shape \(N_i=\lambda_i(2\lambda_i-1)\), the derivative is \(D_Z\lambda_i(4\lambda_i-1)\). For an edge shape \(N_{ij}=4\lambda_i\lambda_j\), it is \(4[\lambda_iD_Z\lambda_j+\lambda_jD_Z\lambda_i]\). Both are homogeneous linear forms in \(\lambda\) on the simplex. With the cylindrical reference radius \(R_{\rm ref}=\sum_k r_k\lambda_k\), the axial Gram matrix is

\[
G_{ij}=\int_T R_{\rm ref}\,(D_ZN_i)(D_ZN_j)\,dR_{\rm ref}\,dZ.
\]

For each product of three barycentric coordinates the unit-simplex integral is \(n_0!n_1!n_2!/5!\), where \(n_0+n_1+n_2=3\). Multiply by the exact positive triangle determinant. This formula also applies to the three active axis-vertex and three active axis-edge elements. It has no denominator containing the radial coordinate. The same six-by-six matrix is applied to radial and axial displacement components independently. Scatter \(G\) on each component’s exact local-to-free map to obtain \(A_0\), and scatter its product with the **full** stored binary reference coordinates to obtain \(b_0\). A fixed coordinate contributes to the reference offset even though it has no free-displacement column. The exact spatial identity is

\[
\partial_I F(u,I)=\alpha(A_0u+b_0),\qquad
\partial_I H(u,I)=\alpha A_0,\qquad \alpha=2\pi g/S_{\rm binary}>0.
\]

The force is the energy gradient on free coordinates. Frozen quadrature weights do not enter these coefficients. The stored midpoint values enter the offset and the authoritative coordinates. The reference triangle determinant comes from its three vertex values. [Exact pointwise target](../../scripts/exact_fe_target_safe.py) · [Independent active moment precedent](../transfer_ambiguity_nonlinear_probe_v1/certification_design/exact_targets_v1/negative_KL_active_energy_v1/independent_actual/replay.py)

## Signed assembly and the norm

The program retains exact per-element Grams, the sparse upper triangle of \(A_0\), and \(b_0\). It writes the exact rational products and absolute row sums with the common factor still removed. Each row of the actual binary \(R\) has a common dyadic denominator. One rational common denominator clears the sparse entries of \(A_0,b_0\). Integer operations therefore form each row of \(RA_0\) and \(Rb_0\) exactly. Absolute values are taken after cancellation across the signed columns of \(RA_0\). This computes the sparse product’s actual infinity norm before multiplying by the upper directed bound for \(\alpha\). The same directed positive \(\alpha\) enclosure at Arb640 scales all signed vectors, including negative values. The delivered quantities are

\[
a \supseteq R(Cu_K+d),\quad k\supseteq RCv,\quad M\ge\|RC\|_\infty,\qquad C=\alpha A_0,\quad d=\alpha b_0.
\]

The whole K–L trial still needs a fresh envelope-domain proof and signed central-dose derivative bounds over that envelope. They remain separate from this state-independent active result. [Continuation consumer and its inputs](../exact_loading_tube_pilot_preparation_v1/PILOT.md)

## Separate short extension on the original K cube

The original independent exact K proof supplies \(\eta_K,Z_K\) and a strict margin \(m=(1-Z_K)r_K-\eta_K>0\) at its unchanged original geometrically certified cube of radius \(r_K=2^{-34}\). A constant predictor permits reusing that identical **state domain** and nominal K derivative bound if the exact dose drift is added. With \(B\ge\|R(Cu_K+d)\|_\infty\) and the computed \(M\), let \(D=B+r_KM\). Choose \(\ell=\min((I_L-I_K)/2,m/(2D))\) for \(D>0\) and half the gap for zero drift. The fixed-cube residual bound becomes \(\eta_K+\ell B\). The derivative bound becomes \(Z_K+\ell M\). The new self-map has strict margin at least \(m/2\), and K’s independently verified same-\(R\) positive symmetric part yields strict local minima on this short branch. That argument identifies K at the left endpoint by uniqueness on its original cube. It can be evaluated and independently replayed from the same active arithmetic. It cannot reach L by construction, and it imports no unproved derivative bound to the wide rectangular K–L envelope. No interval length or branch is claimed by the cold results. [Independent exact K proof](../transfer_ambiguity_nonlinear_probe_v1/certification_design/exact_targets_v1/remaining32_bridge_v1/first_wave_reduction_v1/independent_actual_preparation/actual_replay/max_negative_contrast_change_small_16.json) · [Same binary preconditioner](../transfer_ambiguity_nonlinear_probe_v1/certification_design/exact_targets_v1/remaining32_preconditioner_spd_v1/independent_actual_preparation/actual_replay/max_negative_contrast_change_small_16.json)

## Execute after independent review

From the study directory, a separate owner can rerun the cold checks and then invoke the actual stage only after the source phase admits it.

```powershell
python -B runs/exact_loading_active_operator_preparation_v1/active.py verify-preparation
python -B runs/exact_loading_active_operator_preparation_v1/active.py run --review-sha <INDEPENDENT_PRELAUNCH_SHA256>
```

The write-once `independent_prelaunch.json` must carry status `PASS_INDEPENDENT_KL_ACTIVE_OPERATOR_PRELAUNCH`, `pilot_id` from the protocol, current three `producer_bindings` digests printed by `verify-preparation`, the five `pilot_bindings` digests in the protocol, literal `true` for each `required_review_checks` entry, and `active_operator_only_authorized: true` with `passive_or_geometry_jobs_authorized: false`. The actual command creates `actual/start.json`, `exact_sparse_operator.json`, `result.json`, and `terminal.json`, with one process and no children. A runtime error retains failure and terminal. A successful run emits a producer *candidate* and requires a genuine outer exit observation.

An independent actual replay must reassemble the moments and projection and verify exact signed rows and the active norm with the canonical same binary \(R\), along with the outward common factor and, if reported, original K cube reuse. Its normalized result can carry the required pilot operand status `PASS_INDEPENDENT_KL_EXACT_ACTIVE_OPERATOR_BOUNDS` and authenticate the exact numeric fields `signed_R_active_center_rows_exact`, `signed_R_C_v_rows_exact`, and `R_C_norm_upper_exact` under the same five pilot bindings and identity fields. The pilot consumer additionally requires the independent result’s terminal `KL_LOADING_OPERAND_STAGE_FINISHED` and a separately root-observed zero exit `PASS_ROOT_OBSERVED_KL_LOADING_OPERAND_EXIT`. The producer candidate does not supply this independent proof. [Phase protocol](protocol.json) · [Cold tests](preflight.json)
