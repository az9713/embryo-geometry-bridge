# Cold test of the negative K–L loading order

The two inputs are the fixed negative-endpoint cases small_16 and small_17, corresponding to source records K and L. The cold crosswalk verifies byte equality of the reference coordinates, P2 triangle and surface connectivity, free-coordinate injection, reference activation mask and normalization. The declared material coefficients and exact common gain agree. Each case retains its own exact dose, proof center, positive posterior radius bound and independent exact-spatial root authority. The recentered K proof uses its recentered coordinates. Its archived-observer center is irrelevant to the energy inequality.

For \(g>0\), consider the normalized energy \(E(u,I)=E_{\rm pass}(u)+gI V(u)\). Along a differentiable same-model stationary branch, \(H u'=-g\nabla V\). Thus \(dV/dI=-u'^{T}Hu'/g\le0\) wherever \(H\succeq0\). The test asks only whether the two certified roots obey this necessary ordering at their distinct increasing doses. A positive lower endpoint of \(V_L-V_K\) rules out that smooth stable connection. A nonpositive upper endpoint verifies the necessary condition alone. An interval containing positive and nonpositive numbers leaves the ordering unresolved.

Compute the scalar \(V/(2\pi/S_{\rm binary})\). The factor removed is identical and positive for K and L. For each activated element, the two axial derivatives of the P2 field are affine in material coordinates. Interpolate their corner values by barycentric coordinates and expand
\[
 \frac12\int R_{\rm ref}\!\left|\partial_Z Y\right|^2
 =\frac12\sum_{a,b,c=0}^{2} R_c\, d_a\cdot d_b
          \int\lambda_a\lambda_b\lambda_c.
\]
The moments are exact, nonnegative rationals. No division by reference radius occurs. Corner derivative intervals use the independently certified infinity-norm posterior about the exact binary \(X+Pu_c\). Fixed rim and radial-axis coordinates have zero error. Exact interval squares and products enclose every choice of the shared nodal inputs. Treating dependent intervals separately may enlarge the result but preserves inclusion. Neither a global matrix nor an AD evaluation is needed.

The source, protocol and preflight are frozen before any actual root-energy evaluation. The unit-triangle fixtures have known results \(13/12\) for \(Y=(R+2Z,3Z)\) and \(7/120\) for \(Y=(Z^2,RZ)\). They also exercise signed interval products, reversed orientation and all three possible ordering classifications. The cold preparation read the existing proof, case-contract, box and archive inputs only. It records their 22 current SHA-256 digests.

After an independent cold review, save its exact digest-bound receipt at `independent_prelaunch_review.json`. The required status is `PASS_INDEPENDENT_COLD_NEGATIVE_KL_ACTIVE_ENERGY_TEST`, with `source_sha256`, `protocol_sha256`, `preflight_sha256`, and a nonempty dictionary of literal-true `checks`. Invoke the `run` subcommand with `--protocol-sha256`, `--preflight-sha256` and `--review-sha256` only under root authorization. The output is a write-once candidate in `actual` and requires actual OS observation and independent arithmetic readback before scientific promotion.

[Source-level dose affinity and loading proof](../loading_affine_identity_v1/READBACK.md)
