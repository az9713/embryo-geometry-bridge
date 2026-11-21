# Verified inputs for 32 spatial integrations

The remaining 32 fixed-input states now have independently replayed whole-box geometry and cached derivative bounds on the same closed free-displacement radius \(r=2^{-34}\). The geometry covers all 32 trial boxes. Ten cached derivative bounds were recomputed at the matching radius, and 22 already-reviewed bounds match the same trial cubes.

![Two panels show 32 input bounds, indexed 01–16 for each of two candidates. Every reported whole-box J lower bound exceeds 0.961. Every cached derivative defect upper bound is below 0.309. Open markers identify 10 fresh cached derivative bounds.](geometry_and_cached_defect.png)

The trial boxes use each case’s fixed exact material values, gain times dose and frozen binary preconditioner. The plotted input index labels records. The points do not track a loading history.

## Whole-box geometry

Each case used an admissible fixed-reference quadratic finite-element mesh. The reviewer independently reconstructed all 650 exact nodal coordinate intervals, checked 144 bulk charts and 24 traces, and replayed the Bernstein domain bounds. Across the 32 cases, the reported lower bounds on the deformation Jacobian \(J\) are all at least **0.961291**. The certified lower bound on strong monotonicity is at least **0.300147** in the frozen coordinate scaling. Bounds on the bulk and trace series-domain ratios are **0.038709** and **0.124962** respectively. These are whole domain checks for the stated finite-element model. They place no bound on continuum error.

## Cached derivative bounds

With cached quadrature Hessian \(H_Q\) and the recorded preconditioner \(R\), the stated quantity is a whole-cube upper bound

\[
\bar Z_Q\ge\sup_{u\text{ in cube}}\|I-RH_Q(u)\|_\infty .
\]

| Cached derivative route | Cases | Matching evidence |
|:--|--:|:--|
| Existing exact-radius bound | 22 | Native bound already covers this trial cube. |
| Fresh exact-radius bound | 10 | Seven smaller trial cubes and three enlarged from radius \(2^{-36}\). |

The largest reported \(\bar Z_Q\) is at most **0.308032**, leaving at least **0.691968** below one in this cached derivative bound. The independently reviewed ten-case extension and the original 22-case bounds enter one same-target manifest. Both recentered states (+ small 13 and − small 16) retain their exact shifts \(\sigma\) from the archived observer fields and the stored-coordinate addition gaps \(\xi\). If \(e\) is a verified spatial root error relative to the proof center, the coordinate difference against the originally saved observer field is bounded by

\[
\epsilon_{\text{nodes}}\le 33.5(\sigma+e+\xi)\,\mu\mathrm m .
\]

## Remaining proof obligations

The spatial root test still needs a signed enclosure of the center residual \(RF_I(u_0)\) to obtain its norm bound \(\eta_I\). It also needs a whole-cube norm bound \(b\) on \(R(H_I-H_Q)\). With \(z_I=\bar Z_Q+b\), the sufficient test requires

\[
z_I<1,\qquad \eta_I+z_I r<r .
\]

The current geometry and cached derivative evidence does not supply either spatial defect. Minimum proofs additionally require the applicable positive symmetric-part premise for \(R\). The earlier four exact-spatial certificates remain distinct. Full exact-spatial outcome coverage, source compatibility, loading-path validity and continuum adequacy remain open.

## Evidence

- [32-case geometry replay](../remaining32_geometry_v1/independent_actual_preparation/actual_replay/receipt.json) · [Ten fresh derivative replays](../remaining32_derivative_extension_v1/independent_actual_preparation/actual_replay/receipt.json) · [Joined input authority](../remaining32_bridge_v1/admission/bridge_input_admission.json)
- [Exact target-context review](../remaining32_targets_v1/independent_prelaunch/review.json) · [Earlier four spatial proofs](../four_spatial_checkpoint_v1/SPATIAL.html) · [Earlier four spatial observer bounds](../observer_application_v2/OBSERVERS.html)
