# Independent review of the actual material/gain replay

**Verdict: PASS for the completed seventeen-record nonlinear candidate replay.** The new-material zero state and all seventeen requested input states pass fresh Q8/Q12 force checks, exact stored-P2 geometry checks, nodal/boundary reconstruction and continuation provenance. The primary conditional discrepancy strictly decreases for the saved fields, from approximately 5.379835675 to 5.186387251 micrometres.

This is a directly evaluated improvement of one predeclared **joint material/gain candidate** at one finite-element scenario. It is not a nonlinear optimum, a verified exact-root/continuum statement, an improvement under every observation convention or a completed physical parameter inference.

The independent program and exact endpoint/residual intervals are retained in `evidence/material_cohort_direct_independent_review.json`. All **427 numerical/provenance checks pass** before final packaging checks. The reviewer launched no equilibrium calculation, restarted no job and edited no physical source.

## One frozen proposal and all seventeen records

The sole job completed all seventeen records under frozen driver `448675cc6f59973c07a372f119aea5dd9be0f915398f57a3f9bd7d5cfe91074e` and protocol `1ed458c8f0f7b9c9bacfdf0c04b907d474d06e5607ced5956ecc2d2f1adac953`. The terminal results hash is `638a95bc97055e33a10e4f6c34f4750629b23f53688d19b3c8f9e7ca1ecce2b7`.

The proposed ratios and common gain, displayed approximately, are:

| Quantity | Previous | New |
|---|---:|---:|
| Apical area coefficient | 0.2 | 0.2210341836 |
| Basal area coefficient | 2 | 1.8096748361 |
| Apical baseline tension | 0 | 0 |
| Basal baseline tension | 2 | 2.2103418362 |
| Common activation gain | 8.0629720200 | 7.8126545043 |

The finite logarithmic-volume coefficient remains 5000, shear coefficient remains the force unit, and geometry/discretization remain unchanged. The five-coordinate step is transformed using the declared stored binary64 exponential/multiplication rule; the review independently reconstructs that rule. The displayed decimals are not measured material constants or exact transcendental values.

All seventeen unique source-column records are retained, including I and O. Every requested strength is the exact retained gain-times-dose product; actual solver inputs and their rational binary64 rounding differences match the protocol. No per-record gain or post-result material refit is introduced. The other three convention-specific affine proposals are not separately solved.

## Mechanical and geometric verification

The driver first solves and certifies a **new-material zero** field. Its Q8/Q12 force residuals independently replay at about 1.00e-12 and 5.23e-13. It is not an old-material field relabeled as a valid new-material initializer.

Across all seventeen nonzero targets:

| Independent check | Outcome |
|---|---|
| Fresh Q8 normalized force | Maximum 1.162e-10, below 1e-7 |
| Fresh Q12/surface-16 force on the same field | Maximum 1.161e-10, below 1e-7 |
| Exact determinant/global-embedding certificates | All targets and zero pass |
| Reference/current node and connectivity reconstruction | Bit-for-bit match |
| Stored P2 boundary controls versus actual nodes | Bit-for-bit match |
| Fixed coefficients, axis radius and clamped outer rim | Preserved |
| Exact boundary radial monotonicity | Pass |
| Observer reporting intervals | Resolved on all states |

Exact rational certificate recomputation covers 2592 triangle instances across the eighteen fields. One outward common enclosure for all nonzero targets is

\[
 0.95829100124<J<1.02340430069.
\]

The exact fractions, field hashes and canonical hashes of the independently recomputed geometry certificates are in the receipt. These certificates concern the stored discrete maps. They do not enclose a continuum solution or establish that each floating small-residual field contains an exact root.

The Q12 check is a force evaluation with higher numerical integration at the same coordinates, without reoptimization. Agreement supports this integration check; it is not a rigorous continuum or uniform quadrature-error bound. No energetic stability or global/unique equilibrium proof is supplied.

## Verified continuation and retained attempts

The targets are reached in ascending input order from the new-material zero. All seventeen requested increments are accepted and there are zero rejected attempts. For every retained attempt, the review independently checks parameters, input strength, initial verified-vector hash, final-vector hash, fresh force and local orientation screen. The next initializer is exactly the preceding verified vector, and each final target matches its chain endpoint.

Every accepted target additionally has the full Q12 and exact geometry certification. Original-material fields appear only as the baseline comparison, never as verified new-material warm starts.

The old baseline fields were selected from a different, denser input-continuation history. No uniqueness or branch-identity theorem connects those old and new paths. The verified statement is that the declared solver/continuation procedure produced these new stationary-tolerance fields and their observed loss; it is not a proof of a unique local material branch.

## All four operators on the same new fields

The independent audit reobserves each new field and each corresponding old baseline field under all four declared conventions. Every comparison uses the same seventeen source vectors, common fixed current-radius ROI, reference heights and canonical source-column pairing. The old comparison for alternate conventions is at the **same original primary common gain**, not their separately selected old gains.

The following are approximate micrometre losses; the receipt retains exact rational observer/residual intervals:

| Conditional observation convention | Old same baseline | New primary proposal | Old minus new |
|---|---:|---:|---:|
| Disk / spatial-max components | 5.379836 | 5.186387 | +0.193448 |
| Disk / regional-mean components | 6.467716 | 6.528489 | -0.060772 |
| Line / spatial-max components | 6.183270 | 5.978243 | +0.205027 |
| Line / regional-mean components | 6.183270 | 6.064266 | +0.119004 |

Thus one convention worsens. The primary proposal was selected before these nonlinear results; the audit does not switch conventions afterward to manufacture an improvement.

For an observed value y and predicted interval [l,u], the independent residual bounds are

\[
 \inf_{x\in[l,u]}|x-y|=\max(l-y,y-u,0),\qquad
 \sup_{x\in[l,u]}|x-y|=\max(|l-y|,|u-y|).
\]

Taking maxima over all 68 channels gives the reported discrepancy enclosure. The primary decrease is rigorously enclosed for these stored fields by the outward interval

\[
 \Delta_{\rm old}-\Delta_{\rm new}
 \in[0.193448423,\,0.193448424]\;\mu\mathrm m.
\]

The new worst channel is the apical component of source-column M; the regional-mean height channel of source-column L is close behind. The residual remains above five micrometres, and no experimental discrepancy/noise budget has been established here.

## What the affine prediction learned

The primary affine prediction was approximately 5.183862986 micrometres. The actual new discrepancy is higher by about 0.002524265 micrometres. A direct nonlinear replay was necessary to establish that this particular proposed improvement survives the finite step. The difference at this one step is an empirical surrogate check, not a bound on all future nonlinear remainders.

Gain changes together with the material coefficients, and the old gain was a verified surrogate-selected candidate rather than a certified nonlinear gain optimum. Consequently this result does not isolate the benefit of material-ratio changes or prove improvement beyond reoptimizing gain alone. It demonstrates a lower-loss joint candidate under the primary conditional operator.

## Prelaunch documentary provenance, kept distinct from numerical correctness

The launch job captured review SHA `3fdef9e44a23a572b12e6fc11bec47f5f753eaecb5aefa0f1277609eed90879c`. The final documentary review SHA is `cdcb914af204518ab0cc82f1167d25629699b641058ab7bbde2fabb4baddf5b4`.

The reviewer independently checked that both recovered launch files are byte-identical and match the captured hash. The later file adds only `prelaunch_revision_archive_hashes`; every original key/value, PASS10 check and driver/protocol/candidate/affine binding is unchanged. All appended archive hashes match their retained files. The launch and final documentary hashes are deliberately distinguished, with the recovery addendum retained. No source, numerical policy, parameter or launch was changed by that documentary update.

This resolves the review-file chronology; it is not used as a substitute for the independent force/geometry/observation replay above.

## Limits of the completed component

The component now establishes an actual discrete primary loss decrease, rather than only an affine proposal. It leaves unresolved exact-root enclosures, spatial convergence, physically constrained material/geometry uncertainty, source temporal-peak semantics and larger-region prediction. Reobserving four static conventions is not four independent data sets or an experimental validation.

The earlier continuum/mesh adequacy deficit is not cleared by force convergence at this one grid. No full compatible set, global nonlinear optimum, biological mechanism or confirmed novelty is established. The scientific advance within this workflow is a fully checked local candidate update and a clear account of where its interpretation stops.
