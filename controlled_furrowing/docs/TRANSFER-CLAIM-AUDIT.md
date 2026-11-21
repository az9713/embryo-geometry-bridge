# Audit of a certified ambiguity in small-to-large mechanical transfer

**The general idea is established: different mechanical parameters can fit the same data, and additional experiments can distinguish them.** The sharper candidate is to demonstrate this under the exact finite-height model, complete declared small-region observation contract, and rigorous numerical error bounds. This audit has not constructed that witness pair, established a separating measurement, or confirmed novelty. It identifies the proof obligations and closest overlaps so subsequent computation can test a precise claim.

The audit is bounded: primary local papers and source data were inspected, targeted web searches were executed, and access limits were retained in the [receipt](../evidence/transfer_claim_audit.json). No numerical fitting, larger-outcome-based parameter selection, protocol editing, or external contact was performed.

## The exact candidate claim

Freeze one energy/force law, initial-geometry contract, boundary conditions, dose map, equilibrium/path rule, observation operator, and discrepancy allowances. Let \(\theta\) denote its material and active parameters, and let

\[
\Theta_\Delta=\{\theta\in\Theta_{\rm admissible}:\;
|\mathcal O_{ij}(\theta)-y_{ij}|\le\Delta_j
\quad\text{for every declared small-region record and observable}\}.
\]

The candidate requires two witnesses \(\theta_-,\theta_+\in\Theta_\Delta\), not merely two nearby least-squares minima. They must use the **same** discrepancy vector, allowed geometry/nuisance family, and observation definitions. A tolerance selected afterward to admit two desired witnesses would weaken the claim substantially. A common overall energy-scale change that leaves every deformation unchanged is also insufficient: it is a gauge ambiguity with no transfer consequence.

At minimum the existing 17-record table contains 68 scalar mechanical summaries: maximum shortening, region-averaged shortening, and separate apical and basal maxima. Matching maximum shortening alone is not matching all of them. The three illustrative apical/basal trajectories contain 63 timepoints, but lack a verified crosswalk into all 17 cohort identities. If the claim says **all available small-region dynamics**, those traces and a valid dynamic/time-reference model must also be addressed. A deliberately quasistatic claim should instead name the endpoint summaries it covers and avoid claiming full temporal compatibility. Neither a missing specimen crosswalk nor excluding dynamics from a stated quasistatic claim demonstrates that the data are intrinsically dependent or unusable.

For a prespecified larger activation condition \(d_L\), require validated response intervals \(Q_-(d_L)\) and \(Q_+(d_L)\) with

\[
\operatorname{dist}(Q_-(d_L),Q_+(d_L))=\gamma>0.
\]

The response must be a defined scalar or field functional with units and an observation operator—not an impression from two deformed plots. Examples include a signed profile contrast, a basal/apical displacement partition, or a source-compatible slice observable. A qualitative label whose allowed set is \(\{0,1\}\) for both witnesses cannot supply separation.

This would prove that the stated small-data contract does not determine that larger response within the declared mechanical class. It would not prove that biology is fundamentally unpredictable, that the model is biologically correct, or that every small-data-compatible parameter has been found.

## Closest prior work: what already overlaps?

| Primary work and inspected locator | Established overlap | Difference still requiring proof here |
|---|---|---|
| **Kumar et al. 2024**, [Reverse engineering morphogenesis through Bayesian optimization of physics-based models](https://doi.org/10.1038/s41540-024-00375-z), Results Sec8–9, Fig4, Supplementary Table1 and Fig3 | A finite-height epithelial cross-section model already yields two parameter clusters with similar control shapes: different cell heights and apical/basal contractility. Local Hessian analysis exposes weak parameter directions; extra local cell-length measurements are proposed. The supplement lists actual Control-1 and Control-2 parameter sets, not just a generic warning. | Collagenase-perturbed contours are optimization inputs with newly inferred parameters. The inspected comparison does not establish a small-patch-only, all-summary, no-refit larger-patch separation with rigorous numerical intervals. Do not claim the first equifinal finite-height epithelial model. |
| **Countryman et al. 2025**, [Endogenous OptoRhoGEFs reveal biophysical principles of epithelial tissue furrowing](https://doi.org/10.1038/s41467-025-62483-6), Results Sec7–8, Fig4–5, Methods Sec14–17 | Activation geometry changes furrowing/bending. Basal stiffness and external compression are competing mechanisms; microperforation and oriented ellipses test them experimentally. Similar shortening does not imply similar geometry. | The proposed contribution must be a quantified calibration-to-transfer ambiguity under a complete common contract, not a rediscovery of geometry dependence, basal mechanics, or the value of perturbations. |
| **Boocock et al. 2021**, [Theory of mechanochemical patterning and optimal migration in cell monolayers](https://doi.org/10.1038/s41567-020-01037-7), Methods and Extended Data Fig5–6 | Mechanical and optogenetic perturbations calibrate a mechanochemical monolayer model; explicit identifiability plots compare relaxation parameters using boundary and bulk displacement data. | This is an in-plane mechanochemical migration setting. Its presence already rules out claiming that optogenetic mechanics plus identifiability analysis is new. Final publisher access was incomplete here; indexed primary captions and the author repository were checked. |
| **Ogita et al. 2022**, [Image-based parameter inference for epithelial mechanics](https://doi.org/10.1371/journal.pcbi.1010209), Discussion; Materials and Methods | Shape-based force balance determines only relative parameters without force calibration, cannot by itself decompose active/passive or intrinsic/extrinsic processes, and can suffer correlated parameter terms. Additional force measurements and perturbations are explicitly required. | Their planar junction-level inference is not the present finite-height activation-radius transfer test. A generic force-scale degeneracy or suggestion to measure another variable is already covered. |
| **Oddes and Solav 2023**, [Identifiability of soft tissue constitutive parameters from in-vivo macro-indentation](https://doi.org/10.1016/j.jmbbm.2023.105708), primary abstract | Axisymmetric hyperelastic inverse FE examples show that force-depth data can leave parameters poorly identified, while full-field surface displacement improves identification. | This strongly overlaps the proposed move from scalar endpoints to shape fields. It does not, on the inspected abstract, establish this active-epithelial, all-small-observable transfer certificate. Full methods were not inspected in this audit. |
| **Stegmaier, Skanda and Lebiedz 2013**, [Robust Optimal Design of Experiments for Model Discrimination Using an Interactive Software Tool](https://doi.org/10.1371/journal.pone.0055723), Methods: Model discrimination and Robust optimal designs, Eq4–5 | Models fitting an initial dataset can be separated by designed perturbations. Worst-case parameter variation is incorporated through max–min discrimination and iterative parameter/design updates. | Pair separation and robust optimal measurement selection are established methodology. A nonlinear finite-height mechanics certificate with explicit forward-error control is a specific application requiring demonstration, not a new general design principle. |
| **Karimshoushtari and Novara 2020**, [Design of experiments for nonlinear system identification: A set membership approach](https://doi.org/10.1016/j.automatica.2020.109036), primary abstract | Bounded-error identification and input-constrained experiment design already optimize worst-case model uncertainty using a radius-of-information criterion. | A deterministic compatibility set is not a new replacement for probabilistic inference. Only the source-specific mechanical result and its certified limits could be additional. Full methods were inaccessible here. |
| **Grabowski, Nałęcz-Jawecki and Lipniacki 2023**, [Predictive power of non-identifiable models](https://doi.org/10.1038/s41598-023-37939-8), primary abstract and indexed results | A signaling cascade and mechanical analogue illustrate that useful predictions can survive parameter nonidentifiability; new variables or stimulation protocols change predictive capability. | Parameter ambiguity alone does not prove transfer ambiguity. The larger-response gap must actually be demonstrated; unidentifiable parameters may be irrelevant to that response. Full publisher text was access-limited. |
| **Liu, Maini and Baker 2026**, [Optimal experiment design for practical parameter identifiability and model discrimination](https://doi.org/10.1016/j.mbs.2026.109710), introductory model-discrimination discussion | The primary discussion explicitly relates several data-fitting parameter sets to competing mechanistic hypotheses, followed by design of an input maximizing their output differences. | This closely overlaps the generic proposed workflow. Any novelty claim must concern the specific certified mechanical statement, not that recipe. Access details and inspected coverage are in the receipt. |

These sources do not collectively prove the exact proposed theorem. Equally, failure to find an identical theorem in this bounded search does not establish novelty. The closest direct challenge is Kumar's actual pair of equifinal epithelial parameter clusters, combined with established robust discrimination methodology. Any manuscript should compare against those concrete results, not just cite broad inverse-problem terminology.

## What would count as a rigorous witness?

**1. Nonempty calibration compatibility.** A common monotone peak contract already needs at least 3.131943715 µm uniform maximum-shortening residual for the canonical CSV; stricter concavity needs approximately 3.648762214 µm. These [compatibility floors](CALIBRATION-UNCERTAINTY.md) are lower bounds, not justified measurement-error allowances or evidence that any actual mechanical parameter fits all 68 endpoints. An empty \(\Theta_\Delta\) proves model/contract incompatibility, not ambiguity between two members.

If numerical evaluations satisfy \(|\widehat{\mathcal O}_{ij}-\mathcal O_{ij}|\le\varepsilon_{ij}\), a sufficient membership certificate is

\[
|\widehat{\mathcal O}_{ij}(\theta_\pm)-y_{ij}|+\varepsilon_{ij}(\theta_\pm)\le\Delta_j
\quad\text{for every included }i,j.
\]

An interval straddling the tolerance boundary does not certify membership. Dose uncertainty, specimen-specific growth references, and allowed geometry variation must be given explicit quantifiers; they cannot be silently changed between witnesses.

**2. Physical and mechanical admissibility.** Require positive passive moduli, consistent prestress and initial force balance, admissible activation/control parameters, orientation preservation and noninterpenetration, the declared volume constraint, and a specified stable equilibrium/path. A small Newton residual or positive sampled Jacobians alone do not prove all of these. Multiple stable branches must either be carried or handled by a declared branch-selection rule. If the two witnesses use different initial shapes, say that the ambiguity includes geometry; do not attribute it solely to basal stiffness.

**3. A validated forward-error gap.** If scalar response evaluations have certified error radii \(\varepsilon_-\), \(\varepsilon_+\), then

\[
|\widehat Q_+-\widehat Q_-|>\varepsilon_++\varepsilon_-
\]

is sufficient for a nonzero exact-response gap. Mesh refinements, optimizer convergence and exact geometry certificates address different questions; they do not automatically provide a continuum response-error radius. A result restricted to the finite-dimensional discrete model must say so. Observation-plane, initial-geometry or other allowed nuisance uncertainty must also be propagated into the response intervals before claiming robust separation.

**4. A genuine transfer consequence.** For any single prediction \(p\) based only on the shared calibration data, the triangle inequality gives

\[
\max\{|p-Q_-|,|p-Q_+|\}\ge\tfrac12|Q_+-Q_-|\ge\gamma/2.
\]

This standard two-point optimal-recovery argument would establish an unavoidable worst-case prediction ambiguity **within the assumed compatible class**. It does not assign probabilities or make either witness the true embryo. Two existence witnesses suffice for this narrow counterexample, but do not replace a full identified-set computation required for stronger complete-set predictions.

## Which extra observation is realistically worth testing?

The strongest first candidate is an **uncropped, registered small-region membrane stack that measures the radial apical and basal displacement profiles or basal area expansion**, alongside the preactivation AP/DV surface geometry. Countryman's membrane imaging and segmentation demonstrate an applicable modality. The public composite movie does not establish the necessary annotation accuracy or per-specimen profile itself. The [movie audit](MOVIE-GEOMETRY-FEASIBILITY.md) documents the clipped apical edge and unresolved plane/calibration convention; those display limits must not be converted into new data.

The full shape field can constrain mechanisms hidden by four extrema; measuring initial geometry can independently remove a nuisance explanation. Neither statement guarantees separation for a particular pair. Before recommending an acquisition, evaluate the candidate extra observable \(m_d\) on both retained mechanical witnesses, choose its sampling location/time by their forward predictions alone, and propagate all errors. For predicted centers \(\widehat m_\pm\), numerical radii \(e_\pm\), a measurement error bound \(\eta\), and additional declared model discrepancy radii \(\rho_\pm\), a sufficient discrimination criterion is

\[
|\widehat m_+-\widehat m_-|
>e_++e_-+2\eta+\rho_++\rho_-.
\]

Then the two expanded observation intervals are disjoint. A new measurement could support at most one of them, or reject both. Its required precision is a design output; nominal 2 µm z-spacing is not a validated \(\eta\). No numeric attainable separation has yet been established by this audit.

If the remaining ambiguity is a force-scale gauge, more geometry cannot resolve it; independently calibrated force information is needed, as the epithelial inverse-mechanics literature already explains. Laser recoil without a drag/viscosity calibration is not automatically an absolute stiffness measurement. If the goal is robust discrimination between whole mechanism families \(H_-\) and \(H_+\), optimize the minimum separation over their allowed parameters. Optimizing a gap between two convenient points establishes only pair discrimination. Minimizing over two identical unrestricted copies of \(\Theta_\Delta\) would give zero by choosing the same point twice.

## Contribution and stopping rules

| Statement | Current assessment |
|---|---|
| Nonunique mechanical fits; extra perturbations can help | Established concept; not new here |
| Two source-specific mechanically admissible witnesses fitting every declared small observable | Not constructed or certified by this audit |
| Larger-condition prediction intervals rigorously disjoint | Not established by this audit |
| One practical extra measurement provably separates those witnesses | Candidate modality identified; quantitative separation and precision not established |
| Exact prior publication of the complete finite-height/all-observable/error-certified result | No exact match established in this bounded comparison; originality remains unresolved |
| Full original Gate3/5 requirement satisfied by this literature audit | No; actual compatibility, forward certificates and contribution assessment remain necessary |

Proceed computationally only after the common contract is explicit. If all candidate parameters fail it, report incompatibility. If compatible witnesses yield overlapping larger responses, report that parameter uncertainty has not destroyed this prediction. If a gap survives all declared errors but no feasible extra observation separates it, report unresolved transfer ambiguity. A claim of novelty becomes assessable only after a concrete result exists and the comparison can address its exact hypotheses and conclusions; adding certification alone should not be advertised as a new biological principle.
