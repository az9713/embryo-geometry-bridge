# Uniform stationary roots on one shared cohort cell

**Producer root-stage status: 17 of 17 records now have successful uniform-root receipts on the same shared cell.** Fifteen use the version-two packets; F and K require the separately declared single proof-center refinement described below. Every failed original test remains retained. Source M reuses its immutable independently reviewed certificate; no physical parameter or shared cell was changed.

These are producer-stage statuses. The original cohort replay and the separate recentered-root audit are distinct verification gates. The old same-center observer cannot consume the new F/K packets without an explicitly reviewed center-offset adapter.

The [driver](../scripts/cohort_uniform_roots.py), [protocol](../runs/cohort_uniform_roots_v1/protocol.json), [preflight](../runs/cohort_uniform_roots_v1/preflight.json) and [guard tests](../runs/cohort_uniform_roots_v1/guard_tests.json) precede the expensive launch. A separate independent launch review must bind the exact source, protocol and preflight bytes. The initial prelaunch revision and its missing guard findings are preserved under `prelaunch_revision_0/`.

The first execution stopped after eleven records because its M-reuse guard compared two equal rational doses as literal strings: `0.78208913` and `78208913/100000000`. This implementation error did not change a physical target or invalidate an earlier numerical payload. The original source, failed job, eleven results and caches remain unchanged.

The separately reviewed [completion driver](../scripts/cohort_uniform_roots_completion.py) and [version-two protocol](../runs/cohort_uniform_roots_v2/protocol.json) reuse those eleven payloads exactly, compare M's dose as a Fraction, and calculate only the five remaining records. Its actual process exited successfully, with [terminal results](../runs/cohort_uniform_roots_v2/results.json) reporting 15 certified and two unresolved. The retained F/K tests failed a sufficient self-map inequality; they do not prove root nonexistence.

The independently reviewed [radius-only supplement](../runs/cohort_radius_supplement_v1/protocol.json) also terminated successfully, but all six additional sufficient tests remained unresolved. Its [results](../runs/cohort_radius_supplement_v1/results.json) keep the same parameter cell, predictors and frozen caches. No optimization, cell shrink or new preconditioner preparation was used; Arb point/Hessian assembly was replayed for validation.

| Record | Initial tested radius | Initial \(Z\) | \(Z\) at 2× radius | At 4× | At 8× |
|---|---:|---:|---:|---:|---:|
| F | \(1.73080\times10^{-10}\) | 0.69039 | 1.34310 | 2.64868 | 5.25982 |
| K | \(1.04523\times10^{-10}\) | 0.56532 | 1.04794 | 2.01317 | 3.94384 |

The initial predictor residual upper bounds are approximately \(8.65398\times10^{-11}\) and \(5.22615\times10^{-11}\), respectively. Initially \(Z<1\), but \(\eta+Zr<r\) fails. Increasing these radii makes the retained componentwise derivative bound exceed one. These failures diagnose this proof enclosure; they do not establish nonexistence, instability or biological incompatibility.

The [fresh bound decomposition](COHORT-BOUND-DIAGNOSTICS.md) reproduces both residual bounds exactly. The center-residual term dominates, while the predictor-only derivative bounds are below 0.083. The error-cube enclosure, rather than a demonstrated lack of roots, is the immediate proof obstacle.

## One cell and seventeen exact input laws

Every record uses the **same** five-dimensional actual-M `scale_3` cell: relative half-width \(10^{-8}\) for apical/basal area moduli, basal surface tension and common gain, with apical tension in \([0,2\times10^{-10}]\). No per-record cell shrink or response-dependent selection is permitted. The fixed model retains \(\mu=1\), finite volumetric coefficient \(\kappa=5000\), the 12×6 reference geometry and the actual Duffy-8/surface-12 integration policy.

For source record \(i\) with canonical rational dose \(I_i\), the active coefficient is exactly

\[
 f_i=gI_i.
\]

The earlier floating solver's rounded binary strength is retained separately. It is not substituted for this parameter law. The saved vector supplies a proposed center and numerical preconditioner; the interval proof must validate the exact law, including any resulting residual.

The [common model/observation contract](../runs/cohort_static_observers_v1/model_observation_contract.json) excludes record identity, dose and displacement center. Those belong to per-record metadata. The root stage proves no observation agreement; a separate stage will compose all four conditional static outputs.

## Uniform root statement and the tuple implication

For a record \(i\), with center \(u_{i,0}\) and stored binary predictor \(S_i\), define

\[
 \widehat u_i(\theta)=u_{i,0}+S_i(\theta-\theta_0).
\]

The reviewed kernel encloses the predictor residual and the preconditioned derivative on \(\widehat u_i(\theta)+[-r_i,r_i]^{611}\), uniformly for every \(\theta\) in the shared cell. It includes the explicit numerical tangent defect. If

\[
 Z_i<1,\qquad \eta_i+Z_i r_i<r_i,
\]

Banach's theorem gives a unique stationary root in that record's declared nearby domain for every parameter in the cell. Each field must account for all 611 free coordinates and all 9,504 local quadrature terms.

If all seventeen receipts pass, the finite conjunction gives

\[
 \forall\theta\in C\quad
 \exists!\,(u_1,\ldots,u_{17})\in
 \prod_{i=1}^{17}D_i(\theta):
 \quad F_i(u_i,\theta)=0\quad(i=1,\ldots,17).
\]

Uniqueness is restricted to the product of the selected nearby domains. This is an independent stationary-root tuple, not a common stable or time-dependent path through the experiments. Other stationary branches are not excluded.

## Attempt and retention policy

There is one affine cell attempt per newly prepared record. Its one tested error radius is the larger of \(10^{-11}\) and twice the validated residual upper bound, subject to the declared maximum radius. A failed sufficient test remains **UNRESOLVED**, never OUT. The parameter cell stays unchanged, and the driver continues the other records. A changed source, data file or numerical environment is an integrity failure that stops the job.

Preparation checks require nonempty all-boolean passing tangent verification, exact agreement with the saved free vector/reference/mesh, the frozen integration rule, and the expected dimensions. Thirteen synthetic guard checks reject empty/failed/nonboolean verification, mismatched saved vectors, geometry, integration and term counts, and paths outside the owned workspace/output tree. These are guard tests, not a replacement for the actual root proof.

Each completed cache has a final hash-bound preparation manifest. Interrupted preparation uses a fresh attempt folder; completed per-record results are reused only after their context and artifact hashes match. The Windows job records its real PID and process creation time; the controller checks those read-only and never sends a process signal. A live matching lease blocks duplicate launches. A stale lease requires explicit resume and is archived.

## Limits

This calculation uses only the declared small-region source data and already saved training fields. There is no new equilibrium optimization, held-out response access, gain/material fitting or cell shrink. Root existence alone does not establish an all-68 fit, a biological mechanism, an experimental noise model, energetic stability, a common trajectory, full-body geometry, integration/continuum accuracy, full material-domain coverage or novelty.

## Single proof-center refinement: the actual successful step

The diagnostic result justified a new, frozen [protocol](../runs/cohort_recentered_roots_v1/protocol.json), [prelaunch review](../runs/cohort_recentered_roots_v1/prelaunch_independent_review.json) and [driver](../scripts/recenter_uniform_roots.py) for F/K only. For each record, exactly one numerical proof center was chosen:

\[
 u_1=\operatorname{fl}_{64}\bigl(u_0-\operatorname{mid}(RF(u_0,\theta_0))\bigr).
\]

The exact dyadic Arb midpoint is subtracted from the exact binary old coordinate before a single nearest-even binary64 rounding. A synthetic double-rounding counterexample verifies that the implementation does not first round the correction separately. The midpoint choice is not itself a root certificate.

The original preconditioner \(R\), affine predictor \(S\), parameters, gain-times-dose law, shared cell, reference geometry and quadrature remain unchanged. Every cached array except `u` has an identical fingerprint. A separate physical-potential fingerprint also binds model parameters, base coefficients, dose and cell. Old and new vectors, exact translations, midpoint enclosures and hashes are retained. The auxiliary state is explicitly a numerical proof center with `converged=False`, not a new observation or a claimed nonlinear-solver result.

Fresh Arb assembly and the full affine-cell theorem then gave:

| Record | Center translation \(\|u_1-u_0\|_\infty\), approximately | Fresh \(\eta\) upper | Fresh \(Z\) upper | Normal corrector radius |
|---|---:|---:|---:|---:|
| F | \(8.58448\times10^{-11}\) | \(7.42374\times10^{-13}\) | 0.0754082 | \(10^{-11}\) |
| K | \(5.00616\times10^{-11}\) | \(2.49871\times10^{-12}\) | 0.128929 | \(10^{-11}\) |

Both satisfy the strict self-map and contraction inequalities uniformly over the **full unchanged cell**. The actual process exited successfully; [terminal results](../runs/cohort_recentered_roots_v1/results.json) retain both receipts. This repairs a proof-enclosure limitation rather than introducing a different material fit.

The normal Banach radius is centered at \(u_1+S(\theta-\theta_0)\). The center translation is kept separate; it was not inserted into that radius or its \(\eta,Z\) calculation. Any observer based on the original saved curves must charge the translation and serialization gap explicitly. The previous same-center observer remains an invalid consumer for these packets until a versioned adapter is reviewed. No common physical trajectory, stability theorem, continuum result or full material-domain coverage follows.

## Energy-norm route, not executed for these F/K packets

The [energy-norm study](PARAMETRIC-NORM-FEASIBILITY.md) and [separate weak-field gain certificate](ENERGY-GAIN-CELL.md) demonstrate an alternative to collapsing every local Hessian change into an absolute inverse/row bound. They do not yet certify F or K.

For either unresolved field, freeze a binary metric factor \(W\) and validate \(\|W^TH_0W-I\|_2\le\epsilon_0<1\) against the **actual exact-product target**. Enclose a parameter-dependent energy tube

\[
 u=\widehat u(\theta)+Wz,\qquad\|z\|_2\le\rho.
\]

With \(d=\theta-\theta_0\), a concrete residual bound is

\[
 \eta\le
 \sup_{d\in C-\theta_0}
 \|W^TF_0+W^T(H_0S+B_0)d\|_2+
 K_{\rm pred}\sup_{d\in C-\theta_0}\|W^{-1}Sd\|_2.
\]

Here \(K_{\rm pred}\) bounds the transformed Hessian change along every predictor segment with the full parameter held fixed in the Taylor integral. The tangent metric norms can be bounded using \(S^TH_0S/(1-\epsilon_0)\). Whole-tube local widths must use \(\rho\|(BW)_{j,:}\|_2\), and PSD local majorants must be assembled **before** metric contraction. All widths and global norms require outward validation. The proof still must satisfy \(Z<1\) and \(\eta+Z\rho<\rho\); merely naming an energy norm supplies neither gate.

If that succeeds, its declared uniqueness domain is an energy ball. A coordinate box obtained from row norms of \(W\) can enclose the root for body/observer transfer, but uniqueness in that larger box would require a separate argument. No implementation or success of this F/K alternative is claimed; the actual successful computation used the single proof-center refinement above.
