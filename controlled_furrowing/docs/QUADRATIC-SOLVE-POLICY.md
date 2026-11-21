# Numerical continuation from verified states

The frozen full-P2 energy is unchanged. The earlier target-free pilot retained failed small-gain stationarity checks, particularly at high volume penalty. A later final state could have a small force residual while its overall continuation history still failed. This component changes numerical policy rather than constitutive physics, accepts no failed warm state and retains every attempted increment.

`scripts/robust_quadratic.py` subclasses the frozen full-P2 model. It inherits the energy, gradient, pressure, geometry, reference configuration, boundary conditions and observations. Its only optimizer change is a larger budget: up to 1200 L-BFGS iterations followed by up to 25 safeguarded Newton iterations. The normalized force threshold remains exactly \(10^{-7}\). Positive geometry is still required; neither a small mean volume error nor positive J at quadrature alone substitutes for that gate.

## The accepted-state invariant

The principal API is:

```python
from robust_quadratic import Parameters, SolvePolicy, build_model
model = build_model(Parameters(bulk=5000,
                              basal_surface_tension=2),
                    policy=SolvePolicy())
first = model.advance(1.0, output_directory='private-run/target1')
second = model.advance(5.0, previous=first,
                       output_directory='private-run/target5')
```

The legacy fixed-strength methods are inherited for compatibility; `advance` is the verified-state continuation controller. A supplied previous result must be marked converged, have the identical parameter set, and pass a fresh force and geometry calculation. The state is not trusted just because its file contains a `true` flag. In the absence of a prior result the zero-activation reference is solved and verified first.

Let \((g_k,u_k)\) be the last accepted gain/state and \(g_*\ge g_k\) the next candidate. Each trial starts from an independent copy of \(u_k\). A trial is accepted only if its fixed-strength solve and replay pass the unchanged force/geometry gates. The controller then replaces the accepted pair. On failure it leaves both \(g_k\) and \(u_k\) unchanged and replaces the failed candidate by its midpoint and the original candidate:

\[
 (g_*,\ldots)\longmapsto\left(\frac{g_k+g_*}{2},g_*,\ldots\right).
\]

Thus, by induction, every trial is initialized from a verified state. Failed optimizer vectors are never used as warm starts. The retained initial-vector hashes and accepted-result hashes make this invariant replayable without relying only on prose. A failed state is not promoted merely because a later target happens to converge.

The policy stops explicitly if bisection would require a half-step smaller than 0.001 or if 128 attempts have been made for one requested target. Neither bound guarantees reaching the target. On stopping, the result records the last verified gain, the requested gain and `reached_target=false`; overall `converged` is false even though the last accepted state itself remains verified. There is no numerical success claim for an unreached target. Only nondecreasing activation is currently supported.

## Failed numerical trials are not a physical trajectory

Every actual failed result is retained separately, with force, geometry, vector hash, gain and optimizer history. The accepted path receives its own mechanical potential-work ledger:

\[
 W_k=E(u_k,g_{k+1})-E(u_k,g_k),\qquad
 D_k=E(u_k,g_{k+1})-E(u_{k+1},g_{k+1}).
\]

The telescope \(\Delta E=\sum W_k-\sum D_k\) is checked numerically. These quantities refer to the chosen externally controlled mechanical potential. Failed optimization attempts do not represent time evolution, photons spent or ATP consumption. The controller introduces no biological sensing policy, autonomous regulation or kinetic timescale.

## Target-free test and evidence

The protocol fixes dimensionless \(H=1\), physical conversion \(H=33.5\,\mu\mathrm m\), \(a/H=9.35512037/33.5\), illuminated material depths \(8.5/33.5\) to \(28.5/33.5\), outer radius 3, 12×6 meridional cells, shear modulus 1, apical/basal area moduli 0.2/2, basal tension 2 and apical tension 0. It tests volume coefficients 500 and 5000, with requested gains 1, 5 and 10 in sequence. No measured shortening or larger-region outcome is read or fitted.

Before solving, three random-state objective/gradient comparisons at gains 0/1/5 require bit-for-bit equality with the frozen full-P2 model. This verifies that stronger iteration and adaptive continuation do not change the governing energy. The original failed pilot remains untouched under `runs/quadratic_pilot/`.

`evidence/quadratic_robust_protocol.json` fixes the policy and test before the calculations. `quadratic_robust_verification.json` records the source and all parent hashes, exact energy/gradient comparison and six target outcomes. Complete trial fields, per-target attempt ledgers and final accepted states are retained in `runs/quadratic_robust/`.

All six requested targets were reached:

| Volume coefficient | Target gain | Final force infinity norm | Sampled J range | Attempts / rejected |
|---|---:|---:|---|---:|
| 500 | 1 | 2.73×10⁻¹³ | 0.997548–1.003091 | 1 / 0 |
| 500 | 5 | 1.80×10⁻¹⁰ | 0.979917–1.017614 | 1 / 0 |
| 500 | 10 | 9.09×10⁻¹¹ | 0.966828–1.030380 | 1 / 0 |
| 5000 | 1 | 3.25×10⁻⁹ | 0.999230–1.000903 | 1 / 0 |
| 5000 | 5 | 2.72×10⁻¹¹ | 0.990723–1.004564 | 7 / 3 |
| 5000 | 10 | 4.40×10⁻¹² | 0.980573–1.008819 | 7 / 3 |

These are actual quadrature J evaluations, not exact extrema or conservative coefficient bounds. The remaining local finite-compressibility and continuum-discretization questions are not erased by improved stationarity.

At coefficient 5000, the first attempt from gain 1 to gain 5 failed with force 0.06879. The trial from 1 to 3 also failed. Steps from 1 to 2 and 2 to 3 verified; a subsequent 3 to 5 trial failed, and 3 to 4 then 4 to 5 verified. This is an actual numerical failure/recovery sequence, not a simulated demonstration. Gain 10 required three further rejected attempts. All six rejected optimizer states remain saved. There were eighteen actual attempts and twelve accepted increments in total.

The original numerical batch used controller SHA256 `2e00d7b2b61ecfa9006716ba9a0245af47f2c85ae7b647caf47b66b13532dba9`; that exact source and receipt are preserved under `runs/quadratic_robust/pre_guard_hardening/`. An input-hardening review then added explicit rejection of nonfinite gains/gradients and guaranteed a terminal attempt-ledger flush. It did not change finite-input energy or optimizer arithmetic. Final controller SHA256 is `be6f17e620ddf04533d456fff9d48cd48517f4a0709d921ab510616926b4176b`.

`evidence/quadratic_robust_guard_review.json` verifies nine guard/controller checks and replays all six final states under that final source. The accepted input-hash chains match exactly and the potential-work telescopes agree within 1.43×10⁻¹⁴. The min-step and attempt-limit tests use explicitly labeled injected failures; their fields are kept separate from the actual scientific runs. `quadratic_robust_attempt_replay.json` additionally replays every actual accepted/rejected field, preserving the distinction between a rejected computation and an accepted physical candidate.

This verifies a finite-dimensional numerical controller. It does not establish a continuum error bound, identify material properties, justify an equilibrium-to-transient observation map, validate the larger-region prediction or prove novel biological physics.

![Actual robust continuation attempts, accepted gains and unchanged force gate](../runs/quadratic_robust/solver_policy.png)

The figure displays actual calculations, with accepted attempts in green and rejected attempts in red. It is not a biological trajectory. Image/input hashes are retained in `evidence/quadratic_robust_visual.json`.
