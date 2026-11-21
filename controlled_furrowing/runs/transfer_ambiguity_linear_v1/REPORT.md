# Can small-region measurements hide larger-region shape differences?

**A local diagnostic found a direction worth testing, not a biological discovery.**
The independently checked rounded affine model permits a large-region shape-contrast range of approximately **0.06042515 µm** while every one of its 68 small-region predictions stays within **0.001 µm of the baseline prediction**. Two endpoints can therefore differ by as much as 0.002 µm in a small-region channel; the 0.001 µm constraint applies separately relative to the baseline.

Under the different constraint that each affine small-region prediction stays within the baseline's worst source-data residual allowance, the large contrast range is approximately **0.30181702 µm**. That allowance is approximately **5.37983567 µm**. It is a calculated residual of an imperfect conditional model, **not a measured noise level**.

![Two conditional local affine frontiers; neither curve is a nonlinear identified set.](local_affine_frontier.png)

## Precisely what was calculated

One existing large-region field was used, at the same baseline material and gain as the archived 17-small-record sensitivity calculation. Only the activation radius differs. The selected large diameter is 49.8939753 µm and the optical dose is 0.89767374. Larger experimental response tables were not read. All values refer to one static current-disk/spatial-maximum observation convention.

The five perturbation coordinates are log apical area modulus, log basal area modulus, additive apical surface tension, log basal surface tension and log activation gain. The first, second, fourth and fifth coordinates lie in [−0.1,0.1]; additive apical tension lies in [0,0.02]. These are local search bounds, not an exhaustive physical parameter domain.

The signed shape observable uses the current physical radius:

\[
Q=z_a(0)-z_a(a),\qquad a=24.94698765\,\mu\mathrm m.
\]

On the edge containing the shoulder, write its quadratic radius and height as \(R(t),Z(t)\), with \(R(t)=a\). Differentiation gives

\[
\delta Q=\delta z_a(0)-\sum_jN_j(t)\,\delta Z_j
+\frac{Z'(t)}{R'(t)}\sum_jN_j(t)\,\delta R_j.
\]

The last term accounts for the node motion needed to evaluate height at a **fixed current radius**. The shoulder is strictly inside one edge, its radial derivative is positive, and the symmetry axis remains fixed radially. Exact rational inverse brackets enclose this stored-graph derivative. They do not enclose the error in the mechanical equilibrium or its parameter derivative.

For the discrete force equation \(F(u,\theta)=0\), implicit differentiation at a regular exact root would give

\[
H S+B=0,\qquad H=\partial_uF,\quad B=\partial_\theta F,
\qquad c=(\partial_uQ)S.
\]

Here the matrices and solution are floating diagnostics at the saved state. The five entries of \(c\), in the coordinate order above, are approximately

\[
c=(-0.04576027,\;0.39682220,\;2.33809289,\;-0.55166505,\;-1.77749089).
\]

The archived small-region map is \(r+Jd\), where \(r\) has 68 entries. We solve both signs of the linear objective \(c\cdot d\), subject either to \(\lVert r+Jd\rVert_\infty\le\delta\), or to \(\lVert Jd\rVert_\infty\le\varepsilon\), and the same trust box. The latter compares predictions with the baseline; it does not assert that either prediction fits the experiment.

## Exact certificates, limited to the rounded affine optimization

For a maximization with constraints \(Ad\le b\), any nonnegative multiplier \(\lambda\) yields

\[
c\cdot d\le\lambda\cdot b+
\sum_j\max\{(c-A^T\lambda)_j l_j,
                 (c-A^T\lambda)_j u_j\}.
\]

This follows by adding \(\lambda\cdot(b-Ad)\ge0\) and maximizing the remaining linear form over the box. The producer retains an exactly feasible rational step and these exact rational weak-dual bounds. The independent reviewer reconstructed all 20 primal/dual certificates without calling the producer's optimizer routine. Their largest optimization gap is approximately 1.12×10⁻¹⁵ µm. This is numerical accuracy for the declared affine optimization, not physical prediction accuracy.

| Allowed small prediction change from baseline, µm | Affine large contrast range lower, µm | Upper, µm |
|---:|---:|---:|
| 0.0001 | 0.007711620 | 0.007711620 |
| 0.001 | 0.060425149 | 0.060425149 |
| 0.01 | 0.110337086 | 0.110337086 |
| 0.1 | 0.206601410 | 0.206601410 |

## Evidence and limits

- [Frozen source and input protocol](preparation/protocol.json), [synthetic and preparation checks](preparation/preflight.json), and [independent prelaunch review](preparation/prelaunch_independent_review.json).
- [Actual result with all exact LP witnesses and multipliers](actual/results.json), [saved mechanical matrices](actual/tangents.npz), [worker terminal](actual/terminal.json), and [root's observed OS exit](actual/os_exit_observation.json).
- [Independent actual verification](actual_independent/receipt.json): rebuilt all 68 rows, checked all 20 optimization certificates, replayed the saved linear systems and an adjoint calculation, and checked moving-radius derivatives with multiple perturbation sizes.

The small-data matrix has full numerical rank five. A weak direction exists, not an exact nullspace. Several small-observable derivatives were obtained from central differences of maxima; switching maximizing locations and one-sided physical directions must be tested in the nonlinear model. Exact mean derivatives do not cure that separate issue.

The present work does not establish two source-compatible nonlinear equilibria, a physical compatible-set frontier, continuous stable loading, continuum or nonaxisymmetric error, an experimentally detectable gap, or novelty. Published epithelial inverse-model and optogenetic mechanics studies already motivate identifiability and geometry-dependent responses; see the [specific prior-work audit](../../docs/TRANSFER-CLAIM-AUDIT.html).

## The next falsifiable test

Before any new nonlinear outcomes, freeze two parameter choices from the local frontier. Re-equilibrate every small-dose case and the same larger input without condition-specific refitting. Compare all 68 nonlinear small summaries and the large contrast against the affine predictions. Retain failed fits and numerical failures. A surviving, substantial gap would justify certification of those particular branches; a disappearing gap would refute this local extrapolation. Experimental compatibility still requires a justified source observation and uncertainty contract.

The original five-gate study remains incomplete. The B and L exact-integral root proofs concern separate existing small-case parameter tubes; they do not certify the new large-field tangent or these proposed witnesses.
