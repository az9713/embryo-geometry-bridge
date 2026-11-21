# Independent review of the uniform energy-gain cell

**PASS for the declared fixed-material weak-load family, with explicit limits on independence.** Exact rational acceptance arithmetic and all frozen inputs were checked separately, then the complete unchanged Arb bound-production code was reproduced into a new output directory. This is not described as an independent reimplementation of every local majorant.

The [review code](../scripts/review_energy_gain_cell.py), [71-check receipt](../evidence/energy_gain_cell_independent_review.json), original [energy-cell proof](ENERGY-GAIN-CELL.md), and isolated reproduction under `runs/energy_gain_cell_independent_review/` are retained. No equilibrium solve or physical-source edit occurred. The original successful receipt was unchanged.

## The claim actually verified

The dimensionless relative gain increment satisfies \(|\delta|\le10^{-3}\); the base drive is one, so the drive interval is \([0.999,1.001]\). All material coefficients, geometry and binary quadrature caches are fixed to the old weak-load case. Write

\[
 u=u_0+S\delta+Wz,\qquad \|z\|_2\le\rho=10^{-6},
\]

where S is a fixed binary displacement predictor and W is a fixed nonsingular binary metric matrix. Uniform Arb bounds give the safe rounded summaries

\[
 \eta\le6.100638\times10^{-7},\qquad
 Z\le0.155994806,\qquad
 \eta+Z\rho\le7.660586\times10^{-7}<\rho.
\]

For every admitted gain, the map

\[
 T_\delta(z)=z-W^TF(u_0+S\delta+Wz;\delta)
\]

is a contraction and self-map of the same Euclidean z-ball. It has a unique fixed point there. W is nonsingular, so \(W^TF=0\) implies \(F=0\). The exact rational review independently verifies the strict inequalities, rather than deciding from rounded decimals. The exact transformed Hessian is symmetric and stays within norm Z of the identity; it is positive definite. The interior stationary point is therefore a strict local energy minimum for this finite-quadrature potential.

## Why the predictor residual is bounded uniformly

The energy depends linearly on the active coefficient. At fixed center, \(F(u_0;\delta)=F_0+B\delta\). Therefore

\[
 F(u_0+S\delta;\delta)
 =F_0+(H_0S+B)\delta
 +\int_0^1[H(u_0+tS\delta;\delta)-H_0]S\delta\,dt.
\]

The full gain increment remains in the Hessian on every path subinterval. Only the positional excursion shrinks with t; scaling the gain itself by t here would bound the wrong residual.

Eight exact adjacent t-intervals cover \([0,1]\). On each, a supremum Hessian-variation bound is multiplied by the exact interval length. The independent scalar audit checks the complete partition, every stored rank/majorant combination, the integral upper sum and the final residual/contraction inequalities. The coarser unpartitioned sufficient test had failed; the failed screen remains retained rather than being called nonexistence.

## The energy majorant and its assumptions

For the volumetric outer product, preserve the center and change operators before taking absolute values:

\[
 \|A_0^TD+D^TA_0+D^TD\|_2
 \le2\|A_0\|_2\|D\|_2+\|D\|_2^2.
\]

The remaining symmetric local Hessian differences admit row-diagonal Loewner majorants. Exact positive quadrature weights preserve those bounds under assembly, and congruence by W preserves semidefinite order. Spectral norms are bounded by maximum absolute row sums only for exactly symmetric target matrices. Each of these conditions was checked in the source and input data. The review independently verifies positivity of all 9,504 cached weights and nonsingularity of the triangular binary W, rather than inferring either from a successful plot or floating eigenvalue.

The complete reproduction recomputed all 9,504 local contributions and all nine transformed majorants. Every serialized final bound and every plan bound matched the retained original. This is a full unchanged-source reproduction. Its mathematical interpretation additionally rests on the independent derivation above and the previously independently checked constitutive derivatives; it is not an external expert review.

## Guard finding and scientific boundary

The successful producer's status was chosen from the numerical acceptance inequalities before combining the final source/input booleans. All booleans in this actual packet are true and all bytes match, so this does not invalidate the result. A future source-drifted packet must not be accepted by status alone. This review requires every numerical and metadata check to be true and rechecks the original bytes before and after reproduction. The limitation is retained rather than silently changing the frozen successful source.

No claim here covers actual source record M, the new seventeen material/gain fields, variation of all material ratios, unbounded tails, a physical activation history, global energetic minimization, continuum/integration error or biological validation. A local stationary family over a gain interval is not the completed identified set. The original five-gate objective remains incomplete.
