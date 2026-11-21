# Four strict minima of the cached energy

The two frozen material candidates each have a verified strict local minimum at small record L and at the matched large input. These four minima belong to the normalized binary \(Q8/S12\) quadrature energy. Exact-spatial energy minima remain a separate proof obligation.

## Four positive curvature bounds

The earlier root certificate locates the stationary point \(u_*\) inside a closed \(2^{-36}\) free-coordinate cube. The new minimum proof gives uniform positive curvature on the concentric \(2^{-34}\) cube with the same center, linear free coordinates, exact gain times dose, and frozen binary preconditioner \(R\).

![Four positive lower bounds on the smallest normalized cached Hessian eigenvalue, ranging from 4.0988 to 6.5349 times 10 to the minus four.](curvature_bounds.png)

| Frozen target | Curvature \(m_Q\) (×10⁻⁴) | Root interior margin (×10⁻¹¹) |
|:--|--:|--:|
| Candidate + · small L | 4.0988 | 5.8084 |
| Candidate + · large | 5.9296 | 5.6037 |
| Candidate − · small L | 4.3560 | 5.2856 |
| Candidate − · large | 6.5349 | 5.8158 |

Both columns are lower bounds rounded downward. The interior margins refer to the larger cube in dimensionless free coordinates. The curvature is for the normalized cached energy in those coordinates.

## Why the sign follows

On each admissible cube, positive cached-site analytic margins make the scalar energy \(E_Q\) twice differentiable with symmetric Hessian \(H_Q\). The fresh symmetric-part proof gives \(P=\operatorname{sym}R\succ0\). The larger-cube cached derivative bound gives \(\|I-RH_Q\|_\infty \le \bar Z_Q<1\). With \(A=RH_Q\) and \(S=H_Q^{-1}\), \(A\) is invertible and has eigenvalues with positive real part. The Lyapunov identity yields

\[
\begin{aligned}
AS+SA^\mathsf{T}&=2P,\\
S&=\int_0^\infty e^{-tA}\,2P\,e^{-tA^\mathsf{T}}\,dt\succ0.
\end{aligned}
\]

The Neumann bound and the exact binary norm \(N_R=\|R\|_\infty\) give the uniform quantitative result

\[
\lambda_{\min}(H_Q)\ge m_Q=\frac{1-\bar Z_Q}{N_R}>0 .
\]

Invertibility and contraction alone give no sign. For example, \(R=H_Q=-I\) gives \(RH_Q=I\), yet the quadratic stationary point is a strict maximum.

The root lies strictly inside the larger cube. Every segment from \(u_*\) to \(u\) in that cube stays admissible. Taylor integration with vanishing gradient at \(u_*\) gives

\[
E_Q(u)-E_Q(u_*)\ge \tfrac12 m_Q\|u-u_*\|_2^2 .
\]

This proves one strict minimum per specified cached-energy cube. The four certificates concern isolated fixed targets. Continuous loading, exact spatial integration, continuum accuracy and a measured biological fit require further evidence.

## Proof records

The independently verified receipt resolves the review status of the archived application draft. The exact rational bounds and the same-target proof links are retained for replay. [Independent minimum receipt](independent_readback.json) · [Detailed derivation](PROOF.md) · [Original scalar draft](results.json)

[Cached roots](../../cached_pilots_v1/independent_actual_preparation/actual_replay/receipt.json) · [Whole-cube geometry](../independent_actual/replay/receipt.json) · [Larger cached derivative bounds](../cached_derivative_extension_v1/independent_actual_preparation/actual_replay_v2/receipt.json) · [Actual symmetric-part proofs](../preconditioner_spd_v1/independent_actual_preparation/actual_replay/receipt.json) · [Four-pilot checkpoint](../certification_checkpoint_v1/CHECKPOINT.html)
