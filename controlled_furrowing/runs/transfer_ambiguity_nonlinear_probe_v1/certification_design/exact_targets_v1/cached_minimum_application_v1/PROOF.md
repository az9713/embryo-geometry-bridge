# Coercivity and the local minima of the four cached energies

This is a draft application of the reviewed symmetric-part proof to the four fixed cached-quadrature roots. The four identities are the two frozen candidates at small record L and at the matched large input. Each certificate uses its own exact material point, rational gain times dose, saved binary cache and frozen binary preconditioner. Each energy is the normalized scalar energy differentiated by that cached operator.

For each case, write \(H_Q(u)=\nabla^2 E_Q(u)\). The strictly positive cached determinant, hoop, surface radius, tangent and Jacobian margins certify a valid analytic neighborhood of the closed \(2^{-34}\) free-coordinate cube. The low-Jacobian penalty remains inactive. There are finitely many fixed weighted sites and the free-coordinate injection is linear, so the scalar energy is \(C^2\) and \(H_Q\) is symmetric there. The matching independently reviewed body and trace geometry proves that the entire cube is admissible. The geometry certificate supplies admissibility and the cached derivative replay supplies the Hessian bound.

Let \(P=\operatorname{sym}R\), and use the independently verified \(P\succ0\), exact \(N_R=\|R\|_\infty>0\), and same-cube bound \(\|I-RH_Q(u)\|_\infty\le\bar Z_Q<1\). The applicant takes the larger of the producer and independent replay bounds. Set \(A=RH_Q\). Its eigenvalues lie in the closed disk centered at 1 with radius \(\bar Z_Q\), so \(A\) is invertible and positive stable. The square factors \(R,H_Q\) are invertible. For symmetric \(S=H_Q^{-1}\), the identity

\[
AS+SA^\mathsf T=2P
\]

and decay of \(e^{-tA}\) give

\[
S=\int_0^\infty e^{-tA}\,2P\,e^{-tA^\mathsf T}\,dt\succ0 .
\]

Thus \(H_Q\succ0\). The Neumann norm bound gives

\[
\|H_Q^{-1}\|_\infty \le \frac{N_R}{1-\bar Z_Q},
\qquad
\lambda_{\min}(H_Q)\ge m_Q=\frac{1-\bar Z_Q}{N_R}>0 .
\]

The spectral radius of a symmetric positive definite inverse is its largest eigenvalue and is bounded by its induced infinity norm. A contraction estimate alone would give no sign. For example, \(R=H_Q=-I\) gives \(RH_Q=I\) and zero derivative defect, yet a stationary state of the corresponding quadratic energy is a strict maximum.

The earlier independent root proof gives a stationary root \(u_*\) at distance at most \(e=\eta/(1-z_{\text{small}})\) from the saved center in a \(2^{-36}\) free cube. The script verifies \(e<2^{-36}<2^{-34}\) and checks that the two boxes implement the same center and linear free map. The exact interior margins are \(2^{-36}-e\) and \(2^{-34}-e\). Every line segment in the closed larger cube stays in that cube. Taylor integration from \(u_*\), using \(\nabla E_Q(u_*)=0\), gives the normalized energy gap

\[
E_Q(u)-E_Q(u_*) \ \ge\ \tfrac12 m_Q\|u-u_*\|_2^2 .
\]

The stationary root is a strict local minimizer and the unique minimum on its specified closed cached-energy cube. The exact-spatial energy, continuous loading path, continuum accuracy and a measured fit require separate evidence. This draft computes no matrices, Hessians, forces, AD or equilibrium solves.

[Reviewed symmetric-part method](../preconditioner_spd_v1/METHOD.md) · [Actual independent symmetric-part receipt](../preconditioner_spd_v1/independent_actual_preparation/actual_replay/receipt.json) · [Actual cached derivative replay](../cached_derivative_extension_v1/independent_actual_preparation/actual_replay_v2/receipt.json)
