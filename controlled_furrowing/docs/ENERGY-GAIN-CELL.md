# Validating a gain interval in the energy norm

**The producer calculation passes:** for every relative gain \(\delta\in[-0.001,0.001]\), the specified finite-quadrature energy has a unique stationary point in its predictor-centered energy ball, and that point is a strict local energy minimum. The uniform self-map bound is below the radius by more than 23%. Every material coefficient remains fixed to the original weak-load case. No equilibrium optimization or new experimental fit was used. Independent review of this new certificate is a separate gate.

The [script](../scripts/energy_gain_cell.py) and [receipt](../evidence/energy_gain_cell.json) retain the source/input hashes, all eight predictor-path boxes, the whole error tube, analytic-domain margins and acceptance tests. The [norm feasibility study](PARAMETRIC-NORM-FEASIBILITY.md) explains why the volumetric rank structure matters and preserves the failed constant-remainder screen.

## Actual acceptance values

All 9,504 local terms and the global matrix contractions were evaluated with 160-bit Arb intervals. Exact rational outward endpoints are retained; the table gives safe rounded summaries.

| Quantity | Validated upper bound |
|---|---:|
| Center transformed Hessian defect | \(1.313\times10^{-10}\) |
| Center transformed force norm | \(6.783\times10^{-13}\) |
| Transformed predictor linear defect | \(6.599\times10^{-14}\) |
| Predictor metric norm | 0.017365397 |
| Eight-piece Hessian-variation integral | 0.035130960 |
| Uniform predictor residual \(\eta\) | \(6.100638\times10^{-7}\) |
| Uniform contraction bound \(Z\) | 0.155994806 |
| \(\eta+Z\rho\), with \(\rho=10^{-6}\) | \(7.660586\times10^{-7}\) |
| Root-to-predictor metric distance \(\eta/(1-Z)\) | \(7.228200\times10^{-7}\) |

The self-map margin exceeds \(2.3394\times10^{-7}\). The full error tube retains quadrature-domain lower bounds greater than 0.94284 for the meridional determinant, 0.99236 for hoop stretch and 0.99776 for the three-dimensional Jacobian. Surface radius and tangent norm also stay strictly positive. The cached quadrature weights and energy normalization are positive, as required for the PSD-majorant argument.

The eight successive predictor-path operator upper bounds are approximately 0.008112, 0.015817, 0.023530, 0.031249, 0.038976, 0.046711, 0.054452 and 0.062201. Averaging their exact upper endpoints makes the residual proof pass. Replacing the integral by the last bound had failed in the earlier feasibility study; that failure remains intact.

## The family and its proposed enclosure

The existing binary cache specifies a real-valued smooth energy, with 611 free displacement coordinates, for an original-material 12×6, strength-one state. Its finite volumetric coefficient is \(\kappa=5000\); that is not a measured physiological incompressibility prior. Let \(F(u,\delta)\) be its gradient when the existing illuminated drive is multiplied by \(1+\delta\). All other cached quantities remain fixed.

Let \(u_0\) be the archived binary free vector, \(S\) its archived binary gain tangent, and \(W\) the frozen binary metric from the separately validated norm study. The predictor and its proposed enclosure are

\[
 \widehat u(\delta)=u_0+S\delta,\qquad
 u=\widehat u(\delta)+Wz,\qquad
 \|z\|_2\le\rho=10^{-6}.
\]

This is a parameter-dependent tube, not a single nodal box. The tangent is not assumed exact: the center residual and the defect \(H_0S+b_0\) are reconstructed with Arb and retained in the residual bound.

## Exact predictor-residual identity

The force is affine in the relative gain. Write \(H_0=F_u(u_0,0)\) and \(b_0=F_\delta(u_0,0)\). Integrating along the predictor at fixed \(\delta\) gives

\[
 F(\widehat u(\delta),\delta)
 =F(u_0,0)+(H_0S+b_0)\delta
 +\int_0^1[H(u_0+tS\delta,\delta)-H_0]S\delta\,dt.
\]

Unlike the geometry displacement along this path, the parameter in its Hessian is \(\delta\), not \(t\delta\). Each path box therefore keeps the full gain interval.

For \(t\in[(i-1)/8,i/8]\), the proof bounds the geometry using \(|t\delta|\le ih/8\), with \(h=0.001\), and obtains an operator upper bound \(K_i\). It follows that

\[
 \eta=\|W^TF(u_0,0)\|_2+
 h\|W^T(H_0S+b_0)\|_2+
 h\|W^{-1}S\|_2\frac18\sum_{i=1}^8K_i
\]

bounds the transformed predictor residual uniformly over the gain interval. If \(\epsilon_0\) bounds \(\|W^TH_0W-I\|_2<1\), then the tangent metric norm is enclosed without a floating inverse:

\[
 \|W^{-1}S\|_2^2
 \le\frac{S^TH_0S}{1-\epsilon_0}.
\]

## Whole-tube derivative bounds

At every quadrature point, the linear strain map \(B\) and exact binary \(W,S\) give

\[
 |(B(u-u_0))_j|
 \le h|(BS)_j|+\rho\sqrt{(BWW^TB^T)_{jj}}.
\]

Every quantity here is enclosed with Arb. The pointwise determinant, hoop stretch and surface nondegeneracy margins must stay strictly inside the analytic domain. This is only quadrature-domain positivity; it does not certify the complete three-dimensional body geometry.

The proof preserves the bulk term \(\kappa\nabla\log J\,\nabla\log J^T\) as an outer product, controls its variation using a PSD squared-gradient majorant, and separately bounds the \((\kappa\log J-\mu)\nabla^2\log J\) term. It retains the directional contractions \(\nabla\log J\cdot BS\) before taking absolute values. The [norm study](PARAMETRIC-NORM-FEASIBILITY.md) supplies the detailed inequalities.

All resulting global majorants are transformed by the exact binary \(W\), with full Arb matrix multiplication. For these exactly symmetric enclosed matrices, a maximum absolute row sum bounds the spectral norm. No floating eigenvalue is used in acceptance.

## Acceptance theorem and its limits

For each fixed gain, use

\[
 T_\delta(z)=z-W^TF(\widehat u(\delta)+Wz,\delta).
\]

If the uniform derivative bound is \(Z<1\), and

\[
 \eta+Z\rho<\rho,
\]

then this map is a contraction of the closed Euclidean \(z\)-ball into itself. It supplies a unique stationary point in that gain's predictor-centered ball. Its metric distance from the predictor is at most \(\eta/(1-Z)\). Because the exact transformed Hessian is symmetric and within \(Z\) of the identity, it is positive definite throughout that ball; the stationary point is a strict local minimum of the specified finite-quadrature energy.

A passed derivative bound alone is insufficient. If the residual/self-map inequality fails, the receipt must remain unresolved. Even a successful calculation covers one gain coordinate at one fixed material and one numerical mesh. It does not identify material parameters, cover the full unbounded material domain or all seventeen revised-material states, certify a continuum solution or physical dynamics, resolve source observation semantics, or validate cross-size prediction. It supplies no new-general-concept claim.

The actual final source/input checks all pass. A guard-hardening limitation is retained explicitly: this producer revision selects its status from the numerical inequalities before combining the final unchanged-input/source booleans. Future reuse must check **all** final booleans, rather than trusting the status string alone. No source/input mutation occurred in this run. The successful source and receipt are preserved instead of silently changing that history.
