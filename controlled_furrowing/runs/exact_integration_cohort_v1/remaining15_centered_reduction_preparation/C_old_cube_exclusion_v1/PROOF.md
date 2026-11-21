# Candidate exclusion for record C’s original predictor cube

**Exact comparison passed. Independent review is pending.**

Let \(F_I(u,\theta)\) be the normalized exact-spatial-integration force of the fixed P2 mechanical model and \(F_Q(u,\theta)\) the matching cached quadrature force. Both use the selected source dose \(4089863/12500000\). The original predictor is \(p_\theta=u_0+S(\theta-\theta_0)\) in 611 normalized free coordinates. For every parameter \(\theta\) in the same closed cell, the original cube is \(\|u-p_\theta\|_\infty\le r\), with \(r=10^{-11}\).

The independently replayed signed predictor intervals bound \(R(F_I-F_Q)(p_\theta,\theta)\) over that entire parameter cell. Row 6 in zero-based free-coordinate order lies in \([2.21627353655659e-11,2.94959897880713e-11]\) to display precision, strictly away from zero. Its distance from zero implies \(\|R(F_I-F_Q)(p_\theta,\theta)\|_\infty\ge d_*\) for every \(\theta\). The cached predictor residual has the uniform upper bound \(\|RF_Q(p_\theta,\theta)\|_\infty\le\eta_Q\). The reverse triangle inequality gives
\[
\|RF_I(p_\theta,\theta)\|_\infty\ge d_* - \eta_Q.
\]

Suppose \(u_\star\) is a zero of \(F_I\) in that same cube at the same \(\theta\). The line from \(p_\theta\) to \(u_\star\) stays inside the reviewed domain. The full integral Hessian \(H_I\) includes the retained 24 surfaces, and the reviewed whole-box bound is \(\|I-RH_I\|_\infty\le z\). Hence \(\|RH_I\|_\infty\le1+z\), and integrating the derivative on the line yields
\[
\|RF_I(p_\theta,\theta)\|_\infty
=\|R(F_I(p_\theta,\theta)-F_I(u_\star,\theta))\|_\infty
\le(1+z)\,r.
\]

The exact comparison in [candidate.json](candidate.json) gives \(d_* - \eta_Q\approx 2.16880141206288e-11\), while \((1+z)r\approx 1.08524938827999e-11\). The exact difference is strictly positive, approximately \(1.08355202378289e-11\). Thus the two inequalities are incompatible for every \(\theta\) in the full originally declared closed cell.

This exclusion argument uses a positive residual lower bound. The separate centered sufficient test keeps its unresolved status. The radius is a tolerance in normalized free coordinates. This argument gives no existence or exclusion result elsewhere, no source-data or biological incompatibility, and no rejection of the material law. The candidate uses the prior reviewed force and derivative enclosures and performs no new field computation.
