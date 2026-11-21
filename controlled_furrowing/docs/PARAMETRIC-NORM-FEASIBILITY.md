# A concrete energy norm for parameter-cell validation

**A norm change alone does not certify useful material cells.** This investigation does prove something narrower: the Hessian of the one archived weak-load, finite-quadrature model is positive throughout a displacement-coordinate cube of radius \(10^{-8}\). Preserving the volumetric rank structure also makes a 0.1% gain predictor tube plausible. That tube is still a feasibility result: its complete residual, self-map and directed-rounding gates have not passed here.

The [retained receipt](../evidence/parametric_norm_feasibility.json) separates the initial failed bounds, floating global feasibility screens, and an actual Arb certificate. The [read-only script](../scripts/parametric_norm_feasibility.py) performs no equilibrium optimization and changes no physical model. The target is the original 611-free-coordinate, 12×6, strength-one field used by the [root certificate](ALGEBRAIC-ROOT-CERTIFICATION-PATH.md), with \(\mu=1\) and finite volumetric coefficient \(\kappa=5000\). No result here covers the seventeen revised-material fields or an unbounded material domain.

## 1. What was making the old bound fail

Let \(H_0\) be the exact Hessian at the frozen center. The old componentwise estimate first takes absolute values of local Hessian changes, assembles row sums, then multiplies by \(|R|\), where \(R\) approximates \(H_0^{-1}\). This discards cancellations associated with finite-element compatibility and the stiff volumetric term. At nodal radius \(10^{-8}\), the comparable floating bound is 33.5451, so it cannot prove a contraction.

Choose a binary matrix \(W\) by numerically computing an inverse-transpose Cholesky factor of the archived Hessian. **Its numerical construction is not its proof.** Thereafter every entry of \(W\) is treated as an exact binary constant, and define

\[
 \|x\|_G=\|W^{-1}x\|_2,
 \qquad G=W^{-T}W^{-1}.
\]

The retained \(W\) is triangular with nonzero diagonal, so this is a genuine norm. Independently assembled Arb matrices verify the stronger useful condition

\[
 \|W^T H_0W-I\|_2
 \le \|W^T H_0W-I\|_\infty
 \le 1.312203388023\times10^{-10}.
\]

The norm inequality applies because the exact transformed matrix is symmetric. Thus positive definiteness is established by interval arithmetic, rather than assumed from a successful floating Cholesky factorization.

## 2. A validated positive Hessian on a whole displacement cube

At quadrature term \(q\), let \(B_q\) be the exact linear map from free coordinates to its local deformation variables, \(w_q>0\) its normalized weight, and \(h_q\) its small local Hessian. Arb encloses \(\Delta h_q=h_q(u)-h_q(u_0)\) throughout

\[
 \|u-u_0\|_\infty\le10^{-8}.
\]

Let \(m_{ij}\) bound both \(|\Delta h_{ij}|\) and \(|\Delta h_{ji}|\), and set \(D_{ii}=\sum_jm_{ij}\). Symmetric diagonal dominance gives

\[
 -D_q\preceq\Delta h_q\preceq D_q,
 \qquad
 -M\preceq H(u)-H_0\preceq M,
 \qquad M=\sum_q w_qB_q^TD_qB_q.
\]

All 9,504 terms are included. Full Arb matrix multiplication, with the frozen binary \(W\), proves

\[
 \|W^TMW\|_2\le\|W^TMW\|_\infty
 \le0.072649381762270.
\]

Consequently, everywhere in that cube,

\[
 W^TH(u)W\succeq0.927350618106510\,I.
\]

These displayed decimal constants are safe rounded summaries; the receipt retains exact rational interval endpoints. The earlier root certificate places its unique local stationary point strictly inside this cube. The energy is therefore strictly convex on the cube, and that stationary point is a **strict local energy minimum of this frozen finite-quadrature potential**. This is stronger than stationarity alone. It does not establish physical dynamical stability, a global minimum, branch uniqueness outside the cube, or continuum accuracy.

| Bound on the same \(10^{-8}\) nodal cube | Value | Status |
|---|---:|---|
| Absolute-value inverse/row majorant | 33.5451 | Floating screen; fails |
| Generalized eigenvalue of assembled PSD majorant | 0.0208433 | Floating screen |
| Transformed majorant maximum absolute row sum | 0.072649381762270 | Arb validated |
| Center defect plus transformed majorant | 0.072649381893490 | Arb validated; below one |

The smaller generalized eigenvalue is not used as a certified number.

## 3. Why naive energy balls still fail

A metric ball \(\|W^{-1}e\|_2\le\rho\) gives local widths

\[
 |(B_qe)_i|\le\rho\,\|(B_qW)_{i,:}\|_2.
\]

The maximum nodal conversion constant is approximately 14.05555: \(\|e\|_\infty\le14.05555\rho\). This is a floating diagnostic until its norms are enclosed outward. A metric radius \(10^{-4}\), processed with the same diagonal local majorant, has generalized majorant bound about 120.84. Thus a favorable result on a small nodal cube cannot be promoted to an arbitrary larger energy ball.

Likewise, a 0.1% gain predictor tube of metric radius \(10^{-7}\), treated without volumetric structure, has generalized bound 1.9116 and fails. Those failures are retained.

## 4. Preserve the actual volumetric structure

For the five local bulk variables, write \(\ell=\log J\), \(v=\nabla\ell\), and \(Q=\nabla^2\ell\). The local bulk Hessian is

\[
 h=\mu I+gD+\kappa vv^T+(\kappa\ell-\mu)Q,
 \qquad D=\operatorname{diag}(0,1,0,1,0),
\]

with the local illuminated drive included in \(g\). Its gain Hessian is constant and positive semidefinite; interval subtraction of two full volumetric Hessians is unnecessary for that contribution.

For a displacement predictor \(\widehat u=u_0+S\delta\), first compute \(B_qS\) and contractions with \(v_0,Q_0\), **then** take absolute values. In particular,

\[
 |\ell(q_0+d)-\ell(q_0)|
 \le |v_0^Td|+
 \tfrac12 |d|^T\,\overline{|Q|}\,|d|.
\]

The first term preserves the nearly volume-preserving predictor direction. Replacing it with \(|v_0|^T|d|\) loses that information. Similarly,

\[
 |v(q_0+d)-v_0|
 \le |Q_0d|+\overline{|Q-Q_0|}\,|d|.
\]

For all bulk points, stack the weighted rows \(\sqrt{w_q\kappa}\,v_0^TB_q\) into \(A\), and their variations into \(D_*\). Then the volumetric outer-product variation satisfies

\[
 \|W^T\Delta H_{\rm rank}W\|_2
 \le2\|AW\|_2\|D_*W\|_2+\|D_*W\|_2^2.
\]

If local gradient changes obey \(|\delta v_i|\le d_i\), a PSD majorant for the squared variation is obtained from

\[
 (\delta v^Ty)^2\le
 \left(\sum_i d_i\right)\left(\sum_i d_i y_i^2\right).
\]

The remaining bulk change is bounded separately using

\[
 |\Delta[(\kappa\ell-\mu)Q]|
 \le |\kappa\ell_0-\mu|\,\overline{|Q-Q_0|}
 +\kappa\,\overline{|\ell-\ell_0|}\,\overline{|Q|}.
\]

Surface terms retain their exact local interval Hessian majorants. This decomposition avoids imposing a full \(\kappa\)-scale diagonal perturbation on every soft shear direction.

## 5. Concrete predictor feasibility, with the unpassed gates exposed

The predictor uses the archived gain tangent. Its metric norm is approximately 0.0173654, compared with nodal infinity norm 0.0451285. The following are **floating global norm screens using Arb local derivative enclosures**, not validated material cells:

| Relative gain half-width | Metric tube radius | Naive generalized bound | Rank-aware generalized bound | Rank-aware absolute-row bound |
|---|---:|---:|---:|---:|
| 0.1% | \(10^{-7}\) | 1.9116 | 0.0196733 | 0.0692105 |
| 1% | \(10^{-5}\) | 24.4061 | 0.491516 | 1.60627 |

The first row supports a specific next implementation, rather than the generic advice to change norms. The second row still fails the cheaper absolute-row acceptance test. A floating eigenvalue below one cannot repair that failure.

For a genuine uniform cell, freeze binary \(W,S\) and enclose every contraction and norm outward. In transformed coordinates the map is

\[
 T_\theta(z)=z-W^TF(\widehat u(\theta)+Wz,\theta).
\]

One must prove, uniformly over the declared coefficient cell,

\[
 \eta\ge\sup_\theta\|W^TF(\widehat u(\theta),\theta)\|_2,
 \quad Z\ge\sup_{\theta,\|z\|_2\le\rho}
 \|I-W^TH(\widehat u(\theta)+Wz,\theta)W\|_2,
 \quad \eta+Z\rho<\rho.
\]

Neither a derivative screen \(Z<1\) nor the old one-point root certificate supplies the missing predictor residual \(\eta\). The archived tangent solves the linear system only numerically, so its exact defect must also be retained. For gain half-width \(h\), with gain-force derivative \(b_g\) and constant gain Hessian \(H_g\), a useful residual inequality is

\[
 \eta\le\|W^TF_0\|_2+
 h\|W^T(H_0S_g+b_g)\|_2+
 h\|W^{-1}S_g\|_2\int_0^1 K_{\rm geom}(th)\,dt+
 h^2\|W^TH_gS_g\|_2.
\]

Here \(K_{\rm geom}(th)\) bounds the geometry-only transformed Hessian change along the predictor segment; it excludes the separately handled gain contribution. Subdividing that one-dimensional integral can sharpen the remainder. This inequality is proved by integrating the gradient along the predictor and using the exact affine dependence on gain. Its inputs have **not** all been validated in this investigation, so no useful-width cell is claimed.

A final screen makes the remaining obstruction quantitative. At gain half-width \(h=0.001\), the rank-aware absolute-row bound along the predictor alone is 0.0622008; adding a metric error tube \(\rho=10^{-6}\) gives 0.155995. Bounding the whole Taylor integral by its endpoint bound yields nonlinear residual about \(1.08014\times10^{-6}\), hence

\[
 \eta_{\rm remainder}+Z\rho\approx1.23614\times10^{-6}>\rho.
\]

That conservative screen **fails**, even before adding the small center residual and tangent defect. The receipt retains the failure. It identifies a specific next computation: validated subdivision of the one-dimensional remainder integral, or a sharper validated symmetric norm bound. The floating spectral estimates alone would be more favorable, but are not used to declare a cell.

The archived columns named `log_*` are derivatives at the center. They also equal derivatives with respect to a *relative linear coefficient change* there. The linear predictor above uses that latter coordinate. A uniform logarithmic parameter box must instead enclose the exponential conversion, not silently identify a finite log change with a relative linear change. The zero apical-tension coordinate is additive and physically one-sided.

## 6. Fixed-center coefficient information and the force-scale quotient

At a fixed geometry the Hessian is affine in the mechanical coefficients. The floating generalized absolute norms of its center-scaled coefficient derivatives are:

| Coordinate | Relative Hessian norm | Predictor metric norm |
|---|---:|---:|
| Relative \(\kappa\) | 0.999757 | 0.00414478 |
| Relative apical area modulus | 0.422269 | 0.000680994 |
| Relative basal area modulus | 0.829751 | 0.00142934 |
| Additive apical tension | 14.8063 | 0.0261157 |
| Relative basal tension | 0.960989 | 0.00358246 |
| Relative common gain | 0.136952 | 0.0173654 |

These are numerical diagnostics for this center, not global Lipschitz constants. The additive apical tension is particularly consequential in these coordinates.

Common positive scaling of every force/energy coefficient leaves an unnormalized equilibrium equation unchanged up to a scalar. Fixing \(\mu=1\) therefore removes that gauge and leaves dimensionless ratios; it does not remove the potentially large ratio \(\kappa/\mu\), or bound the remaining ratios. An energy-norm calculation at numerical \(\kappa=5000\) is not a result for an independently measured physiological compressibility and is not a uniform incompressible limit.

## 7. Actionable outcome

The validated positive-Hessian cube supplies a real local stability prerequisite. The rank-aware inequalities identify a concrete route to parameter-dependent tubes: validate the 0.1% gain case first, including its residual; then add coefficient directions with their exact affine derivatives, all failed boxes retained. A fixed-center nodal box and a parameter-dependent predictor tube are different sets and need separate proofs.

There is no justification for declaring that a finite set of local cells covers the full unbounded material domain. That requires either a defensible compact physical domain or proved asymptotic exclusions and a covering argument. No bounded material prior, continuum error bound, experiment-to-observer identification, cross-size validation, or novelty claim is supplied by this norm analysis.

The completed Arb calculation can be replayed without overwriting its receipt:

```powershell
$env:OPENBLAS_NUM_THREADS='1'
python -B controlled_furrowing/scripts/parametric_norm_feasibility.py --replay-nodal
```

Replay recomputes the local interval enclosures and full transformed matrices, then compares the exact stored bounds and frozen metric checksum.
