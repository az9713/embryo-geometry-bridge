# Base-point parameter derivative and predictor variation: preparation only

This frozen component is ready for independent prelaunch review. Preparation passes41 checks. It has computed no actual B-field polynomial moment defect, finite-cache force sum, parameter-variation value, or exact-integral root. The only integrations used in preparation are synthetic analytic toy checks. The actual run remains gated by a nonempty literal-PASS review bound to source, protocol and preflight hashes.

Let \(\Delta(u,\theta)=F_I(u,\theta)-F_Q(u,\theta)\) denote the exact-material-integral versus frozen binary-cache force difference, in the same611 free displacement coordinates. The preconditioner \(R\), exact base displacement \(u_0\), predictor \(S\), physical target, dose and parameter cell \(C\) remain fixed. The base parameters are taken from the independently reviewed SAFE complete point packet: requested exact strength divided by dose supplies the common gain, while the four material coefficients are the exact binary constants used by the point field. The separate binary-rounded convenience strength is not substituted.

## Conditional variation theorem

Assume a proved constant \(b\) satisfies

\[
\sup_{u\in B_{\rm independent},\ \theta\in C}
\lVert R\,\partial_u\Delta(u,\theta)\rVert_\infty\le b.
\]

Here \(B_{\rm independent}\) is the full canonical nodal coordinate box, independently combined with all of \(C\). A bound only on the coupled normal-root family \(u_0+S(\theta-\theta_0)+B_\infty(r)\) does not suffice: at fixed \(\theta\), its segment from \(u_0\) to the predictor may leave that smaller coupled domain. This preparation checks exactly that the independent convex box contains both endpoints of every such segment. All611 predictor-coordinate hulls agree with the fresh geometry factory. Their maximum magnitude \(M\) has exact rational provenance and approximate value \(5.684086827626595\times10^{-10}\).

The force at a fixed displacement is affine in the five varying coefficients. With \(p_\theta=u_0+S(\theta-\theta_0)\),

\[
\begin{aligned}
\Delta(p_\theta,\theta)-\Delta(u_0,\theta_0)
={}&\Delta(p_\theta,\theta)-\Delta(u_0,\theta)\\
&+\sum_j(\theta_j-\theta_{0j})\,\partial_{\theta_j}\Delta(u_0).
\end{aligned}
\]

Applying the fundamental theorem to the first term along the contained segment gives

\[
V\le bM+\sum_j\rho_j A_j,\qquad
\rho_j=\max_{\theta\in C}\lvert\theta_j-\theta_{0j}\rvert,
\quad A_j=\lVert R\,\partial_{\theta_j}\Delta(u_0)\rVert_\infty.
\]

The implementation will assemble signed interval derivatives first and multiply by the fixed exact-binary preconditioner. It will then bound each \(A_j\) and the parameter-only term. It deliberately does not accept a hypothetical bulk-only or unreviewed \(b\) as a completed \(V\) certificate. The full independent-box bulk-plus-surface derivative bound and all later root-reuse conditions remain separate.

## Surface derivatives without repeating square-root integration

At a fixed field, write \(G=J_s^2\). The surface potential is

\[
W_s=K\frac G2+(T-K)\sqrt G+\frac K2.
\]

The constant \(K/2\) has no force. Let \(P\) be the signed integral/cache difference of \(\nabla_u(G/2)\), and \(Q_s\) the corresponding difference of \(\nabla_u\sqrt G\). The already audited base-point edge force obeys

\[
F_0=K_0P+(T_0-K_0)Q_s.
\]

Both actual faces have \(T_0-K_0\ne0\), so interval division by this exact scalar gives

\[
Q_s=\frac{F_0-K_0P}{T_0-K_0},\qquad
\partial_K\Delta F=P-Q_s,\qquad
\partial_T\Delta F=Q_s.
\]

The run needs only polynomial \(P\): the exact SAFE jet supplies \(\nabla_uG/2\), exact material moments give its signed integral-minus-auxiliary-rule term, and an exact rational same-full-six-node cache correction gives the auxiliary-minus-cache term at \(K=T=1\). Normalized old cache maps include their exact radius/length scalings. The OLD positive normalized binary rule weight is applied once. Every actual point packet is source/hash bound; no packet midpoint replaces its interval. Existing square-root/tail uncertainty therefore remains in the recovered parameter derivatives. An equality \(K_0=T_0\) is rejected rather than divided through.

## Gain derivative

The active bulk term is \(\sigma(b^2+d^2)/2\), with \(\sigma=g\,I\,\chi_e\), exact dose \(I\), and the unchanged elementwise mask \(\chi_e\in\{0,1\}\). Therefore

\[
\partial_g\Delta F_e
=I\chi_e\,\Delta\nabla_u\frac{b^2+d^2}{2}.
\]

This is a polynomial derivative. Inactive elements contribute exactly zero. The24 active elements use signed exact polynomial moments plus exact rational coefficient maps evaluated on the SAME full nodal field for the finite auxiliary/cache correction. No constitutive logarithmic or square-root series is needed for these new derivatives.

## Guards and retained history

Synthetic checks include symbolic material/gain partials, both signs of \(T-K\), zero-denominator rejection, exact matched/mismatched polynomial cache examples, predictor hulls, ordered/nonfinite inputs, and false/truthy/mismatched review rejection. Base-point audit provenance, all24 surface maps, all611 endpoint hulls, common cell, cache and predictor hashes are frozen. Former unsafe packets appear only through historical lineage bindings; they are not restored as proof premises.

Preparation attempt0 is retained separately. It halted before actual field integration because it initially interpreted the audit's `exact_reduction_sha256` as the producer point packet. That key names the independent reduction. The corrected guard binds the producer through its audit path/hash entry and separately freezes the independent reduction. No reviewed source or packet was changed.

## Conditional stability lemma sanity check

For the parent's separate root-transfer theorem, suppose a frozen reference Hessian \(H_0\) is symmetric positive definite, \(\lVert I-RH_0\rVert<1\), and the new symmetric Hessian satisfies \(\lVert I-RH_I\rVert\le z<1\). The convex homotopy \(H_t=(1-t)H_0+tH_I\) remains nonsingular by the norm bound and Neumann invertibility. Its symmetric inertia cannot change, so \(H_I\) is positive definite. No symmetry assumption on \(R\) is required. Furthermore

\[
\lVert H_I^{-1}\rVert_\infty\le\frac{\lVert R\rVert_\infty}{1-z},\qquad
\lambda_{\min}(H_I)\ge\frac{1-z}{\lVert R\rVert_\infty}.
\]

This is a conditional finite-dimensional statement, not a completed new numerical stability certificate. It does not imply a physical trajectory, a global energy minimum, PDE stability, or source-matched embryo dynamics.
