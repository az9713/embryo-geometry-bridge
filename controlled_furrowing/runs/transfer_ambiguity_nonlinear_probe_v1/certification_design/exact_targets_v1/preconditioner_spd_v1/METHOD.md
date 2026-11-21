# Exact symmetric part of the frozen preconditioner

The target matrix is \(P=(R+R^\mathsf T)/2\), formed from the exact binary64 entries of the admitted preconditioner. Its positive definiteness is tested with an untrusted binary64 proposal \(K\). If
\[
\|I-KPK^\mathsf T\|_\infty<1
\]
holds for the exact matrices, the symmetric congruence has strictly positive eigenvalues. Its nonsingularity forces square \(K\) to be invertible, so \(P\) is SPD. The two outward evaluations use different product associations, \((KP)K^\mathsf T\) and \(K(PK^\mathsf T)\). The actual \(R\) never changes. A floating Cholesky failure or an insufficient interval bound leaves the result unresolved.

For later use only, suppose a separately reviewed, real symmetric Hessian \(H\) on the matched admissible cube has \(\|I-RH\|_\infty\le z<1\). Then \(A=RH\) is invertible and has spectrum in the open right half-plane. With \(S=H^{-1}=S^\mathsf T\), the identity
\[
AS+SA^\mathsf T=R+R^\mathsf T=2P
\]
follows from \(H=H^\mathsf T\). Because \(e^{-tA}\) decays, integrating the derivative of \(e^{-tA}S e^{-tA^\mathsf T}\) yields
\[
S=\int_0^\infty e^{-tA}(2P)e^{-tA^\mathsf T}\,dt\succ0 .
\]
Hence \(H\succ0\). The Neumann bound gives
\[
\|H^{-1}\|_\infty \le \frac{N_R}{1-z},\qquad
\lambda_{\min}(H)\ge\frac{1-z}{N_R},\qquad N_R\ge\|R\|_\infty .
\]
The last step uses that the spectral radius of the SPD inverse is at most its induced infinity norm. This is the classical linear-algebra argument already used for the earlier cohort. It requires the same \(R\), same symmetry, and valid whole-cube derivative bounds. This matrix phase evaluates no Hessians and certifies neither exact-spatial root existence nor a minimum. A later minimum proof also needs a fresh independently reviewed interior stationary root, \(C^2\) energy on its admissible neighborhood, and an appropriate normalized Hessian bound.

The cold preparation runs only synthetic matrices. Four serial production matrix checks require an independent prelaunch receipt, and any returned certificate still requires a separate independent actual readback. The native cached root and derivative receipts authenticate the matrix/target identities. They do not supply a Hessian sign.
