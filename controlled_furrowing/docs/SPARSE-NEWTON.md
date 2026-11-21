# Sparse automatic differentiation and the 24×12 result

**Status: the sparse solver component and single 24×12 continuation calculation are completed and verified within their discrete scope. The matched 16×8→24×12 comparison still exposes substantial shape sensitivity; continuum accuracy is unresolved.**

The activated mesh study exposed unresolved shape changes. Increasing resolution with a fully vectorized dense AD Hessian would multiply its seed dimension with every quadrature intermediate. This component reduces that computational burden while preserving the declared high-quadrature mechanical energy and its gradient. It changes the numerical Newton policy, and does not claim identical optimizer trajectories, new physics or a continuum accuracy certificate.

`scripts/sparse_quadratic.py` is a separate subclass of the fixed high-quadrature model. It has **no energy or gradient override**. The actual integration remains Duffy Gauss order 8 in each transformed-triangle direction, surface Gauss order 12, and twelve triangles per energy chunk. The full pointwise logarithmic-volume penalty, active potential, surface elasticity and baseline tension remain unchanged.

## Established method, locally verified implementation

Graph coloring for compressed derivative evaluation is established numerical analysis. [Coleman and Moré (1983)](https://epubs.siam.org/doi/abs/10.1137/0720013?journalCode=sjnaam) formulate sparse Jacobian estimation through coloring; their [1982 Hessian technical report](https://ecommons.cornell.edu/entities/publication/6c14f4fa-13f2-450b-b19d-362ef1e2e9a5) treats recovery from matrix-vector products. [Gebremedhin, Manne and Pothen (2005), Section 3.5.1](https://cscapes.cs.purdue.edu/coloringpage/abstracts/sirev.pdf), describes greedy distance-2 coloring within a broader derivative framework. These sources support the established algorithmic context; no originality claim is made. Query/access receipts are retained in `evidence/sparse_quadratic_literature.json`.

The proof below explains the exact recovery rule used here, so the reader need not trust a software label or infer correctness from a speed comparison.

## Why the physical Hessian is sparse

Let \(u\) contain the free displacement coordinates and denote its nodal Hessian by \(\mathcal H\), reserving \(H\) for reference thickness. The high-quadrature potential is a sum of local contributions,

\[
 E(u)=\sum_T E_T(u_{I_T})+\sum_e E_e(u_{I_e}),
\]

where a P2 triangle has at most twelve displacement degrees of freedom and a quadratic boundary edge has at most six. If two free coordinates never occur together in any local term, their mixed second derivative vanishes wherever the physical potential is twice differentiable. The union of these element/edge cliques is therefore a conservative structural Hessian pattern.

This locality depends on the actual energy. The projected pressure in the full-penalty model is a detached diagnostic and does not enter the total physical potential. A governing global projected-volume term would require a different sparsity analysis. The implementation checks the present pattern against every entry of a small dense AD Hessian, including all entries predicted to vanish.

## Coloring and exact compressed recovery

For row \(i\), denote the structural column support by \(S_i\). Choose a color map \(\chi\) such that no two columns in any \(S_i\) have the same color. Greedy coloring of the resulting column-conflict graph supplies this distance-2 collision rule. It need not minimize the number of colors; row uniqueness is the necessary correctness property and is checked explicitly.

For color \(c\), define the seed vector \(v^c_j=1\) if \(\chi(j)=c\), and zero otherwise. Then

\[
 (\mathcal H v^c)_i=\sum_{j\in S_i,\,\chi(j)=c}\mathcal H_{ij}.
\]

The sum contains at most one structural entry. Consequently

\[
 \mathcal H_{ij}=(\mathcal H v^{\chi(j)})_i\qquad(j\in S_i).
\]

This is direct recovery in exact arithmetic, without a finite-difference Hessian or a separately coded constitutive tangent. Batched reverse-mode AD differentiates the actual energy gradient. Since the energy Hessian is symmetric on the smooth physical branch, its transposed-vector convention gives the same products. The implementation records raw asymmetry and uses the symmetric average for the linear solve.

Seeds are evaluated in chunks of four colors, bounding the AD batch dimension independently of the total free-coordinate count. The recovered matrix is stored as CSR. Pattern, coloring, their SHA256 and all source hashes are retained.

At 24×12 cells there are 2,375 free displacement coordinates, 52,005 structural entries and 46 colors. Four-color chunks replace a 2,375-seed dense AD batch. The raw dense matrix alone would require 45,125,000 bytes; sparse values plus 64-bit indices/pointers require approximately 851,088 bytes. Those matrix-storage figures do not include AD graphs, sparse LU fill, quadrature caches or the rest of the process. Actual process memory is measured separately.

## A different, safeguarded Newton policy

At a fixed activation strength, the sparse policy solves

\[
 (\mathcal H+\tau I)d=-g.
\]

It first tries \(\tau=0\). Failed factorizations or non-descent directions trigger increasing diagonal shifts. A finite direction must have a strictly negative computed gradient inner product. Its initial step is capped by a 0.25H coordinate-displacement limit; backtracking then requires the Armijo condition and a sufficient positive-orientation check.

The exact Armijo inequality is \(E(u+\alpha d)\le E(u)+10^{-4}\alpha g^Td\). The code includes the explicit \(10^{-16}\) normalized-energy allowance inherited from the numerical scale of these calculations; strict exact-arithmetic descent is not asserted for changes below that floating allowance. Sparse linear residuals, shifts, directional inner products, accepted fractions and failures are recorded. The linear residual is diagnostic rather than a separate convergence gate.

This policy differs from the earlier dense eigenvalue-based regularization. Matching the energy and gradient does not make the optimization path identical, prove a global minimum or exclude selection of a different local stationary branch.

The direction test alone also does not certify that the unshifted Hessian is positive definite. A final field passing force/geometry checks is a stationary discrete equilibrium; local energetic stability would require an additional curvature/inertia check. It is not silently promoted to a certified local minimum.

The outer [verified-state continuation controller](QUADRATIC-SOLVE-POLICY.md) is unchanged. Accepted states still require normalized force below \(10^{-7}\) and positive geometry. A failed trial is retained and discarded as a warm start; activation increments are bisected from the last verified state. A computational resource improvement cannot replace these physical/numerical checks.

## Verification before any 24×12 solve

The executable protocol fixes tolerances before the checks. On a 6×3 mesh the structural pattern has 161 free coordinates, 3,045 entries and 42 colors. The entire recovered matrix was compared with a dense automatic-differentiation Hessian, rather than checking only selected entries.

| Check | Retained outcome |
|---|---|
| Every structural row has unique column colors | PASS |
| Symmetric structural pattern | PASS |
| Full recovered H versus dense AD | Relative difference 1.06×10⁻¹⁶; absolute 4.55×10⁻¹³ |
| Every dense entry outside the pattern | Exactly zero in the test |
| Raw recovery asymmetry | 4.55×10⁻¹³ |
| Three independent AD Hessian-vector comparisons | Relative differences below 5.4×10⁻¹⁶ |
| Three gradient finite-difference H-vector checks | Relative differences below 9.0×10⁻¹⁰ |
| Energy/gradient equality at gains 0, 1 and 5 | Bit-for-bit identical |
| Activated sparse-policy solve | Force/positive-geometry gate passed |

All fifteen checks pass in `evidence/sparse_quadratic_verification.json`. The full dense/recovered matrices, pattern, vector and colors are retained in `runs/sparse_quadratic/small_hessian_comparison.npz`.

An independent AI reviewer repeated fifteen checks on a different actual-small-geometry parameter case at volume coefficient 5000 and activation 10. It independently rebuilt the cliques and verified row colors, all matrix entries, independent AD/finite-difference products and energy/gradient identity. Its full-matrix relative difference was 7.93×10⁻¹⁷; the receipt is `evidence/joint_refinement_sparse_review.json`. This is internal independent implementation review, not outside expert validation or formal verification of the continuum model.

After the live calculation reached a verified 24×12 gain-0.5 field, a separate derivative-only check recovered its sparse Hessian and compared three new products with independently differentiated AD products and gradient finite differences. Relative errors were below 4.2×10⁻¹⁶ and 2.2×10⁻⁹ respectively. The input field, matrix and receipt are retained in `runs/sparse_quadratic/hessian24_at_verified_gain0p5.npz` and `evidence/sparse_quadratic_24_hvp_check.json`. This was not a second optimization job and did not change the running state.

## Resource evidence is distinct from physical accuracy

Two separate Python processes benchmarked only the small-mesh Hessian operation:

| Method | Time | Sampled incremental process RSS |
|---|---:|---:|
| Dense vectorized AD | 1.276 s | 114,323,456 bytes |
| Colored sparse AD, four seeds per chunk | 1.155 s | 66,830,336 bytes |

Separate processes reduce allocator-history confounding. RSS is sampled every 5 ms and includes native Torch allocations, but can miss brief peaks; it is not a strict memory bound. These small-mesh numbers are not extrapolated into a claimed 24×12 speedup. Host/software variability also limits comparisons. Benchmarks and measurement definitions are retained separately in `evidence/sparse_quadratic_benchmark_dense.json` and `sparse_quadratic_benchmark_sparse.json`.

![Sparse automatic differentiation, compressed products and isolated-process memory measurements](../runs/sparse_quadratic/sparse_hessian_verification.png)

The first two panels display the same small Hessian and its colored compressed products. The memory panel concerns those small-mesh operations only. Figure/input hashes are retained in `evidence/sparse_quadratic_visual.json`.

## One declared 24×12 calculation

Only after both verification gates passed was the single 24×12 job launched. Its frozen source SHA256 is `d4833318151c7b94b4c56d40a96158d55a06258619ef080548d489e2be0394f9`. A matching verification receipt is required by the executable; an existing job record prevents duplicate launch.

The calculation fixes volume coefficient 5000, basal/apical tensions 2/0, basal/apical area moduli 2/0.2, \(H=33.5\,\mu\mathrm m\), activation radius 9.35512037 µm, depth 5–25 µm below the top and outer radius \(3H\). It attempts gains 1, 5 and 10 through verified continuation. No response data, refitted gain or changed physical model enters this job.

`evidence/sparse_quadratic_24_protocol.json` retains inputs, integration policy, sparse policy, pattern and all source hashes. Every full accepted/rejected field remains under `runs/sparse_quadratic/24x12_lambda5000/`. `sparse_quadratic_24_results.json` records target completion, forces, actual sampled J, conservative geometric bounds, attempt counts and per-target resource measurements. A separate agent handles the 16×8 high-quadrature fixed-load state; it is not duplicated here.

### Completed targets and retained failures

All three requested targets reached the unchanged force/positive-geometry gates. The single process ran approximately 40.2 minutes across all target attempts; sampled peak process RSS was 668,655,616 bytes (about 0.623 GiB). Sampling every 5 ms can miss brief peaks, so this is observed resource evidence, not a strict memory bound. Per-target increments are relative to an already allocated process and are not independent total-memory requirements.

| Gain | Final normalized force | Maximum column shortening | Actual Q8 sampled J | Attempts / rejected | Time |
|---|---:|---:|---|---:|---:|
| 1 | 2.98×10⁻¹² | 2.33012 µm | [0.997521, 1.001900] | 3 / 1 | 637.7 s |
| 5 | 3.86×10⁻¹² | 7.39158 µm | [0.974139, 1.014963] | 5 / 2 | 1128.2 s |
| 10 | 4.64×10⁻¹² | 10.68165 µm | [0.960950, 1.024151] | 3 / 1 | 646.6 s |

The accepted sequence was 0→0.5→1→2→3→5→7.5→10. Direct trials 0→1, 1→5, 1→3 and 5→10 failed the force gate and remain retained. The independent full-path replay checked all eleven fields, the seven accepted states, four rejected states, last-verified vector hash chains, energy/work accounting and frozen physical-source hashes. Energy and gradient evaluations are bit-for-bit identical to the high-quadrature parent. The receipt is `evidence/sparse_quadratic_24_verification.json` and executable replay is `runs/sparse_quadratic/verify24.py`.

The retained path contains 227 accepted Newton steps and 44 nonzero diagonal-shift trials. These are numerical policy diagnostics, not a count of distinct physical transitions. The largest recorded linear-system residual is 4.01×10⁻¹². A first rejected trial used the entire 25-iteration Newton allowance; continuation reduced activation increments instead of weakening the acceptance gate.

Higher quadrature was checked without reoptimizing. At gains 5 and 10, Duffy order 12 and surface order 16 give force residuals 3.80×10⁻¹² and 4.37×10⁻¹². `evidence/sparse_quadratic_24_crosscheck.json` retains these evaluation-only checks. The independent review of gain 1 is separately scoped; no duplicate optimization was launched.

For the final gain-10 field, exact rational Bernstein calculations on its stored binary coordinates certify every one of the 576 triangles and give the whole-body enclosure

\[
0.95840970995 < J < 1.02823683506.
\]

The exact fractions, boundary test and global embedding certificate are retained in `evidence/sparse_quadratic_24_final_geometry.json`. This encloses the stored discrete map, not the unknown continuum solution. Its coefficient bounds are distinct from sampled J extrema; neither is a validated biological incompressibility tolerance. The initial orientation checks in the solve use a different, more conservative sufficient certificate.

A derivative-only profiler sample on the saved gain-0.5 field recorded about 171 ms of instrumented CPU operator time for one warmed energy/gradient evaluation. Work was spread across tensor arithmetic, reductions and copies; no single determinant-kernel bottleneck was established. `evidence/sparse_quadratic_gradient_profile.json` retains this performance diagnostic. It changes no state and is not a whole-job timing prediction.

### Matched high-quadrature shape comparison

The comparison uses the completed 16×8 high-quadrature field in `runs/joint_refinement/final.json` and the completed 24×12 gain-10 field. Constitutive parameters, activation, geometry and actual Q8/surface-12 policy match; only mesh dimensions differ. No gain was refitted. Their initial numerical guesses and continuation histories differ, so identical local-branch selection is not established.

The shortening changes by only 0.0339406 µm, from 10.6477085 to 10.6816491 µm. Nevertheless exact rational inversion of their P2 boundary curves on the common **current physical radius** gives

\[
\begin{aligned}
\sup_R |a_{24}(R)-a_{16}(R)|&\in[0.670904,\,0.680214]\;\mu\mathrm m,\\
\sup_R |b_{24}(R)-b_{16}(R)|&\in[0.123578,\,0.133086]\;\mu\mathrm m.
\end{aligned}
\]

The displayed endpoints are rounded outward from the exact receipt. These are enclosures of differences between two stored discrete graphs. They are **not** continuum PDE error bounds, and they do not prove mesh convergence. Agreement of a scalar summary can hide appreciable spatial deformation differences. `evidence/sparse_quadratic_24_comparison.json` retains the exact interval arithmetic, common domain, field hashes and sampled display curves; `runs/sparse_quadratic/compare24.py` replays it.

![Matched gain-10 high-quadrature apical and basal curves, with their Eulerian differences](../runs/sparse_quadratic/24_vs16_highQ.png)

The plotting samples illustrate the curves; the interval algorithm, rather than the sampled maximum, supplies the reported sup enclosure. The first visual had a clipped long panel title and is retained as `24_vs16_highQ.initial_clipped.png`; only the title was shortened. The image and correction hashes are in `evidence/sparse_quadratic_24_comparison_visual.json` and `sparse_quadratic_24_visual_correction.json`.

### Numerical curvature diagnostic at the final equilibria

Stationarity and positive orientation do not imply a minimum. A separate evaluation-only diagnostic therefore recovers the **unshifted** energy Hessian at both saved gain-10 fields. Symmetric LDL factorization has positive numerical inertia in all 1071 free coordinates at 16×8 and all 2375 at 24×12, with relative reconstruction residuals below 1.1×10⁻¹⁵. No negative direction was found, so the conditional negative-direction energy replay was not invoked.

The three shift-invert eigenpairs nearest zero were checked against the recovered matrix and, for the smallest pair, an independently differentiated Hessian-vector product. Relative residuals are below 3.5×10⁻⁹. The smallest reported ordinary-coordinate eigenvalues are approximately 6.2831×10⁻⁴ and 2.6126×10⁻⁴. Numerical LDL inertia, rather than a positive near-zero eigenpair alone, supplies the diagnostic sign count. Protocol, factor-block signs, thresholds and matrices/eigenvectors are retained in `evidence/sparse_equilibrium_stability_diagnostic.json` and `runs/sparse_quadratic/stability_{16x8,24x12}.npz`.

These are floating-arithmetic diagnostics of the declared finite-dimensional, quadrature-defined energy. They are not rigorous spectral enclosures or continuum stability proofs. Ordinary Euclidean-DOF eigenvalue magnitudes depend on basis and mesh normalization: the smaller 24×12 value does **not** establish physical softening or stability loss. Inertia signs are invariant under nonsingular basis congruence, whereas physical cross-mesh eigenfrequencies or comparable curvature magnitudes require a declared mass/metric. No global minimum, unique equilibrium or identical branch is established.

Sparse recovery removed a computational obstacle and enabled a finer mechanically coupled calculation. The result strengthens the numerical audit, but leaves the scientific gate open: spatial convergence, material identification, transient observation mapping and biological prediction remain unresolved. It establishes no new biological concept.
