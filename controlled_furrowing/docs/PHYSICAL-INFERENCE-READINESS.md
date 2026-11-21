# What physical inference is currently justified?

**The public data support a conditional inference problem, but not a uniquely specified biological parameter set or a completed physical identified set.** We can state the domains, reject some exact contracts and calculate candidate discrepancy requirements. The missing observation semantics and forward error bounds cannot be replaced by a better optimizer.

This read-only audit uses the [physical domain](PHYSICAL-PARAMETER-DOMAIN.md), [source observation contract](JOINT-OBSERVATION-CONTRACT.md), [material-input audit](MATERIAL-INPUTS.md) and [completed direct replay](JOINT-DIRECT-REVIEW.md). The [receipt](../evidence/physical_inference_readiness.json) retains source hashes, executed rational-arithmetic checks and a replay program. No larger-region outcomes were read, fitted or used to select assumptions. The original [five-gate plan](../PLAN.md) remains unchanged.

## 1. Separate the four domains before fitting

Let \(w\) denote scale-free material parameters, \(\gamma_i\) specimen geometry, \(\eta\) a cohort-wide observation/branch convention, and \(\delta\in[0,\infty)^4\) four discrepancy allowances. The intended object is

\[
\Theta_\delta(\eta)=\left\{(w,\gamma_1,\ldots,\gamma_{17}):
\left|\mathcal O_k[q_i(w,\gamma_i,I_i);\eta]-y_{ik}\right|
\leq\delta_k\quad\text{for every }i,k\right\}.
\]

Here \(q_i\) must be an actual solution/path under the stated mechanical model and branch rule. Substituting independently adjustable four-vectors for \(q_i\), or letting each record use a convenient different observation convention, changes the inference problem.

| Domain | What the source or mathematics establishes | What would be an added assumption |
|---|---|---|
| Material ratios | Absolute common force scale is invisible to quasistatic shape, under the stated energy and loading assumptions. | Finite upper/lower ratio bounds; transferring another tissue's spring constants; fixing basal stiffness or prestress to a measured value. |
| Height and activation support | Methods Sec14 selects embryos at cell lengths 32–35 µm and activates 5–25 µm below the apical surface. | Exact per-record height, uniform height throughout each patch, or a uniform probability distribution on that interval. |
| Small-region diameter | Nominal 18 µm in the experimental description; 18.71024074 µm in the workbook. | Calling these two conventions endpoints of a measured calibration-error interval. |
| Other geometry | Measurements use coverslip-adjacent tissue; the illustrative images constrain some visible features. | Flatness, axisymmetry, an outer domain radius, curvature bounds, a support law or a verified cell-to-continuum scale separation. |
| Input | Sec17 defines normalized fluorescence from the first activation acquisition, relative to a shared depthwise standard across experiments. | A measured stress, a complete force waveform, exact stress proportionality, or measured independent dose-error intervals. |
| Observation error | The source describes segmentation, correction, PIV registration and reference growth. | A numerical annotation bound or a common independent additive noise law. Pixel spacing and acquisition interval are not those bounds. |

The independently audited 17 records are distinct biological samples; absence of raw-file IDs does not imply dependence. Shared fluorescence normalization and shared image-derived summaries do mean that independence of all errors cannot be inferred from the sample count.

## 2. Material domain: an exhaustive representation, not a convenient box

With fixed length unit \(L_0=33.5\) µm, normalize

\[
c=(\mu,\kappa,K_a/L_0,K_b/L_0,T_a/L_0,T_b/L_0,g),
\qquad w=c/\sum_jc_j.
\]

The currently justified constitutive domain is

\[
\mathcal W=\{w\geq0:\ \sum_jw_j=1,\quad
w_\mu>0,\ w_\kappa>0\}.
\]

It has six independent ratios. None of the source audits supplies a measured finite box exhausting them. If \(\kappa/\mu=r\) is fixed as a numerical incompressibility approximation, impose \(w_\kappa=rw_\mu\) and distinguish convergence in that approximation from inference of a physiological compressibility. Fixing it removes one coordinate by assumption, rather than by the four-summary data.

A mathematically exhaustive computational decomposition is available:

\[
\mathcal W_n=\{w\geq0:\sum_jw_j=1,\quad
w_\mu,w_\kappa\geq1/n\},\qquad
\mathcal W=\bigcup_{n\geq2}\mathcal W_n.
\]

Each \(\mathcal W_n\) is compact. Every admissible vector enters some member because its two required positive coordinates have a positive minimum. The decomposition does **not** prove that any finite \(n\) covers the data-compatible set. A finite exhaustive inference needs a data-derived tail exclusion or an analytic bound. A grid, coordinate descent, or several successful fits inside one member provides candidates; it cannot silently discard all larger members. Degenerating stiffness ratios can also destroy uniform forward stability near excluded simplex faces.

An analogous issue remains for geometry: protocol height ranges alone do not bound initial curvature, outer support or reference trajectories. Those variables need measured constraints or a clearly conditional model family. Choosing flat geometry with a common height is a legitimate hypothesis test; it is not a specimen-derived geometry distribution.

## 3. Discrepancy: infer a family rather than manufacture an error bar

There is currently no evidence-based unique \(\delta\). The useful executable object is the family \(\{\Theta_\delta\}_{\delta\geq0}\), or its compatibility frontier. For a candidate define

\[
D_k(w,\gamma,\eta)=\max_i
|\mathcal O_k[q_i]-y_{ik}|.
\]

Its compatible allowances are precisely the upper orthant \(\delta_k\geq D_k\). Optimizing a common maximum residual in µm can summarize this frontier, but does not turn its optimum into measurement noise. The four entries share physical units; they need not share the same measurement or model-discrepancy distribution.

If validated forward intervals \([\ell_{ik},u_{ik}]\) enclose the *intended physical observables*, use

\[
D_k^- =\max_i\max(\ell_{ik}-y_{ik},y_{ik}-u_{ik},0),\qquad
D_k^+=\max_i\max(|\ell_{ik}-y_{ik}|,|u_{ik}-y_{ik}|).
\]

An allowance below a lower requirement excludes a candidate; an allowance above every upper requirement certifies the candidate's observation inequalities. Interval overlap alone does not construct one common mechanically realizable state. Current stored-polynomial observer intervals do not yet enclose the exact equilibrium or continuum response, so this rule presently certifies a **stored-curve proxy**, not a physical member of \(\Theta_\delta\).

The completed direct candidates require the following componentwise allowances for that proxy. Numbers are rounded; exact fractions are in the receipt.

| Convention | Maximum | Regional mean | Apical | Basal |
|---|---:|---:|---:|---:|
| Disk / spatial-max components | 4.453269 | 5.379836 | 5.114619 | 3.136498 |
| Disk / regional-mean components | 5.141405 | 6.180135 | 5.894309 | 3.322522 |
| Line / spatial-max components | 4.369371 | 5.696665 | 5.509162 | 3.198092 |
| Line / regional-mean components | 4.501763 | 6.245135 | 5.954704 | 3.314836 |

All entries are µm, and every row uses all 68 observations. These are sufficient residual requirements for particular stored fields, not four physical uncertainty estimates or globally optimized physical fits.

## 4. Does an exact observation-compatible set exist?

There is no single answer until the temporal contract is stated.

**One common, fixed-geometry quasistatic branch with a scaled input cap:** if each peak is taken over the whole common path \(0\leq A\leq gI\), the nesting of those paths makes peak response nondecreasing in dose. Source records F and I have doses 0.46402059 and 0.47133397, but maximum shortenings 7.64808362 and 1.38419619 µm. Consequently the exact zero-discrepancy set is empty for this contract, for every choice of material ratios. The sharp common uniform maximum-shortening floor is 3.131943715 µm. This exclusion does not apply automatically when geometry, histories, windows or branches vary between specimens.

**Mean of local temporal peaks with zero initial shortening:** this reading is favored by Methods Sec16, although the exact implementation is not published. Record O has regional mean exactly zero and maximum 2.41505646 µm, in directly crosslinked height-summary cells. With continuous spatial fields, positive weights and the zero initial field included, every local temporal peak is nonnegative. A positive peak at one point persists over a neighborhood, whose positive weighted measure makes their mean positive. The literal zero-error observation set is therefore empty under these assumptions, even before imposing mechanics. The [temporal-order classification](TEMPORAL-ORDER-ROBUST-WITNESS.md) retains this record and separates exact witnesses from nonattained limits.

This is not a criticism of the experiment. The zero is a reported source value, not a validated infinitely precise assertion. Rounding, censoring, reference construction or annotation can change the observation contract. Without a measured spatial regularity/resolution bound, a positive peak can become arbitrarily narrow: this record alone gives no positive uniform error floor. Inventing a measured epsilon, deleting O or silently changing the average/max ordering would conceal the issue.

**Sixteen positive-mean records:** exact kinematic temporal constructions exist under explicit conventions. They demonstrate summary feasibility, not one shared constitutive law, drive history or equilibrium family. The remaining heterogeneous-geometry physical zero-error set is **unresolved**; neither the public summaries nor the latest direct calculation proves it nonempty.

**Sufficiently large declared discrepancy in the flat static-reference subfamily:** if the constitutive domain admits \(g=T_a=T_b=0\), a flat identity deformation is an exact global minimizer: every bulk and area-elastic energy contribution is nonnegative and vanishes there. Its four static-reference outputs are zero; the fixed-reference constant path also has zero temporal peaks. This does not supply the missing cellularization/reference dynamics. It is compatible when component allowances exceed the source maxima, respectively 12.73324306, 8.38319196, 13 and 4.36972161 µm. This proves a trivial conditional nonemptiness statement. It does not justify error bars that large or establish biologically meaningful activation mechanics.

## 5. Operations needed for executable physical inference

1. **Freeze a cohort-wide temporal operator.** Specify the tracked coordinate chart, fixed/current spatial weights, local-peak averaging, separate component aggregation, reference growth, baseline registration and observation window. Carry unresolved global alternatives as alternatives, rather than fitting a different convention per row. Retain canonical pairing or one paired-component permutation with its stated quantifier.
2. **Declare a geometry and input contract.** Retain the protocol height range as such. Obtain per-record initial surfaces/registration where possible; otherwise label the restricted geometry family explicitly. Fluorescence-to-force proportionality and its time course remain hypotheses. Input-error bounds need their own evidence.
3. **Declare discrepancy as a family or obtain measured bounds.** Preserve all 17 samples, including the reported zero. If a quantization/censoring operator is introduced, specify its rule and threshold from source evidence; do not choose it merely to admit a favored model.
4. **Specify which mechanical solution counts.** A Newton warm-start path does not prove energetic stability, uniqueness or a physical loading history. Define the chosen stable branch, global-minimum set, or dynamical selection rule, and state whether compatibility requires existence on one branch or agreement across all admitted branches.
5. **Enclose the required forward observables.** Bound equilibrium error using an applicable existence/stability argument, then integration, spatial approximation and temporal supremum error. Force residuals, mesh differences and dense samples are diagnostics until an argument converts them into the required output enclosures. Local-peak-then-average requires tracking local peaks, not merely maximizing a list of regional means.
6. **Cover the declared compatible domain.** Use validated subdivisions/enclosures over the ratio and geometry domains, with explicit unresolved cells and treatment of tails. Candidate fitting is useful for finding an inner witness and proposing tight bounds; it is not the whole set. Freeze larger-region envelopes only after carrying the remaining compatible combinations and numerical errors.

The primary record states that no custom code was generated. A useful source request would seek exact workflow/settings, derived tracked time courses, references, masks and specimen crosswalks; it should not assert that the authors promised missing custom code. Primary microscopy is author-request, but public numerical sources and illustrative movies already support the bounded analyses documented here. No external contact was made.

These are concrete dependencies of the original physical-inference requirement. The missing acquisition details are external limits; actual compatible-set computation and forward certification remain unperformed operations. Neither category excuses replacing the requirement with a single scenario or a sensitivity grid.
