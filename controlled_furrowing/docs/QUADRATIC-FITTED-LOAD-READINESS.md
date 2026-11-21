# Readiness of the actual activated quadratic model for calibration

**Decision:** the robust solver now produces verified discrete equilibria and embedded bodies on the actual small-region geometry. It is ready for explicitly provisional discrete calibration diagnostics. It is **not yet ready for a continuum-accurate mechanical transfer claim**: activated-state mesh, quadrature and near-incompressibility accuracy remain unresolved. The manufactured benchmark's accuracy cannot be transferred to these much more strongly activated fields.

This review reads six actual final states from `runs/quadratic_robust/lambda_{500.0,5000.0}_target_{1.0,5.0,10.0}/final.json`. It reads no biological response data, performs no fit, and uses none of the failed states in `runs/quadratic_pilot` as a calibration result. Receipt: `evidence/quadratic_robust_geometry.json`, including a reproducible replay program, input hashes, all exact triangle coefficients and exact point evaluations.

## 1. Model and source provenance

Every reviewed final state uses 12×6 meridional cells, P2 displacement, the original **pointwise** logarithmic-volume penalty, and the actual small-patch geometry: height scale 33.5 µm, measured diameter 18.71024074 µm, outer radius three heights, and activation 5–25 µm below the apical surface. The declared surface parameters are \(K_b/(\mu H)=2\), \(K_a/(\mu H)=0.2\), \(T_b/(\mu H)=2\), \(T_a=0\). Activation strengths \(f/\mu=1,5,10\) are numerical load cases, not fitted biological forces.

The numerical batch used retained controller source

```text
2e00d7b2b61ecfa9006716ba9a0245af47f2c85ae7b647caf47b66b13532dba9
```

at `runs/quadratic_robust/pre_guard_hardening/robust_quadratic.py`. Final input guards and terminal-ledger persistence use controller

```text
be6f17e620ddf04533d456fff9d48cd48517f4a0709d921ab510616926b4176b
```

The constitutive/energy source remains

```text
5bb64d6b441f6db8ada9da0f1bcb44716c86ed5d816e44f75bd4e9e021fdf410
```

for `scripts/axisymmetric_quadratic.py`. Current parent-source hashes match the batch receipt. The independent reviewer verifies both controller hashes and replays energies and gradients against the frozen full-P2 model; the controller changes numerical policy, not energy or gradient.

## 2. Algebraic acceptance and failed-step accounting

All six requested targets were reached. Replay of every final field passes the predeclared scaled free-force threshold \(10^{-7}\); the largest final residual is approximately \(3.25\times10^{-9}\). Across the six continuation records, **12 accepted increments** pass replayed stationarity. **Six failed attempts** remain retained. At every attempt, the initial-vector hash equals the most recent accepted-vector hash. A rejected optimizer vector is not used as the next state. Recomputed accepted-path energy ledgers close within numerical rounding.

These checks support a verified discrete stationary-state path. They do not prove a global energy minimum, a unique branch, a biological time trajectory, or physical photon/ATP expenditure. The accepted-path work ledger refers to the specified externally ramped mechanical potential. Failed numerical trials are not physical events.

## 3. Exact whole-body discrete geometry

The independent rational verifier checks every one of the **864 quadratic triangles** in the six final bodies. The degree-four correlated Bernstein method bounds

\[
 J=\det F_m\,R/r
\]

by first forming the polynomial numerator \(\det F_m R\), then dividing its Bernstein coefficients by the degree-elevated affine reference-radius coefficients. Exact zero coefficients on the axis are handled explicitly. The complete derivation is in `QUADRATIC-INDEPENDENT-REVIEW.md`; the implementation is `scripts/mixed_geometry_certificate_v2.py`.

All six bodies also pass two independent sufficient global conditions: the existing exact disk/Jordan-boundary certificate, including quadratic radial monotonicity and top–bottom separation, and the exact vertexwise positive-definiteness test for the symmetric meridional deformation gradient, implying strong monotonicity on the convex reference rectangle. Thus each saved meridional body is globally embedded, and its axis-preserving revolution is an injective three-dimensional body. This excludes self-overlap of the represented body; contact with an external structure absent from the model is outside the claim.

| \(\lambda/\mu\) | \(f/\mu\) | Original quadrature \(J\) range | Certified whole-body discrete \(J\) interval | Whole-body \(|J-1|<1\%\)? |
|---:|---:|---:|---:|---|
| 500 | 1 | [0.997548, 1.003091] | [0.992104, 1.009470] | Yes |
| 500 | 5 | [0.979917, 1.017614] | [0.957415, 1.038645] | No; violation proved |
| 500 | 10 | [0.966828, 1.030380] | [0.929214, 1.066934] | No; violation proved |
| 5000 | 1 | [0.999230, 1.000903] | [0.997586, 1.002742] | Yes |
| 5000 | 5 | [0.990723, 1.004564] | [0.980397, 1.010771] | No; violation proved |
| 5000 | 10 | [0.980573, 1.008819] | [0.959310, 1.021443] | No; violation proved |

Displayed numbers approximate exact rational bounds in the receipt. A bound endpoint need not be attained. The original quadrature extrema are observed values at a finite point set, not global extrema or error bounds.

To distinguish conservative bounds from actual violations, we also evaluate \(J\) exactly at triangle vertices, edge midpoints and centroids with positive reference radius. These are attained values of the stored discrete map. At \(\lambda=5000,f=5\), one such point has \(J=1.0107707674\), exceeding 1.01, although the original quadrature maximum was only 1.004564. At \(\lambda=5000,f=10\), exact evaluated values include 0.9832949015 and 1.0207526236. The stronger-load cases therefore fail a 1% diagnostic by explicit rational point inequalities, not just by lack of a sufficiently sharp bound. Those points still do not necessarily realize the global extrema.

One percent is a transparent numerical diagnostic used here to report local distortion; it is not an experimentally justified tissue compressibility limit. The actual biologically acceptable tolerance must be established independently.

## 4. What the penalty comparison says

The maximum activated material-column shortening changes as follows:

| \(f/\mu\) | \(\lambda=500\), shortening/\(H\) | \(\lambda=5000\), shortening/\(H\) |
|---:|---:|---:|
| 1 | 0.06591696 | 0.05955996 |
| 5 | 0.21421944 | 0.20909121 |
| 10 | 0.31734926 | 0.31294742 |

The endpoint changes by about 9.6%, 2.4% and 1.4% respectively. A relatively small change in one endpoint does not certify the complete geometry. The much stronger volume oscillations in the earlier projected-pressure model are avoided by restoring the full pointwise penalty, but substantial local distortions remain at these activated loads.

Changing \(\lambda\) changes the finite-compressibility continuum material law. If that finite law is the scientific model, \(\lambda\) needs independent identification or declared uncertainty. If instead \(\lambda\) is a numerical penalty approximating an incompressibility constraint, convergence of both geometry and local volume distortion as \(\lambda\to\infty\) must be assessed with sufficient displacement resolution. Raising \(\lambda\) alone is not evidence that locking or discretization error has disappeared.

## 5. Next accuracy gate, before a physical calibration freeze

1. Repeat these **actual activated** cases on a finer matched mesh, initially 16×8, retaining the same activation geometry, constitutive terms, boundary, loads and solver thresholds. Compare against 12×6; refine further if differences remain material. Do not re-fit activation gain separately to conceal resolution changes.
2. Compare complete apical and basal curves on a common physical radial domain, signed displacement extrema, material-column shortening, local \(J\) certificates, energy and force residuals. Distinguish material-column height from Eulerian tissue thickness.
3. Raise volume and surface quadrature independently. Converged algebraic equations at one quadrature rule do not establish accurate nonlinear energy integration.
4. Separate finite material compressibility from the intended incompressible limit. Declare a numerical local-volume tolerance, report its scientific justification separately, and track local extrema rather than global volume alone. Retain all tested penalty branches.
5. Establish a calibration-only numerical uncertainty budget and verify sensitivity of fitted gain and secondary small-region displacements to mesh/penalty choices. Only then freeze a physical small-to-large transfer forecast. Static peak-proxy and missing slice/time issues remain separate biological limitations.

A provisional fit to a precisely declared discrete model can still help diagnose whether that model is even compatible with the small-region measurements. It must be labeled as such; it cannot close the original continuum-mechanics or biological-prediction gates. Adding further material complexity before these accuracy checks would make failure attribution harder.

The retained independent review passes **49 checks**, with all input/source hashes and exact coefficients available for replay. This is progress in solving and verifying the actual discrete mechanics. It establishes no new biological concept or confirmed original general theorem.
