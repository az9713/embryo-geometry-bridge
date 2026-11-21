"""Build a local, receipt-linked report without rewriting the frozen dossier."""
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from build_report import convert

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'runs/transfer_ambiguity_linear_v1'
RESULT_SHA = '68ddc210523cb40001c99babb251916ac22c62eecd836d976fac1ab78dc0fe3f'
REVIEW_SHA = '8fc50b42694b3db513e6d2166c6d539fbe67c66c2b3c6895c90df12f3856d072'


def sha(p):
    return sha256(Path(p).read_bytes()).hexdigest()


def main():
    result_path = OUT/'actual/results.json'
    review_path = OUT/'actual_independent/receipt.json'
    assert sha(result_path) == RESULT_SHA and sha(review_path) == REVIEW_SHA
    result = json.loads(result_path.read_text())
    review = json.loads(review_path.read_text())
    assert review['status'].startswith('PASS')
    rows = result['rows']
    a = [x for x in rows if x['mode'] == 'source_residual_allowance']
    b = [x for x in rows if x['mode'] == 'small_prediction_change']
    for row in rows:
        assert Q(row['contrast_range_width_um_exact'][0]) <= Q(row['contrast_range_width_um_exact'][1])
    plt.style.use('dark_background')
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    fig.patch.set_facecolor('#0b111b')
    for ax in axes:
        ax.set_facecolor('#111e30')
        ax.grid(alpha=.18)
        ax.set_ylabel('Large-region contrast range (µm)')
    axes[0].plot([float(Q(x['allowance_um_exact'])) for x in a],
                 [x['contrast_range_width_um_display'][1] for x in a], 'o-', color='#7ed5f4')
    axes[0].set_xlabel('Assumed common calibration residual allowance (µm)')
    axes[0].set_title('Rounded affine calibration constraint')
    axes[1].loglog([float(Q(x['allowance_um_exact'])) for x in b],
                   [x['contrast_range_width_um_display'][1] for x in b], 'o-', color='#ffd080')
    axes[1].set_xlabel('Maximum change from small-region baseline (µm)')
    axes[1].set_title('Rounded affine near-invisibility constraint')
    fig.suptitle('LOCAL LINEAR MODEL ONLY — nonlinear and experimental compatibility unproved', fontsize=11)
    image = OUT/'local_affine_frontier.png'
    assert not image.exists()
    fig.savefig(image, dpi=180)
    plt.close(fig)
    table = '\n'.join('| '+str(float(Q(x['allowance_um_exact'])))+' | '
                      +format(x['contrast_range_width_um_display'][0], '.9f')+' | '
                      +format(x['contrast_range_width_um_display'][1], '.9f')+' |' for x in b)
    text = r'''# Can small-region measurements hide larger-region shape differences?

**A local diagnostic found a direction worth testing, not a biological discovery.**
The independently checked rounded affine model permits a large-region shape-contrast range of approximately **0.06042515 µm** while every one of its 68 small-region predictions stays within **0.001 µm of the baseline prediction**. Two endpoints can therefore differ by as much as 0.002 µm in a small-region channel; the 0.001 µm constraint applies separately relative to the baseline.

Under the different constraint that each affine small-region prediction stays within the baseline's worst source-data residual allowance, the large contrast range is approximately **0.30181702 µm**. That allowance is approximately **5.37983567 µm**. It is a calculated residual of an imperfect conditional model, **not a measured noise level**.

![Two conditional local affine frontiers; neither curve is a nonlinear identified set.](local_affine_frontier.png)

## Precisely what was calculated

One existing large-region field was used, at the same baseline material and gain as the archived 17-small-record sensitivity calculation. Only the activation radius differs. The selected large diameter is 49.8939753 µm and the optical dose is 0.89767374. Larger experimental response tables were not read. All values refer to one static current-disk/spatial-maximum observation convention.

The five perturbation coordinates are log apical area modulus, log basal area modulus, additive apical surface tension, log basal surface tension and log activation gain. The first, second, fourth and fifth coordinates lie in [−0.1,0.1]; additive apical tension lies in [0,0.02]. These are local search bounds, not an exhaustive physical parameter domain.

The signed shape observable uses the current physical radius:

\[
Q=z_a(0)-z_a(a),\qquad a=24.94698765\,\mu\mathrm m.
\]

On the edge containing the shoulder, write its quadratic radius and height as \(R(t),Z(t)\), with \(R(t)=a\). Differentiation gives

\[
\delta Q=\delta z_a(0)-\sum_jN_j(t)\,\delta Z_j
+\frac{Z'(t)}{R'(t)}\sum_jN_j(t)\,\delta R_j.
\]

The last term accounts for the node motion needed to evaluate height at a **fixed current radius**. The shoulder is strictly inside one edge, its radial derivative is positive, and the symmetry axis remains fixed radially. Exact rational inverse brackets enclose this stored-graph derivative. They do not enclose the error in the mechanical equilibrium or its parameter derivative.

For the discrete force equation \(F(u,\theta)=0\), implicit differentiation at a regular exact root would give

\[
H S+B=0,\qquad H=\partial_uF,\quad B=\partial_\theta F,
\qquad c=(\partial_uQ)S.
\]

Here the matrices and solution are floating diagnostics at the saved state. The five entries of \(c\), in the coordinate order above, are approximately

\[
c=(-0.04576027,\;0.39682220,\;2.33809289,\;-0.55166505,\;-1.77749089).
\]

The archived small-region map is \(r+Jd\), where \(r\) has 68 entries. We solve both signs of the linear objective \(c\cdot d\), subject either to \(\lVert r+Jd\rVert_\infty\le\delta\), or to \(\lVert Jd\rVert_\infty\le\varepsilon\), and the same trust box. The latter compares predictions with the baseline; it does not assert that either prediction fits the experiment.

## Exact certificates, limited to the rounded affine optimization

For a maximization with constraints \(Ad\le b\), any nonnegative multiplier \(\lambda\) yields

\[
c\cdot d\le\lambda\cdot b+
\sum_j\max\{(c-A^T\lambda)_j l_j,
                 (c-A^T\lambda)_j u_j\}.
\]

This follows by adding \(\lambda\cdot(b-Ad)\ge0\) and maximizing the remaining linear form over the box. The producer retains an exactly feasible rational step and these exact rational weak-dual bounds. The independent reviewer reconstructed all 20 primal/dual certificates without calling the producer's optimizer routine. Their largest optimization gap is approximately 1.12×10⁻¹⁵ µm. This is numerical accuracy for the declared affine optimization, not physical prediction accuracy.

| Allowed small prediction change from baseline, µm | Affine large contrast range lower, µm | Upper, µm |
|---:|---:|---:|
'''+table+r'''

## Evidence and limits

- [Frozen source and input protocol](preparation/protocol.json), [synthetic and preparation checks](preparation/preflight.json), and [independent prelaunch review](preparation/prelaunch_independent_review.json).
- [Actual result with all exact LP witnesses and multipliers](actual/results.json), [saved mechanical matrices](actual/tangents.npz), [worker terminal](actual/terminal.json), and [root's observed OS exit](actual/os_exit_observation.json).
- [Independent actual verification](actual_independent/receipt.json): rebuilt all 68 rows, checked all 20 optimization certificates, replayed the saved linear systems and an adjoint calculation, and checked moving-radius derivatives with multiple perturbation sizes.

The small-data matrix has full numerical rank five. A weak direction exists, not an exact nullspace. Several small-observable derivatives were obtained from central differences of maxima; switching maximizing locations and one-sided physical directions must be tested in the nonlinear model. Exact mean derivatives do not cure that separate issue.

The present work does not establish two source-compatible nonlinear equilibria, a physical compatible-set frontier, continuous stable loading, continuum or nonaxisymmetric error, an experimentally detectable gap, or novelty. Published epithelial inverse-model and optogenetic mechanics studies already motivate identifiability and geometry-dependent responses; see the [specific prior-work audit](../../docs/TRANSFER-CLAIM-AUDIT.html).

## The next falsifiable test

Before any new nonlinear outcomes, freeze two parameter choices from the local frontier. Re-equilibrate every small-dose case and the same larger input without condition-specific refitting. Compare all 68 nonlinear small summaries and the large contrast against the affine predictions. Retain failed fits and numerical failures. A surviving, substantial gap would justify certification of those particular branches; a disappearing gap would refute this local extrapolation. Experimental compatibility still requires a justified source observation and uncertainty contract.

The original five-gate study remains incomplete. The B and L exact-integral root proofs concern separate existing small-case parameter tubes; they do not certify the new large-field tangent or these proposed witnesses.
'''
    md = OUT/'REPORT.md'
    assert not md.exists()
    md.write_text(text, encoding='utf-8')
    css = 'body{background:#0b111b;color:#e7eef7;font:17px/1.65 system-ui;margin:0}main{max-width:1080px;margin:auto;padding:32px 20px 70px}a{color:#8ad4ef}h1,h2{line-height:1.2}h2{margin-top:2em;color:#aedbe8}img{max-width:100%}.math.display{overflow:auto;background:#111e30;padding:16px}table{border-collapse:collapse;display:block;overflow:auto}td,th{padding:9px;border:1px solid #41516b}code{overflow-wrap:anywhere}#math-errors{color:#ff9999}'
    html = '<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Local transfer ambiguity: evidence and limits</title><link rel="icon" type="image/svg+xml" href="data:image/svg+xml,%3Csvg%20xmlns%3D%22http%3A%2F%2Fwww.w3.org%2F2000%2Fsvg%22%20viewBox%3D%220%200%2016%2016%22%3E%3Ccircle%20cx%3D%228%22%20cy%3D%228%22%20r%3D%226%22%20fill%3D%22%237ed5f4%22%2F%3E%3C%2Fsvg%3E"><link rel="stylesheet" href="../../../docs/assets/katex/katex.min.css"><style>'+css+'</style></head><body><main><a href="../../docs/index.html">Controlled-furrowing dossier</a>'+convert(text)+'<p id="math-errors"></p></main><script src="../../../docs/assets/katex/katex.min.js"></script><script>document.querySelectorAll(".math").forEach(e=>{try{katex.render(e.dataset.tex,e,{displayMode:e.classList.contains("display"),throwOnError:true})}catch(x){document.getElementById("math-errors").textContent+=String(x)}});</script></body></html>'
    page = OUT/'REPORT.html'
    assert not page.exists()
    page.write_text(html, encoding='utf-8')
    receipt = {'status': 'BUILT_FROM_INDEPENDENTLY_REVIEWED_RESULTS_BROWSER_PENDING',
               'source_sha256': sha(__file__), 'results_sha256': RESULT_SHA, 'review_sha256': REVIEW_SHA,
               'files': {str(p.relative_to(ROOT)): sha(p) for p in (md, page, image)}}
    (OUT/'report_build_receipt.json').open('x', encoding='utf-8').write(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
