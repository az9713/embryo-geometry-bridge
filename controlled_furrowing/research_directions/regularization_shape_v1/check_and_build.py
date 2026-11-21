"""Build/check a prospective theory plan; do not change scientific targets."""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import csv
import html
import importlib.util
import json
import re
import subprocess

HERE = Path(__file__).resolve().parent
MECH = HERE.parents[1]
PROJECT = MECH.parent


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def main():
    # Verify prior-input identity and read scalar counts only. No raw integrals
    # or old AI process observations are reproduced by this documentary check.
    source = MECH / 'runs/nonaxisymmetric_surface_sign_preparation_v1/independent_actual_readback.json'
    old = json.loads(source.read_text(encoding='utf8'))
    linked = old['verified_source_and_actual_sha256']
    failures = [name for name, expected in linked.items()
                if not (MECH/name).is_file() or digest(MECH/name) != expected]
    assert not failures, failures
    counts = [{'case': c['case_id'], 'windows': len(c['all_windows']),
               'counts': c['independent_counts']} for c in old['case_readback']]
    assert old['checked_windows'] == 672
    for row in counts:
        assert row['windows'] == 336
        assert row['counts'] == {'NONNEGATIVE_ON_WINDOW': 276,
                                'UNRESOLVED': 2, 'STRICT_NEGATIVE_PATCH': 58}
    assert old['physical_wavelength_or_continuum_stability_claimed'] is False

    # Exact examples of a rank-one laminate. They are sanity checks for the
    # written general argument, not a computer proof of relaxation or a model fit.
    examples = []
    for k, tension in [(F(1), F(0)), (F(5), F(1)), (F(7), F(2))]:
        target = 1 - tension/k
        for cos, sin in [(F(3, 5), F(4, 5)), (F(5, 13), F(12, 13)),
                         (F(8, 17), F(15, 17))]:
            j = target*cos
            t = target*sin
            assert j*j + t*t == target*target
            psi = lambda s: k/2*(s-1)**2 + tension*s
            assert psi(j) > psi(target)
            assert k*(target-1)+tension == 0
            # A+=(j e1+t e3 | e2), A-=(j e1-t e3 | e2):
            # rank-one difference, exact mean A, and squared area target^2.
            examples.append({'K': str(k), 'T': str(tension), 'J': str(j),
                             'J_star': str(target), 'normal_slope': str(t),
                             'energy_drop': str(psi(j)-psi(target))})

    # Plot assumptions explicitly; no scientific data are used in this figure.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    xs = [i/400 for i in range(1, 601)]
    k, tension, target = 5., 1., .8
    psi = lambda j: k/2*(j-1)**2+tension*j
    ys = [psi(j) for j in xs]
    relaxed = [psi(max(j, target)) for j in xs]
    stress = [k*(j-1)+tension for j in xs]
    relaxed_stress = [max(t, 0) for t in stress]
    with (HERE/'illustration.csv').open('w', newline='', encoding='utf8') as f:
        writer = csv.writer(f)
        writer.writerow(['J', 'original_density', 'relaxed_density',
                         'original_stress', 'relaxed_stress'])
        writer.writerows(zip(xs, ys, relaxed, stress, relaxed_stress))
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                         'figure.facecolor': '#102033', 'axes.facecolor': '#102033',
                         'text.color': '#edf4ff', 'axes.labelcolor': '#edf4ff',
                         'xtick.color': '#ccd9ec', 'ytick.color': '#ccd9ec',
                         'axes.edgecolor': '#6a7f9b', 'svg.fonttype': 'none'})
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), constrained_layout=True)
    axes[0].plot(xs, ys, color='#f3b766', label='Original surface density', lw=2.5)
    axes[0].plot(xs, relaxed, color='#6bd3d0', ls='--', label='Isolated-cell relaxed density', lw=2.5)
    axes[1].plot(xs, stress, color='#f3b766', label='Original derivative', lw=2.5)
    axes[1].plot(xs, relaxed_stress, color='#6bd3d0', ls='--', label='Relaxed derivative', lw=2.5)
    for ax in axes:
        ax.axvline(target, color='#b0bdd1', ls=':', lw=1)
        ax.set_xlabel('Area ratio J')
        ax.grid(alpha=.13)
        legend = ax.legend(facecolor='#1b3049', edgecolor='#6a7f9b', fontsize=9)
        for label in legend.get_texts():
            label.set_color('#edf4ff')
    axes[0].set_ylabel('Energy per reference area (arbitrary units)')
    axes[1].set_ylabel('Surface energy derivative (arbitrary units)')
    axes[1].axhline(0, color='#b0bdd1', lw=.6)
    fig.suptitle('Analytic illustration: K = 5, T = 1, J* = 0.8; no experimental fit', fontsize=12)
    fig.savefig(HERE/'surface-relaxation.svg', metadata={'Creator': 'Research Project', 'Date': None})
    plt.close(fig)

    renderer = MECH/'scripts/build_report.py'
    spec = importlib.util.spec_from_file_location('_plan_markdown', renderer)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    md = (HERE/'PLAN.md').read_text(encoding='utf8')
    body = module.convert(md)
    figure = ('<figure><img src="surface-relaxation.svg" alt="Analytic area-energy illustration: '
              'the relaxed density is constant on the compressed range, and its derivative is zero there.">'
              '<figcaption>Illustration of the isolated-surface lemma only. This is a different '
              'variational density; it does not validate clipping forces in an old bulk equilibrium. '
              'Formula samples: <a href="illustration.csv">CSV</a>.</figcaption></figure>')
    body = body.replace('<h2>5. The proposed continuum theorem</h2>', figure+'<h2>5. The proposed continuum theorem</h2>')
    css = '''*{box-sizing:border-box}html{color-scheme:dark}body{margin:0;background:#0c1727;color:#eaf1fb;font:17px/1.72 system-ui,sans-serif;overflow-wrap:anywhere}main{max-width:1100px;margin:auto;padding:35px clamp(16px,4vw,48px) 70px}a{color:#8ce2e3}h1{font-size:clamp(2rem,4.7vw,3.2rem);line-height:1.13}h2{color:#aee4f2;margin-top:2.2em;line-height:1.25}table{display:block;max-width:100%;overflow:auto;border-collapse:collapse;font-size:.91em}td,th{padding:12px;border:1px solid #455a73;vertical-align:top}.math.display{max-width:100%;overflow:auto;background:#162a41;padding:14px;border-radius:8px}img{display:block;max-width:100%;height:auto}figure{margin:30px 0}figcaption{font-size:.88em;color:#c2cfe1}.notice{border-left:4px solid #e4ba67;padding:16px;background:#1b2a3d}nav{display:flex;flex-wrap:wrap;gap:12px 24px;border-bottom:1px solid #455a73;padding-bottom:18px}li{margin:.5em 0}code{overflow-wrap:anywhere}@media(max-width:600px){body{font-size:16px}main{padding-top:22px}.katex{font-size:1.04em}}'''
    page = ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width,initial-scale=1">'
            '<title>Next theorem target: coarse shape beyond unresolved wrinkles</title>'
            '<link rel="stylesheet" href="../../../docs/assets/katex/katex.min.css">'
            '<style>'+css+'</style></head><body><main><nav>'
            '<a href="PLAN.md">Plan source</a><a href="EVIDENCE.json">Verification scope</a>'
            '<a href="LITERATURE-AUDIT.json">Originality audit</a>'
            '<a href="../../docs/index.html">Existing dossier</a></nav>'
            '<p class="notice">Research proposal. The main continuum and shape-stability '
            'theorems are unproved. No confirmed novelty or biological forecast is claimed.</p>'
            +body+'</main></body></html>')
    pattern = re.compile(r'(<(span|div) class="math(?: display)?" data-tex="([^"]*)">)\s*</\2>', re.S)
    matches = list(pattern.finditer(page))
    jobs = [{'tex': html.unescape(m[3]), 'display': m[2] == 'div'} for m in matches]
    js = r'''const fs=require('node:fs');const k=require(process.argv[1]);const xs=JSON.parse(fs.readFileSync(0,'utf8'));process.stdout.write(JSON.stringify(xs.map(x=>k.renderToString(x.tex,{displayMode:x.display,throwOnError:true,strict:'error',output:'htmlAndMathml'}))));'''
    result = subprocess.run(['node', '-e', js, str(PROJECT/'docs/assets/katex/katex.min.js')],
                            input=json.dumps(jobs), encoding='utf8', capture_output=True, check=True)
    rendered = iter(json.loads(result.stdout))
    page = pattern.sub(lambda m: m[1]+next(rendered)+'</'+m[2]+'>', page)
    expected_math = len(re.findall(r'\\\((.*?)\\\)|\\\[(.*?)\\\]', md, re.S))
    assert expected_math == len(jobs)
    (HERE/'PLAN.html').write_text(page, encoding='utf8', newline='\n')
    evidence = {'status': 'PASS_PLAN_DOCUMENT_AND_SCOPED_CHECKS',
                'scope': 'Strict math compilation, existing input identity and exact laminate examples; not a proof of the proposed continuum theorem.',
                'source_receipt': str(source.relative_to(PROJECT)).replace('\\', '/'),
                'source_receipt_sha256': digest(source), 'prior_linked_hashes_checked': len(linked),
                'stress_counts_readback': counts, 'exact_laminate_examples': examples,
                'strictly_compiled_math_expressions': len(jobs),
                'written_proof_status': {'isolated_surface_cell': 'DERIVATION_STANDARD_NOT_NOVEL',
                                         'conditional_observable_error': 'ELEMENTARY_CONDITIONAL_LEMMA',
                                         'bulk_boundary_relaxation': 'UNPROVED_RESEARCH_TARGET',
                                         'useful_observable_growth_constant': 'UNPROVED_RESEARCH_TARGET'},
                'scientific_sources_modified': False,
                'file_sha256': {name: digest(HERE/name) for name in
                                ['PLAN.md', 'PLAN.html', 'surface-relaxation.svg', 'illustration.csv',
                                 'LITERATURE-AUDIT.json', 'check_and_build.py']}}
    (HERE/'EVIDENCE.json').write_text(json.dumps(evidence, indent=2)+'\n', encoding='utf8')
    print(json.dumps({'status': evidence['status'], 'prior_linked_hashes_checked': len(linked),
                      'strict_math_expressions': len(jobs), 'exact_laminate_examples': len(examples)}))


if __name__ == '__main__':
    main()
