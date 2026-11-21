"""Build a local scientific checkpoint from the independent rational LP review."""

from collections import Counter
import csv
from fractions import Fraction as Q
from hashlib import sha256
from html import escape
from importlib.util import module_from_spec, spec_from_file_location
import json
import math
import os
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.ticker import FormatStrFormatter, MultipleLocator


HERE = Path(__file__).resolve().parent
RUN = HERE.parent
STUDY = RUN.parents[1]
SITE = STUDY.parent
OTHER = STUDY / 'runs/transfer_ambiguity_nonlinear_probe_v1'
REVIEW = RUN / 'independent_readback'
AUTH = REVIEW / 'actual/receipt.json'
ACCEPTED = REVIEW / 'root_acceptance.json'
OS = REVIEW / 'actual/os_exit_observation.json'
TERMINAL = REVIEW / 'actual/terminal.json'
HELD = RUN / 'ASSESSMENT.md'
AT = RUN / 'results.json'
PROTOCOL = RUN / 'protocol.json'
CONVERTER = STUDY / 'scripts/build_report.py'
HELPER = OTHER / 'certification_design/exact_targets_v1/certification_checkpoint_v1/build_checkpoint.py'
ASSETS = SITE / 'docs/assets/katex'
ICON = OTHER / 'report/favicon.png'
SOURCE_CONTRACT = STUDY / 'docs/JOINT-OBSERVATION-CONTRACT.md'
PHYSICS = STUDY / 'docs/PHYSICAL-PARAMETER-DOMAIN.md'
PLAN = STUDY / 'PLAN.md'
ORIGINALITY = OTHER / 'originality_followup/README.md'
SOURCE_REFRESH = STUDY / 'runs/source_contract_public_refresh_v1/PUBLIC-REFRESH.md'
REVISION = HERE / 'revision_0/revision_manifest.json'
WIDTHS = HERE / 'widths.csv'
TITLE = 'How source linkage changes the affine transfer range'
VARIANTS = ('current_disk__spatial_max', 'current_disk__regional_mean',
            'current_radial_line__spatial_max', 'current_radial_line__regional_mean')
LABELS = ('Disk, component maxima', 'Disk, component means',
          'Radial line, component maxima', 'Radial line, component means')
GROUPS = ('ordered68', 'source_joined34')
PINS = {
    AUTH: 'd3b90c6f2b2336afd4141d5b7b4d5bf67bfe39cd64741aded8003a6bbf735ad2',
    ACCEPTED: 'dc21bc539c1cbd96ecc5c88af05bb671315d2a6fb5fc2ce06a45ec8fae5a70f4',
    OS: 'ae7a89b061b5ffc912a0b93da6c70f2166e4d0db37fb23a5c657f3851c5f1690',
    TERMINAL: '28d989de1d54983829cb1a46fe8fb26aed75c55f995bbb11ef561dabd0f79585',
    AT: '9e2d67bd0eba02b300f9a8936513c046b7e170abdb51422eaff77edcf764d517',
    PROTOCOL: 'e409bd1d039c6cefbf87c67d98c8de433bfbbdc0b9b1caa7bacbf0fcb0da69d3',
    HELD: 'fe18a94d058627890d3fca6bafd3fd0073007e64547f61ddf21d9d1df3f3892d',
    CONVERTER: '90148d53d74d8ea5cc8a665e6a3660abf2db72487a1fa59586caba98f9cc38e2',
    HELPER: '50524d01db5003b0839e873efb37d41bb0552ad0f0c162f73bd668ad0e58d31b',
    ASSETS / 'katex.min.css': '5a0ff74f11887f7260dc16b216906f51352b84ca29f21abb9a14b62a2ca6eb15',
    ASSETS / 'katex.min.js': 'e6bfe5deebd4c7ccd272055bab63bd3ab2c73b907b6e6a22d352740a81381fd4',
    ICON: '961e6e85ec16258652bf3c3f322d336fb598b90bccc49821047babfef3925530',
    SOURCE_CONTRACT: 'd805ff0f9c757d35fdd61283ae1591be3873d73d7daadbaac521020454654d04',
    PHYSICS: '8d6a397c6efc3e7a91d9716630355dfb91acbef51928e75437b4ad4b24fd728f',
    PLAN: '993240af5d30c8e5b72ab31a1278a524314d38ede48486b677b021d84bcde500',
    ORIGINALITY: 'be18025c853273ebed98db4cbb47dac48006236631bee671090d22f279410c43',
    SOURCE_REFRESH: '35e37316e8303d4007464a518859a9e702b2c922aad19103efb4e59fa5c535bc',
    REVISION: 'ad8739400a6571c52f4cf22426d4f13666318925ba9f6dd17c34943cda3e2c6b',
}
BG, INK, MUTED, GRID, CYAN, AMBER = '#122034', '#edf6ff', '#c5d5e4', '#3a526e', '#83ddeb', '#ffcb84'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(p):
    return sha256(p.read_bytes()).hexdigest()


def link(p):
    return Path(os.path.relpath(p, HERE)).as_posix()


def pair(values):
    lo, hi = map(Q, values)
    require(lo <= hi, 'ordered proof bounds')
    return lo, hi


def numbers():
    require(all(sha(path) == h for path, h in PINS.items()), 'frozen source or renderer bytes changed')
    audit, adopted, observed, terminal, prod, protocol = map(
        lambda p: json.loads(p.read_bytes()), (AUTH, ACCEPTED, OS, TERMINAL, AT, PROTOCOL))
    require(audit['status'] == 'PASS_INDEPENDENT_ASSOCIATION_SAFE_AFFINE_READBACK'
            and adopted['status'] == 'PASS_ROOT_ACCEPTED_INDEPENDENT_ASSOCIATION_SAFE_AFFINE_READBACK'
            and adopted['independent_receipt_sha256'] == observed['receipt_sha256'] == terminal['receipt_sha256'] == PINS[AUTH]
            and adopted['independent_true_os_sha256'] == PINS[OS] and observed['terminal_sha256'] == PINS[TERMINAL]
            and type(observed['observed_exit_code']) is int and observed['observed_exit_code'] == 0
            and type(audit['pid']) is int and audit['pid'] > 0 and observed['native_pid'] == audit['pid']
            and audit['all_16_endpoint_certificates_verified'] is True and audit['all_8_ranges_verified'] is True
            and audit['physical_or_experimental_identifiability_proved'] is False
            and audit['full_five_gate_complete'] is False and adopted['conditional_rounded_affine_scope_only'] is True
            and prod['protocol_sha256'] == PINS[PROTOCOL] and Q(protocol['allowance_um_exact']) == Q(1, 1000),
            'matching independent accepted archived-scope evidence')
    mapping = {(r['variant'], r['group']): r for r in prod['rows']}
    widths = {}
    for row in audit['rows']:
        key = row['variant'], row['group']
        lo, hi = pair(row['width_um_exact'])
        require(key not in widths and 0 < lo <= hi and hi - lo < Q(19, 10**16)
                and (lo, hi) == pair(mapping[key]['width_um_exact'])
                and type(row['rank_of_exact_stored_matrix']) is int and row['rank_of_exact_stored_matrix'] == 5
                and row['rows'] == mapping[key]['row_count'] == (68 if row['group'] == 'ordered68' else 34),
                'same independently verified width and exact stored matrix rank')
        widths[key] = lo, hi
    require(set(widths) == {(v, g) for v in VARIANTS for g in GROUPS}
            and all(widths[a, 'source_joined34'] == widths[b, 'source_joined34']
                    for a, b in zip(VARIANTS[::2], VARIANTS[1::2])), 'same two exact repeated 34-row widths')
    growth = {}
    for v in VARIANTS:
        a, b = widths[v, 'ordered68'], widths[v, 'source_joined34']
        bounds = b[0]/a[1]-1, b[1]/a[0]-1
        require(0 < bounds[0] <= bounds[1] and bounds == pair(audit['widening'][v]['ratio_increase_exact']),
                'recomputed exact relative increase')
        growth[v] = bounds
    rows = audit['row_labels_and_dose_evidence']
    counts = Counter(r['dose_association'] for r in rows)
    require(len({r['record_id'] for r in rows}) == 17 and len(rows) == 68
            and audit['verified_specimen_ID_joins'] == 0
            and counts == {'EXACT_PUBLISHED_SCALAR_JOIN': 34, 'COMMON_B_R_ORDER_ASSUMPTION_PRESERVING_COMPONENT_PAIR': 34},
            'source association audit')
    law = audit['source_law']
    require(law['L0_um_exact'] == '67/2' and law['tangent_fixed_kappa_over_mu_exact'] == '5000'
            and law['optional_common_conversion_ratio'].startswith('g*beta/mu'), 'fixed scale and gauge semantics')
    return widths, growth, counts


def mid(bounds):
    return float(sum(bounds) / 2)


def edge_float(q, upward):
    x = float(q)
    if (Q.from_float(x) < q if upward else Q.from_float(x) > q):
        x = math.nextafter(x, math.inf if upward else -math.inf)
    return x


def figures(widths):
    fig, ax = plt.subplots(figsize=(7.7, 4.8), dpi=185, facecolor=BG)
    ax.set_facecolor(BG)
    fig.subplots_adjust(left=.355, right=.96, top=.83, bottom=.18)
    glyphs = {}
    for i, (v, label) in enumerate(zip(VARIANTS, LABELS, strict=True)):
        y = 3 - i
        points = [mid(widths[v, g]) for g in GROUPS]
        ax.plot(points, [y, y], color=GRID, linewidth=2.4, zorder=1)
        for j, (group, color) in enumerate(zip(GROUPS, (AMBER, CYAN), strict=True)):
            lo, hi = widths[v, group]
            a, b = edge_float(lo, False), edge_float(hi, True)
            glyphs[v + '/' + group] = {'y': y, 'width_lower_um_exact': str(lo), 'width_upper_um_exact': str(hi),
                    'plot_endpoint_bounds_binary64': [a, b], 'plot_midpoint_binary64': points[j]}
            line, = ax.plot([a,b], [y,y], color=color, linewidth=2, solid_capstyle='butt')
            require(list(line.get_xdata()) == [a,b] and Q.from_float(a) <= lo <= hi <= Q.from_float(b),
                    'plotted outward width interval drift')
            point, = ax.plot([points[j]], [y], marker='s' if j == 0 else 'o', markersize=9,
                             markerfacecolor=BG if j == 0 else color, markeredgecolor=color,
                             markeredgewidth=2, linestyle='None', zorder=4)
            require(point.get_xdata()[0] == mid((lo, hi)), 'plotted marker center drift')
            ax.text(points[j], y+.22, f'{points[j]:.5f}', color=color, fontsize=10, ha='center')
        ax.text(-.03, y, label.replace(', ', '\n'), transform=ax.get_yaxis_transform(), fontsize=11, color=INK,
                ha='right', va='center', linespacing=1.12)
    ax.text(.0, 1.11, '□  68 ordered rows', color=AMBER, fontsize=12, transform=ax.transAxes)
    ax.text(.53, 1.11, '●  34 dose-linked rows', color=CYAN, fontsize=12, transform=ax.transAxes)
    ax.set_xlim(.053,.101)
    ax.set_ylim(-.44,3.49)
    ax.set_yticks([])
    ax.xaxis.set_major_locator(MultipleLocator(.01))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.spines[['top','right','left']].set_visible(False)
    ax.spines['bottom'].set_color(GRID)
    ax.tick_params(axis='x', colors=MUTED, labelsize=11)
    ax.set_xlabel('Width of modeled large-region contrast (µm)', color=MUTED, fontsize=12)
    ax.grid(axis='x', linewidth=.7, color=GRID, alpha=.65)
    with (HERE / 'affine_widths.png').open('xb') as f:
        fig.savefig(f, format='png', facecolor=BG)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.7, 3.5), dpi=185, facecolor=BG)
    fig.subplots_adjust(0, 0, 1, 1)
    ax.set(xlim=(0,1), ylim=(0,1))
    ax.axis('off')

    def node(x, y, w, h, name, color, content):
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.012,rounding_size=.018',
                                    edgecolor=color, facecolor='#182a42', linewidth=1.5))
        ax.text(x+.022, y+h-.11, name, color=color, fontsize=14, fontweight='bold', va='center')
        ax.text(x+.022, y+.10, content, color=INK, fontsize=12, linespacing=1.3, va='center')

    node(.045,.56,.402,.35,'Dose-linked summaries',CYAN,'Fig. 2c/d and Fig. 4h\n17 doses × 2 total channels')
    node(.55,.56,.402,.35,'Intact component pairs',AMBER,'Fig. 2e\n17 pairs: apical and basal')
    ax.plot([.25,.25],[.56,.43],color=CYAN,linewidth=2)
    ax.plot([.75,.75],[.56,.43],color=AMBER,linewidth=2, linestyle=(0,(4,3)))
    ax.text(.25,.34,'Repeated scalar values\nverify the dose joins',color=CYAN,fontsize=12,ha='center',va='center',
            linespacing=1.3)
    ax.text(.75,.34,'Common B–R order\nlinks pairs to doses',color=AMBER,fontsize=12,ha='center',va='center',
            linespacing=1.3)
    ax.text(.50,.105,'0 verified specimen-ID joins across panels',color=MUTED,fontsize=13,ha='center')
    with (HERE / 'source_associations.png').open('xb') as f:
        fig.savefig(f, format='png', facecolor=BG)
    plt.close(fig)
    return glyphs


def csv_file(widths, growth):
    header = ['observer_variant', 'channel_set', 'channels', 'width_lo_um_exact', 'width_hi_um_exact',
              'width_midpoint_um_display', 'increase_fraction_lo_exact', 'increase_fraction_hi_exact']
    with WIDTHS.open('x', encoding='utf-8', newline='') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(header)
        for v in VARIANTS:
            for g in GROUPS:
                lo, hi = widths[v,g]
                w.writerow([v,g,68 if g == 'ordered68' else 34,str(lo),str(hi),f'{mid((lo,hi)):.12f}',*map(str,growth[v])])
    with WIDTHS.open(encoding='utf-8', newline='') as f:
        rows = list(csv.DictReader(f))
    require(len(rows) == 8, 'all exact CSV rows')
    for r in rows:
        v, g = r['observer_variant'],r['channel_set']
        require((Q(r['width_lo_um_exact']), Q(r['width_hi_um_exact'])) == widths[v,g]
                and (Q(r['increase_fraction_lo_exact']),Q(r['increase_fraction_hi_exact'])) == growth[v],
                'wrote changed exact CSV')


def markdown(widths, growth):
    table = '\n'.join(f'| {label} | {mid(widths[v,"ordered68"]):.8f} | '
                      f'{mid(widths[v,"source_joined34"]):.8f} | {100*mid(growth[v]):.2f}% |'
                      for v, label in zip(VARIANTS,LABELS,strict=True))
    disk, line = VARIANTS[0],VARIANTS[2]
    return fr"""# {TITLE}

The modeled large-region contrast has a **{100*mid(growth[disk]):.2f}% wider affine range** under the primary disk-weighted observer when calibration keeps only the 34 dose-linked channels. Under the radial-line alternative with component maxima, the range widens **{100*mid(growth[line]):.2f}%**. The comparison uses independently checked rational linear-programming bounds on archived tangents.

![Four paired width markers. The two disk observers share the 34-row width 0.08471 micrometres, and the two radial-line observers share 0.09363 micrometres. Every 68-row width is smaller.](affine_widths.png)

| Static observer | 68 ordered rows (µm) | 34 dose-linked rows (µm) | Increase |
|:--|--:|--:|--:|
{table}

For each spatial weighting, the two 34-row cases are exactly identical. Changing the component aggregation affects only the excluded component channels. Markers show the midpoint of each certified range-width interval. The exact rational brackets, all shorter than \(1.9\times10^{{-15}}\) µm, are narrower than a plotted pixel. The display rounding applies only to the archived affine optimization. [Exact width and ratio CSV](widths.csv) · [Independent mathematical replay]({link(AUTH)}) · [Acceptance and process evidence]({link(ACCEPTED)})

## What the workbook joins

![Duplicate scalar entries link 17 doses to 17 maximum and 17 mean summaries. Fig. 2e pairs apical and basal within columns and the dose association assumes the same B through R order. No verified specimen identifier crosses these panels.](source_associations.png)

The 34-row group retains the maximum and total regional mean for each dose. Adding the 17 intact apical–basal pairs yields the 68-row group under common B–R column order. The published scalar joins give numerical associations only. The specimen crosswalk, registered boundary profiles and transported region, weights and aggregation order, per-specimen timing and growth references, uncertainty, and a matching large-region measurement remain open. [Observation contract]({link(SOURCE_CONTRACT)}) · [Public data readback]({link(SOURCE_REFRESH)})

## What the affine range measures

For archived small-channel derivative \(J\), large-contrast derivative \(c\), and parameter step \(d\) in the same trust box \(B\), the two endpoint optimizations are

\[
\max_{{d\in B,\;\lVert Jd\rVert_\infty\le\epsilon}}\;\pm c^T d,
\qquad \epsilon=0.001\,\mathrm{{\mu m}}.
\]

Here \(Q\) is the modeled apical center height minus its value at the fixed current radius of 24.94698765 µm. The allowance constrains change from the model baseline in each small channel. Its interpretation is numerical and conditional. The exact stored matrices all have rank five in the archived coordinates. That rank is specific to the five selected parameters, fixed geometry, selected observation conventions, and stored coefficients. A finite perturbation can move the spatial maximum or leave the local linear regime. No experimental error bound or physical identified set follows.

## Exact force scale and conditional dose scale

Multiply all seven force coefficients by \(\alpha>0\). On the positive physical branch, the model energy \(E\), first variation \(DE\), and second variation \(D^2E\) all scale by \(\alpha\). Stationary shapes and their local minimum character are preserved at fixed geometry, activation support, input and displacement constraints, with no independently fixed load. Branch predictions also require scale-invariant selection. The numerical trial penalty is excluded where inactive. With \(L_0=33.5\,\mathrm{{\mu m}}\), the force groups available to shape inference are

\[
\frac{{\kappa}}{{\mu}},\quad \frac{{K_s}}{{\mu L_0}},\quad
\frac{{T_s}}{{\mu L_0}},\quad \frac{{g}}{{\mu}}
\qquad(s=a,b).
\]

The archived tangent fixes \(\kappa/\mu=5000\). The modeled activation enters as \(gI\). If a common unknown conversion \(\beta\) links recorded dose to modeled dose, there is an additional exact invariance

\[
\begin{{aligned}}
I &= \beta I_{{\mathrm{{rec}}}},\\
(g,\beta) &\mapsto(g/a,a\beta),\qquad a>0.
\end{{aligned}}
\]

Under that assumption shape can constrain the effective ratio \(g\beta/\mu\). The archived calculation takes the dose scale as given. [Constitutive law and scope]({link(PHYSICS)})

## Strongest next test

After the fixed-field exact-spatial certificates, bound the 34-channel transfer over one predeclared physical parameter cell. Carry a verified nonlinear remainder and allow the maximizing radius to switch. That would test whether the affine ambiguity persists for the modeled branch. Matching the actual experimental transfer still requires a registered larger-region observation with the same operator and an error budget. Compare simpler models under that same operator before assessing detectability or a scientific contribution. Prior work establishes the general same-law ambiguity construction. [Five-gate plan]({link(PLAN)}) · [Originality follow-up]({link(ORIGINALITY)}) · [Full analytical assessment]({link(HELD)})
"""


CSS = """
*{box-sizing:border-box}
:root{color-scheme:dark;--bg:#0b1320;--panel:#122034;--fg:#edf6ff;--dim:#c5d5e4;--line:#3a526e;--a:#92e4f0;--gold:#ffcb84}
body{margin:0;background:var(--bg);color:var(--fg);font:17px/1.65 system-ui,sans-serif}
main{max-width:900px;min-width:0;margin:auto;padding:24px clamp(16px,5vw,44px) 64px}
nav{display:flex;gap:7px 20px;flex-wrap:wrap;border-bottom:1px solid var(--line);padding:0 0 12px;font-size:.84rem}
a{color:var(--a);text-underline-offset:.2em}
h1,h2{line-height:1.12;text-wrap:balance}
h1{max-width:21ch;font-size:clamp(2.05rem,6vw,3.4rem);letter-spacing:-.035em;margin:28px 0 16px}
h2{font-size:clamp(1.39rem,4vw,1.9rem);color:#d1eff7;letter-spacing:-.02em;margin:1.75em 0 .5em}
h1+p{background:var(--panel);border-left:3px solid var(--gold);padding:13px 16px;border-radius:0 6px 6px 0}
p{margin:.8em 0;overflow-wrap:break-word}
img{display:block;background:var(--panel);width:100%;height:auto;border:1px solid var(--line);border-radius:8px;margin:20px 0}
table{position:relative;display:block;max-width:100%;overflow-x:auto;border-collapse:collapse;margin:1em 0;font-size:.90em;font-variant-numeric:tabular-nums}
th,td{border:1px solid var(--line);padding:8px 9px;text-align:left;vertical-align:top}
td:nth-child(n+2){white-space:nowrap}
thead{background:var(--panel)}
.math.display{max-width:100%;position:relative;overflow-x:auto;background:var(--panel);border-radius:7px;padding:5px 9px}
@media(max-width:500px){body{font-size:16px}main{padding-top:16px}th,td{padding:7px 6px;font-size:.88em}}
"""


def render(md):
    spec = spec_from_file_location('_archived_static_katex', HELPER)
    helper = module_from_spec(spec)
    spec.loader.exec_module(helper)
    sys.path.insert(0, str(CONVERTER.parent))
    import build_report
    require(Path(build_report.__file__).resolve() == CONVERTER, 'wrong Markdown converter')
    body, count = helper.compile_math(build_report.convert(md))
    html = ('<!doctype html>\n<html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width,initial-scale=1">'
            f'<title>{escape(TITLE)}</title><link rel="icon" type="image/png" href="{link(ICON)}">'
            f'<link rel="stylesheet" href="{link(ASSETS/"katex.min.css")}"><style>{CSS}</style></head>'
            '<body><main><nav aria-label="Report resources"><a href="ASSESSMENT.md">Markdown</a>'
            '<a href="widths.csv">Exact width data</a>'
            f'<a href="{link(AUTH)}">Independent replay</a></nav>{body}</main></body></html>\n')
    parser = helper.LocalLinks()
    parser.feed(html)
    paths = {}
    for v in parser.paths:
        target = (HERE / v).resolve()
        require(target.is_relative_to(SITE) and (target.is_file() or target == HERE / 'ASSESSMENT.md'),
                f'broken local evidence link: {v}')
        paths[v] = target.relative_to(SITE).as_posix()
    return html, count, paths


def write(path, value):
    with path.open('x', encoding='utf-8', newline='\n') as f:
        f.write(value)


def main():
    widths, growth, source_counts = numbers()
    csv_file(widths, growth)
    plotted = figures(widths)
    md = markdown(widths, growth)
    html, count, paths = render(md)
    write(HERE / 'ASSESSMENT.md', md)
    write(HERE / 'ASSESSMENT.html', html)
    files = ('ASSESSMENT.md','ASSESSMENT.html','widths.csv','affine_widths.png','source_associations.png')
    prior = json.loads(REVISION.read_bytes())['sha256']
    require(all(sha(HERE/name) == prior[name] for name in files[2:]), 'responsive revision changed data or figures')
    receipt = {'status': 'PASS_BUILD_ASSOCIATION_SAFE_AFFINE_CHECKPOINT', 'source_sha256': sha(Path(__file__)),
               'input_sha256': {p.relative_to(SITE).as_posix(): h for p,h in PINS.items()},
               'file_sha256': {name: sha(HERE / name) for name in files}, 'rendered_math_count': count,
               'local_links': paths, 'plotted_width_data_verified': plotted,
               'exact_repeated_34_width_groups': [['disk maxima','disk means'],['line maxima','line means']],
               'source_association_counts': dict(source_counts), 'verified_specimen_ID_joins': 0,
               'new_numerical_or_biological_claim': False, 'figures_and_csv_unchanged_in_mobile_math_revision': True}
    output = HERE / 'build_receipt.json'
    write(output, json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'build_receipt_sha256': sha(output), 'math': count, 'files': receipt['file_sha256']}))


if __name__ == '__main__':
    main()
