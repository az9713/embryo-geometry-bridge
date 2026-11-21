"""Build the standalone nonlinear transfer report from hash-verified receipts."""
from fractions import Fraction as Q
from hashlib import sha256
from html.parser import HTMLParser
import json
import os
from pathlib import Path
from urllib.parse import urlparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MaxNLocator
from PIL import Image, ImageDraw

from build_report import convert


ROOT = Path(__file__).resolve().parents[1]
RUN = 'runs/transfer_ambiguity_nonlinear_probe_v1'
OUT = ROOT / RUN / 'report'
PRIMARY = f'{RUN}/actual/results.json'
REVIEW = f'{RUN}/independent_actual_preparation/actual_replay/receipt.json'
SWEEP = f'{RUN}/observer_robustness/actual/results.json'
SWEEP_REVIEW = f'{RUN}/observer_robustness/independent_actual/receipt.json'
MANIFEST = f'{RUN}/field_manifest.json'
PROTOCOL = f'{RUN}/preparation/protocol.json'
NOVELTY = f'{RUN}/originality_followup/receipt.json'
ASSETS = ROOT.parent / 'docs/assets/katex'
PINS = {
    PRIMARY: '2ddc2fb3f0da426551bf74b6e2aac7bb489dfe848f8d9e73cfc42fddecbe99d6',
    REVIEW: 'c795ae78ba0b53a01c3d65f466d5fa2f1d233c2a68c4275c04260c3cb5120ebd',
    SWEEP: '53a17296ddae08b4d1fe6b7e9b20ad0a377bf8dc962142bbfb9eb9877a381647',
    SWEEP_REVIEW: '703b574bc51a083b10c7328d9622862e939ddabe791b18ce9e1817f07718f69f',
    MANIFEST: '308ca4c05e476ba6a14cbc7b8d34018983b8f11263e2c224488672b35eb38ac0',
    PROTOCOL: '9d5d38bb5f9785d40bf843dd225eb93b1d5ecbe96f20942aec0023765d119d94',
    NOVELTY: '6797ad5d05f71535f13edc519276e57bc5176d7d5b0e6273b1aa91eba48653d0',
    f'{RUN}/actual/terminal.json': '84aafea79a60878e139ed034c199e4875e91bac21fa650354f12f944e87fd0b9',
    f'{RUN}/originality_followup/README.md': 'be18025c853273ebed98db4cbb47dac48006236631bee671090d22f279410c43',
    f'{RUN}/certification_design/README.md': '4f9ee03c3723f74b5cddecbb237424c75daed6962b331a77862866d24d3a3793',
    'runs/transfer_ambiguity_linear_v1/actual/results.json': '68ddc210523cb40001c99babb251916ac22c62eecd836d976fac1ab78dc0fe3f',
    'runs/transfer_ambiguity_linear_v1/actual_independent/receipt.json': '8fc50b42694b3db513e6d2166c6d539fbe67c66c2b3c6895c90df12f3856d072',
    'scripts/build_report.py': '32ea3c5b638fec45a58f886b8bc4678754510097a06396018988ee60566a94b8',
    '../docs/assets/katex/katex.min.css': '5a0ff74f11887f7260dc16b216906f51352b84ca29f21abb9a14b62a2ca6eb15',
    '../docs/assets/katex/katex.min.js': 'e6bfe5deebd4c7ccd272055bab63bd3ab2c73b907b6e6a22d352740a81381fd4',
}
VARIANTS = ('current_disk__spatial_max', 'current_disk__regional_mean',
            'current_radial_line__spatial_max', 'current_radial_line__regional_mean')
VARIANT_LABELS = ('Disk · component maxima', 'Disk · component means',
                  'Radial line · component maxima', 'Radial line · component means')
ENDPOINTS = ('max_contrast_change', 'max_negative_contrast_change')
BG, PANEL, TEXT, DIM = '#0b1320', '#101d2e', '#eaf3fc', '#a6b8ce'
SKY, GOLD, LILAC, GRID = '#77dcf4', '#ffca84', '#c9adff', '#3b526e'
LITERATURE = {
    'Liu, Maini and Baker (2026)': 'https://doi.org/10.1016/j.mbs.2026.109710',
    'Liu et al. author preprint': 'https://arxiv.org/abs/2506.11311',
    'Kumar et al. (2024)': 'https://www.nature.com/articles/s41540-024-00375-z.pdf',
    'Countryman et al. (2025)': 'https://www.nature.com/articles/s41467-025-62483-6',
    'Foucart and Liao': 'https://arxiv.org/abs/2111.02601',
    'DeVore et al.': 'https://foucart.github.io/publi/DFPW_recovery_final.pdf',
}


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def exact_pair(pair):
    lo, hi = map(Q, pair)
    require(lo <= hi, 'unordered exact interval')
    return lo, hi


def serialized_pair(pair):
    def endpoint(key):
        v = pair[key]
        return Q(v['numerator'], v['denominator'])
    return exact_pair((endpoint('lower'), endpoint('upper')))


def difference(a, b):
    return a[0] - b[1], a[1] - b[0]


def midpoint(pair):
    return (pair[0] + pair[1]) / 2


def open_pinned(path, digest):
    payload = (ROOT / path).read_bytes()
    require(sha256(payload).hexdigest() == digest, f'changed artifact {path}')
    return json.loads(payload)


def reviewed_values():
    for path, expected in PINS.items():
        require(sha(ROOT / path) == expected, f'changed input {path}')
    actual, review, sweep, sweep_review, manifest, protocol, novelty = [
        open_pinned(path, PINS[path]) for path in (PRIMARY, REVIEW, SWEEP, SWEEP_REVIEW, MANIFEST, PROTOCOL, NOVELTY)]
    terminal_path = f'{RUN}/actual/terminal.json'
    terminal = open_pinned(terminal_path, PINS[terminal_path])
    require(terminal['status'] == 'FINISHED' and terminal['results_sha256'] == PINS[PRIMARY]
            and len(terminal['workers']) == 2 and all(w['observed_exit_code'] == 0 for w in terminal['workers']),
            'nonlinear job lacks complete OS exits')
    require(manifest['status'] == 'COMPLETE_44_PRIMARY_CHECKPOINTS' and len(manifest['fields']) == 44
            and manifest['primary_results_sha256'] == PINS[PRIMARY]
            and manifest['primary_terminal_sha256'] == PINS[terminal_path]
            and manifest['independent_primary_review_sha256'] == PINS[REVIEW]
            and manifest['independent_replay_observed_os_exit_code'] == 0,
            'field manifest does not match the reviewed run')
    require(review['status'] == 'PASS_INDEPENDENT_ACTUAL_NONLINEAR_STORED_FIELD_REPLAY'
            and review['actual_results_sha256'] == PINS[PRIMARY] and review['checkpoints'] == 44
            and review['independently_recomputed_small_residuals'] == 136
            and actual['verified_outcome_fields'] == 36 and actual['small_observer_intervals'] == 136,
            'primary replay failed or incomplete')
    require(actual['both_endpoints_satisfy_declared_allowance'] is True
            and review['both_endpoints_satisfy_conditional_stored_allowance'] is True
            and actual['large_Q_stored_intervals_disjoint'] is True
            and all(actual[key] is False for key in ('source_witness_certified', 'exact_integral_FE_root_certified',
                                                     'continuum_accuracy_certified', 'heldout_response_tables_read')),
            'primary decision or scope drift')
    require(sweep_review['status'] == 'PASS_INDEPENDENT_SECONDARY_STORED_OBSERVER_SWEEP'
            and bool(sweep_review['checks']) and all(x is True for x in sweep_review['checks'].values())
            and sweep_review['primary_exact_result_kept'] is True and sweep['primary_replay_exactly_matches'] is True
            and sweep['additional_equilibria'] == sweep['additional_force_evaluations'] == 0
            and sweep['small_constraint_evaluations'] == 544
            and sweep_review['independent_interval_exact_matches'] == 544
            and sweep_review['robust_across_conventions_at_frozen_delta'] is False
            and novelty['novelty_confirmed'] is False, 'unexpected secondary or originality outcome')
    require(sweep_review['input_and_output_sha256'][SWEEP] == PINS[SWEEP],
            'secondary review does not bind the observer result')
    require(tuple(case['endpoint_id'] for case in actual['endpoints']) == ENDPOINTS
            and tuple(case['endpoint_id'] for case in review['endpoints']) == ENDPOINTS
            and tuple(case['endpoint_id'] for case in protocol['endpoints']) == ENDPOINTS,
            'endpoint label order drift')
    delta = Q(protocol['selection']['conditional_source_allowance_um_exact'])
    affine = Q(review['half_step_affine_Q_difference_um_exact'])
    require(delta == Q(sweep_review['source_configured_delta_um_exact'])
            == Q(sweep['frozen_primary_source_allowance_um_exact'])
            and affine == Q(protocol['selection']['half_pair_affine_large_Q_range_um_exact']),
            'changed affine comparison or declared allowance')
    q = [exact_pair(case['large_Q_um_exact']) for case in actual['endpoints']]
    gap = difference(*q)
    offset = difference(gap, (affine, affine))
    required = [exact_pair(case['small_residual']['required_allowance_um_exact']) for case in actual['endpoints']]
    require(all(case['small_residual']['all68_residuals']['constraints'] == 68 for case in actual['endpoints'])
            and required == [exact_pair(case['required_allowance_um_exact']) for case in review['endpoints']]
            and q == [exact_pair(case['large_Q_um_exact']) for case in review['endpoints']]
            and gap == exact_pair(review['stored_large_Q_signed_difference_um_exact'])
            == exact_pair(actual['large_Q_positive_direction_minus_negative_direction_um_exact'])
            and offset == exact_pair(review['nonlinear_minus_affine_signed_difference_um_exact']),
            'independent and primary exact quantities disagree')
    common = [exact_pair(sweep['common_pair_requirement_by_variant'][v]['required_allowance_um_exact'])
              for v in VARIANTS]
    require(common[0] == exact_pair(actual['both_endpoints_required_allowance_um_exact'])
            and common[0][1] < delta < min(p[0] for p in common[1:])
            and exact_pair(sweep_review['all_four_required_common_allowance_um_exact']) == (max(p[0] for p in common),
                                                                                           max(p[1] for p in common)),
            'primary differs from secondary re-evaluation or exact outcome changed')
    profile_rows = [row for row in manifest['fields'] if row['role'] == 'large_final']
    require([row['endpoint_id'] for row in profile_rows] == list(ENDPOINTS), 'wrong large shapes')
    profiles = []
    chosen = {}
    shoulder = Q(protocol['large_observable']['shoulder_radius_um_exact'])
    for index, row in enumerate(profile_rows):
        for kind in ('field', 'verification'):
            path, digest = row[kind], row[f'{kind}_sha256']
            require(review['all_hashes'][path] == digest and sha(ROOT / path) == digest,
                    'large field or curve certificate hash drift')
            chosen[path] = digest
        proof = open_pinned(row['verification'], row['verification_sha256'])['surface']
        require(proof['usable_stored_graph'] is True and exact_pair(proof['Q_um_exact']) == q[index],
                'unverified large stored profile')
        sample = {Q(x['radius_um_exact']): serialized_pair(x['apical_height_um'])
                  for x in proof['sampled_full_current_profiles']}
        require(Q(proof['shoulder_radius_um_exact']) == shoulder and difference(sample[0], sample[shoulder]) == q[index],
                'profile samples do not produce exact Q at the shared shoulder')
        profiles.append(sample)
    samples = [r for branch in review['endpoints'] for r in branch['checkpoints']]
    return {
        'protocol': protocol, 'delta': delta, 'affine': affine, 'q': q,
        'gap': gap, 'offset': offset, 'required': required, 'common': common,
        'margin': delta - common[0][1], 'profiles': profiles, 'shoulder': shoulder,
        'chosen_fields_and_curve_certificates': chosen,
        'checks': {'Q8': max(r['Q8'] for r in samples), 'Q12': max(r['Q12'] for r in samples),
                   'J': min(Q(r['J_lower_exact']) for r in samples),
                   'recomputed_force_evaluations': review['numerical_force_evaluations'],
                   'element_geometry_checks': review['stored_exact_element_checks']},
    }


def start_plot(rows=1, figsize=(8.0, 4.0)):
    fig, axes = plt.subplots(rows, 1, figsize=figsize, layout='constrained', squeeze=False)
    fig.patch.set_facecolor(BG)
    for ax in axes[:, 0]:
        ax.set_facecolor(PANEL)
        ax.tick_params(colors=DIM, labelsize=11, length=3)
        ax.xaxis.label.set_color(TEXT)
        ax.yaxis.label.set_color(TEXT)
        ax.grid(color=GRID, alpha=.35, linewidth=.65)
        for spine in ax.spines.values():
            spine.set_color(GRID)
        ax.set_axisbelow(True)
    return fig, axes[:, 0]


def write_plot(fig, filename):
    path = OUT / filename
    with path.open('xb') as handle:
        fig.savefig(handle, format='png', dpi=180, facecolor=BG)
    plt.close(fig)


def create_figures(v):
    fig, (ax, marks) = start_plot(2, (8.2, 6.4))
    signed = float(midpoint(v['gap']))
    affine = float(v['affine'])
    ax.set_title('Signed large-patch contrast gap', loc='left', color=TEXT, pad=12, fontsize=15)
    ax.scatter((signed, affine), (1, 0), marker='o', c=(SKY, LILAC), s=125, zorder=3)
    ax.hlines((1, 0), .15472, (signed, affine), colors=(SKY, LILAC), alpha=.5, linewidth=1.3)
    ax.set_yticks((1, 0), ('Nonlinear stored fields', 'Frozen affine half-steps'))
    ax.set(xlim=(.15472, .15496), ylim=(-.52, 1.54), xlabel='Q₊ − Q₋ (µm) · detail scale')
    ax.xaxis.set_major_locator(MaxNLocator(5))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.5f'))
    ax.annotate(f'{signed:.9f}', (signed, 1), xytext=(11, 3), textcoords='offset points', color=SKY, fontsize=11)
    ax.annotate(f'{affine:.9f}', (affine, 0), xytext=(-10, 4), textcoords='offset points',
                color=LILAC, ha='right', fontsize=11)
    marks.set_title('Primary small-region requirement', loc='left', color=TEXT, pad=12, fontsize=15)
    delta = float(v['delta'])
    req = [float(x[1]) for x in v['required']]
    marks.axvspan(5.258, delta, color=SKY, alpha=.07)
    marks.axvline(delta, color=GOLD, linestyle='--', linewidth=2)
    for y, value in zip((1, 0), req):
        marks.hlines(y, 5.258, value, color=SKY, alpha=.55)
        marks.scatter(value, y, color=SKY, s=95, zorder=3)
        marks.annotate(f'{value:.6f}', (value, y), textcoords='offset points', xytext=(-9, 10),
                       ha='right' if y == 0 else 'left', color=SKY, fontsize=11)
    marks.annotate(f'δ = {delta:.6f}', (delta, 1.26), xytext=(-5, 0), textcoords='offset points',
                   ha='right', color=GOLD, fontsize=10)
    marks.set_yticks((1, 0), ('Endpoint +', 'Endpoint −'))
    marks.set(xlim=(5.258, 5.427), ylim=(-.52, 1.65), xlabel='Maximum absolute residual (µm)')
    write_plot(fig, 'contrast_and_allowance.png')

    fig, (ax,) = start_plot(1, (8.2, 4.25))
    excess = [float(midpoint(pair) - v['delta']) for pair in v['common']]
    positions = list(range(4))
    ax.barh(positions, excess, color=(SKY, GOLD, GOLD, GOLD), height=.42)
    ax.scatter(excess, positions, color=(SKY, GOLD, GOLD, GOLD), s=30, zorder=3)
    ax.axvline(0, color=TEXT, linewidth=1.6)
    ax.set_yticks(positions, VARIANT_LABELS)
    ax.invert_yaxis()
    ax.set(xlim=(-.19, 1.39), xlabel='Common required allowance minus frozen δ (µm)')
    ax.set_title('Same 34 small fields, four observer conventions', loc='left', color=TEXT, pad=13, fontsize=15)
    for i, x in enumerate(excess):
        ax.annotate(f'{x:+.6f}', (x, i), xytext=(8, 0), textcoords='offset points',
                    color=SKY if i == 0 else GOLD, va='center', fontsize=10)
    write_plot(fig, 'observer_dependence.png')

    fig, (ax, diff) = start_plot(2, (8.2, 6.4))
    plus, minus = v['profiles']
    a = v['shoulder']
    step_radius = sorted(r for r in set(plus) & set(minus) if r <= 54)
    require(step_radius[0] == 0 and a in step_radius, 'shared radius registration not represented')
    x = [float(r) for r in step_radius]
    for sample, label, color in ((plus, 'Endpoint +', SKY), (minus, 'Endpoint −', GOLD)):
        h = [float(midpoint(sample[r])) for r in step_radius]
        ax.plot(x, h, color=color, label=label, linewidth=2.3)
        ax.scatter((0, float(a)), (h[0], float(midpoint(sample[a]))), color=color, s=32, zorder=3)
    ax.set_title('Stored apical profiles under the matched large input', loc='left', color=TEXT, pad=12, fontsize=15)
    ax.set_ylabel('Apical height (µm)')
    ax.legend(facecolor=PANEL, edgecolor=GRID, labelcolor=TEXT, loc='lower right', fontsize=10)
    hz = [float(midpoint(difference(plus[r], minus[r]))) for r in step_radius]
    diff.plot(x, hz, color=LILAC, linewidth=2.2)
    diff.axhline(0, color=GRID)
    diff.set(title='', xlabel='Current radius (µm)', ylabel='Height difference (µm)')
    for panel in (ax, diff):
        panel.axvline(float(a), color=DIM, linestyle='--', linewidth=1)
        panel.set_xlim(-.7, 54)
    diff.set_title('Endpoint + minus endpoint − at the same current radius', loc='left', color=TEXT, pad=12, fontsize=13)
    for r in (Q(0), a):
        z = float(midpoint(difference(plus[r], minus[r])))
        diff.scatter(float(r), z, color=LILAC, s=38, zorder=3)
        diff.annotate(f'{z:+.5f}', (float(r), z), xytext=(7, 4), textcoords='offset points',
                      color=TEXT, fontsize=10)
    write_plot(fig, 'large_current_radius_profiles.png')


def md_report(v):
    p = v['protocol']
    values = {
        'GAP': f"{float(midpoint(v['gap'])):.12f}", 'AFFINE': f"{float(v['affine']):.12f}",
        'OFFSET': f"{float(midpoint(v['offset'])):+.12f}", 'DELTA': f"{float(v['delta']):.11f}",
        'MARGIN': f"{float(v['margin']):.11f}", 'SHOULDER': f"{float(v['shoulder']):.8f}",
        'DIAMETER': f"{float(Q(p['matched_input']['large_reported_diameter_um_exact'])):.7f}",
        'DOSE': f"{float(Q(p['matched_input']['matched_dose_exact'])):.8f}",
        'BASE': f"{float(Q(p['selection']['affine_baseline_source_residual_um_exact'])):.11f}",
        'Q8': f"{v['checks']['Q8']:.3e}", 'Q12': f"{v['checks']['Q12']:.3e}", 'J': f"{float(v['checks']['J']):.6f}",
        'INTERVALWIDTH': f"{float(v['gap'][1] - v['gap'][0]):.2e}",
        'MECHCOUNT': str(v['checks']['recomputed_force_evaluations']),
        'GEOMCOUNT': f"{v['checks']['element_geometry_checks']:,}",
        'ALLFOUR': f"{float(max(pair[1] for pair in v['common'])):.10f}",
    }
    result_table = []
    for name, la, pa in zip(('+', '−'), v['q'], v['required']):
        result_table.append(f'| {name} | {float(midpoint(la)):+.11f} | {float(pa[1]):.11f} | '
                            f'{float(v["delta"] - pa[1]):.11f} |')
    values['RESULT_TABLE'] = '\n'.join(result_table)
    values['OBSERVER_TABLE'] = '\n'.join(
        f'| {label} | {float(interval[1]):.10f} | {float(interval[1] - v["delta"]):+.10f} | '
        f'{("Both satisfy" if interval[1] <= v["delta"] else "At least one exceeds")} |'
        for label, interval in zip(VARIANT_LABELS, v['common']))
    values['STEP_TABLE'] = '\n'.join(
        '| ' + name + ' | ' + ' | '.join(f'{float(Q(z)):+.8f}' for z in case['half_step_exact']) + ' |'
        for name, case in zip(('+','−'), p['endpoints']))
    doc = r'''# Nonlinear transfer ambiguity in one epithelial model

Two frozen parameter choices satisfy all 68 small-region entries per choice under the declared **disk-weighted, spatial-maximum** observer. At the same large optical input, their saved apical curves differ in signed center-to-shoulder contrast by **@@GAP@@ µm**. The result supports a narrow claim about two numerical fields under one explicit observation and residual contract.

**Observer choice changes the conclusion.** The common allowance is **@@DELTA@@ µm**, and the smaller remaining margin is **@@MARGIN@@ µm**. Each of the other three predeclared static observer conventions needs a larger allowance for the same pair. The allowance is a conditional model discrepancy threshold derived from an affine residual. Experimental noise, timing, reference and registration remain separately unresolved. [Primary result](../actual/results.json) · [Independent primary replay](../independent_actual_preparation/actual_replay/receipt.json) · [Independent observer sweep](../observer_robustness/independent_actual/receipt.json)

![Two zoomed panels compare the signed nonlinear and affine gaps, and put both endpoint residual requirements below the frozen allowance.](contrast_and_allowance.png)

## Intervention and parameter selection

Both candidates use the same axisymmetric P2 finite-element law and geometric conventions. They have distinct material coefficients and one common gain *within each candidate* across its 17 small doses and its large intervention. The matched large optical dose is **@@DOSE@@**, the diameter is **@@DIAMETER@@ µm**, and the current-radius shoulder is **@@SHOULDER@@ µm**. Optical dose is held common between candidates, while their mechanical activation strength is the exact candidate gain times that dose, rounded once for the numerical solver.

The frozen protocol selected two ends of the local affine contrast frontier at \(\delta=D_{\mathrm{aff},0}+0.01\,\mu\mathrm m\), where \(D_{\mathrm{aff},0}\approx @@BASE@@\,\mu\mathrm m\) is the archived baseline maximum residual. Each exact selected perturbation was halved before any new nonlinear calculation. Labels \(+\) and \(-\) identify the affine contrast directions. [Affine derivation](../../transfer_ambiguity_linear_v1/REPORT.html) · [Frozen nonlinear protocol](../preparation/protocol.json)

| Endpoint | \(\Delta\log K_a\) | \(\Delta\log K_b\) | \(\Delta T_a\) | \(\Delta\log T_b\) | \(\Delta\log g\) |
|---|---:|---:|---:|---:|---:|
@@STEP_TABLE@@

The small cohort was computed in ascending dose order from a fresh zero equilibrium for each candidate. The large cases each began from a separate new zero and used one-third, two-thirds and full dose checkpoints. This procedure produced 34 small final fields, two large final fields and eight supplementary loading checkpoints. No large experimental response entered selection or scoring. [Whole field manifest](../field_manifest.json)

## Exact comparisons of the saved curves

For record \(i\) and static observable \(k\), let \(y_{ik}\) be the recorded small-region value and \(I_{ik,j}^{s}\) the exact interval for the saved curve of candidate \(j\in\{+,-\}\). Compute

\[
D^{s}_{j}=\max_{i,k}\ \sup_{x\in I_{ik,j}^{s}} |x-y_{ik}|,\qquad
D^{s}_{\mathrm{pair}}=\max(D^{s}_{+},D^{s}_{-}).
\]

The supremum of absolute residual on each closed interval is its larger endpoint residual. The decision uses all 68 rows per candidate and exact rational endpoints. The table displays rounded midpoints of the large contrast intervals and rounded upper bounds on the required allowances. The margin column is \(\delta-D^{s}_{j}\).

| Endpoint | Stored large \(Q\), µm | Required allowance upper bound, µm | Remaining margin, µm |
|---|---:|---:|---:|
@@RESULT_TABLE@@

For each large saved curve, write \(z_a(R)\) for apical height at current physical radius \(R\). The signed statistic is

\[
Q=z_a(0)-z_a(a),\qquad a=@@SHOULDER@@\,\mu\mathrm m .
\]

A negative \(Q\) indicates a center below the shoulder. The actual positive signed gap persists at **@@GAP@@ µm**, compared with the frozen affine half-step figure **@@AFFINE@@ µm**. The actual minus affine gap is **@@OFFSET@@ µm**. This is agreement at the two selected endpoints. The experiment supplies no bound on the affine approximation over the rest of the parameter domain. [Exact gap and difference intervals](../independent_actual_preparation/actual_replay/receipt.json)

![Stored apical profiles and their pointwise height difference, both on current radius with the common shoulder marked.](large_current_radius_profiles.png)

The graph uses midpoints of the retained exact curve evaluations. If \(\Delta z(R)=z_{a,+}(R)-z_{a,-}(R)\), then \(Q_+-Q_-=\Delta z(0)-\Delta z(a)\). The lower panel illustrates that the gap is a difference between two registered heights, each of which also changes across the pair.

## Dependence on the static observer

All four conventions use the same saved small fields, recorded values, flat reference, current-radius region and total maximum shortening. For the regional mean, with \(0\le R\le R_0\), the normalized disk and radial-line weights are \(2R/R_0^2\) and \(1/R_0\), respectively. The two separate apical and basal component values are either spatial maxima or regional means. Separate maxima can occur at different radii. The maximum of the total shortening therefore has its own calculation.

| Static observer | Required common allowance, µm | Excess over frozen \(\delta\), µm | Pair at frozen \(\delta\) |
|---|---:|---:|---|
@@OBSERVER_TABLE@@

![Signed excess above the frozen allowance for four observers. Only disk weighting with component maxima has negative excess.](observer_dependence.png)

An independent sweep reproduced all 544 scalar comparisons from the same 34 fields. The all-four common requirement is **@@ALLFOUR@@ µm** in this sweep. Passing the primary observer provides limited evidence when the source-to-model observation rule remains uncertain. The alternatives here affect static weighting and aggregation. The unresolved source contract also includes timing of reported extrema, the reference for growth and the mapping of the region of interest. [Exact secondary results](../observer_robustness/actual/results.json) · [Independent secondary audit](../observer_robustness/independent_actual/receipt.json) · [Observation contract](../../../docs/JOINT-OBSERVATION-CONTRACT.html)

## Numerical evidence and the missing error accounts

The independent replay verified all 44 checkpoints, recomputed @@MECHCOUNT@@ numerical force evaluations and checked @@GEOMCOUNT@@ discrete element geometries. The largest replayed scaled force infinity norms were **@@Q8@@** for Q8 and **@@Q12@@** for Q12, under the \(10^{-7}\) gate. The smallest exact stored lower Jacobian bound was **@@J@@**. The stored fields passed exact positive orientation and embedding checks. Q8 and Q12 refer to the main and cross-check integration policies, respectively. They assess residual at the saved coordinates in two quadrature rules.

The stored signed gap interval has width about **@@INTERVALWIDTH@@ µm**. This resolves evaluation of the serialized quadratic curves. To bound a stationary point of the exactly integrated finite-element problem, a separate proof must enclose both the signed residual and a valid inverse or contraction on a neighborhood. An observer transfer theorem can then turn a certified nodal displacement bound, including binary addition gaps, into an error on the plotted contrast and all 136 small comparisons. The linked design allocates targets from the stored slacks and includes strict radial-derivative margins. Those targets constrain future certificates; actual root-position, continuum and source errors remain open for this pair. [Pointwise certification design](../certification_design/README.md) · [Nodal observer theorem](../../../docs/NODAL-OBSERVER-ERROR-TRANSFER.html)

Practical discrimination requires a registered measurement of the same \(Q\) plus a budget for the experimental, geometric, temporal and model errors. For two certified response intervals separated by \(\gamma>0\), any single scalar estimate has conditional worst-case absolute error at least \(\gamma/2\). This follows from the triangle inequality. Additional deterministic error radii \(\eta_+,\eta_-\) preserve separation if \(\eta_++\eta_-<\gamma\). Here the stored gap establishes neither a bound on those radii nor a necessary instrument specification.

## Originality and scope

The general same-law discrimination construction is established. Liu, Maini and Baker develop interventions that distinguish candidate parameter sets and give a logistic-model example. Kumar et al. find multiple epithelial parameter clusters that reproduce similar control morphology. Countryman et al. already establish activation-size-dependent furrowing and examine basal mechanics. Deterministic, bounded-error estimation of a quantity of interest has an existing optimal-recovery theory. [Liu et al. (2026)](https://doi.org/10.1016/j.mbs.2026.109710) · [Kumar et al. (2024)](https://www.nature.com/articles/s41540-024-00375-z.pdf) · [Countryman et al. (2025)](https://www.nature.com/articles/s41467-025-62483-6) · [Foucart and Liao](https://arxiv.org/abs/2111.02601) · [DeVore et al.](https://foucart.github.io/publi/DFPW_recovery_final.pdf)

The narrow contribution is a reproducible, predeclared nonlinear pair for this discrete epithelial model with documented dependence on the observation convention. A bounded literature search cannot establish priority. A source-compatible mechanical witness, exact equilibrium error, a continuous loading branch and practically attainable measurement discrimination remain open. [Originality follow-up and search locators](../originality_followup/README.md) · [Originality audit receipt](../originality_followup/receipt.json) · [Build provenance](report_build_receipt.json)
'''
    for key, text in values.items():
        doc = doc.replace(f'@@{key}@@', text)
    require('@@' not in doc, 'unexpanded report placeholders')
    return doc


STYLE = """
*{box-sizing:border-box}body{margin:0;background:#0b1320;color:#eaf3fc;font:17px/1.72 system-ui,sans-serif}
main{max-width:1040px;margin:auto;padding:34px clamp(16px,4vw,44px) 78px}
nav{display:flex;flex-wrap:wrap;gap:8px 24px;border-bottom:1px solid #344963;padding-bottom:16px;font-size:.9em}
a{color:#84dbf5;text-underline-offset:.19em}h1,h2{line-height:1.16}
h1{font-size:clamp(2.1rem,5.3vw,3.65rem);max-width:21ch;margin:30px 0 20px}
h2{color:#b9e8f5;margin-top:2.25em;font-size:clamp(1.35rem,3vw,1.92rem)}
p{margin:1em 0;overflow-wrap:anywhere}h1+p{font-size:1.09em}
h1+p+p{border-left:4px solid #ffca84;background:#142236;padding:16px 18px}
img{width:100%;height:auto;display:block;border:1px solid #344963;border-radius:10px;margin:1.6em 0}
table{position:relative;border-collapse:collapse;display:block;max-width:100%;overflow-x:auto;font-size:.92em;font-variant-numeric:tabular-nums;margin:1.25em 0}
th,td{padding:9px 12px;border:1px solid #344963;vertical-align:top}th{background:#142236;color:#eaf3fc}
.math.display{overflow-x:auto;background:#142236;border-radius:6px;padding:10px}.katex-display{margin:.4em 0}
#math-errors:empty{display:none}#math-errors{color:#ffd2a0}
@media(max-width:520px){body{font-size:16px}main{padding-top:20px}th,td{padding:7px 9px}}
"""


def create_html(markdown_text):
    def link(path):
        return Path(os.path.relpath(path, OUT)).as_posix()
    return (
        '<!doctype html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1">'
        '<title>Nonlinear transfer ambiguity in one epithelial model</title>'
        '<link rel="icon" type="image/png" href="favicon.png" sizes="32x32">'
        f'<link rel="stylesheet" href="{link(ASSETS / "katex.min.css")}"><style>{STYLE}</style>'
        '</head><body><main><nav aria-label="Related documents">'
        f'<a href="{link(ROOT / "docs/index.html")}">Controlled-furrowing dossier</a>'
        '<a href="../../transfer_ambiguity_linear_v1/REPORT.html">Affine study</a>'
        '<a href="REPORT.md">Markdown</a></nav>'
        + convert(markdown_text)
        + '<p id="math-errors" role="status"></p></main>'
        f'<script src="{link(ASSETS / "katex.min.js")}"></script>'
        '<script>document.querySelectorAll(".math").forEach(e=>{try{'
        'katex.render(e.dataset.tex,e,{displayMode:e.classList.contains("display"),throwOnError:true});'
        '}catch(x){document.getElementById("math-errors").textContent+=String(x)+"\\n";}});</script></body></html>'
    )


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.paths = []
        self.math_count = 0

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        self.math_count += int('data-tex' in a)
        for k in ('src', 'href'):
            if k in a:
                self.paths.append(a[k])


def main():
    v = reviewed_values()
    OUT.mkdir(exist_ok=False)
    create_figures(v)
    favicon = Image.new('RGBA', (32, 32), BG)
    draw = ImageDraw.Draw(favicon)
    draw.line([(4, 21), (9, 19), (16, 10), (22, 13), (28, 9)], fill=SKY, width=3)
    draw.line([(4, 25), (9, 22), (16, 15), (22, 17), (28, 12)], fill=GOLD, width=2)
    with (OUT / 'favicon.png').open('xb') as handle:
        favicon.save(handle, format='PNG')
    text = md_report(v)
    page = create_html(text)
    for name, content in (('REPORT.md', text), ('REPORT.html', page)):
        with (OUT / name).open('x', encoding='utf-8', newline='\n') as handle:
            handle.write(content)
    links = Links()
    links.feed(page)
    files = {p.relative_to(ROOT).as_posix(): sha(p) for p in sorted(OUT.iterdir())}
    exact_quantities = {
        'delta_um_exact': str(v['delta']), 'signed_gap_um_exact': list(map(str, v['gap'])),
        'nonlinear_minus_affine_gap_um_exact': list(map(str, v['offset'])),
        'affine_gap_um_exact': str(v['affine']), 'minimum_primary_margin_um_exact': str(v['margin']),
        'large_Q_by_endpoint_um_exact': [list(map(str, x)) for x in v['q']],
        'primary_required_by_endpoint_um_exact': [list(map(str, x)) for x in v['required']],
        'same_field_observer_requirements_um_exact': {k: list(map(str, x)) for k, x in zip(VARIANTS, v['common'])},
    }
    receipt = {
        'status': 'BUILT_FROM_PINNED_INDEPENDENTLY_REVIEWED_RESULTS',
        'source_sha256': sha(Path(__file__)), 'source': Path(__file__).relative_to(ROOT).as_posix(),
        'input_hashes': {**PINS, **v['chosen_fields_and_curve_certificates']},
        'report_and_figure_hashes': files, 'primary_source_links': LITERATURE,
        'figure_data_and_report_rounding_source': exact_quantities,
        'static_checks': {'local_asset_links_verified': True, 'math_spans_detected': links.math_count,
                          'profile_samples_from_reviewed_current_radius_receipts': True,
                          'browser_review_pending': True},
        'additional_equilibria_computed': 0, 'heldout_large_responses_read': False,
    }
    pending_link = (OUT / 'report_build_receipt.json').resolve()
    for href in links.paths:
        url = urlparse(href)
        if not url.scheme and url.path:
            path = (OUT / url.path).resolve()
            require(path.is_file() or path == pending_link, f'broken relative link: {href}')
    require(links.math_count >= 10, 'report lost math markup')
    with pending_link.open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(receipt, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({'status': receipt['status'], 'html': str(OUT / 'REPORT.html'),
                      'build_receipt_sha256': sha(pending_link), 'files': files}))


if __name__ == '__main__':
    main()

