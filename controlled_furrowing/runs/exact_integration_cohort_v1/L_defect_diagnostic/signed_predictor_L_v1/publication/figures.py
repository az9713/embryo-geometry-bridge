"""Render two certificate figures by exact arithmetic on retained L reviews."""

from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, NullLocator


OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[4]
COHORT = 'runs/exact_integration_cohort_v1'
DIAGNOSTIC = f'{COHORT}/L_defect_diagnostic'
SIGNED = f'{DIAGNOSTIC}/signed_predictor_L_v1'
FILES = {
    'direct': (f'{COHORT}/reduction_preparation/optoGEF2_small_col_L_reduction_review.json',
               '5a55e0d139edb39206ab3999681d9858d6b0162f73c7e5103d7a5bed3ebc4514'),
    'scaled': (f'{DIAGNOSTIC}/scaled_offaxis_L_v1/reduction_preparation/actual_independent_review.json',
               '54dee8659601c387712c16732ab6aed90fbf313c1e3b0f77930b3e12b19220e6'),
    'centered': (f'{DIAGNOSTIC}/centered_bulk_L_v1/reduction_preparation/actual_independent_review.json',
                 'cce6f01147b6ed436e3f3a8620e94bbd9df84a8fc1aa45e9bbed9cba8c1f8c4a'),
    'signed': (f'{SIGNED}/actual_independent_review.json',
               'ea1ebcb2f1788eaa3fd491cca6d519e8f40568e9cdc870919a9b32f9f4adc0a9'),
    'independent_256': (f'{SIGNED}/independent_actual/results.json',
                        'bdfa3b85849c47b7d3a80424c8f407711d6531ab7a43b9e5f117abd439a7a5d3'),
    'readback': (f'{SIGNED}/final_mathematical_readback/receipt.json',
                 '524a6bf025b9efff519a7d16d0df57d3901b9f35635ba57d0ff05cccf021f079'),
}
BACKGROUND = '#111921'
PANEL = '#192632'
WHITE = '#f1f5fa'
MUTED = '#afc0ce'
GRID = '#3c4d5a'
AMBER = '#efb476'
CYAN = '#63d3d5'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def all_true(checks):
    return bool(checks) and all(value is True for value in checks.values())


def checked_inputs():
    receipts = {}
    for label, (filename, expected) in FILES.items():
        path = ROOT / filename
        actual = sha(path)
        if actual != expected:
            raise ValueError(f'Changed archive binding for {filename}: {actual}')
        receipts[label] = json.loads(path.read_text(encoding='utf-8'))
    readback = receipts['readback']
    for path, expected in readback['sha256'].items():
        if sha(ROOT / path) != expected:
            raise ValueError(f'Changed final readback premise: {path}')
    return receipts


def archived_bounds(row):
    return {'eta': Q(row['eta_upper_exact']), 'z': Q(row['z_upper_exact']),
            'lhs': Q(row['selfmap_lhs_exact']), 'r': Q(row['normal_radius_exact'])}


def proof_metrics(receipts):
    reviews = receipts['direct'], receipts['scaled']['result'], receipts['centered']['decision']
    independent = receipts['independent_256']
    reviewed = receipts['signed']
    readback = receipts['readback']
    stages = {name: archived_bounds(row) for name, row in zip(('direct', 'scaled', 'centered'), reviews)}
    signed = {'eta': Q(independent['independent_eta_I_exact']), 'z': Q(independent['z_exact']),
              'lhs': Q(independent['lhs_exact']), 'r': Q(independent['r_exact'])}
    stages['signed'] = signed
    for stage in stages.values():
        stage['ratio'] = stage['lhs'] / stage['r']
    eta, z, r = (signed[key] for key in ('eta', 'z', 'r'))
    budget = Q(independent['strict_signed_force_budget_exact'])
    error = Q(independent['root_error_upper_exact'])
    z_cached = Q(independent['old_Z_Q_exact'])
    b = Q(independent['b_exact'])
    norm = Q(independent['R_infinity_norm_exact'])
    curvature = Q(independent['strict_local_minimum_curvature']['normalized_curvature_lower_exact'])
    checks = {
        'same_record_and_passed_final_reviews': (
            reviewed['status'] == 'PASS_INDEPENDENT_SIGNED_ORIGINAL_CUBE_ROOT_CERTIFICATE'
            and all_true(reviewed['checks']) and all_true(readback['checks'])
            and readback['status'] == 'PASS_FINAL_POSTRUN_MATHEMATICAL_READBACK'
            and reviewed['record_id'] == readback['record_id'] == independent['record_id'] == 'optoGEF2_small_col_L'),
        'same_closed_cell': reviewed['same_closed_C_sha256'] == readback['parameter_cell_sha256']
        == independent['same_closed_C_sha256'],
        'independent_256_bit_source_and_exact_metrics_match_review': (
            reviewed['independent_actual_results_sha256'] == FILES['independent_256'][1]
            and independent['precision_bits'] == 256 and
            {k: Q(reviewed[v]) for k, v in {'eta': 'uniform_eta_I_upper_exact', 'z': 'same_box_z_exact',
                                            'lhs': 'uniform_selfmap_lhs_exact', 'r': 'original_radius_exact'}.items()}
            == {k: signed[k] for k in ('eta', 'z', 'lhs', 'r')}),
        'same_radius_and_correct_self_map_arithmetic': all(stage['r'] == r and stage['r'] > 0
            and stage['lhs'] == stage['eta'] + stage['z'] * stage['r'] for stage in stages.values()),
        'only_final_test_passes': (all(row['sufficient_conditions_pass'] is False for row in reviews)
            and min(stages[k]['ratio'] for k in ('direct', 'scaled', 'centered')) > 1 > signed['ratio']
            and 0 <= z < 1 and reviewed['sufficient_test_pass'] is True),
        'same_derivative_budget_after_scaled_stage': stages['scaled']['z'] == stages['centered']['z'] == z
        and z == z_cached + b,
        'strict_force_budget_between_two_centered_residual_bounds': eta < budget < stages['centered']['eta']
        and budget == r * (1 - z),
        'reviewed_error_and_curvature_reproduce': (
            error == eta / (1 - z) == Q(reviewed['root_error_upper_exact']) < r
            and curvature == (1 - z) / norm == Q(reviewed['uniform_curvature_lower_exact']) > 0),
        'scope_and_all_signed_rows_retained': (
            len(independent['independent_full_RF_predictor_intervals']) == 611
            and len(independent['independent_full_signed_anchors']) == 611
            and len(independent['shared_parameter_offsets']) == 5
            and reviewed['site_domain_and_normalization']['site_counts'] == {'bulk': 9216, 'surface': 288}
            and reviewed['root_certified'] is True and reviewed['strict_local_minimum_certified'] is True
            and readback['full_five_gate_goal_completed'] is False),
    }
    if not all_true(checks):
        raise ValueError(f'Certificate figure inputs failed derivation: {checks}')
    return stages, {'r': r, 'strict_eta_budget': budget, 'root_error': error, 'curvature': curvature,
                    'selfmap_margin': r - signed['lhs'], 'Z_Q': z_cached, 'b': b}, checks


def style(ax):
    ax.set_facecolor(PANEL)
    ax.tick_params(length=0, pad=8, colors=MUTED)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.grid(axis='x', linewidth=0.8, color=GRID)
    ax.set_axisbelow(True)


def stage_plot(stages):
    fig, ax = plt.subplots(figsize=(11.5, 5.6), facecolor=BACKGROUND)
    fig.subplots_adjust(left=0.29, right=0.94, top=0.75, bottom=0.22)
    fig.text(0.065, 0.91, 'The L self-map test crosses the threshold', color=WHITE, fontsize=22, weight='bold')
    fig.text(0.065, 0.842, 'Independent upper bounds on the same original predictor cube',
             color=MUTED, fontsize=11.7)
    style(ax)
    ax.set_xscale('log')
    ax.set_xlim(0.52, 27000)
    ax.set_ylim(3.55, -0.55)
    ax.set_xticks([1, 10, 100, 1000, 10000], ['1', '10', '100', '1,000', '10,000'])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.axvspan(0.52, 1, color=CYAN, alpha=0.08)
    ax.axvline(1, color=WHITE, alpha=0.8, ls='--', lw=1.1)
    ax.set_yticks(range(4), ['Direct integration defect', 'Scaled volume coordinates',
                             'Centered bulk, separate norms', 'Signed full predictor residual'])
    ax.set_xlabel(r'Self-map bound divided by radius  $(\eta+zr)/r$', color=MUTED, labelpad=11)
    formats = (',.3f', '.5f', '.5f', '.6f')
    for index, (stage, precision) in enumerate(zip(stages.values(), formats)):
        signed = index == 3
        ratio = float(stage['ratio'])
        color = CYAN if signed else AMBER
        ax.scatter(ratio, index, color=color, s=100, zorder=4)
        ax.annotate(f'{ratio:{precision}}', (ratio, index), xytext=(0, 11), textcoords='offset points',
                    ha='center', color=color, fontsize=11, weight='bold')
    fig.text(0.065, 0.057, 'A ratio below 1 gives strict self-mapping. Above 1, the upper-bound test leaves existence open.',
             color=MUTED, fontsize=10)
    path = OUT / 'L_self_map_stages.png'
    fig.savefig(path, dpi=163, facecolor=BACKGROUND)
    plt.close(fig)
    return path


def force_plot(stages, metrics):
    fig, ax = plt.subplots(figsize=(11.5, 4.9), facecolor=BACKGROUND)
    fig.subplots_adjust(left=0.285, right=0.93, top=0.71, bottom=0.26)
    fig.text(0.065, 0.90, 'Preserving signs makes the force bound fit', color=WHITE, fontsize=22, weight='bold')
    fig.text(0.065, 0.826, 'Same parameter cell, radius, and whole-cube derivative bound', color=MUTED, fontsize=11.7)
    style(ax)
    scale = 10**11
    budget = float(metrics['strict_eta_budget'] * scale)
    ax.set_xlim(0, 5.55)
    ax.set_ylim(1.65, -0.65)
    ax.set_yticks([0, 1], ['Separate force norms', 'Signed rowwise sum'])
    ax.axvspan(0, budget, color=CYAN, alpha=0.08)
    ax.axvline(budget, color=WHITE, ls='--', lw=1.1)
    ax.annotate(f'Strict budget  {budget:.4f}', (budget, 1.63), xytext=(0, -3), textcoords='offset points',
                ha='center', color=WHITE, fontsize=10)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:g}'))
    ax.set_xlabel(r'Uniform bound on $\|RF_I(p(\theta),\theta)\|_\infty$  $(10^{-11})$', color=MUTED, labelpad=10)
    for index, (stage, color) in enumerate((('centered', AMBER), ('signed', CYAN))):
        value = float(stages[stage]['eta'] * scale)
        ax.barh(index, value, color=color, height=0.48, alpha=0.83)
        ax.annotate(f'{value:.4f}', (value, index), xytext=(8, 0), textcoords='offset points',
                    va='center', color=color, fontsize=11, weight='bold')
    fig.text(0.065, 0.078, 'The signed bound combines cached and integration-defect terms before taking absolute values.',
             color=MUTED, fontsize=10)
    fig.text(0.065, 0.038, 'The bars are rigorous upper bounds in normalized nodal coordinates, displayed to four decimals.',
             color=MUTED, fontsize=10)
    path = OUT / 'L_signed_force_budget.png'
    fig.savefig(path, dpi=163, facecolor=BACKGROUND)
    plt.close(fig)
    return path


def main():
    receipts = checked_inputs()
    stages, metrics, checks = proof_metrics(receipts)
    source = Path(__file__).resolve()
    protocol = {
        'source': source.relative_to(ROOT).as_posix(), 'source_sha256': sha(source),
        'inputs': {path: expected for path, expected in FILES.values()},
        'rehash_all_final_readback_bindings': True,
        'stage_figure': 'Compute (eta+zr)/r from each reviewed rational result on the fixed original cube. Plot its logarithmically scaled position against threshold 1.',
        'force_figure': 'Compare the centered separate-norm eta and independent 256-bit signed eta with r(1-z). Convert exact rational values only for display.',
        'scope': 'Postprocessing of bounded proof quantities. No execution of a force integral, cache sum or equilibrium job.',
    }
    path = OUT / 'protocol.json'
    dump(path, protocol)
    plots = (stage_plot(stages), force_plot(stages, metrics))
    report = {
        'status': 'PASS_EXACT_L_FIGURE_READBACK', 'source_sha256': sha(source), 'protocol_sha256': sha(path),
        'input_sha256': protocol['inputs'], 'final_readback_bindings_rehashed': receipts['readback']['sha256'],
        'checks': checks, 'matplotlib': matplotlib.__version__,
        'stages': {name: {key: str(value) for key, value in stage.items()} for name, stage in stages.items()},
        'metrics': {key: str(value) for key, value in metrics.items()},
        'figures': {p.relative_to(ROOT).as_posix(): sha(p) for p in plots},
    }
    dump(OUT / 'receipt.json', report)
    print(json.dumps({'receipt_sha256': sha(OUT / 'receipt.json'), 'figures': report['figures']}))


if __name__ == '__main__':
    main()
