"""Read-only reduction of the completed transfer-pair pilot and retained profiles."""

from datetime import datetime
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


OUT = Path(__file__).resolve().parent
RUN = OUT.parent
ROOT = RUN.parents[1]
RESULT_SHA = '6a2fcb2a0c4b7a51301d025a7725247ec81a999764b7e85a0e97b5f6e139be05'
CASES = ('baseline_small', 'baseline_large', 'joint_small', 'joint_large')
MATERIALS = ('baseline', 'joint')
COLORS = {'baseline': '#70d8e6', 'joint': '#f0b67b', 'apical': '#7fdaee', 'basal': '#edae84'}
SURFACE_LABELS = {'apical': 'Apical', 'basal': 'Basal'}
MATERIAL_LABELS = {'baseline': 'Gain-only calibration', 'joint': 'Joint material/gain calibration'}
BACKGROUND = '#111921'
PANEL = '#18232e'
TEXT = '#eef5fb'
MUTED = '#a9bbc9'
GRID = '#3a4b58'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def rel(path):
    return path.relative_to(ROOT).as_posix()


def load_bound(path, expected):
    require(sha(path) == expected, f'hash changed: {path}')
    return read(path)


def strict_checks(checks):
    return bool(checks) and all(value is True for value in checks.values())


def frac_interval(serial):
    require(serial['lower_closed'] and serial['upper_closed'], 'closed interval required')
    return tuple(Q(serial[key]['numerator'], serial[key]['denominator']) for key in ('lower', 'upper'))


def interval_json(bounds):
    a, b = bounds
    require(a <= b, 'ordered bounds required')
    return {'exact': [str(a), str(b)], 'display_um': [float(a), float(b)],
            'width_um_display': float(b - a)}


def subtract(a, b):
    return a[0] - b[1], a[1] - b[0]


def midpoint(bounds):
    return float(sum(bounds) / 2)


def coefficients(y):
    begin, mid, end = y
    return begin, -3 * begin + 4 * mid - end, 2 * (begin - 2 * mid + end)


def polynomial(coef, t):
    a, b, c = coef
    return a + t * (b + t * c)


def surface_controls(field, height_scale):
    reference = field['reference_nodes']
    nodes = field['current_nodes']
    edge_incidence = {}
    for a, b, c, ab, bc, ca in field['triangles']:
        for i, j, midpoint_id in ((a, b, ab), (b, c, bc), (c, a, ca)):
            edge_incidence.setdefault(tuple(sorted((i, j))), []).append(midpoint_id)
    surfaces = {}
    for name, height, old_key in (('apical', field['parameters']['height'], 'top'), ('basal', 0, 'bottom')):
        edges = []
        for (i, j), midpoints in edge_incidence.items():
            if len(midpoints) != 1 or reference[i][1] != height or reference[j][1] != height:
                continue
            i, j = sorted((i, j), key=lambda n: reference[n][0])
            points = [nodes[n] for n in (i, midpoints[0], j)]
            edges.append((reference[i][0], points))
        pieces = [points for _, points in sorted(edges)]
        require(pieces == field[f'boundary_{old_key}_quadratic_controls'], 'raw connectivity/control mismatch')
        surfaces[name] = [[(height_scale * Q(r), height_scale * Q(z)) for r, z in piece] for piece in pieces]
    return surfaces


def exact_graph(controls):
    require(bool(controls), 'empty trace')
    require(all(a[-1] == b[0] for a, b in zip(controls, controls[1:])), 'trace endpoints do not join')
    segments = []
    for piece in controls:
        radius, height = tuple(coefficients(coord) for coord in zip(*piece))
        end = polynomial(radius, Q(1))
        derivative = min(radius[1], radius[1] + 2 * radius[2])
        require(derivative > 0, 'saved graph lacks a strictly positive radial derivative')
        segments.append((radius, end, height, derivative))
    return segments


def height_bounds(segments, radius, bisections):
    """Enclose z at a current radius by exact polynomial inversion on one segment."""
    for radial, end, height, _ in segments:
        if not radial[0] <= radius <= end:
            continue
        a, b = Q(0), Q(1)
        if radius == radial[0]:
            return height[0], height[0]
        if radius == end:
            z = polynomial(height, b)
            return z, z
        for _ in range(bisections):
            mid = (a + b) / 2
            value = polynomial(radial, mid)
            if value == radius:
                a = b = mid
                break
            if value < radius:
                a = mid
            else:
                b = mid
        z = [polynomial(height, a), polynomial(height, b)]
        if height[2]:
            critical = -height[1] / (2 * height[2])
            if a <= critical <= b:
                z.append(polynomial(height, critical))
        return min(z), max(z)
    raise ValueError('query radius outside the stored graph')


def check_samples(profile, graphs, shoulder, radius_limit, bisections):
    samples = profile['sampled_full_current_profiles']
    radii = [Q(row['radius_um_exact']) for row in samples]
    require(radii[0] == 0 and radii[-1] == radius_limit and shoulder in radii
            and all(a < b for a, b in zip(radii, radii[1:])), 'incomplete or unordered profile samples')
    for row, radius in zip(samples, radii):
        require(row['radius_um'] == float(radius), 'display radius mismatch')
        for name, graph in graphs.items():
            observed = frac_interval(row[f'{name}_height_um'])
            require(observed == height_bounds(graph, radius, bisections), 'retained profile sample differs from P2 graph')
    return samples


def verify(production, result, approved):
    master = read(RUN / 'job.json')
    require(master['results_sha256'] == RESULT_SHA and master['status'] == result['status']
            == 'COMPLETE_FOUR_CONDITIONAL_FIELDS', 'complete terminal supervisor result required')
    require(master['bindings'] == result['bindings'] and master['terminal_cases'] == result['terminal_case_receipts'],
            'supervisor and results inventories differ')
    bindings = result['bindings']
    review = load_bound(approved, bindings['prelaunch_review_sha256'])
    require(review['status'] == 'PASS' and strict_checks(review['checks']), 'independent prelaunch receipt does not pass')
    for key in ('driver_sha256', 'protocol_sha256', 'preflight_sha256'):
        require(review[key] == bindings[key], 'prelaunch binding mismatch')
    require(sha(ROOT / 'scripts/transfer_pair_pilot.py') == bindings['driver_sha256']
            and sha(RUN / 'preparation/preflight.json') == bindings['preflight_sha256'], 'driver/preflight bytes changed')
    for kind in ('source_hashes', 'input_hashes'):
        require(production[kind] == result[kind], f'{kind} differing inventories')
        for path, expected in production[kind].items():
            require(sha(ROOT / path) == expected, f'current {kind} mismatch: {path}')
    for flag in ('source_outcomes_read_or_scored', 'root_validated', 'source_compatible_witness', 'novelty_established'):
        require(result[flag] is False, f'unsupported promotion: {flag}')

    shoulder = Q(production['observable']['shoulder_radius_um_exact'])
    length = Q(production['observable']['length_unit_um_exact'])
    bisections = production['observable']['root_bisections']
    cases = {case['case_id']: case for case in production['cases']}
    jobs = {job['case_id']: job for job in master['terminal_cases']}
    require(tuple(cases) == tuple(jobs) == CASES and len(master['terminal_cases']) == len(CASES)
            and len({job['pid'] for job in jobs.values()}) == len(CASES), 'wrong terminal case/PID inventory')
    metadata = {}
    values = {}
    profiles = {}
    for case_id, case in cases.items():
        directory = RUN / 'jobs' / case_id
        job = jobs[case_id]
        require(read(directory / 'job.json') == job and job['bindings'] == bindings, 'terminal sidecar differs')
        require(job['started'] is True and job['terminal'] is True and job['status'] == 'PROCESS_EXITED'
                and job['exit_code'] == 0 and job['pid'] > 0, 'actual zero exit code required')
        worker_path = ROOT / job['worker_receipt']
        worker = load_bound(worker_path, job['worker_receipt_sha256'])
        require(worker['worker_pid'] == job['pid'] and worker['status'] == job['worker_status']
                == 'COMPLETE_CONDITIONAL_STORED_FIELD', 'OS/worker status or PID mismatch')
        require(worker['case'] == case and worker['bindings'] == bindings
                and worker['input_hashes'] == production['input_hashes']
                and worker['source_hashes'] == production['source_hashes'], 'case/result provenance mismatch')
        started, initialized, finished, exited = map(datetime.fromisoformat, (
            job['started_at_utc'], worker['started_at_utc'], worker['worker_finished_at_utc'], job['ended_at_utc']))
        require(started <= initialized <= finished <= exited and job['elapsed_seconds'] > 0,
                'inconsistent terminal event timestamps')
        require(len(worker['checkpoints']) == len(case['schedule']) == 4, 'missing checkpoint')
        retained = []
        for index, (report, step) in enumerate(zip(worker['checkpoints'], case['schedule'])):
            require(report['index'] == index and all(report[k] == value for k, value in step.items()),
                    'checkpoint schedule differs')
            state = load_bound(ROOT / report['field'], report['field_sha256'])
            checks = load_bound(ROOT / report['verification'], report['verification_sha256'])
            require(state['parameters'] == case['parameters'] and state['reached_target'] is True
                    and state['converged'] is report['converged'] is True
                    and Q.from_float(state['activation_strength']) == Q(step['solver_strength_binary64_exact']),
                    'unexpected field parameters, strength or convergence')
            require(checks['status'] == report['numerical_diagnostic_status'] == 'PASS'
                    and checks['checks'] == report['checks'] and strict_checks(checks['checks']), 'diagnostic gate failed')
            require(checks['root_validated'] is False and worker['root_validated'] is False,
                    'incorrect root promotion')
            q8 = checks['Q8_mechanical_evaluation']['scaled_force_linf']
            q12 = checks['Q12_cross_evaluation']['scaled_force_linf']
            tolerance = production['numerical_policy']['solve']['force_tolerance']
            require(q8 == report['Q8_residual'] and q12 == report['Q12_residual'] and 0 <= max(q8, q12) < tolerance,
                    'normalized force residual binding or threshold failed')
            geometry = checks['geometry']
            require(geometry['all_exact_J_positive'] and geometry['global_embedding_certified']
                    and Q(geometry['J_lower_exact']) > 0 and geometry['J_lower_exact'] == report['exact_J_lower']
                    and geometry['J_upper_exact'] == report['exact_J_upper'], 'stored-body geometry gate failed')
            require(strict_checks(checks['surface']['checks']) and checks['surface']['usable_stored_graph'],
                    'stored surface graph gate failed')
            retained.append({'index': index, 'field': report['field'], 'field_sha256': report['field_sha256'],
                             'verification': report['verification'], 'verification_sha256': report['verification_sha256'],
                             'Q8_normalized_force_linf': q8, 'Q12_normalized_force_linf': q12,
                             'J_exact': [geometry['J_lower_exact'], geometry['J_upper_exact']],
                             'attempts': report['attempts'], 'rejected_attempts': report['rejected_attempts']})
        report = worker['checkpoints'][-1]
        require(report == worker['last_retained_field'] == result['diagnostic_fields_accepted'][case_id],
                'final field inventory differs')
        geometry = checks['geometry']
        surface = checks['surface']
        profile = load_bound(ROOT / worker['profiles'], worker['profiles_sha256'])
        require(profile['case_id'] == case_id and profile['field_sha256'] == report['field_sha256']
                and profile['planned_final_checkpoint'] and profile['diagnostic_gate_passed'], 'profile binding differs')
        require(all(profile[key] == value for key, value in surface.items()), 'profile and final verification differ')
        controls = surface_controls(state, length)
        for name, pieces in controls.items():
            require(pieces == [[tuple(map(Q, point)) for point in piece]
                               for piece in profile['controls_um'][name]], 'physical trace controls differ')
        graphs = {name: exact_graph(pieces) for name, pieces in controls.items()}
        derivatives = {name: min(segment[3] for segment in graph) for name, graph in graphs.items()}
        require({name: str(value) for name, value in derivatives.items()}
                == profile['minimum_radius_derivative_um_per_local_parameter_exact'], 'radial derivative bound differs')
        radius_limit = Q(case['parameters']['outer_radius']) * length
        require(all(graph[0][0][0] == 0 and graph[-1][1] == radius_limit for graph in graphs.values()),
                'surface radial domain differs from model')
        samples = check_samples(profile, graphs, shoulder, radius_limit, bisections)
        graph = graphs['apical']
        center = height_bounds(graph, Q(0), bisections)
        shoulder_height = height_bounds(graph, shoulder, bisections)
        q = subtract(center, shoulder_height)
        require(center == frac_interval(profile['apical_center_um'])
                and shoulder_height == frac_interval(profile['apical_shoulder_um'])
                and q == frac_interval(profile['Q_um'])
                and [str(x) for x in q] == profile['Q_um_exact'] == report['Q_um_exact']
                == result['signed_Q_um'][case_id]['exact'], 'exact-center or Q replay mismatch')
        values[case_id] = q
        profiles[case_id] = samples
        metadata[case_id] = {'pid': job['pid'], 'process_observed_started': job['started'],
                             'process_observed_exited': job['terminal'], 'OS_exit_code': job['exit_code'],
                             'terminal_sidecar': rel(directory / 'job.json'), 'terminal_sidecar_sha256': sha(directory / 'job.json'),
                             'worker_receipt': rel(worker_path), 'worker_receipt_sha256': sha(worker_path),
                             'profiles': worker['profiles'], 'profiles_sha256': worker['profiles_sha256'],
                             'samples': len(samples), 'surface_edge_counts': {k: len(v) for k, v in controls.items()},
                             'radius_domain_um_exact': ['0', str(radius_limit)], 'plot_connects_retained_samples': True,
                             'minimum_radial_derivatives_exact': {k: str(v) for k, v in derivatives.items()},
                             'minimum_graph_thickness_um_lower_exact': profile['minimum_current_thickness_um_lower_exact'],
                             'exact_global_embedding': geometry['global_embedding_certified'],
                             'exact_boundary_status': geometry['boundary']['status'],
                             'checkpoints': retained, 'root_error_enclosure': None}
    return metadata, values, profiles, master


def contrast_arithmetic(values, result):
    effects = {name: subtract(values[f'{name}_large'], values[f'{name}_small']) for name in MATERIALS}
    gaps = {diameter: subtract(values[f'joint_{diameter}'], values[f'baseline_{diameter}'])
            for diameter in ('small', 'large')}
    difference_of_differences = subtract(effects['joint'], effects['baseline'])
    require(difference_of_differences == subtract(gaps['large'], gaps['small']), 'contrast identity mismatch')
    for name, value in effects.items():
        require([str(v) for v in value] == result['signed_diameter_change_large_minus_small_Q_um'][name]['exact'],
                'saved diameter effects differ')
    large_gap = gaps['large']
    pair = result['large_diameter_gap_test']
    require((large_gap[1] < 0 or large_gap[0] > 0)
            and pair['status'] == 'DISJOINT_STORED_CURVE_INTERVALS'
            and [str(v) for v in large_gap] == pair['joint_minus_baseline_large_Q_um']['exact'],
            'invalid archived large diameter separation')
    return effects, gaps, difference_of_differences


def style(ax):
    ax.set_facecolor(PANEL)
    ax.tick_params(colors=MUTED, length=0, pad=8, labelsize=9)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.grid(axis='y', color=GRID, linewidth=0.65, alpha=0.72)
    ax.set_axisbelow(True)
    ax.xaxis.label.set_color(MUTED)
    ax.yaxis.label.set_color(MUTED)


def profile_plot(profiles, metadata, production):
    fig, axes = plt.subplots(2, 2, figsize=(11.4, 7.9), sharex=True, sharey=True, facecolor=BACKGROUND)
    fig.subplots_adjust(left=0.10, right=0.967, bottom=0.16, top=0.78, wspace=0.16, hspace=0.36)
    selected = production['input_selection']
    dose = float(Q(selected['matched_dose_exact']))
    shoulder = float(Q(production['observable']['shoulder_radius_um_exact']))
    fig.text(0.10, 0.935, 'Final apical and basal profiles', color=TEXT, size=23, weight='bold')
    fig.text(0.10, 0.886, f'Shared optical dose {dose:.8f}  •  Displayed over the full stored radial domain', color=MUTED, size=12)
    legend = [Line2D([0], [0], color=COLORS[name], lw=2.7, label=label)
              for name, label in SURFACE_LABELS.items()]
    legend.append(Line2D([0], [0], color=MUTED, lw=1.1, ls='--', label=f'Shoulder radius {shoulder:.3f} µm'))
    fig.legend(handles=legend, loc='upper left', bbox_to_anchor=(0.09, 0.855), ncols=3,
               frameon=False, labelcolor=TEXT, fontsize=10)
    lookup = {case['case_id']: case for case in production['cases']}
    for ax, case_id in zip(axes.flat, CASES):
        style(ax)
        case = lookup[case_id]
        samples = profiles[case_id]
        for name in ('apical', 'basal'):
            heights = [midpoint(frac_interval(row[f'{name}_height_um'])) for row in samples]
            ax.plot([row['radius_um'] for row in samples], heights, color=COLORS[name], lw=2.45)
        ax.axvline(shoulder, color=MUTED, linestyle='--', lw=1.0, alpha=0.8)
        ax.set_xlim(0, 100.5)
        ax.set_ylim(-1.5, 36)
        ax.set_xticks([0, 25, 50, 75, 100])
        diameter = float(Q(case['diameter_um_exact']))
        label = MATERIAL_LABELS[case['candidate']]
        ax.set_title(f'{label}  ·  {diameter:.2f} µm patch', loc='left', fontsize=11, color=TEXT, weight='bold', pad=11)
        row = metadata[case_id]['checkpoints'][-1]
        ax.text(0.98, 0.50, f"Force ∞-norm   Q8 {row['Q8_normalized_force_linf']:.2e}  ·  Q12 {row['Q12_normalized_force_linf']:.2e}",
                transform=ax.transAxes, ha='right', color=MUTED, size=8.7)
    fig.supxlabel('Current radius R (µm)', color=TEXT, fontsize=11, y=0.095)
    fig.supylabel('Height z (µm)', color=TEXT, fontsize=11, x=0.026)
    fig.text(0.10, 0.057, 'All four fields pass stored-field Jacobian, embedding, and radial graph checks.', color=MUTED, fontsize=9.4)
    fig.text(0.10, 0.028, 'Lines connect retained samples. Full-profile supremum and equilibrium-position error remain unresolved.',
             color=MUTED, fontsize=9.4)
    path = OUT / 'current_radius_profiles.png'
    fig.savefig(path, dpi=162, facecolor=BACKGROUND)
    plt.close(fig)
    return path


def contrast_plot(values, effects, gaps, difference_of_differences, production):
    fig, (ax, effect_ax) = plt.subplots(1, 2, figsize=(11.4, 5.6), facecolor=BACKGROUND,
                                       gridspec_kw={'width_ratios': [1.15, 1]})
    fig.subplots_adjust(left=0.10, right=0.96, bottom=0.23, top=0.69, wspace=0.28)
    fig.text(0.10, 0.925, 'Center-to-shoulder shape contrast', color=TEXT, fontsize=23, weight='bold')
    fig.text(0.10, 0.858, 'Q = apical height at the center minus apical height at the fixed shoulder', color=MUTED, fontsize=11.5)
    legend = [Line2D([0], [0], color=COLORS[name], marker='o', lw=2.5, label=MATERIAL_LABELS[name]) for name in MATERIALS]
    fig.legend(handles=legend, loc='upper left', bbox_to_anchor=(0.092, 0.812), ncols=2,
              frameon=False, labelcolor=TEXT, fontsize=10)

    style(ax)
    style(effect_ax)
    ax.set_title('Signed Q at matched dose', color=TEXT, fontsize=12, loc='left', pad=14, weight='bold')
    ax.set_xlim(-0.26, 1.32)
    ax.set_ylim(-8.1, -3.50)
    selection = production['input_selection']
    diameters = [float(Q(selection[key])) for key in ('matched_small_reported_diameter_um_exact',
                                                      'large_reported_diameter_um_exact')]
    ax.set_xticks([0, 1], [f'{d:.2f} µm' for d in diameters])
    ax.set_ylabel('Q (µm)', labelpad=9)
    ax.set_xlabel('Reported patch diameter', labelpad=9)
    for name in MATERIALS:
        ys = [midpoint(values[f'{name}_{level}']) for level in ('small', 'large')]
        color = COLORS[name]
        other = 'joint' if name == 'baseline' else 'baseline'
        ax.plot([0, 1], ys, marker='o', ms=6, lw=2.2, color=color)
        for x, (level, y) in enumerate(zip(('small', 'large'), ys)):
            vertical = 8 if y > midpoint(values[f'{other}_{level}']) else -13
            ax.annotate(f'{y:.4f}', (x, y), xytext=(8 if x else -8, vertical), textcoords='offset points',
                        ha='left' if x else 'right', color=color, size=9)

    effect_ax.set_title('Diameter effect on shape', color=TEXT, fontsize=12, loc='left', pad=14, weight='bold')
    effect_ax.set_ylim(0, 3.74)
    effect_ax.set_xlim(-0.60, 1.60)
    effect_ax.set_xticks([0, 1], ['Gain-only', 'Joint'])
    effect_ax.set_ylabel('Large Q minus small Q (µm)', labelpad=8)
    for i, name in enumerate(MATERIALS):
        val = midpoint(effects[name])
        effect_ax.bar(i, val, width=0.51, color=COLORS[name], alpha=0.77)
        effect_ax.text(i, val + 0.09, f'+{val:.4f}', ha='center', color=COLORS[name], fontsize=10, weight='bold')
    fig.text(0.10, 0.109, f'Joint minus gain-only at the larger diameter  {midpoint(gaps["large"]):+.5f} µm',
             color=TEXT, fontsize=11)
    fig.text(0.10, 0.064, f'Joint minus gain-only diameter effect  {midpoint(difference_of_differences):+.5f} µm',
             color=TEXT, fontsize=11)
    fig.text(0.10, 0.023, 'Stored-curve interval widths are below display precision. True-equilibrium error is unresolved.',
             color=MUTED, fontsize=9.4)
    path = OUT / 'shape_contrast_and_diameter_effect.png'
    fig.savefig(path, dpi=162, facecolor=BACKGROUND)
    plt.close(fig)
    return path


def main():
    result = load_bound(RUN / 'results.json', RESULT_SHA)
    production_path = RUN / 'preparation/protocol.json'
    production = load_bound(production_path, result['bindings']['protocol_sha256'])
    approved = RUN / 'prelaunch_independent_review.json'
    protocol = {'result_sha256': RESULT_SHA, 'analyzer_sha256': sha(Path(__file__)),
                'producer_protocol_sha256': sha(production_path), 'prelaunch_sha256': sha(approved),
                'provenance': 'Frozen inputs, source hashes, OS exit records, and all 16 archived checkpoint file pairs.',
                'calculation': 'Exact rational height inversion reimplemented from P2 nodal connectivity. Signed Q = z_apical(0) - z_apical(shoulder). Derived intervals use endpoint subtraction.',
                'profile_plot': 'Plot the midpoint of each retained interval at every stored current-radius sample over the full domain. Independently verify each sample against raw P2 controls. The connecting segments provide visual interpolation.',
                'diagnostic_scope': 'The predeclared normalized force and exact stored-geometry/graph checks are separate from equilibrium-position accuracy. Equilibrium-position error, full-profile supremum, integration error and continuum error are unresolved.',
                'heldout_outcomes_read': False}
    metadata, values, profiles, master = verify(production, result, approved)
    effects, gaps, dd = contrast_arithmetic(values, result)
    saved_memo = 'runs/exact_integration_cohort_v1/scientific_next_action/receipt.json'
    saved_pair = load_bound(ROOT / saved_memo, production['input_hashes'][saved_memo])['screened_candidate_pair']
    calibration = {name: max(Q(row['sufficient_upper_um_exact']) for row in saved_pair[key].values())
                   for name, key in (('baseline', 'baseline_stored_curve_requirements'),
                                     ('joint', 'joint_stored_curve_requirements'))}
    figures = [profile_plot(profiles, metadata, production), contrast_plot(values, effects, gaps, dd, production)]
    save(OUT / 'protocol.json', protocol)
    save(OUT / 'receipt.json', {
        'status': 'COMPLETE_CONDITIONAL_STORED_CURVE_POSTPROCESSING',
        'source_sha256': sha(Path(__file__)), 'protocol_sha256': sha(OUT / 'protocol.json'),
        'results_sha256': RESULT_SHA, 'supervisor_receipt_sha256': sha(RUN / 'job.json'), 'supervisor_pid': master['pid'],
        'checked_current_source_hashes': production['source_hashes'], 'checked_current_input_hashes': production['input_hashes'],
        'cases': metadata, 'shoulder_radius_um_exact': production['observable']['shoulder_radius_um_exact'],
        'signed_Q_um': {k: interval_json(v) for k, v in values.items()},
        'diameter_effect_large_minus_small_Q_um': {k: interval_json(v) for k, v in effects.items()},
        'joint_minus_baseline_Q_um': {k: interval_json(v) for k, v in gaps.items()},
        'difference_of_differences_joint_minus_baseline_effect_um': interval_json(dd),
        'large_pair_disjoint_stored_curve_intervals': True,
        'conditional_68_row_calibration_residual_context': {
            'rows_per_candidate': 68, 'bound_source': saved_memo,
            'largest_sufficient_stored_curve_discrepancy_requirement_each_candidate': {
                key: {'exact_um': str(value), 'display_um': float(value)} for key, value in calibration.items()},
            'measurement_noise_estimate': False,
            'scope': 'These residual requirements use a different static endpoint observable and supply context only. They do not calibrate the uncertainty of this shape contrast.'},
        'normalized_force_limit': production['numerical_policy']['solve']['force_tolerance'],
        'root_validated': False, 'root_position_error_enclosure_um': None,
        'full_profile_supremum_enclosure_um': None, 'source_compatible_witness': False,
        'larger_source_outcomes_read_or_scored': False, 'novelty_established': False,
        'visualizations': {p.name: {'path': rel(p), 'sha256': sha(p)} for p in figures},
    })
    print(json.dumps({'status': 'COMPLETE', 'receipt_sha256': sha(OUT / 'receipt.json'),
                      'plots': {path.name: sha(path) for path in figures}}))


if __name__ == '__main__':
    main()
