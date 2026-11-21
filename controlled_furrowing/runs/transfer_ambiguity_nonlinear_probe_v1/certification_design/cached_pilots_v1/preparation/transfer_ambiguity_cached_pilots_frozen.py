"""Prepare four pointwise cached-quadrature pilots at frozen exact optical drive."""
import argparse
import hashlib
import importlib.metadata
import json
import os
from fractions import Fraction as Q
from pathlib import Path
import subprocess
import sys
import time
import traceback
from unittest.mock import patch

import transfer_ambiguity_error_budgets as budgets
import transfer_ambiguity_observer_sweep as sweep


SOURCE = Path(__file__).resolve()
ROOT = SOURCE.parents[1]
BASE = f'{sweep.PRIMARY}/certification_design'
BUDGET_RUN = f'{BASE}/budgets_v1'
OUT = ROOT / BASE / 'cached_pilots_v1'
PREP = OUT / 'preparation'
ACTUAL = OUT / 'actual'
BUDGET_RESULT = f'{BUDGET_RUN}/actual/results.json'
BUDGET_TERMINAL = f'{BUDGET_RUN}/actual/terminal.json'
BUDGET_REVIEW = f'{BUDGET_RUN}/independent_actual/receipt.json'
CORE = 'scripts/algebraic_root_feasibility.py'
CORE_REVIEW = 'evidence/algebraic_root_independent_review.json'
BUDGET_SOURCE = 'scripts/transfer_ambiguity_error_budgets.py'
SWEEP_SOURCE = 'scripts/transfer_ambiguity_observer_sweep.py'
PINS = {
    BUDGET_RESULT: '689948766db950ee4ddf399e6731aea08890cc7a1f84c5615d19e7ded68b0872',
    BUDGET_TERMINAL: '83abadac599cba377034a522aa5c160de176b4b8b6b91f6223883415495d9928',
    BUDGET_REVIEW: '978a2a4095bb8735103a5ee2c11ad266a6b0e052e22af2b7a36ef29f1f38265c',
    CORE: '7321d38bf17cf354cc03be2871f3fbaa45bcebbcbfca7be90d0f867b80078acf',
    CORE_REVIEW: 'f327f2fab76f9fdcab3a83882adbb1a4ae3f0cf1e5c269e7335fb909691085d5',
    BUDGET_SOURCE: '410a4944798b1f42f0360cbdf236d9da68f8de2615d93f860bccdd8be54d8859',
    SWEEP_SOURCE: '9eacb1de680cefb84d90d603cc1aca4dd6ef4892372a91cffa3480651c5874a8',
    **budgets.ACTUAL_PINS,
}
PRECISION = 256
RADII = tuple(Q(1, 2**k) for k in (20, 24, 28, 32, 36))
NFREE = 611
TERMS = {'bulk': 144 * 64, 'surface': 24 * 12}
SCOPE = ('Pointwise local-root feasibility for the exact binary-cache quadrature potential at one exact requested activation strength '
         'per saved field. The derivative bound covers the whole requested free-coordinate cube at cached quadrature sites. '
         'Exact spatial integration, whole-body orientation and embedding, energy minimum and PDE implications require separate proofs.')
require, sha, pinned, read, save = sweep.require, sweep.sha, sweep.pinned, sweep.read, sweep.save


def exact_drive(step, material):
    gain, dose = Q(material['gain_exact']), Q(step['dose_exact'])
    strength = gain * dose
    rounded = Q.from_float(float(strength))
    require(gain > 0 and dose > 0 and strength == Q(step['requested_strength_exact'])
            and rounded == Q(step['solver_strength_binary64_exact'])
            and rounded - strength == Q(step['binary64_rounding_difference_exact']),
            'requested strength does not match exact gain times dose and its single binary64 conversion')
    return {'optical_dose_exact': str(dose), 'common_gain_exact': str(gain),
            'activation_strength_exact': str(strength),
            'activation_strength_binary_exact': str(rounded),
            'binary_input_minus_exact_strength': str(rounded - strength)}


def admitted_targets():
    for name, digest in PINS.items():
        require(sha(ROOT / name) == digest, f'frozen target or source changed: {name}')
    for name, module in ((BUDGET_SOURCE, budgets), (SWEEP_SOURCE, sweep)):
        require(Path(module.__file__).resolve() == ROOT / name, 'bound target helper module alias drift')
    old = pinned(ROOT / CORE_REVIEW, PINS[CORE_REVIEW])
    require(old['status'] == 'PASS' and sweep.all_true(old['checks']) and old['source_sha256'] == PINS[CORE],
            'prior review of reusable cached-root algebra missing')
    primary = pinned(ROOT / sweep.PRIMARY_PROT, sweep.PINS[sweep.PRIMARY_PROT])
    sources = {**sweep.PINS, **primary['source_hashes']}
    for name, digest in sources.items():
        require(sha(ROOT / name) == digest, f'frozen model or prerequisite source changed: {name}')
    protocol = {'frozen_44_checkpoint_plan': sweep.plan_checkpoints(primary)}
    comparisons, native_review, selected = budgets.admit(
        protocol, PINS[budgets.NATIVE_REVIEW], PINS[budgets.MANIFEST])
    require(comparisons['heldout_response_tables_read'] is False and len(selected) == 36,
            'wrong primary stored-field authority')

    outcome = pinned(ROOT / BUDGET_RESULT, PINS[BUDGET_RESULT])
    bterminal = pinned(ROOT / BUDGET_TERMINAL, PINS[BUDGET_TERMINAL])
    review = pinned(ROOT / BUDGET_REVIEW, PINS[BUDGET_REVIEW])
    observed = review['observed_budget_terminal']
    require(bterminal['status'] == 'FINISHED' and bterminal['results_sha256'] == PINS[BUDGET_RESULT]
            and bterminal['bindings'] == outcome['bindings'] and observed['pid'] == bterminal['pid']
            and type(observed['exit_code']) is int and observed['exit_code'] == 0,
            'independently observed budget terminal missing')
    require(review['status'] == 'PASS_INDEPENDENT_ACTUAL_CONDITIONAL_OBSERVER_BUDGETS'
            and sweep.all_true(review['all_literal_checks'])
            and review['root_error_bounds_proved'] == outcome['exact_root_error_bounds_proved'] == 0,
            'observer budget prerequisite is incomplete')
    for name in (BUDGET_RESULT, BUDGET_TERMINAL, budgets.NATIVE_REVIEW, budgets.MANIFEST,
                 sweep.PRIMARY_RESULT, sweep.PRIMARY_TERMINAL, BUDGET_SOURCE):
        require(review['source_and_input_hashes'][name] == PINS[name], f'mixed budget review authority: {name}')
    require(outcome['status'] == 'COMPLETE_POSITIVE_CONDITIONAL_36_FIELD_BUDGETS'
            and len(outcome['field_budgets']) == 36 and len(outcome['four_pilot_selection']) == 4,
            'exactly four positive-budget pilots required')

    case_by_id = {x['endpoint_id']: x for x in primary['endpoints']}
    full = {(x['endpoint_id'], x['checkpoint']): x for x in outcome['field_budgets']}
    manifest = {(x['endpoint_id'], x['checkpoint']): x for x in selected}
    require(len(full) == len(manifest) == 36 and set(full) == set(manifest), 'mixed budget field coverage')
    expected = [(eid, stage) for eid in sweep.ENDPOINTS for stage in ('small_17', 'large_03')]
    require([(x['endpoint_id'], x['checkpoint']) for x in outcome['four_pilot_selection']] == expected,
            'four pilots drifted from independently selected two L fields and two large final fields')
    pilots = []
    for reduced, independently_replayed, identity in zip(
            outcome['four_pilot_selection'], review['four_pilot_fields'], expected, strict=True):
        field = full[identity]
        observed_field = manifest[identity]
        goal = Q(reduced['e_goal_dimensionless_exact'])
        shared = ('endpoint_id', 'checkpoint', 'record_id', 'field', 'field_sha256', 'e_goal_dimensionless_exact')
        require(all(field[k] == reduced[k] == independently_replayed[k] for k in shared)
                and goal > 0 and any(r < goal for r in RADII),
                'four reviewed pilot keys or positive exact radius budget do not match')
        require(all(field[k] == observed_field[k] for k in observed_field)
                and native_review['all_hashes'][field['field']] == field['field_sha256']
                and native_review['all_hashes'][field['verification']] == field['verification_sha256'],
                'pilot field or proof differs from whole primary manifest')
        small = field['role'] == 'small_record'
        require((small and field['record_id'] == 'optoGEF2_small_col_L')
                or (field['role'] == 'large_final' and field['record_id'] == primary['matched_input']['large_record_id']),
                'unexpected curated record identity')
        case = case_by_id[identity[0]]
        step = field['step']
        material = case['materials']['small' if small else 'large']
        drive = exact_drive(step, material)
        center = field['center_and_model']
        require(len(center['free_map']) == len(set(center['free_map'])) == NFREE
                and len(center['active_elements']) == 24
                and drive['activation_strength_exact'] == center['exact_target_strength']
                and drive['activation_strength_binary_exact'] == center['solver_binary64_strength_exact']
                and drive['binary_input_minus_exact_strength'] == center['binary64_strength_rounding_error_exact']
                and Q(center['node_addition_gap_dimensionless_exact']) == Q(independently_replayed['binary_gap_exact'])
                and material['parameters'] == center['material_parameters']
                and material['parameter_storage_exact'] == center['material_parameter_storage_exact'],
                'exact drive, binary geometry or material differs from reviewed target')
        pilots.append({**reduced, 'pilot_id': f'{identity[0]}_{identity[1]}', 'step': step,
                       'material': material, 'exact_drive': drive, 'free_map': center['free_map'],
                       'active_elements': center['active_elements'],
                       'xi_dimensionless_exact': center['node_addition_gap_dimensionless_exact'],
                       'verification': field['verification'], 'verification_sha256': field['verification_sha256']})
    return primary, sources, pilots


def prepare_protocol():
    primary, sources, pilots = admitted_targets()
    packages = {name: importlib.metadata.version(name) for name in ('python-flint', 'numpy', 'scipy', 'torch')}
    require(all(packages[name] == primary['runtime'][name] for name in ('numpy', 'scipy', 'torch')),
            'binary-model numerical package versions differ from the frozen primary run')
    return {'status': 'FROZEN_FOUR_CACHED_QUADRATURE_PILOTS', 'scope': SCOPE,
            'source_sha256': sha(SOURCE), 'input_hashes': PINS,
            'source_and_prerequisite_hashes': sources, 'precision_bits': PRECISION, 'torch_cpu_threads': 1,
            'ordered_pilots': pilots, 'numerical_policy': primary['numerical_policy'],
            'runtime_packages': packages,
            'expected_site_counts': TERMS, 'descending_dyadic_radii_exact': [str(r) for r in RADII],
            'preconditioner': 'Fresh sparse Q8 Hessian at saved binary vector and nearest binary64 exact target drive. '
                              'Sparse LU inverse proposal, then all rigorous products use its exact stored binary entries.',
            'proof_target': 'Exact X+injection(u), exact binary shape/weights/normalization/parameters, exact rational gain*dose.',
            'strict_decision': 'r < e_goal, z < 1, eta + z*r < r, and eta/(1-z) < e_goal, with positive cached analytic domain.',
            'local_root_argument': 'Uniform z < 1 implies invertibility of R H and R, and the strict inward bound gives '
                                   'a unique local fixed point whose force is zero for this binary quadrature potential.',
            'execution': 'One fresh worker process per pilot, original order, serial true OS exit validation.',
            'rounded_solver_input_used_as_exact_proof_drive': False, 'new_nonlinear_equilibrium_solves': 0,
            'archived_fixture_cache_or_preconditioner_reused': False}


def upper_rational(core, ball):
    return Q(core.pack(ball)['upper_exact'])


def guard_radius(radius, goal):
    return 0 < radius < goal


def contraction(eta, z, radius, goal):
    require(eta >= 0 and z >= 0 and goal > 0 and radius > 0, 'invalid nonnegative strict decision data')
    lhs = eta + z * radius
    posterior = eta / (1 - z) if z < 1 else None
    strict = bool(guard_radius(radius, goal) and z < 1 and lhs < radius and posterior < goal)
    return {'eta_upper_exact': str(eta), 'z_upper_exact': str(z),
            'eta_plus_zr_upper_exact': str(lhs),
            'strict_selfmap_gap_lower_exact': str(radius - lhs),
            'posterior_error_upper_exact': None if posterior is None else str(posterior),
            'strict_selfmap_contraction_and_observer_radius': strict}


def point_assembly(core, model, u, requested_exact, n, expected_counts):
    require(type(requested_exact) is Q and requested_exact > 0, 'proof drive must be a positive exact Fraction')
    counts = dict.fromkeys(expected_counts, 0)
    def terms():
        for term in core.local_terms(model, u, requested_exact):
            counts[term['kind']] += 1
            yield term
    assembled = core.assemble_point(terms(), n)
    require(counts == expected_counts and len(assembled[2]) == sum(expected_counts.values()),
            'local cached site coverage incomplete or constrained contribution omitted unexpectedly')
    return assembled, counts


def model_input(pilot, policy):
    import numpy as np
    from transfer_pair_pilot import AxisymmetricModel, IntegrationPolicy, Parameters, SolvePolicy, SparsePolicy
    state = pinned(ROOT / pilot['field'], pilot['field_sha256'])
    # The independent manifest also authenticated this sibling; keep the actual worker bound to both.
    require(sha(ROOT / pilot['verification']) == pilot['verification_sha256'], 'stored field certification changed')
    material = pilot['material']
    params = Parameters(**material['parameters'])
    model = AxisymmetricModel(params, SolvePolicy(**policy['solve']),
                              IntegrationPolicy(**policy['solve_quadrature']), SparsePolicy(**policy['sparse']))
    _, geometry = budgets.model_receipt(model, state, pilot['step'], material, policy)
    require(geometry['free_map'] == pilot['free_map'] and len(model.free) == NFREE
            and geometry['active_elements'] == pilot['active_elements']
            and geometry['node_addition_gap_dimensionless_exact'] == pilot['xi_dimensionless_exact']
            and model.energy_scale == params.mu * model.reference_volume
            and model.reference_volume == float(model.integration.numpy().sum())
            and np.isfinite(model.energy_scale) and model.energy_scale > 0,
            'saved model maps, stored coordinates or prescribed quadrature normalization differ')
    u = np.asarray(state['free_displacement'], dtype='<f8')
    require(u.shape == (NFREE,) and np.isfinite(u).all(), 'finite 611 saved coordinates required')
    return model, u


def sparse_proposal(model, u, drive, folder):
    import numpy as np
    from scipy.sparse.linalg import splu
    proposal_drive = float(Q(drive['activation_strength_exact']))
    require(Q.from_float(proposal_drive) == Q(drive['activation_strength_binary_exact']),
            'numerical Hessian proposal strength binary64 relationship drift')
    H, diagnostics = model.recover_hessian(u, proposal_drive)
    require(H.shape == (NFREE, NFREE) and np.isfinite(H.data).all()
            and np.isfinite(diagnostics['gradient']).all(), 'nonfinite or mismatched fresh sparse Hessian')
    path = folder / 'fresh_sparse_hessian.npz'
    with path.open('xb') as handle:
        np.savez_compressed(handle, H_data=H.data, H_indices=H.indices, H_indptr=H.indptr,
                            H_shape=np.array(H.shape), free_indices=model.free, free_vector=u,
                            gradient=diagnostics['gradient'])
    inverse = splu(H.tocsc()).solve(np.eye(NFREE))
    require(inverse.shape == (NFREE, NFREE) and np.isfinite(inverse).all(), 'nonfinite inverse proposal')
    return inverse, {'hessian_archive': path.relative_to(ROOT).as_posix(), 'hessian_archive_sha256': sha(path),
                     'sparse_pattern': diagnostics['pattern'],
                     'matrix_nonzero_count': H.nnz, 'floating_raw_asymmetry_linf': diagnostics['raw_asymmetry_linf'],
                     'numeric_proposal_only': True}


def preconditioned_point(core, g, H, inverse):
    n = len(inverse)
    R = core.arb_mat(n, n, [core.A(float(x)) for x in inverse.ravel()])
    abs_R = core.arb_mat(n, n, [core.A(abs(float(x))) for x in inverse.ravel()])
    I = core.arb_mat(n, n)
    for i in range(n):
        I[i, i] = 1
    defect = I - R * H
    corrected = R * g
    z0 = max(upper_rational(core, sum((core.upper(defect[i, j]) for j in range(n)), core.arb(0)))
             for i in range(n))
    eta = max(upper_rational(core, core.upper(corrected[i, 0])) for i in range(n))
    rnorm = max(upper_rational(core, sum((abs_R[i, j] for j in range(n)), core.arb(0))) for i in range(n))
    rows = {'preconditioned_center_gradient_by_row': [core.pack(corrected[i, 0]) for i in range(n)],
            'preconditioner_infinity_norm_upper_exact': str(rnorm),
            'eta_upper_exact': str(eta), 'center_Z0_upper_exact': str(z0)}
    return eta, z0, abs_R, rows


def test_radius(core, cache, abs_R, eta, z0, radius, goal):
    row = {'radius_requested_exact': str(radius)}
    if not guard_radius(radius, goal):
        return {**row, 'status': 'REJECTED_RADIUS_OUTSIDE_OBSERVER_BUDGET',
                'strict_selfmap_contraction_and_observer_radius': False}
    try:
        variation, min_j, domain = core.variation_rows(cache, NFREE, core.A(radius))
    except ValueError as error:
        return {**row, 'status': 'UNRESOLVED_QUADRATURE_ANALYTIC_DOMAIN', 'reason': str(error),
                'strict_selfmap_contraction_and_observer_radius': False}
    require(domain['term_counts'] == TERMS and set(domain['strict_lower_margins']) == {
        'bulk_det_min', 'bulk_hoop_min', 'surface_radius_min', 'surface_tangent_squared_min',
        'normalized_weight_min'}, 'whole-ball quadrature domain receipt incomplete')
    if (Q(core.pack(min_j)['lower_exact']) <= Q(1, 1000000)
            or any(Q(x['lower_exact']) <= 0 for x in domain['strict_lower_margins'].values())):
        return {**row, 'status': 'UNRESOLVED_QUADRATURE_ANALYTIC_DOMAIN',
                'quad_min_J_enclosure': core.pack(min_j), 'quadrature_analytic_domain': domain,
                'strict_selfmap_contraction_and_observer_radius': False}
    propagated = abs_R * variation
    v = max(upper_rational(core, core.upper(propagated[i, 0])) for i in range(NFREE))
    bounds = contraction(eta, z0 + v, radius, goal)
    return {**row, **bounds, 'preconditioned_variation_upper_exact': str(v),
            'quad_min_J_enclosure': core.pack(min_j), 'quadrature_analytic_domain': domain,
            'status': 'STRICT_CACHED_QUADRATURE_ROOT_CANDIDATE'
            if bounds['strict_selfmap_contraction_and_observer_radius'] else 'UNRESOLVED_STRICT_INEQUALITY'}


def calculate(pilot, protocol, folder, bindings):
    import algebraic_root_feasibility as core
    import numpy as np
    import torch
    require(Path(core.__file__).resolve() == ROOT / CORE and sha(ROOT / CORE) == PINS[CORE],
            'resolved cached algebraic kernel differs from frozen review')
    torch.set_num_threads(1)
    core.ctx.prec = PRECISION
    core.A.cache_clear()
    require(core.ctx.prec == PRECISION and torch.get_num_threads() == 1, 'worker precision or thread limit drift')
    model, u = model_input(pilot, protocol['numerical_policy'])
    for name in (CORE, *protocol['source_and_prerequisite_hashes']):
        if not name.startswith('scripts/'):
            continue
        module = sys.modules.get(Path(name).stem)
        if module is not None:
            require(Path(module.__file__).resolve() == ROOT / name, f'runtime module alias drift: {name}')
    drive = pilot['exact_drive']
    require(drive == exact_drive(pilot['step'], pilot['material']), 'coherent drive reconstruction failed')
    strength = Q(drive['activation_strength_exact'])
    inverse, float_meta = sparse_proposal(model, u, drive, folder)
    packet = {'status': 'FROZEN_FRESH_EXACT_BINARY_QUADRATURE_CACHE', 'scope': SCOPE,
              'pilot_id': pilot['pilot_id'], 'record_id': pilot['record_id'],
              'input_state': pilot['field'], 'input_state_sha256': pilot['field_sha256'],
              'field_verification_sha256': pilot['verification_sha256'], 'parameters': pilot['material']['parameters'],
              'material_parameter_storage_exact': pilot['material']['parameter_storage_exact'],
              'integration_policy': protocol['numerical_policy']['solve_quadrature'], 'free_dof': NFREE,
              'normalization_binary_exact': str(Q.from_float(float(model.energy_scale))),
              'precision_bits': PRECISION, **drive, **float_meta,
              'drive_consumer_contract': 'Read activation_strength_exact as Fraction for analytical proof. '
                                         'Historical fixture replay expects rounded strength and is outside this packet.',
              'frozen_binary_cache': core.freeze_cache(model, u, inverse), 'bindings': bindings}
    packet_path = folder / 'cache_packet.json'
    packet_sha = save(packet_path, packet)
    (g, H, cache, min_j), counts = point_assembly(core, model, u, strength, NFREE, TERMS)
    require(Q(core.pack(min_j)['lower_exact']) > Q(1, 1000000),
            'exact point cached quadrature analytic threshold failed')
    eta, z0, abs_R, row_data = preconditioned_point(core, g, H, inverse)
    center = {**row_data, 'expected_and_rebuilt_site_counts': counts, 'point_quadrature_J_min': core.pack(min_j)}
    center_sha = save(folder / 'point_enclosure.json', center)
    goal = Q(pilot['e_goal_dimensionless_exact'])
    attempts = []
    chosen = None
    for index, radius in enumerate(RADII):
        row = test_radius(core, cache, abs_R, eta, z0, radius, goal)
        packet_name = f'radius_{index:02d}.json'
        artifact_sha = save(folder / packet_name, row)
        attempts.append({'path': packet_name, 'sha256': artifact_sha, **row})
        print(json.dumps({'pilot': pilot['pilot_id'], 'radius': str(radius), 'status': row['status']}), flush=True)
        if row['strict_selfmap_contraction_and_observer_radius']:
            chosen = row
            break
    return {'status': 'CANDIDATE_CERTIFIED_CACHED_QUADRATURE_LOCAL_ROOT' if chosen else
            'UNRESOLVED_CACHED_CONTRACTION_OR_BUDGET', 'scope': SCOPE,
            'pilot': pilot, 'bindings': bindings, 'backend': {
                'arb_precision_bits': int(core.ctx.prec), 'python_flint': importlib.metadata.version('python-flint'),
                'numpy': np.__version__, 'torch_threads': torch.get_num_threads()},
            'cache_packet': packet_path.relative_to(ROOT).as_posix(), 'cache_packet_sha256': packet_sha,
            'embedded_binary_cache_sha256': packet['frozen_binary_cache']['sha256'],
            'point_enclosure_sha256': center_sha, 'eta_upper_exact': str(eta), 'center_Z0_upper_exact': str(z0),
            'radius_attempts': attempts, 'accepted_radius': chosen, 'new_nonlinear_equilibrium_solves': 0,
            'exact_spatial_integration_established': False, 'whole_body_geometry_established': False,
            'local_minimum_established': False, 'independent_actual_result_review_required': True}


def prelaunch(review_path, review_sha):
    expected = {'source_sha256': sha(SOURCE), 'protocol_sha256': sha(PREP / 'protocol.json'),
                'preflight_sha256': sha(PREP / 'preflight.json'), 'budget_actual_sha256': PINS[BUDGET_RESULT],
                'budget_independent_actual_review_sha256': PINS[BUDGET_REVIEW],
                'whole_44_field_manifest_sha256': PINS[budgets.MANIFEST]}
    review = pinned(review_path, review_sha)
    require(sweep.good_review(review) and all(review.get(k) == v for k, v in expected.items()),
            'independent cached-pilot prelaunch receipt does not bind this exact launch')
    receipt = pinned(PREP / 'preflight.json', expected['preflight_sha256'])
    require(receipt['status'] == 'PASS' and sweep.all_true(receipt['checks'])
            and receipt['source_sha256'] == expected['source_sha256']
            and receipt['protocol_sha256'] == expected['protocol_sha256']
            and receipt['actual_preconditioners_recovered'] == receipt['actual_cached_sites_assembled'] == 0
            and sha(PREP / 'transfer_ambiguity_cached_pilots_frozen.py') == expected['source_sha256'],
            'current precision, source or preflight snapshot drift')
    protocol = pinned(PREP / 'protocol.json', expected['protocol_sha256'])
    require(protocol == prepare_protocol(), 'past review does not match the current frozen model/pilot plan')
    return protocol, {**expected, 'prelaunch_review_sha256': review_sha}


def worker(review_path, review_sha, pilot_id, parent_pid):
    protocol, bindings = prelaunch(review_path, review_sha)
    launch = read(ACTUAL / 'started.json')
    require(os.getppid() == parent_pid == launch['pid'] and launch['bindings'] == bindings
            and pilot_id in launch['ordered_pilot_ids'], 'exact authorized supervisor and target required')
    pilot, = [p for p in protocol['ordered_pilots'] if p['pilot_id'] == pilot_id]
    path = ACTUAL / pilot_id
    require(path.is_dir(), 'worker output directory must be created by supervisor')
    root = {'pid': os.getpid(), 'supervisor_pid': parent_pid, 'pilot_id': pilot_id, 'bindings': bindings}
    save(path / 'worker_started.json', {**root, 'utc': sweep.stamp()})
    terminal = {**root, 'status': 'FAILED'}
    start = time.perf_counter()
    try:
        result = calculate(pilot, protocol, path, bindings)
        terminal['results_sha256'] = save(path / 'results.json', result)
        terminal['status'] = 'FINISHED'
    except BaseException as e:
        terminal['error'] = {'type': type(e).__name__, 'message': str(e)}
        save(path / 'failure.json', {**terminal['error'], 'traceback': traceback.format_exc()})
        raise
    finally:
        save(path / 'terminal.json', {**terminal, 'elapsed_seconds': time.perf_counter() - start,
                                       'utc': sweep.stamp()})


def admit_child(exit_code, terminal, identity):
    require(type(exit_code) is int and exit_code == 0 and terminal.get('status') == 'FINISHED'
            and all(terminal.get(k) == v for k, v in identity.items()), 'failed or mixed worker terminal')


def launch(review_path, review_sha):
    protocol, bindings = prelaunch(review_path, review_sha)
    ACTUAL.mkdir(exist_ok=False)
    pid = os.getpid()
    pilot_ids = [p['pilot_id'] for p in protocol['ordered_pilots']]
    save(ACTUAL / 'started.json', {'pid': pid, 'bindings': bindings,
                                   'ordered_pilot_ids': pilot_ids, 'utc': sweep.stamp()})
    terminal = {'status': 'FAILED', 'pid': pid, 'bindings': bindings, 'workers': []}
    results = []
    try:
        for pilot_id in pilot_ids:
            folder = ACTUAL / pilot_id
            folder.mkdir(exist_ok=False)
            command = [sys.executable, '-B', str(SOURCE), 'worker', '--review', str(review_path),
                       '--review-sha256', review_sha, '--pilot-id', pilot_id, '--supervisor-pid', str(pid)]
            with (folder / 'stdout.log').open('x', encoding='utf-8') as log:
                child = subprocess.Popen(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
                identity = {'pid': child.pid, 'supervisor_pid': pid, 'pilot_id': pilot_id, 'bindings': bindings}
                save(folder / 'child_started.json', {**identity, 'command': command, 'utc': sweep.stamp()})
                code = child.wait()
            worker_terminal = folder / 'terminal.json'
            observed = {**identity, 'observed_exit_code': code,
                        'terminal_sha256': sha(worker_terminal) if worker_terminal.is_file() else None,
                        'utc': sweep.stamp()}
            h = save(folder / 'observed_os_exit.json', observed)
            terminal['workers'].append({**observed, 'os_exit_receipt_sha256': h})
            require(observed['terminal_sha256'] is not None, 'worker failed before terminal')
            done = pinned(worker_terminal, observed['terminal_sha256'])
            admit_child(code, done, identity)
            started = read(folder / 'worker_started.json')
            require(all(started.get(k) == v for k, v in identity.items()), 'worker start identity mismatch')
            result = pinned(folder / 'results.json', done['results_sha256'])
            require(result['pilot']['pilot_id'] == pilot_id and result['bindings'] == bindings
                    and result['status'] in ('CANDIDATE_CERTIFIED_CACHED_QUADRATURE_LOCAL_ROOT',
                                             'UNRESOLVED_CACHED_CONTRACTION_OR_BUDGET'),
                    'worker result mixed with a different pilot or scientific status')
            results.append({'pilot_id': pilot_id, 'status': result['status'], 'results_sha256': done['results_sha256'],
                            'cache_packet': result['cache_packet'], 'cache_packet_sha256': result['cache_packet_sha256']})
        require(len(results) == 4, 'partial cohort')
        terminal['results_sha256'] = save(ACTUAL / 'results.json', {
            'status': 'COMPLETE_FOUR_CACHED_QUADRATURE_PILOT_ATTEMPTS', 'scope': SCOPE,
            'bindings': bindings, 'pilots': results, 'independent_actual_review_required': True})
        terminal['status'] = 'FINISHED'
    except BaseException as e:
        terminal['error'] = {'type': type(e).__name__, 'message': str(e)}
        raise
    finally:
        save(ACTUAL / 'terminal.json', {**terminal, 'utc': sweep.stamp()})
        print(json.dumps(terminal), flush=True)


def rejects(action):
    try:
        action()
    except ValueError:
        return True
    return False


def fixture(protocol):
    import algebraic_root_feasibility as core
    import numpy as np
    from types import SimpleNamespace
    calls = []
    def local_terms(model, u, strength):
        calls.append((model, tuple(u), strength, type(strength)))
        yield {'kind': 'bulk'}
    proof = SimpleNamespace(local_terms=local_terms, assemble_point=lambda terms, n: (
        'gradient', 'Hessian', list(terms), n))
    exact = Q(1, 10)
    point_assembly(proof, 'new_pilot', [2], exact, 1, {'bulk': 1})
    too_small = contraction(Q(1, 5), Q(1, 2), Q(1, 3), Q(1))
    equality = contraction(Q(1, 8), Q(1, 2), Q(1, 4), Q(1, 2))
    strictly = contraction(Q(1, 16), Q(1, 2), Q(1, 4), Q(1, 2))
    goal_eq = contraction(Q(1, 16), Q(1, 2), Q(1, 4), Q(1, 4))
    z_eq = contraction(Q(0), Q(1), Q(1, 4), Q(1, 2))
    expected = {'pid': 19, 'supervisor_pid': 1, 'pilot_id': 'test', 'bindings': {'test_pin': 'digest'}}
    term = {**expected, 'status': 'FINISHED'}
    admit_child(0, term, expected)
    p = protocol['ordered_pilots'][0]
    corrupted = {**p['step'], 'requested_strength_exact': p['exact_drive']['activation_strength_binary_exact']}
    prior_precision = core.ctx.prec
    core.ctx.prec = PRECISION
    core.A.cache_clear()
    try:
        g = core.arb_mat(2, 1, [core.A(Q(-1, 16)), core.A(Q(1, 32))])
        H = core.arb_mat(2, 2, [core.A(Q(-1, 2)), 0, 0, core.A(Q(1, 4))])
        eta, z0, _, point = preconditioned_point(core, g, H, np.diag([-2., 4.]))
    finally:
        core.ctx.prec = prior_precision
        core.A.cache_clear()
    return {
        'four_coherent_field_and_budget_sources_without_actual_cached_sites': (
            len(protocol['ordered_pilots']) == 4
            and [p['checkpoint'] for p in protocol['ordered_pilots']] == ['small_17', 'large_03'] * 2
            and len({p['field_sha256'] for p in protocol['ordered_pilots']}) == 4
            and PINS[budgets.NATIVE_REVIEW] == budgets.ACTUAL_PINS[budgets.NATIVE_REVIEW]),
        'exact_Fraction_passed_to_cached_local_terms_without_a_float_conversion': (
            calls == [('new_pilot', (2,), exact, Q)] and exact != Q.from_float(float(exact))
            and rejects(lambda: point_assembly(proof, 'new_pilot', [2], .1, 1, {'bulk': 1}))),
        'reject_binary_solver_drive_substitution_for_actual_fraction': rejects(
            lambda: exact_drive(corrupted, p['material'])),
        'reject_cross_pilot_exact_material_gain': rejects(
            lambda: exact_drive(p['step'], protocol['ordered_pilots'][-1]['material'])),
        'strict_radius_and_posterior_equality_rules': (
            not any(x['strict_selfmap_contraction_and_observer_radius'] for x in (too_small, equality, goal_eq, z_eq))
            and strictly['strict_selfmap_contraction_and_observer_radius']
            and Q(strictly['posterior_error_upper_exact']) == Q(1, 8)),
        'strict_descending_frozen_radii': all(a > b > 0 for a, b in zip(RADII, RADII[1:])),
        'tiny_exact_signed_preconditioner_fixture': (
            eta == Q(1, 8) and z0 == 0 and Q(point['preconditioner_infinity_norm_upper_exact']) == 4),
        'reject_changed_free_scope_and_bad_child_OS_exit': (
            rejects(lambda: point_assembly(proof, 'new_pilot', [2], exact, 1, {'bulk': 2}))
            and rejects(lambda: admit_child(1, term, expected))
            and rejects(lambda: admit_child(False, term, expected))
            and rejects(lambda: admit_child(0, {**term, 'pilot_id': 'other'}, expected))),
    }


def prepare():
    PREP.mkdir(parents=True, exist_ok=False)
    with (PREP / 'transfer_ambiguity_cached_pilots_frozen.py').open('xb') as handle:
        handle.write(SOURCE.read_bytes())
    try:
        original = Path.open
        def forbid_payload(path, *args, **kwargs):
            resolved = path.resolve()
            if resolved.is_relative_to(ROOT / sweep.PRIMARY_ACTUAL):
                require(resolved.name not in ('field.json', 'verification.json'),
                        'preflight cannot open actual field or verification payload')
            return original(path, *args, **kwargs)
        with (patch.object(Path, 'open', forbid_payload),
              patch.object(subprocess, 'Popen', side_effect=AssertionError('preparation cannot launch worker'))):
            protocol = prepare_protocol()
            h = save(PREP / 'protocol.json', protocol)
            checks = fixture(protocol)
        checks['no_actual_output_or_process_preparation'] = not ACTUAL.exists()
        status = 'PASS' if sweep.all_true(checks) else 'FAIL'
        check_hash = save(PREP / 'preflight.json', {
            'status': status, 'checks': checks, 'source_sha256': sha(SOURCE), 'protocol_sha256': h,
            'actual_primary_field_payloads_opened': 0, 'actual_preconditioners_recovered': 0,
            'actual_cached_sites_assembled': 0, 'new_nonlinear_work': 0})
        require(status == 'PASS', 'cached-pilot preflight failed')
        print(json.dumps({'status': status, 'source_sha256': sha(SOURCE),
                          'protocol_sha256': h, 'preflight_sha256': check_hash}), flush=True)
    except BaseException as error:
        save(PREP / 'failure.json', {'type': type(error).__name__, 'message': str(error),
                                     'traceback': traceback.format_exc()})
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    command = parser.add_subparsers(dest='command', required=True)
    command.add_parser('prepare')
    for tag in ('launch', 'worker'):
        sub = command.add_parser(tag)
        sub.add_argument('--review', required=True, type=Path)
        sub.add_argument('--review-sha256', required=True)
        if tag == 'worker':
            sub.add_argument('--pilot-id', required=True)
            sub.add_argument('--supervisor-pid', required=True, type=int)
    args = parser.parse_args()
    if args.command == 'prepare':
        prepare()
    elif args.command == 'launch':
        launch(args.review.resolve(), args.review_sha256)
    else:
        worker(args.review.resolve(), args.review_sha256, args.pilot_id, args.supervisor_pid)


if __name__ == '__main__':
    main()

