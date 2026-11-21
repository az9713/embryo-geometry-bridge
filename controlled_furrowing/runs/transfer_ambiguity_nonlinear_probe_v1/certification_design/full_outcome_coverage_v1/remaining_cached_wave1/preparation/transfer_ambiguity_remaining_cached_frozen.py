"""Prepare and execute two fresh cached-quadrature root attempts for the frozen M states."""

import argparse
from copy import deepcopy
from fractions import Fraction as Q
from hashlib import sha256
import importlib.metadata
import json
from math import isfinite
import os
from pathlib import Path
import subprocess
import sys
import traceback


SOURCE = Path(__file__).resolve()
ROOT = SOURCE.parents[1]
RUN = ROOT/'runs/transfer_ambiguity_nonlinear_probe_v1'
DESIGN = RUN/'certification_design'
INVENTORY = DESIGN/'full_outcome_coverage_v1'
OUT = INVENTORY/'remaining_cached_wave1'
PREP = OUT/'preparation'
ACTUAL = OUT/'actual'
RECORD = 'optoGEF2_small_col_M'
CHECKPOINT = 'small_09'
DOSE = Q(78208913, 100000000)
ENDPOINTS = ('max_contrast_change', 'max_negative_contrast_change')
CASES = tuple(f'{endpoint}_{CHECKPOINT}' for endpoint in ENDPOINTS)
OBSERVER = 'current_disk__spatial_max'
OBS = ('maximum_shortening_um', 'region_averaged_shortening_um',
       'apical_max_shortening_um', 'basal_max_shortening_um')
ORDER = ('apical_area_modulus', 'basal_area_modulus', 'apical_surface_tension',
         'basal_surface_tension', 'common_gain')
NFREE = 611
TERMS = {'bulk': 144*64, 'surface': 24*12}
RADII = tuple(Q(1, 2**k) for k in (20, 24, 28, 32, 36))
PRECISION = 256
CP_SOURCE = ROOT/'scripts/transfer_ambiguity_cached_pilots.py'
CORE_SOURCE = ROOT/'scripts/algebraic_root_feasibility.py'
CORE_REVIEW = ROOT/'evidence/algebraic_root_independent_review.json'
PRIMARY = RUN/'preparation/protocol.json'
FIELD_MANIFEST = RUN/'field_manifest.json'
PRIMARY_REVIEW = RUN/'independent_actual_preparation/actual_replay/receipt.json'
BUDGET = DESIGN/'budgets_v1/actual/results.json'
BUDGET_TERMINAL = DESIGN/'budgets_v1/actual/terminal.json'
BUDGET_REVIEW = DESIGN/'budgets_v1/independent_actual/receipt.json'
LEDGER = INVENTORY/'receipt.json'
PINS = {
    LEDGER: '235e8336c735b1bf1a75889018ae6a18586fd0905a65d38cb683c5ed236df19a',
    PRIMARY: '9d5d38bb5f9785d40bf843dd225eb93b1d5ecbe96f20942aec0023765d119d94',
    FIELD_MANIFEST: '308ca4c05e476ba6a14cbc7b8d34018983b8f11263e2c224488672b35eb38ac0',
    PRIMARY_REVIEW: 'c795ae78ba0b53a01c3d65f466d5fa2f1d233c2a68c4275c04260c3cb5120ebd',
    BUDGET: '689948766db950ee4ddf399e6731aea08890cc7a1f84c5615d19e7ded68b0872',
    BUDGET_TERMINAL: '83abadac599cba377034a522aa5c160de176b4b8b6b91f6223883415495d9928',
    BUDGET_REVIEW: '978a2a4095bb8735103a5ee2c11ad266a6b0e052e22af2b7a36ef29f1f38265c',
    CP_SOURCE: '1b22ebc307c0cee51bf4193ff3999242b1cfe278d07c10cc63619f2795489c58',
    CORE_SOURCE: '7321d38bf17cf354cc03be2871f3fbaa45bcebbcbfca7be90d0f867b80078acf',
    CORE_REVIEW: 'f327f2fab76f9fdcab3a83882adbb1a4ae3f0cf1e5c269e7335fb909691085d5',
    ROOT/'scripts/transfer_ambiguity_error_budgets.py':
        '410a4944798b1f42f0360cbdf236d9da68f8de2615d93f860bccdd8be54d8859',
    ROOT/'scripts/transfer_ambiguity_observer_sweep.py':
        '9eacb1de680cefb84d90d603cc1aca4dd6ef4892372a91cffa3480651c5874a8',
    ROOT/'evidence/nodal_observer_error_independent_review.json':
        '335626b91c72ac5ec37b1b1dd95b4825d3dc847df02d9cf993e5ed6daff5f35a',
    ROOT/'scripts/nodal_observer_error_transfer.py':
        '001769f23ba386fe6bde4530c12e011f72f8733bb1a4a2923182796e54da35de',
    ROOT/'scripts/certify_discrete_geometry.py':
        'ecf329847b8d00f8e33e4ea6ac4f83f3af0c2569528f387a82f93f17b8e5d5ff',
}
SCOPE = ('Two fixed-parameter pointwise binary Q8/S12 cache root attempts at the exact requested drives '
         'and saved M centers. The derivative bound covers the same closed free cube at cached sites. '
         'Spatial integration, whole-body geometry, stability, observer transfer and PDE conclusions are separate.')


def relative(path):
    return path.relative_to(ROOT).as_posix()


def require(ok, message):
    if not ok:
        raise ValueError(message)


def unique(rows, predicate):
    rows = [r for r in rows if predicate(r)]
    require(len(rows) == 1, 'exactly one matching authority row required')
    return rows[0]


def all_true(checks):
    return isinstance(checks, dict) and bool(checks) and all(x is True for x in checks.values())


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def read(path, digest):
    raw = path.read_bytes()
    require(sha256(raw).hexdigest() == digest, f'changed frozen bytes: {path}')
    return json.loads(raw)


def write(path, value):
    with path.open('x', encoding='utf-8', newline='\n') as handle:
        handle.write(json.dumps(value, indent=2, allow_nan=False)+'\n')
    return sha(path)


def authority_inputs():
    p, m, n, b, end, bn, l = [read(path, PINS[path]) for path in (
        PRIMARY, FIELD_MANIFEST, PRIMARY_REVIEW, BUDGET, BUDGET_TERMINAL, BUDGET_REVIEW, LEDGER)]
    require(p['status'] == 'FROZEN_OUTCOME_FREE_NONLINEAR_DIAGNOSTIC'
            and p['observable_variant'] == OBSERVER
            and n['status'] == 'PASS_INDEPENDENT_ACTUAL_NONLINEAR_STORED_FIELD_REPLAY'
            and n['checkpoints'] == len(m['fields']) == 44
            and m['independent_primary_review_sha256'] == PINS[PRIMARY_REVIEW]
            and type(m['independent_replay_observed_os_exit_code']) is int
            and m['independent_replay_observed_os_exit_code'] == 0
            and end['status'] == 'FINISHED' and end['results_sha256'] == PINS[BUDGET]
            and end['bindings'] == b['bindings'] and end['pid'] == bn['observed_budget_terminal']['pid']
            and type(bn['observed_budget_terminal']['exit_code']) is int
            and bn['observed_budget_terminal']['exit_code'] == 0
            and b['status'] == 'COMPLETE_POSITIVE_CONDITIONAL_36_FIELD_BUDGETS'
            and bn['status'] == 'PASS_INDEPENDENT_ACTUAL_CONDITIONAL_OBSERVER_BUDGETS'
            and all_true(bn['all_literal_checks']) and b['exact_root_error_bounds_proved'] == 0
            and l['status'] == 'PASS_READ_ONLY_36_FIELD_COVERAGE_INVENTORY'
            and tuple(l['recommendation_first_pair_case_ids']) == CASES
            and len(l['outcome_states']) == len(b['field_budgets']) == len(bn['fields']) == 36,
            'frozen complete field and budget authority is required')
    for path in (PRIMARY, FIELD_MANIFEST, PRIMARY_REVIEW, BUDGET, BUDGET_TERMINAL, BUDGET_REVIEW):
        require(l['readback_sha256'][relative(path)] == PINS[path], 'mixed inventory authority input')
    for path in (FIELD_MANIFEST, PRIMARY_REVIEW, BUDGET, BUDGET_TERMINAL):
        require(bn['source_and_input_hashes'][relative(path)] == PINS[path], 'mixed independent budget input')
    for name in ('actual/results.json','actual/terminal.json'):
        path = RUN/name
        expected = m['primary_results_sha256' if name.endswith('results.json') else 'primary_terminal_sha256']
        require(n['all_hashes'][relative(path)] == bn['source_and_input_hashes'][relative(path)] == expected
                and sha(path) == expected, 'primary true terminal binding differs')
    require(Q(p['selection']['conditional_source_allowance_um_exact'])
            == Q(b['conditional_source_allowance_um_exact']) == Q(l['conditional_allowance_um_exact']) > 0,
            'observer allowance changed')
    return p, m, n, b, bn, l


def case_inputs(case_id, authorities):
    p, m, n, b, bn, l = authorities
    require(case_id in CASES, 'case is outside the two frozen M states')
    endpoint = ENDPOINTS[CASES.index(case_id)]
    def same(row):
        return row['endpoint_id'] == endpoint and row['checkpoint'] == CHECKPOINT
    ep = unique(p['endpoints'], lambda row: row['endpoint_id'] == endpoint)
    f = unique(m['fields'], same)
    r = unique(b['field_budgets'], same)
    independent = unique(bn['fields'], same)
    ledger = unique(l['outcome_states'], same)
    return ep, f, r, independent, ledger


def case_admission(case_id, ep, field, row, independent, ledger, state, verification, delta):
    require(case_id in CASES, 'case is outside the two frozen M states')
    endpoint = ENDPOINTS[CASES.index(case_id)]
    require(ep['endpoint_id'] == endpoint and ledger['case_id'] == case_id
            and all(r['endpoint_id'] == endpoint and r['checkpoint'] == CHECKPOINT
                    for r in (field, row, independent, ledger))
            and all(r['role'] == 'small_record' for r in (field, row, ledger))
            and row['record_id'] == independent['record_id'] == ledger['record_id'] == RECORD
            and ledger['fresh_cache_and_preconditioner_needed'] is True
            and ledger['pilot_cached_and_geometry_coverage'] is None,
            'mixed candidate/checkpoint or reused cache authority')
    for key in ('field', 'field_sha256', 'verification', 'verification_sha256'):
        require(field[key] == row[key] == ledger[key], 'saved field or verification identity changed')
    require(independent['field'] == field['field'] and independent['field_sha256'] == field['field_sha256']
            and field['field'] == f'{relative(RUN)}/actual/{endpoint}/{CHECKPOINT}/field.json'
            and field['verification'] == f'{relative(RUN)}/actual/{endpoint}/{CHECKPOINT}/verification.json',
            'wrong independent field binding')
    step = unique(ep['small_records'], lambda r: r['record_id'] == RECORD)
    material = ep['materials']['small']
    center = row['center_and_model']
    require(step == row['step'] == field['step'] and Q(step['dose_exact']) == Q(ledger['dose_exact']) == DOSE,
            'different exact physical dose or step')
    params = material['parameters']
    require(params == state['parameters'] == center['material_parameters']
            and material['parameter_storage_exact'] == center['material_parameter_storage_exact']
            and Q(params['mu']) == 1 and Q(params['bulk']) == 5000,
            'saved material parameters differ from the frozen independently reviewed row')
    gain = Q(material['gain_exact'])
    drive = gain*DOSE
    rounded = Q(float(drive))
    require(drive > 0 and drive == Q(step['requested_strength_exact'])
            == Q(center['exact_target_strength']) == Q(ledger['requested_drive_exact'])
            and rounded == Q(step['solver_strength_binary64_exact']) == Q(field['actual_strength_binary64_exact'])
            == Q(state['activation_strength']) == Q(center['solver_binary64_strength_exact'])
            == Q(ledger['solver_binary_drive_exact']) and rounded != drive
            and rounded-drive == Q(step['binary64_rounding_difference_exact'])
            == Q(ledger['solver_minus_exact_drive_exact'])
            == Q(center['binary64_strength_rounding_error_exact']),
            'exact drive and its separate binary64 value must be jointly consistent')
    phys = {k: Q(params[k]) for k in ORDER[:-1]}
    require(phys['apical_area_modulus'] > 0 and phys['basal_area_modulus'] > 0 and min(phys.values()) >= 0,
            'invalid fixed material coefficients')
    phys['common_gain'] = gain
    free = center['free_map']
    require(len(free) == len(set(free)) == NFREE and all(type(x) is int and 0 <= x < 650 for x in free)
            and len(state['free_displacement']) == NFREE and all(isfinite(x) for x in state['free_displacement'])
            and len(state['reference_nodes']) == len(state['current_nodes']) == 325
            and len(state['triangles']) == len(state['activation_elements']) == 144,
            'saved mesh and unique free coordinate layout required')
    active = [i for i,on in enumerate(state['activation_elements']) if on]
    require(len(active) == 24 and active == center['active_elements']
            and all(not set(face['radial_dofs']) & set(free) for face in center['fixed_trace_endpoints'].values()),
            'active mask or fixed radial endpoint layout changed')
    goal, xi = Q(row['e_goal_dimensionless_exact']), Q(center['node_addition_gap_dimensionless_exact'])
    require(row['status'] == 'POSITIVE_CONDITIONAL_ERROR_TARGET'
            and goal == Q(independent['e_goal_dimensionless_exact']) == Q(ledger['observer_goal_exact'])
            and any(r < goal for r in RADII)
            and 0 <= xi == Q(independent['binary_gap_exact']) == Q(ledger['node_addition_gap_exact'])
            and verification['variant'] == OBSERVER and set(row['stored_residuals']) == set(OBS)
            and set(verification['observer_intervals_um_exact']) == set(OBS),
            'positive reviewed observer budget and same saved convention required')
    for k in OBS:
        v = row['stored_residuals'][k]
        obs = Q(step['observed_exact'][k])
        a, b = map(Q, verification['observer_intervals_um_exact'][k])
        margin = delta-max(abs(a-obs),abs(b-obs))
        require(a <= b and obs == Q(v['observed_um_exact'])
                and (a,b) == tuple(map(Q,v['stored_interval_um_exact']))
                and margin == Q(v['margin_um_exact']) == Q(ledger['stored_primary_observer_margins_um_exact'][k]) > 0,
                'saved observer interval, constraint or allowance is mixed')
    # The two authorities supply the independent budget. Also verify the budget algebra here.
    caps = [Q(row['per_observer_allocations_um_exact'][k])/Q(row['conditional_observer_coefficient_um_exact'][k]) for k in OBS]
    require(all(Q(row['per_observer_allocations_um_exact'][k]) == Q(row['stored_residuals'][k]['margin_um_exact'])/2 for k in OBS)
            and goal+xi == min(*caps, Q(row['monotonicity_total_dimensionless_cap_exact'])),
            'observer budget differs from its strict accounting identity')
    return {'case_id': case_id, 'endpoint_id': endpoint, 'checkpoint': CHECKPOINT, 'record_id': RECORD,
            'field': field['field'], 'field_sha256': field['field_sha256'],
            'verification': field['verification'], 'verification_sha256': field['verification_sha256'],
            'step': step, 'material': material, 'free_map': free, 'active_elements': active,
            'xi_dimensionless_exact': str(xi), 'e_goal_dimensionless_exact': str(goal),
            'closed_parameter_cell_exact': {k:[str(x),str(x)] for k,x in phys.items()},
            'parameter_order': list(ORDER), 'observer_variant': OBSERVER,
            'frozen_conditional_allowance_um_exact': str(delta),
            'drive': {'optical_dose_exact': str(DOSE), 'common_gain_exact': str(gain),
                      'activation_strength_exact': str(drive), 'activation_strength_binary_exact': str(rounded),
                      'binary_input_minus_exact_strength': str(rounded-drive)}}


def prepare_protocol():
    for path,digest in PINS.items():
        require(sha(path) == digest, f'frozen input changed: {relative(path)}')
    authorities = authority_inputs()
    p, _, native, _, _, ledger = authorities
    core_review = read(CORE_REVIEW,PINS[CORE_REVIEW])
    require(core_review['status'] == 'PASS' and all_true(core_review['checks'])
            and core_review['source_sha256'] == PINS[CORE_SOURCE], 'cached algebra review changed')
    packages = {name:importlib.metadata.version(name) for name in ('numpy','scipy','torch','python-flint')}
    require(all(packages[k] == p['runtime'][k] for k in ('numpy','scipy','torch')),
            'same core binary package versions required')
    source_hashes = {**p['source_hashes'], **{relative(path):d for path,d in PINS.items() if path.suffix == '.py'}}
    for name,h in source_hashes.items():
        require(sha(ROOT/name) == h, f'pinned source closure changed: {name}')
    inputs = {relative(path):d for path,d in PINS.items() if path.suffix != '.py'}
    cases = []
    delta = Q(p['selection']['conditional_source_allowance_um_exact'])
    for case_id in CASES:
        records = case_inputs(case_id,authorities)
        field = records[1]
        file_data = []
        for k in ('field','verification'):
            path,digest = field[k],field[k+'_sha256']
            require(native['all_hashes'][path] == ledger['readback_sha256'][path] == digest,
                    'selected payload lacks independent byte authority')
            file_data.append(read(ROOT/path,digest))
            inputs[path] = digest
        cases.append(case_admission(case_id,*records,*file_data,delta))
    require(len({c['field_sha256'] for c in cases}) == 2, 'duplicate saved centers')
    return {'status':'FROZEN_TWO_REMAINING_M_CACHED_QUADRATURE_TARGETS','scope':SCOPE,
            'source_sha256':sha(SOURCE),'input_sha256':dict(sorted(inputs.items())),
            'source_closure_sha256':dict(sorted(source_hashes.items())),
            'ordered_case_ids':list(CASES),'cases':cases,'numerical_policy':p['numerical_policy'],
            'runtime_packages':packages,'precision_bits':PRECISION,'torch_cpu_threads':1,
            'finite_radius_ladder_exact':list(map(str,RADII)),'site_counts':TERMS,
            'decision':'Strict r < independent e_goal, z < 1, eta + zr < r and eta/(1-z) < e_goal with positive cached analytic domain.',
            'preconditioner':'Recover Q8/S12 Hessian at the saved vector and nearest binary64 drive only to propose a fresh binary inverse.',
            'exact_cached_force_drive':'Exact protocol gain times exact decimal dose, as a Fraction.',
            'execution':'One process at a time in declared order, observed OS exits and exclusively created outputs.',
            'requires_independent_prelaunch_and_actual_reviews':True,
            'new_nonlinear_equilibrium_solves':0,'old_pilot_case_reused':False}


def rejects(action):
    try:
        action()
    except ValueError:
        return True
    return False


def admit_cached_worker(code, finished, started, expected):
    require(type(code) is int and code == 0 and finished['status'] == 'FINISHED'
            and started == expected and all(finished[k] == v for k,v in expected.items()),
            'worker true OS exit or identity changed')


def fixtures(protocol):
    source = authority_inputs()
    args = list(case_inputs(CASES[0],source))
    field = args[1]
    args.extend(read(ROOT/field[k],field[k+'_sha256']) for k in ('field','verification'))
    args.append(Q(protocol['cases'][0]['frozen_conditional_allowance_um_exact']))
    def mutate(index,key,value):
        fixture = deepcopy(args)
        fixture[index][key] = value
        return rejects(lambda:case_admission(CASES[0],*fixture))
    valid = case_admission(CASES[0],*args)
    mixed_field = case_inputs(CASES[1],source)[1]
    alternate = deepcopy(args[0]['materials'])
    alternate['small']['gain_exact'] = '1'
    worker_id = {'pid':3,'supervisor_pid':2,'case_id':CASES[0],'bindings':{'test':'mock'}}
    worker_terminal = {**worker_id,'status':'FINISHED'}
    admit_cached_worker(0,worker_terminal,worker_id,worker_id)
    return {'admit_two_explicit_distinct_M_states':valid == protocol['cases'][0]
                and all(x['record_id'] == RECORD for x in protocol['cases'])
                and len({x['field_sha256'] for x in protocol['cases']}) == 2,
            'reject_case_outside_preselected_pair':rejects(lambda:case_inputs('max_contrast_change_small_17',source)),
            'reject_other_endpoint_metadata_or_field':(rejects(lambda:case_admission(CASES[1],*args))
                and mutate(1,'field',mixed_field['field']) and mutate(1,'field_sha256',mixed_field['field_sha256'])),
            'reject_changed_raw_bound_file_bytes':rejects(lambda:read(ROOT/field['field'],mixed_field['field_sha256'])),
            'reject_changed_exact_gain_or_rounded_drive_substitution':(mutate(0,'materials',alternate)
                and mutate(4,'requested_drive_exact',valid['drive']['activation_strength_binary_exact'])),
            'reject_changed_budget_or_observer_convention':(mutate(2,'e_goal_dimensionless_exact','1')
                and mutate(3,'e_goal_dimensionless_exact','0')
                and mutate(6,'variant','current_radial_line__regional_mean')),
            'fixed_positive_exact_budget_and_degenerate_cell':all(
                any(r<Q(x['e_goal_dimensionless_exact']) for r in RADII)
                and all(lo==hi for lo,hi in x['closed_parameter_cell_exact'].values()) for x in protocol['cases']),
            'reject_boolean_nonzero_exit_or_other_worker_identity':(
                rejects(lambda:admit_cached_worker(False,worker_terminal,worker_id,worker_id))
                and rejects(lambda:admit_cached_worker(1,worker_terminal,worker_id,worker_id))
                and rejects(lambda:admit_cached_worker(0,{**worker_terminal,'case_id':CASES[1]},worker_id,worker_id))),
            'no_actual_outputs':not ACTUAL.exists()}


def prepare():
    PREP.mkdir(parents=True,exist_ok=False)
    with (PREP/'transfer_ambiguity_remaining_cached_frozen.py').open('xb') as f:
        f.write(SOURCE.read_bytes())
    try:
        protocol = prepare_protocol()
        checks = fixtures(protocol)
        require(all_true(checks),'remaining-outcome cached preparation failed')
        p = write(PREP/'protocol.json',protocol)
        f = write(PREP/'preflight.json',{'status':'PASS','source_sha256':sha(SOURCE),'protocol_sha256':p,
                     'checks':checks,'actual_model_instances_constructed':0,'actual_preconditioners_recovered':0,
                     'actual_binary_caches_frozen':0,'actual_cached_AD_sites_assembled':0,'new_equilibrium_solves':0,
                     'initial_dry_run':'Rejected a mistyped budget terminal SHA before admitting any inputs. The corrected pin matches the independent budget review.'})
        print(json.dumps({'source_sha256':sha(SOURCE),'protocol_sha256':p,'preflight_sha256':f}),flush=True)
    except BaseException as error:
        write(PREP/'failure.json',{'error':str(error),'traceback':traceback.format_exc()})
        raise


def prelaunch(path, digest):
    r = read(path,digest)
    expected = {'source_sha256':sha(SOURCE),'protocol_sha256':sha(PREP/'protocol.json'),
                'preflight_sha256':sha(PREP/'preflight.json'),'coverage_ledger_sha256':PINS[LEDGER],
                'whole_44_field_manifest_sha256':PINS[FIELD_MANIFEST],
                'budget_independent_actual_review_sha256':PINS[BUDGET_REVIEW]}
    require(r['status'] == 'PASS_INDEPENDENT_REMAINING_M_CACHED_PRELAUNCH' and all_true(r['checks'])
            and all(r[k] == v for k,v in expected.items()) and tuple(r['ordered_case_ids']) == CASES,
            'independent prelaunch review must bind this exact driver, pair and protocol')
    prep = read(PREP/'preflight.json',expected['preflight_sha256'])
    require(prep['status'] == 'PASS' and all_true(prep['checks'])
            and prep['source_sha256'] == expected['source_sha256'] and prep['protocol_sha256'] == expected['protocol_sha256']
            and all(type(prep[k]) is int and prep[k] == 0 for k in ('actual_model_instances_constructed',
                'actual_preconditioners_recovered','actual_binary_caches_frozen','actual_cached_AD_sites_assembled','new_equilibrium_solves'))
            and sha(PREP/'transfer_ambiguity_remaining_cached_frozen.py') == expected['source_sha256'],
            'mixed or non-cold preflight')
    protocol = read(PREP/'protocol.json',expected['protocol_sha256'])
    require(protocol == prepare_protocol(),'frozen remaining-outcome protocol no longer matches inputs')
    return protocol,{**expected,'independent_prelaunch_review_sha256':digest}


def actual_case(case,protocol,folder,bindings):
    # prelaunch authenticated all source bytes before importing the numerical stack.
    import algebraic_root_feasibility as core
    import torch
    import transfer_ambiguity_cached_pilots as cp
    import transfer_pair_pilot as model_source
    for name in protocol['source_closure_sha256']:
        module = sys.modules.get(Path(name).stem)
        if module is not None:
            require(Path(module.__file__).resolve()==ROOT/name,f'runtime source alias drift: {name}')
    require(Path(core.__file__).resolve()==CORE_SOURCE and Path(cp.__file__).resolve()==CP_SOURCE
            and Path(model_source.__file__).resolve()==ROOT/'scripts/transfer_pair_pilot.py'
            and cp.NFREE==NFREE and cp.TERMS==TERMS and cp.RADII==RADII,
            'local model/cache kernel scope changed')
    torch.set_num_threads(1)
    core.ctx.prec = PRECISION
    core.A.cache_clear()
    require(torch.get_num_threads()==1 and core.ctx.prec==PRECISION,'floating package policy changed')
    model,u = cp.model_input(case,protocol['numerical_policy'])
    drive = case['drive']
    strength = Q(drive['activation_strength_exact'])
    require(drive==cp.exact_drive(case['step'],case['material'])
            and strength==DOSE*Q(case['closed_parameter_cell_exact']['common_gain'][0]),
            'fresh packet drive does not match independent case')
    R,proposal = cp.sparse_proposal(model,u,drive,folder)
    payload = core.freeze_cache(model,u,R)
    packet = {'status':'FROZEN_FRESH_EXACT_BINARY_QUADRATURE_CACHE','scope':SCOPE,
              'case_id':case['case_id'],'record_id':case['record_id'],'input_state':case['field'],
              'input_state_sha256':case['field_sha256'],'field_verification_sha256':case['verification_sha256'],
              'parameters':case['material']['parameters'],'material_parameter_storage_exact':case['material']['parameter_storage_exact'],
              'closed_parameter_cell_exact':case['closed_parameter_cell_exact'],'parameter_order':case['parameter_order'],
              'integration_policy':protocol['numerical_policy']['solve_quadrature'],'free_dof':NFREE,
              'normalization_binary_exact':str(Q(float(model.energy_scale))),
              'precision_bits':PRECISION,**drive,**proposal,'frozen_binary_cache':payload,'bindings':bindings}
    path=folder/'cache_packet.json'
    ph=write(path,packet)
    (g,H,sites,min_j),counts=cp.point_assembly(core,model,u,strength,NFREE,TERMS)
    require(Q(core.pack(min_j)['lower_exact'])>Q(1,10**6),'point cached analytic threshold failed')
    eta,z0,abs_R,rows=cp.preconditioned_point(core,g,H,R)
    enclosure=write(folder/'point_enclosure.json',{'case_id':case['case_id'],'cache_packet_sha256':ph,**rows,
                   'expected_and_rebuilt_site_counts':counts,'point_quadrature_J_min':core.pack(min_j),'bindings':bindings})
    attempts=[]
    accepted=None
    goal=Q(case['e_goal_dimensionless_exact'])
    for i,r in enumerate(RADII):
        test=cp.test_radius(core,sites,abs_R,eta,z0,r,goal)
        filename=f'radius_{i:02d}.json'
        h=write(folder/filename,{'case_id':case['case_id'],'cache_packet_sha256':ph,**test})
        attempts.append({'path':filename,'sha256':h,**test})
        print(json.dumps({'case_id':case['case_id'],'radius':str(r),'status':test['status']}),flush=True)
        if test['strict_selfmap_contraction_and_observer_radius'] is True:
            accepted=test
            break
    return {'status':'CANDIDATE_CACHED_ROOT_AWAITING_INDEPENDENT_REVIEW' if accepted else
            'UNRESOLVED_CACHED_CONTRACTION_OR_BUDGET','scope':SCOPE,'case':case,'bindings':bindings,
            'cache_packet':relative(path),'cache_packet_sha256':ph,'embedded_binary_cache_sha256':payload['sha256'],
            'point_enclosure_sha256':enclosure,'eta_upper_exact':str(eta),'center_Z0_upper_exact':str(z0),
            'radius_attempts':attempts,'accepted_radius':accepted,'backend':{
                'arb_precision_bits':int(core.ctx.prec),'torch_threads':torch.get_num_threads(),'packages':protocol['runtime_packages']},
            'new_nonlinear_equilibrium_solves':0,'exact_spatial_integration_established':False,
            'whole_body_geometry_established':False,'local_minimum_established':False,'independent_actual_review_required':True}


def worker(review,digest,case_id,parent_pid):
    protocol,bindings=prelaunch(review,digest)
    opening=json.loads((ACTUAL/'started.json').read_bytes())
    require(os.getppid()==parent_pid==opening['pid'] and opening['bindings']==bindings
            and opening['ordered_case_ids']==list(CASES),'current live serial parent required')
    case=unique(protocol['cases'],lambda x:x['case_id']==case_id)
    folder=ACTUAL/case_id
    identity={'pid':os.getpid(),'supervisor_pid':parent_pid,'case_id':case_id,'bindings':bindings}
    write(folder/'worker_started.json',identity)
    terminal={**identity,'status':'FAILED'}
    try:
        terminal['results_sha256']=write(folder/'results.json',actual_case(case,protocol,folder,bindings))
        terminal['status']='FINISHED'
    except BaseException as error:
        terminal['error']=str(error)
        write(folder/'failure.json',{'error':str(error),'traceback':traceback.format_exc()})
        raise
    finally:
        write(folder/'terminal.json',terminal)


def launch(review,digest):
    protocol,bindings=prelaunch(review,digest)
    ACTUAL.mkdir(exist_ok=False)
    pid=os.getpid()
    write(ACTUAL/'started.json',{'pid':pid,'bindings':bindings,'ordered_case_ids':list(CASES)})
    terminal={'pid':pid,'bindings':bindings,'status':'FAILED','workers':[]}
    results=[]
    try:
        for case in protocol['cases']:
            case_id=case['case_id']
            folder=ACTUAL/case_id
            folder.mkdir(exist_ok=False)
            command=[sys.executable,'-B',str(SOURCE),'worker','--review',str(review),'--review-sha256',digest,
                     '--case-id',case_id,'--supervisor-pid',str(pid)]
            with (folder/'stdout.log').open('x',encoding='utf-8') as log:
                child=subprocess.Popen(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
                identity={'pid':child.pid,'supervisor_pid':pid,'case_id':case_id,'bindings':bindings}
                write(folder/'child_started.json',{**identity,'command':command})
                code=child.wait()
            path=folder/'terminal.json'
            observed={**identity,'observed_exit_code':code,'terminal_sha256':sha(path) if path.is_file() else None}
            oh=write(folder/'observed_os_exit.json',observed)
            terminal['workers'].append({**observed,'os_exit_receipt_sha256':oh})
            require(type(code) is int and code==0 and observed['terminal_sha256'] is not None,'cached worker failed')
            done=read(path,observed['terminal_sha256'])
            started=json.loads((folder/'worker_started.json').read_bytes())
            admit_cached_worker(code,done,started,identity)
            result=read(folder/'results.json',done['results_sha256'])
            require(result['case']==case and result['bindings']==bindings
                    and result['status'] in ('CANDIDATE_CACHED_ROOT_AWAITING_INDEPENDENT_REVIEW',
                                             'UNRESOLVED_CACHED_CONTRACTION_OR_BUDGET'),
                    'mixed cached case result')
            cache_path=folder/'cache_packet.json'
            require(result['cache_packet']==relative(cache_path)
                    and sha(cache_path)==result['cache_packet_sha256'],'cached packet bytes or ownership changed')
            results.append({'case_id':case_id,'status':result['status'],'results_sha256':done['results_sha256'],
                            'cache_packet':relative(cache_path),'cache_packet_sha256':result['cache_packet_sha256']})
        require(tuple(x['case_id'] for x in results)==CASES,'partial cached M pair')
        terminal['results_sha256']=write(ACTUAL/'results.json',{'status':'COMPLETE_TWO_M_CACHED_ATTEMPTS_UNREVIEWED',
                                         'scope':SCOPE,'cases':results,'bindings':bindings,
                                         'independent_actual_review_required':True})
        terminal['status']='FINISHED'
    except BaseException as error:
        terminal['error']={'type':type(error).__name__,'message':str(error)}
        raise
    finally:
        write(ACTUAL/'terminal.json',terminal)
        print(json.dumps(terminal),flush=True)


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    sub=cli.add_subparsers(dest='stage',required=True)
    sub.add_parser('prepare')
    for name in ('launch','worker'):
        parser=sub.add_parser(name)
        parser.add_argument('--review',required=True,type=lambda value:Path(value).resolve())
        parser.add_argument('--review-sha256',required=True)
        if name=='worker':
            parser.add_argument('--case-id',required=True,choices=CASES)
            parser.add_argument('--supervisor-pid',required=True,type=int)
    args=cli.parse_args()
    if args.stage=='prepare':
        prepare()
    elif args.stage=='launch':
        launch(args.review,args.review_sha256)
    else:
        worker(args.review,args.review_sha256,args.case_id,args.supervisor_pid)
