"""Prepare and execute bounded fresh cached-quadrature root attempts for the remaining 30 outcomes."""

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
OUT = INVENTORY/'remaining_cached_30_v1'
PREP = OUT/'preparation'
ACTUAL = OUT/'actual'
ENDPOINTS = ('max_contrast_change', 'max_negative_contrast_change')
WAVE_IDS = tuple(f'wave_{i:02d}' for i in range(1,9))
OBSERVER = 'current_disk__spatial_max'
OBS = ('maximum_shortening_um', 'region_averaged_shortening_um',
       'apical_max_shortening_um', 'basal_max_shortening_um')
ORDER = ('apical_area_modulus', 'basal_area_modulus', 'apical_surface_tension',
         'basal_surface_tension', 'common_gain')
NFREE = 611
TERMS = {'bulk': 144*64, 'surface': 24*12}
RADII = tuple(Q(1, 2**k) for k in (20, 24, 28, 32, 34, 36))
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
VERIFIED_M = INVENTORY/'remaining_cached_wave1/VERIFIED-CACHED-COVERAGE.json'
PRIOR_FOUR = DESIGN/'cached_pilots_v1/independent_actual_preparation/actual_replay/receipt.json'
FROZEN_SOURCE = 'transfer_ambiguity_remaining30_cached_frozen.py'
REVIEW_STATUS = 'PASS_INDEPENDENT_REMAINING30_CACHED_WAVE_PRELAUNCH'
PINS = {
    LEDGER: '235e8336c735b1bf1a75889018ae6a18586fd0905a65d38cb683c5ed236df19a',
    VERIFIED_M: 'd5591fe4cfa66ef6e1d25eeaafde6df1000352c94225dd680584dab3018901b1',
    PRIOR_FOUR: '4e891afd79dd2f73508d1f878ddc1113ef107bdff8d9e0f673894f4d2f10ae17',
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
SCOPE = ('Fixed-parameter pointwise binary Q8/S12 cache root attempts at the exact requested drives '
         'and saved centers of 30 small outcomes. The derivative bound covers the same closed free cube at cached sites. '
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
            and len(l['outcome_states']) == len(b['field_budgets']) == len(bn['fields']) == 36,
            'frozen complete field, budget and coverage authority is required')
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


def schedule(ledger, first_four, updated):
    states = {r['case_id']:r for r in ledger['outcome_states']}
    require(len(states)==len(ledger['outcome_states'])==36
            and sum(x['role']=='small_record' for x in states.values())==34,
            'complete unique 34 small and two large ledger states required')
    require(first_four['status']=='PASS_INDEPENDENT_FOUR_CACHED_LOCAL_ROOTS'
            and first_four['root_certificates_cached_quadrature_only']==4
            and first_four['root_certificates_exact_spatial_integral']==0
            and updated['status']=='PASS_VERIFIED_TWO_ADDITIONAL_M_CACHED_LOCAL_ROOTS'
            and all_true(updated['checks']) and updated['new_cached_local_root_count']==2
            and updated['prior_other_cached_local_root_count']==4 and updated['exact_spatial_root_count']==0,
            'independently verified coverage evidence required')
    original = {}
    for r in first_four['field_results']:
        k = r['pilot_id']
        require(k in states and k not in original and r['field_sha256']==states[k]['field_sha256']
                and r['status']=='PASS_INDEPENDENT_CACHED_LOCAL_ROOT'
                and states[k]['pilot_cached_and_geometry_coverage'] is not None,
                'prior cached proof and ledger state differ')
        original[k] = r['receipt_sha256']
    explicit_m = {f'{e}_small_09' for e in ENDPOINTS}
    m = {}
    for r in updated['case_results']:
        k = r['case_id']
        require(k in explicit_m and k not in m and k in states and r['record_id']=='optoGEF2_small_col_M'
                and all(r[f]==states[k][f] for f in ('record_id','field','field_sha256'))
                and Q(r['exact_requested_drive'])==Q(states[k]['requested_drive_exact'])
                and 0<=Q(r['posterior_error_upper_exact'])<Q(r['radius_exact'])<Q(r['observer_goal_exact'])
                ==Q(states[k]['observer_goal_exact']),
                'M proof cannot exclude a different saved field or target')
        proof=ROOT/r['native_proof_file']
        require(proof.resolve().is_relative_to(INVENTORY.resolve()) and sha(proof)==r['native_proof_sha256'],
                'verified M case proof bytes changed')
        m[k] = r['native_proof_sha256']
    require(len(original)==4 and set(m)==explicit_m and not (set(original)&set(m)),
            'four prior pilot and two distinct new cases required')
    existing={**original,**m}
    pending={k for k in states if k not in existing}
    waves={}
    seen=[]
    require(len(ledger['recommended_batch_waves'])==len(WAVE_IDS),'unexpected number of frozen waves')
    for wave_id,group in zip(WAVE_IDS,ledger['recommended_batch_waves'],strict=True):
        rows=[]
        for case_id in group['case_ids']:
            require(case_id in states and case_id not in original,'foreign or completed original pilot in schedule')
            if case_id in pending:
                row=states[case_id]
                require(row['role']=='small_record' and row['record_id'] in group['matched_record_ids'],
                        'wave entries must belong to its declared matched small records')
                rows.append(row)
        require(len(rows) in (2,4),'wave must contain one or two complete matched pairs')
        for record in {r['record_id'] for r in rows}:
            pair=[r for r in rows if r['record_id']==record]
            require(tuple(r['endpoint_id'] for r in pair)==ENDPOINTS
                    and len({(r['checkpoint'],Q(r['dose_exact'])) for r in pair})==1,
                    'pair must retain both frozen endpoints at one exact optical dose')
        waves[wave_id]=rows
        seen.extend(r['case_id'] for r in rows)
    require(len(seen)==len(set(seen))==len(pending)==30 and set(seen)==pending,
            '30 unproved cases must occur once in the fixed waves')
    return waves,existing


def case_inputs(ledger, authorities):
    p,m,_,b,bn,_=authorities
    ep=unique(p['endpoints'],lambda r:r['endpoint_id']==ledger['endpoint_id'])
    def at(row):
        return row['endpoint_id']==ledger['endpoint_id'] and row['checkpoint']==ledger['checkpoint']
    return ep,unique(m['fields'],at),unique(b['field_budgets'],at),unique(bn['fields'],at),ledger


def case_admission(case_id, ep, field, row, independent, ledger, state, verification, delta):
    endpoint,checkpoint = ledger['endpoint_id'],ledger['checkpoint']
    record = row['record_id']
    require(endpoint in ENDPOINTS and checkpoint in {f'small_{i:02d}' for i in range(1,18)}
            and ep['endpoint_id'] == endpoint and ledger['case_id'] == case_id == f'{endpoint}_{checkpoint}'
            and all(r['endpoint_id'] == endpoint and r['checkpoint'] == checkpoint
                    for r in (field, row, independent, ledger))
            and all(r['role'] == 'small_record' for r in (field, row, ledger))
            and independent['record_id'] == ledger['record_id'] == record
            and ledger['fresh_cache_and_preconditioner_needed'] is True
            and ledger['pilot_cached_and_geometry_coverage'] is None,
            'mixed candidate/checkpoint or reused cache authority')
    for key in ('field', 'field_sha256', 'verification', 'verification_sha256'):
        require(field[key] == row[key] == ledger[key], 'saved field or verification identity changed')
    require(independent['field'] == field['field'] and independent['field_sha256'] == field['field_sha256']
            and field['field'] == f'{relative(RUN)}/actual/{endpoint}/{checkpoint}/field.json'
            and field['verification'] == f'{relative(RUN)}/actual/{endpoint}/{checkpoint}/verification.json',
            'wrong independent field binding')
    step = ep['small_records'][int(checkpoint[-2:])-1]
    dose = Q(step['dose_exact'])
    require(step['record_id']==record and dose>0,'frozen record and positive dose required')
    material = ep['materials']['small']
    center = row['center_and_model']
    require(step == row['step'] == field['step'] and Q(step['dose_exact']) == Q(ledger['dose_exact']) == dose,
            'different exact physical dose or step')
    params = material['parameters']
    require(params == state['parameters'] == center['material_parameters']
            and material['parameter_storage_exact'] == center['material_parameter_storage_exact']
            and Q(params['mu']) == 1 and Q(params['bulk']) == 5000,
            'saved material parameters differ from the frozen independently reviewed row')
    gain = Q(material['gain_exact'])
    drive = gain*dose
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
    return {'case_id': case_id, 'endpoint_id': endpoint, 'checkpoint': checkpoint, 'record_id': record,
            'field': field['field'], 'field_sha256': field['field_sha256'],
            'verification': field['verification'], 'verification_sha256': field['verification_sha256'],
            'step': step, 'material': material, 'free_map': free, 'active_elements': active,
            'xi_dimensionless_exact': str(xi), 'e_goal_dimensionless_exact': str(goal),
            'closed_parameter_cell_exact': {k:[str(x),str(x)] for k,x in phys.items()},
            'parameter_order': list(ORDER), 'observer_variant': OBSERVER,
            'frozen_conditional_allowance_um_exact': str(delta),
            'drive': {'optical_dose_exact': str(dose), 'common_gain_exact': str(gain),
                      'activation_strength_exact': str(drive), 'activation_strength_binary_exact': str(rounded),
                      'binary_input_minus_exact_strength': str(rounded-drive)}}


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


def prepare_protocol():
    for path,digest in PINS.items():
        require(sha(path)==digest,f'frozen input changed: {relative(path)}')
    authorities=authority_inputs()
    p,_,native,_,_,ledger=authorities
    core_review=read(CORE_REVIEW,PINS[CORE_REVIEW])
    require(core_review['status']=='PASS' and all_true(core_review['checks'])
            and core_review['source_sha256']==PINS[CORE_SOURCE],'cached algebra review changed')
    packages={name:importlib.metadata.version(name) for name in ('numpy','scipy','torch','python-flint')}
    require(all(packages[k]==p['runtime'][k] for k in ('numpy','scipy','torch')),
            'same core binary package versions required')
    sources={**p['source_hashes'],**{relative(path):h for path,h in PINS.items() if path.suffix=='.py'}}
    for name,h in sources.items():
        require(sha(ROOT/name)==h,f'pinned source closure changed: {name}')
    inputs={relative(path):h for path,h in PINS.items() if path.suffix!='.py'}
    waves,existing=schedule(ledger,read(PRIOR_FOUR,PINS[PRIOR_FOUR]),read(VERIFIED_M,PINS[VERIFIED_M]))
    cases={}
    delta=Q(p['selection']['conditional_source_allowance_um_exact'])
    for wave,rows in waves.items():
        for row in rows:
            case_id=row['case_id']
            records=case_inputs(row,authorities)
            field=records[1]
            payload=[]
            for k in ('field','verification'):
                path,h=field[k],field[k+'_sha256']
                require(native['all_hashes'][path]==ledger['readback_sha256'][path]==h,
                        'selected payload lacks independent byte authority')
                payload.append(read(ROOT/path,h))
                inputs[path]=h
            cases[case_id]=case_admission(case_id,*records,*payload,delta)
    require(len(cases)==len({r['field_sha256'] for r in cases.values()})==30,'duplicate saved states')
    return {'status':'FROZEN_REMAINING30_SMALL_OUTCOME_CACHED_ROOT_TARGETS','scope':SCOPE,
            'source_sha256':sha(SOURCE),'input_sha256':dict(sorted(inputs.items())),
            'source_closure_sha256':dict(sorted(sources.items())),
            'waves':{wave:[r['case_id'] for r in rows] for wave,rows in waves.items()},'cases':cases,
            'covered_cached_cases_excluded':existing,
            'numerical_coverage':{'pending_cached_cases':30,'associated_primary_small_inequalities':120,
                                  'exact_spatial_coverage_still_missing_outside_original_four':32,
                                  'original_four_exact_spatial_roots_assumed':False,
                                  'automatic_observer_transfer':False},
            'numerical_policy':p['numerical_policy'],'runtime_packages':packages,
            'precision_bits':PRECISION,'torch_cpu_threads':1,'site_counts':TERMS,
            'finite_radius_ladder_exact':list(map(str,RADII)),
            'radius_execution':('Skip any radius at or above the independently reviewed case-specific observer goal. '
                                'Try in descending order until strict cached acceptance or exhaustion. '
                                'After the first acceptance, retain one next smaller declared radius as a separately attempted replay backup, if available.'),
            'exact_cached_force_drive':'Exact immutable gain times exact dose. The nearest binary drive is used only for the fresh numerical inverse proposal.',
            'strict_cached_decision':'Positive analytic domain, r < observer goal, z < 1, eta + z*r < r, and eta/(1-z) < observer goal.',
            'execution':'Each disjoint fixed wave requires its own external prelaunch review. One worker at a time under that wave supervisor.',
            'requires_independent_actual_review_for_every_current_point':True,
            'new_equilibrium_solves':0,'prior_case_cache_or_R_reused':False}


def run_radius_ladder(goal,evaluate):
    require(type(goal) is Q and goal>0,'positive exact independent observer goal required')
    attempts=[]
    primary=backup=None
    for index,radius in enumerate(RADII):
        if radius>=goal:
            attempts.append({'radius_requested_exact':str(radius),'status':'SKIPPED_OUTSIDE_OBSERVER_GOAL',
                             'strict_selfmap_contraction_and_observer_radius':False})
            continue
        value=evaluate(index,radius)
        require(Q(value['radius_requested_exact'])==radius
                and type(value['strict_selfmap_contraction_and_observer_radius']) is bool,
                'different radius or nonliteral decision returned from cached bound')
        attempts.append(value)
        if primary is not None:
            backup=value
            break
        if value['strict_selfmap_contraction_and_observer_radius'] is True:
            primary=value
    return attempts,primary,backup


def wave_document(protocol,wave_id,protocol_sha):
    require(wave_id in WAVE_IDS,'unknown frozen wave')
    ids=protocol['waves'][wave_id]
    require(len(ids)==len(set(ids)) and len(ids) in (2,4),'complete disjoint wave required')
    fields=('case_id','record_id','checkpoint','endpoint_id','field','field_sha256','verification_sha256',
            'e_goal_dimensionless_exact','drive')
    return {'status':'FROZEN_REMAINING30_CACHED_WAVE_MANIFEST','wave_id':wave_id,'source_sha256':sha(SOURCE),
            'protocol_sha256':protocol_sha,'ordered_case_ids':ids,
            'case_anchors':[{k:protocol['cases'][i][k] for k in fields} for i in ids],
            'producer_actual_work_in_preparation':0}


def admit_wave_review(review,wave_id,bindings,ids):
    require(review['status']==REVIEW_STATUS and all_true(review['checks'])
            and review['wave_id']==wave_id and review['ordered_case_ids']==ids
            and all(review[k]==v for k,v in bindings.items()),
            'separate independently reviewed prelaunch for this exact frozen wave required')


def fixtures(protocol):
    authorities=authority_inputs()
    prior=read(PRIOR_FOUR,PINS[PRIOR_FOUR])
    m=read(VERIFIED_M,PINS[VERIFIED_M])
    ledger=authorities[-1]
    waves,_=schedule(ledger,prior,m)
    first=waves[WAVE_IDS[0]][0]
    second=waves[WAVE_IDS[0]][1]
    args=list(case_inputs(first,authorities))
    field=args[1]
    args.extend(read(ROOT/field[k],field[k+'_sha256']) for k in ('field','verification'))
    args.append(Q(protocol['cases'][first['case_id']]['frozen_conditional_allowance_um_exact']))
    def altered(index,key,value):
        fixture=deepcopy(args)
        fixture[index][key]=value
        return rejects(lambda:case_admission(first['case_id'],*fixture))
    def radius_fixture(_,r):
        return {'radius_requested_exact':str(r),
                'strict_selfmap_contraction_and_observer_radius':r==Q(1,2**28)}
    attempts,accepted,backup=run_radius_ladder(Q(1,2**24),radius_fixture)
    bad=deepcopy(ledger)
    bad['recommended_batch_waves'][1]['case_ids'][0]=first['case_id']
    wrong_m=deepcopy(m)
    wrong_m['case_results'][0]['field_sha256']='0'*64
    execution_id={'case_id':first['case_id'],'pid':8,'supervisor_pid':7}
    terminal={**execution_id,'status':'FINISHED'}
    admit_cached_worker(0,terminal,execution_id,execution_id)
    wave_id=WAVE_IDS[0]
    review_hashes={'wave_manifest_sha256':'a'*64}
    approved={'status':REVIEW_STATUS,'checks':{'fixture':True},'wave_id':wave_id,
              'ordered_case_ids':protocol['waves'][wave_id],**review_hashes}
    admit_wave_review(approved,wave_id,review_hashes,protocol['waves'][wave_id])
    return {'unique_30_inputs_with_wave_sizes_2_plus_seven_times_4':
                [len(v) for v in protocol['waves'].values()]==[2]+[4]*7
                and set(protocol['cases']).isdisjoint(protocol['covered_cached_cases_excluded']),
            'reject_duplicated_or_substituted_wave_case':rejects(lambda:schedule(bad,prior,m)),
            'reject_m_exclusion_with_mismatched_actual_field':rejects(lambda:schedule(ledger,prior,wrong_m)),
            'reject_cross_record_step_or_field':(rejects(lambda:case_admission(second['case_id'],*args))
                and altered(1,'field_sha256',second['field_sha256']) and altered(4,'dose_exact','1')),
            'reject_rounded_strength_and_changed_independent_budget_or_operator':(
                altered(4,'requested_drive_exact',first['solver_binary_drive_exact'])
                and altered(3,'e_goal_dimensionless_exact','0')
                and altered(6,'variant','current_disk__regional_mean')),
            'strict_goal_rejects_equal_radius_and_keeps_one_narrower_backup':(
                len(attempts)==4 and [x['status'] for x in attempts[:2]]==['SKIPPED_OUTSIDE_OBSERVER_GOAL']*2
                and Q(accepted['radius_requested_exact'])==Q(1,2**28)
                and Q(backup['radius_requested_exact'])==Q(1,2**32)),
            'explicit_smaller_2_to_minus_34_fallback_predeclared':Q(1,2**34) in RADII,
            'reject_nonliteral_process_exit_or_mixed_case_worker':(
                rejects(lambda:admit_cached_worker(False,terminal,execution_id,execution_id))
                and rejects(lambda:admit_cached_worker(1,terminal,execution_id,execution_id))
                and rejects(lambda:admit_cached_worker(0,{**terminal,'case_id':second['case_id']},execution_id,execution_id))),
            'reject_other_wave_review_or_changed_wave_manifest_hash':(
                rejects(lambda:admit_wave_review(approved,WAVE_IDS[1],review_hashes,protocol['waves'][WAVE_IDS[1]]))
                and rejects(lambda:admit_wave_review(approved,wave_id,{'wave_manifest_sha256':'b'*64},protocol['waves'][wave_id]))),
            'no_actual_outputs_created':not ACTUAL.exists()}


def prepare():
    PREP.mkdir(parents=True,exist_ok=False)
    with (PREP/FROZEN_SOURCE).open('xb') as handle:
        handle.write(SOURCE.read_bytes())
    try:
        protocol=prepare_protocol()
        checks=fixtures(protocol)
        require(all_true(checks),'cold remaining-30 preparation failed')
        ph=write(PREP/'protocol.json',protocol)
        files=[]
        for wid in WAVE_IDS:
            path=PREP/(wid+'.json')
            wave=wave_document(protocol,wid,ph)
            h=write(path,wave)
            files.append({'wave_id':wid,'file':path.name,'sha256':h,'ordered_case_ids':wave['ordered_case_ids']})
        ih=write(PREP/'wave_index.json',{'status':'FROZEN_EIGHT_DISJOINT_CACHED_WAVES',
                 'source_sha256':sha(SOURCE),'protocol_sha256':ph,'waves':files,'total_unique_cases':30})
        pf=write(PREP/'preflight.json',{'status':'PASS','source_sha256':sha(SOURCE),'protocol_sha256':ph,
                 'wave_index_sha256':ih,'checks':checks,'models_constructed':0,'preconditioners_recovered':0,
                 'binary_caches_frozen':0,'actual_AD_sites_assembled':0,'new_equilibrium_solves':0})
        print(json.dumps({'source_sha256':sha(SOURCE),'protocol_sha256':ph,'wave_index_sha256':ih,
                          'preflight_sha256':pf,'waves':files}),flush=True)
    except BaseException as err:
        write(PREP/'failure.json',{'error':str(err),'traceback':traceback.format_exc()})
        raise


def prelaunch(wave_id,path,digest):
    require(wave_id in WAVE_IDS,'unknown frozen wave')
    r=read(path,digest)
    ph=sha(PREP/'protocol.json')
    ih=sha(PREP/'wave_index.json')
    wave_file=PREP/(wave_id+'.json')
    wh=sha(wave_file)
    bindings={'source_sha256':sha(SOURCE),'protocol_sha256':ph,'preflight_sha256':sha(PREP/'preflight.json'),
              'wave_index_sha256':ih,'wave_manifest_sha256':wh,'coverage_ledger_sha256':PINS[LEDGER],
              'verified_two_M_coverage_sha256':PINS[VERIFIED_M],'first_four_cached_coverage_sha256':PINS[PRIOR_FOUR]}
    wave=read(wave_file,wh)
    admit_wave_review(r,wave_id,bindings,wave['ordered_case_ids'])
    pf=read(PREP/'preflight.json',bindings['preflight_sha256'])
    require(pf['status']=='PASS' and all_true(pf['checks'])
            and pf['source_sha256']==bindings['source_sha256'] and pf['protocol_sha256']==ph and pf['wave_index_sha256']==ih
            and sha(PREP/FROZEN_SOURCE)==bindings['source_sha256']
            and all(type(pf[k]) is int and pf[k]==0 for k in ('models_constructed','preconditioners_recovered',
                'binary_caches_frozen','actual_AD_sites_assembled','new_equilibrium_solves')),
            'current cold source and preflight mismatch')
    protocol=read(PREP/'protocol.json',ph)
    require(protocol==prepare_protocol(),'current input identities differ from frozen 30-outcome protocol')
    index=read(PREP/'wave_index.json',ih)
    require(index['status']=='FROZEN_EIGHT_DISJOINT_CACHED_WAVES'
            and index['source_sha256']==bindings['source_sha256'] and index['protocol_sha256']==ph
            and type(index['total_unique_cases']) is int and index['total_unique_cases']==30
            and [x['wave_id'] for x in index['waves']]==list(WAVE_IDS)
            and all(x['ordered_case_ids']==protocol['waves'][x['wave_id']]
                    and x['file']==x['wave_id']+'.json' and sha(PREP/x['file'])==x['sha256'] for x in index['waves']),
            'whole immutable wave index is required')
    require(wave==wave_document(protocol,wave_id,ph),'review wave IDs or target anchors differ')
    return protocol,wave,{**bindings,'independent_wave_review_sha256':digest}


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
            and cp.NFREE==NFREE and cp.TERMS==TERMS,
            'local model/cache kernel scope changed')
    torch.set_num_threads(1)
    core.ctx.prec = PRECISION
    core.A.cache_clear()
    require(torch.get_num_threads()==1 and core.ctx.prec==PRECISION,'floating package policy changed')
    model,u = cp.model_input(case,protocol['numerical_policy'])
    drive = case['drive']
    strength = Q(drive['activation_strength_exact'])
    require(drive==cp.exact_drive(case['step'],case['material'])
            and strength==Q(drive['optical_dose_exact'])*Q(case['closed_parameter_cell_exact']['common_gain'][0]),
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
    goal=Q(case['e_goal_dimensionless_exact'])
    def evaluate_radius(i,r):
        require(0<r<goal,'positive radius strictly below this case observer goal required')
        test=cp.test_radius(core,sites,abs_R,eta,z0,r,goal)
        filename=f'radius_{i:02d}.json'
        h=write(folder/filename,{'case_id':case['case_id'],'cache_packet_sha256':ph,**test})
        print(json.dumps({'case_id':case['case_id'],'radius':str(r),'status':test['status']}),flush=True)
        return {'path':filename,'sha256':h,**test}
    attempts,accepted,backup=run_radius_ladder(goal,evaluate_radius)
    return {'status':'CANDIDATE_CACHED_ROOT_AWAITING_INDEPENDENT_REVIEW' if accepted else
            'UNRESOLVED_CACHED_CONTRACTION_OR_BUDGET','scope':SCOPE,'case':case,'bindings':bindings,
            'cache_packet':relative(path),'cache_packet_sha256':ph,'embedded_binary_cache_sha256':payload['sha256'],
            'point_enclosure_sha256':enclosure,'eta_upper_exact':str(eta),'center_Z0_upper_exact':str(z0),
            'radius_attempts':attempts,'accepted_radius':accepted,'narrower_predeclared_replay_backup':backup,'backend':{
                'arb_precision_bits':int(core.ctx.prec),'torch_threads':torch.get_num_threads(),'packages':protocol['runtime_packages']},
            'new_nonlinear_equilibrium_solves':0,'exact_spatial_integration_established':False,
            'whole_body_geometry_established':False,'local_minimum_established':False,'independent_actual_review_required':True}


def worker(wave_id,review,digest,case_id,parent_pid):
    protocol,wave,bindings=prelaunch(wave_id,review,digest)
    parent=ACTUAL/wave_id
    opening=json.loads((parent/'started.json').read_bytes())
    require(type(parent_pid) is int and os.getppid()==parent_pid==opening['pid'] and opening['bindings']==bindings
            and opening['wave_id']==wave_id and opening['ordered_case_ids']==wave['ordered_case_ids']
            and case_id in wave['ordered_case_ids'],'current live serial wave parent and exact approved case required')
    case=protocol['cases'][case_id]
    folder=parent/case_id
    identity={'pid':os.getpid(),'supervisor_pid':parent_pid,'wave_id':wave_id,'case_id':case_id,'bindings':bindings}
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


def launch(wave_id,review,digest):
    protocol,wave,bindings=prelaunch(wave_id,review,digest)
    folder=ACTUAL/wave_id
    folder.mkdir(parents=True,exist_ok=False)
    pid=os.getpid()
    mark={'pid':pid,'wave_id':wave_id,'bindings':bindings,'ordered_case_ids':wave['ordered_case_ids']}
    write(folder/'started.json',mark)
    terminal={**mark,'status':'FAILED','workers':[]}
    results=[]
    try:
        for case_id in wave['ordered_case_ids']:
            case=protocol['cases'][case_id]
            target=folder/case_id
            target.mkdir(exist_ok=False)
            command=[sys.executable,'-B',str(SOURCE),'worker','--wave-id',wave_id,'--review',str(review),'--review-sha256',digest,
                     '--case-id',case_id,'--supervisor-pid',str(pid)]
            with (target/'stdout.log').open('x',encoding='utf-8') as log:
                child=subprocess.Popen(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
                identity={'pid':child.pid,'supervisor_pid':pid,'wave_id':wave_id,'case_id':case_id,'bindings':bindings}
                write(target/'child_started.json',{**identity,'command':command})
                code=child.wait()
            path=target/'terminal.json'
            observed={**identity,'observed_exit_code':code,'terminal_sha256':sha(path) if path.is_file() else None}
            oh=write(target/'observed_os_exit.json',observed)
            terminal['workers'].append({**observed,'os_exit_receipt_sha256':oh})
            require(type(code) is int and code==0 and observed['terminal_sha256'] is not None,
                    'cached worker failed before a complete attempt')
            done=read(path,observed['terminal_sha256'])
            started=json.loads((target/'worker_started.json').read_bytes())
            admit_cached_worker(code,done,started,identity)
            result=read(target/'results.json',done['results_sha256'])
            require(result['case']==case and result['bindings']==bindings
                    and result['status'] in ('CANDIDATE_CACHED_ROOT_AWAITING_INDEPENDENT_REVIEW',
                                             'UNRESOLVED_CACHED_CONTRACTION_OR_BUDGET'),
                    'mixed cached case result')
            packet=target/'cache_packet.json'
            require(result['cache_packet']==relative(packet) and sha(packet)==result['cache_packet_sha256'],
                    'cache packet bytes or ownership changed')
            results.append({'case_id':case_id,'status':result['status'],'results_sha256':done['results_sha256'],
                            'cache_packet':relative(packet),'cache_packet_sha256':result['cache_packet_sha256']})
        require([x['case_id'] for x in results]==wave['ordered_case_ids'],'incomplete cached wave')
        terminal['results_sha256']=write(folder/'results.json',{'status':'COMPLETE_CACHED_WAVE_ATTEMPTS_UNREVIEWED',
                'wave_id':wave_id,'scope':SCOPE,'wave_manifest_sha256':bindings['wave_manifest_sha256'],
                'cases':results,'bindings':bindings,'independent_actual_review_required':True,
                'certified_cached_roots_by_this_supervisor':0})
        terminal['status']='FINISHED'
    except BaseException as error:
        terminal['error']={'type':type(error).__name__,'message':str(error)}
        raise
    finally:
        write(folder/'terminal.json',terminal)
        print(json.dumps(terminal),flush=True)


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    sub=cli.add_subparsers(dest='stage',required=True)
    sub.add_parser('prepare')
    for name in ('launch','worker'):
        parser=sub.add_parser(name)
        parser.add_argument('--wave-id',required=True,choices=WAVE_IDS)
        parser.add_argument('--review',required=True,type=lambda value:Path(value).resolve())
        parser.add_argument('--review-sha256',required=True)
        if name=='worker':
            parser.add_argument('--case-id',required=True)
            parser.add_argument('--supervisor-pid',required=True,type=int)
    args=cli.parse_args()
    if args.stage=='prepare':
        prepare()
    elif args.stage=='launch':
        launch(args.wave_id,args.review,args.review_sha256)
    else:
        worker(args.wave_id,args.review,args.review_sha256,args.case_id,args.supervisor_pid)


