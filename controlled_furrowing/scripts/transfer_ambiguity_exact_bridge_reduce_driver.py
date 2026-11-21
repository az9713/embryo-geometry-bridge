"""Admit and reduce one reviewed singleton spatial bridge after complete replay."""

import argparse
from copy import deepcopy
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
from itertools import combinations_with_replacement
import json
from math import isfinite
import os
from pathlib import Path
import re
import sys
import traceback


SOURCE = Path(__file__).resolve()
ROOT = SOURCE.parents[1]
EXACT = ROOT/'runs/transfer_ambiguity_nonlinear_probe_v1/certification_design/exact_targets_v1'
BRIDGE = EXACT/'bridge_v1'
REPLAY = BRIDGE/'independent_actual_preparation'
OUT = EXACT/'bridge_reduction_v1'
PREP = OUT/'caller_preparation'
ACTUAL = OUT/'actual'
KERNEL_SOURCE = ROOT/'scripts/transfer_ambiguity_exact_bridge_reduction.py'
KERNEL_REVIEW = OUT/'independent_review/final_binding.json'
BRIDGE_SOURCE = ROOT/'scripts/transfer_ambiguity_exact_bridge.py'
BRIDGE_PROTOCOL = BRIDGE/'preparation_v2/protocol.json'
BRIDGE_PREFLIGHT = BRIDGE/'preparation_v2/preflight.json'
REPLAY_SOURCE = REPLAY/'replay.py'
REPLAY_PROTOCOL = REPLAY/'protocol.json'
REPLAY_PREFLIGHT = REPLAY/'preflight.json'
DERIV_ADMISSION = EXACT/'cached_derivative_extension_v1/independent_actual_preparation/bridge_admission.json'
PINS = {
    KERNEL_SOURCE:'d31ca58afe2e1b935d723513fb9fa18f2f978885f8f6ec32e48766f8b15bf0e8',
    KERNEL_REVIEW:'d6e71bcc68e6293fda36e04c3f2b506a7c34c8516492beeda699ad97c433d4e2',
    BRIDGE_SOURCE:'58ae92f849414f185c4e28bbb2ec3b9d05361391c4d6e3fa71f519162e30e5e7',
    BRIDGE_PROTOCOL:'4d38e999565b452a276777b66690d9b8ed4df8f3576cd8ed26d8757b48900581',
    BRIDGE_PREFLIGHT:'8d7acca719a3f9466abc828c2732c5cb825a4012d2a7c1c30445cd306ba91484',
    REPLAY_SOURCE:'e504a67f913968888805c0385c6b3fa0501d6b3e7f5345362b453eb1dde50074',
    REPLAY/'replay_math.py':'e7bb66066d54e131fd009d2bf48e94f3f1fb7bc1955023f66bc96cbe3c6d78e9',
    REPLAY_PROTOCOL:'2ec903e8088596b72d3f59e57ce6700247218b2138fbe9e78edeb1cb046cb9a1',
    REPLAY_PREFLIGHT:'effd86b721146884ff4561e6bd24eabe3f32a260d768fd09218ef0432eb5442b',
    REPLAY/'ROOT-SOURCE-READBACK.md':'60b665f94b38d19a1565be88af9a5941435052c6c063fd8f36e178d98777ebf6',
    DERIV_ADMISSION:'269229312793ec101cc301086d634dd4de89c62a5050ccfc60f185cf761dbc7d',
}
PILOTS = tuple(f'{p}_{s}' for p in ('max_contrast_change','max_negative_contrast_change')
               for s in ('small_17','large_03'))
N = 611
PRECISION = 640
RADIUS = Q(1,2**34)
POINT = 'EXACT_CENTER_POINT_ONLY'
BOX = 'ENTIRE_CLOSED_FIXED_PARAMETER_CUBE'
ORDER = tuple((k,i) for k,n in (('bulk',144),('surface',24)) for i in range(n))
NATIVE_STATUS = 'PASS_INDEPENDENT_ACTUAL_SINGLETON_BRIDGE_COMPONENT_REPLAY'
OS_STATUS = 'PASS_ROOT_OBSERVED_SINGLETON_BRIDGE_REPLAY_TERMINAL'
REVIEW_STATUS = 'PASS_INDEPENDENT_SINGLETON_BRIDGE_REDUCTION_CALLER_PRELAUNCH'
NATIVE_CHECKS = {'true_os_complete_producer_and_both_child_terminals',
                 'same_exact_point_and_whole_box_all_168_bindings',
                 'all_point_cache_force_reconstructed_and_contained',
                 'all_point_force_and_whole_H_contained'}
FAILURE_KEYS = {'direct_point_cache','point_force','whole_box_H'}
SCOPE = ('Single admitted fixed-material radius-2^-34 exact-spatial root arithmetic candidate. '
         'The point full force uses signed cached and spatial-defect rows before norms. '
         'The derivative bound applies on the entire matched embedded cube. '
         'Requires independent actual reduction review. No minimum, source or continuum claim.')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def rel(path):
    return path.relative_to(ROOT).as_posix()


def checkout_path(name):
    require(type(name) is str, 'relative path string required')
    path = (ROOT/name).resolve()
    require(path.is_relative_to(ROOT) and rel(path) == name, 'noncanonical study path')
    return path


def digest_text(value):
    return type(value) is str and re.fullmatch('[0-9a-f]{64}',value) is not None


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def pinned(path, digest):
    require(digest_text(digest), f'missing verified SHA-256 for {path}')
    raw = path.read_bytes()
    require(sha256(raw).hexdigest() == digest, f'changed authority file: {path}')
    return json.loads(raw)


def save(path, value):
    with path.open('x',encoding='utf-8',newline='\n') as f:
        json.dump(value,f,indent=2,allow_nan=False)
        f.write('\n')
    return sha(path)


def yes(checks):
    return type(checks) is dict and bool(checks) and all(x is True for x in checks.values())


def positive_int(x):
    return type(x) is int and x > 0


def true_zero(x):
    return type(x) is int and x == 0


def rational(x):
    require(type(x) in (Q,int,str),'finite explicit rational required')
    return Q(x)


def interval(x):
    require(type(x) in (tuple,list) and len(x)==2,'two exact endpoints required')
    lo,hi=map(rational,x)
    require(lo<=hi,'reversed exact interval')
    return lo,hi


def authority_protocol():
    paths = {rel(p):h for p,h in PINS.items()}
    for path,h in PINS.items():
        require(sha(path) == h, f'changed cold source or proof: {rel(path)}')
    root_review = pinned(KERNEL_REVIEW,PINS[KERNEL_REVIEW])
    require(root_review['status']=='PASS' and yes(root_review['checks'])
            and root_review['applies_only_to_pure_kernel'] is True
            and root_review['actual_packet_authority_admitted'] is False
            and root_review['bound_sha256'][rel(KERNEL_SOURCE)]==PINS[KERNEL_SOURCE],
            'current independent pure reduction review required')
    native = pinned(REPLAY_PROTOCOL,PINS[REPLAY_PROTOCOL])
    pf = pinned(REPLAY_PREFLIGHT,PINS[REPLAY_PREFLIGHT])
    require(native['status']=='FROZEN_INDEPENDENT_SINGLETON_COMPONENT_REPLAY'
            and native['source_sha256']==PINS[REPLAY_SOURCE]
            and native['math_source_sha256']==PINS[REPLAY/'replay_math.py']
            and native['precision_bits']==PRECISION and native['counts']=={'bulk':144,'surface':24}
            and rational(native['radius_exact'])==RADIUS
            and pf['status']=='PASS' and pf['source_sha256']==PINS[REPLAY_SOURCE]
            and pf['protocol_sha256']==PINS[REPLAY_PROTOCOL]
            and true_zero(pf['actual_packet_payloads_opened']) and true_zero(pf['actual_moment_force_or_H_work']),
            'current native replay source and preparation required')
    bridge = pinned(BRIDGE_PROTOCOL,PINS[BRIDGE_PROTOCOL])
    bp = pinned(BRIDGE_PREFLIGHT,PINS[BRIDGE_PREFLIGHT])
    require(bridge['source_sha256']==PINS[BRIDGE_SOURCE]
            and bridge['point_scope']==POINT and bridge['whole_hessian_scope']==BOX
            and bridge['precision_bits']==PRECISION and rational(bridge['whole_H_radius_exact'])==RADIUS
            and tuple(bridge['pilot_order'])==PILOTS and bp['status']=='PASS' and yes(bp['checks']),
            'same bound producer phase and scope required')
    for source in (root_review['bound_sha256'],native['source_and_premise_hashes']):
        for name,h in source.items():
            require(paths.get(name,h)==h and sha(checkout_path(name))==h,f'changed source or prerequisite: {name}')
            paths[name]=h
    for name,h in bridge['input_sha256'].items():
        # The fixed producer protocol pins case payloads. Only source bytes are needed for cold preparation.
        require(paths.get(name,h)==h and digest_text(h),'mixed producer prerequisite pins')
        if name.endswith('.py'):
            require(sha(checkout_path(name))==h,f'changed producer source: {name}')
        paths[name]=h
    return {'status':'FROZEN_SINGLETON_BRIDGE_REDUCTION_CALLER','scope':SCOPE,
            'source_sha256':sha(SOURCE),'source_and_prerequisite_sha256':dict(sorted(paths.items())),
            'pilot_choices':list(PILOTS),'precision_bits':PRECISION,'radius_exact':str(RADIUS),
            'point_force_scope':POINT,'whole_hessian_scope':BOX,'packet_counts':{'bulk':144,'surface':24},
            'native_status_required':NATIVE_STATUS,'native_checks_required':sorted(NATIVE_CHECKS),
            'root_observed_native_OS_status_required':OS_STATUS,
            'required_external_admission_sha256s':['caller_prelaunch_review','original_root_producer_admission',
                                                  'root_observed_native_replay_terminal'],
            'reduction_operand':'Producer signed intervals, after independently replayed point F and whole H containment.',
            'requires_independent_actual_reduction_review':True,'actual_packet_payloads_opened':0,'actual_matrix_work':0}


def native_phase(observed,native,terminal,started,pilot,admission_sha,receipt_sha,terminal_sha):
    require(pilot in PILOTS and observed['pilot_id']==native['pilot_id']==pilot
            and observed['status']==OS_STATUS and positive_int(observed['unified_exec_session'])
            and true_zero(observed['observed_exit_code']) and positive_int(observed['pid'])
            and observed['terminal_sha256']==terminal_sha and observed['receipt_sha256']==receipt_sha,
            'externally observed native replay OS completion required')
    phase = {'pid':observed['pid'],'source_sha256':PINS[REPLAY_SOURCE],
             'protocol_sha256':PINS[REPLAY_PROTOCOL],'preflight_sha256':PINS[REPLAY_PREFLIGHT],
             'root_observation_sha256':admission_sha}
    require(observed['root_observation_sha256']==admission_sha
            and all(all(d[k]==v for k,v in phase.items()) for d in (native,terminal,started))
            and all(positive_int(d['pid']) for d in (native,terminal,started))
            and terminal['status']=='FINISHED' and true_zero(terminal['exit_code'])
            and terminal['receipt_sha256']==receipt_sha and started['status']=='RUNNING'
            and native['status']==NATIVE_STATUS and set(native['checks'])==NATIVE_CHECKS and yes(native['checks'])
            and rational(native['radius_exact'])==RADIUS,
            'unrelated, failed, changed or incomplete native replay')


def check_order(rows):
    require(type(rows) is list and len(rows)==len(ORDER) and all(type(r['slot']) is int for r in rows)
            and tuple((r['kind'],r['slot']) for r in rows)==ORDER,'complete typed 144+24 ordering required')


def typed_map(seq):
    require(type(seq) is list and seq and all(type(x) is int and 0<=x<N for x in seq)
            and len(set(seq))==len(seq),'unique literal integer free projection required')
    return seq


def catalogues(pilot,native,plan,manifest,packet_hashes,plan_hashes):
    require(native['checked_packet_sha256']==packet_hashes and native['plan_packet_sha256']==plan_hashes
            and len(packet_hashes)==len(plan_hashes)==168,'native receipt did not bind the exact entire producer package')
    pdir = BRIDGE/'plans'/pilot
    run = BRIDGE/'actual'/pilot
    out = REPLAY/'actual_replay'/pilot
    for seq in (native['components'],plan['components'],manifest['files']):
        check_order(seq)
    entries=[]
    for row,p,f in zip(native['components'],plan['components'],manifest['files'],strict=True):
        kind,slot = row['kind'],row['slot']
        filename = f'{kind}_{slot:03d}.json'
        replay_file = out/filename
        plan_file = pdir/filename
        prod_file = run/kind/filename
        require(row['replay_file']==rel(replay_file) and digest_text(row['replay_sha256'])
                and row['packet_file']==f['file']==rel(prod_file)
                and row['packet_sha256']==f['sha256']==packet_hashes.get(rel(prod_file))
                and row['plan_file_sha256']==p['sha256']==plan_hashes.get(rel(plan_file))
                and p['file']==filename and row['point_force_scope']==POINT and row['whole_box_H_scope']==BOX,
                'mixed replay report path, producer packet, plan proof, digest or scope')
        typed_map(row['global_indices'])
        entries.append((row,prod_file,plan_file,replay_file))
    return entries


def upper_intervals(rows,m):
    order = list(combinations_with_replacement(range(m),2))
    require(type(rows) is list and len(rows)==len(order)
            and all(type(row) is list and len(row)==4 and type(row[0]) is type(row[1]) is int for row in rows)
            and [(i,j) for i,j,*_ in rows]==order,'complete ordered literal upper-H indexes required')
    return [interval((lo,hi)) for i,j,lo,hi in rows]


def contained(outer,inner):
    require(len(outer)==len(inner) and all(a<=x<=y<=b for (a,b),(x,y) in zip(outer,inner,strict=True)),
            'producer signed interval does not contain independent replay')


def component_operands(row,independent,producer,component,route,record_id,links,proof_digest):
    require(all(independent[k]==row[k] for k in ('kind','slot','global_indices','point_force_scope',
                'whole_box_H_scope','plan_file_sha256','packet_file','packet_sha256'))
            and type(independent['slot']) is type(producer['slot']) is type(component['slot']) is int
            and producer['kind']==component['kind']==row['kind']
            and producer['slot']==component['slot']==row['slot']
            and producer['record_id']==record_id and producer['bindings']==links and producer['root_asserted'] is False
            and type(producer['precision_bits']) is int and producer['precision_bits']==PRECISION
            and producer['plan_file_sha256']==row['plan_file_sha256']
            and producer['proof_sha256']==proof_digest,
            'review, route and component identity differ')
    ids = typed_map(row['global_indices'])
    require(all(typed_map(d['global_indices'])==ids for d in (independent,producer,component,route)),
            'independently matched free projection required')
    f,h = producer['point_force'],producer['whole_box_hessian']
    require(f['scope']==POINT and h['scope']==BOX and rational(h['radius_exact'])==RADIUS,
            'point force and whole-box Hessian scopes must remain separate')
    checks = independent['containment_failures']
    require(type(checks) is dict and set(checks)==FAILURE_KEYS and all(x==[] for x in checks.values()),
            'all three independent containment checks must pass')
    m = len(ids)
    require(all(type(independent[k]) is int for k in ('checked_force_rows','checked_upper_H_entries','checked_cache_sites'))
            and independent['checked_force_rows']==m and independent['checked_upper_H_entries']==m*(m+1)//2
            and independent['checked_cache_sites']==(64 if row['kind']=='bulk' else 12),
            'complete reported local rows and cache sites required')
    fkey = 'F_I_minus_F_cache_intervals_exact'
    hkey = 'H_I_minus_H_cache_upper_exact'
    require(len(f[fkey])==len(f['Qstar_minus_cache_correction_exact'])==m,'complete producer point force and cache correction required')
    list(map(interval,f['Qstar_minus_cache_correction_exact']))
    contained(list(map(interval,f[fkey])),list(map(interval,independent['independent_point_F_exact'])))
    contained(upper_intervals(h[hkey],m),upper_intervals(independent['independent_whole_box_H_exact'],m))
    expected = {'kind':row['kind'],'slot':row['slot'],'global_indices':ids}
    return expected,producer


def loader(path,module_name):
    require(sha(path)==PINS[path],f'changed executable source: {rel(path)}')
    spec = importlib.util.spec_from_file_location(module_name,path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    require(Path(module.__file__).resolve()==path,'resolved local source alias differs')
    return module


def review_gate(path, digest):
    binding = {'source_sha256':sha(SOURCE),'protocol_sha256':sha(PREP/'protocol.json'),
               'preflight_sha256':sha(PREP/'preflight.json'),'pure_kernel_review_sha256':PINS[KERNEL_REVIEW],
               'native_replay_source_sha256':PINS[REPLAY_SOURCE],
               'native_replay_protocol_sha256':PINS[REPLAY_PROTOCOL],'native_replay_preflight_sha256':PINS[REPLAY_PREFLIGHT]}
    r = pinned(path,digest)
    require(r['status']==REVIEW_STATUS and yes(r['checks'])
            and all(r[k]==v for k,v in binding.items()),
            'independent reduction caller prelaunch review required')
    pf = pinned(PREP/'preflight.json',binding['preflight_sha256'])
    require(pf['status']=='PASS' and yes(pf['checks'])
            and all(pf[k]==binding[k] for k in ('source_sha256','protocol_sha256'))
            and true_zero(pf['actual_packet_payloads_opened']) and true_zero(pf['actual_matrix_work'])
            and sha(PREP/'transfer_ambiguity_exact_bridge_reduce_driver_frozen.py')==binding['source_sha256'],
            'unmodified purely cold caller preflight required')
    protocol = pinned(PREP/'protocol.json',binding['protocol_sha256'])
    require(authority_protocol()==protocol,'current source and authority disagree with frozen protocol')
    return protocol,{**binding,'prelaunch_review_sha256':digest}


def authorize(args,bindings):
    pilot = args.pilot_id
    report = REPLAY/'actual_replay'/pilot
    observed = pinned(args.replay_observation,args.replay_observation_sha256)
    require(observed['status']==OS_STATUS and observed['pilot_id']==pilot
            and true_zero(observed['observed_exit_code']) and positive_int(observed['pid'])
            and positive_int(observed['unified_exec_session'])
            and observed['root_observation_sha256']==args.admission_sha256,
            'successful external replay process observation for this producer run required')
    rsha,tsha = observed['receipt_sha256'],observed['terminal_sha256']
    receipt = pinned(report/'receipt.json',rsha)
    terminal = pinned(report/'terminal.json',tsha)
    started_sha = sha(report/'started.json')
    started = pinned(report/'started.json',started_sha)
    native_phase(observed,receipt,terminal,started,pilot,args.admission_sha256,rsha,tsha)
    replay = loader(REPLAY_SOURCE,'_spatial_component_replay_admission')
    # The reviewed helper may hash complete payload bytes, but only after the native replay and exit pass.
    a,p,l,m,ps,fs = replay.admitted(args.admission,args.admission_sha256,pilot)
    require(a['producer_protocol_file']==rel(BRIDGE_PROTOCOL)
            and a['producer_preflight_file']==rel(BRIDGE_PREFLIGHT)
            and a['producer_frozen_source_file']==rel(BRIDGE/'preparation_v2/transfer_ambiguity_exact_bridge_frozen.py')
            and l['source_sha256']==receipt['producer_source_sha256']==PINS[BRIDGE_SOURCE]
            and l['protocol_sha256']==PINS[BRIDGE_PROTOCOL] and l['preflight_sha256']==PINS[BRIDGE_PREFLIGHT]
            and l['plan_sha256']==receipt['plan_sha256'] and l['plan_review_sha256']==receipt['plan_review_sha256']
            and l['derivative_admission_sha256']==receipt['derivative_admission_sha256']==PINS[DERIV_ADMISSION]
            and positive_int(a['supervisor_pid']) and positive_int(receipt['producer_observed_pid'])
            and positive_int(receipt['producer_observed_exec_session'])
            and receipt['producer_observed_pid']==a['supervisor_pid']
            and receipt['producer_observed_exec_session']==a['unified_exec_session']
            and true_zero(receipt['producer_observed_os_exit_code'])
            and receipt['producer_manifest_sha256']==a['sha256'][rel(BRIDGE/'actual'/pilot/'manifest.json')]
            and receipt['producer_terminal_sha256']==a['sha256'][rel(BRIDGE/'actual'/pilot/'terminal.json')],
            'same native-replayed producer run, independent plan and source gates required')
    producer_terminal = pinned(BRIDGE/'actual'/pilot/'terminal.json',receipt['producer_terminal_sha256'])
    require(positive_int(producer_terminal['pid']) and producer_terminal['pid']==a['supervisor_pid']
            and all(positive_int(w['pid']) and positive_int(w['supervisor_pid'])
                    and true_zero(w['observed_exit_code']) and w['supervisor_pid']==a['supervisor_pid']
                    for w in producer_terminal['workers']), 'literal producer subprocess PIDs and exits required')
    entries = catalogues(pilot,receipt,p,m,ps,fs)
    require(all(sha(rv)==e['replay_sha256'] for e,_,_,rv in entries),
            'complete independent replay report bytes must match before any component payload parse')
    return receipt,p,l,entries,{**bindings,
            'root_producer_admission_sha256':args.admission_sha256,'root_observed_replay_exit_sha256':args.replay_observation_sha256,
            'root_producer_admission_file':rel(args.admission) if args.admission.is_relative_to(ROOT) else str(args.admission),
            'root_observed_replay_exit_file':rel(args.replay_observation) if args.replay_observation.is_relative_to(ROOT) else str(args.replay_observation),
            'native_replay_receipt_sha256':rsha,'native_replay_terminal_sha256':tsha,
            'native_replay_started_sha256':started_sha,'native_replay_pid':observed['pid'],
            'producer_manifest_sha256':receipt['producer_manifest_sha256'],
            'producer_terminal_sha256':receipt['producer_terminal_sha256'],'plan_sha256':receipt['plan_sha256']}


def radius_bridge(plan,native,case,box,derivative):
    require(all(rational(v)==RADIUS for v in (plan['whole_box_radius_exact'],native['radius_exact'],
                case['new_derivative_radius_exact'],box['corrector_radius_exact']))
            and plan['point_membership']==POINT and plan['H_membership']==BOX
            and plan['derivative']==derivative,
            'whole-H radius and ZQ certificate must equal the admitted cube')
    z = rational(derivative['ZQ_upper_exact'])
    require(0<=z<1 and z==rational(native['reviewed_cached_ZQ_upper_exact']),
            'scope-matched independent plus producer cached derivative maximum required')
    return z


def context(pilot,native,plan,links,protocol):
    # The fixed source closure was authenticated by the cold protocol before these imports.
    import transfer_ambiguity_exact_bridge as bridge
    import transfer_ambiguity_exact_targets as ex
    import transfer_ambiguity_cached_derivative_extension as dx
    for name,digest in protocol['source_and_prerequisite_sha256'].items():
        if name.startswith('scripts/') and name.endswith('.py'):
            module = sys.modules.get(Path(name).stem)
            if module is not None:
                path = Path(module.__file__).resolve()
                require(path==ROOT/name and sha(path)==digest,f'reduction runtime source mismatch: {name}')
    ex.bootstrap()
    require(bridge.EX is ex and bridge.DX is dx and dx.EX is ex and ex.E.core.ctx.prec==PRECISION,
            'same reviewed factory, derivative authority and Arb precision required')
    _,cases,manifest,data = dx.authority()
    case, = [c for c in cases if c['pilot_id']==pilot]
    ctx,_,entry = ex.open_context(pilot,manifest,data)
    box = pinned(checkout_path(case['box']),case['box_sha256'])
    contract = pinned(checkout_path(entry['model_contract']),entry['model_contract_sha256'])
    derivative = bridge.derivative_inputs(DERIV_ADMISSION,PINS[DERIV_ADMISSION],cases)[pilot]
    z = radius_bridge(plan,native,case,box,derivative)
    require(box==ex.exact_node_box(ctx,RADIUS) and box['root_membership_assumed'] is False
            and case['box_sha256']==entry['box_1_sha256']==plan['context_box_sha256']==native['box_sha256']
            and bridge.digest(bridge.exact_point_inside(ctx.target.nodes(),box))==plan['point_box_sha256']
            and plan['record_id']==case['record_id']==ctx.record_id==native['record_id']==contract['record_id']
            and plan['pilot_id']==ctx.pilot_id==native['pilot_id']==pilot
            and plan['bindings']['derivative_admission_sha256']==links['derivative_admission_sha256']==PINS[DERIV_ADMISSION]
            and case['model_contract_sha256']==entry['model_contract_sha256']==ctx.row['model_observation_contract_sha256']
            and case['cache_packet_sha256']==contract['point_cache_sha256']==ctx.row['point_cache_sha256']==native['cache_packet_sha256']
            and case['embedded_cache_sha256']==contract['frozen_binary_cache_sha256']==ctx.row['frozen_binary_cache_sha256']
            and contract['parameter_cell']==ctx.row['parameter_cell']
            and contract['parameter_cell_sha256']==ctx.row['parameter_cell_sha256']==native['parameter_cell_sha256']
            and all(lo==hi==ctx.base[k] for k,(lo,hi) in ctx.cell.items())
            and ctx.drive==ctx.dose*ctx.base['common_gain']==rational(case['actual_exact_drive'])
            ==rational(contract['exact_requested_drive']['activation_strength_exact']),
            'production, replay and independent cache must bind the exact same drive, state and complete larger cube')
    proof = pinned(checkout_path(case['verified_point']),case['verified_point_sha256'])
    eta,_,_ = dx.validate_center(proof,case)
    require(native['signed_cached_anchor_file']==case['verified_point']==contract['cached_root_independent_receipt']
            and native['signed_cached_anchor_sha256']==case['verified_point_sha256']==contract['cached_root_independent_receipt_sha256'],
            'same independently reviewed point signed preconditioned cached force required')
    anchor = [interval((d['lower_exact'],d['upper_exact'])) for d in proof['signed_center_force_intervals']]
    require(len(anchor)==N and max(max(map(abs,row)) for row in anchor)<=eta,
            'complete signed cached point force enclosure required')
    expected_r = case['preconditioner_binary_sha256']
    require(expected_r==contract['binary_preconditioner_sha256']==native['binary_preconditioner_sha256']
            ==sha256(ctx.R.tobytes(order='C')).hexdigest() and ctx.R.shape==(N,N) and ctx.R.dtype.str=='<f8',
            'exact same frozen binary preconditioner required')
    routes = pinned(checkout_path(entry['routes']),entry['routes_sha256'])
    require(len(routes['bulk'])==144 and len(routes['surface'])==24,'complete factory routes required')
    return ctx,case,entry,routes,anchor,z,bridge.digest,eta


def binary_fraction(x):
    require(type(x) is float and isfinite(x),'finite binary64 scalar required')
    return Q.from_float(x)


def run(args):
    protocol,bindings = review_gate(args.review,args.review_sha256)
    native,plan,links,entries,bindings = authorize(args,bindings)
    folder = ACTUAL/args.pilot_id
    folder.mkdir(parents=True,exist_ok=False)
    mark = {'pid':os.getpid(),'pilot_id':args.pilot_id,'bindings':bindings}
    save(folder/'started.json',{**mark,'status':'RUNNING'})
    terminal = {**mark,'status':'FAILED','exit_code':1}
    try:
        ctx,case,entry,routes,anchor,z,digest,eta_q = context(args.pilot_id,native,plan,links,protocol)
        expected_maps,packets = [],[]
        for row,path,plan_path,replay_path in entries:
            kind,slot = row['kind'],row['slot']
            proof = pinned(plan_path,row['plan_file_sha256'])
            actual = pinned(path,row['packet_sha256'])
            report = pinned(replay_path,row['replay_sha256'])
            expected,packet = component_operands(row,report,actual,proof,routes[kind][slot],ctx.record_id,links,digest(proof))
            expected_maps.append(expected)
            packets.append(packet)
        kernel = loader(KERNEL_SOURCE,'_singleton_spatial_signed_reduction')
        kernel.ctx.prec = PRECISION
        require(kernel.RADIUS==RADIUS and kernel.POINT==POINT and kernel.BOX==BOX,'pure kernel scope changed')
        force,upper_H = kernel.assemble(packets,expected_maps,N)
        # NumPy list decoding preserves the exact binary64 value. Never use a decimal rendering of R.
        R = [list(map(binary_fraction,row)) for row in ctx.R.tolist()]
        operands = {'status':'ADMITTED_SIGNED_SPATIAL_REDUCTION_OPERANDS','pilot_id':args.pilot_id,'record_id':ctx.record_id,
                    'scope':SCOPE,'bindings':bindings,'radius_exact':str(RADIUS),'ZQ_upper_exact':str(z),
                    'closed_cell_exact':ctx.row['parameter_cell'],'parameter_cell_sha256':ctx.row['parameter_cell_sha256'],
                    'same_exact_drive':str(ctx.drive),'dose_exact':str(ctx.dose),'binary_R_sha256':case['preconditioner_binary_sha256'],
                    'binary_R_source':case['cache_packet'],'cache_packet_sha256':case['cache_packet_sha256'],
                    'model_contract_sha256':entry['model_contract_sha256'],'geometry_box_sha256':case['box_sha256'],
                    'cached_signed_anchor_sha256':case['verified_point_sha256'],'cached_eta_upper_exact':str(eta_q),
                    'signed_R_F_cache_exact':[list(map(str,r)) for r in anchor],
                    'signed_F_defect_exact':[list(map(str,r)) for r in force],
                    'whole_box_signed_H_defect_upper_exact':[[i,j,*map(str,v)] for (i,j),v in sorted(upper_H.items())],
                    'maps':expected_maps,'packet_sha256':native['checked_packet_sha256'],
                    'independent_component_reports_sha256':{e['replay_file']:e['replay_sha256'] for e,_,_,_ in entries}}
        oh = save(folder/'assembled_operands.json',operands)
        result = kernel.contract(R,force,upper_H,anchor,z,radius=RADIUS)
        egoal = rational(case['observer_goal_exact'])
        posterior = result['posterior_error_upper_exact']
        rh = save(folder/'results.json',{'status':'SPATIAL_REDUCTION_CANDIDATE_AWAITING_INDEPENDENT_REVIEW','scope':SCOPE,**mark,
                    'assembled_operands_sha256':oh,'arithmetic':result,'observer_design_goal_exact':str(egoal),
                    'posterior_below_observer_design_goal':posterior is not None and rational(posterior)<egoal,
                    'independent_actual_reduction_review_required':True,'exact_spatial_root_certified':False,
                    'force_H_AD_or_equilibrium_jobs_run':0})
        terminal.update(status='FINISHED',exit_code=0,results_sha256=rh,assembled_operands_sha256=oh)
        print(json.dumps({'result_sha256':rh,'classification':result['classification'],'pilot_id':args.pilot_id}),flush=True)
    except BaseException as error:
        save(folder/'failure.json',{'error':str(error),'traceback':traceback.format_exc()})
        raise
    finally:
        save(folder/'terminal.json',terminal)


def rejects(callback):
    try:
        callback()
    except ValueError:
        return True
    return False


def fixtures():
    pilot = PILOTS[0]
    base = 'a'*64
    replay_digest,terminal_digest = 'b'*64,'c'*64
    mark = {'pid':3,'source_sha256':PINS[REPLAY_SOURCE],'protocol_sha256':PINS[REPLAY_PROTOCOL],
            'preflight_sha256':PINS[REPLAY_PREFLIGHT],'root_observation_sha256':base}
    observed = {'status':OS_STATUS,'pilot_id':pilot,'observed_exit_code':0,'unified_exec_session':2,
                'receipt_sha256':replay_digest,'terminal_sha256':terminal_digest,**mark}
    native = {'status':NATIVE_STATUS,'pilot_id':pilot,'checks':dict.fromkeys(sorted(NATIVE_CHECKS),True),
              'radius_exact':str(RADIUS),**mark}
    ended = {'status':'FINISHED','exit_code':0,'receipt_sha256':replay_digest,**mark}
    started = {'status':'RUNNING',**mark}
    def phase(os_doc=observed,n=native,end=ended):
        return native_phase(os_doc,n,end,started,pilot,base,replay_digest,terminal_digest)
    phase()
    components,files,rows=[],[],[]
    ps,hs={},{}
    for k,i in ORDER:
        filename=f'{k}_{i:03d}.json'
        packet=rel(BRIDGE/'actual'/pilot/k/filename)
        proof=rel(BRIDGE/'plans'/pilot/filename)
        independent=rel(REPLAY/'actual_replay'/pilot/filename)
        ident={'kind':k,'slot':i}
        files.append({**ident,'file':packet,'sha256':base})
        components.append({**ident,'file':filename,'sha256':base})
        rows.append({**ident,'global_indices':[1],'point_force_scope':POINT,'whole_box_H_scope':BOX,
                     'packet_file':packet,'packet_sha256':base,'plan_file_sha256':base,
                     'replay_file':independent,'replay_sha256':replay_digest})
        ps[packet],hs[proof]=base,base
    catalog={**native,'checked_packet_sha256':ps,'plan_packet_sha256':hs,'components':rows}
    plan={'components':components}
    manifest={'files':files}
    def load(doc=catalog,p=plan):
        return catalogues(pilot,doc,p,manifest,ps,hs)
    malformed=deepcopy(catalog)
    malformed['components'][0]['slot']=False
    redirected=deepcopy(catalog)
    redirected['components'][0]['replay_file']=rows[-1]['replay_file']
    r=rows[0]
    links={'fixture_id':base}
    selection={'kind':'bulk','slot':0,'global_indices':[1]}
    report={**r,'checked_force_rows':1,'checked_upper_H_entries':1,'checked_cache_sites':64,
            'containment_failures':dict.fromkeys(sorted(FAILURE_KEYS),[]),
            'independent_point_F_exact':[['1','2']],'independent_whole_box_H_exact':[[0,0,'-1','1']]}
    producer={**selection,'bindings':links,'record_id':'fixture','precision_bits':640,'root_asserted':False,
              'plan_file_sha256':base,'proof_sha256':base,
              'point_force':{'scope':POINT,'F_I_minus_F_cache_intervals_exact':[['-1','3']],
                             'Qstar_minus_cache_correction_exact':[['-1','1']]},
              'whole_box_hessian':{'scope':BOX,'radius_exact':str(RADIUS),
                                   'H_I_minus_H_cache_upper_exact':[[0,0,'-2','2']]}}
    def local(n=report,p=producer,route=selection):
        return component_operands(r,n,p,selection,route,'fixture',links,base)
    overlap=deepcopy(producer)
    overlap['point_force']['F_I_minus_F_cache_intervals_exact']=[['3/2','3']]
    bad_index=deepcopy(producer)
    bad_index['whole_box_hessian']['H_I_minus_H_cache_upper_exact']=[[False,0,'-2','2']]
    derivative={'ZQ_upper_exact':'1/3'}
    radius_plan={'whole_box_radius_exact':str(RADIUS),'point_membership':POINT,'H_membership':BOX,'derivative':derivative}
    radius_native={'radius_exact':str(RADIUS),'reviewed_cached_ZQ_upper_exact':'1/3'}
    radius_case={'new_derivative_radius_exact':str(RADIUS)}
    radius_box={'corrector_radius_exact':str(RADIUS)}
    return {'complete_typed_168_catalogues_and_exact_perpacket_paths':len(load())==168,
            'reject_partial_duplicate_bool_slot_and_redirected_native_component':(
                rejects(lambda:load(p={'components':components[:-1]}))
                and rejects(lambda:load(p={'components':[components[1],*components[1:]]}))
                and rejects(lambda:load(malformed)) and rejects(lambda:load(redirected))),
            'reject_unrelated_native_status_nonzero_exit_wrong_pid_or_phase':(
                rejects(lambda:phase(n={**native,'status':'PASS_INDEPENDENT_UNRELATED'}))
                and rejects(lambda:phase(n={**native,'protocol_sha256':base}))
                and rejects(lambda:phase(os_doc={**observed,'observed_exit_code':False}))
                and rejects(lambda:phase(os_doc={**observed,'observed_exit_code':2}))
                and rejects(lambda:phase(os_doc={**observed,'pid':4}))),
            'reject_native_HOLD_or_failed_local_direct_cache_check':(
                rejects(lambda:phase(end={**ended,'status':'FINISHED_WITH_HOLD'}))
                and rejects(lambda:local(n={**report,'containment_failures':{**report['containment_failures'],'direct_point_cache':[0]}}))),
            'independent_local_signed_enclosures_must_be_contained_and_producer_arrays_are_operands':(
                local()==(selection,producer) and rejects(lambda:local(p=overlap))),
            'reject_noninteger_local_hessian_or_global_route_ids':(
                rejects(lambda:local(p=bad_index)) and rejects(lambda:local(route={**selection,'global_indices':[1.0]}))),
            'fixed_exact_radius_and_matching_ZQ_required':(
                radius_bridge(radius_plan,radius_native,radius_case,radius_box,derivative)==Q(1,3)
                and rejects(lambda:radius_bridge(radius_plan,radius_native,radius_case,
                                                 {'corrector_radius_exact':str(2*RADIUS)},derivative))
                and rejects(lambda:radius_bridge(radius_plan,{**radius_native,'reviewed_cached_ZQ_upper_exact':'1/4'},
                                                 radius_case,radius_box,derivative))),
            'preconditioner_input_is_an_exact_binary64_fraction':binary_fraction(0.1)==Q(3602879701896397,2**55)
                and rejects(lambda:binary_fraction('0.1')) and rejects(lambda:binary_fraction(False))}


def prepare():
    PREP.mkdir(parents=True,exist_ok=False)
    with (PREP/'transfer_ambiguity_exact_bridge_reduce_driver_frozen.py').open('xb') as f:
        f.write(SOURCE.read_bytes())
    try:
        protocol=authority_protocol()
        checks=fixtures()
        require(yes(checks),'cold spatial reduction caller fixtures failed')
        ph=save(PREP/'protocol.json',protocol)
        fh=save(PREP/'preflight.json',{'status':'PASS','source_sha256':sha(SOURCE),'protocol_sha256':ph,
                'checks':checks,'actual_packet_payloads_opened':0,'actual_matrix_work':0,
                'test_scope':'Synthetic metadata and two scalar enclosure families, no execution of the kernel or actual replay.',
                'actual_root_claimed':False})
        print(json.dumps({'source_sha256':sha(SOURCE),'protocol_sha256':ph,'preflight_sha256':fh}),flush=True)
    except BaseException as error:
        save(PREP/'failure.json',{'error':str(error),'traceback':traceback.format_exc()})
        raise


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    sub=cli.add_subparsers(dest='stage',required=True)
    sub.add_parser('prepare')
    reduction=sub.add_parser('run')
    reduction.add_argument('--pilot-id',choices=PILOTS,required=True)
    for name in ('review','admission','replay-observation'):
        reduction.add_argument('--'+name,required=True,type=lambda v:Path(v).resolve())
        reduction.add_argument('--'+name+'-sha256',required=True)
    args=cli.parse_args()
    if args.stage=='prepare':
        prepare()
    else:
        run(args)
