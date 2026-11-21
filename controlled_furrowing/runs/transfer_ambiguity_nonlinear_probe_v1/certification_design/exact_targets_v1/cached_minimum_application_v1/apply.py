"""Apply already reviewed symmetric-part and cached derivative bounds to four cached roots."""

import argparse
from base64 import b64decode
from copy import deepcopy
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import re


HERE = Path(__file__).resolve().parent
EXACT = HERE.parent
DESIGN = EXACT.parent
ROOT = EXACT.parents[3]
CACHED = DESIGN/'cached_pilots_v1/independent_actual_preparation/actual_replay'
SPD = EXACT/'preconditioner_spd_v1/independent_actual_preparation/actual_replay'
DERIV = EXACT/'cached_derivative_extension_v1'
D_PREP = DERIV/'independent_actual_preparation'
D_REPLAY = D_PREP/'actual_replay_v2'
GEO = EXACT/'independent_actual/replay'
TARGETS = EXACT/'preparation/target_manifest.json'
THEOREM = EXACT/'preconditioner_spd_v1/METHOD.md'
SHORT = Q(1,2**36)
LARGE = Q(1,2**34)
N = 611
SITE_COUNTS = {'bulk':9216,'surface':288}
MARGIN_KEYS = {'bulk_det_min','bulk_hoop_min','bulk_J_min','normalized_weight_min',
               'surface_radius_min','surface_tangent_squared_min','surface_J_min'}
PIDS = tuple(f'{p}_{s}' for p in ('max_contrast_change','max_negative_contrast_change')
             for s in ('small_17','large_03'))
PINS = {
    CACHED/'receipt.json':'4e891afd79dd2f73508d1f878ddc1113ef107bdff8d9e0f673894f4d2f10ae17',
    CACHED/'terminal.json':'08da2c2f5754d3f8bf8b392157228f52ed316ac2e9d6af53fc8745181b658a8a',
    CACHED/'os_exit_observation.json':'65858f55c1281d3953b775f147c7475a8604891f37cf56ed0fdd47697e6df6e3',
    SPD/'receipt.json':'769a65a128aaa79018c050c7bd61d38999166621742946cd0f6a13876854ba83',
    SPD/'terminal.json':'d3555a7323fd0759a8473d3e734f8c1ee5401911c75246ec873ed5511103be89',
    SPD/'os_exit_observation.json':'b5a4a1ab78dc56c7a1147bb75b657dc141d3c2ed4a08ea31074bb5c5a2a0a3a0',
    D_PREP/'bridge_admission.json':'269229312793ec101cc301086d634dd4de89c62a5050ccfc60f185cf761dbc7d',
    D_REPLAY/'receipt.json':'f766b8b41142ce1f82f3e672f423c32485ae436487d9418c42d907b6f27e8c74',
    D_REPLAY/'terminal.json':'3804066eff55adfdce1c31a454b766b010b4bbfbcfb64b0c360c844f95a21781',
    D_REPLAY/'os_exit_observation.json':'be63043a9bafc6369f89ac2dcfe83674541e753cf221e923af18813ca2db3bc0',
    GEO/'receipt.json':'884e0bc76a27465c86b62a202cab491d764aedf030623592238e09adee86d77d',
    GEO/'terminal.json':'5011dc016786a4c3fd2a184cd4c620ad10e9d374d25d8d7ae02aac2d06d08c14',
    GEO/'os_exit_observation.json':'c2be4be8b403e170e20ec54c34f17a69b4aec4299e6f9f1ab7ecb3d4afa337d9',
    TARGETS:'cb25fbeaef753cea700c187cd8b377b7dcf4aa1fe2b0b246ed17134a07a89bc3',
    THEOREM:'0d6ca85eb910c616e6ccac15d67d63b37082d86818c1797b0304a744bf3e5750',
}
STATUSES = {
    CACHED: ('PASS_INDEPENDENT_FOUR_CACHED_LOCAL_ROOTS','PASS_OBSERVED_SERIAL_INDEPENDENT_REPLAY_TRUE_OS_EXIT'),
    SPD: ('PASS_INDEPENDENT_ACTUAL_FOUR_EXACT_SYMMETRIC_PRECONDITIONER_PARTS_SPD','PASS_OBSERVED_INDEPENDENT_ACTUAL_SPD_TRUE_OS_EXIT'),
    D_REPLAY: ('PASS_INDEPENDENT_FOUR_LARGER_CACHED_DERIVATIVES','PASS_OBSERVED_INDEPENDENT_DERIVATIVE_REPLAY_TRUE_OS_EXIT'),
    GEO: ('PASS_INDEPENDENT_ACTUAL_EIGHT_SINGLETON_GEOMETRY_BOXES','PASS_OBSERVED_INDEPENDENT_GEOMETRY_REPLAY_TRUE_OS_EXIT'),
}
SCOPE=('Application of a classical symmetric-Hessian Lyapunov argument to four already reviewed binary '
       'Q8/S12 cached-quadrature roots at four fixed existing targets. Normalized cached energy and '
       'audited closed free-coordinate cubes. No exact-spatial root or minimum, continuous loading, '
       'continuum, measured fit or novelty claim.')


def require(test,msg):
    if not test:
        raise ValueError(msg)


def yes(data):
    return type(data) is dict and bool(data) and all(v is True for v in data.values())


def name(path):
    return path.relative_to(ROOT).as_posix()


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def digest(obj):
    return sha256(json.dumps(obj,separators=(',',':'),sort_keys=True).encode()).hexdigest()


def pinned(path,h,links,parse=True):
    path=path.resolve()
    require(path.is_relative_to(ROOT) and type(h) is str and re.fullmatch('[0-9a-f]{64}',h),
            'explicit in-study authority path and hash required')
    raw=path.read_bytes()
    nm=name(path)
    require(sha256(raw).hexdigest()==h==links.setdefault(nm,h),f'changed frozen input: {nm}')
    return json.loads(raw) if parse else None


def stated(file,h,links,parse=True):
    path=(ROOT/file).resolve()
    require(path.is_relative_to(ROOT) and name(path)==file,'canonical original in-study path required')
    return pinned(path,h,links,parse)


def save(path,obj):
    with path.open('x',encoding='utf-8',newline='\n') as f:
        json.dump(obj,f,indent=2,allow_nan=False)
        f.write('\n')
    return sha(path)


def exact(x):
    require(type(x) in (Q,str,int),'explicit finite rational required')
    return Q(x)


def zero(x):
    return type(x) is int and x==0


def ordered(rows,k='pilot_id'):
    require(tuple(row[k] for row in rows)==PIDS,'four approved ordered cases required')
    return {r[k]:r for r in rows}


def one(rows,predicate):
    match=[r for r in rows if predicate(r)]
    require(len(match)==1,'one matching proof record required')
    return match[0]


def terminal(folder,l):
    rec=pinned(folder/'receipt.json',PINS[folder/'receipt.json'],l)
    end=pinned(folder/'terminal.json',PINS[folder/'terminal.json'],l)
    os=pinned(folder/'os_exit_observation.json',PINS[folder/'os_exit_observation.json'],l)
    pid=rec['pid'] if folder==SPD else rec['parents' if folder==CACHED else 'source_and_phase']['pid']
    s,t=STATUSES[folder]
    require(rec['status']==s and yes(rec['checks']) and end['status']=='FINISHED' and os['status']==t
            and all(type(x) is int and x==pid and x>0 for x in (pid,end['pid'],os['pid']))
            and zero(end['exit_code']) and zero(os['observed_exit_code'])
            and end['receipt_sha256']==os['receipt_sha256']==PINS[folder/'receipt.json']
            and os['terminal_sha256']==PINS[folder/'terminal.json'],
            'observed independent replay phase and exact receipt must match')
    return rec


def same_node_cubes(small,large):
    shared=('record_id','parameter_cell_sha256','model_observation_contract_sha256','point_cache_sha256',
            'reference_nodes_exact','current_centers_exact','triangles','free','predictor_hulls_exact')
    require(all(small[k]==large[k] for k in shared),'different center, model or free coordinates between cubes')
    free=large['free']
    require(len(free)==len(set(free))==N and all(type(x) is int and 0<=x<650 for x in free),
            'same complete free coordinate layout required')
    active=set(free)
    for radius,box in ((SHORT,small),(LARGE,large)):
        require(exact(box['corrector_radius_exact'])==radius and box['root_membership_assumed'] is False
                and box['any_spatial_derivative_bound_assumed'] is False,'exact same closed radii required')
        expected=[radius if k in active else Q(0) for k in range(650)]
        radii=[exact(v) for row in box['coordinate_radii_exact'] for v in row]
        require(radii==expected,'coordinate hulls do not implement the prescribed free cube')
        centers=[exact(c) for row in box['current_centers_exact'] for c in row]
        hulls=[tuple(map(exact,v)) for row in box['coordinate_hulls_exact'] for v in row]
        require(len(centers)==len(hulls)==650
                and all(lo==c-r and hi==c+r for c,r,(lo,hi) in zip(centers,radii,hulls,strict=True)),
                'wrong exact nodal box endpoints')


def scalar_bounds(point,der,calc,spd):
    eta=exact(point['independent_eta_upper_exact'])
    z_small=exact(point['independent_z_upper_exact'])
    e=exact(point['independent_posterior_error_upper_exact'])
    z1,z2=exact(der['reviewed_ZQ_upper_exact']),exact(calc['ZQ_upper_exact'])
    z=max(z1,z2)
    norm=exact(spd['exact_R_infinity_norm'])
    q=exact(spd['max_reviewed_q_upper_exact'])
    require(point['status']=='PASS_INDEPENDENT_CACHED_LOCAL_ROOT'
            and exact(point['radius_exact'])==SHORT and 0<=z_small<1 and eta>=0
            and eta+z_small*SHORT==exact(point['independent_selfmap_lhs_upper_exact'])<SHORT
            and 0<=e==eta/(1-z_small)<SHORT<LARGE,
            'one existing strictly interior root in the smaller nested cube required')
    require(der['status']=='PASS_INDEPENDENT_LARGER_RADIUS_CACHED_DERIVATIVE'
            and exact(der['radius_exact'])==exact(calc['whole_box_cache_radius_exact'])==LARGE
            and exact(der['producer_ZQ_upper_exact'])==z2 and 0<=min(z1,z2)<=z<1
            and calc['whole_box_controls_only_cached_quadrature'] is True
            and exact(der['prior_point_eta_upper_exact'])==eta,
            'only an admitted same-center larger cached derivative bound with z<1 is sufficient')
    require(spd['status']=='PASS_INDEPENDENT_EXACT_SYMMETRIC_PART_SPD'
            and spd['exact_target_symmetric_part_spd'] is True
            and 0<=q==max(exact(spd['fresh_congruence_q_upper_exact']),exact(spd['producer_congruence_q_upper_exact']))<1
            and exact(spd['strict_congruence_margin_lower_exact'])==1-q
            and 0<norm==exact(point['preconditioner_infinity_norm_exact'])
            ==exact(der['preconditioner_infinity_norm_exact']),
            'same strictly positive Sym(R) and exact binary infinity norm required')
    m=(1-z)/norm
    return {'radius_cached_root_exact':str(SHORT),'radius_coercive_cube_exact':str(LARGE),
            'root_distance_upper_exact':str(e),'root_interior_margin_small_cube_exact':str(SHORT-e),
            'root_interior_margin_coercive_cube_exact':str(LARGE-e),
            'cached_ZQ_on_coercive_cube_upper_exact':str(z),'exact_R_infinity_norm':str(norm),
            'Sym_R_congruence_defect_upper_exact':str(q),'normalized_cached_H_eigenlower_exact':str(m),
            'normalized_cached_energy_quadratic_gap_coefficient_exact':str(m/2),
            'display':{'curvature_lower':float(m),'gap_coefficient':float(m/2),
                       'interior_margin_small':float(SHORT-e),'interior_margin_large':float(LARGE-e)}}


def exact_identity(pilot,entry,contract,packet,box):
    drive=contract['exact_requested_drive']
    dose=exact(drive['optical_dose_exact'])
    gain=exact(drive['common_gain_exact'])
    g=exact(drive['activation_strength_exact'])
    scale=exact(contract['energy_scale_binary_exact'])
    cell=contract['parameter_cell']
    digest_cell=digest(cell)
    require(dose>0 and gain>0 and g==dose*gain==exact(packet['activation_strength_exact'])
            and g!=exact(packet['activation_strength_binary_exact'])
            and 0<scale==exact(packet['normalization_binary_exact'])
            and all(exact(lo)==exact(hi) for lo,hi in cell.values())
            and tuple(map(exact,cell['common_gain']))==(gain,gain)
            and digest_cell==contract['parameter_cell_sha256']==box['parameter_cell_sha256']
            and all(drive[k]==packet[k] for k in drive)
            and contract['pilot_id']==packet['pilot_id']==entry['pilot_id']==pilot
            and contract['record_id']==box['record_id']
            and contract['point_cache_sha256']==box['point_cache_sha256']
            and contract['source_field_sha256']==packet['input_state_sha256']
            and contract['frozen_binary_cache_sha256']==packet['frozen_binary_cache']['sha256']
            and entry['model_contract_sha256']==box['model_observation_contract_sha256'],
            'exact common model, dose, drive, cell, cache and normalization are required')
    return {'pilot_id':pilot,'record_id':contract['record_id'],'field_sha256':contract['source_field_sha256'],
            'cache_packet_sha256':contract['point_cache_sha256'],'embedded_cache_sha256':contract['frozen_binary_cache_sha256'],
            'R_binary_sha256':contract['binary_preconditioner_sha256'],
            'model_contract_sha256':entry['model_contract_sha256'],'parameter_cell_sha256':digest_cell,
            'dose_exact':str(dose),'gain_exact':str(gain),'drive_exact':str(g),'energy_scale_binary_exact':str(scale)}


def same_case_links(expected,spd,der,producer,point,box_digest,point_sha,der_sha):
    calc=producer['case']
    require(spd['target_identity']==expected and spd['record_id']==expected['record_id']
            and spd['cached_point_proof_sha256']==calc['verified_point_sha256']==point_sha
            and spd['derivative_case_sha256']==der_sha
            and spd['larger_radius_box_sha256']==der['box_sha256']==calc['box_sha256']==box_digest
            and spd['cache_packet_sha256']==der['cache_packet_sha256']==point['cache_packet_sha256']
            ==calc['cache_packet_sha256']==expected['cache_packet_sha256']
            and point['cache_array_contents_sha256']==calc['embedded_cache_sha256']==expected['embedded_cache_sha256']
            and calc['preconditioner_binary_sha256']==expected['R_binary_sha256']
            and point['field_sha256']==calc['field_sha256']==expected['field_sha256']
            and point['pilot_id']==calc['pilot_id']==der['pilot_id']==expected['pilot_id']
            and calc['record_id']==expected['record_id']
            and exact(point['exact_strength'])==exact(der['proof_drive_exact'])
            ==exact(calc['actual_exact_drive'])==exact(expected['drive_exact']),
            'same point/derivative/SPD case, cache, binary R, cube and drive required')


def analyticity_and_geometry(der,geo,geo_small,expected,pilot):
    domain=der['whole_box_domains']
    lower=domain['strict_lower_margins']
    require(der['certified_site_counts']==domain['counts']==SITE_COUNTS
            and set(lower)==MARGIN_KEYS and all(0<exact(x['lower_exact'])<=exact(x['upper_exact']) for x in lower.values())
            and min(exact(lower[k]['lower_exact']) for k in ('bulk_J_min','surface_J_min'))>Q(1,10**6),
            'strict cached-site analytic domain with inactive low-J penalty required')
    for r,g in ((SHORT,geo_small),(LARGE,geo)):
        require(g['status']=='PASS_FRESH_SINGLETON_WHOLE_GEOMETRY_ONLY'
                and g['pilot_id']==pilot and g['record_id']==expected['record_id']
                and g['proof_center_field_sha256']==expected['field_sha256']
                and exact(g['radius_exact'])==r
                and 0<exact(g['strong_monotonicity_lower_exact'])
                and 0<exact(g['body_J_bounds_exact'][0])<=exact(g['body_J_bounds_exact'][1]),
                'whole-cube embedding and exact positive geometry for the same field required')
    return {k:item['lower_exact'] for k,item in lower.items()}


def input_rows():
    links={}
    pinned(THEOREM,PINS[THEOREM],links,False)
    p=terminal(CACHED,links)
    s=terminal(SPD,links)
    d=terminal(D_REPLAY,links)
    g=terminal(GEO,links)
    w=pinned(D_PREP/'bridge_admission.json',PINS[D_PREP/'bridge_admission.json'],links)
    m=pinned(TARGETS,PINS[TARGETS],links)
    require(p['root_certificates_cached_quadrature_only']==4 and zero(p['root_certificates_exact_spatial_integral'])
            and s['fresh_actual_matrix_congruences']==4 and s['root_or_minimum_certified'] is False
            and zero(d['spatial_force_hessian_or_solve_jobs']) and zero(g['exact_spatial_root_count'])
            and zero(g['spatial_derivative_or_tail_bound_count'])
            and w['status']=='PASS_REVIEWED_LARGER_CACHED_DERIVATIVE_BRIDGE_INPUTS' and yes(w['checks'])
            and w['native_independent_review_file']==name(D_REPLAY/'receipt.json')
            and w['native_independent_review_sha256']==PINS[D_REPLAY/'receipt.json']
            and w['native_independent_terminal_sha256']==PINS[D_REPLAY/'terminal.json']
            and w['native_independent_os_exit_sha256']==PINS[D_REPLAY/'os_exit_observation.json']
            and exact(w['radius_exact'])==LARGE and w['reviewed_cases']==d['reviewed_cases']
            and m['status']=='FOUR_SINGLETON_TARGETS_PREPARED'
            and tuple(m['pilot_order'])==PIDS and tuple(map(exact,m['radii_exact']))==(SHORT,LARGE),
            'exact scope of four actual root, SPD, geometry and derivative proofs required')
    root_rows,spd_rows,der_rows,targets=map(ordered,(p['field_results'],s['ordered_cases'],d['reviewed_cases'],m['targets']))
    require(len(g['geometry_records'])==8,'complete two-radius geometry proof required')
    d_top=DERIV/'actual/results.json'
    d_end=DERIV/'actual/terminal.json'
    dv,dt=(pinned(path,w[key],links) for path,key in ((d_top,'derivative_results_sha256'),(d_end,'derivative_terminal_sha256')))
    require(all(d['producer_and_input_sha256'][name(path)]==links[name(path)] for path in (d_top,d_end))
            and dt['status']=='FINISHED' and dt['results_sha256']==w['derivative_results_sha256']
            and dv['bindings']==dt['bindings'] and type(dt['pid']) is int and dt['pid']==d['producer_exit']['supervisor_pid']
            and zero(d['producer_exit']['observed_os_exit_code']),'native derivative replay must bind identical actual producer')
    producer_rows=ordered(dv['case_results'])
    bundles=[]
    for pilot in PIDS:
        er,rr,ss,dd,cc=(pool[pilot] for pool in (targets,root_rows,spd_rows,der_rows,producer_rows))
        contract=stated(er['model_contract'],er['model_contract_sha256'],links)
        bp=stated(er['box_0'],er['box_0_sha256'],links)
        bx=stated(er['box_1'],er['box_1_sha256'],links)
        same_node_cubes(bp,bx)
        packet=stated(contract['point_cache'],contract['point_cache_sha256'],links)
        require(sha256(b64decode(packet['frozen_binary_cache']['data'],validate=True)).hexdigest()
                ==packet['frozen_binary_cache']['sha256'],'embedded cached packet bytes changed')
        expected=exact_identity(pilot,er,contract,packet,bx)
        rp=name(CACHED/(pilot+'.json'))
        require(contract['cached_root_independent_receipt']==rp
                and contract['cached_root_independent_receipt_sha256']==rr['receipt_sha256']
                and rr['status']=='PASS_INDEPENDENT_CACHED_LOCAL_ROOT'
                and rr['field_sha256']==expected['field_sha256'],'same original independent cached point receipt required')
        point=stated(rp,rr['receipt_sha256'],links)
        require(point['record_id']==rr['record_id']==contract['physical_source_record_id']==packet['record_id']
                and point['rounded_strength_used_in_proof'] is False
                and tuple(map(exact,(contract['cached_root_radius_exact'],point['radius_exact'])))==(SHORT,SHORT)
                and exact(contract['cached_root_error_exact'])==exact(point['independent_posterior_error_upper_exact']),
                'same smaller original cached local root required')
        require(ss['case_receipt_file']==name(SPD/(pilot+'.json'))
                and dd['native_case_file']==name(D_REPLAY/(pilot+'.json'))
                and ss['status']=='PASS_INDEPENDENT_EXACT_SYMMETRIC_PART_SPD'
                and exact(dd['radius_exact'])==LARGE,
                'explicit same-case native SPD and large derivative receipts required')
        spd=stated(ss['case_receipt_file'],ss['case_receipt_sha256'],links)
        der=stated(dd['native_case_file'],dd['native_case_sha256'],links)
        require(all(spd[k]==ss[k] for k in ('target_identity','fresh_congruence_q_upper_exact','max_reviewed_q_upper_exact',
                                           'exact_R_infinity_norm','cache_packet_sha256'))
                and all(spd['replay_phase'][k]==s[k] for k in ('pid','source_sha256','preflight_sha256','protocol_sha256'))
                and dd['reviewed_ZQ_upper_exact']==der['reviewed_ZQ_upper_exact'],
                'native summary differs from original SPD or derivative case')
        calc_file=DERIV/'actual'/pilot/'results.json'
        require(d['producer_and_input_sha256'][name(calc_file)]==cc['results_sha256']==dd['producer_result_sha256']
                ==der['producer_result_sha256'],'native larger derivative review did not bind exact producer inputs')
        calc=pinned(calc_file,cc['results_sha256'],links)
        require(calc['status']=='FRESH_LARGER_RADIUS_CACHED_DERIVATIVE_UNREVIEWED'
                and calc['ZQ_upper_exact']==cc['ZQ_upper_exact'] and calc['bindings']==dt['bindings']
                and calc['exact_spatial_derivatives_computed'] is False,
                'same cached derivative producer domain required')
        same_case_links(expected,spd,der,calc,point,er['box_1_sha256'],rr['receipt_sha256'],dd['native_case_sha256'])
        solved=[]
        for radius,key in ((SHORT,'box_0_sha256'),(LARGE,'box_1_sha256')):
            gr=one(g['geometry_records'],lambda r:r['pilot_id']==pilot and exact(r['radius_exact'])==radius)
            gp=stated(gr['receipt'],gr['sha256'],links)
            require(gp['box_sha256']==gr['box_sha256']==er[key]
                    and gp['body_geometry_sha256']==gr['body_geometry_sha256'],
                    'independent whole-body geometry belongs to another cube')
            solved.append((gr,gp))
        positivity=analyticity_and_geometry(der,solved[1][1],solved[0][1],expected,pilot)
        bundles.append((expected,point,der,calc,spd,positivity,{'cached_point_sha256':rr['receipt_sha256'],
                    'native_SPD_case_sha256':ss['case_receipt_sha256'],'native_large_cached_ZQ_case_sha256':dd['native_case_sha256'],
                    'producer_large_cached_ZQ_case_sha256':cc['results_sha256'],
                    'small_geometry_sha256':solved[0][0]['sha256'],'large_geometry_sha256':solved[1][0]['sha256'],
                    'small_box_sha256':er['box_0_sha256'],'large_box_sha256':er['box_1_sha256']}))
    return dict(sorted(links.items())),bundles


def rejects(f):
    try:
        f()
    except ValueError:
        return True
    return False


def fixtures(bundles):
    identity,point,der,calc,spd,_,hashes=bundles[0]
    def changed(source,key,v):
        return {**source,key:v}
    def with_link(idn=identity,matrix=spd,solver=calc,root=point):
        return same_case_links(idn,matrix,der,solver,root,hashes['large_box_sha256'],
                               hashes['cached_point_sha256'],hashes['native_large_cached_ZQ_case_sha256'])
    with_link()
    def bounds(p=point,d=der,c=calc,s=spd):
        return scalar_bounds(p,d,c,s)
    values=bounds()
    check_m=exact(values['normalized_cached_H_eigenlower_exact'])
    check_z=exact(values['cached_ZQ_on_coercive_cube_upper_exact'])
    check_r=exact(values['exact_R_infinity_norm'])
    mixed=deepcopy(calc)
    mixed['case']['preconditioner_binary_sha256']=bundles[1][0]['R_binary_sha256']
    norm_other=exact(bundles[1][4]['exact_R_infinity_norm'])
    return {'exact_cached_scalar_margin_and_gap_identities':
             check_m>0 and check_m*check_r+check_z==1
             and exact(values['normalized_cached_energy_quadratic_gap_coefficient_exact'])==check_m/2,
            'reject_different_binary_R_or_signed_infinity_norm':rejects(lambda:with_link(solver=mixed))
             and rejects(lambda:bounds(s=changed(spd,'exact_R_infinity_norm',str(norm_other)))),
            'reject_different_cache_field_model_or_root_source':(
             rejects(lambda:with_link(matrix=changed(spd,'cache_packet_sha256','0'*64)))
             and rejects(lambda:with_link(matrix=changed(spd,'cached_point_proof_sha256','0'*64)))
             and rejects(lambda:with_link(root=changed(point,'field_sha256','0'*64)))
             and rejects(lambda:with_link(idn=changed(identity,'model_contract_sha256','0'*64)))
             and rejects(lambda:pinned(THEOREM,'0'*64,{},False))),
            'reject_different_radius_or_nonpositive_curvature_domain':(
             rejects(lambda:bounds(d=changed(der,'radius_exact',str(SHORT))))
             and rejects(lambda:bounds(d=changed(der,'reviewed_ZQ_upper_exact','1')))
             and rejects(lambda:bounds(c=changed(calc,'whole_box_controls_only_cached_quadrature',False)))),
            'reject_missing_SPD_guard_or_congruence_equality':(
             rejects(lambda:bounds(s=changed(spd,'exact_target_symmetric_part_spd',False)))
             and rejects(lambda:bounds(s=changed(spd,'max_reviewed_q_upper_exact','1')))),
            'reject_posterior_root_bound_at_small_cube_boundary':rejects(
             lambda:bounds(p=changed(point,'independent_posterior_error_upper_exact',str(SHORT)))),
            'contraction_alone_counterexample_preserved':(
             (-1)*(-1)==1 and (-1)<0 and 'R=H_Q=-I' in (HERE/'PROOF.md').read_text(encoding='utf-8'))}


def protocol(links):
    return {'status':'FROZEN_SCALAR_CACHE_MINIMUM_APPLICATION','scope':SCOPE,'source_sha256':sha(Path(__file__)),
            'proof_sha256':sha(HERE/'PROOF.md'),'theorem_source_sha256':PINS[THEOREM],
            'source_and_input_sha256':links,'ordered_pilot_ids':list(PIDS),
            'radius_cached_root_exact':str(SHORT),'radius_coercive_cube_exact':str(LARGE),
            'normalized_curvature':'(1-max(pinned producer ZQ, independent native ZQ))/exact binary R infinity norm',
            'energy_gap':'m_Q/2 times squared Euclidean displacement from the existing stationary root',
            'C2_premise':'Positive aligned cached analytic site domains on each entire closed cube, fixed linear free coordinates.',
            'application_requires_independent_review_before_promotion':True,'new_matrix_AD_force_or_solver_work':0}


def prepare():
    links,bundles=input_rows()
    p=protocol(links)
    checks=fixtures(bundles)
    require(yes(checks),'cached minimum application preflight failed')
    with (HERE/'apply_frozen.py').open('xb') as f:
        f.write(Path(__file__).read_bytes())
    ph=save(HERE/'protocol.json',p)
    pf=save(HERE/'preflight.json',{'status':'PASS','source_sha256':p['source_sha256'],'protocol_sha256':ph,
          'checks':checks,'new_production_matrix_AD_force_or_solver_work':0})
    print(json.dumps({'source_sha256':p['source_sha256'],'protocol_sha256':ph,'preflight_sha256':pf}),flush=True)


def apply(ph,fh):
    p=pinned(HERE/'protocol.json',ph,{})
    pf=pinned(HERE/'preflight.json',fh,{})
    links,bundles=input_rows()
    require(p==protocol(links) and pf['status']=='PASS' and yes(pf['checks'])
            and pf['protocol_sha256']==ph and pf['source_sha256']==p['source_sha256']
            and p['source_sha256']==sha(HERE/'apply_frozen.py')
            and zero(pf['new_production_matrix_AD_force_or_solver_work']),
            'frozen source, math scope, preflight or input bytes changed')
    rows=[]
    for identity,point,der,calc,spd,positivity,hashes in bundles:
        b=scalar_bounds(point,der,calc,spd)
        rows.append({'target_identity':identity,'premise_sha256':hashes,'analytic_cached_site_lower_margins':positivity,
                     'cached_site_counts':SITE_COUNTS,**b})
    r={'status':'DRAFT_FOUR_STRICT_CACHED_QUADRATURE_MINIMA_AWAITING_INDEPENDENT_REVIEW','scope':SCOPE,
       'source_sha256':p['source_sha256'],'proof_sha256':p['proof_sha256'],'protocol_sha256':ph,'preflight_sha256':fh,
       'source_and_input_sha256':links,'candidate_count':4,'rows':rows,
       'new_matrix_AD_force_or_solver_work':0,'exact_spatial_root_or_minimum_certified':False,
       'full_five_gate_completion_claimed':False}
    rh=save(HERE/'results.json',r)
    reread=pinned(HERE/'results.json',rh,{})
    for row in reread['rows']:
        m,z,nr,e=map(exact,(row['normalized_cached_H_eigenlower_exact'],
                            row['cached_ZQ_on_coercive_cube_upper_exact'],row['exact_R_infinity_norm'],
                            row['root_distance_upper_exact']))
        require(0<m and 0<nr and z+m*nr==1 and 0<=e<SHORT<LARGE
                and exact(row['normalized_cached_energy_quadratic_gap_coefficient_exact'])*2==m
                and e+exact(row['root_interior_margin_small_cube_exact'])==SHORT
                and e+exact(row['root_interior_margin_coercive_cube_exact'])==LARGE,
                'serialized exact output arithmetic drift')
    rb=save(HERE/'readback.json',{'status':'AUTHOR_READBACK_OF_SCALAR_DRAFT_PASS','results_sha256':rh,
             'proof_sha256':p['proof_sha256'],'scope':SCOPE,'exact_scalar_identity_readbacks':4,
             'independent_actual_application_review_required':True,'no_new_production_numerics':True})
    print(json.dumps({'status':r['status'],'results_sha256':rh,'readback_sha256':rb,
                      'display':[row['display'] for row in rows]}),flush=True)


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    cmd=cli.add_subparsers(dest='stage',required=True)
    cmd.add_parser('prepare')
    a=cmd.add_parser('apply')
    a.add_argument('--protocol-sha256',required=True)
    a.add_argument('--preflight-sha256',required=True)
    args=cli.parse_args()
    if args.stage=='prepare':
        prepare()
    else:
        apply(args.protocol_sha256,args.preflight_sha256)
