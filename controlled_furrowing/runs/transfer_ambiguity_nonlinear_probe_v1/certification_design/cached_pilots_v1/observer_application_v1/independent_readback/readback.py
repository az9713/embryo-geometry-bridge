"""Independent exact arithmetic readback of the cached-root observer application."""

import base64
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
import io
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
APP = HERE.parent
CACHE = APP.parent
DESIGN = CACHE.parent
RUN = DESIGN.parent
ROOT = RUN.parents[1]
EXACT = DESIGN / 'exact_targets_v1'
CACHED_REPLAY = CACHE / 'independent_actual_preparation/actual_replay'
GEOMETRY_REPLAY = EXACT / 'independent_actual/replay'
L = Q(67, 2)
RADIUS = Q(1, 2**36)
VARIANTS = ('current_disk__spatial_max', 'current_disk__regional_mean',
            'current_radial_line__spatial_max', 'current_radial_line__regional_mean')
OBS = ('maximum_shortening_um', 'region_averaged_shortening_um',
       'apical_max_shortening_um', 'basal_max_shortening_um')
ENDPOINTS = ('max_contrast_change', 'max_negative_contrast_change')
PINS = {
    APP/'receipt.json': '0081a5939496f0c7153360d82bf63384dbe172d1da07b69e06e568a1f6e2b5a4',
    APP/'results.json': '21dc32c0b8bf073ffa9bf7a3c163cca107698e0f847304d53483caf071db1a40',
    APP/'protocol.json': '28a3e85bbce2d2821d8a5e38399667065eec4952f2252ead3492718cdd293299',
    APP/'preflight.json': '7e6c0fe2110c672fa3b0f15d6852a8254f71c8e74699e7943fa0a2aa9a88c90b',
    CACHED_REPLAY/'receipt.json': '4e891afd79dd2f73508d1f878ddc1113ef107bdff8d9e0f673894f4d2f10ae17',
    GEOMETRY_REPLAY/'receipt.json': '884e0bc76a27465c86b62a202cab491d764aedf030623592238e09adee86d77d',
    RUN/'observer_robustness/independent_actual/receipt.json': '703b574bc51a083b10c7328d9622862e939ddabe791b18ce9e1817f07718f69f',
    RUN/'observer_robustness/actual/results.json': '53a17296ddae08b4d1fe6b7e9b20ad0a377bf8dc962142bbfb9eb9877a381647',
    RUN/'independent_actual_preparation/actual_replay/receipt.json': 'c795ae78ba0b53a01c3d65f466d5fa2f1d233c2a68c4275c04260c3cb5120ebd',
    RUN/'field_manifest.json': '308ca4c05e476ba6a14cbc7b8d34018983b8f11263e2c224488672b35eb38ac0',
    RUN/'preparation/protocol.json': '9d5d38bb5f9785d40bf843dd225eb93b1d5ecbe96f20942aec0023765d119d94',
    EXACT/'preparation/target_manifest.json': 'cb25fbeaef753cea700c187cd8b377b7dcf4aa1fe2b0b246ed17134a07a89bc3',
    DESIGN/'budgets_v1/actual/results.json': '689948766db950ee4ddf399e6731aea08890cc7a1f84c5615d19e7ded68b0872',
    DESIGN/'budgets_v1/independent_actual/receipt.json': '978a2a4095bb8735103a5ee2c11ad266a6b0e052e22af2b7a36ef29f1f38265c',
    ROOT/'evidence/nodal_observer_error_independent_review.json': '335626b91c72ac5ec37b1b1dd95b4825d3dc847df02d9cf993e5ed6daff5f35a',
}
SCOPE = ('Independent arithmetic application to four reviewed cached-quadrature roots. Two small L states, '
         '32 static comparisons, and two large states. The prior root, geometry, stored observer and '
         'nodal-to-graph theorem are reused. Exact-spatial FE, stability, PDE, source compatibility, '
         'the remaining 32 small states and full-goal certification remain separate.')


def check(ok, message):
    if not ok:
        raise ValueError(message)


def sha(raw):
    return sha256(raw).hexdigest()


def bind(path, expected, links, parse=True):
    path = path.resolve()
    check(path.is_relative_to(ROOT), 'path outside study')
    name = path.relative_to(ROOT).as_posix()
    raw = path.read_bytes()
    check(sha(raw) == expected == links.get(name, expected), f'changed bound file {name}')
    links[name] = expected
    return json.loads(raw) if parse else None


def one(rows, predicate):
    found = [r for r in rows if predicate(r)]
    check(len(found) == 1, 'unique matching authority row required')
    return found[0]


def interval(vals):
    a,b=map(Q,vals)
    check(a<=b,'ordered exact endpoints required')
    return a,b


def as_text(vals):
    return list(map(str,vals))


def box_q(rows):
    return [tuple(map(Q,row)) for row in rows]


def physical_endpoints(reference,triangles,z):
    edges={}
    for a,b,c,ab,bc,ca in triangles:
        for i,j,m in ((a,b,ab),(b,c,bc),(c,a,ca)):
            if reference[i][1] == reference[j][1] == z:
                key=tuple(sorted((i,j),key=lambda k:reference[k][0]))
                check(edges.get(key,m)==m,'conforming P2 surface nodes required')
                edges[key]=m
    ordered=sorted((i,m,j) for (i,j),m in edges.items())
    ordered.sort(key=lambda entry:reference[entry[0]][0])
    check(len(ordered)==12 and all(a[-1]==b[0] for a,b in zip(ordered,ordered[1:])),
          'expected twelve conforming ordered trace edges')
    return ordered


def current_coordinates(state,packet,proof,contract,box):
    p=packet['frozen_binary_cache']
    raw=base64.b64decode(p['data'],validate=True)
    check(sha(raw)==p['sha256']==proof['cache_array_contents_sha256']==contract['frozen_binary_cache_sha256'],
          'independently reviewed frozen cache bytes required')
    with np.load(io.BytesIO(raw),allow_pickle=False) as npz:
        X=box_q(npz['reference'].tolist())
        tris=npz['triangles'].tolist()
        ix=list(map(int,npz['free']))
        u=list(map(Q,npz['u'].tolist()))
        R=npz['preconditioner']
        apertures=[npz[f'surface_{f}_edges'].tolist() for f in range(2)]
        check(R.shape==(611,611) and np.isfinite(R).all() and sha(R.tobytes(order='C'))==contract['binary_preconditioner_sha256'],
              'same independently reviewed preconditioner required')
    check(X==box_q(state['reference_nodes'])==box_q(box['reference_nodes_exact'])
          and tris==state['triangles']==box['triangles']
          and u==list(map(Q,state['free_displacement'])) and ix==box['free'] and len(ix)==611
          and len(X)==325 and len(tris)==144,'cached, serialized and reviewed geometry data must coincide')
    rim=Q(state['parameters']['outer_radius'])
    free=[2*i+d for i,(r,_) in enumerate(X) for d in range(2) if r!=rim and (r!=0 or d==1)]
    check(ix==free,'dereferenced ordered axis/rim-free mask must match')
    increments=dict(zip(ix,u,strict=True))
    Y=[tuple(x+increments.get(2*i+d,0) for d,x in enumerate(node)) for i,node in enumerate(X)]
    stored=box_q(state['current_nodes'])
    check(Y==box_q(box['current_centers_exact']) and len(stored)==325
          and [tuple(Q(float(c)) for c in row) for row in Y]==stored,
          'saved current nodes must be nearest binary representation of the exact center')
    worst=Q(0)
    for i,(xr,yr,sr,rr,hr) in enumerate(zip(X,Y,stored,box['coordinate_radii_exact'],
                                           box['coordinate_hulls_exact'],strict=True)):
        for j,(x,y,s,r,h) in enumerate(zip(xr,yr,sr,rr,hr,strict=True)):
            active=2*i+j in increments
            expected=RADIUS if active else Q(0)
            check(Q(r)==expected and interval(h)==(y-expected,y+expected)
                  and (active or x==y==s),'one or more of 650 exact coordinate hulls drifted')
            worst=max(worst,abs(y-s))
    height=Q(state['parameters']['height'])
    faces={}
    dofs={}
    for name,z,cache in zip(('bottom','top'),(Q(0),height),apertures,strict=True):
        tri_trace=physical_endpoints(X,tris,z)
        check([list(edge) for edge in tri_trace]==cache,'quadratic trace differs from cached integration data')
        i,j=tri_trace[0][0],tri_trace[-1][-1]
        check(2*i not in increments and 2*j not in increments and Y[i][0]==0 and Y[j][0]==rim,
              'fixed axis and rim radius required by observer theorem')
        dofs[name]=[2*i,2*j]
        faces[name]=[tuple(tuple(L*v for v in stored[k]) for k in edge) for edge in tri_trace]
        check(faces[name]==[tuple(tuple(L*Q(v) for v in point) for point in edge)
                            for edge in state[f'boundary_{name}_quadratic_controls']],
              'verified observer graph controls differ from the serialized full field')
    return faces,dofs,worst,L*height


def face_error(controls,epsilon,roi):
    # Independently compute endpoint derivatives: exact P2 derivative is affine in t.
    signs=(-1,1)
    derivative=[]
    for edge in controls:
        p,m,q=edge
        # For Lagrange t={0,1/2,1}, d/dt = (q-p) + (2t-1)*2*(p+q-2m).
        d=[tuple((b-a)+s*2*(a+b-2*c) for a,b,c in zip(p,q,m,strict=True)) for s in signs]
        derivative.append(d)
    radial_min=min(d[0] for edge in derivative for d in edge)
    check(controls[0][0][0]<=roi[0]<roi[1]<=controls[-1][-1][0]
          and radial_min>8*epsilon,'reviewed strictly increasing common radius requirement failed')
    slope=max(max(abs(d[1]) for d in edge)/min(d[0] for d in edge) for edge in derivative)
    err=Q(5,4)*(1+slope)*epsilon
    return radial_min,slope,err


def required_bound(saved,error,observed,delta):
    p,q=saved
    L0,U0=p-error,q+error
    v=(max(Q(0),L0-observed,observed-U0),max(abs(U0-observed),abs(L0-observed)))
    return (L0,U0),v,'IN' if v[1]<=delta else 'OUT' if v[0]>delta else 'UNRESOLVED'


def compare_pilot(candidate,docs,links,delta):
    pid=candidate['pilot_id']
    proof_top=docs[CACHED_REPLAY/'receipt.json']
    geo_top=docs[GEOMETRY_REPLAY/'receipt.json']
    proof_summary=one(proof_top['field_results'],lambda e:e['pilot_id']==pid)
    field=one(docs[RUN/'field_manifest.json']['fields'],lambda e:e['field']==proof_summary['field'])
    stage=field['endpoint_id']
    target=one(docs[EXACT/'preparation/target_manifest.json']['targets'],lambda r:r['pilot_id']==pid)
    contract=bind(ROOT/target['model_contract'],target['model_contract_sha256'],links)
    proof=bind(ROOT/contract['cached_root_independent_receipt'],proof_summary['receipt_sha256'],links)
    cache_digest=contract['point_cache_sha256']
    check(proof['status']=='PASS_INDEPENDENT_CACHED_LOCAL_ROOT' and proof['pilot_id']==pid
          and candidate['cached_root_receipt_sha256']==proof_summary['receipt_sha256']
          and proof['cache_packet_sha256']==cache_digest==proof_top['producer_and_payload_sha256'][contract['point_cache']]
          and proof['field_sha256']==field['field_sha256']==contract['source_field_sha256']
          and field['endpoint_id']==candidate['endpoint_id'] and field['role']==candidate['role']
          and proof['record_id']==candidate['record_id'] and proof['rounded_strength_used_in_proof'] is False,
          'root/candidate field authority drift')
    native_primary=docs[RUN/'independent_actual_preparation/actual_replay/receipt.json']
    for k in ('field','verification'):
        check(native_primary['all_hashes'][field[k]]==field[k+'_sha256'],'field lacks independent nonlinear binding')
    state=bind(ROOT/field['field'],field['field_sha256'],links)
    verification=bind(ROOT/field['verification'],field['verification_sha256'],links)
    packet=bind(ROOT/contract['point_cache'],cache_digest,links)
    step=field['step']
    dose,gain=Q(step['dose_exact']),Q(packet['common_gain_exact'])
    drive=dose*gain
    check(drive==Q(step['requested_strength_exact'])==Q(proof['exact_strength'])==Q(candidate['exact_drive'])
          and drive!=Q(state['activation_strength'])==Q(step['solver_strength_binary64_exact'])
          and dose==Q(candidate['exact_dose']) and gain==Q(candidate['exact_gain'])
          and packet['parameters']==state['parameters']
          and tuple(map(Q,contract['parameter_cell']['common_gain']))==(gain,gain),'exact nonrounded physical drive mismatch')
    eta,z0,v,z,r=[Q(proof[k]) for k in ('independent_eta_upper_exact','independent_Z0_upper_exact',
             'independent_variation_upper_exact','independent_z_upper_exact','radius_exact')]
    check(r==RADIUS and z==z0+v and 0<=z<1 and 0<=eta and eta+z*r<r
          and Q(proof['independent_selfmap_lhs_upper_exact'])==eta+z*r,
          'strict independent cached-radius contraction must hold')
    err=eta/(1-z)
    check(err==Q(proof['independent_posterior_error_upper_exact'])<r
          and [Q(candidate[k]) for k in ('eta_exact','z_exact','certified_radius_exact','posterior_radius_exact')]==[eta,z,r,err],
          'recomputed exact posterior error does not match')
    geom=one(geo_top['geometry_records'],lambda e:e['pilot_id']==pid and Q(e['radius_exact'])==r)
    geom_proof=bind(ROOT/geom['receipt'],geom['sha256'],links)
    check(geom_proof['status']=='PASS_FRESH_SINGLETON_WHOLE_GEOMETRY_ONLY'
          and geom_proof['proof_center_field_sha256']==field['field_sha256']
          and geom['sha256']==candidate['geometry_receipt_sha256']
          and geom_proof['box_sha256']==geom['box_sha256']==candidate['box_sha256']==target['box_0_sha256']
          and Q(geom_proof['strong_monotonicity_lower_exact'])>0
          and interval(geom_proof['body_J_bounds_exact'])[0]>0,
          'current r=2^-36 exact same-center complete geometry premise required')
    box=bind(ROOT/target['box_0'],target['box_0_sha256'],links)
    faces,dofs,xi,height=current_coordinates(state,packet,proof,contract,box)
    epsilon=L*(err+xi)
    subject=one(docs[DESIGN/'budgets_v1/actual/results.json']['field_budgets'],lambda x:x['field']==field['field'])
    check(xi==Q(subject['center_and_model']['node_addition_gap_dimensionless_exact'])
          ==Q(candidate['serialized_X_plus_u_gap_exact']) and epsilon==Q(candidate['physical_nodal_error_um_exact'])
          and dofs==candidate['fixed_radial_endpoint_dofs'],'serialized addition error or fixed endpoint readback mismatch')
    pair_protocol=docs[RUN/'preparation/protocol.json']
    is_small=field['role']=='small_record'
    roi=(Q(0),L*Q(state['parameters']['patch_radius']) if is_small
         else Q(pair_protocol['large_observable']['shoulder_radius_um_exact']))
    check(roi==interval(candidate['ROI_um_exact']),'fixed common observer ROI changed')
    graph={}
    for label,index in (('top','apical'),('bottom','basal')):
        radial,slope,e=face_error(faces[label],epsilon,roi)
        actual=candidate['graph_transfer'][index]
        check(radial-8*epsilon==Q(actual['radial_derivative_margin_exact'])
              and slope==Q(actual['base_graph_Lipschitz_upper_exact'])
              and e==Q(actual['Eulerian_vertical_error_um_exact'])
              and interval(actual['ROI_um_exact'])==roi,'independent per-face coefficient or error bound mismatch')
        graph[index]={'stored_radial_local_lower':str(radial),'stored_graph_slope_upper':str(slope),
                      'C_um_per_dimensionless_exact':str(L*Q(5,4)*(1+slope)),'error_um_exact':str(e)}
    ea,eb=(Q(graph[k]['error_um_exact']) for k in ('apical','basal'))
    record={'pilot_id':pid,'endpoint_id':stage,'role':field['role'],'field_sha256':field['field_sha256'],
            'posterior_exact':str(err),'binary_addition_gap_exact':str(xi),'nodal_epsilon_um_exact':str(epsilon),
            'fixed_R_trace_dofs':dofs,'per_face':graph,'ROI_um_exact':as_text(roi)}
    if not is_small:
        check(field['role']=='large_final','unexpected target role')
        stored=interval(verification['surface']['Q_um_exact'])
        endpoint=one(native_primary['endpoints'],lambda x:x['endpoint_id']==stage)
        check(stored==interval(endpoint['large_Q_um_exact'])==interval(candidate['stored_Q_um_exact'])
              and Q(verification['surface']['shoulder_radius_um_exact'])==roi[1],'stored large Q/shoulder changed')
        root=(stored[0]-epsilon-ea,stored[1]+epsilon+ea)
        check(root==interval(candidate['cached_root_Q_um_exact']) and epsilon+ea==Q(candidate['Q_error_upper_um_exact']),
              'large fixed-axis and current-radius inflation mismatch')
        return {**record,'stored_Q_um_exact':as_text(stored),'root_Q_um_exact':as_text(root)},[]
    all_sweep=docs[RUN/'observer_robustness/actual/results.json']
    verified_sweep=docs[RUN/'observer_robustness/independent_actual/receipt.json']
    small_path=f"{(RUN/'observer_robustness/actual').relative_to(ROOT).as_posix()}/{stage}/{field['checkpoint']}.json"
    meta=one(all_sweep['field_observer_receipts'],lambda x:x['path']==small_path)
    check(verified_sweep['input_and_output_sha256'][small_path]==meta['sha256'],'sweep data lack independent binding')
    sweep_field=bind(ROOT/small_path,meta['sha256'],links)
    check(all(sweep_field[k]==field[k] for k in ('field','field_sha256','verification','verification_sha256','step')),
          'small stored observer intervals from another field')
    errors=dict(zip(OBS,(ea+eb,ea+eb,ea,eb),strict=True))
    def from_packed(x):
        return Q(x['numerator'],x['denominator'])
    oc=sweep_field['observer_contract']
    check(tuple(map(from_packed,oc['current_radius_ROI_um']))==roi
          and from_packed(oc['reference_apical_um'])==height and from_packed(oc['reference_basal_um'])==0,
          'same static registered ROI and flat reference required')
    comparisons=[]
    for variant in VARIANTS:
        for name in OBS:
            old=one(candidate['root_observer_comparisons'],lambda w:w['variant']==variant and w['observable']==name)
            saved=interval(sweep_field['variant_observables_um_exact'][variant][name])
            observed=Q(step['observed_exact'][name])
            error=errors[name]
            root,absolute,verdict=required_bound(saved,error,observed,delta)
            check(saved==interval(old['saved_um_exact']) and observed==Q(old['recorded_value_um_exact'])
                  and error==Q(old['observer_error_upper_um_exact']) and root==interval(old['cached_root_um_exact'])
                  and absolute==interval(old['absolute_residual_um_exact']) and verdict==old['status_at_frozen_delta'],
                  f'exact static observer comparison mismatch: {stage}, {variant}, {name}')
            comparisons.append({'endpoint_id':stage,'variant':variant,'observable':name,
                                'root_interval_um_exact':as_text(root),'residual_um_exact':as_text(absolute),
                                'status':verdict})
    return record,comparisons


def main():
    links={}
    docs={p:bind(p,h,links) for p,h in PINS.items()}
    app=docs[APP/'receipt.json']
    result=docs[APP/'results.json']
    protocol=docs[APP/'protocol.json']
    prep=docs[APP/'preflight.json']
    check(app['results_sha256']==PINS[APP/'results.json'] and app['protocol_sha256']==PINS[APP/'protocol.json']
          and app['preflight_sha256']==PINS[APP/'preflight.json'] and app['cached_root_review_sha256']==PINS[CACHED_REPLAY/'receipt.json']
          and app['geometry_review_sha256']==PINS[GEOMETRY_REPLAY/'receipt.json']
          and app['exact_spatial_integral_roots_claimed'] is False and app['new_AD_force_Hessian_or_solve'] is False
          and prep['status']=='PASS_PREPARATION_ONLY','candidate receipt authority mismatch')
    bind(APP/'apply.py',app['source_sha256'],links,parse=False)
    bind(APP/'PROOF.md',app['proof_markdown_sha256'],links,parse=False)
    for source in (protocol['inputs_sha256'],result['inputs_sha256']):
        for name,digest in source.items():
            bind(ROOT/name,digest,links,parse=False)
    for d in (CACHED_REPLAY,GEOMETRY_REPLAY):
        OSfile=d/'os_exit_observation.json'
        obs=bind(OSfile,protocol['inputs_sha256'][OSfile.relative_to(ROOT).as_posix()],links)
        tfile=d/'terminal.json'
        terminal=bind(tfile,obs['terminal_sha256'],links)
        check(terminal['status']=='FINISHED' and type(obs['observed_exit_code']) is int
              and obs['observed_exit_code']==terminal['exit_code']==0 and obs['pid']==terminal['pid']
              and terminal['receipt_sha256']==obs['receipt_sha256']==PINS[d/'receipt.json'],
              'prior independent true OS exit chain differs')
    primary=docs[RUN/'preparation/protocol.json']
    delta=Q(primary['selection']['conditional_source_allowance_um_exact'])
    check(delta==Q(result['frozen_model_discrepancy_allowance_um_exact'])==Q(app['same_saved_allowance_exact']),
          'unchanged frozen conditional discrepancy required')
    expected=tuple(f'{p}_{s}' for p in ENDPOINTS for s in ('small_17','large_03'))
    check(tuple(p['pilot_id'] for p in result['four_pilots'])==expected,'all four ordered piloted states required')
    records,rows=[],[]
    for candidate in result['four_pilots']:
        record,local=compare_pilot(candidate,docs,links,delta)
        records.append(record)
        rows.extend(local)
    counts=dict(Counter(r['status'] for r in rows))
    check(len(rows)==result['L_static_observer_comparisons']==32 and counts==result['comparison_counts']=={'IN':28,'OUT':4}
          and all(r['variant'].startswith('current_radial_line__') and r['observable']=='region_averaged_shortening_um'
                  for r in rows if r['status']=='OUT'),'wrong all32 static classifications')
    large=[r for r in records if r['role']=='large_final']
    def difference(xs):
        a,b=map(interval,xs)
        return a[0]-b[1],a[1]-b[0]
    old_gap=difference(r['stored_Q_um_exact'] for r in large)
    root_gap=difference(r['root_Q_um_exact'] for r in large)
    shown=interval(json.loads(app['positive_cached_quadrature_gap_outward_um'], parse_float=Q))
    check(old_gap==interval(result['stored_large_signed_gap_um_exact'])
          and 0<shown[0]<=root_gap[0]<=root_gap[1]<=shown[1]
          and root_gap==interval(result['cached_root_large_signed_gap_um_exact']),
          'cached root large signed gap or outward decimal bounds differ')
    with (HERE/'receipt.json').open('x',encoding='utf-8',newline='\n') as f:
        json.dump({'status':'PASS_INDEPENDENT_CACHED_ROOT_OBSERVER_APPLICATION_READBACK','scope':SCOPE,
                   'source_sha256':sha(Path(__file__).read_bytes()),'input_sha256':links,
                   'cached_operator':'exact-real frozen binary-cache quadrature at exact nonrounded physical strengths',
                   'exact_large_gap_um':as_text(root_gap),'comparison_counts':counts,
                   'conditional_delta_um_exact':str(delta),'four_checked_fields':records,'all_32_comparisons':rows,
                   'reused_reviewed_mathematical_premises':['cached independent root enclosures',
                       'independent same-box body embedding','reviewed P2 nodal-to-graph theorem',
                       'independently reviewed saved observer vectors'],
                   'producer_application_functions_called':False,'new_force_or_AD':False,
                   'new_geometry_certifications':False,'exact_spatial_roots_certified':0,
                   'tracked_small_root_fields':2,'full_goal_claimed':False},f,indent=2)
        f.write('\n')
    print(json.dumps({'status':'PASS','receipt_sha256':sha((HERE/'receipt.json').read_bytes()),
                      'counts':counts,'large_gap_display':[float(v) for v in root_gap]}),flush=True)


if __name__=='__main__':
    main()

