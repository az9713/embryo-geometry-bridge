"""Independent exact rational observer readback from the four spatial proof centers."""
import argparse
from base64 import b64decode
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import re

import numpy as np


HERE=Path(__file__).resolve().parent
APP=HERE.parent
EXACT=APP.parent
DESIGN=EXACT.parent
RUN=DESIGN.parent
ROOT=RUN.parents[1]
SPATIAL=EXACT/'bridge_reduction_v1/independent_actual_preparation'
PRIMARY=RUN/'independent_actual_preparation/actual_replay/receipt.json'
SECONDARY=RUN/'observer_robustness/independent_actual/receipt.json'
OLD=DESIGN/'cached_pilots_v1/observer_application_v1/independent_readback/receipt.json'
REVIEW=APP/'independent_source_review/receipt.json'
SOURCE=ROOT/'scripts/transfer_ambiguity_exact_four_observers.py'
PHASE={'source_sha256':'4779f97e37fdb7b4bd6cb66dc5a7934ebd6adf12f5c007ea88e0c26244845cb1',
       'protocol_sha256':'080729ec73de6d92d1f2e4ceeaa0ac748bf252a3c56295b19eeccb8338bb3a14',
       'preflight_sha256':'43f518b1c2fcefee7812214e02c8bcd1a7995189603783cd14582e44d5feb3d9',
       'independent_source_review_sha256':'2411ec6b9bc0b762af3b0ec617138c8748e487eb114693c58ce2133a4937c684'}
RADIUS=Q(1,2**34)
LENGTH=Q(67,2)
ENDPOINTS=('max_contrast_change','max_negative_contrast_change')
CASES=tuple(p+'_'+role for p in ENDPOINTS for role in ('small_17','large_03'))
VARIANTS=('current_disk__spatial_max','current_disk__regional_mean',
          'current_radial_line__spatial_max','current_radial_line__regional_mean')
OBS=('maximum_shortening_um','region_averaged_shortening_um','apical_max_shortening_um','basal_max_shortening_um')


def check(condition,msg):
    if not condition:
        raise ValueError(msg)


def name(path):
    return path.relative_to(ROOT).as_posix()


def digest(raw):
    return sha256(raw).hexdigest()


def read(path,h):
    raw=path.read_bytes()
    check(digest(raw)==h,f'changed pinned bytes: {path}')
    return json.loads(raw)


def exact(x):
    check(type(x) in (str,int,Q),'strict exact rational required')
    return Q(x)


def pair(w):
    a,b=map(exact,w)
    check(a<=b,'invalid directed interval')
    return a,b


def one(rows,test):
    row,=(r for r in rows if test(r))
    return row


def admission(os,expected_result):
    session=os.get('unified_exec_session')
    chunk=os.get('unified_exec_chunk_id')
    process_id=(type(session) is int and session>0 or type(chunk) is str and re.fullmatch('[0-9a-f]{6,64}',chunk))
    check(os['status']=='PASS_ROOT_OBSERVED_EXACT_FOUR_OBSERVER_EXIT' and process_id
          and type(os['observed_exit_code']) is int and os['observed_exit_code']==0
          and os['results_sha256']==expected_result and all(os[k]==v for k,v in PHASE.items()),
          'external successful producer OS observation required')


def build_traces(field,box,packet,identity,contract):
    payload=packet['frozen_binary_cache']
    raw=b64decode(payload['data'],validate=True)
    check(len(raw)==payload['bytes'] and digest(raw)==payload['sha256']==identity['embedded_cache_sha256'],
          'binary point-cache mismatch')
    with np.load(BytesIO(raw),allow_pickle=False) as state:
        R=state['preconditioner']
        check(R.dtype.str=='<f8' and R.shape==(611,611) and np.isfinite(R).all()
              and digest(R.tobytes(order='C'))==identity['R_binary_sha256'], 'same exact binary R required')
        X=[tuple(map(Q,row)) for row in state['reference'].tolist()]
        U=tuple(map(Q,state['u'].tolist()))
        ix=state['free'].tolist()
        tri=state['triangles'].tolist()
        trace_cache=[{tuple(row) for row in state[f'surface_{face}_edges'].tolist()} for face in range(2)]
    params=field['parameters']; rim=Q(params['outer_radius']); height=Q(params['height'])
    stable_free=[2*i+d for i,(r,_) in enumerate(X) for d in (0,1) if r!=rim and (r!=0 or d==1)]
    check(len(stable_free)==len(ix)==len(U)==611 and all(type(k) is int for k in ix)
          and ix==stable_free==box['free'] and len(tri)==144 and tri==box['triangles']==field['triangles']
          and len(X)==325 and X==[tuple(map(Q,row)) for row in field['reference_nodes']]
          ==[tuple(map(exact,row)) for row in box['reference_nodes_exact']]
          and U==tuple(map(Q,field['free_displacement'])) and exact(box['corrector_radius_exact'])==RADIUS
          and box['root_membership_assumed'] is False and box['any_spatial_derivative_bound_assumed'] is False
          and box['record_id']==identity['record_id'] and box['point_cache_sha256']==identity['cache_packet_sha256']
          and box['model_observation_contract_sha256']==identity['model_contract_sha256']
          and box['parameter_cell_sha256']==contract['parameter_cell_sha256'],'different exact FE center, map or cube')
    inc=dict(zip(ix,U,strict=True))
    S=[tuple(map(Q,row)) for row in field['current_nodes']]
    gap=Q(0)
    for i,(xval,wrapped,ctr,rads,stored) in enumerate(zip(X,box['coordinate_hulls_exact'],box['current_centers_exact'],
             box['coordinate_radii_exact'],S,strict=True)):
        for d,(x,intr,c,rad,s) in enumerate(zip(xval,wrapped,ctr,rads,stored,strict=True)):
            index=2*i+d; y=x+inc.get(index,0); r=RADIUS if index in inc else Q(0)
            check(y==exact(c) and r==exact(rad) and pair(intr)==(y-r,y+r) and Q(float(y))==s,
                  'one of 650 exact center and whole-cube nodal hull links differs')
            gap=max(gap,abs(y-s))
    faces={}
    for label,z,prior in zip(('bottom','top'),(Q(0),height),trace_cache,strict=True):
        edges=set()
        for a,b,c,ab,bc,ca in tri:
            for i,j,m in ((a,b,ab),(b,c,bc),(c,a,ca)):
                if X[i][1]==X[j][1]==z:
                    edges.add((i,m,j) if X[i][0]<X[j][0] else (j,m,i))
        edge_order=sorted(edges,key=lambda trio:X[trio[0]][0])
        check(edges==prior and len(edge_order)==12 and all(a[-1]==b[0] for a,b in zip(edge_order,edge_order[1:]))
              and 2*edge_order[0][0] not in inc and 2*edge_order[-1][-1] not in inc
              and S[edge_order[0][0]][0]==0 and S[edge_order[-1][-1]][0]==rim,
              'fixed-radius endpoints and connected P2 trace needed')
        points=[[S[k] for k in edge] for edge in edge_order]
        checks=[[tuple(map(Q,node)) for node in segment] for segment in field[f'boundary_{label}_quadratic_controls']]
        check(points==checks,'stored observer geometry differs from the same serialized nodal field')
        faces[label]=[tuple(tuple(LENGTH*v for v in p) for p in points_on_edge) for points_on_edge in points]
    return faces,gap,LENGTH*height


def face_bound(trace,eps,roi):
    area_min=[]; lipschitz=[]
    check(trace[0][0][0]<=roi[0]<roi[1]<=trace[-1][-1][0],'fixed-radius ROI not covered')
    for p,m,q in trace:
        a=tuple(-3*v+4*w-x for v,w,x in zip(p,m,q,strict=True))
        b=tuple(v-4*w+3*x for v,w,x in zip(p,m,q,strict=True))
        low=min(a[0],b[0])
        check(low>8*eps,'insufficient strict radial derivative after spatial root error')
        area_min.append(low)
        lipschitz.append(max(abs(a[1]),abs(b[1]))/low)
    radial,slope=min(area_min),max(lipschitz)
    return radial-8*eps,slope,(eps+eps/4)*(1+slope)


def enclose(saved,change):
    check(change>=0,'nonnegative validated observer error required')
    left,right=pair(saved)
    return left-change,right+change


def residual(root,recorded,delta):
    start,end=root
    lower=recorded-end if recorded>end else start-recorded if recorded<start else Q(0)
    upper=max(recorded-start,end-recorded)
    status='OUT' if lower>delta else 'IN' if upper<=delta else 'UNRESOLVED'
    return lower,upper,status


def reproduce(read_doc,hashes,protocol):
    four=read_doc(SPATIAL/'all_four_readback.json')
    tm=read_doc(EXACT/'preparation/target_manifest.json')
    fm=read_doc(RUN/'field_manifest.json')
    geo=read_doc(EXACT/'independent_actual/replay/receipt.json')
    budget=read_doc(DESIGN/'budgets_v1/actual/results.json')
    primary=read_doc(PRIMARY)
    secondary=read_doc(SECONDARY)
    old=read_doc(OLD)
    delta=exact(protocol['allowance_um_exact'])
    check(delta>0 and delta==exact(budget['conditional_source_allowance_um_exact'])
          ==exact(old['conditional_delta_um_exact'])==exact(secondary['source_configured_delta_um_exact'])
          and fm['independent_primary_review_sha256']==hashes[name(PRIMARY)]
          and tuple(r['pilot_id'] for r in four['cases'])==CASES==tuple(tm['pilot_order'])
          and len(four['current_file_sha256'])==3006
          and all(hashes[k]==v for k,v in four['current_file_sha256'].items()),'changed primary scope or proof closure')
    outputs=[]; traces=[]
    for case in four['cases']:
        pid=case['pilot_id']; identity=case['target_identity']; root=case['selfmap']
        native=read_doc(ROOT/case['independent_proof_file'])
        old_p=one(old['four_checked_fields'],lambda p:p['pilot_id']==pid)
        t=one(tm['targets'],lambda e:e['pilot_id']==pid)
        contract=read_doc(ROOT/t['model_contract'])
        fid=one(fm['fields'],lambda e:e['field']==contract['source_field'])
        c=read_doc(ROOT/contract['point_cache'])
        g=one(geo['geometry_records'],lambda e:e['pilot_id']==pid and exact(e['radius_exact'])==RADIUS)
        gp=read_doc(ROOT/g['receipt'])
        check(hashes[t['model_contract']]==t['model_contract_sha256']==identity['model_contract_sha256']
              and identity['field_sha256']==fid['field_sha256']==contract['source_field_sha256']==c['input_state_sha256']
              ==primary['all_hashes'][fid['field']]==hashes[fid['field']]==old_p['field_sha256']
              and fid['verification_sha256']==primary['all_hashes'][fid['verification']]==hashes[fid['verification']]
              ==contract['source_geometry_verification_sha256']==c['field_verification_sha256']
              and hashes[contract['point_cache']]==identity['cache_packet_sha256']==contract['point_cache_sha256']
              and identity['R_binary_sha256']==contract['binary_preconditioner_sha256']
              and identity['energy_scale_binary_exact']==contract['energy_scale_binary_exact']==c['normalization_binary_exact']
              and case['independent_proof_sha256']==hashes[case['independent_proof_file']]
              and native['posterior_below_observer_state_design_goal'] is True
              and native['target_identity']==identity and native['conservative_compared_bounds']==root
              and gp['status']=='PASS_FRESH_SINGLETON_WHOLE_GEOMETRY_ONLY'
              and gp['proof_center_field_sha256']==fid['field_sha256']
              and g['sha256']==hashes[g['receipt']]==native['same_geometry_admissibility_sha256']
              and gp['box_sha256']==g['box_sha256']==hashes[t['box_1']]==t['box_1_sha256'],
              'different field, point cache, R, normalization, or reviewed whole geometry')
        dose=exact(fid['step']['dose_exact']); gain=exact(c['common_gain_exact']); drive=dose*gain
        check(dose==exact(identity['dose_exact'])==exact(c['optical_dose_exact'])
              and gain==exact(identity['gain_exact'])
              and drive==exact(fid['step']['requested_strength_exact'])==exact(identity['drive_exact'])
              ==exact(c['activation_strength_exact'])==exact(contract['exact_requested_drive']['activation_strength_exact']),
              'exact gain times dose differs')
        e,z,r=map(exact,(root['eta_upper_exact'],root['z_upper_exact'],case['radius_exact']))
        check(0<=e and 0<=z<1 and r==RADIUS and root['strict_selfmap'] is True
              and e+z*r==exact(root['selfmap_lhs_upper_exact'])<r, 'strict fixed-input spatial root witness failed')
        posterior=e/(1-z)
        check(posterior==exact(root['posterior_error_upper_exact'])<r,
              'different spatial posterior bound')
        field=read_doc(ROOT/fid['field']); verified=read_doc(ROOT/fid['verification'])
        faces,xi,height=build_traces(field,read_doc(ROOT/t['box_1']),c,identity,contract)
        budget_row=one(budget['field_budgets'],lambda j:j['field']==fid['field'])
        roi=pair(budget_row['fixed_roi_um_exact']); eps=LENGTH*(posterior+xi)
        check(xi==exact(old_p['binary_addition_gap_exact'])
              ==exact(budget_row['center_and_model']['node_addition_gap_dimensionless_exact'])
              and roi==pair(old_p['ROI_um_exact']) and posterior<exact(budget_row['e_goal_dimensionless_exact'])
              and (fid['endpoint_id'],fid['role'])==(old_p['endpoint_id'],old_p['role']),
              'different exact rounding gap, physical role or frozen observer budget')
        (ra,ma,ea),(rb,mb,eb)=(face_bound(faces[face],eps,roi) for face in ('top','bottom'))
        for side,slope in (('apical',ma),('basal',mb)):
            check(slope==exact(old_p['per_face'][side]['stored_graph_slope_upper']),
                  'newly extracted exact graph slope differs from independently reviewed stored curve')
        row={'pilot_id':pid,'role':fid['role'],'endpoint_id':fid['endpoint_id'],
             'spatial_root_proof_sha256':case['independent_proof_sha256'],'posterior_error_exact':str(posterior),
             'xi_exact':str(xi),'physical_nodal_error_um_exact':str(eps),
             'new_radial_margins_exact':{'apical':str(ra),'basal':str(rb)}}
        if fid['role']=='large_final':
            full=one(primary['endpoints'],lambda y:y['endpoint_id']==fid['endpoint_id'])
            q=pair(verified['surface']['Q_um_exact'])
            check(fid['checkpoint']=='large_03' and q==pair(full['large_Q_um_exact'])==pair(old_p['stored_Q_um_exact'])
                  and exact(verified['surface']['shoulder_radius_um_exact'])==roi[1] and roi[0]==0,
                  'different large reference, current radius or stored Q')
            row.update(stored_Q_um_exact=list(map(str,q)),Q_error_um_exact=str(eps+ea),
                       exact_spatial_root_Q_um_exact=list(map(str,enclose(q,eps+ea))))
        else:
            check(fid['role']=='small_record' and fid['checkpoint']=='small_17', 'only the original two L records admitted')
            spath=RUN/'observer_robustness/actual'/fid['endpoint_id']/'small_17.json'
            sweep=read_doc(spath)
            check(hashes[name(spath)]==secondary['input_and_output_sha256'][name(spath)]
                  and all(sweep[k]==fid[k] for k in ('field','field_sha256','verification','verification_sha256','step'))
                  and sweep['record_id']=='optoGEF2_small_col_L','different independently reviewed L observer data')
            oc=sweep['observer_contract']; unbox=lambda x:Q(x['numerator'],x['denominator'])
            check(tuple(map(unbox,oc['current_radius_ROI_um']))==roi and unbox(oc['reference_apical_um'])==height
                  and unbox(oc['reference_basal_um'])==0,'different fixed reference profiles or ROI')
            comparisons=[]
            for v in VARIANTS:
                for o,error in zip(OBS,(ea+eb,ea+eb,ea,eb),strict=True):
                    vals=enclose(sweep['variant_observables_um_exact'][v][o],error)
                    lo,hi,status=residual(vals,exact(fid['step']['observed_exact'][o]),delta)
                    comparisons.append({'variant':v,'observable':o,'root_interval_um_exact':list(map(str,vals)),
                                        'absolute_residual_um_exact':list(map(str,(lo,hi))),'status':status})
            row['comparisons']=comparisons
        outputs.append(row)
        traces.append({'pilot_id':pid,'serialized_addition_gap_exact':str(xi),'fixed_ROI_um_exact':list(map(str,roi)),
                       'graph_slope_upper_exact':{'apical':str(ma),'basal':str(mb)},
                       'vertical_error_um_exact':{'apical':str(ea),'basal':str(eb)},'whole_cube_nodal_hulls_checked':650})
    return outputs,traces


def report(result_sha,exit_sha):
    protocol=read(APP/'protocol.json',PHASE['protocol_sha256'])
    preflight=read(APP/'preflight.json',PHASE['preflight_sha256'])
    review=read(REVIEW,PHASE['independent_source_review_sha256'])
    check(digest(SOURCE.read_bytes())==PHASE['source_sha256']==protocol['source_sha256']==preflight['source_sha256']
          and review['status']=='PASS_INDEPENDENT_EXACT_FOUR_OBSERVER_SOURCE'
          and all(review[k]==PHASE[k] for k in ('source_sha256','protocol_sha256','preflight_sha256')),
          'unchanged externally reviewed observer phase required')
    admission(read(APP/'observed_os_exit.json',exit_sha),result_sha)
    actual=read(APP/'results.json',result_sha)
    hashes=protocol['inputs_sha256']
    extras={name(APP/'protocol.json'):PHASE['protocol_sha256'],name(APP/'preflight.json'):PHASE['preflight_sha256'],
            name(REVIEW):PHASE['independent_source_review_sha256']}
    check(actual['status']=='APPLIED_FOUR_INDEPENDENT_SPATIAL_ROOTS_AWAITING_READBACK'
          and actual['bindings']=={k:PHASE[k] for k in ('source_sha256','protocol_sha256','preflight_sha256')}
          and actual['inputs_sha256']=={**hashes,**extras} and actual['scope']==protocol['scope']
          and actual['full136_calibration_claimed'] is False and actual['source_fit_PDE_loading_or_novelty_claimed'] is False,
          'producer authority differs from the source reviewed phase')
    independent,trace=reproduce(lambda p:read(p,hashes[name(p)]),hashes,protocol)
    for row,prod in zip(independent,actual['cases'],strict=True):
        lhs=dict(row); rhs=dict(prod)
        if 'comparisons' in lhs:
            def keyed(rows):
                return {(x['variant'],x['observable']):x for x in rows}
            a,b=lhs.pop('comparisons'),rhs.pop('comparisons')
            check(len(a)==len(b)==len(keyed(a))==len(keyed(b))==16 and keyed(a)==keyed(b),
                  'one of 16 independent exact L observer rows disagrees')
        check(lhs==rhs,'independent root, observer or geometry scalars disagree')
    two={r['endpoint_id']:pair(r['exact_spatial_root_Q_um_exact']) for r in independent if r['role']=='large_final'}
    p,m=(two[x] for x in ENDPOINTS); gap=p[0]-m[1],p[1]-m[0]
    counts=Counter(x['status'] for w in independent for x in w.get('comparisons',[]))
    check(list(map(str,gap))==actual['large_signed_gap_um_exact'] and actual['strict_positive_gap'] is (gap[0]>0)
          and dict(counts)==actual['small_comparison_counts'] and sum(counts.values())==32,'global signed gap or count mismatch')
    output={'status':'PASS_INDEPENDENT_EXACT_FOUR_SPATIAL_OBSERVER_READBACK',**PHASE,
            'independent_source_sha256':digest(Path(__file__).read_bytes()),'producer_results_sha256':result_sha,
            'producer_observed_os_exit_sha256':exit_sha,'spatial_proof_closure_sha256':hashes[name(SPATIAL/'all_four_readback.json')],
            'independently_matched_spatial_cases':4,'independently_matched_L_alternative_comparisons':32,
            'comparison_counts':dict(counts),'large_signed_gap_um_exact':list(map(str,gap)),'strict_positive_gap':gap[0]>0,
            'conditional_allowance_um_exact':protocol['allowance_um_exact'],'cases':independent,'trace_checks':trace,
            'producer_functions_executed':0,'new_mechanics_or_geometry_certification_jobs':0,
            'all_136_primary_comparisons_claimed':False,'source_compatibility_or_continuum_or_novelty_claimed':False,
            'scope':'Exact arithmetic reconstruction with independent state, P2 trace, signed difference and observed residual calculations. Shared frozen P2 graph theorem, stored-field observer intervals and independently certified spatial roots.'}
    path=HERE/'receipt.json'
    with path.open('x',encoding='utf-8') as f:
        json.dump(output,f,indent=2)
        f.write('\n')
    print(json.dumps({'status':output['status'],'receipt_sha256':digest(path.read_bytes()),
                      'count':dict(counts),'gap_um_display':list(map(float,gap))}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--producer-results-sha256',required=True)
    parser.add_argument('--producer-os-sha256',required=True)
    args=parser.parse_args()
    report(args.producer_results_sha256,args.producer_os_sha256)
