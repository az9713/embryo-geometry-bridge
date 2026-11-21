"""Independent composition arithmetic, with disclosed reviewed-component replays.

The body and baseline static observers reuse previously reviewed implementations;
predictor/rounding radii, boundary extraction, transport and residual arithmetic
are rebuilt separately. Finalization also requires independent root capsules.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,base64,hashlib,io,json,time
import numpy as np
from root_body_geometry import certify_box
from quadratic_endpoint_observer import observe_state
from joint_calibration_contract import OBSERVABLES

ROOT=Path(__file__).resolve().parents[1]
PROD=ROOT/'runs/cohort_static_observers_v1'
RESULT=PROD/'composition_final/results.json'
RUN=ROOT/'runs/cohort_observer_review_v1'
read=lambda p:json.loads(Path(p).read_text(encoding='utf-8'))
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
rel=lambda p:Path(p).relative_to(ROOT).as_posix()
q=lambda v:Q.from_float(float(v)) if isinstance(v,(float,np.floating)) else Q(v)
ORDER=('apical_area_modulus','basal_area_modulus','apical_surface_tension','basal_surface_tension','common_gain')

def save(p,v):
    if p.exists():raise RuntimeError('preserve review artifact '+str(p))
    p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8')

def boundaries(state):
    ref=state['reference_nodes'];edges={}
    for a,b,c,d,e,f in state['triangles']:
        for i,j,m in ((a,b,d),(b,c,e),(c,a,f)):edges[tuple(sorted((i,j)))]=m
    out={}
    for name,level in [('top',q(state['parameters']['height'])),('bottom',Q(0))]:
        selected=[]
        for (i,j),m in edges.items():
            if q(ref[i][1])==q(ref[j][1])==level:
                if q(ref[i][0])>q(ref[j][0]):i,j=j,i
                selected.append((q(ref[i][0]),[i,m,j]))
        out[name]=[nodes for _,nodes in sorted(selected)]
    return out

def arithmetic():
    start=time.time();p=read(PROD/'protocol.json');result=read(RESULT);selection=read(PROD/'provider_selection.json')
    context=read(ROOT/p['root_protocol_file']);contract=read(PROD/'model_observation_contract.json')
    files={rel(RESULT):sha(RESULT),rel(PROD/'protocol.json'):sha(PROD/'protocol.json'),rel(PROD/'provider_selection.json'):sha(PROD/'provider_selection.json'),rel(PROD/'model_observation_contract.json'):sha(PROD/'model_observation_contract.json'),rel(Path(__file__)):sha(__file__)}
    files.update(p['source_hashes']);files.update(p['input_hashes'])
    for key,h in files.items():
        if sha(ROOT/key)!=h:raise RuntimeError('changed frozen input '+key)
    rows={v['record_id']:v for v in context['records']};providers={v['record_id']:v for v in selection['records']}
    audits=[];residuals=[];missing=[]
    from csv import DictReader
    with (ROOT/'data/derived/calibration_small_18um.csv').open(encoding='utf-8-sig',newline='') as f:observations={r['record_id']:r for r in DictReader(f)}
    for entry in result['records']:
        rid=entry['record_id'];row=rows[rid];provider=providers[rid];rfile=ROOT/provider['root_result_file'];root=read(rfile)
        outfile=ROOT/entry['file'];o=read(outfile);files[rel(rfile)]=sha(rfile);files[rel(outfile)]=sha(outfile)
        checks={'observer_record_bytes_bound':sha(outfile)==entry['sha256'],'root_record_bytes_bound':sha(rfile)==provider['root_result_sha256'],
            'same_record_dose_cell_contract':o['record_id']==rid and root['record_id']==rid and Q(root['dose_exact'])==Q(row['dose_exact']) and root['parameter_cell_sha256']==result['parameter_cell_sha256'] and root['model_observation_contract_sha256']==result['common_contract_sha256']}
        if root['status']!='CERTIFIED_UNIFORM_LOCAL_ROOT':
            checks['unresolved_root_never_receives_invented_endpoints']=o['status']=='UNRESOLVED' and 'prediction_intervals_exact' not in o
            missing.extend([[rid,k] for k in OBSERVABLES]);audits.append({'record_id':rid,'status':'EXPECTED_UNRESOLVED','checks':checks,'root_certificate_sha256':sha(rfile),'observer_certificate_sha256':sha(outfile)})
            continue
        state=read(ROOT/row['saved_state']);packet=read(ROOT/root['point_cache_file'])
        files[root['point_cache_file']]=sha(ROOT/root['point_cache_file']);files[root['predictor_file']]=sha(ROOT/root['predictor_file'])
        raw=base64.b64decode(packet['frozen_binary_cache']['data']);assert hashlib.sha256(raw).hexdigest()==packet['frozen_binary_cache']['sha256']
        with np.load(io.BytesIO(raw),allow_pickle=False) as z:reference=z['reference'].copy();free=list(map(int,z['free']));u=z['u'].copy()
        with np.load(ROOT/root['predictor_file'],allow_pickle=False) as z:S=z['S'].copy()
        selected=root['uniform_cell_result'];radius=Q(selected['selected_radius_exact'])
        C=root['parameter_cell'];delta=[(Q(C[k][0])-Q(context['base_parameters_exact'][k]),Q(C[k][1])-Q(context['base_parameters_exact'][k])) for k in ORDER]
        hulls=[]
        for coefficients in S:
            products=[(q(s)*a,q(s)*b) for s,(a,b) in zip(coefficients,delta)]
            hulls.append((sum(min(v) for v in products),sum(max(v) for v in products)))
        lookup={k:i for i,k in enumerate(free)};radii=[];gaps=[]
        for k,v in enumerate(reference.ravel()):
            exact=q(v)+(q(u[lookup[k]]) if k in lookup else Q(0));gap=abs(exact-q(state['current_nodes'][k//2][k%2]));gaps.append(gap)
            radii.append(gap+(max(map(abs,hulls[lookup[k]]))+radius if k in lookup else Q(0)))
        checks['full_corrector_radius_not_posterior']=Q(o['full_corrector_radius_exact'])==radius
        checks['all611_predictor_hulls_exact']=list(map(lambda v:tuple(map(Q,v)),o['predictor_hulls_exact']))==hulls
        checks['all_serialization_gaps_and_nodal_radii_exact']=[Q(v) for pair in o['serialization_gaps_exact'] for v in pair]==gaps and [Q(v) for pair in o['coordinate_radii_normalized_exact'] for v in pair]==radii
        exact_ref=[[q(v) for v in pair] for pair in reference];centers=[[q(v) for v in pair] for pair in state['current_nodes']]
        body=certify_box(exact_ref,centers,state['triangles'],[radii[i:i+2] for i in range(0,len(radii),2)])
        checks['reviewed_whole_body_component_replay']=body==o['body_geometry'] and body['status']=='CERTIFIED_ENTIRE_NODAL_BOX_EMBEDDING'
        edge=boundaries(state);L=Q(contract['physical_length_um_exact']);roi=list(map(Q,contract['roi_um_exact']));errors={}
        for name,nodes in edge.items():
            eps=L*max(radii[2*n+j] for triple in nodes for n in triple for j in (0,1));minimum=None;lip=Q(0)
            for i,m,j in nodes:
                points=[tuple(L*q(v) for v in state['current_nodes'][n]) for n in (i,m,j)]
                linear=tuple(-3*points[0][k]+4*points[1][k]-points[2][k] for k in (0,1));quadratic=tuple(2*points[0][k]-4*points[1][k]+2*points[2][k] for k in (0,1))
                der=min(linear[0],linear[0]+2*quadratic[0]);minimum=der if minimum is None else min(minimum,der)
                if der<=0:raise RuntimeError('nonmonotone base graph')
                lip=max(lip,max(abs(linear[1]),abs(linear[1]+2*quadratic[1]))/der)
            endpoint_nodes=(nodes[0][0],nodes[-1][-1])
            checks[name+'_fixed_radial_endpoints_independently_extracted']=all(2*n not in lookup and radii[2*n]==0 for n in endpoint_nodes)
            E=Q(5,4)*eps*(1+lip);errors[name]=E;transfer=o['surface_error_transfer'][name]
            checks[name+'_transport_exact_and_strict_margin']=minimum>8*eps and E==Q(transfer['Eulerian_vertical_error_um_exact']) and eps==Q(transfer['nodal_coordinate_error_um_exact']) and lip==Q(transfer['base_graph_Lipschitz_upper_exact']) and roi==list(map(Q,transfer['ROI_um_exact']))
        ea,eb=errors['top'],errors['bottom'];widths=dict(zip(OBSERVABLES,(ea+eb,ea+eb,ea,eb)))
        checks['endpoint_transport_assignments_exact']=all(Q(o['endpoint_error_um_exact'][k])==v for k,v in widths.items())
        baseline=observe_state(state,roi_um=roi,length_unit_um=L,reference_apical_um=L,reference_basal_um=0,weighting='current_disk',component_aggregation='spatial_max',tolerance_um=p['static_endpoint_tolerance_um_exact'])
        values=(baseline['four_endpoints']['maximum_Eulerian_height_shortening'],baseline['mean_alternatives']['current_disk']['height_regional_mean_shortening'],baseline['spatial_max_components']['apical_downward'],baseline['spatial_max_components']['basal_upward'])
        checks['reviewed_baseline_observer_component_replayed']=all(tuple(map(Q,o['base_endpoints_exact'][k]))==(v.lo,v.hi) for k,v in zip(OBSERVABLES,values))
        checks['base_reporting_certificate_status_and_actual_width']=baseline['status']=='CERTIFIED_STATIC_INTERVALS' and Q(o['base_maximum_reporting']['actual_interval_width_um_exact'])==baseline['maximum_shortening_certificate']['interval'].hi-baseline['maximum_shortening_certificate']['interval'].lo
        for k,v in zip(OBSERVABLES,values):
            lo,hi=map(Q,o['prediction_intervals_exact'][k]);E=widths[k];y=Q(observations[rid][k])
            checks[k+'_uniform_interval_exact']=lo==v.lo-E and hi==v.hi+E
            residuals.append({'record_id':rid,'observable':k,'observed_exact':str(y),'prediction_interval_exact':[str(lo),str(hi)],'minimum_residual_exact':str(max(lo-y,y-hi,Q(0))),'maximum_residual_exact':str(max(abs(lo-y),abs(hi-y)))})
        audits.append({'record_id':rid,'status':'PASS_ARITHMETIC_CONDITIONAL_ON_ROOT_REVIEW' if all(checks.values()) else 'FAIL','checks':checks,'root_certificate_sha256':sha(rfile),'observer_certificate_sha256':sha(outfile),'branch_domain_id':root['branch_domain_id']})
        print(json.dumps({'record':rid,'status':audits[-1]['status'],'checks':len(checks)}),flush=True)
    lower=max(Q(v['minimum_residual_exact']) for v in residuals);aggregate=result['aggregate']
    final={'all17_preserved':len(audits)==17,'actual60_intervals_8missing':len(residuals)==60 and len(missing)==8,
        'all_independent_application_checks_pass':all(all(v is True for v in a['checks'].values()) for a in audits),
        'necessary_lower_exact_reproduced':lower==Q(aggregate['maximum_residual_lower_exact']),
        'missing_rows_identical':missing==aggregate['missing_rows'],
        'no_sufficient_upper_or_IN':aggregate['maximum_residual_upper_exact'] is None and aggregate['common_delta_family']['IN_static_stationary_root_tuple_if'] is None,
        'all_files_unchanged':all(sha(ROOT/k)==h for k,h in files.items())}
    save(RUN/'arithmetic.json',{'status':'PASS_ARITHMETIC_PENDING_ROOT_REVIEW' if all(final.values()) else 'FAIL','checks':final,'record_audits':audits,'file_hashes':files,'independent_residuals':residuals,
        'lower_exact':str(lower),'missing_rows':missing,'elapsed_seconds':time.time()-start,
        'scope':'Independent application arithmetic; whole-body and baseline static observer are disclosed replays of previously reviewed components. No independent reimplementation of those components, all17 root proof, complete fit, continuum bound or novelty is claimed.'})
    print(json.dumps({'status':read(RUN/'arithmetic.json')['status'],'checks':sum(len(a['checks']) for a in audits),'covered':len(residuals),'missing':len(missing)}))
    if not all(final.values()):raise SystemExit(1)

def finalize():
    arithmetic=read(RUN/'arithmetic.json');rootsfile=ROOT/'evidence/cohort_uniform_roots_independent_review.json';roots=read(rootsfile)
    if arithmetic['status']!='PASS_ARITHMETIC_PENDING_ROOT_REVIEW':raise RuntimeError('passing application arithmetic required')
    if any(sha(ROOT/k)!=h for k,h in arithmetic['file_hashes'].items()):raise RuntimeError('application snapshot changed')
    capsules={v['record_id']:v['audit'] for v in roots['record_audits']};checks={};out=[]
    for audit in arithmetic['record_audits']:
        rid=audit['record_id'];r=capsules[rid]
        if r['certificate_sha256']!=audit['root_certificate_sha256']:raise RuntimeError('independent root audit targets another provider')
        if audit['status']=='EXPECTED_UNRESOLVED':checks[rid+'_still_unresolved']=r['status']=='UNRESOLVED';continue
        checks[rid+'_independent_root_uniform_proof_PASS']=r['status']=='PASS' and bool(r['checks']) and all(v is True for v in r['checks'].values()) and r['uniform_quantifier_verified'] is True
        out.append({**audit,'status':'PASS_CONDITIONAL_UNIFORM_OBSERVER_APPLICATION','independent_root_review_file':rel(rootsfile),'independent_root_review_sha256':sha(rootsfile)})
    checks['exactly15_audited_root_observer_pairs']=len(out)==15
    result={'status':'PASS_PARTIAL_APPLICATION' if all(checks.values()) else 'FAIL','checks':checks,'record_audits':out,
        'arithmetic_review_sha256':sha(RUN/'arithmetic.json'),'independent_root_review_sha256':sha(rootsfile),
        'covered_observables':60,'required_observables':68,'missing_rows':arithmetic['missing_rows'],'necessary_common_discrepancy_lower_exact':arithmetic['lower_exact'],
        'sufficient_common_discrepancy_upper':None,'full_compatible_tuple_certified':False,'review_source_sha256':sha(__file__),
        'scope':'Verified60 conditional static outputs on the selected15 root families; F/K8channels unresolved. Partial application verification is not completion of the full five-gate research objective.'}
    save(ROOT/'evidence/cohort_static_observers_independent_review.json',result)
    print(json.dumps({'status':result['status'],'checks':len(checks),'covered':60,'missing':8}))
    if not all(checks.values()):raise SystemExit(1)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['arithmetic','finalize']);args=ap.parse_args()
    arithmetic() if args.stage=='arithmetic' else finalize()
