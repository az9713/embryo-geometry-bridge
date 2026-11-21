"""Exact168-packet reduction and rigorous negative SELF-MAP-center test."""
from pathlib import Path
from fractions import Fraction as F
import argparse,base64,hashlib,io,json
import numpy as np
HERE=Path(__file__).resolve().parent;CF=HERE.parents[2]
BF=CF/'runs/exact_integration_bridge_v1/bulk_assembly_B/results.json'
SF=CF/'runs/exact_integration_bridge_v1/surface_assembly_B/results.json'
PF=CF/'runs/exact_integration_bridge_v1/complete_point_B/results.json'
EXPECTED={'bulk':'c7938f431321d5b189efcd057fbf6635fbf0ea0d583cf0431f633ceeb7c23dfd','surface':'821bc183209a763703613060412590eb88a400abd7bcf9ca482b838815ceadcb','point':'57bb157ba4ea40d790cb01c3775e03b6dbb7a811c8961f78d32167532af3b98d'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def interval(v):
    x,y=F(v['lower_exact']),F(v['upper_exact'])
    if x>y:raise ValueError('reversed interval')
    return x,y
def add(x,y):return x[0]+y[0],x[1]+y[1]
def scaled(k,x):return(k*x[0],k*x[1])if k>=0 else(k*x[1],k*x[0])
def absmax(x):return max(map(abs,x))
def absmin(x):return F(0)if x[0]<=0<=x[1]else min(map(abs,x))
def encloses(p,x):a,b=interval(p);return a<=x[0]and x[1]<=b
def canon(p):return{key:[str(F(x))for x in pair]for key,pair in p.items()}
def main(surface_review,surface_review_sha):
    if any(sha(file)!=EXPECTED[k]for k,file in[('bulk',BF),('surface',SF),('point',PF)]):raise ValueError('authoritative terminal outputs changed')
    checks={};bindings={};localhashes=[]
    def bound(path,digest=None):
        path=Path(path);h=sha(path)
        if digest is not None and h!=digest:raise ValueError('bound source changed: '+str(path))
        key=path.resolve().relative_to(CF.resolve()).as_posix()
        if key in bindings and bindings[key]!=h:raise ValueError('conflicting hash')
        bindings[key]=h;return read(path)
    bulk=bound(BF);surface=bound(SF);point=bound(PF);bp=bound(BF.parent/'protocol_v2.json');sp=bound(SF.parent/'protocol.json')
    for packet in[bp,sp,point]:
        for k,h in packet['source_and_input_hashes'].items():
            if sha(CF/k)!=h:raise ValueError('frozen component input changed '+k)
            name=(CF/k).resolve().relative_to(CF.resolve()).as_posix()
            if name in bindings and bindings[name]!=h:raise ValueError('incompatible shared target binding')
            bindings[name]=h
    recordfile=CF/'runs/cohort_uniform_roots_v2/records/00_optoGEF2_small_col_B/result.json';record=bound(recordfile);cachefile=CF/record['point_cache_file'];cache=bound(cachefile,record['artifact_hashes'][record['point_cache_file']]);raw=base64.b64decode(cache['frozen_binary_cache']['data'],validate=True)
    if hashlib.sha256(raw).hexdigest()!=cache['frozen_binary_cache']['sha256']:raise ValueError('binary cache corrupt')
    with np.load(io.BytesIO(raw),allow_pickle=False)as z:arrays={k:z[k].copy()for k in z.files}
    # Original independently assembled uniform-residual capsule and its base.
    bundle=bound(CF/'evidence/cohort_uniform_roots_independent_review.json','67c9d5195c36fc11814d08517070fe53ea0f10f67d45911df517bd900e098297');entry=next(v for v in bundle['record_audits']if v['record_id']==record['record_id']);old=bound(CF/entry['review_file'],entry['review_sha256']);oldprotocol=bound(CF/'runs/cohort_uniform_root_review_v1/protocol.json',old['review_protocol_sha256']);baseprotocol=bound(CF/'runs/cohort_uniform_roots_v1/protocol.json',oldprotocol['producer_protocol_sha256'])
    checks['old_independent_capsule_valid_for_exact_selected_root']=old['status']=='PASS'and all(v is True for v in old['checks'].values())and old==entry['audit']and old['certificate_sha256']==sha(recordfile)and old['parameter_cell_sha256']==record['parameter_cell_sha256']
    theta={k:F(cache['parameters'][k])for k in['apical_area_modulus','basal_area_modulus','apical_surface_tension','basal_surface_tension']};theta['common_gain']=F(record['base_strength_exact'])/F(record['dose_exact']);base={k:F(v)for k,v in baseprotocol['base_parameters_exact'].items()}
    checks['same_exact_old_predictor_base_theta']=theta==base=={k:F(v)for k,v in point['base_parameters_exact'].items()}
    checks['theta0_inside_same_C_including_closed_boundary']=all(F(record['parameter_cell'][k][0])<=v<=F(record['parameter_cell'][k][1])for k,v in theta.items())and canon(record['parameter_cell'])==canon(baseprotocol['shared_parameter_cell'])
    state=bound(CF/cache['input_state'],cache['input_state_sha256']);checks['same_exact_saved_u0_reference_and_free_coordinates']=np.array_equal(arrays['u'],np.asarray(state['free_displacement']))and np.array_equal(arrays['reference'],np.asarray(state['reference_nodes']))and np.array_equal(arrays['triangles'],np.asarray(state['triangles']))and len(arrays['free'])==611
    checks['same_positive_binary_normalization_and_preconditioner']=arrays['energy_scale'].shape==(1,)and arrays['energy_scale'][0]>0 and arrays['preconditioner'].shape==(611,611)and sp['point_cache_file']==record['point_cache_file']and bp['record_id']==record['record_id']
    checks['exact_gain_times_dose_not_rounded_solver_strength']=theta['common_gain']*F(record['dose_exact'])==F(record['base_strength_exact'])and cache['uniform_target_warning'].startswith('Use exact shared')
    # Uniform eta at u0 follows because p_theta=u0+S(theta-base), now delta=0.
    eta_ind=F(old['bounds']['eta']['upper_exact']);eta_producer=F(record['uniform_cell_result']['eta_uniform']['upper_exact']);eta=max(eta_ind,eta_producer);radius=F(record['uniform_cell_result']['selected_radius_exact'])
    checks['same_exact_radius_and_conservative_independent_eta']=radius==F(point['old_radius_exact'])and F(old['bounds']['radius']['lower_exact'])<=radius<=F(old['bounds']['radius']['upper_exact'])and F(point['old_uniform_eta_upper_exact'])==eta_producer and eta>=eta_ind and eta>=eta_producer
    # Domain proof is independent and covers all24 coherent closed edges.
    dom=bound(CF/'runs/exact_integration_bridge_v1/surface_assembly_B/prelaunch_independent_review.json','852ea1c24bcfb54e5b2733acd518c3e24d2be5c0274d3c06b2e36287fcab33e5');checks['all24_independent_wholeedge_domains_PASS']=dom['status']=='PASS'and all(v is True for v in dom['checks'].values())and dom['protocol_sha256']==sha(SF.parent/'protocol.json')
    Fsum=[(F(0),F(0))for _ in range(611)];Hsum={};surfaceF=[(F(0),F(0))for _ in range(611)];surfaceH={};seenbulk=set();seensurface=set();cases={(x['face'],x['edge']):x for x in sp['cases']};selectors={x['element']:x for x in bp['selectors']}
    for kind,items in[('bulk',bulk['elements']),('surface',surface['edges'])]:
        for item in items:
            local=bound(CF/item['file'],item['sha256']);ids=local['global_indices'];n=len(ids)
            if len(set(ids))!=n or any(not 0<=k<611 for k in ids)or len(local['force_integral_minus_cache'])!=n:raise ValueError('invalid local free map')
            if kind=='bulk':
                tag=local['element'];seenbulk.add(tag);selected=selectors[tag];good=ids==selected['global_indices']and local['truncation_order']==selected['degree']and F(local['rho_exact'])==F(selected['point_rho_exact']);hkey='H_integral_minus_cache_upper_triangle'
            else:
                tag=(local['face'],local['edge']);seensurface.add(tag);selected=cases[tag];good=ids==selected['global_indices']and local['order']==sp['truncation_order']and local['precision_bits']==sp['precision_bits']and F(local['q0_exact'])==F(selected['domain']['q0_exact'])and F(local['rho_exact'])==F(selected['domain']['rho_exact']);hkey='H_integral_minus_cache_upper'
            ht={(v['i'],v['j']):v['enclosure']for v in local[hkey]}
            if len(ht)!=len(local[hkey])or set(ht)!={(i,j)for i in range(n)for j in range(i,n)}:raise ValueError('local H coverage incomplete')
            checks[f'{kind}_{tag}_hash_policy_domain_free_map']=good
            for k,v in zip(ids,local['force_integral_minus_cache']):
                Fsum[k]=add(Fsum[k],interval(v))
                if kind=='surface':surfaceF[k]=add(surfaceF[k],interval(v))
            for(i,j),v in ht.items():
                key=tuple(sorted((ids[i],ids[j])));Hsum[key]=add(Hsum.get(key,(F(0),F(0))),interval(v))
                if kind=='surface':surfaceH[key]=add(surfaceH.get(key,(F(0),F(0))),interval(v))
            localhashes.append({'kind':kind,'file':item['file'],'sha256':item['sha256']})
    checks['all168_contributions_exactly_144_plus24']=len(localhashes)==168 and len(bulk['elements'])==len(seenbulk)==144 and seenbulk==set(range(144))and len(surface['edges'])==len(seensurface)==24 and seensurface==set(cases)
    checks['surface611_scatter_contains_exact_local_sums']=len(surface['force_integral_minus_cache'])==611 and all(encloses(p,x)for p,x in zip(surface['force_integral_minus_cache'],surfaceF))
    sh={(v['i'],v['j']):v['enclosure']for v in surface['H_upper']};checks['surface_sparse_H_unique_complete_contains_sums']=len(sh)==len(surface['H_upper'])and set(sh)==set(surfaceH)and all(encloses(sh[k],v)for k,v in surfaceH.items())
    checks['complete611_force_contains_all168_exact_sums']=len(point['force_integral_minus_cache'])==611 and all(encloses(p,x)for p,x in zip(point['force_integral_minus_cache'],Fsum))
    checks['combined_H6676_slots_match']=len(Hsum)==point['H_upper_slots']==6676
    R=arrays['preconditioner'];RF=[];Rnorm=F(0)
    for row in R:
        rr=[F.from_float(float(x))for x in row];v=(F(0),F(0))
        for a,b in zip(rr,Fsum):v=add(v,scaled(a,b))
        RF.append(v);Rnorm=max(Rnorm,sum(map(abs,rr),F(0)))
    checks['all611_complete_signed_RF_enclose_exact_reduction']=len(point['preconditioned_force_defect'])==611 and all(encloses(p,v)for p,v in zip(point['preconditioned_force_defect'],RF))
    lo=max(map(absmin,RF));hi=max(map(absmax,RF));claimed_lo=F(point['preconditioned_force_defect_norm_lower_exact']);claimed_hi=interval(point['preconditioned_force_defect_norm_upper'])[1]
    checks['claimed_full_defect_norm_bounds_valid']=claimed_lo<=lo<=hi<=claimed_hi
    rows=[F(0)]*611
    for(i,j),v in Hsum.items():
        rows[i]+=absmax(v)
        if i!=j:rows[j]+=absmax(v)
    hn=max(rows);checks['complete_H_and_RH_norm_upper_bounds']=interval(point['point_H_defect_row_norm_upper'])[1]>=hn and interval(point['point_conservative_preconditioned_H_defect_norm_upper'])[1]>=Rnorm*hn
    independent_map_lower=max(F(0),lo-eta);producer_map_formula=max(F(0),claimed_lo-eta_producer)
    checks['producer_map_lower_formula_exact']=F(point['new_map_center_distance_lower_exact'])==producer_map_formula
    checks['independent_new_integral_map_center_OUTSIDE_OLD_CUBE']=independent_map_lower>radius
    # Independent actualsurface component review capsule, kept separate from
    # complete rational reduction and genericengine mathematical review.
    component=bound(CF/surface_review,surface_review_sha);chome=(CF/surface_review).parent;cp=bound(chome/'protocol.json',component['protocol_sha256'])
    for k,v in cp['input_hashes'].items():
        if sha(CF/k)!=v:raise ValueError('independent surface input changed '+k)
    checks['surface_component_independent_review_PASS']=component['status'].startswith('PASS')and component['all_bindings_current']is True and len(component['outputs'])==6 and cp['precision_bits']==640 and cp['truncation_order']==28
    six=[]
    for output in component['outputs']:
        edge=bound(chome/output['file'],output['sha256']);key=(edge['face'],edge['edge']);matched=next(x for x in surface['edges']if(x['face'],x['edge'])==key);current=read(CF/matched['file']);oldH={(x['i'],x['j']):interval(x['enclosure'])for x in edge['Hessian_enclosures']};newH={(x['i'],x['j']):interval(x['enclosure'])for x in current['H_integral_minus_cache_upper']}
        def meet(a,b):return a[0]<=b[1]and b[0]<=a[1]
        ok=edge['status'].startswith('PASS')and bool(edge['checks'])and all(v is True for v in edge['checks'].values())and edge['global_indices']==current['global_indices']and len(edge['force_enclosures'])==len(current['force_integral_minus_cache'])and all(meet(interval(a),interval(b))for a,b in zip(edge['force_enclosures'],current['force_integral_minus_cache']))and set(oldH)==set(newH)and all(meet(oldH[k],v)for k,v in newH.items())
        checks[f'independent_surface{key}_matches_complete24_point_output']=ok;six.append(key)
    checks['six_distinct_independent_surface_cases_bound_to_current24']=len(set(six))==6
    checks['all_frozen_inputs_still_match']=all(sha(CF/k)==v for k,v in bindings.items())
    exact={'force_endpoint_sums':[[str(a),str(b)]for a,b in Fsum],'H_upper_endpoint_sums':[{'i':i,'j':j,'interval':[str(a),str(b)]}for(i,j),(a,b)in sorted(Hsum.items())],'RF_endpoint_sums':[[str(a),str(b)]for a,b in RF],'RF_norm_lower_exact':str(lo),'RF_norm_upper_exact':str(hi),'conservative_H_row_upper_exact':str(hn),'conservative_RH_upper_exact':str(Rnorm*hn),'producer_eta_exact':str(eta_producer),'independent_eta_exact':str(eta_ind),'eta_used_exact':str(eta),'radius_exact':str(radius),'new_integral_map_center_distance_lower_exact':str(independent_map_lower),'strict_excess_over_radius_exact':str(independent_map_lower-radius)}
    HERE.mkdir(parents=True,exist_ok=True);(HERE/'exact_reduction.json').write_text(json.dumps(exact,indent=2)+'\n')
    result={'status':'PASS_COMPLETE_POINT_DISCREPANCY_AND_NEGATIVE_SELF_MAP_CERTIFICATE'if all(checks.values())else'FAIL_REVIEW','checks':checks,'review_source_sha256':sha(Path(__file__)),'source_sha256':sha(CF/'scripts/exact_point_bridge_B.py'),'producer_point_results_sha256':sha(PF),'source_and_input_hashes':bindings,'exact_reduction_sha256':sha(HERE/'exact_reduction.json'),'local_packets':localhashes,'norm_and_map_certificate':{k:v for k,v in exact.items()if k.endswith('_exact')},'base_parameters_exact':{k:str(v)for k,v in theta.items()},'proof':'At same theta0 in C, p_theta0=u0. Independentlycertified oldeta bounds||R FQ(u0)||. Exact completebulk+surfaceerrorintervals prove||R(FI-FQ)(u0)||>=L. Reverse trianglegives||R FI(u0)||>=L-max(eta_author,eta_independent)>r. Thus NEW map u-R FI(u), with oldR and oldcube, sendsoldcenteroutside.','scope_limits':['Oldquadraturemap/proofs remainvalid; targetequationchanged.','No equilibrium nonexistence, no displacementlowerbound for an integralroot, no proofof where suchroot lies.','No parameter/tube-uniforma or b, newintegralroot/stability/PDE/sourcebiologicalvalidation.','Exact endpointreduction is rationalarithmetic onvalidatedlocalenclosures, not exactevaluationofintegrals.']}
    (HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n');(CF/'evidence/exact_complete_point_B_independent_review.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':result['status'],'checks':len(checks),'failed':[k for k,v in checks.items()if not v],'RF_lower':float(lo),'RF_upper':float(hi),'eta_used':float(eta),'map_center_lower':float(independent_map_lower),'radius':float(radius),'receipt_sha256':sha(CF/'evidence/exact_complete_point_B_independent_review.json')}))
    if not all(checks.values()):raise SystemExit(1)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--surface-review',required=True);p.add_argument('--surface-review-sha',required=True);a=p.parse_args();main(a.surface_review,a.surface_review_sha)
