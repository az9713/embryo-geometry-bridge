"""Independent exact endpoint reduction of SAFE144+24 packets; no integration."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,base64,hashlib,io,json,datetime
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
BF=ROOT/'runs/exact_integration_bridge_v1/bulk_assembly_B_safe_v2/results.json'
SF=ROOT/'runs/exact_integration_bridge_v1/surface_assembly_B_safe_v2/results.json'
SUR=ROOT/'runs/exact_integration_bridge_v1/surface_independent_replay_safe_v2'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def write(p,d):
    if p.exists():raise ValueError('preserve '+str(p))
    p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2)+'\n')
def iv(p):
    lo,hi=Q(p['lower_exact']),Q(p['upper_exact']);assert lo<=hi;return lo,hi
def add(a,b):return a[0]+b[0],a[1]+b[1]
def scalar(a,b):return(a*b[0],a*b[1])if a>=0 else(a*b[1],a*b[0])
def hiabs(a):return max(map(abs,a))
def loabs(a):return Q(0)if a[0]<=0<=a[1]else min(map(abs,a))
def contains(p,b):a=iv(p);return a[0]<=b[0]and b[1]<=a[1]
def overlap(a,b):return a[0]<=b[1]and b[0]<=a[1]
def selfchecks():
    tests={'negative_scalar_endpoint_reversal':scalar(Q(-2),(Q(-1),Q(3)))==(Q(-6),Q(2)),'exact_signed_sum':add((Q(-1),Q(2)),(Q(3),Q(4)))==(Q(2),Q(6)),'zero_crossing_absolute_lower':loabs((Q(-1),Q(3)))==0,'negative_interval_absolute_lower':loabs((Q(-4),Q(-2)))==2,'negative_interval_absolute_upper':hiabs((Q(-4),Q(-2)))==4,'closed_endpoint_overlap':overlap((Q(0),Q(1)),(Q(1),Q(2))),'center_escape_does_not_imply_no_fixedpoint':Q(6,5)>1 and Q(6,5)-Q(4,5)/2==Q(4,5)and -1<=Q(4,5)<=1}
    assert all(v is True for v in tests.values()),tests
    return tests
def main(bulksha,pointfile,pointsha):
    if not bulksha or sha(BF)!=bulksha:raise ValueError('authoritative bulk terminal SHA required')
    if sha(SF)!='024c3070f2f3e27866dbcf771cff6c08cb5648396b16bc709af552e9797028d9':raise ValueError('authoritative SAFE surface result differs')
    checks=selfchecks();bindings={}
    def bound(p,expected=None):
        p=Path(p);h=sha(p)
        if expected is not None:assert h==expected,p
        key=p.resolve().relative_to(ROOT.resolve()).as_posix()
        if key in bindings:assert bindings[key]==h
        bindings[key]=h;return read(p)
    bulk=bound(BF,bulksha);surface=bound(SF);bp=bound(BF.parent/'protocol_v2.json');sp=bound(SF.parent/'protocol.json')
    checks['both_SAFE_producers_true_terminal_labels']=bulk['status']=='SAFE_V2_ALL144_POINT_BULK_ASSEMBLED_UNREVIEWED'and surface['status']=='SAFE_V2_ALL24_POINT_SURFACE_ASSEMBLED_UNREVIEWED'and bulk['all_frozen_bindings_unchanged']is True and surface['all_frozen_bindings_unchanged']is True
    for p in [bp,sp]:
        for name,h in p['source_and_input_hashes'].items():
            assert sha(ROOT/name)==h,name
            if name in bindings:assert bindings[name]==h
            bindings[name]=h
    bpre=bound(BF.parent/'prelaunch_independent_review.json');spre=bound(SF.parent/'prelaunch_independent_review.json')
    checks['both_SAFE_prelaunch_reviews_bound_literal_PASS']=all(v['status']=='PASS'and bool(v['checks'])and all(x is True for x in v['checks'].values())for v in [bpre,spre])and bpre['protocol_sha256']==sha(BF.parent/'protocol_v2.json')and spre['protocol_sha256']==sha(SF.parent/'protocol.json')
    recfile=ROOT/'runs/cohort_uniform_roots_v2/records/00_optoGEF2_small_col_B/result.json';rec=bound(recfile);cache=bound(ROOT/rec['point_cache_file'],rec['artifact_hashes'][rec['point_cache_file']]);blob=base64.b64decode(cache['frozen_binary_cache']['data'],validate=True);assert hashlib.sha256(blob).hexdigest()==cache['frozen_binary_cache']['sha256']
    with np.load(io.BytesIO(blob),allow_pickle=False)as z:a={k:z[k].copy()for k in z.files}
    assert a['preconditioner'].shape==(611,611)and all(np.isfinite(x).all()for x in a.values())
    oldbundle=bound(ROOT/'evidence/cohort_uniform_roots_independent_review.json','67c9d5195c36fc11814d08517070fe53ea0f10f67d45911df517bd900e098297');entry=next(x for x in oldbundle['record_audits']if x['record_id']==rec['record_id']);old=bound(ROOT/entry['review_file'],entry['review_sha256']);oldproto=bound(ROOT/'runs/cohort_uniform_root_review_v1/protocol.json',old['review_protocol_sha256']);baseproto=bound(ROOT/'runs/cohort_uniform_roots_v1/protocol.json',oldproto['producer_protocol_sha256'])
    checks['oldQ_eta_capsule_same_original_B_u0_C']=old['status']=='PASS'and all(v is True for v in old['checks'].values())and old==entry['audit']and old['certificate_sha256']==sha(recfile)and old['parameter_cell_sha256']==rec['parameter_cell_sha256']
    theta={k:Q(cache['parameters'][k])for k in ['apical_area_modulus','basal_area_modulus','apical_surface_tension','basal_surface_tension']};theta['common_gain']=Q(rec['base_strength_exact'])/Q(rec['dose_exact']);base={k:Q(v)for k,v in baseproto['base_parameters_exact'].items()}
    checks['same_exact_base_theta_mu_kappa_gain_dose']=theta==base and Q(cache['parameters']['mu'])==1 and Q(cache['parameters']['bulk'])==5000
    checks['base_theta_in_same_closed_C']=all(Q(rec['parameter_cell'][k][0])<=v<=Q(rec['parameter_cell'][k][1])for k,v in theta.items())and {k:list(map(Q,v))for k,v in rec['parameter_cell'].items()}=={k:list(map(Q,v))for k,v in baseproto['shared_parameter_cell'].items()}
    state=bound(ROOT/cache['input_state'],cache['input_state_sha256']);checks['same_exact_u0_reference_free_sigma_R']=np.array_equal(a['u'],np.asarray(state['free_displacement']))and np.array_equal(a['reference'],np.asarray(state['reference_nodes']))and np.array_equal(a['triangles'],np.asarray(state['triangles']))and a['energy_scale'].shape==(1,)and a['energy_scale'][0]>0 and sp['point_cache_file']==rec['point_cache_file']and bp['record_id']==rec['record_id']
    eta_ind=Q(old['bounds']['eta']['upper_exact']);eta_author=Q(rec['uniform_cell_result']['eta_uniform']['upper_exact']);eta=max(eta_ind,eta_author);radius=Q(rec['uniform_cell_result']['selected_radius_exact'])
    checks['old_eta_radius_same_selected_Q_proof']=Q(old['bounds']['radius']['lower_exact'])<=radius<=Q(old['bounds']['radius']['upper_exact'])and eta>=eta_ind and eta>=eta_author
    totalF=[(Q(0),Q(0))for _ in range(611)];totalH={};localpackets=[];localdict={}
    for kind,output,items,protocol in [('bulk',bulk,bulk['elements'],bp),('surface',surface,surface['edges'],sp)]:
        partF=[(Q(0),Q(0))for _ in range(611)];partH={};seen=set();sels={x['element']:x for x in bp['selectors']}if kind=='bulk'else{(x['face'],x['edge']):x for x in sp['cases']}
        for item in items:
            p=ROOT/item['file'];x=bound(p,item['sha256']);ids=x['global_indices'];n=len(ids);assert len(ids)==len(set(ids))and all(type(i)is int and 0<=i<611 for i in ids)and len(x['force_integral_minus_cache'])==n
            tag=x['element']if kind=='bulk'else(x['face'],x['edge']);assert tag not in seen;seen.add(tag);sel=sels[tag];assert ids==sel['global_indices'];hkey='H_integral_minus_cache_upper_triangle'if kind=='bulk'else'H_integral_minus_cache_upper';slots={(z['i'],z['j']):z['enclosure']for z in x[hkey]};assert len(slots)==len(x[hkey])and set(slots)=={(i,j)for i in range(n)for j in range(i,n)}
            if kind=='bulk':
                good=x['truncation_order']==sel['degree']and Q(x['rho_exact'])==Q(sel['point_rho_exact'])and iv(x['force_max_tail_upper'])[1]<=Q(bp['force_tail_target_exact'])and iv(x['H_max_tail_upper'])[1]<=Q(bp['H_entry_tail_target_exact'])
            else:good=x['order']==28 and x['precision_bits']==512 and Q(x['q0_exact'])==Q(sel['domain']['q0_exact'])and Q(x['rho_exact'])==Q(sel['domain']['rho_exact'])and Q(sel['domain']['G_lower_exact'])>0 and Q(sel['domain']['hoop_lower_exact'])>0
            checks[f'{kind}_{tag}_local_hash_map_domain_tail_policy']=bool(good)
            for i,v in zip(ids,x['force_integral_minus_cache']):partF[i]=add(partF[i],iv(v));totalF[i]=add(totalF[i],iv(v))
            for(i,j),v in slots.items():
                key=tuple(sorted((ids[i],ids[j])));partH[key]=add(partH.get(key,(Q(0),Q(0))),iv(v));totalH[key]=add(totalH.get(key,(Q(0),Q(0))),iv(v))
            localpackets.append({'kind':kind,'file':item['file'],'sha256':item['sha256']});localdict[kind,tag]=x
        checks[kind+'_exact_unique_complete_coverage']=seen==set(range(144))if kind=='bulk'else seen=={(k,e)for k in [0,1]for e in range(12)}
        checks[kind+'_611_force_matches_exact_endpoint_sum']=len(output['force_integral_minus_cache'])==611 and all(contains(v,s)for v,s in zip(output['force_integral_minus_cache'],partF))
        ohkey='bulk_H_upper'if kind=='bulk'else'H_upper';oh={(x['i'],x['j']):x['enclosure']for x in output[ohkey]};checks[kind+'_all_symmetric_H_slots_match_exact_sum']=len(oh)==len(output[ohkey])and set(oh)==set(partH)and all(contains(oh[k],v)for k,v in partH.items())
    checks['all168_packets_and6676_slots']=len(localpackets)==168 and len(totalH)==6676
    RF=[];Rnorm=Q(0)
    for row in a['preconditioner']:
        exact=[Q.from_float(float(v))for v in row];value=(Q(0),Q(0))
        for r,f in zip(exact,totalF):value=add(value,scalar(r,f))
        RF.append(value);Rnorm=max(Rnorm,sum(map(abs,exact),Q(0)))
    lower=max(map(loabs,RF));upper=max(map(hiabs,RF));rows=[Q(0)]*611
    for(i,j),v in totalH.items():
        rows[i]+=hiabs(v)
        if i!=j:rows[j]+=hiabs(v)
    hnorm=max(rows);maplower=max(Q(0),lower-eta);mapoutside=maplower>radius
    # All six newly integrated bulk and surface cases must bind to these exact
    # SAFE producer packets. Old v1 point outputs are never inspected here.
    replay=bound(HERE/'results.json');rproto=bound(HERE/'protocol.json',replay['protocol_sha256']);checks['safe6bulk_review_PASS']=replay['status']=='PASS_SAFE_SIX_POINT_BULK_REPLAYS'and all(v is True for v in replay['checks'].values())and replay['producer_results_sha256']==bulksha and replay['source_sha256']==sha(HERE/'replay.py')
    for k,v in rproto['input_hashes'].items():assert sha(ROOT/k)==v
    for item in replay['elements']:
        x=bound(HERE/item['file'],item['sha256']);prod=localdict['bulk',x['element']];ph={(z['i'],z['j']):iv(z['enclosure'])for z in prod['H_integral_minus_cache_upper_triangle']};ih={(z['i'],z['j']):iv(z['interval'])for z in x['H_upper']};checks[f'safe_bulk_replay{x["element"]}_all_entries_match']=x['global_indices']==prod['global_indices']and len(x['force'])==len(prod['force_integral_minus_cache'])and all(overlap(iv(v),iv(w))for v,w in zip(x['force'],prod['force_integral_minus_cache']))and set(ph)==set(ih)and all(overlap(v,ph[k])for k,v in ih.items())
    sreceipt=bound(ROOT/'evidence/exact_surface_safe_independent_review.json','248ba5000b27f8e29ede2b501b00bf16e5fbdc6624f06a65a8ced339c39c37a2');sreview=bound(SUR/'results.json','fd6fd10fe63827d8cdbdcc8b12d4df2a5537e324d50fa31b05a74673878b3859');sproto=bound(SUR/'protocol.json',sreview['protocol_sha256']);checks['safe6surface_review_PASS']=sreceipt['status'].startswith('PASS')and all(v is True for v in sreceipt['checks'].values())and sreview['all_bindings_current']is True and len(sreview['outputs'])==6 and sproto['precision_bits']==640 and sproto['truncation_order']==28
    for k,v in sproto['input_hashes'].items():assert sha(ROOT/k)==v
    for item in sreview['outputs']:
        x=bound(SUR/item['file'],item['sha256']);tag=(x['face'],x['edge']);prod=localdict['surface',tag];ph={(z['i'],z['j']):iv(z['enclosure'])for z in prod['H_integral_minus_cache_upper']};ih={(z['i'],z['j']):iv(z['enclosure'])for z in x['Hessian_enclosures']};checks[f'safe_surface_replay{tag}_all_entries_match']=x['global_indices']==prod['global_indices']and len(x['force_enclosures'])==len(prod['force_integral_minus_cache'])and all(overlap(iv(v),iv(w))for v,w in zip(x['force_enclosures'],prod['force_integral_minus_cache']))and set(ph)==set(ih)and all(overlap(v,ph[k])for k,v in ih.items())
    if pointfile:
        assert 'safe_v2' in pointfile
        point=bound(ROOT/pointfile,pointsha)
        for k,v in point['source_and_input_hashes'].items():assert sha(ROOT/k)==v
        checks['SAFE_combined_point_exact_base_eta_radius_scope']=point['status']=='SAFE_V2_COMPLETE_POINT_DISCREPANCY_UNREVIEWED'and {k:Q(v)for k,v in point['base_parameters_exact'].items()}==theta and Q(point['old_uniform_eta_upper_exact'])==eta and Q(point['old_radius_exact'])==radius and point['uniform_a_or_b_established']is False and point['new_integral_root_certified']is False
        checks['SAFE_combined_point_all611F_and_RF_enclose_exact_sum']=len(point['force_integral_minus_cache'])==len(point['preconditioned_force_defect'])==611 and all(contains(v,s)for v,s in zip(point['force_integral_minus_cache'],totalF))and all(contains(v,s)for v,s in zip(point['preconditioned_force_defect'],RF));checks['SAFE_combined_point_norm_bounds_valid']=Q(point['preconditioned_force_defect_norm_lower_exact'])<=lower<=upper<=iv(point['preconditioned_force_defect_norm_upper'])[1]and iv(point['point_H_defect_row_norm_upper'])[1]>=hnorm and iv(point['point_conservative_preconditioned_H_defect_norm_upper'])[1]>=Rnorm*hnorm
        checks['SAFE_combined_point_map_formula_and_boolean_valid']=Q(point['new_map_center_distance_lower_exact'])==max(Q(0),Q(point['preconditioned_force_defect_norm_lower_exact'])-eta)and point['old_map_center_test']['center_leaves_old_cube_certified']==(Q(point['new_map_center_distance_lower_exact'])>radius)
    checks['all_frozen_bindings_current']=all(sha(ROOT/k)==v for k,v in bindings.items())
    exact={'F_endpoint_sums':[[str(x),str(y)]for x,y in totalF],'H_upper_endpoint_sums':[{'i':i,'j':j,'interval':list(map(str,v))}for(i,j),v in sorted(totalH.items())],'RF_endpoint_sums':[[str(x),str(y)]for x,y in RF],'RF_norm_lower_exact':str(lower),'RF_norm_upper_exact':str(upper),'H_row_norm_upper_exact':str(hnorm),'Rnorm_exact':str(Rnorm),'conservative_RH_norm_upper_exact':str(Rnorm*hnorm),'eta_author_exact':str(eta_author),'eta_independent_exact':str(eta_ind),'eta_used_exact':str(eta),'old_radius_exact':str(radius),'new_integral_map_center_distance_lower_exact':str(maplower),'new_integral_map_center_outside_old_cube':mapoutside}
    write(HERE/'exact_168_reduction.json',exact)
    result={'status':'PASS_SAFE_V2_COMPLETE_POINT_REDUCTION'if all(checks.values())else'FAIL_SAFE_V2_COMPLETE_POINT_REDUCTION','checks':checks,'review_source_sha256':sha(__file__),'source_and_input_hashes':bindings,'all168_packets':localpackets,'exact_reduction_sha256':sha(HERE/'exact_168_reduction.json'),'certificate':{k:v for k,v in exact.items()if k.endswith('_exact')or k.endswith('_cube')},'old_v1_claims_restored':False,'new_conditional_negative_map_certificate_proved':mapoutside,'proof':'At exact base theta0 in C the old predictor is u0. Unaffected oldQ certificate bounds ||R FQ(u0)|| by max(author,independent eta). New SAFE complete discrepancy intervals give ||R(FI-FQ)(u0)||>=L. Therefore ||R FI(u0)||>=max(0,L-eta). Only if this exceeds r is the NEW map center outside the OLD cube. No root nonexistence or root location follows.','reintegration_coverage':{'bulk':6,'bulk_total':144,'surface':6,'surface_total':24,'remaining_packets':'Reviewed generic safe engines, complete local domain/tail/hash/map checks and exact signed reduction; not independently reintegrated individually.'},'scope':'One original B-center and exact base theta only; no uniform a/b, new integral root/stability/PDE/biological validation or full-five-gate completion.'}
    write(HERE/'reduction_results.json',result);write(ROOT/'evidence/exact_complete_point_B_independent_review_safe_v2.json',result);print(json.dumps({'status':result['status'],'checks':len(checks),'failed':[k for k,v in checks.items()if not v],'RF_lower':float(lower),'RF_upper':float(upper),'eta':float(eta),'maplower':float(maplower),'r':float(radius),'map_outside':mapoutside,'receipt_sha256':sha(ROOT/'evidence/exact_complete_point_B_independent_review_safe_v2.json')}),flush=True)
    assert all(v is True for v in checks.values()),checks
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--bulk-results-sha',required=True);p.add_argument('--point-file');p.add_argument('--point-sha');x=p.parse_args();main(x.bulk_results_sha,x.point_file,x.point_sha)
