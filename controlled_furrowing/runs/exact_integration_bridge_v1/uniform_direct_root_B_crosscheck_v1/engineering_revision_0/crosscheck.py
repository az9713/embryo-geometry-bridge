"""Exact scalar crosscheck; no integration, matrix multiplication, or solve."""
from pathlib import Path
from fractions import Fraction as Q
import json, hashlib, sys, datetime, time
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
BASE=ROOT/'runs/exact_integration_bridge_v1'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p): return json.loads(p.read_text(encoding='utf-8'))
def save(p,x): p.write_text(json.dumps(x,indent=2)+'\n',encoding='utf-8')
def utc(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sup_interval(x):
    if not isinstance(x,(tuple,list)) or len(x)!=2: raise ValueError('interval pair')
    lo,hi=map(Q,x)
    if lo>hi: raise ValueError('ordered interval')
    return max(abs(lo),abs(hi))
def scalar(eta,m,zq,b,d,r):
    eta,m,zq,b,d,r=map(Q,(eta,m,zq,b,d,r))
    if min(eta,m,zq,b,d)<0 or r<=0: raise ValueError('nonnegative bounds and positive radius')
    z=zq+b; e=eta+m+z*d
    ok=z<1 and e+z*r<r
    return dict(z=z,eta=e,lhs=e+z*r,passed=ok,error=e/(1-z) if ok else None)
def selftest():
    checks={}
    checks['zero_crossing_interval_retained']=sup_interval(['-1','2'])==2
    checks['negative_interval_absolute']=sup_interval(['-3','-2'])==3
    checks['strict_selfmap_pass']=scalar(0,Q(1,20),Q(1,4),0,0,1)['passed']
    checks['strict_boundary_fails']=not scalar(0,Q(3,4),Q(1,4),0,0,1)['passed']
    checks['noncontraction_fails']=not scalar(0,0,1,0,0,1)['passed']
    checks['actual_binary_difference_not_rounded_subtraction']=Q.from_float(.3)-Q.from_float(.2)!=Q.from_float(.3-.2)
    for name,fun in [('reverse_interval',lambda:sup_interval([2,1])),('negative_bound',lambda:scalar(-1,0,0,0,0,1)),('zero_radius',lambda:scalar(0,0,0,0,0,0))]:
        try: fun(); checks[name]=False
        except ValueError: checks[name]=True
    assert all(checks.values())
    return checks
def inputs():
    paths=[BASE/'uniform_full_B_physical_axis_independent_reduction_v2'/n for n in ['protocol.json','results.json','exact_reduction.json','delivery.json']]
    paths += [BASE/'cached_Q_common_B_v1'/n for n in ['protocol.json','results.json','delivery.json']]
    paths += [BASE/'retube_B_preparation'/n for n in ['protocol.json','rounding_ledger.json','candidate_vectors.npz']]
    paths += [ROOT/'evidence/cohort_uniform_roots_independent_review.json',ROOT/'runs/cohort_uniform_root_review_v1/records/optoGEF2_small_col_B.json',ROOT/'runs/cohort_uniform_roots_v2/records/00_optoGEF2_small_col_B/result.json',ROOT/'evidence/exact_retube_B_geometry_independent_review.json',ROOT/'evidence/exact_cached_Q_common_B_reduction_review.json']
    reduction=read(paths[1]); paths += [ROOT/p for p in reduction['input_hashes']]
    return {str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in paths}
def prepare():
    assert not (HERE/'protocol.json').exists()
    s=selftest(); save(HERE/'selfchecks.json',{'status':'PASS','checks':s})
    save(HERE/'protocol.json',{'version':1,'utc':utc(),'source_sha256':sha(Path(__file__)),'input_hashes':inputs(),'selfchecks_sha256':sha(HERE/'selfchecks.json'),'arithmetic':'Exact Fraction endpoint bounds, binary floats converted by Fraction.from_float; no floating acceptance','premises':'Reviewed uniform exact-integral-minus-cached-Q force/H boxes on SAME C and common nodal cube; old B cached-Q uniform center residual; fresh cached-Q contraction','formula':'eta_direct=eta_Q+max_i sup|Dbox_i-d_i|+(z_Q+b_absR_H)*||d||inf; strict eta_direct+z*r<r','scope':'One B exact-integration finite-P2 stationary-root family if strict test succeeds; no new integration/solver or independent physical kernel implementation.'})
    print('PREPARED',sha(HERE/'protocol.json'),flush=True)
def run():
    t=time.perf_counter(); p=read(HERE/'protocol.json'); checks={}
    checks['frozen_source_matches']=p['source_sha256']==sha(Path(__file__))
    checks['all_frozen_inputs_match']=all(sha(ROOT/k)==v for k,v in p['input_hashes'].items())
    checks['selfchecks_match']=sha(HERE/'selfchecks.json')==p['selfchecks_sha256']
    assert all(checks.values())
    red=read(BASE/'uniform_full_B_physical_axis_independent_reduction_v2/exact_reduction.json')
    rr=read(BASE/'uniform_full_B_physical_axis_independent_reduction_v2/results.json')
    cq=read(BASE/'cached_Q_common_B_v1/results.json'); cp=read(BASE/'cached_Q_common_B_v1/protocol.json')
    old=read(ROOT/'runs/cohort_uniform_root_review_v1/records/optoGEF2_small_col_B.json')
    if 'audit'in old: old=old['audit']
    # Standalone record has the same capsule structure as the combined receipt.
    if 'bounds'not in old and 'independent_result'in old: old=old['independent_result']
    ledger=read(BASE/'retube_B_preparation/rounding_ledger.json')
    with np.load(BASE/'retube_B_preparation/candidate_vectors.npz',allow_pickle=False) as a:
        d=[Q.from_float(float(x))-Q.from_float(float(y)) for x,y in zip(a['u0'],a['u1'])]
    D=red['signed_RF_endpoint_sums']; assert len(d)==len(D)==611
    mrows=[sup_interval([Q(v[0])-x,Q(v[1])-x]) for v,x in zip(D,d)]
    m=max(mrows); dn=max(map(abs,d)); r=Q(cp['old_radius_exact'])
    eta=Q(old['bounds']['eta']['upper_exact']); zq=Q(cq['ZQ_upper_exact']); b=Q(red['b_absR_H_upper_exact'])
    checks['old_B_review_individually_PASS']=old['status']=='PASS' and all(v is True for v in old['checks'].values())
    checks['old_eta_matches_prepared_ledger']=eta==Q(read(BASE/'retube_B_preparation/protocol.json')['eta_old_uniform_exact'])
    checks['actual_binary_shift_matches_all611_retained_values']=d==list(map(Q,ledger['applied_decrement_exact']))
    checks['same_parameter_cell']=rr['parameter_cell_sha256']==cq['parameter_cell_sha256']==old['parameter_cell_sha256']
    checks['same_common_radius']=rr['common_radius_exact']==cq['common_radius_exact'] and r+dn==Q(cq['common_radius_exact'])
    checks['fresh_Q_same_Z']=zq==Q(red['ZQ_upper_exact'])
    checks['full168_component_checks']=all(v is True for v in rr['checks'].values())
    checks['fresh_Q_checks']=all(v is True for v in cq['checks'].values())
    answer=scalar(eta,m,zq,b,dn,r)
    checks['all_frozen_inputs_still_match']=all(sha(ROOT/k)==v for k,v in p['input_hashes'].items())
    assert all(checks.values())
    out={'status':'CERTIFIED_DIRECT_UNIFORM_RESIDUAL_ROUTE' if answer['passed'] else 'UNRESOLVED_DIRECT_UNIFORM_RESIDUAL_SELFМAP','utc':utc(),'protocol_sha256':sha(HERE/'protocol.json'),'source_sha256':sha(Path(__file__)),'checks':checks,'arithmetic':'Exact rational scalar inequalities; rounded displays never used for acceptance','parameter_cell_sha256':rr['parameter_cell_sha256'],'record_id':cq['record_id'],'common_geometry_domain_id':rr['common_geometry_domain_id'],'old_eta_exact':str(eta),'m_direct_exact':str(m),'m_maximizing_free_index':mrows.index(m),'actual_shift_norm_exact':str(dn),'old_radius_exact':str(r),'common_radius_exact':cq['common_radius_exact'],'ZQ_exact':str(zq),'b_exact_row_majorant':str(b),'Z_integral_exact':str(answer['z']),'eta_direct_exact':str(answer['eta']),'selfmap_lhs_exact':str(answer['lhs']),'strict_selfmap_pass':answer['passed'],'root_error_exact':None if answer['error']is None else str(answer['error']),'all611_m_rows_exact':list(map(str,mrows)),'display':{'eta_Q':float(eta),'m_direct':float(m),'z':float(answer['z']),'eta_direct':float(answer['eta']),'selfmap_lhs':float(answer['lhs']),'radius':float(r),'root_error':None if answer['error']is None else float(answer['error'])},'reused_physical_kernels':True,'no_new_matrix_integration_or_solver':True,'scope':'Separate scalar sufficient test on same B/C/common domain. Failure is inconclusive; success certifies one finite-P2 stationary-root family only. Generic physical kernels and their proof reviews are shared, not wholly independent physics. No PDE, other records, global domains, temporal dynamics, source-aligned biology or novelty claim.','seconds':time.perf_counter()-t}
    save(HERE/'results.json',out); print(out['status'],len(checks),'checks',out['display'],flush=True)
if __name__=='__main__':
    {'prepare':prepare,'run':run}[sys.argv[1]]()
