"""All144 order-selector audit and tiny sparse-assembly tests; no integration."""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import hashlib,json,sys,importlib.util
import numpy as np
from flint import arb,ctx
H=Path(__file__).resolve().parent;C=next(p for p in H.parents if (p/'scripts/exact_fe_target_safe.py').exists());OUT=C/'runs/exact_integration_bridge_v1/bulk_assembly_B_safe_v2'
SOURCE=C/'scripts/exact_bulk_bridge_assembly_safe.py'
EXPECTED='84a36fe274e3606d77adfe668fcdfc975e84996319d7bfc3cf4883d31de94950'
PROTOCOL='92488086a01c04d8e0e307f606484a05a701328d192465b7ea9757a97c64d015'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def independent_rho(f,poly):
    shifted={}
    if f['chart']['kind']=='off_axis':
        lo=min(x[0]for x in f['corners']);hi=max(x[0]for x in f['corners']);rc=(lo+hi)/2;scale=(hi-lo)/2;d=max(0,-min(i for i,j in poly));z=max(abs(x[1]-f['zcenter'])for x in f['corners'])
        for(i,j),v in poly.items():
            for k in range(i+d+1):shifted[k,j]=shifted.get((k,j),F(0))+v*comb(i+d,k)*rc**(i+d-k)*scale**k
        bound=sum(abs(v)*z**j for(i,j),v in shifted.items())/lo**d
    else:
        for(i,j),v in poly.items():
            for a in range(i+1):
                for b in range(j+1):shifted[a,b]=shifted.get((a,b),F(0))+v*comb(i,a)*comb(j,b)*F(1,2)**(i-a+j-b)
        bound=sum(abs(v)*F(1,2)**(i+j)for(i,j),v in shifted.items())
    return min(F(21,500),bound),bound
def main():
    common=read(H.parent/'common_review.json')
    assert common['status']=='PASS' and all(v is True for v in common['checks'].values())
    if sha(SOURCE)!=EXPECTED or sha(OUT/'protocol_v2.json')!=PROTOCOL:raise ValueError('frozen wrapper changed')
    helper=C/'runs/exact_integration_bridge_v1/bulk_pilot_independent_replay/replay.py'
    if sha(helper)!='82f19f5aaf7a2ae986768fb2450c882e6b625cd6a03f5e0a6846ff215b927686':raise ValueError('independent jet helper changed')
    spec=importlib.util.spec_from_file_location('previous_independent_point_replay',helper);I=importlib.util.module_from_spec(spec);sys.modules[spec.name]=I;spec.loader.exec_module(I)
    import exact_fe_target_safe as SAFE_T
    I.T=SAFE_T  # Every independently formed jet now uses exact-zero-only arithmetic.
    ctx.prec=512;I.C.A.cache_clear();p=read(OUT/'protocol_v2.json');pre=read(OUT/'preflight_v2.json');inputs={};checks={};details=[]
    for name,digest in p['source_and_input_hashes'].items():
        if sha(C/name)!=digest:raise ValueError('frozen input differs '+name)
        key=(C/name).resolve().relative_to(C.resolve()).as_posix()
        if key in inputs and inputs[key]!=digest:raise ValueError('conflicting expected hash')
        inputs[key]=digest
    checks['all_current_source_protocol_preflight_inputs_bound']=pre['source_sha256']==EXPECTED and pre['protocol_sha256']==PROTOCOL and pre['status']=='PASS_SAFE_V2_PREPARATION_ONLY'
    checks['all144_selector_IDs_unique_and_complete']=len(p['selectors'])==144 and {x['element']for x in p['selectors']}==set(range(144))
    checks['point_scope_and_no_root_claim']='original saved B center' in p['scope']and p['root_transfer_proved']is False
    record=read(C/'runs/cohort_uniform_roots_v2/records/00_optoGEF2_small_col_B/result.json');target=I.T.Target(read(C/record['point_cache_file']));rule=target.auxiliary_rule()
    for row in p['selectors']:
        i=row['element'];f=target.element(i);jet=I.independent_jet(f['q'],f['B']);rho,bound=independent_rho(f,I.T.sub(jet['value'],I.T.ONE))
        # Independent simplex integral of affine r, with the same reviewed
        # normalized positive pointmass arithmetic. No force/H integral.
        rm=f['det']*sum(v[0]for v in f['corners'])/6
        massI=target.normalization()*I.C.A(rm);massQ=I.C.A(sum(rule['volume_weights'][i]));mass=F(str(I.C.upper(massI+massQ).fmpq()))
        maxg=max(I.majorant(x,f)for x in jet['gradient']);maxh=max(I.majorant(x,f)for x in jet['hessian_upper'].values());chosen=None
        for n in range(4,25):
            t1,t2=I.M.volumetric_derivative_tail_bounds(n,rho,1,5000);fb=mass*t1*maxg;hb=mass*(t2*maxg**2+t1*maxh)
            if fb<=F(1,10**22)and hb<=F(1,10**17):chosen=n;break
        good=chosen==row['degree']and rho==F(row['point_rho_exact'])and bound==F(row['independent_coefficient_majorant_exact'])and row['global_indices']==f['global_indices']
        # Last-bit upper enclosures can differ between valid arithmetic
        # associations; compare the declared exact tail values separately.
        stated_ok=F(row['force_tail_budget_upper_exact'])<=F(1,10**22)and F(row['H_entry_tail_budget_upper_exact'])<=F(1,10**17)
        checks[f'element_{i:03d}_independent_rho_smallest_degree_and_free_map']=good and stated_ok
        details.append({'element':i,'degree':chosen,'rho_exact':str(rho),'majorant_exact':str(bound),'force_tail_upper_exact':str(fb),'H_entry_tail_upper_exact':str(hb),'declared_bound_exactly_reproduced':fb==F(row['force_tail_budget_upper_exact'])and hb==F(row['H_entry_tail_budget_upper_exact']),'selector_PASS':bool(good and stated_ok)})
    # Sparse mapping with UNSORTED element-free indices, shared nodes, and both
    # positive/negative entries; compare with an independent dense double loop.
    blocks=[([3,0,2],[F(1),F(2),F(-1)],[[F(2),F(-1),F(1,2)],[F(-1),F(3),F(2)],[F(1,2),F(2),F(4)]]),([2,1],[F(-3),F(4)],[[F(5),F(-2)],[F(-2),F(6)]])]
    sparse={};g=[F(0)]*4;dense=np.full((4,4),F(0),dtype=object);wantg=[F(0)]*4
    for ids,force,hess in blocks:
        for k,v in zip(ids,force):g[k]+=v
        for i in range(len(ids)):
            for j in range(i,len(ids)):
                key=tuple(sorted((ids[i],ids[j])));sparse[key]=sparse.get(key,F(0))+hess[i][j]
        for i,gi in enumerate(ids):
            wantg[gi]+=force[i]
            for j,gj in enumerate(ids):dense[gi,gj]+=hess[i][j]
    actual=np.full((4,4),F(0),dtype=object)
    for(i,j),v in sparse.items():actual[i,j]=v;actual[j,i]=v
    checks['unsorted_sparse_upper_projection_matches_dense_all_entries']=np.array_equal(actual,dense)
    checks['shared_node_force_scatter_matches_dense']=g==wantg
    rows=[F(0)]*4
    for(i,j),v in sparse.items():
        rows[i]+=abs(v)
        if i!=j:rows[j]+=abs(v)
    checks['symmetric_row_norm_counts_offdiagonal_once']=rows==[sum(map(abs,row),F(0))for row in dense]
    R=[[F((-1)**(i+j)*(i+2),j+3)for j in range(4)]for i in range(4)];Rg=[sum(x*y for x,y in zip(row,g))for row in R]
    Rnorm=max(sum(map(abs,row),F(0))for row in R);Hnorm=max(rows);RH=[[sum(R[i][k]*dense[k,j]for k in range(4))for j in range(4)]for i in range(4)]
    checks['conservative_RH_row_norm_submultiplicative']=max(sum(map(abs,row),F(0))for row in RH)<=Rnorm*Hnorm
    checks['signed_RF_does_not_absolute_before_projection']=Rg==[F(-13,5),F(39,10),F(-26,5),F(13,2)]
    lo,hi=F(-2,3),F(5,7);ball=I.C.A(lo).union(I.C.A(hi));checks['endpoint_packet_union_outward']=F(str(ball.lower().fmpq()))<=lo and hi<=F(str(ball.upper().fmpq()))
    checks['all_predeclared_inputs_still_match']=all(sha(C/k)==v for k,v in inputs.items())and sha(SOURCE)==EXPECTED and sha(OUT/'protocol_v2.json')==PROTOCOL
    out={'status':'PASS'if all(checks.values())else'FAIL','checks':{**checks,**common['checks']},'common_review_sha256':sha(H.parent/'common_review.json'),'source_sha256':EXPECTED,'protocol_sha256':PROTOCOL,'preflight_sha256':sha(OUT/'preflight_v2.json'),'review_source_sha256':sha(Path(__file__)),'independent_jet_helper_sha256':sha(helper),'input_hashes':inputs,'selectors':details,'actual_144_force_or_H_integrals_run':False,'toy_RF_exact':list(map(str,Rg)),'scope':'Independent all144 POINT rho/smallest-order selector recomputation, plus unsorted sparse assembly/matrix-norm toy proofs. Reuses reviewed exacttarget/moment/tail components. Not actual144integration or a uniformtube/root proof.','nonblocking_caveat':'Local pilot moment-cache size diagnostics remain cumulative acrossinstances, not per-elementunique; no effect on bounds.'}
    assert not (H/'review.json').exists() and not (OUT/'prelaunch_independent_review.json').exists()
    (H/'review.json').write_text(json.dumps(out,indent=2)+'\n');(OUT/'prelaunch_independent_review.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'checks':len(checks),'failed':[k for k,v in checks.items()if not v],'receipt_sha256':sha(OUT/'prelaunch_independent_review.json')}))
    if not all(checks.values()):raise SystemExit(1)
if __name__=='__main__':main()
