"""Independent whole-edge rational domain certificates; no integrals."""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import hashlib,json,sys
from flint import arb,ctx
H=Path(__file__).resolve().parent;C=H.parents[2];OUT=C/'runs/exact_integration_bridge_v1/surface_assembly_B'
SOURCE=C/'scripts/exact_surface_bridge_assembly.py';EXPECTED='3d6c73ed55b710fe3ffa5b34a6834e0c3d3601c737c502348aaf75c248ae95ab';PROTOCOL='2bba6d5e11ead452ce195a035effb612ff726d38039998cfa957a368b53068c3';PREFLIGHT='1787b0f3409aa40b7e6525269775c6c78bbab55de57acbd33df258abe15f22f3'
sys.path.insert(0,str(C/'scripts'));import exact_surface_bridge_pilot as P
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def bernstein(poly,lo,hi):
    if any(j or i<0 for i,j in poly):raise ValueError('ordinary univariate polynomial required')
    n=max((i for i,j in poly),default=0);power=[F(0)]*(n+1)
    for(i,j),c in poly.items():
        for k in range(i+1):power[k]+=c*comb(i,k)*lo**(i-k)*(hi-lo)**k
    return[sum(power[j]*F(comb(k,j),comb(n,j))for j in range(k+1))for k in range(n+1)]
def prove(poly,lo,hi,depth=0):
    if min(bernstein(poly,lo,hi))>=0:return True,1
    if depth==8:return False,1
    mid=(lo+hi)/2;a,na=prove(poly,lo,mid,depth+1);b,nb=prove(poly,mid,hi,depth+1);return a and b,na+nb
def main():
    if sha(SOURCE)!=EXPECTED or sha(OUT/'protocol.json')!=PROTOCOL or sha(OUT/'preflight.json')!=PREFLIGHT:raise ValueError('frozen assembly differs')
    ctx.prec=512;P.core.A.cache_clear();p=read(OUT/'protocol.json');pre=read(OUT/'preflight.json');inputs={};checks={};details=[]
    for k,v in p['source_and_input_hashes'].items():
        if sha(C/k)!=v:raise ValueError('input changed '+k)
        key=(C/k).resolve().relative_to(C.resolve()).as_posix()
        if key in inputs and inputs[key]!=v:raise ValueError('conflicting canonical hash')
        inputs[key]=v
    checks['all24_cases_exactly_once']=[(x['face'],x['edge'])for x in p['cases']]==[(face,e)for face in[0,1]for e in range(12)]
    checks['N28_precision512_pointonly']=p['truncation_order']==28 and p['precision_bits']==512 and 'ONE saved B-center point' in p['scope']and p['root_transfer_proved']is False
    checks['prepared_status_and_checks_bound']=pre['status']=='PASS_PREPARATION_ONLY'and all(v is True for v in pre['checks'].values())and pre['protocol_sha256']==PROTOCOL and pre['source_sha256']==EXPECTED
    target=P.T.Target(read(C/p['point_cache_file']))
    for case in p['cases']:
        face,edge=case['face'],case['edge'];f=target.surface(face,edge);dom=case['domain'];leaves=dom['leaves'];lo=F(dom['G_lower_exact']);hi=F(dom['G_upper_exact']);rho=F(dom['rho_exact']);q0=F(dom['q0_exact']);hoop=F(dom['hoop_lower_exact']);cover=bool(leaves)and F(leaves[0]['s_exact'][0])==0 and F(leaves[-1]['s_exact'][1])==1 and all(F(x['s_exact'][1])==F(y['s_exact'][0])for x,y in zip(leaves,leaves[1:]));valid=cover and 0<lo<=hi and hoop>0 and q0==(lo+hi)/2 and rho==(hi-lo)/(hi+lo)and 0<=rho<1 and f['global_indices']==case['global_indices']
        a,b,c=f['q_numerators'];den=f['hoop_denominator'];den2=P.T.mul(den,den);Gnum=P.T.mul(P.T.mul(a,a),P.T.add(P.T.mul(b,b),P.T.mul(c,c)));nproof=0
        for leaf in leaves:
            ll,hh=map(F,leaf['s_exact']);gl,gh=map(F,leaf['G_bounds_exact']);hl=F(leaf['q_bounds_exact'][0][0]);valid &= 0<=ll<hh<=1 and gl>=lo and gh<=hi and hl>=hoop and min(bernstein(den,ll,hh))>0
            polynomials=[P.T.sub(a,P.T.scale(den,hl)),P.T.sub(Gnum,P.T.scale(den2,gl)),P.T.sub(P.T.scale(den2,gh),Gnum)]
            for polynomial in polynomials:
                ok,number=prove(polynomial,ll,hh);valid &= ok;nproof+=number
        checks[f'face{face}_edge{edge:02d}_independent_Bernstein_hoop_G_relative_domain']=bool(valid)
        details.append({'face':face,'edge':edge,'original_leaf_count':len(leaves),'Bernstein_inequality_subinterval_certificates':nproof,'q0_exact':str(q0),'rho_exact':str(rho),'PASS':bool(valid)})
    # Same reviewed bulk-wrapper scatter pattern: unsorted local maps and shared
    # nodes. Surface material directions/weights stay in the reviewed edgeengine.
    ids1=[3,0,2];ids2=[2,1];blocks=[(ids1,[1,-2,3],[[2,-1,4],[-1,3,5],[4,5,7]]),(ids2,[-4,6],[[11,-3],[-3,13]])];g=[0]*4;upper={};dense=[[0]*4 for _ in range(4)];gd=[0]*4
    for ids,gg,hh in blocks:
        for k,v in zip(ids,gg):g[k]+=v
        for i in range(len(ids)):
            for j in range(i,len(ids)):
                key=tuple(sorted((ids[i],ids[j])));upper[key]=upper.get(key,0)+hh[i][j]
        for i,gi in enumerate(ids):
            gd[gi]+=gg[i]
            for j,gj in enumerate(ids):dense[gi][gj]+=hh[i][j]
    checks['sparse_force_shared_nodes_signed']=g==gd
    checks['sparse_symmetric_upper_unsorted_maps']=all(upper.get(tuple(sorted((i,j))),0)==dense[i][j]for i in range(4)for j in range(4))
    val=P.A(F(-2,7)).union(P.A(F(3,11)));checks['serialized_interval_unpack_outward']=F(str(val.lower().fmpq()))<=F(-2,7)and F(3,11)<=F(str(val.upper().fmpq()))
    checks['literal_review_validator_positive_and_negative_controls']=False
    good={'status':'PASS','checks':{'ok':True},'source_sha256':EXPECTED,'protocol_sha256':PROTOCOL,'preflight_sha256':PREFLIGHT};P.validate_prelaunch(good,EXPECTED,PROTOCOL,PREFLIGHT)
    try:P.validate_prelaunch({**good,'checks':{'ok':1}},EXPECTED,PROTOCOL,PREFLIGHT)
    except ValueError:checks['literal_review_validator_positive_and_negative_controls']=True
    checks['sources_and_all_inputs_still_frozen']=sha(SOURCE)==EXPECTED and all(sha(C/k)==v for k,v in inputs.items())
    result={'status':'PASS'if all(checks.values())else'FAIL','checks':checks,'source_sha256':EXPECTED,'protocol_sha256':PROTOCOL,'preflight_sha256':PREFLIGHT,'review_source_sha256':sha(Path(__file__)),'input_hashes':inputs,'independent_domain_certificates':details,'actual_surface_integrals_run_by_reviewer':False,'proof':'ExactrationalBernstein coefficients certify numerator/denominator polynomial inequalities on a complete closedleafcover. Positivehoop andG implyvalidsqrtbranch; q0/rho formulas thenboundrelativeG onwholeedge.','scope':'PrelaunchALL24savedBcenterPOINTsurfacewrapper only. Reviewededgekernelreused, independentactual6edgereplayhandledseparately. No parameter/root tube/equilibriumtransfer/PDEclaim.'}
    (H/'review.json').write_text(json.dumps(result,indent=2)+'\n');(OUT/'prelaunch_independent_review.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':result['status'],'checks':len(checks),'failed':[k for k,v in checks.items()if not v],'receipt_sha256':sha(OUT/'prelaunch_independent_review.json')}))
    if not all(checks.values()):raise SystemExit(1)
if __name__=='__main__':main()
