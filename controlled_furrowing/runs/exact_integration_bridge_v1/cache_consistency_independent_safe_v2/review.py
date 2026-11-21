"""Independent third-order interval AD and two local uniform consistency checks."""
from pathlib import Path
from fractions import Fraction as F
from math import factorial
import hashlib,importlib.util,itertools,json,sys
import numpy as np
from flint import arb,arb_mat,ctx
HERE=Path(__file__).resolve().parent;CF=HERE.parents[2];OUT=CF/'runs/exact_integration_bridge_v1/cache_consistency_preparation'
sys.path.insert(0,str(CF/'scripts'));import exact_cache_consistency as P
EXPECTED='d363567cf2efbf34b2b3b39fb5282dc7dcabced1b3b96a9717fbcf382d5f17d6'
ZERO=(0,0,0,0,0)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def clean(d):
    out={}
    for key,value in d.items():
        if isinstance(value,arb):
            if not value.is_finite():raise ValueError('finite AD coefficient required')
            exact_zero=value.is_zero()
        else:exact_zero=F(value)==0
        if not exact_zero:out[key]=value
    return out
def add(a,b):
    out=dict(a)
    for k,v in b.items():out[k]=out.get(k,arb(0))+v
    return clean(out)
def scale(a,b):return clean({k:v*b for k,v in a.items()})
def mul(a,b):
    out={}
    for x,v in a.items():
        for y,w in b.items():
            z=tuple(i+j for i,j in zip(x,y))
            if sum(z)<=3:out[z]=out.get(z,arb(0))+v*w
    return clean(out)
def constant(v):return{ZERO:P.A(v)}
def variable(v,i):
    key=list(ZERO);key[i]=1;return{ZERO:P.A(v),tuple(key):arb(1)}
def log(a):
    base=a.get(ZERO,arb(0))
    if not base>0:raise ValueError('ADlog entirebox mustbe positive')
    d=dict(a);d.pop(ZERO,None);out=constant(base.log());power={ZERO:arb(1)}
    for n in range(1,4):power=mul(power,d);out=add(out,scale(power,(-1)**(n+1)/(P.A(n)*base**n)))
    return out
def energy_jet(q,mu,kappa,drive):
    x=[variable(v,i)for i,v in enumerate(q)];J=mul(add(mul(x[0],x[3]),scale(mul(x[1],x[2]),-1)),x[4]);ell=log(J)
    quad={}
    for v in x:quad=add(quad,mul(v,v))
    W=add(scale(quad,P.A(mu)/2),add(scale(ell,-P.A(mu)),scale(mul(ell,ell),P.A(kappa)/2)))
    return add(W,scale(add(mul(x[1],x[1]),mul(x[3],x[3])),P.A(drive)/2))
def derivative(jet,indices):
    multi=tuple(indices.count(i)for i in range(5));factor=1
    for n in multi:factor*=factorial(n)
    return jet.get(multi,arb(0))*factor
def matrix(rows):return arb_mat(len(rows),len(rows[0]),[v for row in rows for v in row])
def hessian(q,drive):
    jet=energy_jet(q,1,5000,drive);return matrix([[derivative(jet,[i,j])for j in range(5)]for i in range(5)])
def pack(x):return P.core.pack(x)
def overlaps_packet(x,p):return F(str(x.lower().fmpq()))<=F(p['upper_exact'])and F(p['lower_exact'])<=F(str(x.upper().fmpq()))
def adversarial_tests():
    checks={};x=arb(0,1);half=P.A(F(1,2))
    checks['safe_zero_crossing_constant_retained']=ZERO in clean({ZERO:x}) and clean({ZERO:x})[ZERO].contains(half)
    checks['safe_exact_zero_removed']=clean({ZERO:arb(0)})=={}
    checks['safe_fraction_exact_zero_removed']=clean({ZERO:F(0)})=={}
    checks['safe_nonzero_fraction_retained']=clean({ZERO:F(2)})[ZERO]==2
    checks['safe_uncertain_addition_retained']=ZERO in add({ZERO:arb(1)},{ZERO:arb(-1,1)}) and add({ZERO:arb(1)},{ZERO:arb(-1,1)})[ZERO].contains(half)
    checks['safe_exact_cancellation_removed']=add({ZERO:arb(1)},{ZERO:arb(-1)})=={}
    checks['safe_uncertain_scaling_retained']=scale({ZERO:x},2)[ZERO].contains(P.A(1))
    checks['safe_uncertain_multiplication_retained']=mul({ZERO:x},{ZERO:arb(1)})[ZERO].contains(half)
    power={ZERO:arb(1)}
    for n in range(1,7):
        power=mul(power,{ZERO:x});checks['safe_uncertain_power_'+str(n)]=ZERO in power and power[ZERO].contains(P.A(F(1,2)**n))
    key=(1,0,0,0,0)
    checks['safe_uncertain_derivative_coefficient_retained']=derivative(clean({key:x}),[0]).contains(half)
    for label,value in [('nan',arb('nan')),('positive_inf',arb('+inf')),('negative_inf',arb('-inf'))]:
        try:clean({ZERO:value});checks['safe_nonfinite_rejected_'+label]=False
        except ValueError:checks['safe_nonfinite_rejected_'+label]=True
    assert all(v is True for v in checks.values()),checks
    return checks

def main():
    safe_protocol=read(HERE/'protocol.json')
    assert sha(Path(__file__))==safe_protocol['review_source_sha256']
    assert all(sha(CF/k)==v for k,v in safe_protocol['input_hashes'].items())
    if sha(Path(P.__file__))!=EXPECTED or sha(OUT/'protocol.json')!='0c7076d2b62653bbc67b12056a86cfea8f9d04844a6c7f4f685ca3d04ed574cc'or sha(OUT/'preflight.json')!='8841ab3f7576483bd515a43ad4b598a0cc5867812b77da72b0e4a285e7ec21a1':raise ValueError('frozen component differs')
    ctx.prec=384;P.core.A.cache_clear();protocol=read(OUT/'protocol.json');pre=read(OUT/'preflight.json');checks=adversarial_tests();bindings={}
    for name,d in protocol['source_and_input_hashes'].items():
        if sha(CF/name)!=d:raise ValueError('frozen input differs')
        bindings[name]=d
    row=protocol['selected_record'];target=P.target_core.Target(read(CF/row['point_cache_file']))
    with np.load(CF/row['predictor_file'],allow_pickle=False)as z:S=z['S'].copy()
    helper=CF/'runs/exact_integration_bridge_v1/bulk_pilot_independent_replay/replay.py';spec=importlib.util.spec_from_file_location('independent_point_maps_for_uniform_review',helper);I=importlib.util.module_from_spec(spec);sys.modules[spec.name]=I;spec.loader.exec_module(I)
    if sha(helper)!='82f19f5aaf7a2ae986768fb2450c882e6b625cd6a03f5e0a6846ff215b927686':raise ValueError('independent maphelper differs')
    # All125 thirds via direct interval AD of the ENERGY, not source chain rule.
    point=[F(11,10),F(1,20),F(-1,25),F(9,10),F(21,20)];jet=energy_jet(list(map(P.A,point)),1,5000,F(3,7));tensor,_=P.bulk_third_tensor(point,1,5000)
    for i,j,k in itertools.product(range(5),repeat=3):checks[f'direct_energy_AD_third_{i}{j}{k}']=derivative(jet,[i,j,k]).overlaps(tensor[i][j][k])
    # Independent exact endpoint construction for the OLD entire root cube.
    deltas=[(F(row['parameter_cell'][k][0])-F(row['base_parameters_exact'][k]),F(row['parameter_cell'][k][1])-F(row['base_parameters_exact'][k]))for k in P.ORDER];radius=F(row['selected_corrector_radius_exact']);Ybounds=[[ (F.from_float(float(v)),F.from_float(float(v)))for v in pair]for pair in target.a['reference']]
    for i,k in enumerate(target.a['free']):
        terms=[(F.from_float(float(v))*lo,F.from_float(float(v))*hi)for v,(lo,hi)in zip(S[i],deltas)];base=F.from_float(float(target.a['u'][i]))+F.from_float(float(target.a['reference'][int(k)//2,int(k)%2]));Ybounds[int(k)//2][int(k)%2]=(base+sum(min(v)for v in terms)-radius,base+sum(max(v)for v in terms)+radius)
    fixtures=[]
    for item in pre['two_actual_local_uniform_fixtures']:
        file=OUT/item['file']
        if sha(file)!=item['sha256']:raise ValueError('actual local fixture changed')
        old=read(file);element,iq=old['element'],old['quadrature_index'];ids=list(map(int,target.a['triangles'][element]));corners=[[F.from_float(float(v))for v in p]for p in target.a['reference'][ids[:3]]];lam=target.auxiliary_rule()['volume_barycentric'][iq];N,D,rr=I.direct_basis(corners,lam);star=[[F(0)]*12 for _ in range(5)];cache=[[F(0)]*12 for _ in range(5)]
        for k in range(6):
            star[0][2*k]=D[k][0];star[1][2*k]=D[k][1];star[2][2*k+1]=D[k][0];star[3][2*k+1]=D[k][1];star[4][2*k]=N[k]/rr
            dx=F.from_float(float(target.a['shape_grad'][element,iq,k,0]));dz=F.from_float(float(target.a['shape_grad'][element,iq,k,1]));cache[0][2*k]=dx;cache[1][2*k]=dz;cache[2][2*k+1]=dx;cache[3][2*k+1]=dz;cache[4][2*k]=F.from_float(float(target.a['shape'][iq,k]))/F.from_float(float(target.a['qr'][element,iq]))
        maps=P.exact_bulk_maps(target,element,iq);checks[f'fixture{element}_{iq}_independent_fullmaps']=maps['star_full']==star and maps['cache_full']==cache
        selected=[k for k in range(12)if 2*ids[k//2]+k%2 in target.index];localbounds=[Ybounds[node][c]for node in ids for c in range(2)];Y=[P.A(a).union(P.A(b))for a,b in localbounds];delta=[[a-b for a,b in zip(x,y)]for x,y in zip(star,cache)]
        qc=[sum((P.A(v)*y for v,y in zip(b,Y)),arb(0))for b in cache];dq=[sum((P.A(v)*y for v,y in zip(b,Y)),arb(0))for b in delta];line=[a+arb(0).union(b)for a,b in zip(qc,dq)];driveLo,driveHi=[F(x)*F(row['dose_exact'])*F.from_float(float(target.a['active'][element,0]))for x in row['parameter_cell']['common_gain']];drive=P.A(driveLo).union(P.A(driveHi))
        # Separate four-term expansion plus independently differentiated W.
        J=(line[0]*line[3]-line[1]*line[2])*line[4];checks[f'fixture{element}_{iq}_convex_AD_domain']=J>P.A(F(1,1000000))and line[4]>0 and(line[0]*line[3]-line[1]*line[2])>0
        ad=energy_jet(line,1,5000,drive);third=[[[derivative(ad,[i,j,k])for k in range(5)]for j in range(5)]for i in range(5)];dH=matrix([[sum((third[i][j][k]*dq[k]for k in range(5)),arb(0))for j in range(5)]for i in range(5)]);Hc=hessian(qc,drive);Bs=matrix([[P.A(b[k])for k in selected]for b in star]);Bc=matrix([[P.A(b[k])for k in selected]for b in cache]);dB=matrix([[P.A(b[k])for k in selected]for b in delta]);bound=dB.transpose()*Hc*Bc+Bc.transpose()*Hc*dB+dB.transpose()*Hc*dB+Bs.transpose()*dH*Bs;n=len(selected)
        checks[f'fixture{element}_{iq}_all_uniform_H_entries_overlap']=all(overlaps_packet(bound[i,j],old['unweighted_local_free_H_difference'][i][j])for i in range(n)for j in range(n))
        witnesses=0
        for mode in range(4):
            values=[lo if mode==0 else hi if mode==1 else(lo+hi)/2 if mode==2 else(lo if k%2 else hi)for k,(lo,hi)in enumerate(localbounds)]
            for active in[driveLo,driveHi]:
                qs=[P.A(sum(v*y for v,y in zip(b,values)))for b in star];qcache=[P.A(sum(v*y for v,y in zip(b,values)))for b in cache];actual=Bs.transpose()*hessian(qs,active)*Bs-Bc.transpose()*hessian(qcache,active)*Bc
                checks[f'fixture{element}_{iq}_direct_point_witness{witnesses}']=all(actual[i,j].overlaps(bound[i,j])and overlaps_packet(actual[i,j],old['unweighted_local_free_H_difference'][i][j])for i in range(n)for j in range(n));witnesses+=1
        fixtures.append({'element':element,'iq':iq,'local_dofs':n,'independent_uniform_H_difference':[[pack(bound[i,j])for j in range(n)]for i in range(n)],'delta_q':[pack(v)for v in dq],'convex_line_J':pack(J),'point_witnesses':witnesses,'old_fixture_sha256':item['sha256']})
    # Reproduce, do not silently suppress, the nonblocking integer-index guard.
    rule_typeerror=False
    try:P.exact_bulk_maps(target,0,True)
    except TypeError:rule_typeerror=True
    checks['boolean_rule_bypasses_explicit_guard_then_TypeError']=rule_typeerror
    element_typeerror=False
    try:P.exact_bulk_maps(target,True,0)
    except TypeError:element_typeerror=True
    checks['boolean_element_bypasses_explicit_guard_then_TypeError']=element_typeerror
    checks['source_and_original_inputs_still_frozen']=sha(Path(P.__file__))==EXPECTED and all(sha(CF/k)==v for k,v in bindings.items())
    result={'status':'PASS_SAFE_V2_UNIFORM_LOCAL_BULK_CACHE_COMPONENT_WITH_GUARD_FINDING'if all(checks.values())else'FAIL_COMPONENT_REVIEW','checks':checks,'review_protocol_sha256':sha(HERE/'protocol.json'),'supersedes_flawed_independent_review_sha256':safe_protocol['old_flawed_receipt_sha256'],'unsafe_old_T_polynomial_operations_used_with_Arb':False,'source_sha256':EXPECTED,'protocol_sha256':sha(OUT/'protocol.json'),'preflight_sha256':sha(OUT/'preflight.json'),'review_source_sha256':sha(Path(__file__)),'independent_map_helper_sha256':sha(helper),'precision_bits':384,'fixtures':fixtures,'input_hashes':bindings,'independent_method':'Generic third-order intervalTaylorAD ofscalarenergy; independentbarycentric maps/fullcubeendpointconstruction; four-term matrixexpansion ratherthan producerthree-term grouping; sixteen directpointwitnesses.','finding':{'kind':'NONBLOCKING_API_GUARD','description':'BothBooleanindices bypass explicitintegerguard andlaterNumPyindexing triggersshape-relatedTypeError; neitherreturnedvalidmaps. Actualdeclaredintegerfixtures unaffected; hardenfutureguard withoutinvalidatingfrozenvalidinputmath.'},'scope':'Two UNWEIGHTED local auxiliary sites over supplied old-root nodalcube/common-driveinterval. No spatialintegration, weightedsum/full144cachebound,surface/uniformexact-integralroottransfer.'}
    assert all(sha(CF/k)==v for k,v in safe_protocol['input_hashes'].items())
    assert not (HERE/'results.json').exists() and not (CF/'evidence/exact_cache_consistency_independent_review_v2.json').exists()
    (HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n');(CF/'evidence/exact_cache_consistency_independent_review_v2.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':result['status'],'checks':len(checks),'failed':[k for k,v in checks.items()if not v],'receipt_sha256':sha(CF/'evidence/exact_cache_consistency_independent_review_v2.json')}))
    if not all(checks.values()):raise SystemExit(1)
if __name__=='__main__':main()
