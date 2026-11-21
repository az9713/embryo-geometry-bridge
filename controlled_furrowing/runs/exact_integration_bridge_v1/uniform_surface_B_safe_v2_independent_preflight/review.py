"""Independent preparation-only interval domains/maps and shared-K/T review.
NO actual edge integrals or weighted cache sums are computed here.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import hashlib,json,sys,datetime
from flint import arb
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];sys.path.insert(0,str(ROOT/'scripts'))
import exact_uniform_surface_bridge_B_safe as U
OUT=ROOT/'runs/exact_integration_bridge_v1/uniform_surface_B_safe_v2'
EXPECTED_SOURCE='61cf47d2b1a9e7ee46f9fa680c24b4fdf7b4e31e4889379c56951ab6a7455ee9'
EXPECTED_PROTOCOL='91d65ac402c763f3eadd12dc3bee18ab76a38ab3aef1c5b19f90a825d69badfa'
EXPECTED_PREFLIGHT='d0d49afb6f362ade7cee86e3ea45774bf918cd5c6016ff87af7eabe55ca67d0f'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def write(p,d):
    if p.exists():raise ValueError('preserve '+str(p))
    p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2)+'\n')
def add(a,b):return a[0]+b[0],a[1]+b[1]
def scale(k,a):return(k*a[0],k*a[1])if k>=0 else(k*a[1],k*a[0])
def divide(a,b):
    assert b[0]>0
    v=[x/y for x in a for y in b];return min(v),max(v)
def square(a):return(F(0)if a[0]<=0<=a[1]else min(x*x for x in a),max(x*x for x in a))
def contains(a,b):return a[0]<=b[0]<=b[1]<=a[1]
def hull(v):return F(str(v.lower().fmpq())),F(str(v.upper().fmpq()))
def dot(row,Y):
    x=(F(0),F(0))
    for c,y in zip(row,Y):x=add(x,scale(c,y))
    return x
def bernstein(coeffs,lo,hi):
    # Direct single linear transform, rather than source two-stage power array.
    n=max(coeffs,default=0);out=[]
    for i in range(n+1):
        b=(F(0),F(0))
        for k,c in coeffs.items():
            factor=sum((F(comb(k,j)*comb(i,j),comb(n,j))*lo**(k-j)*(hi-lo)**j for j in range(min(k,i)+1)),F(0))
            b=add(b,scale(factor,c))
        out.append(b)
    return out
def raw_numerators(reference,Y):
    lo,hi=reference[0][0],reference[2][0];L=hi-lo;assert L>0
    phi=[[F(1),F(-3),F(2)],[F(0),F(4),F(-4)],[F(0),F(-1),F(2)]]
    radial=[dot([phi[j][i]for j in range(3)],[Y[j][0]for j in range(3)])for i in range(3)];vertical=[dot([phi[j][i]for j in range(3)],[Y[j][1]for j in range(3)])for i in range(3)]
    if lo==0:
        assert radial[0]==(0,0)
        a={i:scale(1/L,x)for i,x in enumerate(radial[1:])};den={0:(F(1),F(1))}
    else:a=dict(enumerate(radial));den={0:(lo,lo),1:(L,L)}
    return[a,{i:scale(F(i+1)/L,x)for i,x in enumerate(radial[1:])},{i:scale(F(i+1)/L,x)for i,x in enumerate(vertical[1:])}],den
def actual_maps(target,face,edge,iq):
    a=target.a;ids=list(map(int,a[f'surface_{face}_edges'][edge]));ref=[[U.T.rational(v)for v in a['reference'][j]]for j in ids];lo,hi=ref[0][0],ref[2][0];L=hi-lo
    oldN=list(map(U.T.rational,a[f'surface_{face}_N'][iq]));oldD=list(map(U.T.rational,a[f'surface_{face}_dN'][iq]));s=oldN[2]+oldN[1]/2;assert 0<s<1
    Ns=[(1-s)*(1-2*s),4*s*(1-s),s*(2*s-1)];Ds=[-3+4*s,4-8*s,-1+4*s];R=lo+L*s;r0=U.T.rational(a[f'surface_{face}_r0'][edge,iq]);l0=U.T.rational(a[f'surface_{face}_l0'][edge]);assert R>0 and r0>0 and l0>0
    D=[L/(r0*l0),1/L,1/L];star=[[F(0)]*6 for _ in range(3)];cache=[[F(0)]*6 for _ in range(3)]
    for j in range(3):
        star[0][2*j]=Ns[j]/R;star[1][2*j]=Ds[j]/L;star[2][2*j+1]=Ds[j]/L;cache[0][2*j]=D[0]*oldN[j];cache[1][2*j]=D[1]*oldD[j];cache[2][2*j+1]=D[2]*oldD[j]
    selected=[k for k in range(6)if 2*ids[k//2]+k%2 in target.index];indices=[target.index[2*ids[k//2]+k%2]for k in selected];weight=U.T.rational(a[f'surface_{face}_weights'][edge,iq])/target.energy_scale
    return ids,star,cache,D,selected,indices,weight
def main():
    assert sha(Path(U.__file__))==EXPECTED_SOURCE and sha(OUT/'protocol.json')==EXPECTED_PROTOCOL and sha(OUT/'preflight.json')==EXPECTED_PREFLIGHT
    U.set_precision(640);p=read(OUT/'protocol.json');pre=read(OUT/'preflight.json');checks={};details=[];site_details=[]
    def ck(k,v):
        checks[k]=bool(v)
        if checks[k]is not True:raise AssertionError(k)
    for k,v in p['source_and_input_hashes'].items():assert sha(ROOT/k)==v,k
    ck('all48_preparation_checks_literal_true',len(pre['checks'])==48 and all(v is True for v in pre['checks'].values()))
    ck('all24_case_keys_exactly_once',[(c['face'],c['edge'])for c in p['cases']]==[(f,e)for f in [0,1]for e in range(12)])
    ck('same640_N28_uniform_scope',p['precision_bits']==640 and p['truncation_order']==28 and p['root_membership_assumed']is False)
    target=U.T.Target(read(ROOT/p['point_cache_file']));box=read(ROOT/p['common_geometry_box_file']);nodes=U.common_nodes(target,box)
    ck('same_common_box_and_cell',sha(ROOT/p['common_geometry_box_file'])==p['common_geometry_box_sha256']and box['parameter_cell_sha256']==p['parameter_cell_sha256']==p['selected_record']['parameter_cell_sha256'])
    Yall=[[(F(c)-F(r),F(c)+F(r))for c,r in zip(cc,rr)]for cc,rr in zip(box['current_centers_exact'],box['coordinate_radii_exact'])]
    ck('all_axis_nodes_exact_zero_and_other_radial_boxes_positive',all((Yall[j][0]==(0,0)and nodes[j][0].is_zero())if target.a['reference'][j,0]==0 else Yall[j][0][0]>0 and nodes[j][0]>0 for j in range(325)))
    ck('safe_arithmetic_target_and_helper_hashes',sha(Path(U.T.__file__))=='889b61cfa794d432021fbd4670b781169033d4b8451554b18e7f7a569d9d3954'and sha(Path(U.JET.__file__))=='1439e20ec127d39565467d0788173b1148b5ab47f659b3b39486b81aaf1ee373')
    local=read(ROOT/'evidence/exact_surface_cache_consistency_independent_review.json');ck('local_cache_derivative_review_current',local['status']=='PASS_INDEPENDENT_LOCAL_SURFACE_CACHE_CONSISTENCY_REVIEW'and all(v is True for v in local['checks'].values())and local['source_sha256']==sha(Path(U.L.__file__)))
    for case in p['cases']:
        face,edge=case['face'],case['edge'];ids=list(map(int,target.a[f'surface_{face}_edges'][edge]));ref=[[U.T.rational(v)for v in target.a['reference'][j]]for j in ids];assert len({z[1]for z in ref})==1
        f=target.surface(face,edge,nodes=nodes);raw,rawden=raw_numerators(ref,[Yall[j]for j in ids]);coeff=[{i:hull(v)if isinstance(v,arb)else(F(v),F(v))for(i,j),v in poly.items()}for poly in f['q_numerators']];den={i:(F(v),F(v))for(i,j),v in f['hoop_denominator'].items()}
        ck(f'edge{face}_{edge}_raw_P2_interval_coefficients_enclosed',all(contains(cc.get(k,(0,0)),v)for cc,rr in zip(coeff,raw)for k,v in rr.items())and all(den.get(k,(0,0))==v for k,v in rawden.items()))
        dom=case['domain'];leaves=dom['leaves'];cover=bool(leaves)and F(leaves[0]['s_exact'][0])==0 and F(leaves[-1]['s_exact'][1])==1 and all(F(a['s_exact'][1])==F(b['s_exact'][0])for a,b in zip(leaves,leaves[1:]));valid=cover
        lower=[];upper=[];hoops=[]
        for leaf in leaves:
            lo,hi=map(F,leaf['s_exact']);bs=[bernstein(c,lo,hi)for c in coeff];bd=bernstein(den,lo,hi);db=(min(x[0]for x in bd),max(x[1]for x in bd));assert db[0]>0
            qr=[(min(x[0]for x in row),max(x[1]for x in row))for row in bs];qr[0]=divide(qr[0],db);sq=list(map(square,qr));G=(sq[0][0]*(sq[1][0]+sq[2][0]),sq[0][1]*(sq[1][1]+sq[2][1]));lower.append(G[0]);upper.append(G[1]);hoops.append(qr[0][0])
            valid=valid and [list(map(str,x))for x in qr]==leaf['q_bounds_exact']and list(map(str,G))==leaf['G_bounds_exact']and [[list(map(str,x))for x in row]for row in bs]==leaf['q_numerator_Bernstein_coefficients_exact']and [list(map(str,x))for x in bd]==leaf['hoop_denominator_Bernstein_coefficients_exact']
        gl,gh=min(lower),max(upper);rho=(gh-gl)/(gh+gl);valid=valid and gl>0 and min(hoops)>0 and(F(dom['q0_exact'])==(gl+gh)/2)and F(dom['rho_exact'])==rho<1 and list(map(str,(gl,gh)))==dom['G_bounds_exact']and F(dom['hoop_lower_exact'])==min(hoops)
        ck(f'edge{face}_{edge}_all_closed_leaf_Bernstein_domains_exactly_reproduced',valid)
        prefix='basal'if face==0 else'apical';Klo,Khi=map(F,p['selected_record']['parameter_cell'][prefix+'_area_modulus']);Tlo,Thi=map(F,p['selected_record']['parameter_cell'][prefix+'_surface_tension']);ck(f'edge{face}_{edge}_positive_K_nonnegative_T_same_C',0<Klo<=Khi and 0<=Tlo<=Thi)
        for iq,frozen in enumerate(case['cache_site_domains']):
            ids2,star,cache,D,selected,indices,weight=actual_maps(target,face,edge,iq);assert ids2==ids
            yy=[Yall[j][c]for j in ids for c in [0,1]];delta=[[a-b for a,b in zip(s,c)]for s,c in zip(star,cache)];qc=[dot(row,yy)for row in cache];dq=[dot(row,yy)for row in delta];line=[(v[0]+min(F(0),dd[0]),v[1]+max(F(0),dd[1]))for v,dd in zip(qc,dq)];n=add(square(line[1]),square(line[2]))
            good=frozen['quadrature_index']==iq and frozen['nodes']==ids and frozen['selected_local']==selected and frozen['global_indices']==indices==case['global_indices']and frozen['normalization_D_exact']==list(map(str,D))and frozen['auxiliary_normalized_weight_exact']==str(weight)and weight>0 and frozen['cache_q_exact']==[list(map(str,x))for x in qc]and frozen['delta_q_exact']==[list(map(str,x))for x in dq]and frozen['convex_line_exact']==[list(map(str,x))for x in line]and line[0][0]>0 and n[0]>0 and F(frozen['tangent_norm_squared_lower_exact'])==n[0]
            ck(f'edge{face}_{edge}_site{iq}_independent_fullY_delta_and_old_normalization',good)
            site_details.append({'face':face,'edge':edge,'site':iq,'same_exact_normalization_and_fullY_line':True})
        ck(f'edge{face}_{edge}_twelve_sites_and_free_projection',len(case['cache_site_domains'])==12 and f['global_indices']==case['global_indices'])
        details.append({'face':face,'edge':edge,'leaf_count':len(leaves),'rho_exact':str(rho),'direct_P2_and_interval_Bernstein_proof':True})
    # Independent shared-K/T endpoint tests of the uniform series and its
    # two scalar derivative tails. This is synthetic algebra, not field I/Q.
    K=U.iball((F(1),F(3)));T=U.iball((F(0),F(4)));polyA,polyB=U.uniform_series({(1,0):F(1,100)},F(4),4,K,T)
    ck('shared_KT_zero_crossing_prefactor_retained',bool(polyA)and bool(polyB)and any(v.contains(0)and not v.is_zero()for v in polyB.values()))
    for k in [F(1),F(3)]:
        for t in [F(0),F(4)]:
            for s in [F(0),F(1,2),F(1)]:
                e=s/100;G=4*(1+U.A(e));aa=U.A(k/2)+U.A(t-k)/(2*G.sqrt());bb=-U.A(t-k)/(4*G*G.sqrt());trA=U.A(k/2)+U.A((t-k)/4)*sum((U.A(c*e**n)for n,c in enumerate(U.S.binomial_coefficients(F(-1,2),4))),arb(0));trB=-U.A((t-k)/32)*sum((U.A(c*e**n)for n,c in enumerate(U.S.binomial_coefficients(F(-3,2),4))),arb(0));ta=U.A(abs(t-k)*U.S.binomial_tail(F(-1,2),4,F(1,100))/4);tb=U.A(abs(t-k)*U.S.binomial_tail(F(-3,2),4,F(1,100))/32)
                ck(f'shared_KT_corner{k}_{t}_{s}_coefficient_enclosure_and_tails',U.T.evaluate(polyA,s,0).contains(trA)and U.T.evaluate(polyB,s,0).contains(trB)and U.core.upper(aa-trA)<=ta and U.core.upper(bb-trB)<=tb)
    ck('no_actual_integrals_or_weighted_cache_sums_started',not(OUT/'started.json').exists()and not(OUT/'results.json').exists())
    ck('all_source_input_hashes_still_frozen',all(sha(ROOT/k)==v for k,v in p['source_and_input_hashes'].items())and sha(Path(U.__file__))==EXPECTED_SOURCE)
    result={'status':'PASS','checks':checks,'source_sha256':EXPECTED_SOURCE,'protocol_sha256':EXPECTED_PROTOCOL,'preflight_sha256':EXPECTED_PREFLIGHT,'review_source_sha256':sha(__file__),'source_and_input_hashes':p['source_and_input_hashes'],'whole_edge_details':details,'local_site_details':site_details,'independent_method':'Direct rawP2 interval coefficient hulls; single-transform exact rational intervalBernstein coefficients; independently normalized star/cache full maps and SAME-Y coefficientwise deltas at all288sites; sharedmaterial endpoint polynomial/tail tests. Reviewed safe field/moment/jet/local-third components reused explicitly.','scope':'Prelaunch uniformSURFACEcomponent only. No actual24edge integration or weighted cache sums, fullbulk+surface/error/root/PDE/biology certificate.'}
    assert all(v is True for v in checks.values())
    write(HERE/'review.json',result);write(OUT/'prelaunch_independent_review.json',result);print(json.dumps({'status':'PASS','checks':len(checks),'leaves':sum(x['leaf_count']for x in details),'sites':len(site_details),'receipt_sha256':sha(OUT/'prelaunch_independent_review.json')}),flush=True)
if __name__=='__main__':main()
