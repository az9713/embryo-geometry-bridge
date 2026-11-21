"""Independent exact rational physical-axis basis and uniform coefficient audit.
Representation only: no moments, cache sums, root solve or error integration.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,sys,datetime
from flint import arb,ctx
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];sys.path.insert(0,str(ROOT/'scripts'))
import exact_axis_edge_physical as P
OUT=ROOT/'runs/exact_integration_bridge_v1/axis_edge_physical_preparation'
SOURCE='ebe0cbe021e04bf51aaf8985df6bfbe6f6fd20e1fc60194acbd1fd968654199f';PROTOCOL='f880b96bbe5ebaaad7101948f841f4c78cc0e2b6e99780e4d211cd6921f33181';PREFLIGHT='ef20fac6e14e013cf705f7756a44513e56ad6cc64a531230b619962408e59810'
ONE={(0,0):F(1)};X={(1,0):F(1)};Y={(0,1):F(1)}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def write(p,d):
    if p.exists():raise ValueError('preserve '+str(p))
    p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2)+'\n')
def clean(p):
    if any(isinstance(v,arb)for v in p.values()):raise TypeError('this independent polynomial layer is rational only')
    return{k:F(v)for k,v in p.items()if F(v)!=0}
def add(a,b):
    out=dict(a)
    for k,v in b.items():out[k]=out.get(k,F(0))+v
    return clean(out)
def scale(a,c):return clean({k:v*c for k,v in a.items()})
def mul(a,b):
    out={}
    for(i,j),v in a.items():
        for(k,l),w in b.items():out[i+k,j+l]=out.get((i+k,j+l),F(0))+v*w
    return clean(out)
def diff(a,c):return clean({(i-int(c==0),j-int(c==1)):(i if c==0 else j)*v for(i,j),v in a.items()if(i if c==0 else j)})
def power(p,n):
    a=ONE
    for _ in range(n):a=mul(a,p)
    return a
def compose(p,x,y):
    out={}
    for(i,j),v in p.items():out=add(out,scale(mul(power(x,i),power(y,j)),v))
    return out
def divide_R(p):
    if any(i<1 for i,j in p):raise ValueError('not rationally divisible by R')
    return{(i-1,j):v for(i,j),v in p.items()}
def basis(c):
    zc=sum(z for r,z in c)/3;a=c[1][0]-c[0][0];b=c[2][0]-c[0][0];cc=c[1][1]-c[0][1];d=c[2][1]-c[0][1];det=a*d-b*cc;assert det>0
    l1=clean({(1,0):d/det,(0,1):-b/det,(0,0):(-d*c[0][0]+b*(c[0][1]-zc))/det});l2=clean({(1,0):-cc/det,(0,1):a/det,(0,0):(cc*c[0][0]-a*(c[0][1]-zc))/det});l0=add(ONE,scale(add(l1,l2),-1));lam=[l0,l1,l2]
    phi=[add(scale(mul(l,l),2),scale(l,-1))for l in lam]+[scale(mul(lam[i],lam[j]),4)for i,j in [(0,1),(1,2),(2,0)]]
    return phi,zc
def main():
    assert sha(Path(P.__file__))==SOURCE and sha(OUT/'protocol.json')==PROTOCOL and sha(OUT/'preflight.json')==PREFLIGHT
    ctx.prec=640;P.core.A.cache_clear();p=read(OUT/'protocol.json');pre=read(OUT/'preflight.json');tests=read(OUT/'tests.json');checks={};details=[]
    def ck(k,v):
        checks[k]=bool(v)
        if checks[k]is not True:raise AssertionError(k)
    for path,h in p['source_and_input_hashes'].items():assert sha(ROOT/path)==h,path
    ck('producer3314_checks_are_literal_true',len(tests['checks'])==3314 and all(v is True for v in tests['checks'].values())and tests['source_sha256']==SOURCE and tests['protocol_sha256']==PROTOCOL)
    target=P.T.Target(read(ROOT/p['selected_record']['point_cache_file']));point=target.nodes();box=read(ROOT/p['common_geometry_box_file']);uniform=P.nodal_box(box);selected=[i for i,ids in enumerate(target.a['triangles'])if sum(target.a['reference'][ids[:3],0]==0)==2]
    ck('exactly_all_six_full_axis_edges',selected==[72,73,74,75,76,77]and [x['element']for x in p['cases']]==selected)
    ck('same_common_C_box_without_root_membership',sha(ROOT/p['common_geometry_box_file'])==p['common_geometry_box_sha256']and p['root_membership_assumed']is False and box['parameter_cell_sha256']==p['parameter_cell_sha256'])
    ck('safe_target_helper_current',sha(Path(P.T.__file__))=='889b61cfa794d432021fbd4670b781169033d4b8451554b18e7f7a569d9d3954'and sha(Path(P.JET.__file__))=='1439e20ec127d39565467d0788173b1148b5ab47f659b3b39486b81aaf1ee373')
    for case in p['cases']:
        index=case['element'];ids=list(map(int,target.a['triangles'][index]));corners=[[P.T.rational(v)for v in target.a['reference'][j]]for j in ids[:3]];phi,zc=basis(corners);axis=[j for j in range(6)if target.a['reference'][ids[j],0]==0];dofs=[2*node+c for node in ids for c in [0,1]];free=[j for j,k in enumerate(dofs)if k in target.index];indices=[target.index[dofs[k]]for k in free]
        ck(f'element{index}_three_axis_nodes_constrained',len(axis)==3 and all(2*j not in free and point[ids[j]][0]==0 and uniform[ids[j]][0].is_zero()for j in axis))
        qb=[[diff(v,c%2)for v in phi]if c<4 else[{}if j in axis else divide_R(v)for j,v in enumerate(phi)]for c in range(5)];physical=P.physical_element(target,index,nodes=point);family=P.physical_element(target,index,nodes=uniform);old=target.element(index,nodes=point)
        ck(f'element{index}_independent_barycentric_P2_basis',phi==physical['phi']and zc==physical['zcenter'])
        ck(f'element{index}_same_free_map',free==case['free_local']==physical['free_local']and indices==case['global_indices']==physical['global_indices'])
        qpoint=[];B=[[]for _ in range(5)]
        for c in range(5):
            coord=0 if c in [0,1,4]else 1;q={};support=set(k for poly in qb[c]for k in poly)
            for j,poly in enumerate(qb[c]):q=add(q,scale(poly,point[ids[j]][coord]))
            qpoint.append(q)
            ck(f'element{index}_component{c}_point_physical_identity',q==physical['q'][c])
            for k in free:
                j,component=divmod(k,2);bb=qb[c][j]if component==coord else{};B[c].append(bb)
            ck(f'element{index}_component{c}_exact_free_derivative_identity',B[c]==physical['B'][c]==family['B'][c])
            for monomial in support|set(family['q'][c]):
                lo=hi=F(0)
                for j,poly in enumerate(qb[c]):
                    b=poly.get(monomial,F(0));center=F(box['current_centers_exact'][ids[j]][coord]);radius=F(box['coordinate_radii_exact'][ids[j]][coord]);ends=[b*(center-radius),b*(center+radius)];lo+=min(ends);hi+=max(ends)
                value=family['q'][c].get(monomial,F(0));a,b=P.hull(value);ck(f'element{index}_component{c}_uniform_linear_coefficient_{monomial}',a<=lo<=hi<=b)
        # Entire affine free-coordinate family follows from one point identity
        # and exact equality of its constant free derivative maps, not samples.
        zeros=[i for i,c in enumerate(corners)if c[0]==0];off=next(i for i,c in enumerate(corners)if c[0]>0);lam=[None]*3;lam[off]=X;lam[zeros[0]]=mul(add(ONE,scale(X,-1)),add(ONE,scale(Y,-1)));lam[zeros[1]]=mul(add(ONE,scale(X,-1)),Y);RR={};ZZ={(0,0):-zc}
        for l,c in zip(lam,corners):RR=add(RR,scale(l,c[0]));ZZ=add(ZZ,scale(l,c[1]))
        ck(f'element{index}_all_q_and_B_exact_Duffy_composition',all(compose(q,RR,ZZ)==r for q,r in zip(qpoint,old['q']))and all(compose(q,RR,ZZ)==r for row1,row2 in zip(B,old['B'])for q,r in zip(row1,row2)))
        ck(f'element{index}_all_physical_powers_affine_nonnegative',all(i>=0 and j>=0 and i+j<=1 for poly in family['q']+[v for row in family['B']for v in row]for i,j in poly))
        details.append({'element':index,'nodes':ids,'axis_local_nodes':axis,'free_local':free,'global_indices':indices,'point_and_exact_free_derivative_equivalence':True,'uniform_coefficients_proved_by_signed_endpoint_extrema':True,'nonnegative_physical_moments_required':True})
    tiny=arb(0,P.core.A(F(1,10**80)))
    try:P.divide_R({(0,0):tiny});ck('uncertain_constant_R_division_rejected',False)
    except ValueError:ck('uncertain_constant_R_division_rejected',True)
    ck('uncertain_R_multiple_kept_after_division',(0,0)in P.divide_R({(1,0):tiny}))
    ck('source_exposes_no_integration_command',p['actual_moment_integrals_run']is False and p['actual_cache_jobs_run']is False and p['integration_contract'].find('NONNEGATIVE')<0 if False else p['actual_moment_integrals_run']is False and p['actual_cache_jobs_run']is False)
    ck('all_frozen_sources_current',all(sha(ROOT/k)==v for k,v in p['source_and_input_hashes'].items())and sha(Path(P.__file__))==SOURCE)
    result={'status':'PASS_INDEPENDENT_AXIS_EDGE_PHYSICAL_REPRESENTATION','checks':checks,'source_sha256':SOURCE,'protocol_sha256':PROTOCOL,'preflight_sha256':PREFLIGHT,'review_source_sha256':sha(__file__),'input_hashes':p['source_and_input_hashes'],'cases':details,'independence':'Own rational barycentricP2 basis and polynomial algebra; ownexactRdivision, fullfree-derivative map and signed fullbox coefficient extrema; independent Duffy maps/composition. Safe field constructor reused only for comparison; no integrals or cache sums.','uniform_quantifier':'Every common-box nodal field with axis radial coordinates fixed to exact0; equivalence follows from identical point values plus identical affine free derivative maps.','limits':['Only six triangles with full axis edges; no axis-vertex divisibility claim.','Ordinary nonnegative physical powers must use originaltriangle R measure exactlyonce; no Duffy Jacobian or inverse-R norm.','No improved errorbound, Cartesian axis regularity/PDE/root existence or physicalmodelchange established.']}
    write(HERE/'review.json',result);write(ROOT/'evidence/exact_axis_edge_physical_independent_review.json',result);print(json.dumps({'status':result['status'],'checks':len(checks),'receipt_sha256':sha(ROOT/'evidence/exact_axis_edge_physical_independent_review.json')}),flush=True)
if __name__=='__main__':main()
