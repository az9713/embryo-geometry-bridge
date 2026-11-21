"""Cheap local algebra/support study; never integrates a full field or solves."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import hashlib,json,sys,time,math
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'scripts'))
import exact_fe_target as T
OUT=Path(__file__).parent
EXPECTED='4fd89ddf23c59050a92f18f7b4b6e69ef8e5d08e57ba0c1949cf0b61b2c869a4'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def norm(p):return sum(map(abs,p.values()),F(0))
def summarize(p):
    if not p:return {'support':0,'radial_range':None,'axial_range':None,'max_total_degree':None,'coefficient_max_bits':0}
    return {'support':len(p),'radial_range':[min(i for i,j in p),max(i for i,j in p)],'axial_range':[min(j for i,j in p),max(j for i,j in p)],'max_total_degree':max(i+j for i,j in p),'coefficient_max_bits':max(max(c.numerator.bit_length(),c.denominator.bit_length())for c in p.values())}
def groups(ps):
    return {'count':len(ps),'nonzero_count':sum(bool(p)for p in ps),'support_sum':sum(map(len,ps)),'support_max':max(map(len,ps),default=0),'union_support':summarize({k:F(1)for p in ps for k in p}),'support_histogram':dict(sorted(Counter(map(len,ps)).items()))}
def minkowski(a,b):return {(i+k,j+l)for i,j in a for k,l in b}
def support_growth(base,other=None):
    out=[];s={(0,0)};other=set(other or {(0,0)})
    for n in range(17):
        out.append({'power':n,'power_support_upper':len(s),'with_factor_support_upper':len(minkowski(s,other))})
        s=minkowski(s,set(base))
    return out

def jacobian_derivatives(q,B):
    a,b,c,d,h=q;det=T.sub(T.mul(a,d),T.mul(b,c));J=T.mul(det,h);n=len(B[0]);du=[]
    for k in range(n):du.append(T.sub(T.add(T.mul(B[0][k],d),T.mul(a,B[3][k])),T.add(T.mul(B[1][k],c),T.mul(b,B[2][k]))))
    Ju=[T.add(T.mul(v,h),T.mul(det,B[4][k]))for k,v in enumerate(du)];Juv=[]
    for k in range(n):
        for l in range(k,n):
            d2=T.sub(T.add(T.mul(B[0][k],B[3][l]),T.mul(B[0][l],B[3][k])),T.add(T.mul(B[1][k],B[2][l]),T.mul(B[1][l],B[2][k])))
            Juv.append(T.add(T.mul(d2,h),T.add(T.mul(du[k],B[4][l]),T.mul(du[l],B[4][k]))))
    return J,Ju,Juv

def surface_derivatives(q,B):
    a,b,c=q;n=len(B[0]);v=T.add(T.mul(b,b),T.mul(c,c));a2=T.mul(a,a);G=T.mul(a2,v);du=[];d2=[]
    for k in range(n):
        vk=T.scale(T.add(T.mul(b,B[1][k]),T.mul(c,B[2][k])),2)
        ak=T.scale(T.mul(a,B[0][k]),2)
        du.append(T.add(T.mul(ak,v),T.mul(a2,vk)))
    for k in range(n):
        for l in range(k,n):
            ak=T.scale(T.mul(a,B[0][k]),2);al=T.scale(T.mul(a,B[0][l]),2)
            vk=T.scale(T.add(T.mul(b,B[1][k]),T.mul(c,B[2][k])),2);vl=T.scale(T.add(T.mul(b,B[1][l]),T.mul(c,B[2][l])),2)
            vv=T.scale(T.add(T.mul(B[1][k],B[1][l]),T.mul(B[2][k],B[2][l])),2)
            aa=T.scale(T.mul(B[0][k],B[0][l]),2)
            d2.append(T.add(T.add(T.mul(aa,v),T.mul(ak,vl)),T.add(T.mul(al,vk),T.mul(a2,vv))))
    return G,du,d2

def radial_cancellation(J,f):
    z=T.sub(J,T.ONE);co=f['corners'];rc=(min(x[0]for x in co)+max(x[0]for x in co))/2;hr=(max(x[0]for x in co)-min(x[0]for x in co))/2;hz=max(abs(x[1]-f['z_center'])for x in co);rlo=min(x[0]for x in co);rhi=max(x[0]for x in co)
    # Physical-R Laurent coefficient enclosure on the entire containing box.
    physical=sum(abs(c)*(rhi**i if i>=0 else rlo**i)*hz**j for (i,j),c in z.items())
    # J-1=P(R,zeta)/R; shift only its polynomial numerator before series powers.
    P=T.mul(z,T.X)
    if any(i<0 for i,j in P):raise ValueError('unexpected denominator beyond R')
    shifted=T.compose(P,{(0,0):rc,(1,0):hr},{(0,1):hz})
    local=norm(shifted)/rlo
    return {'Rcenter_exact':str(rc),'Rhalfwidth_exact':str(hr),'zeta_halfwidth_exact':str(hz),'physical_Laurent_absolute_sum_bound_exact':str(physical),'local_numerator_L1_over_Rmin_bound_exact':str(local),'coefficient_bound_reduction_factor':float(physical/local)if local else None,'max_numerator_degree':summarize(P)['max_total_degree'],'shifted_numerator_support':len(shifted),'order16_radial_backexpansion_cancellation_budget_bits':math.ceil(16*max(i for i,j in P)*math.log2(1+float(rc/hr))),'note':'Coefficient L1 bounds on a containing box, not attained strain ranges or certified series-tail acceptance. Negative R powers were NOT approximated by a finite shift.'}

def main():
    if (OUT/'results.json').exists():raise ValueError('write-once results exist')
    started=time.monotonic();src=ROOT/'scripts/exact_fe_target.py'
    if sha(src)!=EXPECTED:raise ValueError('frozen exact target source changed')
    protocol=T.read(ROOT/'runs/exact_integration_bridge_v1/target_preparation/protocol.json');row=protocol['records'][0];packet=T.read(ROOT/row['point_cache_file']);target=T.Target(packet)
    selected={}
    for i,ids in enumerate(target.a['triangles']):
        count=sum(target.a['reference'][ids[:3],0]==0);kind={0:'off_axis',1:'axis_vertex',2:'axis_edge'}[count]
        selected.setdefault(kind,i)
    selected['off_axis_outer']=max((i for i,ids in enumerate(target.a['triangles'])if all(target.a['reference'][ids[:3],0]>0)),key=lambda i:sum(target.a['reference'][target.a['triangles'][i,:3],0]))
    bulk=[]
    for kind,ti in selected.items():
        f=target.element(ti);J,Ju,Juv=jacobian_derivatives(f['q'],f['B']);z=T.sub(J,T.ONE);JuJu=[T.mul(a,b)for i,a in enumerate(Ju)for b in Ju[i:]]
        factors={k for p in Ju+Juv+JuJu for k in p};growth=support_growth(z,factors)
        bulk.append({'kind':kind,'element':ti,'free_local_dofs':len(f['free_local']),'q':[summarize(p)for p in f['q']],'B':groups([p for row in f['B']for p in row]),'J':summarize(J),'J_minus1':summarize(z),'J_u':groups(Ju),'J_uv_uppertriangle':groups(Juv),'J_u_times_J_v_uppertriangle':groups(JuJu),'series_support_upper_by_order':growth,'series16_naive_per_derivative_product_support_upper':growth[16]['with_factor_support_upper'],'radial_coordinate_cancellation':radial_cancellation(J,f)if kind.startswith('off_axis')else None})
    surfaces=[]
    for k in(0,1):
        for e in(0,1,11):
            f=target.surface(k,e);G,Gu,Guv=surface_derivatives(f['q_physical_R'],f['B_physical_R']);z=T.sub(G,T.ONE);GuGu=[T.mul(a,b)for i,a in enumerate(Gu)for b in Gu[i:]];factors={key for p in Gu+Guv+GuGu for key in p}
            surfaces.append({'surface':k,'edge':e,'Rlo_exact':str(f['Rlo']),'Rhi_exact':str(f['Rhi']),'free_local_dofs':len(f['free_local']),'q':[summarize(p)for p in f['q_physical_R']],'G_Js_squared':summarize(G),'G_u':groups(Gu),'G_uv_uppertriangle':groups(Guv),'G_u_times_G_v_uppertriangle':groups(GuGu),'series_support_upper_by_order':support_growth(z,factors)})
    inputs={'scripts/exact_fe_target.py':EXPECTED,'runs/exact_integration_bridge_v1/target_preparation/protocol.json':sha(ROOT/'runs/exact_integration_bridge_v1/target_preparation/protocol.json'),row['point_cache_file']:sha(ROOT/row['point_cache_file'])}
    result={'status':'PASS_LOCAL_ALGEBRA_DESIGN_PROBE','source_sha256':sha(__file__),'input_hashes':inputs,'record_id':row['record_id'],'bulk':bulk,'surfaces':surfaces,'checks':{'exact_target_source_unchanged':sha(src)==EXPECTED,'all3_bulk_chart_classes_probed':set(selected)>={'off_axis','axis_edge','axis_vertex'},'both_surface_families_first_interior_outer_probed':len(surfaces)==6,'no_actual_series_coefficients_above_power0_formed':True,'no_field_integration_or_solve':True},'runtime_seconds':time.monotonic()-started,'scope':'Actual B-center local polynomial/sparsity/support and coefficient-conditioning diagnostics. Minkowski supports are upper bounds, not actual nonzero power16 supports. No power16 coefficient expansion, physical range, validated nonlinear remainder, full field integration, equilibrium solve or exact-target root bridge.'}
    (OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print(json.dumps({'status':result['status'],'bulk':[(x['kind'],x['J']['support'],x['J_u']['support_sum'],x['J_uv_uppertriangle']['support_sum'],x['series16_naive_per_derivative_product_support_upper'])for x in bulk],'runtime_seconds':result['runtime_seconds'],'results_sha256':sha(OUT/'results.json')}),flush=True)
if __name__=='__main__':main()
