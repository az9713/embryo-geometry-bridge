"""Uniform local BULK exact-Qstar/cache Hessian consistency primitives.
Small coefficient differences are formed BEFORE evaluating the same nodal box.
No whole-field integral, 144-element cache sum, solve or root transfer is run.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations_with_replacement,product
import argparse,datetime,hashlib,json,time
import numpy as np
from flint import arb,arb_mat,ctx
import exact_fe_target as target_core
import algebraic_root_feasibility as core
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'runs/exact_integration_bridge_v1/cache_consistency_preparation'
PREP='runs/exact_integration_bridge_v1/target_preparation/protocol.json'
ORDER=target_core.ORDER

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf8'))
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    if p.exists():raise ValueError('write-once output exists: '+str(p))
    p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf8')
def A(v):
    if isinstance(v,arb):
        if not v.is_finite():raise ValueError('finite interval required')
        return v
    return core.A(target_core.rational(v))
def q(v):return target_core.rational(v)
def mat(rows):return arb_mat(len(rows),len(rows[0]),[A(v)for row in rows for v in row])
def maxabs(values):return max((core.upper(v)for v in values),default=arb(0))
def pack_matrix(value):return [[core.pack(value[i,j])for j in range(value.ncols())]for i in range(value.nrows())]

def jacobian_tensors(strain):
    if len(strain)!=5:raise ValueError('five bulk strain coordinates required')
    a,b,c,d,h=map(A,strain);D=a*d-b*c;J=D*h
    if not D>0 or not h>0 or not J>A(F(1,1000000)):raise ValueError('whole strain box must stay in positive-J analytic/no-penalty region')
    firstD=[d,-c,-b,a];secondD=lambda i,j:1 if(i,j)in((0,3),(3,0))else-1 if(i,j)in((1,2),(2,1))else 0
    first=[v*h for v in firstD]+[D];second=[[arb(0)for j in range(5)]for i in range(5)];third=[[[arb(0)for k in range(5)]for j in range(5)]for i in range(5)]
    for i in range(5):
        for j in range(5):
            if i<4 and j<4:second[i][j]=A(secondD(i,j))*h
            elif i==4 and j<4:second[i][j]=firstD[j]
            elif j==4 and i<4:second[i][j]=firstD[i]
            for k in range(5):
                indices=(i,j,k)
                if indices.count(4)==1:
                    remaining=[v for v in indices if v!=4];third[i][j][k]=A(secondD(*remaining))
    return J,first,second,third

def bulk_third_tensor(strain,mu,kappa):
    """Third derivative of W on a convex strain box; active quadratic vanishes.
    Scalar phi_third=(2*kappa*logJ-3*kappa-2*mu)/J**3.
    """
    mu,kappa=A(mu),A(kappa)
    if not mu>0 or not kappa>=0:raise ValueError('positive mu/nonnegative bulk modulus required')
    J,g,h,t=jacobian_tensors(strain);log=J.log();p1=(kappa*log-mu)/J;p2=(kappa+mu-kappa*log)/J**2;p3=(2*kappa*log-3*kappa-2*mu)/J**3
    cache={};tensor=[[[arb(0)for k in range(5)]for j in range(5)]for i in range(5)]
    for i in range(5):
        for j in range(5):
            for k in range(5):
                key=tuple(sorted((i,j,k)))
                if key not in cache:
                    u,v,w=key;cache[key]=p3*g[u]*g[v]*g[w]+p2*(h[u][v]*g[w]+h[u][w]*g[v]+h[v][w]*g[u])+p1*t[u][v][w]
                tensor[i][j][k]=cache[key]
    return tensor,{'J':J,'phi_prime':p1,'phi_second':p2,'phi_third':p3}

def consistency_bound(Bstar_full,Bcache_full,same_nodes,selected_local,mu,kappa,drive):
    """Uniform UNWEIGHTED local-free Hessian difference enclosure.
    theta is shared on both sides; delta_q=delta_B_full times SAME Y.
    """
    if len(Bstar_full)!=5 or len(Bcache_full)!=5 or not same_nodes:raise ValueError('nonempty full-coordinate five-row maps required')
    nfull=len(same_nodes)
    if any(len(row)!=nfull for row in Bstar_full+Bcache_full):raise ValueError('full maps must share the same nodal vector')
    selected_local=list(selected_local)
    if not selected_local or len(set(selected_local))!=len(selected_local)or any(isinstance(i,bool)or not isinstance(i,(int,np.integer))or not 0<=i<nfull for i in selected_local):raise ValueError('unique valid local free selection required')
    Bs=[[q(v)for v in row]for row in Bstar_full];Bc=[[q(v)for v in row]for row in Bcache_full];dB=[[s-c for s,c in zip(sr,cr)]for sr,cr in zip(Bs,Bc)];Y=list(map(A,same_nodes))
    qc=[sum((A(v)*y for v,y in zip(row,Y)),arb(0))for row in Bc];dq=[sum((A(v)*y for v,y in zip(row,Y)),arb(0))for row in dB]
    line=[v+arb(0).union(w)for v,w in zip(qc,dq)]
    tensor,domain=bulk_third_tensor(line,mu,kappa)
    mu,kappa,drive=map(A,(mu,kappa,drive))
    if not drive>=0:raise ValueError('nonnegative shared active drive required')
    _,hc,_=core.bulk_local(qc,mu,kappa,drive);Hc=mat(hc)
    dH=arb_mat(5,5,[sum((tensor[i][j][k]*dq[k]for k in range(5)),arb(0))for i in range(5)for j in range(5)])
    star=mat([[row[i]for i in selected_local]for row in Bs]);cache=mat([[row[i]for i in selected_local]for row in Bc]);delta=mat([[row[i]for i in selected_local]for row in dB])
    first=delta.transpose()*Hc*cache;second=star.transpose()*Hc*delta;third=star.transpose()*dH*star
    return {'difference':first+second+third,'term_deltaBT_Hc_Bc':first,'term_BstarT_Hc_deltaB':second,'term_BstarT_deltaH_Bstar':third,'local_deltaH':dH,'delta_q':dq,'cache_q':qc,'convex_line_box':line,'domain':domain,'delta_B_full':dB,'star_B_free':star,'cache_B_free':cache,'selected_local':selected_local,'scope':'All same-Y/shared-parameter inputs in supplied boxes; ONE local auxiliary point, unweighted. No full-cache/field/root transfer.'}

def exact_bulk_maps(target,element,quadrature_index):
    """Exact full12-column star/cache maps at an actual derived Qstar site."""
    a=target.a
    if not isinstance(element,int)or not 0<=element<len(a['triangles'])or not isinstance(quadrature_index,int)or not 0<=quadrature_index<len(a['shape']):raise ValueError('valid exact integer element/rule index required')
    rule=target.auxiliary_rule();lam=rule['volume_barycentric'][quadrature_index];ids=a['triangles'][element];g=target_core.material_basis(a['reference'][ids[:3]])
    R=sum(l*x[0]for l,x in zip(lam,g['corners']));zeta=sum(l*x[1]for l,x in zip(lam,g['corners']))-g['z_center'];oldR=q(a['qr'][element,quadrature_index])
    if R<=0 or oldR<=0:raise ValueError('positive interior reference radii required')
    star=[[F(0)for _ in range(12)]for _ in range(5)];cache=[[F(0)for _ in range(12)]for _ in range(5)]
    for j,p in enumerate(g['phi']):
        gr=[target_core.evaluate(target_core.diff(p,k),R,zeta)for k in(0,1)];old=[q(v)for v in a['shape_grad'][element,quadrature_index,j]]
        for B,v in[(star,gr),(cache,old)]:B[0][2*j]=v[0];B[1][2*j]=v[1];B[2][2*j+1]=v[0];B[3][2*j+1]=v[1]
        star[4][2*j]=target_core.evaluate(p,R,zeta)/R;cache[4][2*j]=q(a['shape'][quadrature_index,j])/oldR
    selected=[k for k in range(12)if 2*int(ids[k//2])+k%2 in target.index]
    return {'star_full':star,'cache_full':cache,'nodes':list(map(int,ids)),'selected_local':selected,'global_indices':[target.index[2*int(ids[k//2])+k%2]for k in selected],'Qstar_physical_exact':[R,zeta+g['z_center']],'auxiliary_normalized_weight_exact':rule['volume_weights'][element][quadrature_index]}


def synthetic_tests():
    checks={}
    def ck(k,v):
        checks[k]=bool(v)
        if not v:raise AssertionError(k)
    import sympy as sy
    a,b,c,d,h,mu,kappa,drive=sy.symbols('a b c d h mu kappa drive');xs=[a,b,c,d,h];D=a*d-b*c;J=D*h;ell=sy.log(J)
    js=[sy.diff(J,x)for x in xs];j2=[[sy.diff(J,x,y)for y in xs]for x in xs];j3=[[[sy.diff(J,x,y,z)for z in xs]for y in xs]for x in xs]
    p1=(kappa*ell-mu)/J;p2=(kappa+mu-kappa*ell)/J**2;p3=(2*kappa*ell-3*kappa-2*mu)/J**3
    scalar=sy.symbols('j',positive=True);phi=-mu*sy.log(scalar)+kappa*sy.log(scalar)**2/2
    for n,pn in[(1,p1),(2,p2),(3,p3)]:ck('scalar_phi_derivative_'+str(n),sy.simplify(sy.diff(phi,scalar,n).subs(scalar,J)-pn)==0)
    W=mu/2*(sum(x*x for x in xs)-3)-mu*ell+kappa*ell**2/2+drive/2*(b*b+d*d)
    for i,j,k in combinations_with_replacement(range(5),3):
        expected=p3*js[i]*js[j]*js[k]+p2*(j2[i][j]*js[k]+j2[i][k]*js[j]+j2[j][k]*js[i])+p1*j3[i][j][k]
        ck(f'Sympy_full_W_third_{i}{j}{k}',sy.cancel(sy.diff(W,xs[i],xs[j],xs[k])-expected)==0)
    point=[F(11,10),F(1,20),F(-1,25),F(9,10),F(21,20)];jp,jg,jh,jt=jacobian_tensors(point);subs=dict(zip(xs,map(lambda v:sy.Rational(v.numerator,v.denominator),point)))
    for i in range(5):
        ck(f'J_first_{i}',jg[i].contains(A(F(str(js[i].subs(subs))))))
        for j in range(5):
            ck(f'J_second_{i}{j}',jh[i][j].contains(A(F(str(j2[i][j].subs(subs))))))
            for k in range(5):ck(f'J_third_{i}{j}{k}',jt[i][j][k].contains(A(F(str(j3[i][j][k].subs(subs))))))
    tensor,_=bulk_third_tensor(point,F(1),F(7));ck('third_tensor_permutation_envelopes_identical',all((tensor[i][j][k].lower().fmpq(),tensor[i][j][k].upper().fmpq())==(tensor[k][j][i].lower().fmpq(),tensor[k][j][i].upper().fmpq())==(tensor[j][i][k].lower().fmpq(),tensor[j][i][k].upper().fmpq())for i in range(5)for j in range(5)for k in range(5)))
    for name,strain in [('negative_det',[F(-1),0,0,1,1]),('zero_hoop',[1,0,0,1,0]),('penalty_boundary',[1,0,0,1,F(1,1000000)]),('line_crosses_det_zero',[A(1).union(A(-1)),0,0,1,1])]:
        try:bulk_third_tensor(strain,1,7);ck('analytic_guard_'+name,False)
        except ValueError:checks['analytic_guard_'+name]=True
    I=[[F(int(i==j))for j in range(5)]for i in range(5)];zero=consistency_bound(I,I,point,list(range(5)),1,7,A(0).union(A(2)))
    ck('identical_coefficient_maps_zero_deltaq_and_global_bound',all(v==0 for v in zero['delta_q'])and all(v==0 for v in zero['difference'].entries()))
    eps=F(1,2**40);star=[[v+eps*F((i+2*j)%5-2)for j,v in enumerate(row)]for i,row in enumerate(I)];fixtures=[]
    for number,(base,kap)in enumerate([(point,F(7)),([F(6,5),F(-1,10),F(1,20),F(9,10),F(19,20)],F(5000))]):
        radius=F(1,1000);box=[A(v-radius).union(A(v+radius))for v in base];bound=consistency_bound(star,I,box,list(range(5)),1,kap,A(0).union(A(2)));Dout=bound['difference'];witnesses=0
        for signs in product((-1,1),repeat=5):
            yy=[A(v+s*radius)for v,s in zip(base,signs)]
            for active in(0,2):
                qc=yy;qs=[sum((A(v)*y for v,y in zip(row,yy)),arb(0))for row in star];_,hc,_=core.bulk_local(qc,A(1),A(kap),A(active));_,hs,_=core.bulk_local(qs,A(1),A(kap),A(active));sm=mat(star);direct=sm.transpose()*mat(hs)*sm-mat(hc)
                ck(f'interval_fixture{number}_vertex{witnesses}_enclosed',all(Dout[i,j].contains(direct[i,j])for i in range(5)for j in range(5)));witnesses+=1
        qcbox=box;qsbox=[sum((A(v)*y for v,y in zip(row,box)),arb(0))for row in star];naive=[s-c for s,c in zip(qsbox,qcbox)]
        ck('coefficientwise_deltaq_small_fixture'+str(number),maxabs(bound['delta_q'])<A(F(1,10**8))and maxabs(naive)>A(F(1,1000)))
        fixtures.append({'fixture':number,'bulk_modulus_exact':str(kap),'nodal_halfwidth_exact':str(radius),'delta_B_scale_exact':str(eps),'interval_parameter_drive_exact':['0','2'],'endpoint_witnesses':witnesses,'signed_decomposition_max_upper':core.pack(maxabs(Dout.entries())),'coefficientwise_deltaq_max_upper':core.pack(maxabs(bound['delta_q'])),'naive_deltaq_max_upper':core.pack(maxabs(naive))})
    return checks,fixtures

def select_B():
    prep=read(ROOT/PREP);row=next(r for r in prep['records']if r['record_id']=='optoGEF2_small_col_B');root=read(ROOT/row['root_certificate_file']);bindings={}
    for collection in(prep['source_hashes'],prep['input_hashes']):
        for path,digest in collection.items():
            target_core.bind(bindings,path,digest)
            if sha(ROOT/path)!=digest:raise ValueError('changed target input '+path)
    target_core.bind(bindings,PREP,sha(ROOT/PREP));target_core.bind(bindings,'scripts/exact_cache_consistency.py',sha(__file__))
    if root['artifact_hashes'][row['point_cache_file']]!=row['point_cache_sha256']or root['artifact_hashes'][row['predictor_file']]!=row['predictor_sha256']:raise ValueError('authoritative old cache/predictor binding required')
    target=target_core.Target(read(ROOT/row['point_cache_file']))
    with np.load(ROOT/row['predictor_file'],allow_pickle=False)as z:S=z['S'].copy()
    return prep,row,target,S,bindings

def point_fixture(target,S,row,element,iq):
    maps=exact_bulk_maps(target,element,iq);tube=target.node_tube(S,row['parameter_cell'],row['base_parameters_exact'],F(row['selected_corrector_radius_exact']));Y=[tube[node][c]for node in maps['nodes']for c in(0,1)]
    dose=F(row['dose_exact']);glo,ghi=map(F,row['parameter_cell']['common_gain']);active=q(target.a['active'][element,0]);drive=A(glo*dose*active).union(A(ghi*dose*active))
    bound=consistency_bound(maps['star_full'],maps['cache_full'],Y,maps['selected_local'],1,5000,drive)
    return {'status':'CERTIFIED_UNIFORM_LOCAL_CACHE_CONSISTENCY_UNREVIEWED','element':element,'quadrature_index':iq,'record_id':row['record_id'],'parameter_cell_sha256':row['parameter_cell_sha256'],'branch_domain_id':row['branch_domain_id'],'same_nodal_tube':'exact selected proof center plus affine predictor hull plus full normal corrector cube; fixed reference coordinates unchanged','global_indices':maps['global_indices'],'delta_B_max_abs_exact':str(max(abs(x)for row in bound['delta_B_full']for x in row)),'delta_q':[core.pack(v)for v in bound['delta_q']],'convex_line_J':core.pack(bound['domain']['J']),'local_deltaH_max_upper':core.pack(maxabs(bound['local_deltaH'].entries())),'unweighted_local_free_H_difference':pack_matrix(bound['difference']),'signed_term_maxima':{k:core.pack(maxabs(bound[k].entries()))for k in('term_deltaBT_Hc_Bc','term_BstarT_Hc_deltaB','term_BstarT_deltaH_Bstar')},'auxiliary_normalized_weight_exact_NOT_applied':str(maps['auxiliary_normalized_weight_exact']),'scope':'Only ONE auxiliary point, uniform over this supplied old-root nodal cube and shared drive interval. No weights applied, no quadrature sum, whole-field H bound, surface bound or exact-integral root transfer.'}

def prepare():
    if(OUT/'protocol.json').exists():raise ValueError('write-once preparation exists')
    ctx.prec=256;checks,toys=synthetic_tests();prep,row,target,S,bindings=select_B()
    protocol={'version':1,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_and_input_hashes':bindings,'selected_record':row,'precision_bits':256,'local_actual_fixture_points':[[0,0],[6,63]],'formula':'BstarT Hstar Bstar-BcT Hc Bc = deltaBT Hc Bc + BstarT Hc deltaB + BstarT deltaH Bstar. deltaq=deltaBfull*SAME Y. deltaH_ij=int_0^1 sum_k Wthird_ijk(qc+t*deltaq)*deltaq_k dt.','scalar_phi_third':'(2*kappa*log(J)-3*kappa-2*mu)/J^3','J':'(a*d-b*c)*h; exact first, second and third tensors.','analytic_guard':'Entire convex strain box: detFm>0, hoop>0 and J>1e-6, hence no trial penalty.','coefficient_types':'Canonical/cached B coefficients are exact rational binary/canonical values; interval subtraction of BstarY and BcacheY is forbidden for proof deltaq.','same_shared_parameter_policy':'mu/kappa/drive identical on both sides; all allowed values in their supplied boxes. Active drive quadratic contributes zero to third derivative.','scope':'BULK ONLY, two actual quadrature-site local fixtures. No whole144/full611 cache job, surface contribution, integration, solve or root transfer.'}
    write(OUT/'protocol.json',protocol);write(OUT/'synthetic_tests.json',{'status':'PASS','checks':checks,'interval_parameter_fixtures':toys,'source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json')})
    actual=[]
    for element,iq in protocol['local_actual_fixture_points']:
        result=point_fixture(target,S,row,element,iq);name=f'local_element{element}_q{iq}.json';write(OUT/name,result);actual.append({'file':name,'sha256':sha(OUT/name),'status':result['status']})
    for k,h in bindings.items():
        if sha(ROOT/k)!=h:raise ValueError('bound source/input changed during local fixtures '+k)
    preflight={'status':'PASS_LOCAL_PREPARATION_UNREVIEWED','source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json'),'synthetic_tests_sha256':sha(OUT/'synthetic_tests.json'),'synthetic_checks':len(checks),'all_synthetic_checks_literal_true':all(v is True for v in checks.values()),'two_actual_local_uniform_fixtures':actual,'all_frozen_inputs_current':True,'full_cache_job_launched':False,'full_uniform_Hstar_minus_Hcache_proved':False,'surface_variant_implemented':False,'root_transfer_proved':False}
    write(OUT/'preflight.json',preflight);print(json.dumps({'status':preflight['status'],'source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json'),'preflight_sha256':sha(OUT/'preflight.json'),'synthetic_checks':len(checks),'actual_local_points':len(actual)}),flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('command',choices=['test','prepare']);args=ap.parse_args();ctx.prec=256
    if args.command=='prepare':prepare()
    else:
        c,f=synthetic_tests();print(json.dumps({'status':'PASS','checks':len(c),'fixtures':f}),flush=True)
