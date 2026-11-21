"""Independent four-element POINT replay. Never an assembled/uniform bridge."""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import argparse,datetime,hashlib,json,sys,time
import numpy as np
from flint import arb,ctx
HERE=Path(__file__).resolve().parent;CF=HERE.parents[2]
sys.path.insert(0,str(CF/'scripts'));sys.path.insert(0,str(CF/'runs/exact_integration_bridge_v1/bridge_design_probe'))
import exact_fe_target as T
import exact_fe_moments as M
import polynomial_jet as DUAL
import algebraic_root_feasibility as C
ROOTOUT=CF/'runs/exact_integration_bridge_v1/bulk_pilot'
EXPECTED_DRIVER='91a4f9a41b69c7b9e9ae5d5e0021087a7dfa4cf9e80625273bb4945d0ac3d0b3'
EXPECTED_PROTOCOL='746b260f5fbe7d049958681533165784e4847a15ca7d302ae9cdf5712cda95f0'
EXPECTED_RESULTS='4f8a821eab8d894119632240990f1ecafb7c6e7feb290a5768b726cf97ba8ec0'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def write(p,x):
    if Path(p).exists():raise ValueError('immutable output exists: '+str(p))
    Path(p).parent.mkdir(parents=True,exist_ok=True);Path(p).write_text(json.dumps(x,indent=2)+'\n')
def rational(x):return F.from_float(float(x))if isinstance(x,(float,np.floating))else F(x)
def product(items):
    p=T.ONE
    for x in items:p=T.mul(p,x)
    return p
def independent_jet(q,B):
    # Different construction: derivatives of two signed cubic monomials,
    # rather than differentiating a meridional determinant and multiplying h.
    terms=[(1,(0,3,4)),(-1,(1,2,4))];n=len(B[0]);value={};g=[{}for _ in range(n)];H={(i,j):{}for i in range(n)for j in range(i,n)}
    for sign,indices in terms:
        value=T.add(value,T.scale(product(q[k]for k in indices),sign))
        for i in range(n):
            for a in range(3):g[i]=T.add(g[i],T.scale(product([B[indices[a]][i]]+[q[indices[b]]for b in range(3)if b!=a]),sign))
        for i,j in H:
            for a in range(3):
                for b in range(3):
                    if a==b:continue
                    H[i,j]=T.add(H[i,j],T.scale(product([B[indices[a]][i],B[indices[b]][j]]+[q[indices[c]]for c in range(3)if c not in(a,b)]),sign))
    return {'value':value,'gradient':g,'hessian_upper':H,'local_dofs':n}
def coefficients(n):
    log=[F(0)]+[F((-1)**(k+1),k)for k in range(1,n+2)]
    first=[F(-(-1)**k)+5000*sum((log[j]*(-1)**(k-j)for j in range(k+1)),F(0))for k in range(n+2)]
    return first[:n+1],[(k+1)*first[k+1]for k in range(n+1)]
def series(e,n):
    first,second=coefficients(n);power=T.ONE;a={};b={};ee={k:C.A(v)for k,v in e.items()}
    for k in range(n+1):
        a=T.add(a,T.scale(power,C.A(first[k])));b=T.add(b,T.scale(power,C.A(second[k])))
        if k<n:power=T.mul(power,ee)
    return a,b
def majorant(p,f):
    if f['chart']['kind']!='off_axis':return sum(map(abs,p.values()),F(0))
    rlo=min(v[0]for v in f['corners']);rhi=max(v[0]for v in f['corners']);z=max(abs(v[1]-f['zcenter'])for v in f['corners'])
    return sum((abs(c)*max(rlo**i,rhi**i)*z**j for(i,j),c in p.items()),F(0))
class Functional:
    def __init__(self,target,f,index,rule):
        self.target=target;self.f=f;self.points=[];self.lambdas=rule['volume_barycentric'];verts=f['corners'];zero=[i for i,x in enumerate(verts)if x[0]==0];off=[i for i in range(3)if i not in zero]
        for lam,w in zip(self.lambdas,rule['volume_weights'][index]):
            if f['chart']['kind']=='axis_vertex':
                ss=lam[off[0]]+lam[off[1]];xy=(ss,lam[off[1]]/ss)
            elif f['chart']['kind']=='axis_edge':xy=(lam[off[0]],lam[zero[1]]/(1-lam[off[0]]))
            else:xy=(sum(l*v[0]for l,v in zip(lam,verts)),sum(l*v[1]for l,v in zip(lam,verts))-f['zcenter'])
            self.points.append((xy,w))
        self.balls=[(C.A(x),C.A(y),C.A(w))for(x,y),w in self.points]
        self.massI=self.integral(0,0);self.massQ=C.A(sum(w for _,w in self.points))
        if not self.massI>0 or not self.massQ>0:raise ValueError('positive two-mass measures required')
    @lru_cache(None)
    def integral(self,i,j):
        if self.f['chart']['kind']=='off_axis':return M.normalized_bulk_moment(self.f['corners'],i,j,self.target.energy_scale,self.f['zcenter'],640)
        return self.target.normalization()*C.A(M.axis_chart_moment(self.f['chart']['kind'],i,j,self.f['chart']['R1'],self.f['det']))
    @lru_cache(None)
    def node_power(self,k,d,e):return self.balls[k][d]**e
    @lru_cache(None)
    def defect(self,i,j):
        discrete=sum((w*self.node_power(k,0,i)*self.node_power(k,1,j)for k,(_,_,w)in enumerate(self.balls)),arb(0))
        return self.integral(i,j)-discrete
    def apply(self,p):return sum((T.cmul(c,self.defect(*e))for e,c in p.items()),arb(0))
def direct_basis(corners,lam):
    a=corners[1][0]-corners[0][0];b=corners[2][0]-corners[0][0];c=corners[1][1]-corners[0][1];d=corners[2][1]-corners[0][1];det=a*d-b*c
    dl=[((c-d)/det,(b-a)/det),(d/det,-b/det),(-c/det,a/det)]
    shape=[x*(2*x-1)for x in lam]+[4*lam[i]*lam[j]for i,j in[(0,1),(1,2),(2,0)]]
    grad=[[(4*lam[i]-1)*dl[i][k]for k in range(2)]for i in range(3)]+[[4*(lam[i]*dl[j][k]+lam[j]*dl[i][k])for k in range(2)]for i,j in[(0,1),(1,2),(2,0)]]
    return shape,grad,sum(x*y[0]for x,y in zip(lam,corners))
def correction(target,f,index,functional,drive):
    arrays=target.a;nodes=target.nodes();ys=[nodes[int(k)]for k in arrays['triangles'][index]];n=len(f['free_local']);g=[arb(0)for _ in range(n)];H={(i,j):arb(0)for i in range(n)for j in range(i,n)};basis_checks=0
    for iq,(lam,(xy,w))in enumerate(zip(functional.lambdas,functional.points)):
        N,D,rr=direct_basis(f['corners'],lam);star=[[F(0)for _ in range(12)]for _ in range(5)];old=[[F(0)for _ in range(12)]for _ in range(5)]
        for k in range(6):
            star[0][2*k]=D[k][0];star[1][2*k]=D[k][1];star[2][2*k+1]=D[k][0];star[3][2*k+1]=D[k][1];star[4][2*k]=N[k]/rr
            x=rational(arrays['shape_grad'][index,iq,k,0]);z=rational(arrays['shape_grad'][index,iq,k,1]);old[0][2*k]=x;old[1][2*k]=z;old[2][2*k+1]=x;old[3][2*k+1]=z;old[4][2*k]=rational(arrays['shape'][iq,k])/rational(arrays['qr'][index,iq])
        flat=[v for pair in ys for v in pair];qs=[sum(x*y for x,y in zip(row,flat))for row in star]
        if qs!=[T.evaluate(p,*xy)for p in f['q']]:raise ValueError('independent canonical-point field mismatch')
        for k in range(5):
            for i,local in enumerate(f['free_local']):
                if star[k][local]!=T.evaluate(f['B'][k][i],*xy):raise ValueError('independent canonical-point free basis mismatch')
                basis_checks+=1
        qc=[sum(x*y for x,y in zip(row,flat))for row in old]
        gs,hs,_=C.bulk_local(list(map(C.A,qs)),C.A(1),C.A(5000),C.A(drive));gc,hc,_=C.bulk_local(list(map(C.A,qc)),C.A(1),C.A(5000),C.A(drive))
        for i,li in enumerate(f['free_local']):
            g[i]+=C.A(w)*sum((C.A(star[k][li])*gs[k]-C.A(old[k][li])*gc[k]for k in range(5)),arb(0))
            for j in range(i,n):
                lj=f['free_local'][j];H[i,j]+=C.A(w)*sum((C.A(star[k][li])*hs[k][l]*C.A(star[l][lj])-C.A(old[k][li])*hc[k][l]*C.A(old[l][lj])for k in range(5)for l in range(5)),arb(0))
    return g,H,basis_checks
def bounds(x):return F(str(x.lower().fmpq())),F(str(x.upper().fmpq()))
def overlaps(x,p):
    lo,hi=bounds(x);return lo<=F(p['upper_exact'])and F(p['lower_exact'])<=hi
def selftests():
    checks={};q=[T.constant(F(2)),T.constant(F(1,3)),T.constant(F(-1,7)),T.constant(F(3,2)),T.constant(F(6,5))];B=[[T.constant(F(i+1,j+3))for j in range(2)]for i in range(5)]
    ours=independent_jet(q,B);reference=DUAL.bulk_jacobian_jet(q,B)
    checks['cubic_monomial_jet_matches_reviewed_determinant_jet']=ours==reference
    checks['convolution_coefficients_order2']=coefficients(2)==([F(-1),F(5001),F(-7501)],[F(5001),F(-15002),F(27503)])
    corners=[[F(1),F(0)],[F(2),F(0)],[F(1),F(1)]];N,D,r=direct_basis(corners,[F(1,3)]*3)
    checks['direct_canonical_basis_partition']=sum(N)==1 and all(sum(x[k]for x in D)==0 for k in range(2))and r==F(4,3)
    checks['component_source_current']=sha(CF/'scripts/exact_fe_moments.py')=='52bd1e4b54f4ae03ca66592bf4ab295fef30713e5fe171e6fc3e955acf2b9925'
    if not all(checks.values()):raise ValueError(checks)
    return checks
def prepare():
    ctx.prec=640;C.A.cache_clear();tests=selftests()
    if sha(CF/'scripts/exact_bulk_bridge_pilot.py')!=EXPECTED_DRIVER or sha(ROOTOUT/'protocol_v2.json')!=EXPECTED_PROTOCOL:raise ValueError('producer source/protocol changed')
    if sha(ROOTOUT/'results.json')!=EXPECTED_RESULTS:raise ValueError('authoritative producer terminal receipt differs')
    terminal=read(ROOTOUT/'results.json')
    if terminal['status']!='BOUNDED_ELEMENT_PILOT_COMPLETE_UNREVIEWED'or len(terminal['elements'])!=4:raise ValueError('producer not fully terminal')
    producer=read(ROOTOUT/'protocol_v2.json');inputs={}
    for collection in['source_and_input_hashes','target_preparation_input_hashes']:
        for name,digest in producer[collection].items():
            path=CF/name
            if sha(path)!=digest:raise ValueError('producer input changed')
            canonical=path.resolve().relative_to(CF.resolve()).as_posix()
            if canonical in inputs and inputs[canonical]!=digest:raise ValueError('conflicting binding')
            inputs[canonical]=digest
    if {x['file']for x in terminal['elements']}!={f'element_{i:03d}.json'for i in[0,6,72,66]}:raise ValueError('wrong four producer element packets')
    for x in terminal['elements']:
        file=ROOTOUT/x['file']
        if sha(file)!=x['sha256']:raise ValueError('producer element differs')
        inputs[file.relative_to(CF).as_posix()]=x['sha256']
    for file in[Path(__file__),HERE/'METHOD.md',ROOTOUT/'results.json',ROOTOUT/'protocol_v2.json',CF/'scripts/exact_fe_moments.py',CF/'scripts/exact_fe_target.py',Path(DUAL.__file__),Path(C.__file__)]:inputs[file.resolve().relative_to(CF.resolve()).as_posix()]=sha(file)
    protocol={'source_sha256':sha(Path(__file__)),'input_hashes':inputs,'elements':[0,6,72,66],'order':18,'precision_bits':640,'rho_exact':'21/500','producer_driver_sha256':EXPECTED_DRIVER,'producer_protocol_sha256':EXPECTED_PROTOCOL,'scope':'Same four POINT bulk elements only; independently constructed cubic jet/convolution/canonical point maps; reused reviewed exactmoment/tail/constitutive/dual kernels.','selfchecks':tests,'actual_replay_started':False}
    write(HERE/'protocol.json',protocol);write(HERE/'selfchecks.json',{'status':'PASS','checks':tests,'protocol_sha256':sha(HERE/'protocol.json')});print(json.dumps({'status':'PREPARED_ONLY','source_sha256':sha(Path(__file__)),'protocol_sha256':sha(HERE/'protocol.json'),'selfchecks':tests}),flush=True)
def run():
    protocol=read(HERE/'protocol.json')
    for name,d in protocol['input_hashes'].items():
        if sha(CF/name)!=d:raise ValueError('frozen replay input changed '+name)
    ctx.prec=640;C.A.cache_clear();selftests();write(HERE/'started.json',{'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'protocol_sha256':sha(HERE/'protocol.json')})
    prod=read(ROOTOUT/'protocol_v2.json');record=read(CF/'runs/cohort_uniform_roots_v2/records/00_optoGEF2_small_col_B/result.json');target=T.Target(read(CF/record['point_cache_file']));rule=target.auxiliary_rule();checks={};rows=[]
    for index in protocol['elements']:
        start=time.monotonic();f=target.element(index);jet=independent_jet(f['q'],f['B']);n=jet['local_dofs'];functional=Functional(target,f,index,rule);e=T.sub(jet['value'],T.ONE);a,b=series(e,18);vol=DUAL.contract_volumetric_defect(jet,a,b,functional.defect)
        checks[f'element{index}_direct_gradient_moment_vs_dual']=functional.apply(T.mul(a,jet['gradient'][0])).overlaps(vol['gradient'][0])
        direct=T.add(T.mul(a,jet['hessian_upper'][0,0]),T.mul(b,T.mul(jet['gradient'][0],jet['gradient'][0])))
        checks[f'element{index}_direct_H00_moment_vs_dual']=functional.apply(direct).overlaps(vol['hessian_upper'][0,0])
        drive=F(record['base_strength_exact'])*rational(target.a['active'][index,0]);diag=[F(1),1+drive,F(1),1+drive,F(1)]
        t1,t2=M.volumetric_derivative_tail_bounds(18,F(21,500),1,5000);mass=C.upper(functional.massI+functional.massQ);jg=[majorant(p,f)for p in jet['gradient']];gf=[];hf={}
        for i in range(n):
            p={}
            for k in range(5):p=T.add(p,T.scale(T.mul(f['q'][k],f['B'][k][i]),diag[k]))
            center=functional.apply(p)+vol['gradient'][i];tail=mass*C.A(t1*jg[i]);gf.append(center+(-tail).union(tail))
        for i,j in jet['hessian_upper']:
            p={}
            for k in range(5):p=T.add(p,T.scale(T.mul(f['B'][k][i],f['B'][k][j]),diag[k]))
            center=functional.apply(p)+vol['hessian_upper'][i,j];tail=mass*C.A(t2*jg[i]*jg[j]+t1*majorant(jet['hessian_upper'][i,j],f));hf[i,j]=center+(-tail).union(tail)
        cg,ch,bcount=correction(target,f,index,functional,drive);gf=[v+cg[i]for i,v in enumerate(gf)];hf={k:v+ch[k]for k,v in hf.items()};old=read(ROOTOUT/f'element_{index:03d}.json')
        checks[f'element{index}_same_chart_global_free_projection']=old['chart']==f['chart']['kind']and old['global_indices']==f['global_indices']and old['local_dofs']==n
        checks[f'element{index}_all_force_intervals_overlap']=len(old['force_integral_minus_cache'])==n and all(overlaps(v,p)for v,p in zip(gf,old['force_integral_minus_cache']))
        oldH={(x['i'],x['j']):x['enclosure']for x in old['H_integral_minus_cache_upper_triangle']}
        checks[f'element{index}_all_Hessian_intervals_overlap']=set(oldH)==set(hf)and all(overlaps(v,oldH[k])for k,v in hf.items())
        out={'element':index,'order':18,'precision_bits':640,'chart':f['chart']['kind'],'global_indices':f['global_indices'],'local_dofs':n,'canonical_basis_point_entry_checks':bcount,'independent_force':[C.pack(v)for v in gf],'independent_Hessian_upper':[{'i':i,'j':j,'interval':C.pack(v)}for(i,j),v in hf.items()],'seconds':time.monotonic()-start,'max_force_upper':C.pack(max(C.upper(v)for v in gf)),'max_Hessian_entry_upper':C.pack(max(C.upper(v)for v in hf.values()))}
        file=HERE/f'element_{index:03d}.json';write(file,out);rows.append({'file':file.name,'sha256':sha(file),'force_entries':n,'Hessian_entries':len(hf),'canonical_basis_point_entry_checks':bcount});print(json.dumps({'element':index,'seconds':out['seconds'],'checks_so_far':checks}),flush=True)
    checks['all_frozen_inputs_still_match']=all(sha(CF/k)==v for k,v in protocol['input_hashes'].items())
    result={'status':'PASS_INDEPENDENT_POINT_BULK_PILOT'if all(checks.values())else'FAIL_INDEPENDENT_POINT_BULK_PILOT','checks':checks,'protocol_sha256':sha(HERE/'protocol.json'),'review_source_sha256':sha(Path(__file__)),'elements':rows,'producer_driver_sha256':EXPECTED_DRIVER,'producer_protocol_sha256':EXPECTED_PROTOCOL,'precision_bits':640,'order':18,'reviewed_reuse':['ExactFEtarget polynomial algebra/field construction','Exactrational/logarithmicmoments and derivative-tail majorants','Constitutive five-variable gradient/Hessian','Signedmoment-dual contraction'],'independent_components':['Signed cubicmonomial Jacobianjet','Log/geometricseries coefficient convolution','Direct barycentric canonical strain/freebasismaps atallpoints','Manual different-association pointconstitutive lifts','Direct moment-product checks forgradient0/H00','Higherorder18 andprecision640 with intervaloverlap comparison'],'scope':'Four bulk elements at saved recordBcenter only. No surfaces/globalassembly/tubeuniformity/roottransfer/PDE/sourcevalidation.','root_transfer_proved':False}
    write(HERE/'results.json',result);write(CF/'evidence/exact_bulk_pilot_independent_review.json',result);print(json.dumps({'status':result['status'],'receipt_sha256':sha(CF/'evidence/exact_bulk_pilot_independent_review.json')}),flush=True)
    if not all(checks.values()):raise SystemExit(1)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=['prepare','run']);a=p.parse_args();prepare()if a.command=='prepare'else run()
