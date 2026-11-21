"""Independent interval-AD review of one frozen binary-quadrature root.

No equilibrium optimization. Rebuilds strain maps and second derivatives from
the embedded cache; does not call the producer assembler or row-bound routine.
"""
from pathlib import Path
from fractions import Fraction as Q
from functools import lru_cache
import argparse,base64,hashlib,io,json,os,time
import numpy as np
from flint import arb,arb_mat,fmpq,ctx

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'runs/algebraic_root_independent_review'
CERT=ROOT/'evidence/algebraic_root_certification_feasibility.json'
RECEIPT=ROOT/'evidence/algebraic_root_independent_review.json'
EXPECTED_PRODUCER='7321d38bf17cf354cc03be2871f3fbaa45bcebbcbfca7be90d0f867b80078acf'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text())
def save(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')

@lru_cache(maxsize=300000)
def A(x):
    q=Q.from_float(float(x)) if isinstance(x,(float,np.floating)) else Q(x)
    return arb(fmpq(q.numerator,q.denominator))
def UP(x):return x.abs_upper()
def pack(x):return {'ball':str(x),'lower_exact':str(x.lower().fmpq()),'upper_exact':str(x.upper().fmpq())}
def overlap(a,b):return bool(a.lower()<=b.upper() and b.lower()<=a.upper())
def contains(a,b):return bool(a.lower()<=b.lower() and b.upper()<=a.upper())

class Jet:
    """Second-order forward AD with Arb values, gradients and Hessians."""
    def __init__(self,value,g,H):self.v=value;self.g=g;self.H=H;self.n=len(g)
    @classmethod
    def constant(cls,value,n):return cls(value if isinstance(value,arb) else A(value),[arb(0)for _ in range(n)],[[arb(0)for _ in range(n)]for _ in range(n)])
    @classmethod
    def variable(cls,value,n,index):
        out=cls.constant(value,n);out.g[index]=arb(1);return out
    def coerce(self,x):return x if isinstance(x,Jet) else Jet.constant(x,self.n)
    def __add__(self,x):
        x=self.coerce(x);return Jet(self.v+x.v,[a+b for a,b in zip(self.g,x.g)],[[a+b for a,b in zip(aa,bb)]for aa,bb in zip(self.H,x.H)])
    __radd__=__add__
    def __neg__(self):return Jet(-self.v,[-x for x in self.g],[[-x for x in row]for row in self.H])
    def __sub__(self,x):return self+-self.coerce(x)
    def __rsub__(self,x):return self.coerce(x)+-self
    def __mul__(self,x):
        x=self.coerce(x);return Jet(self.v*x.v,[self.g[i]*x.v+self.v*x.g[i]for i in range(self.n)],
            [[self.H[i][j]*x.v+self.v*x.H[i][j]+self.g[i]*x.g[j]+x.g[i]*self.g[j]for j in range(self.n)]for i in range(self.n)])
    __rmul__=__mul__
    def unary(self,value,first,second):
        return Jet(value,[first*x for x in self.g],[[first*self.H[i][j]+second*self.g[i]*self.g[j]for j in range(self.n)]for i in range(self.n)])
    def inverse(self):return self.unary(1/self.v,-1/(self.v*self.v),2/(self.v*self.v*self.v))
    def __truediv__(self,x):return self*self.coerce(x).inverse()
    def log(self):
        if not self.v>0:raise ValueError('positive logarithm argument required')
        return self.unary(self.v.log(),1/self.v,-1/(self.v*self.v))
    def sqrt(self):
        if not self.v>0:raise ValueError('positive square-root argument required')
        value=self.v.sqrt();return self.unary(value,1/(2*value),-1/(4*value*value*value))

def bulk_AD(values,mu,kappa,drive):
    q=[Jet.variable(v,5,i)for i,v in enumerate(values)];a,b,c,d,h=q
    det=a*d-b*c;J=det*h
    if not det.v>0 or not h.v>0 or not J.v>A(Q(1,1000000)):raise ValueError('independent ball outside analytic no-penalty region')
    L=J.log();energy=mu/2*(sum(v*v for v in q)-3-2*L)+kappa/2*L*L+drive/2*(b*b+d*d)
    return energy.g,energy.H,J.v,{'bulk_det_min':det.v,'bulk_hoop_min':h.v,'bulk_J_min':J.v}

def surface_AD(values,K,T,scale):
    R,x,z=[Jet.variable(v,3,i)for i,v in enumerate(values)]
    length2=x*x+z*z
    if not R.v>0 or not length2.v>0:raise ValueError('independent surface ball degenerate')
    js=scale*R*length2.sqrt();energy=K/2*(js-1)*(js-1)+T*js
    return energy.g,energy.H,js.v,{'surface_radius_min':R.v,'surface_tangent_squared_min':length2.v,'surface_J_min':js.v}

def local_tests():
    import algebraic_root_feasibility as producer
    if sha(Path(producer.__file__))!=EXPECTED_PRODUCER:raise RuntimeError('producer source changed')
    ctx.prec=256;A.cache_clear();producer.A.cache_clear();checks={}
    bulk_points=[['1','1/5','1/10','11/10','9/10'],['9/10','-1/7','1/8','11/10','21/20'],['6/5','1/11','-1/9','9/10','19/20']]
    surface_points=[['2','3','-4'],['1/3','-2','1/5'],['7/5','1/10','-3/2']]
    negative=False
    for kind,points,fn,manual,params in [('bulk',bulk_points,bulk_AD,producer.bulk_local,[A(1),A(5000),A(Q(3,7))]),
            ('surface',surface_points,surface_AD,producer.surface_local,[A(Q(2,3)),A(Q(7,5)),A(Q(3,4))])]:
        for case,values in enumerate(points):
            center=list(map(A,values));g,H,J,domain=fn(center,*params);gm,Hm,Jm=manual(center,*params)
            checks[f'{kind}_{case}_point_gradient_all']=all(overlap(a,b)for a,b in zip(g,gm))
            checks[f'{kind}_{case}_point_Hessian_all']=all(overlap(a,b)for aa,bb in zip(H,Hm)for a,b in zip(aa,bb))
            checks[f'{kind}_{case}_point_J']=overlap(J,Jm)
            negative=negative or any(H[i][j].upper()<0 for i in range(len(H))for j in range(len(H))if i!=j)
            radius=A(Q(1,10**8));box=[v+arb(0,radius.upper())for v in center]
            bg,bH,bJ,_=fn(box,*params);mg,mH,mJ=manual(box,*params)
            for sign in (-1,1):
                point=[v+sign*radius/A(i+2)for i,v in enumerate(center)];pg,pH,pJ,_=fn(point,*params)
                checks[f'{kind}_{case}_ball_point{sign}_AD_gradient']=all(contains(a,b)for a,b in zip(bg,pg))
                checks[f'{kind}_{case}_ball_point{sign}_AD_Hessian']=all(contains(a,b)for aa,bb in zip(bH,pH)for a,b in zip(aa,bb))
                checks[f'{kind}_{case}_ball_point{sign}_producer_Hessian']=all(contains(a,b)for aa,bb in zip(mH,pH)for a,b in zip(aa,bb))
            checks[f'{kind}_{case}_ball_Hessian_overlap']=all(overlap(a,b)for aa,bb in zip(bH,mH)for a,b in zip(aa,bb))
    checks['negative_offdiagonal_coupling_exercised']=bool(negative)
    for key,call in [('negative_J',lambda:bulk_AD([A(1),A(0),A(0),A(-1),A(1)],A(1),A(5),A(1))),
                     ('zero_surface_tangent',lambda:surface_AD([A(1),A(0),A(0)],A(1),A(1),A(1)))]:
        try:call();checks[key+'_guard']=False
        except ValueError:checks[key+'_guard']=True
    save(OUT/'synthetic.json',{'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'source_sha256':sha(__file__),
        'producer_source_sha256':EXPECTED_PRODUCER,'scope':'Separate interval forward-AD derivative implementation on signed/off-diagonal point and ball cases; not yet actual root replay'})
    print(json.dumps({'stage':'independent_local_AD','checks':len(checks),'pass':all(checks.values())}),flush=True)
    if not all(checks.values()):raise RuntimeError('independent local AD checks failed')

def terms_from_arrays(arrays,parameters,strength):
    """Independent exact linear strain maps including all fixed-node terms."""
    free=np.asarray(arrays['free']);lookup={int(x):i for i,x in enumerate(free)}
    coords=[A(float(x))for x in arrays['reference'].ravel()]
    for i,index in enumerate(free):coords[int(index)]=coords[int(index)]+A(float(arrays['u'][i]))
    denominator=A(float(arrays['energy_scale'][0]))
    if not denominator>0:raise ValueError('normalization must be positive')
    for triangle,nodes in enumerate(arrays['triangles']):
        dofs=[2*int(node)+coordinate for node in nodes for coordinate in (0,1)]
        selected=[(j,lookup[dof])for j,dof in enumerate(dofs)if dof in lookup]
        for point in range(len(arrays['shape'])):
            radius=A(float(arrays['qr'][triangle,point]))
            if not radius>0:raise ValueError('positive frozen quadrature radius required')
            M=[[arb(0)for _ in dofs]for _ in range(5)]
            for node in range(6):
                dr=A(float(arrays['shape_grad'][triangle,point,node,0]));dz=A(float(arrays['shape_grad'][triangle,point,node,1]))
                M[0][2*node]=dr;M[1][2*node]=dz;M[2][2*node+1]=dr;M[3][2*node+1]=dz
                M[4][2*node]=A(float(arrays['shape'][point,node]))/radius
            q=[sum((v*coords[index]for v,index in zip(row,dofs)),arb(0))for row in M]
            active=arrays['active'][triangle,0 if arrays['active'].shape[1]==1 else point]
            yield {'kind':'bulk','q':q,'B':[[row[j]for j,_ in selected]for row in M],
                'indices':[i for _,i in selected],'weight':A(float(arrays['volume_weights'][triangle,point]))/denominator,
                'parameters':(A(parameters['mu']),A(parameters['bulk']),A(float(active))*A(strength))}
    for surface,K,T in [(0,parameters['basal_area_modulus'],parameters['basal_surface_tension']),
                        (1,parameters['apical_area_modulus'],parameters['apical_surface_tension'])]:
        key=lambda name:arrays[f'surface_{surface}_{name}']
        for edge,nodes in enumerate(key('edges')):
            dofs=[2*int(node)+coordinate for node in nodes for coordinate in (0,1)]
            selected=[(j,lookup[dof])for j,dof in enumerate(dofs)if dof in lookup]
            for point in range(len(key('N'))):
                radius=A(float(key('r0')[edge,point]));length=A(float(key('l0')[edge]))
                if not radius>0 or not length>0:raise ValueError('positive frozen surface reference required')
                M=[[arb(0)for _ in dofs]for _ in range(3)]
                for node in range(3):
                    M[0][2*node]=A(float(key('N')[point,node]));M[1][2*node]=A(float(key('dN')[point,node]));M[2][2*node+1]=M[1][2*node]
                yield {'kind':'surface','q':[sum((v*coords[index]for v,index in zip(row,dofs)),arb(0))for row in M],
                    'B':[[row[j]for j,_ in selected]for row in M],'indices':[i for _,i in selected],
                    'weight':A(float(key('weights')[edge,point]))/denominator,'parameters':(A(K),A(T),1/(radius*length))}

def assemble(terms,n):
    H=arb_mat(n,n);g=[arb(0)for _ in range(n)];cache=[];counts={'bulk':0,'surface':0}
    for index,term in enumerate(terms):
        if not term['weight']>0:raise ValueError('positive quadrature weights required')
        fn=bulk_AD if term['kind']=='bulk' else surface_AD;gradient,hessian,J,domain=fn(term['q'],*term['parameters'])
        counts[term['kind']]+=1;dim=len(term['indices'])
        if not dim:continue
        B=arb_mat(len(term['B']),dim,[x for row in term['B']for x in row]);HH=arb_mat(len(hessian),len(hessian),[x for row in hessian for x in row]);gg=arb_mat(len(gradient),1,gradient)
        localH=B.transpose()*HH*B;localg=B.transpose()*gg
        for a,i in enumerate(term['indices']):
            g[i]+=term['weight']*localg[a,0]
            for b,j in enumerate(term['indices']):H[i,j]+=term['weight']*localH[a,b]
        term['Hcenter']=hessian;term['absB']=[[UP(x)for x in row]for row in term['B']]
        term['rowsumB']=[sum(row,arb(0)).upper()for row in term['absB']];cache.append(term)
        if index%1536==0:print(json.dumps({'phase':'independent_AD_point_assembly','terms':index}),flush=True)
    return arb_mat(n,1,g),H,cache,counts

def variation(cache,n,radius):
    rows=[arb(0)for _ in range(n)];margins={};counts={'bulk':0,'surface':0}
    for index,term in enumerate(cache):
        q=[center+arb(0,(radius*coefficient).upper())for center,coefficient in zip(term['q'],term['rowsumB'])]
        fn=bulk_AD if term['kind']=='bulk' else surface_AD;gradient,H,J,domain=fn(q,*term['parameters'])
        counts[term['kind']]+=1;domain['normalized_weight_min']=term['weight']
        for name,v in domain.items():
            low=v.lower();margins[name]=low if name not in margins else min(margins[name],low)
        delta=[[UP(a-b)for a,b in zip(aa,bb)]for aa,bb in zip(H,term['Hcenter'])]
        projected=[sum((d*s for d,s in zip(row,term['rowsumB'])),arb(0))for row in delta]
        for j,i in enumerate(term['indices']):rows[i]+=term['weight']*sum((term['absB'][a][j]*projected[a]for a in range(len(q))),arb(0))
        if index%1536==0:print(json.dumps({'phase':'independent_AD_ball_assembly','terms':index}),flush=True)
    return arb_mat(n,1,[row.upper()for row in rows]),{'counts':counts,'strict_lower_margins':{k:pack(v)for k,v in margins.items()}}

def verify(precision):
    import torch
    from sparse_quadratic import AxisymmetricModel,Parameters
    ctx.prec=precision;A.cache_clear();torch.set_num_threads(1)
    test=read(OUT/'synthetic.json')
    if test['status']!='PASS' or test['source_sha256']!=sha(__file__):raise RuntimeError('matching independent synthetic checks required')
    producer_path=ROOT/'scripts/algebraic_root_feasibility.py'
    if sha(producer_path)!=EXPECTED_PRODUCER:raise RuntimeError('producer source changed')
    if (OUT/'protocol.json').exists():raise RuntimeError('independent replay already recorded; do not duplicate')
    cert_bytes=CERT.read_bytes();certificate=json.loads(cert_bytes);certificate_hash=hashlib.sha256(cert_bytes).hexdigest()
    if certificate['status']!='CERTIFIED_LOCAL_ROOT_BINARY_QUADRATURE':raise RuntimeError('producer certificate not successful')
    (OUT/'certificate_snapshot.json').write_bytes(cert_bytes)
    checks={'producer_source_bound':certificate['script_sha256']==EXPECTED_PRODUCER,
        'all_physical_sources_unchanged':all(sha(ROOT/'scripts'/name)==value for name,value in certificate['source_hashes'].items()),
        'actual_state_hash':sha(ROOT/certificate['input_state'])==certificate['input_state_sha256'],
        'matrix_archive_hash':sha(ROOT/certificate['matrix_archive'])==certificate['matrix_archive_sha256']}
    raw=base64.b64decode(certificate['frozen_binary_cache']['data']);checks['cache_payload_checksum']=hashlib.sha256(raw).hexdigest()==certificate['frozen_binary_cache']['sha256']
    with np.load(io.BytesIO(raw),allow_pickle=False)as z:arrays={k:z[k].copy()for k in z.files}
    n=len(arrays['u']);checks['actual611_DOF']=n==certificate['free_dof']==611
    checks['finite_preconditioner_square']=arrays['preconditioner'].shape==(n,n) and bool(np.all(np.isfinite(arrays['preconditioner'])))
    state=read(ROOT/certificate['input_state']);model=AxisymmetricModel(Parameters(**certificate['parameters']))
    live={'reference':model.reference,'free':model.free,'triangles':model.triangles,'shape':model.shape.numpy(),
        'shape_grad':model.shape_grad.numpy(),'qr':model.qr.numpy(),'volume_weights':model.integration.numpy(),
        'active':model.active.numpy(),'energy_scale':np.array([model.energy_scale]),'u':np.asarray(state['free_displacement'])}
    for k,surface in enumerate(model.surfaces):
        for name,value in zip(('edges','N','dN','r0','l0','weights'),surface[:6]):live[f'surface_{k}_{name}']=value.numpy()
    checks['all_frozen_cache_arrays_match_actual_model']=all(np.array_equal(value,arrays[key])for key,value in live.items())
    with np.load(ROOT/certificate['matrix_archive'],allow_pickle=False)as z:
        checks['archive_vector_and_free_dof_match']=np.array_equal(z['free_vector'],arrays['u']) and np.array_equal(z['free_indices'],arrays['free'])
    checks['state_parameters_strength_match']=state['parameters']==certificate['parameters'] and Q.from_float(state['activation_strength'])==Q(certificate['activation_strength_binary_exact'])
    successful=next(row for row in certificate['radius_attempts']if row['strict_selfmap_and_contraction']);rq=Q(successful['radius_requested_exact']);radius=A(rq)
    eta0=Q(certificate['preconditioned_residual_eta']['upper_exact']);Zold=Q(successful['Z_total']['upper_exact'])
    checks['producer_exact_rational_Banach_inequalities']=Zold<1 and eta0+Zold*rq<rq
    checks['producer_error_bound_covers_exact_ratio']=Q(successful['root_distance_upper']['upper_exact'])>=eta0/(1-Zold)
    protocol={'source_sha256':sha(__file__),'producer_source_sha256':EXPECTED_PRODUCER,'certificate_sha256':certificate_hash,
        'precision_bits':precision,'radius_exact':str(rq),'target':certificate['target'],'cache_sha256':certificate['frozen_binary_cache']['sha256'],
        'independence':'New second-order forward interval AD, frozen-array strain reconstruction, full point assembly and full Hessian-variation row bound. Producer assembler/variation routines not called.',
        'no_equilibrium_optimization':True,'preflight_checks':checks}
    save(OUT/'protocol.json',protocol);save(OUT/'job.json',{'status':'RUNNING','pid':os.getpid(),'protocol_sha256':sha(OUT/'protocol.json')})
    if not all(checks.values()):raise RuntimeError('preflight checks failed')
    started=time.monotonic();strength=float(Q(certificate['activation_strength_binary_exact']))
    g,H,cache,counts=assemble(terms_from_arrays(arrays,certificate['parameters'],strength),n)
    checks['all9504_terms_independently_assembled']=counts=={'bulk':9216,'surface':288} and len(cache)==9504
    R=arrays['preconditioner'];RM=arb_mat(n,n,[A(float(x))for x in R.ravel()]);absR=arb_mat(n,n,[A(abs(float(x)))for x in R.ravel()])
    defect=arb_mat(n,n)
    for i in range(n):defect[i,i]=1
    defect-=RM*H;Z0=max(sum((UP(defect[i,j])for j in range(n)),arb(0)).upper()for i in range(n))
    correction=RM*g;eta=max(UP(correction[i,0])for i in range(n))
    vr,domain=variation(cache,n,radius);mapped=absR*vr;Z=Z0+max(UP(mapped[i,0])for i in range(n))
    etaQ=Q(str(eta.upper().fmpq()));ZQ=Q(str(Z.upper().fmpq()))
    checks['independent_exact_rational_contraction']=ZQ<1
    checks['independent_exact_rational_strict_selfmap']=etaQ+ZQ*rq<rq
    checks['analytic_ball_no_penalty_and_positive_surfaces']=all(Q(v['lower_exact'])>0 for v in domain['strict_lower_margins'].values()) and Q(domain['strict_lower_margins']['bulk_J_min']['lower_exact'])>Q(1,1000000)
    # Separate point consistency diagnostic against original Torch AD; binary
    # node addition is rounded there, unlike this exact X+u mathematical target.
    _,tg=model.objective(arrays['u'],strength);TH,meta=model.recover_hessian(arrays['u'],strength)
    gm=np.array([float(g[i,0].mid())for i in range(n)]);Hm=np.array([[float(H[i,j].mid())for j in range(n)]for i in range(n)])
    gdiff=float(np.max(abs(gm-tg)));habs=float(np.max(abs(Hm-TH.toarray())));hrel=habs/max(1.,float(np.max(abs(TH.data))))
    checks['original_AD_point_gradient_consistency']=gdiff<1e-8
    checks['original_AD_point_Hessian_consistency']=hrel<1e-10
    checks['higher_than_producer_precision']=precision>certificate['precision_bits']
    checks['source_and_certificate_unchanged_at_finish']=sha(producer_path)==EXPECTED_PRODUCER and sha(CERT)==certificate_hash
    bound=etaQ/(1-ZQ)if ZQ<1 else None
    result={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'certificate_sha256':certificate_hash,
        'source_sha256':EXPECTED_PRODUCER,'reviewer_source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json'),
        'actual_free_DOF':n,'term_counts':counts,'precision_bits':precision,'radius_exact':str(rq),
        'independent_eta':pack(eta),'independent_center_defect_Z0':pack(Z0),'independent_total_Z':pack(Z),
        'strict_selfmap_gap_exact':str(rq-etaQ-ZQ*rq),'root_distance_upper_exact':str(bound)if bound is not None else None,
        'root_distance_upper_display':float(bound)if bound is not None else None,'analytic_ball_domain':domain,
        'point_AD_comparison':{'gradient_absolute_linf_difference':gdiff,'Hessian_absolute_max_difference':habs,'Hessian_relative_max_difference':hrel,
            'scope':'Numerical mapping/rounding diagnostic, not equality of rounded Torch evaluations and exact-real cache potential'},
        'cache_sha256':certificate['frozen_binary_cache']['sha256'],'seconds':time.monotonic()-started,
        'scope_limits':['One local unique root of the exact-real finite binary-cache potential in this ball only',
            'No energetic stability/global branch or parameter-cell theorem','No quadrature-integration/PDE or source-observation/biological validation',
            'No certificate for new material-cohort fields or all17records','No equilibrium optimization performed']}
    save(OUT/'higher_precision_replay.json',result);save(RECEIPT,result)
    save(OUT/'job.json',{'status':'FINISHED' if all(checks.values())else'FINISHED_WITH_UNRESOLVED_CHECKS','pid':os.getpid(),'receipt_sha256':sha(RECEIPT)})
    print(json.dumps({'status':result['status'],'checks':len(checks),'eta':str(eta),'Z':str(Z),'error':result['root_distance_upper_display'],'seconds':result['seconds']}),flush=True)
    if not all(checks.values()):raise RuntimeError('independent root review has unresolved checks')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['test','verify']);parser.add_argument('--precision',type=int,default=256);args=parser.parse_args()
    if args.stage=='test':local_tests()
    else:verify(args.precision)
