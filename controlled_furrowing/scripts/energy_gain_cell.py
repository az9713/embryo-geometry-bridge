"""Arb energy-norm uniform gain-cell attempt; no equilibrium optimization.

The mathematical family uses the existing frozen binary-quadrature cache.
All norm contractions, predictor widths and acceptance tests here use Arb.
Failed attempts are retained. One gain interval is not a full material domain.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse, base64, hashlib, io, json, time, traceback
import numpy as np
from flint import arb,arb_mat,ctx
import algebraic_root_feasibility as source
from parametric_norm_feasibility import log_data

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'evidence/energy_gain_cell.json'
NORM=ROOT/'evidence/parametric_norm_feasibility.json'
CERT=ROOT/'evidence/algebraic_root_certification_feasibility.json'
LOCKED_NORM_SOURCE='8614a4099e169b9a1c3e72292fd1a57d80ad7f79f3550f4bc9b197cc9cf65192'


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(r):OUT.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n',encoding='utf8')
def upper(x):return source.upper(x)
def norm2(v):return sum((upper(v[i,0])**2 for i in range(v.nrows())),arb(0)).sqrt().upper()
def rootupper(x):return upper(x).sqrt().upper()
def bound(x):return source.pack(x)


def diagonal_majorant(M):
    n=M.nrows();D=arb_mat(n,n)
    for i in range(n):D[i,i]=sum((max(upper(M[i,j]),upper(M[j,i])) for j in range(n)),arb(0)).upper()
    return D


def norm_upper(M):
    """For an exactly symmetric real matrix enclosed by M, spectral<=abs rows."""
    return max(sum((upper(M[i,j]) for j in range(M.ncols())),arb(0)).upper() for i in range(M.nrows()))


def add_local(globalM,localM,indices,weight):
    for i,gi in enumerate(indices):
        for j,gj in enumerate(indices):globalM[gi,gj]+=weight*localM[i,j]


def attempt(halfwidth=Q(1,1000),rho=Q(1,10**6),pieces=8):
    if OUT.exists():raise RuntimeError('preserve existing gain-cell attempt')
    started=time.time();ctx.prec=160
    c=json.loads(CERT.read_text(encoding='utf8'));nr=json.loads(NORM.read_text(encoding='utf8'))
    if sha(Path(__file__).with_name('parametric_norm_feasibility.py'))!=LOCKED_NORM_SOURCE:raise RuntimeError('frozen norm source changed')
    proof=nr['validated_nodal_energy']
    if proof['status']!='CERTIFIED_POSITIVE_HESSIAN_ON_NODAL_CUBE' or not all(proof['checks'].values()):raise RuntimeError('actual metric proof required')
    if any(sha(ROOT/k)!=v for k,v in nr['input_hashes'].items()):raise RuntimeError('metric inputs changed')
    if sha(source.__file__)!=c['script_sha256']:raise RuntimeError('point source changed')
    raw=base64.b64decode(proof['W_cache']['data'])
    if hashlib.sha256(raw).hexdigest()!=proof['W_cache']['sha256']:raise RuntimeError('metric cache changed')
    with np.load(io.BytesIO(raw),allow_pickle=False) as z:W=z['W'].copy()
    model,u,_=source.restore_cache(c);n=len(u)
    with np.load(ROOT/c['matrix_archive'],allow_pickle=False) as z:S=z['displacement_tangents'][:,5].copy()
    if W.shape!=(n,n) or len(S)!=n:raise RuntimeError('metric/tangent dimension mismatch')
    WM=arb_mat(n,n,[source.A(float(v)) for v in W.ravel()]);SM=arb_mat(n,1,[source.A(float(v)) for v in S])
    metric_cov=WM*WM.transpose()
    H0=arb_mat(n,n);g0=arb_mat(n,1);Bg=arb_mat(n,1);Hg=arb_mat(n,n);A0=arb_mat(n,n)
    plans=[{'name':f'predictor_path_piece_{i}','shape_halfwidth':halfwidth*Q(i,pieces),'error_radius':Q(0),
            'time_interval':[str(Q(i-1,pieces)),str(Q(i,pieces))]} for i in range(1,pieces+1)]
    plans.append({'name':'full_error_tube','shape_halfwidth':halfwidth,'error_radius':rho})
    for p in plans:p.update(D=arb_mat(n,n),Q=arb_mat(n,n),surface=arb_mat(n,n),margins={})
    record={'status':'RUNNING_VALIDATED_ENERGY_GAIN_CELL','script_sha256':sha(__file__),'no_equilibrium_optimization':True,
        'input_hashes':{CERT.relative_to(ROOT).as_posix():sha(CERT),NORM.relative_to(ROOT).as_posix():sha(NORM),c['matrix_archive']:sha(ROOT/c['matrix_archive'])},
        'source_hashes':{'parametric_norm_feasibility.py':LOCKED_NORM_SOURCE,'algebraic_root_feasibility.py':sha(source.__file__),**c['source_hashes']},
        'gain_relative_halfwidth_exact':str(halfwidth),'metric_error_radius_exact':str(rho),'predictor_integral_pieces':pieces,
        'metric_cache_sha256':proof['W_cache']['sha256'],'model_cache_sha256':c['frozen_binary_cache']['sha256'],
        'predictor_column':'archived binary log_common_gain tangent, used as derivative with respect to relative linear gain at center',
        'fixed_parameters':c['parameters'],'precision_bits':ctx.prec,
        'scope':'One relative gain interval, fixed original material and finite binary-quadrature cache; no continuum/full-Theta/all17 certificate.'}
    save(record);count=0
    for term in source.local_terms(model,u,float(Q(c['activation_strength_binary_exact']))):
        ids=term['indices'];d=len(ids);k=len(term['B']);BM=arb_mat(k,d,[v for row in term['B'] for v in row])
        subP=arb_mat(d,d,[metric_cov[i,j] for i in ids for j in ids]);Gram=BM*subP*BM.transpose()
        subS=arb_mat(d,1,[SM[i,0] for i in ids]);direction=BM*subS
        rownorm=[rootupper(Gram[i,i]) for i in range(k)]
        fun=source.bulk_local if term['kind']=='bulk' else source.surface_local
        grad,hc,_=fun(term['q'],*term['parameters']);HC=arb_mat(k,k,[x for row in hc for x in row]);gc=arb_mat(k,1,grad)
        add_local(H0,BM.transpose()*HC*BM,ids,term['weight'])
        gl=BM.transpose()*gc
        for i,gi in enumerate(ids):g0[gi,0]+=term['weight']*gl[i,0]
        if term['kind']=='bulk':
            ell,v,Q0=log_data(term['q']);VM=arb_mat(1,k,v);QM=arb_mat(k,k,[x for row in Q0 for x in row])
            kap=term['parameters'][1];mu=term['parameters'][0];drive=term['parameters'][2]
            av=VM*BM;add_local(A0,av.transpose()*av,ids,term['weight']*kap)
            Gdiag=arb_mat(k,k);Gdiag[1,1]=drive;Gdiag[3,3]=drive
            add_local(Hg,BM.transpose()*Gdiag*BM,ids,term['weight'])
            ag=BM.transpose()*arb_mat(k,1,[0,drive*term['q'][1],0,drive*term['q'][3],0])
            for i,gi in enumerate(ids):Bg[gi,0]+=term['weight']*ag[i,0]
            vd=upper((VM*direction)[0,0]);vnorm=rootupper((VM*Gram*VM.transpose())[0,0])
            Qdirection=QM*direction;QGram=QM*Gram*QM.transpose()
            Qrownorm=[rootupper(QGram[i,i]) for i in range(k)]
        for p in plans:
            sh=source.A(p['shape_halfwidth']);er=source.A(p['error_radius'])
            widths=[(sh*upper(direction[i,0])+er*rownorm[i]).upper() for i in range(k)]
            qb=[q+arb(0,width) for q,width in zip(term['q'],widths)]
            if term['kind']=='bulk':
                det=qb[0]*qb[3]-qb[1]*qb[2];J=det*qb[4]
                if not det>0 or not qb[4]>0 or not J>source.A(Q(1,10**6)):raise ValueError('tube leaves analytic positive-J domain')
                _,_,Qb=log_data(qb)
                Qabs=[[upper(x) for x in row] for row in Qb]
                Qdelta=[[upper(Qb[i][j]-Q0[i][j]) for j in range(k)] for i in range(k)]
                dell=(sh*vd+er*vnorm+sum((widths[i]*Qabs[i][j]*widths[j]/2 for i in range(k) for j in range(k)),arb(0))).upper()
                dv=[(sh*upper(Qdirection[i,0])+er*Qrownorm[i]+sum((Qdelta[i][j]*widths[j] for j in range(k)),arb(0))).upper() for i in range(k)]
                sumdv=sum(dv,arb(0)).upper();DM=arb_mat(k,k)
                for i in range(k):DM[i,i]=(kap*sumdv*dv[i]).upper()
                add_local(p['D'],BM.transpose()*DM*BM,ids,term['weight'])
                cf=upper(kap*ell-mu)
                change=arb_mat(k,k,[cf*Qdelta[i][j]+kap*dell*Qabs[i][j] for i in range(k) for j in range(k)])
                add_local(p['Q'],BM.transpose()*diagonal_majorant(change)*BM,ids,term['weight'])
                margins={'det':det.lower(),'hoop':qb[4].lower(),'J':J.lower()}
            else:
                _,hb,_=source.surface_local(qb,*term['parameters'])
                change=arb_mat(k,k,[hb[i][j]-hc[i][j] for i in range(k) for j in range(k)])
                add_local(p['surface'],BM.transpose()*diagonal_majorant(change)*BM,ids,term['weight'])
                margins={'surface_R':qb[0].lower(),'surface_tangent_norm2':(qb[1]*qb[1]+qb[2]*qb[2]).lower()}
            for name,value in margins.items():p['margins'][name]=value if name not in p['margins'] else min(p['margins'][name],value)
        count+=1
        if count%1536==0:
            record['progress']={'phase':'Arb local majorants','terms':count,'elapsed_seconds':time.time()-started};save(record);print(json.dumps(record['progress']),flush=True)
    transformed=WM.transpose()*H0*WM;defect=transformed-arb_mat(n,n,[1 if i==j else 0 for i in range(n) for j in range(n)])
    eps0=norm_upper(defect);a=norm_upper(WM.transpose()*A0*WM);gain=source.A(halfwidth)*norm_upper(WM.transpose()*Hg*WM)
    rows=[]
    for p in plans:
        print(json.dumps({'phase':'Arb transformed majorants','plan':p['name'],'elapsed_seconds':time.time()-started}),flush=True)
        dd=norm_upper(WM.transpose()*p['D']*WM);qq=norm_upper(WM.transpose()*p['Q']*WM);ss=norm_upper(WM.transpose()*p['surface']*WM)
        rank=2*(a*dd).sqrt()+dd;total=rank+qq+ss+gain
        rows.append({'name':p['name'],'shape_halfwidth_exact':str(p['shape_halfwidth']),'metric_error_radius_exact':str(p['error_radius']),
            'time_interval':p.get('time_interval'),'rank_center_squared_norm_upper':bound(a),'rank_change_squared_norm_upper':bound(dd),
            'volumetric_rank_variation_bound':bound(rank),'bulk_logQ_variation_bound':bound(qq),'surface_variation_bound':bound(ss),
            'gain_parameter_variation_bound':bound(gain),'total_transformed_Hessian_variation_bound':bound(total),
            'analytic_domain_strict_lower_margins':{k:bound(v) for k,v in p['margins'].items()}})
        record['completed_plan_bounds']=rows;save(record)
    eta0=norm2(WM.transpose()*g0);linear=norm2(WM.transpose()*(H0*SM+Bg))
    snorm=((SM.transpose()*H0*SM)[0,0]/(1-eps0)).sqrt().upper()
    integral=sum((source.A(Q(row['total_transformed_Hessian_variation_bound']['upper_exact']))/pieces for row in rows[:-1]),arb(0)).upper()
    eta=eta0+source.A(halfwidth)*linear+source.A(halfwidth)*snorm*integral
    Z=eps0+source.A(Q(rows[-1]['total_transformed_Hessian_variation_bound']['upper_exact']))
    lhs=eta+Z*source.A(rho);passed=bool(Z<1 and lhs<source.A(rho))
    record.update(status='CERTIFIED_UNIFORM_ENERGY_GAIN_CELL' if passed else 'UNRESOLVED_ENERGY_GAIN_CELL',
        terms=count,center_metric_defect=bound(eps0),center_residual_energy_norm=bound(eta0),predictor_linear_defect_energy_norm=bound(linear),
        predictor_metric_norm_upper=bound(snorm),path_integral_H_variation_upper=bound(integral),uniform_predictor_eta=bound(eta),uniform_contraction_Z=bound(Z),
        eta_plus_Zrho=bound(lhs),metric_radius=bound(source.A(rho)),selfmap_margin=bound(source.A(rho)-lhs),
        checks={'all_9504_terms':count==9504,'center_metric_positive':bool(eps0<1),'strict_contraction':bool(Z<1),'strict_selfmap':bool(lhs<source.A(rho)),
                'all_input_hashes_unchanged':all(sha(ROOT/k)==v for k,v in record['input_hashes'].items()),
                'all_sources_unchanged':all(sha(Path(__file__).with_name(k))==v for k,v in record['source_hashes'].items()),'no_equilibrium_solve':True},
        root_error_about_predictor_metric_upper=bound(eta/(1-Z)) if Z<1 else None,
        completed_elapsed_seconds=time.time()-started,
        conclusion='For every relative gain in the declared interval, a unique stationary point exists in its predictor-centered energy ball; strict local energetic minimum follows from the symmetric transformed Hessian bound.' if passed else 'No uniform-cell existence claim. Failed acceptance bounds retained.',
        limitations=['one fixed original material and weak-load finite binary-quadrature family','quadrature positivity does not certify whole-body geometry','no physical dynamics/global branch/continuum/fullTheta or all17 result'])
    save(record);print(json.dumps({k:v for k,v in record.items() if k not in ('completed_plan_bounds','source_hashes','fixed_parameters')},indent=2),flush=True)
    return passed


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',action='store_true');args=ap.parse_args()
    if not args.run:raise SystemExit('Explicit --run required; no equilibrium optimization')
    try:attempt()
    except Exception as err:
        if OUT.exists():
            r=json.loads(OUT.read_text(encoding='utf8'));r['status']='FAILED_ENERGY_CELL_ATTEMPT';r['error']=repr(err);r['traceback']=traceback.format_exc();save(r)
        raise
