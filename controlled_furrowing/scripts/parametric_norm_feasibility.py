"""Read-only norm/decomposition feasibility at one frozen weak-load field.

No equilibrium solve. Local interval Hessian majorants are rigorous Arb
enclosures; global spectral/trace contractions below are floating feasibility
diagnostics unless explicitly marked otherwise. Never a full-Theta certificate.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse, base64, hashlib, io, json, time
import numpy as np
from scipy.linalg import cholesky, solve_triangular, eigvalsh
from scipy.sparse import csr_matrix
from flint import arb, arb_mat, ctx
import algebraic_root_feasibility as source

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'evidence/parametric_norm_feasibility.json'


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def upper_float(x):
    """Outward binary upper of a nonnegative Arb bound."""
    q=Q(str(source.upper(x).fmpq()));v=float(q)
    return np.nextafter(v,np.inf) if Q.from_float(v)<q else v


def log_data(q):
    a,b,c,d,h=q;det=a*d-b*c
    first=[d,-c,-b,a];v=[x/det for x in first]+[1/h]
    second=[[arb(0) for _ in range(5)] for _ in range(5)]
    for i in range(4):
        for j in range(4):
            cross=1 if (i,j) in ((0,3),(3,0)) else -1 if (i,j) in ((1,2),(2,1)) else 0
            second[i][j]=cross/det-first[i]*first[j]/det**2
    second[4][4]=-1/h**2
    return det.log()+h.log(),v,second


def structured(residual_screen=False):
    """Additional screen preserving the volumetric outer-product structure."""
    if not OUT.exists():raise RuntimeError('retain initial screen first')
    record=json.loads(OUT.read_text(encoding='utf8'))
    result_key='gain_residual_screen' if residual_screen else 'structured_rank_screen'
    if result_key in record:raise RuntimeError('preserve completed structured screen')
    started=time.time();ctx.prec=128
    cert=json.loads((ROOT/'evidence/algebraic_root_certification_feasibility.json').read_text(encoding='utf8'))
    model,u,R=source.restore_cache(cert)
    with np.load(ROOT/cert['matrix_archive'],allow_pickle=False) as z:
        H=csr_matrix((z['H_data'],z['H_indices'],z['H_indptr']),shape=tuple(z['H_shape'])).toarray();S=z['displacement_tangents'].copy()
    n=len(u);L=cholesky(H,lower=True);W=solve_triangular(L.T,np.eye(n))
    if residual_screen:
        plans=[{'name':'gain_predictor_point1percent_zero_tube','nodal_radius':0.,'energy_radius':0.,'relative_gain_halfwidth':.001},
               {'name':'gain_predictor_point1percent_tube_1e-6','nodal_radius':0.,'energy_radius':1e-6,'relative_gain_halfwidth':.001}]
    else:plans=[{k:p[k] for k in ('name','nodal_radius','energy_radius','relative_gain_halfwidth')} for p in record['plans']]
    for p in plans:p.update(D=np.zeros((n,n)),Qchange=np.zeros((n,n)),surface=np.zeros((n,n)),gain=np.zeros((n,n)),max_log_variation=0.)
    A0=np.zeros((n,n));count=0
    for term in source.local_terms(model,u,float(Q(cert['activation_strength_binary_exact']))):
        ids=np.array(term['indices'],int);ix=np.ix_(ids,ids)
        B=np.array([[float(x.mid()) for x in row] for row in term['B']]);BW=B@W[ids,:]
        d=B@S[ids,5];w=float(term['weight'].mid());kind=term['kind']
        if kind=='bulk':
            ell,v,Q0=log_data(term['q']);vf=np.array([float(x.mid()) for x in v]);Qf=np.array([[float(x.mid()) for x in row] for row in Q0])
            mu,kap,drive=map(lambda x:float(x.mid()),term['parameters']);av=vf@B
            A0[ix]+=w*kap*np.outer(av,av)
        else:_,H0,_=source.surface_local(term['q'],*term['parameters'])
        for p in plans:
            radii=p['nodal_radius']*np.abs(B).sum(axis=1)+p['energy_radius']*np.linalg.norm(BW,axis=1)+p['relative_gain_halfwidth']*np.abs(d)
            qb=[q+arb(0,source.A(float(np.nextafter(t*(1+1e-12),np.inf)))) for q,t in zip(term['q'],radii)]
            if kind=='bulk':
                _,vb,Qb=log_data(qb)
                Qabs=np.array([[upper_float(x) for x in row] for row in Qb]);Qdelta=np.array([[upper_float(Qb[a][b]-Q0[a][b]) for b in range(5)] for a in range(5)])
                # Preserve directional cancellations in the first derivative.
                ellfirst=p['nodal_radius']*np.abs(av).sum()+p['energy_radius']*np.linalg.norm(av@W[ids,:])+p['relative_gain_halfwidth']*abs(vf@d)
                dell=ellfirst+.5*float(radii@Qabs@radii)
                # Integral Taylor form for grad(log J); center contractions first.
                Q_B=Qf@B;Q_BW=Qf@BW
                dv=p['nodal_radius']*np.abs(Q_B).sum(axis=1)+p['energy_radius']*np.linalg.norm(Q_BW,axis=1)+p['relative_gain_halfwidth']*np.abs(Qf@d)+Qdelta@radii
                # (dv.y)^2 <= sum(dv_i)*sum(dv_i*y_i^2).
                dd=kap*dv.sum()*dv
                p['D'][ix]+=w*(B.T@(dd[:,None]*B))
                c0=upper_float(term['parameters'][1]*ell-term['parameters'][0])
                rem=c0*Qdelta+kap*dell*Qabs
                rd=np.maximum(rem,rem.T).sum(axis=1)
                p['Qchange'][ix]+=w*(B.T@(rd[:,None]*B))
                dg=np.array([0.,drive,0.,drive,0.])*p['relative_gain_halfwidth']
                p['gain'][ix]+=w*(B.T@(dg[:,None]*B))
                p['max_log_variation']=max(p['max_log_variation'],dell)
            else:
                _,Hb,_=source.surface_local(qb,*term['parameters'])
                delta=np.array([[upper_float(Hb[a][b]-H0[a][b]) for b in range(3)] for a in range(3)])
                sd=np.maximum(delta,delta.T).sum(axis=1)
                p['surface'][ix]+=w*(B.T@(sd[:,None]*B))
        count+=1
        if count%1536==0:print(json.dumps({'phase':'rank_structure','terms':count,'elapsed':time.time()-started}),flush=True)
    def spectral(C):
        X=W.T@C@W;return float(eigvalsh(X,subset_by_index=[n-1,n-1])[0]),float(np.abs(X).sum(axis=1).max())
    a,a_row=spectral(A0);rows=[]
    for p in plans:
        d,drow=spectral(p['D']);q,qrow=spectral(p['Qchange']);s,srow=spectral(p['surface']);g,grow=spectral(p['gain'])
        rank=2*np.sqrt(a*d)+d;rankrow=2*np.sqrt(a_row*drow)+drow
        rows.append({k:v for k,v in p.items() if not isinstance(v,np.ndarray)}|{
            'rank_center_norm_squared_float':a,'rank_delta_norm_squared_majorant_float':d,
            'bulk_rank_variation_bound_float':float(rank),'bulk_logQ_variation_bound_float':q,
            'surface_variation_bound_float':s,'gain_coefficient_variation_bound_float':g,
            'total_rank_structured_bound_float':float(rank+q+s+g),
            'absrow_only_rank_structured_bound_float':float(rankrow+qrow+srow+grow)})
    result={'status':'FLOATING_FEASIBILITY_ONLY_NOT_CERTIFIED_CELL','source_sha256':sha(__file__),'plans':rows,
        'formula':'||Delta H||metric <= 2||A0W||||DW||+||DW||^2+||Qremainder||metric+surface+gain; first-order logJ and gradlogJ predictor contractions performed before abs.',
        'elapsed_seconds':time.time()-started,'caveat':'Arb local interval derivatives but matrix assembly, Cholesky metric contractions and BW/BS norms floating; these numbers do not establish a validated parameter cell.'}
    if residual_screen:
        # This is only an inexpensive operator-based predictor-residual screen.
        # Validation additionally needs exact center residual/linear-defect norms.
        snorm=float(np.linalg.norm(L.T@S[:,5]));h=.001;rho=1e-6
        eta=rows[0]['absrow_only_rank_structured_bound_float']*h*snorm
        Z=rows[1]['absrow_only_rank_structured_bound_float']
        result['residual_screen']={'gain_halfwidth':h,'error_metric_radius':rho,'predictor_energy_length_float':h*snorm,
            'nonlinear_remainder_eta_upper_float':eta,'eta_plus_Zrho_float':eta+Z*rho,
            'unvalidated_operator_screen_pass':eta+Z*rho<rho,
            'omitted_until_validation':'center F0 and numerical tangent linear-system defect, outward BW/BS/global norm bounds; neither actual uniform cell nor root existence asserted'}
    record[result_key]=result
    OUT.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n',encoding='utf8');print(json.dumps(result,indent=2),flush=True)


def validate_nodal_energy(replay=False):
    """Actual Arb SPD validation on a fixed1e-8 nodal cube, not a material cell."""
    record=json.loads(OUT.read_text(encoding='utf8'))
    if 'validated_nodal_energy' in record and not replay:raise RuntimeError('retain completed validation')
    if replay and 'validated_nodal_energy' not in record:raise RuntimeError('no retained validation to replay')
    started=time.time();ctx.prec=128
    certpath=ROOT/'evidence/algebraic_root_certification_feasibility.json';cert=json.loads(certpath.read_text(encoding='utf8'))
    model,u,_=source.restore_cache(cert);n=len(u)
    with np.load(ROOT/cert['matrix_archive'],allow_pickle=False) as z:
        Hfloat=csr_matrix((z['H_data'],z['H_indices'],z['H_indptr']),shape=tuple(z['H_shape'])).toarray()
    W=solve_triangular(cholesky(Hfloat,lower=True).T,np.eye(n))
    if not np.array_equal(W,np.triu(W)) or np.any(np.diag(W)==0):raise RuntimeError('invalid exact binary metric')
    H=arb_mat(n,n);M=arb_mat(n,n);radius=source.A(Q(1,10**8));terms=0
    for term in source.local_terms(model,u,float(Q(cert['activation_strength_binary_exact']))):
        fun=source.bulk_local if term['kind']=='bulk' else source.surface_local
        _,hc,_=fun(term['q'],*term['parameters']);k=len(hc);indices=term['indices'];d=len(indices)
        widths=[sum((source.upper(v) for v in row),arb(0))*radius for row in term['B']]
        qb=[q+arb(0,w.upper()) for q,w in zip(term['q'],widths)]
        _,hb,_=fun(qb,*term['parameters'])
        delta=[[max(source.upper(hb[a][b]-hc[a][b]),source.upper(hb[b][a]-hc[b][a])) for b in range(k)] for a in range(k)]
        D=[sum(row,arb(0)).upper() for row in delta]
        BM=arb_mat(k,d,[v for row in term['B'] for v in row]);HM=arb_mat(k,k,[v for row in hc for v in row])
        DM=arb_mat(k,k)
        for a in range(k):DM[a,a]=D[a]
        localH=BM.transpose()*HM*BM;localM=BM.transpose()*DM*BM;weight=term['weight']
        for a,i in enumerate(indices):
            for b,j in enumerate(indices):H[i,j]+=weight*localH[a,b];M[i,j]+=weight*localM[a,b]
        terms+=1
        if terms%1536==0:print(json.dumps({'phase':'validated_metric_assembly','terms':terms,'elapsed':time.time()-started}),flush=True)
    WM=arb_mat(n,n,[source.A(float(v)) for v in W.ravel()])
    print(json.dumps({'phase':'validated_metric_matrix_products','elapsed':time.time()-started}),flush=True)
    center=WM.transpose()*H*WM;majorant=WM.transpose()*M*WM
    center_defect=arb(0);variation=arb(0)
    for i in range(n):
        row=sum((source.upper(center[i,j]-(1 if i==j else 0)) for j in range(n)),arb(0)).upper()
        center_defect=max(center_defect,row)
        row=sum((source.upper(majorant[i,j]) for j in range(n)),arb(0)).upper();variation=max(variation,row)
    total=center_defect+variation;margin=1-total
    buff=io.BytesIO();np.savez_compressed(buff,W=W);raw=buff.getvalue()
    root_success=[a for a in cert['radius_attempts'] if a['strict_selfmap_and_contraction']]
    existing_inside=any(Q(a['radius_requested_exact'])<Q(1,10**8) for a in root_success)
    checks={'all_9504_terms':terms==9504,'binary_W_triangular_nonsingular':bool(np.array_equal(W,np.triu(W)) and np.all(np.diag(W)!=0)),
        'center_metric_defect_below_one':bool(center_defect<1),'entire_cube_positive_metric_Hessian':bool(total<1),
        'existing_certified_root_strictly_inside_cube':existing_inside,
        'producer_source_and_inputs_current':sha(source.__file__)==cert['script_sha256'] and all(sha(ROOT/p)==v for p,v in record['input_hashes'].items())}
    result={'status':'CERTIFIED_POSITIVE_HESSIAN_ON_NODAL_CUBE' if all(checks.values()) else 'FAILED_ENERGY_VALIDATION','checks':checks,
        'script_sha256':sha(__file__),'norm':'||W^-1 delta_u||2, W exact frozen binary matrix; G=W^-T W^-1',
        'coordinate_cube_radius_exact':'1/100000000','center_metric_defect':source.pack(center_defect),'metric_variation_bound':source.pack(variation),
        'combined_bound':source.pack(total),'positive_coercivity_margin':source.pack(margin),
        'W_cache':{'encoding':'base64 np.savez_compressed exact binary64','sha256':hashlib.sha256(raw).hexdigest(),'data':base64.b64encode(raw).decode(),'shape':list(W.shape)},
        'precision_bits':ctx.prec,'terms':terms,'elapsed_seconds':time.time()-started,
        'scope':'Strict convexity and strict local energetic minimality of the already certified one weak-field binary-quadrature root within the declared nodal cube. Not physical dynamics, global branch, parameter cell, quadrature integration, continuum or biological validation.'}
    if replay:
        old=record['validated_nodal_energy']
        keys=('center_metric_defect','metric_variation_bound','combined_bound','positive_coercivity_margin')
        if not all(result[k]==old[k] for k in keys) or result['W_cache']['sha256']!=old['W_cache']['sha256']:
            raise RuntimeError('read-only validation replay differs from retained proof')
        result['read_only_replay_exact_bounds_match']=True
    else:
        record['validated_nodal_energy']=result;OUT.write_text(json.dumps(record,indent=2)+'\n',encoding='utf8')
    print(json.dumps({k:v for k,v in result.items() if k!='W_cache'},indent=2),flush=True)
    if not all(checks.values()):raise RuntimeError('Retained validation failure')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run',action='store_true');parser.add_argument('--structured',action='store_true');parser.add_argument('--validate-nodal',action='store_true');parser.add_argument('--residual-screen',action='store_true');parser.add_argument('--replay-nodal',action='store_true')
    args=parser.parse_args()
    if args.validate_nodal:return validate_nodal_energy()
    if args.replay_nodal:return validate_nodal_energy(replay=True)
    if args.residual_screen:return structured(residual_screen=True)
    if args.structured:return structured()
    if not args.run:raise SystemExit('Explicit --run for read-only feasibility calculation')
    if OUT.exists():raise RuntimeError('Preserve existing feasibility receipt')
    started=time.time();ctx.prec=128
    certpath=ROOT/'evidence/algebraic_root_certification_feasibility.json'
    cert=json.loads(certpath.read_text(encoding='utf8'))
    archive=ROOT/cert['matrix_archive']
    if sha(archive)!=cert['matrix_archive_sha256']:raise RuntimeError('stale matrix archive')
    model,u,R=source.restore_cache(cert)
    with np.load(archive,allow_pickle=False) as saved:
        H=csr_matrix((saved['H_data'],saved['H_indices'],saved['H_indptr']),shape=tuple(saved['H_shape'])).toarray()
        S=saved['displacement_tangents'].copy();Bparameter=saved['parameter_force_derivatives'].copy()
    n=len(u);L=cholesky((H+H.T)/2,lower=True);W=solve_triangular(L.T,np.eye(n))
    # Exact binary W is nonsingular because it is triangular with nonzero diagonal.
    # The metric is defined from W, not by assuming numerical Cholesky was exact.
    wrowmax=float(np.linalg.norm(W,axis=1).max())
    columns=['relative_kappa','relative_apical_area','relative_basal_area','additive_apical_tension','relative_basal_tension','relative_common_gain']
    plans=[
        {'name':'nodal_cube_1e-8','nodal_radius':1e-8,'energy_radius':0.,'relative_gain_halfwidth':0.},
        {'name':'energy_ball_1e-4','nodal_radius':0.,'energy_radius':1e-4,'relative_gain_halfwidth':0.},
        {'name':'gain_predictor_1percent_tube_1e-5','nodal_radius':0.,'energy_radius':1e-5,'relative_gain_halfwidth':.01},
        {'name':'gain_predictor_point1percent_tube_1e-7','nodal_radius':0.,'energy_radius':1e-7,'relative_gain_halfwidth':.001},
    ]
    for p in plans:p.update(E=np.zeros((n,n)),parameter_E=np.zeros((n,n)),local_trace_upper=0.,terms=0)
    Hassembled=np.zeros_like(H);component=[np.zeros_like(H) for _ in columns]
    minJ=1.;kinds={'bulk':0,'surface':0};strength=float(Q(cert['activation_strength_binary_exact']))
    for term in source.local_terms(model,u,strength):
        kind=term['kind'];kinds[kind]+=1;ids=np.array(term['indices'],int)
        fun=source.bulk_local if kind=='bulk' else source.surface_local
        _,hc,J=fun(term['q'],*term['parameters']);minJ=min(minJ,float(J.lower())) if kind=='bulk' else minJ
        B=np.array([[float(x.mid()) for x in row] for row in term['B']]);weight=float(term['weight'].mid())
        BW=B@W[ids,:];rowabs=np.abs(B).sum(axis=1);row2=np.linalg.norm(BW,axis=1)
        hcfloat=np.array([[float(x.mid()) for x in row] for row in hc]);ix=np.ix_(ids,ids)
        Hassembled[ix]+=weight*(B.T@hcfloat@B)
        # Hessians are affine in physical coefficients at fixed geometry.
        if kind=='bulk':
            mu,kappa,drive=term['parameters']
            candidates=[(0,(mu,2*kappa,drive)),(5,(mu,kappa,2*drive))]
        else:
            K,T,scale=term['parameters'];isbasal=bool(T>0)
            candidates=[(2 if isbasal else 1,(2*K,T,scale)),(4 if isbasal else 3,(K,2*T if isbasal else T+1,scale))]
        for col,pars in candidates:
            _,hj,_=fun(term['q'],*pars)
            mat=np.array([[float((hj[a][b]-hc[a][b]).mid()) for b in range(len(hc))] for a in range(len(hc))])
            component[col][ix]+=weight*(B.T@mat@B)
        for plan in plans:
            direction=B@S[ids,5]
            # Float radii multiplied by outward norm inflation. This is only a
            # feasibility screen; final norm certificates must enclose BW too.
            qradii=plan['nodal_radius']*rowabs+plan['energy_radius']*row2+plan['relative_gain_halfwidth']*np.abs(direction)
            qball=[v+arb(0,source.A(float(np.nextafter(d*(1+1e-12),np.inf)))) for v,d in zip(term['q'],qradii)]
            _,hb,_=fun(qball,*term['parameters'])
            change=np.array([[upper_float(hb[a][b]-hc[a][b]) for b in range(len(hc))] for a in range(len(hc))])
            change=np.maximum(change,change.T);D=change.sum(axis=1)
            # Symmetric diagonal-dominance gives -diag(D)<=Delta h<=diag(D).
            plan['E'][ix]+=weight*(B.T@(D[:,None]*B))
            plan['local_trace_upper']+=weight*float(np.dot(D,row2*row2))
            if kind=='bulk' and plan['relative_gain_halfwidth']:
                # Activation is quadratic: its Hessian is geometry independent.
                drive=upper_float(term['parameters'][2])
                dg=np.array([0.,drive,0.,drive,0.])*plan['relative_gain_halfwidth']
                plan['parameter_E'][ix]+=weight*(B.T@(dg[:,None]*B))
            plan['terms']+=1
        if sum(kinds.values())%1536==0:print(json.dumps({'phase':'local_majorants','terms':sum(kinds.values()),'elapsed':time.time()-started}),flush=True)
    rows=[]
    for plan in plans:
        transformed=W.T@plan['E']@W;parameter=W.T@plan['parameter_E']@W
        geom=float(eigvalsh(transformed,subset_by_index=[n-1,n-1])[0]);param=float(eigvalsh(parameter,subset_by_index=[n-1,n-1])[0])
        rows.append({k:v for k,v in plan.items() if k not in ('E','parameter_E')}|{
            'geometry_majorant_generalized_lambda_max_float':geom,
            'geometry_majorant_transformed_absrows_float':float(np.abs(transformed).sum(axis=1).max()),
            'geometry_majorant_trace_float':float(np.trace(transformed)),
            'relative_gain_H_majorant_lambda_max_float':param,
            'raw_nodal_inf_preconditioned_abs_majorant_float':float((np.abs(R)@np.abs(plan['E']).sum(axis=1)).max()),
            'total_local_majorant_lambda_float':geom+param,
            'gain_predictor_energy_displacement_upper_float':plan['relative_gain_halfwidth']*float(np.linalg.norm(L.T@S[:,5])),
            'energy_ball_implied_nodal_error_upper_float':plan['energy_radius']*wrowmax,
        })
    components=[]
    for label,C in zip(columns,component):
        vals=eigvalsh(W.T@C@W)
        components.append({'coordinate':label,'generalized_eigenvalue_range_float':[float(vals[0]),float(vals[-1])],
                           'absolute_relative_norm_float':float(np.max(np.abs(vals))),
                           'energy_tangent_norm_float':float(np.linalg.norm(L.T@S[:,len(components)])),
                           'nodal_tangent_linf_float':float(np.max(np.abs(S[:,len(components)])))})
    checks={'source_and_matrix_archive_bound':sha(Path(source.__file__))==cert['script_sha256'] and sha(archive)==cert['matrix_archive_sha256'],
            'exact_binary_metric_nonsingular_triangular':bool(np.all(np.diag(W)!=0) and np.array_equal(W,np.triu(W))),
            'archived_H_symmetric':bool(np.array_equal(H,H.T)),
            'independent_local_H_matches_archive':bool(np.max(np.abs(Hassembled-H))<1e-8),
            'all_six_tangent_linear_systems_match':bool(np.max(np.abs(H@S+Bparameter))<1e-8),
            'local_Arb_domains_positive':minJ>0,
            'no_equilibrium_solve_performed':True}
    result={'status':'COMPLETE_FLOATING_GLOBAL_NORM_FEASIBILITY_WITH_ARB_LOCAL_ENCLOSURES','checks':checks,
        'script_sha256':sha(__file__),'input_hashes':{certpath.relative_to(ROOT).as_posix():sha(certpath),archive.relative_to(ROOT).as_posix():sha(archive)},
        'frozen_cache_sha256':cert['frozen_binary_cache']['sha256'],'free_dofs':n,'term_counts':kinds,'precision_bits':ctx.prec,
        'metric_definition':'W=binary inverse transpose of a floating Cholesky factor; norm ||W^-1 delta_u||2. W is exactly triangular/nonsingular; numerical Cholesky is not assumed exact.',
        'metric_center_identity_defect_absrows_float':float(np.abs(np.eye(n)-W.T@H@W).sum(axis=1).max()),
        'max_metric_to_nodal_linf_constant_float':wrowmax,'independent_H_linf_error_float':float(np.max(np.abs(Hassembled-H))),
        'point_min_quadrature_J_float':minJ,'plans':rows,'coefficient_components':components,
        'scope':'Fixed finite-binary-quadrature weak state only. Local Arb enclosures do not certify floating matrix contractions/norms/parameter tubes. No equilibrium optimization, material identification, continuum or full-domain certificate.',
        'elapsed_seconds':time.time()-started}
    OUT.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf8')
    print(json.dumps({'status':result['status'],'checks':checks,'plans':rows,'elapsed':result['elapsed_seconds']},indent=2),flush=True)
    if not all(checks.values()):raise RuntimeError('Feasibility integrity check failed')


if __name__=='__main__':main()
