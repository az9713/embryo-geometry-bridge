"""Independent exact sym(R) positivity and inherited uniform Hessian coercivity.
Frozen K is only a proposal. No producer numerical-certifier call, no H assembly,
equilibrium solve, new fitting, physical source edit or prior receipt rewrite.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,base64,hashlib,importlib.metadata,io,json,time
import numpy as np
from flint import arb,arb_mat,fmpq,ctx
ROOT=Path(__file__).resolve().parents[1]
TARGET=ROOT/'runs/cohort_energy_stability_v1'
RUN=ROOT/'runs/cohort_energy_stability_review_v1'
OUT=ROOT/'evidence/cohort_energy_stability_independent_review.json'
JOIN=ROOT/'evidence/shared_cell_cohort_inference.json'
THEORY=ROOT/'evidence/cohort_energy_stability_theory_review.json'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
rel=lambda p:Path(p).resolve().relative_to(ROOT.resolve()).as_posix()
yes=lambda d:isinstance(d,dict) and bool(d) and all(v is True for v in d.values())

def path(name):
    p=(ROOT/name).resolve()
    if not p.is_relative_to(ROOT.resolve()):raise ValueError('outside study child')
    return p

def write(p,value):
    p=Path(p)
    if p.exists():raise RuntimeError('preserve immutable independent artifact '+str(p))
    p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf8')

def rational(v):
    if isinstance(v,(bool,np.bool_)):raise ValueError('boolean not numerical bound')
    if isinstance(v,(float,np.floating)):
        if not np.isfinite(v):raise ValueError('nonfinite binary value')
        return Q(*float(v).as_integer_ratio())
    return Q(v)

def A(v):
    q=rational(v);return arb(fmpq(q.numerator,q.denominator))

def mat(v,n=None):
    v=np.asarray(v)
    if v.dtype!=np.dtype('float64') or v.ndim!=2 or not len(v) or v.shape[0]!=v.shape[1] or (n is not None and v.shape!=(n,n)) or not np.isfinite(v).all():raise ValueError('finite square binary64 matrix required')
    return v

def decode(cache):
    p=cache['frozen_binary_cache'];raw=base64.b64decode(p['data'],validate=True)
    if hashlib.sha256(raw).hexdigest()!=p['sha256']:raise ValueError('binary cache checksum')
    with np.load(io.BytesIO(raw),allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}

def norm_R_integer(R):
    """Independent exact dyadic-row sums via ONE common power-of-two denominator."""
    R=mat(R);ratios=[[float(v).as_integer_ratio() for v in row] for row in R];D=max(d for row in ratios for n,d in row)
    totals=[sum(abs(n)*(D//d) for n,d in row) for row in ratios]
    return Q(max(totals),D)

def positive_m(Z,N):
    Z,N=rational(Z),rational(N)
    if not 0<=Z<1 or N<=0:raise ValueError('full-tube Z<1 and positive norm required')
    return (1-Z)/N

def exact_P(R):
    R=mat(R);n=len(R)
    # Entries added AFTER exact binary conversion; never fl64(R+R.T)/2.
    return arb_mat(n,n,[(A(R[i,j])+A(R[j,i]))/arb(2) for i in range(n) for j in range(n)])

def congruence_bounds(R,K,precision):
    ctx.prec=precision;R=mat(R);K=mat(K,len(R));n=len(R);P=exact_P(R);T=arb_mat(n,n,[A(v) for v in K.ravel()])
    # Different association from the producer's (K*P)*K.T: independently P*K.T first.
    B=T*(P*T.transpose());rows=[]
    for i in range(n):
        s=arb(0)
        for j in range(n):s+=(arb(int(i==j))-B[i,j]).abs_upper()
        rows.append(Q(str(s.upper().fmpq())))
    return max(rows),rows

def guards():
    ctx.prec=192;checks={};I=np.eye(2,dtype=np.float64)
    q,_=congruence_bounds(np.array([[2.,5.],[-5.,2.]]),.5*I,192);checks['nonsymmetric_R_exact_SPDpart_certified']=q<1
    q,_=congruence_bounds(np.diag([1.,-1.]),I,192);checks['indefinite_part_not_certified']=q>=1
    q,_=congruence_bounds(np.array([[1.,4.],[0.,1.]]),I,192);checks['positive_eigenvalues_not_enough']=q>=1
    q,_=congruence_bounds(I,np.diag([1.,0.]),192);checks['singular_K_boundary_rejected']=q==1
    q,_=congruence_bounds(np.diag([2.**-1000,1.]),np.diag([2.**500,1.]),192);checks['tiny_exact_binary_SPD_no_float_eigenclaim']=q==0
    tricky=np.array([[2.,1.],[2.**-53,2.]]);P=exact_P(tricky);exact=Q(1,2)+Q(1,2**54)
    checks['P_offdiag_formed_exactly_before_rounding']=P[0,1].contains(A(exact)) and not P[0,1].contains(A(Q.from_float(float(((tricky+tricky.T)/2)[0,1]))))
    checks['independent_integer_dyadic_Rnorm']=norm_R_integer(np.array([[.5,-.25],[1.,0.]]))==1
    checks['exact_rational_coercivity']=positive_m(Q(1,4),Q(3,2))==Q(1,2)
    checks['analytic_H_diag_coercivity_case']=positive_m(0,norm_R_integer(np.diag([.5,.25])))==2
    for name,v in [('negativeZ',Q(-1,10)),('boundaryZ',1),('aboveZ',Q(11,10)),('boolZ',True)]:
        try:positive_m(v,1);checks['reject_'+name]=False
        except ValueError:checks['reject_'+name]=True
    for name,v in [('boolmatrix',np.eye(2,dtype=bool)),('NaN',np.array([[np.nan]])),('complex',np.eye(2,dtype=complex)),('nonsquare',np.zeros((2,3)))]:
        try:mat(v);checks['reject_'+name]=False
        except ValueError:checks['reject_'+name]=True
    try:positive_m(0,0);checks['reject_zero_norm']=False
    except ValueError:checks['reject_zero_norm']=True
    files={}
    try:bind(files,rel(__file__),'0'*64);checks['reject_stale_expected_hash_before_merge']=False
    except ValueError:checks['reject_stale_expected_hash_before_merge']=not files
    files={rel(__file__):'f'*64}
    try:bind(files,rel(__file__),sha(__file__));checks['reject_conflicting_existing_expected_hash']=False
    except ValueError:checks['reject_conflicting_existing_expected_hash']=files[rel(__file__)]=='f'*64
    return checks

def bind(files,name,expected):
    # Check ORIGINAL expected digest BEFORE merge; never overwrite it with currentbytes.
    if sha(path(name))!=expected:raise ValueError('stale original binding '+name)
    if name in files and files[name]!=expected:raise ValueError('conflicting original binding '+name)
    files[name]=expected

def selected(prov,joined,files):
    rid=prov['record_id'];rootfile=prov['root_certificate_file'];bind(files,rootfile,prov['root_certificate_sha256']);root=read(path(rootfile))
    reviewfile=prov['root_review_file'];bind(files,reviewfile,prov['root_review_sha256']);bundle=read(path(reviewfile));matches=[a for a in bundle['record_audits'] if a['record_id']==rid]
    if len(matches)!=1:raise ValueError('one independent capsule perrecord')
    wrapper=matches[0];capsule=wrapper['audit'];bind(files,wrapper['review_file'],wrapper['review_sha256'])
    if read(path(wrapper['review_file']))!=capsule or capsule['status']!='PASS' or not yes(capsule['checks']) or capsule['uniform_quantifier_verified'] is not True:raise ValueError('specific root capsule not literal PASS')
    expected={'certificate_sha256':prov['root_certificate_sha256'],'parameter_cell_sha256':joined['parameter_cell_sha256'],'model_observation_contract_sha256':joined['model_observation_contract_sha256'],'branch_domain_id':root['branch_domain_id'],'branch_law_semantics':'independent_stationary_root_per_record','branch_scope':'selected_unique_root_in_declared_domain'}
    if any(capsule[k]!=v for k,v in expected.items()):raise ValueError('capsule source/cell/domain mismatch')
    if root['status']!='CERTIFIED_UNIFORM_LOCAL_ROOT' or root['parameter_cell']!=joined['parameter_cell'] or root['model_observation_contract_sha256']!=joined['model_observation_contract_sha256']:raise ValueError('root family differs')
    for name,digest in root['artifact_hashes'].items():bind(files,name,digest)
    if rid.endswith('_M'):
        bind(files,capsule['reuse_review_file'],capsule['reuse_review_sha256']);nested=read(path(capsule['reuse_review_file']));audit=nested['root_audit']
        if audit['status']!='PASS' or not yes(audit['checks']) or audit['uniform_quantifier_verified'] is not True:raise ValueError('nested M independent root audit')
        oldfile=root['reuse_bindings']['root_certificate'];bind(files,oldfile,audit['certificate_sha256']);oldroot=read(path(oldfile))
        if oldroot['uniform_cell_result']!=root['uniform_cell_result'] or audit['branch_domain_id']!=root['branch_domain_id'] or audit['parameter_cell_sha256']!=joined['parameter_cell_sha256']:raise ValueError('M mechanical tube not identical')
        for name,digest in nested['reviewed_hashes'].items():bind(files,name,digest)
        bounds=nested['independent_bounds'];Z=Q(bounds['Z_uniform']['upper_exact']);source={'file':capsule['reuse_review_file'],'object_path':'independent_bounds.Z_uniform'}
    else:
        bounds=capsule['bounds'];Z=Q(bounds['Z']['upper_exact']);source={'file':wrapper['review_file'],'object_path':'bounds.Z'}
    r=Q(root['uniform_cell_result']['selected_radius_exact']);j=Q(bounds['quadrature_J_lower']['lower_exact'])
    if not 0<=Z<1 or r<=0 or not Q(bounds['radius']['lower_exact'])<=r<=Q(bounds['radius']['upper_exact']) or Q(bounds['eta_plus_Zr']['upper_exact'])>=r or j<=Q(1,10**6):raise ValueError('full-tube strict root/domain premise not bound')
    cachefile=root['point_cache_file']
    if cachefile not in files:raise ValueError('cache MUST have immutable artifact/review expecteddigest before currentread')
    cache=read(path(cachefile));arrays=decode(cache);R=mat(arrays['preconditioner'],611)
    return {'record_id':rid,'root_certificate_file':rootfile,'root_certificate_sha256':prov['root_certificate_sha256'],'root_capsule_file':wrapper['review_file'],'root_capsule_sha256':wrapper['review_sha256'],'point_cache_file':cachefile,'point_cache_sha256':files[cachefile],'frozen_binary_cache_sha256':cache['frozen_binary_cache']['sha256'],'R_binary_sha256':hashlib.sha256(np.ascontiguousarray(R).tobytes()).hexdigest(),'independent_Z_source':source,'independent_full_tube_Z_upper_exact':str(Z),'exact_selected_corrector_radius':str(r),'quadrature_J_lower_exact':str(j),'branch_domain_id':root['branch_domain_id'],'parameter_cell_sha256':joined['parameter_cell_sha256'],'model_observation_contract_sha256':joined['model_observation_contract_sha256']}

def prepare():
    if (RUN/'protocol.json').exists():raise RuntimeError('preserve frozen independentprotocol')
    terminal=read(TARGET/'terminal.json');producer=read(TARGET/'protocol.json');result=read(TARGET/'results.json');joined=read(JOIN);theory=read(THEORY)
    if terminal['status']!='FINISHED' or terminal['results_sha256']!=sha(TARGET/'results.json'):raise ValueError('producer true terminal required')
    if not yes(theory['checks']):raise ValueError('independent theory applicability review required')
    files={}
    for name,digest in producer['source_hashes'].items():bind(files,name,digest)
    for name,digest in producer['input_hashes'].items():bind(files,name,digest)
    for q in [TARGET/'protocol.json',TARGET/'results.json',TARGET/'terminal.json',TARGET/'prelaunch_independent_review.json',JOIN,THEORY]:bind(files,rel(q),sha(q))
    records=[];outputs={r['record_id']:r for r in result['records']};prow={r['record_id']:r for r in producer['records']};trows={r['record_id']:r for r in theory['selected_records']}
    for prov in joined['proof_provenance']:
        row=selected(prov,joined,files);entry=outputs[row['record_id']];bind(files,entry['file'],entry['sha256']);answer=read(path(entry['file']));reported=prow[row['record_id']]
        if answer['record_id']!=row['record_id'] or answer['root_binding']!=reported:raise ValueError('producer recordbinding mismatch')
        for key in ['root_certificate_file','root_certificate_sha256','point_cache_file','point_cache_sha256','frozen_binary_cache_sha256','independent_full_tube_Z_upper_exact','exact_selected_corrector_radius','quadrature_J_lower_exact','branch_domain_id','parameter_cell_sha256','model_observation_contract_sha256']:
            if reported[key]!=row[key]:raise ValueError('independently rebuilt producercontext mismatch '+key)
        if reported['preconditioner_binary_sha256']!=row['R_binary_sha256']:raise ValueError('actualRbytes different')
        tr=trows[row['record_id']]
        if tr['certificate_sha256']!=row['root_certificate_sha256'] or Q(tr['uniform_Z_upper_exact'])!=Q(row['independent_full_tube_Z_upper_exact']) or tr['branch_domain_id']!=row['branch_domain_id']:raise ValueError('independent theory/root premise differs')
        row['producer_result_file']=entry['file'];row['producer_result_sha256']=entry['sha256']
        if answer['status']=='CERTIFIED_UNIFORM_STRICT_LOCAL_ENERGY_MINIMUM':
            bind(files,answer['K_proposal_file'],answer['K_proposal_sha256']);row['K_proposal_file']=answer['K_proposal_file'];row['K_proposal_sha256']=answer['K_proposal_sha256']
            with np.load(path(row['K_proposal_file']),allow_pickle=False) as z:mat(z['K'],611)
        else:row['producer_unresolved_status']=answer['status']
        records.append(row)
    checks={'all17_selected_records':len(records)==17 and len({r['record_id'] for r in records})==17,'all_current_original_expected_digests_validated_before_merge':True,'all_independent_fulltube_Z_radius_C_cache_bindings':True,'M_nested_source_handled':any(r['record_id'].endswith('_M') and r['independent_Z_source']['object_path']=='independent_bounds.Z_uniform' for r in records),'all_small_analytic_guards_PASS':yes(guards()),'producer_true_terminal':True,'no_actual611_congruence_before_review':not (RUN/'records').exists()}
    write(RUN/'protocol.json',{'review_source_sha256':sha(__file__),'precision_bits':192,'file_hashes':files,'producer_protocol_sha256':sha(TARGET/'protocol.json'),'records':records,'python_flint_version':importlib.metadata.version('python-flint'),'parameter_cell_sha256':joined['parameter_cell_sha256'],'model_observation_contract_sha256':joined['model_observation_contract_sha256'],'method':'Independent exact-entry P construction, different-association Arb192 K*(P*K.T), exactdyadic integer rownorm, rational m. No producer numericalcertifier/no Hassembly.','scope':'Uniform discrete strictlocalenergyminima only, notdynamicalstability/commonpath/globalminimum/continuum/sourcevalidation ornovelty.'})
    write(RUN/'guard_tests.json',{'status':'PASS','checks':guards(),'source_sha256':sha(__file__),'no_actual611_work':True})
    write(RUN/'preflight.json',{'status':'PASS','checks':checks,'review_source_sha256':sha(__file__),'protocol_sha256':sha(RUN/'protocol.json'),'guard_tests_sha256':sha(RUN/'guard_tests.json'),'actual611_congruences_computed':0})
    print(json.dumps({'status':'PREFLIGHT_PASS_ONLY','source_sha256':sha(__file__),'protocol_sha256':sha(RUN/'protocol.json'),'preflight_sha256':sha(RUN/'preflight.json')}),flush=True)

def frozen(p):
    if sha(__file__)!=p['review_source_sha256'] or importlib.metadata.version('python-flint')!=p['python_flint_version']:raise ValueError('review source/backend changed')
    for name,digest in p['file_hashes'].items():
        if sha(path(name))!=digest:raise ValueError('frozen input changed '+name)
    return True

def run():
    p=read(RUN/'protocol.json');gate=read(RUN/'prelaunch_review.json')
    if gate.get('status')!='PASS' or not yes(gate.get('checks')) or gate.get('review_source_sha256')!=sha(__file__) or gate.get('protocol_sha256')!=sha(RUN/'protocol.json') or gate.get('preflight_sha256')!=sha(RUN/'preflight.json'):raise ValueError('concrete independent review prelaunch required')
    if OUT.exists():raise RuntimeError('preserve finished stabilityreview')
    start=time.time();records=[]
    for row in p['records']:
        frozen(p);beg=time.time();producer=read(path(row['producer_result_file']))
        if 'K_proposal_file' not in row:
            audit={'record_id':row['record_id'],'status':'UNRESOLVED','checks':{'producer_has_frozen_K':False},'producer_result_sha256':row['producer_result_sha256']}
        else:
            arrays=decode(read(path(row['point_cache_file'])));R=mat(arrays['preconditioner'],611)
            with np.load(path(row['K_proposal_file']),allow_pickle=False) as z:K=mat(z['K'],611)
            q,rows=congruence_bounds(R,K,p['precision_bits']);N=norm_R_integer(R);Z=Q(row['independent_full_tube_Z_upper_exact']);m=positive_m(Z,N)
            checks={'dimension611_exact_binary_K_R':R.shape==K.shape==(611,611),'actualRbytes_match_reviewed_rootcache':hashlib.sha256(np.ascontiguousarray(R).tobytes()).hexdigest()==row['R_binary_sha256'],'independent_exact_sympart_congruence_norm_below1':q<1,'exact_Rnorm_reproduces_producer_rational':N==Q(producer['exact_R_infinity_norm']),'exact_m_reproduces_positive_rational_bound':m>0 and m==Q(producer['uniform_free_coordinate_Hessian_lower_bound_exact']),'same_independent_fulltube_Z':Z==Q(producer['independent_full_tube_Z_upper_exact']),'all_current_frozen_inputs_still_match':frozen(p),'no_producer_certifier_or_Hessian_assembly_called':True}
            audit={'record_id':row['record_id'],'status':'PASS' if yes(checks) else 'UNRESOLVED','checks':checks,'producer_result_sha256':row['producer_result_sha256'],'root_certificate_sha256':row['root_certificate_sha256'],'parameter_cell_sha256':row['parameter_cell_sha256'],'model_observation_contract_sha256':row['model_observation_contract_sha256'],'branch_domain_id':row['branch_domain_id'],'precision_bits':p['precision_bits'],'independent_q_infinity_upper_exact':str(q),'all611_row_upper_bounds_exact':list(map(str,rows)),'exact_R_infinity_norm':str(N),'independent_fulltube_Z_upper_exact':str(Z),'uniform_Hessian_lower_bound_exact':str(m),'normal_root_radius_exact':row['exact_selected_corrector_radius'],'K_proposal_sha256':row['K_proposal_sha256'],'elapsed_seconds':time.time()-beg}
        file=RUN/'records'/(row['record_id']+'.json');write(file,audit);records.append({'record_id':row['record_id'],'review_file':rel(file),'review_sha256':sha(file),'audit':audit});print(json.dumps({'record':row['record_id'],'status':audit['status'],'seconds':time.time()-start}),flush=True)
    checks={'all17_independently_reviewed':len(records)==17,'all17_exact_sympart_and_uniform_m_PASS':all(r['audit']['status']=='PASS' for r in records),'all_frozen_bindings_current':frozen(p),'no_Hessian_or_equilibrium_jobs':True,'classical_lemma_and_exact_declared_tube_scope':True}
    write(OUT,{'status':'PASS' if yes(checks) else 'PARTIAL_UNRESOLVED','checks':checks,'record_audits':records,'review_source_sha256':sha(__file__),'protocol_sha256':sha(RUN/'protocol.json'),'preflight_sha256':sha(RUN/'preflight.json'),'prelaunch_review_sha256':sha(RUN/'prelaunch_review.json'),'parameter_cell_sha256':p['parameter_cell_sha256'],'model_observation_contract_sha256':p['model_observation_contract_sha256'],'precision_bits':p['precision_bits'],'scope':p['scope'],'independent_implementation':'Own exact-entry P, different-association intervalmatrix product and integer dyadic norm. FrozenK is onlya proposal; no producer certify_sympart call. Existingwhole-tube root/Hessian bounds and independentlyreviewed scalar-energy applicability reused.','original_five_gate_goal_complete':False,'confirmed_novelty':False,'elapsed_seconds':time.time()-start});print(json.dumps({'status':read(OUT)['status'],'receipt_sha256':sha(OUT),'seconds':time.time()-start}),flush=True)
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['guards','prepare','run']);args=parser.parse_args()
    if args.stage=='guards':print(json.dumps({'status':'PASS' if yes(guards()) else 'FAIL','checks':guards(),'source_sha256':sha(__file__)}))
    elif args.stage=='prepare':prepare()
    else:run()
