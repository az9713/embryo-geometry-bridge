"""Certify sym(R)>0, then uniform local energy minima for17 frozen roots.

No equilibrium, Hessian assembly, parameter fitting, or physical target edits.
The exact target is P_ij=(binary_Rij+binary_Rji)/2, NEVER rounded float sym(R).
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,base64,datetime,hashlib,importlib.metadata,io,json,os,time,traceback
import numpy as np
import scipy
from scipy.linalg import cholesky,solve_triangular
from flint import arb,arb_mat,fmpq,ctx

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'runs/cohort_energy_stability_v1'
JOIN=ROOT/'evidence/shared_cell_cohort_inference.json'
JOINPROTO=ROOT/'runs/shared_cell_cohort_v1/protocol.json'
THEORY=ROOT/'evidence/cohort_energy_stability_theory_review.json'
IDS={'optoGEF2_small_col_'+chr(i) for i in range(ord('B'),ord('R')+1)}

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf8'))
def rel(p):return Path(p).resolve().relative_to(ROOT.resolve()).as_posix()
def path(name):
    p=(ROOT/name).resolve();p.relative_to(ROOT.resolve());return p
def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def q(x):
    if isinstance(x,(bool,np.bool_)):raise ValueError('boolean is not a bound')
    return Q.from_float(float(x)) if isinstance(x,(float,np.floating)) else Q(x)
def A(x):
    v=q(x);return arb(fmpq(v.numerator,v.denominator))
def pack(x):return {'ball':str(x),'lower_exact':str(x.lower().fmpq()),'upper_exact':str(x.upper().fmpq())}
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);data=json.dumps(v,indent=2,allow_nan=False)+'\n'
    if p.exists():
        if p.read_text(encoding='utf8')!=data:raise ValueError('immutable output already exists: '+str(p))
    else:p.write_text(data,encoding='utf8')
def check_hashes(bindings):
    for name,digest in bindings.items():
        if sha(path(name))!=digest:raise ValueError('bound input changed: '+name)
def merge_expected(destination,incoming):
    """Never let a current-byte snapshot replace authoritative provenance."""
    for name,digest in incoming.items():
        name=rel(path(name))
        if not isinstance(digest,str) or len(digest)!=64 or any(c not in '0123456789abcdef' for c in digest):raise ValueError('literal SHA256 digest required')
        if name in destination and destination[name]!=digest:raise ValueError('conflicting expected hashes: '+name)
        destination[name]=digest
    return destination
def combined_hashes(*maps):
    answer={}
    for bindings in maps:merge_expected(answer,bindings)
    return answer
def verify_expected_file(name,digest,bindings):
    if sha(path(name))!=digest:raise ValueError('authoritative file hash differs BEFORE decoding: '+name)
    merge_expected(bindings,{name:digest})
def require_checks(c):
    if not isinstance(c,dict) or not c or any(v is not True for v in c.values()):raise ValueError('nonempty literal-True checks required')
def matrix(value,n=None):
    a=np.asarray(value)
    if a.dtype!=np.dtype('float64') or a.ndim!=2 or a.shape[0]!=a.shape[1] or not a.shape[0] or (n is not None and a.shape!=(n,n)) or not np.isfinite(a).all():raise ValueError('finite square exact-binary64 matrix required')
    return a
def decode(packet):
    blob=packet['frozen_binary_cache'];raw=base64.b64decode(blob['data'],validate=True)
    if hashlib.sha256(raw).hexdigest()!=blob['sha256']:raise ValueError('frozen cache checksum differs')
    with np.load(io.BytesIO(raw),allow_pickle=False) as z:arrays={k:z[k].copy() for k in z.files}
    if any(not np.isfinite(a).all() for a in arrays.values()):raise ValueError('nonfinite cache array')
    matrix(arrays['preconditioner'],611)
    if arrays['u'].shape!=(611,) or arrays['free'].shape!=(611,) or len(set(map(int,arrays['free'])))!=611:raise ValueError('611 independent linear free coordinates required')
    return arrays

def upper_Z(value):
    z=q(value)
    if not 0<=z<1:raise ValueError('independently audited FULL-TUBE upper Z must be in[0,1)')
    return z

def certify_sympart(R,K=None,precision_bits=160):
    """Exact-binary target, rigorous Arb residual; floatCholesky is only proposal."""
    R=matrix(R);n=len(R);ctx.prec=precision_bits
    if K is None:
        try:
            with np.errstate(over='raise',invalid='raise',divide='raise'):
                P_float=(R+R.T)/2
                L=cholesky(P_float,lower=True,check_finite=True)
                K=solve_triangular(L,np.eye(n,dtype=np.float64),lower=True,check_finite=True)
            matrix(K,n)
        except (np.linalg.LinAlgError,ValueError,FloatingPointError) as error:
            return {'status':'UNRESOLVED_FLOATING_PROPOSAL','reason':repr(error),'positive_or_negative_eigenvalue_inference':False},None
    else:K=matrix(K,n)
    Rb=arb_mat(n,n,[A(v) for v in R.ravel()]);Kb=arb_mat(n,n,[A(v) for v in K.ravel()])
    # Addition/division happen AFTER exact binary entries enter Arb. P_float
    # above never defines the certified matrix.
    Pb=(Rb+Rb.transpose())/arb(2)
    B=Kb*Pb*Kb.transpose();row_bounds=[]
    for i in range(n):
        row=sum(((arb(int(i==j))-B[i,j]).abs_upper() for j in range(n)),arb(0))
        row_bounds.append(Q(str(row.upper().fmpq())))
    bound=max(row_bounds)
    result={'status':'CERTIFIED_EXACT_SYMMETRIC_PART_SPD' if bound<1 else 'UNRESOLVED_SYMMETRIC_PART_SUFFICIENT_BOUND',
        'exact_target':'Pij=(exact_binary_Rij+exact_binary_Rji)/2','precision_bits':precision_bits,'dimension':n,
        'q_infinity_upper_exact':str(bound),'row_absolute_sum_upper_exact':list(map(str,row_bounds)),
        'criterion':'||I-K P K^T||infinity<1; Pexact symmetric =>K P K^T SPD =>squareK invertible =>P SPD',
        'floating_Cholesky_and_K_are_untrusted_proposals_only':True,'R_symmetry_assumed':False,
        'R_binary_sha256':hashlib.sha256(np.ascontiguousarray(R).tobytes()).hexdigest()}
    return result,K

def exact_R_norm(R):
    R=matrix(R)
    return max(sum((abs(q(v)) for v in row),Q(0)) for row in R)

def validate_capsule(capsule,provenance,root,cell,contract):
    if capsule.get('status')!='PASS':raise ValueError('specific independent root capsule must PASS')
    require_checks(capsule.get('checks'))
    expected={'certificate_sha256':provenance['root_certificate_sha256'],'parameter_cell_sha256':cell,
        'model_observation_contract_sha256':contract,'branch_domain_id':root['branch_domain_id'],
        'branch_law_semantics':'independent_stationary_root_per_record','branch_scope':'selected_unique_root_in_declared_domain','uniform_quantifier_verified':True}
    if any(capsule.get(k)!=v for k,v in expected.items()):raise ValueError('independent capsule root/cell/domain/quantifier mismatch')

def selected_record(provenance,joined,inputs):
    file=provenance['root_certificate_file'];root=read(path(file))
    if sha(path(file))!=provenance['root_certificate_sha256'] or root['status']!='CERTIFIED_UNIFORM_LOCAL_ROOT':raise ValueError('selected root bytes/status differ')
    reviewfile=provenance['root_review_file'];bundle=read(path(reviewfile))
    if sha(path(reviewfile))!=provenance['root_review_sha256']:raise ValueError('selected independent review changed')
    wrapper=next(a for a in bundle['record_audits'] if a['record_id']==provenance['record_id']);capsule=wrapper['audit']
    capsulefile=wrapper['review_file']
    if sha(path(capsulefile))!=wrapper['review_sha256'] or read(path(capsulefile))!=capsule:raise ValueError('independent capsule byte association differs')
    validate_capsule(capsule,provenance,root,joined['parameter_cell_sha256'],joined['model_observation_contract_sha256'])
    proof=root['uniform_cell_result'];radius=q(proof['selected_radius_exact'])
    if radius<=0 or proof['status']!='CERTIFIED_UNIFORM_LOCAL_ROOT' or proof['cell_bounds_exact']!=joined['parameter_cell']:raise ValueError('selected full root tube differs')
    merge_expected(inputs,{file:provenance['root_certificate_sha256'],reviewfile:provenance['root_review_sha256'],capsulefile:wrapper['review_sha256']})
    merge_expected(inputs,root.get('artifact_hashes',{}))
    if provenance['record_id'].endswith('_M'):
        nestedfile=capsule['reuse_review_file'];nested=read(path(nestedfile))
        if sha(path(nestedfile))!=capsule['reuse_review_sha256']:raise ValueError('immutable M nested review changed')
        audit=nested['root_audit'];require_checks(audit['checks'])
        if audit['status']!='PASS' or audit['uniform_quantifier_verified'] is not True or audit['parameter_cell_sha256']!=joined['parameter_cell_sha256'] or audit['branch_domain_id']!=root['branch_domain_id']:raise ValueError('immutable M root audit scope mismatch')
        oldfile=root['reuse_bindings']['root_certificate'];oldroot=read(path(oldfile))
        if sha(path(oldfile))!=audit['certificate_sha256'] or oldroot['uniform_cell_result']!=proof or q(oldroot['uniform_cell_result']['selected_radius_exact'])!=radius:raise ValueError('M exact mechanical tube equivalence failed')
        bounds=nested['independent_bounds'];zpack=bounds['Z_uniform'];rpack=bounds['radius'];jpack=bounds['quadrature_J_lower']
        if audit['checks'].get('independent_analytic_domain_J_positive') is not True:raise ValueError('M analytic-domain review missing')
        merge_expected(inputs,{nestedfile:capsule['reuse_review_sha256'],oldfile:audit['certificate_sha256']})
        merge_expected(inputs,nested['reviewed_hashes'])
        bound_source={'file':nestedfile,'object_path':'independent_bounds.Z_uniform','root_audit_object_path':'root_audit'}
    else:
        bounds=capsule['bounds'];zpack=bounds['Z'];rpack=bounds['radius'];jpack=bounds['quadrature_J_lower']
        if capsule['checks'].get('all_analytic_domain_J_gt_1e6_inverse') is not True:raise ValueError('independent analytic-domain review missing')
        bound_source={'file':capsulefile,'object_path':'bounds.Z','root_audit_object_path':'$'}
    Z=upper_Z(zpack['upper_exact'])
    if not q(rpack['lower_exact'])<=radius<=q(rpack['upper_exact']):raise ValueError('independent radius enclosure does not cover exact selected radius')
    if q(jpack['lower_exact'])<=Q(1,1000000):raise ValueError('scalar C2 penalty-inactive analytic domain not certified')
    if q(bounds['eta_plus_Zr']['upper_exact'])>=radius:raise ValueError('independent full tube strict self-map missing')
    cachefile=root['point_cache_file']
    if cachefile not in root.get('artifact_hashes',{}):raise ValueError('authoritative selected-root cache hash required')
    cache_expected=root['artifact_hashes'][cachefile]
    verify_expected_file(cachefile,cache_expected,inputs)
    packet=read(path(cachefile));arrays=decode(packet)
    sources={'scripts/'+name:digest for name,digest in root['source_hashes'].items()}
    row={'record_id':provenance['record_id'],'root_certificate_file':file,'root_certificate_sha256':sha(path(file)),
        'independent_root_capsule_file':capsulefile,'independent_root_capsule_sha256':sha(path(capsulefile)),
        'independent_Z_source':bound_source,'independent_full_tube_Z_upper_exact':str(Z),
        'exact_selected_corrector_radius':str(radius),'independent_radius_enclosure':rpack,
        'quadrature_J_lower_exact':jpack['lower_exact'],'point_cache_file':cachefile,'point_cache_sha256':cache_expected,
        'frozen_binary_cache_sha256':packet['frozen_binary_cache']['sha256'],
        'preconditioner_binary_sha256':hashlib.sha256(np.ascontiguousarray(arrays['preconditioner']).tobytes()).hexdigest(),
        'dimension':611,'branch_domain_id':root['branch_domain_id'],'parameter_cell_sha256':joined['parameter_cell_sha256'],
        'model_observation_contract_sha256':joined['model_observation_contract_sha256'],
        'C2_scalar_energy_and_linear_free_constraints':'Independent analytic-domain/root proof plus separately audited theory applicability; no new H assembly.',
        'output_folder':rel(OUT/'records'/provenance['record_id'])}
    return row,sources

def guards():
    checks={};I=np.eye(2,dtype=np.float64)
    cert,_=certify_sympart(np.array([[2.,5.],[-5.,2.]]));checks['nonsymmetric_R_SPD_part_accepted']=cert['status']=='CERTIFIED_EXACT_SYMMETRIC_PART_SPD'
    for label,R in [('symmetric_indefinite',np.diag([1.,-1.])),('positive_eigenvalues_but_indefinite_sympart',np.array([[1.,4.],[0.,1.]]))]:
        c,_=certify_sympart(R,I);checks[label+'_not_certified']=c['status']!='CERTIFIED_EXACT_SYMMETRIC_PART_SPD'
    c,_=certify_sympart(np.diag([2.**-1020,1.]));checks['near_singular_exact_binary_SPD_accepted']=c['status']=='CERTIFIED_EXACT_SYMMETRIC_PART_SPD'
    c,_=certify_sympart(I,np.diag([1.,0.]));checks['rank_deficient_K_q_boundary_rejected']=q(c['q_infinity_upper_exact'])==1 and c['status']!='CERTIFIED_EXACT_SYMMETRIC_PART_SPD'
    for label,p,accepted in [('below',2.**-52,True),('equal',0.,False),('above',-2.**-52,False)]:
        c,_=certify_sympart(np.diag([p,1.]),I);checks['exact_q_'+label+'_one']=(c['status']=='CERTIFIED_EXACT_SYMMETRIC_PART_SPD') is accepted
    c,_=certify_sympart(np.diag([1.,-1.]));checks['failed_Cholesky_only_unresolved_no_inertia_claim']=c['status']=='UNRESOLVED_FLOATING_PROPOSAL' and c['positive_or_negative_eigenvalue_inference'] is False
    checks['exact_R_infinity_norm']=exact_R_norm(np.array([[.5,-.25],[1.,0.]]))==1
    for label,bad in [('NaN',np.array([[np.nan]])),('nonsquare',np.zeros((2,3))),('bool',np.eye(2,dtype=bool)),('complex',np.eye(2,dtype=complex))]:
        try:matrix(bad);checks['reject_'+label]=False
        except ValueError:checks['reject_'+label]=True
    for label,z in [('one','1'),('above','1001/1000'),('negative','-1/1000'),('boolean',True)]:
        try:upper_Z(z);checks['reject_independent_full_Z_'+label]=False
        except ValueError:checks['reject_independent_full_Z_'+label]=True
    h0=hashlib.sha256(b'authoritative original cache').hexdigest();h1=hashlib.sha256(b'tampered cache').hexdigest()
    fixture=OUT/'guard_fixtures/tampered_cache.json';fixture.parent.mkdir(parents=True,exist_ok=True)
    if fixture.exists() and fixture.read_bytes()!=b'tampered cache':raise ValueError('guard fixture changed')
    if not fixture.exists():fixture.write_bytes(b'tampered cache')
    try:verify_expected_file(rel(fixture),h0,{});checks['tampered_cache_rejected_before_decode']=False
    except ValueError:checks['tampered_cache_rejected_before_decode']=True
    try:merge_expected({rel(fixture):h0},{rel(fixture):h1});checks['current_snapshot_cannot_override_authoritative_hash']=False
    except ValueError:checks['current_snapshot_cannot_override_authoritative_hash']=True
    try:combined_hashes({rel(fixture):h0},{rel(fixture):h1});checks['source_input_binding_conflict_rejected']=False
    except ValueError:checks['source_input_binding_conflict_rejected']=True
    checks['identical_duplicate_binding_allowed']=combined_hashes({rel(fixture):h0},{rel(fixture):h0})=={rel(fixture):h0}
    result={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'source_sha256':sha(__file__),
        'scope':'Small analytical/adversarial matrices only; no actual611-dimensional certificate calculated.'}
    write(OUT/'guard_tests.json',result)
    if not all(checks.values()):raise ValueError('stability guards failed')
    return result

def versions():return {'numpy':np.__version__,'scipy':scipy.__version__,'python-flint':importlib.metadata.version('python-flint')}

def freeze():
    test=read(OUT/'guard_tests.json');require_checks(test['checks'])
    if test['status']!='PASS' or test['source_sha256']!=sha(__file__):raise ValueError('current-source passing guards required')
    joined=read(JOIN);prior=read(JOINPROTO);theory=read(THEORY);require_checks(theory['checks'])
    if sha(JOIN)!= '1f8c4dc1ad8924f8fda13e6748d5e73aac426b5771b02404803d27c1955aa677':raise ValueError('selected cohort changed')
    if len(joined['proof_provenance'])!=17 or {r['record_id'] for r in joined['proof_provenance']}!=IDS:raise ValueError('exact17 selected root families required')
    inputs=combined_hashes(prior['file_hashes'],theory['input_hashes'])
    for f in (JOIN,JOINPROTO,THEORY,OUT/'guard_tests.json',ROOT/'evidence/shared_cell_cohort_independent_review.json'):merge_expected(inputs,{rel(f):sha(f)})
    records=[];sources={rel(__file__):sha(__file__)}
    for provenance in joined['proof_provenance']:
        row,local_sources=selected_record(provenance,joined,inputs);merge_expected(sources,local_sources)
        t=next(x for x in theory['selected_records'] if x['record_id']==row['record_id'])
        if t['certificate_sha256']!=row['root_certificate_sha256'] or t['branch_domain_id']!=row['branch_domain_id'] or q(t['uniform_Z_upper_exact'])!=q(row['independent_full_tube_Z_upper_exact']):raise ValueError('separate theory/source applicability differs')
        records.append(row)
    check_hashes(combined_hashes(sources,inputs))
    protocol={'version':1,'created_utc':stamp(),'source_hashes':sources,'input_hashes':inputs,'environment_versions':versions(),
        'records':records,'precision_bits':160,'parameter_cell':joined['parameter_cell'],'parameter_cell_sha256':joined['parameter_cell_sha256'],
        'model_observation_contract_sha256':joined['model_observation_contract_sha256'],
        'P_target':'P_ij=(exact_binary_Rij+exact_binary_Rji)/2 formed inArb; neverroundedfloating symP',
        'K_policy':'Floating scipy lower Cholesky of fl64((R+R.T)/2), solve L*K=I; K proposal then frozen binaryNPZ. No K invertibility assumption.',
        'SPD_gate':'exact rational q_upper for Arb160 ||I-K P K^T||infinity <1',
        'failed_proposal_or_bound':'UNRESOLVED; no sign/eigenvalue/nonexistence inference; retain failure and continue17',
        'lemma':'Independent fulltube ||I-RH||infinity<=Z<1 gives RH positive-stable. Exactsym(R)>0 and symmetricC2 Hessian imply H>0 by Lyapunov identity A^T H+H A=2H sym(R)H, A=RH.',
        'coercivity_formula':'m=(1-independent_uniform_Z_upper)/exact_binary_R_infinity_norm; H^-1=(RH)^-1 R and Neumann bound imply lambda_min(H)>=m>0 after HSPD.',
        'coercivity_units':'Normalized free-coordinate energy curvature only; not a physical relaxation time, force or continuum modulus.',
        'target':'Exact scalar C2 finite binaryquadrature energy on fixed independent linearfreecoordinates; existing reviewed inactivepenalty/positiveanalyticdomain.',
        'root_domains':'Exactly the17 previously selected fullpredictor+NORMALcorrector tubes, including recentered F/K; no new root or H assembly.',
        'scope':'Uniform strict local energy minima if all17 SPD gates pass. No globalminimum, PDE, dynamical/path stability, sourcebiology validation, datafit or novelty claim.',
        'prelaunch_required_fields':['status=PASS','checks all literal true','source_sha256','protocol_sha256','preflight_sha256']}
    write(OUT/'protocol.json',protocol)
    checks={'all17_cache_preconditioners_finite_square611':len(records)==17,'all17_specific_independent_full_Z_bound':all(q(r['independent_full_tube_Z_upper_exact'])<1 for r in records),
        'M_nested_independent_bounds_and_root_audit_bound':any(r['record_id'].endswith('_M') and r['independent_Z_source']['object_path']=='independent_bounds.Z_uniform' for r in records),
        'actual_recentered_F_K_selected':all('cohort_recentered_roots_v1' in r['point_cache_file'] for r in records if r['record_id'].endswith(('_F','_K'))),
        'all_guard_checks_pass':all(test['checks'].values()),'no_actual611_SPD_certificate_calculated_before_review':not (OUT/'records').exists(),
        'all_frozen_sources_inputs_current':True,'C2_scalar_energy_and_linear_free_constraints_theory_review_bound':True}
    write(OUT/'preflight.json',{'status':'PASS','checks':checks,'source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json'),
        'guard_tests_sha256':sha(OUT/'guard_tests.json'),'actual611_certifications':0,'runtime_expectation':'Unmeasured estimate1–10minutes total17 sequential611Arb160 two-matrix products; actual perrecord timings retained. No actual611benchmark before review.'})
    return protocol

def run():
    p=read(OUT/'protocol.json');reviewfile=OUT/'prelaunch_independent_review.json';review=read(reviewfile);require_checks(review.get('checks'))
    if review.get('status')!='PASS' or review.get('source_sha256')!=sha(__file__) or review.get('protocol_sha256')!=sha(OUT/'protocol.json') or review.get('preflight_sha256')!=sha(OUT/'preflight.json'):raise ValueError('matching concrete independent prelaunch review required')
    if p['environment_versions']!=versions():raise ValueError('frozen library versions changed')
    check_hashes(combined_hashes(p['source_hashes'],p['input_hashes']));jobpath=OUT/'job.json'
    if jobpath.exists():raise ValueError('existing job retained; do not restart')
    job={'status':'RUNNING','pid':os.getpid(),'started_utc':stamp(),'source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json'),'prelaunch_review_sha256':sha(reviewfile)}
    write(jobpath,job);start=time.time();summaries=[]
    try:
        for row in p['records']:
            check_hashes(combined_hashes(p['source_hashes'],p['input_hashes']));beg=time.time()
            verify_expected_file(row['point_cache_file'],row['point_cache_sha256'],{})
            arrays=decode(read(path(row['point_cache_file'])));R=arrays['preconditioner']
            certificate,K=certify_sympart(R,precision_bits=p['precision_bits']);folder=path(row['output_folder']);extra={}
            if K is not None:
                folder.mkdir(parents=True,exist_ok=True);buf=io.BytesIO();np.savez_compressed(buf,K=K);proposal=folder/'K_proposal.npz'
                if proposal.exists():raise ValueError('existing proposal retained')
                proposal.write_bytes(buf.getvalue());extra={'K_proposal_file':rel(proposal),'K_proposal_sha256':sha(proposal)}
            answer={'record_id':row['record_id'],'status':'UNRESOLVED','protocol_sha256':sha(OUT/'protocol.json'),'root_binding':row,
                'symmetric_part_certificate':certificate,**extra,'no_new_Hessian_or_equilibrium_assembly':True}
            if certificate['status']=='CERTIFIED_EXACT_SYMMETRIC_PART_SPD':
                norm=exact_R_norm(R);Z=q(row['independent_full_tube_Z_upper_exact']);m=(1-Z)/norm
                if norm<=0 or m<=0:raise ValueError('positive exact coercivity bound required')
                answer.update(status='CERTIFIED_UNIFORM_STRICT_LOCAL_ENERGY_MINIMUM',exact_R_infinity_norm=str(norm),
                    independent_full_tube_Z_upper_exact=str(Z),uniform_free_coordinate_Hessian_lower_bound_exact=str(m),
                    whole_tube_Hessian_positive=True,strict_local_minimum_at_selected_root=True,
                    dynamical_or_common_path_or_globalminimum_or_PDE_claim=False)
            answer['elapsed_seconds']=time.time()-beg;file=folder/'result.json';write(file,answer)
            summaries.append({'record_id':row['record_id'],'status':answer['status'],'file':rel(file),'sha256':sha(file)})
            print(json.dumps({'record_id':row['record_id'],'status':answer['status'],'seconds':answer['elapsed_seconds']}),flush=True)
        check_hashes(combined_hashes(p['source_hashes'],p['input_hashes']));complete=all(r['status']=='CERTIFIED_UNIFORM_STRICT_LOCAL_ENERGY_MINIMUM' for r in summaries)
        result={'status':'CERTIFIED_ALL17_UNIFORM_LOCAL_ENERGY_MINIMA' if complete else 'PARTIAL_UNRESOLVED','records':summaries,
            'protocol_sha256':sha(OUT/'protocol.json'),'parameter_cell_sha256':p['parameter_cell_sha256'],
            'elapsed_seconds':time.time()-start,'independent_mathematical_replay_required':True,'scope':p['scope']}
        write(OUT/'results.json',result)
        # A separate terminal receipt preserves the immutable launch receipt.
        write(OUT/'terminal.json',{'status':'FINISHED','pid':os.getpid(),'finished_utc':stamp(),'results_sha256':sha(OUT/'results.json')})
        return result
    except BaseException:
        write(OUT/'terminal.json',{'status':'FAILED','pid':os.getpid(),'failed_utc':stamp(),'traceback':traceback.format_exc()});raise

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['test','freeze','run']);args=parser.parse_args()
    result=guards() if args.stage=='test' else freeze() if args.stage=='freeze' else run()
    print(json.dumps({'stage':args.stage,'status':result.get('status','FROZEN')}))
