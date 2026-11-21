"""Independent cohort uniform-root replay; saved states only, no equilibrium solve.

Uses reviewed exact-cache local strain formulas, but separately assembles all
forces, Hessians, five parameter columns, affine Taylor remainders and uniform
Jacobian bounds. Never imports or calls the producer's build_point/evaluate_cell.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,hashlib,importlib.metadata,json,os,time
import numpy as np
from flint import arb,arb_mat,ctx
import algebraic_root_feasibility as local

ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'runs/cohort_uniform_root_review_v1'
PRODUCER=ROOT/'runs/cohort_uniform_roots_v1'
OUT=ROOT/'evidence/cohort_uniform_roots_independent_review.json'
NAMES=('apical_area_modulus','basal_area_modulus','apical_surface_tension','basal_surface_tension','common_gain')
OLD_M_REVIEW_SHA='4b03db70583be266a8d4ae581f19c8db0099bc028ddd6af82a5c990a91d8e600'
EXPECTED_CELL='64f93fefa034f4c004b93451d98b2fd04e0072855be8c5e3bba9babb8f394dbf'
EXPECTED_CONTRACT='2498493bf50049295eb2ad0d2360b1defaf387ee8db8828ceb576f3d7e942f82'
QUANTIFIER='for_all_parameters_in_closed_cell_exists_unique_root_in_declared_branch_domain'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf8'))
def relative(p):return Path(p).resolve().relative_to(ROOT.resolve()).as_posix()
def path(s):
    if not isinstance(s,str):raise ValueError('string relative path required')
    q=(ROOT/s).resolve()
    if not q.is_relative_to(ROOT.resolve()):raise ValueError('outside study')
    return q

def write(p,v):
    p=Path(p)
    if p.exists():raise RuntimeError('immutable review file already exists: '+str(p))
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf8')

def cell_hash(cell):
    canonical={k:[str(Q(a)),str(Q(b))] for k,(a,b) in sorted(cell.items())}
    return hashlib.sha256(json.dumps(canonical,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def alltrue(checks):return isinstance(checks,dict) and bool(checks) and all(v is True for v in checks.values())
def A(v):return local.A(v)
def up(v):return local.upper(v)
def pack(v):return local.pack(v)
def max_upper(values):return max(v.upper() for v in values)
def infnorm(v):return max(up(v[i,0]) for i in range(v.nrows()))
def rowsum(v):return max(sum((up(v[i,j]) for j in range(v.ncols())),arb(0)).upper() for i in range(v.nrows()))

def strict_acceptance(eta,Z,r):
    if not isinstance(r,arb) or not bool(r>0):raise ValueError('positive exact radius required')
    return bool(Z<1 and eta>=0 and Z>=0 and eta+Z*r<r)

def guard_tests():
    ctx.prec=192;local.A.cache_clear();checks={}
    checks['valid_strict_test']=strict_acceptance(A(Q(1,10)),A(Q(1,5)),A(1))
    checks['selfmap_equality_rejected']=not strict_acceptance(A(Q(4,5)),A(Q(1,5)),A(1))
    checks['Z_equality_rejected']=not strict_acceptance(A(0),A(1),A(1))
    checks['negative_eta_rejected']=not strict_acceptance(A(-1),A(Q(1,5)),A(1))
    checks['negative_Z_rejected']=not strict_acceptance(A(0),A(-1),A(1))
    try:strict_acceptance(A(0),A(0),A(0));checks['zero_radius_rejected']=False
    except ValueError:checks['zero_radius_rejected']=True
    checks['overlapping_interval_max_uses_all_upper_endpoints']=max_upper([arb('1 +/- 0.25'),arb('1.1 +/- 0.3')])>=arb('1.4')
    checks['overlapping_interval_max_uses_all_upper_endpoints']=bool(checks['overlapping_interval_max_uses_all_upper_endpoints'])
    checks['nonboolean_checks_rejected']=not alltrue({'x':1}) and not alltrue({})
    checks['dose_fraction_normalization']=Q('0.78208913')==Q('78208913/100000000')
    checks['exact_cell_hash_fraction_normalization']=cell_hash({'a':['0.5','1']})==cell_hash({'a':['1/2','1']})
    try:path('../outside.json');checks['outside_path_rejected']=False
    except ValueError:checks['outside_path_rejected']=True
    # Exact parameter derivatives at a signed nondegenerate surface point.
    import sympy as sp
    RR,xx,zz,KK,TT,ss=sp.symbols('RR xx zz KK TT ss',real=True)
    js=ss*RR*sp.sqrt(xx*xx+zz*zz);potential=KK*(js-1)**2/2+TT*js
    for i,q in enumerate((RR,xx,zz)):
        checks['surface_K_column_'+str(i)]=sp.simplify(sp.diff(potential,q,KK)-(js-1)*sp.diff(js,q))==0
        checks['surface_T_column_'+str(i)]=sp.simplify(sp.diff(potential,q,TT)-sp.diff(js,q))==0
    aa,bb,cc,dd,hh,gg,II,mask=sp.symbols('aa bb cc dd hh gg II mask',real=True)
    eactive=gg*II*mask*(bb*bb+dd*dd)/2
    for i,q in enumerate((aa,bb,cc,dd,hh)):
        expected=II*mask*q if i in (1,3) else 0
        checks['bulk_gain_column_'+str(i)]=sp.simplify(sp.diff(eactive,q,gg)-expected)==0
    # The row majorant is valid even before the approximate inverse is applied.
    BB=sp.Matrix([[1,-2],[3,1]]);DH=sp.Matrix([[1,-1],[2,0]]);RRm=sp.Matrix([[1,1],[0,1]])
    exact=RRm*BB.T*DH*BB
    bm=BB.applyfunc(abs);localmajor=DH.applyfunc(abs)*sp.Matrix([sum(bm[i,j] for j in range(2)) for i in range(2)])
    bound=RRm.applyfunc(abs)*bm.T*localmajor
    checks['projection_absolute_row_majorant_exact_example']=all(sum(abs(exact[i,j]) for j in range(2))<=bound[i] for i in range(2))
    return checks

def bindings(protocol):
    checks={'review_source_unchanged':sha(__file__)==protocol['review_source_sha256'],
            'protocol_producer_unchanged':sha(PRODUCER/'protocol.json')==protocol['producer_protocol_sha256'],
            'all_frozen_sources_inputs_unchanged':all(sha(path(k))==v for k,v in protocol['file_hashes'].items()),
            'python_flint_version_unchanged':importlib.metadata.version('python-flint')==protocol['python_flint_version']}
    if not alltrue(checks):raise RuntimeError('frozen binding failure '+str([k for k,v in checks.items() if not v]))
    return checks

def prepare(manifest_path):
    """Freeze a final per-record provider manifest AFTER producer terminal state."""
    if (RUN/'protocol.json').exists():raise RuntimeError('preserve existing review protocol')
    manifest=read(manifest_path);producer=read(PRODUCER/'protocol.json');terminal_file=path(manifest.get('terminal_job_file',relative(PRODUCER/'job.json')));job=read(terminal_file)
    if job.get('status')!='FINISHED':raise RuntimeError('producer must be truly terminal before final review manifest')
    if manifest.get('producer_protocol_sha256')!=sha(PRODUCER/'protocol.json'):raise ValueError('manifest producer mismatch')
    rows=manifest.get('records')
    if not isinstance(rows,list) or len(rows)!=17:raise ValueError('exactly17 manifest rows required')
    checks={'all17_canonical_records': [v['record_id'] for v in rows]==[v['record_id'] for v in producer['records']],
            'same_expected_cell':producer['shared_parameter_cell_sha256']==EXPECTED_CELL==cell_hash(producer['shared_parameter_cell']),
            'same_expected_common_contract':producer['common_model_observation_contract_sha256']==EXPECTED_CONTRACT,
            'all_guard_tests':alltrue(guard_tests()),'producer_terminal_FINISHED':job['status']=='FINISHED'}
    files={relative(manifest_path):sha(manifest_path),relative(PRODUCER/'protocol.json'):sha(PRODUCER/'protocol.json'),relative(terminal_file):sha(terminal_file)}
    for k,v in producer['input_hashes'].items():
        if sha(path(k))!=v:raise RuntimeError('producer frozen input changed')
        files[k]=v
    for k,v in producer['source_hashes'].items():
        q=ROOT/'scripts'/k
        if sha(q)!=v:raise RuntimeError('producer frozen source changed')
        files[relative(q)]=v
    for k,v in producer['common_contract_source_hashes'].items():
        q=ROOT/'scripts'/k
        if sha(q)!=v:raise RuntimeError('common contract source changed')
        files[relative(q)]=v
    for candidate,original in zip(rows,producer['records']):
        q=path(candidate['certificate_file']);result=read(q)
        if result['record_id']!=original['record_id'] or Q(result['dose_exact'])!=Q(original['dose_exact']):raise ValueError('record or dose mismatch')
        if result['parameter_cell_sha256']!=EXPECTED_CELL or cell_hash(result['parameter_cell'])!=EXPECTED_CELL:raise ValueError('shared cell mismatch')
        if result['model_observation_contract_sha256']!=EXPECTED_CONTRACT:raise ValueError('observation contract mismatch')
        provider_path=path(candidate['provider_protocol_file']);provider=read(provider_path)
        if result['protocol_sha256']!=sha(provider_path):raise ValueError('selected provider protocol mismatch')
        files[relative(provider_path)]=sha(provider_path)
        if not isinstance(provider.get('source_hashes'),dict) or not provider['source_hashes']:raise ValueError('provider nonempty frozen source closure required')
        for source_name,source_hash in provider['source_hashes'].items():
            source_path=path(source_name) if '/' in source_name or '\\' in source_name else ROOT/'scripts'/source_name
            if sha(source_path)!=source_hash:raise ValueError('selected provider frozen source changed')
            files[relative(source_path)]=source_hash
        for input_name,input_hash in provider.get('input_hashes',{}).items():
            input_path=path(input_name)
            if sha(input_path)!=input_hash:raise ValueError('selected provider input changed')
            files[relative(input_path)]=input_hash
        if result['saved_state_sha256']!=original['saved_state_sha256'] or result['saved_state']!=original['saved_state']:raise ValueError('source state mismatch')
        if result['branch_law_semantics']!='independent_stationary_root_per_record' or result['branch_domain_id']!=original['branch_domain_id']:raise ValueError('branch mismatch')
        if result['quantifier']!=QUANTIFIER:raise ValueError('quantifier mismatch')
        if sha(q)!=candidate['certificate_sha256']:raise ValueError('manifest result hash mismatch')
        files[relative(q)]=sha(q)
        for k,v in result.get('artifact_hashes',{}).items():
            if sha(path(k))!=v:raise ValueError('per-record cache/input changed')
            files[k]=v
        for k in ('point_cache_file','predictor_file'):
            pp=path(result[k]);files[relative(pp)]=sha(pp)
        cache=read(path(result['point_cache_file']));mp=path(cache['matrix_archive']);files[relative(mp)]=sha(mp)
        if sha(mp)!=cache['matrix_archive_sha256']:raise ValueError('matrix archive hash mismatch')
        # A supplement is permitted only through its own frozen prelaunch receipt.
        for dependency in candidate.get('additional_frozen_files',[]):
            dp=path(dependency['file'])
            if sha(dp)!=dependency['sha256']:raise ValueError('supplement binding mismatch')
            files[relative(dp)]=sha(dp)
        reviewpath=path(candidate['provider_prelaunch_review']);review=read(reviewpath)
        if review.get('status')!='PASS' or not alltrue(review.get('checks')) or review.get('protocol_sha256')!=sha(provider_path) or review.get('source_sha256') not in provider['source_hashes'].values():raise ValueError('bound provider prelaunch review required')
        files[relative(reviewpath)]=sha(reviewpath)
    for dependency in manifest.get('historical_files',[]):
        hp=path(dependency['file'])
        if sha(hp)!=dependency['sha256']:raise ValueError('historical failure/retry binding mismatch')
        files[relative(hp)]=sha(hp)
    checks['all17_source_dose_cell_branch_and_artifact_bindings']=True
    if not alltrue(checks):raise RuntimeError(checks)
    protocol={'status':'FROZEN_INDEPENDENT_REPLAY_PROTOCOL','precision_bits':192,'review_source_sha256':sha(__file__),
              'producer_protocol_sha256':sha(PRODUCER/'protocol.json'),'file_hashes':files,'python_flint_version':importlib.metadata.version('python-flint'),
              'manifest':relative(manifest_path),'records':rows,'shared_parameter_cell_sha256':EXPECTED_CELL,'model_observation_contract_sha256':EXPECTED_CONTRACT,
              'method':'Fresh analytic global assembly from reviewed frozen-cache local algebra, independent columns and Taylor/row bounds; no producer point/cell routines.',
              'acceptance':'Strict independent Z<1 and eta+Z*r<r at EXACT producer FULL corrector radius. No arbitrary crossprecision majorant equality.',
              'actual_M_reuse_review_sha256':OLD_M_REVIEW_SHA,'no_equilibrium_optimization':True,'read_large_outcomes':False,
              'scope':'Uniform finite-quadrature stationary branches only; no observer/global stability/PDE/full identified set claim.'}
    write(RUN/'protocol.json',protocol);write(RUN/'preflight.json',{'status':'PASS','checks':checks,'protocol_sha256':sha(RUN/'protocol.json'),'review_source_sha256':sha(__file__),'new_hessian_calculations':0})
    print(json.dumps({'status':'PREFLIGHT_PASS_ONLY','protocol_sha256':sha(RUN/'protocol.json'),'preflight_sha256':sha(RUN/'preflight.json'),'checks':len(checks)}),flush=True)


def selected_cell_equal(result):
    return cell_hash(result['uniform_cell_result']['cell_bounds_exact'])==cell_hash(result['parameter_cell'])

def review_record(original,result,precision):
    started=time.time();ctx.prec=precision;local.A.cache_clear();checks={}
    state=read(path(original['saved_state']));cache=read(path(result['point_cache_file']))
    model,u,R=local.restore_cache(cache);n=len(u);dose=Q(original['dose_exact'])
    basep=read(PRODUCER/'protocol.json')['base_parameters_exact'];base=[Q(basep[k]) for k in NAMES]
    with np.load(path(result['predictor_file']),allow_pickle=False) as zz:S=zz['S'].copy()
    with np.load(path(cache['matrix_archive']),allow_pickle=False) as zz:
        storedS=zz['displacement_tangents'][:,[1,2,3,4,5]].copy();storedfree=zz['free_indices'].copy();storedu=zz['free_vector'].copy()
    storedS/=np.array([float(base[0]),float(base[1]),1.,float(base[3]),float(base[4])])[None,:]
    checks['binary_predictor_reproduces_archived_tangent_rule']=np.array_equal(S,storedS)
    checks['saved_state_center_free_map_reference_reproduced']=np.array_equal(u,storedu) and np.array_equal(model.free,storedfree) and np.array_equal(u,np.asarray(state['free_displacement'])) and np.array_equal(model.reference,np.asarray(state['reference_nodes']))
    producer=read(PRODUCER/'protocol.json');checks['cache_exact_source_state_and_saved_strength_bound']=cache['input_state']==original['saved_state'] and cache['input_state_sha256']==original['saved_state_sha256'] and Q.from_float(float(state['activation_strength']))==Q(original['solver_strength_binary_exact'])
    checks['material_geometry_exact_source_record_bound']=state['parameters']==producer['fixed_full_parameters']==cache['parameters'] and Q(original['base_gain_exact'])==base[4] and Q(original['exact_requested_strength'])==dose*base[4]
    checks['reference_mesh_reproduced']=np.array_equal(model.triangles,np.asarray(state['triangles']))
    checks['producer_declares_certified_uniform_root']=result['status']=='CERTIFIED_UNIFORM_LOCAL_ROOT'
    checks['cell_matches_result_and_expected_hash']=cell_hash(result['parameter_cell'])==EXPECTED_CELL and selected_cell_equal(result)
    checks['finite_square_preconditioner_and_predictor']=R.shape==(n,n) and S.shape==(n,5) and np.isfinite(R).all() and np.isfinite(S).all()
    selected=result['uniform_cell_result'];cell=result['parameter_cell']
    radius=Q(selected['selected_radius_exact']) if 'selected_radius_exact' in selected else Q(selected['radius_attempts'][-1]['radius_exact'])
    if radius<=0:raise ValueError('nonpositive corrector radius')
    RM=arb_mat(n,n,[A(float(v)) for v in R.ravel()]);AR=arb_mat(n,n,[A(abs(float(v))) for v in R.ravel()]);SM=arb_mat(n,5,[A(float(v)) for v in S.ravel()])
    g=arb_mat(n,1);H=arb_mat(n,n);B0=arb_mat(n,5);terms=[];surface_count=bulk_count=0
    basal_count=len(model.surfaces[0][0])*len(model.surfaces[0][1])
    for term in local.local_terms(model,u,base[4]*dose):
        ids=term['indices'];nd=len(ids);k=len(term['q']);BM=arb_mat(k,nd,[v for row in term['B'] for v in row])
        fun=local.bulk_local if term['kind']=='bulk' else local.surface_local
        grad,h,_=fun(term['q'],*term['parameters']);HM=arb_mat(k,k,[v for row in h for v in row]);GM=arb_mat(k,1,grad)
        localH=BM.transpose()*HM*BM;localg=BM.transpose()*GM;C=arb_mat(k,5)
        if term['kind']=='bulk':
            element,iq=divmod(bulk_count,len(model.shape));active=model.active.numpy();mask=A(float(active[element,0 if active.shape[1]==1 else iq]));bulk_count+=1
            C[1,4]=mask*A(dose)*term['q'][1];C[3,4]=mask*A(dose)*term['q'][3];term['mask']=mask
        else:
            basal=surface_count<basal_count;surface_count+=1;term['basal']=basal
            rr,x,z=term['q'];scale=term['parameters'][2];length=(x*x+z*z).sqrt();js=scale*rr*length
            djs=[scale*length,scale*rr*x/length,scale*rr*z/length];ki,ti=(1,3) if basal else (0,2)
            for a in range(k):C[a,ki]=(js-1)*djs[a];C[a,ti]=djs[a]
        localB=BM.transpose()*C
        for a,i in enumerate(ids):
            g[i,0]+=term['weight']*localg[a,0]
            for b,j in enumerate(ids):H[i,j]+=term['weight']*localH[a,b]
            for j in range(5):B0[i,j]+=term['weight']*localB[a,j]
        term['HC']=h;term['P']=BM*arb_mat(nd,5,[A(float(S[i,j])) for i in ids for j in range(5)])
        term['absB']=[[up(v) for v in row] for row in term['B']];term['rowsB']=[sum(row,arb(0)).upper() for row in term['absB']]
        if not bool(term['weight']>0):raise ValueError('positive normalized quadrature weights required')
        terms.append(term)
        if len(terms)%3072==0:print(json.dumps({'record':original['record_id'],'phase':'independent_center','terms':len(terms),'seconds':time.time()-started}),flush=True)
    theta=[A(Q(cell[k][0])).union(A(Q(cell[k][1]))) for k in NAMES]
    delta=[theta[i]-A(base[i]) for i in range(5)];hull=[v.union(arb(0)) for v in delta];r=A(radius)
    remrows=[arb(0) for _ in range(n)];drows=[arb(0) for _ in range(n)];minJ=None
    for term in terms:
        k=len(term['q']);offset=[sum((term['P'][a,j]*hull[j] for j in range(5)),arb(0)) for a in range(k)]
        for remainder,rows in ((True,remrows),(False,drows)):
            qb=[term['q'][a]+offset[a]+(arb(0,(r*term['rowsB'][a]).upper()) if not remainder else 0) for a in range(k)]
            if term['kind']=='bulk':
                _,hb,J=local.bulk_local(qb,term['parameters'][0],term['parameters'][1],term['mask']*A(dose)*theta[4]);minJ=J.lower() if minJ is None else min(minJ,J.lower())
            else:
                ki,ti=(1,3) if term['basal'] else (0,2);_,hb,_=local.surface_local(qb,theta[ki],theta[ti],term['parameters'][2])
            width=[up(v) for v in offset] if remainder else term['rowsB']
            change=[[up(hb[a][b]-term['HC'][a][b]) for b in range(k)] for a in range(k)]
            localchange=[sum((change[a][b]*width[b] for b in range(k)),arb(0)) for a in range(k)]
            for j,i in enumerate(term['indices']):rows[i]+=term['weight']*sum((term['absB'][a][j]*localchange[a] for a in range(k)),arb(0))
    identity=arb_mat(n,n,[1 if i==j else 0 for i in range(n) for j in range(n)]);Z0=rowsum(identity-RM*H)
    linear=RM*g+RM*(H*SM+B0)*arb_mat(5,1,delta);remainder=AR*arb_mat(n,1,[x.upper() for x in remrows])
    eta=max_upper(up(linear[i,0])+up(remainder[i,0]) for i in range(n));variation=AR*arb_mat(n,1,[x.upper() for x in drows]);Z=Z0+infnorm(variation);lhs=eta+Z*r
    checks.update({'all9504terms_611free_dofs':len(terms)==9504 and n==611,'all9216bulk_288surface_terms':bulk_count==9216 and surface_count==288,
                   'preconditioner_nonsingular_from_center_defect':bool(Z0<1),'uniform_contraction':bool(Z<1),'uniform_strict_selfmap':bool(lhs<r),
                   'all_analytic_domain_J_gt_1e6_inverse':bool(minJ>A(Q(1,1000000))),'strict_Banach_acceptance':strict_acceptance(eta,Z,r),
                   'no_producer_point_cell_or_optimization_called':True})
    return {'status':'PASS' if alltrue(checks) else 'UNRESOLVED','checks':checks,'precision_bits':precision,'record_id':original['record_id'],
            'free_dofs':n,'local_terms':len(terms),'certificate_sha256':None,'parameter_cell_sha256':cell_hash(cell),
            'model_observation_contract_sha256':EXPECTED_CONTRACT,'uniform_quantifier_verified':alltrue(checks),
            'branch_law_semantics':'independent_stationary_root_per_record','branch_scope':result['branch_scope'],'branch_domain_id':result['branch_domain_id'],
            'bounds':{'Z0':pack(Z0),'eta':pack(eta),'Z':pack(Z),'radius':pack(r),'eta_plus_Zr':pack(lhs),'quadrature_J_lower':pack(minJ)},
            'scope':'Finite-quadrature local stationary-root family only; no observer/PDE/energetic stability/common path assertion.','elapsed_seconds':time.time()-started}


def reuse_M(original,result):
    reviewpath=ROOT/'evidence/parametric_actual_M_independent_review.json';old=read(reviewpath);prior=read(ROOT/'runs/parametric_root_cells_v1/actual_M/root_certificate.json')
    oldprotocol=read(ROOT/'runs/parametric_root_cells_v1/actual_M/protocol.json')
    checks={'immutable_actual_M_review_exact_hash':sha(reviewpath)==OLD_M_REVIEW_SHA,'immutable_actual_M_independent_root_audit_PASS':old['root_audit']['status']=='PASS' and alltrue(old['root_audit']['checks']),
            'uniform_result_identical_to_reviewed_M':result['uniform_cell_result']==prior['uniform_cell_result'],
            'identical_immutable_M_cache_and_predictor':sha(path(result['point_cache_file']))==oldprotocol['point_certificate_sha256'] and sha(path(result['predictor_file']))==oldprotocol['predictor_sha256'],
            'immutable_M_source_cache_review_binding':sha(ROOT/'runs/parametric_root_cells_v1/actual_M/root_certificate.json')==old['root_audit']['certificate_sha256'] and sha(ROOT/'runs/parametric_root_cells_v1/actual_M/protocol.json')==old['reviewed_hashes']['runs/parametric_root_cells_v1/actual_M/protocol.json'],
            'same_mechanical_source_dose_cell':Q(result['dose_exact'])==Q(prior['dose_exact']) and cell_hash(result['parameter_cell'])==cell_hash(prior['parameter_cell']) and result['saved_state_sha256']==sha(path(original['saved_state'])),
            'cohort_source_branch_wrapper_explicit':result['branch_law_semantics']=='independent_stationary_root_per_record' and result['quantifier']==QUANTIFIER and result['model_observation_contract_sha256']==EXPECTED_CONTRACT}
    return {'status':'PASS' if alltrue(checks) else 'UNRESOLVED','checks':checks,'record_id':original['record_id'],'certificate_sha256':None,
            'parameter_cell_sha256':EXPECTED_CELL,'model_observation_contract_sha256':EXPECTED_CONTRACT,'uniform_quantifier_verified':alltrue(checks),
            'branch_law_semantics':result['branch_law_semantics'],'branch_scope':result['branch_scope'],'branch_domain_id':result['branch_domain_id'],
            'reuse_review_file':relative(reviewpath),'reuse_review_sha256':sha(reviewpath),'no_fresh_M_assembly':True,
            'scope':'Reuses immutable independently reviewed identical mechanical root family; new common observation wrapper is not itself an observer proof.'}


def run():
    p=read(RUN/'protocol.json');pre=read(RUN/'preflight.json');approval=read(RUN/'prelaunch_review.json')
    if approval.get('status')!='PASS' or not alltrue(approval.get('checks')) or approval.get('review_source_sha256')!=sha(__file__) or approval.get('protocol_sha256')!=sha(RUN/'protocol.json') or approval.get('preflight_sha256')!=sha(RUN/'preflight.json'):raise RuntimeError('parent prelaunch review required')
    if OUT.exists():raise RuntimeError('preserve completed cohort review')
    bindings(p);ctx.prec=p['precision_bits'];producer=read(PRODUCER/'protocol.json');started=time.time();capsules=[]
    for selected,original in zip(p['records'],producer['records']):
        bindings(p);resultpath=path(selected['certificate_file']);result=read(resultpath);destination=RUN/'records'/(original['record_id']+'.json')
        if destination.exists():
            review=read(destination)
            if review['review_protocol_sha256']!=sha(RUN/'protocol.json') or review['certificate_sha256']!=sha(resultpath):raise RuntimeError('stale review-prefix artifact')
        else:
            review=reuse_M(original,result) if original['mode']=='reuse_immutable_actual_M' else review_record(original,result,p['precision_bits'])
            review['certificate_sha256']=sha(resultpath);review['review_protocol_sha256']=sha(RUN/'protocol.json');review['all_frozen_inputs_unchanged']=alltrue(bindings(p));write(destination,review)
        capsules.append({'record_id':original['record_id'],'review_file':relative(destination),'review_sha256':sha(destination),'audit':review})
        progress={'status':'RUNNING_INDEPENDENT_REPLAY','completed':len(capsules),'records':17,'PASS':sum(z['audit']['status']=='PASS' for z in capsules),'latest_record':original['record_id']}
        (RUN/'progress.json').write_text(json.dumps(progress,indent=2)+'\n',encoding='utf8');print(json.dumps(progress),flush=True)
    checks={'all17_records_reviewed':len(capsules)==17,'one_immutable_M_review_reused':sum(z['audit'].get('no_fresh_M_assembly',False) for z in capsules)==1,
            'all17_independent_acceptance_PASS':all(z['audit']['status']=='PASS' and alltrue(z['audit']['checks']) for z in capsules),
            'all_frozen_bindings_still_match':alltrue(bindings(p)),'no_new_equilibrium_solve':True,'no_PDE_observer_stability_or_identified_set_claim':True}
    write(OUT,{'status':'PASS' if alltrue(checks) else 'PARTIAL_UNRESOLVED','checks':checks,'record_audits':capsules,'review_source_sha256':sha(__file__),
               'protocol_sha256':sha(RUN/'protocol.json'),'preflight_sha256':sha(RUN/'preflight.json'),'prelaunch_review_sha256':sha(RUN/'prelaunch_review.json'),
               'parameter_cell_sha256':EXPECTED_CELL,'model_observation_contract_sha256':EXPECTED_CONTRACT,'precision_bits':p['precision_bits'],
               'scope':'All17 independent uniform finite-quadrature stationary-root acceptance only. Observers/body geometry/completeness/source matching remain separate.','elapsed_seconds':time.time()-started})
    print(json.dumps({'status':read(OUT)['status'],'receipt_sha256':sha(OUT),'seconds':time.time()-started}),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['guards','prepare','run']);parser.add_argument('--manifest',type=Path);args=parser.parse_args()
    if args.action=='guards':print(json.dumps({'status':'PASS' if alltrue(guard_tests()) else 'FAIL','checks':guard_tests(),'source_sha256':sha(__file__)}))
    elif args.action=='prepare':
        if args.manifest is None:parser.error('--manifest required')
        prepare(args.manifest.resolve())
    else:run()
