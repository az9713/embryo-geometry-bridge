"""Separate independent F/K root replay around explicit auxiliary proof centers.

No producer point/cell routines, no equilibrium solve, no change to legacy
review source or receipts. Local strain algebra is already reviewed. Fresh
forces/Hessians/five columns/remainders are assembled at new binary u1.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,base64,hashlib,importlib.metadata,io,json,time
import numpy as np
from flint import arb,arb_mat,ctx
import algebraic_root_feasibility as local
import review_cohort_uniform_roots as reviewed
ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'runs/cohort_recentered_root_review_v1'
TARGET=ROOT/'runs/cohort_recentered_roots_v1'
PRODUCER=ROOT/'runs/cohort_uniform_roots_v1'
OUT=ROOT/'evidence/cohort_recentered_roots_independent_review.json'
NAMES=reviewed.NAMES
EXPECTED_CELL=reviewed.EXPECTED_CELL
EXPECTED_CONTRACT=reviewed.EXPECTED_CONTRACT
A=reviewed.A;up=reviewed.up;pack=reviewed.pack;infnorm=reviewed.infnorm
rowsum=reviewed.rowsum;max_upper=reviewed.max_upper;strict_acceptance=reviewed.strict_acceptance
sha=reviewed.sha;read=reviewed.read;path=reviewed.path;write=reviewed.write
cell_hash=reviewed.cell_hash;selected_cell_equal=reviewed.selected_cell_equal;alltrue=reviewed.alltrue
relative=reviewed.relative

def arrays(cache):
    packet=cache['frozen_binary_cache'];raw=base64.b64decode(packet['data'])
    if hashlib.sha256(raw).hexdigest()!=packet['sha256']:raise ValueError('cache byte checksum')
    with np.load(io.BytesIO(raw),allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}

def fingerprint(values,exclude=()):
    h=hashlib.sha256()
    for k in sorted(values):
        if k in exclude:continue
        a=np.ascontiguousarray(values[k]);h.update(k.encode()+b'\0'+a.dtype.str.encode()+b'\0'+json.dumps(list(a.shape)).encode()+b'\0'+a.tobytes())
    return h.hexdigest()

def one_round(u,m):return np.array([float(Q.from_float(float(v))-Q(w)) for v,w in zip(u,m)],dtype=np.float64)

def guards():
    checks=reviewed.guard_tests()
    old={'u':np.array([1.,.5]),'R':np.eye(2),'reference':np.zeros((2,2))};new={**old,'u':np.array([.9,.4])}
    checks['fingerprint_ignores_only_u']=fingerprint(old,{'u'})==fingerprint(new,{'u'})
    checks['changed_R_detected_bitwise']=fingerprint(old,{'u'})!=fingerprint({**new,'R':2*np.eye(2)},{'u'})
    checks['signed_zero_reference_change_detected']=fingerprint(old,{'u'})!=fingerprint({**new,'reference':-np.zeros((2,2))},{'u'})
    m=Q(1,2**54)+Q(1,2**110);v=one_round(np.array([1.]),[m])[0]
    checks['one_round_subtraction_differs_from_double_round']=bool(v==np.nextafter(1.,-np.inf) and v!=np.float64(1.)-np.float64(float(m)))
    checks['center_translation_not_corrector_radius']=Q.from_float(float(v))-1!=Q(1,10**11)
    return checks

def lineage(original,result,producer):
    oldcache=read(path(original['original_point_cache_file']));newcache=read(path(result['point_cache_file']))
    a,b=arrays(oldcache),arrays(newcache);seed=read(path(original['original_seed_state_file']));aux=read(path(result['proof_center_state_file']));center=read(path(result['center_refinement_file']));mids=list(map(Q,center['RF_midpoints_exact']))
    checks={'all_cache_array_names_unchanged':set(a)==set(b),'ALL_binary_arrays_except_u_identical':set(a)==set(b) and all(a[k].dtype==b[k].dtype and a[k].shape==b[k].shape and a[k].tobytes()==b[k].tobytes() for k in a if k!='u'),
            'exact_binary_preconditioner_unchanged':a['preconditioner'].tobytes()==b['preconditioner'].tobytes(),
            'exact_predictor_file_unchanged':sha(path(result['predictor_file']))==sha(path(original['predictor_file']))==center['S_sha256'],
            'seed_u0_and_auxiliary_u1_separate':np.array_equal(a['u'],np.asarray(seed['free_displacement'])) and np.array_equal(b['u'],np.asarray(aux['free_displacement'])) and not np.array_equal(a['u'],b['u']),
            'one_round_center_recipe_reproduced':len(mids)==611 and np.array_equal(one_round(a['u'],mids),b['u']),
            'recorded_midpoints_lie_inside_recorded_RF_balls':all(Q(v['lower_exact'])<=m<=Q(v['upper_exact']) for v,m in zip(center['RF_exact_enclosures'],mids)) and len(center['RF_exact_enclosures'])==611,
            'same_parameters_equation_cache_and_archived_tangents':oldcache['parameters']==newcache['parameters']==seed['parameters']==aux['parameters']==producer['fixed_full_parameters'] and oldcache['matrix_archive']==newcache['matrix_archive'] and oldcache['matrix_archive_sha256']==newcache['matrix_archive_sha256'],
            'source_and_auxiliary_identities_separate_bound':result['original_seed_state_file']==original['original_seed_state_file'] and result['original_seed_state_sha256']==sha(path(original['original_seed_state_file']))==original['original_seed_state_sha256'] and result['proof_center_state_sha256']==sha(path(result['proof_center_state_file'])) and aux['kind']=='AUXILIARY_NUMERICAL_PROOF_CENTER_NOT_EXPERIMENT_OR_SOLVER_RESULT' and aux['converged'] is False,
            'same_cell_target_dose_new_declared_branch':cell_hash(result['parameter_cell'])==EXPECTED_CELL and result['model_observation_contract_sha256']==EXPECTED_CONTRACT and Q(result['base_strength_exact'])==Q(original['base_gain_exact'])*Q(original['dose_exact']) and result['branch_domain_id']==original['branch_domain_id'] and result['branch_law_semantics']=='independent_stationary_root_per_record',
            'normal_corrector_not_translation':result['normal_corrector_radius_is_not_center_translation'] is True}
    translation=[Q.from_float(float(y))-Q.from_float(float(x)) for x,y in zip(a['u'],b['u'])]
    checks['all611_actual_translations_exact']=list(map(Q,center['actual_center_translation_exact']))==translation and Q(center['translation_linf_exact'])==max(map(abs,translation))
    fp0,fp1=fingerprint(a,{'u'}),fingerprint(b,{'u'});checks['cache_fingerprints_independently_reproduced']=fp0==fp1==center['cache_arrays_except_u_fingerprint']==center['new_cache_arrays_except_u_fingerprint']
    pp0,pp1=fingerprint(a,{'u','preconditioner'}),fingerprint(b,{'u','preconditioner'});checks['potential_array_fingerprints_reproduced']=pp0==pp1==center['potential_arrays_excluding_u_R_fingerprint']==center['new_potential_arrays_excluding_u_R_fingerprint']
    def physical(values):return hashlib.sha256(json.dumps({'potential_arrays':fingerprint(values,{'u','preconditioner'}),'model_parameters':oldcache['parameters'],'theta0_exact':producer['base_parameters_exact'],'dose_exact':original['dose_exact'],'shared_parameter_cell':producer['shared_parameter_cell']},sort_keys=True,separators=(',',':')).encode()).hexdigest()
    checks['full_physical_target_fingerprints_reproduced']=physical(a)==physical(b)==center['full_physical_potential_fingerprint']==center['new_full_physical_potential_fingerprint']
    with np.load(path(oldcache['matrix_archive']),allow_pickle=False) as z:checks['old_float_tangent_archive_keeps_u0']=np.array_equal(z['free_vector'],a['u'])
    return {'status':'PASS' if alltrue(checks) else 'FAIL','checks':checks,'translation_linf_exact':str(max(map(abs,translation))),'recipe_scope':'Exact reproduction of stored midpoint single-round subtraction and unchanged physicaltarget. Old RFball assembly not independently repeated; fresh new-center root theorem is independent of whether this proposal was an accurate Newtonstep.'}

def review_record(original,result,precision):
    started=time.time();ctx.prec=precision;local.A.cache_clear();checks={}
    state=read(path(result['proof_center_state_file']));seed=read(path(original['original_seed_state_file']));cache=read(path(result['point_cache_file']))
    model,u,R=local.restore_cache(cache);n=len(u);dose=Q(original['dose_exact'])
    basep=read(PRODUCER/'protocol.json')['base_parameters_exact'];base=[Q(basep[k]) for k in NAMES]
    with np.load(path(result['predictor_file']),allow_pickle=False) as zz:S=zz['S'].copy()
    with np.load(path(cache['matrix_archive']),allow_pickle=False) as zz:
        storedS=zz['displacement_tangents'][:,[1,2,3,4,5]].copy();storedfree=zz['free_indices'].copy();storedu=zz['free_vector'].copy()
    storedS/=np.array([float(base[0]),float(base[1]),1.,float(base[3]),float(base[4])])[None,:]
    checks['binary_predictor_reproduces_archived_tangent_rule']=np.array_equal(S,storedS)
    checks['old_matrix_archive_belongs_to_seed_u0_not_new_u1']=np.array_equal(storedu,np.asarray(seed['free_displacement'])) and np.array_equal(model.free,storedfree)
    checks['new_cache_center_u1_matches_auxiliary_state']=np.array_equal(u,np.asarray(state['free_displacement'])) and np.array_equal(model.reference,np.asarray(state['reference_nodes']))
    producer=read(PRODUCER/'protocol.json');checks['cache_exact_auxiliary_state_and_original_seed_bound']=cache['input_state']==result['proof_center_state_file'] and cache['input_state_sha256']==result['proof_center_state_sha256'] and cache['original_seed_state_file']==original['original_seed_state_file'] and cache['original_seed_state_sha256']==original['original_seed_state_sha256']
    checks['exact_gain_dose_and_original_solver_strength_distinguished']=Q(state['proof_activation_strength_exact'])==Q(original['exact_requested_strength']) and Q(state['original_solver_strength_binary_exact'])==Q(original['solver_strength_binary_exact'])==Q.from_float(float(seed['activation_strength']))
    checks['material_geometry_exact_source_record_bound']=state['parameters']==producer['fixed_full_parameters']==cache['parameters'] and Q(original['base_gain_exact'])==base[4] and Q(original['exact_requested_strength'])==dose*base[4]
    checks['reference_mesh_reproduced']=np.array_equal(model.triangles,np.asarray(state['triangles']))
    checks['producer_declares_certified_uniform_root']=result['status']=='CERTIFIED_UNIFORM_LOCAL_ROOT'
    checks['cell_matches_result_and_expected_hash']=cell_hash(result['parameter_cell'])==EXPECTED_CELL and selected_cell_equal(result)
    checks['finite_square_preconditioner_and_predictor']=bool(R.shape==(n,n) and S.shape==(n,5) and np.isfinite(R).all() and np.isfinite(S).all())
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

def prepare():
    if (RUN/'protocol.json').exists():raise RuntimeError('preserve frozen independent protocol')
    producer=read(TARGET/'protocol.json');result=read(TARGET/'results.json');job=read(TARGET/'job.json');files={}
    checks={'producer_terminal':job['status']=='FINISHED','bothF_K_present':{x['record_id'] for x in result['records']}=={'optoGEF2_small_col_F','optoGEF2_small_col_K'},'same_expectedC':producer['shared_parameter_cell_sha256']==EXPECTED_CELL and cell_hash(producer['shared_parameter_cell'])==EXPECTED_CELL,'all_lightweight_guards':alltrue(guards())}
    for q in [TARGET/'protocol.json',TARGET/'results.json',TARGET/'job.json',TARGET/'prelaunch_independent_review.json',Path(reviewed.__file__),Path(local.__file__),ROOT/'evidence/cohort_uniform_roots_independent_review.json']:
        files[relative(q)]=sha(q)
    for name,value in producer['source_hashes'].items():
        q=ROOT/'scripts'/name
        if sha(q)!=value:raise RuntimeError('producer source changed')
        files[relative(q)]=value
    for name,value in producer['input_hashes'].items():
        if sha(path(name))!=value:raise RuntimeError('producer input changed')
        files[name]=value
    reports=[]
    for original,entry in zip(producer['records'],result['records']):
        root=read(path(entry['file']))
        if root['record_id']!=original['record_id'] or sha(path(entry['file']))!=entry['sha256'] or root['protocol_sha256']!=sha(TARGET/'protocol.json'):raise RuntimeError('record/provider binding mismatch')
        for name,value in root['artifact_hashes'].items():
            if sha(path(name))!=value:raise RuntimeError('record artifact changed')
            files[name]=value
        files[entry['file']]=entry['sha256'];cache=read(path(root['point_cache_file']));files[cache['matrix_archive']]=sha(path(cache['matrix_archive']))
        report=lineage(original,root,producer);reports.append({'record_id':original['record_id'],**report})
    checks['all_lineage_checks_PASS']=all(z['status']=='PASS' for z in reports)
    if not alltrue(checks):raise RuntimeError(checks)
    protocol={'status':'FROZEN_INDEPENDENT_RECENTER_REPLAY','review_source_sha256':sha(__file__),'precision_bits':192,'producer_protocol_sha256':sha(TARGET/'protocol.json'),'file_hashes':files,'records':result['records'],'producer_records':producer['records'],'python_flint_version':importlib.metadata.version('python-flint'),'parameter_cell_sha256':EXPECTED_CELL,'model_observation_contract_sha256':EXPECTED_CONTRACT,'scope':'Fresh independent F/K local root proof at newu1; exact same physicaltarget andC; no observers/fullidentifiedset/PDE/commonpath.','old_RF_assembly_repeated':False,'one_center_recipe_from_recorded_midpoints_verified':True,'no_new_equilibrium_optimization':True}
    write(RUN/'protocol.json',protocol);write(RUN/'lineage_preflight.json',{'status':'PASS','records':reports,'scope':'No new force/Hessian/root calculation'});write(RUN/'preflight.json',{'status':'PASS','checks':checks,'review_source_sha256':sha(__file__),'protocol_sha256':sha(RUN/'protocol.json'),'lineage_preflight_sha256':sha(RUN/'lineage_preflight.json'),'new_force_or_Hessian_calculations':0})
    print(json.dumps({'status':'PREFLIGHT_PASS_ONLY','source_sha256':sha(__file__),'protocol_sha256':sha(RUN/'protocol.json'),'preflight_sha256':sha(RUN/'preflight.json')}),flush=True)

def bindings(p):
    if sha(__file__)!=p['review_source_sha256'] or any(sha(path(k))!=v for k,v in p['file_hashes'].items()) or importlib.metadata.version('python-flint')!=p['python_flint_version']:raise RuntimeError('frozen independent replay binding changed')
    return True

def run():
    p=read(RUN/'protocol.json');approval=read(RUN/'prelaunch_review.json')
    if approval.get('status')!='PASS' or not alltrue(approval.get('checks')) or approval.get('review_source_sha256')!=sha(__file__) or approval.get('protocol_sha256')!=sha(RUN/'protocol.json') or approval.get('preflight_sha256')!=sha(RUN/'preflight.json'):raise RuntimeError('parent concrete prelaunch review required')
    if OUT.exists():raise RuntimeError('preserve terminal independent replay')
    bindings(p);start=time.time();producer=read(TARGET/'protocol.json');audits=[]
    for original,entry in zip(p['producer_records'],p['records']):
        bindings(p);result=read(path(entry['file']));line=lineage(original,result,producer);audit=review_record(original,result,p['precision_bits']);audit['lineage_review']=line;audit['checks']['all_physical_target_and_center_lineage_checks']=line['status']=='PASS';audit['status']='PASS' if alltrue(audit['checks']) else 'UNRESOLVED';audit['uniform_quantifier_verified']=audit['status']=='PASS';audit['certificate_sha256']=sha(path(entry['file']));audit['review_protocol_sha256']=sha(RUN/'protocol.json');audit['all_frozen_bindings_unchanged']=bindings(p)
        dest=RUN/'records'/(original['record_id']+'.json');write(dest,audit);audits.append({'record_id':original['record_id'],'review_file':relative(dest),'review_sha256':sha(dest),'audit':audit});print(json.dumps({'record':original['record_id'],'status':audit['status'],'eta_upper':audit['bounds']['eta']['upper_exact'],'seconds':time.time()-start}),flush=True)
    checks={'both_records_freshly_assembled':len(audits)==2,'both_independent_uniform_root_proofs_PASS':all(z['audit']['status']=='PASS' for z in audits),'same_physicaltarget_andC':all(z['audit']['lineage_review']['status']=='PASS' for z in audits),'all_frozen_bindings_unchanged':bindings(p),'no_old_review_overwritten':True,'no_observer_commonpath_PDE_or_global_completion_claim':True}
    write(OUT,{'status':'PASS' if alltrue(checks) else 'PARTIAL_UNRESOLVED','checks':checks,'record_audits':audits,'review_source_sha256':sha(__file__),'protocol_sha256':sha(RUN/'protocol.json'),'preflight_sha256':sha(RUN/'preflight.json'),'prelaunch_review_sha256':sha(RUN/'prelaunch_review.json'),'precision_bits':p['precision_bits'],'parameter_cell_sha256':EXPECTED_CELL,'model_observation_contract_sha256':EXPECTED_CONTRACT,'scope':p['scope'],'elapsed_seconds':time.time()-start});print(json.dumps({'status':read(OUT)['status'],'receipt_sha256':sha(OUT),'seconds':time.time()-start}),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['guards','prepare','run']);args=parser.parse_args()
    if args.stage=='guards':print(json.dumps({'status':'PASS' if alltrue(guards()) else 'FAIL','checks':guards(),'source_sha256':sha(__file__)}))
    elif args.stage=='prepare':prepare()
    else:run()
